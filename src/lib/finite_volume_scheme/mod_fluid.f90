module mod_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic
  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data, TINY_MACH
  use mod_nondimensionalization, only: rho_0, v_0, p_0, e_0
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_floating_point_utils, only: near_zero
  use mod_surrogate, only: surrogate
  use mod_units
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_surrogate, only: surrogate
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_flux_tensor, only: operator(.dot.), H => flux_tensor_t
  use mod_time_integrator_factory, only: time_integrator_factory

  implicit none

  private
  public :: fluid_t, new_fluid

  type, extends(integrand_t) :: fluid_t
    real(rk), dimension(:, :, :), allocatable :: conserved_vars !< ((rho, rho*u, rho*v, rho*E), i, j); Conserved quantities
  contains
    procedure, public :: initialize
    procedure, private :: initialize_from_ini
    procedure, private :: initialize_from_hdf5
    procedure, public :: t => time_derivative
    procedure, public :: get_sound_speed
    procedure, public :: residual_smoother
    procedure, public :: get_max_sound_speed
    procedure, public :: get_primitive_vars
    procedure, public :: sanity_check
    procedure, nopass, private :: flux_edges
    procedure, pass(lhs), public :: type_plus_type => add_fluid
    procedure, pass(lhs), public :: type_minus_type => subtract_fluid
    procedure, pass(lhs), public :: type_mul_real => fluid_mul_real
    procedure, pass(rhs), public :: real_mul_type => real_mul_fluid
    procedure, pass(lhs), public :: assign => assign_fluid
    procedure, public :: force_finalization
    procedure, private, nopass :: add_fields
    final :: finalize
  end type fluid_t

contains

  function new_fluid(input, finite_volume_scheme) result(fluid)
    !< Fluid constructor
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    type(fluid_t), pointer :: fluid

    allocate(fluid)
    call fluid%initialize(input, finite_volume_scheme)
  end function new_fluid

  subroutine initialize(self, input, finite_volume_scheme)
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    class(strategy), pointer :: time_integrator => null()

    integer(ik) :: alloc_status

    alloc_status = 0
    call debug_print('Initializing fluid_t', __FILE__, __LINE__)

    associate(imin=>finite_volume_scheme%grid%ilo_bc_cell, &
              imax=>finite_volume_scheme%grid%ihi_bc_cell, &
              jmin=>finite_volume_scheme%grid%jlo_bc_cell, &
              jmax=>finite_volume_scheme%grid%jhi_bc_cell)

      allocate(self%conserved_vars(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to allocate fluid_t%conserved_vars"
      end if
      self%conserved_vars = 0.0_rk
    end associate

    time_integrator => time_integrator_factory(input)
    allocate(self%time_integrator, source=time_integrator, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%time_integrator"
    deallocate(time_integrator)

    if(input%read_init_cond_from_file .or. input%restart_from_file) then
      call self%initialize_from_hdf5(input, finite_volume_scheme)
    else
      call self%initialize_from_ini(input)
    end if

  end subroutine initialize

  subroutine initialize_from_hdf5(self, input, finite_volume_scheme)
    !< Initialize from an .hdf5 file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the hdf5 file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    type(hdf5_file) :: h5
    logical :: file_exists
    character(:), allocatable :: filename
    character(32) :: str_buff = ''

    real(rk), dimension(:, :), allocatable :: density
    real(rk), dimension(:, :), allocatable :: x_velocity
    real(rk), dimension(:, :), allocatable :: y_velocity
    real(rk), dimension(:, :), allocatable :: pressure

    call debug_print('Initializing fluid_t from hdf5', __FILE__, __LINE__)

    if(input%restart_from_file) then
      filename = trim(input%restart_file)
    else
      filename = trim(input%initial_condition_file)
    end if

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      write(*, '(a)') 'Error in finite_volume_scheme_t%initialize_from_hdf5(); file not found: "'//filename//'"'
      error stop 'Error in finite_volume_scheme_t%initialize_from_hdf5(); file not found, exiting...'
    end if

    call h5%initialize(filename=filename, status='old', action='r')

    call h5%get('/density', density)
    if(input%restart_from_file) then
      call h5%readattr('/density', 'units', str_buff)
      select case(trim(str_buff))
      case('g/cc')
      case default
        error stop "Unknown density units in .h5 file. Acceptable units are 'g/cc'."
      end select
    end if

    call h5%get('/x_velocity', x_velocity)
    call h5%get('/y_velocity', y_velocity)
    if(input%restart_from_file) then
      call h5%readattr('/x_velocity', 'units', str_buff)
      select case(trim(str_buff))
      case('km/s')
        x_velocity = x_velocity * km_per_sec_to_cm_per_sec
        y_velocity = y_velocity * km_per_sec_to_cm_per_sec
      case('cm/s')
      case default
        error stop "Unknown velocity units in .h5 file. Acceptable units are 'km/s' and 'cm/s'."
      end select
    end if

    call h5%get('/pressure', pressure)
    if(input%restart_from_file) then
      call h5%readattr('/pressure', 'units', str_buff)
      select case(trim(str_buff))
      case('barye')
      case('Mbar')
        pressure = pressure * mega_bar_to_barye
      case default
        error stop "Unknown density units in .h5 file. Acceptable units are 'g/cc'."
      end select
    end if

    call h5%finalize()

    ! Non-dimensionalize
    density = density / rho_0
    x_velocity = x_velocity / v_0
    y_velocity = y_velocity / v_0
    pressure = pressure / p_0

    if(any(near_zero(pressure))) then
      error stop "Some (or all) of the pressure array is ~0 in fluid_t%initialize_from_hdf5"
    end if

    if(any(near_zero(density))) then
      error stop "Some (or all) of the density array is ~0 in fluid_t%initialize_from_hdf5"
    end if

    associate(imin=>finite_volume_scheme%grid%ilo_bc_cell, imax=>finite_volume_scheme%grid%ihi_bc_cell, &
              jmin=>finite_volume_scheme%grid%jlo_bc_cell, jmax=>finite_volume_scheme%grid%jhi_bc_cell)

      self%conserved_vars(1, imin:imax, jmin:jmax) = density
      self%conserved_vars(2, imin:imax, jmin:jmax) = density * x_velocity
      self%conserved_vars(3, imin:imax, jmin:jmax) = density * y_velocity
      self%conserved_vars(4, imin:imax, jmin:jmax) = density &
                                                     * eos%calculate_total_energy(pressure=pressure, &
                                                                                  density=density, &
                                                                                  x_velocity=x_velocity, &
                                                                                  y_velocity=y_velocity)
    end associate

    write(*, '(a)') 'Initial fluid stats'
    write(*, '(a)') '==================================================='
    write(*, '(a, f0.4)') 'EOS Gamma:                     ', eos%get_gamma()
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max density    [non-dim]: ', minval(density), maxval(density)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max x-velocity [non-dim]: ', minval(x_velocity), maxval(x_velocity)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max y-velocity [non-dim]: ', minval(x_velocity), maxval(y_velocity)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max pressure   [non-dim]: ', minval(pressure), maxval(pressure)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max density    [dim]:     ', minval(density) * rho_0, maxval(density) * rho_0
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max x-velocity [dim]:     ', minval(x_velocity) * v_0, maxval(x_velocity) * v_0
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max y-velocity [dim]:     ', minval(y_velocity) * v_0, maxval(y_velocity) * v_0
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max pressure   [dim]:     ', minval(pressure) * p_0, maxval(pressure) * p_0
    write(*, '(a)') '==================================================='
    write(*, *)

  end subroutine initialize_from_hdf5

  subroutine initialize_from_ini(self, input)
    !< Initialize from an .ini file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the .ini file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    real(rk) :: density, x_velocity, y_velocity, pressure

    ! Non-dimensionalize
    density = input%init_density / rho_0
    x_velocity = input%init_x_velocity / v_0
    y_velocity = input%init_y_velocity / v_0
    pressure = input%init_pressure / p_0

    call debug_print('Initializing fluid_t from .ini', __FILE__, __LINE__)
    write(*, '(a,4(f0.3, 1x))') 'Initializing with [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure
    self%conserved_vars(1, :, :) = density
    self%conserved_vars(2, :, :) = density * x_velocity
    self%conserved_vars(3, :, :) = density * y_velocity
    self%conserved_vars(4, :, :) = density * eos%calculate_total_energy(pressure=pressure, &
                                                                        density=density, &
                                                                        x_velocity=x_velocity, &
                                                                        y_velocity=y_velocity)

    if(near_zero(input%init_pressure)) then
      error stop "Some (or all) of the pressure array is ~0 in fluid_t%initialize_from_hdf5"
    end if

    if(near_zero(input%init_density)) then
      error stop "Some (or all) of the density array is ~0 in fluid_t%initialize_from_hdf5"
    end if

  end subroutine initialize_from_ini

  subroutine force_finalization(self)
    class(fluid_t), intent(inout) :: self

    call debug_print('Running fluid_t%force_finalization()', __FILE__, __LINE__)
    if(allocated(self%conserved_vars)) deallocate(self%conserved_vars)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine force_finalization

  subroutine finalize(self)
    type(fluid_t), intent(inout) :: self

    call debug_print('Running fluid_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%conserved_vars)) deallocate(self%conserved_vars)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine finalize

  function time_derivative(self, fv) result(d_dt)
    !< Implementation of dU/dt

    class(fluid_t), intent(in) :: self
    class(finite_volume_scheme_t), intent(inout) :: fv

    ! Locals
    class(integrand_t), allocatable :: d_dt !< dU/dt (integrand_t to satisfy parent interface)
    type(fluid_t), allocatable :: local_d_dt !< dU/dt
    integer(ik) :: alloc_status, error_code, i, j

    real(rk), dimension(:, :, :), allocatable :: primitive_vars
    !< ((rho, u, v, p), i, j); Primitive variables at each cell center

    real(rk), dimension(:, :, :), allocatable :: evolved_corner_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each corner

    ! Indexing the midpoints is a pain, so they're split by the up/down edges and left/right edges
    real(rk), dimension(:, :, :), allocatable :: evolved_downup_midpoints_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each midpoint defined by vectors that go up/down
    !< (edges 2 (right edge) and 4 (left edge))

    real(rk), dimension(:, :, :), allocatable :: evolved_leftright_midpoints_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each midpoint defined by vectors that go left/right
    !< (edges 1 (bottom edge) and 3 (top edge))

    real(rk), dimension(:, :, :, :, :), allocatable, target :: reconstructed_state
    !< (((rho, u, v, p)), point, node/midpoint, i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints. Note, this DOES repeat nodes, since corners and midpoints are
    !< shared by neighboring cells, but each point has its own reconstructed value based on the parent cell's state

    integer(ik), dimension(3) :: bounds
    integer(ik), dimension(5) :: recon_bounds

    call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

    associate(imin=>fv%grid%ilo_bc_cell, imax=>fv%grid%ihi_bc_cell, &
              jmin=>fv%grid%jlo_bc_cell, jmax=>fv%grid%jhi_bc_cell)

      allocate(reconstructed_state(4, 4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate reconstructed_state in fluid_t%time_derivative()"
    end associate

    associate(imin_node=>fv%grid%ilo_node, imax_node=>fv%grid%ihi_node, &
              jmin_node=>fv%grid%jlo_node, jmax_node=>fv%grid%jhi_node, &
              imin_cell=>fv%grid%ilo_cell, imax_cell=>fv%grid%ihi_cell, &
              jmin_cell=>fv%grid%jlo_cell, jmax_cell=>fv%grid%jhi_cell)

      ! corners
      allocate(evolved_corner_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_corner_state in fluid_t%time_derivative()"

      ! left/right midpoints
      allocate(evolved_leftright_midpoints_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_leftright_midpoints_state in fluid_t%time_derivative()"

      ! down/up midpoints
      allocate(evolved_downup_midpoints_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_downup_midpoints_state"

    end associate

    allocate(local_d_dt, source=self)
    call local_d_dt%get_primitive_vars(primitive_vars, fv)

    ! Now we can reconstruct the entire domain
    bounds = lbound(primitive_vars)
    call fv%reconstruction_operator%set_primitive_vars_pointer(primitive_vars=primitive_vars, &
                                                               lbounds=bounds)
    call fv%reconstruct(primitive_vars=primitive_vars, &
                        cell_lbounds=bounds, &
                        reconstructed_state=reconstructed_state)

    recon_bounds = lbound(reconstructed_state)
    call fv%evolution_operator%set_reconstructed_state_pointer( &
      reconstructed_state_target=reconstructed_state, &
      lbounds=recon_bounds)

    ! Apply the reconstructed state to the ghost layers
    call fv%apply_reconstructed_state_bc(reconstructed_state, lbounds=recon_bounds)
    call fv%apply_cell_gradient_bc()

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
    bounds = lbound(evolved_downup_midpoints_state)
    call fv%evolution_operator%evolve(evolved_state=evolved_downup_midpoints_state, &
                                      location='down/up midpoint', &
                                      lbounds=bounds, &
                                      error_code=error_code)

    if(error_code /= 0) then
      fv%error_code = error_code
    end if

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
    bounds = lbound(evolved_leftright_midpoints_state)
    call fv%evolution_operator%evolve(evolved_state=evolved_leftright_midpoints_state, &
                                      location='left/right midpoint', &
                                      lbounds=bounds, &
                                      error_code=error_code)

    if(error_code /= 0) then
      fv%error_code = error_code
    end if

    ! Evolve, i.e. E0(R_omega), at all corner nodes
    call debug_print('Evolving corner nodes', __FILE__, __LINE__)
    bounds = lbound(evolved_corner_state)
    call fv%evolution_operator%evolve(evolved_state=evolved_corner_state, &
                                      location='corner', &
                                      lbounds=bounds, &
                                      error_code=error_code)

    if(error_code /= 0) then
      fv%error_code = error_code
    end if

    nullify(fv%reconstruction_operator%primitive_vars)
    nullify(fv%evolution_operator%reconstructed_state)

    call self%flux_edges(grid=fv%grid, &
                         evolved_corner_state=evolved_corner_state, &
                         evolved_leftright_midpoints_state=evolved_leftright_midpoints_state, &
                         evolved_downup_midpoints_state=evolved_downup_midpoints_state, &
                         dU_dt=local_d_dt%conserved_vars)

    call move_alloc(local_d_dt, d_dt)
    call d_dt%set_temp(calling_function='fluid_t%time_derivative (d_dt)', line=__LINE__)

    ! Now deallocate everything
    deallocate(primitive_vars, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate primitive_vars in fluid_t%time_derivative()"

    deallocate(evolved_corner_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate evolved_corner_state in fluid_t%time_derivative()"

    deallocate(evolved_downup_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate evolved_downup_midpoints_state in fluid_t%time_derivative()"

    deallocate(evolved_leftright_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate evolved_leftright_midpoints_state in fluid_t%time_derivative()"

    deallocate(reconstructed_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate reconstructed_state"

  end function time_derivative

  subroutine flux_edges(grid, evolved_corner_state, evolved_leftright_midpoints_state, &
                        evolved_downup_midpoints_state, dU_dt)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    class(grid_t), intent(in) :: grid
    real(rk), dimension(:, grid%ilo_node:, grid%jlo_node:), intent(in) :: evolved_corner_state
    real(rk), dimension(:, grid%ilo_cell:, grid%jlo_node:), intent(in) :: evolved_leftright_midpoints_state
    real(rk), dimension(:, grid%ilo_node:, grid%jlo_cell:), intent(in) :: evolved_downup_midpoints_state
    real(rk), dimension(:, grid%ilo_bc_cell:, grid%jlo_bc_cell:), intent(out) :: dU_dt

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j, edge, xy
    real(rk), dimension(4) :: top_left_corner, top_right_corner, bottom_left_corner, bottom_right_corner
    real(rk), dimension(4) :: bottom_midpoint, right_midpoint, top_midpoint, left_midpoint
    real(rk), dimension(2, 4) :: n_hat
    real(rk), dimension(4) :: delta_l
    real(rk), dimension(4) :: bottom_flux
    real(rk), dimension(4) :: right_flux
    real(rk), dimension(4) :: top_flux
    real(rk), dimension(4) :: left_flux

    ilo = grid%ilo_cell
    ihi = grid%ihi_cell
    jlo = grid%jlo_cell
    jhi = grid%jhi_cell

    bottom_flux = 0.0_rk
    right_flux = 0.0_rk
    top_flux = 0.0_rk
    left_flux = 0.0_rk

    top_left_corner = 0.0_rk
    top_right_corner = 0.0_rk
    bottom_left_corner = 0.0_rk
    bottom_right_corner = 0.0_rk
    bottom_midpoint = 0.0_rk
    right_midpoint = 0.0_rk
    top_midpoint = 0.0_rk
    left_midpoint = 0.0_rk

    call debug_print('Running fluid_t%flux_edges()', __FILE__, __LINE__)

    do j = jlo, jhi
      do i = ilo, ihi

        top_left_corner = evolved_corner_state(:, i, j + 1)
        bottom_left_corner = evolved_corner_state(:, i, j)
        top_right_corner = evolved_corner_state(:, i + 1, j + 1)
        bottom_right_corner = evolved_corner_state(:, i + 1, j)
        bottom_midpoint = evolved_leftright_midpoints_state(:, i, j)
        top_midpoint = evolved_leftright_midpoints_state(:, i, j + 1)
        right_midpoint = evolved_downup_midpoints_state(:, i + 1, j)
        left_midpoint = evolved_downup_midpoints_state(:, i, j)
        delta_l = grid%cell_edge_lengths(:, i, j)

        do edge = 1, 4
          do xy = 1, 2
            n_hat(xy, edge) = grid%cell_edge_norm_vectors(xy, edge, i, j)
          end do
        end do

        ! Edge 1 (bottom)
        bottom_flux = ( &
                      ((H(bottom_left_corner) + &
                        4.0_rk * H(bottom_midpoint) + &
                        H(bottom_right_corner)) .dot.n_hat(:, 1)) * (delta_l(1) / 6.0_rk))

        ! Edge 2 (right)
        right_flux = ( &
                     ((H(bottom_right_corner) + &
                       4.0_rk * H(right_midpoint) + &
                       H(top_right_corner)) .dot.n_hat(:, 2)) * (delta_l(2) / 6.0_rk))

        ! Edge 3 (top)
        top_flux = ( &
                   ((H(top_right_corner) + &
                     4.0_rk * H(top_midpoint) + &
                     H(top_left_corner)) .dot.n_hat(:, 3)) * (delta_l(3) / 6.0_rk))

        ! Edge 4 (left)
        left_flux = ( &
                    ((H(top_left_corner) + &
                      4.0_rk * H(left_midpoint) + &
                      H(bottom_left_corner)) .dot.n_hat(:, 4)) * (delta_l(4) / 6.0_rk))

        dU_dt(:, i, j) = (-1.0_rk / grid%cell_volume(i, j)) * (bottom_flux + right_flux + top_flux + left_flux)
      end do ! i
    end do ! j
  end subroutine flux_edges

  subroutine residual_smoother(self)
    class(fluid_t), intent(inout) :: self

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j
    real(rk) :: eps
    real(rk) :: u, v
    real(rk), dimension(:, :), allocatable :: sound_speed

    ! call self%get_sound_speed(sound_speed)

    ! ilo = lbound(self%conserved_vars, dim=2)
    ! ihi = ubound(self%conserved_vars, dim=2)
    ! jlo = lbound(self%conserved_vars, dim=3)
    ! jhi = ubound(self%conserved_vars, dim=3)

    ! do j = jlo, jhi
    !   do i = ilo, ihi
    !     u = self%conserved_vars(2, i, j) / self%conserved_vars(1, i, j)
    !     v = self%conserved_vars(3, i, j) / self%conserved_vars(1, i, j)

    !     if(abs(u / sound_speed(i, j)) < TINY_MACH) then
    !       self%conserved_vars(2, i, j) = 0.0_rk
    !     end if

    !     if(abs(v / sound_speed(i, j)) < TINY_MACH) then
    !       self%conserved_vars(3, i, j) = 0.0_rk
    !     end if
    !   end do
    ! end do
    ! deallocate(sound_speed)

  end subroutine

  function subtract_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: difference
    type(fluid_t), allocatable :: local_difference

    call debug_print('Running fluid_t%subtract_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(fluid_t)
      allocate(local_difference, source=lhs)
      ! local_difference%conserved_vars = lhs%conserved_vars + rhs%conserved_vars
      call add_fields(a=lhs%conserved_vars, b=rhs%conserved_vars, c=local_difference%conserved_vars) ! c=a-b
    class default
      error stop 'fluid_t%subtract_fluid: unsupported rhs class'
    end select

    call move_alloc(local_difference, difference)
    call difference%set_temp(calling_function='subtract_fluid (difference)', line=__LINE__)
  end function subtract_fluid

  function add_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: sum
    type(fluid_t), allocatable :: local_sum

    call debug_print('Running fluid_t%add_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(fluid_t)
      allocate(local_sum, source=lhs)
      ! local_sum%conserved_vars = lhs%conserved_vars + rhs%conserved_vars
      call add_fields(a=lhs%conserved_vars, b=rhs%conserved_vars, c=local_sum%conserved_vars) ! c=a+b
    class default
      error stop 'fluid_t%add_fluid: unsupported rhs class'
    end select

    call local_sum%set_temp(calling_function='fluid_t%add_fluid(local_sum)', line=__LINE__)
    call move_alloc(local_sum, sum)
    call sum%set_temp(calling_function='fluid_t%add_fluid(sum)', line=__LINE__)
  end function add_fluid

  pure subroutine add_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a + b
    real(rk), dimension(:, :, :), intent(in) :: a
    real(rk), dimension(:, :, :), intent(in) :: b
    real(rk), dimension(:, :, :), intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo, ihi, jlo, jhi, klo, khi

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)
    klo = lbound(a, dim=3)
    khi = ubound(a, dim=3)

    do concurrent(k=klo:khi)
      do concurrent(j=jlo:jhi)
        do i = ilo, ihi
          c(i, j, k) = a(i, j, k) + b(i, j, k)
        end do
      end do
    end do
  end subroutine add_fields

  pure subroutine subtract_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a - b
    real(rk), dimension(:, :, :), intent(in) :: a
    real(rk), dimension(:, :, :), intent(in) :: b
    real(rk), dimension(:, :, :), intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo, ihi, jlo, jhi, klo, khi

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)
    klo = lbound(a, dim=3)
    khi = ubound(a, dim=3)

    do concurrent(k=klo:khi)
      do concurrent(j=jlo:jhi)
        do i = ilo, ihi
          c(i, j, k) = a(i, j, k) - b(i, j, k)
        end do
      end do
    end do
  end subroutine subtract_fields

  pure subroutine mult_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :, :), intent(in) :: a
    real(rk), dimension(:, :, :), intent(in) :: b
    real(rk), dimension(:, :, :), intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo, ihi, jlo, jhi, klo, khi

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)
    klo = lbound(a, dim=3)
    khi = ubound(a, dim=3)

    do concurrent(k=klo:khi)
      do concurrent(j=jlo:jhi)
        do i = ilo, ihi
          c(i, j, k) = a(i, j, k) * b(i, j, k)
        end do
      end do
    end do
  end subroutine mult_fields

  pure subroutine mult_field_by_real(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :, :), intent(in) :: a
    real(rk), intent(in) :: b
    real(rk), dimension(:, :, :), intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo, ihi, jlo, jhi, klo, khi

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)
    klo = lbound(a, dim=3)
    khi = ubound(a, dim=3)

    do concurrent(k=klo:khi)
      do concurrent(j=jlo:jhi)
        do i = ilo, ihi
          c(i, j, k) = a(i, j, k) * b
        end do
      end do
    end do
  end subroutine mult_field_by_real

  function fluid_mul_real(lhs, rhs) result(product)
    !< Implementation of the fluid * real operation
    class(fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(integrand_t), allocatable :: product

    class(fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running fluid_t%fluid_mul_real()', __FILE__, __LINE__)

    allocate(local_product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fluid_t%fluid_mul_real"

    ! local_product%conserved_vars = lhs%conserved_vars * rhs
    call mult_field_by_real(a=lhs%conserved_vars, b=rhs, c=local_product%conserved_vars)

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='fluid_mul_real (product)', line=__LINE__)
  end function fluid_mul_real

  function real_mul_fluid(lhs, rhs) result(product)
    !< Implementation of the real * fluid operation
    class(fluid_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    class(integrand_t), allocatable :: product

    type(fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running fluid_t%real_mul_fluid()', __FILE__, __LINE__)

    allocate(local_product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fluid_t%real_mul_fluid"

    local_product%time_integrator = rhs%time_integrator
    ! local_product%conserved_vars = rhs%conserved_vars * lhs
    call mult_field_by_real(a=rhs%conserved_vars, b=lhs, c=local_product%conserved_vars)

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='real_mul_fluid (product)', line=__LINE__)
  end function real_mul_fluid

  subroutine assign_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(fluid_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs
    integer(ik) :: alloc_status

    alloc_status = 0
    call debug_print('Running fluid_t%assign_fluid', __FILE__, __LINE__)

    call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
    select type(rhs)
    class is(fluid_t)
      lhs%time_integrator = rhs%time_integrator
      lhs%conserved_vars = rhs%conserved_vars
    class default
      error stop 'Error in fluid_t%assign_fluid: unsupported class'
    end select

    call rhs%clean_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
  end subroutine assign_fluid

  subroutine get_sound_speed(self, sound_speed)
    !< Calculate the sound speed for the entire domain
    class(fluid_t), intent(in) :: self
    real(rk), dimension(:, :), intent(out), allocatable :: sound_speed
    integer(ik) :: ilo, ihi, jlo, jhi

    call debug_print('Running fluid_t%calculate_sound_speed()', __FILE__, __LINE__)

    ilo = lbound(self%conserved_vars, dim=2)
    ihi = ubound(self%conserved_vars, dim=2)
    jlo = lbound(self%conserved_vars, dim=3)
    jhi = ubound(self%conserved_vars, dim=3)

    if(.not. allocated(sound_speed)) allocate(sound_speed(ilo:ihi, jlo:jhi))

    call eos%sound_speed_from_conserved(conserved_vars=self%conserved_vars, sound_speed=sound_speed)
  end subroutine get_sound_speed

  real(rk) function get_max_sound_speed(self) result(max_cs)
    !< Find the maximum sound speed in the domain
    class(fluid_t), intent(in) :: self
    real(rk), dimension(:, :), allocatable :: sound_speed
    integer(ik) :: ilo, ihi, jlo, jhi
    max_cs = 0.0_rk
    call debug_print('Running fluid_t%get_max_sound_speed()', __FILE__, __LINE__)

    call self%get_sound_speed(sound_speed)
    ilo = lbound(sound_speed, dim=1) + 1
    ihi = ubound(sound_speed, dim=1) - 1
    jlo = lbound(sound_speed, dim=2) + 1
    jhi = ubound(sound_speed, dim=2) - 1

    max_cs = maxval(sound_speed(ilo:ihi, jlo:jhi))
    deallocate(sound_speed)
  end function get_max_sound_speed

  subroutine get_primitive_vars(self, primitive_vars, fv)
    !< Convert the current conserved_vars [rho, rho u, rho v, rho E] to primitive [rho, u, v, p]. Note, the
    !< work is done in the EOS module due to the energy to pressure conversion
    class(fluid_t), intent(in) :: self
    class(finite_volume_scheme_t), intent(inout) :: fv
    real(rk), dimension(:, :, :), allocatable, intent(out) :: primitive_vars
    integer(ik) :: i
    integer(ik), dimension(3) :: bounds

    call debug_print('Running fluid_t%get_primitive_vars()', __FILE__, __LINE__)
    allocate(primitive_vars, mold=self%conserved_vars)

    call eos%conserved_to_primitive(self%conserved_vars, primitive_vars)
    bounds = lbound(primitive_vars)
    call fv%apply_primitive_vars_bc(primitive_vars, bounds)
  end subroutine get_primitive_vars

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(fluid_t), intent(in) :: self
    integer(ik) :: i, j, ilo, jlo, ihi, jhi
    logical :: invalid_numbers, negative_numbers
    integer(ik), intent(out) :: error_code

    error_code = 0

    negative_numbers = .false.
    invalid_numbers = .false.

    ilo = lbound(self%conserved_vars, dim=2)
    ihi = ubound(self%conserved_vars, dim=2)
    jlo = lbound(self%conserved_vars, dim=3)
    jhi = ubound(self%conserved_vars, dim=3)
    do j = jlo, jhi
      do i = ilo, ihi
        if(.not. any(ieee_is_finite(self%conserved_vars(:, i, j)))) then
          write(std_err, '(2(a,i0), a, 4(es10.3))') 'Infinite numbers in conserved_vars(:, ', i, ', ', j, ')', &
            self%conserved_vars(:, i, j)
          invalid_numbers = .true.
        end if

        if(any(ieee_is_nan(self%conserved_vars(:, i, j)))) then
          write(std_err, '(2(a,i0), a, 4(es10.3))') 'NaNs in conserved_vars(:, ', i, ', ', j, ')', &
            self%conserved_vars(:, i, j)
          invalid_numbers = .true.
        end if
      end do
    end do

    if(minval(self%conserved_vars(1, :, :)) < 0.0_rk) then
      write(std_err, '(a, 2(i0, 1x), a, es10.3)') "Error: Negative density at fluid_t%conserved_vars(1,i,j): (", &
        minloc(self%conserved_vars(1, :, :)), ") density = ", minval(self%conserved_vars(1, :, :))
      negative_numbers = .true.
    end if

    if(minval(self%conserved_vars(4, :, :)) < 0.0_rk) then
      write(std_err, '(a, 2(i0, 1x))') "Error: Negative rho E (density * total energy) at fluid_t%conserved_vars(4,i,j): ", &
        minloc(self%conserved_vars(4, :, :))
      negative_numbers = .true.
    end if

    if(invalid_numbers .or. negative_numbers) then
      write(std_out, '(a)') "Invalid or negative numbers in the conserved variables [rho, rho u, rho v, rho E]"
      write(std_err, '(a)') "Invalid or negative numbers in the conserved variables [rho, rho u, rho v, rho E]"
      error_code = 1
    end if

  end subroutine sanity_check

end module mod_fluid
