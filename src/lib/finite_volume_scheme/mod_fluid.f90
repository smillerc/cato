module mod_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use mod_globals, only: debug_print
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_floating_point_utils, only: near_zero
  use mod_surrogate, only: surrogate
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
    procedure, public :: get_max_sound_speed
    procedure, public :: get_primitive_vars
    procedure, nopass, private :: flux_edges
    procedure, pass(lhs), public :: type_plus_type => add_fluid
    procedure, pass(lhs), public :: type_minus_type => subtract_fluid
    procedure, pass(lhs), public :: type_mul_real => fluid_mul_real
    procedure, pass(rhs), public :: real_mul_type => real_mul_fluid
    procedure, pass(lhs), public :: assign => assign_fluid
    procedure, public :: force_finalization
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

      ! allocate(self%primitive_vars(4, imin:imax, jmin:jmax), stat=alloc_status)
      ! if(alloc_status /= 0) then
      !   error stop "Unable to allocate fluid_t%primitive_vars"
      ! end if
      ! self%primitive_vars = 0.0_rk

      ! allocate(self%sound_speed(imin:imax, jmin:jmax), stat=alloc_status)
      ! if(alloc_status /= 0) then
      !   error stop "Unable to allocate fluid_t%sound_speed"
      ! end if
      ! self%sound_speed = 0.0_rk
    end associate

    time_integrator => time_integrator_factory(input)
    allocate(self%time_integrator, source=time_integrator, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%time_integrator"
    deallocate(time_integrator)

    if(input%read_init_cond_from_file) then
      call self%initialize_from_hdf5(input, finite_volume_scheme)
    else
      call self%initialize_from_ini(input)
    end if

    ! call self%update_primitive_vars()
    ! call self%calculate_sound_speed()
    ! self%primitives_updated = .true.

    write(*, '(a)') 'Initial fluid stats'
    write(*, '(a)') '======================================'
    write(*, '(a, f0.3)') 'EOS Gamma:                ', eos%get_gamma()
    ! write(*, '(a, 2(es10.3, 1x))') 'Min/Max Sound Speed:      ', &
    !   minval(self%sound_speed), maxval(self%sound_speed)
    ! write(*, '(a, 2(es10.3, 1x))') 'Min/Max Density:          ', &
    !   minval(self%primitive_vars(1, :, :)), maxval(self%primitive_vars(1, :, :))
    ! write(*, '(a, 2(es10.3, 1x))') 'Min/Max X-Velocity:       ', &
    !   minval(self%primitive_vars(2, :, :)), maxval(self%primitive_vars(2, :, :))
    ! write(*, '(a, 2(es10.3, 1x))') 'Min/Max Y-Velocity:       ', &
    !   minval(self%primitive_vars(3, :, :)), maxval(self%primitive_vars(3, :, :))
    ! write(*, '(a, 2(es10.3, 1x))') 'Min/Max Pressure:         ', &
    !   minval(self%primitive_vars(4, :, :)), maxval(self%primitive_vars(4, :, :))
    ! write(*, *)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max rho:              ', &
      minval(self%conserved_vars(1, :, :)), maxval(self%conserved_vars(1, :, :))
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max rho u:            ', &
      minval(self%conserved_vars(2, :, :)), maxval(self%conserved_vars(2, :, :))
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max rho v:            ', &
      minval(self%conserved_vars(3, :, :)), maxval(self%conserved_vars(3, :, :))
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max rho E:            ', &
      minval(self%conserved_vars(4, :, :)), maxval(self%conserved_vars(4, :, :))
    write(*, *)

  end subroutine initialize

  subroutine initialize_from_hdf5(self, input, finite_volume_scheme)
    !< Initialize from an .hdf5 file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the hdf5 file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    type(hdf5_file) :: h5

    real(rk), dimension(:, :), allocatable :: density
    real(rk), dimension(:, :), allocatable :: x_velocity
    real(rk), dimension(:, :), allocatable :: y_velocity
    real(rk), dimension(:, :), allocatable :: pressure

    call debug_print('Initializing fluid_t from hdf5', __FILE__, __LINE__)
    call h5%initialize(filename=input%initial_condition_file, status='old', action='r')
    call h5%get('/density', density)
    call h5%get('/x_velocity', x_velocity)
    call h5%get('/y_velocity', y_velocity)
    call h5%get('/pressure', pressure)
    call h5%finalize()

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

  end subroutine initialize_from_hdf5

  subroutine initialize_from_ini(self, input)
    !< Initialize from an .ini file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the .ini file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    call debug_print('Initializing fluid_t from .ini', __FILE__, __LINE__)
    write(*, '(a,4(f0.3, 1x))') 'Initializing with [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure
    self%conserved_vars(1, :, :) = input%init_density
    self%conserved_vars(2, :, :) = input%init_density * input%init_x_velocity
    self%conserved_vars(3, :, :) = input%init_density * input%init_y_velocity
    self%conserved_vars(4, :, :) = input%init_density * eos%calculate_total_energy(pressure=input%init_pressure, &
                                                                                   density=input%init_density, &
                                                                                   x_velocity=input%init_x_velocity, &
                                                                                   y_velocity=input%init_y_velocity)

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
    integer(ik) :: alloc_status

    real(rk), dimension(:, :, :), allocatable :: primitive_vars
    !< ((rho, u, v, p), i, j); Primitive variables at each cell center

    real(rk), dimension(:, :, :), allocatable :: evolved_corner_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each corner

    real(rk), dimension(:, :, :), allocatable :: corner_reference_state
    !< ((rho, u, v, p), i, j); Reference state (tilde) at each corner

    ! Indexing the midpoints is a pain, so they're split by the up/down edges and left/right edges
    real(rk), dimension(:, :, :), allocatable :: evolved_downup_midpoints_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each midpoint defined by vectors that go up/down
    !< (edges 2 (right edge) and 4 (left edge))

    real(rk), dimension(:, :, :), allocatable :: evolved_leftright_midpoints_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each midpoint defined by vectors that go left/right
    !< (edges 1 (bottom edge) and 3 (top edge))

    real(rk), dimension(:, :, :), allocatable :: downup_midpoints_reference_state
    !< ((rho, u, v, p), i, j); Reference state (tilde) at each midpoint defined by vectors that go up/down
    !< (edges 2 (right edge) and 4 (left edge))

    real(rk), dimension(:, :, :), allocatable :: leftright_midpoints_reference_state
    !< ((rho, u, v, p), i, j); Reference state (tilde) at each midpoint defined by vectors that go left/right
    !< (edges 1 (bottom edge) and 3 (top edge))

    real(rk), dimension(:, :, :, :, :), allocatable, target :: reconstructed_state
    !< (((rho, u, v, p)), point, node/midpoint, i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints. Note, this DOES repeat nodes, since corners and midpoints are
    !< shared by neighboring cells, but each point has its own reconstructed value based on the parent cell's state

    call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

    associate(imin=>fv%grid%ilo_bc_cell, imax=>fv%grid%ihi_bc_cell, &
              jmin=>fv%grid%jlo_bc_cell, jmax=>fv%grid%jhi_bc_cell)

      allocate(reconstructed_state(4, 4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate reconstructed_state in fluid_t%time_derivative()"
      reconstructed_state = 0.0_rk
    end associate

    associate(imin_node=>fv%grid%ilo_node, imax_node=>fv%grid%ihi_node, &
              jmin_node=>fv%grid%jlo_node, jmax_node=>fv%grid%jhi_node, &
              imin_cell=>fv%grid%ilo_cell, imax_cell=>fv%grid%ihi_cell, &
              jmin_cell=>fv%grid%jlo_cell, jmax_cell=>fv%grid%jhi_cell)

      ! corners
      allocate(evolved_corner_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_corner_state in fluid_t%time_derivative()"
      evolved_corner_state = 0.0_rk

      allocate(corner_reference_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate corner_reference_state in fluid_t%time_derivative()"
      corner_reference_state = 0.0_rk

      ! left/right midpoints
      allocate(evolved_leftright_midpoints_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_leftright_midpoints_state in fluid_t%time_derivative()"
      evolved_leftright_midpoints_state = 0.0_rk

      allocate(leftright_midpoints_reference_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate leftright_midpoints_reference_state in fluid_t%time_derivative()"
      leftright_midpoints_reference_state = 0.0_rk

      ! down/up midpoints
      allocate(evolved_downup_midpoints_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_downup_midpoints_state"
      evolved_downup_midpoints_state = 0.0_rk

      allocate(downup_midpoints_reference_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate downup_midpoints_reference_state"
      downup_midpoints_reference_state = 0.0_rk
    end associate

    allocate(local_d_dt, source=self)
    ! local_d_dt%conserved_vars = self%conserved_vars

    allocate(primitive_vars, mold=self%conserved_vars)
    call local_d_dt%get_primitive_vars(primitive_vars)

    call fv%apply_source_terms(primitive_vars, lbound(primitive_vars))

    ! First put primitive vars in ghost layers
    call fv%apply_primitive_vars_bc(primitive_vars, lbound(primitive_vars))

    ! Reference state is a neighbor average (which is why edges needed ghost primitive variables)
    call fv%calculate_reference_state(primitive_vars=primitive_vars, &
                                      leftright_midpoints_reference_state=leftright_midpoints_reference_state, &
                                      downup_midpoints_reference_state=downup_midpoints_reference_state, &
                                      corner_reference_state=corner_reference_state, &
                                      cell_lbounds=lbound(primitive_vars))

    ! Now we can reconstruct the entire domain
    call fv%reconstruction_operator%set_primitive_vars_pointer(primitive_vars=primitive_vars, &
                                                               lbounds=lbound(primitive_vars))

    call fv%reconstruct(primitive_vars=primitive_vars, &
                        cell_lbounds=lbound(primitive_vars), &
                        reconstructed_state=reconstructed_state)

    call fv%evolution_operator%set_reconstructed_state_pointer( &
      reconstructed_state_target=reconstructed_state, &
      lbounds=lbound(reconstructed_state))

    ! Apply the reconstructed state to the ghost layers
    call fv%apply_reconstructed_state_bc(reconstructed_state, lbounds=lbound(reconstructed_state))

    call fv%apply_cell_gradient_bc()

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
    call fv%evolution_operator%evolve_leftright_midpoints(reference_state=leftright_midpoints_reference_state, &
                                                          evolved_state=evolved_leftright_midpoints_state, &
                                                          lbounds=lbound(leftright_midpoints_reference_state))

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
    call fv%evolution_operator%evolve_downup_midpoints(reference_state=downup_midpoints_reference_state, &
                                                       evolved_state=evolved_downup_midpoints_state, &
                                                       lbounds=lbound(evolved_downup_midpoints_state))

    ! Evolve, i.e. E0(R_omega), at all corner nodes
    call debug_print('Evolving corner nodes', __FILE__, __LINE__)
    call fv%evolution_operator%evolve_corners(reference_state=corner_reference_state, &
                                              evolved_state=evolved_corner_state, &
                                              lbounds=lbound(evolved_corner_state))

    nullify(fv%reconstruction_operator%primitive_vars)
    nullify(fv%evolution_operator%reconstructed_state)

    call self%flux_edges(grid=fv%grid, &
                         evolved_corner_state=evolved_corner_state, &
                         evolved_leftright_midpoints_state=evolved_leftright_midpoints_state, &
                         evolved_downup_midpoints_state=evolved_downup_midpoints_state, &
                         new_conserved_vars=local_d_dt%conserved_vars)

    call move_alloc(local_d_dt, d_dt)
    call d_dt%set_temp(calling_function='fluid_t%time_derivative (d_dt)', line=__LINE__)

    ! Now deallocate everything
    deallocate(primitive_vars, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate primitive_vars in fluid_t%time_derivative()"

    deallocate(evolved_corner_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate evolved_corner_state in fluid_t%time_derivative()"

    deallocate(corner_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate corner_reference_state in fluid_t%time_derivative()"

    deallocate(evolved_downup_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate evolved_downup_midpoints_state in fluid_t%time_derivative()"

    deallocate(evolved_leftright_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate evolved_leftright_midpoints_state in fluid_t%time_derivative()"

    deallocate(downup_midpoints_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate downup_midpoints_reference_state in fluid_t%time_derivative()"

    deallocate(leftright_midpoints_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate leftright_midpoints_reference_state in fluid_t%time_derivative()"

    deallocate(reconstructed_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate reconstructed_state"

  end function time_derivative

  subroutine flux_edges(grid, evolved_corner_state, evolved_leftright_midpoints_state, &
                        evolved_downup_midpoints_state, new_conserved_vars)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    class(grid_t), intent(in) :: grid
    real(rk), dimension(:, grid%ilo_node:, grid%jlo_node:), intent(in) :: evolved_corner_state
    real(rk), dimension(:, grid%ilo_cell:, grid%jlo_node:), intent(in) :: evolved_leftright_midpoints_state
    real(rk), dimension(:, grid%ilo_node:, grid%jlo_cell:), intent(in) :: evolved_downup_midpoints_state
    real(rk), dimension(:, grid%ilo_bc_cell:, grid%jlo_bc_cell:), intent(out) :: new_conserved_vars

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j
    real(rk), dimension(4) :: top_left_corner, top_right_corner, bottom_left_corner, bottom_right_corner
    real(rk), dimension(4) :: bottom_midpoint, right_midpoint, top_midpoint, left_midpoint
    real(rk), dimension(2, 4) :: n_hat
    real(rk), dimension(4) :: delta_l
    real(rk), dimension(4) :: edge_flux_1
    real(rk), dimension(4) :: edge_flux_2
    real(rk), dimension(4) :: edge_flux_3
    real(rk), dimension(4) :: edge_flux_4

    ilo = grid%ilo_cell
    ihi = grid%ihi_cell
    jlo = grid%jlo_cell
    jhi = grid%jhi_cell

    edge_flux_1 = 0.0_rk
    edge_flux_2 = 0.0_rk
    edge_flux_3 = 0.0_rk
    edge_flux_4 = 0.0_rk

    top_left_corner = 0.0_rk
    top_right_corner = 0.0_rk
    bottom_left_corner = 0.0_rk
    bottom_right_corner = 0.0_rk
    bottom_midpoint = 0.0_rk
    right_midpoint = 0.0_rk
    top_midpoint = 0.0_rk
    left_midpoint = 0.0_rk

    call debug_print('Running fluid_t%flux_edges()', __FILE__, __LINE__)

    new_conserved_vars(:, grid%ilo_bc_cell, :) = 0.0_rk
    new_conserved_vars(:, grid%ihi_bc_cell, :) = 0.0_rk
    new_conserved_vars(:, :, grid%jlo_bc_cell) = 0.0_rk
    new_conserved_vars(:, :, grid%jhi_bc_cell) = 0.0_rk

    do concurrent(j=jlo:jhi)
      do concurrent(i=ilo:ihi)

        top_left_corner = evolved_corner_state(:, i, j + 1)
        bottom_left_corner = evolved_corner_state(:, i, j)
        top_right_corner = evolved_corner_state(:, i + 1, j + 1)
        bottom_right_corner = evolved_corner_state(:, i + 1, j)
        bottom_midpoint = evolved_leftright_midpoints_state(:, i, j)
        top_midpoint = evolved_leftright_midpoints_state(:, i, j + 1)
        right_midpoint = evolved_downup_midpoints_state(:, i + 1, j)
        left_midpoint = evolved_downup_midpoints_state(:, i, j)
        delta_l = grid%cell_edge_lengths(:, i, j)
        n_hat = grid%cell_edge_norm_vectors(:, :, i, j)

        ! Edge 1 (bottom)
        edge_flux_1 = ( &
                      ((H(bottom_left_corner) + &
                        4.0_rk * H(bottom_midpoint) + &
                        H(bottom_right_corner)) .dot.n_hat(:, 1)) * (delta_l(1) / 6.0_rk))

        ! associate(E0_R_omega_k1=>evolved_corner_state(:, i, j), &
        !           E0_R_omega_kc=>evolved_leftright_midpoints_state(:, i, j), &
        !           E0_R_omega_k2=>evolved_corner_state(:, i + 1, j), &
        !           n_hat=>grid%cell_edge_norm_vectors(:, 1, i, j), &
        !           delta_l=>grid%cell_edge_lengths(1, i, j))

        !   ! Eq. 13, for edge 1
        !   edge_flux = edge_flux + &
        !               ( &
        !               ((H(E0_R_omega_k1) + &
        !                 4.0_rk * H(E0_R_omega_kc) + &
        !                 H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        ! end associate

        ! Edge 2 (right)
        edge_flux_2 = ( &
                      ((H(bottom_right_corner) + &
                        4.0_rk * H(right_midpoint) + &
                        H(top_right_corner)) .dot.n_hat(:, 2)) * (delta_l(2) / 6.0_rk))
        ! associate(E0_R_omega_k1=>evolved_corner_state(:, i + 1, j), &
        !           E0_R_omega_kc=>evolved_downup_midpoints_state(:, i + 1, j), &
        !           E0_R_omega_k2=>evolved_corner_state(:, i + 1, j + 1), &
        !           n_hat=>grid%cell_edge_norm_vectors(:, 2, i, j), &
        !           delta_l=>grid%cell_edge_lengths(2, i, j))

        !   ! Eq. 13, for edge 2
        !   edge_flux = edge_flux + &
        !               ( &
        !               ((H(E0_R_omega_k1) + &
        !                 4.0_rk * H(E0_R_omega_kc) + &
        !                 H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        ! end associate

        ! Edge 3 (top)
        edge_flux_3 = ( &
                      ((H(top_right_corner) + &
                        4.0_rk * H(top_midpoint) + &
                        H(top_left_corner)) .dot.n_hat(:, 3)) * (delta_l(3) / 6.0_rk))

        ! associate(E0_R_omega_k1=>evolved_corner_state(:, i + 1, j + 1), &
        !           E0_R_omega_kc=>evolved_leftright_midpoints_state(:, i, j + 1), &
        !           E0_R_omega_k2=>evolved_corner_state(:, i, j + 1), &
        !           n_hat=>grid%cell_edge_norm_vectors(:, 3, i, j), &
        !           delta_l=>grid%cell_edge_lengths(3, i, j))

        !   ! Eq. 13, for edge 3
        !   edge_flux = edge_flux + &
        !               ( &
        !               ((H(E0_R_omega_k1) + &
        !                 4.0_rk * H(E0_R_omega_kc) + &
        !                 H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        ! end associate

        ! Edge 4 (left)
        edge_flux_4 = ( &
                      ((H(top_left_corner) + &
                        4.0_rk * H(left_midpoint) + &
                        H(bottom_left_corner)) .dot.n_hat(:, 4)) * (delta_l(4) / 6.0_rk))

        ! associate(E0_R_omega_k1=>evolved_corner_state(:, i, j + 1), &
        !           E0_R_omega_kc=>evolved_downup_midpoints_state(:, i, j), &
        !           E0_R_omega_k2=>evolved_corner_state(:, i, j), &
        !           n_hat=>grid%cell_edge_norm_vectors(:, 4, i, j), &
        !           delta_l=>grid%cell_edge_lengths(4, i, j))

        !   ! Eq. 13, for edge 4
        !   edge_flux = edge_flux + &
        !               ( &
        !               ((H(E0_R_omega_k1) + &
        !                 4.0_rk * H(E0_R_omega_kc) + &
        !                 H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        ! end associate

        new_conserved_vars(:, i, j) = (-1.0_rk / grid%cell_volume(i, j)) * (edge_flux_1 + edge_flux_2 + edge_flux_3 + edge_flux_4)
      end do ! i
    end do ! j
    ! new_conserved_vars(3,:,:) = 0.0_rk
  end subroutine flux_edges

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
      local_difference%conserved_vars = lhs%conserved_vars - rhs%conserved_vars
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
      local_sum%conserved_vars = lhs%conserved_vars + rhs%conserved_vars
    class default
      error stop 'fluid_t%add_fluid: unsupported rhs class'
    end select

    call local_sum%set_temp(calling_function='fluid_t%add_fluid(local_sum)', line=__LINE__)
    call move_alloc(local_sum, sum)
    call sum%set_temp(calling_function='fluid_t%add_fluid(sum)', line=__LINE__)
  end function add_fluid

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

    local_product%conserved_vars = lhs%conserved_vars * rhs

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
    local_product%conserved_vars = rhs%conserved_vars * lhs

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

    ilo = lbound(self%conserved_vars, dim=2)
    ihi = ubound(self%conserved_vars, dim=2)
    jlo = lbound(self%conserved_vars, dim=3)
    jhi = ubound(self%conserved_vars, dim=3)

    allocate(sound_speed(ilo:ihi, jlo:jhi))

    call debug_print('Running fluid_t%calculate_sound_speed()', __FILE__, __LINE__)

    call eos%sound_speed_from_conserved(conserved_vars=self%conserved_vars, &
                                        sound_speed=sound_speed)
  end subroutine get_sound_speed

  real(rk) function get_max_sound_speed(self) result(max_cs)
    !< Find the maximum sound speed in the domain
    class(fluid_t), intent(in) :: self
    real(rk), dimension(:, :), allocatable :: sound_speed

    call debug_print('Running fluid_t%get_max_sound_speed()', __FILE__, __LINE__)

    call self%get_sound_speed(sound_speed)
    max_cs = maxval(sound_speed)
  end function get_max_sound_speed

  pure subroutine get_primitive_vars(self, primitive_vars)
    !< Convert the current conserved_vars [rho, rho u, rho v, rho E] to primitive [rho, u, v, p]. Note, the
    !< work is done in the EOS module due to the energy to pressure conversion
    class(fluid_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(out) :: primitive_vars

    call eos%conserved_to_primitive(self%conserved_vars, primitive_vars)
  end subroutine get_primitive_vars

end module mod_fluid
