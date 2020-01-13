module mod_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

  use mod_globals, only: debug_print
  use mod_parallel
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_floating_point_utils, only: near_zero
  use mod_surrogate, only: surrogate
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_surrogate, only: surrogate
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use hdf5_interface, only: hdf5_file
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_flux_tensor, only: operator(.dot.), H => flux_tensor_t
  use mod_time_integrator_factory, only: time_integrator_factory

  implicit none

  private
  public :: fluid_t, new_fluid

  type, extends(integrand_t) :: fluid_t
    real(rk), dimension(:, :, :), allocatable :: conserved_vars
    integer(ik) :: lb(2)         !< Lower bounds
    integer(ik) :: ub(2)         !< Upper bounds
    integer(ik) :: dims(2)       !< Data dimensions
    integer(ik) :: neighbors(4)  !< Neighbor images
    integer(ik) :: edge_size
    ! real(rk) :: timestep
    ! real(rk) :: time
    ! integer(ik) :: iteration
  contains
    procedure, public :: initialize
    procedure, private :: initialize_from_ini
    procedure, private :: initialize_from_hdf5
    procedure, public :: t => time_derivative
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
  end

  subroutine initialize(self, input, finite_volume_scheme)
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    class(strategy), pointer :: time_integrator => null()
    integer(ik) :: edge_size, indices(4)
    integer(ik) :: alloc_status

    alloc_status = 0
    call debug_print('Initializing fluid_t', __FILE__, __LINE__)

    associate(imin=>finite_volume_scheme%grid%ilo_bc_cell, &
              imax=>finite_volume_scheme%grid%ihi_bc_cell, &
              jmin=>finite_volume_scheme%grid%jlo_bc_cell, &
              jmax=>finite_volume_scheme%grid%jhi_bc_cell, &
              ni=>finite_volume_scheme%grid%ni_cell, &
              nj=>finite_volume_scheme%grid%nj_cell)

      indices = tile_indices([ni, nj])
      self%lb = indices([1, 3])
      self%ub = indices([2, 4])
      self%neighbors = tile_neighbors_2d(periodic=.true.)
      self%edge_size = max(self%ub(1) - self%lb(1) + 1, &
                           self%ub(2) - self%lb(2) + 1)
      call co_max(self%edge_size)
      allocate(self%conserved_vars(4, imin:imax, jmin:jmax), stat=alloc_status)
      ! ((rho,u,v,p),i,j) Conserved variables for each cell
      if(alloc_status /= 0) then
        error stop "Unable to allocate fluid_t%conserved_vars"
      end if
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

    write(*, '(a)') 'Initial fluid stats'
    write(*, '(a)') '==================='
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max Density:  ', minval(self%conserved_vars(1, :, :)), maxval(self%conserved_vars(1, :, :))
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max U:        ', minval(self%conserved_vars(2, :, :)), maxval(self%conserved_vars(2, :, :))
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max V:        ', minval(self%conserved_vars(3, :, :)), maxval(self%conserved_vars(3, :, :))
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max Pressure: ', minval(self%conserved_vars(4, :, :)), maxval(self%conserved_vars(4, :, :))
    write(*, *)
  end subroutine

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

      self%conserved_vars(1, imin:imax, jmin:jmax) = density    ! (1:imax,1:jmax)
      self%conserved_vars(2, imin:imax, jmin:jmax) = x_velocity ! (1:imax,1:jmax)
      self%conserved_vars(3, imin:imax, jmin:jmax) = y_velocity ! (1:imax,1:jmax)
      self%conserved_vars(4, imin:imax, jmin:jmax) = pressure   ! (1:imax,1:jmax)
    end associate

  end subroutine initialize_from_hdf5

  subroutine initialize_from_ini(self, input)
    !< Initialize from an .ini file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the .ini file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    call debug_print('Initializing fluid_t from .ini', __FILE__, __LINE__)
    write(*, '(a,4(f0.3, 1x))') 'Initializing fluid_t%conserved_vars to [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure
    self%conserved_vars(1, :, :) = input%init_density
    self%conserved_vars(2, :, :) = input%init_x_velocity
    self%conserved_vars(3, :, :) = input%init_y_velocity
    self%conserved_vars(4, :, :) = input%init_pressure

    if(near_zero(input%init_pressure)) then
      error stop "Some (or all) of the pressure array is ~0 in fluid_t%initialize_from_hdf5"
    end if

    if(near_zero(input%init_density)) then
      error stop "Some (or all) of the density array is ~0 in fluid_t%initialize_from_hdf5"
    end if

  end subroutine initialize_from_ini

  subroutine force_finalization(self)
    class(fluid_t), intent(inout) :: self

    call debug_print('Calling fluid_t%force_finalization()', __FILE__, __LINE__)
    if(allocated(self%conserved_vars)) deallocate(self%conserved_vars)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine

  subroutine finalize(self)
    type(fluid_t), intent(inout) :: self

    call debug_print('Calling fluid_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%conserved_vars)) deallocate(self%conserved_vars)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine

  function time_derivative(self, finite_volume_scheme) result(d_dt)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    class(fluid_t), intent(in) :: self
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme

    class(integrand_t), allocatable :: d_dt
    type(fluid_t), allocatable :: local_d_dt

    real(rk), dimension(4) :: edge_flux
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j

    call debug_print('Calculating d/dt, e.g. integrating edge fluxes', __FILE__, __LINE__)
    allocate(local_d_dt, source=self)

    ilo = finite_volume_scheme%grid%ilo_cell
    jlo = finite_volume_scheme%grid%jlo_cell
    ihi = finite_volume_scheme%grid%ihi_cell
    jhi = finite_volume_scheme%grid%jhi_cell

    associate(grid=>finite_volume_scheme%grid, &
              corner_state=>finite_volume_scheme%evolved_corner_state, &
              leftright_midpoints_state=>finite_volume_scheme%evolved_leftright_midpoints_state, &
              downup_midpoints_state=>finite_volume_scheme%evolved_downup_midpoints_state)

      ! do concurrent(j=jlo:jhi)
      !   do concurrent(i=ilo:ihi)
      do j = jlo, jhi
        do i = ilo, ihi

          edge_flux = 0.0_rk
          ! Edge 1 (bottom)
          associate(E0_R_omega_k1=>corner_state(:, i, j), &
                    E0_R_omega_kc=>leftright_midpoints_state(:, i, j), &
                    E0_R_omega_k2=>corner_state(:, i + 1, j), &
                    n_hat=>grid%cell_edge_norm_vectors(:, 1, i, j), &
                    delta_l=>grid%cell_edge_lengths(1, i, j))

            ! Eq. 13, for edge 1
            edge_flux = edge_flux + &
                        ( &
                        ((H(E0_R_omega_k1) + &
                          4.0_rk * H(E0_R_omega_kc) + &
                          H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

          end associate

          ! Edge 2 (right)
          associate(E0_R_omega_k1=>corner_state(:, i + 1, j), &
                    E0_R_omega_kc=>downup_midpoints_state(:, i + 1, j), &
                    E0_R_omega_k2=>corner_state(:, i + 1, j + 1), &
                    n_hat=>grid%cell_edge_norm_vectors(:, 2, i, j), &
                    delta_l=>grid%cell_edge_lengths(2, i, j))

            ! Eq. 13, for edge 2
            edge_flux = edge_flux + &
                        ( &
                        ((H(E0_R_omega_k1) + &
                          4.0_rk * H(E0_R_omega_kc) + &
                          H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

          end associate

          ! Edge 3 (top)
          associate(E0_R_omega_k1=>corner_state(:, i + 1, j + 1), &
                    E0_R_omega_kc=>leftright_midpoints_state(:, i, j + 1), &
                    E0_R_omega_k2=>corner_state(:, i, j + 1), &
                    n_hat=>grid%cell_edge_norm_vectors(:, 3, i, j), &
                    delta_l=>grid%cell_edge_lengths(3, i, j))

            ! Eq. 13, for edge 3
            edge_flux = edge_flux + &
                        ( &
                        ((H(E0_R_omega_k1) + &
                          4.0_rk * H(E0_R_omega_kc) + &
                          H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

          end associate

          ! Edge 4 (left)
          associate(E0_R_omega_k1=>corner_state(:, i, j + 1), &
                    E0_R_omega_kc=>downup_midpoints_state(:, i, j), &
                    E0_R_omega_k2=>corner_state(:, i, j), &
                    n_hat=>grid%cell_edge_norm_vectors(:, 4, i, j), &
                    delta_l=>grid%cell_edge_lengths(4, i, j))

            ! Eq. 13, for edge 4
            edge_flux = edge_flux + &
                        ( &
                        ((H(E0_R_omega_k1) + &
                          4.0_rk * H(E0_R_omega_kc) + &
                          H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

          end associate

          local_d_dt%conserved_vars(:, i, j) = (-1.0_rk / grid%cell_volume(i, j)) * edge_flux

        end do ! i
      end do ! j
    end associate

    call move_alloc(local_d_dt, d_dt)
    call d_dt%set_temp(calling_function='time_derivative (d_dt)', line=__LINE__)

  end function time_derivative

  function subtract_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: difference
    type(fluid_t), allocatable :: local_difference

    call debug_print('Calling fluid%subtract_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(fluid_t)
      allocate(local_difference, source=lhs)
      local_difference%conserved_vars = lhs%conserved_vars - rhs%conserved_vars
    class default
      error stop 'fluid%subtract_fluid: unsupported rhs class'
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

    call debug_print('Calling fluid%add_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(fluid_t)
      allocate(local_sum, source=lhs)
      local_sum%conserved_vars = lhs%conserved_vars + rhs%conserved_vars
    class default
      error stop 'fluid%add_fluid: unsupported rhs class'
    end select

    call move_alloc(local_sum, sum)
    call sum%set_temp(calling_function='add_fluid(sum)', line=__LINE__)
  end function add_fluid

  function fluid_mul_real(lhs, rhs) result(product)
    !< Implementation of the fluid * real operation
    class(fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(integrand_t), allocatable :: product

    type(fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Calling fluid%fluid_mul_real()', __FILE__, __LINE__)

    allocate(local_product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fluid%fluid_mul_real"

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

    call debug_print('Calling fluid%real_mul_fluid()', __FILE__, __LINE__)

    allocate(local_product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fluid%fluid_mul_real"

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
    call debug_print('Calling assign_fluid', __FILE__, __LINE__)

    call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
    select type(rhs)
    class is(fluid_t)
      lhs%time_integrator = rhs%time_integrator
      lhs%conserved_vars = rhs%conserved_vars
    class default
      error stop 'fluid%assign_fluid: unsupported class'
    end select

    call rhs%clean_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
  end subroutine assign_fluid

end module mod_fluid
