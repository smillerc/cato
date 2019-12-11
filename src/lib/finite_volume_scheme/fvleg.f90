module mod_fvleg

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_reconstruction_factory, only: reconstruction_factory_t
  use mod_local_evo_operator, only: local_evo_operator_t
  use mod_second_order_reconstruction, only: second_order_reconstruction_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_flux_tensor, only: H => flux_tensor_t
  use mod_grid_factory, only: grid_factory
  use mod_bc_factory, only: bc_factory
  use mod_grid, only: grid_t
  use hdf5_interface, only: hdf5_file
  use mod_time_integrator_factory, only: time_integrator_factory

  ! use mod_first_order_reconstruction, only: first_order_reconstruction_t

  implicit none
  private
  public :: fvleg_t, new_fvleg

  type, extends(finite_volume_scheme_t) :: fvleg_t
    !< Implementation of the finite volume local evolution Galerkin (FVLEG) scheme type

  contains
    procedure, public :: initialize
    procedure, public :: reconstruct => reconstruct_fvleg
    procedure, public :: evolve_domain
    procedure, public :: apply_conserved_vars_bc
    procedure, public :: apply_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc
    procedure, public :: t => time_derivative
    procedure, private :: integrate_fluxes
    procedure, private :: duplicate
    procedure, private :: initialize_from_hdf5
    procedure, private :: initialize_from_ini
    procedure, pass(lhs), public :: type_plus_type => add_fvleg
    procedure, pass(lhs), public :: type_minus_type => subtract_fvleg
    procedure, pass(lhs), public :: type_mul_real => fvleg_mul_real
    procedure, pass(rhs), public :: real_mul_type => real_mul_fvleg
    procedure, pass(lhs), public :: assign => assign_fvleg
    final :: finalize
  end type

  interface new_fvleg
    module procedure :: initialize
  end interface

contains

  subroutine initialize(self, input)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(fvleg_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    type(reconstruction_factory_t) :: recon_factory
    integer(ik) :: alloc_status

    self%title = input%title
    self%grid = grid_factory(input)

    self%time_integrator = time_integrator_factory(input)

    ! Set boundary conditions
    self%bc_plus_x = bc_factory(bc_type=input%plus_x_bc, location='+x')
    self%bc_plus_y = bc_factory(bc_type=input%plus_y_bc, location='+y')
    self%bc_minus_x = bc_factory(bc_type=input%minus_x_bc, location='-x')
    self%bc_minus_y = bc_factory(bc_type=input%minus_y_bc, location='-y')

    associate(imin=>self%grid%ilo_bc_cell, imax=>self%grid%ihi_bc_cell, &
              jmin=>self%grid%jlo_bc_cell, jmax=>self%grid%jhi_bc_cell)

      allocate(self%conserved_vars(4, imin:imax, jmin:jmax), stat=alloc_status)
      ! ((rho,u,v,p),i,j) Conserved variables for each cell
      if(alloc_status /= 0) then
        write(*, *) "Allocation status: ", alloc_status
        error stop "Unable to allocate fvleg_t%conserved_vars"
      end if

      allocate(self%reconstructed_state(4, 4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      ! ((rho, u ,v, p), point, node/midpoint, i, j); this is a cell-based value, so imax=ni-1, etc
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%reconstructed_state"
    end associate

    recon_factory = reconstruction_factory_t(input)
    self%reconstruction_operator = recon_factory%create_reconstruction(grid=self%grid)
    self%evolution_operator = local_evo_operator_t(input, self%grid, self%reconstruction_operator)

    associate(imin=>self%grid%ilo_node, imax=>self%grid%ihi_node, &
              jmin=>self%grid%jlo_node, jmax=>self%grid%jhi_node)

      allocate(self%evolved_corner_state(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_corner_state"
      ! ((rho,u,v,p), i, j); Reconstructed U at each corner

      allocate(self%corner_reference_state(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%corner_reference_state"
      ! ((rho, u ,v, p), i, j); Reference state (tilde) at each corner

      allocate(self%evolved_downup_midpoints_state(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_downup_midpoints_state"
      ! ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the up/down edges (edges 2 and 4)

      allocate(self%evolved_leftright_midpoints_state(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_leftright_midpoints_state"
      ! ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges (edges 1 and 3)

      allocate(self%downup_midpoints_reference_state(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%downup_midpoints_reference_state"
      ! ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the up/down edges (edges 2 and 4)

      allocate(self%leftright_midpoints_reference_state(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%leftright_midpoints_reference_state"
      ! ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges (edges 1 and 3)
    end associate

    if(input%read_init_cond_from_file) then
      call self%initialize_from_hdf5(input)
    else
      call self%initialize_from_ini(input)
    end if

    self%initiated = .true.
  end subroutine initialize

  subroutine initialize_from_hdf5(self, input)
    !< Initialize from an .hdf5 file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the hdf5 file
    class(fvleg_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    type(hdf5_file) :: h5

    real(rk), dimension(:, :), allocatable :: density
    real(rk), dimension(:, :), allocatable :: x_velocity
    real(rk), dimension(:, :), allocatable :: y_velocity
    real(rk), dimension(:, :), allocatable :: pressure

    integer(ik) :: i, j
    call h5%initialize(filename=input%initial_condition_file, status='old', action='r')
    call h5%get('/density', density)
    call h5%get('/x_velocity', x_velocity)
    call h5%get('/y_velocity', y_velocity)
    call h5%get('/pressure', pressure)
    call h5%finalize()

    if (minval(pressure) < 0.0_rk) then
      error stop "Pressure < 0 in initialize_from_hdf5"
    end if

    associate(imin=>self%grid%ilo_bc_cell, imax=>self%grid%ihi_bc_cell, &
              jmin=>self%grid%jlo_bc_cell, jmax=>self%grid%jhi_bc_cell)

      self%conserved_vars(1, imin:imax, jmin:jmax) = density    ! (1:imax,1:jmax)
      self%conserved_vars(2, imin:imax, jmin:jmax) = x_velocity ! (1:imax,1:jmax)
      self%conserved_vars(3, imin:imax, jmin:jmax) = y_velocity ! (1:imax,1:jmax)
      self%conserved_vars(4, imin:imax, jmin:jmax) = pressure   ! (1:imax,1:jmax)
    end associate
      
      
    end subroutine initialize_from_hdf5

  subroutine initialize_from_ini(self, input)
    !< Initialize from an .ini file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the .ini file
    class(fvleg_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    write(*, '(a,4(f0.3, 1x))') 'Initializing fvleg_t%conserved_vars to [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure
    self%conserved_vars(1, :, :) = input%init_density
    self%conserved_vars(2, :, :) = input%init_x_velocity
    self%conserved_vars(3, :, :) = input%init_y_velocity
    self%conserved_vars(4, :, :) = input%init_pressure

  end subroutine initialize_from_ini

  subroutine finalize(self)
    !< Implementation of the class cleanup
    type(fvleg_t), intent(inout) :: self
    integer(ik) :: alloc_status

    alloc_status = 0

    print *, 'Finalizing fvleg_t'

    if(allocated(self%evolution_operator)) then
      deallocate(self%evolution_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolution_operator"
    end if

    if(allocated(self%grid)) then
      deallocate(self%grid, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%grid"
    end if

    if(allocated(self%reconstruction_operator)) then
      deallocate(self%reconstruction_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstruction_operator"
    end if

    if(allocated(self%conserved_vars)) then
      deallocate(self%conserved_vars, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%conserved_vars"
    end if

    if(allocated(self%reconstructed_state)) then
      deallocate(self%reconstructed_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstructed_state"
    end if

    if(allocated(self%evolved_corner_state)) then
      deallocate(self%evolved_corner_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_corner_state"
    end if

    if(allocated(self%corner_reference_state)) then
      deallocate(self%corner_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%corner_reference_state"
    end if

    if(allocated(self%evolved_downup_midpoints_state)) then
      deallocate(self%evolved_downup_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_downup_midpoints_state"
    end if

    if(allocated(self%evolved_leftright_midpoints_state)) then
      deallocate(self%evolved_leftright_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_leftright_midpoints_state"
    end if

    if(allocated(self%downup_midpoints_reference_state)) then
      deallocate(self%downup_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%downup_midpoints_reference_state"
    end if

    if(allocated(self%leftright_midpoints_reference_state)) then
      deallocate(self%leftright_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%leftright_midpoints_reference_state"
    end if

  end subroutine finalize

  subroutine reconstruct_fvleg(self)
    !< Implementation of the FVLEG reconstruction.
    !< This reconstructs the entire grid at all the nodes/midpoints
    class(fvleg_t), intent(inout) :: self

    print *, 'Reconstructing ...'
    call self%reconstruction_operator%reconstruct_domain(self%conserved_vars, self%reconstructed_state)
  end subroutine reconstruct_fvleg

  subroutine apply_conserved_vars_bc(self)
    !< Apply the boundary conditions
    class(fvleg_t), intent(inout) :: self

    print *, 'Applying conserbed var bc ...'
    call self%bc_plus_x%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
    call self%bc_plus_y%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
    call self%bc_minus_x%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
    call self%bc_minus_y%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
  end subroutine apply_conserved_vars_bc

  subroutine apply_reconstructed_state_bc(self)
    !< Apply the boundary conditions
    class(fvleg_t), intent(inout) :: self

    print *, 'Applying bc reconstructed state ...'
    call self%bc_plus_x%apply_reconstructed_state_bc (reconstructed_state=self%reconstructed_state)
    call self%bc_plus_y%apply_reconstructed_state_bc (reconstructed_state=self%reconstructed_state)
    call self%bc_minus_x%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
    call self%bc_minus_y%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
  end subroutine apply_reconstructed_state_bc

  subroutine apply_cell_gradient_bc(self, r_omega)
    !< Apply the boundary conditions
    class(fvleg_t), intent(inout) :: self
    class(abstract_reconstruction_t), intent(in) :: r_omega
      select type(r_omega)
      class is (second_order_reconstruction_t)
        print *, 'Applying bc cell gradient state ...'
        call self%bc_plus_x%apply_cell_gradient_bc (cell_gradient=self%reconstruction_operator%cell_gradient)
        call self%bc_plus_y%apply_cell_gradient_bc (cell_gradient=self%reconstruction_operator%cell_gradient)
        call self%bc_minus_x%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
        call self%bc_minus_y%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
      end select
      
  end subroutine apply_cell_gradient_bc

  function time_derivative(self) result(dU_dt)
    ! //TODO: make pure
    !< Implementation of dU/dt = 1/Omega_ij Sum(1,k) Integral(H . n dl)
    class(fvleg_t), intent(in) :: self
    class(integrand_t), allocatable :: dU_dt
    type(fvleg_t), allocatable :: local_dU_dt
    ! integer(ik) :: alloc_status

    print *, 'Finding dU/dt...'
    allocate(fvleg_t :: local_dU_dt)
    local_dU_dt = self%duplicate()

    call local_dU_dt%calculate_reference_state()
    call local_dU_dt%apply_conserved_vars_bc()
    call local_dU_dt%reconstruct()                    ! reconstruct only the real domain // TODO: make sure init applies ghost conditions at the start?
    call local_dU_dt%apply_cell_gradient_bc(self%reconstruction_operator)
    call local_dU_dt%apply_reconstructed_state_bc()   ! ghost cells now pick up reconstructed and conserved states
    call local_dU_dt%evolve_domain()                  ! now evolve the flow (with b.c. and reconstructed state info)

    ! Set the RHS of Eq. 3 in the main text
    local_dU_dt%conserved_vars = local_dU_dt%integrate_fluxes()

    call move_alloc(local_dU_dt, dU_dt)

  end function time_derivative

  subroutine evolve_domain(self)
    class(fvleg_t), intent(inout) :: self

    print *, 'Evolving the domain...'

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    print *, 'Left/right midpoints...'
    call self%evolution_operator%evolve_leftright_midpoints(conserved_vars=self%conserved_vars, &
                                                            reconstructed_state=self%reconstructed_state, &
                                                            reference_state=self%leftright_midpoints_reference_state, &
                                                            evolved_state=self%evolved_leftright_midpoints_state)

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    print *, 'Down/up midpoints ..'
    call self%evolution_operator%evolve_downup_midpoints(conserved_vars=self%conserved_vars, &
                                                         reconstructed_state=self%reconstructed_state, &
                                                         reference_state=self%downup_midpoints_reference_state, &
                                                         evolved_state=self%evolved_downup_midpoints_state)

    ! Evolve, i.e. E0(R_omega), at all corner nodes
    print *, 'Corner nodes...'
    call self%evolution_operator%evolve_corners(conserved_vars=self%conserved_vars, &
                                                reconstructed_state=self%reconstructed_state, &
                                                reference_state=self%corner_reference_state, &
                                                evolved_state=self%evolved_corner_state)
    print *, 'Done'
  end subroutine

  function integrate_fluxes(self) result(rhs)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    use mod_flux_tensor, only: operator(.dot.)
    class(fvleg_t), intent(in) :: self
    real(rk), dimension(4) :: edge_flux

    integer(ik) :: ilo, ihi, jlo, jhi
    ! real(rk), dimension(size(self%conserved_vars, dim=1), &
    !                     size(self%conserved_vars, dim=2), &
    !                     size(self%conserved_vars, dim=3)) :: rhs !< RHS of Eq 3 in the main text

    real(rk), dimension(:, :, :), allocatable :: rhs

    ! type(flux_tensor_t) :: H1, H2, H3
    integer(ik) :: i, j, k

    ilo = self%grid%ilo_cell
    jlo = self%grid%jlo_cell
    ihi = self%grid%ihi_cell
    jhi = self%grid%jhi_cell

    print *, 'Integrating fluxes...'
    allocate(rhs, mold=self%conserved_vars)

    do j = jlo, jhi
      do i = ilo, ihi
    ! do concurrent(j=self%grid%jlo_cell:self%grid%jhi_cell)
    !   do concurrent(i=self%grid%ilo_cell:self%grid%ihi_cell)

        edge_flux = 0.0_rk
        print*, 'i, j:', i, j
        ! Edge 1 (bottom)
        print*, 'bottom'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i, j), &
                  E0_R_omega_kc=>self%evolved_leftright_midpoints_state(:, i, j), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i + 1, j), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 1, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(1, i, j))

          ! Eq. 13, for edge 1
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        ! Edge 2 (right)
        print*, 'right'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i + 1, j), &
                  E0_R_omega_kc=>self%evolved_downup_midpoints_state(:, i + 1, j), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i + 1, j + 1), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 2, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(2, i, j))

          ! Eq. 13, for edge 2
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        ! Edge 3 (top)
        print*, 'top'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i + 1, j + 1), &
                  E0_R_omega_kc=>self%evolved_leftright_midpoints_state(:, i, j + 1), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i, j + 1), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 3, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(3, i, j))

          ! Eq. 13, for edge 3
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        ! Edge 4 (left)
        print*, 'left'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i, j + 1), &
                  E0_R_omega_kc=>self%evolved_downup_midpoints_state(:, i, j), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i, j), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 4, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(4, i, j))

          ! Eq. 13, for edge 4
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        rhs(:, i, j) = (-1.0_rk / self%grid%cell_volume(i, j)) * edge_flux

      end do ! i
    end do ! j

  end function integrate_fluxes

  function subtract_fvleg(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fvleg type
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: difference
    type(fvleg_t), allocatable :: local_difference

    print *, 'Calling subtract_fvleg'
    select type(rhs)
    class is(fvleg_t)
      local_difference = lhs%duplicate()
      local_difference%conserved_vars = lhs%conserved_vars - rhs%conserved_vars
    class default
      error stop 'fvleg_t%subtract_fvleg: unsupported rhs class'
    end select

    call move_alloc(local_difference, difference)
  end function subtract_fvleg

  function add_fvleg(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fvleg type
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: sum
    type(fvleg_t), allocatable :: local_sum

    print *, 'Calling add_fvleg'
    select type(rhs)
    class is(fvleg_t)
      allocate(fvleg_t :: local_sum)
      local_sum = lhs%duplicate()
      local_sum%conserved_vars = lhs%conserved_vars + rhs%conserved_vars

    class default
      error stop 'fvleg_t%add_fvleg: unsupported rhs class'
    end select
    call move_alloc(local_sum, sum)
  end function add_fvleg

  function fvleg_mul_real(lhs, rhs) result(product)
    !< Implementation of the fvleg_t * real operation
    class(fvleg_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(integrand_t), allocatable :: product
    class(fvleg_t), allocatable :: local_product

    allocate(fvleg_t :: local_product)
    local_product = lhs%duplicate()
    local_product%conserved_vars = lhs%conserved_vars * rhs

    call move_alloc(local_product, product)

  end function fvleg_mul_real

  function real_mul_fvleg(lhs, rhs) result(product)
    !< Implementation of the real * fvleg_t operation
    class(fvleg_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    class(integrand_t), allocatable :: product
    class(fvleg_t), allocatable :: local_product

    allocate(fvleg_t :: local_product)
    local_product = rhs%duplicate()
    local_product%conserved_vars = rhs%conserved_vars * lhs
    call move_alloc(local_product, product)

  end function real_mul_fvleg

  subroutine assign_fvleg(lhs, rhs)
    !< Implementation of the (=) operator for the fvleg type. e.g. lhs = rhs
    class(fvleg_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs
    integer(ik) :: alloc_status

    alloc_status = 0
    print *, 'Calling assign_fvleg'
    select type(rhs)
    class is(fvleg_t)

      allocate(lhs%grid, source=rhs%grid, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate lhs%grid from rhs%grid in fvleg_t%assign"
      lhs%grid = rhs%grid

      lhs%bc_plus_x = rhs%bc_plus_x
      lhs%bc_plus_y = rhs%bc_plus_y
      lhs%bc_minus_x = rhs%bc_minus_x
      lhs%bc_minus_y = rhs%bc_minus_y

      if(.not. allocated(lhs%evolution_operator)) then
        allocate(lhs%evolution_operator, source=rhs%evolution_operator, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%evolution_operator from rhs%evolution_operator in fvleg_t%assign"
        end if
      end if
      lhs%evolution_operator = rhs%evolution_operator

      if(.not. allocated(lhs%time_integrator)) then
        allocate(lhs%time_integrator, source=rhs%time_integrator, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%time_integrator from rhs%time_integrator in fvleg_t%assign"
        endif
      end if
      lhs%time_integrator = rhs%time_integrator

      if(.not. allocated(lhs%reconstruction_operator)) then
        allocate(lhs%reconstruction_operator, source=rhs%reconstruction_operator, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%reconstruction_operator "// &
            "from rhs%reconstruction_operator in fvleg_t%assign"
        endif
      end if
      lhs%reconstruction_operator = rhs%reconstruction_operator

      if(.not. allocated(lhs%conserved_vars)) then
        allocate(lhs%conserved_vars, mold=rhs%conserved_vars, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%conserved_vars from "// &
            "rhs%conserved_vars in fvleg_t%assign"
        endif
      end if
      lhs%conserved_vars = rhs%conserved_vars

      if(.not. allocated(lhs%evolved_corner_state)) then
        allocate(lhs%evolved_corner_state, &
                 mold=rhs%evolved_corner_state, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%evolved_corner_state "// &
            "from rhs%evolved_corner_state in fvleg_t%assign"
        endif
      endif
      lhs%evolved_corner_state = rhs%evolved_corner_state

      if(.not. allocated(lhs%corner_reference_state)) then
        allocate(lhs%corner_reference_state, &
                 mold=rhs%corner_reference_state, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%corner_reference_state "// &
            "from rhs%corner_reference_state in fvleg_t%assign"
        endif
      endif
      lhs%corner_reference_state = rhs%corner_reference_state

      if(.not. allocated(lhs%evolved_downup_midpoints_state)) then
        allocate(lhs%evolved_downup_midpoints_state, &
                 mold=rhs%evolved_downup_midpoints_state, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%evolved_downup_midpoints_state "// &
            "from rhs%evolved_downup_midpoints_state in fvleg_t%assign"
        end if
      endif
      lhs%evolved_downup_midpoints_state = rhs%evolved_downup_midpoints_state

      if(.not. allocated(lhs%evolved_leftright_midpoints_state)) then
        allocate(lhs%evolved_leftright_midpoints_state, &
                 mold=rhs%evolved_leftright_midpoints_state, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%evolved_leftright_midpoints_state "// &
            "from rhs%evolved_leftright_midpoints_state in fvleg_t%assign"
        end if
      endif
      lhs%evolved_leftright_midpoints_state = rhs%evolved_leftright_midpoints_state

      if(.not. allocated(lhs%downup_midpoints_reference_state)) then
        allocate(lhs%downup_midpoints_reference_state, &
                 mold=rhs%downup_midpoints_reference_state, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%downup_midpoints_reference_state "// &
            "from rhs%downup_midpoints_reference_state in fvleg_t%assign"
        end if
      endif
      lhs%downup_midpoints_reference_state = rhs%downup_midpoints_reference_state

      if(.not. allocated(lhs%leftright_midpoints_reference_state)) then
        allocate(lhs%leftright_midpoints_reference_state, &
                 mold=rhs%leftright_midpoints_reference_state, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%leftright_midpoints_reference_state "// &
            "from rhs%leftright_midpoints_reference_state in fvleg_t%assign"
        end if
      endif
      lhs%leftright_midpoints_reference_state = rhs%leftright_midpoints_reference_state

      if(.not. allocated(lhs%reconstructed_state)) then
        allocate(lhs%reconstructed_state, mold=rhs%reconstructed_state, stat=alloc_status)
        if(alloc_status /= 0) then
          error stop "Unable to allocate lhs%reconstructed_state from "// &
            "rhs%reconstructed_state in fvleg_t%assign"
        end if
      endif
      lhs%reconstructed_state = rhs%reconstructed_state

    class default
      error stop 'fvleg_t%assign_fvleg: unsupported class'
    end select
  end subroutine assign_fvleg

  function duplicate(self) result(duplicate_fvleg)
    !< Duplicate the fvleg_t type
    class(fvleg_t), intent(in) :: self
    class(fvleg_t), allocatable :: duplicate_fvleg

    integer(ik) :: alloc_status

    allocate(fvleg_t :: duplicate_fvleg, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg"

    duplicate_fvleg%grid = self%grid
    duplicate_fvleg%evolution_operator = self%evolution_operator
    duplicate_fvleg%reconstruction_operator = self%reconstruction_operator
    duplicate_fvleg%time_integrator = self%time_integrator

    duplicate_fvleg%bc_plus_x = self%bc_plus_x
    duplicate_fvleg%bc_plus_y = self%bc_plus_y
    duplicate_fvleg%bc_minus_x = self%bc_minus_x
    duplicate_fvleg%bc_minus_y = self%bc_minus_y

    allocate(duplicate_fvleg%conserved_vars, mold=self%conserved_vars, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%conserved_vars"
    duplicate_fvleg%conserved_vars = self%conserved_vars

    allocate(duplicate_fvleg%evolved_corner_state, mold=self%evolved_corner_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%evolved_corner_state"
    duplicate_fvleg%evolved_corner_state = 0.0_rk

    allocate(duplicate_fvleg%corner_reference_state, mold=self%corner_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%corner_reference_state"
    duplicate_fvleg%corner_reference_state = 0.0_rk

    allocate(duplicate_fvleg%evolved_downup_midpoints_state, mold=self%evolved_downup_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%evolved_downup_midpoints_state"
    duplicate_fvleg%evolved_downup_midpoints_state = 0.0_rk

    allocate(duplicate_fvleg%evolved_leftright_midpoints_state, mold=self%evolved_leftright_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%evolved_leftright_midpoints_state"
    duplicate_fvleg%evolved_leftright_midpoints_state = 0.0_rk

    allocate(duplicate_fvleg%downup_midpoints_reference_state, mold=self%downup_midpoints_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%downup_midpoints_reference_state"
    duplicate_fvleg%downup_midpoints_reference_state = 0.0_rk

    allocate(duplicate_fvleg%leftright_midpoints_reference_state, mold=self%leftright_midpoints_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%leftright_midpoints_reference_state"
    duplicate_fvleg%leftright_midpoints_reference_state = 0.0_rk

    allocate(duplicate_fvleg%reconstructed_state, mold=self%reconstructed_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate duplicate_fvleg%reconstructed_state"
    duplicate_fvleg%reconstructed_state = 0.0_rk

  end function

end module mod_fvleg
