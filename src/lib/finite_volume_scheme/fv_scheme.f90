module mod_finite_volume_schemes

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use mod_globals, only: debug_print
  use mod_floating_point_utils, only: near_zero
  use mod_input, only: input_t
  use mod_reconstruction_factory, only: reconstruction_factory
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_second_order_reconstruction, only: second_order_reconstruction_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_flux_tensor, only: H => flux_tensor_t
  use mod_grid_factory, only: grid_factory
  use mod_bc_factory, only: bc_factory
  use mod_evo_operator_factory, only: evo_operator_factory
  use mod_grid, only: grid_t
  use hdf5_interface, only: hdf5_file
  use mod_boundary_conditions, only: boundary_condition_t

  implicit none
  private
  public :: finite_volume_scheme_t, make_fv_scheme

  type :: finite_volume_scheme_t
    !< Abstract representation of the finite volume scheme. This is essentially a puppeteer that
    !< manages the reconstruction, grid, and evolution operator to ultimately calculate the conserved
    !< state variables of each finite cell. The reconstruction, grid, and evolution implementations are passed
    !< on to decendents like the FVLEG scheme.

    character(len=32) :: title = ''
    integer(ik) :: timestep = 0
    real(rk) :: delta_t = 0.0_rk
    real(rk) :: time = 0.0_rk
    logical :: initiated = .false.

    class(abstract_reconstruction_t), allocatable :: reconstruction_operator
    !< R_Omega reconstruction operator used to reconstruct the corners/midpoints based on the cell
    !< average (and gradient if high(er) order reconstruction used)

    ! real(rk), dimension(:, :, :), allocatable :: conserved_vars
    ! !< ((rho, u ,v, p), i, j); Conserved variables for each cell

    class(abstract_evo_operator_t), allocatable :: evolution_operator
    !< Evolution operator to construct (rho, u, v, p) at each corner and midpoint

    class(grid_t), allocatable :: grid
    !< Grid class to hold geometric information (edge lengths, volumes, etc.)

    class(boundary_condition_t), allocatable :: bc_plus_x
    class(boundary_condition_t), allocatable :: bc_plus_y
    class(boundary_condition_t), allocatable :: bc_minus_x
    class(boundary_condition_t), allocatable :: bc_minus_y

    ! Corner/midpoint index convention
    ! --------------------------------
    !
    !   C----M----C----M----C
    !   |         |         |
    !   O    x    O    x    O
    !   |         |         |
    !   C----M----C----M----C
    !   |         |         |
    !   O    x    O    x    O
    !   |         |         |
    !   C----M----C----M----C

    ! C: corner, M: left/right midpoint, O: up/down midpoint, x: cell

    ! Since the reconstructed state at the corners (C) and midpoints (M) are reused by multiple cells,
    ! the datastructures are set up for maximum reuse.
    ! If they were indexed via cell, each cell would duplicate information since they share corners and midpoints

    real(rk), dimension(:, :, :), allocatable :: evolved_corner_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each corner

    real(rk), dimension(:, :, :), allocatable :: corner_reference_state
    !< ((rho, u ,v, p), i, j); Reference state (tilde) at each corner

    ! Indexing the midpoints is a pain, so they're split by the up/down edges and left/right edges

    real(rk), dimension(:, :, :), allocatable :: evolved_downup_midpoints_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the up/down edges (edges 2 and 4)

    real(rk), dimension(:, :, :), allocatable :: evolved_leftright_midpoints_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges (edges 1 and 3)

    real(rk), dimension(:, :, :), allocatable :: downup_midpoints_reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the down/up edges (edges 2 and 4)

    real(rk), dimension(:, :, :), allocatable :: leftright_midpoints_reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges (edges 1 and 3)

    real(rk), dimension(:, :, :, :, :), allocatable :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints. Note, this DOES repeat nodes, since corners and midpoints are
    !< shared by neighboring cells, but each point has its own reconstructed value based on the parent cell's state

  contains
    procedure, public :: calculate_reference_state
    procedure, public :: initialize
    procedure, public :: reconstruct
    procedure, public :: apply_conserved_vars_bc
    procedure, public :: evolve_domain
    procedure, public :: apply_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc
    procedure, public :: apply_source_terms
    final :: finalize
  end type

  interface make_fv_scheme
    module procedure :: constructor
  end interface

contains

  function constructor(input) result(fv)
    class(input_t), intent(in) :: input
    type(finite_volume_scheme_t), pointer :: fv

    allocate(fv)
    call fv%initialize(input)

  end function

  subroutine initialize(self, input)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(finite_volume_scheme_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    class(boundary_condition_t), pointer :: bc => null()
    class(grid_t), pointer :: grid => null()
    class(abstract_reconstruction_t), pointer :: r_omega => null()
    class(abstract_evo_operator_t), pointer :: E0 => null()

    integer(ik) :: alloc_status
    alloc_status = 0

    call debug_print('Initializing finite_volume_scheme_t', __FILE__, __LINE__)

    self%title = trim(input%title)

    grid => grid_factory(input)
    allocate(self%grid, source=grid, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%grid"
    deallocate(grid)

    ! Set boundary conditions
    bc => bc_factory(bc_type=input%plus_x_bc, location='+x', input=input)
    allocate(self%bc_plus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_plus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%plus_y_bc, location='+y', input=input)
    allocate(self%bc_plus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_plus_y"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_x_bc, location='-x', input=input)
    allocate(self%bc_minus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_minus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_y_bc, location='-y', input=input)
    allocate(self%bc_minus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_minus_y"
    deallocate(bc)

    write(*, '(a)') "Boundary Conditions"
    write(*, '(a)') "==================="
    write(*, '(3(a),i0,a)') "+x: ", trim(self%bc_plus_x%name), ' (priority = ', self%bc_plus_x%priority, ')'
    write(*, '(3(a),i0,a)') "-x: ", trim(self%bc_minus_x%name), ' (priority = ', self%bc_minus_x%priority, ')'
    write(*, '(3(a),i0,a)') "+y: ", trim(self%bc_plus_y%name), ' (priority = ', self%bc_plus_y%priority, ')'
    write(*, '(3(a),i0,a)') "-y: ", trim(self%bc_minus_y%name), ' (priority = ', self%bc_minus_y%priority, ')'
    write(*, *)

    associate(imin=>self%grid%ilo_bc_cell, imax=>self%grid%ihi_bc_cell, &
              jmin=>self%grid%jlo_bc_cell, jmax=>self%grid%jhi_bc_cell)

      allocate(self%reconstructed_state(4, 4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      ! ((rho, u ,v, p), point, node/midpoint, i, j); this is a cell-based value, so imax=ni-1, etc
      if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%reconstructed_state"
    end associate

    ! r_omega => reconstruction_factory(input=input, grid_target=self%grid, &
    !                                   conserved_vars_target=self%conserved_vars, lbounds=lbound(self%conserved_vars))
    r_omega => reconstruction_factory(input=input, grid_target=self%grid)
    allocate(self%reconstruction_operator, source=r_omega, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%reconstruction_operator"
    deallocate(r_omega)

    call debug_print('Making an E0 operator', __FILE__, __LINE__)
    E0 => evo_operator_factory(input=input, grid_target=self%grid, &
                               recon_operator_target=self%reconstruction_operator, &
                               reconstructed_state_target=self%reconstructed_state, lbounds=lbound(self%reconstructed_state))

    allocate(self%evolution_operator, source=E0, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%evolution_operator"
    deallocate(E0)

    associate(imin_node=>self%grid%ilo_node, imax_node=>self%grid%ihi_node, &
              jmin_node=>self%grid%jlo_node, jmax_node=>self%grid%jhi_node, &
              imin_cell=>self%grid%ilo_node, imax_cell=>self%grid%ihi_node, &
              jmin_cell=>self%grid%jlo_node, jmax_cell=>self%grid%jhi_node)

      ! corners
      if(.not. allocated(self%evolved_corner_state)) then
        ! ((rho,u,v,p), i, j); Reconstructed U at each corner
        allocate(self%evolved_corner_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%evolved_corner_state"
      end if

      if(.not. allocated(self%corner_reference_state)) then
        ! ((rho, u ,v, p), i, j); Reference state (tilde) at each corner
        allocate(self%corner_reference_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%corner_reference_state"
      endif

      ! left/right midpoints
      if(.not. allocated(self%evolved_leftright_midpoints_state)) then
        ! ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges (edges 1 and 3)
        allocate(self%evolved_leftright_midpoints_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%evolved_leftright_midpoints_state"
      end if

      if(.not. allocated(self%leftright_midpoints_reference_state)) then
        ! ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges (edges 1 and 3)
        allocate(self%leftright_midpoints_reference_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%leftright_midpoints_reference_state"
      end if

      ! down/up midpoints
      if(.not. allocated(self%evolved_downup_midpoints_state)) then
        ! ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the down/up edges (edges 2 and 4)
        allocate(self%evolved_downup_midpoints_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%evolved_downup_midpoints_state"
      end if

      if(.not. allocated(self%downup_midpoints_reference_state)) then
        ! ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the down/up edges (edges 2 and 4)
        allocate(self%downup_midpoints_reference_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%downup_midpoints_reference_state"
      end if

    end associate

    ! if(input%read_init_cond_from_file) then
    !   call self%initialize_from_hdf5(input)
    ! else
    !   call self%initialize_from_ini(input)
    ! end if

    self%initiated = .true.
  end subroutine initialize

  subroutine finalize(self)
    !< Implementation of the class cleanup
    type(finite_volume_scheme_t), intent(inout) :: self
    integer(ik) :: alloc_status

    alloc_status = 0

    call debug_print('Calling finite_volume_scheme_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%evolution_operator)) then
      deallocate(self%evolution_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%evolution_operator"
    end if

    if(allocated(self%grid)) then
      deallocate(self%grid, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%grid"
    end if

    if(allocated(self%reconstruction_operator)) then
      deallocate(self%reconstruction_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%reconstruction_operator"
    end if

    if(allocated(self%bc_plus_x)) then
      deallocate(self%bc_plus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_plus_x"
    end if

    if(allocated(self%bc_plus_y)) then
      deallocate(self%bc_plus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_plus_y"
    end if

    if(allocated(self%bc_minus_x)) then
      deallocate(self%bc_minus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_minus_x"
    end if

    if(allocated(self%bc_minus_y)) then
      deallocate(self%bc_minus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_minus_y"
    end if

    if(allocated(self%reconstructed_state)) then
      deallocate(self%reconstructed_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%reconstructed_state"
    end if

    if(allocated(self%evolved_corner_state)) then
      deallocate(self%evolved_corner_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%evolved_corner_state"
    end if

    if(allocated(self%corner_reference_state)) then
      deallocate(self%corner_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%corner_reference_state"
    end if

    if(allocated(self%evolved_downup_midpoints_state)) then
      deallocate(self%evolved_downup_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%evolved_downup_midpoints_state"
    end if

    if(allocated(self%evolved_leftright_midpoints_state)) then
      deallocate(self%evolved_leftright_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%evolved_leftright_midpoints_state"
    end if

    if(allocated(self%downup_midpoints_reference_state)) then
      deallocate(self%downup_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%downup_midpoints_reference_state"
    end if

    if(allocated(self%leftright_midpoints_reference_state)) then
      deallocate(self%leftright_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%leftright_midpoints_reference_state"
    end if

  end subroutine finalize

  subroutine apply_source_terms(self)
    class(finite_volume_scheme_t), intent(inout) :: self
  end subroutine

  subroutine reconstruct(self, conserved_vars, lbounds)
    !< Implementation of the FVLEG reconstruction.
    !< This reconstructs the entire grid at all the nodes/midpoints
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(in) :: conserved_vars

    call debug_print('Reconstructing the domain', __FILE__, __LINE__)
    call self%reconstruction_operator%set_conserved_vars_pointer(conserved_vars=conserved_vars, &
                                                                 lbounds=lbounds)
    call self%reconstruction_operator%set_grid_pointer(self%grid)
    call self%reconstruction_operator%reconstruct_domain(reconstructed_domain=self%reconstructed_state, &
                                                         lbounds=lbound(self%reconstructed_state))

  end subroutine reconstruct

  subroutine apply_conserved_vars_bc(self, conserved_vars, lbounds)
    !< Apply the boundary conditions
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: conserved_vars

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Calling apply_conserved_var_bc', __FILE__, __LINE__)

    max_priority_bc = max(self%bc_plus_x%priority, self%bc_plus_y%priority, &
                          self%bc_minus_x%priority, self%bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(self%bc_plus_x%priority == priority) then
        call self%bc_plus_x%apply_conserved_var_bc(conserved_vars=conserved_vars)
      end if

      if(self%bc_plus_y%priority == priority) then
        call self%bc_plus_y%apply_conserved_var_bc(conserved_vars=conserved_vars)
      end if

      if(self%bc_minus_x%priority == priority) then
        call self%bc_minus_x%apply_conserved_var_bc(conserved_vars=conserved_vars)
      end if

      if(self%bc_minus_y%priority == priority) then
        call self%bc_minus_y%apply_conserved_var_bc(conserved_vars=conserved_vars)
      end if

    end do

  end subroutine apply_conserved_vars_bc

  subroutine apply_reconstructed_state_bc(self)
    !< Apply the boundary conditions
    class(finite_volume_scheme_t), intent(inout) :: self

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Calling apply_reconstructed_state_bc', __FILE__, __LINE__)

    max_priority_bc = max(self%bc_plus_x%priority, self%bc_plus_y%priority, &
                          self%bc_minus_x%priority, self%bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(self%bc_plus_x%priority == priority) then
        call self%bc_plus_x%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
      end if

      if(self%bc_plus_y%priority == priority) then
        call self%bc_plus_y%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
      end if

      if(self%bc_minus_x%priority == priority) then
        call self%bc_minus_x%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
      end if

      if(self%bc_minus_y%priority == priority) then
        call self%bc_minus_y%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
      end if

    end do
  end subroutine apply_reconstructed_state_bc

  subroutine apply_cell_gradient_bc(self)
    !< Apply the boundary conditions
    class(finite_volume_scheme_t), intent(inout) :: self

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    max_priority_bc = max(self%bc_plus_x%priority, self%bc_plus_y%priority, &
                          self%bc_minus_x%priority, self%bc_minus_y%priority)

    if(self%reconstruction_operator%order > 1) then
      call debug_print('Calling apply_cell_gradient_bc', __FILE__, __LINE__)
      do priority = max_priority_bc, 0, -1

        if(self%bc_plus_x%priority == priority) then
          call self%bc_plus_x%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
        end if

        if(self%bc_plus_y%priority == priority) then
          call self%bc_plus_y%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
        end if

        if(self%bc_minus_x%priority == priority) then
          call self%bc_minus_x%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
        end if

        if(self%bc_minus_y%priority == priority) then
          call self%bc_minus_y%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
        end if

      end do
    end if

  end subroutine apply_cell_gradient_bc

  subroutine calculate_reference_state(self, conserved_vars, lbounds)
    !< Calculate the reference state at each corner/midpoint. This is just an average of
    !< the neighboring cells
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: conserved_vars

    integer(ik) :: i, j
    integer(ik) :: ilo, jlo, ihi, jhi

    ! left/right midpoints -> needs to average cells above and below
    ilo = lbound(self%leftright_midpoints_reference_state, dim=2)
    ihi = ubound(self%leftright_midpoints_reference_state, dim=2)
    jlo = lbound(self%leftright_midpoints_reference_state, dim=3)
    jhi = ubound(self%leftright_midpoints_reference_state, dim=3)
    do concurrent(j=jlo:jhi)
      do concurrent(i=ilo:ihi)
        associate(U_tilde=>self%leftright_midpoints_reference_state, U=>conserved_vars)
          ! U_tilde(:, i, j) = 0.5_rk * (U(:, i, j) + U(:, i, j - 1))
          U_tilde(:, i, j) = max(U(:, i, j), U(:, i, j - 1))
          ! debug_write(*,*) U(:, i, j), U(:, i, j - 1)
        end associate
      end do
    end do

    ! up/down midpoints -> needs to average cells right and left
    ilo = lbound(self%downup_midpoints_reference_state, dim=2)
    ihi = ubound(self%downup_midpoints_reference_state, dim=2)
    jlo = lbound(self%downup_midpoints_reference_state, dim=3)
    jhi = ubound(self%downup_midpoints_reference_state, dim=3)
    do concurrent(j=jlo:jhi)
      do concurrent(i=ilo:ihi)
        associate(U_tilde=>self%downup_midpoints_reference_state, U=>conserved_vars)
          ! U_tilde(:, i, j) = 0.5_rk * (U(:, i - 1, j) + U(:, i, j))
          U_tilde(:, i, j) = max(U(:, i - 1, j), U(:, i, j))
          ! debug_write(*,*) U(:, i - 1, j), U(:, i, j)
        end associate
      end do
    end do

    ! Corners
    ilo = lbound(self%corner_reference_state, dim=2)
    ihi = ubound(self%corner_reference_state, dim=2)
    jlo = lbound(self%corner_reference_state, dim=3)
    jhi = ubound(self%corner_reference_state, dim=3)
    do concurrent(j=jlo:jhi)
      do concurrent(i=ilo:ihi)
        associate(U_tilde=>self%corner_reference_state, U=>conserved_vars)
          ! U_tilde(:, i, j) = 0.25_rk * (U(:, i, j) + U(:, i - 1, j) + &
          !                               U(:, i, j - 1) + U(:, i - 1, j - 1))
          U_tilde(:, i, j) = max(U(:, i, j), U(:, i - 1, j), &
                                 U(:, i, j - 1), U(:, i - 1, j - 1))

          ! debug_write(*,*)
        end associate
      end do
    end do

  end subroutine calculate_reference_state

  subroutine evolve_domain(self)
    class(finite_volume_scheme_t), intent(inout) :: self

    call debug_print('Evolving the domain', __FILE__, __LINE__)

    ! Set the pointers for the evolution operator (makes the code w/in it easier)
    call self%evolution_operator%set_grid_pointer(grid_target=self%grid)
    call self%evolution_operator%set_reconstruction_operator_pointer(operator_target=self%reconstruction_operator)
    call self%evolution_operator%set_reconstructed_state_pointer(reconstructed_state_target=self%reconstructed_state, &
                                                                 lbounds=lbound(self%reconstructed_state))

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
    call self%evolution_operator%evolve_leftright_midpoints(reference_state=self%leftright_midpoints_reference_state, &
                                                            evolved_state=self%evolved_leftright_midpoints_state)

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
    call self%evolution_operator%evolve_downup_midpoints(reference_state=self%downup_midpoints_reference_state, &
                                                         evolved_state=self%evolved_downup_midpoints_state)

    ! Evolve, i.e. E0(R_omega), at all corner nodes
    call debug_print('Evolving corner nodes', __FILE__, __LINE__)
    call self%evolution_operator%evolve_corners(reference_state=self%corner_reference_state, &
                                                evolved_state=self%evolved_corner_state)

    ! call self%evolution_operator%nullify_pointer_members()
  end subroutine evolve_domain

end module mod_finite_volume_schemes
