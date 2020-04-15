module mod_local_evo_operator
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, &
                                                                              std_err => error_unit, &
                                                                              std_out => output_unit
  use mod_globals, only: debug_print, print_recon_data

  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: near_zero, equal
  use math_constants, only: pi, rad2deg
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
  ! use mod_corner_mach_cone, only: corner_mach_cone_t, new_corner_cone
  ! use mod_midpoint_mach_cone, only: midpoint_mach_cone_t, new_midpoint_cone
  use mod_mach_cone_collection, only: mach_cone_collection_t
  use mod_grid, only: grid_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_reconstruction_factory, only: reconstruction_factory

  implicit none

  private
  public :: local_evo_operator_t

  type, extends(abstract_evo_operator_t) :: local_evo_operator_t
    !< Local Evolution Operator (E0)
  contains
    procedure, public :: initialize
    procedure, public :: evolve
    procedure, public, nopass :: e0_operator
    ! procedure, public :: evolve_leftright_midpoints
    ! procedure, public :: evolve_downup_midpoints
    ! procedure, public :: evolve_corners
    ! procedure, public :: e0_operator_corner
    ! procedure, public :: e0_operator_midpoint
    procedure, public :: copy
    final :: finalize
  end type local_evo_operator_t

contains

  subroutine initialize(self, input, grid_target, recon_operator_target)
    !< Constructor for the FVLEG operator
    class(local_evo_operator_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    class(grid_t), intent(in), target :: grid_target
    class(abstract_reconstruction_t), intent(in), target :: recon_operator_target

    self%name = "FVLEG"
    call self%set_tau(input%tau)
    self%grid => grid_target
    self%reconstruction_operator => recon_operator_target

  end subroutine

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self

    call debug_print('Running local_evo_operator_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%name)) deallocate(self%name)
    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%reconstructed_state)) nullify(self%reconstructed_state)
    if(associated(self%reconstruction_operator)) nullify(self%reconstruction_operator)
  end subroutine finalize

  subroutine copy(out_evo, in_evo)
    class(abstract_evo_operator_t), intent(in) :: in_evo
    class(local_evo_operator_t), intent(inout) :: out_evo

    call debug_print('Running local_evo_operator_t%copy()', __FILE__, __LINE__)

    if(associated(out_evo%grid)) nullify(out_evo%grid)
    out_evo%grid => in_evo%grid

    if(associated(out_evo%reconstructed_state)) nullify(out_evo%reconstructed_state)
    out_evo%reconstructed_state => in_evo%reconstructed_state

    if(associated(out_evo%reconstruction_operator)) nullify(out_evo%reconstruction_operator)
    out_evo%reconstruction_operator => in_evo%reconstruction_operator

    if(allocated(out_evo%name)) deallocate(out_evo%name)
    allocate(out_evo%name, source=in_evo%name)

  end subroutine

  subroutine evolve(self, location, evolved_state, lbounds, error_code)
    !< Create Mach cones and evolve the state at all of the left/right midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the left and right.

    class(local_evo_operator_t), intent(inout) :: self

    integer(ik), dimension(3), intent(in) :: lbounds
    character(len=*), intent(in) :: location !< Mach cone location ['corner', 'left/right midpoint', or 'down/up midpoint']
    integer(ik), intent(out) :: error_code
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(out) :: evolved_state

    ! Locals
    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik), parameter :: N = 1 !< used to select the corners in the reconstructed_state
    integer(ik), parameter :: M = 2 !< used to select the midpoint in the reconstructed_state

    integer(ik), dimension(:, :, :, :), allocatable :: neighbor_cell_indices
    !< ((i,j), cell_id, i, j); neighbor cell (i,j) sets for each location

    real(rk), dimension(:, :, :, :), allocatable :: reconstructed_state
    !< ((rho, u, v, p), cell_id, i, j); the reconstructed state of the corner with respect to each cell

    error_code = 0

    ilo = lbound(evolved_state, dim=2)
    ihi = ubound(evolved_state, dim=2)
    jlo = lbound(evolved_state, dim=3)
    jhi = ubound(evolved_state, dim=3)

    select case(trim(location))
    case('corner')
      call debug_print('In local_evo_operator_t%evolve() -> location=corner', __FILE__, __LINE__)
      !       cell 4                   cell 3
      !      (i-1,j)                   (i,j)
      !  N4----M3----N3    P3   N4----M3----N3
      !  |            |    |    |            |
      !  M4    C4    M2    |    M4    C3    M2
      !  |            |    |    |            |
      !  N1----M1----N2    |    N1----M1----N2
      !                    |
      !  P4----------------O-----------------P2
      !                    |
      !  N4----M3----N3    |    N4----M3----N3
      !  |            |    |    |            |
      !  M4    C1    M2    |    M4    C2    M2
      !  |            |    |    |            |
      !  N1----M1----N2    P1   N1----M1----N2
      !      cell 1                  cell 2
      !     (i-1,j-1)               (i,j-1)

      ! For corners, the edge vectors go down (0-P1), right(O-P2), up(O-P3), right (O-P4)
      ! O is the corner point shared by the 4 cells

      allocate(reconstructed_state(4, 4, ilo:ihi, jlo:jhi))
      allocate(neighbor_cell_indices(2, 4, ilo:ihi, jlo:jhi))

      do j = jlo, jhi
        do i = ilo, ihi
          neighbor_cell_indices(:, 1, i, j) = [i - 1, j - 1] ! lower left
          neighbor_cell_indices(:, 2, i, j) = [i, j - 1]     ! lower right
          neighbor_cell_indices(:, 3, i, j) = [i, j]         ! upper right
          neighbor_cell_indices(:, 4, i, j) = [i - 1, j]     ! upper left
        end do
      end do

      ! if(.not. allocated(self%corner_mach_cones)) then
      !   allocate(self%corner_mach_cones)
      ! end if

      do j = jlo, jhi
        do i = ilo, ihi
          ! reconstructed_state indexing; ((rho, u ,v, p), point, node/midpoint, i, j)

          ! Cell 1: lower left cell -> corner is in the upper right (N3) of its parent cell
          reconstructed_state(:, 1, i, j) = self%reconstructed_state(:, 3, N, i - 1, j - 1)

          ! Cell 2: lower right cell -> corner is in the upper left (N4) of its parent cell
          reconstructed_state(:, 2, i, j) = self%reconstructed_state(:, 4, N, i, j - 1)

          ! Cell 3: upper right cell-> corner is in the lower left (N1) of its parent cell
          reconstructed_state(:, 3, i, j) = self%reconstructed_state(:, 1, N, i, j)

          ! Cell 4: upper left cell -> corner is in the lower right (N2) of its parent cell
          reconstructed_state(:, 4, i, j) = self%reconstructed_state(:, 2, N, i - 1, j)

        end do
      end do

      call self%corner_mach_cones%initialize(edge_vectors=self%grid%corner_edge_vectors, &
                                             vector_scaling=self%grid%corner_edge_vectors_scale, &
                                             reconstructed_state=reconstructed_state, &
                                             cell_indices=neighbor_cell_indices, &
                                             cone_location=trim(location))
      call self%e0_operator(cones=self%corner_mach_cones, primitive_vars=evolved_state)

    case('left/right midpoint')
      call debug_print('In local_evo_operator_t%evolve() -> location=left/right midpoint', __FILE__, __LINE__)
      !       cell 1
      !      (i-1,j)
      !  N4----M3----N3
      !  |            |
      !  M4    C1     M2
      !  |            |
      !  N1----M1----N2
      !
      !  P1----O----P2  edge vector
      !
      !  N4----M3----N3
      !  |           |
      !  M4    C2    M2
      !  |           |
      !  N1----M1----N2
      !       cell 2
      !      (i,j-1)
      !
      ! For left/right midpoints, the edge vectors go left (O-P1) then right (O-P2).
      ! The neighboring cells are up (i, j+1) and down (i, j-1)

      allocate(reconstructed_state(4, 2, ilo:ihi, jlo:jhi))
      allocate(neighbor_cell_indices(2, 2, ilo:ihi, jlo:jhi))

      do j = jlo, jhi
        do i = ilo, ihi
          neighbor_cell_indices(:, 1, i, j) = [i, j]     ! above
          neighbor_cell_indices(:, 2, i, j) = [i, j - 1] ! below
        end do
      end do

      do j = jlo, jhi
        do i = ilo, ihi
          ! reconstructed_state indexing; ((rho, u ,v, p), point, node/midpoint, i, j)

          ! Cell 1 -> use M1 from the cell above
          reconstructed_state(:, 1, i, j) = self%reconstructed_state(:, 1, M, i, j)

          ! Cell 2 -> use M3 from the cell below
          reconstructed_state(:, 2, i, j) = self%reconstructed_state(:, 3, M, i, j - 1)
        end do
      end do

      ! if(.not. allocated(self%leftright_midpoint_mach_cones)) then
      !   allocate(self%leftright_midpoint_mach_cones)
      ! end if

      call self%leftright_midpoint_mach_cones%initialize(edge_vectors=self%grid%leftright_midpoint_edge_vectors, &
                                                         vector_scaling=self%grid%leftright_midpoint_edge_vectors_scale, &
                                                         reconstructed_state=reconstructed_state, &
                                                         cell_indices=neighbor_cell_indices, &
                                                         cone_location=trim(location))
      call self%e0_operator(cones=self%leftright_midpoint_mach_cones, primitive_vars=evolved_state)
    case('down/up midpoint')
      call debug_print('In local_evo_operator_t%evolve() -> location=down/up midpoint', __FILE__, __LINE__)
      !       cell 1          edge          cell 2
      !      (i-1,j)         vector          (i,j)
      !  N4----M3----N3       P2       N4----M3----N3
      !  |           |        |       |            |
      !  M4    C1    M2       O       M4    C2     M2
      !  |           |        |       |            |
      !  N1----M1----N2      P1       N1----M1----N2
      !
      ! For down/up midpoints, the edge vectors go down (O-P1) then up (O-P2).
      ! The neighboring cells are left (i, j) and right (i-1, j)

      allocate(reconstructed_state(4, 2, ilo:ihi, jlo:jhi))
      allocate(neighbor_cell_indices(2, 2, ilo:ihi, jlo:jhi))

      do j = jlo, jhi
        do i = ilo, ihi
          neighbor_cell_indices(:, 1, i, j) = [i - 1, j] ! left
          neighbor_cell_indices(:, 2, i, j) = [i, j]     ! right
        end do
      end do

      do j = jlo, jhi
        do i = ilo, ihi
          ! reconstructed_state indexing; ((rho, u ,v, p), point, node/midpoint, i, j)

          ! Cell 1: use M2 from the left cell
          reconstructed_state(:, 1, i, j) = self%reconstructed_state(:, 2, M, i - 1, j)

          ! Cell 2: cell to the right -> use M4 from the right cell
          reconstructed_state(:, 2, i, j) = self%reconstructed_state(:, 4, M, i, j)
        end do
      end do

      ! if(.not. allocated(self%downup_midpoint_mach_cones)) then
      !   allocate(self%downup_midpoint_mach_cones)
      ! end if

      call self%downup_midpoint_mach_cones%initialize(edge_vectors=self%grid%downup_midpoint_edge_vectors, &
                                                      vector_scaling=self%grid%downup_midpoint_edge_vectors_scale, &
                                                      reconstructed_state=reconstructed_state, &
                                                      cell_indices=neighbor_cell_indices, &
                                                      cone_location=trim(location))
      call self%e0_operator(cones=self%downup_midpoint_mach_cones, primitive_vars=evolved_state)
    case default
      error stop "Unsupported location in local_evo_operator_t%evolve()"
    end select

    if(allocated(neighbor_cell_indices)) deallocate(neighbor_cell_indices)
    if(allocated(reconstructed_state)) deallocate(reconstructed_state)
  end subroutine

  subroutine e0_operator(cones, primitive_vars)
    class(mach_cone_collection_t), intent(in) :: cones
    real(rk), dimension(:, :, :), intent(inout) :: primitive_vars

    associate(a_tilde=>cones%reference_sound_speed, &
              rho_a_tilde_sq=>cones%reference_density * cones%reference_sound_speed**2, &
              dtheta=>cones%dtheta, &
              sin_dtheta=>cones%sin_dtheta, &
              cos_dtheta=>cones%cos_dtheta, &
              sin_d2theta=>cones%sin_d2theta, &
              cos_d2theta=>cones%cos_d2theta, &
              rho_p_prime=>cones%density_p_prime, &
              pressure_p_prime=>cones%pressure_p_prime, &
              density=>primitive_vars(1, :, :), &
              u=>primitive_vars(2, :, :), &
              v=>primitive_vars(3, :, :), &
              pressure=>primitive_vars(4, :, :), &
              u_i_star=>cones%u, &
              v_i_star=>cones%v, &
              p_i_star=>cones%p)

      ! Note: These are all the normalized versions, e.g. p_i*, u_i*, of the equations
      pressure = sum(p_i_star * dtheta - u_i_star * sin_dtheta - v_i_star * cos_dtheta, dim=1) * (rho_a_tilde_sq / (2.0_rk * pi))

      ! Note: the a_tilde**2 term is not combined on purpose, so that floating point subtraction is more accurate
      density = rho_p_prime + (pressure / a_tilde**2) - (pressure_p_prime / a_tilde**2)

      u = (a_tilde / pi) * sum(-p_i_star * sin_dtheta &
                               + u_i_star * ((dtheta / 2.0_rk) + (sin_d2theta / 4.0_rk)) &
                               - v_i_star * (cos_d2theta / 4.0_rk), dim=1)

      v = (a_tilde / pi) * sum(p_i_star * cos_dtheta &
                               - u_i_star * (cos_d2theta / 4.0_rk) &
                               + v_i_star * ((dtheta / 2.0_rk) - (sin_d2theta / 4.0_rk)), dim=1) / pi
    end associate

  end subroutine

end module mod_local_evo_operator
