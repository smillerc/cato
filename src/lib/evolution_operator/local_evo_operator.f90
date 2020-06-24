module mod_local_evo_operator
#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */
  use, intrinsic :: iso_fortran_env, only: &
    ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use mod_globals, only: debug_print, print_recon_data, track_single_cell_cone, enable_debug_cone_output, &
                         du_midpoint_unit, lr_midpoint_unit, corner_unit, single_cell_cone_unit

  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: near_zero, equal, neumaier_sum
  use math_constants, only: pi, rad2deg
  use mod_input, only: input_t
  use mod_mach_cone_collection, only: mach_cone_collection_t

  use mod_grid, only: grid_t, C1, M1, C2, M2, C3, M3, C4, M4
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_reconstruction_factory, only: reconstruction_factory

  implicit none

  private
  public :: local_evo_operator_t

  type :: local_evo_operator_t
    !< Local Evolution Operator (E0)
    !< This class provides the framework for any evolution operators. For example, the local evolution
    !< operator will inherit from this class. Most of the attributes are the same, just that the
    !< implementation varies between those who inherit this class

    class(grid_t), pointer :: grid => null()
    !< pointer to the grid object

    real(rk), dimension(:, :, :), pointer :: reconstructed_rho => null()
    !< (point 1:8, i, j); pointer to reconstructed density
    real(rk), dimension(:, :, :), pointer :: reconstructed_u => null()
    !< (point 1:8, i, j); pointer to reconstructed x-velocity
    real(rk), dimension(:, :, :), pointer :: reconstructed_v => null()
    !< (point 1:8, i, j); pointer to reconstructed y-velocity
    real(rk), dimension(:, :, :), pointer :: reconstructed_p => null()
    !< (point 1:8, i, j); pointer to reconstructed pressure

    class(abstract_reconstruction_t), pointer :: reconstruction_operator => null()
    !< pointer to the R_Omega operator used to provide values at the P' location

    character(:), allocatable :: name  !< Name of the evolution operator
    real(rk) :: tau !< time increment to evolve
    integer(ik) :: iteration = 0
    real(rk) :: time = 0.0_rk
    real(rk) :: time_step = 0.0_rk
  contains
    procedure, public :: initialize
    procedure, public :: evolve
    procedure, public, nopass :: e0_operator
    final :: finalize
  end type local_evo_operator_t

contains

  subroutine initialize(self, dt, grid_target, recon_operator_target)
    !< Constructor for the FVLEG operator
    class(local_evo_operator_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid_target
    real(rk), intent(in) :: dt !< timestep
    class(abstract_reconstruction_t), intent(in), target :: recon_operator_target

    self%name = "FVLEG"
    self%time_step = dt
    self%grid => grid_target
    self%reconstruction_operator => recon_operator_target

  end subroutine initialize

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self

    call debug_print('Running local_evo_operator_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)

    if(associated(self%reconstructed_rho)) nullify(self%reconstructed_rho)
    if(associated(self%reconstructed_u)) nullify(self%reconstructed_u)
    if(associated(self%reconstructed_v)) nullify(self%reconstructed_v)
    if(associated(self%reconstructed_p)) nullify(self%reconstructed_p)
    if(associated(self%reconstruction_operator)) nullify(self%reconstruction_operator)

    ! if (allocated(self%leftright_midpoint_neighbors)) deallocate(self%leftright_midpoint_neighbors)
    ! if (allocated(self%downup_midpoint_neighbors)) deallocate(self%downup_midpoint_neighbors)
    ! if (allocated(self%corner_neighbors)) deallocate(self%corner_neighbors)

    if(allocated(self%name)) deallocate(self%name)

  end subroutine finalize

  subroutine evolve(self, dt, location, evolved_rho, evolved_u, evolved_v, evolved_p, lbounds, error_code)
    !< Create Mach cones and evolve the state at all of the left/right midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the left and right.

    class(local_evo_operator_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    character(len=*), intent(in) :: location !< Mach cone location ['corner', 'left/right midpoint', or 'down/up midpoint']
    integer(ik), intent(out) :: error_code
    real(rk), intent(in) :: dt !< time step
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: evolved_rho
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: evolved_u
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: evolved_v
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: evolved_p

    ! Locals
    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi

    integer(ik), dimension(:, :, :, :), allocatable :: leftright_midpoint_neighbors
    !< ((i,j), (cell 1:2), i, j); (i,j) indices of the cells on either side of the left/right vector midpoint

    integer(ik), dimension(:, :, :, :), allocatable :: downup_midpoint_neighbors
    !< ((i,j), (cell 1:2), i, j); (i,j) indices of the cells on either side of the down/up vector midpoint

    integer(ik), dimension(:, :, :, :), allocatable :: corner_neighbors
    !< ((i,j), (cell 1:4), i, j); (i,j) indices of the cells surrounding the corner vector set

    type(mach_cone_collection_t) :: leftright_midpoint_mach_cones
    type(mach_cone_collection_t) :: downup_midpoint_mach_cones
    type(mach_cone_collection_t) :: corner_mach_cones

    real(rk) :: tau !< time increment

    real(rk), dimension(:, :, :), allocatable :: reconstructed_rho
    !< (point 1:8, i, j); the reconstructed density of the corners/midpoints with respect to each cell

    real(rk), dimension(:, :, :), allocatable :: reconstructed_u
    !< (point 1:8, i, j); the reconstructed x-velocity of the corners/midpoints with respect to each cell

    real(rk), dimension(:, :, :), allocatable :: reconstructed_v
    !< (point 1:8, i, j); the reconstructed y-velocity of the corners/midpoints with respect to each cell

    real(rk), dimension(:, :, :), allocatable :: reconstructed_p
    !< (point 1:8, i, j); the reconstructed pressure of the corners/midpoints with respect to each cell

    error_code = 0
    self%time_step = dt
    tau = dt / 2.0_rk
    ilo = lbound(evolved_rho, dim=1)
    ihi = ubound(evolved_rho, dim=1)
    jlo = lbound(evolved_rho, dim=2)
    jhi = ubound(evolved_rho, dim=2)

    select case(trim(location))
    case('corner')
      call debug_print('In local_evo_operator_t%evolve() -> location=corner', __FILE__, __LINE__)
      !       cell 4                   cell 3
      !      (i-1,j)                   (i,j)
      !  C4----M3----C3    P3   C4----M3----C3
      !  |            |    |    |            |
      !  M4    X4    M2    |    M4    X3    M2
      !  |            |    |    |            |
      !  C1----M1----C2    |    C1----M1----C2
      !                    |
      !  P4----------------O-----------------P2
      !                    |
      !  C4----M3----C3    |    C4----M3----C3
      !  |            |    |    |            |
      !  M4    X1    M2    |    M4    X2    M2
      !  |            |    |    |            |
      !  C1----M1----C2    P1   C1----M1----C2
      !      cell 1                  cell 2
      !     (i-1,j-1)               (i,j-1)

      ! For corners, the edge vectors go down (0-P1), right(O-P2), up(O-P3), right (O-P4)
      ! O is the corner point shared by the 4 cells
      if(.not. allocated(corner_neighbors)) then
        allocate(corner_neighbors(2, 4, ilo:ihi, jlo:jhi))
        do j = jlo, jhi
          do i = ilo, ihi
            corner_neighbors(:, 1, i, j) = [i - 1, j - 1] ! lower left
            corner_neighbors(:, 2, i, j) = [i, j - 1]     ! lower right
            corner_neighbors(:, 3, i, j) = [i, j]         ! upper right
            corner_neighbors(:, 4, i, j) = [i - 1, j]     ! upper left
          end do
        end do
      end if

      allocate(reconstructed_rho(4, ilo:ihi, jlo:jhi))
      allocate(reconstructed_u(4, ilo:ihi, jlo:jhi))
      allocate(reconstructed_v(4, ilo:ihi, jlo:jhi))
      allocate(reconstructed_p(4, ilo:ihi, jlo:jhi))

      !$omp parallel default(none) &
      !$omp shared(self) &
      !$omp shared(reconstructed_rho, reconstructed_u, reconstructed_v, reconstructed_p) &
      !$omp firstprivate(ilo,ihi,jlo,jhi) &
      !$omp private(i,j)
      !$omp do
      do j = jlo, jhi
        do i = ilo, ihi
          ! reconstructed_state indexing; ((rho, u ,v, p), point, node/midpoint, i, j)

          ! Cell 1: lower left cell -> corner is in the upper right (C3) of its parent cell
          reconstructed_rho(1, i, j) = self%reconstructed_rho(C3, i - 1, j - 1)
          reconstructed_u(1, i, j) = self%reconstructed_u(C3, i - 1, j - 1)
          reconstructed_v(1, i, j) = self%reconstructed_v(C3, i - 1, j - 1)
          reconstructed_p(1, i, j) = self%reconstructed_p(C3, i - 1, j - 1)

          ! Cell 2: lower right cell -> corner is in the upper left (C4) of its parent cell
          reconstructed_rho(2, i, j) = self%reconstructed_rho(C4, i, j - 1)
          reconstructed_u(2, i, j) = self%reconstructed_u(C4, i, j - 1)
          reconstructed_v(2, i, j) = self%reconstructed_v(C4, i, j - 1)
          reconstructed_p(2, i, j) = self%reconstructed_p(C4, i, j - 1)

          ! Cell 3: upper right cell-> corner is in the lower left (C1) of its parent cell
          reconstructed_rho(3, i, j) = self%reconstructed_rho(C1, i, j)
          reconstructed_u(3, i, j) = self%reconstructed_u(C1, i, j)
          reconstructed_v(3, i, j) = self%reconstructed_v(C1, i, j)
          reconstructed_p(3, i, j) = self%reconstructed_p(C1, i, j)

          ! Cell 4: upper left cell -> corner is in the lower right (C2) of its parent cell
          reconstructed_rho(4, i, j) = self%reconstructed_rho(C2, i - 1, j)
          reconstructed_u(4, i, j) = self%reconstructed_u(C2, i - 1, j)
          reconstructed_v(4, i, j) = self%reconstructed_v(C2, i - 1, j)
          reconstructed_p(4, i, j) = self%reconstructed_p(C2, i - 1, j)

        end do
      end do
      !$omp end do
      !$omp end parallel

      call corner_mach_cones%initialize(tau=tau, edge_vectors=self%grid%corner_edge_vectors, &
                                        reconstructed_rho=reconstructed_rho, &
                                        reconstructed_u=reconstructed_u, &
                                        reconstructed_v=reconstructed_v, &
                                        reconstructed_p=reconstructed_p, &
                                        cell_indices=corner_neighbors, &
                                        cone_location=trim(location), &
                                        lbounds=lbound(reconstructed_rho), &
                                        reconstruction_operator=self%reconstruction_operator)
      call self%e0_operator(cones=corner_mach_cones, &
                            rho=evolved_rho, u=evolved_u, v=evolved_v, p=evolved_p)

    case('left/right midpoint')
      call debug_print('In local_evo_operator_t%evolve() -> location=left/right midpoint', __FILE__, __LINE__)
      !       cell 1
      !      (i-1,j)
      !  C4----M3----C3
      !  |            |
      !  M4    X1     M2
      !  |            |
      !  C1----M1----C2
      !
      !  P1----O----P2  edge vector
      !
      !  C4----M3----C3
      !  |           |
      !  M4    X2    M2
      !  |           |
      !  C1----M1----C2
      !       cell 2
      !      (i,j-1)
      !
      ! For left/right midpoints, the edge vectors go left (O-P1) then right (O-P2).
      ! The neighboring cells are up (i, j+1) and down (i, j-1)

      if(.not. allocated(leftright_midpoint_neighbors)) then
        allocate(leftright_midpoint_neighbors(2, 2, ilo:ihi, jlo:jhi))
        do j = jlo, jhi
          do i = ilo, ihi
            leftright_midpoint_neighbors(:, 1, i, j) = [i, j]     ! above
            leftright_midpoint_neighbors(:, 2, i, j) = [i, j - 1] ! below
          end do
        end do
      end if

      allocate(reconstructed_rho(2, ilo:ihi, jlo:jhi))
      allocate(reconstructed_u(2, ilo:ihi, jlo:jhi))
      allocate(reconstructed_v(2, ilo:ihi, jlo:jhi))
      allocate(reconstructed_p(2, ilo:ihi, jlo:jhi))

      !$omp parallel default(none) &
      !$omp shared(self) &
      !$omp shared(reconstructed_rho, reconstructed_u, reconstructed_v, reconstructed_p) &
      !$omp firstprivate(ilo,ihi,jlo,jhi) &
      !$omp private(i,j)
      !$omp do
      do j = jlo, jhi
        do i = ilo, ihi
          ! Cell X1 -> use M1 from the cell above
          reconstructed_rho(1, i, j) = self%reconstructed_rho(M1, i, j)
          reconstructed_u(1, i, j) = self%reconstructed_u(M1, i, j)
          reconstructed_v(1, i, j) = self%reconstructed_v(M1, i, j)
          reconstructed_p(1, i, j) = self%reconstructed_p(M1, i, j)

          ! Cell X2 -> use M3 from the cell below
          reconstructed_rho(2, i, j) = self%reconstructed_rho(M3, i, j - 1)
          reconstructed_u(2, i, j) = self%reconstructed_u(M3, i, j - 1)
          reconstructed_v(2, i, j) = self%reconstructed_v(M3, i, j - 1)
          reconstructed_p(2, i, j) = self%reconstructed_p(M3, i, j - 1)
        end do
      end do
      !$omp end do
      !$omp end parallel

      call leftright_midpoint_mach_cones%initialize(tau=tau, edge_vectors=self%grid%leftright_midpoint_edge_vectors, &
                                                    reconstructed_rho=reconstructed_rho, &
                                                    reconstructed_u=reconstructed_u, &
                                                    reconstructed_v=reconstructed_v, &
                                                    reconstructed_p=reconstructed_p, &
                                                    cell_indices=leftright_midpoint_neighbors, &
                                                    cone_location=trim(location), &
                                                    lbounds=lbound(reconstructed_rho), &
                                                    reconstruction_operator=self%reconstruction_operator)
      call self%e0_operator(cones=leftright_midpoint_mach_cones, &
                            rho=evolved_rho, u=evolved_u, v=evolved_v, p=evolved_p)

    case('down/up midpoint')
      call debug_print('In local_evo_operator_t%evolve() -> location=down/up midpoint', __FILE__, __LINE__)
      !       cell 1          edge          cell 2
      !      (i-1,j)         vector          (i,j)
      !  C4----M3----C3       P2       C4----M3----C3
      !  |           |        |       |            |
      !  M4    X1    M2       O       M4    X2     M2
      !  |           |        |       |            |
      !  C1----M1----C2      P1       C1----M1----C2
      !
      ! For down/up midpoints, the edge vectors go down (O-P1) then up (O-P2).
      ! The neighboring cells are left (i, j) and right (i-1, j)
      if(.not. allocated(downup_midpoint_neighbors)) then
        allocate(downup_midpoint_neighbors(2, 2, ilo:ihi, jlo:jhi))
        do j = jlo, jhi
          do i = ilo, ihi
            downup_midpoint_neighbors(:, 1, i, j) = [i - 1, j] ! left
            downup_midpoint_neighbors(:, 2, i, j) = [i, j]     ! right
          end do
        end do
      end if

      allocate(reconstructed_rho(2, ilo:ihi, jlo:jhi))
      allocate(reconstructed_u(2, ilo:ihi, jlo:jhi))
      allocate(reconstructed_v(2, ilo:ihi, jlo:jhi))
      allocate(reconstructed_p(2, ilo:ihi, jlo:jhi))

      !$omp parallel default(none) &
      !$omp shared(self) &
      !$omp shared(reconstructed_rho, reconstructed_u, reconstructed_v, reconstructed_p) &
      !$omp firstprivate(ilo,ihi,jlo,jhi) &
      !$omp private(i,j)
      !$omp do
      do j = jlo, jhi
        do i = ilo, ihi
          ! Cell X1: cell to the left -> use M2
          reconstructed_rho(1, i, j) = self%reconstructed_rho(M2, i - 1, j)
          reconstructed_u(1, i, j) = self%reconstructed_u(M2, i - 1, j)
          reconstructed_v(1, i, j) = self%reconstructed_v(M2, i - 1, j)
          reconstructed_p(1, i, j) = self%reconstructed_p(M2, i - 1, j)

          ! Cell X2: cell to the right -> use M4
          reconstructed_rho(2, i, j) = self%reconstructed_rho(M4, i, j)
          reconstructed_u(2, i, j) = self%reconstructed_u(M4, i, j)
          reconstructed_v(2, i, j) = self%reconstructed_v(M4, i, j)
          reconstructed_p(2, i, j) = self%reconstructed_p(M4, i, j)
        end do
      end do
      !$omp end do
      !$omp end parallel

      call downup_midpoint_mach_cones%initialize(tau=tau, edge_vectors=self%grid%downup_midpoint_edge_vectors, &
                                                 reconstructed_rho=reconstructed_rho, &
                                                 reconstructed_u=reconstructed_u, &
                                                 reconstructed_v=reconstructed_v, &
                                                 reconstructed_p=reconstructed_p, &
                                                 cell_indices=downup_midpoint_neighbors, &
                                                 cone_location=trim(location), &
                                                 lbounds=lbound(reconstructed_rho), &
                                                 reconstruction_operator=self%reconstruction_operator)
      call self%e0_operator(cones=downup_midpoint_mach_cones, &
                            rho=evolved_rho, u=evolved_u, v=evolved_v, p=evolved_p)

    case default
      error stop "Unsupported location in local_evo_operator_t%evolve()"
    end select

    if(allocated(reconstructed_rho)) deallocate(reconstructed_rho)
    if(allocated(reconstructed_u)) deallocate(reconstructed_u)
    if(allocated(reconstructed_v)) deallocate(reconstructed_v)
    if(allocated(reconstructed_p)) deallocate(reconstructed_p)
  end subroutine evolve

  subroutine e0_operator(cones, rho, u, v, p)
    type(mach_cone_collection_t), intent(in) :: cones
    real(rk), dimension(:, :), intent(inout) :: rho
    real(rk), dimension(:, :), intent(inout) :: u
    real(rk), dimension(:, :), intent(inout) :: v
    real(rk), dimension(:, :), intent(inout) :: p

    real(rk), allocatable, dimension(:, :) :: rho_a_tilde
    integer(ik) :: i, j, c
    real(rk) :: t1, t2, t3
    real(rk), dimension(3) :: t_tot

    call debug_print('Running local_evo_operator_t%e0_operator()', __FILE__, __LINE__)

    allocate(rho_a_tilde, mold=cones%reference_density)

    associate(a_tilde=>cones%reference_sound_speed, &
              dtheta=>cones%dtheta, &
              sin_dtheta=>cones%sin_dtheta, &
              cos_dtheta=>cones%cos_dtheta, &
              sin_d2theta=>cones%sin_d2theta, &
              cos_d2theta=>cones%cos_d2theta, &
              rho_p_prime=>cones%density_p_prime, &
              pressure_p_prime=>cones%pressure_p_prime, &
              u_i=>cones%u, &
              v_i=>cones%v, &
              p_i=>cones%p)

      !$omp parallel default(none) &
      !$omp shared(cones, p, rho, u, v, rho_a_tilde) &
      !$omp private(i,j) &
      !$omp private(t1, t2, t3, t_tot)
      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          rho_a_tilde(i, j) = cones%reference_density(i, j) * cones%reference_sound_speed(i, j)
        end do
      end do
      !$omp end do

      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          t1 = neumaier_sum(p_i(:, i, j) * dtheta(:, i, j))
          t2 = neumaier_sum(-rho_a_tilde(i, j) * u_i(:, i, j) * sin_dtheta(:, i, j))
          t3 = neumaier_sum(rho_a_tilde(i, j) * v_i(:, i, j) * cos_dtheta(:, i, j))
          t_tot = [t1, t2, t3]
          p(i, j) = neumaier_sum(t_tot) / (2.0_rk * pi)
        end do
      end do
      !$omp end do

      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          t_tot = [rho_p_prime(i, j),(p(i, j) / a_tilde(i, j)**2), -pressure_p_prime(i, j) / a_tilde(i, j)**2]
          rho(i, j) = neumaier_sum(t_tot)
        end do
      end do
      !$omp end do

      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          t1 = neumaier_sum(-p_i(:, i, j) * sin_dtheta(:, i, j)) / rho_a_tilde(i, j)
          if(abs(t1) < 1e-5_rk * a_tilde(i, j)) t1 = 0.0_rk

          t2 = neumaier_sum(u_i(:, i, j) * ((dtheta(:, i, j) / 2.0_rk) + (sin_d2theta(:, i, j) / 4.0_rk)))
          if(abs(t2) < 1e-5_rk * a_tilde(i, j)) t2 = 0.0_rk

          t3 = neumaier_sum(-v_i(:, i, j) * cos_d2theta(:, i, j)) / 4.0_rk
          if(abs(t3) < 1e-5_rk * a_tilde(i, j)) t3 = 0.0_rk

          t_tot = [t1, t2, t3]
          u(i, j) = neumaier_sum(t_tot) / pi
        end do
      end do
      !$omp end do

      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          t1 = neumaier_sum(p_i(:, i, j) * cos_dtheta(:, i, j)) / rho_a_tilde(i, j)
          if(abs(t1) < 1e-5_rk * a_tilde(i, j)) t1 = 0.0_rk

          t2 = neumaier_sum(-u_i(:, i, j) * cos_d2theta(:, i, j)) / 4.0_rk
          if(abs(t2) < 1e-5_rk * a_tilde(i, j)) t2 = 0.0_rk

          t3 = neumaier_sum(v_i(:, i, j) * ((dtheta(:, i, j) / 2.0_rk) - (sin_d2theta(:, i, j) / 4.0_rk)))
          if(abs(t3) < 1e-5_rk * a_tilde(i, j)) t3 = 0.0_rk

          t_tot = [t1, t2, t3]
          v(i, j) = neumaier_sum(t_tot) / pi
        end do
      end do
      !$omp end do
      !$omp end parallel
    end associate
    deallocate(rho_a_tilde)
  end subroutine e0_operator

end module mod_local_evo_operator
