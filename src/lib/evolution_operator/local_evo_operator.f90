module mod_local_evo_operator
#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */
  use, intrinsic :: iso_fortran_env, only: &
    ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use mod_globals, only: debug_print, print_recon_data, track_single_cell_cone, enable_debug_cone_output, &
                         du_midpoint_unit, lr_midpoint_unit, corner_unit, single_cell_cone_unit

  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: near_zero, equal
  use math_constants, only: pi, rad2deg
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
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

    real(rk) :: tau !< time increment
    real(rk), dimension(:, :, :, :), allocatable :: reconstructed_state
    !< ((rho, u, v, p), cell_id, i, j); the reconstructed state of the corner with respect to each cell

    error_code = 0
    tau = self%time_step / 2.0_rk
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
      if(.not. allocated(self%corner_neighbors)) then
        allocate(self%corner_neighbors(2, 4, ilo:ihi, jlo:jhi))
        do j = jlo, jhi
          do i = ilo, ihi
            self%corner_neighbors(:, 1, i, j) = [i - 1, j - 1] ! lower left
            self%corner_neighbors(:, 2, i, j) = [i, j - 1]     ! lower right
            self%corner_neighbors(:, 3, i, j) = [i, j]         ! upper right
            self%corner_neighbors(:, 4, i, j) = [i - 1, j]     ! upper left
          end do
        end do
      end if

      allocate(reconstructed_state(4, 4, ilo:ihi, jlo:jhi))
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

      call self%corner_mach_cones%initialize(tau=tau, edge_vectors=self%grid%corner_edge_vectors, &
                                             reconstructed_state=reconstructed_state, &
                                             cell_indices=self%corner_neighbors, &
                                             cone_location=trim(location))
      call self%e0_operator(cones=self%corner_mach_cones, primitive_vars=evolved_state)

      ! !$omp parallel private(corner_cone, evolved_primitive_vars) shared(tau, reconstructed_state)
      ! !$omp do
      ! do j = jlo, jhi
      !   do i = ilo, ihi
      !     corner_cone = new_corner_cone(tau=tau, &
      !                                   edge_vectors=self%grid%corner_edge_vectors(:, :, i, j), &
      !                                   reconstructed_state=reconstructed_state(:, :, i, j), &
      !                                   cell_indices=self%corner_neighbors(:, :, i, j), &
      !                                   cone_location='corner')
      !     call self%e0_operator_corner(corner_cone, evolved_primitive_vars, 'evolve_corners', error_code)
      !     evolved_state(:, i, j) = evolved_primitive_vars
      !   end do
      ! end do
      ! !$omp end do
      ! !$omp end parallel

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

      if(.not. allocated(self%leftright_midpoint_neighbors)) then
        allocate(self%leftright_midpoint_neighbors(2, 2, ilo:ihi, jlo:jhi))
        do j = jlo, jhi
          do i = ilo, ihi
            self%leftright_midpoint_neighbors(:, 1, i, j) = [i, j]     ! above
            self%leftright_midpoint_neighbors(:, 2, i, j) = [i, j - 1] ! below
          end do
        end do
      end if

      allocate(reconstructed_state(4, 2, ilo:ihi, jlo:jhi))
      do j = jlo, jhi
        do i = ilo, ihi
          ! reconstructed_state indexing; ((rho, u ,v, p), point, node/midpoint, i, j)

          ! Cell 1 -> use M1 from the cell above
          reconstructed_state(:, 1, i, j) = self%reconstructed_state(:, 1, M, i, j)

          ! Cell 2 -> use M3 from the cell below
          reconstructed_state(:, 2, i, j) = self%reconstructed_state(:, 3, M, i, j - 1)
        end do
      end do

      call self%leftright_midpoint_mach_cones%initialize(tau=tau, edge_vectors=self%grid%leftright_midpoint_edge_vectors, &
                                                         reconstructed_state=reconstructed_state, &
                                                         cell_indices=self%leftright_midpoint_neighbors, &
                                                         cone_location=trim(location))
      call self%e0_operator(cones=self%leftright_midpoint_mach_cones, primitive_vars=evolved_state)
      ! !$omp parallel private(midpoint_cone, evolved_primitive_vars) shared(tau, reconstructed_state)
      ! !$omp do
      ! do j = jlo, jhi
      !   do i = ilo, ihi
      !     midpoint_cone = new_midpoint_cone(tau=tau, &
      !                                       edge_vectors=self%grid%leftright_midpoint_edge_vectors(:, :, i, j), &
      !                                       reconstructed_state=reconstructed_state(:, :, i, j), &
      !                                       cell_indices=self%leftright_midpoint_neighbors(:, :, i, j), &
      !                                       cone_location='left/right midpoint')
      !     call self%e0_operator_midpoint(midpoint_cone, evolved_primitive_vars, 'evolve_leftright_midpoints', error_code)
      !     evolved_state(:, i, j) = evolved_primitive_vars
      !   end do
      ! end do
      ! !$omp end do
      ! !$omp end parallel

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
      if(.not. allocated(self%downup_midpoint_neighbors)) then
        allocate(self%downup_midpoint_neighbors(2, 2, ilo:ihi, jlo:jhi))
        do j = jlo, jhi
          do i = ilo, ihi
            self%downup_midpoint_neighbors(:, 1, i, j) = [i - 1, j] ! left
            self%downup_midpoint_neighbors(:, 2, i, j) = [i, j]     ! right
          end do
        end do
      end if

      allocate(reconstructed_state(4, 2, ilo:ihi, jlo:jhi))
      do j = jlo, jhi
        do i = ilo, ihi
          ! reconstructed_state indexing; ((rho, u ,v, p), point, node/midpoint, i, j)

          ! Cell 1: use M2 from the left cell
          reconstructed_state(:, 1, i, j) = self%reconstructed_state(:, 2, M, i - 1, j)

          ! Cell 2: cell to the right -> use M4 from the right cell
          reconstructed_state(:, 2, i, j) = self%reconstructed_state(:, 4, M, i, j)
        end do
      end do

      call self%downup_midpoint_mach_cones%initialize(tau=tau, edge_vectors=self%grid%downup_midpoint_edge_vectors, &
                                                      reconstructed_state=reconstructed_state, &
                                                      cell_indices=self%downup_midpoint_neighbors, &
                                                      cone_location=trim(location))
      call self%e0_operator(cones=self%downup_midpoint_mach_cones, primitive_vars=evolved_state)

      ! !$omp parallel private(midpoint_cone, evolved_primitive_vars) shared(tau, reconstructed_state)
      ! !$omp do
      ! do j = jlo, jhi
      !   do i = ilo, ihi
      !     midpoint_cone = new_midpoint_cone(tau=tau, &
      !                                       edge_vectors=self%grid%downup_midpoint_edge_vectors(:, :, i, j), &
      !                                       reconstructed_state=reconstructed_state(:, :, i, j), &
      !                                       cell_indices=self%downup_midpoint_neighbors(:, :, i, j), &
      !                                       cone_location='down/up midpoint')
      !     call self%e0_operator_midpoint(midpoint_cone, evolved_primitive_vars, 'evolve_downup_midpoints', error_code)
      !     evolved_state(:, i, j) = evolved_primitive_vars
      !   end do
      ! end do
      ! !$omp end do
      ! !$omp end parallel

    case default
      error stop "Unsupported location in local_evo_operator_t%evolve()"
    end select

    if(allocated(reconstructed_state)) deallocate(reconstructed_state)
  end subroutine

  subroutine e0_operator(cones, primitive_vars)
    class(mach_cone_collection_t), intent(in) :: cones
    real(rk), dimension(:, :, :), intent(out) :: primitive_vars

    real(rk), allocatable, dimension(:, :) :: rho_a_tilde
    real(rk), allocatable, dimension(:, :) :: density
    real(rk), allocatable, dimension(:, :) :: x_vel
    real(rk), allocatable, dimension(:, :) :: y_vel
    real(rk), allocatable, dimension(:, :) :: pressure
    integer(ik) :: i, j

    allocate(rho_a_tilde, mold=cones%reference_density)
    allocate(density, mold=cones%density_p_prime)
    allocate(x_vel, mold=cones%density_p_prime)
    allocate(y_vel, mold=cones%density_p_prime)
    allocate(pressure, mold=cones%density_p_prime)

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

      !$omp parallel default(shared) private(i,j)
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
          pressure(i, j) = sum(p_i(:, i, j) * dtheta(:, i, j) - &
                               rho_a_tilde(i, j) * u_i(:, i, j) * sin_dtheta(:, i, j) + &
                               rho_a_tilde(i, j) * v_i(:, i, j) * cos_dtheta(:, i, j)) / (2.0_rk * pi)
        end do
      end do
      !$omp end do

      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          density(i, j) = rho_p_prime(i, j) + (pressure(i, j) / a_tilde(i, j)**2) - (pressure_p_prime(i, j) / a_tilde(i, j)**2)
        end do
      end do
      !$omp end do

      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          x_vel(i, j) = sum((-p_i(:, i, j) / rho_a_tilde(i, j)) * sin_dtheta(:, i, j) &
                            + u_i(:, i, j) * ((dtheta(:, i, j) / 2.0_rk) + (sin_d2theta(:, i, j) / 4.0_rk)) &
                            - v_i(:, i, j) * (cos_d2theta(:, i, j) / 4.0_rk)) / pi

          y_vel(i, j) = sum((p_i(:, i, j) / rho_a_tilde(i, j)) * cos_dtheta(:, i, j) &
                            - u_i(:, i, j) * (cos_d2theta(:, i, j) / 4.0_rk) &
                            + v_i(:, i, j) * ((dtheta(:, i, j) / 2.0_rk) - (sin_d2theta(:, i, j) / 4.0_rk))) / pi
        end do
      end do
      !$omp end do

      !$omp do
      do j = 1, cones%nj
        do i = 1, cones%ni
          primitive_vars(1, i, j) = density(i, j)
          primitive_vars(2, i, j) = x_vel(i, j)
          primitive_vars(3, i, j) = y_vel(i, j)
          primitive_vars(4, i, j) = pressure(i, j)
        end do
      end do
      !$omp end do

      !$omp end parallel
    end associate

    ! call cones(i=500,j=1)
    ! print*, "rho: ", density(500,1)
    ! print*, "u  : ", x_vel(500,1)
    ! print*, "v  : ", y_vel(500,1)
    ! print*, "p  : ", pressure(500,1)
    ! print*
    ! error stop

    deallocate(rho_a_tilde)
    deallocate(density)
    deallocate(x_vel)
    deallocate(y_vel)
    deallocate(pressure)

  end subroutine

end module mod_local_evo_operator
