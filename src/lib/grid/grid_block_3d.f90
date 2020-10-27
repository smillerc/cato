module mod_grid_block_3d

  use iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_grid_block, only: grid_block_t
  implicit none (type, external)

  type, extends(grid_block_t) :: grid_block_3d_t
    private

    ! Bounds information
    integer(ik), dimension(3), public :: dims = 0         !< (i, j); dimensions
    integer(ik), dimension(3), public :: lbounds = 0      !< (i, j); lower bounds (not including halo cells)
    integer(ik), dimension(3), public :: ubounds = 0      !< (i, j); upper bounds (not including halo cells)
    integer(ik), dimension(3), public :: lbounds_halo = 0 !< (i, j); lower bounds (including halo cells)
    integer(ik), dimension(3), public :: ubounds_halo = 0 !< (i, j); upper bounds (including halo cells)

    integer(ik), public :: ilo = 0 !< Lower bound for i
    integer(ik), public :: ihi = 0 !< Upper bound for i
    integer(ik), public :: jlo = 0 !< Lower bound for j
    integer(ik), public :: jhi = 0 !< Upper bound for j
    integer(ik), public :: klo = 0 !< Lower bound for k
    integer(ik), public :: khi = 0 !< Upper bound for k
    integer(ik), public :: ilo_halo = 0 !< Lower halo bound for i
    integer(ik), public :: ihi_halo = 0 !< Upper halo bound for i
    integer(ik), public :: jlo_halo = 0 !< Lower halo bound for j
    integer(ik), public :: jhi_halo = 0 !< Upper halo bound for j
    integer(ik), public :: klo_halo = 0 !< Lower halo bound for k
    integer(ik), public :: khi_halo = 0 !< Upper halo bound for k

    ! Boundary condition checks
    logical, public :: on_ihi_bc = .true. !< does this grid live on the ihi boundary?
    logical, public :: on_ilo_bc = .true. !< does this grid live on the ilo boundary?
    logical, public :: on_jhi_bc = .true. !< does this grid live on the jhi boundary?
    logical, public :: on_jlo_bc = .true. !< does this grid live on the jlo boundary?
    logical, public :: on_khi_bc = .true. !< does this grid live on the khi boundary?
    logical, public :: on_klo_bc = .true. !< does this grid live on the klo boundary?

    real(rk), dimension(:, :, :), allocatable, public :: cell_volume !< (i, j, k); volume of each cell
    real(rk), dimension(:, :, :), allocatable, public :: cell_dx     !< (i, j, k); dx spacing of each cell
    real(rk), dimension(:, :, :), allocatable, public :: cell_dy     !< (i, j, k); dy spacing of each cell
    real(rk), dimension(:, :, :), allocatable, public :: node_x      !< (i, j, k); x location of each node
    real(rk), dimension(:, :, :), allocatable, public :: node_y      !< (i, j, k); y location of each node
    real(rk), dimension(:, :, :), allocatable, public :: centroid_x  !< (i, j, k); x location of the cell centroid
    real(rk), dimension(:, :, :), allocatable, public :: centroid_y  !< (i, j, k); y location of the cell centroid
    real(rk), dimension(:, :, :, :), allocatable, public :: edge_lengths         !< ((edge_1:edge_n), i, j); length of each edge
    real(rk), dimension(:, :, :, :, :), allocatable, public :: edge_norm_vectors !< ((x,y), edge, i, j); normal direction vector of each face
  contains
    procedure :: initialize => init_3d_block
    final :: finalize_3d_block
  end type grid_block_3d_t
contains
  
  subroutine init_3d_block(self)
    class(grid_block_3d_t), intent(inout) :: self
  end subroutine

  subroutine finalize_3d_block(self)
    type(grid_block_3d_t), intent(inout) :: self
  end subroutine
end module mod_grid_block_3d
