module mod_grid_block
  !< Summary: Provide the implementation of a domain-decomposed grid block
  !< Date: 10/07/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_out => output_unit, std_error => error_unit

  implicit none

  type :: block_master_t
    !< Master class that knows the extents and metadata of the entire domain
    integer(ik) :: n_blocks = 1 !< How many total blocks?
    integer(ik), dimension(3) :: domain_shape = [1, 1, 1] !< (i blocks, j blocks, k blocks); # blocks in each dimension, e.g. 1x3x1 blocks
    integer(ik) :: dimensionality = 0 !< 1D, 2D, or 3D
    integer(ik), dimension(3) :: lbounds = [0, 0, 0]     !< (i, j, k); Lower bounds for the domain not including boundary ghost cells
    integer(ik), dimension(3) :: ubounds = [0, 0, 0]     !< (i, j, k); Upper bounds for the domain not including boundary ghost cells
    integer(ik), dimension(3) :: lbounds_bc = [0, 0, 0]  !< (i, j, k); Lower bounds for the domain including boundary ghost cells
    integer(ik), dimension(3) :: ubounds_bc = [0, 0, 0]  !< (i, j, k); Upper bounds for the domain including boundary ghost cells
  end type

  type, abstract :: grid_block_t
    private
    integer(ik), public :: n_halo_cells = 0  !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    logical :: is_axisymmetric = .false.     !< Is this an axisymmetric grid? (need for volume computation)
  contains
    procedure(initialize), deferred :: initialize
  end type grid_block_t

  abstract interface
    subroutine initialize(self)
      import :: grid_block_t
      class(grid_block_t), intent(inout) :: self
    end subroutine
  end interface

  type, extends(grid_block_t) :: grid_block_1d_t
    private

    ! Bounds information
    integer(ik), dimension(1), public :: dims = 0         !< (i, j); dimensions
    integer(ik), dimension(1), public :: lbounds = 0      !< (i, j); lower bounds (not including halo cells)
    integer(ik), dimension(1), public :: ubounds = 0      !< (i, j); upper bounds (not including halo cells)
    integer(ik), dimension(1), public :: lbounds_halo = 0 !< (i, j); lower bounds (including halo cells)
    integer(ik), dimension(1), public :: ubounds_halo = 0 !< (i, j); upper bounds (including halo cells)

    integer(ik), public :: ilo = 0 !< Lower bound for i
    integer(ik), public :: ihi = 0 !< Upper bound for i
    integer(ik), public :: ilo_halo = 0 !< Lower halo bound for i
    integer(ik), public :: ihi_halo = 0 !< Upper halo bound for i

    ! Boundary condition checks
    logical, public :: on_ihi_bc = .true. !< does this grid live on the ihi boundary?
    logical, public :: on_ilo_bc = .true. !< does this grid live on the ilo boundary?

    real(rk), dimension(:), allocatable, public :: cell_volume !< (i); volume of each cell
    real(rk), dimension(:), allocatable, public :: cell_dx     !< (i); dx spacing of each cell
    real(rk), dimension(:), allocatable, public :: cell_dy     !< (i); dy spacing of each cell
    real(rk), dimension(:), allocatable, public :: node_x      !< (i); x location of each node
    real(rk), dimension(:), allocatable, public :: node_y      !< (i); y location of each node
    real(rk), dimension(:), allocatable, public :: centroid_x  !< (i); x location of the cell centroid
    real(rk), dimension(:), allocatable, public :: centroid_y  !< (i); y location of the cell centroid
    real(rk), dimension(:, :), allocatable, public :: edge_lengths         !< ((edge_1:edge_n), i); length of each edge
    real(rk), dimension(:, :, :), allocatable, public :: edge_norm_vectors !< ((x,y), edge, i); normal direction vector of each face
  contains
    procedure :: initialize => init_1d_block
    final :: finalize_1d_block
  end type grid_block_1d_t

  type, extends(grid_block_t) :: grid_block_2d_t
    private

    ! Bounds information
    integer(ik), dimension(2), public :: dims = 0         !< (i, j); dimensions
    integer(ik), dimension(2), public :: lbounds = 0      !< (i, j); lower bounds (not including halo cells)
    integer(ik), dimension(2), public :: ubounds = 0      !< (i, j); upper bounds (not including halo cells)
    integer(ik), dimension(2), public :: lbounds_halo = 0 !< (i, j); lower bounds (including halo cells)
    integer(ik), dimension(2), public :: ubounds_halo = 0 !< (i, j); upper bounds (including halo cells)

    integer(ik), public :: ilo = 0 !< Lower bound for i
    integer(ik), public :: ihi = 0 !< Upper bound for i
    integer(ik), public :: jlo = 0 !< Lower bound for j
    integer(ik), public :: jhi = 0 !< Upper bound for j
    integer(ik), public :: ilo_halo = 0 !< Lower halo bound for i
    integer(ik), public :: ihi_halo = 0 !< Upper halo bound for i
    integer(ik), public :: jlo_halo = 0 !< Lower halo bound for j
    integer(ik), public :: jhi_halo = 0 !< Upper halo bound for j

    ! Boundary condition checks
    logical, public :: on_ihi_bc = .true. !< does this grid live on the ihi boundary?
    logical, public :: on_ilo_bc = .true. !< does this grid live on the ilo boundary?
    logical, public :: on_jhi_bc = .true. !< does this grid live on the jhi boundary?
    logical, public :: on_jlo_bc = .true. !< does this grid live on the jlo boundary?

    real(rk), dimension(:, :), allocatable, public :: cell_volume !< (i, j); volume of each cell
    real(rk), dimension(:, :), allocatable, public :: cell_dx     !< (i, j); dx spacing of each cell
    real(rk), dimension(:, :), allocatable, public :: cell_dy     !< (i, j); dy spacing of each cell
    real(rk), dimension(:, :), allocatable, public :: node_x      !< (i, j); x location of each node
    real(rk), dimension(:, :), allocatable, public :: node_y      !< (i, j); y location of each node
    real(rk), dimension(:, :), allocatable, public :: centroid_x  !< (i, j); x location of the cell centroid
    real(rk), dimension(:, :), allocatable, public :: centroid_y  !< (i, j); y location of the cell centroid
    real(rk), dimension(:, :, :), allocatable, public :: edge_lengths         !< ((edge_1:edge_n), i, j); length of each edge
    real(rk), dimension(:, :, :, :), allocatable, public :: edge_norm_vectors !< ((x,y), edge, i, j); normal direction vector of each face
  contains
    procedure :: initialize => init_2d_block
    final :: finalize_2d_block
  end type grid_block_2d_t

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

  subroutine init_1d_block(self)
    class(grid_block_1d_t), intent(inout) :: self
  end subroutine

  subroutine init_2d_block(self)
    class(grid_block_2d_t), intent(inout) :: self
  end subroutine

  subroutine init_3d_block(self)
    class(grid_block_3d_t), intent(inout) :: self
  end subroutine

  subroutine finalize_1d_block(self)
    type(grid_block_1d_t), intent(inout) :: self
  end subroutine

  subroutine finalize_2d_block(self)
    type(grid_block_2d_t), intent(inout) :: self
  end subroutine

  subroutine finalize_3d_block(self)
    type(grid_block_3d_t), intent(inout) :: self
  end subroutine

end module mod_grid_block
