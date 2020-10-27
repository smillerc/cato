module mod_grid_block
  !< Summary: Provide the implementation of a domain-decomposed grid block
  !< Date: 10/07/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_out => output_unit, std_error => error_unit
  use mod_input, only: input_t

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

    integer(ik) :: host_image_id = 0  !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    integer(ik) :: n_halo_cells = 0   !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    logical :: is_axisymmetric = .false. !< Is this an axisymmetric grid? (need for volume computation)
    logical :: is_uniform = .false.      !< Are all the cells the same size/shape?

    ! Bounds information
    integer(ik), dimension(:), allocatable :: global_node_dims !< (i, j, k); # of array indices in each dimension
    integer(ik), dimension(:), allocatable :: global_cell_dims !< (i, j, k); # of array indices in each dimension

    integer(ik), dimension(:), allocatable :: domain_node_shape !< (i, j, k); # of array indices in each dimension
    integer(ik), dimension(:), allocatable :: domain_cell_shape !< (i, j, k); # of array indices in each dimension

    integer(ik), dimension(:), allocatable :: global_node_dims_halo !< (i, j, k); # of array indices in each dimension
    integer(ik), dimension(:), allocatable :: global_cell_dims_halo !< (i, j, k); # of array indices in each dimension

    integer(ik), dimension(:), allocatable :: domain_node_shape_halo !< (i, j, k); # of array indices in each dimension
    integer(ik), dimension(:), allocatable :: domain_cell_shape_halo !< (i, j, k); # of array indices in each dimension

    integer(ik), dimension(:), allocatable :: node_lbounds !< (i, j, k); lower bounds (not including halo cells)
    integer(ik), dimension(:), allocatable :: node_ubounds !< (i, j, k); upper bounds (not including halo cells)

    integer(ik), dimension(:), allocatable :: node_lbounds_halo !< (i, j, k); lower bounds (including halo cells)
    integer(ik), dimension(:), allocatable :: node_ubounds_halo !< (i, j, k); upper bounds (including halo cells)

    integer(ik), dimension(:), allocatable :: cell_lbounds !< (i, j, k); lower bounds (not including halo cells)
    integer(ik), dimension(:), allocatable :: cell_ubounds !< (i, j, k); upper bounds (not including halo cells)

    integer(ik), dimension(:), allocatable :: cell_lbounds_halo !< (i, j, k); lower bounds (including halo cells)
    integer(ik), dimension(:), allocatable :: cell_ubounds_halo !< (i, j, k); upper bounds (including halo cells)

    ! Boundary condition checks
    logical, dimension(:,:), allocatable :: on_bc !< ((lo,hi),(i,j,k)) does this grid live on the ihi boundary?

  contains
    procedure(initialize), deferred :: initialize
  end type grid_block_t

  abstract interface
    subroutine initialize(self, input)
      import :: grid_block_t, input_t
      class(grid_block_t), intent(inout) :: self
      class(input_t), intent(in) :: input
    end subroutine
  end interface
contains


end module mod_grid_block
