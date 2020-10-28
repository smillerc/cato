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
    logical :: is_axisymmetric = .false. !< Is this an axisymmetric grid? (need for volume computation)
    logical :: is_uniform = .false.      !< Are all the cells the same size/shape?
    integer(ik) :: host_image_id = 0     !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    integer(ik) :: n_halo_cells = 0      !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    integer(ik) :: total_cells = 0       !<
    integer(ik), dimension(3) :: block_dims = 0   !< (ni,nj,nk); local block cell dimensions
    integer(ik), dimension(3) :: global_dims = 0   !< (ni,nj,nk); global cell dimensions
    integer(ik), dimension(3) :: lbounds = 0   !< (i,j,k); block lower cell bounds excluding halo cells
    integer(ik), dimension(3) :: ubounds = 0   !< (i,j,k); block upper cell bounds excluding halo cells
    integer(ik), dimension(3) :: lbounds_global = 0 !< (i,j,k); block lower cell bounds excluding halo cells
    integer(ik), dimension(3) :: ubounds_global = 0 !< (i,j,k); block upper cell bounds excluding halo cells
    integer(ik), dimension(3) :: lbounds_halo = 0   !< (i,j,k); block lower cell bounds including halo cells
    integer(ik), dimension(3) :: ubounds_halo = 0   !< (i,j,k); block upper cell bounds including halo cells

    logical :: on_ilo_bc = .false.
    logical :: on_ihi_bc = .false.
    logical :: on_jlo_bc = .false.
    logical :: on_jhi_bc = .false.
    logical :: on_klo_bc = .false.
    logical :: on_khi_bc = .false.

    integer(ik), dimension(:), allocatable :: neighbors !< (direction); neighbor image
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
