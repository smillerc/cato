module mod_first_order_reconstruction
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use mod_globals, only: debug_print, n_ghost_layers
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t

  implicit none

  private
  public :: first_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: first_order_reconstruction_t
  contains
    procedure, public :: initialize => init_first_order
    procedure, public :: reconstruct
    ! procedure, public :: copy
    final :: finalize
  end type

contains

  subroutine init_first_order(self, input, grid_target)
    !< Construct the first_order_reconstruction_t type

    class(first_order_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_reconstruction_t', __FILE__, __LINE__)

    self%order = 1
    self%name = 'cell_average_reconstruction'

    self%grid => grid_target

  end subroutine init_first_order

  subroutine finalize(self)
    !< Finalize the first_order_reconstruction_t type
    type(first_order_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running first_order_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    ! if(associated(self%rho)) nullify(self%rho)
    ! if(associated(self%u)) nullify(self%u)
    ! if(associated(self%v)) nullify(self%v)
    ! if(associated(self%p)) nullify(self%p)

    ! if(allocated(self%cell_gradient)) deallocate(self%cell_gradient)

  end subroutine finalize

  ! subroutine copy(out_recon, in_recon)
  !   class(abstract_reconstruction_t), intent(in) :: in_recon
  !   class(first_order_reconstruction_t), intent(inout) :: out_recon

  !   call debug_print('Running first_order_reconstruction_t%copy()', __FILE__, __LINE__)

  !   if(associated(out_recon%grid)) nullify(out_recon%grid)
  !   out_recon%grid => in_recon%grid

  !   if(associated(out_recon%rho)) nullify(out_recon%primitive_vars)
  !   out_recon%primitive_vars => in_recon%primitive_vars

  !   if(associated(out_recon%u)) nullify(out_recon%primitive_vars)
  !   out_recon%primitive_vars => in_recon%primitive_vars

  !   if(associated(out_recon%primitive_vars)) nullify(out_recon%primitive_vars)
  !   out_recon%primitive_vars => in_recon%primitive_vars

  !   if(associated(out_recon%primitive_vars)) nullify(out_recon%primitive_vars)
  !   out_recon%primitive_vars => in_recon%primitive_vars

  !   if(allocated(out_recon%name)) deallocate(out_recon%name)
  !   allocate(out_recon%name, source=in_recon%name)

  ! end subroutine

  ! pure function reconstruct_point(self, xy, cell_ij) result(U_bar)
  !   !< Reconstruct the value of the primitive variables (U) at location (x,y)
  !   !< withing a cell (i,j)

  !   class(first_order_reconstruction_t), intent(in) :: self
  !   real(rk), dimension(2), intent(in) :: xy !< where should U_bar be reconstructed at?
  !   real(rk), dimension(4) :: U_bar  !< U_bar = reconstructed [rho, u, v, p]
  !   integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
  !   integer(ik) :: i, j

  !   ! i = cell_ij(1); j = cell_ij(2)
  !   ! U_bar = self%primitive_vars(:, i, j)

  ! end function reconstruct_point

  subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds)
    !< Reconstruct the entire domain. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient

    class(first_order_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: primitive_var !< (i,j); cell primitive variable to reconstruct
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(out) :: reconstructed_var
    !< ((corner1:midpoint4), i, j); reconstructed variable, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    integer(ik) :: i, j, p
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc

    if(.not. associated(self%grid)) error stop "Grid not associated"

    ! Bounds do not include ghost cells. Ghost cells get their
    ! reconstructed values and gradients from the boundary conditions
    ilo_bc = lbound(primitive_var, dim=1)
    ihi_bc = ubound(primitive_var, dim=1)
    jlo_bc = lbound(primitive_var, dim=2)
    jhi_bc = ubound(primitive_var, dim=2)

    ilo = ilo_bc + n_ghost_layers
    ihi = ihi_bc - n_ghost_layers
    jlo = jlo_bc + n_ghost_layers
    jhi = jhi_bc - n_ghost_layers

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(reconstructed_var, primitive_var)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        do p = 1, 8
          reconstructed_var(p, i, j) = primitive_var(i, j)
        end do
      end do
    end do
    !$omp end do simd
    !$omp end parallel

    self%domain_has_been_reconstructed = .true.
  end subroutine reconstruct

end module mod_first_order_reconstruction
