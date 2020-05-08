module mod_piecewise_constant_reconstruction
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use mod_globals, only: debug_print, n_ghost_layers
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t

  implicit none

  private
  public :: piecewise_constant_reconstruction_t

  type, extends(abstract_reconstruction_t) :: piecewise_constant_reconstruction_t
  contains
    procedure, public :: initialize => init_first_order
    procedure, public :: reconstruct
    final :: finalize
  end type

contains

  subroutine init_first_order(self, input, grid_target)
    !< Construct the piecewise_constant_reconstruction_t type

    class(piecewise_constant_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_reconstruction_t', __FILE__, __LINE__)

    self%order = 1
    self%name = 'cell_average_reconstruction'

    self%grid => grid_target

  end subroutine init_first_order

  subroutine finalize(self)
    !< Finalize the piecewise_constant_reconstruction_t type
    type(piecewise_constant_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running piecewise_constant_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    ! if(associated(self%rho)) nullify(self%rho)
    ! if(associated(self%u)) nullify(self%u)
    ! if(associated(self%v)) nullify(self%v)
    ! if(associated(self%p)) nullify(self%p)

    ! if(allocated(self%cell_gradient)) deallocate(self%cell_gradient)

  end subroutine finalize

  subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds)
    !< Reconstruct the entire domain. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient

    class(piecewise_constant_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), contiguous :: primitive_var !< (i,j); cell primitive variable to reconstruct
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(out), contiguous :: reconstructed_var
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

  end subroutine reconstruct

end module mod_piecewise_constant_reconstruction
