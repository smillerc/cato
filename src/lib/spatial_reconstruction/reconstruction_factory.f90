module mod_reconstruction_factory

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_input, only: input_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_piecewise_constant_reconstruction, only: piecewise_constant_reconstruction_t
  use mod_piecewise_linear_reconstruction, only: piecewise_linear_reconstruction_t

  use mod_grid, only: grid_t

  implicit none

  private
  public :: reconstruction_factory

contains

  function reconstruction_factory(input, grid_target) result(operator)
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target
    class(abstract_reconstruction_t), pointer :: operator

    character(len=:), allocatable :: recon_type

    recon_type = trim(input%reconstruction_type)
    call debug_print('Making a "'//recon_type//'" reconstruction operator', __FILE__, __LINE__)

    select case(recon_type)
    case('piecewise_linear')
      allocate(piecewise_linear_reconstruction_t :: operator)
      call operator%initialize(input=input, grid_target=grid_target)
    case('piecewise_constant', 'cell_average')
      allocate(piecewise_constant_reconstruction_t :: operator)
      call operator%initialize(input=input, grid_target=grid_target)
    case default
      error stop "Error in reconstruction_factory_t, unknown reconstruction type"
    end select
  end function

end module mod_reconstruction_factory
