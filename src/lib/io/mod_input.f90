module mod_input

  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: input_t

  type :: input_t
    integer(ik) :: ni
    integer(ik) :: nj

    real(rk) :: xmin
    real(rk) :: xmax
    real(rk) :: ymin
    real(rk) :: ymax

  contains
    procedure, public :: initialize
  end type input_t

contains
  subroutine initialize(self, ni, nj, xmin, xmax, ymin, ymax)
    class(input_t), intent(inout) :: self
    integer(ik), intent(in) :: ni, nj
    real(rk), intent(in) :: xmin, xmax, ymin, ymax

    self%ni = ni
    self%nj = nj
    self%xmin = xmin
    self%xmax = xmax
    self%ymin = ymin
    self%ymax = ymax

  end subroutine initialize
end module mod_input
