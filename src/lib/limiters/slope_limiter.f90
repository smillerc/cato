module mod_slope_limiter
  use iso_fortran_env, only: int32, real64

  implicit none

  private
  public :: slope_limiter_t

  type, abstract :: slope_limiter_t
    character(:), allocatable :: name
    integer(int32) :: error_code
    character(:), allocatable :: error_message
  contains
    procedure(initialize), deferred, public :: initialize
    procedure(limit), deferred, nopass, public :: limit
  end type

  abstract interface
    subroutine initialize(self, name)
      import :: slope_limiter_t
      class(slope_limiter_t), intent(inout) :: self
      character(*) :: name
    end subroutine

    pure function limit(a, b) result(slope)
      import :: real64
      real(real64), intent(in) :: a, b
    end function
  end interface

  ! type, extends(slope_limiter_t) :: sun_ren_2009_limiter
  ! end type

contains

  ! pure function limit(a, b) result(slope)
  !   !< Slope limiter based on Equation 10 in DOI: 10.1016/j.jcp.2009.04.001
  !   !< Finds the equation of $$L(a,b) = \frac{max(ab,0)(a+b)}{a^2+b^2}$$

  !   real(real64) :: slope
  !   real(real64), intent(in) :: a, b

  !   real(real64) :: max_ab

  !   max_ab = max(a * b, 0.0_real64)

  !   if(max_ab > 0) then
  !     slope = max_ab * (a + b) / (a**2 + b**2)
  !   else
  !     slope = 0.0_real64
  !   end if

  ! end function limit

end module mod_slope_limiter
