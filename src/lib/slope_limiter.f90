module slope_limiter
  use iso_fortran_env, only: int32, real64

  implicit none

  private
  public :: limit
contains

  pure function limit(a, b) result(slope)
    !< Slope limiter based on Equation 10 in DOI: 10.1016/j.jcp.2009.04.001
    !< Finds the equation of $$L(a,b) = \frac{max(ab,0)(a+b)}{a^2+b^2}$$

    real(real64) :: slope
    real(real64), intent(in) :: a, b

    real(real64) :: max_ab

    max_ab = max(a * b, 0.0_real64)

    if(max_ab > 0) then
      slope = max_ab * (a + b) / (a**2 + b**2)
    else
      slope = 0.0_real64
    end if

  end function limit

end module slope_limiter
