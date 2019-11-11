module mod_slope_limiter
  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: slope_limiter_t

  type :: slope_limiter_t
    character(:), allocatable :: name
    procedure(limit), pointer, nopass :: limit
  end type

  interface slope_limiter_t
    module procedure constructor
  end interface

  interface
    pure function limit(a, b) result(slope)
      import :: rk
      real(rk) :: slope
      real(rk), intent(in) :: a, b
    end function
  end interface

contains

  type(slope_limiter_t) pure function constructor(name) result(limiter)
    character(len=*), intent(in) :: name

    limiter%name = trim(name)

    select case(trim(name))
    case('sun_ren_09')
      limiter%limit => sun_ren_09_limit
    case default
      limiter%limit => sun_ren_09_limit
      limiter%name = 'sun_ren_09'
    end select

  end function

  pure function sun_ren_09_limit(a, b) result(slope)
    !< Slope limiter based on Equation 10 in DOI: 10.1016/j.jcp.2009.04.001
    !< Finds the equation of $$L(a,b) = \frac{max(ab,0)(a+b)}{a^2+b^2}$$

    real(rk) :: slope
    real(rk), intent(in) :: a, b

    real(rk) :: max_ab

    max_ab = max(a * b, 0.0_rk)

    if(max_ab > 0) then
      slope = max_ab * (a + b) / (a**2 + b**2)
    else
      slope = 0.0_rk
    end if

  end function sun_ren_09_limit

end module mod_slope_limiter
