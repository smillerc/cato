module mod_slope_limiter
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_floating_point_utils, only: equal

  implicit none

  private
  public :: slope_limiter_t

  type :: slope_limiter_t
    character(:), allocatable :: name
    procedure(limit), pointer, nopass :: limit => null()
  contains
    final :: finalize
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
    case('upwind')
      limiter%limit => upwind_limit
    case('minmod')
      limiter%limit => minmod_limit
    case('weno')
      limiter%limit => weno_limit
    case default
      limiter%limit => upwind_limit
      limiter%name = 'upwind'
    end select

  end function

  subroutine finalize(self)
    type(slope_limiter_t), intent(inout) :: self
    if(associated(self%limit)) nullify(self%limit)
    if(allocated(self%name)) deallocate(self%name)
  end subroutine finalize

  pure function upwind_limit(a, b) result(slope)
    !< Slope limiter based on Equation 10 in DOI: 10.1016/j.jcp.2009.04.001
    !< Finds the equation of $$L(a,b) = \frac{max(ab,0)(a+b)}{a^2+b^2}$$

    real(rk) :: slope
    real(rk), intent(in) :: a, b

    real(rk) :: denom

    denom = a**2 + b**2

    if(denom > 0.0_rk) then
      slope = (max(a * b, 0.0_rk) * (a + b)) / denom
    else
      slope = 0.0_rk
    end if

  end function upwind_limit

  pure function minmod_limit(a, b) result(slope)
    real(rk) :: slope
    real(rk), intent(in) :: a, b

    if(a * b > 0.0_rk) then
      if(abs(a) < abs(b)) then
        slope = a
      else if(abs(b) < abs(a)) then
        slope = b
      else if(equal(abs(a), abs(b))) then
        slope = a
      end if
    else
      slope = 0.0_rk
    end if

  end function

  pure function weno_limit(a, b) result(slope)
    !< Basic WENO slope limiter

    real(rk) :: slope
    real(rk), intent(in) :: a, b
    real(rk) :: omega_1, omega_2  !< weights
    real(rk) :: epsilon

    epsilon = 1e-6_rk
    omega_1 = (epsilon + a**2)**(-2)
    omega_2 = (epsilon + b**2)**(-2)

    slope = (omega_1 * a + omega_2 * b) / (omega_1 + omega_2)

  end function

end module mod_slope_limiter
