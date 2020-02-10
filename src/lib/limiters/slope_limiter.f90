module mod_slope_limiter
  !< Summary: Define the various slope limiters
  !< Note: The slope limiters in the form listed below is based on [1] and the symmetric forms
  !<       of the slope (not flux) limiters. Slope and flux limiters are very similar, but
  !<       have some subtle differences. See the reference for more info. The formulas are based on Eq 10.
  !< References:
  !< [1] M. Berger, M. Aftosmis, S. Muman, "Analysis of Slope Limiters on Irregular Grids",
  !<     43rd AIAA Aerospace Sciences Meeting and Exhibit (2005), https://doi.org/10.2514/6.2005-490

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_floating_point_utils, only: equal
  use math_constants, only: pi

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
    pure real(rk) function limit(f) result(slope)
      import :: rk
      real(rk), intent(in) :: f
    end function
  end interface

contains

  type(slope_limiter_t) pure function constructor(name) result(limiter)
    !< Slope limiter factory/constructor
    character(len=*), intent(in) :: name

    limiter%name = trim(name)

    select case(trim(name))
    case('barth_jesperson')
      limiter%limit => barth_jespersen
    case('minmod')
      limiter%limit => minmod
    case('van_leer')
      limiter%limit => van_leer
    case('sine')
      limiter%limit => sine
    case('van_albada')
      limiter%limit => van_albada
    case default
      limiter%limit => minmod
      limiter%name = 'minmod'
    end select
  end function

  subroutine finalize(self)
    type(slope_limiter_t), intent(inout) :: self
    if(associated(self%limit)) nullify(self%limit)
    if(allocated(self%name)) deallocate(self%name)
  end subroutine finalize

  pure real(rk) function barth_jespersen(f) result(slope)
    !< Barth Jesperson slope limiter. See Eq. 10 in [1]
    real(rk), intent(in) :: f !< f  = (u_i - u_{i-1}) / (u_{i+1} - u_{i-1})

    slope = min(1.0_rk, 4.0_rk * f, 4.0_rk * (1.0_rk - f))
  end function barth_jespersen

  pure real(rk) function van_leer(f) result(slope)
    !< van Leer slope limiter. See Eq. 10 in [1]
    real(rk), intent(in) :: f !< f  = (u_i - u_{i-1}) / (u_{i+1} - u_{i-1})

    slope = 4.0_rk * f * (1.0_rk - f)
  end function van_leer

  pure real(rk) function sine(f) result(slope)
    !< Sine slope limiter. See Eq. 10 in [1]
    real(rk), intent(in) :: f !< f  = (u_i - u_{i-1}) / (u_{i+1} - u_{i-1})

    slope = sin(pi * f)
  end function sine

  pure real(rk) function van_albada(f) result(slope)
    !< van Albada slope limiter. See Eq. 10 in [1]
    real(rk), intent(in) :: f !< f  = (u_i - u_{i-1}) / (u_{i+1} - u_{i-1})

    slope = (2.0_rk * f * (1.0_rk - f)) / (f**2 + (1.0_rk - f)**2)
  end function van_albada

  pure real(rk) function minmod(f) result(slope)
    !< Min-mod slope limiter. See Eq. 10 in [1]
    real(rk), intent(in) :: f !< f  = (u_i - u_{i-1}) / (u_{i+1} - u_{i-1})

    slope = min(2.0_rk * f, 2 * (1.0_rk - f))
  end function minmod

end module mod_slope_limiter
