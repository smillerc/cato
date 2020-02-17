module mod_slope_limiter
  !< Summary: Define the various slope limiters
  !< Note: The slope limiters in the form listed below is based on [1] the slope (not flux)
  !<       limiters. Slope and flux limiters are very similar, but have some subtle differences.
  !<       See the reference for more info. The formulas are based on Eq 10.
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
    real(rk) function limit(R) result(phi_lim)
      import :: rk
      real(rk), intent(in) :: R
    end function
  end interface

contains

  type(slope_limiter_t) function constructor(name) result(limiter)
    !< Slope limiter factory/constructor
    character(len=*), intent(in) :: name

    limiter%name = trim(name)

    select case(trim(name))
    case('none')
      limiter%limit => none
    case('barth_jesperson')
      limiter%limit => barth_jespersen
    case('minmod')
      limiter%limit => minmod
    case('van_leer')
      limiter%limit => van_leer
    case('van_albada')
      limiter%limit => van_albada
    case default
      error stop "Error in slope_limiter_t%constructor(): Unknown slope limiter name"
    end select
  end function

  subroutine finalize(self)
    type(slope_limiter_t), intent(inout) :: self
    if(associated(self%limit)) nullify(self%limit)
    if(allocated(self%name)) deallocate(self%name)
  end subroutine finalize

  real(rk) function none(R) result(phi_lim)
    !< Unlimited slope
    real(rk), intent(in) :: R !< R  = (u_{i+1} - u_i) / (u_i - u_{i-1})
    phi_lim = 1.0_rk
  end function none

  real(rk) function barth_jespersen(R) result(phi_lim)
    !< Barth Jesperson slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: R !< R  = (u_{i+1} - u_i) / (u_i - u_{i-1})
    if(R <= 0.0_rk) then
      phi_lim = 0.0_rk
    else
      phi_lim = min(1.0_rk, 4.0_rk / (R + 1.0_rk),(4.0_rk * R) / (R + 1.0_rk))
    end if
  end function barth_jespersen

  pure real(rk) function van_leer(R) result(phi_lim)
    !< van Leer slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: R !< R  = (u_{i+1} - u_i) / (u_i - u_{i-1})
    if(R <= 0.0_rk) then
      phi_lim = 0.0_rk
    else
      phi_lim = 4.0_rk * R / (R + 1.0_rk)**2
    end if
  end function van_leer

  pure real(rk) function van_albada(R) result(phi_lim)
    !< van Albada slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: R !< R  = (u_{i+1} - u_i) / (u_i - u_{i-1})
    if(R <= 0.0_rk) then
      phi_lim = 0.0_rk
    else
      phi_lim = (2.0_rk * R) / (R**2 + 1.0_rk)
    end if
  end function van_albada

  pure real(rk) function minmod(R) result(phi_lim)
    !< Min-mod slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: R !< R  = (u_{i+1} - u_i) / (u_i - u_{i-1})
    if(R <= 0.0_rk) then
      phi_lim = 0.0_rk
    else
      phi_lim = min(2.0_rk / (1.0_rk + R),(2.0_rk * R) / (1.0_rk + R))
    end if
  end function minmod

end module mod_slope_limiter
