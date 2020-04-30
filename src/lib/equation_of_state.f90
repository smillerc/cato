module mod_eos
  !> Summary:  Short module description
  !> Date: 04/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic

#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */

  use mod_input, only: input_t
  implicit none

  private
  public :: set_equation_of_state, eos, universal_gas_constant

  real(rk), parameter :: universal_gas_constant = 8.31446261815324e7_rk !< R in ergs/(mol K)

  type :: eos_t
    private
    real(rk) :: gamma = 5.0_rk / 3.0_rk
    real(rk) :: c_p = (3.0_rk / 2.0_rk) * universal_gas_constant !< Specific heat capacity at constant pressure
    real(rk) :: c_v = (5.0_rk / 2.0_rk) * universal_gas_constant !< Specific heat capacity at constant volume
    real(rk) :: R = universal_gas_constant  !< Specific gas constant
  contains
    procedure, public :: get_gamma
    procedure, public :: set_gamma
    procedure, public :: get_cp
    procedure, public :: get_cv
    procedure, public :: total_energy
    procedure, public :: sound_speed
    procedure, public :: temperature
    procedure, public :: conserved_to_primitive
    procedure, public :: primitive_to_conserved
  end type eos_t

  interface new_eos
    module procedure :: constructor
  end interface

  type(eos_t) :: eos !< singleton equation of state object

contains
  type(eos_t) function constructor(input) result(eq_of_state)
    !< Constructor for the equation of state type
    class(input_t), intent(in) :: input
    call eq_of_state%set_gamma(input%polytropic_index)
  end function constructor

  subroutine set_gamma(self, gamma)
    class(eos_t), intent(inout) :: self
    real(rk), intent(in) :: gamma
    self%gamma = gamma
  end subroutine set_gamma

  real(rk) pure function get_gamma(self) result(gamma)
    class(eos_t), intent(in) :: self
    gamma = self%gamma
  end function get_gamma

  real(rk) pure function get_cp(self) result(cp)
    class(eos_t), intent(in) :: self
    cp = self%c_p
  end function get_cp

  real(rk) pure function get_cv(self) result(cv)
    class(eos_t), intent(in) :: self
    cv = self%c_v
  end function get_cv

  subroutine set_equation_of_state(input)
    class(input_t), intent(in) :: input
    eos = constructor(input)
  end subroutine set_equation_of_state

  elemental real(rk) function total_energy(self, rho, u, v, p) result(E)
    !< Calculate total energy

    !$omp declare simd
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho   !< density
    real(rk), intent(in) :: u    !< x-velocity
    real(rk), intent(in) :: v    !< y-velocity
    real(rk), intent(in) :: p    !< pressure

    E = (p / (rho * (self%gamma - 1.0_rk))) + ((u**2 + v**2) / 2.0_rk)
  end function total_energy

  elemental real(rk) function sound_speed(self, p, rho)
    !< Calculate sound speed

    !$omp declare simd
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: p   !< pressure
    real(rk), intent(in) :: rho !< density

    sound_speed = sqrt(self%gamma * p / rho)
  end function sound_speed

  elemental real(rk) function temperature(self, pressure, density)
    !< Calculate sound speed

    !$omp declare simd
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: density

    temperature = pressure / (density * self%R)
  end function temperature

  elemental subroutine conserved_to_primitive(self, rho, rho_u, rho_v, rho_E, u, v, p)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    !$omp declare simd
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho   !< density
    real(rk), intent(in) :: rho_u !< density * x-velocity
    real(rk), intent(in) :: rho_v !< density * y-velocity
    real(rk), intent(in) :: rho_E !< density * total energy
    real(rk), intent(out) :: u    !< x-velocity
    real(rk), intent(out) :: v    !< y-velocity
    real(rk), intent(out) :: p    !< pressure

    u = rho_u / rho
    v = rho_v / rho
    p = rho * (self%gamma - 1.0_rk) * ((rho_E / rho) - ((u**2 + v**2) / 2.0_rk))

  end subroutine conserved_to_primitive

  elemental subroutine primitive_to_conserved(self, rho, u, v, p, rho_u, rho_v, rho_E)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    !$omp declare simd
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho
    real(rk), intent(in) :: u
    real(rk), intent(in) :: v
    real(rk), intent(in) :: p
    real(rk), intent(out) :: rho_u
    real(rk), intent(out) :: rho_v
    real(rk), intent(out) :: rho_E

    rho_u = u * rho
    rho_v = v * rho
    rho_E = (p / (self%gamma - 1.0_rk)) + rho * ((u**2 + v**2) / 2.0_rk)

  end subroutine primitive_to_conserved

end module mod_eos
