#undef pure
module mod_eos

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  implicit none

  private
  public :: set_equation_of_state, eos

  type :: eos_t
    private
    real(rk) :: gamma = 5.0_rk / 3.0_rk
  contains
    procedure, public :: get_gamma
    procedure, public :: set_gamma
    procedure, public :: calc_pressure_from_dens_and_internal_energy
    procedure, public :: calc_internal_energy_from_press_and_dens
    procedure, public :: total_energy
    procedure, public :: sound_speed
    procedure, public :: energy_to_pressure
    procedure, public :: conserved_to_primitive
    procedure, public :: primitive_to_conserved
    procedure, public :: calc_density_from_isentropic_press
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
  end function

  subroutine set_gamma(self, gamma)
    class(eos_t), intent(inout) :: self
    real(rk), intent(in) :: gamma
    self%gamma = gamma
  end subroutine

  real(rk) pure function get_gamma(self) result(gamma)
    class(eos_t), intent(in) :: self
    gamma = self%gamma
  end function

  subroutine set_equation_of_state(input)
    class(input_t), intent(in) :: input
    eos = constructor(input)
  end subroutine

  elemental function calc_pressure_from_dens_and_internal_energy(self, density, internal_energy) result(pressure)
    !< Ideal gas law -> P = (1-gamma)*rho*e
    class(eos_t), intent(in) :: self
    real(rk) :: pressure
    real(rk), intent(in) :: density
    real(rk), intent(in) :: internal_energy

    pressure = (1.0_rk - self%gamma) * density * internal_energy
  end function

  elemental function calc_internal_energy_from_press_and_dens(self, pressure, density) result(internal_energy)
    !< Calculate internal energy
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: density
    real(rk) :: internal_energy

    internal_energy = (pressure / density) / (self%gamma - 1.0_rk)
  end function

  elemental real(rk) function total_energy(self, pressure, density, x_velocity, y_velocity)
    !< Calculate total energy
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: density
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: x_velocity
    real(rk), intent(in) :: y_velocity

    total_energy = (pressure / (self%gamma - 1.0_rk)) + 0.5_rk * density * (x_velocity**2 + y_velocity**2)
  end function

  elemental real(rk) function sound_speed(self, pressure, density)
    !< Calculate sound speed
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: density
    ! real(rk) :: sound_speed

    if(pressure < 0) error stop "Pressure is < 0 in eos_t%calc_sound_speed"
    if(density < 0) error stop "Density is < 0 in eos_t%calc_sound_speed"
    sound_speed = sqrt(self%gamma * pressure / density)
  end function

  elemental function calc_density_from_isentropic_press(self, P_1, P_2, rho_1) result(rho_2)
    !< Calculate the isentropic density based on a given pressure ratio and initial density
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: P_1
    real(rk), intent(in) :: P_2
    real(rk), intent(in) :: rho_1
    real(rk) :: rho_2

    if(P_1 <= 0.0_rk) error stop "P_1 <= 0"
    if(P_2 <= 0.0_rk) error stop "P_2 <= 0"
    if(rho_1 <= 0.0_rk) error stop "rho_1 <= 0"
    rho_2 = rho_1 * (P_2 / P_1)**(self%gamma - 1)
  end function

  pure function conserved_to_primitive(self, conserved_vars) result(primitive_vars)
    class(eos_t), intent(in) :: self
    real(rk), dimension(4), intent(in) :: conserved_vars !< [rho, rho*u, rho*v, E]
    real(rk), dimension(4) :: primitive_vars !< [rho, u, v, p]

    real(rk) :: rho
    real(rk) :: pressure
    real(rk) :: u
    real(rk) :: v

    rho = conserved_vars(1)
    u = conserved_vars(2) / rho
    v = conserved_vars(3) / rho

    associate(e=>conserved_vars(4), gamma=>self%get_gamma())
      pressure = (gamma - 1.0_rk) * (e - 0.5_rk * rho * (u**2 + v**2))
    end associate

    primitive_vars = [rho, u, v, pressure]

  end function

  pure function primitive_to_conserved(self, primitive_vars) result(conserved_vars)
    class(eos_t), intent(in) :: self
    real(rk), dimension(4), intent(in) :: primitive_vars !< [rho, u, v, p]
    real(rk), dimension(4) :: conserved_vars !< [rho, rho*u, rho*v, E]
    real(rk) :: e  !< total energy

    associate(gamma=>self%get_gamma(), &
              rho=>primitive_vars(1), &
              u=>primitive_vars(2), &
              v=>primitive_vars(3), &
              p=>primitive_vars(4))
      e = (p / (gamma - 1.0_rk)) + 0.5_rk * rho * (u**2 + v**2)
      conserved_vars = [rho, rho * u, rho * v, e]
    end associate

  end function

  elemental real(rk) function energy_to_pressure(self, energy, rho, u, v) result(pressure)
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: energy
    real(rk), intent(in) :: rho
    real(rk), intent(in) :: u
    real(rk), intent(in) :: v

    pressure = (self%gamma - 1.0_rk) * (energy - 0.5_rk * rho * (u**2 + v**2))
  end function

end module mod_eos
