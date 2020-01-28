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
    ! procedure, public :: calc_pressure_from_dens_and_internal_energy
    ! procedure, public :: calc_internal_energy_from_press_and_dens
    procedure, public :: calculate_total_energy
    procedure, public :: calc_sound_speed
    procedure, public :: calc_sound_speed_from_primitive
    procedure, public :: total_energy_to_pressure
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

  subroutine set_equation_of_state(input)
    class(input_t), intent(in) :: input
    eos = constructor(input)
  end subroutine set_equation_of_state

  ! elemental function calc_pressure_from_dens_and_internal_energy(self, density, internal_energy) result(pressure)
  !   !< Ideal gas law -> P = (1-gamma)*rho*e
  !   class(eos_t), intent(in) :: self
  !   real(rk) :: pressure
  !   real(rk), intent(in) :: density
  !   real(rk), intent(in) :: internal_energy

  !   pressure = (1.0_rk - self%gamma) * density * internal_energy
  ! end function calc_pressure_from_dens_and_internal_energy

  ! elemental function calc_internal_energy_from_press_and_dens(self, pressure, density) result(internal_energy)
  !   !< Calculate internal energy
  !   class(eos_t), intent(in) :: self
  !   real(rk), intent(in) :: pressure
  !   real(rk), intent(in) :: density
  !   real(rk) :: internal_energy

  !   internal_energy = (pressure / density) / (self%gamma - 1.0_rk)
  ! end function calc_internal_energy_from_press_and_dens

  elemental real(rk) function calculate_total_energy(self, pressure, density, x_velocity, y_velocity) result(total_energy)
    !< Calculate total energy
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: density
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: x_velocity
    real(rk), intent(in) :: y_velocity

    ! total_energy = (pressure / (self%gamma - 1.0_rk)) + (density / 2.0_rk) * (x_velocity**2 + y_velocity**2)
    total_energy = (pressure / (density * (self%gamma - 1.0_rk))) + ((x_velocity**2 + y_velocity**2) / 2.0_rk)
  end function calculate_total_energy

  elemental real(rk) function calc_sound_speed(self, pressure, density) result(sound_speed)
    !< Calculate sound speed
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: density

    if(pressure < 0) error stop "Pressure is < 0 in eos_t%calc_sound_speed"
    if(density < 0) error stop "Density is < 0 in eos_t%calc_sound_speed"
    sound_speed = sqrt(self%gamma * pressure / density)
  end function calc_sound_speed

  subroutine calc_sound_speed_from_primitive(self, primitive_vars, sound_speed)
    !< Calculate sound speed
    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(in) :: primitive_vars
    real(rk), dimension(:, :), intent(out) :: sound_speed

    if(minval(primitive_vars(4, :, :)) < 0.0_rk) then
      write(*, *) "Pressure is < 0 in eos_t%calc_sound_speed_from_primitive() at (i,j)", minloc(primitive_vars(4, :, :))
      error stop
    end if

    if(minval(primitive_vars(1, :, :)) < 0.0_rk) then
      write(*, *) "Density is < 0 in eos_t%calc_sound_speed_from_primitive() at (i,j)", minloc(primitive_vars(1, :, :))
    end if

    sound_speed = sqrt(self%gamma * primitive_vars(4, :, :) / primitive_vars(1, :, :))
  end subroutine calc_sound_speed_from_primitive

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
  end function calc_density_from_isentropic_press

  elemental real(rk) function total_energy_to_pressure(self, total_energy, rho, u, v) result(pressure)
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: total_energy
    real(rk), intent(in) :: rho
    real(rk), intent(in) :: u
    real(rk), intent(in) :: v

    ! pressure = (self%gamma - 1.0_rk) * (total_energy - (rho / 2.0_rk) * (u**2 + v**2))
    pressure = rho * (self%gamma - 1.0_rk) * (total_energy - ((u**2 + v**2) / 2.0_rk))
  end function total_energy_to_pressure

end module mod_eos
