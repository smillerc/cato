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
    procedure, public :: calculate_total_energy
    procedure, public :: sound_speed
    procedure, public :: sound_speed_from_primitive
    procedure, public :: sound_speed_from_conserved
    procedure, public :: conserved_to_primitive
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

  elemental real(rk) function sound_speed(self, pressure, density)
    !< Calculate sound speed
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: density

    ! if(pressure < 0) error stop "Pressure is < 0 in eos_t%sound_speed"
    ! if(density < 0) error stop "Density is < 0 in eos_t%sound_speed"
    sound_speed = sqrt(self%gamma * abs(pressure / density))
  end function sound_speed

  subroutine sound_speed_from_primitive(self, primitive_vars, sound_speed)
    !< Calculate sound speed using the primitive variables [rho, u, v, p]
    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(in) :: primitive_vars
    real(rk), dimension(:, :), intent(out) :: sound_speed

    ! if(minval(primitive_vars(4, :, :)) < 0.0_rk) then
    !   print*, 'minloc(primitive_vars(4, :, :))', minloc(primitive_vars(4, :, :))
    !   error stop "Pressure is < 0 in eos_t%calc_sound_speed_from_primitive()"
    ! end if

    ! if(minval(primitive_vars(1, :, :)) < 0.0_rk) then
    !   print*, 'minloc(primitive_vars(1, :, :))', minloc(primitive_vars(1, :, :))
    !   error stop "Density is < 0 in eos_t%calc_sound_speed_from_primitive()"
    ! end if

    sound_speed = sqrt(self%gamma * abs(primitive_vars(4, :, :) / primitive_vars(1, :, :)))
  end subroutine sound_speed_from_primitive

  subroutine sound_speed_from_conserved(self, conserved_vars, sound_speed)
    !< Calculate the sound speed using the conserved variables [rho, rho u, rho v, rho E]
    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(in) :: conserved_vars
    real(rk), dimension(:, :), intent(out) :: sound_speed
    real(rk), dimension(:, :, :), allocatable :: primitive_vars

    allocate(primitive_vars, mold=conserved_vars)
    call self%conserved_to_primitive(conserved_vars, primitive_vars)
    call self%sound_speed_from_primitive(primitive_vars, sound_speed)
    deallocate(primitive_vars)
  end subroutine

  elemental function calc_density_from_isentropic_press(self, P_1, P_2, rho_1) result(rho_2)
    !< Calculate the isentropic density based on a given pressure ratio and initial density
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: P_1
    real(rk), intent(in) :: P_2
    real(rk), intent(in) :: rho_1
    real(rk) :: rho_2

    if(P_1 <= 0.0_rk) error stop "Error in eos_t%calc_density_from_isentropic_press(): P_1 <= 0"
    if(P_2 <= 0.0_rk) error stop "Error in eos_t%calc_density_from_isentropic_press(): P_2 <= 0"
    if(rho_1 <= 0.0_rk) error stop "Error in eos_t%calc_density_from_isentropic_press(): rho_1 <= 0"
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

  subroutine conserved_to_primitive(self, conserved_vars, primitive_vars)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(in) :: conserved_vars
    real(rk), dimension(:, :, :), intent(out) :: primitive_vars

    ! rho
    primitive_vars(1, :, :) = conserved_vars(1, :, :)

    ! u
    primitive_vars(2, :, :) = conserved_vars(2, :, :) / conserved_vars(1, :, :)

    ! v
    primitive_vars(3, :, :) = conserved_vars(3, :, :) / conserved_vars(1, :, :)

    ! pressure
    associate(gamma=>self%gamma, E=>conserved_vars(4, :, :) / conserved_vars(1, :, :), &
              rho=>conserved_vars(1, :, :), &
              u=>primitive_vars(2, :, :), &
              v=>primitive_vars(3, :, :))

      primitive_vars(4, :, :) = rho * (gamma - 1.0_rk) * (E - ((u**2 + v**2) / 2.0_rk))
    end associate

  end subroutine conserved_to_primitive

end module mod_eos
