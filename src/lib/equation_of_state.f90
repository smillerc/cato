module mod_eos
  !> Summary:  Short module description
  !> Date: 04/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_globals, only: n_ghost_layers

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
    procedure, private :: total_energy_2d
    procedure, private :: sound_speed_2d
    procedure, private :: temperature_2d
    procedure, private :: conserved_to_primitive_2d
    procedure, private :: primitive_to_conserved_2d
    procedure, private :: total_energy_0d
    procedure, private :: sound_speed_0d
    procedure, private :: temperature_0d
    procedure, private :: conserved_to_primitive_0d
    procedure, private :: primitive_to_conserved_0d

    generic, public :: total_energy => total_energy_2d, total_energy_0d
    generic, public :: sound_speed => sound_speed_0d, sound_speed_2d
    generic, public :: temperature => temperature_0d, temperature_2d
    generic, public :: conserved_to_primitive => conserved_to_primitive_0d, conserved_to_primitive_2d
    generic, public :: primitive_to_conserved => primitive_to_conserved_0d, primitive_to_conserved_2d

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

  impure subroutine total_energy_2d(self, rho, u, v, p, E)
    !< Calculate total energy
    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: rho   !< density
    real(rk), dimension(:, :), intent(in) :: u    !< x-velocity
    real(rk), dimension(:, :), intent(in) :: v    !< y-velocity
    real(rk), dimension(:, :), intent(in) :: p    !< pressure
    real(rk), dimension(:, :), intent(inout) :: E    !< total energy

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: gamma_m_one

    gamma_m_one = self%gamma - 1.0_rk
    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    !$omp parallel default(none) &
    !$omp shared(rho, u, v, p, E) &
    !$omp firstprivate(ilo, ihi, jlo, jhi, gamma_m_one) &
    !$omp private(i, j)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        E(i, j) = (p(i, j) / (rho(i, j) * gamma_m_one)) + ((u(i, j)**2 + v(i, j)**2) / 2.0_rk)
      end do
    end do
    !$omp end do simd
    !$omp end parallel
  end subroutine total_energy_2d

  impure subroutine sound_speed_2d(self, p, rho, cs)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: p   !< pressure
    real(rk), dimension(:, :), intent(in) :: rho !< density
    real(rk), dimension(:, :), intent(inout) :: cs !< sound speed

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: gamma

    gamma = self%gamma
    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    !$omp parallel default(none) &
    !$omp shared(rho, p, cs) &
    !$omp firstprivate(ilo, ihi, jlo, jhi, gamma) &
    !$omp private(i, j)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        cs(i, j) = sqrt(gamma * p(i, j) / rho(i, j))
      end do
    end do
    !$omp end do simd
    !$omp end parallel
  end subroutine sound_speed_2d

  impure subroutine temperature_2d(self, p, rho, t)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: p
    real(rk), dimension(:, :), intent(in) :: rho
    real(rk), dimension(:, :), intent(inout) :: t

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: R

    R = self%R
    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    !$omp parallel default(none) &
    !$omp shared(rho, p, t) &
    !$omp firstprivate(ilo, ihi, jlo, jhi, R) &
    !$omp private(i, j)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        t(i, j) = p(i, j) / (rho(i, j) * R)
      end do
    end do
    !$omp end do simd
    !$omp end parallel
  end subroutine temperature_2d

  impure subroutine conserved_to_primitive_2d(self, rho, rho_u, rho_v, rho_E, u, v, p)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: rho   !< density
    real(rk), dimension(:, :), intent(in) :: rho_u !< density * x-velocity
    real(rk), dimension(:, :), intent(in) :: rho_v !< density * y-velocity
    real(rk), dimension(:, :), intent(in) :: rho_E !< density * total energy
    real(rk), dimension(:, :), intent(inout) :: u    !< x-velocity
    real(rk), dimension(:, :), intent(inout) :: v    !< y-velocity
    real(rk), dimension(:, :), intent(inout) :: p    !< pressure

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: gamma_m_one

    gamma_m_one = self%gamma - 1.0_rk

    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    !$omp parallel default(none) &
    !$omp shared(rho, u, v, p, rho_u, rho_v, rho_E) &
    !$omp firstprivate(ilo, ihi, jlo, jhi, gamma_m_one) &
    !$omp private(i, j)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        u(i, j) = rho_u(i, j) / rho(i, j)
        v(i, j) = rho_v(i, j) / rho(i, j)
        p(i, j) = rho(i, j) * gamma_m_one * ((rho_E(i, j) / rho(i, j)) - &
                                             ((u(i, j)**2 + v(i, j)**2) / 2.0_rk))
      end do
    end do
    !$omp end do simd
    !$omp end parallel

  end subroutine conserved_to_primitive_2d

  impure subroutine primitive_to_conserved_2d(self, rho, u, v, p, rho_u, rho_v, rho_E)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: rho
    real(rk), dimension(:, :), intent(in) :: u
    real(rk), dimension(:, :), intent(in) :: v
    real(rk), dimension(:, :), intent(in) :: p
    real(rk), dimension(:, :), intent(inout) :: rho_u
    real(rk), dimension(:, :), intent(inout) :: rho_v
    real(rk), dimension(:, :), intent(inout) :: rho_E

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: gamma_m_one

    gamma_m_one = self%gamma - 1.0_rk

    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    !$omp parallel default(none) &
    !$omp shared(rho, u, v, p, rho_u, rho_v, rho_E) &
    !$omp firstprivate(ilo, ihi, jlo, jhi, gamma_m_one) &
    !$omp private(i, j)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        rho_u(i, j) = u(i, j) * rho(i, j)
        rho_v(i, j) = v(i, j) * rho(i, j)
        rho_E(i, j) = (p(i, j) / gamma_m_one) + &
                      rho(i, j) * ((u(i, j)**2 + v(i, j)**2) / 2.0_rk)
      end do
    end do
    !$omp end do simd
    !$omp end parallel

  end subroutine primitive_to_conserved_2d

  pure subroutine total_energy_0d(self, rho, u, v, p, E)
    !< Calculate total energy
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho   !< density
    real(rk), intent(in) :: u    !< x-velocity
    real(rk), intent(in) :: v    !< y-velocity
    real(rk), intent(in) :: p    !< pressure
    real(rk), intent(inout) :: E    !< total energy

    ! Locals
    real(rk) :: gamma_m_one

    gamma_m_one = self%gamma - 1.0_rk

    E = (p / (rho * gamma_m_one)) + ((u**2 + v**2) / 2.0_rk)

  end subroutine total_energy_0d

  pure subroutine sound_speed_0d(self, p, rho, cs)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: p   !< pressure
    real(rk), intent(in) :: rho !< density
    real(rk), intent(inout) :: cs !< sound speed

    cs = sqrt(self%gamma * p / rho)

  end subroutine sound_speed_0d

  pure subroutine temperature_0d(self, p, rho, t)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: p
    real(rk), intent(in) :: rho
    real(rk), intent(inout) :: t

    t = p / (rho * self%R)

  end subroutine temperature_0d

  pure subroutine conserved_to_primitive_0d(self, rho, rho_u, rho_v, rho_E, u, v, p)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho   !< density
    real(rk), intent(in) :: rho_u !< density * x-velocity
    real(rk), intent(in) :: rho_v !< density * y-velocity
    real(rk), intent(in) :: rho_E !< density * total energy
    real(rk), intent(inout) :: u    !< x-velocity
    real(rk), intent(inout) :: v    !< y-velocity
    real(rk), intent(inout) :: p    !< pressure

    ! Locals
    real(rk) :: gamma_m_one

    gamma_m_one = self%gamma - 1.0_rk

    u = rho_u / rho
    v = rho_v / rho
    p = rho * gamma_m_one * ((rho_E / rho) - ((u**2 + v**2) / 2.0_rk))

  end subroutine conserved_to_primitive_0d

  pure subroutine primitive_to_conserved_0d(self, rho, u, v, p, rho_u, rho_v, rho_E)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho
    real(rk), intent(in) :: u
    real(rk), intent(in) :: v
    real(rk), intent(in) :: p
    real(rk), intent(inout) :: rho_u
    real(rk), intent(inout) :: rho_v
    real(rk), intent(inout) :: rho_E

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: gamma_m_one

    gamma_m_one = self%gamma - 1.0_rk

    rho_u = u * rho
    rho_v = v * rho
    rho_E = (p / gamma_m_one) + rho * ((u**2 + v**2) / 2.0_rk)

  end subroutine primitive_to_conserved_0d

end module mod_eos
