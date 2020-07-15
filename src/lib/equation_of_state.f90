! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_eos
  !< Summary:  Short module description
  !< Date: 04/22/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
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
    procedure, private :: total_enthalpy_2d
    procedure, private :: sound_speed_2d
    procedure, private :: temperature_2d
    procedure, private :: conserved_to_primitive_2d
    procedure, private :: primitive_to_conserved_2d
    procedure, private :: total_energy_0d
    procedure, private :: total_enthalpy_0d
    procedure, private :: sound_speed_0d
    procedure, private :: temperature_0d
    procedure, private :: conserved_to_primitive_0d
    procedure, private :: primitive_to_conserved_0d

    generic, public :: total_energy => total_energy_2d, total_energy_0d
    generic, public :: total_enthalpy => total_enthalpy_2d, total_enthalpy_0d
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
    real(rk), dimension(:, :), contiguous, intent(in) :: rho   !< density
    real(rk), dimension(:, :), contiguous, intent(in) :: u    !< x-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: v    !< y-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: p    !< pressure
    real(rk), dimension(:, :), contiguous, intent(inout) :: E    !< total energy

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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(E, rho, u, v, p:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        E(i, j) = (p(i, j) / (rho(i, j) * gamma_m_one)) + ((u(i, j)**2 + v(i, j)**2) / 2.0_rk)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine total_energy_2d

  impure subroutine total_enthalpy_2d(self, rho, u, v, p, H)
    !< Calculate total enthalpy
    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), contiguous, intent(in) :: rho  !< density
    real(rk), dimension(:, :), contiguous, intent(in) :: u    !< x-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: v    !< y-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: p    !< pressure
    real(rk), dimension(:, :), contiguous, intent(inout) :: H !< total enthalpy

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: gamma

    gamma = self%gamma
    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    !$omp parallel default(none) &
    !$omp shared(rho, u, v, p, H) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, gamma)
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(H, rho, u, v, p:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        H(i, j) = (p(i, j) / rho(i, j)) * (gamma / (gamma - 1.0_rk)) + ((u(i, j)**2 + v(i, j)**2) / 2.0_rk)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine total_enthalpy_2d

  impure subroutine sound_speed_2d(self, p, rho, cs)
    !< Calculate sound speed
    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), contiguous, intent(in) :: p     !< pressure
    real(rk), dimension(:, :), contiguous, intent(in) :: rho   !< density
    real(rk), dimension(:, :), contiguous, intent(inout) :: cs !< sound speed

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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(cs,p,rho:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        cs(i, j) = sqrt(gamma * abs(p(i, j) / rho(i, j)))
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine sound_speed_2d

  impure subroutine temperature_2d(self, p, rho, t)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), contiguous, intent(in) :: p
    real(rk), dimension(:, :), contiguous, intent(in) :: rho
    real(rk), dimension(:, :), contiguous, intent(inout) :: t

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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(r, p, rho:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        t(i, j) = p(i, j) / (rho(i, j) * R)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine temperature_2d

  impure subroutine conserved_to_primitive_2d(self, rho, rho_u, rho_v, rho_E, u, v, p)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), contiguous, intent(in) :: rho   !< density
    real(rk), dimension(:, :), contiguous, intent(in) :: rho_u !< density * x-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: rho_v !< density * y-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: rho_E !< density * total energy
    real(rk), dimension(:, :), contiguous, intent(inout) :: u  !< x-velocity
    real(rk), dimension(:, :), contiguous, intent(inout) :: v  !< y-velocity
    real(rk), dimension(:, :), contiguous, intent(inout) :: p  !< pressure

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: gamma_m_one, vel

    gamma_m_one = self%gamma - 1.0_rk

    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    !$omp parallel default(none) &
    !$omp shared(rho, u, v, p, rho_u, rho_v, rho_E) &
    !$omp firstprivate(ilo, ihi, jlo, jhi, gamma_m_one) &
    !$omp private(i, j, vel)
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(rho, u, v, p, rho_u, rho_v, rho_E:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        u(i, j) = rho_u(i, j) / rho(i, j)
        if(abs(u(i, j)) < 1e-30_rk) u(i, j) = 0.0_rk

        v(i, j) = rho_v(i, j) / rho(i, j)
        if(abs(v(i, j)) < 1e-30_rk) v(i, j) = 0.0_rk

        vel = 0.5_rk * (u(i, j)**2 + v(i, j)**2)
        if(vel < epsilon(1.0_rk)) vel = 0.0_rk

        p(i, j) = rho(i, j) * gamma_m_one * ((rho_E(i, j) / rho(i, j)) - vel)
      end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine conserved_to_primitive_2d

  impure subroutine primitive_to_conserved_2d(self, rho, u, v, p, rho_u, rho_v, rho_E)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    real(rk), dimension(:, :), contiguous, intent(in) :: rho      !< density
    real(rk), dimension(:, :), contiguous, intent(in) :: u        !< x-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: v        !< v-velocity
    real(rk), dimension(:, :), contiguous, intent(in) :: p        !< pressure
    real(rk), dimension(:, :), contiguous, intent(inout) :: rho_u !< density * x-velocity
    real(rk), dimension(:, :), contiguous, intent(inout) :: rho_v !< density * y-velocity
    real(rk), dimension(:, :), contiguous, intent(inout) :: rho_E !< density * total energy

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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(rho, u, v, p, rho_u, rho_v, rho_E:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        rho_u(i, j) = u(i, j) * rho(i, j)
        rho_v(i, j) = v(i, j) * rho(i, j)
        rho_E(i, j) = (p(i, j) / gamma_m_one) + &
                      rho(i, j) * ((u(i, j)**2 + v(i, j)**2) / 2.0_rk)
      end do
    end do
    !$omp end do
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

  pure subroutine total_enthalpy_0d(self, rho, u, v, p, H)
    !< Calculate total enthalpy
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho  !< density
    real(rk), intent(in) :: u    !< x-velocity
    real(rk), intent(in) :: v    !< y-velocity
    real(rk), intent(in) :: p    !< pressure
    real(rk), intent(inout) :: H !< total enthalpy

    H = (p / rho) * (self%gamma / (self%gamma - 1.0_rk)) + ((u**2 + v**2) / 2.0_rk)

  end subroutine total_enthalpy_0d

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
