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
  !< Summary:  Provide basic equation of state (EOS) functions
  !< Date: 04/22/2020
  !< Author: Sam Miller

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_field, only: field_2d_t
  use mod_globals, only: n_ghost_layers
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
    procedure, private :: total_energy_field_2d, total_energy_real
    procedure, private :: total_enthalpy_field_2d, total_enthalpy_real
    procedure, private :: sound_speed_field_2d, sound_speed_real
    procedure, private :: temperature_field_2d, temperature_real
    procedure, private :: conserved_to_primitive_field_2d, conserved_to_primitive_real
    procedure, private :: primitive_to_conserved_field_2d, primitive_to_conserved_real

    generic, public :: total_energy => total_energy_field_2d, total_energy_real
    generic, public :: total_enthalpy => total_enthalpy_field_2d, total_enthalpy_real
    generic, public :: sound_speed => sound_speed_real, sound_speed_field_2d
    generic, public :: temperature => temperature_real, temperature_field_2d
    generic, public :: conserved_to_primitive => conserved_to_primitive_real, conserved_to_primitive_field_2d
    generic, public :: primitive_to_conserved => primitive_to_conserved_real, primitive_to_conserved_field_2d

  endtype eos_t

  interface new_eos
    module procedure :: constructor
  endinterface

  type(eos_t) :: eos !< singleton equation of state object

contains
  type(eos_t) function constructor(input) result(eq_of_state)
    !< Constructor for the equation of state type
    class(input_t), intent(in) :: input
    call eq_of_state%set_gamma(input%polytropic_index)
  endfunction constructor

  subroutine set_gamma(self, gamma)
    class(eos_t), intent(inout) :: self
    real(rk), intent(in) :: gamma
    self%gamma = gamma
  endsubroutine set_gamma

  real(rk) pure function get_gamma(self) result(gamma)
    class(eos_t), intent(in) :: self
    gamma = self%gamma
  endfunction get_gamma

  real(rk) pure function get_cp(self) result(cp)
    class(eos_t), intent(in) :: self
    cp = self%c_p
  endfunction get_cp

  real(rk) pure function get_cv(self) result(cv)
    class(eos_t), intent(in) :: self
    cv = self%c_v
  endfunction get_cv

  subroutine set_equation_of_state(input)
    class(input_t), intent(in) :: input
    eos = constructor(input)
  endsubroutine set_equation_of_state

  elemental subroutine total_energy_real(self, rho, u, v, p, E)
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

    E = (p / (rho * gamma_m_one)) + 0.5_rk * (u * u + v * v)

  endsubroutine total_energy_real

  elemental subroutine total_enthalpy_real(self, rho, u, v, p, H)
    !< Calculate total enthalpy
    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: rho  !< density
    real(rk), intent(in) :: u    !< x-velocity
    real(rk), intent(in) :: v    !< y-velocity
    real(rk), intent(in) :: p    !< pressure
    real(rk), intent(inout) :: H !< total enthalpy

    H = (p / rho) * (self%gamma / (self%gamma - 1.0_rk)) + 0.5_rk * (u * u + v * v)

  endsubroutine total_enthalpy_real

  elemental subroutine sound_speed_real(self, p, rho, cs)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: p   !< pressure
    real(rk), intent(in) :: rho !< density
    real(rk), intent(inout) :: cs !< sound speed

    cs = sqrt(self%gamma * p / rho)

  endsubroutine sound_speed_real

  elemental subroutine temperature_real(self, p, rho, t)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    real(rk), intent(in) :: p
    real(rk), intent(in) :: rho
    real(rk), intent(inout) :: t

    t = p / (rho * self%R)

  endsubroutine temperature_real

  elemental subroutine conserved_to_primitive_real(self, rho, rho_u, rho_v, rho_E, u, v, p)
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
    p = rho * gamma_m_one * ((rho_E / rho) - 0.5_rk * (u * u + v * v))

  endsubroutine conserved_to_primitive_real

  elemental subroutine primitive_to_conserved_real(self, rho, u, v, p, rho_u, rho_v, rho_E)
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
    real(rk) :: gamma_m_one

    gamma_m_one = self%gamma - 1.0_rk

    rho_u = u * rho
    rho_v = v * rho
    rho_E = (p / gamma_m_one) + rho * 0.5_rk * (u * u + v * v)

  endsubroutine primitive_to_conserved_real

  subroutine total_energy_field_2d(self, rho, u, v, p, E)
    !< Calculate total energy
    class(eos_t), intent(in) :: self
    class(field_2d_t), intent(in) :: rho   !< density
    class(field_2d_t), intent(in) :: u    !< x-velocity
    class(field_2d_t), intent(in) :: v    !< y-velocity
    class(field_2d_t), intent(in) :: p    !< pressure
    class(field_2d_t), intent(inout) :: E    !< total energy

    ! Locals
    real(rk) :: gamma_m_one
    integer(ik) :: i, j

    gamma_m_one = self%gamma - 1.0_rk

    ! E = (p / (rho * gamma_m_one)) + (0.5_rk * (u * u + v * v))

    do j = rho%jlo, rho%jhi
      !$omp simd
      do i = rho%ilo, rho%ihi
        E%data(i,j) = (p%data(i,j) / (rho%data(i,j) * gamma_m_one)) + (0.5_rk * (u%data(i,j) * u%data(i,j) + v%data(i,j) * v%data(i,j)))
      enddo
      !$omp end simd
    enddo

    E%name = 'E'

  endsubroutine total_energy_field_2d

  subroutine total_enthalpy_field_2d(self, rho, u, v, p, H)
    !< Calculate total enthalpy
    class(eos_t), intent(in) :: self
    class(field_2d_t), intent(in) :: rho  !< density
    class(field_2d_t), intent(in) :: u    !< x-velocity
    class(field_2d_t), intent(in) :: v    !< y-velocity
    class(field_2d_t), intent(in) :: p    !< pressure
    class(field_2d_t), intent(inout) :: H !< total enthalpy

    real(rk) :: gamma_term
    integer(ik) :: i, j

    gamma_term = self%gamma / (self%gamma - 1.0_rk)

    ! H = (p / rho) * (self%gamma / (self%gamma - 1.0_rk)) + (0.5_rk * (u * u + v * v))

    do j = rho%jlo, rho%jhi
      !$omp simd
      do i = rho%ilo, rho%ihi
        H%data(i,j) = (p%data(i,j) / rho%data(i,j)) * gamma_term + (0.5_rk * (u%data(i,j) * u%data(i,j) + v%data(i,j) * v%data(i,j)))
      enddo
      !$omp end simd
    enddo

    H%name = 'H'

  endsubroutine total_enthalpy_field_2d

  subroutine sound_speed_field_2d(self, p, rho, cs)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    class(field_2d_t), intent(in) :: p   !< pressure
    class(field_2d_t), intent(in) :: rho !< density
    class(field_2d_t), intent(inout) :: cs !< sound speed
    
    integer(ik) :: i, j
    real(rk) :: sqrt_gamma

    sqrt_gamma = sqrt(self%gamma)

    do j = rho%jlo, rho%jhi
      !$omp simd
      do i = rho%ilo, rho%ihi
        cs%data(i,j) = sqrt_gamma * sqrt(p%data(i,j)/ rho%data(i,j))
      enddo
      !$omp end simd
    enddo

    cs%name = 'cs'
  endsubroutine sound_speed_field_2d

  subroutine temperature_field_2d(self, p, rho, t)
    !< Calculate sound speed

    class(eos_t), intent(in) :: self
    class(field_2d_t), intent(in) :: p
    class(field_2d_t), intent(in) :: rho
    class(field_2d_t), intent(inout) :: t

    t = p / (rho * self%R)
    t%name = 'T'

  endsubroutine temperature_field_2d

  subroutine conserved_to_primitive_field_2d(self, rho, rho_u, rho_v, rho_E, u, v, p)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    class(field_2d_t), intent(in) :: rho   !< density
    class(field_2d_t), intent(in) :: rho_u !< density * x-velocity
    class(field_2d_t), intent(in) :: rho_v !< density * y-velocity
    class(field_2d_t), intent(in) :: rho_E !< density * total energy
    class(field_2d_t), intent(inout) :: u    !< x-velocity
    class(field_2d_t), intent(inout) :: v    !< y-velocity
    class(field_2d_t), intent(inout) :: p    !< pressure

    ! Locals
    real(rk) :: gamma_m_one
    integer(ik) :: i, j

    gamma_m_one = self%gamma - 1.0_rk

    ! u = rho_u / rho
    ! v = rho_v / rho
    
    do j = rho%jlo, rho%jhi
      !$omp simd
      do i = rho%ilo, rho%ihi
        u%data(i,j) = rho_u%data(i,j) / rho%data(i,j)
        v%data(i,j) = rho_v%data(i,j) / rho%data(i,j)
      enddo
      !$omp end simd
    enddo

    ! p = rho * gamma_m_one * ((rho_E / rho) - (0.5_rk * (u * u + v * v)))
    do j = rho%jlo, rho%jhi
      !$omp simd
      do i = rho%ilo, rho%ihi
        p%data(i,j) = rho%data(i,j) * gamma_m_one * ((rho_E%data(i,j) / rho%data(i,j)) - &
                      (0.5_rk * (u%data(i,j) * u%data(i,j) + v%data(i,j) * v%data(i,j))))
      enddo
      !$omp end simd
    enddo

    u%name = 'u'
    v%name = 'v'
    p%name = 'p'

  endsubroutine conserved_to_primitive_field_2d

  subroutine primitive_to_conserved_field_2d(self, rho, u, v, p, rho_u, rho_v, rho_E)
    !< Convert conserved quantities [rho, rho u, rho v, rho E] into primitive [rho, u, v, p]. This
    !< is in the EOS class due to requirement of converting energy into pressure

    class(eos_t), intent(in) :: self
    class(field_2d_t), intent(in) :: rho
    class(field_2d_t), intent(in) :: u
    class(field_2d_t), intent(in) :: v
    class(field_2d_t), intent(in) :: p
    class(field_2d_t), intent(inout) :: rho_u
    class(field_2d_t), intent(inout) :: rho_v
    class(field_2d_t), intent(inout) :: rho_E

    ! Locals
    real(rk) :: inv_gamma_m_one
    integer(ik) :: i, j

    inv_gamma_m_one = 1.0_rk  / (self%gamma - 1.0_rk)

    do j = rho%jlo, rho%jhi
      !$omp simd
      do i = rho%ilo, rho%ihi
        rho_u%data(i, j) = u%data(i, j) * rho%data(i, j)
        rho_v%data(i, j) = v%data(i, j) * rho%data(i, j)
        rho_E%data(i, j) = (p%data(i, j) * inv_gamma_m_one) + rho%data(i, j) * \
                            0.5_rk * (u%data(i, j) * u%data(i, j) + v%data(i, j) * v%data(i, j))        
      enddo
      !$omp end simd
    enddo
    
    rho_v%name = 'rho_v'
    rho_E%name = 'rho_E'
    rho_u%name = 'rho_u'
  endsubroutine primitive_to_conserved_field_2d

endmodule mod_eos
