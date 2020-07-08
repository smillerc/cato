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

module mod_flux_tensor
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_vector, only: vector_t
  use mod_eos, only: eos

  implicit none

  private
  public :: flux_tensor_t, operator(.dot.) !, operator(*), operator(-), operator(+)

  type flux_tensor_t
    real(rk), dimension(4, 2) :: state  !< ((i,j), conserved_quantites)
  contains
    procedure, pass(lhs), private :: assign_from_flux
    procedure, pass(flux) :: real_mul_flux
    procedure, private :: flux_mul_real
    procedure, private :: flux_div_real
    procedure, private :: flux_add_real
    procedure, private :: flux_sub_real
    procedure, private :: flux_add_flux
    procedure, private :: flux_sub_flux

    generic :: assignment(=) => assign_from_flux
    generic :: operator(+) => flux_add_flux, flux_add_real
    generic :: operator(-) => flux_sub_flux, flux_sub_real
    generic :: operator(*) => flux_mul_real, real_mul_flux
    generic :: operator(/) => flux_div_real
    ! generic :: operator(.dot.) => flux_dot_vector
  end type

  interface flux_tensor_t
    module procedure :: constructor
  end interface

  interface operator(.dot.)
    procedure :: flux_dot_vector
  end interface

contains

  type(flux_tensor_t) pure function constructor(primitive_variables) result(flux_tensor)
    !< Implementation of the flux tensor (H) construction. This requires the primitive variables
    real(rk), dimension(4), intent(in) :: primitive_variables !< [rho, u, v, p]

    real(rk) :: total_energy !< total energy

    ! The flux tensor is H = Fi + Gj
    associate(rho=>primitive_variables(1), &
              u=>primitive_variables(2), &
              v=>primitive_variables(3), &
              p=>primitive_variables(4), &
              H=>flux_tensor%state)

      call eos%total_energy(p=p, rho=rho, u=u, v=v, E=total_energy)

      ! F
      H(:, 1) = [rho * u, &
                 rho * u**2 + p, &
                 rho * u * v, &
                 u * (rho * total_energy + p)]

      ! G
      H(:, 2) = [rho * v, &
                 rho * u * v, &
                 rho * v**2 + p, &
                 v * (rho * total_energy + p)]
    end associate

  end function

  pure subroutine assign_from_flux(lhs, rhs)
    class(flux_tensor_t), intent(inout) :: lhs
    class(flux_tensor_t), intent(in) :: rhs
    lhs%state = rhs%state
  end subroutine assign_from_flux

  type(flux_tensor_t) pure function flux_add_flux(self, other_flux) result(output)
    class(flux_tensor_t), intent(in) :: self
    class(flux_tensor_t), intent(in) :: other_flux
    output%state = self%state + other_flux%state
  end function

  type(flux_tensor_t) pure function flux_add_real(self, scalar) result(output)
    class(flux_tensor_t), intent(in) :: self
    real(rk), intent(in) :: scalar
    output%state = self%state + scalar
  end function

  type(flux_tensor_t) pure function flux_sub_flux(self, other_flux) result(output)
    class(flux_tensor_t), intent(in) :: self
    class(flux_tensor_t), intent(in) :: other_flux
    output%state = self%state - other_flux%state
  end function

  type(flux_tensor_t) pure function flux_sub_real(self, scalar) result(output)
    class(flux_tensor_t), intent(in) :: self
    real(rk), intent(in) :: scalar
    output%state = self%state - scalar
  end function

  type(flux_tensor_t) pure function flux_mul_real(self, scalar) result(output)
    class(flux_tensor_t), intent(in) :: self
    real(rk), intent(in) :: scalar
    output%state = self%state * scalar
  end function

  type(flux_tensor_t) pure function real_mul_flux(scalar, flux) result(output)
    class(flux_tensor_t), intent(in) :: flux
    real(rk), intent(in) :: scalar
    output%state = flux%state * scalar
  end function

  type(flux_tensor_t) pure function flux_div_real(self, scalar) result(output)
    class(flux_tensor_t), intent(in) :: self
    real(rk), intent(in) :: scalar
    output%state = self%state / scalar
  end function

  pure function flux_dot_vector(lhs, vec) result(output)
    !< Implementation of the dot product between a flux tensor and a vector, e.g. H . v
    type(flux_tensor_t), intent(in) :: lhs  !< Left-hand side of the dot product
    real(rk), dimension(2), intent(in) :: vec  !< Right-hand side of the dot product
    real(rk), dimension(4) :: output

    ! The flux tensor is H = Fi + Gj
    output = lhs%state(:, 1) * vec(1) + lhs%state(:, 2) * vec(2)
  end function

end module mod_flux_tensor
