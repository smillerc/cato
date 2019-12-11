module mod_flux_tensor
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_vector, only: vector_t
  use mod_eos, only: eos

  implicit none

  private
  public :: flux_tensor_t, operator(.dot.) !, operator(*), operator(-), operator(+)

  type flux_tensor_t
    real(rk), dimension(2, 4) :: state
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

  type(flux_tensor_t) function constructor(primitive_variables) result(flux_tensor)
    !< Implementation of the flux tensor (H) construction. This requires the primitive variables
    real(rk), dimension(4), intent(in) :: primitive_variables !< [rho, u, v, p]

    real(rk) :: total_energy, internal_energy

    associate(rho=>primitive_variables(1), &
              u=>primitive_variables(2), v=>primitive_variables(3), &
              p=>primitive_variables(4), gamma=>eos%get_gamma())

      print*, 'H:', p, rho, gamma
      internal_energy = (p / rho) / (gamma - 1)
      total_energy = rho * (internal_energy + 0.5_rk * (u**2 + v**2))
    end associate

    ! The flux tensor is H = Fi + Gj
    associate(rho=>primitive_variables(1), &
              u=>primitive_variables(2), v=>primitive_variables(3), &
              p=>primitive_variables(4), &
              E=>total_energy, H=>flux_tensor%state)
      ! F
      H(1, :) = [rho * u, rho * u**2 + p, rho * u * v,(E + p) * u]

      ! G
      H(2, :) = [rho * v, rho * u * v, rho * v**2 + p,(E + p) * v]
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
    output(1) = dot_product(lhs%state(:, 1), vec)
    output(2) = dot_product(lhs%state(:, 2), vec)
    output(3) = dot_product(lhs%state(:, 3), vec)
    output(4) = dot_product(lhs%state(:, 4), vec)
  end function

  ! type(flux_tensor_t) pure function vector_dot_flux(vec, rhs) result(H)
  !   !< Implementation of the dot product between a flux tensor and a vector, e.g. v . H
  !   type(vector_t), intent(in) :: vec  !< Right-hand side of the dot product
  !   type(flux_tensor_t), intent(in) :: rhs  !< left-hand side of the dot product
  !   ! H%state = dot_product([vec%x, vec%y], rhs%state)
  ! end function

  ! type(flux_tensor_t) pure function flux_add_flux(lhs, rhs) result(flux_tensor)
  !   type(flux_tensor_t), intent(in) :: lhs  !< Left-hand side
  !   type(flux_tensor_t), intent(in) :: rhs  !< Right-hand side
  !   flux_tensor%state = lhs%state + rhs%state
  ! end function

  ! type(flux_tensor_t) pure function flux_minus_flux(lhs, rhs) result(flux_tensor)
  !   type(flux_tensor_t), intent(in) :: lhs  !< Left-hand side
  !   type(flux_tensor_t), intent(in) :: rhs  !< Right-hand side
  !   flux_tensor%state = lhs%state - rhs%state
  ! end function

  ! type(flux_tensor_t) pure function flux_multiply_scalar_rhs(scalar, rhs) result(flux_tensor)
  !   type(flux_tensor_t), intent(in) :: rhs  !< Right-hand side
  !   real(rk), intent(in) :: scalar
  !   flux_tensor%state = rhs%state * scalar
  ! end function

  ! type(flux_tensor_t) pure function flux_multiply_scalar_lhs(lhs, scalar) result(flux_tensor)
  !   type(flux_tensor_t), intent(in) :: lhs  !< Left-hand side
  !   real(rk), intent(in) :: scalar
  !   flux_tensor%state = lhs%state * scalar
  ! end function

end module mod_flux_tensor
