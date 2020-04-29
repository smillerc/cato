module mod_flux_tensor
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_vector, only: vector_t
  use mod_eos, only: eos

  implicit none

  private
  public :: flux_tensor_t, operator(.dot.) !, operator(*), operator(-), operator(+)

contains

  pure subroutine get_fluxes(primitive_variables, fluxes, lbounds)
    !< Implementation of the flux tensor (H) construction. This requires the primitive variables
    real(rk), dimension(:, :, :), intent(in) :: primitive_variables !< [rho, u, v, p]
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(:, :, :, :), allocatable, intent(out) :: fluxes

    real(rk) :: total_energy !< total energy

    ! The flux tensor is H = Fi + Gj
    associate(rho=>primitive_variables(1), &
              u=>primitive_variables(2), &
              v=>primitive_variables(3), &
              p=>primitive_variables(4), &
              H=>flux_tensor%state)

      total_energy = eos%calculate_total_energy(pressure=p, density=rho, x_velocity=u, y_velocity=v)

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

  end subroutine

  pure function flux_dot_vector(lhs, vec) result(output)
    !< Implementation of the dot product between a flux tensor and a vector, e.g. H . v
    type(flux_tensor_t), intent(in) :: lhs  !< Left-hand side of the dot product
    real(rk), dimension(2), intent(in) :: vec  !< Right-hand side of the dot product
    real(rk), dimension(4) :: output

    ! The flux tensor is H = Fi + Gj
    output = lhs%state(:, 1) * vec(1) + lhs%state(:, 2) * vec(2)
  end function

end module mod_flux_tensor
