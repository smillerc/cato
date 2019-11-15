module mod_fvleg

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_reconstruction_factory, only: reconstruction_factory_t
  use mod_local_evo_operator, only: local_evo_operator_t
  use mod_second_order_reconstruction, only: second_order_reconstruction_t
  use mod_first_order_reconstruction, only: first_order_reconstruction_t

  implicit none
  private
  public :: fvleg_t

  type, extends(finite_volume_scheme_t) :: fvleg_t
    !< Implementation of the finite volume local evolution Galerkin (fvleg) scheme type
    class(local_evo_operator_t), allocatable :: local_evolution_operator
  contains
    procedure, public :: reconstruct => reconstruct_fvleg
    procedure, public :: eval_fluxes => eval_fvleg_fluxes
    procedure, public :: t => time_derivative

    procedure, public :: add => add_fvleg
    procedure, public :: multiply => multiply_fvleg
    procedure, public :: assign => assign_fvleg
  end type

contains

  function constructor(input) result(scheme)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(input_t), intent(in) :: input
    class(fvleg_t), allocatable :: scheme
    type(reconstruction_factory_t) :: recon_factory

    recon_factory = reconstruction_factory_t()
    allocate(fvleg_t :: scheme)

    scheme%conserved_vars = conserved_vars_t(input%ni, input%nj)
    scheme%reconstruction_operator = recon_factory%create_reconstruction(input)
  end function

  subroutine reconstruct_fvleg(self)
    class(fvleg_t), intent(inout) :: self
  end subroutine

  pure function get_fluxes(density, x_velocity, y_velocity, pressure, gamma) result(H)
    !< Create the H tensor based on the given state variables
    real(rk), intent(in) :: density
    real(rk), intent(in) :: x_velocity
    real(rk), intent(in) :: y_velocity
    real(rk), intent(in) :: pressure
    real(rk), intent(in) :: gamma
    real(rk), dimension(2, 4) :: H

    real(rk) :: total_energy, internal_energy

    internal_energy = (pressure / density) / (gamma - 1)
    total_energy = density * (internal_energy + 0.5_rk * (x_velocity**2 + y_velocity**2))

    ! The flux tensor is H = Fi + Gj
    associate(rho=>density, u=>x_velocity, v=>y_velocity, p=>pressure, E=>total_energy)
      ! F
      H(1) = [rho * u, rho * u**2 + p, rho * u * v,(E + p) * u]

      ! G
      H(2) = [rho * v, rho * u * v, rho * v**2 + p,(E + p) * v]
    end associate
  end function

  pure function eval_fvleg_fluxes(self) result(H)

    !< Evaluate the fluxes in the equation dU/dt + dF/dx + dG/dy = 0, where H = Fi + Gj
    class(fvleg_t), intent(inout) :: self
    real(rk), dimension(2, 4) :: H
    real(rk), dimension(4) :: F, G

  end function

  function time_derivative(self) result(dU_dt)
    !< Implementation of dU/dt
    class(fvleg_t), intent(in) :: self
    class(integrand_t), allocatable :: dU_dt
    type(fvleg_t), allocatable :: local_dU_dt

    ! allocate(local_dU_dt)

    ! select type(self%reconstruction_operator)
    ! class is (first_order_reconstruction_t)
    !   dU_dt = self%first_order_edge_flux()
    ! class is (second_order_reconstruction_t)
    !   dU_dt = self%second_order_edge_flux()
    ! end select

  end function time_derivative

  pure function first_order_edge_flux(self) result(H)
    class(fvleg_t), intent(inout) :: self
    real(rk), dimension(2, 4) :: H

  end function first_order_edge_flux

  pure function second_order_edge_flux(self) result(H)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    class(fvleg_t), intent(inout) :: self
    real(rk), dimension(2, 4) :: H

    integer(ik) :: i, j, k

    do j = self%grid%jlo, self%grid%jhi
      do i = self%grid%ilo, self%grid%ihi

        r_omega_u_bar => self%reconstruction_operator
        edge_flux = 0.0_rk

        do k = 1, 4  ! edge
          associate(H=>self%get_fluxes, E0=>self%local_evolution_operator%evolve, &
                    n_hat=>self%grid%cell_edge_norm_vectors(i, j, k, :, :), &
                    delta_l=>self%grid%cell_edge_lengths(i, j, k))
            edge_flux = edge_flux + &
                        H(E0(r_omega_u_bar, edge=k, loc='k,1'))) + &
                        4 * H(E0(r_omega_u_bar, edge=k, loc='k,c'))) + &
                        H(E0(r_omega_u_bar, edge=k, loc='k,2'))) .dot.n_hat * delta_l / 6.0_rk
          end associate
        end do
      end do
    end do

  end function second_order_edge_flux

  function add_fvleg(lhs, rhs) result(operator_result)
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs
    class(integrand_t), allocatable :: operator_result
  end function add_fvleg

  function multiply_fvleg(lhs, rhs) result(operator_result)
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), allocatable :: operator_result
    real(rk), intent(in) :: rhs
  end function multiply_fvleg

  subroutine assign_fvleg(lhs, rhs)
    class(fvleg_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs
  end subroutine assign_fvleg

end module mod_fvleg
