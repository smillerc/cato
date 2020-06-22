module mod_ausm_fluid
  !> Summary: Provide the basis for the fluid class that uses the family of AUSM solvers
  !> Date: 06/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data
  use mod_nondimensionalization, only: scale_factors_set, rho_0, v_0, p_0, e_0, t_0
  use mod_base_fluid, only: base_fluid_t

  implicit none

  private
  public :: base_fluid_t, new_fluid

  logical, parameter :: filter_small_mach = .false.

  type, extends(base_fluid_t) :: ausm_fluid_t
  contains
    procedure :: force_finalization
    procedure, public :: t => time_derivative
    procedure, nopass, private :: flux_edges
    final :: finalize
  end type ausm_fluid_t

contains

  subroutine time_derivative()
  end subroutine time_derivative

  subroutine flux_edges()
  end subroutine flux_edges

  subroutine force_finalization(self)
    type(ausm_fluid_t), intent(inout) :: self
    call self%finalize()
  end subroutine force_finalization

  subroutine finalize(self)
    type(ausm_fluid_t), intent(inout) :: self

    call debug_print('Running fluid_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%rho)) deallocate(self%rho)
    if(allocated(self%u)) deallocate(self%u)
    if(allocated(self%v)) deallocate(self%v)
    if(allocated(self%p)) deallocate(self%p)
    if(allocated(self%rho_u)) deallocate(self%rho_u)
    if(allocated(self%rho_v)) deallocate(self%rho_v)
    if(allocated(self%rho_E)) deallocate(self%rho_E)
    if(allocated(self%cs)) deallocate(self%cs)
    if(allocated(self%mach_u)) deallocate(self%mach_u)
    if(allocated(self%mach_v)) deallocate(self%mach_v)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine finalize

end module mod_ausm_fluid
