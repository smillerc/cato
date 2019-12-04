module mod_periodic_bc

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: periodic_bc_t

  type, extends(boundary_condition_t) :: periodic_bc_t
  contains
    procedure, public :: initialize => init_periodic_bc
    procedure, public :: apply_bc => apply_periodic_bc
    procedure, private :: apply_plus_x
    procedure, private :: apply_plus_y
    procedure, private :: apply_minus_x
    procedure, private :: apply_minus_y

  end type periodic_bc_t

contains

  subroutine init_periodic_bc(self, location)
    class(periodic_bc_t), intent(inout) :: self
    ! class(input_t), intent(in) :: input
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)

    self%name = 'periodic'
    self%location = location
  end subroutine init_periodic_bc

  subroutine apply_periodic_bc(self, conserved_vars, reconstructed_state)
    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, :, :), intent(inout) :: reconstructed_state

    select case(self%location)
    case('+x')
      call self%apply_plus_x(conserved_vars, reconstructed_state)
    case('-x')
      call self%apply_minus_x(conserved_vars, reconstructed_state)
    case('+y')
      call self%apply_plus_y(conserved_vars, reconstructed_state)
    case('-y')
      call self%apply_minus_y(conserved_vars, reconstructed_state)
    case default
      error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_bc()"
    end select

  end subroutine apply_periodic_bc

  subroutine apply_plus_x(self, conserved_vars, reconstructed_state)
    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, :, :), intent(inout) :: reconstructed_state

    integer(ik) :: imin     !< Min i index (w/in real domain)
    integer(ik) :: imax     !< Max i index (w/in real domain)
    integer(ik) :: jmin     !< Min j index (w/in real domain)
    integer(ik) :: jmax     !< Max j index (w/in real domain)
    integer(ik) :: imin_bc  !< Min i index (including boundary)
    integer(ik) :: imax_bc  !< Max i index (including boundary)
    integer(ik) :: jmin_bc  !< Min j index (including boundary)
    integer(ik) :: jmax_bc  !< Max j index (including boundary)

    ! Top Right Corner
    ! Bottom Right Corner
    ! Middle Right Section
  end subroutine

  subroutine apply_plus_y(self, conserved_vars, reconstructed_state)
    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, :, :), intent(inout) :: reconstructed_state

    integer(ik) :: imin     !< Min i index (w/in real domain)
    integer(ik) :: imax     !< Max i index (w/in real domain)
    integer(ik) :: jmin     !< Min j index (w/in real domain)
    integer(ik) :: jmax     !< Max j index (w/in real domain)
    integer(ik) :: imin_bc  !< Min i index (including boundary)
    integer(ik) :: imax_bc  !< Max i index (including boundary)
    integer(ik) :: jmin_bc  !< Min j index (including boundary)
    integer(ik) :: jmax_bc  !< Max j index (including boundary)

    ! Top Left Corner
    ! Top Right Corner
    ! Top Middle Section
  end subroutine

  subroutine apply_minus_x(self, conserved_vars, reconstructed_state)
    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, :, :), intent(inout) :: reconstructed_state

    integer(ik) :: imin     !< Min i index (w/in real domain)
    integer(ik) :: imax     !< Max i index (w/in real domain)
    integer(ik) :: jmin     !< Min j index (w/in real domain)
    integer(ik) :: jmax     !< Max j index (w/in real domain)
    integer(ik) :: imin_bc  !< Min i index (including boundary)
    integer(ik) :: imax_bc  !< Max i index (including boundary)
    integer(ik) :: jmin_bc  !< Min j index (including boundary)
    integer(ik) :: jmax_bc  !< Max j index (including boundary)

    ! Top Left Corner
    ! Bottom Left Corner
    ! Middle Left Section
  end subroutine

  subroutine apply_minus_y(self, conserved_vars, reconstructed_state)
    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, :, :), intent(inout) :: reconstructed_state

    integer(ik) :: imin     !< Min i index (w/in real domain)
    integer(ik) :: imax     !< Max i index (w/in real domain)
    integer(ik) :: jmin     !< Min j index (w/in real domain)
    integer(ik) :: jmax     !< Max j index (w/in real domain)
    integer(ik) :: imin_bc  !< Min i index (including boundary)
    integer(ik) :: imax_bc  !< Max i index (including boundary)
    integer(ik) :: jmin_bc  !< Min j index (including boundary)
    integer(ik) :: jmax_bc  !< Max j index (including boundary)

    ! Bottom Left Corner
    ! Bottom Right Corner
    ! Bottom Middle Section
  end subroutine
end module mod_periodic_bc
