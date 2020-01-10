module mod_pressure_input_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t
  use linear_interpolation_module, only: linear_interp_1d
  use mod_eos, only: eos
  implicit none

  private
  public :: pressure_input_bc_t, pressure_input_bc_constructor

  type, extends(boundary_condition_t) :: pressure_input_bc_t
    logical :: constant_pressure = .true.
    real(rk) :: pressure_input = 1.0_rk

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_pressure_input
    character(len=100) :: input_filename
  contains
    procedure, private :: read_pressure_input
    procedure, public :: apply_conserved_var_bc => apply_pressure_input_conserved_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_pressure_input_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc => apply_pressure_input_cell_gradient_bc
    procedure, public :: copy => copy_pressure_input_bc
  end type
contains

  function pressure_input_bc_constructor(location, input) result(bc)
    type(pressure_input_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input

    allocate(bc)
    bc%name = 'pressure_input'
    bc%location = location
  end function pressure_input_bc_constructor

  subroutine read_pressure_input(self)
    !< Read the pressure input file and initialize the linear iterpolated object. This
    !< makes it easy to get the input pressure at any given time
    class(pressure_input_bc_t), intent(inout) :: self

    real(rk), dimension(:), allocatable :: time_sec
    real(rk), dimension(:), allocatable :: pressure_barye
    integer(ik) :: input_unit, io_status, interp_status
    integer(ik) :: line, nlines

    open(newunit=input_unit, file=trim(self%input_filename))
    nlines = 0
    do
      read(input_unit, *, iostat=io_status)
      if(io_status /= 0) exit
      nlines = nlines + 1
    end do
    close(input_unit)

    open(newunit=input_unit, file=trim(self%input_filename))
    allocate(time_sec(nlines))
    allocate(pressure_barye(nlines))

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time_sec(line), pressure_barye(line)
      if(io_status /= 0) exit
    end do
    close(input_unit)

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_pressure_input%initialize(time_sec, pressure_barye, interp_status)
    if(interp_status /= 0) error stop "Error initializing pressure_input_bc_t%temporal_pressure_input"

    deallocate(time_sec)
    deallocate(pressure_barye)
  end subroutine

  subroutine copy_pressure_input_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(pressure_input_bc_t), intent(inout) :: out_bc
  end subroutine

  subroutine apply_pressure_input_conserved_var_bc(self, conserved_vars)
    !< Apply pressure_input boundary conditions to the conserved state vector field
    class(pressure_input_bc_t), intent(inout) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    !< ((rho, u ,v, p), i, j); Conserved variables for each cell

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index

    integer(ik) :: interp_stat
    real(rk) :: boundary_pressure

    left_ghost = lbound(conserved_vars, dim=2)
    right_ghost = ubound(conserved_vars, dim=2)
    bottom_ghost = lbound(conserved_vars, dim=3)
    top_ghost = ubound(conserved_vars, dim=3)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    ! Since we know the pressure/density in the fluid edge cell,
    ! and we also know what we want the pressure to be on the boundary (ghost)
    ! We use the isentropic relation rho2/rho2 = (P2/P1)^(gamma-1) to find the density
    ! in the ghost zone
    !           |
    !  Fluid    |   Ghost
    !  P1,      |    P2
    !  rho_1    |  rho2
    !           |

    if(self%constant_pressure) then
      boundary_pressure = self%pressure_input
    else
      call self%temporal_pressure_input%evaluate(x=self%get_time(), f=boundary_pressure, istat=interp_stat)
 if(interp_stat /= 0) error stop "Unable to interpolate pressure within pressure_input_bc_t%apply_pressure_input_conserved_var_bc()"
    end if

    select case(self%location)
    case('+x')

      conserved_vars(1, right_ghost, :) = eos%calc_density_from_isentropic_press(p_1=conserved_vars(4, right, :), &
                                                                                 rho_1=conserved_vars(1, right, :), &
                                                                                 p_2=boundary_pressure)

      conserved_vars(2, right_ghost, :) = 0.0_rk ! or should it be the same as the interior cell?
      conserved_vars(3, right_ghost, :) = 0.0_rk
      conserved_vars(4, right_ghost, :) = boundary_pressure
    case('-x')
      ! conserved_vars(:, left_ghost, :) = conserved_vars(:, right, bottom)
    case('+y')
      ! conserved_vars(:, :, top_ghost) = conserved_vars(:, left:right, bottom)
    case('-y')
      ! conserved_vars(:, :, bottom_ghost) = conserved_vars(:, left:right, top)
    case default
      error stop "Unsupported location to apply the bc at in pressure_input_bc_t%apply_pressure_input_cell_gradient_bc()"
    end select

  end subroutine apply_pressure_input_conserved_var_bc

  subroutine apply_pressure_input_reconstructed_state_bc(self, reconstructed_state)
    !< Apply pressure_input boundary conditions to the reconstructed state vector field

    class(pressure_input_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

    ! integer(ik) :: left         !< Min i real cell index
    ! integer(ik) :: right        !< Max i real cell index
    ! integer(ik) :: bottom       !< Min j real cell index
    ! integer(ik) :: top          !< Max j real cell index
    ! integer(ik) :: left_ghost   !< Min i ghost cell index
    ! integer(ik) :: right_ghost  !< Max i ghost cell index
    ! integer(ik) :: bottom_ghost !< Min j ghost cell index
    ! integer(ik) :: top_ghost    !< Max j ghost cell index

    ! left_ghost = lbound(reconstructed_state, dim=4)
    ! right_ghost = ubound(reconstructed_state, dim=4)
    ! bottom_ghost = lbound(reconstructed_state, dim=5)
    ! top_ghost = ubound(reconstructed_state, dim=5)
    ! left = left_ghost + 1
    ! right = right_ghost - 1
    ! bottom = bottom_ghost + 1
    ! top = top_ghost - 1

    ! select case(self%location)
    ! case('+x')
    !   reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)
    !   reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)
    !   reconstructed_state(:, :, :, right_ghost, bottom:top) = reconstructed_state(:, :, :, left, bottom:top)
    ! case('-x')
    !   reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)
    !   reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)
    !   reconstructed_state(:, :, :, left_ghost, bottom:top) = reconstructed_state(:, :, :, right, bottom:top)
    ! case('+y')
    !   reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)
    !   reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)
    !   reconstructed_state(:, :, :, left:right, top_ghost) = reconstructed_state(:, :, :, left:right, bottom)
    ! case('-y')
    !   reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)
    !   reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)
    !   reconstructed_state(:, :, :, left:right, bottom_ghost) = reconstructed_state(:, :, :, left:right, top)
    ! case default
    !   error stop "Unsupported location to apply the bc at in pressure_input_bc_t%apply_pressure_input_reconstructed_state_bc()"
    ! end select

  end subroutine apply_pressure_input_reconstructed_state_bc

  subroutine apply_pressure_input_cell_gradient_bc(self, cell_gradient)
    !< Apply pressure_input boundary conditions to the reconstructed state vector field

    class(pressure_input_bc_t), intent(in) :: self
    real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

    ! integer(ik) :: left         !< Min i real cell index
    ! integer(ik) :: right        !< Max i real cell index
    ! integer(ik) :: bottom       !< Min j real cell index
    ! integer(ik) :: top          !< Max j real cell index
    ! integer(ik) :: left_ghost   !< Min i ghost cell index
    ! integer(ik) :: right_ghost  !< Max i ghost cell index
    ! integer(ik) :: bottom_ghost !< Min j ghost cell index
    ! integer(ik) :: top_ghost    !< Max j ghost cell index

    ! left_ghost = lbound(cell_gradient, dim=3)
    ! right_ghost = ubound(cell_gradient, dim=3)
    ! bottom_ghost = lbound(cell_gradient, dim=4)
    ! top_ghost = ubound(cell_gradient, dim=4)
    ! left = left_ghost + 1
    ! right = right_ghost - 1
    ! bottom = bottom_ghost + 1
    ! top = top_ghost - 1

    ! select case(self%location)
    ! case('+x')
    !   cell_gradient(:, :, right_ghost, top_ghost) = cell_gradient(:, :, left, bottom)
    !   cell_gradient(:, :, right_ghost, bottom_ghost) = cell_gradient(:, :, left, top)
    !   cell_gradient(:, :, right_ghost, bottom:top) = cell_gradient(:, :, left, bottom:top)
    ! case('-x')
    !   cell_gradient(:, :, left_ghost, top_ghost) = cell_gradient(:, :, right, bottom)
    !   cell_gradient(:, :, left_ghost, bottom_ghost) = cell_gradient(:, :, right, top)
    !   cell_gradient(:, :, left_ghost, bottom:top) = cell_gradient(:, :, right, bottom:top)
    ! case('+y')
    !   cell_gradient(:, :, left_ghost, top_ghost) = cell_gradient(:, :, right, bottom)
    !   cell_gradient(:, :, right_ghost, top_ghost) = cell_gradient(:, :, left, bottom)
    !   cell_gradient(:, :, left:right, top_ghost) = cell_gradient(:, :, left:right, bottom)
    ! case('-y')
    !   cell_gradient(:, :, left_ghost, bottom_ghost) = cell_gradient(:, :, right, top)
    !   cell_gradient(:, :, right_ghost, bottom_ghost) = cell_gradient(:, :, left, top)
    !   cell_gradient(:, :, left:right, bottom_ghost) = cell_gradient(:, :, left:right, top)
    ! case default
    !   error stop "Unsupported location to apply the bc at in pressure_input_bc_t%apply_pressure_input_cell_gradient_bc()"
    ! end select

  end subroutine apply_pressure_input_cell_gradient_bc
end module mod_pressure_input_bc
