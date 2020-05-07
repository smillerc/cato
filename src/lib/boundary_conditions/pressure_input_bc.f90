module mod_pressure_input_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_nondimensionalization, only: p_0, rho_0, t_0
  use mod_input, only: input_t
  use mod_inlet_outlet, only: subsonic_inlet, subsonic_outlet, supersonic_inlet, supersonic_outlet
  use linear_interpolation_module, only: linear_interp_1d
  use mod_eos, only: eos
  implicit none

  private
  public :: pressure_input_bc_t, pressure_input_bc_constructor

  type, extends(boundary_condition_t) :: pressure_input_bc_t
    logical :: constant_pressure = .false.
    real(rk) :: pressure_input = 0.0_rk
    real(rk) :: temperature_input = 273.0_rk !K
    real(rk) :: density_input = 1e-3_rk !K
    real(rk) :: scale_factor = 1.0_rk !< scaling factor (e.g. 1x, 2x) to scale the input up/down

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_pressure_input
    character(len=100) :: input_filename

    real(rk), dimension(:), allocatable :: edge_rho !< boundary (ghost/edge/etc) density
    real(rk), dimension(:), allocatable :: edge_u   !< boundary (ghost/edge/etc) x-velocity
    real(rk), dimension(:), allocatable :: edge_v   !< boundary (ghost/edge/etc) y-velocity
    real(rk), dimension(:), allocatable :: edge_p   !< boundary (ghost/edge/etc) pressure

  contains
    procedure, private :: read_pressure_input
    procedure, private :: get_desired_pressure
    procedure, public :: apply_primitive_var_bc => apply_pressure_input_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_pressure_input_reconstructed_state_bc
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
    bc%density_input = 1e-3_rk / rho_0
    bc%constant_pressure = input%apply_constant_bc_pressure
    bc%scale_factor = input%bc_pressure_scale_factor
    if(.not. bc%constant_pressure) then
      bc%input_filename = trim(input%bc_pressure_input_file)
      call bc%read_pressure_input()
    else
      bc%pressure_input = input%constant_bc_pressure_value / p_0 ! non-dimensionalize
    end if

  end function pressure_input_bc_constructor

  subroutine read_pressure_input(self)
    !< Read the pressure input file and initialize the linear iterpolated object. This
    !< makes it easy to get the input pressure at any given time
    class(pressure_input_bc_t), intent(inout) :: self

    real(rk), dimension(:), allocatable :: time_sec
    real(rk), dimension(:), allocatable :: pressure_barye
    character(len=300) :: line_buffer
    logical :: has_header_line
    logical :: file_exists
    integer(ik) :: input_unit, io_status, interp_status
    integer(ik) :: line, nlines

    file_exists = .false.
    has_header_line = .false.

    inquire(file=trim(self%input_filename), exist=file_exists)

    if(.not. file_exists) then
      error stop 'Error in pressure_input_bc_t%read_pressure_input(); pressure input file not found, exiting...'
    end if

    open(newunit=input_unit, file=trim(self%input_filename))
    nlines = 0
    do
      read(input_unit, *, iostat=io_status) line_buffer
      if(nlines == 0) then
        line_buffer = adjustl(line_buffer)
        if(line_buffer(1:1) == '#') has_header_line = .true.
      end if

      if(io_status /= 0) exit
      nlines = nlines + 1
    end do
    close(input_unit)

    open(newunit=input_unit, file=trim(self%input_filename))
    if(has_header_line) then
      read(input_unit, *, iostat=io_status) ! skip the first line
      nlines = nlines - 1
    end if

    allocate(time_sec(nlines))
    allocate(pressure_barye(nlines))
    time_sec = 0.0_rk
    pressure_barye = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time_sec(line), pressure_barye(line)
      if(io_status /= 0) exit
    end do
    close(input_unit)

    self%max_time = maxval(time_sec)

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_pressure_input%initialize(time_sec, pressure_barye, interp_status)
    if(interp_status /= 0) then
      write(*, *) time_sec
      write(*, *) pressure_barye
      write(*, *) interp_status
      error stop "Error initializing pressure_input_bc_t%temporal_pressure_input"
    end if

    deallocate(time_sec)
    deallocate(pressure_barye)
  end subroutine

  subroutine copy_pressure_input_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(pressure_input_bc_t), intent(inout) :: out_bc

    call debug_print('Running pressure_input_bc_t%copy_pressure_input_bc()', __FILE__, __LINE__)

    out_bc%name = in_bc%name
    out_bc%location = in_bc%location
    select type(in_bc)
    class is(pressure_input_bc_t)
      out_bc%constant_pressure = in_bc%constant_pressure
      out_bc%pressure_input = in_bc%pressure_input
      out_bc%temporal_pressure_input = in_bc%temporal_pressure_input
      out_bc%input_filename = in_bc%input_filename
    class default
      error stop 'pressure_input_bc_t%copy_pressure_input_bc: unsupported in_bc class'
    end select

  end subroutine

  real(rk) function get_desired_pressure(self) result(desired_pressure)
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik) :: interp_stat

    if(self%constant_pressure) then
      desired_pressure = self%pressure_input
    else
      call self%temporal_pressure_input%evaluate(x=self%get_time() * t_0, f=desired_pressure, istat=interp_stat)
      if(interp_stat /= 0) then
        error stop "Unable to interpolate pressure within pressure_input_bc_t%get_desired_pressure()"
      end if
    end if

    ! apply scale factor (for user convienence) and non-dimensional scale factor
    desired_pressure = desired_pressure * self%scale_factor / p_0
  end function get_desired_pressure

  subroutine apply_pressure_input_primitive_var_bc(self, rho, u, v, p, lbounds)
    !< Apply pressure_input boundary conditions to the conserved state vector field
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    logical :: inflow
    logical :: outflow

    real(rk) :: desired_boundary_pressure, boundary_density, gamma
    real(rk) :: mach_u, mach_v, mach, cs
    real(rk), dimension(:), allocatable :: edge_pressure
    real(rk), dimension(:), allocatable :: domain_rho
    real(rk), dimension(:), allocatable :: domain_u
    real(rk), dimension(:), allocatable :: domain_v
    real(rk), dimension(:), allocatable :: domain_p
    real(rk), dimension(4) :: boundary_prim_vars, domain_prim_vars
    real(rk), dimension(2) :: boundary_norm

    gamma = eos%get_gamma()
    inflow = .false.
    outflow = .false.

    left_ghost = lbound(rho, dim=1)
    right_ghost = ubound(rho, dim=1)
    bottom_ghost = lbound(rho, dim=2)
    top_ghost = ubound(rho, dim=2)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    desired_boundary_pressure = self%get_desired_pressure()

    if(allocated(self%edge_rho)) deallocate(self%edge_rho)
    if(allocated(self%edge_u)) deallocate(self%edge_u)
    if(allocated(self%edge_v)) deallocate(self%edge_v)
    if(allocated(self%edge_p)) deallocate(self%edge_p)

    select case(self%location)
    case('+x', '-x')
      allocate(self%edge_rho(bottom_ghost:top_ghost))
      allocate(self%edge_u(bottom_ghost:top_ghost))
      allocate(self%edge_v(bottom_ghost:top_ghost))
      allocate(self%edge_p(bottom_ghost:top_ghost))
      allocate(domain_rho(bottom_ghost:top_ghost))
      allocate(domain_u(bottom_ghost:top_ghost))
      allocate(domain_v(bottom_ghost:top_ghost))
      allocate(domain_p(bottom_ghost:top_ghost))
    case('+y', '-y')
      allocate(self%edge_rho(left_ghost:right_ghost))
      allocate(self%edge_u(left_ghost:right_ghost))
      allocate(self%edge_v(left_ghost:right_ghost))
      allocate(self%edge_p(left_ghost:right_ghost))
      allocate(domain_rho(left_ghost:right_ghost))
      allocate(domain_u(left_ghost:right_ghost))
      allocate(domain_v(left_ghost:right_ghost))
      allocate(domain_p(left_ghost:right_ghost))
    end select

    self%edge_rho = 0.0_rk
    self%edge_u = 0.0_rk
    self%edge_v = 0.0_rk
    self%edge_p = 0.0_rk

    domain_rho = 0.0_rk
    domain_u = 0.0_rk
    domain_v = 0.0_rk
    domain_p = 0.0_rk

    select case(self%location)
    case('+x')
      call debug_print('Running pressure_input_bc_t%apply_pressure_input_primitive_var_bc() +x', __FILE__, __LINE__)

      domain_rho = sum(rho(right, bottom_ghost:top_ghost)) / size(rho(right, bottom_ghost:top_ghost))
      domain_u = sum(u(right, bottom_ghost:top_ghost)) / size(u(right, bottom_ghost:top_ghost))
      domain_v = 0.0_rk
      domain_p = sum(p(right, bottom_ghost:top_ghost)) / size(p(right, bottom_ghost:top_ghost))
      boundary_norm = [1.0_rk, 0.0_rk]

      do j = bottom, top
        associate(rho_d=>domain_rho(j), u_d=>domain_u(j), &
                  v_d=>domain_v(j), p_d=>domain_p(j))
          domain_prim_vars = [rho_d, u_d, v_d, p_d]

          call eos%sound_speed(p=p_d, rho=rho_d, cs=cs)
          mach_u = u_d / cs

          if(mach_u > 0.0_rk) then ! outlet
            if(mach_u > 1.0_rk) then

              boundary_prim_vars = supersonic_outlet(domain_prim_vars=domain_prim_vars)
            else
              boundary_prim_vars = subsonic_outlet(domain_prim_vars=domain_prim_vars, &
                                                   exit_pressure=desired_boundary_pressure, &
                                                   boundary_norm=boundary_norm)
            end if
          else ! inlet
            if(abs(mach_u) > 1.0_rk) then
              error stop 'supersonic inlet'
            else
              boundary_prim_vars = subsonic_inlet(domain_prim_vars=domain_prim_vars, &
                                                  boundary_norm=boundary_norm, &
                                                  inlet_density=self%density_input, &
                                                  inlet_total_press=desired_boundary_pressure, &
                                                  inlet_flow_angle=0.0_rk)
            end if
          end if
        end associate
        self%edge_rho(j) = boundary_prim_vars(1)
        self%edge_u(j) = boundary_prim_vars(2)
        self%edge_v(j) = boundary_prim_vars(3)
        self%edge_p(j) = boundary_prim_vars(4)
      end do

      rho(right_ghost, bottom:top) = self%edge_rho(bottom:top)
      u(right_ghost, bottom:top) = self%edge_u(bottom:top)
      v(right_ghost, bottom:top) = self%edge_v(bottom:top)
      p(right_ghost, bottom:top) = self%edge_p(bottom:top)
    case default
      error stop "Unsupported location to apply the bc at in "// &
        "pressure_input_bc_t%apply_pressure_input_primitive_var_bc()"
    end select

    if(allocated(edge_pressure)) deallocate(edge_pressure)

  end subroutine apply_pressure_input_primitive_var_bc

  subroutine apply_pressure_input_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply pressure boundary conditions to the reconstructed state vector field

    class(pressure_input_bc_t), intent(inout) :: self

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index
    integer(ik) :: n, p

    left_ghost = lbound(recon_rho, dim=2)
    right_ghost = ubound(recon_rho, dim=2)
    bottom_ghost = lbound(recon_rho, dim=3)
    top_ghost = ubound(recon_rho, dim=3)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    select case(self%location)
    case('+x')
      call debug_print('Running pressure_input_bc_t%apply_pressure_input_reconstructed_state_bc() +x', __FILE__, __LINE__)
      do p = 1, 8
        recon_rho(p, right_ghost, :) = self%edge_rho
        recon_u(p, right_ghost, :) = self%edge_u
        recon_v(p, right_ghost, :) = self%edge_v
        recon_p(p, right_ghost, :) = self%edge_p
      end do

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
    case default
      error stop "Unsupported location to apply the bc at in pressure_input_bc_t%apply_pressure_input_reconstructed_state_bc()"
    end select

    if(allocated(self%edge_rho)) deallocate(self%edge_rho)
    if(allocated(self%edge_u)) deallocate(self%edge_u)
    if(allocated(self%edge_v)) deallocate(self%edge_v)
    if(allocated(self%edge_p)) deallocate(self%edge_p)
  end subroutine apply_pressure_input_reconstructed_state_bc

end module mod_pressure_input_bc
