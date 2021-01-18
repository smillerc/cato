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

module mod_pressure_input_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_error, only: error_msg
  use mod_globals, only: enable_debug_print, debug_print
  use mod_field, only: field_2d_t
  use mod_grid_block_2d, only: grid_block_2d_t
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
    logical :: constant_density = .false.
    real(rk) :: pressure_input = 0.0_rk
    real(rk) :: temperature_input = 273.0_rk ! K
    real(rk) :: density_input = 5e-3_rk !< Boundary density g/cc
    real(rk) :: scale_factor = 1.0_rk !< scaling factor (e.g. 1x, 2x) to scale the input up/down

    ! Temporal inputs
    class(linear_interp_1d), allocatable :: temporal_pressure_input
    class(linear_interp_1d), allocatable :: temporal_density_input
    character(len=100) :: input_filename

    real(rk), dimension(:), allocatable :: edge_rho !< boundary (ghost/edge/etc) density
    real(rk), dimension(:), allocatable :: edge_u   !< boundary (ghost/edge/etc) x-velocity
    real(rk), dimension(:), allocatable :: edge_v   !< boundary (ghost/edge/etc) y-velocity
    real(rk), dimension(:), allocatable :: edge_p   !< boundary (ghost/edge/etc) pressure

  contains
    procedure, private :: read_pressure_input
    procedure, private :: get_desired_pressure
    procedure, private :: get_desired_density
    procedure, public :: apply => apply_pressure_input_primitive_var_bc
    final :: finalize
  endtype
contains

  function pressure_input_bc_constructor(location, input, grid, time) result(bc)
    type(pressure_input_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid
    real(rk), intent(in) :: time

    if(enable_debug_print) call debug_print('Running pressure_input_bc_t%pressure_input_bc_constructor() ', __FILE__, __LINE__)

    allocate(bc)
    bc%name = 'pressure_input'
    bc%location = location
    bc%priority = 0 !< currently the lowest (last to go) priority
    call bc%set_indices(grid)

    bc%time = time
    bc%density_input = input%bc_density / rho_0
    bc%constant_pressure = input%apply_constant_bc_pressure
    bc%scale_factor = input%bc_pressure_scale_factor
    if(.not. bc%constant_pressure) then
      bc%input_filename = trim(input%bc_pressure_input_file)
      call bc%read_pressure_input()
    else
      bc%pressure_input = input%constant_bc_pressure_value / p_0 ! non-dimensionalize
    endif

  endfunction pressure_input_bc_constructor

  subroutine read_pressure_input(self)
    !< Read the pressure input file and initialize the linear iterpolated object. This
    !< makes it easy to get the input pressure at any given time
    class(pressure_input_bc_t), intent(inout) :: self

    real(rk), dimension(:), allocatable :: time_sec       !< User input for time in seconds
    real(rk), dimension(:), allocatable :: pressure_barye !< User input for pressure in barye
    real(rk), dimension(:), allocatable :: density_gcc    !< User input for density in g/cc
    character(len=300) :: line_buffer
    logical :: has_header_line
    logical :: file_exists
    integer(ik) :: input_unit, io_status, interp_status
    integer(ik) :: line, nlines

    if(enable_debug_print) call debug_print('Running pressure_input_bc_t%read_pressure_input() ', __FILE__, __LINE__)

    file_exists = .false.
    has_header_line = .false.

    inquire(file=trim(self%input_filename), exist=file_exists)

    if(.not. file_exists) then
      call error_msg(module_name='mod_pressure_input_bc', class_name='pressure_input_bc_t', procedure_name='read_pressure_input', &
                     message="Pressure input file '"//trim(self%input_filename)//"' not found", &
                     file_name=__FILE__, line_number=__LINE__)
    endif

    ! Open and determine how many lines there are
    open(newunit=input_unit, file=trim(self%input_filename))
    nlines = 0
    do
      read(input_unit, *, iostat=io_status) line_buffer
      if(nlines == 0) then
        line_buffer = adjustl(line_buffer)
        if(line_buffer(1:1) == '#') has_header_line = .true.
      endif

      if(io_status /= 0) exit
      nlines = nlines + 1
    enddo
    close(input_unit)

    ! Re-open and now allocate the arrays to the proper size based on the # of lines in the file
    open(newunit=input_unit, file=trim(self%input_filename))
    if(has_header_line) then
      read(input_unit, *, iostat=io_status) ! skip the first line
      nlines = nlines - 1
    endif

    allocate(time_sec(nlines))
    allocate(pressure_barye(nlines))
    allocate(density_gcc(nlines))
    time_sec = 0.0_rk
    pressure_barye = 0.0_rk
    density_gcc = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time_sec(line), pressure_barye(line), density_gcc(line)
      if(io_status /= 0) exit
    enddo
    close(input_unit)

    ! Apply scaling and non-dimensionalization
    time_sec = time_sec / t_0
    pressure_barye = (pressure_barye / p_0) * self%scale_factor
    density_gcc = (density_gcc / rho_0) * self%scale_factor

    self%max_time = maxval(time_sec)

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    allocate(self%temporal_pressure_input)
    call self%temporal_pressure_input%initialize(time_sec, pressure_barye, interp_status)
    if(interp_status /= 0) then
      write(*, *) time_sec
      write(*, *) pressure_barye
      write(*, *) density_gcc
      write(*, *) interp_status
      error stop "Error initializing pressure_input_bc_t%temporal_pressure_input"
    endif

    allocate(self%temporal_density_input)
    call self%temporal_density_input%initialize(time_sec, density_gcc, interp_status)
    if(interp_status /= 0) then
      write(*, *) time_sec
      write(*, *) pressure_barye
      write(*, *) density_gcc
      write(*, *) interp_status
      error stop "Error initializing pressure_input_bc_t%temporal_density_input"
    endif

    deallocate(time_sec)
    deallocate(pressure_barye)
    deallocate(density_gcc)
  endsubroutine read_pressure_input

  real(rk) function get_desired_pressure(self) result(desired_pressure)
    !< Get the pressure specified in the input file or external file (for time-dependent cases).
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik) :: interp_stat

    if(enable_debug_print) call debug_print('Running pressure_input_bc_t%get_desired_pressure() ', __FILE__, __LINE__)

    if(self%constant_pressure) then
      desired_pressure = self%pressure_input
    else
      call self%temporal_pressure_input%evaluate(x=self%get_time(), f=desired_pressure, istat=interp_stat)
      if(interp_stat /= 0) then
        error stop "Unable to interpolate pressure within pressure_input_bc_t%get_desired_pressure()"
      endif
    endif

    if(desired_pressure < tiny(1.0_rk)) then
      error stop "Error in pressure_input_bc_t%get_desired_pressure: desired_pressure < tiny(1.0_rk)"
    endif

  endfunction get_desired_pressure

  real(rk) function get_desired_density(self) result(desired_density)
    !< Get the density specified in the input file or external file (for time-dependent cases).
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik) :: interp_stat

    if(enable_debug_print) call debug_print('Running pressure_input_bc_t%get_desired_density() ', __FILE__, __LINE__)

    if(self%constant_pressure) then
      desired_density = self%density_input
    else
      call self%temporal_density_input%evaluate(x=self%get_time(), f=desired_density, istat=interp_stat)

      if(interp_stat /= 0) then
        error stop "Unable to interpolate density within pressure_input_bc_t%get_desired_density()"
      endif
    endif

    if(desired_density < tiny(1.0_rk)) then
      error stop "Error in pressure_input_bc_t%get_desired_density: desired_density < tiny(1.0_rk)"
    endif
  endfunction get_desired_density

  subroutine apply_pressure_input_primitive_var_bc(self, rho, u, v, p)
    !< Apply pressure_input boundary conditions to the primitive variables
    class(pressure_input_bc_t), intent(inout) :: self
    class(field_2d_t), intent(inout) :: rho
    class(field_2d_t), intent(inout) :: u
    class(field_2d_t), intent(inout) :: v
    class(field_2d_t), intent(inout) :: p

    integer(ik) :: i, j
    logical :: inflow
    logical :: outflow

    real(rk) :: desired_boundary_density  !< user-specified boundary value for density
    real(rk) :: desired_boundary_pressure !< user-specified boundary value for pressure
    real(rk) :: gamma, boundary_cs
    real(rk) :: mach_u, cs
    real(rk), dimension(:), allocatable :: edge_pressure
    real(rk), dimension(:), allocatable :: domain_rho
    real(rk), dimension(:), allocatable :: domain_u
    real(rk), dimension(:), allocatable :: domain_v
    real(rk), dimension(:), allocatable :: domain_p
    real(rk), dimension(4) :: boundary_prim_vars, domain_prim_vars
    real(rk), dimension(2) :: boundary_norm

    if(enable_debug_print) then
      call debug_print('Running pressure_input_bc_t%apply_pressure_input_primitive_var_bc() ', __FILE__, __LINE__)
    endif

    if(rho%on_ihi_bc .or. rho%on_ilo_bc .or. &
       rho%on_jhi_bc .or. rho%on_jlo_bc) then

      gamma = eos%get_gamma()
      inflow = .false.
      outflow = .false.

      boundary_prim_vars = 0.0_rk
      desired_boundary_pressure = self%get_desired_pressure()
      desired_boundary_density = self%get_desired_density()

      if(allocated(self%edge_rho)) deallocate(self%edge_rho)
      if(allocated(self%edge_u)) deallocate(self%edge_u)
      if(allocated(self%edge_v)) deallocate(self%edge_v)
      if(allocated(self%edge_p)) deallocate(self%edge_p)
    endif

    associate(left => self%ilo, right => self%ihi, &
              ! bottom => self%jlo, &
              ! top => self%jhi, &
              bottom => minval(self%jlo_ghost), &
              top => maxval(self%jhi_ghost), &
              left_ghost => minval(self%ilo_ghost), &
              right_ghost => maxval(self%ihi_ghost), &
              bottom_ghost => minval(self%jlo_ghost), &
              top_ghost => maxval(self%jhi_ghost), &
              n => self%n_ghost_layers)

      ! print*, 'pressure bc edges'
      ! print*, self%ilo, self%ihi, self%jlo, self%jhi
      ! print*, self%ilo_ghost, self%ihi_ghost, self%jlo_ghost, self%jhi_ghost
      ! print*

      select case(self%location)
      case('+x', '-x')
        if(rho%on_ihi_bc .or. rho%on_ilo_bc) then
          allocate(self%edge_rho(bottom_ghost:top_ghost))
          allocate(self%edge_u(bottom_ghost:top_ghost))
          allocate(self%edge_v(bottom_ghost:top_ghost))
          allocate(self%edge_p(bottom_ghost:top_ghost))
          allocate(domain_rho(bottom_ghost:top_ghost))
          allocate(domain_u(bottom_ghost:top_ghost))
          allocate(domain_v(bottom_ghost:top_ghost))
          allocate(domain_p(bottom_ghost:top_ghost))

          self%edge_rho = 0.0_rk
          self%edge_u = 0.0_rk
          self%edge_v = 0.0_rk
          self%edge_p = 0.0_rk
          domain_rho = 0.0_rk
          domain_u = 0.0_rk
          domain_v = 0.0_rk
          domain_p = 0.0_rk
        endif
      case('+y', '-y')
        if(rho%on_jhi_bc .or. rho%on_jlo_bc) then
          allocate(self%edge_rho(left_ghost:right_ghost))
          allocate(self%edge_u(left_ghost:right_ghost))
          allocate(self%edge_v(left_ghost:right_ghost))
          allocate(self%edge_p(left_ghost:right_ghost))
          allocate(domain_rho(left_ghost:right_ghost))
          allocate(domain_u(left_ghost:right_ghost))
          allocate(domain_v(left_ghost:right_ghost))
          allocate(domain_p(left_ghost:right_ghost))
          self%edge_rho = 0.0_rk
          self%edge_u = 0.0_rk
          self%edge_v = 0.0_rk
          self%edge_p = 0.0_rk
          domain_rho = 0.0_rk
          domain_u = 0.0_rk
          domain_v = 0.0_rk
          domain_p = 0.0_rk
        endif
      endselect

      select case(self%location)
      case('+x')
        if(rho%on_ihi_bc) then
          if(enable_debug_print) then
            call debug_print('Running pressure_input_bc_t%apply_pressure_input_primitive_var_bc() +x', &
                             __FILE__, __LINE__)
          endif

          domain_rho = rho%data(right, bottom_ghost:top_ghost)
          domain_u = u%data(right, bottom_ghost:top_ghost)
          domain_v = 0.0_rk
          domain_p = p%data(right, bottom_ghost:top_ghost)
          boundary_norm = [1.0_rk, 0.0_rk] ! outward normal vector at +x

          do j = bottom, top
            associate(rho_d => domain_rho(j), u_d => domain_u(j), &
                      v_d => domain_v(j), p_d => domain_p(j))
              domain_prim_vars = [rho_d, u_d, v_d, p_d]

              call eos%sound_speed(p=p_d, rho=rho_d, cs=cs)
              call eos%sound_speed(p=desired_boundary_pressure, rho=desired_boundary_density, cs=boundary_cs)
              mach_u = u_d / cs

              if(mach_u > 0.0_rk) then ! outlet
                if(mach_u <= 1.0_rk .or. maxval(domain_p(bottom:top)) < desired_boundary_pressure) then
                  boundary_prim_vars = subsonic_outlet(domain_prim_vars=domain_prim_vars, &
                                                       exit_pressure=desired_boundary_pressure, &
                                                       boundary_norm=boundary_norm)
                else
                  boundary_prim_vars = supersonic_outlet(domain_prim_vars=domain_prim_vars)
                endif
              else ! inlet
                if(abs(mach_u) > 1.0_rk) then
                  boundary_prim_vars = [desired_boundary_density, u_d, v_d, desired_boundary_pressure]
                else
                  boundary_prim_vars = subsonic_inlet(domain_prim_vars=domain_prim_vars, &
                                                      boundary_norm=boundary_norm, &
                                                      inlet_density=desired_boundary_density, &
                                                      inlet_total_press=desired_boundary_pressure, &
                                                      inlet_flow_angle=0.0_rk)
                endif
              endif
            endassociate
            self%edge_rho(j) = boundary_prim_vars(1)
            self%edge_u(j) = boundary_prim_vars(2)
            self%edge_v(j) = 0.0_rk
            self%edge_p(j) = boundary_prim_vars(4)

            ! write(*, '(4(es10.3), a, 4(es10.3))') domain_prim_vars, ":", boundary_prim_vars
          enddo

          do i = 1, self%n_ghost_layers
            rho%data(self%ihi_ghost(i), bottom:top) = self%edge_rho(bottom:top)
            u%data(self%ihi_ghost(i), bottom:top) = self%edge_u(bottom:top)
            v%data(self%ihi_ghost(i), bottom:top) = self%edge_v(bottom:top)
            p%data(self%ihi_ghost(i), bottom:top) = self%edge_p(bottom:top)
          enddo

!           write(*, '(a, 50(es12.3))') 'rho', rho%data(:,1)
!           write(*, '(a, 50(es12.3))') 'u  ', u%data(:,1)
!           write(*, '(a, 50(es12.3))') 'v  ', v%data(:,1)
!           write(*, '(a, 50(es12.3))') 'p  ', p%data(:,1)
! print*
!           write(*, '(a, 50(es12.3))') 'rho', rho%data(22,:)
!           write(*, '(a, 50(es12.3))') 'u  ',   u%data(22,:)
!           write(*, '(a, 50(es12.3))') 'v  ',   v%data(22,:)
!           write(*, '(a, 50(es12.3))') 'p  ',   p%data(22,:)
        endif ! on_ihi_bc
      case default
        call error_msg(module_name='mod_pressure_bc', class_name='pressure_input_bc_t', &
                       procedure_name='apply_pressure_input_primitive_var_bc', &
                       message="Unsupported pressure_bc location '"//self%location//"'", &
                       file_name=__FILE__, line_number=__LINE__)
      endselect

    endassociate
    if(allocated(edge_pressure)) deallocate(edge_pressure)

  endsubroutine apply_pressure_input_primitive_var_bc

  subroutine finalize(self)
    type(pressure_input_bc_t), intent(inout) :: self
    if(enable_debug_print) call debug_print('Running pressure_input_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)

    if(allocated(self%edge_rho)) deallocate(self%edge_rho)
    if(allocated(self%edge_u)) deallocate(self%edge_u)
    if(allocated(self%edge_v)) deallocate(self%edge_v)
    if(allocated(self%edge_p)) deallocate(self%edge_p)

    if(allocated(self%temporal_density_input)) deallocate(self%temporal_density_input)
    if(allocated(self%temporal_pressure_input)) deallocate(self%temporal_pressure_input)
  endsubroutine finalize
endmodule mod_pressure_input_bc
