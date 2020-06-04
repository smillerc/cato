module mod_contour_writer

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, real32
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_fluid, only: fluid_t
  use mod_units
  use mod_nondimensionalization, only: rho_0, v_0, p_0, t_0, l_0, e_0
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_input, only: input_t
  use mod_functional, only: operator(.reverse.)
  use mod_globals, only: compiler_flags_str, compiler_version_str, git_hash, git_ref, &
                         git_local_changes, cato_version, &
                         compile_host, compile_os, build_type, set_global_options, globals_set

  implicit none

  private
  public :: contour_writer_t

  type :: contour_writer_t
    !< Type that manages writing out data to hdf5
    private
    type(hdf5_file) :: hdf5_file
    character(len=:), allocatable :: format !< xdmf or just plain hdf5
    character(len=:), allocatable :: hdf5_filename
    character(len=:), allocatable :: xdmf_filename
    character(len=:), allocatable :: results_folder
    logical, private :: plot_64bit = .false.
    logical, private :: plot_ghost_cells = .false.
    logical, private :: plot_reconstruction_states = .false.
    logical, private :: plot_reference_states = .false.
    logical, private :: plot_evolved_states = .false.

    integer(ik) :: ilo_node, ihi_node
    integer(ik) :: jlo_node, jhi_node

    integer(ik) :: ilo_cell, ihi_cell
    integer(ik) :: jlo_cell, jhi_cell

  contains
    procedure, public :: write_contour
    procedure, private :: write_xdmf
    procedure, private :: write_hdf5
    final :: finalize
  end type

  interface contour_writer_t
    module procedure :: constructor
  end interface
contains

  type(contour_writer_t) function constructor(input) result(writer)
    class(input_t), intent(in) :: input

    integer(ik) :: dt(8)
    character(len=32) :: str_buff = ''
    integer(ik) :: cstat, estat
    character(100) :: cmsg

    call date_and_time(values=dt)

    if(input%append_date_to_result_folder) then
      ! make the results folder look like: 'results_2019_08_12-21_50' -> results_YYYY_MM_DD-HH_mm
      write(str_buff, '(a, i4, 5(a, i2.2))') &
        'results_', dt(1), '_', dt(2), '_', dt(3), '-', dt(5), '_', dt(6)
    else
      str_buff = 'results'
    end if

    writer%results_folder = trim(str_buff)
    call execute_command_line('mkdir -p '//writer%results_folder, exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
    if(cstat > 0) then
      write(*, '(a)') "Unable to make results folder: "//trim(cmsg)
      error stop
    else if(cstat < 0) then
      write(*, '(a)') "Unable to make results folder"
      error stop "Unable to make results folder"
    end if

    writer%format = input%contour_io_format

    ! Turn on debug plotting?
    writer%plot_reconstruction_states = input%plot_reconstruction_states
    writer%plot_reference_states = input%plot_reference_states
    writer%plot_evolved_states = input%plot_evolved_states
    writer%plot_ghost_cells = input%plot_ghost_cells
    writer%plot_64bit = input%plot_64bit

    write(*, '(a)') "Contour Writer Info"
    write(*, '(a)') "==================="
    write(*, '(a, a)') 'Output folder: ', writer%results_folder
    ! write(*, '(a, l1)') 'plot_reconstruction_states: ', writer%plot_reconstruction_states
    ! write(*, '(a, l1)') 'plot_reference_states: ', writer%plot_reference_states
    ! write(*, '(a, l1)') 'plot_evolved_states: ', writer%plot_evolved_states
    ! write(*, '(a, l1)') 'plot_ghost_cells: ', writer%plot_ghost_cells
    ! write(*, '(a, l1)') 'plot_64bit: ', writer%plot_64bit
    write(*, *)

  end function

  subroutine finalize(self)
    type(contour_writer_t), intent(inout) :: self
    if(allocated(self%format)) deallocate(self%format)
    if(allocated(self%hdf5_filename)) deallocate(self%hdf5_filename)
    if(allocated(self%xdmf_filename)) deallocate(self%xdmf_filename)
  end subroutine

  subroutine write_contour(self, fluid, fv_scheme, time, iteration)
    class(contour_writer_t), intent(inout) :: self
    class(fluid_t), intent(in) :: fluid
    class(finite_volume_scheme_t), intent(in) :: fv_scheme
    integer(ik), intent(in) :: iteration
    real(rk), intent(in) :: time
    character(50) :: char_buff

    write(char_buff, '(a,i0.7)') 'step_', iteration
    self%hdf5_filename = trim(char_buff)//'.h5'
    self%xdmf_filename = trim(char_buff)//'.xdmf'

    self%ilo_cell = fv_scheme%grid%ilo_bc_cell
    self%ihi_cell = fv_scheme%grid%ihi_bc_cell
    self%jlo_cell = fv_scheme%grid%jlo_bc_cell
    self%jhi_cell = fv_scheme%grid%jhi_bc_cell

    self%ilo_node = fv_scheme%grid%ilo_bc_node
    self%ihi_node = fv_scheme%grid%ihi_bc_node
    self%jlo_node = fv_scheme%grid%jlo_bc_node
    self%jhi_node = fv_scheme%grid%jhi_bc_node

    write(*, '(a,a)') "Saving contour file: "//self%hdf5_filename
    select case(self%format)
    case('xdmf')
      call self%write_hdf5(fluid, fv_scheme, time, iteration)
      call self%write_xdmf(fluid, fv_scheme, time, iteration)
    case('hdf5', 'h5')
      call self%write_hdf5(fluid, fv_scheme, time, iteration)
    case default
      print *, 'Contour format:', self%format
      error stop "Unsupported I/O contour format"
    end select

  end subroutine

  subroutine write_hdf5(self, fluid, fv_scheme, time, iteration)
    class(contour_writer_t), intent(inout) :: self
    class(finite_volume_scheme_t), intent(in) :: fv_scheme
    class(fluid_t), intent(in) :: fluid
    integer(ik), intent(in) :: iteration
    real(rk), intent(in) :: time
    character(32) :: dataset_name

    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: time_w_dims, delta_t_w_dims
    real(rk), dimension(:, :), allocatable :: io_data_buffer
    integer(ik), dimension(:, :), allocatable :: int_data_buffer

    time_w_dims = time * io_time_units * t_0
    delta_t_w_dims = fv_scheme%delta_t * t_0

    call self%hdf5_file%initialize(filename=self%results_folder//'/'//self%hdf5_filename, &
                                   status='new', action='w', comp_lvl=6)

    ! Header info
    call self%hdf5_file%add('/title', fv_scheme%title)
    call self%hdf5_file%add('/iteration', iteration)
    call self%hdf5_file%writeattr('/iteration', 'description', 'Iteration Count')
    call self%hdf5_file%writeattr('/iteration', 'units', 'dimensionless')

    call self%hdf5_file%add('/time', time_w_dims)
    call self%hdf5_file%writeattr('/time', 'description', 'Simulation Time')
    call self%hdf5_file%writeattr('/time', 'units', io_time_label)

    call self%hdf5_file%add('/delta_t', delta_t_w_dims)
    call self%hdf5_file%writeattr('/delta_t', 'description', 'Simulation Timestep')
    call self%hdf5_file%writeattr('/delta_t', 'units', 'seconds')

    ! Version info
    if(.not. globals_set) call set_global_options()

    call self%hdf5_file%add('/cato_info', 'CATO info')
    call self%hdf5_file%writeattr('/cato_info', 'compiler_flags', compiler_flags_str)
    call self%hdf5_file%writeattr('/cato_info', 'compiler_version', compiler_version_str)
    call self%hdf5_file%writeattr('/cato_info', 'git_hast', git_hash)
    call self%hdf5_file%writeattr('/cato_info', 'git_ref', git_ref)
    call self%hdf5_file%writeattr('/cato_info', 'git_changes', git_local_changes)
    call self%hdf5_file%writeattr('/cato_info', 'version', cato_version)
    call self%hdf5_file%writeattr('/cato_info', 'compile_hostname', compile_host)
    call self%hdf5_file%writeattr('/cato_info', 'compile_os', compile_os)
    call self%hdf5_file%writeattr('/cato_info', 'build_type', build_type)

    ! Node Data
    ilo = self%ilo_node
    ihi = self%ihi_node
    jlo = self%jlo_node
    jhi = self%jhi_node
    allocate(io_data_buffer(ilo:ihi, jlo:jhi))

    dataset_name = '/x'
    io_data_buffer = fv_scheme%grid%node_x(ilo:ihi, jlo:jhi) * l_0
    io_data_buffer = io_data_buffer * io_length_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'X Coordinate')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_length_label))

    dataset_name = '/y'
    io_data_buffer = fv_scheme%grid%node_y(ilo:ihi, jlo:jhi) * l_0
    io_data_buffer = io_data_buffer * io_length_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Y Coordinate')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_length_label))

    deallocate(io_data_buffer)

    ! Cell Data
    ilo = self%ilo_cell
    ihi = self%ihi_cell
    jlo = self%jlo_cell
    jhi = self%jhi_cell
    allocate(io_data_buffer(ilo:ihi, jlo:jhi))
    allocate(int_data_buffer(ilo:ihi, jlo:jhi))

    ! Write a simple flag to tag ghost cells
    dataset_name = '/ghost_cell'
    int_data_buffer = 0
    associate(ilo_g=>fv_scheme%grid%ilo_bc_cell, &
              ihi_g=>fv_scheme%grid%ihi_bc_cell, &
              jlo_g=>fv_scheme%grid%jlo_bc_cell, &
              jhi_g=>fv_scheme%grid%jhi_bc_cell)

      int_data_buffer(ilo_g, :) = 1
      int_data_buffer(ihi_g, :) = 1
      int_data_buffer(:, jlo_g) = 1
      int_data_buffer(:, jhi_g) = 1
    end associate
    call self%hdf5_file%add(trim(dataset_name), int_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Ghost Cell [0=no, 1=yes]')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', 'dimensionless')

    ! Indexing
    dataset_name = '/i'
    io_data_buffer = 0
    do i = ilo, ihi
      io_data_buffer(i, :) = i
    end do

    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell i Index')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', 'dimensionless')

    dataset_name = '/j'
    io_data_buffer = 0
    do j = jlo, jhi
      io_data_buffer(:, j) = j
    end do

    io_data_buffer = io_data_buffer * io_density_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell j Index')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', 'dimensionless')

    deallocate(int_data_buffer)

    ! Primitive Variables
    dataset_name = '/density'
    io_data_buffer = fluid%rho * rho_0
    io_data_buffer = io_data_buffer * io_density_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell Density')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_density_label))

    dataset_name = '/x_velocity'
    io_data_buffer = fluid%u * v_0
    io_data_buffer = io_data_buffer * io_velocity_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell X Velocity')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_velocity_label))

    dataset_name = '/y_velocity'
    io_data_buffer = fluid%v * v_0
    io_data_buffer = io_data_buffer * io_velocity_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell Y Velocity')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_velocity_label))

    dataset_name = '/pressure'
    io_data_buffer = fluid%p * p_0
    io_data_buffer = io_data_buffer * io_pressure_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell Pressure')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_pressure_label))

    ! dataset_name = '/total_energy'
    ! associate(rhoE=>fluid%rho_E, &
    !           rho=>fluid%rho)
    !   io_data_buffer = (rhoE / rho) * e_0
    ! end associate
    ! io_data_buffer = io_data_buffer * io_pressure_units
    ! call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    ! call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell Total Energy')
    ! call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_energy_label))

    dataset_name = '/sound_speed'
    io_data_buffer = fluid%cs * v_0 * io_velocity_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell Sound Speed')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_velocity_label))

    dataset_name = '/temperature'
    associate(p=>fluid%p, &
              rho=>fluid%rho)
      call eos%temperature(p=p, rho=rho, t=io_data_buffer)
    end associate
    io_data_buffer = io_data_buffer * io_temperature_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell Temperature')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_temperature_label))

    ! Volume
    dataset_name = '/volume'
    io_data_buffer = fv_scheme%grid%cell_volume(:, :) * l_0**2
    io_data_buffer = io_data_buffer * io_volume_units
    call self%hdf5_file%add(trim(dataset_name), io_data_buffer)
    call self%hdf5_file%writeattr(trim(dataset_name), 'description', 'Cell Volume')
    call self%hdf5_file%writeattr(trim(dataset_name), 'units', trim(io_velocity_label))

    if(allocated(int_data_buffer)) deallocate(int_data_buffer)
    if(allocated(io_data_buffer)) deallocate(io_data_buffer)
    call self%hdf5_file%finalize()
  end subroutine

  subroutine write_xdmf(self, fluid, fv_scheme, time, iteration)
    class(contour_writer_t), intent(inout) :: self
    class(finite_volume_scheme_t), intent(in) :: fv_scheme
    class(fluid_t), intent(in) :: fluid
    integer(ik), intent(in) :: iteration
    real(rk), intent(in) :: time
    integer(ik) :: xdmf_unit
    character(50) :: char_buff
    character(10) :: unit_label = ''
    character(:), allocatable :: cell_shape, node_shape

    open(newunit=xdmf_unit, file=self%results_folder//'/'//self%xdmf_filename, status='replace')

    write(char_buff, '(2(i0,1x))') .reverse. &
      shape(fluid%rho(self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    cell_shape = trim(char_buff)

    write(char_buff, '(2(i0,1x))') .reverse. &
      shape(fv_scheme%grid%node_x(self%ilo_node:self%ihi_node, self%jlo_node:self%jhi_node))
    node_shape = trim(char_buff)

    write(xdmf_unit, '(a)') '<?xml version="1.0" ?>'
    write(xdmf_unit, '(a)') '<Xdmf version="2.2">'
    write(xdmf_unit, '(a)') '  <Domain>'
    write(xdmf_unit, '(a)') '    <Grid GridType="Uniform" Name="grid">'
    write(xdmf_unit, '(a,g0.3,a)') '      <Time Value="', time * io_time_units * t_0, ' '//trim(io_time_label)//'"/>'
    write(xdmf_unit, '(a)') '      <Topology NumberOfElements="'//node_shape//'" TopologyType="2DSMesh"/>'

    write(xdmf_unit, '(a)') '      <Geometry GeometryType="X_Y">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//node_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/x</DataItem>'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//node_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/y</DataItem>'
    write(xdmf_unit, '(a)') '      </Geometry>'

    unit_label = "["//trim(io_density_label)//"]"
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Density '//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/density</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    unit_label = "["//trim(io_velocity_label)//"]"
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Sound Speed '//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/sound_speed</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    unit_label = "["//trim(io_temperature_label)//"]"
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Temperature '//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/temperature</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    unit_label = "["//trim(io_pressure_label)//"]"
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Pressure '//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/pressure</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    ! unit_label = "["//trim(io_energy_label)//"]"
    ! write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Total Energy '//trim(unit_label)//'">'
    ! write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
    !   '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/total_energy</DataItem>'
    ! write(xdmf_unit, '(a)') '      </Attribute>'

    ! Velocity
    unit_label = "["//trim(io_velocity_label)//"]"
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Vector" Center="Cell" Name="Velocity '// &
      trim(unit_label)//'" Dimensions="'//cell_shape//' 2">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape//' 2" ItemType="Function" Function="JOIN($0, $1)">'
    write(xdmf_unit, '(a)') '          <DataItem DataType="Float" Dimensions="'//cell_shape// &
      '" Format="HDF" Precision="4">'//self%hdf5_filename//':/x_velocity</DataItem>'
    write(xdmf_unit, '(a)') '          <DataItem DataType="Float" Dimensions="'//cell_shape// &
      '" Format="HDF" Precision="4">'//self%hdf5_filename//':/y_velocity</DataItem>'
    write(xdmf_unit, '(a)') '        </DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    ! Mach Number
    unit_label = ""
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Vector" Center="Cell" Name="Mach" Dimensions="'//cell_shape//' 2">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape//' 2" ItemType="Function" Function="JOIN(ABS($0/$2), ABS($1/$2))">'
    write(xdmf_unit, '(a)') '          <DataItem DataType="Float" Dimensions="'//cell_shape// &
      '" Format="HDF" Precision="4">'//self%hdf5_filename//':/x_velocity</DataItem>'
    write(xdmf_unit, '(a)') '          <DataItem DataType="Float" Dimensions="'//cell_shape// &
      '" Format="HDF" Precision="4">'//self%hdf5_filename//':/y_velocity</DataItem>'
    write(xdmf_unit, '(a)') '          <DataItem DataType="Float" Dimensions="'//cell_shape// &
      '" Format="HDF" Precision="4">'//self%hdf5_filename//':/sound_speed</DataItem>'
    write(xdmf_unit, '(a)') '        </DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    unit_label = "["//trim(io_volume_label)//"]"
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Volume '//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/volume</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    unit_label = ""
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Ghost Cell '//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/ghost_cell</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    unit_label = ""
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="i'//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/i</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    unit_label = ""
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="j'//trim(unit_label)//'">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/j</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '    </Grid>'
    write(xdmf_unit, '(a)') '  </Domain>'
    write(xdmf_unit, '(a)') '</Xdmf>'

    close(xdmf_unit)

    deallocate(cell_shape)
    deallocate(node_shape)
  end subroutine
end module mod_contour_writer
