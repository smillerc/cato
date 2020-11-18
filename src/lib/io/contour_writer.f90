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

module mod_contour_writer

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, real32
  use mod_master_puppeteer, only: master_puppeteer_t
  use mod_globals, only: n_ghost_layers, debug_print
  use mod_fluid, only: fluid_t
  use mod_grid_block, only: grid_block_t
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_error, only: error_msg
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
    integer(ik), private :: compression_level = 6

    integer(ik) :: ilo_node, ihi_node
    integer(ik) :: jlo_node, jhi_node

    integer(ik) :: ilo_cell, ihi_cell
    integer(ik) :: jlo_cell, jhi_cell

  contains
    procedure, public :: write_contour
    procedure, private :: write_xdmf
    procedure, private :: write_hdf5
    procedure, private :: write_2d_integer_data
    procedure, private :: write_2d_real_data
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
    integer(ik) :: cstat
    integer(ik) :: estat
    character(100) :: cmsg

    cmsg = ''
    cstat = 0
    estat = 0

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

    if(this_image() == 1) then
      write(*, '(a)') "Contour Writer Info"
      write(*, '(a)') "==================="
      write(*, '(a, a)') 'Output folder: ', writer%results_folder
      write(*, '(a, l1)') 'plot_64bit: ', writer%plot_64bit
      write(*, '(a, l1)') 'plot_ghost_cells: ', writer%plot_ghost_cells
      write(*, *)
    endif

  end function

  subroutine finalize(self)
    type(contour_writer_t), intent(inout) :: self
    if(allocated(self%format)) deallocate(self%format)
    if(allocated(self%hdf5_filename)) deallocate(self%hdf5_filename)
    if(allocated(self%xdmf_filename)) deallocate(self%xdmf_filename)
  end subroutine

  subroutine write_contour(self, master, iteration)
    class(contour_writer_t), intent(inout) :: self
    class(master_puppeteer_t), intent(in) :: master
    integer(ik), intent(in) :: iteration
    real(rk) :: time
    character(50) :: char_buff

    time = master%time
    
    write(char_buff, '(a,i0.7)') 'step_', iteration
    self%hdf5_filename = trim(char_buff)//'.h5'
    self%xdmf_filename = trim(char_buff)//'.xdmf'

    if(self%plot_ghost_cells) then
      self%ilo_cell = master%grid%lbounds_halo(1)
      self%ihi_cell = master%grid%ubounds_halo(1)
      self%jlo_cell = master%grid%lbounds_halo(2)
      self%jhi_cell = master%grid%ubounds_halo(2)

      self%ilo_node = master%grid%lbounds_halo(1)
      self%ihi_node = master%grid%ubounds_halo(1) + 1
      self%jlo_node = master%grid%lbounds_halo(2)
      self%jhi_node = master%grid%ubounds_halo(2) + 1
    else
      self%ilo_cell = master%grid%lbounds(1)
      self%ihi_cell = master%grid%ubounds(1)
      self%jlo_cell = master%grid%lbounds(2)
      self%jhi_cell = master%grid%ubounds(2)

      self%ilo_node = master%grid%lbounds(1)
      self%ihi_node = master%grid%ubounds(1) + 1
      self%jlo_node = master%grid%lbounds(2)
      self%jhi_node = master%grid%ubounds(2) + 1
    end if

    if(this_image() == 1) write(*, '(a,a)') "Saving contour file: "//self%hdf5_filename
    select case(self%format)
    case('xdmf')
      call self%write_hdf5(master, time, iteration)
      call self%write_xdmf(master, time, iteration)
    case('hdf5', 'h5')
      call self%write_hdf5(master, time, iteration)
    case default
      call error_msg(module_name='mod_contour_writer', class_name='contour_writer_t', procedure_name='write_contour', &
                     message="Unsupported I/O contour format'"//self%format//"' (must be .xdmf .h5 or .hdf5)", &
                     file_name=__FILE__, line_number=__LINE__)
    end select

  end subroutine write_contour

  subroutine write_hdf5(self, master, time, iteration)
    class(contour_writer_t), intent(inout) :: self
    class(master_puppeteer_t), intent(in) :: master
    integer(ik), intent(in) :: iteration
    real(rk), intent(in) :: time
    character(32) :: dataset_name

    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk) :: time_w_dims, delta_t_w_dims
    real(rk), dimension(:, :), allocatable :: io_data_buffer
    integer(ik), dimension(:, :), allocatable :: int_data_buffer

    real(rk), allocatable :: int_gather_coarray(:, :)[:]
    

    time_w_dims = time * io_time_units * t_0
    delta_t_w_dims = master%dt * t_0

    if(this_image() == 1) then
      call self%hdf5_file%initialize(filename=self%results_folder//'/'//self%hdf5_filename, &
                                     status='new', action='w', comp_lvl=self%compression_level)

      ! Header info
      call self%hdf5_file%add('/title', master%title)
      call self%hdf5_file%add('/iteration', iteration)
      call self%hdf5_file%writeattr('/iteration', 'description', 'Iteration Count')
      call self%hdf5_file%writeattr('/iteration', 'units', 'dimensionless')

      call self%hdf5_file%add('/time', time_w_dims)
      call self%hdf5_file%writeattr('/time', 'description', 'Simulation Time')
      call self%hdf5_file%writeattr('/time', 'units', io_time_label)

      call self%hdf5_file%add('/delta_t', delta_t_w_dims)
      call self%hdf5_file%writeattr('/delta_t', 'description', 'Simulation Timestep')
      call self%hdf5_file%writeattr('/delta_t', 'units', 'seconds')

      call self%hdf5_file%add('/n_ghost_layers', master%grid%n_halo_cells)

      ! Version info
      if(.not. globals_set) call set_global_options()

      call self%hdf5_file%add('/cato_info', 'CATO info')
      call self%hdf5_file%writeattr('/cato_info', 'compiler_flags', compiler_flags_str)
      call self%hdf5_file%writeattr('/cato_info', 'compiler_version', compiler_version_str)
      call self%hdf5_file%writeattr('/cato_info', 'git_hash', git_hash)
      call self%hdf5_file%writeattr('/cato_info', 'git_ref', git_ref)
      call self%hdf5_file%writeattr('/cato_info', 'git_changes', git_local_changes)
      call self%hdf5_file%writeattr('/cato_info', 'version', cato_version)
      call self%hdf5_file%writeattr('/cato_info', 'compile_hostname', compile_host)
      call self%hdf5_file%writeattr('/cato_info', 'compile_os', compile_os)
      call self%hdf5_file%writeattr('/cato_info', 'build_type', build_type)
    end if

    ! Node Data
    dataset_name = '/x'
    
    
    io_data_buffer = master%grid%gather(var='x', image=1)
    if (this_image() == 1) then
      call self%write_2d_real_data(data=io_data_buffer * l_0 * io_length_units, name='/x', description='X Coordinate', units=trim(io_length_label))
    end if

    dataset_name = '/y'
    io_data_buffer = master%grid%gather(var='y', image=1)
    if (this_image() == 1) then
      call self%write_2d_real_data(data=io_data_buffer * l_0 * io_length_units, name='/y', description='Y Coordinate', units=trim(io_length_label))
    end if
    
    if (allocated(io_data_buffer)) deallocate(io_data_buffer)
    
    ! Cell Data
    ilo = 1
    ihi = master%grid%global_dims(1)
    jlo = 1
    jhi = master%grid%global_dims(2)
    if (this_image() == 1) allocate(int_data_buffer(ilo:ihi, jlo:jhi))

    ! ! if(self%plot_ghost_cells) then
    ! !   ! Write a simple flag to tag ghost cells
    ! !   dataset_name = '/ghost_cell'
    ! !   int_data_buffer = 1
    ! !   associate(ilo_r => master%grid%cell_lbounds(1), ihi_r => master%grid%cell_ubounds(1), &
    ! !             jlo_r => master%grid%cell_lbounds(2), jhi_r => master%grid%cell_ubounds(2))
    ! !     int_data_buffer(ilo_r:ihi_r, jlo_r:jhi_r) = 0
    ! !   end associate
    ! !   call self%write_2d_integer_data(data=int_data_buffer, name='/ghost_cell', &
    ! !                                   description='Ghost Cell [0=no, 1=yes]', units='dimensionless')
    ! ! end if

    ! Indexing
    dataset_name = '/image_id'
    allocate(int_gather_coarray(master%fluid%rho%global_dims(1), master%fluid%rho%global_dims(2))[*])
    associate(ilo => master%fluid%rho%lbounds(1), ihi => master%fluid%rho%ubounds(1), &
             jlo => master%fluid%rho%lbounds(2), jhi => master%fluid%rho%ubounds(2))

      int_gather_coarray(ilo:ihi, jlo:jhi)[1] = master%fluid%rho%host_image_id
      sync all
    end associate
    
    if(this_image() == 1) then
        int_data_buffer = int_gather_coarray
        call self%write_2d_integer_data(data=int_data_buffer, name='/image_id', &
                                        description='Coarray Image Index', units='dimensionless')
    endif
    deallocate(int_gather_coarray)

    if (this_image() == 1) then
      dataset_name = '/i'
      int_data_buffer = 0
      do i = ilo, ihi
        int_data_buffer(i, :) = i
      end do

      call self%write_2d_integer_data(data=int_data_buffer, name='/i', &
                                    description='Cell i Index', units='dimensionless')
    endif

    
    if (this_image() == 1) then
      dataset_name = '/j'
      int_data_buffer = 0
      do j = jlo, jhi
        int_data_buffer(:, j) = j
      end do
      call self%write_2d_integer_data(data=int_data_buffer, name='/j', &
                                    description='Cell j Index', units='dimensionless')
    endif

    if (allocated(int_data_buffer)) deallocate(int_data_buffer)

    ! call self%write_2d_integer_data(data=master%fluid%continuous_sensor(ilo:ihi, jlo:jhi), name='/continuity_sensor', &
    !                            description='Continuity Sensor [0=continuous, 1=linear discontinuity, 2=non-linear discontinuity]', &
    !                                 units='dimensionless')

    ! Primitive Variables
    dataset_name = '/density'
    io_data_buffer = master%fluid%rho%gather(image=1)
    if (this_image() == 1) call self%write_2d_real_data(data=io_data_buffer * rho_0 * io_density_units, name='/density', &
                                 description='Cell Density', units=trim(io_density_label))

    dataset_name = '/x_velocity'
    io_data_buffer = master%fluid%u%gather(image=1)
    if (this_image() == 1) call self%write_2d_real_data(data=io_data_buffer * v_0 * io_velocity_units, name='/x_velocity', &
                                 description='Cell X Velocity', units=trim(io_velocity_label))

    dataset_name = '/y_velocity'
    io_data_buffer = master%fluid%v%gather(image=1)
    if (this_image() == 1) call self%write_2d_real_data(data=io_data_buffer * v_0 * io_velocity_units, name='/y_velocity', &
                                 description='Cell Y Velocity', units=trim(io_velocity_label))

    dataset_name = '/pressure'
    io_data_buffer = master%fluid%p%gather(image=1)
    if (this_image() == 1) call self%write_2d_real_data(data=io_data_buffer * p_0 * io_pressure_units, name='/pressure', &
                                 description='Cell Pressure', units=trim(io_pressure_label))

    dataset_name = '/sound_speed'
    io_data_buffer = master%fluid%cs%gather(image=1)
    if (this_image() == 1) call self%write_2d_real_data(data=io_data_buffer * v_0 * io_velocity_units, name='/sound_speed', &
                                 description='Cell Sound Speed', units=trim(io_velocity_label))
    ! Volume
    dataset_name = '/volume'
    io_data_buffer = master%grid%gather(var='volume', image=1)
    
    if (this_image() == 1) call self%write_2d_real_data(data=io_data_buffer * (l_0*l_0) * io_volume_units, name='/volume', &
                                 description='Cell Volume', units=trim(io_volume_label))

    if(allocated(int_data_buffer)) deallocate(int_data_buffer)
    if(allocated(io_data_buffer)) deallocate(io_data_buffer)
    if(this_image() == 1) call self%hdf5_file%finalize()
  end subroutine write_hdf5

  subroutine write_xdmf(self, master, time, iteration)
    class(contour_writer_t), intent(inout) :: self
    class(master_puppeteer_t), intent(in) :: master
    integer(ik), intent(in) :: iteration
    real(rk), intent(in) :: time
    integer(ik) :: xdmf_unit
    character(50) :: char_buff
    character(10) :: unit_label = ''
    character(:), allocatable :: cell_shape, node_shape

    open(newunit=xdmf_unit, file=self%results_folder//'/'//self%xdmf_filename, status='replace')

    write(char_buff, '(2(i0,1x))') master%fluid%rho%global_dims(2), master%fluid%rho%global_dims(1)
    cell_shape = trim(char_buff)

    write(char_buff, '(2(i0,1x))') master%fluid%rho%global_dims(2) + 1, master%fluid%rho%global_dims(1) + 1
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

    ! unit_label = "["//trim(io_temperature_label)//"]"
    ! write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Temperature '//trim(unit_label)//'">'
    ! write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
    !   '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/temperature</DataItem>'
    ! write(xdmf_unit, '(a)') '      </Attribute>'

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
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape &
      //' 2" ItemType="Function" Function="JOIN(ABS($0/$2), ABS($1/$2))">'
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

    ! if(self%plot_ghost_cells) then
    !   write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Ghost Cell">'
    !   write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
    !     '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/ghost_cell</DataItem>'
    !   write(xdmf_unit, '(a)') '      </Attribute>'
    ! endif

    ! write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Continuity Sensor">'
    ! write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
    !   '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/continuity_sensor</DataItem>'
    ! write(xdmf_unit, '(a)') '      </Attribute>'
    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="image_id">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/image_id</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="i">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/i</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="j">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Int" Precision="2">'//self%hdf5_filename//':/j</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '    </Grid>'
    write(xdmf_unit, '(a)') '  </Domain>'
    write(xdmf_unit, '(a)') '</Xdmf>'

    close(xdmf_unit)

    if(allocated(cell_shape)) deallocate(cell_shape)
    if(allocated(node_shape)) deallocate(node_shape)
  end subroutine write_xdmf

  subroutine write_2d_integer_data(self, data, name, description, units)
    !< Helper subroutine to write data to the hdf5 file.
    class(contour_writer_t), intent(inout) :: self
    integer(ik), dimension(:, :), intent(in) :: data
    character(len=*), intent(in) :: name         !< dataset name, e.g. '/density'
    character(len=*), intent(in) :: description  !< dataset description
    character(len=*), intent(in) :: units        !< dataset units (if any)

    call self%hdf5_file%add(name, data)
    call self%hdf5_file%writeattr(name, 'description', trim(description))
    call self%hdf5_file%writeattr(name, 'units', trim(units))
  end subroutine write_2d_integer_data

  subroutine write_2d_real_data(self, data, name, description, units)
    !< Helper subroutine to write data to the hdf5 file.
    class(contour_writer_t), intent(inout) :: self
    real(rk), dimension(:, :), intent(in) :: data
    character(len=*), intent(in) :: name         !< dataset name, e.g. '/density'
    character(len=*), intent(in) :: description  !< dataset description
    character(len=*), intent(in) :: units        !< dataset units (if any)

    integer(ik):: ilo, ihi, jlo, jhi, i, j
    real(real32), dimension(:, :), allocatable :: single_prec_data

    ilo = lbound(data, dim=1)
    ihi = ubound(data, dim=1)
    jlo = lbound(data, dim=2)
    jhi = ubound(data, dim=2)

    if(self%plot_64bit) then
      call self%hdf5_file%add(name, data)
    else

      allocate(single_prec_data(ilo:ihi, jlo:jhi))

      ! Conversion checks
      do j = lbound(data, dim=2), ubound(data, dim=2)
        do i = lbound(data, dim=1), ubound(data, dim=1)
          if(abs(data(i, j)) > huge(1.0_real32)) then
            if(data(i, j) > 0.0_rk) single_prec_data(i, j) = huge(1.0_real32)
            if(data(i, j) < 0.0_rk) single_prec_data(i, j) = -huge(1.0_real32)
          else if(abs(data(i, j)) < tiny(1.0_real32)) then
            single_prec_data(i, j) = 0.0_rk
          else
            single_prec_data(i, j) = real(data(i, j), real32)
          end if
        end do
      end do

      call self%hdf5_file%add(name, single_prec_data)

      deallocate(single_prec_data)
    end if

    call self%hdf5_file%writeattr(name, 'description', trim(description))
    call self%hdf5_file%writeattr(name, 'units', trim(units))
  end subroutine write_2d_real_data
end module mod_contour_writer
