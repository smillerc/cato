module mod_contour_writer

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, real32
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_fluid, only: fluid_t
  use hdf5_interface, only: hdf5_file
  use mod_input, only: input_t
  use mod_functional, only: operator(.reverse.)
  use mod_globals, only: compiler_flags_str, compiler_version_str, git_hash, git_ref, &
                         git_local_changes, fvleg_2d_version, &
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

    writer%format = input%contour_io_format

    ! Turn on debug plotting?
    writer%plot_reconstruction_states = input%plot_reconstruction_states
    writer%plot_reference_states = input%plot_reference_states
    writer%plot_evolved_states = input%plot_evolved_states
    writer%plot_ghost_cells = input%plot_ghost_cells
    writer%plot_64bit = input%plot_64bit

    write(*, '(a)') "Contour Writer Inputs"
    write(*, '(a)') "====================="
    write(*, '(a, l1)') 'plot_reconstruction_states: ', writer%plot_reconstruction_states
    write(*, '(a, l1)') 'plot_reference_states: ', writer%plot_reference_states
    write(*, '(a, l1)') 'plot_evolved_states: ', writer%plot_evolved_states
    write(*, '(a, l1)') 'plot_ghost_cells: ', writer%plot_ghost_cells
    write(*, '(a, l1)') 'plot_64bit: ', writer%plot_64bit
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

    if(self%plot_ghost_cells) then
      self%ilo_cell = fv_scheme%grid%ilo_bc_cell
      self%ihi_cell = fv_scheme%grid%ihi_bc_cell
      self%jlo_cell = fv_scheme%grid%jlo_bc_cell
      self%jhi_cell = fv_scheme%grid%jhi_bc_cell

      self%ilo_node = fv_scheme%grid%ilo_bc_node
      self%ihi_node = fv_scheme%grid%ihi_bc_node
      self%jlo_node = fv_scheme%grid%jlo_bc_node
      self%jhi_node = fv_scheme%grid%jhi_bc_node

    else
      self%ilo_cell = fv_scheme%grid%ilo_cell
      self%ihi_cell = fv_scheme%grid%ihi_cell
      self%jlo_cell = fv_scheme%grid%jlo_cell
      self%jhi_cell = fv_scheme%grid%jhi_cell

      self%ilo_node = fv_scheme%grid%ilo_node
      self%ihi_node = fv_scheme%grid%ihi_node
      self%jlo_node = fv_scheme%grid%jlo_node
      self%jhi_node = fv_scheme%grid%jhi_node
    end if

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
    character(50) :: char_buff

    ! sync all
    ! if(this_image() == 1) then

    call self%hdf5_file%initialize(filename=self%hdf5_filename, &
                                   status='new', action='w', comp_lvl=6)

    ! Header info
    call self%hdf5_file%add('/title', fv_scheme%title)

    call self%hdf5_file%add('/iteration', iteration)
    call self%hdf5_file%writeattr('/iteration', 'description', 'Iteration Count')
    call self%hdf5_file%writeattr('/iteration', 'units', 'dimensionless')

    call self%hdf5_file%add('/time', time)
    call self%hdf5_file%writeattr('/time', 'description', 'Simulation Time')
    call self%hdf5_file%writeattr('/time', 'units', 'seconds')

    call self%hdf5_file%add('/delta_t', fv_scheme%delta_t)
    call self%hdf5_file%writeattr('/delta_t', 'description', 'Simulation Timestep')
    call self%hdf5_file%writeattr('/delta_t', 'units', 'seconds')

    ! Version info
    if(.not. globals_set) call set_global_options()
    call self%hdf5_file%writeattr('/', 'compiler_flags', compiler_flags_str)
    call self%hdf5_file%writeattr('/', 'compiler_version', compiler_version_str)
    call self%hdf5_file%writeattr('/', 'git_hast', git_hash)
    call self%hdf5_file%writeattr('/', 'git_ref', git_ref)
    call self%hdf5_file%writeattr('/', 'git_changes', git_local_changes)
    call self%hdf5_file%writeattr('/', 'version', fvleg_2d_version)
    call self%hdf5_file%writeattr('/', 'compile_hostname', compile_host)
    call self%hdf5_file%writeattr('/', 'compile_os', compile_os)
    call self%hdf5_file%writeattr('/', 'build_type', build_type)

    ! Grid
    if(self%plot_64bit) then
      call self%hdf5_file%add('/x', fv_scheme%grid%node_x(self%ilo_node:self%ihi_node, self%jlo_node:self%jhi_node))
    else
      call self%hdf5_file%add('/x', real(fv_scheme%grid%node_x(self%ilo_node:self%ihi_node, self%jlo_node:self%jhi_node), real32))
    end if
    call self%hdf5_file%writeattr('/x', 'description', 'X Coordinate')
    call self%hdf5_file%writeattr('/x', 'units', 'cm')

    if(self%plot_64bit) then
      call self%hdf5_file%add('/y', fv_scheme%grid%node_y(self%ilo_node:self%ihi_node, self%jlo_node:self%jhi_node))
    else
      call self%hdf5_file%add('/y', real(fv_scheme%grid%node_y(self%ilo_node:self%ihi_node, self%jlo_node:self%jhi_node), real32))
    end if
    call self%hdf5_file%writeattr('/y', 'description', 'Y Coordinate')
    call self%hdf5_file%writeattr('/y', 'units', 'cm')

    ! Conserved Variables
    if(self%plot_64bit) then
      call self%hdf5_file%add('/density', fluid%conserved_vars(1, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    else
call self%hdf5_file%add('/density', real(fluid%conserved_vars(1, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell), real32))
    end if
    call self%hdf5_file%writeattr('/density', 'description', 'Cell Density')
    call self%hdf5_file%writeattr('/density', 'units', 'g/cc')

    if(self%plot_64bit) then
      call self%hdf5_file%add('/x_velocity', fluid%conserved_vars(2, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    else
      call self%hdf5_file%add('/x_velocity', real(fluid%conserved_vars(2, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell), real32))
    end if
    call self%hdf5_file%writeattr('/x_velocity', 'description', 'Cell X Velocity')
    call self%hdf5_file%writeattr('/x_velocity', 'units', 'cm/s')

    if(self%plot_64bit) then
      call self%hdf5_file%add('/y_velocity', fluid%conserved_vars(3, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    else
      call self%hdf5_file%add('/y_velocity', real(fluid%conserved_vars(3, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell), real32))
    end if
    call self%hdf5_file%writeattr('/y_velocity', 'description', 'Cell Y Velocity')
    call self%hdf5_file%writeattr('/y_velocity', 'units', 'cm/s')

    if(self%plot_64bit) then
      call self%hdf5_file%add('/pressure', fluid%conserved_vars(4, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    else
      call self%hdf5_file%add('/pressure', real(fluid%conserved_vars(4, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell), real32))
    end if
    call self%hdf5_file%writeattr('/pressure', 'description', 'Cell Pressure')
    call self%hdf5_file%writeattr('/pressure', 'units', 'barye')

    ! Volume
    if(self%plot_64bit) then
      call self%hdf5_file%add('/volume', fv_scheme%grid%cell_volume(self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    else
      call self%hdf5_file%add('/volume', real(fv_scheme%grid%cell_volume(self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell), real32))
    end if
    call self%hdf5_file%writeattr('/volume', 'description', 'Cell Volume')
    call self%hdf5_file%writeattr('/volume', 'units', 'cc')

    ! Source Terms (if any)

    if(self%plot_evolved_states) then
      call self%hdf5_file%add('/evolved_corner_state', fv_scheme%evolved_corner_state)
      call self%hdf5_file%writeattr('/evolved_corner_state', 'indices', '((rho, u, v, p), i, j)')

      call self%hdf5_file%add('/evolved_downup_midpoints_state', fv_scheme%evolved_downup_midpoints_state)
      call self%hdf5_file%writeattr('/evolved_downup_midpoints_state', 'indices', '((rho, u, v, p), i, j)')

      call self%hdf5_file%add('/evolved_leftright_midpoints_state', fv_scheme%evolved_leftright_midpoints_state)
      call self%hdf5_file%writeattr('/evolved_leftright_midpoints_state', 'indices', '((rho, u, v, p), i, j)')
    end if

    if(self%plot_reconstruction_states) then
      call self%hdf5_file%add('/reconstructed_state', fv_scheme%reconstructed_state(:,:,:,self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
      call self%hdf5_file%writeattr('/reconstructed_state', 'indices', '((rho, u ,v, p), point, node/midpoint, i, j)')
    end if

    if(self%plot_reference_states) then
      call self%hdf5_file%add('/corner_reference_state', fv_scheme%corner_reference_state)
      call self%hdf5_file%writeattr('/corner_reference_state', 'indices', '((rho, u, v, p), i, j)')

      call self%hdf5_file%add('/downup_midpoints_reference_state', fv_scheme%downup_midpoints_reference_state)
      call self%hdf5_file%writeattr('/downup_midpoints_reference_state', 'indices', '((rho, u, v, p), i, j)')

      call self%hdf5_file%add('/leftright_midpoints_reference_state', fv_scheme%leftright_midpoints_reference_state)
      call self%hdf5_file%writeattr('/leftright_midpoints_reference_state', 'indices', '((rho, u, v, p), i, j)')
    end if

    ! Inputs
    call self%hdf5_file%finalize()
    ! end if
  end subroutine

  subroutine write_xdmf(self, fluid, fv_scheme, time, iteration)
    class(contour_writer_t), intent(inout) :: self
    class(finite_volume_scheme_t), intent(in) :: fv_scheme
    class(fluid_t), intent(in) :: fluid
    integer(ik), intent(in) :: iteration
    real(rk), intent(in) :: time
    integer(ik) :: xdmf_unit
    character(50) :: char_buff
    character(:), allocatable :: cell_shape, node_shape

    open(newunit=xdmf_unit, file=self%xdmf_filename, status='replace')

    write(char_buff, '(2(i0,1x))') .reverse.shape(fluid%conserved_vars(1, self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    cell_shape = trim(char_buff)

    write(char_buff, '(2(i0,1x))') .reverse.shape(fv_scheme%grid%node_x(self%ilo_cell:self%ihi_cell, self%jlo_cell:self%jhi_cell))
    node_shape = trim(char_buff)

    write(xdmf_unit, '(a)') '<?xml version="1.0" ?>'
    write(xdmf_unit, '(a)') '<Xdmf version="2.2">'
    write(xdmf_unit, '(a)') '  <Domain>'
    write(xdmf_unit, '(a)') '    <Grid GridType="Uniform" Name="grid">'
    write(xdmf_unit, '(a,g0.3,a)') '      <Time Value="', time, ' second"/>'
    write(xdmf_unit, '(a)') '      <Topology NumberOfElements="'//node_shape//'" TopologyType="2DSMesh"/>'

    write(xdmf_unit, '(a)') '      <Geometry GeometryType="X_Y">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//node_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/x</DataItem>'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//node_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/y</DataItem>'
    write(xdmf_unit, '(a)') '      </Geometry>'

    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Density [g/cc]">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/density</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="X Velocity [cm/s]">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/x_velocity</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Y Velocity [cm/s]">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/y_velocity</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Pressure [barye]">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/pressure</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '      <Attribute AttributeType="Scalar" Center="Cell" Name="Volume [cc]">'
    write(xdmf_unit, '(a)') '        <DataItem Dimensions="'//cell_shape// &
      '" Format="HDF" NumberType="Float" Precision="4">'//self%hdf5_filename//':/volume</DataItem>'
    write(xdmf_unit, '(a)') '      </Attribute>'

    write(xdmf_unit, '(a)') '    </Grid>'
    write(xdmf_unit, '(a)') '  </Domain>'
    write(xdmf_unit, '(a)') '</Xdmf>'
    close(xdmf_unit)

    deallocate(cell_shape)
    deallocate(node_shape)
  end subroutine
end module mod_contour_writer
