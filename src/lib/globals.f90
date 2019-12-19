module mod_globals
  use iso_fortran_env, only: compiler_options, compiler_version, std_out => output_unit

  implicit none

  logical, parameter :: enable_debug_print = .true.
  logical, parameter :: enable_file_and_line_stats = .false.

  logical :: globals_set = .false.
  character(:), allocatable :: compiler_flags_str
  character(:), allocatable :: compiler_version_str
  character(:), allocatable :: git_hash
  character(:), allocatable :: git_ref
  character(:), allocatable :: git_local_changes
  character(:), allocatable :: fvleg_2d_version
  character(:), allocatable :: compile_host
  character(:), allocatable :: compile_os
  character(:), allocatable :: build_type

contains

  subroutine set_global_options()
    if(.not. globals_set) then
      globals_set = .true.
      include 'version.h'
      compiler_flags_str = compiler_options()
      compiler_version_str = compiler_version()
    end if
  end subroutine set_global_options

  subroutine print_version_stats()
    call set_global_options()
    ! if(this_image() == 1) then
    write(std_out, '(a,a)')
    write(std_out, '(a,a)') "Version: ", fvleg_2d_version
    write(std_out, '(a,a)') "Build Type: ", build_type
    write(std_out, '(a,a)') "Compile OS: ", compile_os
    write(std_out, '(a,a)') "Compiler Flags: ", compiler_flags_str
    write(std_out, '(a,a)') "Compiler Version: ", compiler_version_str
    write(std_out, '(a,a)') "Git Hash: ", git_hash
    write(std_out, '(a,a)') "Git Ref: ", git_ref
    write(std_out, '(a,a)') "Commited Changes (Clean/Dirty): ", git_local_changes
    write(std_out, '(a,a)') "Compile Host: ", compile_host
    write(std_out, '(a,a)')
    ! end if
  end subroutine

  subroutine debug_print(str, file, line_number)
    character(len=*), intent(in) :: str
    character(len=*), intent(in), optional :: file
    integer, intent(in), optional :: line_number

    if(enable_debug_print) then
      if(this_image() == 1) then

        if(enable_file_and_line_stats) then
          if(present(file) .and. present(line_number)) then
            write(std_out, '(a,a,i0,3(a))') file, ':', line_number, ' "', str, '"'
          else
            write(std_out, '(a)') str
          end if
        else
          write(std_out, '(a)') str
        end if

      end if
    end if
  end subroutine

end module mod_globals
