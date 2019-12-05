module mod_globals
  use iso_fortran_env, only: compiler_options, compiler_version, output_unit

  implicit none

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
    include 'version.h'
    compiler_flags_str = compiler_options()
    compiler_version_str = compiler_version()
  end subroutine set_global_options

  subroutine print_version_stats()
    call set_global_options()
    if(this_image() == 1) then
      write(output_unit, '(a,a)')
      write(output_unit, '(a,a)') "Version: ", fvleg_2d_version
      write(output_unit, '(a,a)') "Build Type: ", build_type
      write(output_unit, '(a,a)') "Compile OS: ", compile_os
      write(output_unit, '(a,a)') "Compiler Flags: ", compiler_flags_str
      write(output_unit, '(a,a)') "Compiler Version: ", compiler_version_str
      write(output_unit, '(a,a)') "Git Hash: ", git_hash
      write(output_unit, '(a,a)') "Git Ref: ", git_ref
      write(output_unit, '(a,a)') "Commited Changes (Clean/Dirty): ", git_local_changes
      write(output_unit, '(a,a)') "Compile Host: ", compile_host
      write(output_unit, '(a,a)')
    end if
  end subroutine

end module mod_globals
