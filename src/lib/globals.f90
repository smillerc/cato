module mod_globals
  use iso_fortran_env, only: compiler_options, compiler_version, std_out => output_unit, &
                             ik => int32, rk => real64

  implicit none

  logical, parameter :: enable_debug_print = .false.
  logical, parameter :: enable_file_and_line_stats = .false.

  real(rk), parameter :: TINY_DIST = 5e-16_rk

  logical :: globals_set = .false.

  character(len=5), protected :: global_dimensionality
  logical, protected :: is_1d = .false.
  logical, protected :: is_1d_in_x = .false.
  logical, protected :: is_1d_in_y = .false.
  logical, protected :: is_2d = .true.

  character(:), allocatable :: compiler_flags_str
  character(:), allocatable :: compiler_version_str
  character(:), allocatable :: git_hash
  character(:), allocatable :: git_ref
  character(:), allocatable :: git_local_changes
  character(:), allocatable :: cato_version
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

  subroutine set_domain_dimensionality(dimensionality)
    !< Set the global dimensionality of the domain. This helps the code
    !< make shortcuts elsewhere
    character(len=4), intent(in) :: dimensionality

    select case(trim(dimensionality))
    case('1D_X')
      is_1d = .true.
      is_1d_in_x = .true.
      is_1d_in_y = .false.
      is_2d = .false.
    case('1D_Y')
      is_1d = .true.
      is_1d_in_x = .false.
      is_1d_in_y = .true.
      is_2d = .false.
    case('2D_XY')
      is_1d = .false.
      is_1d_in_x = .true.
      is_1d_in_y = .false.
      is_2d = .true.
    end select
  end subroutine set_domain_dimensionality

  subroutine print_version_stats()
    call set_global_options()
    ! if(this_image() == 1) then
    write(std_out, '(a,a)')
    write(std_out, '(a,a)') "Version: ", cato_version
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

  subroutine print_evolved_cell_data(primitive_variable, i, j, corners, updown_midoints, leftright_midpoints, primitive_vars)

    integer(ik), intent(in) :: i, j
    real(rk), intent(in), dimension(:, :, :) :: corners
    real(rk), intent(in), dimension(:, :, :) :: updown_midoints
    real(rk), intent(in), dimension(:, :, :) :: leftright_midpoints
    real(rk), intent(in), dimension(:, 0:, 0:) :: primitive_vars
    character(len=*) :: primitive_variable
    integer(ik) :: l

    select case(trim(primitive_variable))
    case('rho', 'density')
      l = 1
    case('u', 'x velocity')
      l = 2
    case('l', 'y velocity')
      l = 3
    case('p', 'pressure')
      l = 4
    case default
      error stop 'Invalid primitive variable name'
    end select

    associate(cell_ave=>primitive_vars(l, i, j), &
              C1=>corners(l, i, j), &
              C2=>corners(l, i + 1, j), &
              C3=>corners(l, i + 1, j + 1), &
              C4=>corners(l, i, j + 1), &
              M1=>leftright_midpoints(l, i, j), &
              M2=>updown_midoints(l, i + 1, j), &
              M3=>leftright_midpoints(l, i, j + 1), &
              M4=>updown_midoints(l, i, j))

      write(*, '(a)') "*********************************************"
      write(*, '(a)') "Evolved state for: "//trim(primitive_variable)
      write(*, '(a)') "*********************************************"
      write(*, *)
      write(*, '(16x,es11.3)') M3
      write(*, '(2(es11.3,a))') C4, " C4--------M3--------C3", C3
      write(*, *) "           |                    |"
      write(*, *) "           |                    |"
      write(*, '((12x,a,es11.3,a))') "|    ", cell_ave, "     |"
      write(*, '(2(es11.3,a))') M4, " M4       (i,j)      M2", M2
      write(*, *) "           |                    |"
      write(*, *) "           |                    |"
      write(*, *) "           |                    |"
      write(*, '(2(es11.3,a))') C1, " C1--------M1--------C2", C2
      write(*, '(16x,es11.3)') M1
      write(*, *)
    end associate

    associate(V=>primitive_vars)
      write(*, *) "Cell neighbors: ", trim(primitive_variable)
      write(*, *) "   (i-1,j+1)  |    (i,j+1)   |   (i+1,j+1)"
      write(*, '(2(es12.3,a),es14.3)') V(l, i - 1, j + 1), "   |", V(l, i, j + 1), "  |", V(l, i + 1, j + 1)
      write(*, *) "              |              |"
      write(*, *) "--------------------------------------------"
      write(*, '(a,i0,",",i0,a)') &
        "   (i-1,j)    |    (", i, j, ")     |   (i+1,j)  "
      write(*, '(2(es12.3,a),es14.3)') V(l, i - 1, j), "   |", V(l, i, j), "  |", V(l, i + 1, j)
      write(*, *) "              |              |"
      write(*, *) "--------------------------------------------"
      write(*, *) "   (i-1,j-1)  |    (i,j-1)   |   (i+1,j-1)"
      write(*, '(2(es12.3,a),es14.3)') V(l, i - 1, j - 1), "   |", V(l, i, j - 1), "  |", V(l, i + 1, j - 1)
      write(*, *) "              |              |"
    end associate

  end subroutine print_evolved_cell_data

  subroutine print_recon_data(primitive_variable, i, j, reconstructed_domain, primitive_vars)

    integer(ik), intent(in) :: i, j
    real(rk), intent(in), dimension(:, :, :, 0:, 0:) :: reconstructed_domain
    real(rk), intent(in), dimension(:, 0:, 0:), optional :: primitive_vars
    character(len=*) :: primitive_variable
    integer(ik) :: l
    real(rk) :: cell_ave

    select case(trim(primitive_variable))
    case('rho', 'density')
      l = 1
    case('u', 'x velocity')
      l = 2
    case('l', 'y velocity')
      l = 3
    case('p', 'pressure')
      l = 4
    case default
      error stop 'Invalid primitive variable name'
    end select

    if(present(primitive_vars)) then
      cell_ave = primitive_vars(l, i, j)
    end if

    associate(C1=>reconstructed_domain(l, 1, 1, i, j), &
              C2=>reconstructed_domain(l, 2, 1, i, j), &
              C3=>reconstructed_domain(l, 3, 1, i, j), &
              C4=>reconstructed_domain(l, 4, 1, i, j), &
              M1=>reconstructed_domain(l, 1, 2, i, j), &
              M2=>reconstructed_domain(l, 2, 2, i, j), &
              M3=>reconstructed_domain(l, 3, 2, i, j), &
              M4=>reconstructed_domain(l, 4, 2, i, j))

      write(*, '(a)') "*********************************************"
      write(*, '(a)') "Reconstructed state for: "//trim(primitive_variable)
      write(*, '(a)') "*********************************************"
      write(*, *)
      write(*, '(16x,es11.3)') M3
      write(*, '(2(es11.3,a))') C4, " C4--------M3--------C3", C3
      write(*, *) "           |                    |"
      write(*, *) "           |                    |"
      write(*, '((12x,a,es11.3,a))') "|    ", cell_ave, "     |"
      write(*, '(es11.3,a,i6,",",i6,a,es11.3)') &
        M4, " M4 (", i, j, ")  M2", M2
      write(*, *) "           |                    |"
      write(*, *) "           |                    |"
      write(*, *) "           |                    |"
      write(*, '(2(es11.3,a))') C1, " C1--------M1--------C2", C2
      write(*, '(16x,es11.3)') M1
      write(*, *)
    end associate

    if(present(primitive_vars)) then
      associate(V=>primitive_vars)
        write(*, *) "Cell neighbors: ", trim(primitive_variable)
        write(*, *) "   (i-1,j+1)  |    (i,j+1)   |   (i+1,j+1)"
        write(*, '(2(es12.3,a),es14.3)') V(l, i - 1, j + 1), "   |", V(l, i, j + 1), "  |", V(l, i + 1, j + 1)
        write(*, *) "              |              |"
        write(*, *) "--------------------------------------------"
        write(*, '(a,i0,",",i0,a)') &
          "   (i-1,j)     |   (", i, j, ")    |   (i+1,j)  "
        write(*, '(2(es12.3,a),es14.3)') V(l, i - 1, j), "   |", V(l, i, j), "  |", V(l, i + 1, j)
        write(*, *) "              |              |"
        write(*, *) "--------------------------------------------"
        write(*, *) "   (i-1,j-1)  |    (i,j-1)   |   (i+1,j-1)"
        write(*, '(2(es12.3,a),es14.3)') V(l, i - 1, j - 1), "   |", V(l, i, j - 1), "  |", V(l, i + 1, j - 1)
        write(*, *) "              |              |"
      end associate
    end if

  end subroutine print_recon_data
end module mod_globals
