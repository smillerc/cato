module test_field_mod
  use iso_fortran_env
  use mod_field
  use mod_parallel
  use caf_testing, only: assert_equal
  implicit none(type, external)

  integer, parameter :: ni = 8
  integer, parameter :: nj = 8
contains

  subroutine test_halo_exchange()
    type(field_2d_t) :: field
    integer :: i, j

    field = field_2d(name='rho', long_name='density', &
                     descrip='Cell-centered mass density', &
                     units='g/cc', global_dims=[8, 8], n_halo_cells=2)

    call assert_equal(desired=4, actual=field%domain_shape(1), file=__FILE__, line=__LINE__)
    call assert_equal(desired=4, actual=field%domain_shape(2), file=__FILE__, line=__LINE__)
    call assert_equal(desired=2, actual=field%n_halo_cells, file=__FILE__, line=__LINE__)

    if(this_image() == 1) then
      field%data = 1.0_real64
    end if
    if(this_image() == 2) then
      field%data = 2.0_real64
    end if
    if(this_image() == 3) then
      field%data = 3.0_real64
    end if
    if(this_image() == 4) then
      field%data = 4.0_real64
    end if
    call field%zero_out_halo

    associate(ilo => field%lbounds(1), ihi => field%ubounds(1), &
              jlo => field%lbounds(2), jhi => field%ubounds(2), &
              ilo_halo => field%lbounds_halo(1), ihi_halo => field%ubounds_halo(1), &
              jlo_halo => field%lbounds_halo(2), jhi_halo => field%ubounds_halo(2), &
              nh => field%n_halo_cells, &
              neighbors => field%neighbors)

      if(this_image() == 2) then
        do j = jhi_halo, jlo_halo, -1
          write(*, '( 100(f4.1, 1x))') field%data(:, j)
        end do
      end if
    end associate

    call field%sync_edges()
    associate(ilo => field%lbounds(1), ihi => field%ubounds(1), &
              jlo => field%lbounds(2), jhi => field%ubounds(2), &
              ilo_halo => field%lbounds_halo(1), ihi_halo => field%ubounds_halo(1), &
              jlo_halo => field%lbounds_halo(2), jhi_halo => field%ubounds_halo(2), &
              nh => field%n_halo_cells, &
              neighbors => field%neighbors)

      if(this_image() == 2) then
        print *
        do j = jhi_halo, jlo_halo, -1
          write(*, '( 100(f4.1, 1x))') field%data(:, j)
        end do
      end if
    end associate

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success"
  end subroutine test_halo_exchange

  subroutine test_mult_arithmetic()
    type(field_2d_t) :: field

    real(real64), dimension(:, :), allocatable :: x_2d
    real(real64), dimension(ni, nj) :: global
    real(real64) :: x_1d, field_sum

    integer :: i, j

    field = field_2d(name='rho', long_name='density', descrip='Cell-centered mass density', &
                     units='g/cc', global_dims=[ni, nj], n_halo_cells=2)

    associate(ilo => field%lbounds(1), ihi => field%ubounds(1), &
              jlo => field%lbounds(2), jhi => field%ubounds(2))
      allocate(x_2d(ilo:ihi, jlo:jhi))
    end associate

    x_2d = 10.0_real64
    x_1d = 20.0_real64
    
    select case(this_image())
    case (1)
      field%data = 1.0_real64
    case (2)
      field%data = 2.0_real64
    case (3)
      field%data = 3.0_real64
    case (4)
      field%data = 4.0_real64
    end select
    
    select case(this_image())
    case (1)
      field_sum = field%sum()
      call assert_equal(desired=16.0_real64, actual=field_sum, file=__FILE__, line=__LINE__)
    case (2)
      field_sum = field%sum()
      call assert_equal(desired=32.0_real64, actual=field_sum, file=__FILE__, line=__LINE__)
    case (3)
      field_sum = field%sum()
      call assert_equal(desired=48.0_real64, actual=field_sum, file=__FILE__, line=__LINE__)
    case (4)
      field_sum = field%sum()
      call assert_equal(desired=64.0_real64, actual=field_sum, file=__FILE__, line=__LINE__)
    end select
    
    field = field * x_2d
    select case(this_image())
    case (1)
      call assert_equal(desired=160.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (2)
      call assert_equal(desired=320.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (3)
      call assert_equal(desired=480.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (4)
      call assert_equal(desired=640.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    end select

    field = field / 10.0_real64
    select case(this_image())
    case (1)
      call assert_equal(desired=16.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (2)
      call assert_equal(desired=32.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (3)
      call assert_equal(desired=48.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (4)
      call assert_equal(desired=64.0_real64, actual=field%sum(), file=__FILE__, line=__LINE__)
    end select

    sync all

    global = field%gather(image=1)
    if (this_image() == 1) then
      call assert_equal(desired=160.0_real64, actual=sum(global), file=__FILE__, line=__LINE__)
    end if

    write(*, '(a, i0, a)') "Image: ", this_image(), " Success"
  end subroutine test_mult_arithmetic
end module test_field_mod

program test_field
  use test_field_mod
  implicit none(type, external)

  ! sync all
  if(this_image() == 1) print*, new_line('') // "Running test_halo_exchange" // new_line('') 
  call test_halo_exchange()

  if (num_images() /= 4) then
    error stop "test_field is designed to run with 4 images, currently num_images() /= 4"
  end if

  sync all
  
  if(this_image() == 1) print*, new_line('') // "Running test_mult_arithmetic" //new_line('')
  call test_mult_arithmetic()
end program test_field
