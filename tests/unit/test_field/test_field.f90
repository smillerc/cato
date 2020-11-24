module test_field_mod
  use iso_fortran_env, only: rk => real64, ik => int32
  use mod_field
  use mod_parallel
  use caf_testing, only: assert_equal
  implicit none

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

    select case(this_image())
    case(1)
      field%data = 0.0_rk
      field%data(1:4,4) = [5.0_rk, 4.0_rk, 3.0_rk, 2.0_rk]
      field%data(1:4,3) = [9.0_rk, 8.0_rk, 7.0_rk, 6.0_rk]
      field%data(1:4,2) = [5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
      field%data(1:4,1) = [1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk]
    case(2)
      field%data = 0.0_rk
      field%data(5:8,4) = [4.0_rk, 3.0_rk, 2.0_rk, 1.0_rk]
      field%data(5:8,3) = [8.0_rk, 7.0_rk, 6.0_rk, 5.0_rk]
      field%data(5:8,2) = [6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk]
      field%data(5:8,1) = [2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]
    case(3)
      field%data = 0.0_rk
      field%data(1:4,8) = [3.0_rk, 2.0_rk, 1.0_rk, 9.0_rk]
      field%data(1:4,7) = [7.0_rk, 6.0_rk, 5.0_rk, 4.0_rk]
      field%data(1:4,6) = [7.0_rk, 8.0_rk, 9.0_rk, 8.0_rk]
      field%data(1:4,5) = [3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
    case(4)
      field%data = 0.0_rk
      field%data(5:8,8) = [2.0_rk, 1.0_rk, 2.0_rk, 3.0_rk]
      field%data(5:8,7) = [6.0_rk, 5.0_rk, 4.0_rk, 3.0_rk]
      field%data(5:8,6) = [8.0_rk, 9.0_rk, 8.0_rk, 7.0_rk]
      field%data(5:8,5) = [4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk]
    end select

    call field%zero_out_halo

    !      Domain used for halo exchange testing (8x8 grid)
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 0 | 0 | 0 | 0 |   | 0 | 0 | 0 | 0 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 0 | 0 | 0 | 0 |   | 0 | 0 | 0 | 0 || 0 | 0 |
    !    |===|===||===|===|===|===|   |===|===|===|===||===|===|
    !    | 0 | 0 || 3 | 2 | 1 | 9 |   | 2 | 1 | 2 | 3 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 7 | 6 | 5 | 4 |   | 6 | 5 | 4 | 3 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 7 | 8 | 9 | 8 |   | 8 | 9 | 8 | 7 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 3 | 4 | 5 | 6 |   | 4 | 5 | 6 | 7 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|

    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 5 | 4 | 3 | 2 |   | 4 | 3 | 2 | 1 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 9 | 8 | 7 | 6 |   | 8 | 7 | 6 | 5 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 5 | 6 | 7 | 8 |   | 6 | 7 | 8 | 9 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 1 | 2 | 3 | 4 |   | 2 | 3 | 4 | 5 || 0 | 0 | 
    !    |===|===||===|===|===|===|   |===|===|===|===||===|===|
    !    | 0 | 0 || 0 | 0 | 0 | 0 |   | 0 | 0 | 0 | 0 || 0 | 0 |
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|
    !    | 0 | 0 || 0 | 0 | 0 | 0 |   | 0 | 0 | 0 | 0 || 0 | 0 | 
    !    |---|---||---|---|---|---|   |---|---|---|---||---|---|

    call field%sync_edges()
    associate(ilo => field%lbounds(1), ihi => field%ubounds(1), &
              jlo => field%lbounds(2), jhi => field%ubounds(2), &
              ilo_halo => field%lbounds_halo(1), ihi_halo => field%ubounds_halo(1), &
              jlo_halo => field%lbounds_halo(2), jhi_halo => field%ubounds_halo(2), &
              nh => field%n_halo_cells, &
              neighbors => field%neighbors)

      select case(this_image())
      case(1)
        print*, "Image: 1"
        do j = jhi_halo, jlo_halo, -1
          write(*, '( 100(f4.1, 1x))') field%data(:, j)
        end do
        
        ! jhi
        call assert_equal([0.0_rk, 0.0_rk, 7.0_rk, 8.0_rk, 9.0_rk, 8.0_rk, 8.0_rk, 9.0_rk], &
                          field%data(ilo_halo:ihi_halo, jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal([0.0_rk, 0.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 4.0_rk, 5.0_rk], &
                          field%data(ilo_halo:ihi_halo, jhi_halo - 1), file=__FILE__, line=__LINE__)

        ! jlo
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jlo_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jlo_halo + 1), file=__FILE__, line=__LINE__)
      
        ! ihi
        call assert_equal([0.0_rk, 0.0_rk, 3.0_rk, 7.0_rk, 7.0_rk, 3.0_rk, 5.0_rk, 9.0_rk], &
                          field%data(ihi_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal([0.0_rk, 0.0_rk, 2.0_rk, 6.0_rk, 8.0_rk, 4.0_rk, 4.0_rk, 8.0_rk], &
                          field%data(ihi_halo - 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)

        ! ilo
        call assert_equal(0.0_rk, field%data(ilo_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ilo_halo + 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)

      case(2)
        ! jhi
        call assert_equal([9.0_rk, 8.0_rk, 8.0_rk, 9.0_rk, 8.0_rk, 7.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ilo_halo:ihi_halo, jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal([5.0_rk, 6.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ilo_halo:ihi_halo, jhi_halo - 1), file=__FILE__, line=__LINE__)

        ! ihi
        call assert_equal(0.0_rk, field%data(ihi_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ihi_halo - 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        
        ! jlo
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jlo_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jlo_halo + 1), file=__FILE__, line=__LINE__)
        
        ! ilo
        call assert_equal([0.0_rk, 0.0_rk, 3.0_rk, 7.0_rk, 7.0_rk, 3.0_rk, 5.0_rk, 9.0_rk], &
                          field%data(ilo_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal([0.0_rk, 0.0_rk, 4.0_rk, 8.0_rk, 6.0_rk, 2.0_rk, 6.0_rk, 8.0_rk], &
                          field%data(ilo_halo + 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)

      case(3)
        ! jhi
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jhi_halo - 1), file=__FILE__, line=__LINE__)

        ! ihi
        call assert_equal([7.0_rk, 3.0_rk, 5.0_rk, 9.0_rk, 5.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ihi_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal([8.0_rk, 4.0_rk, 4.0_rk, 8.0_rk, 6.0_rk, 2.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ihi_halo - 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        
        ! jlo
        call assert_equal([0.0_rk, 0.0_rk, 9.0_rk, 8.0_rk, 7.0_rk, 6.0_rk, 8.0_rk, 7.0_rk], &
                          field%data(ilo_halo:ihi_halo, jlo_halo), file=__FILE__, line=__LINE__)
        call assert_equal([0.0_rk, 0.0_rk, 5.0_rk, 4.0_rk, 3.0_rk, 2.0_rk, 4.0_rk, 3.0_rk], &
                          field%data(ilo_halo:ihi_halo, jlo_halo + 1), file=__FILE__, line=__LINE__)
        
        ! ilo
        call assert_equal(0.0_rk, field%data(ilo_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ilo_halo + 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
      case(4)
        ! jhi
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ilo_halo:ihi_halo, jhi_halo - 1), file=__FILE__, line=__LINE__)

        ! ihi
        call assert_equal(0.0_rk, field%data(ihi_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, field%data(ihi_halo - 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        
        ! jlo
        call assert_equal([7.0_rk, 6.0_rk, 8.0_rk, 7.0_rk, 6.0_rk, 5.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ilo_halo:ihi_halo, jlo_halo), file=__FILE__, line=__LINE__)
        call assert_equal([3.0_rk, 2.0_rk, 4.0_rk, 3.0_rk, 2.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ilo_halo:ihi_halo, jlo_halo + 1), file=__FILE__, line=__LINE__)
        
        ! ilo
        call assert_equal([7.0_rk, 3.0_rk, 5.0_rk, 9.0_rk, 5.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ilo_halo, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
        call assert_equal([6.0_rk, 2.0_rk, 6.0_rk, 8.0_rk, 4.0_rk, 9.0_rk, 0.0_rk, 0.0_rk], &
                          field%data(ilo_halo + 1, jlo_halo:jhi_halo), file=__FILE__, line=__LINE__)
      end select

    end associate

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success"
  end subroutine test_halo_exchange

  subroutine test_mult_arithmetic()
    type(field_2d_t) :: field

    real(rk), dimension(:, :), allocatable :: x_2d
    real(rk), dimension(ni, nj) :: global
    real(rk) :: x_1d, field_sum

    integer :: i, j

    field = field_2d(name='rho', long_name='density', descrip='Cell-centered mass density', &
                     units='g/cc', global_dims=[ni, nj], n_halo_cells=2)

    associate(ilo => field%lbounds(1), ihi => field%ubounds(1), &
              jlo => field%lbounds(2), jhi => field%ubounds(2))
      allocate(x_2d(ilo:ihi, jlo:jhi))
    end associate

    x_2d = 10.0_rk
    x_1d = 20.0_rk
    
    select case(this_image())
    case (1)
      field%data = 1.0_rk
    case (2)
      field%data = 2.0_rk
    case (3)
      field%data = 3.0_rk
    case (4)
      field%data = 4.0_rk
    end select
    
    select case(this_image())
    case (1)
      field_sum = field%sum()
      call assert_equal(desired=16.0_rk, actual=field_sum, file=__FILE__, line=__LINE__)
    case (2)
      field_sum = field%sum()
      call assert_equal(desired=32.0_rk, actual=field_sum, file=__FILE__, line=__LINE__)
    case (3)
      field_sum = field%sum()
      call assert_equal(desired=48.0_rk, actual=field_sum, file=__FILE__, line=__LINE__)
    case (4)
      field_sum = field%sum()
      call assert_equal(desired=64.0_rk, actual=field_sum, file=__FILE__, line=__LINE__)
    end select
    
    field = field * x_2d
    select case(this_image())
    case (1)
      call assert_equal(desired=160.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (2)
      call assert_equal(desired=320.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (3)
      call assert_equal(desired=480.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (4)
      call assert_equal(desired=640.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    end select

    field = field / 10.0_rk
    select case(this_image())
    case (1)
      call assert_equal(desired=16.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (2)
      call assert_equal(desired=32.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (3)
      call assert_equal(desired=48.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    case (4)
      call assert_equal(desired=64.0_rk, actual=field%sum(), file=__FILE__, line=__LINE__)
    end select

    sync all

    global = field%gather(image=1)
    if (this_image() == 1) then
      call assert_equal(desired=160.0_rk, actual=sum(global), file=__FILE__, line=__LINE__)
    end if

    write(*, '(a, i0, a)') "Image: ", this_image(), " Success"
  end subroutine test_mult_arithmetic
end module test_field_mod

program test_field
  use test_field_mod
  implicit none

  ! sync all
  if(this_image() == 1) print*, new_line('') // "Running test_halo_exchange" // new_line('') 
  call test_halo_exchange()

  ! if (num_images() /= 4) then
  !   error stop "test_field is designed to run with 4 images, currently num_images() /= 4"
  ! end if

  ! sync all
  
  ! if(this_image() == 1) print*, new_line('') // "Running test_mult_arithmetic" //new_line('')
  ! call test_mult_arithmetic()
end program test_field
