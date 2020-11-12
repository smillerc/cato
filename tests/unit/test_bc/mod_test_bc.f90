module mod_test_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  use mod_bc_factory, only: bc_factory
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_grid_factory, only: grid_factory
  use mod_grid_block, only: grid_block_t
  use mod_grid_block_2d, only: grid_block_2d_t, new_2d_grid_block
  use mod_field, only: field_2d_t, field_2d
  use caf_testing

  implicit none

  integer(ik), parameter :: n_ghost_layers = 2
  integer(ik), parameter :: ni_nodes = 5
  integer(ik), parameter :: nj_nodes = 5
  integer(ik), parameter :: ni_cells = ni_nodes - 1
  integer(ik), parameter :: nj_cells = nj_nodes - 1
  integer(ik), parameter :: ilo_c = -n_ghost_layers + 1       !< low index for all cells
  integer(ik), parameter :: ihi_c = ni_cells + n_ghost_layers !<
  integer(ik), parameter :: jlo_c = -n_ghost_layers + 1       !<
  integer(ik), parameter :: jhi_c = nj_cells + n_ghost_layers !<

  real(rk), dimension(ilo_c:ihi_c) :: top_row = 0.0_rk
  real(rk), dimension(ilo_c:ihi_c) :: bottom_row = 0.0_rk
  real(rk), dimension(jlo_c:jhi_c) :: right_col = 0.0_rk
  real(rk), dimension(jlo_c:jhi_c) :: left_col = 0.0_rk

  integer(ik) :: alloc_status

  type(input_t) :: input
  class(grid_block_2d_t), pointer :: grid
  type(field_2d_t) :: rho
  type(field_2d_t) :: u
  type(field_2d_t) :: v
  type(field_2d_t) :: p
  integer(ik), parameter :: bottom_edge = 1 !< cell bottom edge index
  integer(ik), parameter :: right_edge = 2  !< cell right edge index
  integer(ik), parameter :: top_edge = 3    !< cell top edge index
  integer(ik), parameter :: left_edge = 4   !< cell left edge index

    !   Domain used for BC testing (read in from simple.h5)
    !    |---|---||---|---|---|---||---|---|
    ! 6  | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 | <- j=6; aka top_ghost (grid%ubounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    ! 5  | 5 | 2 || 1 | 5 | 5 | 2 || 1 | 5 |
    !    |===|===||===|===|===|===||===|===|
    ! 4  | 7 | 3 || 4 | 7 | 7 | 3 || 4 | 7 | <- j=4; aka top (grid%ubounds(2))
    !    |---|---||---|---|---|---||---|---|
    ! 3  | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 |
    !    |---|---||---|---|---|---||---|---|
    ! 2  | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 |
    !    |---|---||---|---|---|---||---|---|
    ! 1  | 5 | 2 || 1 | 5 | 5 | 2 || 1 | 5 | <- j=1; aka bottom (grid%lbounds(2))
    !    |===|===||===|===|===|===||===|===|
    ! 0  | 7 | 3 || 4 | 7 | 7 | 3 || 4 | 7 |
    !    |---|---||---|---|---|---||---|---|
    ! -1 | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 | <- j=-1; aka bottom_ghost (grid%lbounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    !     -1   0    1   2   3   4    5   6

    ! i0; aka left_ghost, grid%lbounds_halo(1)
    ! i1; aka left, grid%lbounds(1)
    ! i4; aka right, grid%ubounds(1)
    ! i5; aka right_ghost grid%ubounds_halo(1)

contains

  subroutine startup()
    input = input_t(spatial_reconstruction='minmod', flux_solver='AUSMPW+')
    input%initial_condition_file = 'simple.h5'

    grid => new_2d_grid_block(input)

    rho = field_2d(name='', long_name='', descrip='', units='', &
                   global_dims=grid%global_dims, n_halo_cells=n_ghost_layers)
    call rho%read_from_h5(filename='simple.h5', dataset='/density')

    u = field_2d(name='', long_name='', descrip='', units='', &
                 global_dims=grid%global_dims, n_halo_cells=n_ghost_layers)
    call u%read_from_h5(filename='simple.h5', dataset='/x_velocity')

    v = field_2d(name='', long_name='', descrip='', units='', &
                 global_dims=grid%global_dims, n_halo_cells=n_ghost_layers)
    call v%read_from_h5(filename='simple.h5', dataset='/y_velocity')

    p = field_2d(name='', long_name='', descrip='', units='', &
                 global_dims=grid%global_dims, n_halo_cells=n_ghost_layers)
    call p%read_from_h5(filename='simple.h5', dataset='/pressure')

  end subroutine startup

  subroutine cleanup()
    deallocate(grid)
  end subroutine cleanup

  subroutine test_x_periodic()
    !< Test the periodic boundary condition
    class(boundary_condition_t), pointer :: bc_minus_x => null()
    class(boundary_condition_t), pointer :: bc_plus_x => null()
    

    sync all

    if(this_image() == 1) print*, new_line('') // "Running test_x_periodic" // new_line('')
    input%plus_x_bc = 'periodic'
    input%plus_y_bc = 'periodic'
    input%minus_x_bc = 'zero_gradient'
    input%minus_y_bc = 'zero_gradient'

    bc_minus_x => bc_factory(bc_type='periodic', location='-x', input=input, grid=grid, time=0.0_rk)
    bc_plus_x => bc_factory(bc_type='periodic', location='+x', input=input, grid=grid, time=0.0_rk)
    call assert_equal('-x', bc_minus_x%location, file=__FILE__, line=__LINE__)
    call assert_equal('+x', bc_plus_x%location, file=__FILE__, line=__LINE__)

    call rho%zero_out_halo()
    call u%zero_out_halo()
    call v%zero_out_halo()
    call p%zero_out_halo()

    call bc_minus_x%apply(rho=rho, u=u, v=v, p=p)
    call bc_plus_x%apply(rho=rho, u=u, v=v, p=p)

    ! Now test the results
    associate(left => grid%lbounds(1), right => grid%ubounds(1), &
              bottom => grid%lbounds(2), top => grid%ubounds(2), &
              left_ghost => grid%lbounds_halo(1), right_ghost => grid%ubounds_halo(1), &
              bottom_ghost => grid%lbounds_halo(2), top_ghost => grid%ubounds_halo(2))

      select case(num_images())
      case(4)
        select case(this_image())
        case(1)
          call assert_equal([5.0_rk, 9.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([2.0_rk, 6.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
        case(2)
          call assert_equal([5.0_rk, 9.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 8.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
        case(3)
          call assert_equal([9.0_rk, 7.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([6.0_rk, 3.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
        case(4)
          call assert_equal([9.0_rk, 7.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 4.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
        end select
      case(1)
        ! right ghost layer
        right_col = [0.0_rk, 0.0_rk, 1.0_rk, 8.0_rk, 8.0_rk, 4.0_rk, 0.0_rk, 0.0_rk]
        call assert_equal(right_col, rho%data(right_ghost - 1, :), file=__FILE__, line=__LINE__)

        right_col = [0.0_rk, 0.0_rk, 5.0_rk, 9.0_rk, 9.0_rk, 7.0_rk, 0.0_rk, 0.0_rk]
        call assert_equal(right_col, rho%data(right_ghost, :), file=__FILE__, line=__LINE__)
      case default
        error stop "Test failed. This is only configured to test for 1 or 4 images"
      end select
    end associate

    deallocate(bc_plus_x)
    deallocate(bc_minus_x)

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success" // " [test_plus_x_periodic]"
  end subroutine test_x_periodic

  subroutine test_y_periodic()
    !< Test the periodic boundary condition
    class(boundary_condition_t), pointer :: bc_plus_y => null()
    class(boundary_condition_t), pointer :: bc_minus_y => null()
    
    

    sync all
    if(this_image() == 1) print*, new_line('') // "Running test_y_periodic" // new_line('')
    
    input%plus_x_bc = 'periodic'
    input%plus_y_bc = 'periodic'
    input%minus_x_bc = 'zero_gradient'
    input%minus_y_bc = 'zero_gradient'

    bc_minus_y => bc_factory(bc_type='periodic', location='-y', input=input, grid=grid, time=0.0_rk)
    bc_plus_y => bc_factory(bc_type='periodic', location='+y', input=input, grid=grid, time=0.0_rk)
    call assert_equal('-y', bc_minus_y%location, file=__FILE__, line=__LINE__)
    call assert_equal('+y', bc_plus_y%location, file=__FILE__, line=__LINE__)

    call rho%zero_out_halo()
    call u%zero_out_halo()
    call v%zero_out_halo()
    call p%zero_out_halo()

    call bc_minus_y%apply(rho=rho, u=u, v=v, p=p)
    call bc_plus_y%apply(rho=rho, u=u, v=v, p=p)
    sync all

    ! Now test the results
    associate(left => grid%lbounds(1), right => grid%ubounds(1), &
              bottom => grid%lbounds(2), top => grid%ubounds(2), &
              left_ghost => grid%lbounds_halo(1), right_ghost => grid%ubounds_halo(1), &
              bottom_ghost => grid%lbounds_halo(2), top_ghost => grid%ubounds_halo(2))

      select case(num_images())
      case(4)
        select case(this_image())
        case(1)
          call assert_equal([4.0_rk, 7.0_rk], rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 9.0_rk], rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
        case(2)
          call assert_equal([7.0_rk, 3.0_rk], rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk], rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
        case(3)
          call assert_equal([8.0_rk, 9.0_rk], rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 5.0_rk], rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        case(4)
          call assert_equal([9.0_rk, 6.0_rk], rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 2.0_rk], rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        end select
      case(1)
        ! top ghost layer
        top_row = [0.0_rk, 0.0_rk, 8.0_rk, 9.0_rk, 9.0_rk, 6.0_rk, 0.0_rk, 0.0_rk]
        call assert_equal(top_row, rho%data(:, top_ghost), file=__FILE__, line=__LINE__)
  
        top_row = [0.0_rk, 0.0_rk, 1.0_rk, 5.0_rk, 5.0_rk, 2.0_rk, 0.0_rk, 0.0_rk]
        call assert_equal(top_row, rho%data(:, top_ghost - 1), file=__FILE__, line=__LINE__)
      case default
        error stop "Test failed. This is only configured to test for 1 or 4 images"
      end select
    end associate

    deallocate(bc_minus_y)
    deallocate(bc_plus_y)

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success" // " [test_plus_y_periodic]"
  end subroutine test_y_periodic

  subroutine test_periodic_all()
    !< Test the periodic boundary condition
    class(boundary_condition_t), pointer :: bc_plus_x => null()
    class(boundary_condition_t), pointer :: bc_plus_y => null()
    class(boundary_condition_t), pointer :: bc_minus_x => null()
    class(boundary_condition_t), pointer :: bc_minus_y => null()
    
    integer :: j

    sync all
    if(this_image() == 1) print*, new_line('') // "Running test_periodic_all" // new_line('')
    input%plus_x_bc = 'periodic'
    input%plus_y_bc = 'periodic'
    input%minus_x_bc = 'periodic'
    input%minus_y_bc = 'periodic'
    
    bc_plus_x => bc_factory(bc_type='periodic', location='+x', input=input, grid=grid, time=0.0_rk)
    bc_plus_y => bc_factory(bc_type='periodic', location='+y', input=input, grid=grid, time=0.0_rk)
    bc_minus_x => bc_factory(bc_type='periodic', location='-x', input=input, grid=grid, time=0.0_rk)
    bc_minus_y => bc_factory(bc_type='periodic', location='-y', input=input, grid=grid, time=0.0_rk)

    call assert_equal('+x', bc_plus_x%location, file=__FILE__, line=__LINE__)
    call assert_equal('+y', bc_plus_y%location, file=__FILE__, line=__LINE__)
    call assert_equal('-x', bc_minus_x%location, file=__FILE__, line=__LINE__)
    call assert_equal('-y', bc_minus_y%location, file=__FILE__, line=__LINE__)

    call rho%zero_out_halo()
    call u%zero_out_halo()
    call v%zero_out_halo()
    call p%zero_out_halo()

    call bc_plus_x%apply(rho=rho, u=u, v=v, p=p)
    call bc_plus_y%apply(rho=rho, u=u, v=v, p=p)
    call bc_minus_x%apply(rho=rho, u=u, v=v, p=p)
    call bc_minus_y%apply(rho=rho, u=u, v=v, p=p)

    !    Expected domain after periodic BCs are applied
    !    |---|---||---|---|---|---||---|---|
    ! 6  | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 | <- j=6; aka top_ghost (grid%ubounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    ! 5  | 5 | 2 || 1 | 5 | 5 | 2 || 1 | 5 |
    !    |===|===||===|===|===|===||===|===|
    ! 4  | 7 | 3 || 4 | 7 | 7 | 3 || 4 | 7 | <- j=4; aka top (grid%ubounds(2))
    !    |---|---||---|---|---|---||---|---|
    ! 3  | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 |
    !    |---|---||---|---|---|---||---|---|
    ! 2  | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 |
    !    |---|---||---|---|---|---||---|---|
    ! 1  | 5 | 2 || 1 | 5 | 5 | 2 || 1 | 5 | <- j=1; aka bottom (grid%lbounds(2))
    !    |===|===||===|===|===|===||===|===|
    ! 0  | 7 | 3 || 4 | 7 | 7 | 3 || 4 | 7 |
    !    |---|---||---|---|---|---||---|---|
    ! -1 | 9 | 6 || 8 | 9 | 9 | 6 || 8 | 9 | <- j=-1; aka bottom_ghost (grid%lbounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    !     -1   0    1   2   3   4    5   6

    ! Now test the results
    associate(left => grid%lbounds(1), right => grid%ubounds(1), &
              bottom => grid%lbounds(2), top => grid%ubounds(2), &
              left_ghost => grid%lbounds_halo(1), right_ghost => grid%ubounds_halo(1), &
              bottom_ghost => grid%lbounds_halo(2), top_ghost => grid%ubounds_halo(2))

      select case(num_images())
      case(4)
        select case(this_image())
        case(1) ! lower left image
          call assert_equal([5.0_rk, 9.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([2.0_rk, 6.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([7.0_rk, 3.0_rk, 4.0_rk, 7.0_rk], rho%data(left_ghost:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk], rho%data(left_ghost:right, bottom_ghost), file=__FILE__, line=__LINE__)
        case(2) ! lower right image
          call assert_equal([5.0_rk, 9.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 8.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([7.0_rk, 3.0_rk, 4.0_rk, 7.0_rk], rho%data(left:right_ghost, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk], rho%data(left:right_ghost, bottom_ghost), file=__FILE__, line=__LINE__)
        case(3) ! upper left image
          call assert_equal([9.0_rk, 7.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([6.0_rk, 3.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk], rho%data(left_ghost:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 2.0_rk, 1.0_rk, 5.0_rk], rho%data(left_ghost:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        case(4) ! upper right image
          call assert_equal([9.0_rk, 7.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 4.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk], rho%data(left:right_ghost, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 2.0_rk, 1.0_rk, 5.0_rk], rho%data(left:right_ghost, top_ghost - 1), file=__FILE__, line=__LINE__)
        end select
      case(1)
        print *, "After BC's applied"
        do j = top_ghost, bottom_ghost, -1
          write(*, '(i3, a, 10(f6.0))') j, " : ", rho%data(:, j)
        end do
        print *

        ! bottom ghost layers
        bottom_row = [7.0_rk, 3.0_rk, 4.0_rk, 7.0_rk, 7.0_rk, 3.0_rk, 4.0_rk, 7.0_rk]
        call assert_equal(bottom_row, rho%data(:, bottom_ghost + 1), file=__FILE__, line=__LINE__)

        bottom_row = [9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk, 9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk]
        call assert_equal(bottom_row, rho%data(:, bottom_ghost), file=__FILE__, line=__LINE__)

        ! top ghost layer
        top_row = [9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk, 9.0_rk, 6.0_rk, 8.0_rk, 9.0_rk]
        call assert_equal(top_row, rho%data(:, top_ghost), file=__FILE__, line=__LINE__)

        top_row = [5.0_rk, 2.0_rk, 1.0_rk, 5.0_rk, 5.0_rk, 2.0_rk, 1.0_rk, 5.0_rk]
        call assert_equal(top_row, rho%data(:, top_ghost - 1), file=__FILE__, line=__LINE__)

        ! right ghost layer
        right_col = [8.0_rk, 4.0_rk, 1.0_rk, 8.0_rk, 8.0_rk, 4.0_rk, 1.0_rk, 8.0_rk]
        call assert_equal(right_col, rho%data(right_ghost - 1, :), file=__FILE__, line=__LINE__)

        right_col = [9.0_rk, 7.0_rk, 5.0_rk, 9.0_rk, 9.0_rk, 7.0_rk, 5.0_rk, 9.0_rk]
        call assert_equal(right_col, rho%data(right_ghost, :), file=__FILE__, line=__LINE__)

        ! left ghost layer
        left_col = [6.0_rk, 3.0_rk, 2.0_rk, 6.0_rk, 6.0_rk, 3.0_rk, 2.0_rk, 6.0_rk]
        call assert_equal(left_col, rho%data(left_ghost + 1, :), file=__FILE__, line=__LINE__)

        left_col = [9.0_rk, 7.0_rk, 5.0_rk, 9.0_rk, 9.0_rk, 7.0_rk, 5.0_rk, 9.0_rk]
        call assert_equal(left_col, rho%data(left_ghost, :), file=__FILE__, line=__LINE__)
      case default
        error stop "Test failed. This is only configured to test for 1 or 4 images"
      end select

    end associate

    deallocate(bc_plus_x)
    deallocate(bc_plus_y)
    deallocate(bc_minus_x)
    deallocate(bc_minus_y)

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success" // " [test_periodic_all]"
  end subroutine test_periodic_all

  subroutine test_symmetry()
    !< Test the symmetry boundary condition
    class(boundary_condition_t), pointer :: bc_plus_x => null()
    class(boundary_condition_t), pointer :: bc_plus_y => null()
    class(boundary_condition_t), pointer :: bc_minus_x => null()
    class(boundary_condition_t), pointer :: bc_minus_y => null()
    integer(ik) :: top_ghost, left_ghost, right_ghost, bottom_ghost
    
    integer :: j

    real(rk), dimension(4) ::    top_sym_row = 0.0_rk
    real(rk), dimension(4) :: bottom_sym_row = 0.0_rk
    real(rk), dimension(4) ::  right_sym_col = 0.0_rk
    real(rk), dimension(4) ::   left_sym_col = 0.0_rk

    sync all
    if(this_image() == 1) print*, new_line('') // "Running test_symmetry" // new_line('')

    bottom_ghost = grid%lbounds_halo(2)
    top_ghost = grid%ubounds_halo(2)

    left_ghost = grid%lbounds_halo(1)
    right_ghost = grid%ubounds_halo(1)

    bc_plus_x => bc_factory(bc_type='symmetry', location='+x', input=input, grid=grid, time=0.0_rk)
    bc_plus_y => bc_factory(bc_type='symmetry', location='+y', input=input, grid=grid, time=0.0_rk)
    bc_minus_x => bc_factory(bc_type='symmetry', location='-x', input=input, grid=grid, time=0.0_rk)
    bc_minus_y => bc_factory(bc_type='symmetry', location='-y', input=input, grid=grid, time=0.0_rk)

    call assert_equal('+x', bc_plus_x%location, file=__FILE__, line=__LINE__)
    call assert_equal('+y', bc_plus_y%location, file=__FILE__, line=__LINE__)
    call assert_equal('-x', bc_minus_x%location, file=__FILE__, line=__LINE__)
    call assert_equal('-y', bc_minus_y%location, file=__FILE__, line=__LINE__)



    !   Domain cartoon with symmetry applied to all the boundaries
    !   |---|---||---|---|---|---||---|---|
    !   | - | - || 8 | 9 | 9 | 6 || - | - | <- j=6; aka top_ghost (grid%ubounds_halo(2))
    !   |---|---||---|---|---|---||---|---|
    !   | - | - || 4 | 7 | 7 | 3 || - | - |
    !   |===|===||===|===|===|===||===|===|
    !   | 7 | 4 || 4 | 7 | 7 | 3 || 3 | 7 | <- j=4; aka top (grid%ubounds(2))
    !   |---|---||---|---|---|---||---|---|
    !   | 9 | 8 || 8 | 9 | 9 | 6 || 6 | 9 |
    !   |---|---||---|---|---|---||---|---|
    !   | 9 | 8 || 8 | 9 | 9 | 6 || 6 | 9 |
    !   |---|---||---|---|---|---||---|---|
    !   | 5 | 1 || 1 | 5 | 5 | 2 || 2 | 5 | <- j=1; aka bottom (grid%lbounds(2))
    !   |===|===||===|===|===|===||===|===|
    !   | - | - || 1 | 5 | 5 | 2 || - | - |
    !   |---|---||---|---|---|---||---|---|
    !   | - | - || 8 | 9 | 9 | 6 || - | - | <- j=-1; aka bottom_ghost (grid%lbounds_halo(2))
    !   |---|---||---|---|---|---||---|---|
    !    i0   i1  i2  i3  i4   i5

    ! i0; aka left_ghost, grid%lbounds_halo(1)
    ! i1; aka left, grid%lbounds(1)
    ! i4; aka right, grid%ubounds(1)
    ! i5; aka right_ghost grid%ubounds_halo(1)

    call rho%zero_out_halo()
    call u%zero_out_halo()
    call v%zero_out_halo()
    call p%zero_out_halo()

    call bc_plus_x%apply(rho=rho, u=u, v=v, p=p)
    call bc_plus_y%apply(rho=rho, u=u, v=v, p=p)
    call bc_minus_x%apply(rho=rho, u=u, v=v, p=p)
    call bc_minus_y%apply(rho=rho, u=u, v=v, p=p)

    ! Now test the results
    associate(left => grid%lbounds(1), right => grid%ubounds(1), &
              bottom => grid%lbounds(2), top => grid%ubounds(2), &
              left_ghost => grid%lbounds_halo(1), right_ghost => grid%ubounds_halo(1), &
              bottom_ghost => grid%lbounds_halo(2), top_ghost => grid%ubounds_halo(2))

      select case(num_images())
      case(4)
        select case(this_image())
        case(1) ! lower left image
          ! jlo side
          call assert_equal([1.0_rk, 5.0_rk], rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 9.0_rk], rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 5.0_rk], u%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 9.0_rk], u%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
          call assert_equal(-[1.0_rk, 5.0_rk], v%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal(-[8.0_rk, 9.0_rk], v%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 5.0_rk], p%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 9.0_rk], p%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)

          ! ilo side
          call assert_equal([5.0_rk, 9.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 8.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[5.0_rk, 9.0_rk], u%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[1.0_rk, 8.0_rk], u%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 9.0_rk], v%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 8.0_rk], v%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 9.0_rk], p%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 8.0_rk], p%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
        case(2) ! lower right image
          ! ihi side
          call assert_equal([5.0_rk, 9.0_rk],  rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([2.0_rk, 6.0_rk],  rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[5.0_rk, 9.0_rk], u%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[2.0_rk, 6.0_rk], u%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 9.0_rk],  v%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([2.0_rk, 6.0_rk],  v%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 9.0_rk],  p%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([2.0_rk, 6.0_rk],  p%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)

          ! jlo side
          call assert_equal([5.0_rk, 2.0_rk],  rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk],  rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 2.0_rk],  u%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk],  u%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
          call assert_equal(-[5.0_rk, 2.0_rk], v%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal(-[9.0_rk, 6.0_rk], v%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 2.0_rk],  p%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk],  p%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)

        case(3) ! upper left image
          ! ilo side
          call assert_equal([9.0_rk, 7.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 4.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[9.0_rk, 7.0_rk], u%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[8.0_rk, 4.0_rk], u%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 7.0_rk], v%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 4.0_rk], v%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 7.0_rk], p%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 4.0_rk], p%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)

          ! jhi side
          call assert_equal([8.0_rk, 9.0_rk], rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([4.0_rk, 7.0_rk], rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 9.0_rk], u%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([4.0_rk, 7.0_rk], u%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
          call assert_equal(-[8.0_rk, 9.0_rk], v%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal(-[4.0_rk, 7.0_rk], v%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 9.0_rk], p%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([4.0_rk, 7.0_rk], p%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)

        case(4) ! upper right image

          ! jhi side
          call assert_equal([9.0_rk, 6.0_rk], rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([7.0_rk, 3.0_rk], rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk], u%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([7.0_rk, 3.0_rk], u%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
          call assert_equal(-[9.0_rk, 6.0_rk], v%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal(-[7.0_rk, 3.0_rk], v%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 6.0_rk], p%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([7.0_rk, 3.0_rk], p%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)

          ! ihi side
          call assert_equal([9.0_rk, 7.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([6.0_rk, 3.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[9.0_rk, 7.0_rk], u%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(-[6.0_rk, 3.0_rk], u%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 7.0_rk], v%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([6.0_rk, 3.0_rk], v%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([9.0_rk, 7.0_rk], p%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([6.0_rk, 3.0_rk], p%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
        end select
      case(1)

        ! Print out in all it's glory
        print *, "After BC's applied"
        do j = top_ghost, bottom_ghost, -1
          write(*, '(i3, a, 10(f6.0))') j, " : ", u%data(:, j)
        end do
        print *


        ! -- Bottom --
        bottom_sym_row = [1.0_rk, 5.0_rk, 5.0_rk, 2.0_rk]
        call assert_equal(bottom_sym_row, rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
        call assert_equal(bottom_sym_row, u%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
        call assert_equal(-bottom_sym_row, v%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
        call assert_equal(bottom_sym_row, p%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)

        ! 2nd ghost layer (down 1 layer)
        bottom_sym_row = [8.0_rk, 9.0_rk, 9.0_rk, 6.0_rk]
        call assert_equal(bottom_sym_row, rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
        call assert_equal(bottom_sym_row, u%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
        call assert_equal(-bottom_sym_row, v%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
        call assert_equal(bottom_sym_row, p%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)

        ! Top --
        top_sym_row = [8.0_rk, 9.0_rk, 9.0_rk, 6.0_rk]
        call assert_equal(top_sym_row, rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
        call assert_equal(top_sym_row, u%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
        call assert_equal(-top_sym_row, v%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
        call assert_equal(top_sym_row, p%data(left:right, top_ghost), file=__FILE__, line=__LINE__)

        ! 2nd ghost layer (down 1 layer)
        top_sym_row = [4.0_rk, 7.0_rk, 7.0_rk, 3.0_rk]
        call assert_equal(top_sym_row, rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        call assert_equal(top_sym_row, u%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        call assert_equal(-top_sym_row, v%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        call assert_equal(top_sym_row, p%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)

        ! -- Right --
        right_sym_col = [2.0_rk, 6.0_rk, 6.0_rk, 3.0_rk]
        call assert_equal(right_sym_col, rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(-right_sym_col, u%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(right_sym_col, v%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(right_sym_col, p%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)

        ! 2nd ghost layer
        right_sym_col = [5.0_rk, 9.0_rk, 9.0_rk, 7.0_rk]
        call assert_equal(right_sym_col, rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(-right_sym_col, u%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(right_sym_col, v%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(right_sym_col, p%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)

        ! -- Left --
        left_sym_col = [1.0_rk, 8.0_rk, 8.0_rk, 4.0_rk]
        call assert_equal(left_sym_col, rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(-left_sym_col, u%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(left_sym_col, v%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(left_sym_col, p%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)

        ! 2nd ghost layer
        left_sym_col = [5.0_rk, 9.0_rk, 9.0_rk, 7.0_rk]
        call assert_equal(left_sym_col, rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(-left_sym_col, u%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(left_sym_col, v%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal(left_sym_col, p%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
      case default
        error stop "Test failed. This is only configured to test for 1 or 4 images"
      end select
    end associate

    deallocate(bc_plus_x)
    deallocate(bc_plus_y)
    deallocate(bc_minus_x)
    deallocate(bc_minus_y)

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success" // " [test_symmetry]"
  end subroutine test_symmetry

  subroutine test_zero_gradient()
    !< Test the periodic boundary condition
    class(boundary_condition_t), pointer :: bc_plus_x => null()
    class(boundary_condition_t), pointer :: bc_plus_y => null()
    class(boundary_condition_t), pointer :: bc_minus_x => null()
    class(boundary_condition_t), pointer :: bc_minus_y => null()
    
    

    if(this_image() == 1) print*, new_line('') // "Running test_zero_gradient" // new_line('')
    input%plus_x_bc = 'zero_gradient'
    input%plus_y_bc = 'zero_gradient'
    input%minus_x_bc = 'zero_gradient'
    input%minus_y_bc = 'zero_gradient'
    
    bc_plus_x => bc_factory(bc_type='zero_gradient', location='+x', input=input, grid=grid, time=0.0_rk)
    bc_plus_y => bc_factory(bc_type='zero_gradient', location='+y', input=input, grid=grid, time=0.0_rk)
    bc_minus_x => bc_factory(bc_type='zero_gradient', location='-x', input=input, grid=grid, time=0.0_rk)
    bc_minus_y => bc_factory(bc_type='zero_gradient', location='-y', input=input, grid=grid, time=0.0_rk)

    call assert_equal('+x', bc_plus_x%location, file=__FILE__, line=__LINE__)
    call assert_equal('+y', bc_plus_y%location, file=__FILE__, line=__LINE__)
    call assert_equal('-x', bc_minus_x%location, file=__FILE__, line=__LINE__)
    call assert_equal('-y', bc_minus_y%location, file=__FILE__, line=__LINE__)

    call rho%zero_out_halo()
    call u%zero_out_halo()
    call v%zero_out_halo()
    call p%zero_out_halo()

    call bc_plus_x%apply(rho=rho, u=u, v=v, p=p)
    call bc_plus_y%apply(rho=rho, u=u, v=v, p=p)
    call bc_minus_x%apply(rho=rho, u=u, v=v, p=p)
    call bc_minus_y%apply(rho=rho, u=u, v=v, p=p)

    !    Expected domain after periodic BCs are applied
    !    |---|---||---|---|---|---||---|---|
    ! 6  | 0 | 0 || 4 | 7 | 7 | 3 || 0 | 0 | <- j=6; aka top_ghost (grid%ubounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    ! 5  | 0 | 0 || 4 | 7 | 7 | 3 || 0 | 0 |
    !    |===|===||===|===|===|===||===|===|
    ! 4  | 4 | 4 || 4 | 7 | 7 | 3 || 3 | 3 | <- j=4; aka top (grid%ubounds(2))
    !    |---|---||---|---|---|---||---|---|
    ! 3  | 8 | 8 || 8 | 9 | 9 | 6 || 6 | 6 |
    !    |---|---||---|---|---|---||---|---|
    ! 2  | 8 | 8 || 8 | 9 | 9 | 6 || 6 | 6 |
    !    |---|---||---|---|---|---||---|---|
    ! 1  | 1 | 1 || 1 | 5 | 5 | 2 || 2 | 2 | <- j=1; aka bottom (grid%lbounds(2))
    !    |===|===||===|===|===|===||===|===|
    ! 0  | 0 | 0 || 1 | 5 | 5 | 2 || 0 | 0 |
    !    |---|---||---|---|---|---||---|---|
    ! -1 | 0 | 0 || 1 | 5 | 5 | 2 || 0 | 0 | <- j=-1; aka bottom_ghost (grid%lbounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    !     -1   0    1   2   3   4    5   6

    ! Now test the results
    associate(left => grid%lbounds(1), right => grid%ubounds(1), &
              bottom => grid%lbounds(2), top => grid%ubounds(2), &
              left_ghost => grid%lbounds_halo(1), right_ghost => grid%ubounds_halo(1), &
              bottom_ghost => grid%lbounds_halo(2), top_ghost => grid%ubounds_halo(2))

      select case(num_images())
      case(4)
        select case(this_image())
        case(1) ! lower left image
          call assert_equal([1.0_rk, 8.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 8.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 5.0_rk], rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([1.0_rk, 5.0_rk], rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
        case(2) ! lower right image
          call assert_equal([2.0_rk, 6.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([2.0_rk, 6.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 2.0_rk], rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
          call assert_equal([5.0_rk, 2.0_rk], rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)
        case(3) ! upper left image
          call assert_equal([8.0_rk, 4.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([8.0_rk, 4.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([4.0_rk, 7.0_rk], rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([4.0_rk, 7.0_rk], rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        case(4) ! upper right image
          call assert_equal([6.0_rk, 3.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([6.0_rk, 3.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal([7.0_rk, 3.0_rk], rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
          call assert_equal([7.0_rk, 3.0_rk], rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)
        end select
      case(1)
        ! bottom ghost layers
        call assert_equal([1.0_rk, 5.0_rk, 5.0_rk, 2.0_rk], rho%data(left:right, bottom_ghost + 1), file=__FILE__, line=__LINE__)
        call assert_equal([1.0_rk, 5.0_rk, 5.0_rk, 2.0_rk], rho%data(left:right, bottom_ghost), file=__FILE__, line=__LINE__)

        ! top ghost layer
        call assert_equal([4.0_rk, 7.0_rk, 7.0_rk, 3.0_rk], rho%data(left:right, top_ghost), file=__FILE__, line=__LINE__)
        call assert_equal([4.0_rk, 7.0_rk, 7.0_rk, 3.0_rk], rho%data(left:right, top_ghost - 1), file=__FILE__, line=__LINE__)

        ! right ghost layer
        call assert_equal([2.0_rk, 6.0_rk, 6.0_rk, 3.0_rk], rho%data(right_ghost - 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal([2.0_rk, 6.0_rk, 6.0_rk, 3.0_rk], rho%data(right_ghost, bottom:top), file=__FILE__, line=__LINE__)

        ! left ghost layer
        call assert_equal([1.0_rk, 8.0_rk, 8.0_rk, 4.0_rk], rho%data(left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
        call assert_equal([1.0_rk, 8.0_rk, 8.0_rk, 4.0_rk], rho%data(left_ghost, bottom:top), file=__FILE__, line=__LINE__)
      case default
        error stop "Test failed. This is only configured to test for 1 or 4 images"
      end select

    end associate

    deallocate(bc_plus_x)
    deallocate(bc_plus_y)
    deallocate(bc_minus_x)
    deallocate(bc_minus_y)

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success" // " [test_periodic_all]"
  end subroutine test_zero_gradient

  subroutine test_pressure()
    !< Test the periodic boundary condition
    class(boundary_condition_t), pointer :: bc_plus_x => null()
    
    integer :: j

    sync all
    if(this_image() == 1) print*, new_line('') // "Running test_pressure" // new_line('')
    input%plus_x_bc = 'pressure_input'
    
    input%bc_density = 1.0_rk
    input%apply_constant_bc_pressure = .true.
    input%bc_pressure_scale_factor = 1.0_rk
    input%constant_bc_pressure_value = 4.0_rk

    bc_plus_x => bc_factory(bc_type='pressure_input', location='+x', input=input, grid=grid, time=0.0_rk)

    call assert_equal('+x', bc_plus_x%location, file=__FILE__, line=__LINE__)

    call rho%zero_out_halo()
    call u%zero_out_halo()
    call v%zero_out_halo()
    call p%zero_out_halo()

    call bc_plus_x%apply(rho=rho, u=u, v=v, p=p)

    !    Expected domain after periodic BCs are applied
    !    The +x only works for now... rho is 1 and p is 4 for inputs (showing p below)
    !    |---|---||---|---|---|---||---|---|
    ! 6  | 0 | 0 || 0 | 0 | 0 | 0 || 0 | 0 | <- j=6; aka top_ghost (grid%ubounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    ! 5  | 0 | 0 || 0 | 0 | 0 | 0 || 0 | 0 |
    !    |===|===||===|===|===|===||===|===|
    ! 4  | 0 | 0 || 4 | 7 | 7 | 3 || 4 | 4 | <- j=4; aka top (grid%ubounds(2))
    !    |---|---||---|---|---|---||---|---|
    ! 3  | 0 | 0 || 8 | 9 | 9 | 6 || 4 | 4 |
    !    |---|---||---|---|---|---||---|---|
    ! 2  | 0 | 0 || 8 | 9 | 9 | 6 || 4 | 4 |
    !    |---|---||---|---|---|---||---|---|
    ! 1  | 0 | 0 || 1 | 5 | 5 | 2 || 4 | 4 | <- j=1; aka bottom (grid%lbounds(2))
    !    |===|===||===|===|===|===||===|===|
    ! 0  | 0 | 0 || 0 | 0 | 0 | 0 || 0 | 0 |
    !    |---|---||---|---|---|---||---|---|
    ! -1 | 0 | 0 || 0 | 0 | 0 | 0 || 0 | 0 | <- j=-1; aka bottom_ghost (grid%lbounds_halo(2))
    !    |---|---||---|---|---|---||---|---|
    !     -1   0    1   2   3   4    5   6

    ! Now test the results
    associate(left => grid%lbounds(1), right => grid%ubounds(1), &
              bottom => grid%lbounds(2), top => grid%ubounds(2), &
              left_ghost => grid%lbounds_halo(1), right_ghost => grid%ubounds_halo(1), &
              bottom_ghost => grid%lbounds_halo(2), top_ghost => grid%ubounds_halo(2))

      select case(num_images())
      case(4)
        select case(this_image())
        case(1) ! lower left image
          call assert_equal(0.0_rk, rho%data(left_ghost:left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, rho%data(left:right, bottom_ghost:bottom_ghost + 1), file=__FILE__, line=__LINE__)
        case(2) ! lower right image
          call assert_equal(0.0_rk, rho%data(right_ghost - 1:right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, rho%data(left:right, bottom_ghost:bottom_ghost + 1), file=__FILE__, line=__LINE__)
        case(3) ! upper left image
          call assert_equal(0.0_rk, rho%data(left_ghost:left_ghost + 1, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, rho%data(left:right, top_ghost - 1:top_ghost), file=__FILE__, line=__LINE__)
        case(4) ! upper right image
          call assert_equal(0.0_rk, rho%data(right_ghost - 1:right_ghost, bottom:top), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, rho%data(left:right, top_ghost - 1:top_ghost), file=__FILE__, line=__LINE__)
        end select
      case(1)
        ! Print out in all it's glory
        print *, "After pressure BC applied"
        do j = top_ghost, bottom_ghost, -1
          write(*, '(i3, a, 10(f6.0))') j, " : ", p%data(:, j)
        end do
        print *

        ! jlo
        call assert_equal(0.0_rk, rho%data(:, bottom_ghost:bottom_ghost + 1), file=__FILE__, line=__LINE__)
        
        ! jhi
        call assert_equal(0.0_rk, rho%data(:, top_ghost - 1:top_ghost), file=__FILE__, line=__LINE__)
        
        ! ihi
        call assert_equal(0.0_rk, rho%data(right_ghost - 1:right_ghost, :), file=__FILE__, line=__LINE__)
        
        ! ilo
        call assert_equal(0.0_rk, rho%data(left_ghost:left_ghost + 1, :), file=__FILE__, line=__LINE__)
      case default
        error stop "Test failed. This is only configured to test for 1 or 4 images"
      end select

    end associate

    deallocate(bc_plus_x)

    sync all
    write(*, '(a, i0, a)') "Image: ", this_image(), " Success" // " [test_periodic_all]"
  end subroutine test_pressure

  function gather(f, image)
    !< Performs a gather of field data to image.
    class(field_2d_t), intent(in) :: f
    integer(ik), intent(in) :: image
    real(rk) :: gather(f%global_dims(1), f%global_dims(2))

    real(rk), allocatable :: gather_coarray(:, :)[:]
    allocate(gather_coarray(f%global_dims(1), f%global_dims(2))[*])

    associate(is => f%lbounds(1), ie => f%ubounds(1), &
              js => f%lbounds(2), je => f%ubounds(2))
      gather_coarray(is:ie, js:je)[image] = f%data(is:ie, js:je)
      sync all
      if(this_image() == image) gather = gather_coarray
    end associate
    deallocate(gather_coarray)
  end function gather
end module mod_test_bc