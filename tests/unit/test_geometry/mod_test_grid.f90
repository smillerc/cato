module mod_test_grid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, stdout => output_unit
  use mod_grid_block, only: grid_block_t
  use mod_grid_block_2d, only: grid_block_2d_t, new_2d_grid_block
  use mod_input, only: input_t
  use caf_testing
  implicit none

  type(input_t) :: input
  class(grid_block_2d_t), pointer :: grid
  integer(ik), parameter :: bottom = 1
  integer(ik), parameter :: right = 2
  integer(ik), parameter :: top = 3
  integer(ik), parameter :: left = 4

    !   Domain used for testing (x and y coordinates of the nodes are shown).
    !   H -> halo cell, I2 -> CAF image #2
    !
    !  4 |--------|--------||--------|--------|--------|--------||--------|--------|
    !    |   H    |   H    ||    H   |    H   |   H    |   H    ||   H    |   H    |
    !  3 |--------|--------||--------|--------|--------|--------||--------|--------|
    !    |   H    |   H    ||    H   |    H   |   H    |   H    ||   H    |   H    |
    !  2 |========|========||========|========|========|========||========|========|
    !    |   H    |   H    ||   I3   |   I3   |   I4   |   I4   ||   H    |   H    |
    !  1 |--------|--------||--------|--------|--------|--------||--------|--------|
    !    |   H    |   H    ||   I3   |   I3   |   I4   |   I4   ||   H    |   H    |
    !  0 |--------|--------||--------|--------|--------|--------||--------|--------|
    !    |   H    |   H    ||   I1   |   I1   |   I2   |   I2   ||   H    |   H    |
    ! -1 |--------|--------||--------|--------|--------|--------||--------|--------|
    !    |   H    |   H    ||   I1   |   I1   |   I2   |   I2   ||   H    |   H    |
    ! -2 |========|========||========|========|========|========||========|========|
    !    |   H    |   H    ||    H   |    H   |   H    |   H    ||   H    |   H    |
    ! -3 |--------|--------||--------|--------|--------|--------||--------|--------|
    !    |   H    |   H    ||    H   |    H   |   H    |   H    ||   H    |   H    |
    ! -4 |--------|--------||--------|--------|--------|--------||--------|--------|
    !   -8       -6        -4       -2        0        2        4         6        8
    !                                         X

contains

  subroutine startup()
    input = input_t(spatial_reconstruction='minmod', flux_solver='AUSMPW+')
    input%initial_condition_file = 'simple.h5'

    grid => new_2d_grid_block(input)
  end subroutine startup

  subroutine cleanup()
    deallocate(grid)
  end subroutine cleanup

  subroutine test_xy_coords()

    integer(ik) :: i, j

    sync all
    if(this_image() == 1) print*, new_line('') // "Running test_xy_coords" // new_line('')

    associate(left => grid%lbounds(1), right => grid%ubounds(1), &
              bottom => grid%lbounds(2), top => grid%ubounds(2), &
              left_ghost => grid%lbounds_halo(1), right_ghost => grid%ubounds_halo(1), &
              bottom_ghost => grid%lbounds_halo(2), top_ghost => grid%ubounds_halo(2))

      select case(num_images())
      case(4)
        select case(this_image())
        case(1)
          do j = lbound(grid%node_x, dim=2), ubound(grid%node_x, dim=2)
            call assert_equal([-8.0_rk, -6.0_rk, -4.0_rk, -2.0_rk, 0.0_rk, 2.0_rk, 4.0_rk], &
                            grid%node_x(:, j), file=__FILE__, line=__LINE__)
          end do

          call assert_equal(2.0_rk, grid%node_y(:, 5), file=__FILE__, line=__LINE__)
          call assert_equal(1.0_rk, grid%node_y(:, 4), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, grid%node_y(:, 3), file=__FILE__, line=__LINE__)
          call assert_equal(-1.0_rk, grid%node_y(:, 2), file=__FILE__, line=__LINE__)
          call assert_equal(-2.0_rk, grid%node_y(:, 1), file=__FILE__, line=__LINE__)
          call assert_equal(-3.0_rk, grid%node_y(:, 0), file=__FILE__, line=__LINE__)
          call assert_equal(-4.0_rk, grid%node_y(:, -1), file=__FILE__, line=__LINE__)
        case(2)
          do j = lbound(grid%node_x, dim=2), ubound(grid%node_x, dim=2)
            call assert_equal([-4.0_rk, -2.0_rk, 0.0_rk, 2.0_rk, 4.0_rk, 6.0_rk, 8.0_rk], &
                            grid%node_x(:, j), file=__FILE__, line=__LINE__)
          end do

          call assert_equal(2.0_rk, grid%node_y(:, 5), file=__FILE__, line=__LINE__)
          call assert_equal(1.0_rk, grid%node_y(:, 4), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, grid%node_y(:, 3), file=__FILE__, line=__LINE__)
          call assert_equal(-1.0_rk, grid%node_y(:, 2), file=__FILE__, line=__LINE__)
          call assert_equal(-2.0_rk, grid%node_y(:, 1), file=__FILE__, line=__LINE__)
          call assert_equal(-3.0_rk, grid%node_y(:, 0), file=__FILE__, line=__LINE__)
          call assert_equal(-4.0_rk, grid%node_y(:, -1), file=__FILE__, line=__LINE__)
        case(3)
          do j = lbound(grid%node_x, dim=2), ubound(grid%node_x, dim=2)
            call assert_equal([-8.0_rk, -6.0_rk, -4.0_rk, -2.0_rk, 0.0_rk, 2.0_rk, 4.0_rk], &
                            grid%node_x(:, j), file=__FILE__, line=__LINE__)
          end do

          call assert_equal(4.0_rk, grid%node_y(:, 7), file=__FILE__, line=__LINE__)
          call assert_equal(3.0_rk, grid%node_y(:, 6), file=__FILE__, line=__LINE__)
          call assert_equal(2.0_rk, grid%node_y(:, 5), file=__FILE__, line=__LINE__)
          call assert_equal(1.0_rk, grid%node_y(:, 4), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, grid%node_y(:, 3), file=__FILE__, line=__LINE__)
          call assert_equal(-1.0_rk, grid%node_y(:, 2), file=__FILE__, line=__LINE__)
          call assert_equal(-2.0_rk, grid%node_y(:, 1), file=__FILE__, line=__LINE__)
        case(4)
          do j = lbound(grid%node_x, dim=2), ubound(grid%node_x, dim=2)
            call assert_equal([-4.0_rk, -2.0_rk, 0.0_rk, 2.0_rk, 4.0_rk, 6.0_rk, 8.0_rk], &
                            grid%node_x(:, j), file=__FILE__, line=__LINE__)
          end do

          call assert_equal(4.0_rk, grid%node_y(:, 7), file=__FILE__, line=__LINE__)
          call assert_equal(3.0_rk, grid%node_y(:, 6), file=__FILE__, line=__LINE__)
          call assert_equal(2.0_rk, grid%node_y(:, 5), file=__FILE__, line=__LINE__)
          call assert_equal(1.0_rk, grid%node_y(:, 4), file=__FILE__, line=__LINE__)
          call assert_equal(0.0_rk, grid%node_y(:, 3), file=__FILE__, line=__LINE__)
          call assert_equal(-1.0_rk, grid%node_y(:, 2), file=__FILE__, line=__LINE__)
          call assert_equal(-2.0_rk, grid%node_y(:, 1), file=__FILE__, line=__LINE__)
        end select
      case(1)
        ! Single-image testing is a bit easier
        do j = lbound(grid%node_x, dim=2), ubound(grid%node_x, dim=2)
          call assert_equal([-8.0_rk, -6.0_rk, -4.0_rk, -2.0_rk, 0.0_rk, 2.0_rk, 4.0_rk, 6.0_rk, 8.0_rk], &
                            grid%node_x(:, j), file=__FILE__, line=__LINE__)
        end do

        call assert_equal(4.0_rk, grid%node_y(:, 7), file=__FILE__, line=__LINE__)
        call assert_equal(3.0_rk, grid%node_y(:, 6), file=__FILE__, line=__LINE__)
        call assert_equal(2.0_rk, grid%node_y(:, 5), file=__FILE__, line=__LINE__)
        call assert_equal(1.0_rk, grid%node_y(:, 4), file=__FILE__, line=__LINE__)
        call assert_equal(0.0_rk, grid%node_y(:, 3), file=__FILE__, line=__LINE__)
        call assert_equal(-1.0_rk, grid%node_y(:, 2), file=__FILE__, line=__LINE__)
        call assert_equal(-2.0_rk, grid%node_y(:, 1), file=__FILE__, line=__LINE__)
        call assert_equal(-3.0_rk, grid%node_y(:, 0), file=__FILE__, line=__LINE__)
        call assert_equal(-4.0_rk, grid%node_y(:, -1), file=__FILE__, line=__LINE__)

      case default
        error stop "Test failed. This is only configured to test for 1 or 4 images"
      end select

    end associate
    
    ! edge_norm_vectors((x,y), edge, i, j)
    call assert_equal(0.0_rk,  grid%edge_norm_vectors(1,bottom,:,:), file=__FILE__, line=__LINE__)
    call assert_equal(-1.0_rk,  grid%edge_norm_vectors(2,bottom,:,:), file=__FILE__, line=__LINE__)

    call assert_equal(1.0_rk,  grid%edge_norm_vectors(1,right,:,:), file=__FILE__, line=__LINE__)
    call assert_equal(0.0_rk,  grid%edge_norm_vectors(2,right,:,:), file=__FILE__, line=__LINE__)

    call assert_equal(0.0_rk,  grid%edge_norm_vectors(1,top,:,:), file=__FILE__, line=__LINE__)
    call assert_equal(1.0_rk,  grid%edge_norm_vectors(2,top,:,:), file=__FILE__, line=__LINE__)

    call assert_equal(-1.0_rk,  grid%edge_norm_vectors(1,left,:,:), file=__FILE__, line=__LINE__)
    call assert_equal(0.0_rk,  grid%edge_norm_vectors(2,left,:,:), file=__FILE__, line=__LINE__)

    call assert_equal(2.0_rk, grid%volume, file=__FILE__, line=__LINE__)

    call assert_equal(2.0_rk, grid%edge_lengths(1,:,:), file=__FILE__, line=__LINE__)
    call assert_equal(1.0_rk, grid%edge_lengths(2,:,:), file=__FILE__, line=__LINE__)
    call assert_equal(2.0_rk, grid%edge_lengths(3,:,:), file=__FILE__, line=__LINE__)
    call assert_equal(1.0_rk, grid%edge_lengths(4,:,:), file=__FILE__, line=__LINE__)

    write(*, '(a, i0, a)') "Image: ", this_image(), " Success" // " [test_xy_coords]"
  end subroutine test_xy_coords

  ! subroutine test_corner_edge_vector_query()

  !   class(grid_block_t), pointer :: grid
  !   type(input_t) :: input

  !   real(rk), dimension(2, 0:4) :: actual_corner_vectors
  !   real(rk), dimension(2, 0:4) :: grid_corner_vectors

  !   print *, 'Testing test_corner_edge_vector_query()'
  !   input = input_t(spatial_reconstruction='MUSCL', &
  !                   ni_nodes=ni_nodes, nj_nodes=nj_nodes, &
  !                   read_init_cond_from_file=.false., &
  !                   xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

  !   grid => grid_factory(input)

  !   ! grid_corner_vectors = grid%get_corner_vectors(cell_ij=[1, 1], corner='lower-left')
  !   grid_corner_vectors = grid%corner_edge_vectors(:, :, 1, 1)
  !   actual_corner_vectors(:, 0) = [-2.0_rk, -2.0_rk]
  !   actual_corner_vectors(:, 1) = [-2.0_rk, -4.0_rk]
  !   actual_corner_vectors(:, 2) = [-1.0_rk, -2.0_rk]
  !   actual_corner_vectors(:, 3) = [-2.0_rk, 0.0_rk]
  !   actual_corner_vectors(:, 4) = [-3.0_rk, -2.0_rk]

  !   @assertEqual(actual_corner_vectors, grid_corner_vectors)
  !   deallocate(grid)
  ! end subroutine test_corner_edge_vector_query

  ! subroutine test_leftright_midpoint_edge_vector_query()

  !   class(grid_block_t), pointer :: grid
  !   type(input_t) :: input

  !   real(rk), dimension(2, 0:2) :: actual_midpoint_vectors !< ((x,y), (tail,head), (vector1:vector2))
  !   real(rk), dimension(2, 0:2) :: grid_midpoint_vectors !< ((x,y), (tail,head), (vector1:vector2))

  !   print *, 'Testing test_leftright_midpoint_edge_vector_query()'
  !   input = input_t(spatial_reconstruction='MUSCL', &
  !                   ni_nodes=ni_nodes, nj_nodes=nj_nodes, &
  !                   read_init_cond_from_file=.false., &
  !                   xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

  !   grid => grid_factory(input)

  !   ! grid_midpoint_vectors = grid%get_midpoint_vectors(cell_ij=[1, 1], edge='bottom')
  !   grid_midpoint_vectors = grid%leftright_midpoint_edge_vectors(:, :, 1, 1)

  !   actual_midpoint_vectors(:, 0) = [-1.5_rk, -2.0_rk] ! common origin
  !   actual_midpoint_vectors(:, 1) = [-2.0_rk, -2.0_rk] ! vector 1
  !   actual_midpoint_vectors(:, 2) = [-1.0_rk, -2.0_rk]  ! vector 2

  !   @assertEqual(actual_midpoint_vectors, grid_midpoint_vectors)
  !   deallocate(grid)
  ! end subroutine test_leftright_midpoint_edge_vector_query

  ! subroutine test_downup_midpoint_edge_vector_query()

  !   class(grid_block_t), pointer :: grid
  !   type(input_t) :: input

  !   real(rk), dimension(2, 0:2) :: actual_midpoint_vectors !< ((x,y), (tail,head), (vector1:vector2))
  !   real(rk), dimension(2, 0:2) :: grid_midpoint_vectors !< ((x,y), (tail,head), (vector1:vector2))

  !   print *, 'Testing test_downup_midpoint_edge_vector_query()'
  !   input = input_t(spatial_reconstruction='MUSCL', &
  !                   ni_nodes=ni_nodes, nj_nodes=nj_nodes, &
  !                   read_init_cond_from_file=.false., &
  !                   xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

  !   grid => grid_factory(input)

  !   ! grid_midpoint_vectors = grid%get_midpoint_vectors(cell_ij=[1, 1], edge='left')
  !   grid_midpoint_vectors = grid%downup_midpoint_edge_vectors(:, :, 1, 1)

  !   actual_midpoint_vectors(:, 0) = [-2.0_rk, -1.0_rk] ! common origin
  !   actual_midpoint_vectors(:, 1) = [-2.0_rk, -2.0_rk] ! vector 1
  !   actual_midpoint_vectors(:, 2) = [-2.0_rk, 0.0_rk]  ! vector 2
  !   @assertEqual(actual_midpoint_vectors, grid_midpoint_vectors)
  !   deallocate(grid)
  ! end subroutine test_downup_midpoint_edge_vector_query

end module mod_test_grid
