module test_parallel
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_parallel
  use funit

  implicit none

contains

  @test
  subroutine test_tiling()

    @assertEqual([1, 1], num_tiles(1))
    @assertEqual([2, 1], num_tiles(2))
    @assertEqual([3, 1], num_tiles(3))
    @assertEqual([2, 2], num_tiles(4))
    @assertEqual([5, 1], num_tiles(5))
    @assertEqual([3, 2], num_tiles(6))

  end subroutine

end module test_parallel
