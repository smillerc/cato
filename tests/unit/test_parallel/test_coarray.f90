module partition_mod
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_parallel

  implicit none

contains

  subroutine test_tiling()
    integer(ik) :: neighbors(4)
    neighbors = tile_neighbors_2d(periodic=.false.)
    if(this_image() == 1) write(*, '(a, i3)') 'Number of Images: ', num_images()

    write(*, '(a,i3, a, 4(i3, 1x))') 'image: ', this_image(), ' neighbors: ', neighbors
    sync all
  end subroutine

end module

program test_coarray_partitioning
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use partition_mod
  implicit none

  if(this_image() == 1) print *, "Testing coarray partitioning..."

  call test_tiling()

  if(this_image() == 1) print *, "Success"

end program test_coarray_partitioning
