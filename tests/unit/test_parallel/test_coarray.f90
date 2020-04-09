module partition_mod
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_parallel

  implicit none

contains

  subroutine test_tiling()
    integer(ik) :: neighbors(4)

    if(this_image() == 1) print*, "Calling tile_neighbors_2d(periodic=.false.)"
    if(this_image() == 1) write(*, '(a, i3)') 'Number of Images: ', num_images()

    neighbors = tile_neighbors_2d(periodic=.false.)
    sync all
    write(*, '(a,i3, a, 4(i3, 1x))') 'image: ', this_image(), ' neighbors: ', neighbors
    sync all
    
    if(this_image() == 1) print*, "Calling tile_neighbors_2d(periodic=.true.)"

    neighbors = tile_neighbors_2d(periodic=.true.)
    write(*, '(a,i3, a, 4(i3, 1x))') 'image: ', this_image(), ' neighbors: ', neighbors
  end subroutine

  subroutine test_local_global_indexing()

    integer(ik) :: indices(4)
    integer(ik) :: neighbors(4)
    integer(ik) :: lower_bounds(2), upper_bounds(2)
    integer(ik) :: ilo, ihi, jlo, jhi

    logical :: img_ilo_bc, img_ihi_bc, img_jlo_bc, img_jhi_bc

    integer(ik) :: ilo_img, ihi_img, jlo_img, jhi_img
    integer(ik) :: dims(2)

    if(this_image() == 1) print*, "Calling test_local_global_indexing()"
    if(this_image() == 1) write(*, '(a, i3)') 'Number of Images: ', num_images()

    img_ilo_bc = .false. 
    img_ihi_bc = .false. 
    img_jlo_bc = .false. 
    img_jhi_bc = .false.

    ! Global indices ranges
    ilo = 0
    ihi = 200
    jlo = 0
    jhi = 100
    dims = [ihi+1, jhi+1]
    indices = tile_indices(dims) - 1
    ! lower_bounds = indices([1, 3])
    ! upper_bounds = indices([2, 4])

    ilo_img = indices(1)
    jlo_img = indices(3)
    ihi_img = indices(2)
    jhi_img = indices(4)

    if (ilo_img == ilo) img_ilo_bc = .true.
    if (ihi_img == ihi) img_ihi_bc = .true.
    if (jlo_img == jlo) img_jlo_bc = .true.
    if (jhi_img == jhi) img_jhi_bc = .true.

    neighbors = tile_neighbors_2d(periodic=.true.)

    ! write(*,'(a, i0, a, 4(i4))') 'image: ', this_image(), ' neighbors: ', neighbors
    write(*,'(a, i0, 4(a, i0))') 'image: ', this_image(), &
    ' i = ', ilo_img, ":", ihi_img, &
    ' j = ', jlo_img, ":", jhi_img


    write(*,'(a, i0, 4(a, 1x, l, 1x))') 'image: ', this_image(), &
    " img_ilo_bc: ", img_ilo_bc,&
    " img_ihi_bc: ", img_ihi_bc,&
    " img_jlo_bc: ", img_jlo_bc,&
    " img_jhi_bc: ", img_jhi_bc
    print*

    
  end subroutine

end module

program test_coarray_partitioning
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use partition_mod
  implicit none

  if(this_image() == 1) print *, "Testing coarray partitioning..."

  call test_tiling()
  sync all
  call test_local_global_indexing()
  sync all
  if(this_image() == 1) print *, "Success"

end program test_coarray_partitioning
