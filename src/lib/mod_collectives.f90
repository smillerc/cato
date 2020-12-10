module collectives
  use iso_fortran_env
  implicit none

contains

  real(real64) function sum_to_root(x, root) result(x_sum)
    !< Broadcast the sum of x to all images using a binary-tree algorithm.
    real(real64), intent(in) :: x                !< data to perform the operation on
    integer(int32), intent(in), optional :: root !< desired image to send result to, defaults to 1
    real(real64), save :: y[*]                   !< (image); coarray buffer to perform operation on
    integer(int32) :: img_idx                    !< image loop index
    logical :: is_valid_image                    !< ensure only valid images (w/in the tree) do the operation

    x_sum = 0.0_real64

    ! Allocate the coarray buffer and copy over the data
    !if (.not. allocated(y)) allocate(y[*])
    y = x
    sync all

    img_idx = 1
    do while(img_idx < num_images())
      is_valid_image = this_image() + img_idx <= num_images() .and. &
                       mod(this_image() - 1, 2 * img_idx) == 0

      ! Sum up the values if we're in a valid image index
      if(is_valid_image) then
        y = y + y[this_image() + img_idx]
      endif

      img_idx = 2 * img_idx
      sync all
    enddo

    ! Now broadcast to the image of choice
    if(present(root)) then
      if(this_image() == root) x_sum = y[1]
    else
      if(this_image() == 1) x_sum = y
    endif

  endfunction sum_to_root

  real(real64) function max_to_root(x, root) result(x_max)
    !< Broadcast the maximum value of x to all images
    real(real64), intent(in) :: x                !< data to perform the operation on
    integer(int32), intent(in), optional :: root !< desired image to send result to, defaults to 1
    real(real64), save :: y[*]                   !< (image); coarray buffer to perform operation on
    integer(int32) :: img_idx                    !< image loop index
    logical :: is_valid_image                    !< ensure only valid images (w/in the tree) do the operation

    x_max = 0.0_real64

    ! Allocate the coarray buffer and copy over the data
    !if (.not. allocated(y)) allocate(y[*])
    y = x
    sync all

    img_idx = 1
    do while(img_idx < num_images())
      is_valid_image = this_image() + img_idx <= num_images() .and. &
                       mod(this_image() - 1, 2 * img_idx) == 0

      ! Sum up the values if we're in a valid image index
      if(is_valid_image) then
        y = max(y, y[this_image() + img_idx])
      endif

      img_idx = 2 * img_idx
      sync all
    enddo

    ! Now broadcast to the image of choice
    if(present(root)) then
      if(this_image() == root) x_max = y[1]
    else
      if(this_image() == 1) x_max = y
    endif

  endfunction max_to_root

  real(real64) function root_to_all(x, root) result(s)
    !< Broadcast the value of x to all images
    real(real64), intent(in) :: x                !< data to perform the operation on
    integer(int32), intent(in), optional :: root !< desired image to send result to, defaults to 1
    real(real64), save :: y[*]                   !< (image); coarray buffer to perform operation on
    integer(int32) :: img_idx                    !< image loop index
    logical :: is_valid_image                    !< ensure only valid images (w/in the tree) do the operation

    !if (.not. allocated(y)) allocate(y[*])

    if(present(root)) then
      if(this_image() == root) y[1] = x
    else
      if(this_image() == 1) y = x
    endif
    sync all

    img_idx = 1
    do while(img_idx < num_images())
      img_idx = 2 * img_idx
    enddo

    do while(img_idx > 0)
      is_valid_image = this_image() + img_idx <= num_images() .and. &
                       mod(this_image() - 1, 2 * img_idx) == 0

      if(is_valid_image) y[this_image() + img_idx] = y
      img_idx = img_idx / 2
      sync all
    enddo
    s = y
  endfunction root_to_all

  real(real64) function sum_to_all(x) result(x_sum)
    !< Broadcast the global sum to all images
    real(real64), intent(in) :: x !< data to perform the operation on
    x_sum = root_to_all(sum_to_root(x))
  endfunction sum_to_all

  real(real64) function max_to_all(x) result(x_max)
    !< Broadcast the global maximum to all images
    real(real64), intent(in) :: x !< data to perform the operation on
    x_max = root_to_all(max_to_root(x))
  endfunction max_to_all

  real(real64) function min_to_all(x) result(x_min)
    !< Broadcast the global minimum to all images
    real(real64), intent(in) :: x !< data to perform the operation on
    x_min = -max_to_all(-x)
  endfunction min_to_all
endmodule collectives
