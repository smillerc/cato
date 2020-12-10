program test_grid
  use mod_test_grid
  implicit none

  call startup() 
  call test_xy_coords()
  call cleanup()
end program test_grid
