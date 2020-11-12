program test_bc
  use mod_test_bc
  implicit none

  call startup() 
  call test_symmetry()
  call test_x_periodic()
  call test_y_periodic()
  call test_periodic_all()

  call cleanup()
end program test_bc
