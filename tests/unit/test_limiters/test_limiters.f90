module test_mod
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_slope_limiter, only: slope_limiter_t
  use mod_flux_limiter, only: flux_limiter_t
  use caf_testing, only: assert_equal
  implicit none

contains

subroutine test_basic_flux_limiters()

    type(flux_limiter_t) :: limiter

    limiter = flux_limiter_t(name='minmod')
    call assert_equal(0.0_rk, limiter%limit(-1.0_rk))
    call assert_equal(0.0_rk, limiter%limit(0.0_rk))
    call assert_equal(1.0_rk, limiter%limit(1.0_rk))
    call assert_equal(1.0_rk, limiter%limit(4.0_rk))
    ! call assert_equal('minmod', limiter%name)

    limiter = flux_limiter_t(name='van_leer')
    call assert_equal(0.0_rk, limiter%limit(-1.0_rk))
    call assert_equal(0.0_rk, limiter%limit(0.0_rk))
    call assert_equal(1.0_rk, limiter%limit(1.0_rk))
    ! call assert_equal('van_leer', limiter%name)

    limiter = flux_limiter_t(name='superbee')
    call assert_equal(0.0_rk, limiter%limit(-1.0_rk))
    call assert_equal(0.0_rk, limiter%limit(0.0_rk))
    call assert_equal(1.0_rk, limiter%limit(0.5_rk))
    call assert_equal(1.0_rk, limiter%limit(1.0_rk))
    call assert_equal(2.0_rk, limiter%limit(2.0_rk))
    call assert_equal(2.0_rk, limiter%limit(4.0_rk))
    ! call assert_equal('superbee', limiter%name)
  end subroutine test_basic_flux_limiters

end module test_mod

program test_limiters
  use test_mod
  implicit none(type, external)

  if(this_image() == 1) print*, new_line('') // "Running test_basic_flux_limiters" // new_line('') 
  call test_basic_flux_limiters()
end program test_limiters
