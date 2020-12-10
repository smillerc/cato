! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_nondimensionalization
  !< Summary: Provide scaling factors to non-dimensionalize the Euler equations
  !< Date: 04/25/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !      [1]

  use iso_fortran_env, only: rk => real64

  implicit none

  !TODO: Scale so that the cell size is always 1!
  real(rk), protected :: t_0 = 1.0_rk   !< scale factor for time
  real(rk), protected :: rho_0 = 1.0_rk !< scale factor for density
  real(rk), protected :: l_0 = 1.0_rk   !< scale factor for length
  real(rk), protected :: v_0 = 1.0_rk   !< scale factor for velocity
  real(rk), protected :: p_0 = 1.0_rk   !< scale factor for pressure
  real(rk), protected :: e_0 = 1.0_rk   !< scale factor for energy (same as pressure)

  logical, parameter :: auto_scale_by_grid = .true.
  logical, protected :: length_scale_set = .false.
  logical, protected :: scale_factors_set = .false.

contains

  subroutine set_length_scale(length_scale)
    real(rk), intent(in) :: length_scale

    if(length_scale < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_length_scale(), length_scale < 0"
    endif

    l_0 = length_scale
    length_scale_set = .true.
  endsubroutine

  subroutine set_scale_factors(density_scale, pressure_scale)
    !< Set the non-dimensional scale factors based on the provided scales. This
    !< allows the user to set these protected factors outside of the module.

    real(rk), intent(in) :: density_scale
    real(rk), intent(in) :: pressure_scale

    if(.not. length_scale_set) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), "// &
        "the length scale needs to be set first (via the grid)"
    endif

    if(density_scale < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), density_scale < 0"
    endif

    if(pressure_scale < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), pressure_scale < 0"
    endif

    rho_0 = density_scale
    p_0 = pressure_scale
    v_0 = sqrt(p_0 / rho_0)
    e_0 = p_0
    t_0 = l_0 / v_0

    if(this_image() == 1) then
      print *
      write(*, '(a)') "Scale factors (CATO is non-dimensional)"
      write(*, '(a)') "========================================"
      write(*, '(a, es10.3)') "Time scale factor (t_0):      ", t_0
      write(*, '(a, es10.3)') "Density scale factor (rho_0): ", rho_0
      write(*, '(a, es10.3)') "Length scale factor (l_0):    ", l_0
      write(*, '(a, es10.3)') "Pressure scale factor (p_0):  ", p_0
      write(*, '(a, es10.3)') "Velocity scale factor (v_0):  ", v_0
      write(*, '(a, es10.3)') "Energy scale factor (e_0):    ", e_0
      write(*, '(a)') "========================================"
      print *
    endif

    scale_factors_set = .true.

  endsubroutine set_scale_factors

  ! subroutine set_scale_factors(time_scale, length_scale, density_scale)
  !   !< Set the non-dimensional scale factors based on the provided scales. This
  !   !< allows the user to set these protected factors outside of the module.

  !   real(rk), intent(in) :: time_scale
  !   real(rk), intent(in) :: length_scale
  !   real(rk), intent(in) :: density_scale

  !   if(time_scale < 0.0_rk) then
  !     error stop "Error in mod_nondimensionalization::set_scale_factors(), time_scale < 0"
  !   end if

  !   if(length_scale < 0.0_rk) then
  !     error stop "Error in mod_nondimensionalization::set_scale_factors(), length_scale < 0"
  !   end if

  !   if(density_scale < 0.0_rk) then
  !     error stop "Error in mod_nondimensionalization::set_scale_factors(), density_scale < 0"
  !   end if

  !   t_0 = time_scale
  !   rho_0 = density_scale
  !   l_0 = length_scale
  !   v_0 = l_0 / t_0
  !   p_0 = rho_0 * v_0**2
  !   e_0 = p_0

  !   print *
  !   write(*, '(a)') "Scale factors (CATO is non-dimensional)"
  !   write(*, '(a)') "========================================"
  !   write(*, '(a, es10.3)') "Time scale factor (t_0):      ", t_0
  !   write(*, '(a, es10.3)') "Density scale factor (rho_0): ", rho_0
  !   write(*, '(a, es10.3)') "Length scale factor (l_0):    ", l_0
  !   write(*, '(a, es10.3)') "Pressure scale factor (p_0):  ", p_0
  !   write(*, '(a, es10.3)') "Velocity scale factor (v_0):  ", v_0
  !   write(*, '(a, es10.3)') "Energy scale factor (e_0):    ", e_0
  !   write(*, '(a)') "========================================"
  !   print *

  ! end subroutine set_scale_factors
endmodule mod_nondimensionalization
