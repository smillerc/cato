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
  real(rk), protected :: t_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: t_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: t_to_dim = 1.0_rk     !< scale factor for time

  real(rk), protected :: density_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: density_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: density_to_dim = 1.0_rk     !< scale factor for time

  real(rk), protected :: vel_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: vel_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: vel_to_dim = 1.0_rk     !< scale factor for time

  real(rk), protected :: len_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: len_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: len_to_dim = 1.0_rk     !< scale factor for time

  real(rk), protected :: press_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: press_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: press_to_dim = 1.0_rk     !< scale factor for time

  real(rk), protected :: energy_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: energy_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: energy_to_dim = 1.0_rk     !< scale factor for time

  real(rk), protected :: pow_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: pow_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: pow_to_dim = 1.0_rk     !< scale factor for time

  real(rk), protected :: force_ref = 1.0_rk        !< reference time in sec
  real(rk), protected :: force_to_nondim = 1.0_rk  !< factor to convert time from dimensional to non-dimensional
  real(rk), protected :: force_to_dim = 1.0_rk     !< scale factor for time

  ! real(rk), protected :: t_0 = 1.0_rk   !< scale factor for length
  ! real(rk), protected :: rho_0 = 1.0_rk   !< scale factor for length
  ! real(rk), protected :: l_0 = 1.0_rk   !< scale factor for length
  ! real(rk), protected :: v_0 = 1.0_rk   !< scale factor for velocity
  ! real(rk), protected :: p_0 = 1.0_rk   !< scale factor for pressure
  ! real(rk), protected :: e_0 = 1.0_rk   !< scale factor for energy (same as pressure)

  logical, protected :: apply_nondimensionalization = .true.
  ! logical, parameter :: auto_scale_by_grid = .true.
  ! logical, protected :: length_scale_set = .false.
  logical, protected :: scale_factors_set = .false.

contains

  subroutine set_nondim_flag(apply_dim)
    logical, intent(in) :: apply_dim
    apply_nondimensionalization = apply_dim
  endsubroutine

  ! subroutine set_length_scale(length_scale)
  !   real(rk), intent(in) :: length_scale

  !   if(apply_nondimensionalization) then
  !     if(length_scale < 0.0_rk) then
  !       error stop "Error in mod_nondimensionalization::set_length_scale(), length_scale < 0"
  !     endif

  !     l_0 = length_scale
  !   endif
  !   length_scale_set = .true.

  ! endsubroutine

  subroutine set_refrence_quantities(ref_length, ref_velocity, ref_density)
    real(rk), intent(in) :: ref_length   !< reference
    real(rk), intent(in) :: ref_velocity !< reference
    real(rk), intent(in) :: ref_density  !< reference

    t_ref = ref_length / ref_velocity
    len_ref = ref_length
    density_ref = ref_density
    vel_ref = ref_velocity
    press_ref = density_ref * vel_ref**2
    energy_ref = density_ref * vel_ref**2
    pow_ref = density_ref * vel_ref**2 / t_ref

    ! Factors used to convert from dim to non-dim and vice versa
    ! e.g. t (that the code uses) can be converted to real units by => t * t_to_dim
    !      and to convert from units to non-dim vesion => t * t_to_nondim
    len_to_nondim = 1.0_rk / len_ref
    len_to_dim = len_ref

    t_to_nondim = 1.0_rk / t_ref
    t_to_dim = t_ref

    density_to_nondim = 1.0_rk / density_ref
    density_to_dim = density_ref

    vel_to_nondim = 1.0_rk / vel_ref
    vel_to_dim = vel_ref

    press_to_nondim = 1.0_rk / press_ref
    press_to_dim = press_ref

    energy_to_nondim = 1.0_rk / energy_ref
    energy_to_dim = energy_ref

    pow_to_nondim = 1.0_rk / pow_ref
    pow_to_dim = pow_ref

    if(density_to_nondim < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), density_scale < 0"
    endif

    if(press_to_nondim < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), pressure_scale < 0"
    endif

    if(this_image() == 1) then
      print *
      write(*, '(a)') "Scale factors (CATO is non-dimensional)"
      write(*, '(a)') "============================================"

      write(*, '(a, es10.3)') "Reference length        (input): ", len_ref
      write(*, '(a, es10.3)') "Reference density       (input): ", density_ref
      write(*, '(a, es10.3)') "Reference velocity      (input): ", vel_ref
      write(*, '(a, es10.3)') "Reference time     (calculated): ", t_ref
      write(*, '(a, es10.3)') "Reference energy   (calculated): ", energy_ref
      write(*, '(a, es10.3)') "Reference pressure (calculated): ", press_ref
      print *

      write(*, '(a, es10.3)') "Time scale factor (t_tilde)      : ", t_to_nondim
      write(*, '(a, es10.3)') "Density scale factor (rho_tilde) : ", density_to_nondim
      write(*, '(a, es10.3)') "Length scale factor (l_tilde)    : ", len_to_nondim
      write(*, '(a, es10.3)') "Pressure scale factor (p_tilde)  : ", press_to_nondim
      write(*, '(a, es10.3)') "Velocity scale factor (v_tilde)  : ", vel_to_nondim
      write(*, '(a, es10.3)') "Energy scale factor (e_tilde)    : ", energy_to_nondim
      write(*, '(a)') "============================================"
      print *
    endif

    scale_factors_set = .true.

  endsubroutine

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
