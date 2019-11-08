module mod_conserved_vars
  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: conserved_vars_t

  type :: conserved_vars_t
    real(rk), allocatable, dimension(:, :) :: density
    real(rk), allocatable, dimension(:, :) :: x_velocity
    real(rk), allocatable, dimension(:, :) :: y_velocity
    real(rk), allocatable, dimension(:, :) :: pressure
  end type

contains

end module mod_conserved_vars
