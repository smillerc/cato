module mod_conserved_vars
  use iso_fortran_env, only: int32, real64

  implicit none

  private
  public :: conserved_vars_t

  type :: conserved_vars_t
    real(real64), allocatable, dimension(:, :) :: density
    real(real64), allocatable, dimension(:, :) :: x_velocity
    real(real64), allocatable, dimension(:, :) :: y_velocity
    real(real64), allocatable, dimension(:, :) :: pressure
  end type

contains

end module mod_conserved_vars
