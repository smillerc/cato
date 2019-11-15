module mod_flux_tensor
  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: flux_tensor_t

  type flux_tensor_t
    real(rk), dimension(2, 4) :: state
  end type

  interface operator(.dot.)
    module procedure flux_dot_product
  end interface

contains

  type(flux_tensor_t) pure function constructor(primitive_variables) result(H)
    real(rk), dimension(4) :: primitive_variables
  end function

  pure function flux_dot_product result(H)
    real(rk), dimension(4) :: H
  end function

end module mod_flux_tensor
