module mod_vector_2d
  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  type vector_2d_t
    real(rk), dimension(2) :: x
    real(rk), dimension(2) :: y
    real(rk), dimension(2) :: length
  contains
    procedure, public :: initialize
    procedure, public ::
  end type

  interface operator(.unitnorm.)
    module procedure unit_normal
  end interface

  interface operator(.cross.)
    module procedure cross_product
  end interface

  interface operator(.dot.)
    module procedure d_product
  end interface

contains

  pure type(vector_2d_t) function vector(x, y) result(vec)
    real(rk), intent(in), dimension(2) :: x
    real(rk), intent(in), dimension(2) :: y

    vec
  end function

end module mod_vector_2d

