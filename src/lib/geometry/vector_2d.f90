module mod_vector_2d
  use iso_fortran_env, only: int32, real64

  implicit none

  type vector_2d_t
    real(real64), dimension(2) :: x
    real(real64), dimension(2) :: y
    real(real64), dimension(2) :: length
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
    real(real64), intent(in), dimension(2) :: x
    real(real64), intent(in), dimension(2) :: y

    vec
  end function

end module mod_vector_2d

