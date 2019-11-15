module mod_vector

  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between

  type vector_t
    real(rk) :: x
    real(rk) :: y
    real(rk) :: length
  contains
  end type

  interface vector_t
    module procedure :: constructor, constructor_from_2d
  end interface

  interface operator(.unitnorm.)
    module procedure unit_normal
  end interface

  interface operator(.cross.)
    module procedure vec_cross_product
  end interface

  interface operator(.dot.)
    module procedure vec_dot_product
  end interface

contains

  type(vector_t) pure function constructor(x, y) result(vec)
    !< Construcgtor from an (x,y) pair
    real(rk), intent(in) :: x
    real(rk), intent(in) :: y

    vec%x = x
    vec%y = y
    vec%length = norm2([vec%x, vec%y])
  end function

  type(vector_t) pure function constructor_from_2d(x, y) result(vec)
    !< Constructor from a set of (x,y) pairs
    real(rk), intent(in), dimension(2) :: x
    real(rk), intent(in), dimension(2) :: y

    vec%x = x(2) - x(1)
    vec%y = y(2) - y(1)
    vec%length = norm2([vec%x, vec%y])
  end function

  type(vector_t) pure function unit_normal(vec) result(norm_vec)
    !< Creat the .unitnorm. operator for the vector_t type
    class(vector_t), intent(in) :: vec
    norm_vec = vector_t(x=vec%x / vec%length, y=vec%y / vec%length)
  end function

  real(rk) pure function vec_dot_product(vec1, vec2) result(scalar_dot_product)
    !< Create the .dot. operator for the vector_t type
    class(vector_t), intent(in) :: vec1, vec2
    scalar_dot_product = dot_product([vec1%x, vec1%y],[[vec2%x, vec2%y]])
  end function

  real(rk) pure function vec_cross_product(vec1, vec2) result(cross_product)
    !< Create the .cross. operator for the vector_t type. Since these vectors are only 2d, then
    !< the cross product is a scalar quantity
    class(vector_t), intent(in) :: vec1, vec2

    cross_product = vec1%x * vec2%y - vec2%x * vec1%x
  end function

  real(rk) pure function angle_between(vec1, vec2) result(angle)
    class(vector_t), intent(in) :: vec1, vec2

    angle = acos(vec1.dot.vec2) / (vec1%length * vec2%length)
  end function

end module mod_vector
