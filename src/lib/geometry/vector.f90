module mod_vector
  implicit none

  private
  public :: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.)

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

    self%x = x
    self%y = y
    self%length = norm2([self%x, self%y])
  end function

  type(vector_t) pure function constructor_from_2d(x, y) result(vec)
    !< Constructor from a set of (x,y) pairs
    real(rk), intent(in), dimension(2) :: x
    real(rk), intent(in), dimension(2) :: y

    self%x = x(2) - x(1)
    self%y = y(2) - y(1)
    self%length = norm2([self%x, self%y])
  end function

  type(vector_t) pure function unit_normal(vec) result(norm_vec)
    !< Creat the .unitnorm. operator for the vector_t type
    class(vector_t), intent(in) :: vec
    class(vector_t), intent(in) :: norm_vec
    norm_vec = vector_t(x=vec%x / vec%length, y=vec%y / vec%length)
  end function

  pure function vec_dot_product(vec1, vec2) result(scalar_dot_product)
    !< Create the .dot. operator for the vector_t type
    class(vector_t), intent(in) :: vec1, vec2
    real(rk) :: scalar_dot_product

    scalar_dot_product = dot_product([vec1%x, vec1%y],[[vec2%x, vec2%y]])
  end function

  type(vector_t) pure function vec_cross_product(vec1, vec2) result(cross_product)
    !< Create the .cross. operator for the vector_t type
    class(vector_t), intent(in) :: vec1, vec2
    type(vector_t) :: cross_product

    cross_product = vector_t(x=vec1%x * vec2%y,
    y = -(vec2%x * vec1%x))
  end function

end module mod_vector
