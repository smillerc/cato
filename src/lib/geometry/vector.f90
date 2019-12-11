module mod_vector

  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between

  type vector_t
    real(rk) :: x = 0.0_rk
    real(rk) :: y = 0.0_rk
    real(rk) :: length = 0.0_rk
  contains
    private
    procedure, pass :: write => write_vector
    generic, public :: write(formatted) => write
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
    ! //TODO: fix this so no temp array copy warning issued
    !< Constructor from a set of (x,y) pairs
    real(rk), intent(in), dimension(2) :: x !< (x1,x2)
    real(rk), intent(in), dimension(2) :: y !< (y1,y2)

    vec%x = x(2) - x(1)
    vec%y = y(2) - y(1)
    vec%length = norm2([vec%x, vec%y])
  end function

  subroutine write_vector(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of `write(*,*) vector_t`

    class(vector_t), intent(in) :: self  !< vector class
    integer, intent(in) :: unit  !< input/output unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write(unit, '(3(a,f0.4))', iostat=iostat) 'Vector(x,y): (', self%x, ', ', self%y, ') Length: ', self%length
  end subroutine

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

    cross_product = vec1%x * vec2%y - vec1%y * vec2%x
  end function

  real(rk) pure function angle_between(vec1, vec2) result(angle)
    class(vector_t), intent(in) :: vec1, vec2

    angle = acos(vec1.dot.vec2) / (vec1%length * vec2%length)
  end function

end module mod_vector
