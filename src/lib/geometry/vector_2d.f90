module mod_vector_2d
  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: vector_2d_t, operator(.unitnorm.)

  type vector_2d_t
    real(rk), dimension(2) :: x
    real(rk), dimension(2) :: y
    real(rk) :: length
  contains
    procedure, public :: initialize
    ! procedure, private :: unit_normal
  end type

  ! constructor
  ! interface vector_2d_t
  !   module procedure make_vector
  ! end interface

  interface operator(.unitnorm.)
    module procedure unit_normal
  end interface

  ! interface operator(.cross.)
  !   module procedure cross_product
  ! end interface

  ! interface operator(.dot.)
  !   module procedure d_product
  ! end interface

contains

  subroutine initialize(self, x_coords, y_coords)
    class(vector_2d_t), intent(inout) :: self
    !< Constructor for vector_2d_t
    real(rk), intent(in), dimension(2) :: x_coords
    real(rk), intent(in), dimension(2) :: y_coords

    self%x = x_coords
    self%y = y_coords
    self%length = sqrt((self%x(2) - self%x(1))**2 + (self%y(2) - self%y(1))**2)
  end subroutine

  ! pure type(vector_2d_t) function make_vector(x_coords, y_coords) result(vec)
  !   !< Constructor for vector_2d_t
  !   real(rk), intent(in), dimension(2) :: x_coords
  !   real(rk), intent(in), dimension(2) :: y_coords

  !   vec%x = x_coords
  !   vec%y = y_coords
  !   vec%length = sqrt((vec%x(2) - vec%x(1))**2 + (vec%y(2) - vec%y(1))**2)

  ! end function

  type(vector_2d_t) function unit_normal(vec) result(norm_vec)
    !< Normalize a vector_2d_t to unit length
    class(vector_2d_t), intent(in) :: vec

    real(rk), dimension(2) :: vec_at_origin

    vec_at_origin(:) = [vec%x(2) - vec%x(1), vec%y(2) - vec%y(1)]

    vec_at_origin = vec_at_origin / norm2(vec_at_origin)
    call norm_vec%initialize(x_coords=[vec%x(1), vec%x(1) + vec_at_origin(1)], y_coords=[vec%y(1), vec%y(1) + vec_at_origin(2)])

  end function

end module mod_vector_2d
