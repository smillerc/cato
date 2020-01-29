module mod_source
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  type, abstract :: source_t
    character(len=:), allocatable :: source_type
    integer(ik) :: ilo = 0 !< Index to apply source term at
    integer(ik) :: jlo = 0 !< Index to apply source term at
    integer(ik) :: ihi = 0 !< Index to apply source term at
    integer(ik) :: jhi = 0 !< Index to apply source term at
  contains
    procedure(apply_source), deferred :: apply_source
    procedure(copy_source), public, deferred :: copy
    generic :: assignment(=) => copy
  end type

  abstract interface
    subroutine copy_source(out_source, in_source)
      import :: source_t
      class(source_t), intent(in) :: in_source
      class(source_t), intent(inout) :: out_source
    end subroutine copy_source

    subroutine apply_source(self, primitive_vars, lbounds, time)
      import :: source_t, rk, ik
      class(source_t), intent(inout) :: self
      integer(ik), dimension(3), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: primitive_vars
      real(rk), intent(in) :: time
    end subroutine apply_source
  end interface
contains

end module mod_source
