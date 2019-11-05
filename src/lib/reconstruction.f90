module reconstruction
  use iso_fortran_env, only: int32, real64

  use slope_limiter, only: limit

  implicit none

  ! type, abstract :: reconstruction_base_t
  ! end type

  ! type, extends(reconstruction_base_t) :: piecewise_linear_reconstruction_t
  ! end type

contains

  ! pure function reconstruct_piecewise() result(W)
  !   !< Follows the reconstruction scheme of DOI: 10.1016/j.jcp.2009.04.001, although it
  !   !< inherits it from DOI: 10.1016/j.jcp.2006.03.018, where the explaination is a bit clearer,
  !   !< and without typos
  ! end function

end module reconstruction
