module mod_surrogate
  use mod_hermetic, only: hermetic
  implicit none
  private
  public :: surrogate

  type, abstract, extends(hermetic) :: surrogate
    !< The surrogate type only serves to stand in for types used in other places.
    !< This is particularily useful to avoid circular references
  end type
end module
