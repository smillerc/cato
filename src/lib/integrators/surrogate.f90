module mod_surrogate
  implicit none
  private
  public :: surrogate

  type, abstract :: surrogate
    !< The surrogate type only serves to stand in for types used in other places.
    !< This is particularily useful to avoid circular references
  end type
end module
