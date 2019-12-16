# Contributing Guide

## Design

1. Design to the interface
2. The code should comment itself as much as possible
3. Make the code look like the math as much as possible

## Coding Conventions

Refering to 2 and 3 above; Do things like the following

```fortran
associate(rho => self%density, &
          p => self%pressure, &
          gamma => self%polytropic_index, &
          cs => self%speed_of_sound)
  cs = sqrt(gamma*P/rho)
end associate
```


#### Naming Conventions
Use `snake_case` naming exclusively. Fortran is a case insensitive language, so the compiler does not differentiate between `camelCase` or `PascalCase` variables.




#### Modules

All modules should have the prefix `mod_`, as in `mod_grid`. I prefer to keep one module per file as a general rule of thumb, but this may not always be the case.

Modules should have the following structure

```fortran
module mod_new_type
  !< Summary: What this module does

  use, intrinsic iso_fortran_env, only: ik => int32, rk => real64
  use mod_parent_type, only : parent_type_t
  use mod_input, only : input_t

  implicit none

  private
  public :: new_type_t

  type, extends(parent_type_t) :: new_type_t
    !< Class to do stuff
    private
  contains
    procedure, public :: initialize
    final :: finalize
  end type new_type_t

  interface new_type_t
    module procedure :: constructor
  end interface

contains

  type(new_type_t) pure function constructor(input)
    type(input_t), intent(in) :: input
  end function

  pure subroutine finalize(self)
    type(new_type_t), intent(inout) :: self
  end function

end mod_new_type
```

### Derived Types (aka Objects/Classes)
