! MIT License
! Copyright (c) 2020 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_local_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic

  use mod_error, only: ALL_OK, NEG_DENSITY, NEG_PRESSURE, NANS_FOUND, error_msg
  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data, n_ghost_layers
  use mod_nondimensionalization, only: scale_factors_set, rho_0, v_0, p_0, e_0, t_0
  use mod_floating_point_utils, only: near_zero, nearly_equal, neumaier_sum, neumaier_sum_2, neumaier_sum_3, neumaier_sum_4
  use mod_functional, only: operator(.sort.)
  use mod_units
  ! use mod_distinguisher, only: distinguish
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_flux_solver, only: flux_solver_t
  ! use mod_ausm_plus_solver, only: ausm_plus_solver_t
  ! use mod_fvleg_solver, only: fvleg_solver_t
  ! use mod_m_ausmpw_plus_solver, only: m_ausmpw_plus_solver_t
  use mod_ausmpw_plus_solver, only: ausmpw_plus_solver_t
  ! use mod_slau_solver, only: slau_solver_t

  implicit none

  private
  public :: local_fluid_t, new_fluid

  logical, parameter :: filter_small_mach = .false.

  type :: local_fluid_t
    !< Fluid solver physics package

    private ! make all private by default

    real(rk), dimension(:, :), allocatable, public :: rho    !< (i, j); density (conserved & primitive)
    real(rk), dimension(:, :), allocatable, public :: rho_u  !< (i, j); density * x-velocity (conserved)
    real(rk), dimension(:, :), allocatable, public :: rho_v  !< (i, j); density * y-velocity (conserved)
    real(rk), dimension(:, :), allocatable, public :: rho_E  !< (i, j); density * total energy (conserved)
    real(rk), dimension(:, :), allocatable, public :: u      !< (i, j); x-velocity (primitive)
    real(rk), dimension(:, :), allocatable, public :: v      !< (i, j); y-velocity (primitive)
    real(rk), dimension(:, :), allocatable, public :: p      !< (i, j); pressure (primitive)
    real(rk), dimension(:, :), allocatable, public :: cs     !< (i, j); sound speed
    real(rk), dimension(:, :), allocatable, public :: mach_u !< (i, j); mach number in the x-direction
    real(rk), dimension(:, :), allocatable, public :: mach_v !< (i, j); mach number in the y-direction

    ! Intel compiler alignment hints
    !dir$ attributes align:__ALIGNBYTES__ :: rho
    !dir$ attributes align:__ALIGNBYTES__ :: rho_u
    !dir$ attributes align:__ALIGNBYTES__ :: rho_v
    !dir$ attributes align:__ALIGNBYTES__ :: rho_E
    class(flux_solver_t), allocatable :: solver !< solver scheme used to flux quantities at cell interfaces

    ! Time variables
    character(len=10) :: time_integration_scheme = 'ssp_rk2'
    real(rk) :: time = 0.0_rk !< current simulation time
    real(rk) :: dt = 0.0_rk   !< time step
    integer(ik) :: iteration = 0 !< current iteration number

    logical, public :: prim_vars_updated = .false.
    logical :: smooth_residuals = .true.

    ! Residual history
    character(len=32) :: residual_hist_file = 'residual_hist.csv'
    logical :: residual_hist_header_written = .false.
  contains
    ! Private methods
    private
    procedure :: sanity_check
    procedure, nopass :: add_fields
    procedure, nopass :: mult_fields
    procedure, nopass :: subtract_fields

    ! Operators
    procedure, pass(lhs), public :: add_fluid
    procedure, pass(lhs), public :: subtract_fluid
    procedure, pass(lhs), public :: fluid_mul_real
    procedure, pass(rhs), public :: real_mul_fluid
    procedure, pass(lhs), public :: assign_fluid

    ! Public methods
    procedure, public :: initialize
    procedure, public :: force_finalization

    ! Finalizer
    final :: finalize

    ! Map operators to corresponding procedures
    generic :: operator(+) => add_fluid
    generic :: operator(-) => subtract_fluid
    generic :: operator(*) => real_mul_fluid, fluid_mul_real
    generic :: assignment(=) => assign_fluid
  end type local_fluid_t

contains

  function new_fluid(input, grid) result(fluid)
    !< Fluid constructor
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid
    type(local_fluid_t), pointer :: fluid

    allocate(fluid)
    call fluid%initialize(input, grid)
  end function new_fluid

  subroutine initialize(self, input, grid)
    class(local_fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid
    class(flux_solver_t), pointer :: solver => null()

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi, io
  end subroutine initialize

  subroutine force_finalization(self)
    class(local_fluid_t), intent(inout) :: self

    call debug_print('Running local_fluid_t%force_finalization()', __FILE__, __LINE__)
    if(allocated(self%rho)) deallocate(self%rho)
    if(allocated(self%u)) deallocate(self%u)
    if(allocated(self%v)) deallocate(self%v)
    if(allocated(self%p)) deallocate(self%p)
    if(allocated(self%rho_u)) deallocate(self%rho_u)
    if(allocated(self%rho_v)) deallocate(self%rho_v)
    if(allocated(self%rho_E)) deallocate(self%rho_E)
    if(allocated(self%cs)) deallocate(self%cs)
    if(allocated(self%mach_u)) deallocate(self%mach_u)
    if(allocated(self%mach_v)) deallocate(self%mach_v)
    ! if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine force_finalization

  subroutine finalize(self)
    type(local_fluid_t), intent(inout) :: self

    call debug_print('Running local_fluid_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%rho)) deallocate(self%rho)
    if(allocated(self%u)) deallocate(self%u)
    if(allocated(self%v)) deallocate(self%v)
    if(allocated(self%p)) deallocate(self%p)
    if(allocated(self%rho_u)) deallocate(self%rho_u)
    if(allocated(self%rho_v)) deallocate(self%rho_v)
    if(allocated(self%rho_E)) deallocate(self%rho_E)
    if(allocated(self%cs)) deallocate(self%cs)
    if(allocated(self%mach_u)) deallocate(self%mach_u)
    if(allocated(self%mach_v)) deallocate(self%mach_v)
    if(allocated(self%solver)) deallocate(self%solver)
  end subroutine finalize

  function subtract_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(local_fluid_t), intent(in) :: lhs
    class(local_fluid_t), intent(in) :: rhs
    type(local_fluid_t), allocatable :: difference

    call debug_print('Running local_fluid_t%subtract_fluid()', __FILE__, __LINE__)

    allocate(difference, source=lhs)
    call subtract_fields(a=lhs%rho, b=rhs%rho, c=difference%rho) ! c=a+b
    call subtract_fields(a=lhs%rho_u, b=rhs%rho_u, c=difference%rho_u) ! c=a+b
    call subtract_fields(a=lhs%rho_v, b=rhs%rho_v, c=difference%rho_v) ! c=a+b
    call subtract_fields(a=lhs%rho_E, b=rhs%rho_E, c=difference%rho_E) ! c=a+b
    difference%prim_vars_updated = .false.
  end function subtract_fluid

  function add_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(local_fluid_t), intent(in) :: lhs
    class(local_fluid_t), intent(in) :: rhs
    type(local_fluid_t), allocatable :: sum

    call debug_print('Running local_fluid_t%add_fluid()', __FILE__, __LINE__)

    allocate(sum, source=lhs)
    call add_fields(a=lhs%rho, b=rhs%rho, c=sum%rho) ! c=a+b
    call add_fields(a=lhs%rho_u, b=rhs%rho_u, c=sum%rho_u) ! c=a+b
    call add_fields(a=lhs%rho_v, b=rhs%rho_v, c=sum%rho_v) ! c=a+b
    call add_fields(a=lhs%rho_E, b=rhs%rho_E, c=sum%rho_E) ! c=a+b
    sum%prim_vars_updated = .false.

    ! call sum%set_temp(calling_function='local_fluid_t%add_fluid(sum)', line=__LINE__)
    ! call move_alloc(sum, sum)
    ! call sum%set_temp(calling_function='local_fluid_t%add_fluid(sum)', line=__LINE__)
  end function add_fluid

  subroutine add_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a + b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), dimension(:, :), contiguous, intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    real(rk) :: diff, threshold
    integer(ik) :: i, j
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, diff, threshold) &
    !$omp shared(a,b,c)
    !$omp do
    do j = jlo, jhi
      !$omp simd
      do i = ilo, ihi
        c(i, j) = a(i, j) + b(i, j)
        diff = abs(c(i, j) - a(i, j))
        threshold = abs(a(i, j) + b(i, j)) * 1e-9_rk
        if(diff < threshold) c(i, j) = a(i, j)
      end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine add_fields

  subroutine subtract_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a - b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), dimension(:, :), contiguous, intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    real(rk) :: diff, threshold
    integer(ik) :: i, j
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, diff, threshold) &
    !$omp shared(a,b,c)
    !$omp do
    do j = jlo, jhi
      !$omp simd
      do i = ilo, ihi
        c(i, j) = a(i, j) - b(i, j)
        diff = abs(c(i, j) - a(i, j))
        threshold = abs(a(i, j) + b(i, j)) * 1e-9_rk
        if(diff < threshold) then
          c(i, j) = a(i, j)
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine subtract_fields

  subroutine mult_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), dimension(:, :), contiguous, intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    integer(ik) :: i, j
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(a,b,c)
    !$omp do
    do j = jlo, jhi

      !$omp simd

      do i = ilo, ihi
        c(i, j) = a(i, j) * b(i, j)
      end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine mult_fields

  subroutine mult_field_by_real(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    integer(ik) :: i, j
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(a,b,c)
    !$omp do
    do j = jlo, jhi

      !$omp simd
      do i = ilo, ihi
        c(i, j) = a(i, j) * b
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine mult_field_by_real

  function fluid_mul_real(lhs, rhs) result(product)
    !< Implementation of the fluid * real operation
    class(local_fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    type(local_fluid_t), allocatable :: product

    integer(ik) :: alloc_status

    call debug_print('Running local_fluid_t%fluid_mul_real()', __FILE__, __LINE__)

    allocate(product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) then
      call error_msg(module='mod_fluid', class='local_fluid_t', procedure='fluid_mul_real', &
                     message="Unable to allocate lhs%solver", file_name=__FILE__, line_number=__LINE__)
    end if

    call mult_field_by_real(a=lhs%rho, b=rhs, c=product%rho)
    call mult_field_by_real(a=lhs%rho_u, b=rhs, c=product%rho_u)
    call mult_field_by_real(a=lhs%rho_v, b=rhs, c=product%rho_v)
    call mult_field_by_real(a=lhs%rho_E, b=rhs, c=product%rho_E)
    product%prim_vars_updated = .false.
  end function fluid_mul_real

  function real_mul_fluid(lhs, rhs) result(product)
    !< Implementation of the real * fluid operation
    class(local_fluid_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    type(local_fluid_t), allocatable :: product
    integer(ik) :: alloc_status

    call debug_print('Running local_fluid_t%real_mul_fluid()', __FILE__, __LINE__)

    allocate(product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) then
      call error_msg(module='mod_fluid', class='local_fluid_t', procedure='real_mul_fluid', &
                     message="Unable to allocate product", file_name=__FILE__, line_number=__LINE__)
    end if

    call mult_field_by_real(a=rhs%rho, b=lhs, c=product%rho)
    call mult_field_by_real(a=rhs%rho_u, b=lhs, c=product%rho_u)
    call mult_field_by_real(a=rhs%rho_v, b=lhs, c=product%rho_v)
    call mult_field_by_real(a=rhs%rho_E, b=lhs, c=product%rho_E)
    product%prim_vars_updated = .false.
  end function real_mul_fluid

  subroutine assign_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(local_fluid_t), intent(inout) :: lhs
    type(local_fluid_t), intent(in) :: rhs
    integer(ik) :: error_code
    integer(ik) :: alloc_status

    call debug_print('Running local_fluid_t%assign_fluid()', __FILE__, __LINE__)
    lhs%residual_hist_header_written = rhs%residual_hist_header_written
    lhs%rho = rhs%rho
    lhs%rho_u = rhs%rho_u
    lhs%rho_v = rhs%rho_v
    lhs%rho_E = rhs%rho_E
    if(allocated(lhs%solver)) deallocate(lhs%solver)

    allocate(lhs%solver, source=rhs%solver, stat=alloc_status)
    if(alloc_status /= 0) then
      call error_msg(module='mod_fluid', class='local_fluid_t', procedure='assign_fluid', &
                     message="Unable to allocate lhs%solver", file_name=__FILE__, line_number=__LINE__)
    end if

    lhs%time_integration_scheme = rhs%time_integration_scheme
    lhs%time = rhs%time
    lhs%dt = rhs%dt
    lhs%iteration = rhs%iteration
    lhs%prim_vars_updated = rhs%prim_vars_updated
    lhs%smooth_residuals = rhs%smooth_residuals
    lhs%residual_hist_file = rhs%residual_hist_file
    lhs%residual_hist_header_written = rhs%residual_hist_header_written
  end subroutine assign_fluid

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(local_fluid_t), intent(in) :: self
    integer(ik) :: i, j, ilo, jlo, ihi, jhi
    logical :: invalid_numbers, negative_numbers
    character(len=64) :: err_message = ''
    integer(ik), intent(out) :: error_code

    call debug_print('Running local_fluid_t%sanity_check()', __FILE__, __LINE__)

    error_code = ALL_OK

    negative_numbers = .false.
    invalid_numbers = .false.

    ilo = lbound(self%rho, dim=1) + n_ghost_layers
    ihi = ubound(self%rho, dim=1) - n_ghost_layers
    jlo = lbound(self%rho, dim=2) + n_ghost_layers
    jhi = ubound(self%rho, dim=2) - n_ghost_layers

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(self, error_code, err_message)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(self%rho(i, j) < 0.0_rk) then
          write(err_message, '(a, i0, ", ", i0, a)') "Negative density at (", i, j, ")"
          call error_msg(module='mod_fluid', class='local_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NEG_DENSITY
        end if
      end do
    end do
    !$omp end do nowait

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(self%p(i, j) < 0.0_rk) then
          write(err_message, '(a, i0, ", ", i0, a)') "Negative pressure found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='local_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NEG_PRESSURE
        end if
      end do
    end do
    !$omp end do nowait

    ! NaN checks

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%rho(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN density found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='local_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%u(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN x-velocity found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='local_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%v(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN y-velocity found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='local_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%p(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN pressure found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='local_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait

    !$omp end parallel
  end subroutine sanity_check

end module mod_local_fluid
