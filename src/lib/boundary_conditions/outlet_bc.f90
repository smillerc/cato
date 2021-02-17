! MIT License
! Copyright (c) 2019 Sam Miller
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

module mod_outlet_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_error, only: error_msg
  use mod_globals, only: enable_debug_print, debug_print
  use mod_field, only: field_2d_t
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_nondimensionalization

  implicit none

  private
  public :: outlet_bc_t, outlet_bc_constructor

  type, extends(boundary_condition_t) :: outlet_bc_t
    real(rk) :: outlet_pressure = 0.0_rk
  contains
    procedure, public :: apply => apply_outlet_primitive_var_bc
    final :: finalize
  endtype
contains

  function outlet_bc_constructor(location, input, grid) result(bc)
    type(outlet_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid

    allocate(bc)
    bc%name = 'outlet'
    bc%location = location
    bc%priority = 0
    bc%outlet_pressure = input%constant_bc_pressure_value * press_to_nondim

    if (bc%outlet_pressure <= 0.0_rk) then
      call error_msg(module_name='mod_outlet_bc', class_name='outlet_bc_t', &
                    procedure_name='outlet_bc_constructor', &
                    message="Outlet ambient pressure is <= 0!", &
                    file_name=__FILE__, line_number=__LINE__)
    endif

    call bc%set_indices(grid) 
  endfunction outlet_bc_constructor

  subroutine apply_outlet_primitive_var_bc(self, rho, u, v, p)

    class(outlet_bc_t), intent(inout) :: self
    class(field_2d_t), intent(inout) :: rho  !< density
    class(field_2d_t), intent(inout) :: u    !< x-velocity
    class(field_2d_t), intent(inout) :: v    !< y-velocity
    class(field_2d_t), intent(inout) :: p    !< pressure

    real(rk) :: cs, gamma, vel, nx, ny, mach

    integer(ik) :: i, j, ihi, ihi_halo

    gamma = eos%get_gamma()

    select case(self%location)
    case('+x')

      if(rho%on_ihi_bc) then
        if(enable_debug_print) then
          call debug_print('Running outlet_bc_t%apply_outlet_primitive_var_bc() +x', &
                           __FILE__, __LINE__)
        endif
      
        ihi = rho%ihi 
        ihi_halo = rho%ihi_halo

        do i = ihi+1, ihi_halo
          do j = rho%jlo_halo, rho%jhi_halo

            associate(p_d => p%data(ihi, j), p_b => self%outlet_pressure, rho_0 => rho%data(ihi, j))
              cs = sqrt(gamma * p_d / rho_0)
              vel = sqrt(u%data(ihi,j)**2 + v%data(ihi,j)**2)
              ! print*, i, j, 'cs', cs, 'vel', vel, 'p_d', p_d, 'p_b', p_b, 'rho_0', rho_0
              mach = vel / cs

              nx = u%data(ihi, j) / (vel + 1e-30_rk)
              ny = v%data(ihi, j) / (vel + 1e-30_rk)

              ! print*, 'nx, ny', nx, ny, 'outlet mach', mach

              if (mach > 1.0_rk) then
                ! Supersonic extrapolation
                rho%data(i,j) = rho%data(ihi, j)
                u%data  (i,j) =   u%data(ihi, j)
                v%data  (i,j) =   v%data(ihi, j)
                p%data  (i,j) =   p%data(ihi, j)

              else !if (mach >= 0.0_rk .and. mach < 1.0_rk) then
                ! Subsonic extrapolation
                rho%data(i,j) = rho%data(ihi, j) + (p_b - p_d) / cs**2
                u%data  (i,j) =   u%data(ihi, j) + nx * (p_d - p_b)/(rho_0 * cs)
                v%data  (i,j) =   v%data(ihi, j) + ny * (p_d - p_b)/(rho_0 * cs)
                p%data  (i,j) =   self%outlet_pressure
              ! else
              !   call error_msg(module_name='mod_outlet_bc', class_name='outlet_bc_t', &
              !                       procedure_name='apply_outlet_primitive_var_bc', &
              !                       message="Inflow at the +x outflow BC!", &
              !                       file_name=__FILE__, line_number=__LINE__)
              end if
            end associate
          end do

          
        enddo
      endif
    case('-x')
      if(rho%on_ilo_bc) then
        if(enable_debug_print) then
          call debug_print('Running outlet_bc_t%apply_outlet_primitive_var_bc() -x', &
                           __FILE__, __LINE__)
        endif
        do i = 1, self%n_ghost_layers
          rho%data(self%ilo_ghost(i), :) = rho%data(self%ilo, :)
          u%data(self%ilo_ghost(i), :) = u%data(self%ilo, :)
          v%data(self%ilo_ghost(i), :) = v%data(self%ilo, :)
          p%data(self%ilo_ghost(i), :) = p%data(self%ilo, :)
        enddo
      endif
    case('+y')
      if(rho%on_jhi_bc) then
        if(enable_debug_print) then
          call debug_print('Running outlet_bc_t%apply_outlet_primitive_var_bc() +y', &
                           __FILE__, __LINE__)
        endif
        do i = 1, self%n_ghost_layers
          rho%data(:, self%jhi_ghost(i)) = rho%data(:, self%jhi)
          u%data(:, self%jhi_ghost(i)) = u%data(:, self%jhi)
          v%data(:, self%jhi_ghost(i)) = v%data(:, self%jhi)
          p%data(:, self%jhi_ghost(i)) = p%data(:, self%jhi)
        enddo
      endif
    case('-y')
      if(rho%on_jlo_bc) then
        if(enable_debug_print) then
          call debug_print('Running outlet_bc_t%apply_outlet_primitive_var_bc() -y', &
                           __FILE__, __LINE__)
        endif
        do i = 1, self%n_ghost_layers
          rho%data(:, self%jlo_ghost(i)) = rho%data(:, self%jlo)
          u%data(:, self%jlo_ghost(i)) = u%data(:, self%jlo)
          v%data(:, self%jlo_ghost(i)) = v%data(:, self%jlo)
          p%data(:, self%jlo_ghost(i)) = p%data(:, self%jlo)
        enddo
      endif
    case default
      error stop "Unsupported location to apply the bc at in outlet_bc_t%apply_outlet_cell_gradient_bc()"
    endselect

  endsubroutine apply_outlet_primitive_var_bc

  subroutine finalize(self)
    type(outlet_bc_t), intent(inout) :: self
    if(enable_debug_print) call debug_print('Running outlet_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  endsubroutine finalize

endmodule mod_outlet_bc
