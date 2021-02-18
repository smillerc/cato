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

module mod_energy_source
  !< Define the class for injecting energy into the domain

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  use mod_globals, only: debug_print
  use mod_field, only: field_2d_t, field_2d
  use mod_units
  use math_constants, only: universal_gas_const
  use mod_source, only: source_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use linear_interpolation_module, only: linear_interp_1d

  implicit none

  private
  public :: energy_source_t, new_energy_source

  type, extends(source_t) :: energy_source_t
  contains
    procedure :: integrate
    final :: finalize
  endtype

contains

  function new_energy_source(input) result(energy_source)
    type(energy_source_t), pointer :: energy_source
    class(input_t), intent(in) :: input
    allocate(energy_source)

    energy_source%source_type = 'energy'
    energy_source%constant_source = input%apply_constant_source
    if(.not. energy_source%constant_source) then
      energy_source%input_filename = trim(input%source_file)
      call energy_source%read_input_file()
    else
      energy_source%constant_source_value = input%constant_source_value
    endif

    open(newunit=energy_source%io_unit, file='input_source_value.dat')
    write(energy_source%io_unit, '(a)') 'Time [sec] Energy [erg]'
  endfunction

  subroutine finalize(self)
    type(energy_source_t), intent(inout) :: self
    logical :: is_open = .false.

    ! if(allocated(self%data)) deallocate(self%data)
    if(allocated(self%source_type)) deallocate(self%source_type)
    if(allocated(self%source_geometry)) deallocate(self%source_geometry)
    if(allocated(self%input_filename)) deallocate(self%input_filename)

    inquire(unit=self%io_unit, opened=is_open)
    if(is_open) close(self%io_unit)

  endsubroutine finalize

  type(field_2d_t) function integrate(self, time, dt, density) result(d_dt)
    !< Create the source field to be passed to the fluid class or others
    class(energy_source_t), intent(inout) :: self
    class(field_2d_t), intent(in) :: density
    real(rk), intent(in) :: dt, time
  end function integrate


endmodule mod_energy_source
