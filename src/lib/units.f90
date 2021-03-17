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

module mod_units
  !< Summary:  Provide capability to handle units, nondimensionalization, etc.
  !< Author: Sam Miller

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  ! Defaulting to CGS, so everything is initially 1
  real(rk), protected :: io_velocity_units = 1.0_rk
  real(rk), protected :: io_time_units = 1.0_rk
  real(rk), protected :: io_length_units = 1.0_rk
  real(rk), protected :: io_pressure_units = 1.0_rk
  real(rk), protected :: io_mass_units = 1.0_rk
  real(rk), protected :: io_density_units = 1.0_rk
  real(rk), protected :: io_volume_units = 1.0_rk
  real(rk), protected :: io_temperature_units = 1.0_rk
  real(rk), protected :: io_energy_units = 1.0_rk
  real(rk), protected :: io_energy_density_units = 1.0_rk
  real(rk), protected :: io_force_units = 1.0_rk
  real(rk), protected :: io_power_units = 1.0_rk

  character(len=10), protected :: io_velocity_label = 'cm/s'
  character(len=10), protected :: io_time_label = 's'
  character(len=10), protected :: io_length_label = 'cm'
  character(len=10), protected :: io_pressure_label = 'barye'
  character(len=10), protected :: io_mass_label = 'g'
  character(len=10), protected :: io_density_label = 'g/cc'
  character(len=10), protected :: io_volume_label = 'cc'
  character(len=10), protected :: io_temperature_label = 'K'
  character(len=10), protected :: io_energy_label = 'erg'
  character(len=10), protected :: io_energy_density_label = 'erg/cc'
  character(len=10), protected :: io_force_label = 'dyne'
  character(len=10), protected :: io_power_label = 'erg/s'

  ! Time conversion
  real(rk), parameter :: pico_seconds = 1e12_rk !< sec to ps conversion
  real(rk), parameter :: nano_seconds = 1e9_rk !< sec to ps conversion

  real(rk), parameter :: ns_to_s = 1e-9_rk !< ns to sec conversion

  ! Velocity
  real(rk), parameter :: km_per_sec = 1e-5_rk !< cm/s to km/s conversion
  real(rk), parameter :: um_per_ns = 1e-5_rk !< cm/s to um/ns conversion
  real(rk), parameter :: km_per_sec_to_cm_per_sec = 1e5_rk !< km/s to cm/s conversion

  ! Distance
  real(rk), parameter :: microns = 1e4_rk !< cm to microns conversion
  real(rk), parameter :: cm_to_um = 1e4_rk !< cm to um conversion
  real(rk), parameter :: um_to_cm = 1e-4_rk !< um to cm conversion

  ! Pressure
  real(rk), parameter :: mega_bar = 1e-12_rk !< barye to Megabar conversion
  real(rk), parameter :: mega_bar_to_barye = 1e12_rk !< Mbar to barye convers8ihn

  real(rk), parameter :: atm = 9.86923e-7_rk !< barye to atm conversion
  real(rk), parameter :: giga_pascal = 1e-10_rk !< barye to GPa conversion

  ! Temperature
  real(rk), parameter :: kelvin_to_electron_volt = 1.160451812e-4_rk

  ! Energy
  real(rk), parameter :: joule = 1e-7_rk ! erg to joules

  ! Power
  real(rk), parameter :: watt = 1e-10_rk ! erg/s to watts

  character(len=3), protected :: unit_system
  ! ICF: length = microns, mass = g, time = ns, pressure = Mbar
  ! CGS: length = cm, mass = g, time = sec, pressure = barye

contains

  subroutine set_output_unit_system(system)
    character(len=*), intent(in) :: system

    select case(trim(system))
    case("cgs")
      if(this_image() == 1) write(*, '(a)') "Setting I/O unit system to CGS conventions (g/cc, cm, K, barye, erg)"
      unit_system = "cgs"
      io_velocity_units = 1.0_rk
      io_velocity_label = 'cm/s'
      io_time_units = 1.0_rk
      io_time_label = 's'
      io_length_units = 1.0_rk
      io_length_label = 'cm'
      io_pressure_units = 1.0_rk
      io_pressure_label = 'barye'
      io_mass_units = 1.0_rk
      io_mass_label = 'g'
      io_density_units = 1.0_rk
      io_density_label = 'g/cc'
      io_temperature_units = 1.0_rk
      io_temperature_label = 'K'
      io_volume_units = 1.0_rk
      io_volume_label = 'cc'
      io_energy_density_units = 1.0_rk
      io_energy_density_label = 'erg/cc'
      
    case("icf")
      if(this_image() == 1) write(*, '(a)') "Setting I/O unit system to ICF conventions (g/cc, microns, eV, Mbar, erg)"
      unit_system = "icf"
      io_velocity_units = km_per_sec
      io_velocity_label = 'km/s'
      io_time_units = nano_seconds
      io_time_label = 'ns'
      io_length_units = microns
      io_length_label = 'um'
      io_pressure_units = mega_bar
      io_pressure_label = 'Mbar'
      io_mass_units = 1.0_rk
      io_mass_label = 'g'
      io_density_units = 1.0_rk
      io_density_label = 'g/cc'
      io_temperature_units = kelvin_to_electron_volt
      io_temperature_label = 'eV'
      io_volume_units = 1.0_rk
      io_volume_label = 'cc'
      io_energy_density_units = 1.0_rk
      io_energy_density_label = 'erg/cc'

    endselect
  endsubroutine set_output_unit_system

endmodule mod_units
