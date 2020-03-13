module mod_units
  !< Convert from CGS to some other unit

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

  character(len=10), protected :: io_velocity_label = 'cm/s'
  character(len=10), protected :: io_time_label = 's'
  character(len=10), protected :: io_length_label = 'cm'
  character(len=10), protected :: io_pressure_label = 'barye'
  character(len=10), protected :: io_mass_label = 'g'
  character(len=10), protected :: io_density_label = 'g/cc'
  character(len=10), protected :: io_volume_label = 'cc'
  character(len=10), protected :: io_temperature_label = 'K'

  ! Time conversion
  real(rk), parameter :: pico_seconds = 1e12_rk !< sec to ps conversion
  real(rk), parameter :: nano_seconds = 1e9_rk !< sec to ps conversion

  ! Vel
  real(rk), parameter :: km_per_sec = 1e-5_rk !< cm/s to km/s conversion
  real(rk), parameter :: um_per_ns = 1e-5_rk !< cm/s to um/ns conversion

  ! Distance
  real(rk), parameter :: microns = 1e4_rk !< cm to microns conversion

  ! Pressure
  real(rk), parameter :: megabar = 1e-12_rk !< barye to megabar conversion

  ! Temperature
  real(rk), parameter :: kelvin_to_electron_volt = 1.160451812e4_rk

  character(len=3), protected :: unit_system
  ! ICF: length = microns, mass = g, time = ns, pressure = Mbar
  ! CGS: length = cm, mass = g, time = sec, pressure = barye

contains

  subroutine set_output_unit_system(system)
    character(len=*), intent(in) :: system

    select case(trim(system))
    case("cgs")
      write(*, '(a)') "Setting I/O unit system to CGS conventions (g/cc, cm, K, barye)"
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

    case("icf")
      write(*, '(a)') "Setting I/O unit system to ICF conventions (g/cc, microns, eV, Mbar)"
      unit_system = "icf"
      io_velocity_units = km_per_sec
      io_velocity_label = 'km/s'
      io_time_units = nano_seconds
      io_time_label = 'ns'
      io_length_units = microns
      io_length_label = 'um'
      io_pressure_units = megabar
      io_pressure_label = 'Mbar'
      io_mass_units = 1.0_rk
      io_mass_label = 'g'
      io_density_units = 1.0_rk
      io_density_label = 'g/cc'
      io_temperature_units = kelvin_to_electron_volt
      io_temperature_label = 'eV'
      io_volume_units = 1.0_rk
      io_volume_label = 'cc'

    end select
  end subroutine set_output_unit_system

end module mod_units