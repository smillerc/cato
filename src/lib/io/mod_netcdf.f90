! !++++++++++++++++++++++++++++++++++++++++++++++++++++
! ! Example of writing simple profile data conforming
! ! to CF Conventions version 1.6, Appendix H.3.
! !
! ! Contributed by Frank Kreienkamp
! !++++++++++++++++++++++++++++++++++++++++++++++++++++

! program example_h3_cf1_6
!   use netcdf
!   implicit none

!   integer(kind=4), parameter :: iObs = 13, iLenStringName = 40
!   integer(kind=4) :: iNcid, i, iStationNo
!   integer(kind=4) :: iDimStation_ID, iDimTime_ID, iDimLenStringName_ID
!   integer(kind=4) :: iVarLON_ID, iVarLAT_ID, iVarAlt_ID, iVarStationName_ID, iVarStationInfo_ID, iVarStationElevation_ID
!   integer(kind=4) :: iVarTime_ID, iVarHumidity_ID, iVarTemp_ID
!   integer(kind=4) :: iStart(NF90_MAX_VAR_DIMS), iCount(NF90_MAX_VAR_DIMS)

!   real(4) :: rTemp(1), rField(iObs)

!   character(len=40) :: cStationName

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! Create file.
!   call check(nf90_create("example_h3_cf1_6.nc", NF90_CLOBBER, iNcid))
!   !call check(nf90_create("example_h3_cf1_6.nc", NF90_NETCDF4, iNcid))

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! Define the dimensions. The Station-record dimension is defined to have
!   ! unlimited length - it can grow as needed.
!   call check(nf90_def_dim(iNcid, "station", NF90_UNLIMITED, iDimStation_ID))
!   call check(nf90_def_dim(iNcid, "obs", iObs, iDimTime_ID))
!   call check(nf90_def_dim(iNcid, "NoCharStationName", iLenStringName, iDimLenStringName_ID))

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! Define of variables.
!   call check(nf90_def_var(iNcid, "lon", NF90_REAL4,(/iDimStation_ID/), iVarLON_ID))
!   call check(nf90_def_var(iNcid, "lat", NF90_REAL4,(/iDimStation_ID/), iVarLAT_ID))
!   call check(nf90_def_var(iNcid, "alt", NF90_REAL4,(/iDimStation_ID/), iVarAlt_ID))
!   call check(nf90_def_var(iNcid, "station_name", NF90_CHAR,(/iDimLenStringName_ID, iDimStation_ID/), iVarStationName_ID))
!   call check(nf90_def_var(iNcid, "station_info", NF90_INT4,(/iDimStation_ID/), iVarStationInfo_ID))
!   call check(nf90_def_var(iNcid, "station_elevation", NF90_REAL4,(/iDimStation_ID/), iVarStationElevation_ID))
!   call check(nf90_def_var(iNcid, "time", NF90_REAL8,(/iDimTime_ID, iDimStation_ID/), iVarTime_ID))
!   call check(nf90_def_var(iNcid, "humidity", NF90_REAL4,(/iDimTime_ID, iDimStation_ID/), iVarHumidity_ID))
!   call check(nf90_def_var(iNcid, "temp", NF90_REAL4,(/iDimTime_ID, iDimStation_ID/), iVarTemp_ID))

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! Assign attributes to the variables.
!   call check(nf90_put_att(iNcid, iVarLON_ID, "standard_name", "longitude"))
!   call check(nf90_put_att(iNcid, iVarLON_ID, "long_name", "station_longitude"))
!   call check(nf90_put_att(iNcid, iVarLON_ID, "units", "degrees_easth"))

!   call check(nf90_put_att(iNcid, iVarLAT_ID, "standard_name", "latitude"))
!   call check(nf90_put_att(iNcid, iVarLAT_ID, "long_name", "station_latitude"))
!   call check(nf90_put_att(iNcid, iVarLAT_ID, "units", "degrees_north"))

!   call check(nf90_put_att(iNcid, iVarAlt_ID, "long_name", "height"))
!   call check(nf90_put_att(iNcid, iVarAlt_ID, "standard_name", "vertical distance abovethe surface"))
!   call check(nf90_put_att(iNcid, iVarAlt_ID, "units", "m"))
!   call check(nf90_put_att(iNcid, iVarAlt_ID, "positve", "up"))
!   call check(nf90_put_att(iNcid, iVarAlt_ID, "axis", "Z"))

!   call check(nf90_put_att(iNcid, iVarStationName_ID, "long_name", "station name"))
!   call check(nf90_put_att(iNcid, iVarStationName_ID, "cf_role", "timeseries_id"))

!   call check(nf90_put_att(iNcid, iVarStationInfo_ID, "long_name", "any kind of station info"))

!   call check(nf90_put_att(iNcid, iVarStationElevation_ID, "long_name", "height above the geoid"))
!   call check(nf90_put_att(iNcid, iVarStationElevation_ID, "standard_name", "surface_altitude"))
!   call check(nf90_put_att(iNcid, iVarStationElevation_ID, "units", "m"))

!   call check(nf90_put_att(iNcid, iVarTime_ID, "long_name", "time"))
!   call check(nf90_put_att(iNcid, iVarTime_ID, "standard_name", "time of measurement"))
!   call check(nf90_put_att(iNcid, iVarTime_ID, "units", "days since 1970-01-01 00:00:00"))
!   call check(nf90_put_att(iNcid, iVarTime_ID, "missing_value", -999.9))

!   call check(nf90_put_att(iNcid, iVarHumidity_ID, "standard_name", "specific_humidity"))
!   call check(nf90_put_att(iNcid, iVarHumidity_ID, "coordinates", "time lat lon alt"))
!   call check(nf90_put_att(iNcid, iVarHumidity_ID, "_FillValue", -999.9))

!   call check(nf90_put_att(iNcid, iVarTemp_ID, "standard_name", "air_temperature"))
!   call check(nf90_put_att(iNcid, iVarTemp_ID, "units", "Celsius"))
!   call check(nf90_put_att(iNcid, iVarTemp_ID, "coordinates", "time lat lon alt"))
!   call check(nf90_put_att(iNcid, iVarTemp_ID, "_FillValue", -999.9))

!   call check(nf90_put_att(iNcid, nf90_global, "featureType", "timeSeries"))

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! End define mode.
!   call check(nf90_enddef(iNcid))

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! Adding 1. station
!   iStationNo = 1
!   istart = 0
!   iCount = 0
!   istart(1) = iStationNo
!   iCount(1) = 1
!   rTemp(1) = 13.32
!   call check(nf90_put_var(iNcid, iVarLON_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   rTemp(1) = 52.45
!   call check(nf90_put_var(iNcid, iVarLAT_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   rTemp(1) = 2
!   call check(nf90_put_var(iNcid, iVarAlt_ID, rTemp(:), istart(1:1), iCount(1:1)))

!   rTemp(1) = 10381
!   call check(nf90_put_var(iNcid, iVarStationInfo_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   rTemp(1) = 51
!   call check(nf90_put_var(iNcid, iVarStationElevation_ID, rTemp(:), istart(1:1), iCount(1:1)))

!   istart = 0
!   iCount = 0
!   istart(1) = 1
!   iCount(1) = iLenStringName
!   istart(2) = iStationNo
!   iCount(2) = 1
!   cStationName = ""
!   cStationName = "Berlin-Dahlem (FU)"//CHAR(0)
!   call check(nf90_put_var(iNcid, iVarStationName_ID, trim(cStationName), istart(1:2), iCount(1:2)))

!   istart(1) = 1
!   iCount(1) = iObs
!   istart(2) = iStationNo
!   iCount(2) = 1
!   do i = 1, iObs
!     rField(i) = i
!   end do ! i
!   call check(nf90_put_var(iNcid, iVarTime_ID, rField(:), istart(1:2), iCount(1:2)))

!   do i = 1, iObs
!     rField(i) = 70.+i
!   end do ! i
!   call check(nf90_put_var(iNcid, iVarHumidity_ID, rField(:), istart(1:2), iCount(1:2)))

!   do i = 1, iObs
!     rField(i) = 25.-i
!   end do ! i
!   call check(nf90_put_var(iNcid, iVarTemp_ID, rField(:), istart(1:2), iCount(1:2)))

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! Adding 2. station
!   iStationNo = 2
!   istart = 0
!   iCount = 0
!   istart(1) = iStationNo
!   iCount(1) = 1
!   rTemp(1) = 13.07
!   call check(nf90_put_var(iNcid, iVarLON_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   rTemp(1) = 52.38
!   call check(nf90_put_var(iNcid, iVarLAT_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   rTemp(1) = 2
!   call check(nf90_put_var(iNcid, iVarAlt_ID, rTemp(:), istart(1:1), iCount(1:1)))

!   rTemp(1) = 10379
!   call check(nf90_put_var(iNcid, iVarStationInfo_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   rTemp(1) = 81
!   call check(nf90_put_var(iNcid, iVarStationElevation_ID, rTemp(:), istart(1:1), iCount(1:1)))

!   istart(1) = 1
!   iCount(1) = iLenStringName
!   istart(2) = iStationNo
!   iCount(2) = 1
!   cStationName = ""
!   cStationName = "Potsdam"//CHAR(0)
!   call check(nf90_put_var(iNcid, iVarStationName_ID, trim(cStationName), istart(1:2), iCount(1:2)))

!   istart(1) = 1
!   iCount(1) = iObs
!   istart(2) = iStationNo
!   iCount(2) = 1
!   do i = 1, iObs
!     rField(i) = i
!   end do ! i
!   call check(nf90_put_var(iNcid, iVarTime_ID, rField(:), istart(1:2), iCount(1:2)))

!   do i = 1, iObs
!     rField(i) = 70.-i
!   end do ! i
!   call check(nf90_put_var(iNcid, iVarHumidity_ID, rField(:), istart(1:2), iCount(1:2)))

!   do i = 1, iObs
!     rField(i) = 15.+i
!   end do ! i
!   call check(nf90_put_var(iNcid, iVarTemp_ID, rField(:), istart(1:2), iCount(1:2)))

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! reading station 1
!   iStationNo = 1
!   istart = 0
!   iCount = 0
!   istart(1) = iStationNo
!   iCount(1) = 1
!   call check(nf90_get_var(iNcid, iVarLON_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, " lon: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarLAT_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, ", lat: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarAlt_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, ", alt: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarStationInfo_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f10.2)') "Station ", iStationNo, ", StationsInfo: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarStationElevation_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, ", StationElevation: ", rTemp(:)

!   istart(1) = 1
!   iCount(1) = iLenStringName
!   istart(2) = iStationNo
!   iCount(2) = 1
!   cStationName = ""
!   call check(nf90_get_var(iNcid, iVarStationName_ID, cStationName, istart(1:2), iCount(1:2)))
!   write(6, '(1a,1i1,1a,1a)') "Station ", iStationNo, ", StationName: ", trim(cStationName)

!   istart(1) = 1
!   iCount(1) = iObs
!   istart(2) = iStationNo
!   iCount(2) = 1
!   call check(nf90_get_var(iNcid, iVarTime_ID, rField(:), istart(1:2), iCount(1:2)))
!   write(6, '(1a,1i1,1a,13f7.2)') "Station ", iStationNo, ", Time: ", rField(:)
!   call check(nf90_get_var(iNcid, iVarHumidity_ID, rField(:), istart(1:2), iCount(1:2)))
!   write(6, '(1a,1i1,1a,13f7.2)') "Station ", iStationNo, ", Humidity: ", rField(:)
!   call check(nf90_get_var(iNcid, iVarTemp_ID, rField(:), istart(1:2), iCount(1:2)))
!   write(6, '(1a,1i1,1a,13f7.2)') "Station ", iStationNo, ", Temp: ", rField(:)

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! reading station 2
!   iStationNo = 2
!   istart = 0
!   iCount = 0
!   istart(1) = iStationNo
!   iCount(1) = 1
!   call check(nf90_get_var(iNcid, iVarLON_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, " lon: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarLAT_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, ", lat: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarAlt_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, ", alt: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarStationInfo_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f10.2)') "Station ", iStationNo, ", StationsInfo: ", rTemp(:)
!   call check(nf90_get_var(iNcid, iVarStationElevation_ID, rTemp(:), istart(1:1), iCount(1:1)))
!   write(6, '(1a,1i1,1a,1f7.2)') "Station ", iStationNo, ", StationElevation: ", rTemp(:)

!   istart(1) = 1
!   iCount(1) = iLenStringName
!   istart(2) = iStationNo
!   iCount(2) = 1
!   cStationName = ""
!   call check(nf90_get_var(iNcid, iVarStationName_ID, cStationName, istart(1:2), iCount(1:2)))
!   !write(6,'(1a,1i1,1a,1a)') "Station ",iStationNo,", StationName: ",trim(cStationName)
!   !write(6,'(1a,1i1,1a,1a)') "Station ",iStationNo,", StationName: ",cStationName(1:iCount(1))
!   write(6, '(1a,1i1,1a,1a)') "Station ", iStationNo, ", StationName: ", cStationName(1:(index(cStationName, char(0)) - 1))

!   istart(1) = 1
!   iCount(1) = iObs
!   istart(2) = iStationNo
!   iCount(2) = 1
!   call check(nf90_get_var(iNcid, iVarTime_ID, rField(:), istart(1:2), iCount(1:2)))
!   write(6, '(1a,1i1,1a,13f7.2)') "Station ", iStationNo, ", Time: ", rField(:)
!   call check(nf90_get_var(iNcid, iVarHumidity_ID, rField(:), istart(1:2), iCount(1:2)))
!   write(6, '(1a,1i1,1a,13f7.2)') "Station ", iStationNo, ", Humidity: ", rField(:)
!   call check(nf90_get_var(iNcid, iVarTemp_ID, rField(:), istart(1:2), iCount(1:2)))
!   write(6, '(1a,1i1,1a,13f7.2)') "Station ", iStationNo, ", Temp: ", rField(:)

!   !++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! Close file.
!   call check(nf90_close(iNcid))

! end program example_h3_cf1_6

! subroutine check(status)
!   use netcdf
!   IMPLICIT NONE
!   !---------------------------------------------------------------------
!   integer(4), intent(in) :: status
!   !
!   !- End of header -----------------------------------------------------
!   !---------------------------------------------------------------------
!   ! Subroutine Body
!   !---------------------------------------------------------------------

!   if(status /= nf90_noerr) then
!     print *, trim(nf90_strerror(status))
!     stop "Stopped"
!   end if
! end subroutine check

! !     This is part of the netCDF package.
! !     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! !     See COPYRIGHT file for conditions of use.

! !     This is a very simple example which writes a 2D array of
! !     sample data. To handle this in netCDF we create two shared
! !     dimensions, "x" and "y", and a netCDF variable, called "data".

! !     This example demonstrates the netCDF Fortran 90 API. This is part
! !     of the netCDF tutorial, which can be found at:
! !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! !     Full documentation of the netCDF Fortran 90 API can be found at:
! !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! !     $Id: simple_xy_wr.f90,v 1.7 2006/12/09 18:44:58 russ Exp $

! program simple_xy_wr
!   use netcdf
!   implicit none

!   ! This is the name of the data file we will create.
!   character(len=*), parameter :: FILE_NAME = "simple_xy.nc"

!   ! We are writing 2D data, a 6 x 12 grid.
!   integer, parameter :: NDIMS = 2
!   integer, parameter :: NX = 6, NY = 12

!   ! When we create netCDF files, variables and dimensions, we get back
!   ! an ID for each one.
!   integer :: ncid, varid, dimids(NDIMS)
!   integer :: x_dimid, y_dimid

!   ! This is the data array we will write. It will just be filled with
!   ! a progression of integers for this example.
!   integer :: data_out(NY, NX)

!   ! Loop indexes, and error handling.
!   integer :: x, y

!   ! Create some pretend data. If this wasn't an example program, we
!   ! would have some real data to write, for example, model output.
!   do x = 1, NX
!     do y = 1, NY
!       data_out(y, x) = (x - 1) * NY + (y - 1)
!     end do
!   end do

!   ! Always check the return code of every netCDF function call. In
!   ! this example program, wrapping netCDF calls with "call check()"
!   ! makes sure that any return which is not equal to nf90_noerr (0)
!   ! will print a netCDF error message and exit.

!   ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
!   ! overwrite this file, if it already exists.
!   call check(nf90_create(FILE_NAME, NF90_CLOBBER, ncid))

!   ! Define the dimensions. NetCDF will hand back an ID for each.
!   call check(nf90_def_dim(ncid, "x", NX, x_dimid))
!   call check(nf90_def_dim(ncid, "y", NY, y_dimid))

!   ! The dimids array is used to pass the IDs of the dimensions of
!   ! the variables. Note that in fortran arrays are stored in
!   ! column-major format.
!   dimids = (/y_dimid, x_dimid/)

!   ! Define the variable. The type of the variable in this case is
!   ! NF90_INT (4-byte integer).
!   call check(nf90_def_var(ncid, "data", NF90_INT, dimids, varid))

!   ! End define mode. This tells netCDF we are done defining metadata.
!   call check(nf90_enddef(ncid))

!   ! Write the pretend data to the file. Although netCDF supports
!   ! reading and writing subsets of data, in this case we write all the
!   ! data in one operation.
!   call check(nf90_put_var(ncid, varid, data_out))

!   ! Close the file. This frees up any internal netCDF resources
!   ! associated with the file, and flushes any buffers.
!   call check(nf90_close(ncid))

!   print *, "*** SUCCESS writing example file simple_xy.nc! "

! contains
!   subroutine check(status)
!     integer, intent(in) :: status

!     if(status /= nf90_noerr) then
!       print *, trim(nf90_strerror(status))
!       stop "Stopped"
!     end if
!   end subroutine check
! end program simple_xy_wr
