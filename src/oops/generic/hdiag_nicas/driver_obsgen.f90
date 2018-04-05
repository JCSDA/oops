!----------------------------------------------------------------------
! Module: driver_obsgen
!> Purpose: generate random observations locations
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_obsgen

use tools_const, only: pi
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type
use type_obsop, only: obsop_type
use type_rng, only: rng

implicit none

private
public :: run_obsgen

contains

!----------------------------------------------------------------------
! Subroutine: run_obsgen
!> Purpose: generate random observations locations
!----------------------------------------------------------------------
subroutine run_obsgen(nam,geom,obsop)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(obsop_type),intent(inout) :: obsop !< Observation operator data

! Local variables
integer :: info,iobs
real(kind_real) :: lat,lon
logical :: readobs
character(len=1024) :: filename

if (nam%new_obsop) then
   ! Define number of observations
   filename = trim(nam%datadir)//'/'//trim(nam%prefix)//'_obs_in.dat'
   inquire(file=trim(filename),exist=readobs)
   if (readobs) then
      ! Read observation network
      obsop%nobs = 0
      open(unit=100,file=trim(filename),status='old')
      do
         read(100,*,iostat=info) lat,lon
         if (info/=0) exit
         obsop%nobs = obsop%nobs+1
      end do
      close(unit=100)
      obsop%nobs = min(obsop%nobs,nam%nobs)
      if (obsop%nobs<1) call msgerror('no observation in the file provided')
   else
      ! Generate random observation network
      obsop%nobs = nam%nobs
   end if

   ! Allocation
   allocate(obsop%lonobs(obsop%nobs))
   allocate(obsop%latobs(obsop%nobs))

   ! Get observation locations
   if (mpl%main) then
      if (readobs) then
         ! Read observation network
         open(unit=100,file=trim(nam%datadir)//'/'//trim(nam%prefix)//'_obs_in.dat',status='old')
         iobs = 0
         do while (iobs<obsop%nobs)
            read(100,*) lat,lon
            iobs = iobs+1
            obsop%lonobs(iobs) = lon
            obsop%latobs(iobs) = lat
         end do
         close(unit=100)
      else
         ! Generate random observation network
         if (abs(maxval(geom%area)-4.0*pi)<1.0e-1) then
            ! Limited-area domain
            call rng%rand_real(minval(geom%lon),maxval(geom%lon),obsop%lonobs)
            call rng%rand_real(minval(geom%lat),maxval(geom%lat),obsop%latobs)
         else
            ! Global domain
            call rng%rand_real(-pi,pi,obsop%lonobs)
            call rng%rand_real(-1.0_kind_real,1.0_kind_real,obsop%latobs)
            obsop%latobs = 0.5*pi-acos(obsop%latobs)
         end if
      end if
   end if

   ! Broadcast data
   call mpl%bcast(obsop%lonobs,mpl%ioproc)
   call mpl%bcast(obsop%latobs,mpl%ioproc)

   ! Print results
   write(mpl%unit,'(a7,a,i8)') '','Number of observations: ',obsop%nobs
end if

end subroutine run_obsgen

end module driver_obsgen
