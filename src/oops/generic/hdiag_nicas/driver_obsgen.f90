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
use type_geom, only: geomtype
use type_mpl, only: mpl,mpl_bcast
use type_nam, only: namtype
use type_odata, only: odatatype
use type_randgen, only: rand_real

implicit none

logical,parameter :: readobs = .false. !< Read observations
logical,parameter :: allobs = .true.   !< All observation are used

private
public :: run_obsgen

contains

!----------------------------------------------------------------------
! Subroutine: run_obsgen
!> Purpose: generate random observations locations
!----------------------------------------------------------------------
subroutine run_obsgen(nam,geom,odata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam        !< Namelist
type(geomtype),intent(in) :: geom      !< Geometry
type(odatatype),intent(inout) :: odata !< Observation operator data

! Local variables
integer :: info,iobs,active
real(kind_real) :: lat,lon

if (nam%new_obsop) then
   ! Define number of observations
   if (readobs) then
      ! Read observation network
      odata%nobs = 0
      open(unit=100,file=trim(nam%datadir)//'/'//trim(nam%prefix)//'_obs_in.dat',status='old')
      do
         read(100,*,iostat=info) lat,lon,active
         if (info/=0) exit
         if (allobs.or.(active==1)) odata%nobs = odata%nobs+1
      end do
      close(unit=100)
      odata%nobs = min(odata%nobs,nam%nobs)
      if (odata%nobs<1) call msgerror('no active observation in the file provided')
   else
      ! Generate random observation network
      odata%nobs = nam%nobs
   end if

   ! Allocation
   allocate(odata%lonobs(odata%nobs))
   allocate(odata%latobs(odata%nobs))

   ! Get observation locations
   if (mpl%main) then
      if (readobs) then
         ! Read observation network
         open(unit=100,file=trim(nam%datadir)//'/'//trim(nam%prefix)//'_obs_in.dat',status='old')
         iobs = 0
         do while (iobs<odata%nobs)
            read(100,*) lat,lon,active
            if (allobs.or.(active==1)) then
               iobs = iobs+1
               odata%lonobs(iobs) = lon
               odata%latobs(iobs) = lat
            end if
         end do
         close(unit=100)
      else
         ! Generate random observation network
         if (.true.) then
            ! Limited-area domain
            call rand_real(minval(geom%lon),maxval(geom%lon),odata%lonobs)
            call rand_real(minval(geom%lat),maxval(geom%lat),odata%latobs)
         else
            ! Global domain
            call rand_real(-pi,pi,odata%lonobs)
            call rand_real(-1.0_kind_real,1.0_kind_real,odata%latobs)
            odata%latobs = 0.5*pi-acos(odata%latobs)
         end if
      end if
   end if

   ! Broadcast data
   call mpl_bcast(odata%lonobs,mpl%ioproc)
   call mpl_bcast(odata%latobs,mpl%ioproc)

   ! Print results
   write(mpl%unit,'(a7,a,i8)') '','Number of observations: ',odata%nobs
end if

end subroutine run_obsgen

end module driver_obsgen
