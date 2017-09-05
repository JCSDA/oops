!----------------------------------------------------------------------
! Module: type_odata
!> Purpose: observation operator data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_odata

use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_com, only: comtype
use type_geom, only: geomtype
use type_linop, only: linoptype
use type_mpl, only: mpl,mpl_recv,mpl_send

implicit none

! Observation operator data derived type
type odatatype
   ! Geometry
   type(geomtype),pointer :: geom                  !< Geometry

   ! Number of points
   integer :: nc0
   integer :: nl0
   integer :: nobs

   ! Interpolation data
   type(linoptype) :: interp

   ! MPI distribution
   integer,allocatable :: iobs_to_iproc(:)
   integer,allocatable :: iobs_to_iobsa(:)
   integer,allocatable :: iproc_to_nobsa(:)
end type odatatype

! Local observation operator data derived type
type odataloctype
   ! Number of points
   integer :: nc0a
   integer :: nc0b
   integer :: nl0
   integer :: nobsa

   ! Interpolation data
   type(linoptype) :: interp

   ! Communication data
   type(comtype) :: com
end type odataloctype

private
public :: odatatype,odataloctype,yobs_com_gl,yobs_com_lg

contains

!----------------------------------------------------------------------
! Subroutine: yobs_com_gl
!> Purpose: communicate full field from global to local distribution
!----------------------------------------------------------------------
subroutine yobs_com_gl(odata,yobs)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata       !< Sampling data
real(kind_real),allocatable,intent(inout) :: yobs(:,:)        !< Field

! Local variables
integer :: iobs,iobsa,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      allocate(sbuf(odata%iproc_to_nobsa(iproc)*odata%nl0))

      ! Initialization
      call msr(sbuf)

      ! Prepare buffer
      do iobs=1,odata%nobs
         jproc = odata%iobs_to_iproc(iobs)
         if (jproc==iproc) then
            iobsa = odata%iobs_to_iobsa(iobs)
            sbuf((iobsa-1)*odata%nl0+1:iobsa*odata%nl0) = yobs(iobs,:)
         end if
      end do

      if (iproc==mpl%ioproc) then
         ! Allocation
         allocate(rbuf(odata%iproc_to_nobsa(iproc)*odata%nl0))

         ! Copy data
         rbuf = sbuf
      else
         ! Send data to iproc
         call mpl_send(odata%iproc_to_nobsa(iproc)*odata%nl0,sbuf,iproc,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Allocation
   allocate(rbuf(odata%iproc_to_nobsa(mpl%myproc)*odata%nl0))

   ! Receive data from ioproc
   call mpl_recv(odata%iproc_to_nobsa(mpl%myproc)*odata%nl0,rbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Rellocation
if (allocated(yobs)) deallocate(yobs)
allocate(yobs(odata%iproc_to_nobsa(mpl%myproc),odata%nl0))

! Copy from buffer
do iobsa=1,odata%iproc_to_nobsa(mpl%myproc)
   yobs(iobsa,:) = rbuf((iobsa-1)*odata%nl0+1:iobsa*odata%nl0)
end do

end subroutine yobs_com_gl

!----------------------------------------------------------------------
! Subroutine: yobs_com_lg
!> Purpose: communicate full field from local to global distribution
!----------------------------------------------------------------------
subroutine yobs_com_lg(odata,yobs)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata       !< Sampling data
real(kind_real),allocatable,intent(inout) :: yobs(:,:)        !< Field

! Local variables
integer :: iobs,iobsa,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Allocation
allocate(sbuf(odata%iproc_to_nobsa(mpl%myproc)*odata%nl0))

! Prepare buffer
do iobsa=1,odata%iproc_to_nobsa(mpl%myproc)
   sbuf((iobsa-1)*odata%nl0+1:iobsa*odata%nl0) = yobs(iobsa,:)
end do

! Release memory
deallocate(yobs)

! Communication
if (mpl%main) then
   ! Allocation
   allocate(yobs(odata%nobs,odata%nl0))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(rbuf(odata%iproc_to_nobsa(iproc)*odata%nl0))

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Receive data from iproc
         call mpl_recv(odata%iproc_to_nobsa(iproc)*odata%nl0,rbuf,iproc,mpl%tag)
      end if

      ! Copy from buffer
      do iobs=1,odata%nobs
         jproc = odata%iobs_to_iproc(iobs)
         if (jproc==iproc) then
            iobsa = odata%iobs_to_iobsa(iobs)
            yobs(iobs,:) = rbuf((iobsa-1)*odata%nl0+1:iobsa*odata%nl0)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Sending data to iproc
   call mpl_send(odata%iproc_to_nobsa(mpl%myproc)*odata%nl0,sbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

end subroutine yobs_com_lg

end module type_odata
