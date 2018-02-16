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
use type_nam, only: namtype

implicit none

! Observation operator data derived type
type odatatype
   ! Namelist
   type(namtype),pointer :: nam             !< Namelist

   ! Geometry
   type(geomtype),pointer :: geom           !< Geometry

   ! Observations
   integer :: nobs                          !< Number of observations
   real(kind_real),allocatable :: lonobs(:) !< Observations longitudes
   real(kind_real),allocatable :: latobs(:) !< Observations latitudes

   ! Interpolation data
   type(linoptype) :: hfull                 !< Full interpolation data

   ! MPI distribution
   integer,allocatable :: obs_to_proc(:)    !< Observation to processor
   integer,allocatable :: obs_to_obsa(:)    !< Observation to local observation
   integer,allocatable :: proc_to_nobsa(:)  !< Processor to local number of observations
   integer,allocatable :: c0b_to_c0(:)      !< Subset Sc0, halo B to global
   integer,allocatable :: c0_to_c0b(:)      !< Subset Sc0, global to halo B
   integer,allocatable :: c0a_to_c0b(:)     !< Subset Sc0, halo A to halo B

   ! Required data to apply an observation operator

   ! Number of points
   integer :: nc0a                          !< Halo A size
   integer :: nc0b                          !< Halo B size

   ! Number of observations
   integer :: nobsa                         !< Local number of observations

   ! Interpolation data
   type(linoptype) :: h                     !< Interpolation data

   ! Communication data
   type(comtype) :: com                     !< Communication data
end type odatatype

private
public :: odatatype,yobs_com_gl,yobs_com_lg

contains

!----------------------------------------------------------------------
! Subroutine: yobs_com_gl
!> Purpose: communicate full field from global to local distribution
!----------------------------------------------------------------------
subroutine yobs_com_gl(odata,yobs)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata                    !< Observation operator data
real(kind_real),allocatable,intent(inout) :: yobs(:,:) !< Observations

! Local variables
integer :: iobs,iobsa,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Associate
associate(geom=>odata%geom)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      allocate(sbuf(odata%proc_to_nobsa(iproc)*geom%nl0))

      ! Initialization
      call msr(sbuf)

      ! Prepare buffer
      do iobs=1,odata%nobs
         jproc = odata%obs_to_proc(iobs)
         if (jproc==iproc) then
            iobsa = odata%obs_to_obsa(iobs)
            sbuf((iobsa-1)*geom%nl0+1:iobsa*geom%nl0) = yobs(iobs,:)
         end if
      end do

      if (iproc==mpl%ioproc) then
         ! Allocation
         allocate(rbuf(odata%proc_to_nobsa(iproc)*geom%nl0))

         ! Copy data
         rbuf = sbuf
      else
         ! Send data to iproc
         call mpl_send(odata%proc_to_nobsa(iproc)*geom%nl0,sbuf,iproc,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Allocation
   allocate(rbuf(odata%proc_to_nobsa(mpl%myproc)*geom%nl0))

   ! Receive data from ioproc
   call mpl_recv(odata%proc_to_nobsa(mpl%myproc)*geom%nl0,rbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Rellocation
if (allocated(yobs)) deallocate(yobs)
allocate(yobs(odata%proc_to_nobsa(mpl%myproc),geom%nl0))

! Copy from buffer
do iobsa=1,odata%proc_to_nobsa(mpl%myproc)
   yobs(iobsa,:) = rbuf((iobsa-1)*geom%nl0+1:iobsa*geom%nl0)
end do

! End associate
end associate

end subroutine yobs_com_gl

!----------------------------------------------------------------------
! Subroutine: yobs_com_lg
!> Purpose: communicate full field from local to global distribution
!----------------------------------------------------------------------
subroutine yobs_com_lg(odata,yobs)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata                    !< Observation operator data
real(kind_real),allocatable,intent(inout) :: yobs(:,:) !< Observations

! Local variables
integer :: iobs,iobsa,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Associate
associate(geom=>odata%geom)

! Allocation
allocate(sbuf(odata%proc_to_nobsa(mpl%myproc)*geom%nl0))

! Prepare buffer
do iobsa=1,odata%proc_to_nobsa(mpl%myproc)
   sbuf((iobsa-1)*geom%nl0+1:iobsa*geom%nl0) = yobs(iobsa,:)
end do

! Release memory
deallocate(yobs)

! Communication
if (mpl%main) then
   ! Allocation
   allocate(yobs(odata%nobs,geom%nl0))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(rbuf(odata%proc_to_nobsa(iproc)*geom%nl0))

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Receive data from iproc
         call mpl_recv(odata%proc_to_nobsa(iproc)*geom%nl0,rbuf,iproc,mpl%tag)
      end if

      ! Copy from buffer
      do iobs=1,odata%nobs
         jproc = odata%obs_to_proc(iobs)
         if (jproc==iproc) then
            iobsa = odata%obs_to_obsa(iobs)
            yobs(iobs,:) = rbuf((iobsa-1)*geom%nl0+1:iobsa*geom%nl0)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Sending data to iproc
   call mpl_send(odata%proc_to_nobsa(mpl%myproc)*geom%nl0,sbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! End associate
end associate

end subroutine yobs_com_lg

end module type_odata
