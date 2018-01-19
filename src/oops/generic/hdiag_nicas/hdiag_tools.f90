!----------------------------------------------------------------------
! Module: hdiag_tools.f90
!> Purpose: diagnostics tools routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_tools

use netcdf
use omp_lib
use tools_const, only: req,sphere_dist,rad2deg,gc99
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isnotmsi,isnotmsr,isanynotmsr,isallnotmsr
use tools_nc, only: ncfloat,ncerr
use tools_qsort, only: qsort
use type_hdata, only: hdatatype
use type_linop, only: apply_linop
use type_mpl, only: mpl,mpl_send,mpl_recv
implicit none

interface diag_com_lg
  module procedure diag_com_lg_single
  module procedure diag_com_lg
end interface

private
public :: diag_filter,diag_interpolation,diag_com_lg

contains

!----------------------------------------------------------------------
! Subroutine: diag_filter
!> Purpose: filter diagnostics
!----------------------------------------------------------------------
subroutine diag_filter(hdata,il0,filter_type,r,diag)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata              !< HDIAG data
integer,intent(in) :: il0                        !< Level
character(len=*),intent(in) :: filter_type       !< Filter type
real(kind_real),intent(in) :: r                  !< Filter support radius
real(kind_real),intent(inout) :: diag(hdata%nc2) !< Filtered diagnostics

! Local variables
integer :: ic2,jc2,nc2eff
integer,allocatable :: order(:)
real(kind_real) :: diag_tmp(hdata%nc2),distnorm,norm,wgt
real(kind_real),allocatable :: diag_eff(:),diag_eff_dist(:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Copy diagnostics
diag_tmp = diag

!$omp parallel do schedule(static) private(ic2,diag_eff,diag_eff_dist,nc2eff,jc2,distnorm,norm,order,wgt)
do ic2=1,hdata%nc2
   if (hdata%c1l0_log(hdata%c2_to_c1(ic2),il0)) then
      ! Allocation
      allocate(diag_eff(hdata%nc2))
      allocate(diag_eff_dist(hdata%nc2))

      ! Build diag_eff of valid points
      nc2eff = 0
      jc2 = 1
      do while (hdata%nn_c2_dist(jc2,ic2,min(il0,geom%nl0i))<r)
         ! Check the point validity
         if (isnotmsr(diag_tmp(hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i))))) then
            nc2eff = nc2eff+1
            diag_eff(nc2eff) = diag_tmp(hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i)))
            diag_eff_dist(nc2eff) = hdata%nn_c2_dist(jc2,ic2,min(il0,geom%nl0i))
         end if
         jc2 = jc2+1
         if (jc2>hdata%nc2) exit
      end do

      ! Apply filter
      if (nc2eff>0) then
         select case (trim(filter_type))
         case ('average')
            ! Compute average
            diag(ic2) = sum(diag_eff(1:nc2eff))/float(nc2eff)
         case ('gc99')
            ! Gaspari-Cohn (1999) kernel
            diag(ic2) = 0.0
            norm = 0.0
            do jc2=1,nc2eff
               distnorm = diag_eff_dist(jc2)/r
               wgt = gc99(distnorm)
               diag(ic2) = diag(ic2)+wgt*diag_eff(jc2)
               norm = norm+wgt
            end do
            diag(ic2) = diag(ic2)/norm
         case ('median')
            ! Compute median
            allocate(order(nc2eff))
            call qsort(nc2eff,diag_eff(1:nc2eff),order)
            if (mod(nc2eff,2)==0) then
               diag(ic2) = 0.5*(diag_eff(nc2eff/2)+diag_eff(nc2eff/2+1))
            else
               diag(ic2) = diag_eff((nc2eff+1)/2)
            end if
            deallocate(order)
         case default
            ! Wrong filter
            call msgerror('wrong filter type')
         end select
      else
         call msr(diag(ic2))
      end if

      ! Release memory
      deallocate(diag_eff)
      deallocate(diag_eff_dist)
   end if
end do
!$omp end parallel do

! End associate
end associate

end subroutine diag_filter

!----------------------------------------------------------------------
! Subroutine: diag_interpolation
!> Purpose: compute interpolation from diagnostic points to grid-points, on multiple levels
!----------------------------------------------------------------------
subroutine diag_interpolation(hdata,fld_c2,fld)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                               !< HDIAG data
real(kind_real),intent(in) :: fld_c2(hdata%nc2,hdata%geom%nl0)    !< Diagnostic
real(kind_real),intent(out) :: fld(hdata%geom%nc0,hdata%geom%nl0) !< Interpolated diagnostic

! Local variables
integer :: il0

! Associate
associate(geom=>hdata%geom)

do il0=1,geom%nl0
   ! Compute interpolation on a single level
   call apply_linop(hdata%h(min(il0,geom%nl0i)),fld_c2(:,il0),fld(:,il0))
end do

! End associate
end associate

end subroutine diag_interpolation

!----------------------------------------------------------------------
! Subroutine: diag_com_lg_single
!> Purpose: communicate diagnostic from local to global distribution, single level
!----------------------------------------------------------------------
subroutine diag_com_lg_single(hdata,diag)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                  !< HDIAG data
real(kind_real),allocatable,intent(inout) :: diag(:) !< Diagnostic

! Local variables
integer :: ic2,ic2a,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Allocation
allocate(sbuf(hdata%nc2a))

! Prepare buffer
sbuf = diag

! Release memory
deallocate(diag)

! Communication
if (mpl%main) then
   ! Allocation
   allocate(diag(hdata%nc2))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(rbuf(hdata%proc_to_nc2a(iproc)))

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Receive data from iproc
         call mpl_recv(hdata%proc_to_nc2a(iproc),rbuf,iproc,mpl%tag)
      end if

      ! Copy from buffer
      do ic2=1,hdata%nc2
         jproc = hdata%c2_to_proc(ic2)
         if (jproc==iproc) then
            ic2a = hdata%c2_to_c2a(ic2)
            diag(ic2) = rbuf(ic2a)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Sending data to iproc
   call mpl_send(hdata%nc2a,sbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! End associate
end associate

end subroutine diag_com_lg_single

!----------------------------------------------------------------------
! Subroutine: diag_com_lg
!> Purpose: communicate diagnostic from local to global distribution
!----------------------------------------------------------------------
subroutine diag_com_lg(hdata,diag)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                    !< HDIAG data
real(kind_real),allocatable,intent(inout) :: diag(:,:) !< Diagnostic

! Local variables
integer :: ic2,ic2a,iproc,jproc
real(kind_real),allocatable :: diag_loc(:,:),sbuf(:),rbuf(:)
logical,allocatable :: mask_unpack(:,:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Allocation
allocate(sbuf(hdata%nc2a*geom%nl0))

! Prepare buffer
sbuf = pack(diag,mask=.true.)

! Release memory
deallocate(diag)

! Communication
if (mpl%main) then
   ! Allocation
   allocate(diag(hdata%nc2,geom%nl0))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(diag_loc(hdata%proc_to_nc2a(iproc),geom%nl0))
      allocate(mask_unpack(hdata%proc_to_nc2a(iproc),geom%nl0))
      allocate(rbuf(hdata%proc_to_nc2a(iproc)*geom%nl0))

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Receive data from iproc
         call mpl_recv(hdata%proc_to_nc2a(iproc)*geom%nl0,rbuf,iproc,mpl%tag)
      end if

      ! Copy from buffer
      mask_unpack = .true.
      diag_loc = unpack(rbuf,mask=mask_unpack,field=diag_loc)
      do ic2=1,hdata%nc2
         jproc = hdata%c2_to_proc(ic2)
         if (jproc==iproc) then
            ic2a = hdata%c2_to_c2a(ic2)
            diag(ic2,:) = diag_loc(ic2a,:)
         end if
      end do

      ! Release memory
      deallocate(diag_loc)
      deallocate(mask_unpack)
      deallocate(rbuf)
   end do
else
   ! Sending data to iproc
   call mpl_send(hdata%nc2a*geom%nl0,sbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! End associate
end associate

end subroutine diag_com_lg

end module hdiag_tools
