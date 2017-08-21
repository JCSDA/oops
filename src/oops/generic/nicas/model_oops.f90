!----------------------------------------------------------------------
! Module: model_oops.f90
!> Purpose: OOPS required routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_oops

use module_namelist, only: nam
use netcdf
use tools_const, only: pi,deg2rad,rad2deg,req
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msi,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast,mpl_allreduce_sum
use type_ndata, only: ndatatype,ndata_alloc

implicit none

private
public :: model_oops_coord,model_oops_write

contains

!----------------------------------------------------------------------
! Subroutine: model_oops_coord
!> Purpose: load OOPS coordinates
!----------------------------------------------------------------------
subroutine model_oops_coord(lats,lons,areas,levs,mask3d,mask2d,glbind,ndata)

implicit none

! Passed variables
real(kind_real),intent(in) :: lats(:)
real(kind_real),intent(in) :: lons(:)
real(kind_real),intent(in) :: areas(:)
real(kind_real),intent(in) :: levs(:)
integer,intent(in) :: mask3d(:)
integer,intent(in) :: mask2d(:)
integer,intent(in) :: glbind(:)
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: nc0a,nc0ag(mpl%nproc),ic0,ic0a,il0,offset,iproc
integer,allocatable :: glbindg(:),order(:)
logical,allocatable :: lmask3d(:,:)

! TODO: change that one day
ndata%nl0 = nam%nl

! Local number of nodes and levels
nc0a = size(lats)
ndata%nlev = size(levs)

! Check TODO: check dimensions consistency of all arrays
if (any(nam%levs>ndata%nlev)) call msgerror('not enough levels in model_oops')

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         nc0ag(iproc) = nc0a
      else
         ! Receive data on ioproc
         call mpl_recv(nc0ag(iproc),iproc,mpl%tag)
      end if
   end do
else
   ! Send data to ioproc
   call mpl_send(nc0a,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl_bcast(nc0ag,mpl%ioproc)

! Global number of nodes
ndata%nc0 = sum(nc0ag)

! Allocation
call ndata_alloc(ndata)
allocate(ndata%ic0_to_iproc(ndata%nc0))
allocate(ndata%ic0_to_ic0a(ndata%nc0))
allocate(lmask3d(nc0a,ndata%nl0))
allocate(glbindg(ndata%nc0))

! Define local index and MPI task
ndata%nc0amax = maxval(nc0ag)
ic0 = 0
do iproc=1,mpl%nproc
   do ic0a=1,nc0ag(iproc)
      ic0 = ic0+1
      ndata%ic0_to_iproc(ic0) = iproc
      ndata%ic0_to_ic0a(ic0) = ic0a
   end do
end do

! Convert 3d mask
do il0=1,ndata%nl0
   offset = (nam%levs(il0)-1)*nc0a
   do ic0a=1,nc0a
      if (mask3d(offset+ic0a)==0) then
         lmask3d(ic0a,il0) = .false.
      elseif (mask3d(offset+ic0a)==1) then
         lmask3d(ic0a,il0) = .true.
      else
         call msgerror('wrong 3d mask value in model_oops_coord')
      end if
   end do
end do

! Communication and reordering
if (mpl%main) then
   ! Allocation
   offset = 0
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         ndata%lon(offset+1:offset+nc0ag(iproc)) = lons*deg2rad
         ndata%lat(offset+1:offset+nc0ag(iproc)) = lats*deg2rad
         do il0=1,ndata%nl0
            ndata%mask(offset+1:offset+nc0ag(iproc),il0) = lmask3d(:,il0)
         end do
         glbindg(offset+1:offset+nc0ag(iproc)) = glbind
      else
         ! Receive data on ioproc
         call mpl_recv(nc0ag(iproc),ndata%lon(offset+1:offset+nc0ag(iproc)),iproc,mpl%tag)
         call mpl_recv(nc0ag(iproc),ndata%lat(offset+1:offset+nc0ag(iproc)),iproc,mpl%tag+1)
         do il0=1,ndata%nl0
            call mpl_recv(nc0ag(iproc),ndata%mask(offset+1:offset+nc0ag(iproc),il0),iproc,mpl%tag+1+il0)
         end do
         call mpl_recv(nc0ag(iproc),glbindg(offset+1:offset+nc0ag(iproc)),iproc,mpl%tag+2+ndata%nl0)
      end if

      !  Update offset
      offset = offset+nc0ag(iproc)
   end do
else
   ! Send data to ioproc
   call mpl_send(nc0a,lons*deg2rad,mpl%ioproc,mpl%tag)
   call mpl_send(nc0a,lats*deg2rad,mpl%ioproc,mpl%tag+1)
   do il0=1,ndata%nl0
      call mpl_send(nc0a,lmask3d(:,il0),mpl%ioproc,mpl%tag+1+il0)
   end do
   call mpl_send(nc0a,glbind,mpl%ioproc,mpl%tag+2+ndata%nl0)
end if
mpl%tag = mpl%tag+3+ndata%nl0

! Broadcast data
call mpl_bcast(ndata%lon,mpl%ioproc)
call mpl_bcast(ndata%lat,mpl%ioproc)
call mpl_bcast(ndata%mask,mpl%ioproc)
call mpl_bcast(glbindg,mpl%ioproc)

! Reorder data
if (all(glbindg>0)) then
   ! Allocation
   allocate(order(ndata%nc0))

   ! Sort glbindg
   call qsort(ndata%nc0,glbindg,order)

   ! Check glbindg
   do ic0=2,ndata%nc0
      if (glbindg(ic0)<=glbindg(ic0-1)) call msgerror('wrong glbindg in model_oops')
   end do

   ! Reorder columns
   ndata%lon = ndata%lon(order)
   ndata%lat = ndata%lat(order)
   do il0=1,ndata%nl0
      ndata%mask(:,il0) = ndata%mask(order,il0)
   end do
   ndata%ic0_to_iproc = ndata%ic0_to_iproc(order)
   ndata%ic0_to_ic0a = ndata%ic0_to_ic0a(order)
end if

! Normalized area
do il0=1,ndata%nl0
   call mpl_allreduce_sum(sum(areas,mask=lmask3d(:,il0))/req**2,ndata%area(il0))
end do

! Vertical unit
ndata%vunit = levs

end subroutine model_oops_coord

!----------------------------------------------------------------------
! Subroutine: model_oops_write
!> Purpose: write OOPS field
!----------------------------------------------------------------------
subroutine model_oops_write(ncid,varname,ndata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                             !< NetCDF file ID
character(len=*),intent(in) :: varname                 !< Variable name
type(ndatatype),intent(in) :: ndata                    !< Sampling data
real(kind_real),intent(in) :: fld(ndata%nc0,ndata%nl0) !< Written field

! Local variables
integer :: il0,ierr
integer :: nc0_id,nlev_id,nt_id,fld_id,lon_id,lat_id
character(len=1024) :: subr = 'model_oops_write'

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'nc0',nc0_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nc0',ndata%nc0,nc0_id))
   ierr = nf90_inq_dimid(ncid,'nlev',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nlev',ndata%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlev_id,nc0_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,ndata%nl0
   if (isanynotmsr(fld(:,il0))) then
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld(:,il0),(/il0,1/),(/1,ndata%nc0/)))
   end if
end do

! Write coordinates
ierr = nf90_inq_varid(ncid,'lon',lon_id)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
   call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
   call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
   call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
   call ncerr(subr,nf90_enddef(ncid))
   call ncerr(subr,nf90_put_var(ncid,lon_id,ndata%lon*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,lat_id,ndata%lat*rad2deg))
end if

end subroutine model_oops_write

end module model_oops
