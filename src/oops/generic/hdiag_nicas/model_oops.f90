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

use netcdf
use tools_const, only: pi,rad2deg,req
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msi,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use tools_qsort, only: qsort
use type_geom, only: geomtype,geom_alloc
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast,mpl_allreduce_sum
use type_nam, only: namtype

implicit none

private
public :: model_oops_coord,model_oops_write

contains

!----------------------------------------------------------------------
! Subroutine: model_oops_coord
!> Purpose: load OOPS coordinates
!----------------------------------------------------------------------
subroutine model_oops_coord(nam,geom,lats,lons,areas,vunit,imask,glbind)

implicit none

! Passed variables
type(namtype),intent(inout) :: nam               !< Namelist
type(geomtype),intent(inout) :: geom             !< Geometry
real(kind_real),intent(in) :: lats(geom%nc0a)    !< Latitudes
real(kind_real),intent(in) :: lons(geom%nc0a)    !< Longitudes
real(kind_real),intent(in) :: areas(geom%nc0a)   !< Areas
real(kind_real),intent(in) :: vunit(geom%nlev)   !< Vertical unit
integer,intent(in) :: imask(geom%nc0a*geom%nlev) !< Mask
integer,intent(in) :: glbind(geom%nc0a)          !< Global index

! Local variables
integer :: nc0ag(mpl%nproc),ic0,ic0a,il0,offset,iproc
integer,allocatable :: glbindg(:),order(:)
logical,allocatable :: lmask(:,:)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         nc0ag(iproc) = geom%nc0a
      else
         ! Receive data on ioproc
         call mpl_recv(nc0ag(iproc),iproc,mpl%tag)
      end if
   end do
else
   ! Send data to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl_bcast(nc0ag,mpl%ioproc)

! Number of levels
geom%nl0 = geom%nlev

! Global number of nodes
geom%nc0 = sum(nc0ag)

! Print summary
write(mpl%unit,'(a7,a)') '','Distribution summary:'
do iproc=1,mpl%nproc
   write(mpl%unit,'(a10,a,i3,a,i8,a)') '','Proc #',iproc,': ',nc0ag(iproc),' grid-points'
end do
write(mpl%unit,'(a10,a,i8,a)') '','Total: ',geom%nc0,' grid-points'

! Allocation
call geom_alloc(geom)
allocate(geom%ic0_to_iproc(geom%nc0))
allocate(geom%ic0_to_ic0a(geom%nc0))
allocate(lmask(geom%nc0a,geom%nl0))
allocate(glbindg(geom%nc0))

! Define local index and MPI task
ic0 = 0
do iproc=1,mpl%nproc
   do ic0a=1,nc0ag(iproc)
      ic0 = ic0+1
      geom%ic0_to_iproc(ic0) = iproc
      geom%ic0_to_ic0a(ic0) = ic0a
   end do
end do

! Convert mask
do il0=1,geom%nl0
   offset = (nam%levs(il0)-1)*geom%nc0a
   do ic0a=1,geom%nc0a
      if (imask(offset+ic0a)==0) then
         lmask(ic0a,il0) = .false.
      elseif (imask(offset+ic0a)==1) then
         lmask(ic0a,il0) = .true.
      else
         call msgerror('wrong mask value in model_oops_coord')
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
         geom%lon(offset+1:offset+nc0ag(iproc)) = lons
         geom%lat(offset+1:offset+nc0ag(iproc)) = lats
         do il0=1,geom%nl0
            geom%mask(offset+1:offset+nc0ag(iproc),il0) = lmask(:,il0)
         end do
         glbindg(offset+1:offset+nc0ag(iproc)) = glbind
      else
         ! Receive data on ioproc
         call mpl_recv(nc0ag(iproc),geom%lon(offset+1:offset+nc0ag(iproc)),iproc,mpl%tag)
         call mpl_recv(nc0ag(iproc),geom%lat(offset+1:offset+nc0ag(iproc)),iproc,mpl%tag+1)
         do il0=1,geom%nl0
            call mpl_recv(nc0ag(iproc),geom%mask(offset+1:offset+nc0ag(iproc),il0),iproc,mpl%tag+1+il0)
         end do
         call mpl_recv(nc0ag(iproc),glbindg(offset+1:offset+nc0ag(iproc)),iproc,mpl%tag+2+geom%nl0)
      end if

      !  Update offset
      offset = offset+nc0ag(iproc)
   end do
else
   ! Send data to ioproc
   call mpl_send(geom%nc0a,lons,mpl%ioproc,mpl%tag)
   call mpl_send(geom%nc0a,lats,mpl%ioproc,mpl%tag+1)
   do il0=1,geom%nl0
      call mpl_send(geom%nc0a,lmask(:,il0),mpl%ioproc,mpl%tag+1+il0)
   end do
   call mpl_send(geom%nc0a,glbind,mpl%ioproc,mpl%tag+2+geom%nl0)
end if
mpl%tag = mpl%tag+3+geom%nl0

! Broadcast data
call mpl_bcast(geom%lon,mpl%ioproc)
call mpl_bcast(geom%lat,mpl%ioproc)
call mpl_bcast(geom%mask,mpl%ioproc)
call mpl_bcast(glbindg,mpl%ioproc)

! Reorder data
if (all(glbindg>0)) then
   ! Allocation
   allocate(order(geom%nc0))

   ! Sort glbindg
   call qsort(geom%nc0,glbindg,order)

   ! Check glbindg
   do ic0=2,geom%nc0
      if (glbindg(ic0)<=glbindg(ic0-1)) call msgerror('wrong glbindg in model_oops')
   end do

   ! Reorder columns
   geom%lon = geom%lon(order)
   geom%lat = geom%lat(order)
   do il0=1,geom%nl0
      geom%mask(:,il0) = geom%mask(order,il0)
   end do
   geom%ic0_to_iproc = geom%ic0_to_iproc(order)
   geom%ic0_to_ic0a = geom%ic0_to_ic0a(order)
end if

! Normalized area
do il0=1,geom%nl0
   call mpl_allreduce_sum(sum(areas,mask=lmask(:,il0))/req**2,geom%area(il0))
end do

! Vertical unit
geom%vunit = vunit

end subroutine model_oops_coord

!----------------------------------------------------------------------
! Subroutine: model_oops_write
!> Purpose: write OOPS field
!----------------------------------------------------------------------
subroutine model_oops_write(geom,ncid,varname,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                    !< Geometry
integer,intent(in) :: ncid                           !< NetCDF file ID
character(len=*),intent(in) :: varname               !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0,geom%nl0) !< Written field

! Local variables
integer :: il0,ierr
integer :: nc0_id,nlev_id,fld_id,lon_id,lat_id
character(len=1024) :: subr = 'model_oops_write'

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'nc0',nc0_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
   ierr = nf90_inq_dimid(ncid,'nlev',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nlev',geom%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlev_id,nc0_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,geom%nl0
   if (isanynotmsr(fld(:,il0))) then
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld(:,il0),(/il0,1/),(/1,geom%nc0/)))
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
   call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
end if

end subroutine model_oops_write

end module model_oops
