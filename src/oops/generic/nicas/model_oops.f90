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
use tools_const, only: pi,deg2rad,rad2deg
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msi,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_esmf, only: esmf_create_field
use type_mpl, only: mpl
use type_sdata, only: sdatatype,sdata_alloc
implicit none

private
public :: model_oops_coord,model_oops_write

contains

!----------------------------------------------------------------------
! Subroutine: model_oops_coord
!> Purpose: load OOPS coordinates
!----------------------------------------------------------------------
subroutine model_oops_coord(dims,lats,lons,levs,area,mask,sdata)

implicit none

! Passed variables
integer,intent(in) :: dims(:)
real(kind_real),intent(in) :: lats(:)
real(kind_real),intent(in) :: lons(:)
real(kind_real),intent(in) :: levs(:)
real(kind_real),intent(in) :: area(:)
integer,intent(in) :: mask(:)
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ic0,il0,offset

! TODO: change that one day
sdata%nl0 = nam%nl

! Number of nodes
sdata%nc0 = size(lats)
sdata%nlev = size(levs)
if (sdata%nc0/=product(dims)) call msgerror ('product(dims) should be equal to nc0')

write(mpl%unit,*) 'in nicas'
write(mpl%unit,*) sdata%nc0,sdata%nlev,sdata%nl0
write(mpl%unit,*) 

! Normalized area
allocate(sdata%area(sdata%nl0))
sdata%area = area(nam%levs)

! Pack
call sdata_alloc(sdata)
sdata%lon = lons*deg2rad
sdata%lat = lats*deg2rad
do il0=1,sdata%nl0
   offset = (nam%levs(il0)-1)*sdata%nc0
   write(mpl%unit,*) 'test',il0,offset,offset+1,offset+sdata%nc0,size(mask)
   do ic0=1,sdata%nc0
      if (mask(offset+ic0)==0) then
         sdata%mask(ic0,il0) = .false.
      elseif (mask(offset+ic0)==1) then
         sdata%mask(ic0,il0) = .true.
      else
         call msgerror('wrong mask value in model_oops_coord')
      end if
   end do
end do

! Vertical unit
sdata%vunit = levs

! ESMF field for interpolations
if (size(dims)==1) then
   ! Create ESMF field from mesh
   call esmf_create_field(sdata,sdata%nc0,sdata%lon,sdata%lat,any(sdata%mask,dim=2),sdata%c0field)
elseif (size(dims)==2) then
   ! Create ESMF field from grid
   sdata%nlon = dims(1)
   sdata%nlat = dims(2)
   call esmf_create_field(sdata)
else
   ! Not yet implemented
   call msgerror('wrong dims size in create_nicas')
end if

end subroutine model_oops_coord

!----------------------------------------------------------------------
! Subroutine: model_oops_write
!> Purpose: write OOPS field
!----------------------------------------------------------------------
subroutine model_oops_write(ncid,varname,sdata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                             !< NetCDF file ID
character(len=*),intent(in) :: varname                 !< Variable name
type(sdatatype),intent(in) :: sdata                    !< Sampling data
real(kind_real),intent(in) :: fld(sdata%nc0,sdata%nl0) !< Written field

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
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nc0',sdata%nc0,nc0_id))
   ierr = nf90_inq_dimid(ncid,'nlev',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nlev',sdata%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlev_id,nc0_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,sdata%nl0
   if (isanynotmsr(fld(:,il0))) then
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld(:,il0),(/il0,1/),(/1,sdata%nc0/)))
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
   call ncerr(subr,nf90_put_var(ncid,lon_id,sdata%lon*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,lat_id,sdata%lat*rad2deg))
end if

end subroutine model_oops_write

end module model_oops
