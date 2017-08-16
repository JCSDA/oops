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
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndata_alloc

implicit none

private
public :: model_oops_coord,model_oops_write

contains

!----------------------------------------------------------------------
! Subroutine: model_oops_coord
!> Purpose: load OOPS coordinates
!----------------------------------------------------------------------
subroutine model_oops_coord(lats,lons,levs,mask,ndata)

implicit none

! Passed variables
real(kind_real),intent(in) :: lats(:)
real(kind_real),intent(in) :: lons(:)
real(kind_real),intent(in) :: levs(:)
integer,intent(in) :: mask(:)
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ic0,il0,offset

! TODO: change that one day
ndata%nl0 = nam%nl

! Number of nodes
ndata%nc0 = size(lats)
ndata%nlev = size(levs)

! Pack
call ndata_alloc(ndata)
ndata%lon = lons*deg2rad
ndata%lat = lats*deg2rad
do il0=1,ndata%nl0
   offset = (nam%levs(il0)-1)*ndata%nc0
   do ic0=1,ndata%nc0
      if (mask(offset+ic0)==0) then
         ndata%mask(ic0,il0) = .false.
      elseif (mask(offset+ic0)==1) then
         ndata%mask(ic0,il0) = .true.
      else
         call msgerror('wrong mask value in model_oops_coord')
      end if
   end do
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
