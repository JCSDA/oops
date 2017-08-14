!----------------------------------------------------------------------
! Module: module_gem.f90
!> Purpose: GEM model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_gem

use module_namelist, only: nam
use netcdf
use tools_const, only: pi,req,deg2rad,ps
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_esmf, only: esmf_create_field
use type_ndata, only: ndatatype,ndata_alloc

implicit none

private
public :: model_gem_coord,model_gem_read,model_gem_write

contains

!----------------------------------------------------------------------
! Subroutine: model_gem_coord
!> Purpose: get GEM coordinates
!----------------------------------------------------------------------
subroutine model_gem_coord(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,a_id,b_id
real(kind=8),allocatable :: lon(:,:),lat(:,:),a(:),b(:)
character(len=1024) :: subr = 'model_gem_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'lon',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'lat',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=ndata%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=ndata%nlat))
ndata%nc0 = ndata%nlon*ndata%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'lev',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=ndata%nlev))

! Allocation
allocate(lon(ndata%nlon,ndata%nlat))
allocate(lat(ndata%nlon,ndata%nlat))
allocate(a(ndata%nlev))
allocate(b(ndata%nlev))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'lon',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'lat',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'ap',a_id))
call ncerr(subr,nf90_inq_varid(ncid,'b',b_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon(:,1)))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat(1,:)))
call ncerr(subr,nf90_get_var(ncid,a_id,a))
call ncerr(subr,nf90_get_var(ncid,b_id,b))
call ncerr(subr,nf90_close(ncid))

! Compute normalized area
allocate(ndata%area(ndata%nl0))
ndata%area = 4.0*pi

! Convert to radian
lon(:,1) = lon(:,1)*deg2rad
lat(1,:) = lat(1,:)*deg2rad

! Fill arrays
do ilat=1,ndata%nlat
   lon(:,ilat) = lon(:,1)
end do
do ilon=1,ndata%nlon
   lat(ilon,:) = lat(1,:)
end do

! Pack
call ndata_alloc(ndata)
ndata%lon = pack(real(lon,kind_real),mask=.true.)
ndata%lat = pack(real(lat,kind_real),mask=.true.)
ndata%mask = .true.

! Vertical unit
if (nam%logpres) then
   ndata%vunit = log(a(nam%levs(1:ndata%nl0))+b(nam%levs(1:ndata%nl0))*ps)
else
   ndata%vunit = float(nam%levs(1:ndata%nl0))
end if

! Create ESMF field from grid
call esmf_create_field(ndata)

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(a)
deallocate(b)

end subroutine model_gem_coord

!----------------------------------------------------------------------
! Subroutine: model_gem_read
!> Purpose: read GEM field
!----------------------------------------------------------------------
subroutine model_gem_read(ncid,varname,ndata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                              !< NetCDF file ID
character(len=*),intent(in) :: varname                  !< Variable name
type(ndatatype),intent(in) :: ndata                     !< Sampling data
real(kind_real),intent(out) :: fld(ndata%nc0,ndata%nl0) !< Read field

! Local variables
integer :: il0,xt
integer :: fld_id
integer,allocatable :: fld_loc_int(:,:)
real(kind_real) :: add_offset,scale_factor
real(kind_real),allocatable :: fld_loc(:,:)
character(len=1024) :: subr = 'model_gem_read'

! Initialize field
call msr(fld)

! Get variable id
call ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

! Check variable type
call ncerr(subr,nf90_inquire_variable(ncid,fld_id,xtype=xt))
if (xt==nf90_short) then
   allocate(fld_loc_int(ndata%nlon,ndata%nlat))
elseif (xt==nf90_double) then
   allocate(fld_loc(ndata%nlon,ndata%nlat))
else
   call msgerror('wrong variable type')
end if

! Read variable
do il0=1,ndata%nl0
   if (xt==nf90_short) then
      call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc_int,(/1,1,nam%levs(il0)/),(/ndata%nlon,ndata%nlat,1/)))
      call ncerr(subr,nf90_get_att(ncid,fld_id,'add_offset',add_offset))
      call ncerr(subr,nf90_get_att(ncid,fld_id,'scale_factor',scale_factor))
      fld(:,il0) = pack(add_offset+scale_factor*float(fld_loc_int),mask=.true.)
   elseif (xt==nf90_double) then
      call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%levs(il0)/),(/ndata%nlon,ndata%nlat,1/)))
      fld(:,il0) = pack(fld_loc,mask=.true.)
   else
      call msgerror('wrong netcdf variable type')
   end if
end do

end subroutine model_gem_read

!----------------------------------------------------------------------
! Subroutine: model_gem_write
!> Purpose: write GEM field
!----------------------------------------------------------------------
subroutine model_gem_write(ncid,varname,ndata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                             !< NetCDF file ID
character(len=*),intent(in) :: varname                 !< Variable name
type(ndatatype),intent(in) :: ndata                    !< Sampling data
real(kind_real),intent(in) :: fld(ndata%nc0,ndata%nl0) !< Written field

! Local variables
integer :: il0,ierr
integer :: nlon_id,nlat_id,nlev_id,fld_id
real(kind_real) :: fld_loc(ndata%nlon,ndata%nlat)
logical :: mask_unpack(ndata%nlon,ndata%nlat)
character(len=1024) :: subr = 'model_gem_write'

! Initialization
mask_unpack = .true.

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'lon',nlon_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'lon',ndata%nlon,nlon_id))
   ierr = nf90_inq_dimid(ncid,'lat',nlat_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'lat',ndata%nlat,nlat_id))
   ierr = nf90_inq_dimid(ncid,'lev',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'lev',ndata%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,ndata%nl0
   if (isanynotmsr(fld(:,il0))) then
      call msr(fld_loc)
      fld_loc = unpack(fld(:,il0),mask=mask_unpack,field=fld_loc)
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/ndata%nlon,ndata%nlat,1/)))
   end if
end do

end subroutine model_gem_write

end module model_gem
