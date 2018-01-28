!----------------------------------------------------------------------
! Module: module_aro.f90
!> Purpose: AROME model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_aro

use netcdf
use tools_const, only: deg2rad,rad2deg,req,ps
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geomtype,geom_alloc
use type_nam, only: namtype

implicit none

private
public :: model_aro_coord,model_aro_read,model_aro_write

contains

!----------------------------------------------------------------------
! Subroutine: model_aro_coord
!> Purpose: load AROME coordinates
!----------------------------------------------------------------------
subroutine model_aro_coord(nam,geom)

implicit none

! Passed variables
type(namtype),intent(in) :: nam      !< Namelist
type(geomtype),intent(inout) :: geom !< Geometry

! Local variables
integer :: ncid,nlon_id,nlat_id,nlev_id,pp_id,lon_id,lat_id,cmask_id,a_id,b_id
real(kind_real) :: dx,dy
real(kind=8),allocatable :: lon(:,:),lat(:,:),cmask(:,:),a(:),b(:)
character(len=1024) :: subr = 'model_aro_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'X',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'Y',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
geom%nc0 = geom%nlon*geom%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'Z',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon,geom%nlat))
allocate(lat(geom%nlon,geom%nlat))
allocate(cmask(geom%nlon,geom%nlat))
allocate(a(geom%nlev+1))
allocate(b(geom%nlev+1))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'longitude',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'latitude',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'cmask',cmask_id))
call ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_A',a_id))
call ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_B',b_id))
call ncerr(subr,nf90_inq_varid(ncid,'Projection_parameters',pp_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call ncerr(subr,nf90_get_var(ncid,cmask_id,cmask))
call ncerr(subr,nf90_get_var(ncid,a_id,a))
call ncerr(subr,nf90_get_var(ncid,b_id,b))
call ncerr(subr,nf90_get_att(ncid,pp_id,'x_resolution',dx))
call ncerr(subr,nf90_get_att(ncid,pp_id,'y_resolution',dy))
call ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*real(deg2rad,kind=8)
lat = lat*real(deg2rad,kind=8)

! Pack
call geom_alloc(geom)
geom%lon = pack(real(lon,kind_real),mask=.true.)
geom%lat = pack(real(lat,kind_real),mask=.true.)
geom%mask = .true.

! Compute normalized area
geom%area = float(geom%nlon*geom%nlat)*dx*dy/req**2

! Vertical unit
if (nam%logpres) then
   geom%vunit(1:nam%nl) = log(0.5*(a(1:nam%nl)+a(2:nam%nl+1))+0.5*(b(1:nam%nl)+b(2:nam%nl+1))*ps)
   if (geom%nl0>nam%nl) geom%vunit(geom%nl0) = log(ps)
else
   geom%vunit = float(nam%levs(1:geom%nl0))
end if

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(cmask)
deallocate(a)
deallocate(b)

end subroutine model_aro_coord

!----------------------------------------------------------------------
! Subroutine: model_aro_read
!> Purpose: read AROME field
!----------------------------------------------------------------------
subroutine model_aro_read(nam,geom,ncid,its,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                              !< Namelist
type(geomtype),intent(in) :: geom                            !< Geometry
integer,intent(in) :: ncid                                   !< NetCDF file ID
integer,intent(in) :: its                                    !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0,geom%nl0,nam%nv) !< Read field

! Local variables
integer :: iv,il0,dum
integer :: fld_id
real(kind_real) :: fld_loc(geom%nlon,geom%nlat)
character(len=3) :: ilchar
character(len=1024) :: subr = 'model_aro_read'

! Initialize field
call msr(fld)

do iv=1,nam%nv
   ! 3d variable
   do il0=1,nam%nl
      ! Get id
      write(ilchar,'(i3.3)') nam%levs(il0)
      call ncerr(subr,nf90_inq_varid(ncid,'S'//ilchar//trim(nam%varname(iv)),fld_id))

      ! Read data
      call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc))
      fld(:,il0,iv) = pack(real(fld_loc,kind_real),mask=.true.)
   end do

   if (trim(nam%addvar2d(iv))/='') then
      ! 2d variable

      ! Get id
      call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

      ! Read data
      call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc))
      fld(:,geom%nl0,iv) = pack(real(fld_loc,kind_real),mask=.true.)

      ! Variable change for surface pressure
      if (trim(nam%addvar2d(iv))=='SURFPRESSION') fld(:,geom%nl0,iv) = exp(fld(:,geom%nl0,iv))
   end if
end do

! Use timeslot to avoid warning
dum = its

end subroutine model_aro_read

!----------------------------------------------------------------------
! Subroutine: model_aro_write
!> Purpose: write AROME field
!----------------------------------------------------------------------
subroutine model_aro_write(geom,ncid,varname,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                    !< Geometry
integer,intent(in) :: ncid                           !< NetCDF file ID
character(len=*),intent(in) :: varname               !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0,geom%nl0) !< Written field

! Local variables
integer :: il0,info
integer :: nlon_id,nlat_id,nlev_id,fld_id,lon_id,lat_id
real(kind_real) :: fld_loc(geom%nlon,geom%nlat)
logical :: mask_unpack(geom%nlon,geom%nlat)
character(len=1024) :: subr = 'model_aro_write'

! Initialization
mask_unpack = .true.

! Get variable id
info = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (info/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   info = nf90_inq_dimid(ncid,'X',nlon_id)
   if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'X',geom%nlon,nlon_id))
   info = nf90_inq_dimid(ncid,'Y',nlat_id)
   if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Y',geom%nlat,nlat_id))
   info = nf90_inq_dimid(ncid,'Z',nlev_id)
   if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Z',geom%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,geom%nl0
   if (isanynotmsr(fld(:,il0))) then
      call msr(fld_loc)
      fld_loc = unpack(fld(:,il0),mask=mask_unpack,field=fld_loc)
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/geom%nlon,geom%nlat,1/)))
   end if
end do

! Write coordinates
info = nf90_inq_varid(ncid,'longitude',lon_id)
if (info/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   info = nf90_inq_dimid(ncid,'X',nlon_id)
   if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'X',geom%nlon,nlon_id))
   info = nf90_inq_dimid(ncid,'Y',nlat_id)
   if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Y',geom%nlat,nlat_id))
   call ncerr(subr,nf90_def_var(ncid,'longitude',ncfloat,(/nlon_id,nlat_id/),lon_id))
   call ncerr(subr,nf90_def_var(ncid,'latitude',ncfloat,(/nlon_id,nlat_id/),lat_id))
   call ncerr(subr,nf90_enddef(ncid))
   fld_loc = unpack(geom%lon*rad2deg,mask=mask_unpack,field=fld_loc)
   call ncerr(subr,nf90_put_var(ncid,lon_id,fld_loc))
   fld_loc = unpack(geom%lat*rad2deg,mask=mask_unpack,field=fld_loc)
   call ncerr(subr,nf90_put_var(ncid,lat_id,fld_loc))
end if

end subroutine model_aro_write

end module model_aro
