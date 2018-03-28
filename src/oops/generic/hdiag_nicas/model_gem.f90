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

use netcdf
use tools_const, only: pi,req,deg2rad,rad2deg,ps
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: model_gem_coord,model_gem_read,model_gem_write

contains

!----------------------------------------------------------------------
! Subroutine: model_gem_coord
!> Purpose: get GEM coordinates
!----------------------------------------------------------------------
subroutine model_gem_coord(nam,geom)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,a_id,b_id
real(kind=8),allocatable :: lon(:,:),lat(:,:),a(:),b(:)
character(len=1024) :: subr = 'model_gem_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'lon',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'lat',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
geom%nc0 = geom%nlon*geom%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'lev',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon,geom%nlat))
allocate(lat(geom%nlon,geom%nlat))
allocate(geom%rgmask(geom%nlon,geom%nlat))
allocate(a(geom%nlev))
allocate(b(geom%nlev))

! Initialization
geom%rgmask = .true.

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

! Convert to radian
lon(:,1) = lon(:,1)*deg2rad
lat(1,:) = lat(1,:)*deg2rad

! Fill arrays
do ilat=1,geom%nlat
   lon(:,ilat) = lon(:,1)
end do
do ilon=1,geom%nlon
   lat(ilon,:) = lat(1,:)
end do

! Pack
call geom%alloc
geom%lon = pack(real(lon,kind_real),mask=.true.)
geom%lat = pack(real(lat,kind_real),mask=.true.)
geom%mask = .true.

! Compute normalized area
geom%area = 4.0*pi

! Vertical unit
if (nam%logpres) then
   geom%vunit(1:nam%nl) = log(a(nam%levs(1:nam%nl))+b(nam%levs(1:nam%nl))*ps)
   if (geom%nl0>nam%nl) geom%vunit(geom%nl0) = log(ps)
else
   geom%vunit = float(nam%levs(1:geom%nl0))
end if

! Not redundant grid
geom%redgrid = .false.

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
subroutine model_gem_read(nam,geom,ncid,its,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                             !< Namelist
type(geom_type),intent(in) :: geom                           !< Geometry
integer,intent(in) :: ncid                                   !< NetCDF file ID
integer,intent(in) :: its                                    !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0,xt,dum
integer :: fld_id
integer,allocatable :: fld_loc_int(:,:)
real(kind_real) :: add_offset,scale_factor
real(kind_real) :: fld_glb(geom%nc0,geom%nl0)
real(kind_real),allocatable :: fld_loc(:,:)
character(len=1024) :: subr = 'model_gem_read'

! Initialize field
call msr(fld)

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable

      ! Get variable id
      call ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))

      ! Check variable type
      call ncerr(subr,nf90_inquire_variable(ncid,fld_id,xtype=xt))
      if (xt==nf90_short) then
         allocate(fld_loc_int(geom%nlon,geom%nlat))
      elseif (xt==nf90_double) then
         allocate(fld_loc(geom%nlon,geom%nlat))
      else
         call msgerror('wrong variable type')
      end if

      ! 3d variable
      do il0=1,nam%nl
         if (xt==nf90_short) then
            call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc_int,(/1,1,nam%levs(il0)/),(/geom%nlon,geom%nlat,1/)))
            call ncerr(subr,nf90_get_att(ncid,fld_id,'add_offset',add_offset))
            call ncerr(subr,nf90_get_att(ncid,fld_id,'scale_factor',scale_factor))
            fld_glb(:,il0) = pack(add_offset+scale_factor*float(fld_loc_int),mask=.true.)
         elseif (xt==nf90_double) then
            call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%levs(il0)/),(/geom%nlon,geom%nlat,1/)))
            fld_glb(:,il0) = pack(fld_loc,mask=.true.)
         else
            call msgerror('wrong netcdf variable type')
         end if
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         if (xt==nf90_short) then
            call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc_int,(/1,1/),(/geom%nlon,geom%nlat/)))
            call ncerr(subr,nf90_get_att(ncid,fld_id,'add_offset',add_offset))
            call ncerr(subr,nf90_get_att(ncid,fld_id,'scale_factor',scale_factor))
            fld_glb(:,geom%nl0) = pack(add_offset+scale_factor*float(fld_loc_int),mask=.true.)
         elseif (xt==nf90_double) then
            call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1/),(/geom%nlon,geom%nlat/)))
            fld_glb(:,geom%nl0) = pack(fld_loc,mask=.true.)
         else
            call msgerror('wrong netcdf variable type')
         end if
      end if
   end if

   ! Split over processors
   call geom%fld_com_gl(fld_glb,fld(:,:,iv))
end do

! Use timeslot to avoid warning
dum = its

end subroutine model_gem_read

!----------------------------------------------------------------------
! Subroutine: model_gem_write
!> Purpose: write GEM field
!----------------------------------------------------------------------
subroutine model_gem_write(geom,ncid,varname,fld)

implicit none

! Passed variables
type(geom_type),intent(in) :: geom                    !< Geometry
integer,intent(in) :: ncid                            !< NetCDF file ID
character(len=*),intent(in) :: varname                !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
integer :: il0,info
integer :: nlon_id,nlat_id,nlev_id,fld_id,lon_id,lat_id
real(kind_real) :: fld_loc(geom%nlon,geom%nlat),fld_glb(geom%nc0,geom%nl0)
character(len=1024) :: subr = 'model_gem_write'

! Local to global
call geom%fld_com_lg(fld,fld_glb)

if (mpl%main) then
   ! Get variable id
   info = nf90_inq_varid(ncid,trim(varname),fld_id)

   ! Define dimensions and variable if necessary
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      info = nf90_inq_dimid(ncid,'lon',nlon_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'lon',geom%nlon,nlon_id))
      info = nf90_inq_dimid(ncid,'lat',nlat_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'lat',geom%nlat,nlat_id))
      info = nf90_inq_dimid(ncid,'lev',nlev_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'lev',geom%nl0,nlev_id))
      call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
      call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
   end if

   ! Write data
   do il0=1,geom%nl0
      if (isanynotmsr(fld(:,il0))) then
         call msr(fld_loc)
         fld_loc = unpack(fld_glb(:,il0),geom%rgmask,fld_loc)
         call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/geom%nlon,geom%nlat,1/)))
      end if
   end do

   ! Write coordinates
   info = nf90_inq_varid(ncid,'lon',lon_id)
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nlon_id,nlat_id/),lon_id))
      call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nlon_id,nlat_id/),lat_id))
      call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
      call msr(fld_loc)
      fld_loc = unpack(geom%lon*rad2deg,geom%rgmask,fld_loc)
      call ncerr(subr,nf90_put_var(ncid,lon_id,fld_loc))
      fld_loc = unpack(geom%lat*rad2deg,geom%rgmask,fld_loc)
      call ncerr(subr,nf90_put_var(ncid,lat_id,fld_loc))
   end if
end if

end subroutine model_gem_write

end module model_gem
