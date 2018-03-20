!----------------------------------------------------------------------
! Module: module_wrf.f90
!> Purpose: WRF model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_wrf

use netcdf
use tools_const, only: deg2rad,rad2deg,req,ps
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: model_wrf_coord,model_wrf_read,model_wrf_write

contains

!----------------------------------------------------------------------
! Subroutine: model_wrf_coord
!> Purpose: get WRF coordinates
!----------------------------------------------------------------------
subroutine model_wrf_coord(nam,geom)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,pres_id
real(kind=4) :: dx,dy
real(kind=4),allocatable :: lon(:,:),lat(:,:),pres(:)
character(len=1024) :: subr = 'model_wrf_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'west_east',nlon_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call ncerr(subr,nf90_inq_dimid(ncid,'south_north',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
geom%nc0 = geom%nlon*geom%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'bottom_top',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon,geom%nlat))
allocate(lat(geom%nlon,geom%nlat))
allocate(geom%rgmask(geom%nlon,geom%nlat))
allocate(pres(geom%nlev))

! Initialization
geom%rgmask = .true.

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'XLONG',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'XLAT',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'PB',pres_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon,(/1,1,1/),(/geom%nlon,geom%nlat,1/)))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat,(/1,1,1/),(/geom%nlon,geom%nlat,1/)))
call ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'DX',dx))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'DY',dy))
call ncerr(subr,nf90_close(ncid))

! Conversion to radian
lon = lon*real(deg2rad,kind=4)
lat = lat*real(deg2rad,kind=4)

! Pack
call geom%alloc
geom%lon = pack(real(lon,kind_real),mask=.true.)
geom%lat = pack(real(lat,kind_real),mask=.true.)
geom%mask = .true.

! Compute normalized area
geom%area = float(geom%nlon*geom%nlat)*dx*dy/req**2

! Vertical unit
if (nam%logpres) then
   geom%vunit(1:nam%nl) = log(pres(nam%levs(1:nam%nl)))
   if (geom%nl0>nam%nl) geom%vunit(geom%nl0) = log(ps)
else
   geom%vunit = float(nam%levs(1:geom%nl0))
end if

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(pres)

end subroutine model_wrf_coord

!----------------------------------------------------------------------
! Subroutine: model_wrf_read
!> Purpose: read WRF field
!----------------------------------------------------------------------
subroutine model_wrf_read(nam,geom,ncid,its,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                              !< Namelist
type(geom_type),intent(in) :: geom                            !< Geometry
integer,intent(in) :: ncid                                    !< NetCDF file ID
integer,intent(in) :: its                                     !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0
integer :: fld_id
real(kind=4),allocatable :: fld_loc(:,:)
real(kind_real) :: fld_glb(geom%nc0,geom%nl0)
character(len=1024) :: subr = 'model_wrf_read'

! Initialize field
call msr(fld)

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable
   
      ! Allocation
      select case (trim(nam%varname(iv)))
      case ('U')
         allocate(fld_loc(geom%nlon+1,geom%nlat))
      case ('V')
         allocate(fld_loc(geom%nlon,geom%nlat+1))
      case default
         allocate(fld_loc(geom%nlon,geom%nlat))
      end select
   
      ! Get variable id
      call ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))
   
      ! 3d variable
      do il0=1,nam%nl
         select case (trim(nam%varname(iv)))
         case ('U')
            call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%levs(il0),nam%timeslot(its)/),(/geom%nlon+1,geom%nlat,1,1/)))
            fld_glb(:,il0) = pack(real(0.5*(fld_loc(1:geom%nlon,:)+fld_loc(2:geom%nlon+1,:)),kind_real),mask=.true.)
         case ('V')
            call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%levs(il0),nam%timeslot(its)/),(/geom%nlon,geom%nlat+1,1,1/)))
            fld_glb(:,il0) = pack(real(0.5*(fld_loc(:,1:geom%nlat)+fld_loc(:,2:geom%nlat+1)),kind_real),mask=.true.)
         case default
            call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%levs(il0),nam%timeslot(its)/),(/geom%nlon,geom%nlat,1,1/)))
            fld_glb(:,il0) = pack(real(fld_loc,kind_real),mask=.true.)
         end select
      end do
   
      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable
   
         ! Get id
         call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))
   
         ! Read data
         call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%timeslot(its)/),(/geom%nlon,geom%nlat,1/)))
         fld_glb(:,geom%nl0) = pack(real(fld_loc,kind_real),.true.)
      end if
   
      ! Release memory
      deallocate(fld_loc)
   end if

   ! Split over processors
   call geom%fld_com_gl(fld_glb,fld(:,:,iv))
end do

end subroutine model_wrf_read

!----------------------------------------------------------------------
! Subroutine: model_wrf_write
!> Purpose: write WRF field
!----------------------------------------------------------------------
subroutine model_wrf_write(geom,ncid,varname,fld)

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
character(len=1024) :: subr = 'model_wrf_write'

! Local to global
call geom%fld_com_lg(fld,fld_glb)

if (mpl%main) then
   ! Get variable id
   info = nf90_inq_varid(ncid,trim(varname),fld_id)
   
   ! Define dimensions and variable if necessary
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      info = nf90_inq_dimid(ncid,'west_east',nlon_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'west_east',geom%nlon,nlon_id))
      info = nf90_inq_dimid(ncid,'south_north',nlat_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'south_north',geom%nlat,nlat_id))
      info = nf90_inq_dimid(ncid,'bottom_top',nlev_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'bottom_top',geom%nl0,nlev_id))
      call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
      call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
   end if
   
   ! Write data
   do il0=1,geom%nl0
      if (isanynotmsr(fld_glb(:,il0))) then
         call msr(fld_loc)
         fld_loc = unpack(fld_glb(:,il0),geom%rgmask,fld_loc)
         call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/geom%nlon,geom%nlat,1/)))
      end if
   end do

   ! Write coordinates
   info = nf90_inq_varid(ncid,'XLONG',lon_id)
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      call ncerr(subr,nf90_def_var(ncid,'XLONG',ncfloat,(/nlon_id,nlat_id/),lon_id))
      call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_def_var(ncid,'XLAT',ncfloat,(/nlon_id,nlat_id/),lat_id))
      call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
      call msr(fld_loc)
      fld_loc = unpack(geom%lon*rad2deg,geom%rgmask,fld_loc)
      call ncerr(subr,nf90_put_var(ncid,lon_id,fld_loc))
      fld_loc = unpack(geom%lat*rad2deg,geom%rgmask,fld_loc)
      call ncerr(subr,nf90_put_var(ncid,lat_id,fld_loc))
   end if
end if

end subroutine model_wrf_write

end module model_wrf
