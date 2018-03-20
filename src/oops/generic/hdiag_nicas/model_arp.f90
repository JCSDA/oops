!----------------------------------------------------------------------
! Module: module_arp.f90
!> Purpose: ARPEGE model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_arp

use netcdf
use tools_const, only: pi,req,deg2rad,ps,rad2deg
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msi,msr,isanynotmsr,isnotmsr,isnotmsi
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: model_arp_coord,model_arp_read,model_arp_write

contains

!----------------------------------------------------------------------
! Subroutine: model_arp_coord
!> Purpose: get ARPEGE coordinates
!----------------------------------------------------------------------
subroutine model_arp_coord(nam,geom)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,a_id,b_id
real(kind=8),allocatable :: lon(:,:),lat(:,:),a(:),b(:)

character(len=1024) :: subr = 'model_arp_coord'

! Open file and get dimensions
call msi(geom%nlon)
call msi(geom%nlat)
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'longitude',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'latitude',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
call ncerr(subr,nf90_inq_dimid(ncid,'Z',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon,geom%nlat))
allocate(lat(geom%nlon,geom%nlat))
allocate(geom%rgmask(geom%nlon,geom%nlat))
allocate(a(geom%nlev+1))
allocate(b(geom%nlev+1))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'longitude',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'latitude',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_A',a_id))
call ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_B',b_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call ncerr(subr,nf90_get_var(ncid,a_id,a))
call ncerr(subr,nf90_get_var(ncid,b_id,b))
call ncerr(subr,nf90_close(ncid))

! Define mask for the reduced Gaussian grid
do ilon=1,geom%nlon
   do ilat=1,geom%nlat
      if (lon(ilon,ilat)>-1000.0) then
         ! Valid point
         geom%rgmask(ilon,ilat) = .true.

         ! Convert to radian
         lon(ilon,ilat) = lon(ilon,ilat)*real(deg2rad,kind=8)
         lat(ilon,ilat) = lat(ilon,ilat)*real(deg2rad,kind=8)
      else
         ! Missing point
         geom%rgmask(ilon,ilat) = .false.
         call msr(lon(ilon,ilat))
         call msr(lat(ilon,ilat))
      end if
   end do
end do
geom%nc0 = count(geom%rgmask)

! Pack
call geom%alloc
geom%lon = pack(real(lon,kind_real),mask=geom%rgmask)
geom%lat = pack(real(lat,kind_real),mask=geom%rgmask)
geom%mask = .true.

! Compute normalized area
geom%area = 4.0*pi

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
deallocate(a)
deallocate(b)

end subroutine model_arp_coord

!----------------------------------------------------------------------
! Subroutine: model_arp_read
!> Purpose: read ARPEGE field
!----------------------------------------------------------------------
subroutine model_arp_read(nam,geom,ncid,its,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                              !< Namelist
type(geom_type),intent(in) :: geom                            !< Geometry
integer,intent(in) :: ncid                                    !< NetCDF file ID
integer,intent(in) :: its                                     !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0,dum
integer :: fld_id
real(kind_real) :: fld_loc(geom%nlon,geom%nlat),fld_glb(geom%nc0,geom%nl0)
character(len=3) :: ilchar
character(len=1024) :: subr = 'model_arp_read'

! Initialize field
call msr(fld)

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable
      do il0=1,nam%nl
         ! Get id
         write(ilchar,'(i3.3)') nam%levs(il0)
         call ncerr(subr,nf90_inq_varid(ncid,'S'//ilchar//trim(nam%varname(iv)),fld_id))

         ! Read data
         call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc))
         fld_glb(:,il0) = pack(real(fld_loc,kind_real),mask=geom%rgmask)
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc))
         fld_glb(:,geom%nl0) = pack(real(fld_loc,kind_real),mask=geom%rgmask)

         ! Variable change for surface pressure
         if (trim(nam%addvar2d(iv))=='SURFPRESSION') fld_glb(:,geom%nl0) = exp(fld_glb(:,geom%nl0))
      end if
   end if

   ! Split over processors
   call geom%fld_com_gl(fld_glb,fld(:,:,iv))
end do

! Use timeslot to avoid warning
dum = its

end subroutine model_arp_read

!----------------------------------------------------------------------
! Subroutine: model_arp_write
!> Purpose: write ARPEGE field
!----------------------------------------------------------------------
subroutine model_arp_write(geom,ncid,varname,fld)

implicit none

! Passed variables
type(geom_type),intent(in) :: geom                    !< Geometry
integer,intent(in) :: ncid                            !< NetCDF file ID
character(len=*),intent(in) :: varname                !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
integer :: il0,info,ilon,ilat
integer :: nlon_id,nlat_id,nlev_id,fld_id,lat_id,lon_id
real(kind_real) :: fld_loc(geom%nlon,geom%nlat),fld_glb(geom%nc0,geom%nl0)
character(len=1024) :: subr = 'model_arp_write'

! Local to global
call geom%fld_com_lg(fld,fld_glb)

if (mpl%main) then
   ! Get variable id
   info = nf90_inq_varid(ncid,trim(varname),fld_id)

   ! Define dimensions and variable if necessary
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      info = nf90_inq_dimid(ncid,'longitude',nlon_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'longitude',geom%nlon,nlon_id))
      info = nf90_inq_dimid(ncid,'latitude',nlat_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'latitude',geom%nlat,nlat_id))
      info = nf90_inq_dimid(ncid,'level',nlev_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'level',geom%nl0,nlev_id))
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
   info = nf90_inq_varid(ncid,'longitude',lon_id)
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      call ncerr(subr,nf90_def_var(ncid,'longitude',ncfloat,(/nlon_id,nlat_id/),lon_id))
      call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_def_var(ncid,'latitude',ncfloat,(/nlon_id,nlat_id/),lat_id))
      call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
      call msr(fld_loc)
      fld_loc = unpack(geom%lon,geom%rgmask,fld_loc)
      do ilon=1,geom%nlon
         do ilat=1,geom%nlat
            if (isnotmsr(fld_loc(ilon,ilat))) fld_loc(ilon,ilat) = fld_loc(ilon,ilat)*rad2deg
         end do
      end do
      call ncerr(subr,nf90_put_var(ncid,lon_id,fld_loc))
      call msr(fld_loc)
      fld_loc = unpack(geom%lat,geom%rgmask,fld_loc)
      do ilon=1,geom%nlon
         do ilat=1,geom%nlat
            if (isnotmsr(fld_loc(ilon,ilat))) fld_loc(ilon,ilat) = fld_loc(ilon,ilat)*rad2deg
         end do
      end do
      call ncerr(subr,nf90_put_var(ncid,lat_id,fld_loc))
   end if
end if

end subroutine model_arp_write

end module model_arp
