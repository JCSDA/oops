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

use module_namelist, only: nam
use netcdf
use tools_const, only: pi,req,deg2rad,ps,rad2deg
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msi,msr,isanynotmsr,isnotmsr,isnotmsi
use tools_nc, only: ncerr,ncfloat
use type_esmf, only: esmf_create_field
use type_sdata, only: sdatatype,sdata_alloc
implicit none

private
public :: model_arp_coord,model_arp_read,model_arp_write

contains

!----------------------------------------------------------------------
! Subroutine: model_arp_coord
!> Purpose: get ARPEGE coordinates
!----------------------------------------------------------------------
subroutine model_arp_coord(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ilon,ilat,il0
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,a_id,b_id
real(kind=8),allocatable :: lon(:,:),lat(:,:),a(:),b(:)

character(len=1024) :: subr = 'model_arp_coord'

! Open file and get dimensions
call msi(sdata%nlon)
call msi(sdata%nlat)
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'longitude',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'latitude',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=sdata%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=sdata%nlat))
call ncerr(subr,nf90_inq_dimid(ncid,'Z',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=sdata%nlev))

! Allocation
allocate(lon(sdata%nlon,sdata%nlat))
allocate(lat(sdata%nlon,sdata%nlat))
allocate(sdata%rgmask(sdata%nlon,sdata%nlat))
allocate(a(sdata%nlev+1))
allocate(b(sdata%nlev+1))

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

! Compute normalized area
allocate(sdata%area(sdata%nl0))
sdata%area = 4.0*pi

! Define mask for the reduced Gaussian grid
do ilon=1,sdata%nlon
   do ilat=1,sdata%nlat
      if (lon(ilon,ilat)>-1000.0) then
         ! Valid point
         sdata%rgmask(ilon,ilat) = .true.

         ! Convert to radian
         lon(ilon,ilat) = lon(ilon,ilat)*real(deg2rad,kind=8)
         lat(ilon,ilat) = lat(ilon,ilat)*real(deg2rad,kind=8)
      else
         ! Missing point
         sdata%rgmask(ilon,ilat) = .false.
         call msr(lon(ilon,ilat))
         call msr(lat(ilon,ilat))
      end if
   end do
end do
sdata%nc0 = count(sdata%rgmask)

! Pack
call sdata_alloc(sdata)
sdata%lon = pack(real(lon,kind_real),mask=sdata%rgmask)
sdata%lat = pack(real(lat,kind_real),mask=sdata%rgmask)
sdata%mask = .true.

! Vertical unit
if (nam%logpres) then
   do il0=1,sdata%nl0
      if (nam%levs(il0)<=sdata%nlev) then
         sdata%vunit(il0) = log(0.5*(a(nam%levs(il0))+a(nam%levs(il0)+1))+0.5*(b(nam%levs(il0))+b(nam%levs(il0)+1))*ps)
      else
         sdata%vunit(il0) = log(ps)
      end if
   end do
else
   sdata%vunit = float(nam%levs(1:sdata%nl0))
end if

! Create ESMF field from mesh
call esmf_create_field(sdata,sdata%nc0,sdata%lon,sdata%lat,any(sdata%mask,dim=2),sdata%c0field)

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
subroutine model_arp_read(ncid,varname,sdata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                              !< NetCDF file ID
character(len=*),intent(in) :: varname                  !< Variable name
type(sdatatype),intent(in) :: sdata                     !< Sampling data
real(kind_real),intent(out) :: fld(sdata%nc0,sdata%nl0) !< Read field

! Local variables
integer :: il0
integer :: fld_id
real(kind_real) :: fld_loc(sdata%nlon,sdata%nlat,sdata%nl0)
character(len=1024) :: subr = 'model_arp_read'

! Initialize field
call msr(fld)

! Get variable id
call ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

! Read data
do il0=1,sdata%nl0
   call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc(:,:,il0),(/1,1,nam%levs(il0)/),(/sdata%nlon,sdata%nlat,1/)))
end do

! Pack data
do il0=1,sdata%nl0
   fld(:,il0) = pack(real(fld_loc(:,:,il0),kind_real),mask=sdata%rgmask)
end do

end subroutine model_arp_read

!----------------------------------------------------------------------
! Subroutine: model_arp_write
!> Purpose: write ARPEGE field
!----------------------------------------------------------------------
subroutine model_arp_write(ncid,varname,sdata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                             !< NetCDF file ID
character(len=*),intent(in) :: varname                 !< Variable name
type(sdatatype),intent(in) :: sdata                    !< Sampling data
real(kind_real),intent(in) :: fld(sdata%nc0,sdata%nl0) !< Written field

! Local variables
integer :: il0,ierr,ilon,ilat
integer :: nlon_id,nlat_id,nlev_id,fld_id,lat_id,lon_id
real(kind_real) :: fld_loc(sdata%nlon,sdata%nlat)
character(len=1024) :: subr = 'model_arp_write'

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'longitude',nlon_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'longitude',sdata%nlon,nlon_id))
   ierr = nf90_inq_dimid(ncid,'latitude',nlat_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'latitude',sdata%nlat,nlat_id))
   ierr = nf90_inq_dimid(ncid,'level',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'level',sdata%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,sdata%nl0
   if (isanynotmsr(fld(:,il0))) then
      call msr(fld_loc)
      fld_loc = unpack(fld(:,il0),mask=sdata%rgmask,field=fld_loc)
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/sdata%nlon,sdata%nlat,1/)))
   end if
end do

! Write coordinates
ierr = nf90_inq_varid(ncid,'longitude',lon_id)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   call ncerr(subr,nf90_def_var(ncid,'longitude',ncfloat,(/nlon_id,nlat_id/),lon_id))
   call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'latitude',ncfloat,(/nlon_id,nlat_id/),lat_id))
   call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
   call msr(fld_loc)
   fld_loc = unpack(sdata%lon,mask=sdata%rgmask,field=fld_loc)
   do ilon=1,sdata%nlon
      do ilat=1,sdata%nlat
         if (isnotmsr(fld_loc(ilon,ilat))) fld_loc(ilon,ilat) = fld_loc(ilon,ilat)*rad2deg
      end do
   end do
   call ncerr(subr,nf90_put_var(ncid,lon_id,fld_loc))
   call msr(fld_loc)
   fld_loc = unpack(sdata%lat,mask=sdata%rgmask,field=fld_loc)
   do ilon=1,sdata%nlon
      do ilat=1,sdata%nlat
         if (isnotmsr(fld_loc(ilon,ilat))) fld_loc(ilon,ilat) = fld_loc(ilon,ilat)*rad2deg
      end do
   end do
   call ncerr(subr,nf90_put_var(ncid,lat_id,fld_loc))
end if

end subroutine model_arp_write

end module model_arp
