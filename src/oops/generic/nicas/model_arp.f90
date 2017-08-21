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
use type_ndata, only: ndatatype,ndata_alloc

implicit none

private
public :: model_arp_coord,model_arp_read,model_arp_write

contains

!----------------------------------------------------------------------
! Subroutine: model_arp_coord
!> Purpose: get ARPEGE coordinates
!----------------------------------------------------------------------
subroutine model_arp_coord(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ilon,ilat,il0
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,a_id,b_id
real(kind=8),allocatable :: lon(:,:),lat(:,:),a(:),b(:)

character(len=1024) :: subr = 'model_arp_coord'

! Open file and get dimensions
call msi(ndata%nlon)
call msi(ndata%nlat)
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'longitude',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'latitude',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=ndata%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=ndata%nlat))
call ncerr(subr,nf90_inq_dimid(ncid,'Z',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=ndata%nlev))

! Allocation
allocate(lon(ndata%nlon,ndata%nlat))
allocate(lat(ndata%nlon,ndata%nlat))
allocate(ndata%rgmask(ndata%nlon,ndata%nlat))
allocate(a(ndata%nlev+1))
allocate(b(ndata%nlev+1))

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
do ilon=1,ndata%nlon
   do ilat=1,ndata%nlat
      if (lon(ilon,ilat)>-1000.0) then
         ! Valid point
         ndata%rgmask(ilon,ilat) = .true.

         ! Convert to radian
         lon(ilon,ilat) = lon(ilon,ilat)*real(deg2rad,kind=8)
         lat(ilon,ilat) = lat(ilon,ilat)*real(deg2rad,kind=8)
      else
         ! Missing point
         ndata%rgmask(ilon,ilat) = .false.
         call msr(lon(ilon,ilat))
         call msr(lat(ilon,ilat))
      end if
   end do
end do
ndata%nc0 = count(ndata%rgmask)

! Pack
call ndata_alloc(ndata)
ndata%lon = pack(real(lon,kind_real),mask=ndata%rgmask)
ndata%lat = pack(real(lat,kind_real),mask=ndata%rgmask)
ndata%mask = .true.

! Compute normalized area
ndata%area = 4.0*pi

! Vertical unit
if (nam%logpres) then
   do il0=1,ndata%nl0
      if (nam%levs(il0)<=ndata%nlev) then
         ndata%vunit(il0) = log(0.5*(a(nam%levs(il0))+a(nam%levs(il0)+1))+0.5*(b(nam%levs(il0))+b(nam%levs(il0)+1))*ps)
      else
         ndata%vunit(il0) = log(ps)
      end if
   end do
else
   ndata%vunit = float(nam%levs(1:ndata%nl0))
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
subroutine model_arp_read(ncid,varname,ndata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                              !< NetCDF file ID
character(len=*),intent(in) :: varname                  !< Variable name
type(ndatatype),intent(in) :: ndata                     !< Sampling data
real(kind_real),intent(out) :: fld(ndata%nc0,ndata%nl0) !< Read field

! Local variables
integer :: il0
integer :: fld_id
real(kind_real) :: fld_loc(ndata%nlon,ndata%nlat,ndata%nl0)
character(len=1024) :: subr = 'model_arp_read'

! Initialize field
call msr(fld)

! Get variable id
call ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

! Read data
do il0=1,ndata%nl0
   call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc(:,:,il0),(/1,1,nam%levs(il0)/),(/ndata%nlon,ndata%nlat,1/)))
end do

! Pack data
do il0=1,ndata%nl0
   fld(:,il0) = pack(real(fld_loc(:,:,il0),kind_real),mask=ndata%rgmask)
end do

end subroutine model_arp_read

!----------------------------------------------------------------------
! Subroutine: model_arp_write
!> Purpose: write ARPEGE field
!----------------------------------------------------------------------
subroutine model_arp_write(ncid,varname,ndata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                             !< NetCDF file ID
character(len=*),intent(in) :: varname                 !< Variable name
type(ndatatype),intent(in) :: ndata                    !< Sampling data
real(kind_real),intent(in) :: fld(ndata%nc0,ndata%nl0) !< Written field

! Local variables
integer :: il0,ierr,ilon,ilat
integer :: nlon_id,nlat_id,nlev_id,fld_id,lat_id,lon_id
real(kind_real) :: fld_loc(ndata%nlon,ndata%nlat)
character(len=1024) :: subr = 'model_arp_write'

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'longitude',nlon_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'longitude',ndata%nlon,nlon_id))
   ierr = nf90_inq_dimid(ncid,'latitude',nlat_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'latitude',ndata%nlat,nlat_id))
   ierr = nf90_inq_dimid(ncid,'level',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'level',ndata%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,ndata%nl0
   if (isanynotmsr(fld(:,il0))) then
      call msr(fld_loc)
      fld_loc = unpack(fld(:,il0),mask=ndata%rgmask,field=fld_loc)
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/ndata%nlon,ndata%nlat,1/)))
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
   fld_loc = unpack(ndata%lon,mask=ndata%rgmask,field=fld_loc)
   do ilon=1,ndata%nlon
      do ilat=1,ndata%nlat
         if (isnotmsr(fld_loc(ilon,ilat))) fld_loc(ilon,ilat) = fld_loc(ilon,ilat)*rad2deg
      end do
   end do
   call ncerr(subr,nf90_put_var(ncid,lon_id,fld_loc))
   call msr(fld_loc)
   fld_loc = unpack(ndata%lat,mask=ndata%rgmask,field=fld_loc)
   do ilon=1,ndata%nlon
      do ilat=1,ndata%nlat
         if (isnotmsr(fld_loc(ilon,ilat))) fld_loc(ilon,ilat) = fld_loc(ilon,ilat)*rad2deg
      end do
   end do
   call ncerr(subr,nf90_put_var(ncid,lat_id,fld_loc))
end if

end subroutine model_arp_write

end module model_arp
