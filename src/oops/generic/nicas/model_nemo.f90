!----------------------------------------------------------------------
! Module: module_nemo.f90
!> Purpose: NEMO model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_nemo

use module_namelist, only: nam
use netcdf
use tools_const, only: req,deg2rad,sphere_dist
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_ndata, only: ndatatype,ndata_alloc

implicit none

private
public :: model_nemo_coord,model_nemo_read,model_nemo_write

contains

!----------------------------------------------------------------------
! Subroutine: model_nemo_coord
!> Purpose: get NEMO coordinates
!----------------------------------------------------------------------
subroutine model_nemo_coord(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: il0,ic0,jc0,ilat,jlat,klat,ilon,jlon,klon,i
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,tmask_id,e1t_id,e2t_id
integer(kind=1),allocatable :: tmask(:,:,:)
real(kind=4),allocatable :: lon(:,:),lat(:,:),e1t(:,:,:),e2t(:,:,:)
character(len=1024) :: subr = 'model_nemo_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'x',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'y',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=ndata%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=ndata%nlat))
ndata%nc0 = ndata%nlon*ndata%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'z',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=ndata%nlev))

! Allocation
allocate(lon(ndata%nlon,ndata%nlat))
allocate(lat(ndata%nlon,ndata%nlat))
allocate(tmask(ndata%nlon,ndata%nlat,ndata%nl0))
allocate(e1t(ndata%nlon,ndata%nlat,ndata%nl0))
allocate(e2t(ndata%nlon,ndata%nlat,ndata%nl0))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'nav_lon',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'nav_lat',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'tmask',tmask_id))
call ncerr(subr,nf90_inq_varid(ncid,'e1t',e1t_id))
call ncerr(subr,nf90_inq_varid(ncid,'e2t',e2t_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon,(/1,1/),(/ndata%nlon,ndata%nlat/)))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat,(/1,1/),(/ndata%nlon,ndata%nlat/)))
do il0=1,ndata%nl0
   call ncerr(subr,nf90_get_var(ncid,tmask_id,tmask(:,:,il0),(/1,1,nam%levs(il0),1/),(/ndata%nlon,ndata%nlat,1,1/)))
   call ncerr(subr,nf90_get_var(ncid,e1t_id,e1t(:,:,il0),(/1,1,1/),(/ndata%nlon,ndata%nlat,1/)))
   call ncerr(subr,nf90_get_var(ncid,e2t_id,e2t(:,:,il0),(/1,1,1/),(/ndata%nlon,ndata%nlat,1/)))
end do
call ncerr(subr,nf90_close(ncid))

! Compute normalized area
allocate(ndata%area(ndata%nl0))
do il0=1,ndata%nl0
   ndata%area(il0) = sum(e1t(:,:,il0)*e2t(:,:,il0),mask=tmask(:,:,il0)>0.0)/req**2
end do

! Convert to radian
lon = lon*real(deg2rad,kind=4)
lat = lat*real(deg2rad,kind=4)

if (nam%network) then
   ! Find grid neighbors
   allocate(ndata%net_nnb(ndata%nc0))
   ndata%net_nnb = 8
   allocate(ndata%net_inb(8,ndata%nc0))
   do ilat=1,ndata%nlat
      do ilon=1,ndata%nlon
         ic0 = (ilat-1)*ndata%nlon+ilon
         i = 0
         do jlat=ilat-1,ilat+1
            klat = jlat
            if (klat==0) klat = ndata%nlat
            if (klat==ndata%nlat+1) klat = 1
            do jlon=ilon-1,ilon+1
               klon = jlon
               if (klon==0) klon = ndata%nlon
               if (klon==ndata%nlon+1) klon = 1
               if ((jlat/=ilat).or.(jlon/=ilon)) then
                  i = i+1
                  jc0 = (klat-1)*ndata%nlon+klon
                  ndata%net_inb(i,ic0) = jc0
               end if
            end do
         end do
      end do
   end do
end if

! Pack
call ndata_alloc(ndata)
ndata%lon = pack(real(lon,kind_real),mask=.true.)
ndata%lat = pack(real(lat,kind_real),mask=.true.)
do il0=1,ndata%nl0
   ! Land/sea mask
   ndata%mask(:,il0) = pack(tmask(:,:,il0)>0,mask=.true.)
end do

! Vertical unit
ndata%vunit = float(nam%levs(1:ndata%nl0))

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(tmask)

end subroutine model_nemo_coord

!----------------------------------------------------------------------
! Subroutine: model_nemo_read
!> Purpose: read NEMO field
!----------------------------------------------------------------------
subroutine model_nemo_read(ncid,varname,ndata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                              !< NetCDF file ID
character(len=*),intent(in) :: varname                  !< Variable name
type(ndatatype),intent(in) :: ndata                     !< Sampling data
real(kind_real),intent(out) :: fld(ndata%nc0,ndata%nl0) !< Read field

! Local variables
integer :: il0
integer :: fld_id
real(kind=8) :: fld_loc(ndata%nlon,ndata%nlat)
character(len=1024) :: subr = 'model_nemo_read'

! Initialize field
call msr(fld)

! Get variable id
call ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

! Read variable
do il0=1,ndata%nl0
   call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%levs(il0),1/),(/ndata%nlon,ndata%nlat,1,1/)))
   fld(:,il0) = pack(real(fld_loc,kind_real),mask=.true.)
end do

end subroutine model_nemo_read

!----------------------------------------------------------------------
! Subroutine: model_nemo_write
!> Purpose: write NEMO field
!----------------------------------------------------------------------
subroutine model_nemo_write(ncid,varname,ndata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                             !< NetCDF file ID
character(len=*),intent(in) :: varname                 !< Variable name
type(ndatatype),intent(in) :: ndata                    !< Sampling data
real(kind_real),intent(in) :: fld(ndata%nc0,ndata%nl0) !< Written field

! Local variables
integer :: il0,ic0,ierr
integer :: nlon_id,nlat_id,nlev_id,nt_id,fld_id
real(kind_real) :: fld_tmp(ndata%nc0),fld_loc(ndata%nlon,ndata%nlat)
logical :: mask_unpack(ndata%nlon,ndata%nlat)
character(len=1024) :: subr = 'model_nemo_write'

! Initialization
mask_unpack = .true.

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'x',nlon_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'x',ndata%nlon,nlon_id))
   ierr = nf90_inq_dimid(ncid,'y',nlat_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'y',ndata%nlat,nlat_id))
   ierr = nf90_inq_dimid(ncid,'z',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'z',ndata%nl0,nlev_id))
   ierr = nf90_inq_dimid(ncid,'t',nt_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'t',1,nt_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id,nt_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,ndata%nl0
   if (isanynotmsr(fld(:,il0))) then
      call msr(fld_tmp)
      do ic0=1,ndata%nc0
         if (ndata%mask(ic0,il0)) fld_tmp(ic0) = fld(ic0,il0)
      end do
      fld_loc = unpack(fld_tmp,mask=mask_unpack,field=fld_loc)
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0,1/),(/ndata%nlon,ndata%nlat,1,1/)))
   end if
end do

end subroutine model_nemo_write

end module model_nemo
