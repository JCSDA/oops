!----------------------------------------------------------------------
! Module: module_ifs
!> Purpose: IFS model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module model_ifs

use netcdf
use tools_const, only: pi,deg2rad,rad2deg,ps
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msr,isanynotmsr,isallnotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: model_ifs_coord,model_ifs_read

contains

!----------------------------------------------------------------------
! Subroutine: model_ifs_coord
!> Purpose: get IFS coordinates
!----------------------------------------------------------------------
subroutine model_ifs_coord(nam,geom)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: ilon,ilat,ic0
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,pres_id
real(kind=4),allocatable :: lon(:),lat(:),pres(:)
character(len=1024) :: subr = 'model_ifs_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'longitude',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'latitude',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
geom%nmg = geom%nlon*geom%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'level',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon))
allocate(lat(geom%nlat))
allocate(pres(geom%nlev))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'longitude',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'latitude',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'pf',pres_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*real(deg2rad,kind=4)
lat = lat*real(deg2rad,kind=4)

! Not redundant grid
call geom%find_redundant

! Pack
call geom%alloc
ic0 = 0
do ilon=1,geom%nlon
   do ilat=1,geom%nlat
      ic0 = ic0+1
      geom%c0_to_lon(ic0) = ilon
      geom%c0_to_lat(ic0) = ilat
      geom%lon(ic0) = real(lon(ilon),kind_real)
      geom%lat(ic0) = real(lat(ilat),kind_real)
      geom%mask(ic0,:) = .true.
   end do
end do

! Compute normalized area
geom%area = 4.0*pi

! Vertical unit
do ic0=1,geom%nc0
   if (nam%logpres) then
      geom%vunit(ic0,1:nam%nl) = log(pres(nam%levs(1:nam%nl)))
      if (geom%nl0>nam%nl) geom%vunit(ic0,geom%nl0) = log(ps)
   else
      geom%vunit(ic0,:) = real(nam%levs(1:geom%nl0),kind_real)
   end if
end do

! Release memory
deallocate(lon)
deallocate(lat)

end subroutine model_ifs_coord

!----------------------------------------------------------------------
! Subroutine: model_ifs_read
!> Purpose: read IFS field
!----------------------------------------------------------------------
subroutine model_ifs_read(nam,geom,filename,its,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                              !< Namelist
type(geom_type),intent(in) :: geom                            !< Geometry
character(len=*),intent(in) :: filename                       !< File name
integer,intent(in) :: its                                     !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0,ic0,ilon,ilat
integer :: ncid,fld_id
real(kind=4),allocatable :: fld_tmp(:,:,:)
real(kind_real) :: fld_c0(geom%nc0)
character(len=1024) :: subr = 'model_ifs_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(geom%nlon,geom%nlat,geom%nl0))

   ! Open file
   call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
end if

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable

      ! Get variable id
      call ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))

      ! Read data
      do il0=1,nam%nl
         call ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0),(/1,1,nam%levs(il0),nam%timeslot(its)/), &
       & (/geom%nlon,geom%nlat,1,1/)))
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         call ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,geom%nl0),(/1,1,nam%timeslot(its)/),(/geom%nlon,geom%nlat,1/)))
      end if
   end if

   ! Global to local
   do il0=1,geom%nl0
      if (mpl%main) then
         do ic0=1,geom%nc0
            ilon = geom%c0_to_lon(ic0)
            ilat = geom%c0_to_lat(ic0)
            fld_c0(ic0) = real(fld_tmp(ilon,ilat,il0),kind_real)
         end do
      end if
      call mpl%scatterv(geom%proc_to_nc0a,geom%nc0,fld_c0,geom%nc0a,fld(:,il0,iv))
   end do
end do

if (mpl%main) then
   ! Close file
   call ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(fld_tmp)
end if

end subroutine model_ifs_read

end module model_ifs
