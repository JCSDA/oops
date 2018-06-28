!----------------------------------------------------------------------
! Module: module_nemo
!> Purpose: NEMO model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyrimght © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module model_nemo

use netcdf
use tools_const, only: req,deg2rad,rad2deg
use tools_func, only: sphere_dist
use tools_kinds,only: kind_real
use tools_missing, only: msr,isanynotmsr,isnotmsi
use tools_nc, only: ncfloat
use tools_qsort, only: qsort
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

private
public :: model_nemo_coord,model_nemo_read

contains

!----------------------------------------------------------------------
! Subroutine: model_nemo_coord
!> Purpose: get NEMO coordinates
!----------------------------------------------------------------------
subroutine model_nemo_coord(mpl,nam,geom)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: il0,img,ilat,ilon,ic0
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,tmask_id,e1t_id,e2t_id
integer,allocatable :: mg_to_lon(:),mg_to_lat(:)
integer(kind=1),allocatable :: tmask(:,:,:)
real(kind=4),allocatable :: lon(:,:),lat(:,:),e1t(:,:,:),e2t(:,:,:)
real(kind_real),allocatable :: lon_mg(:),lat_mg(:),area_mg(:)
logical,allocatable :: lmask_mg(:,:)
character(len=1024) :: subr = 'model_nemo_coord'

! Open file and get dimensions
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'x',nlon_id))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'y',nlat_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
geom%nmg = geom%nlon*geom%nlat
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'z',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon,geom%nlat))
allocate(lat(geom%nlon,geom%nlat))
allocate(tmask(geom%nlon,geom%nlat,geom%nl0))
allocate(e1t(geom%nlon,geom%nlat,geom%nl0))
allocate(e2t(geom%nlon,geom%nlat,geom%nl0))
allocate(mg_to_lon(geom%nmg))
allocate(mg_to_lat(geom%nmg))
allocate(lon_mg(geom%nmg))
allocate(lat_mg(geom%nmg))
allocate(area_mg(geom%nmg))
allocate(lmask_mg(geom%nmg,geom%nl0))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'nav_lon',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'nav_lat',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'tmask',tmask_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'e1t',e1t_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'e2t',e2t_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
do il0=1,geom%nl0
   call mpl%ncerr(subr,nf90_get_var(ncid,tmask_id,tmask(:,:,il0),(/1,1,nam%levs(il0),1/),(/geom%nlon,geom%nlat,1,1/)))
   call mpl%ncerr(subr,nf90_get_var(ncid,e1t_id,e1t(:,:,il0),(/1,1,1/),(/geom%nlon,geom%nlat,1/)))
   call mpl%ncerr(subr,nf90_get_var(ncid,e2t_id,e2t(:,:,il0),(/1,1,1/),(/geom%nlon,geom%nlat,1/)))
end do
call mpl%ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*real(deg2rad,kind=4)
lat = lat*real(deg2rad,kind=4)

! Redundant grid
img = 0
do ilon=1,geom%nlon
   do ilat=1,geom%nlat
      img = img+1
      mg_to_lon(img) = ilon
      mg_to_lat(img) = ilat
      lon_mg(img) = real(lon(ilon,ilat),kind_real)
      lat_mg(img) = real(lat(ilon,ilat),kind_real)
      area_mg(img) = real(e1t(ilon,ilat,1)*e2t(ilon,ilat,1),kind_real)/req**2
      do il0=1,geom%nl0
         lmask_mg(img,il0) = (tmask(ilon,ilat,il0)>0)
      end do
   end do
end do
call geom%find_redundant(mpl,lon_mg,lat_mg)

! Pack
call geom%alloc
geom%c0_to_lon = mg_to_lon(geom%c0_to_mg)
geom%c0_to_lat = mg_to_lat(geom%c0_to_mg)
geom%lon = lon_mg(geom%c0_to_mg)
geom%lat = lat_mg(geom%c0_to_mg)
do il0=1,geom%nl0
   geom%mask(:,il0) = lmask_mg(geom%c0_to_mg,il0)
   geom%area(il0) = sum(area_mg(geom%c0_to_mg),geom%mask(:,il0))/req**2
end do

! Vertical unit
do ic0=1,geom%nc0
   geom%vunit(ic0,:) = real(nam%levs(1:geom%nl0),kind_real)
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(tmask)

end subroutine model_nemo_coord

!----------------------------------------------------------------------
! Subroutine: model_nemo_read
!> Purpose: read NEMO field
!----------------------------------------------------------------------
subroutine model_nemo_read(mpl,nam,geom,filename,its,fld)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl                              !< MPI data
type(nam_type),intent(in) :: nam                              !< Namelist
type(geom_type),intent(in) :: geom                            !< Geometry
character(len=*),intent(in) :: filename                       !< File name
integer,intent(in) :: its                                     !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0,ic0,ilon,ilat
integer :: ncid,fld_id
real(kind=8) :: fld_tmp2
real(kind=8),allocatable :: fld_tmp(:,:,:)
real(kind_real) :: fld_c0(geom%nc0)
character(len=1024) :: subr = 'model_nemo_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(geom%nlon,geom%nlat,geom%nl0))

   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
end if

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable

      ! Get variable id
      call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))

      ! Read data
      do il0=1,nam%nl
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0),(/1,1,nam%levs(il0),nam%timeslot(its)/), &
       & (/geom%nlon,geom%nlat,1,1/)))
         select case (trim(nam%varname(iv)))
         case ('un')
            do ilat=1,geom%nlat
               do ilon=1,geom%nlon
                  if (ilon==1) then
                     call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/geom%nlon,ilat,nam%levs(il0),nam%timeslot(its)/)))
                  else
                     call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon-1,ilat,nam%levs(il0),nam%timeslot(its)/)))
                  end if
                  fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
               end do
            end do
         case ('vn')
            do ilat=1,geom%nlat
               do ilon=1,geom%nlon
                  if (ilat==1) then
                     call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon,geom%nlat,nam%levs(il0),nam%timeslot(its)/)))
                  else
                     call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon,ilat-1,nam%levs(il0),nam%timeslot(its)/)))
                  end if
                  fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
               end do
            end do
         end select
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,geom%nl0),(/1,1,nam%timeslot(its)/),(/geom%nlon,geom%nlat,1/)))
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
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(fld_tmp)
end if

end subroutine model_nemo_read

end module model_nemo
