!----------------------------------------------------------------------
! Module: module_mpas
!> Purpose: MPAS model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module model_mpas

use netcdf
use tools_const, only: pi,ps
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr,isanynotmsr
use tools_nc, only: ncfloat
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

private
public :: model_mpas_coord,model_mpas_read

contains

!----------------------------------------------------------------------
! Subroutine: model_mpas_coord
!> Purpose: get MPAS coordinates
!----------------------------------------------------------------------
subroutine model_mpas_coord(mpl,nam,geom)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: ncid,ng_id,nlev_id,lon_id,lat_id,pres_id
integer :: ic0
real(kind=4),allocatable :: lon(:),lat(:),pres(:)
character(len=1024) :: subr = 'model_mpas_coord'

! Open file and get dimensions
call msi(geom%nlon)
call msi(geom%nlat)
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nCells',ng_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,ng_id,len=geom%nmg))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nVertLevels',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nmg))
allocate(lat(geom%nmg))
allocate(pres(geom%nlev))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'lonCell',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'latCell',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'pressure_base',pres_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call mpl%ncerr(subr,nf90_close(ncid))

! Not redundant grid
call geom%find_redundant(mpl)

! Pack
call geom%alloc
geom%lon = real(lon,kind_real)
geom%lat = real(lat,kind_real)
geom%mask = .true.

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
deallocate(pres)

end subroutine model_mpas_coord

!----------------------------------------------------------------------
! Subroutine: model_mpas_read
!> Purpose: read MPAS field
!----------------------------------------------------------------------
subroutine model_mpas_read(mpl,nam,geom,filename,its,fld)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl                              !< MPI data
type(nam_type),intent(in) :: nam                              !< Namelist
type(geom_type),intent(in) :: geom                            !< Geometry
character(len=*),intent(in) :: filename                       !< File name
integer,intent(in) :: its                                     !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0,img,ic0
integer :: ncid,fld_id
real(kind=4),allocatable :: fld_tmp(:,:)
real(kind_real) :: fld_c0(geom%nc0)
character(len=1024) :: subr = 'model_mpas_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(geom%nmg,geom%nl0))

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
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,il0),(/1,nam%levs(il0),nam%timeslot(its)/),(/geom%nmg,1,1/)))
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,il0),(/1,nam%timeslot(its)/),(/geom%nmg,1/)))
      end if
   end if

   ! Global to local
   do il0=1,geom%nl0
      if (mpl%main) then
         do ic0=1,geom%nc0
            img = geom%c0_to_mg(ic0)
            fld_c0(ic0) = real(fld_tmp(img,il0),kind_real)
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

end subroutine model_mpas_read

end module model_mpas
