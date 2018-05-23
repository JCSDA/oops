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
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: model_mpas_coord,model_mpas_read

contains

!----------------------------------------------------------------------
! Subroutine: model_mpas_coord
!> Purpose: get MPAS coordinates
!----------------------------------------------------------------------
subroutine model_mpas_coord(nam,geom)

implicit none

! Passed variables
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
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'nCells',ng_id))
call ncerr(subr,nf90_inquire_dimension(ncid,ng_id,len=geom%nmg))
call ncerr(subr,nf90_inq_dimid(ncid,'nVertLevels',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nmg))
allocate(lat(geom%nmg))
allocate(pres(geom%nlev))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'lonCell',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'latCell',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'pressure_base',pres_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call ncerr(subr,nf90_close(ncid))

! Not redundant grid
call geom%find_redundant

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
subroutine model_mpas_read(nam,geom,filename,its,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                              !< Namelist
type(geom_type),intent(in) :: geom                            !< Geometry
character(len=*),intent(in) :: filename                       !< File name
integer,intent(in) :: its                                     !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0,ic0a,ic0,iproc
integer :: ncid,fld_id
real(kind=4) :: fld_tmp
character(len=1024) :: subr = 'model_mpas_read'

! Initialize field
call msr(fld)

do iproc=1,mpl%nproc
   if (mpl%myproc==iproc) then
      ! Open file
      call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))

      do iv=1,nam%nv
         ! 3d variable

         ! Get variable id
         call ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))

         ! 3d variable
         do il0=1,nam%nl
            do ic0a=1,geom%nc0a
               ic0 = geom%c0a_to_c0(ic0a)
               call ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp,(/nam%levs(il0),ic0,nam%timeslot(its)/)))
               fld(ic0a,il0,iv) = real(fld_tmp,kind_real)
            end do
         end do

         if (trim(nam%addvar2d(iv))/='') then
            ! 2d variable

            ! Get id
            call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

            ! Read data
            do ic0a=1,geom%nc0a
               ic0 = geom%c0a_to_c0(ic0a)
               call ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp,(/ic0,nam%timeslot(its)/)))
               fld(ic0a,geom%nl0,iv) = real(fld_tmp,kind_real)
            end do
         end if
      end do

      ! Close file
      call ncerr(subr,nf90_close(ncid))
   end if

   ! Wait
   call mpl%barrier
end do

end subroutine model_mpas_read

end module model_mpas
