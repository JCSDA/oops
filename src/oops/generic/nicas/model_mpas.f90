!----------------------------------------------------------------------
! Module: module_mpas.f90
!> Purpose: MPAS model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_mpas

use module_namelist, only: namtype
use netcdf
use tools_const, only: pi,deg2rad,rad2deg
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msi,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geomtype,geom_alloc

implicit none

private
public :: model_mpas_coord,model_mpas_read,model_mpas_write

contains

!----------------------------------------------------------------------
! Subroutine: model_mpas_coord
!> Purpose: get MPAS coordinates
!----------------------------------------------------------------------
subroutine model_mpas_coord(nam,geom)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Local variables
integer :: ncid,nc0_id,nlev_id,lon_id,lat_id,pres_id
real(kind=4),allocatable :: lon(:),lat(:),pres(:)
character(len=1024) :: subr = 'model_mpas_coord'

! Open file and get dimensions
call msi(geom%nlon)
call msi(geom%nlat)
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'nCells',nc0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc0_id,len=geom%nc0))
call ncerr(subr,nf90_inq_dimid(ncid,'nVertLevels',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nc0))
allocate(lat(geom%nc0))
allocate(pres(geom%nlev))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'lonCell',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'latCell',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'pressure_base',pres_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call ncerr(subr,nf90_close(ncid))

! Pack
call geom_alloc(geom)
geom%lon = real(lon,kind_real)
geom%lat = real(lat,kind_real)
geom%mask = .true.

! Compute normalized area
geom%area = 4.0*pi

! Vertical unit
if (nam%logpres) then
   geom%vunit = log(pres(nam%levs(1:geom%nl0)))
else
   geom%vunit = float(nam%levs(1:geom%nl0))
end if

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(pres)

end subroutine model_mpas_coord

!----------------------------------------------------------------------
! Subroutine: model_mpas_read
!> Purpose: read MPAS field
!----------------------------------------------------------------------
subroutine model_mpas_read(nam,ncid,varname,geom,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
integer,intent(in) :: ncid                              !< NetCDF file ID
character(len=*),intent(in) :: varname                  !< Variable name
type(geomtype),intent(in) :: geom                     !< Sampling data
real(kind_real),intent(out) :: fld(geom%nc0,geom%nl0) !< Read field

! Local variables
integer :: il0
integer :: fld_id
real(kind=4) :: fld_loc(geom%nc0)
character(len=1024) :: subr = 'model_mpas_read'

! Initialize field
call msr(fld)

! Get variable id
call ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

! Read variable
do il0=1,geom%nl0
   call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/nam%levs(il0),1,1/),(/1,geom%nc0,1/)))
   fld(:,il0) = real(fld_loc,kind_real)
end do

end subroutine model_mpas_read

!----------------------------------------------------------------------
! Subroutine: model_mpas_write
!> Purpose: write MPAS field
!----------------------------------------------------------------------
subroutine model_mpas_write(ncid,varname,geom,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                             !< NetCDF file ID
character(len=*),intent(in) :: varname                 !< Variable name
type(geomtype),intent(in) :: geom                    !< Sampling data
real(kind_real),intent(in) :: fld(geom%nc0,geom%nl0) !< Written field

! Local variables
integer :: il0,ierr
integer :: nc0_id,nlev_id,nt_id,fld_id,lon_id,lat_id
character(len=1024) :: subr = 'model_mpas_write'

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'nc0',nc0_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
   ierr = nf90_inq_dimid(ncid,'nVertLevels',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nVertLevels',geom%nl0,nlev_id))
   ierr = nf90_inq_dimid(ncid,'Time',nt_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Time',1,nt_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlev_id,nc0_id,nt_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,geom%nl0
   if (isanynotmsr(fld(:,il0))) then
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld(:,il0),(/il0,1,1/),(/1,geom%nc0,1/)))
   end if
end do

! Write coordinates
ierr = nf90_inq_varid(ncid,'lonCell',lon_id)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   call ncerr(subr,nf90_def_var(ncid,'lonCell',ncfloat,(/nc0_id/),lon_id))
   call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
   call ncerr(subr,nf90_def_var(ncid,'latCell',ncfloat,(/nc0_id/),lat_id))
   call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
   call ncerr(subr,nf90_enddef(ncid))
   call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
end if

end subroutine model_mpas_write

end module model_mpas
