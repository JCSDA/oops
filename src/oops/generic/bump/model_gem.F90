!----------------------------------------------------------------------
! Module: module_gem
!> Purpose: GEM model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module model_gem

use netcdf
use tools_const, only: pi,req,deg2rad,rad2deg,ps
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: model_gem_coord,model_gem_read

contains

!----------------------------------------------------------------------
! Subroutine: model_gem_coord
!> Purpose: get GEM coordinates
!----------------------------------------------------------------------
subroutine model_gem_coord(nam,geom)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: ic0,ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,a_id,b_id
real(kind=8),allocatable :: lon(:),lat(:),a(:),b(:)
character(len=1024) :: subr = 'model_gem_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'lon',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'lat',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
geom%nmg = geom%nlon*geom%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'lev',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon))
allocate(lat(geom%nlat))
allocate(a(geom%nlev))
allocate(b(geom%nlev))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'lon',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'lat',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'ap',a_id))
call ncerr(subr,nf90_inq_varid(ncid,'b',b_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call ncerr(subr,nf90_get_var(ncid,a_id,a))
call ncerr(subr,nf90_get_var(ncid,b_id,b))
call ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*real(deg2rad,kind=8)
lat = lat*real(deg2rad,kind=8)

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
      geom%vunit(ic0,1:nam%nl) = log(a(nam%levs(1:nam%nl))+b(nam%levs(1:nam%nl))*ps)
      if (geom%nl0>nam%nl) geom%vunit(ic0,geom%nl0) = log(ps)
   else
      geom%vunit(ic0,:) = real(nam%levs(1:geom%nl0),kind_real)
   end if
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(a)
deallocate(b)

end subroutine model_gem_coord

!----------------------------------------------------------------------
! Subroutine: model_gem_read
!> Purpose: read GEM field
!----------------------------------------------------------------------
subroutine model_gem_read(nam,geom,filename,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                             !< Namelist
type(geom_type),intent(in) :: geom                           !< Geometry
character(len=*),intent(in) :: filename                       !< File name
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) !< Field

! Local variables
integer :: iproc,iv,il0,xt,ic0a,ic0,ilon,ilat
integer :: ncid,fld_id
integer :: fld_int
real(kind_real) :: add_offset,scale_factor
character(len=1024) :: subr = 'model_gem_read'

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

         ! Check variable type
         call ncerr(subr,nf90_inquire_variable(ncid,fld_id,xtype=xt))
         select case (xt)
         case (nf90_short)
            call ncerr(subr,nf90_get_att(ncid,fld_id,'add_offset',add_offset))
            call ncerr(subr,nf90_get_att(ncid,fld_id,'scale_factor',scale_factor))
         case (nf90_double)
         case default
            call msgerror('wrong variable type')
         end select

         ! 3d variable
         do il0=1,nam%nl
            do ic0a=1,geom%nc0a
               ic0 = geom%c0a_to_c0(ic0a)
               ilon = geom%c0_to_lon(ic0)
               ilat = geom%c0_to_lat(ic0)
               if (xt==nf90_short) then
                  call ncerr(subr,nf90_get_var(ncid,fld_id,fld_int,(/ilon,ilat,nam%levs(il0)/)))
                  fld(ic0a,il0,iv) = add_offset+scale_factor*real(fld_int,kind_real)
               elseif (xt==nf90_double) then
                  call ncerr(subr,nf90_get_var(ncid,fld_id,fld(ic0a,il0,iv),(/ilon,ilat,nam%levs(il0)/)))
               end if
            end do
         end do

         if (trim(nam%addvar2d(iv))/='') then
            ! 2d variable

            ! Get id
            call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

            ! Read data
            do ic0a=1,geom%nc0a
               ic0 = geom%c0a_to_c0(ic0a)
               ilon = geom%c0_to_lon(ic0)
               ilat = geom%c0_to_lat(ic0)
               if (xt==nf90_short) then
                  call ncerr(subr,nf90_get_var(ncid,fld_id,fld_int,(/ilon,ilat/)))
                  fld(ic0a,geom%nl0,iv) = add_offset+scale_factor*real(fld_int,kind_real)
               elseif (xt==nf90_double) then
                  call ncerr(subr,nf90_get_var(ncid,fld_id,fld(ic0a,geom%nl0,iv),(/ilon,ilat/)))
               end if
            end do
         end if
      end do

      ! Close file
      call ncerr(subr,nf90_close(ncid))
   end if

   ! Wait
   call mpl%barrier
end do

end subroutine model_gem_read

end module model_gem
