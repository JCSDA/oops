!----------------------------------------------------------------------
! Module: model_interface.f90
!> Purpose: model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_interface

use model_aro, only: model_aro_coord,model_aro_read,model_aro_write
use model_arp, only: model_arp_coord,model_arp_read,model_arp_write
use model_gem, only: model_gem_coord,model_gem_read,model_gem_write
use model_geos, only: model_geos_coord,model_geos_read,model_geos_write
use model_gfs, only: model_gfs_coord,model_gfs_read,model_gfs_write
use model_ifs, only: model_ifs_coord,model_ifs_read,model_ifs_write
use model_mpas, only: model_mpas_coord,model_mpas_read,model_mpas_write
use model_nemo, only: model_nemo_coord,model_nemo_read,model_nemo_write
use model_oops, only: model_oops_coord,model_oops_write
use model_wrf, only: model_wrf_coord,model_wrf_read,model_wrf_write
use netcdf
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msi,msr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: model_coord,model_read,load_ensemble,model_write

contains

!----------------------------------------------------------------------
! Subroutine: model_coord
!> Purpose: get coordinates
!----------------------------------------------------------------------
subroutine model_coord(nam,geom)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam      !< Namelist variables
type(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: iv

! Number of levels
geom%nl0 = nam%nl
do iv=1,nam%nv
   if (trim(nam%addvar2d(iv))/='') geom%nl0 = nam%nl+1
end do

! Select model
if (trim(nam%model)=='aro') call model_aro_coord(nam,geom)
if (trim(nam%model)=='arp') call model_arp_coord(nam,geom)
if (trim(nam%model)=='gem') call model_gem_coord(nam,geom)
if (trim(nam%model)=='geos') call model_geos_coord(nam,geom)
if (trim(nam%model)=='gfs') call model_gfs_coord(nam,geom)
if (trim(nam%model)=='ifs') call model_ifs_coord(nam,geom)
if (trim(nam%model)=='mpas') call model_mpas_coord(nam,geom)
if (trim(nam%model)=='nemo') call model_nemo_coord(nam,geom)
if (trim(nam%model)=='oops') call msgerror('OOPS should not call model_coord')
if (trim(nam%model)=='wrf') call model_wrf_coord(nam,geom)

end subroutine model_coord

!----------------------------------------------------------------------
! Subroutine: model_read
!> Purpose: read model field
!----------------------------------------------------------------------
subroutine model_read(nam,geom,filename,ie,jsub,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                                      !< Namelist
type(geom_type),intent(in) :: geom                                    !< Geometry
character(len=*),intent(in) :: filename                               !< File name
integer,intent(in) :: ie                                              !< Ensemble member index
integer,intent(in) :: jsub                                            !< Sub-ensemble index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,ncid
character(len=1024) :: fullname
character(len=1024) :: subr = 'model_read'

! Initialization
call msr(fld)

do its=1,nam%nts
   ! Define filename
   select case (trim(nam%model))
   case ('aro','arp','gem','gfs')
      if (jsub==0) then
         write(fullname,'(a,a,i2.2,a,i4.4,a)') trim(filename),'_',nam%timeslot(its),'_',ie,'.nc'
      else
         write(fullname,'(a,a,i2.2,a,i4.4,a,i4.4,a)') trim(filename),'_',nam%timeslot(its),'_',jsub,'_',ie,'.nc'
      end if
   case ('geos','ifs','mpas','nemo','wrf')
      if (jsub==0) then
         write(fullname,'(a,a,i4.4,a)') trim(filename),'_',ie,'.nc'
      else
         write(fullname,'(a,a,i4.4,a,i4.4,a)') trim(filename),'_',jsub,'_',ie,'.nc'
      end if
   case ('oops')
      call msgerror('OOPS should not call model_read')
   end select

   ! Open file
   if (mpl%main) then
      call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(fullname),nf90_nowrite,ncid))
   else
      call msi(ncid)
   end if

   ! Select model
   if (trim(nam%model)=='aro') call model_aro_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='arp') call model_arp_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='gem') call model_gem_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='geos') call model_geos_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='gfs') call model_gfs_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='ifs') call model_ifs_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='mpas') call model_mpas_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='nemo') call model_nemo_read(nam,geom,ncid,its,fld(:,:,:,its))
   if (trim(nam%model)=='wrf') call model_wrf_read(nam,geom,ncid,its,fld(:,:,:,its))

   ! Close file
   if (mpl%main) call ncerr(subr,nf90_close(ncid))
end do

end subroutine model_read

!----------------------------------------------------------------------
! Subroutine: load_ensemble
!> Purpose: load ensemble
!----------------------------------------------------------------------
subroutine load_ensemble(nam,geom,ens1)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                                                   !< Namelist
type(geom_type),intent(in) :: geom                                                 !< Geometry
real(kind_real),intent(out) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: isub,jsub,ie,ietot
real(kind_real),allocatable :: fld(:,:,:,:)

! Initialization
ietot = 0
call msr(ens1)

do isub=1,nam%ens1_nsub
   if (nam%ens1_nsub==1) then
      write(mpl%unit,'(a7,a)',advance='no') '','Full ensemble, member:'
   else
      write(mpl%unit,'(a7,a,i4,a)',advance='no') '','Sub-ensemble ',isub,', member:'
   end if

   do ie=1,nam%ens1_ne/nam%ens1_nsub
      write(mpl%unit,'(i4)',advance='no') nam%ens1_ne_offset+ie

      ! Read member
      allocate(fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))
      if (nam%ens1_nsub==1) then
         jsub = 0
      else
         jsub = isub
      end if
      call model_read(nam,geom,'ens1',nam%ens1_ne_offset+ie,jsub,fld)

      ! Copy
      ietot = ietot+1
      ens1(:,:,:,:,ietot) = fld

      ! Release memory
      deallocate(fld)
   end do
   write(mpl%unit,'(a)') ''
end do

end subroutine load_ensemble

!----------------------------------------------------------------------
! Subroutine: model_write
!> Purpose: write model field
!----------------------------------------------------------------------
subroutine model_write(nam,geom,filename,varname,fld)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam                      !< Namelist
type(geom_type),intent(in) :: geom                    !< Geometry
character(len=*),intent(in) :: filename               !< File name
character(len=*),intent(in) :: varname                !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
integer :: ic0a,ic0,jc0a,jc0,il0,info
integer :: ncid
real(kind_real) :: fld_loc(geom%nc0a,geom%nl0)
character(len=1024) :: subr = 'model_write'

! Apply mask
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      if (geom%mask(ic0,il0)) then
         fld_loc(ic0a,il0) = fld(ic0a,il0)
      else
         call msr(fld_loc(ic0a,il0))
      end if
   end do
end do

if (allocated(geom%mesh%redundant)) then
   ! Copy redundant points
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      jc0 = geom%mesh%redundant(ic0)
      if (isnotmsi(jc0)) then
         jc0a = geom%c0_to_c0a(jc0)
         fld_loc(ic0a,:) = fld_loc(jc0a,:)
      end if
   end do
end if

if (mpl%main) then
   ! Check if the file exists
   info = nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_noclobber,nf90_64bit_offset),ncid)
   if (info==nf90_noerr) then
      ! Write namelist parameters
      call nam%ncwrite(ncid)

      ! Define attribute
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))

      ! End definition mode
      call ncerr(subr,nf90_enddef(ncid))
   else
      ! Open file
      call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_write,ncid))
   end if
else
   call msi(ncid)
end if

! Select model
if (trim(nam%model)=='aro') call model_aro_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='arp') call model_arp_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='gem') call model_gem_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='geos') call model_geos_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='gfs') call model_gfs_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='ifs') call model_ifs_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='mpas') call model_mpas_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='nemo') call model_nemo_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='oops') call model_oops_write(geom,ncid,varname,fld_loc)
if (trim(nam%model)=='wrf') call model_wrf_write(geom,ncid,varname,fld_loc)

! Close file
if (mpl%main) call ncerr(subr,nf90_close(ncid))

end subroutine model_write

end module model_interface
