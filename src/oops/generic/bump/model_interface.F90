!----------------------------------------------------------------------
! Module: model_interface
!> Purpose: model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module model_interface

use model_aro, only: model_aro_coord,model_aro_read
use model_arp, only: model_arp_coord,model_arp_read
use model_gem, only: model_gem_coord,model_gem_read
use model_geos, only: model_geos_coord,model_geos_read
use model_gfs, only: model_gfs_coord,model_gfs_read
use model_ifs, only: model_ifs_coord,model_ifs_read
use model_mpas, only: model_mpas_coord,model_mpas_read
use model_nemo, only: model_nemo_coord,model_nemo_read
use model_wrf, only: model_wrf_coord,model_wrf_read
use netcdf
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr,isnotmsi
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

private
public :: model_coord,model_read

contains

!----------------------------------------------------------------------
! Subroutine: model_coord
!> Purpose: get coordinates
!----------------------------------------------------------------------
subroutine model_coord(mpl,rng,nam,geom)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl      !< MPI data
type(rng_type),intent(inout) :: rng   !< Random number generator
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
if (trim(nam%model)=='aro') call model_aro_coord(mpl,nam,geom)
if (trim(nam%model)=='arp') call model_arp_coord(mpl,nam,geom)
if (trim(nam%model)=='gem') call model_gem_coord(mpl,nam,geom)
if (trim(nam%model)=='geos') call model_geos_coord(mpl,nam,geom)
if (trim(nam%model)=='gfs') call model_gfs_coord(mpl,nam,geom)
if (trim(nam%model)=='ifs') call model_ifs_coord(mpl,nam,geom)
if (trim(nam%model)=='mpas') call model_mpas_coord(mpl,nam,geom)
if (trim(nam%model)=='nemo') call model_nemo_coord(mpl,nam,geom)
if (trim(nam%model)=='online') call mpl%abort('online model should not call model_coord')
if (trim(nam%model)=='wrf') call model_wrf_coord(mpl,nam,geom)

! Define distribution
call geom%define_distribution(mpl,nam,rng)

end subroutine model_coord

!----------------------------------------------------------------------
! Subroutine: model_read
!> Purpose: read model field
!----------------------------------------------------------------------
subroutine model_read(mpl,nam,geom,filename,ie,jsub,fld)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl                                      !< MPI data
type(nam_type),intent(in) :: nam                                      !< Namelist
type(geom_type),intent(in) :: geom                                    !< Geometry
character(len=*),intent(in) :: filename                               !< File name
integer,intent(in) :: ie                                              !< Ensemble member index
integer,intent(in) :: jsub                                            !< Sub-ensemble index
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its
character(len=1024) :: fullname

! Initialization
call msr(fld)

do its=1,nam%nts
   ! Define filename
   if (jsub==0) then
      write(fullname,'(a,a,i2.2,a,i4.4,a)') trim(filename),'_',nam%timeslot(its),'_',ie,'.nc'
   else
      write(fullname,'(a,a,i2.2,a,i4.4,a,i4.4,a)') trim(filename),'_',nam%timeslot(its),'_',jsub,'_',ie,'.nc'
   end if

   ! Select model
   if (trim(nam%model)=='aro') call model_aro_read(mpl,nam,geom,fullname,fld(:,:,:,its))
   if (trim(nam%model)=='arp') call model_arp_read(mpl,nam,geom,fullname,fld(:,:,:,its))
   if (trim(nam%model)=='gem') call model_gem_read(mpl,nam,geom,fullname,fld(:,:,:,its))
   if (trim(nam%model)=='geos') call model_geos_read(mpl,nam,geom,fullname,its,fld(:,:,:,its))
   if (trim(nam%model)=='gfs') call model_gfs_read(mpl,nam,geom,fullname,fld(:,:,:,its))
   if (trim(nam%model)=='ifs') call model_ifs_read(mpl,nam,geom,fullname,its,fld(:,:,:,its))
   if (trim(nam%model)=='mpas') call model_mpas_read(mpl,nam,geom,fullname,its,fld(:,:,:,its))
   if (trim(nam%model)=='nemo') call model_nemo_read(mpl,nam,geom,fullname,its,fld(:,:,:,its))
   if (trim(nam%model)=='wrf') call model_wrf_read(mpl,nam,geom,fullname,its,fld(:,:,:,its))
end do

end subroutine model_read

end module model_interface
