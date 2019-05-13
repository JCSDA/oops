! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_tools_mod

use config_mod
use datetime_mod
use duration_mod
use iso_c_binding
use kinds
use netcdf
use qg_constants_mod

implicit none

private
public :: genfilename,ncerr,baroclinic_instability,large_vortices
! ------------------------------------------------------------------------------
real(kind_real),parameter :: ubot = -2.0_kind_real !< Zonal wind at the surface (m/s)
real(kind_real),parameter :: utop = 58.0_kind_real !< Zonal wind at the top (m/s)
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Generate filename
function genfilename(conf,length,vdate)

implicit none

! Passed variables
type(c_ptr),intent(in) :: conf     !< Configuration
integer,intent(in) :: length       !< Length
type(datetime),intent(in) :: vdate !< Date and time

! Result
character(len=length) :: genfilename

! Local variables
integer :: lenfn
character(len=length) :: fdbdir,expver,typ,validitydate,referencedate,sstep,prefix,mmb
type(datetime) :: rdate
type(duration) :: step

! Get configuration parameters
fdbdir = config_get_string(conf,len(fdbdir),'datadir')
expver = config_get_string(conf,len(expver),'exp')
typ = config_get_string(conf,len(typ),'type')

! Ensemble case
if (typ=='ens') then
  mmb = config_get_string(conf,len(mmb),'member')
  lenfn = len_trim(fdbdir) + 1 + len_trim(expver) + 1 + len_trim(typ) + 1 + len_trim(mmb)
  prefix = trim(fdbdir) // '/' // trim(expver) // '.' // trim(typ) // '.' // trim(mmb)
else
  lenfn = len_trim(fdbdir) + 1 + len_trim(expver) + 1 + len_trim(typ)
  prefix = trim(fdbdir) // '/' // trim(expver) // '.' // trim(typ)
endif

! Forecast / ensemble cases
if ((typ=='fc').or.(typ=='ens')) then
  referencedate = config_get_string(conf,len(referencedate),'date')
  call datetime_to_string(vdate,validitydate)
  call datetime_create(trim(referencedate),rdate)
  call datetime_diff(vdate,rdate,step)
  call duration_to_string(step,sstep)
  lenfn = lenfn+1+len_trim(referencedate)+1+len_trim(sstep)
  genfilename = trim(prefix)//'.'//trim(referencedate)//'.'// trim(sstep)
endif

! Analysis case
if (typ=='an') then
  call datetime_to_string(vdate,validitydate)
  lenfn = lenfn+1+len_trim(validitydate)
  genfilename = trim(prefix)//'.'//trim(validitydate)
endif

! Check filename length
if (lenfn>length) call abor1_ftn('genfilename: filename too long')

end function genfilename
! ------------------------------------------------------------------------------
!> Check NetCDF status
subroutine ncerr(info)

implicit none

! Passed variables
integer,intent(in) :: info !< Info index

! Check status
if (info/=nf90_noerr) call abor1_ftn(trim(nf90_strerror(info)))

end subroutine ncerr
! ------------------------------------------------------------------------------
!> Generate values for baroclinic instability
subroutine baroclinic_instability(x,y,z,var,res)

implicit none

! Passed variables
real(kind_real),intent(in) :: x    !< X value
real(kind_real),intent(in) :: y    !< Y value
real(kind_real),intent(in) :: z    !< Z value
character(len=1),intent(in) :: var !< Variable
real(kind_real),intent(out) :: res !< Results

! Local variable
real(kind_real) :: u

! Define zonal wind
u = ubot+(utop-ubot)*z/domain_depth

select case (var)
case ('x')
  ! Streamfunction
  res = -u*y
case ('q')
  call abor1_ftn('baroclinic_instability: cannot define q')
case ('u')
  ! Zonal wind
  res = u
case ('v')
  ! Meridional wind
  res = 0.0
case default
  call abor1_ftn('baroclinic_instability: wrong variable')
end select

end subroutine baroclinic_instability
! ------------------------------------------------------------------------------
!> Generate values for large vortices
subroutine large_vortices(x,y,z,var,res)

implicit none

! Passed variables
real(kind_real),intent(in) :: x    !< X value
real(kind_real),intent(in) :: y    !< Y value
real(kind_real),intent(in) :: z    !< Z value
character(len=1),intent(in) :: var !< Variable
real(kind_real),intent(out) :: res !< Results

! Local variable
real(kind_real) :: ff

! Define wind speed
ff = ubot+(utop-ubot)*z/domain_depth

select case (var)
case ('x')
  ! Streamfunction
  res = ff*domain_meridional*cos(2.0*pi*x/domain_zonal)*sin(pi*y/domain_meridional)
case ('q')
  call abor1_ftn('large_vortices: cannot define q')
case ('u')
  ! Zonal wind
  res = -ff*pi*cos(2.0*pi*x/domain_zonal)*cos(pi*y/domain_meridional)
case ('v')
  ! Meridional wind
  res = -2.0*ff*pi*domain_meridional/domain_zonal*sin(2.0*pi*x/domain_zonal)*sin(pi*y/domain_meridional)
case default
  call abor1_ftn('large_vortices: wrong variable')
end select

end subroutine large_vortices
! ------------------------------------------------------------------------------
end module qg_tools_mod
