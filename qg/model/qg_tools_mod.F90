! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_tools_mod

use fckit_configuration_module, only: fckit_configuration
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
function genfilename(f_conf,length,vdate,date_cols)
use string_utils

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
integer,intent(in) :: length                   !< Length
type(datetime),intent(in) :: vdate             !< Date and time
logical,intent(in) :: date_cols                !< Date written with colons or not

! Result
character(len=2*length) :: genfilename

! Local variables
integer :: lenfn
character(len=length) :: fdbdir,expver,typ,validitydate,referencedate,sstep,mmb,iter
character(len=2*length) :: prefix
character(len=:),allocatable :: str

type(datetime) :: rdate
type(duration) :: step

! Get configuration parameters
if (f_conf%has("datadir")) then
  call f_conf%get_or_die("datadir",str)
  call swap_name_member(f_conf, str, 6)
  fdbdir = str
else
  fdbdir = "."
endif
call f_conf%get_or_die("exp",str)
call swap_name_member(f_conf, str, 6)
expver = str
call f_conf%get_or_die("type",str)
typ = str

! Ensemble case
if (typ=='ens') then
  call f_conf%get_or_die("member",str)
  mmb = str
  lenfn = len_trim(fdbdir) + 1 + len_trim(expver) + 1 + len_trim(typ) + 1 + len_trim(mmb)
  prefix = trim(fdbdir) // '/' // trim(expver) // '.' // trim(typ) // '.' // trim(mmb)
else
  lenfn = len_trim(fdbdir) + 1 + len_trim(expver) + 1 + len_trim(typ)
  prefix = trim(fdbdir) // '/' // trim(expver) // '.' // trim(typ)
endif

! Forecast / ensemble cases
if ((typ=='fc').or.(typ=='ens')) then
  call f_conf%get_or_die("date",str)
  referencedate = str
  if (date_cols) then
    call datetime_create(trim(referencedate),rdate)
    call datetime_diff(vdate,rdate,step)
    call duration_to_string(step,sstep)
  else
    call datetime_create(trim(referencedate),rdate)
    call datetime_to_string_io(rdate,referencedate)
    call datetime_diff(vdate,rdate,step)
    call duration_to_string(step,sstep)
  endif
  call datetime_delete(rdate)

  lenfn = lenfn+1+len_trim(referencedate)+1+len_trim(sstep)
  genfilename = trim(prefix)//'.'//trim(referencedate)//'.'// trim(sstep)//'.nc'
endif

! Analysis, diagnostic or increment case
if ((typ=='an').or.(typ=='diag').or.(typ=='in')) then
  if (date_cols) then
    call datetime_to_string(vdate,validitydate)
  else
    call datetime_to_string_io(vdate,validitydate)
  endif
  lenfn = lenfn+1+len_trim(validitydate)
  genfilename = trim(prefix)//'.'//trim(validitydate)//'.nc'
endif

if (typ=='krylov') then
  call f_conf%get_or_die("iteration",str)
  iter = str
  if (date_cols) then
    call datetime_to_string(vdate,validitydate)
  else
    call datetime_to_string_io(vdate,validitydate)
  endif
  lenfn = lenfn+1+len_trim(iter)+1+len_trim(validitydate)
  genfilename = trim(prefix)//'.'// trim(iter)//'.'//trim(validitydate)//'.nc'
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
