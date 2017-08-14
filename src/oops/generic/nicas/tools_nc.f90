!----------------------------------------------------------------------
! Module: tools_nc
!> Purpose: NetCDF routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_nc

use netcdf
use tools_display, only: msgerror
use tools_kinds, only: kind_real

implicit none

! NetCDF type for real
integer,parameter :: ncfloat = nf90_double

private
public :: ncfloat
public :: ncerr

contains

!----------------------------------------------------------------------
! Subroutine: ncerr
!> Purpose: handle NetCDF error
!----------------------------------------------------------------------
subroutine ncerr(subr,info)

implicit none

! Passed variables
character(len=*),intent(in) :: subr !< Calling subroutine
integer,intent(in) :: info          !< Info index

! Check status
if (info/=nf90_noerr) call msgerror('in '//trim(subr)//': '//trim(nf90_strerror(info)))

end subroutine ncerr

end module tools_nc
