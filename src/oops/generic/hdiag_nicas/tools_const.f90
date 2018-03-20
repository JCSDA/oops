!----------------------------------------------------------------------
! Module: tools_const
!> Purpose: usual constants
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_const

use tools_kinds, only: kind_real

implicit none

! Constants
real(kind_real),parameter :: pi = acos(-1.0)    !< Pi
real(kind_real),parameter :: deg2rad = pi/180.0 !< Degree to radian
real(kind_real),parameter :: rad2deg = 180.0/pi !< Radian to degree
real(kind_real),parameter :: req = 6371229.0    !< Earth radius (m)
real(kind_real),parameter :: reqkm = 6371.229   !< Earth radius (km)
real(kind_real),parameter :: ps = 101325.0      !< Reference surface pressure

private
public :: pi,deg2rad,rad2deg,req,reqkm,ps

end module tools_const
