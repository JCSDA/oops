!----------------------------------------------------------------------
! Module: tools_const
! Purpose: define usual constants and missing values
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_const

use tools_kinds, only: kind_real

implicit none

! Constants
real(kind_real),parameter :: pi = acos(-1.0_kind_real)    ! Pi
real(kind_real),parameter :: deg2rad = pi/180.0_kind_real ! Degree to radian
real(kind_real),parameter :: rad2deg = 180.0_kind_real/pi ! Radian to degree
real(kind_real),parameter :: req = 6371229.0_kind_real    ! Earth radius
real(kind_real),parameter :: reqkm = 6371.229_kind_real   ! Earth radius
real(kind_real),parameter :: ps = 101325.0_kind_real      ! Reference surface pressure

! Missing values
integer,parameter :: msvali = -999                        ! Integer missing value
real(kind_real),parameter :: msvalr = -999.0_kind_real    ! Real missing value

private
public :: pi,deg2rad,rad2deg,req,reqkm,ps,msvali,msvalr

end module tools_const
