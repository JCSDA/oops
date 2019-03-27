!----------------------------------------------------------------------
! Module: tools_kinds
! Purpose: kinds definition
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_kinds

use iso_c_binding
use netcdf, only: nf90_double

implicit none

integer,parameter :: kind_real = c_double       ! Real kind
integer,parameter :: nc_kind_real = nf90_double ! NetCDF real kind

private
public kind_real,nc_kind_real

end module tools_kinds
