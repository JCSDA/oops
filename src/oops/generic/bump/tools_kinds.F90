!----------------------------------------------------------------------
! Module: tools_kinds
! Purpose: kinds definition
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_kinds

use,intrinsic :: iso_c_binding

implicit none

! Real kind
integer,parameter :: kind_real = c_double ! Real kind

private
public kind_real

end module tools_kinds
