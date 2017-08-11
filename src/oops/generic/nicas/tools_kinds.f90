!----------------------------------------------------------------------
! Module: tools_kinds
!> Purpose: kinds definition
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_kinds

use, intrinsic :: iso_c_binding
implicit none

! Real kind
integer, parameter :: kind_real = c_double !< Real kind

private
public kind_real

end module tools_kinds
