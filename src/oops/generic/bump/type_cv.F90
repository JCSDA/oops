!----------------------------------------------------------------------
! Module: type_cv
!> Purpose: control vector derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_cv

use tools_kinds, only: kind_real
use type_cv_blk, only: cv_blk_type

implicit none

! Control vector derived type
type cv_type
   type(cv_blk_type),allocatable :: blk(:) !< Control variable blocks
end type cv_type

private
public :: cv_type

end module type_cv
