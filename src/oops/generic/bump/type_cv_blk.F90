!----------------------------------------------------------------------
! Module: type_cv_blk
!> Purpose: control vector derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_cv_blk

use tools_kinds, only: kind_real

implicit none

! Control vector block derived type
type cv_blk_type
   integer :: n                            !< Control variable block size
   real(kind_real),allocatable :: alpha(:) !< Control vector block field
end type cv_blk_type

private
public :: cv_blk_type

end module type_cv_blk
