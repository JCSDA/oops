!----------------------------------------------------------------------
! Module: type_cv_blk
! Purpose: control vector derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_cv_blk

use tools_kinds, only: kind_real

implicit none

! Control vector block derived type
type cv_blk_type
   integer :: n                            ! Control variable block size
   real(kind_real),allocatable :: alpha(:) ! Control vector block field
end type cv_blk_type

private
public :: cv_blk_type

end module type_cv_blk
