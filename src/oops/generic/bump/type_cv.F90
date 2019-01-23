!----------------------------------------------------------------------
! Module: type_cv
! Purpose: control vector derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_cv

use tools_kinds, only: kind_real
use type_cv_blk, only: cv_blk_type

implicit none

! Control vector derived type
type cv_type
   integer :: n                            ! Total control variable size
   integer :: nbe                          ! Number of control variable blocks
   type(cv_blk_type),allocatable :: blk(:) ! Control variable blocks
contains
   procedure :: pack => cv_pack
   procedure :: unpack => cv_unpack
end type cv_type

private
public :: cv_type

contains

!----------------------------------------------------------------------
! Subroutine: cv_pack
! Purpose: pack
!----------------------------------------------------------------------
subroutine cv_pack(cv,pcv)

! Passed variables
class(cv_type),intent(in) :: cv          ! Control variable
real(kind_real),intent(out) :: pcv(cv%n) ! Packed control variable

! Local variable
integer :: ib,offset

! Initialization
offset = 0

do ib=1,cv%nbe
   if (cv%blk(ib)%n>0) then
      ! Pack control variable
      pcv(offset+1:offset+cv%blk(ib)%n) = cv%blk(ib)%alpha

      ! Update
      offset = offset+cv%blk(ib)%n
   end if
end do

end subroutine cv_pack

!----------------------------------------------------------------------
! Subroutine: cv_unpack
! Purpose: unpack
!----------------------------------------------------------------------
subroutine cv_unpack(cv,pcv)

! Passed variables
class(cv_type),intent(inout) :: cv      ! Control variable
real(kind_real),intent(in) :: pcv(cv%n) ! Packed control variable

! Local variable
integer :: ib,offset

! Initialization
offset = 0

do ib=1,cv%nbe
   if (cv%blk(ib)%n>0) then
      ! Unpack control variable
      cv%blk(ib)%alpha = pcv(offset+1:offset+cv%blk(ib)%n)

      ! Update
      offset = offset+cv%blk(ib)%n
   end if
end do

end subroutine cv_unpack

end module type_cv
