!----------------------------------------------------------------------
! Module: type_msv
! Purpose: deal with missing values
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_msv

use tools_kinds,only: kind_real

implicit none

type msv_type
   ! Missing values
   integer :: vali         ! Missing value for integers
   real(kind_real) :: valr ! Missing value for reals
contains
   procedure :: init => msv_init
   procedure :: msv_isnoti_0d
   procedure :: msv_isnoti_1d
   procedure :: msv_isnoti_2d
   procedure :: msv_isnoti_3d
   procedure :: msv_isnoti_4d
   procedure :: msv_isnoti_5d
   generic :: isnoti => msv_isnoti_0d,msv_isnoti_1d,msv_isnoti_2d,msv_isnoti_3d,msv_isnoti_4d,msv_isnoti_5d
   procedure :: msv_isi_0d
   procedure :: msv_isi_1d
   procedure :: msv_isi_2d
   procedure :: msv_isi_3d
   procedure :: msv_isi_4d
   procedure :: msv_isi_5d
   generic :: isi => msv_isi_0d,msv_isi_1d,msv_isi_2d,msv_isi_3d,msv_isi_4d,msv_isi_5d
   procedure :: msv_isanynoti_1d
   procedure :: msv_isanynoti_2d
   procedure :: msv_isanynoti_3d
   procedure :: msv_isanynoti_4d
   procedure :: msv_isanynoti_5d
   generic :: isanynoti => msv_isanynoti_1d,msv_isanynoti_2d,msv_isanynoti_3d,msv_isanynoti_4d,msv_isanynoti_5d
   procedure :: msv_isanyi_1d
   procedure :: msv_isanyi_2d
   procedure :: msv_isanyi_3d
   procedure :: msv_isanyi_4d
   procedure :: msv_isanyi_5d
   generic :: isanyi => msv_isanyi_1d,msv_isanyi_2d,msv_isanyi_3d,msv_isanyi_4d,msv_isanyi_5d
   procedure :: msv_isallnoti_1d
   procedure :: msv_isallnoti_2d
   procedure :: msv_isallnoti_3d
   procedure :: msv_isallnoti_4d
   procedure :: msv_isallnoti_5d
   generic :: isallnoti => msv_isallnoti_1d,msv_isallnoti_2d,msv_isallnoti_3d,msv_isallnoti_4d,msv_isallnoti_5d
   procedure :: msv_isalli_1d
   procedure :: msv_isalli_2d
   procedure :: msv_isalli_3d
   procedure :: msv_isalli_4d
   procedure :: msv_isalli_5d
   generic :: isalli => msv_isalli_1d,msv_isalli_2d,msv_isalli_3d,msv_isalli_4d,msv_isalli_5d
   procedure :: msv_isnotr_0d
   procedure :: msv_isnotr_1d
   procedure :: msv_isnotr_2d
   procedure :: msv_isnotr_3d
   procedure :: msv_isnotr_4d
   procedure :: msv_isnotr_5d
   generic :: isnotr => msv_isnotr_0d,msv_isnotr_1d,msv_isnotr_2d,msv_isnotr_3d,msv_isnotr_4d,msv_isnotr_5d
   procedure :: msv_isr_0d
   procedure :: msv_isr_1d
   procedure :: msv_isr_2d
   procedure :: msv_isr_3d
   procedure :: msv_isr_4d
   procedure :: msv_isr_5d
   generic :: isr => msv_isr_0d,msv_isr_1d,msv_isr_2d,msv_isr_3d,msv_isr_4d,msv_isr_5d
   procedure :: msv_isanynotr_1d
   procedure :: msv_isanynotr_2d
   procedure :: msv_isanynotr_3d
   procedure :: msv_isanynotr_4d
   procedure :: msv_isanynotr_5d
   generic :: isanynotr => msv_isanynotr_1d,msv_isanynotr_2d,msv_isanynotr_3d,msv_isanynotr_4d,msv_isanynotr_5d
   procedure :: msv_isanyr_1d
   procedure :: msv_isanyr_2d
   procedure :: msv_isanyr_3d
   procedure :: msv_isanyr_4d
   procedure :: msv_isanyr_5d
   generic :: isanyr => msv_isanyr_1d,msv_isanyr_2d,msv_isanyr_3d,msv_isanyr_4d,msv_isanyr_5d
   procedure :: msv_isallnotr_1d
   procedure :: msv_isallnotr_2d
   procedure :: msv_isallnotr_3d
   procedure :: msv_isallnotr_4d
   procedure :: msv_isallnotr_5d
   generic :: isallnotr => msv_isallnotr_1d,msv_isallnotr_2d,msv_isallnotr_3d,msv_isallnotr_4d,msv_isallnotr_5d
   procedure :: msv_isallr_1d
   procedure :: msv_isallr_2d
   procedure :: msv_isallr_3d
   procedure :: msv_isallr_4d
   procedure :: msv_isallr_5d
   generic :: isallr => msv_isallr_1d,msv_isallr_2d,msv_isallr_3d,msv_isallr_4d,msv_isallr_5d
end type msv_type

private
public :: msv_type

contains

!----------------------------------------------------------------------
! Subroutine: msv_init
! Purpose: initialize missing values object
!----------------------------------------------------------------------
subroutine msv_init(msv,msvali,msvalr)

implicit none

! Passed variables
class(msv_type),intent(inout) :: msv ! Missing values
integer,intent(in) :: msvali         ! Missing value for integers
real(kind_real),intent(in) :: msvalr ! Missing value for reals

! Copy missing values
msv%vali = msvali
msv%valr = msvalr

end subroutine msv_init

!----------------------------------------------------------------------
! Function: msv_isnoti_0d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function msv_isnoti_0d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i           ! Integer

! Returned value
logical :: msv_isnoti_0d

msv_isnoti_0d = abs(i-msv%vali)>0

end function msv_isnoti_0d

!----------------------------------------------------------------------
! Function: msv_isnoti_1d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function msv_isnoti_1d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:) ! Integer

! Returned value
logical :: msv_isnoti_1d(max(size(i),1))

if (size(i)>0) then
   msv_isnoti_1d = abs(i-msv%vali)>0
else
   msv_isnoti_1d = .false.
end if

end function msv_isnoti_1d

!----------------------------------------------------------------------
! Function: msv_isnoti_2d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function msv_isnoti_2d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:) ! Integer

! Returned value
logical :: msv_isnoti_2d(max(size(i,1),1),max(size(i,2),1))

if (size(i)>0) then
   msv_isnoti_2d = abs(i-msv%vali)>0
else
   msv_isnoti_2d = .false.
end if

end function msv_isnoti_2d

!----------------------------------------------------------------------
! Function: msv_isnoti_3d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function msv_isnoti_3d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:) ! Integer

! Returned value
logical :: msv_isnoti_3d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1))

if (size(i)>0) then
   msv_isnoti_3d = abs(i-msv%vali)>0
else
   msv_isnoti_3d = .false.
end if

end function msv_isnoti_3d

!----------------------------------------------------------------------
! Function: msv_isnoti_4d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function msv_isnoti_4d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:,:) ! Integer

! Returned value
logical :: msv_isnoti_4d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1))

if (size(i)>0) then
   msv_isnoti_4d = abs(i-msv%vali)>0
else
   msv_isnoti_4d = .false.
end if

end function msv_isnoti_4d

!----------------------------------------------------------------------
! Function: msv_isnoti_5d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function msv_isnoti_5d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:,:,:) ! Integer

! Returned value
logical :: msv_isnoti_5d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1),max(size(i,5),1))

if (size(i)>0) then
   msv_isnoti_5d = abs(i-msv%vali)>0
else
   msv_isnoti_5d = .false.
end if

end function msv_isnoti_5d

!----------------------------------------------------------------------
! Function: msv_isi_0d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function msv_isi_0d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i           ! Integer

! Returned value
logical :: msv_isi_0d

msv_isi_0d = .not.msv%isnoti(i)

end function msv_isi_0d

!----------------------------------------------------------------------
! Function: msv_isi_1d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function msv_isi_1d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:)        ! Integer

! Returned value
logical :: msv_isi_1d(max(size(i),1))

msv_isi_1d = .not.msv%isnoti(i)

end function msv_isi_1d

!----------------------------------------------------------------------
! Function: msv_isi_2d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function msv_isi_2d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:)      ! Integer

! Returned value
logical :: msv_isi_2d(max(size(i,1),1),max(size(i,2),1))

msv_isi_2d = .not.msv%isnoti(i)

end function msv_isi_2d

!----------------------------------------------------------------------
! Function: msv_isi_3d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function msv_isi_3d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:)    ! Integer

! Returned value
logical :: msv_isi_3d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1))

msv_isi_3d = .not.msv%isnoti(i)

end function msv_isi_3d

!----------------------------------------------------------------------
! Function: msv_isi_4d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function msv_isi_4d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:,:)  ! Integer

! Returned value
logical :: msv_isi_4d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1))

msv_isi_4d = .not.msv%isnoti(i)

end function msv_isi_4d

!----------------------------------------------------------------------
! Function: msv_isi_5d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function msv_isi_5d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
integer,intent(in) :: i(:,:,:,:,:) ! Integer

! Returned value
logical :: msv_isi_5d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1),max(size(i,5),1))

msv_isi_5d = .not.msv%isnoti(i)

end function msv_isi_5d

!----------------------------------------------------------------------
! Function: msv_isanynoti_1d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynoti_1d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:)        ! Integer

msv_isanynoti_1d = any(msv%isnoti(i))

end function msv_isanynoti_1d

!----------------------------------------------------------------------
! Function: msv_isanynoti_2d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynoti_2d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:)      ! Integer

msv_isanynoti_2d = any(msv%isnoti(i))

end function msv_isanynoti_2d

!----------------------------------------------------------------------
! Function: msv_isanynoti_3d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynoti_3d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:)    ! Integer

msv_isanynoti_3d = any(msv%isnoti(i))

end function msv_isanynoti_3d

!----------------------------------------------------------------------
! Function: msv_isanynoti_4d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynoti_4d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:,:)  ! Integer

msv_isanynoti_4d = any(msv%isnoti(i))

end function msv_isanynoti_4d

!----------------------------------------------------------------------
! Function: msv_isanynoti_5d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynoti_5d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
integer,intent(in) :: i(:,:,:,:,:) ! Integer

msv_isanynoti_5d = any(msv%isnoti(i))

end function msv_isanynoti_5d

!----------------------------------------------------------------------
! Function: msv_isanyi_1d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyi_1d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:)        ! Integer

msv_isanyi_1d = .not.msv%isallnoti(i)

end function msv_isanyi_1d

!----------------------------------------------------------------------
! Function: msv_isanyi_2d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyi_2d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:)      ! Integer

msv_isanyi_2d = .not.msv%isallnoti(i)

end function msv_isanyi_2d

!----------------------------------------------------------------------
! Function: msv_isanyi_3d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyi_3d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:)    ! Integer

msv_isanyi_3d = .not.msv%isallnoti(i)

end function msv_isanyi_3d

!----------------------------------------------------------------------
! Function: msv_isanyi_4d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyi_4d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:,:)  ! Integer

msv_isanyi_4d = .not.msv%isallnoti(i)

end function msv_isanyi_4d

!----------------------------------------------------------------------
! Function: msv_isanyi_5d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyi_5d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
integer,intent(in) :: i(:,:,:,:,:) ! Integer

msv_isanyi_5d = .not.msv%isallnoti(i)

end function msv_isanyi_5d

!----------------------------------------------------------------------
! Function: msv_isallnoti_1d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnoti_1d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:)        ! Integer

msv_isallnoti_1d = all(msv%isnoti(i))

end function msv_isallnoti_1d

!----------------------------------------------------------------------
! Function: msv_isallnoti_2d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnoti_2d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:)      ! Integer

msv_isallnoti_2d = all(msv%isnoti(i))

end function msv_isallnoti_2d

!----------------------------------------------------------------------
! Function: msv_isallnoti_3d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnoti_3d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:)    ! Integer

msv_isallnoti_3d = all(msv%isnoti(i))

end function msv_isallnoti_3d

!----------------------------------------------------------------------
! Function: msv_isallnoti_4d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnoti_4d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:,:)  ! Integer

msv_isallnoti_4d = all(msv%isnoti(i))

end function msv_isallnoti_4d

!----------------------------------------------------------------------
! Function: msv_isallnoti_5d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnoti_5d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
integer,intent(in) :: i(:,:,:,:,:) ! Integer

msv_isallnoti_5d = all(msv%isnoti(i))

end function msv_isallnoti_5d

!----------------------------------------------------------------------
! Function: msv_isalli_1d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isalli_1d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:)        ! Integer

msv_isalli_1d = .not.msv%isanynoti(i)

end function msv_isalli_1d

!----------------------------------------------------------------------
! Function: msv_isalli_2d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isalli_2d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:)      ! Integer

msv_isalli_2d = .not.msv%isanynoti(i)

end function msv_isalli_2d

!----------------------------------------------------------------------
! Function: msv_isalli_3d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isalli_3d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:)    ! Integer

msv_isalli_3d = .not.msv%isanynoti(i)

end function msv_isalli_3d

!----------------------------------------------------------------------
! Function: msv_isalli_4d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isalli_4d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
integer,intent(in) :: i(:,:,:,:)  ! Integer

msv_isalli_4d = .not.msv%isanynoti(i)

end function msv_isalli_4d

!----------------------------------------------------------------------
! Function: msv_isalli_5d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function msv_isalli_5d(msv,i)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
integer,intent(in) :: i(:,:,:,:,:) ! Integer

msv_isalli_5d = .not.msv%isanynoti(i)

end function msv_isalli_5d

!----------------------------------------------------------------------
! Function: msv_isnotr_0d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function msv_isnotr_0d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
real(kind_real),intent(in) :: r   ! Real number

! Returned value
logical :: msv_isnotr_0d

msv_isnotr_0d = abs(r-msv%valr)>0.0

end function msv_isnotr_0d

!----------------------------------------------------------------------
! Function: msv_isnotr_1d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function msv_isnotr_1d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
real(kind_real),intent(in) :: r(:) ! Real number

! Returned value
logical :: msv_isnotr_1d(size(r))

if (size(r)>0) then
   msv_isnotr_1d = abs(r-msv%valr)>0.0
else
   msv_isnotr_1d = .false.
end if

end function msv_isnotr_1d

!----------------------------------------------------------------------
! Function: msv_isnotr_2d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function msv_isnotr_2d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv    ! Missing values
real(kind_real),intent(in) :: r(:,:) ! Real number

! Returned value
logical :: msv_isnotr_2d(max(size(r,1),1),max(size(r,2),1))

if (size(r)>0) then
   msv_isnotr_2d = abs(r-msv%valr)>0.0
else
   msv_isnotr_2d = .false.
end if

end function msv_isnotr_2d

!----------------------------------------------------------------------
! Function: msv_isnotr_3d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function msv_isnotr_3d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv      ! Missing values
real(kind_real),intent(in) :: r(:,:,:) ! Real number

! Returned value
logical :: msv_isnotr_3d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1))

if (size(r)>0) then
   msv_isnotr_3d = abs(r-msv%valr)>0.0
else
   msv_isnotr_3d = .false.
end if

end function msv_isnotr_3d

!----------------------------------------------------------------------
! Function: msv_isnotr_4d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function msv_isnotr_4d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv        ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:) ! Real number

! Returned value
logical :: msv_isnotr_4d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1))

if (size(r)>0) then
   msv_isnotr_4d = abs(r-msv%valr)>0.0
else
   msv_isnotr_4d = .false.
end if

end function msv_isnotr_4d

!----------------------------------------------------------------------
! Function: msv_isnotr_5d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function msv_isnotr_5d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv          ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real number

! Returned value
logical :: msv_isnotr_5d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1),max(size(r,5),1))

if (size(r)>0) then
   msv_isnotr_5d = abs(r-msv%valr)>0.0
else
   msv_isnotr_5d = .false.
end if

end function msv_isnotr_5d

!----------------------------------------------------------------------
! Function: msv_isr_0d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function msv_isr_0d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv ! Missing values
real(kind_real),intent(in) :: r   ! Real number

! Returned value
logical :: msv_isr_0d

msv_isr_0d = .not.msv%isnotr(r)

end function msv_isr_0d

!----------------------------------------------------------------------
! Function: msv_isr_1d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function msv_isr_1d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
real(kind_real),intent(in) :: r(:) ! Real number

! Returned value
logical :: msv_isr_1d(size(r))

msv_isr_1d = .not.msv%isnotr(r)

end function msv_isr_1d

!----------------------------------------------------------------------
! Function: msv_isr_2d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function msv_isr_2d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv    ! Missing values
real(kind_real),intent(in) :: r(:,:) ! Real number

! Returned value
logical :: msv_isr_2d(max(size(r,1),1),max(size(r,2),1))

msv_isr_2d = .not.msv%isnotr(r)

end function msv_isr_2d

!----------------------------------------------------------------------
! Function: msv_isr_3d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function msv_isr_3d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv      ! Missing values
real(kind_real),intent(in) :: r(:,:,:) ! Real number

! Returned value
logical :: msv_isr_3d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1))

msv_isr_3d = .not.msv%isnotr(r)

end function msv_isr_3d

!----------------------------------------------------------------------
! Function: msv_isr_4d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function msv_isr_4d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv        ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:) ! Real number

! Returned value
logical :: msv_isr_4d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1))

msv_isr_4d = .not.msv%isnotr(r)

end function msv_isr_4d

!----------------------------------------------------------------------
! Function: msv_isr_5d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function msv_isr_5d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv          ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real number

! Returned value
logical :: msv_isr_5d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1),max(size(r,5),1))

msv_isr_5d = .not.msv%isnotr(r)

end function msv_isr_5d

!----------------------------------------------------------------------
! Function: msv_isanynotr_1d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynotr_1d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
real(kind_real),intent(in) :: r(:) ! Real

msv_isanynotr_1d = any(msv%isnotr(r))

end function msv_isanynotr_1d

!----------------------------------------------------------------------
! Function: msv_isanynotr_2d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynotr_2d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv    ! Missing values
real(kind_real),intent(in) :: r(:,:) ! Real

msv_isanynotr_2d = any(msv%isnotr(r))

end function msv_isanynotr_2d

!----------------------------------------------------------------------
! Function: msv_isanynotr_3d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynotr_3d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv      ! Missing values
real(kind_real),intent(in) :: r(:,:,:) ! Real

msv_isanynotr_3d = any(msv%isnotr(r))

end function msv_isanynotr_3d

!----------------------------------------------------------------------
! Function: msv_isanynotr_4d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynotr_4d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv        ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

msv_isanynotr_4d = any(msv%isnotr(r))

end function msv_isanynotr_4d

!----------------------------------------------------------------------
! Function: msv_isanynotr_5d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isanynotr_5d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv          ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

msv_isanynotr_5d = any(msv%isnotr(r))

end function msv_isanynotr_5d

!----------------------------------------------------------------------
! Function: msv_isanyr_1d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyr_1d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
real(kind_real),intent(in) :: r(:) ! Real

msv_isanyr_1d = .not.msv%isallr(r)

end function msv_isanyr_1d

!----------------------------------------------------------------------
! Function: msv_isanyr_2d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyr_2d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv    ! Missing values
real(kind_real),intent(in) :: r(:,:) ! Real

msv_isanyr_2d = .not.msv%isallr(r)

end function msv_isanyr_2d

!----------------------------------------------------------------------
! Function: msv_isanyr_3d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyr_3d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv      ! Missing values
real(kind_real),intent(in) :: r(:,:,:) ! Real

msv_isanyr_3d = .not.msv%isallr(r)

end function msv_isanyr_3d

!----------------------------------------------------------------------
! Function: msv_isanyr_4d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyr_4d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv        ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

msv_isanyr_4d = .not.msv%isallr(r)

end function msv_isanyr_4d

!----------------------------------------------------------------------
! Function: msv_isanyr_5d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isanyr_5d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv          ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

msv_isanyr_5d = .not.msv%isallr(r)

end function msv_isanyr_5d

!----------------------------------------------------------------------
! Function: msv_isallnotr_1d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnotr_1d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
real(kind_real),intent(in) :: r(:) ! Real

msv_isallnotr_1d = all(msv%isnotr(r))

end function msv_isallnotr_1d

!----------------------------------------------------------------------
! Function: msv_isallnotr_2d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnotr_2d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv    ! Missing values
real(kind_real),intent(in) :: r(:,:) ! Real

msv_isallnotr_2d = all(msv%isnotr(r))

end function msv_isallnotr_2d

!----------------------------------------------------------------------
! Function: msv_isallnotr_3d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnotr_3d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv      ! Missing values
real(kind_real),intent(in) :: r(:,:,:) ! Real

msv_isallnotr_3d = all(msv%isnotr(r))

end function msv_isallnotr_3d

!----------------------------------------------------------------------
! Function: msv_isallnotr_4d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnotr_4d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv        ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

msv_isallnotr_4d = all(msv%isnotr(r))

end function msv_isallnotr_4d

!----------------------------------------------------------------------
! Function: msv_isallnotr_5d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function msv_isallnotr_5d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv          ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

msv_isallnotr_5d = all(msv%isnotr(r))

end function msv_isallnotr_5d

!----------------------------------------------------------------------
! Function: msv_isallr_1d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isallr_1d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv  ! Missing values
real(kind_real),intent(in) :: r(:) ! Real

msv_isallr_1d = .not.msv%isanynotr(r)

end function msv_isallr_1d

!----------------------------------------------------------------------
! Function: msv_isallr_2d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isallr_2d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv    ! Missing values
real(kind_real),intent(in) :: r(:,:) ! Real

msv_isallr_2d = .not.msv%isanynotr(r)

end function msv_isallr_2d

!----------------------------------------------------------------------
! Function: msv_isallr_3d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isallr_3d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv      ! Missing values
real(kind_real),intent(in) :: r(:,:,:) ! Real

msv_isallr_3d = .not.msv%isanynotr(r)

end function msv_isallr_3d

!----------------------------------------------------------------------
! Function: msv_isallr_4d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isallr_4d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv        ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

msv_isallr_4d = .not.msv%isanynotr(r)

end function msv_isallr_4d

!----------------------------------------------------------------------
! Function: msv_isallr_5d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function msv_isallr_5d(msv,r)

implicit none

! Passed variables
class(msv_type),intent(in) :: msv          ! Missing values
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

msv_isallr_5d = .not.msv%isanynotr(r)

end function msv_isallr_5d

end module type_msv
