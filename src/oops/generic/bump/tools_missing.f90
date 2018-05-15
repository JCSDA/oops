!----------------------------------------------------------------------
! Module: tools_missing
!> Purpose: deal with missing values
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module tools_missing

use tools_const, only: msvali,msvalr
use tools_display, only: msgerror
use tools_kinds,only: kind_real

implicit none

interface msi
  module procedure msi_0d
  module procedure msi_1d
  module procedure msi_2d
  module procedure msi_3d
  module procedure msi_4d
  module procedure msi_5d
  module procedure msi_6d
end interface
interface msr
  module procedure msr_0d
  module procedure msr_1d
  module procedure msr_2d
  module procedure msr_3d
  module procedure msr_4d
  module procedure msr_5d
  module procedure msr_6d
end interface
interface isnotmsi
  module procedure isnotmsi_0d
  module procedure isnotmsi_1d
  module procedure isnotmsi_2d
  module procedure isnotmsi_3d
  module procedure isnotmsi_4d
  module procedure isnotmsi_5d
end interface
interface isnotmsr
  module procedure isnotmsr_0d
  module procedure isnotmsr_1d
  module procedure isnotmsr_2d
  module procedure isnotmsr_3d
  module procedure isnotmsr_4d
  module procedure isnotmsr_5d
end interface
interface isanynotmsi
  module procedure isanynotmsi_1d
  module procedure isanynotmsi_2d
  module procedure isanynotmsi_3d
  module procedure isanynotmsi_4d
  module procedure isanynotmsi_5d
end interface
interface isanynotmsr
  module procedure isanynotmsr_1d
  module procedure isanynotmsr_2d
  module procedure isanynotmsr_3d
  module procedure isanynotmsr_4d
  module procedure isanynotmsr_5d
end interface
interface isallnotmsi
  module procedure isallnotmsi_1d
  module procedure isallnotmsi_2d
  module procedure isallnotmsi_3d
  module procedure isallnotmsi_4d
  module procedure isallnotmsi_5d
end interface
interface isallnotmsr
  module procedure isallnotmsr_1d
  module procedure isallnotmsr_2d
  module procedure isallnotmsr_3d
  module procedure isallnotmsr_4d
  module procedure isallnotmsr_5d
end interface

private
public :: msi,msr,isnotmsi,isnotmsr,isanynotmsi,isanynotmsr,isallnotmsi,isallnotmsr

contains

!----------------------------------------------------------------------
! Subroutine: msi_0d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_0d(i)

implicit none

! Passed variables
integer,intent(out),optional :: i !< Integer

i = msvali

end subroutine msi_0d

!----------------------------------------------------------------------
! Subroutine: msi_1d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_1d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:) !< Integer

i = msvali

end subroutine msi_1d

!----------------------------------------------------------------------
! Subroutine: msi_2d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_2d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:) !< Integer

i = msvali

end subroutine msi_2d

!----------------------------------------------------------------------
! Subroutine: msi_3d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_3d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:) !< Integer

i = msvali

end subroutine msi_3d

!----------------------------------------------------------------------
! Subroutine: msi_4d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_4d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:,:) !< Integer

i = msvali

end subroutine msi_4d

!----------------------------------------------------------------------
! Subroutine: msi_5d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_5d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:,:,:) !< Integer

i = msvali

end subroutine msi_5d

!----------------------------------------------------------------------
! Subroutine: msi_6d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_6d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:,:,:,:) !< Integer

i = msvali

end subroutine msi_6d

!----------------------------------------------------------------------
! Subroutine: msr_0d
!> Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_0d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r !< Real number

r = msvalr

end subroutine msr_0d

!----------------------------------------------------------------------
! Subroutine: msr_1d
!> Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:) !< Real number

r = msvalr

end subroutine msr_1d

!----------------------------------------------------------------------
! Subroutine: msr_2d
!> Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:) !< Real number

r = msvalr

end subroutine msr_2d

!----------------------------------------------------------------------
! Subroutine: msr_3d
!> Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:) !< Real number

r = msvalr

end subroutine msr_3d

!----------------------------------------------------------------------
! Subroutine: msr_4d
!> Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:,:) !< Real number

r = msvalr

end subroutine msr_4d

!----------------------------------------------------------------------
! Subroutine: msr_5d
!> Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:,:,:) !< Real number

r = msvalr

end subroutine msr_5d

!----------------------------------------------------------------------
! Subroutine: msr_6d
!> Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_6d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:,:,:,:) !< Real number

r = msvalr

end subroutine msr_6d

!----------------------------------------------------------------------
! Function: isnotmsi_0d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_0d(i)

implicit none

! Passed variables
integer,intent(in) :: i !< Integer

! Returned value
logical :: isnotmsi_0d

isnotmsi_0d = abs(i-msvali)>0

end function isnotmsi_0d

!----------------------------------------------------------------------
! Function: isnotmsi_1d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) !< Integer

! Returned value
logical :: isnotmsi_1d(max(size(i),1))

if (size(i)>0) then
   isnotmsi_1d = abs(i-msvali)>0
else
   call msgerror('size zero for isnotmsi_1d')
end if

end function isnotmsi_1d

!----------------------------------------------------------------------
! Function: isnotmsi_2d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) !< Integer

! Returned value
logical :: isnotmsi_2d(max(size(i,1),1),max(size(i,2),1))

if (size(i)>0) then
   isnotmsi_2d = abs(i-msvali)>0
else
   call msgerror('size zero for isnotmsi_2d')
end if

end function isnotmsi_2d

!----------------------------------------------------------------------
! Function: isnotmsi_3d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) !< Integer

! Returned value
logical :: isnotmsi_3d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1))

if (size(i)>0) then
   isnotmsi_3d = abs(i-msvali)>0
else
   call msgerror('size zero for isnotmsi_3d')
end if

end function isnotmsi_3d

!----------------------------------------------------------------------
! Function: isnotmsi_4d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) !< Integer

! Returned value
logical :: isnotmsi_4d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1))

if (size(i)>0) then
   isnotmsi_4d = abs(i-msvali)>0
else
   call msgerror('size zero for isnotmsi_4d')
end if

end function isnotmsi_4d

!----------------------------------------------------------------------
! Function: isnotmsi_5d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) !< Integer

! Returned value
logical :: isnotmsi_5d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1),max(size(i,5),1))

if (size(i)>0) then
   isnotmsi_5d = abs(i-msvali)>0
else
   call msgerror('size zero for isnotmsi_5d')
end if

end function isnotmsi_5d

!----------------------------------------------------------------------
! Function: isnotmsr_0d
!> Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_0d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r !< Real number

! Returned value
logical :: isnotmsr_0d

isnotmsr_0d = abs(r-msvalr)>0.0

end function isnotmsr_0d

!----------------------------------------------------------------------
! Function: isnotmsr_1d
!> Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) !< Real number

! Returned value
logical :: isnotmsr_1d(size(r))

if (size(r)>0) then
   isnotmsr_1d = abs(r-msvalr)>0.0
else
   call msgerror('size zero for isnotmsr_1d')
end if

end function isnotmsr_1d

!----------------------------------------------------------------------
! Function: isnotmsr_2d
!> Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) !< Real number

! Returned value
logical :: isnotmsr_2d(max(size(r,1),1),max(size(r,2),1))

if (size(r)>0) then
   isnotmsr_2d = abs(r-msvalr)>0.0
else
   call msgerror('size zero for isnotmsr_2d')
end if

end function isnotmsr_2d

!----------------------------------------------------------------------
! Function: isnotmsr_3d
!> Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) !< Real number

! Returned value
logical :: isnotmsr_3d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1))

if (size(r)>0) then
   isnotmsr_3d = abs(r-msvalr)>0.0
else
   call msgerror('size zero for isnotmsr_3d')
end if

end function isnotmsr_3d

!----------------------------------------------------------------------
! Function: isnotmsr_4d
!> Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) !< Real number

! Returned value
logical :: isnotmsr_4d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1))

if (size(r)>0) then
   isnotmsr_4d = abs(r-msvalr)>0.0
else
   call msgerror('size zero for isnotmsr_4d')
end if

end function isnotmsr_4d

!----------------------------------------------------------------------
! Function: isnotmsr_5d
!> Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) !< Real number

! Returned value
logical :: isnotmsr_5d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1),max(size(r,5),1))

if (size(r)>0) then
   isnotmsr_5d = abs(r-msvalr)>0.0
else
   call msgerror('size zero for isnotmsr_5d')
end if

end function isnotmsr_5d

!----------------------------------------------------------------------
! Function: isanynotmsi_1d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) !< Integer

isanynotmsi_1d = any(abs(i-msvali)>0)

end function isanynotmsi_1d

!----------------------------------------------------------------------
! Function: isanynotmsi_2d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) !< Integer

isanynotmsi_2d = any(abs(i-msvali)>0)

end function isanynotmsi_2d

!----------------------------------------------------------------------
! Function: isanynotmsi_3d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) !< Integer

isanynotmsi_3d = any(abs(i-msvali)>0)

end function isanynotmsi_3d

!----------------------------------------------------------------------
! Function: isanynotmsi_4d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) !< Integer

isanynotmsi_4d = any(abs(i-msvali)>0)

end function isanynotmsi_4d

!----------------------------------------------------------------------
! Function: isanynotmsi_5d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) !< Integer

isanynotmsi_5d = any(abs(i-msvali)>0)

end function isanynotmsi_5d

!----------------------------------------------------------------------
! Function: isanynotmsr_1d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) !< Real

isanynotmsr_1d = any(abs(r-msvalr)>0)

end function isanynotmsr_1d

!----------------------------------------------------------------------
! Function: isanynotmsr_2d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) !< Real

isanynotmsr_2d = any(abs(r-msvalr)>0)

end function isanynotmsr_2d

!----------------------------------------------------------------------
! Function: isanynotmsr_3d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) !< Real

isanynotmsr_3d = any(abs(r-msvalr)>0)

end function isanynotmsr_3d

!----------------------------------------------------------------------
! Function: isanynotmsr_4d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) !< Real

isanynotmsr_4d = any(abs(r-msvalr)>0)

end function isanynotmsr_4d

!----------------------------------------------------------------------
! Function: isanynotmsr_5d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) !< Real

isanynotmsr_5d = any(abs(r-msvalr)>0)

end function isanynotmsr_5d

!----------------------------------------------------------------------
! Function: isallnotmsi_1d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) !< Integer

isallnotmsi_1d = all(abs(i-msvali)>0)

end function isallnotmsi_1d

!----------------------------------------------------------------------
! Function: isallnotmsi_2d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) !< Integer

isallnotmsi_2d = all(abs(i-msvali)>0)

end function isallnotmsi_2d

!----------------------------------------------------------------------
! Function: isallnotmsi_3d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) !< Integer

isallnotmsi_3d = all(abs(i-msvali)>0)

end function isallnotmsi_3d

!----------------------------------------------------------------------
! Function: isallnotmsi_4d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) !< Integer

isallnotmsi_4d = all(abs(i-msvali)>0)

end function isallnotmsi_4d

!----------------------------------------------------------------------
! Function: isallnotmsi_5d
!> Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) !< Integer

isallnotmsi_5d = all(abs(i-msvali)>0)

end function isallnotmsi_5d

!----------------------------------------------------------------------
! Function: isallnotmsr_1d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) !< Real

isallnotmsr_1d = all(abs(r-msvalr)>0)

end function isallnotmsr_1d

!----------------------------------------------------------------------
! Function: isallnotmsr_2d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) !< Real

isallnotmsr_2d = all(abs(r-msvalr)>0)

end function isallnotmsr_2d

!----------------------------------------------------------------------
! Function: isallnotmsr_3d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) !< Real

isallnotmsr_3d = all(abs(r-msvalr)>0)

end function isallnotmsr_3d

!----------------------------------------------------------------------
! Function: isallnotmsr_4d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) !< Real

isallnotmsr_4d = all(abs(r-msvalr)>0)

end function isallnotmsr_4d

!----------------------------------------------------------------------
! Function: isallnotmsr_5d
!> Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) !< Real

isallnotmsr_5d = all(abs(r-msvalr)>0)

end function isallnotmsr_5d

end module tools_missing
