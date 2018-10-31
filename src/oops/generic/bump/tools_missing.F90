!----------------------------------------------------------------------
! Module: tools_missing
! Purpose: deal with missing values
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_missing

use tools_const, only: msvali,msvalr
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
interface ismsi
  module procedure ismsi_0d
  module procedure ismsi_1d
  module procedure ismsi_2d
  module procedure ismsi_3d
  module procedure ismsi_4d
  module procedure ismsi_5d
end interface
interface isnotmsr
  module procedure isnotmsr_0d
  module procedure isnotmsr_1d
  module procedure isnotmsr_2d
  module procedure isnotmsr_3d
  module procedure isnotmsr_4d
  module procedure isnotmsr_5d
end interface
interface ismsr
  module procedure ismsr_0d
  module procedure ismsr_1d
  module procedure ismsr_2d
  module procedure ismsr_3d
  module procedure ismsr_4d
  module procedure ismsr_5d
end interface
interface isanynotmsi
  module procedure isanynotmsi_1d
  module procedure isanynotmsi_2d
  module procedure isanynotmsi_3d
  module procedure isanynotmsi_4d
  module procedure isanynotmsi_5d
end interface
interface isanymsi
  module procedure isanymsi_1d
  module procedure isanymsi_2d
  module procedure isanymsi_3d
  module procedure isanymsi_4d
  module procedure isanymsi_5d
end interface
interface isanynotmsr
  module procedure isanynotmsr_1d
  module procedure isanynotmsr_2d
  module procedure isanynotmsr_3d
  module procedure isanynotmsr_4d
  module procedure isanynotmsr_5d
end interface
interface isanymsr
  module procedure isanymsr_1d
  module procedure isanymsr_2d
  module procedure isanymsr_3d
  module procedure isanymsr_4d
  module procedure isanymsr_5d
end interface
interface isallnotmsi
  module procedure isallnotmsi_1d
  module procedure isallnotmsi_2d
  module procedure isallnotmsi_3d
  module procedure isallnotmsi_4d
  module procedure isallnotmsi_5d
end interface
interface isallmsi
  module procedure isallmsi_1d
  module procedure isallmsi_2d
  module procedure isallmsi_3d
  module procedure isallmsi_4d
  module procedure isallmsi_5d
end interface
interface isallnotmsr
  module procedure isallnotmsr_1d
  module procedure isallnotmsr_2d
  module procedure isallnotmsr_3d
  module procedure isallnotmsr_4d
  module procedure isallnotmsr_5d
end interface
interface isallmsr
  module procedure isallmsr_1d
  module procedure isallmsr_2d
  module procedure isallmsr_3d
  module procedure isallmsr_4d
  module procedure isallmsr_5d
end interface

private
public :: msi,msr,isnotmsi,ismsi,isnotmsr,ismsr,isanynotmsi,isanymsi,isanynotmsr,isanymsr,isallnotmsi,isallmsi,isallnotmsr,isallmsr

contains

!----------------------------------------------------------------------
! Subroutine: msi_0d
! Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_0d(i)

implicit none

! Passed variables
integer,intent(out),optional :: i ! Integer

i = msvali

end subroutine msi_0d

!----------------------------------------------------------------------
! Subroutine: msi_1d
! Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_1d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:) ! Integer

i = msvali

end subroutine msi_1d

!----------------------------------------------------------------------
! Subroutine: msi_2d
! Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_2d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:) ! Integer

i = msvali

end subroutine msi_2d

!----------------------------------------------------------------------
! Subroutine: msi_3d
! Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_3d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:) ! Integer

i = msvali

end subroutine msi_3d

!----------------------------------------------------------------------
! Subroutine: msi_4d
! Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_4d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:,:) ! Integer

i = msvali

end subroutine msi_4d

!----------------------------------------------------------------------
! Subroutine: msi_5d
! Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_5d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:,:,:) ! Integer

i = msvali

end subroutine msi_5d

!----------------------------------------------------------------------
! Subroutine: msi_6d
! Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_6d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:,:,:,:) ! Integer

i = msvali

end subroutine msi_6d

!----------------------------------------------------------------------
! Subroutine: msr_0d
! Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_0d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r ! Real number

r = msvalr

end subroutine msr_0d

!----------------------------------------------------------------------
! Subroutine: msr_1d
! Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:) ! Real number

r = msvalr

end subroutine msr_1d

!----------------------------------------------------------------------
! Subroutine: msr_2d
! Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:) ! Real number

r = msvalr

end subroutine msr_2d

!----------------------------------------------------------------------
! Subroutine: msr_3d
! Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:) ! Real number

r = msvalr

end subroutine msr_3d

!----------------------------------------------------------------------
! Subroutine: msr_4d
! Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:,:) ! Real number

r = msvalr

end subroutine msr_4d

!----------------------------------------------------------------------
! Subroutine: msr_5d
! Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:,:,:) ! Real number

r = msvalr

end subroutine msr_5d

!----------------------------------------------------------------------
! Subroutine: msr_6d
! Purpose: set real number to missing value
!----------------------------------------------------------------------
subroutine msr_6d(r)

implicit none

! Passed variables
real(kind_real),intent(out) :: r(:,:,:,:,:,:) ! Real number

r = msvalr

end subroutine msr_6d

!----------------------------------------------------------------------
! Function: isnotmsi_0d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_0d(i)

implicit none

! Passed variables
integer,intent(in) :: i ! Integer

! Returned value
logical :: isnotmsi_0d

isnotmsi_0d = abs(i-msvali)>0

end function isnotmsi_0d

!----------------------------------------------------------------------
! Function: isnotmsi_1d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) ! Integer

! Returned value
logical :: isnotmsi_1d(max(size(i),1))

if (size(i)>0) then
   isnotmsi_1d = abs(i-msvali)>0
else
   isnotmsi_1d = .false.
end if

end function isnotmsi_1d

!----------------------------------------------------------------------
! Function: isnotmsi_2d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) ! Integer

! Returned value
logical :: isnotmsi_2d(max(size(i,1),1),max(size(i,2),1))

if (size(i)>0) then
   isnotmsi_2d = abs(i-msvali)>0
else
   isnotmsi_2d = .false.
end if

end function isnotmsi_2d

!----------------------------------------------------------------------
! Function: isnotmsi_3d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) ! Integer

! Returned value
logical :: isnotmsi_3d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1))

if (size(i)>0) then
   isnotmsi_3d = abs(i-msvali)>0
else
   isnotmsi_3d = .false.
end if

end function isnotmsi_3d

!----------------------------------------------------------------------
! Function: isnotmsi_4d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) ! Integer

! Returned value
logical :: isnotmsi_4d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1))

if (size(i)>0) then
   isnotmsi_4d = abs(i-msvali)>0
else
   isnotmsi_4d = .false.
end if

end function isnotmsi_4d

!----------------------------------------------------------------------
! Function: isnotmsi_5d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
function isnotmsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) ! Integer

! Returned value
logical :: isnotmsi_5d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1),max(size(i,5),1))

if (size(i)>0) then
   isnotmsi_5d = abs(i-msvali)>0
else
   isnotmsi_5d = .false.
end if

end function isnotmsi_5d

!----------------------------------------------------------------------
! Function: ismsi_0d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function ismsi_0d(i)

implicit none

! Passed variables
integer,intent(in) :: i ! Integer

! Returned value
logical :: ismsi_0d

ismsi_0d = .not.isnotmsi(i)

end function ismsi_0d

!----------------------------------------------------------------------
! Function: ismsi_1d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function ismsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) ! Integer

! Returned value
logical :: ismsi_1d(max(size(i),1))

ismsi_1d = .not.isnotmsi(i)

end function ismsi_1d

!----------------------------------------------------------------------
! Function: ismsi_2d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function ismsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) ! Integer

! Returned value
logical :: ismsi_2d(max(size(i,1),1),max(size(i,2),1))

ismsi_2d = .not.isnotmsi(i)

end function ismsi_2d

!----------------------------------------------------------------------
! Function: ismsi_3d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function ismsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) ! Integer

! Returned value
logical :: ismsi_3d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1))

ismsi_3d = .not.isnotmsi(i)

end function ismsi_3d

!----------------------------------------------------------------------
! Function: ismsi_4d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function ismsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) ! Integer

! Returned value
logical :: ismsi_4d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1))

ismsi_4d = .not.isnotmsi(i)

end function ismsi_4d

!----------------------------------------------------------------------
! Function: ismsi_5d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
function ismsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) ! Integer

! Returned value
logical :: ismsi_5d(max(size(i,1),1),max(size(i,2),1),max(size(i,3),1),max(size(i,4),1),max(size(i,5),1))

ismsi_5d = .not.isnotmsi(i)

end function ismsi_5d

!----------------------------------------------------------------------
! Function: isnotmsr_0d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_0d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r ! Real number

! Returned value
logical :: isnotmsr_0d

isnotmsr_0d = abs(r-msvalr)>0.0

end function isnotmsr_0d

!----------------------------------------------------------------------
! Function: isnotmsr_1d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) ! Real number

! Returned value
logical :: isnotmsr_1d(size(r))

if (size(r)>0) then
   isnotmsr_1d = abs(r-msvalr)>0.0
else
   isnotmsr_1d = .false.
end if

end function isnotmsr_1d

!----------------------------------------------------------------------
! Function: isnotmsr_2d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) ! Real number

! Returned value
logical :: isnotmsr_2d(max(size(r,1),1),max(size(r,2),1))

if (size(r)>0) then
   isnotmsr_2d = abs(r-msvalr)>0.0
else
   isnotmsr_2d = .false.
end if

end function isnotmsr_2d

!----------------------------------------------------------------------
! Function: isnotmsr_3d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) ! Real number

! Returned value
logical :: isnotmsr_3d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1))

if (size(r)>0) then
   isnotmsr_3d = abs(r-msvalr)>0.0
else
   isnotmsr_3d = .false.
end if

end function isnotmsr_3d

!----------------------------------------------------------------------
! Function: isnotmsr_4d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) ! Real number

! Returned value
logical :: isnotmsr_4d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1))

if (size(r)>0) then
   isnotmsr_4d = abs(r-msvalr)>0.0
else
   isnotmsr_4d = .false.
end if

end function isnotmsr_4d

!----------------------------------------------------------------------
! Function: isnotmsr_5d
! Purpose: check if an real number is not set to missing value
!----------------------------------------------------------------------
function isnotmsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real number

! Returned value
logical :: isnotmsr_5d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1),max(size(r,5),1))

if (size(r)>0) then
   isnotmsr_5d = abs(r-msvalr)>0.0
else
   isnotmsr_5d = .false.
end if

end function isnotmsr_5d

!----------------------------------------------------------------------
! Function: ismsr_0d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function ismsr_0d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r ! Real number

! Returned value
logical :: ismsr_0d

ismsr_0d = .not.isnotmsr(r)

end function ismsr_0d

!----------------------------------------------------------------------
! Function: ismsr_1d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function ismsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) ! Real number

! Returned value
logical :: ismsr_1d(size(r))

ismsr_1d = .not.isnotmsr(r)

end function ismsr_1d

!----------------------------------------------------------------------
! Function: ismsr_2d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function ismsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) ! Real number

! Returned value
logical :: ismsr_2d(max(size(r,1),1),max(size(r,2),1))

ismsr_2d = .not.isnotmsr(r)

end function ismsr_2d

!----------------------------------------------------------------------
! Function: ismsr_3d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function ismsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) ! Real number

! Returned value
logical :: ismsr_3d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1))

ismsr_3d = .not.isnotmsr(r)

end function ismsr_3d

!----------------------------------------------------------------------
! Function: ismsr_4d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function ismsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) ! Real number

! Returned value
logical :: ismsr_4d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1))

ismsr_4d = .not.isnotmsr(r)

end function ismsr_4d

!----------------------------------------------------------------------
! Function: ismsr_5d
! Purpose: check if an real number is set to missing value
!----------------------------------------------------------------------
function ismsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real number

! Returned value
logical :: ismsr_5d(max(size(r,1),1),max(size(r,2),1),max(size(r,3),1),max(size(r,4),1),max(size(r,5),1))

ismsr_5d = .not.isnotmsr(r)

end function ismsr_5d

!----------------------------------------------------------------------
! Function: isanynotmsi_1d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) ! Integer

isanynotmsi_1d = any(isnotmsi(i))

end function isanynotmsi_1d

!----------------------------------------------------------------------
! Function: isanynotmsi_2d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) ! Integer

isanynotmsi_2d = any(isnotmsi(i))

end function isanynotmsi_2d

!----------------------------------------------------------------------
! Function: isanynotmsi_3d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) ! Integer

isanynotmsi_3d = any(isnotmsi(i))

end function isanynotmsi_3d

!----------------------------------------------------------------------
! Function: isanynotmsi_4d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) ! Integer

isanynotmsi_4d = any(isnotmsi(i))

end function isanynotmsi_4d

!----------------------------------------------------------------------
! Function: isanynotmsi_5d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) ! Integer

isanynotmsi_5d = any(isnotmsi(i))

end function isanynotmsi_5d

!----------------------------------------------------------------------
! Function: isanymsi_1d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isanymsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) ! Integer

isanymsi_1d = .not.isallnotmsi(i)

end function isanymsi_1d

!----------------------------------------------------------------------
! Function: isanymsi_2d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isanymsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) ! Integer

isanymsi_2d = .not.isallnotmsi(i)

end function isanymsi_2d

!----------------------------------------------------------------------
! Function: isanymsi_3d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isanymsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) ! Integer

isanymsi_3d = .not.isallnotmsi(i)

end function isanymsi_3d

!----------------------------------------------------------------------
! Function: isanymsi_4d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isanymsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) ! Integer

isanymsi_4d = .not.isallnotmsi(i)

end function isanymsi_4d

!----------------------------------------------------------------------
! Function: isanymsi_5d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isanymsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) ! Integer

isanymsi_5d = .not.isallnotmsi(i)

end function isanymsi_5d

!----------------------------------------------------------------------
! Function: isanynotmsr_1d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) ! Real

isanynotmsr_1d = any(isnotmsr(r))

end function isanynotmsr_1d

!----------------------------------------------------------------------
! Function: isanynotmsr_2d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) ! Real

isanynotmsr_2d = any(isnotmsr(r))

end function isanynotmsr_2d

!----------------------------------------------------------------------
! Function: isanynotmsr_3d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) ! Real

isanynotmsr_3d = any(isnotmsr(r))

end function isanynotmsr_3d

!----------------------------------------------------------------------
! Function: isanynotmsr_4d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

isanynotmsr_4d = any(isnotmsr(r))

end function isanynotmsr_4d

!----------------------------------------------------------------------
! Function: isanynotmsr_5d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isanynotmsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

isanynotmsr_5d = any(isnotmsr(r))

end function isanynotmsr_5d

!----------------------------------------------------------------------
! Function: isanymsr_1d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isanymsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) ! Real

isanymsr_1d = .not.isallmsr(r)

end function isanymsr_1d

!----------------------------------------------------------------------
! Function: isanymsr_2d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isanymsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) ! Real

isanymsr_2d = .not.isallmsr(r)

end function isanymsr_2d

!----------------------------------------------------------------------
! Function: isanymsr_3d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isanymsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) ! Real

isanymsr_3d = .not.isallmsr(r)

end function isanymsr_3d

!----------------------------------------------------------------------
! Function: isanymsr_4d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isanymsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

isanymsr_4d = .not.isallmsr(r)

end function isanymsr_4d

!----------------------------------------------------------------------
! Function: isanymsr_5d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isanymsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

isanymsr_5d = .not.isallmsr(r)

end function isanymsr_5d

!----------------------------------------------------------------------
! Function: isallnotmsi_1d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) ! Integer

isallnotmsi_1d = all(isnotmsi(i))

end function isallnotmsi_1d

!----------------------------------------------------------------------
! Function: isallnotmsi_2d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) ! Integer

isallnotmsi_2d = all(isnotmsi(i))

end function isallnotmsi_2d

!----------------------------------------------------------------------
! Function: isallnotmsi_3d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) ! Integer

isallnotmsi_3d = all(isnotmsi(i))

end function isallnotmsi_3d

!----------------------------------------------------------------------
! Function: isallnotmsi_4d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) ! Integer

isallnotmsi_4d = all(isnotmsi(i))

end function isallnotmsi_4d

!----------------------------------------------------------------------
! Function: isallnotmsi_5d
! Purpose: check if an integer is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) ! Integer

isallnotmsi_5d = all(isnotmsi(i))

end function isallnotmsi_5d

!----------------------------------------------------------------------
! Function: isallmsi_1d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isallmsi_1d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:) ! Integer

isallmsi_1d = .not.isanynotmsi(i)

end function isallmsi_1d

!----------------------------------------------------------------------
! Function: isallmsi_2d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isallmsi_2d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:) ! Integer

isallmsi_2d = .not.isanynotmsi(i)

end function isallmsi_2d

!----------------------------------------------------------------------
! Function: isallmsi_3d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isallmsi_3d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:) ! Integer

isallmsi_3d = .not.isanynotmsi(i)

end function isallmsi_3d

!----------------------------------------------------------------------
! Function: isallmsi_4d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isallmsi_4d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:) ! Integer

isallmsi_4d = .not.isanynotmsi(i)

end function isallmsi_4d

!----------------------------------------------------------------------
! Function: isallmsi_5d
! Purpose: check if an integer is set to missing value
!----------------------------------------------------------------------
logical function isallmsi_5d(i)

implicit none

! Passed variables
integer,intent(in) :: i(:,:,:,:,:) ! Integer

isallmsi_5d = .not.isanynotmsi(i)

end function isallmsi_5d

!----------------------------------------------------------------------
! Function: isallnotmsr_1d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) ! Real

isallnotmsr_1d = all(isnotmsr(r))

end function isallnotmsr_1d

!----------------------------------------------------------------------
! Function: isallnotmsr_2d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) ! Real

isallnotmsr_2d = all(isnotmsr(r))

end function isallnotmsr_2d

!----------------------------------------------------------------------
! Function: isallnotmsr_3d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) ! Real

isallnotmsr_3d = all(isnotmsr(r))

end function isallnotmsr_3d

!----------------------------------------------------------------------
! Function: isallnotmsr_4d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

isallnotmsr_4d = all(isnotmsr(r))

end function isallnotmsr_4d

!----------------------------------------------------------------------
! Function: isallnotmsr_5d
! Purpose: check if a real is not set to missing value
!----------------------------------------------------------------------
logical function isallnotmsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

isallnotmsr_5d = all(isnotmsr(r))

end function isallnotmsr_5d

!----------------------------------------------------------------------
! Function: isallmsr_1d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isallmsr_1d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:) ! Real

isallmsr_1d = .not.isanynotmsr(r)

end function isallmsr_1d

!----------------------------------------------------------------------
! Function: isallmsr_2d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isallmsr_2d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:) ! Real

isallmsr_2d = .not.isanynotmsr(r)

end function isallmsr_2d

!----------------------------------------------------------------------
! Function: isallmsr_3d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isallmsr_3d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:) ! Real

isallmsr_3d = .not.isanynotmsr(r)

end function isallmsr_3d

!----------------------------------------------------------------------
! Function: isallmsr_4d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isallmsr_4d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:) ! Real

isallmsr_4d = .not.isanynotmsr(r)

end function isallmsr_4d

!----------------------------------------------------------------------
! Function: isallmsr_5d
! Purpose: check if a real is set to missing value
!----------------------------------------------------------------------
logical function isallmsr_5d(r)

implicit none

! Passed variables
real(kind_real),intent(in) :: r(:,:,:,:,:) ! Real

isallmsr_5d = .not.isanynotmsr(r)

end function isallmsr_5d

end module tools_missing
