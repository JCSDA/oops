!----------------------------------------------------------------------
! Module: tools_missing
!> Purpose: deal with missing values
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_missing

use tools_kinds,only: kind_real
implicit none

! Constants
integer,parameter :: msvali = -999            !< Integer missing value
real(kind_real),parameter :: msvalr = -999.0  !< Real missing value

interface msi
  module procedure msi_0d
  module procedure msi_1d
  module procedure msi_2d
  module procedure msi_3d
  module procedure msi_4d
  module procedure msi_5d
end interface
interface msr
  module procedure msr_0d
  module procedure msr_1d
  module procedure msr_2d
  module procedure msr_3d
  module procedure msr_4d
  module procedure msr_5d
end interface
interface isnotmsi
  module procedure isnotmsi_0d
  module procedure isnotmsi_1d
  module procedure isnotmsi_2d
  module procedure isnotmsi_3d
end interface
interface isnotmsr
  module procedure isnotmsr_0d
  module procedure isnotmsr_1d
  module procedure isnotmsr_2d
  module procedure isnotmsr_3d
end interface
interface isanynotmsi
  module procedure isanynotmsi_1d
  module procedure isanynotmsi_2d
end interface
interface isanynotmsr
  module procedure isanynotmsr_1d
  module procedure isanynotmsr_2d
end interface
interface isallnotmsi
  module procedure isallnotmsi_1d
  module procedure isallnotmsi_2d
end interface
interface isallnotmsr
  module procedure isallnotmsr_1d
  module procedure isallnotmsr_2d
end interface

private
public :: msvali,msvalr
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
! Subroutine: msi_3d
!> Purpose: set integer to missing value
!----------------------------------------------------------------------
subroutine msi_5d(i)

implicit none

! Passed variables
integer,intent(out) :: i(:,:,:,:,:) !< Integer

i = msvali

end subroutine msi_5d

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
logical :: isnotmsi_1d(size(i))

isnotmsi_1d = abs(i-msvali)>0

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
logical :: isnotmsi_2d(size(i,1),size(i,2))

isnotmsi_2d = abs(i-msvali)>0

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
logical :: isnotmsi_3d(size(i,1),size(i,2),size(i,3))

isnotmsi_3d = abs(i-msvali)>0

end function isnotmsi_3d

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

isnotmsr_0d = abs(r-msvalr)>tiny(1.0)

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

isnotmsr_1d = abs(r-msvalr)>tiny(1.0)

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
logical :: isnotmsr_2d(size(r,1),size(r,2))

isnotmsr_2d = abs(r-msvalr)>tiny(1.0)

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
logical :: isnotmsr_3d(size(r,1),size(r,2),size(r,3))

isnotmsr_3d = abs(r-msvalr)>tiny(1.0)

end function isnotmsr_3d

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

end module tools_missing
