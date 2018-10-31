!----------------------------------------------------------------------
! Module: tools_repro
! Purpose: reproducibility functions
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_repro

use tools_const, only: pi
use tools_kinds, only: kind_real

implicit none

real(kind_real),parameter :: rth = 1.0e-12 ! Reproducibility threshold

private
public :: rth
public :: eq,inf,infeq,sup,supeq,indist

contains

!----------------------------------------------------------------------
! Function: eq
! Purpose: equal test for reals
!----------------------------------------------------------------------
function eq(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x ! First real
real(kind_real),intent(in) :: y ! Second real

! Returned variable
logical :: eq

eq = abs(x-y)<rth

end function eq

!----------------------------------------------------------------------
! Function: inf
! Purpose: inferior test for reals
!----------------------------------------------------------------------
function inf(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x ! First real
real(kind_real),intent(in) :: y ! Second real

! Returned variable
logical :: inf

inf = (abs(x-y)>rth*abs(x+y)).and.(x<y)

end function inf

!----------------------------------------------------------------------
! Function: infeq
! Purpose: inferior or equal test for reals
!----------------------------------------------------------------------
function infeq(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x ! First real
real(kind_real),intent(in) :: y ! Second real

! Returned variable
logical :: infeq

infeq = inf(x,y).or.eq(x,y)

end function infeq

!----------------------------------------------------------------------
! Function: sup
! Purpose: superior test for reals
!----------------------------------------------------------------------
function sup(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x ! First real
real(kind_real),intent(in) :: y ! Second real

! Returned variable
logical :: sup

sup = (abs(x-y)>rth*abs(x+y)).and.(x>y)

end function sup

!----------------------------------------------------------------------
! Function: supeq
! Purpose: superior or equal test for reals
!----------------------------------------------------------------------
function supeq(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x ! First real
real(kind_real),intent(in) :: y ! Second real

! Returned variable
logical :: supeq

supeq = sup(x,y).or.eq(x,y)

end function supeq

!----------------------------------------------------------------------
! Function: indist
! Purpose: indistiguishability test
!----------------------------------------------------------------------
function indist(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x ! First real
real(kind_real),intent(in) :: y ! Second real

! Returned variable
logical :: indist

indist = abs(x)<rth*abs(y)

end function indist

end module tools_repro
