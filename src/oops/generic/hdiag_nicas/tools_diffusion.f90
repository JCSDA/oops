!----------------------------------------------------------------------
! Module: tools_diffusion.f90
!> Purpose: 3D diffusion function from Weaver and Mirouze (2013)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_diffusion

use tools_display, only: msgerror
use tools_kinds, only: kind_real
implicit none

private
public :: matern

contains

!----------------------------------------------------------------------
! Function: diffusion
!> Purpose: compute the normalized diffusion function from eq. (55) of Mirouze and Weaver (2013), for the 3d case (d = 3)
!----------------------------------------------------------------------
function matern(M,x)

implicit none

! Result
real(kind_real) :: matern

! Passed variables
integer,intent(in) :: M
real(kind_real),intent(in) :: x

! Local variables
integer :: j
real(kind_real) :: xtmp,beta

! Check
if (M<2) call msgerror('M should be larger than 2')
if (mod(M,2)>0) call msgerror('M should be even')

! Initialization
matern = 0.0
beta = 1.0
xtmp = x*sqrt(float(2*M-3))

do j=0,M-3
   ! Update sum
   matern = matern+beta*(xtmp)**(M-2-j)
   
   ! Update beta
   beta = beta*float((j+1+M-2)*(-j+M-2))/float(2*(j+1))
end do

! Last term and normalization
matern = matern/beta+1.0

! Exponential factor
matern = matern*exp(-xtmp)

end function matern

end module tools_diffusion
