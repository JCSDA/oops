!----------------------------------------------------------------------
! Module: tools_minim
!> Purpose: bound constrained minimization routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_minim

use tools_asa047, only: nelmin
use tools_compass_search, only: compass_search
use tools_display, only: msgwarning
use tools_kinds, only: kind_real
use tools_praxis, only: praxis
use type_min, only: mintype
use type_mpl, only: mpl

implicit none

! Minimization parameters
real(kind_real),parameter :: reqmin = 1.0e-8    !< Nelder-Mead parameter
integer,parameter :: konvge = 10                !< Nelder-Mead parameter
integer,parameter :: kcount = 1000              !< Nelder-Mead parameter
real(kind_real),parameter :: delta_tol = 1.0e-3 !< Compass search parameter
integer,parameter :: k_max = 500                !< Compass search parameter
real(kind_real),parameter :: t0 = 1.0e-3        !< PRAXIS parameter

private
public :: minim

contains

!----------------------------------------------------------------------
! subroutine: minim
!> Purpose: minimize ensuring bounds constraints
!----------------------------------------------------------------------
subroutine minim(mindata,func,lprt)

implicit none

! Passed variables
type(mintype),intent(inout) :: mindata         !< Minimization data
interface
   subroutine func(mindata,x,f)                !< Cost function
   use tools_kinds, only: kind_real
   use type_min, only: mintype
   type(mintype),intent(in) :: mindata         !< Minimization data
   real(kind_real),intent(in) :: x(mindata%nx) !< Control vector
   real(kind_real),intent(out) :: f            !< Cost function value
   end subroutine
end interface
logical,intent(in) :: lprt                     !< Print key

! Local variables
integer :: icount,numres,info,ix
real(kind_real) :: guess(mindata%nx),xmin(mindata%nx),y,ynewlo,step(mindata%nx)
real(kind_real) :: delta_init,h0

! Initialization
do ix=1,mindata%nx
   if (abs(mindata%norm(ix))>0.0) then
      guess(ix) = mindata%guess(ix)/mindata%norm(ix)
   else
      guess(ix) = 0.0
   end if
end do

! Initial cost
mindata%f_guess = 0.0
call func(mindata,guess,y)
mindata%f_guess = y

select case (trim(mindata%fit_type))
case ('nelder_mead')
   ! Initialization
   step = 0.1

   ! Nelder-Mead algorithm
   call nelmin(mindata,func,mindata%nx,guess,xmin,ynewlo,reqmin,step,konvge,kcount,icount,numres,info)
case ('compass_search')
   ! Initialization
   delta_init = 0.1

   ! Compass search
   call compass_search(mindata,func,mindata%nx,guess,delta_tol,delta_init,k_max,xmin,ynewlo,icount)
case ('praxis')
   ! Initialization
   h0 = 0.1
   xmin = guess

   ! Praxis
   ynewlo = praxis(mindata,func,t0,h0,mindata%nx,0,xmin)
end select

! Test
if (ynewlo<y) then
   mindata%x = xmin*mindata%norm
   if (lprt) write(mpl%unit,'(a7,a,f6.1,a)') '','Minimizer '//trim(mindata%fit_type)//', cost function decrease:', &
 & abs(ynewlo-y)/y*100.0,'%'
else
   mindata%x = mindata%guess
   if (lprt) call msgwarning('Minimizer '//trim(mindata%fit_type)//' failed')
end if

end subroutine minim

end module tools_minim
