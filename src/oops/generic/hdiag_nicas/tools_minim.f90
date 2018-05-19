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
use type_mdata, only: mdatatype
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
subroutine minim(mdata,func,lprt)

implicit none

! Passed variables
type(mdatatype),intent(inout) :: mdata       !< Minimization data
interface
   subroutine func(mdata,x,f)                !< Cost function
   use tools_kinds, only: kind_real
   use type_mdata, only: mdatatype
   type(mdatatype),intent(in) :: mdata       !< Minimization data
   real(kind_real),intent(in) :: x(mdata%nx) !< Control vector
   real(kind_real),intent(out) :: f          !< Cost function value
   end subroutine
end interface
logical,intent(in) :: lprt                   !< Print key

! Local variables
integer :: icount,numres,info,ix
real(kind_real) :: guess(mdata%nx),xmin(mdata%nx),y,ynewlo,step(mdata%nx)
real(kind_real) :: delta_init,h0

! Initialization
do ix=1,mdata%nx
   if (abs(mdata%norm(ix))>0.0) then
      guess(ix) = mdata%guess(ix)/mdata%norm(ix)
   else
      guess(ix) = 0.0
   end if
end do

! Initial cost
mdata%f_guess = 0.0
call func(mdata,guess,y)
mdata%f_guess = y

select case (trim(mdata%fit_type))
case ('nelder_mead')
   ! Initialization
   step = 0.1

   ! Nelder-Mead algorithm
   call nelmin(mdata,func,mdata%nx,guess,xmin,ynewlo,reqmin,step,konvge,kcount,icount,numres,info)
case ('compass_search')
   ! Initialization
   delta_init = 0.1

   ! Compass search
   call compass_search(mdata,func,mdata%nx,guess,delta_tol,delta_init,k_max,xmin,ynewlo,icount)
case ('praxis')
   ! Initialization
   h0 = 0.1
   xmin = guess

   ! Praxis
   ynewlo = praxis(mdata,func,t0,h0,mdata%nx,0,xmin)
end select

! Test
if (ynewlo<y) then
   mdata%x = xmin*mdata%norm
   if (lprt) write(mpl%unit,'(a7,a,f6.1,a)') '','Minimizer '//trim(mdata%fit_type)//', cost function decrease:', &
 & abs(ynewlo-y)/y*100.0,'%'
else
   mdata%x = mdata%guess
   if (lprt) call msgwarning('Minimizer '//trim(mdata%fit_type)//' failed')
end if

end subroutine minim

end module tools_minim
