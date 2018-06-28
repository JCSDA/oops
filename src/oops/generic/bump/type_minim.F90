!----------------------------------------------------------------------
! Module: type_minim
!> Purpose: minimization data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_minim

use tools_const, only: rth
use tools_func, only: fit_diag,fit_lct
use tools_kinds, only: kind_real
use tools_missing, only: isnotmsr
use type_mpl, only: mpl_type

implicit none

real(kind_real),parameter :: rho = 0.5_kind_real    !< Convergence parameter for the Hooke algorithm
real(kind_real),parameter :: tol = 1.0e-6_kind_real !< Tolerance for the Hooke algorithm
integer,parameter :: itermax = 10                   !< Maximum number of iteration for the Hooke algorithm

! Minimization data derived type
type minim_type
   ! Generic data
   integer :: nx                              !< Control vector size
   integer :: ny                              !< Function output size
   real(kind_real),allocatable :: x(:)        !< Control vector
   real(kind_real),allocatable :: guess(:)    !< Control vector guess
   real(kind_real),allocatable :: norm(:)     !< Control vector norm
   real(kind_real),allocatable :: binf(:)     !< Control vector lower bound
   real(kind_real),allocatable :: bsup(:)     !< Control vector upper bound
   real(kind_real),allocatable :: obs(:)      !< Observation
   character(len=1024) :: cost_function       !< Cost function
   real(kind_real) :: f_guess                 !< Guess cost
   real(kind_real) :: f_min                   !< Minimum cost
   character(len=1024) :: algo                !< Minimization algorithm

   ! Common data
   integer :: nl0                             !< Number of levels
   integer :: nc3                             !< Number of classes

   ! Specific data (fit)
   integer :: nl0r                            !< Reduced number of levels
   logical :: lhomh                           !< Vertically homogenous horizontal support radius key
   logical :: lhomv                           !< Vertically homogenous vertical support radius key
   integer,allocatable :: l0rl0_to_l0(:,:)    !< Reduced level to level
   real(kind_real),allocatable :: distvr(:,:) !< Vertical distance
   real(kind_real),allocatable :: disth(:)    !< Horizontal distance

   ! Specific data (LCT)
   integer :: nscales                         !< Number of LCT scales
   integer,allocatable :: ncomp(:)            !< Number of LCT components
   real(kind_real),allocatable :: dx(:,:)     !< Zonal separation
   real(kind_real),allocatable :: dy(:,:)     !< Meridian separation
   real(kind_real),allocatable :: dz(:)       !< Vertical separation
   logical,allocatable :: dmask(:,:)          !< Mask
contains
   procedure :: compute => minim_compute
   procedure :: cost => minim_cost
   procedure :: cost_fit_diag => minim_cost_fit_diag
   procedure :: cost_fit_lct => minim_cost_fit_lct
   procedure :: hooke => minim_hooke
   procedure :: best_nearby => minim_best_nearby
end type minim_type

private
public :: minim_type

contains

!----------------------------------------------------------------------
! subroutine: minim_compute
!> Purpose: minimize ensuring bounds constraints
!----------------------------------------------------------------------
subroutine minim_compute(minim,mpl,lprt)

implicit none

! Passed variables
class(minim_type),intent(inout) :: minim !< Minimization data
type(mpl_type),intent(in) :: mpl         !< MPI data
logical,intent(in) :: lprt               !< Print key

! Local variables
integer :: ix
real(kind_real) :: guess(minim%nx)

! Check
if (minim%nx<=0) call mpl%abort('nx should be positive to minimize')
if (minim%ny<=0) call mpl%abort('nx should be positive to minimize')

! Initialization
do ix=1,minim%nx
   if (abs(minim%norm(ix))>0.0) then
      guess(ix) = minim%guess(ix)/minim%norm(ix)
   else
      guess(ix) = 0.0
   end if
end do

! Initial cost
call minim%cost(mpl,guess,minim%f_guess)

select case (trim(minim%algo))
case ('hooke')
   ! Hooke algorithm
   call minim%hooke(mpl,guess)
end select

! Final cost
call minim%cost(mpl,minim%x,minim%f_min)

! Test
if (minim%f_min<minim%f_guess) then
   minim%x = minim%x*minim%norm
   if (lprt) then
      write(mpl%unit,'(a13,a,f6.1,a)') '','Minimizer '//trim(minim%algo)//', cost function decrease:', &
    & abs(minim%f_min-minim%f_guess)/minim%f_guess*100.0,'%'
      call flush(mpl%unit)
   end if
else
   minim%x = minim%guess
   if (lprt) call mpl%warning('Minimizer '//trim(minim%algo)//' failed')
end if

end subroutine minim_compute

!----------------------------------------------------------------------
! subroutine: minim_cost
!> Purpose: compute cost function
!----------------------------------------------------------------------
subroutine minim_cost(minim,mpl,x,f)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim     !< Minimization data
type(mpl_type),intent(in) :: mpl          !< MPI data
real(kind_real),intent(in) :: x(minim%nx) !< Control vector
real(kind_real),intent(out) :: f          !< Cost function value

select case (minim%cost_function)
case ('fit_diag')
   call minim%cost_fit_diag(mpl,x,f)
case ('fit_lct')
   call minim%cost_fit_lct(mpl,x,f)
end select

end subroutine minim_cost

!----------------------------------------------------------------------
! Function: minim_cost_fit_diag
!> Purpose: diagnosic fit function cost
!----------------------------------------------------------------------
subroutine minim_cost_fit_diag(minim,mpl,x,f)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim     !< Minimization data
type(mpl_type),intent(in) :: mpl          !< MPI data
real(kind_real),intent(in) :: x(minim%nx) !< Control vector
real(kind_real),intent(out) :: f          !< Cost function value

! Local variables
integer :: offset,ix
real(kind_real) :: fit_rh(minim%nl0),fit_rv(minim%nl0)
real(kind_real) :: fit(minim%nc3,minim%nl0r,minim%nl0)
real(kind_real) :: xtmp(minim%nx),fit_pack(minim%ny),xx

! Renormalize
xtmp = x*minim%norm

! Get data
offset = 0
if (minim%lhomh) then
   fit_rh = xtmp(offset+1)
   offset = offset+1
else
   fit_rh = xtmp(offset+1:offset+minim%nl0)
   offset = offset+minim%nl0
end if
if (minim%lhomv) then
   fit_rv = xtmp(offset+1)
   offset = offset+1
else
   fit_rv = xtmp(offset+1:offset+minim%nl0)
   offset = offset+minim%nl0
end if

! Compute function
call fit_diag(mpl,minim%nc3,minim%nl0r,minim%nl0,minim%l0rl0_to_l0,minim%disth,minim%distvr,fit_rh,fit_rv,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Cost
f = sum((minim%obs-fit_pack)**2,mask=isnotmsr(minim%obs).and.isnotmsr(fit_pack))

! Bound penalty
do ix=1,minim%nx
   if (minim%bsup(ix)>minim%binf(ix)) then
      xx = (xtmp(ix)-minim%binf(ix))/(minim%bsup(ix)-minim%binf(ix))
      if (xx<0.0) then
         f = f+minim%f_guess*xx**2
      elseif (xx>1.0) then
         f = f+minim%f_guess*(xx-1.0)**2
      end if
   end if
end do

end subroutine minim_cost_fit_diag

!----------------------------------------------------------------------
! Function: minim_cost_fit_lct
!> Purpose: LCT fit function cost
!----------------------------------------------------------------------
subroutine minim_cost_fit_lct(minim,mpl,x,f)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim     !< Minimization data
type(mpl_type),intent(in) :: mpl          !< MPI data
real(kind_real),intent(in) :: x(minim%nx) !< Control vector
real(kind_real),intent(out) :: f          !< Cost function value

! Local variables
integer :: ix
real(kind_real) :: fit(minim%nc3,minim%nl0)
real(kind_real) :: xtmp(minim%nx),fit_pack(minim%ny),xx,coef(minim%nscales)

! Renormalize
xtmp = x*minim%norm

! Compute function
if (minim%nscales>1) then
   coef(1:minim%nscales-1) = xtmp(sum(minim%ncomp)+1:sum(minim%ncomp)+minim%nscales-1)
   coef(minim%nscales) = 1.0-sum(coef(1:minim%nscales-1))
else
   coef(1) = 1.0
end if
call fit_lct(mpl,minim%nc3,minim%nl0,minim%dx,minim%dy,minim%dz,minim%dmask,minim%nscales,minim%ncomp, &
 & xtmp(1:sum(minim%ncomp)),coef,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Cost
f = sum((minim%obs-fit_pack)**2,mask=isnotmsr(minim%obs).and.isnotmsr(fit_pack))

! Bound penalty
do ix=1,minim%nx
   if (minim%bsup(ix)>minim%binf(ix)) then
      xx = (xtmp(ix)-minim%binf(ix))/(minim%bsup(ix)-minim%binf(ix))
      if (xx<0.0) then
         f = f+minim%f_guess*xx**2
      elseif (xx>1.0) then
         f = f+minim%f_guess*(xx-1.0)**2
      end if
   end if
end do

end subroutine minim_cost_fit_lct

!----------------------------------------------------------------------
! Subroutine: minim_hooke
!> Purpose: seeks a minimizer of a scalar function of several variables
!> Author: ALGOL original by Arthur Kaupe, C version by Mark Johnson, FORTRAN90 version by John Burkardt
!----------------------------------------------------------------------
subroutine minim_hooke(minim,mpl,guess)

implicit none

! Passed variables
class(minim_type),intent(inout) :: minim      !< Minimization data
type(mpl_type),intent(in) :: mpl              !< MPI data
real(kind_real),intent(in) :: guess(minim%nx) !< Guess

! Local variables
integer :: funevals,i,iters,keep
real(kind_real) :: fbefore,newf,steplength,tmp
real(kind_real) :: delta(minim%nx),newx(minim%nx)

! Initialization
newx = guess
minim%x = guess
do i=1,minim%nx
   if (abs(guess(i))>rth) then
      delta(i) = rho*abs(guess(i))
   else
      delta(i) = rho
   end if
end do
funevals = 0
steplength = rho
iters = 0
call minim%cost(mpl,newx,fbefore)
funevals = funevals + 1
newf = fbefore

! Iterative search
do while ((iters<itermax).and.(abs(tol-steplength)>rth*abs(tol+steplength)).and.(tol<steplength))
   ! Update iteration
   iters = iters + 1

   ! Find best new point, one coordinate at a time
   newx = minim%x
   call minim%best_nearby(mpl,delta,newx,fbefore,funevals,newf)

   ! If we made some improvements, pursue that direction
   keep = 1

   do while ((abs(newf-fbefore)>rth*abs(newf+fbefore)).and.(newf<fbefore).and.(keep==1))
      do i=1,minim%nx
         ! Arrange the sign of delta
         if ((abs(newx(i)-minim%x(i))>rth*abs(newx(i)+minim%x(i))).and.(newx(i)>minim%x(i))) then
            delta(i) = abs(delta(i))
         else
            delta(i) = -abs(delta(i))
         end if

         ! Now, move further in this direction.
         tmp = minim%x(i)
         minim%x(i) = newx(i)
         newx(i) = newx(i)+newx(i)-tmp
      end do

      ! Update
      fbefore = newf
      call minim%best_nearby(mpl,delta,newx,fbefore,funevals,newf)

      ! If the further (optimistic) move was bad...
      if ((abs(newf-fbefore)>rth*abs(newf+fbefore)).and.(fbefore<newf)) exit

      ! Make sure that the differences between the new and the old points
      ! are due to actual displacements; beware of roundoff errors that
      ! might cause NEWF < FBEFORE.
      keep = 0

      do i=1,minim%nx
         if ((abs(0.5*abs(delta(i))-abs(newx(i)-minim%x(i)))>rth*abs(0.5*abs(delta(i))+abs(newx(i)-minim%x(i)))) &
       & .and.(0.5*abs(delta(i))<abs(newx(i)-minim%x(i)))) then
            keep = 1
            exit
         end if
      end do
   end do

   if ((abs(tol-steplength)>rth*abs(tol+steplength)).and.(.not.(tol>steplength)) &
 & .and.((abs(newf-fbefore)<tiny(1.0)).or.(abs(newf-fbefore)>rth*abs(newf+fbefore))).and.(.not.(fbefore>newf))) then
      steplength = steplength*rho
      delta = delta*rho
   end if
end do

end subroutine minim_hooke

!----------------------------------------------------------------------
! Subroutine: minim_best_nearby
!> Purpose: looks for a better nearby point, one coordinate at a time
!> Author: ALGOL original by Arthur Kaupe, C version by Mark Johnson, FORTRAN90 version by John Burkardt
!----------------------------------------------------------------------
subroutine minim_best_nearby(minim,mpl,delta,point,prevbest,funevals,minf)

implicit none

! Passed variables
class(minim_type),intent(inout) :: minim         !< Minimization data
type(mpl_type),intent(in) :: mpl                 !< MPI data
real(kind_real),intent(inout) :: delta(minim%nx) !< Step
real(kind_real),intent(inout) :: point(minim%nx) !< Point
real(kind_real),intent(in) :: prevbest           !< Best existing cost
integer,intent(inout) :: funevals                !< Number of evaluations
real(kind_real),intent(out) :: minf              !< Minimum cost

! Local variables
integer :: i
real(kind_real) :: ftmp
real(kind_real) :: z(minim%nx)

! Initialization
minf = prevbest
z = point

do i=1,minim%nx
   z(i) = point(i)+delta(i)
   call minim%cost(mpl,z,ftmp)
   funevals = funevals+1
   if ((abs(ftmp-minf)>rth*abs(ftmp+minf)).and.(ftmp<minf)) then
      minf = ftmp
   else
      delta(i) = -delta(i)
      z(i) = point(i)+delta(i)
      call minim%cost(mpl,z,ftmp)
      funevals = funevals+1
      if ((abs(ftmp-minf)>rth*abs(ftmp+minf)).and.(ftmp<minf)) then
         minf = ftmp
      else
         z(i) = point(i)
      end if
   end if
end do

! Update
point = z

end subroutine minim_best_nearby

end module type_minim
