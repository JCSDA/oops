!----------------------------------------------------------------------
! Module: type_minim
! Purpose: minimization data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_minim

use tools_fit, only: ver_smooth
use tools_func, only: fit_diag,fit_diag_dble,fit_lct
use tools_kinds, only: kind_real
use tools_repro, only: rth,eq,inf,infeq,sup
use type_mpl, only: mpl_type

implicit none

real(kind_real),parameter :: rho = 0.5_kind_real    ! Convergence parameter for the Hooke algorithm
real(kind_real),parameter :: tol = 1.0e-6_kind_real ! Tolerance for the Hooke algorithm
integer,parameter :: itermax = 10                   ! Maximum number of iteration for the Hooke algorithm

! Minimization data derived type
type minim_type
   ! Generic data
   integer :: nx                              ! Control vector size
   integer :: ny                              ! Function output size
   real(kind_real),allocatable :: x(:)        ! Control vector
   real(kind_real),allocatable :: guess(:)    ! Control vector guess
   real(kind_real),allocatable :: binf(:)     ! Control vector lower bound
   real(kind_real),allocatable :: bsup(:)     ! Control vector upper bound
   real(kind_real),allocatable :: obs(:)      ! Observation
   character(len=1024) :: cost_function       ! Cost function
   real(kind_real) :: f_guess                 ! Guess cost
   real(kind_real) :: f_min                   ! Minimum cost
   character(len=1024) :: algo                ! Minimization algorithm

   ! Common data
   integer :: nl0                             ! Number of levels
   integer :: nc3                             ! Number of classes

   ! Specific data (fit)
   integer :: nl0r                            ! Reduced number of levels
   logical :: lhomh                           ! Vertically homogenous horizontal support radius key
   logical :: lhomv                           ! Vertically homogenous vertical support radius key
   integer,allocatable :: l0rl0_to_l0(:,:)    ! Reduced level to level
   real(kind_real),allocatable :: disth(:)    ! Horizontal distance
   real(kind_real),allocatable :: distv(:,:)  ! Vertical distance

   ! Specific data (LCT)
   integer :: nscales                         ! Number of LCT scales
   real(kind_real),allocatable :: dx(:,:)     ! Zonal separation
   real(kind_real),allocatable :: dy(:,:)     ! Meridian separation
   real(kind_real),allocatable :: dz(:,:)     ! Vertical separation
   logical,allocatable :: dmask(:,:)          ! Mask
contains
   procedure :: compute => minim_compute
   procedure :: cost => minim_cost
   procedure :: cost_fit_diag => minim_cost_fit_diag
   procedure :: cost_fit_diag_dble => minim_cost_fit_diag_dble
   procedure :: cost_fit_lct => minim_cost_fit_lct
   procedure :: hooke => minim_hooke
   procedure :: best_nearby => minim_best_nearby
   procedure :: vt_dir => minim_vt_dir
   procedure :: vt_inv => minim_vt_inv
end type minim_type

private
public :: minim_type

contains

!----------------------------------------------------------------------
! subroutine: minim_compute
! Purpose: minimize ensuring bounds constraints
!----------------------------------------------------------------------
subroutine minim_compute(minim,mpl,lprt)

implicit none

! Passed variables
class(minim_type),intent(inout) :: minim ! Minimization data
type(mpl_type),intent(inout) :: mpl      ! MPI data
logical,intent(in) :: lprt               ! Print key

! Local variables
real(kind_real) :: guess(minim%nx)
character(len=1024),parameter :: subr = 'minim_compute'

! Check
if (minim%nx<=0) call mpl%abort(subr,'nx should be positive to minimize')
if (minim%ny<=0) call mpl%abort(subr,'nx should be positive to minimize')

! Initialization
guess = minim%guess
call minim%vt_inv(mpl,guess)

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
   call minim%vt_dir(minim%x)
   if (lprt) then
      write(mpl%info,'(a13,a,f6.1,a,e9.2,a,e9.2,a)') '','Minimizer '//trim(minim%algo)//', cost function decrease:', &
    & abs(minim%f_min-minim%f_guess)/minim%f_guess*100.0,'% (',minim%f_guess,' to ',minim%f_min,')'
      call mpl%flush
   end if
else
   minim%x = minim%guess
   if (lprt) call mpl%warning(subr,'Minimizer '//trim(minim%algo)//' failed')
end if

end subroutine minim_compute

!----------------------------------------------------------------------
! Subroutine: minim_cost
! Purpose: compute cost function
!----------------------------------------------------------------------
subroutine minim_cost(minim,mpl,x,f)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim     ! Minimization data
type(mpl_type),intent(inout) :: mpl       ! MPI data
real(kind_real),intent(in) :: x(minim%nx) ! Control vector
real(kind_real),intent(out) :: f          ! Cost function value

select case (minim%cost_function)
case ('fit_diag')
   call minim%cost_fit_diag(mpl,x,f)
case ('fit_diag_dble')
   call minim%cost_fit_diag_dble(mpl,x,f)
case ('fit_lct')
   call minim%cost_fit_lct(mpl,x,f)
end select

end subroutine minim_cost

!----------------------------------------------------------------------
! Subroutine: minim_cost_fit_diag
! Purpose: diagnosic fit function cost
!----------------------------------------------------------------------
subroutine minim_cost_fit_diag(minim,mpl,x,f)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim     ! Minimization data
type(mpl_type),intent(inout) :: mpl       ! MPI data
real(kind_real),intent(in) :: x(minim%nx) ! Control vector
real(kind_real),intent(out) :: f          ! Cost function value

! Local variables
integer :: offset,il0
real(kind_real) :: fo,fh,fv,norm,fit_rh_avg,fit_rv_avg
real(kind_real) :: fit_rh(minim%nl0),fit_rv(minim%nl0)
real(kind_real) :: fit(minim%nc3,minim%nl0r,minim%nl0)
real(kind_real) :: xtmp(minim%nx),fit_pack(minim%ny)

! Check control vector size
offset = 0
if (minim%lhomh) then
   offset = offset+1
else
   offset = offset+minim%nl0
end if
if (minim%lhomv) then
   offset = offset+1
else
   offset = offset+minim%nl0
end if

! Renormalize
xtmp = x
call minim%vt_dir(xtmp)

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
call fit_diag(mpl,minim%nc3,minim%nl0r,minim%nl0,minim%l0rl0_to_l0,minim%disth,minim%distv,fit_rh,fit_rv,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Observations penalty
fo = sum((fit_pack-minim%obs)**2,mask=mpl%msv%isnotr(minim%obs).and.mpl%msv%isnotr(fit_pack))
norm = sum(minim%obs**2,mask=mpl%msv%isnotr(minim%obs).and.mpl%msv%isnotr(fit_pack))
if (norm>0.0) fo = fo/norm

! Smoothing penalty
fh = 0.0
fv = 0.0
do il0=2,minim%nl0-1
   fit_rh_avg = 0.5*(fit_rh(il0-1)+fit_rh(il0+1))
   norm = fit_rh_avg**2
   if (norm>0.0) fh = fh+(fit_rh(il0)-fit_rh_avg)**2/norm
   fit_rv_avg = 0.5*(fit_rv(il0-1)+fit_rv(il0+1))
   norm = fit_rv_avg**2
   if (norm>0.0) fv = fv+(fit_rv(il0)-fit_rv_avg)**2/norm
end do

! Full penalty function
f = fo+fh+fv

end subroutine minim_cost_fit_diag

!----------------------------------------------------------------------
! Subroutine: minim_cost_fit_diag_dble
! Purpose: diagnosic fit function cost, double fit
!----------------------------------------------------------------------
subroutine minim_cost_fit_diag_dble(minim,mpl,x,f)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim     ! Minimization data
type(mpl_type),intent(inout) :: mpl       ! MPI data
real(kind_real),intent(in) :: x(minim%nx) ! Control vector
real(kind_real),intent(out) :: f          ! Cost function value

! Local variables
integer :: offset,il0
real(kind_real) :: fo,fh,fv,fr,fc,norm,fit_rh_avg,fit_rv_avg,fit_rv_rfac_avg,fit_rv_coef_avg
real(kind_real) :: fit_rh(minim%nl0),fit_rv(minim%nl0),fit_rv_rfac(minim%nl0),fit_rv_coef(minim%nl0)
real(kind_real) :: fit(minim%nc3,minim%nl0r,minim%nl0)
real(kind_real) :: xtmp(minim%nx),fit_pack(minim%ny)

! Renormalize
xtmp = x
call minim%vt_dir(xtmp)

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
   fit_rv_rfac = xtmp(offset+1)
   offset = offset+1
   fit_rv_coef = xtmp(offset+1)
   offset = offset+1
else
   fit_rv = xtmp(offset+1:offset+minim%nl0)
   offset = offset+minim%nl0
   fit_rv_rfac = xtmp(offset+1:offset+minim%nl0)
   offset = offset+minim%nl0
   fit_rv_coef = xtmp(offset+1:offset+minim%nl0)
   offset = offset+minim%nl0
end if

! Compute function
call fit_diag_dble(mpl,minim%nc3,minim%nl0r,minim%nl0,minim%l0rl0_to_l0,minim%disth,minim%distv,fit_rh,fit_rv, &
 & fit_rv_rfac,fit_rv_coef,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Observations penalty
fo = sum((fit_pack-minim%obs)**2,mask=mpl%msv%isnotr(minim%obs).and.mpl%msv%isnotr(fit_pack))
norm = sum(minim%obs**2,mask=mpl%msv%isnotr(minim%obs).and.mpl%msv%isnotr(fit_pack))
if (norm>0.0) fo = fo/norm

! Smoothing penalty
fh = 0.0
fv = 0.0
fr = 0.0
fc = 0.0
do il0=2,minim%nl0-1
   fit_rh_avg = 0.5*(fit_rh(il0-1)+fit_rh(il0+1))
   norm = fit_rh_avg**2
   if (norm>0.0) fh = fh+(fit_rh(il0)-fit_rh_avg)**2/norm
   fit_rv_avg = 0.5*(fit_rv(il0-1)+fit_rv(il0+1))
   norm = fit_rv_avg**2
   if (norm>0.0) fv = fv+(fit_rv(il0)-fit_rv_avg)**2/norm
   fit_rv_rfac_avg = 0.5*(fit_rv_rfac(il0-1)+fit_rv_rfac(il0+1))
   norm = fit_rv_rfac_avg**2
   if (norm>0.0) fr = fr+(fit_rv_rfac(il0)-fit_rv_rfac_avg)**2/norm
   fit_rv_coef_avg = 0.5*(fit_rv_coef(il0-1)+fit_rv_coef(il0+1))
   norm = fit_rv_coef_avg**2
   if (norm>0.0) fc = fc+(fit_rv_coef(il0)-fit_rv_coef_avg)**2/norm
end do

! Full penalty function
f = fo+fh+fv+fr+fc

end subroutine minim_cost_fit_diag_dble

!----------------------------------------------------------------------
! Function: minim_cost_fit_lct
! Purpose: LCT fit function cost
!----------------------------------------------------------------------
subroutine minim_cost_fit_lct(minim,mpl,x,f)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim     ! Minimization data
type(mpl_type),intent(inout) :: mpl       ! MPI data
real(kind_real),intent(in) :: x(minim%nx) ! Control vector
real(kind_real),intent(out) :: f          ! Cost function value

! Local variables
integer :: iscales,icomp
real(kind_real) :: norm
real(kind_real) :: fit(minim%nc3,minim%nl0),D(4,minim%nscales)
real(kind_real) :: xtmp(minim%nx),fit_pack(minim%ny),coef(minim%nscales)

! Renormalize
xtmp = x
call minim%vt_dir(xtmp)

! Compute function
do iscales=1,minim%nscales
   do icomp=1,4
      D(icomp,iscales) = xtmp((iscales-1)*4+icomp)
   end do
end do
if (minim%nscales>1) then
   coef(1:minim%nscales-1) = xtmp(minim%nscales*4+1:minim%nscales*4+minim%nscales-1)
   coef(minim%nscales) = 1.0-sum(coef(1:minim%nscales-1))
else
   coef(1) = 1.0
end if
call fit_lct(mpl,minim%nc3,minim%nl0,minim%dx,minim%dy,minim%dz,minim%dmask,minim%nscales, &
 & xtmp(1:minim%nscales*4),coef,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Observations penalty
f = sum((fit_pack-minim%obs)**2,mask=mpl%msv%isnotr(minim%obs).and.mpl%msv%isnotr(fit_pack))
norm = sum(minim%obs**2,mask=mpl%msv%isnotr(minim%obs).and.mpl%msv%isnotr(fit_pack))
if (norm>0.0) f = f/norm

end subroutine minim_cost_fit_lct

!----------------------------------------------------------------------
! Subroutine: minim_hooke
! Purpose: seeks a minimizer of a scalar function of several variables
! Author: ALGOL original by Arthur Kaupe, C version by Mark Johnson, FORTRAN90 version by John Burkardt
!----------------------------------------------------------------------
subroutine minim_hooke(minim,mpl,guess)

implicit none

! Passed variables
class(minim_type),intent(inout) :: minim      ! Minimization data
type(mpl_type),intent(inout) :: mpl           ! MPI data
real(kind_real),intent(in) :: guess(minim%nx) ! Guess

! Local variables
integer :: funevals,i,iters,keep
real(kind_real) :: fbefore,newf,steplength,tmp
real(kind_real) :: delta(minim%nx),newx(minim%nx)

! Initialization
newx = guess
minim%x = guess
do i=1,minim%nx
   if (sup(abs(guess(i)),rth)) then
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
do while ((iters<itermax).and.inf(tol,steplength))
   ! Update iteration
   iters = iters + 1

   ! Find best new point, one coordinate at a time
   newx = minim%x
   call minim%best_nearby(mpl,delta,newx,fbefore,funevals,newf)

   ! If we made some improvements, pursue that direction
   keep = 1

   do while (inf(newf,fbefore).and.(keep==1))
      do i=1,minim%nx
         ! Arrange the sign of delta
         if (sup(newx(i),minim%x(i))) then
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
      if (inf(fbefore,newf)) exit

      ! Make sure that the differences between the new and the old points
      ! are due to actual displacements; beware of roundoff errors that
      ! might cause NEWF < FBEFORE.
      keep = 0

      do i=1,minim%nx
         if (inf(0.5*abs(delta(i)),abs(newx(i)-minim%x(i)))) then
            keep = 1
            exit
         end if
      end do
   end do

   if (infeq(tol,steplength).and.infeq(fbefore,newf)) then
      steplength = steplength*rho
      delta = delta*rho
   end if
end do

end subroutine minim_hooke

!----------------------------------------------------------------------
! Subroutine: minim_best_nearby
! Purpose: looks for a better nearby point, one coordinate at a time
! Author: ALGOL original by Arthur Kaupe, C version by Mark Johnson, FORTRAN90 version by John Burkardt
!----------------------------------------------------------------------
subroutine minim_best_nearby(minim,mpl,delta,point,prevbest,funevals,minf)

implicit none

! Passed variables
class(minim_type),intent(inout) :: minim         ! Minimization data
type(mpl_type),intent(inout) :: mpl              ! MPI data
real(kind_real),intent(inout) :: delta(minim%nx) ! Step
real(kind_real),intent(inout) :: point(minim%nx) ! Point
real(kind_real),intent(in) :: prevbest           ! Best existing cost
integer,intent(inout) :: funevals                ! Number of evaluations
real(kind_real),intent(out) :: minf              ! Minimum cost

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
   if (inf(ftmp,minf)) then
      minf = ftmp
   else
      delta(i) = -delta(i)
      z(i) = point(i)+delta(i)
      call minim%cost(mpl,z,ftmp)
      funevals = funevals+1
      if (inf(ftmp,minf)) then
         minf = ftmp
      else
         z(i) = point(i)
      end if
   end if
end do

! Update
point = z

end subroutine minim_best_nearby

!----------------------------------------------------------------------
! Subroutine: vt_dir
! Purpose: direct variable transform
!----------------------------------------------------------------------
subroutine minim_vt_dir(minim,x)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim        ! Minimization data
real(kind_real),intent(inout) :: x(minim%nx) ! Vector

! Linear expansion of the hyperbolic tangent of the variable
x = minim%binf+0.5*(1.0+tanh(x))*(minim%bsup-minim%binf)

end subroutine minim_vt_dir

!----------------------------------------------------------------------
! Subroutine: vt_inv
! Purpose: inverse variable transform
!----------------------------------------------------------------------
subroutine minim_vt_inv(minim,mpl,x)

implicit none

! Passed variables
class(minim_type),intent(in) :: minim        ! Minimization data
type(mpl_type),intent(inout) :: mpl          ! MPI data
real(kind_real),intent(inout) :: x(minim%nx) ! Vector

! Local variables
integer :: ix
character(len=1024),parameter :: subr = 'minim_vt_inv'

! Inverse hyperbolic tangent of the linearly bounded variable
if (any((x<minim%binf).or.(x>minim%bsup))) call mpl%abort(subr,'variable out of bounds in vt_inv')
do ix=1,minim%nx
   if (sup(minim%bsup(ix),minim%binf(ix))) then
      x(ix) = atanh(2.0*(x(ix)-minim%binf(ix))/(minim%bsup(ix)-minim%binf(ix))-1.0)
   else
      x(ix) = 0.0
   end if
end do

end subroutine minim_vt_inv

end module type_minim
