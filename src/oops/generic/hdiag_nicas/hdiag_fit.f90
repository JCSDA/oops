!----------------------------------------------------------------------
! Module: hdiag_fit.f90
!> Purpose: fit routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_fit

use omp_lib
use tools_const, only: gc99,reqkm
use tools_display, only: msgerror,msgwarning,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_fit, only: fast_fit,ver_smooth
use tools_minim, only: minim
use tools_missing, only: msi,msr,isnotmsr
use type_curve, only: curvetype,curve_pack,curve_unpack
use type_geom, only: geomtype
use type_hdata, only: hdatatype
use type_mdata, only: mdatatype
use type_mpl, only: mpl
use type_nam, only: namtype

implicit none

integer,parameter :: nsc = 50 !< Scaling optimization parameter
logical :: prt = .false.      !< Print fit optimization results

interface compute_fit
  module procedure compute_fit
  module procedure compute_fit_local
end interface

private
public :: compute_fit

contains

!----------------------------------------------------------------------
! Subroutine: compute_fit
!> Purpose: compute a semi-positive definite fit of a raw function
!----------------------------------------------------------------------
subroutine compute_fit(hdata,ib,curve)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata    !< HDIAG data
integer,intent(in) :: ib               !< Block index
type(curvetype),intent(inout) :: curve !< Curve

! Local variables
integer :: il0,jl0r,jl0,offset,isc
real(kind_real) :: distvr(hdata%nam%nl0r,hdata%geom%nl0),rawv(hdata%nam%nl0r)
real(kind_real) :: alpha,alpha_opt,mse,mse_opt
real(kind_real) :: fit_rh(hdata%geom%nl0),fit_rv(hdata%geom%nl0),fit(hdata%nam%nc3,hdata%nam%nl0r,hdata%geom%nl0)
type(mdatatype) :: mdata

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Check
if (trim(nam%fit_type)=='none') call msgerror('cannot compute fit if fit_type = none')

! Initialization
call msr(curve%fit_rh)
call msr(curve%fit_rv)
call msr(curve%fit)

! Reduced vertical distance
call msr(distvr)
do il0=1,geom%nl0
   do jl0r=1,nam%nl0r
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
      distvr(jl0r,il0) = geom%distv(jl0,il0)
   end do
end do

! Fast fit
do il0=1,geom%nl0
   ! Get zero separation level
   jl0r = bpar%il0rz(il0,ib)

   ! Horizontal fast fit
   call fast_fit(nam%nc3,1,geom%disth,curve%raw(:,jl0r,il0),curve%fit_rh(il0))

   ! Vertical fast fit
   rawv = curve%raw(1,:,il0)
   call fast_fit(nam%nl0r,jl0r,distvr(:,il0),rawv,curve%fit_rv(il0))
end do
if (nam%lhomh) curve%fit_rh = sum(curve%fit_rh,mask=isnotmsr(curve%fit_rh))/float(count(isnotmsr(curve%fit_rh)))
if (nam%lhomv) curve%fit_rv = sum(curve%fit_rv,mask=isnotmsr(curve%fit_rv))/float(count(isnotmsr(curve%fit_rv)))

! Scaling optimization (brute-force)
mse_opt = huge(1.0)
alpha_opt = 1.0
do isc=1,nsc
   ! Scaling factor
   alpha = 0.5+float(isc-1)/float(nsc-1)*(2.0-0.5)

   ! Scaled radii
   fit_rh = alpha*curve%fit_rh
   fit_rv = alpha*curve%fit_rv

   ! Fit
   call define_fit(nam%nc3,nam%nl0r,geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth,distvr,fit_rh,fit_rv,fit)

   ! MSE
   mse = sum((fit-curve%raw)**2,mask=isnotmsr(curve%raw))
   if (mse<mse_opt) then
      mse_opt = mse
      alpha_opt = alpha
   end if
end do
do il0=1,geom%nl0
   if (isnotmsr(curve%fit_rh(il0))) curve%fit_rh(il0) = alpha_opt*curve%fit_rh(il0)
   if (isnotmsr(curve%fit_rv(il0))) curve%fit_rv(il0) = alpha_opt*curve%fit_rv(il0)
end do
if (prt) write(mpl%unit,'(a7,a,f6.1,a)') '','Scaling optimization, cost function decrease:',abs(mse_opt-mse)/mse*100.0,'%'

select case (trim(nam%fit_type))
case ('nelder_mead','compass_search','praxis')
   ! Allocation
   mdata%nx = 0
   if (nam%lhomh) then
      mdata%nx = mdata%nx+1
   else
      mdata%nx = mdata%nx+geom%nl0
   end if
   if (nam%lhomv) then
      mdata%nx = mdata%nx+1
   else
      mdata%nx = mdata%nx+geom%nl0
   end if
   mdata%ny = nam%nc3*nam%nl0r*geom%nl0
   allocate(mdata%x(mdata%nx))
   allocate(mdata%guess(mdata%nx))
   allocate(mdata%norm(mdata%nx))
   allocate(mdata%binf(mdata%nx))
   allocate(mdata%bsup(mdata%nx))
   allocate(mdata%obs(mdata%ny))
   allocate(mdata%l0rl0_to_l0(nam%nl0r,geom%nl0))
   allocate(mdata%disth(nam%nc3))
   allocate(mdata%distvr(nam%nl0r,geom%nl0))

   ! Fill mdata
   offset = 0
   if (nam%lhomh) then
      mdata%guess(offset+1) = curve%fit_rh(1)
      offset = offset+1
   else
      mdata%guess(offset+1:offset+geom%nl0) = curve%fit_rh
      offset = offset+geom%nl0
   end if
   if (nam%lhomv) then
      mdata%guess(offset+1) = curve%fit_rv(1)
      offset = offset+1
   else
      mdata%guess(offset+1:offset+geom%nl0) = curve%fit_rv
      offset = offset+geom%nl0
   end if
   mdata%norm = mdata%guess
   mdata%binf = 0.5*mdata%guess
   mdata%bsup = 1.5*mdata%guess
   mdata%obs = pack(curve%raw,mask=.true.)
   mdata%fit_type = nam%fit_type
   mdata%nc3 = nam%nc3
   mdata%nl0r = nam%nl0r
   mdata%nl0 = geom%nl0
   mdata%lhomh = nam%lhomh
   mdata%lhomv = nam%lhomv
   mdata%l0rl0_to_l0 = bpar%l0rl0b_to_l0(:,:,ib)
   mdata%disth = geom%disth
   mdata%distvr = distvr

   ! Compute fit
   call minim(mdata,func,prt)

   ! Apply bounds
   mdata%x = max(mdata%binf,min(mdata%x,mdata%bsup))

   ! Copy parameters
   offset = 0
   if (nam%lhomh) then
      curve%fit_rh = mdata%x(offset+1)
      offset = offset+1
   else
      do il0=1,geom%nl0
         if (isnotmsr(curve%fit_rh(il0))) curve%fit_rh(il0) = mdata%x(offset+il0)
      end do
      offset = offset+geom%nl0
   end if
   if (nam%lhomv) then
      curve%fit_rv = mdata%x(offset+1)
      offset = offset+1
   else
      do il0=1,geom%nl0
         if (isnotmsr(curve%fit_rv(il0))) curve%fit_rv(il0) = mdata%x(offset+il0)
      end do
      offset = offset+geom%nl0
   end if

   ! Dummy call to avoid warnings
   call dummy(mdata)
end select

! Smooth vertically
call ver_smooth(geom%nl0,geom%vunit,nam%rvflt,curve%fit_rh)
call ver_smooth(geom%nl0,geom%vunit,nam%rvflt,curve%fit_rv)

! Rebuild fit
call define_fit(nam%nc3,nam%nl0r,geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth,distvr,curve%fit_rh,curve%fit_rv,curve%fit)

! End associate
end associate

end subroutine compute_fit

!----------------------------------------------------------------------
! Subroutine: compute_fit_local
!> Purpose: compute a semi-positive definite fit of a raw function, multiple curves
!----------------------------------------------------------------------
subroutine compute_fit_local(hdata,ib,curve)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                !< HDIAG data
integer,intent(in) :: ib                           !< Block index
type(curvetype),intent(inout) :: curve(hdata%nc2a) !< Curve

! Local variables
integer :: ic2a,progint
logical :: done(hdata%nc2a)

! Initialization
call prog_init(progint,done)

! Loop over points
do ic2a=1,hdata%nc2a
   ! Compute fit
   call compute_fit(hdata,ib,curve(ic2a))

   ! Print progression
   done(ic2a) = .true.
   call prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'

end subroutine compute_fit_local

!----------------------------------------------------------------------
! Subroutine: define_fit
!> Purpose: define the fit
!----------------------------------------------------------------------
subroutine define_fit(nc3,nl0r,nl0,l0rl0_to_l0,disth,distvr,rh,rv,fit)

implicit none

! Passed variables
integer,intent(in) :: nc3                        !< Number of classes
integer,intent(in) :: nl0r                       !< Reduced number of levels
integer,intent(in) :: nl0                        !< Number of levels
integer,intent(in) :: l0rl0_to_l0(nl0r,nl0)      !< Reduced level to level
real(kind_real),intent(in) :: disth(nc3)         !< Horizontal distance
real(kind_real),intent(in) :: distvr(nl0r,nl0)   !< Vertical distance
real(kind_real),intent(in) :: rh(nl0)            !< Horizontal support radius
real(kind_real),intent(in) :: rv(nl0)            !< Vertical support radius
real(kind_real),intent(out) :: fit(nc3,nl0r,nl0) !< Fit

! Local variables
integer :: jl0r,jl0,il0,kl0r,kl0,jc3,kc3,ip,jp,np,np_new
integer,allocatable :: plist(:,:),plist_new(:,:)
real(kind_real) :: rhsq,rvsq,distnorm,disttest
real(kind_real),allocatable :: dist(:,:)
logical :: add_to_front

! Initialization
fit = 0.0

do il0=1,nl0
   ! Allocation
   allocate(plist(nc3*nl0r,2))
   allocate(plist_new(nc3*nl0r,2))
   allocate(dist(nc3,nl0r))

   ! Initialize the front
   np = 1
   call msi(plist)
   plist(1,1) = 1
   do jl0r=1,nl0r
      if (l0rl0_to_l0(jl0r,il0)==il0) plist(1,2) = jl0r
   end do
   dist = 1.0
   dist(plist(1,1),plist(1,2)) = 0.0

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         jc3 = plist(ip,1)
         jl0r = plist(ip,2)
         jl0 = l0rl0_to_l0(jl0r,il0)

         ! Loop over neighbors
         do kc3=max(jc3-1,1),min(jc3+1,nc3)
            do kl0r=max(jl0r-1,1),min(jl0r+1,nl0r)
               kl0 = l0rl0_to_l0(kl0r,il0)
               if (isnotmsr(rh(jl0)).and.isnotmsr(rh(kl0))) then
                  rhsq = 0.5*(rh(jl0)**2+rh(kl0)**2)
               else
                  rhsq = 0.0
               end if
               if (isnotmsr(rv(jl0)).and.isnotmsr(rv(kl0))) then
                  rvsq = 0.5*(rv(jl0)**2+rv(kl0)**2)
               else
                  rvsq = 0.0
               end if
               distnorm = 0.0
               if (rhsq>0.0) then
                  distnorm = distnorm+(disth(kc3)-disth(jc3))**2/rhsq
               elseif (kc3/=jc3) then
                  distnorm = distnorm+0.5*huge(1.0)
               end if
               if (rvsq>0.0) then
                  distnorm = distnorm+distvr(kl0r,jl0)**2/rvsq
               elseif (kl0r/=jl0r) then
                  distnorm = distnorm+0.5*huge(1.0)
               end if
               disttest = dist(jc3,jl0r)+sqrt(distnorm)
               if (disttest<1.0) then
                  ! Point is inside the support
                  if (disttest<dist(kc3,kl0r)) then
                     ! Update distance
                     dist(kc3,kl0r) = disttest

                     ! Check if the point should be added to the front (avoid duplicates)
                     add_to_front = .true.
                     do jp=1,np_new
                        if ((plist_new(jp,1)==kc3).and.(plist_new(jp,2)==kl0r)) then
                           add_to_front = .false.
                           exit
                        end if
                     end do

                     if (add_to_front) then
                        ! Add point to the front
                        np_new = np_new+1
                        plist_new(np_new,1) = kc3
                        plist_new(np_new,2) = kl0r
                     end if
                  end if
               end if
            end do
         end do
      end do

      ! Copy new front
      np = np_new
      plist(1:np,:) = plist_new(1:np,:)
   end do

   do jl0r=1,nl0r
      do jc3=1,nc3
         ! Gaspari-Cohn (1999) function
         distnorm = dist(jc3,jl0r)
         if (distnorm<1.0) fit(jc3,jl0r,il0) = gc99(distnorm)
      end do
   end do

   ! Release memory
   deallocate(plist)
   deallocate(plist_new)
   deallocate(dist)
end do

end subroutine define_fit

!----------------------------------------------------------------------
! Function: func
!> Purpose: fit function cost
!----------------------------------------------------------------------
subroutine func(mdata,x,f)

implicit none

! Passed variables
type(mdatatype),intent(in) :: mdata       !< Minimization data
real(kind_real),intent(in) :: x(mdata%nx) !< Control vector
real(kind_real),intent(out) :: f          !< Cost function value

! Local variables
integer :: offset,ix
real(kind_real) :: fit_rh(mdata%nl0),fit_rv(mdata%nl0)
real(kind_real) :: fit(mdata%nc3,mdata%nl0r,mdata%nl0)
real(kind_real) :: xtmp(mdata%nx),fit_pack(mdata%ny),xx

! Renormalize
xtmp = x*mdata%norm

! Get data
offset = 0
if (mdata%lhomh) then
   fit_rh = xtmp(offset+1)
   offset = offset+1
else
   fit_rh = xtmp(offset+1:offset+mdata%nl0)
   offset = offset+mdata%nl0
end if
if (mdata%lhomv) then
   fit_rv = xtmp(offset+1)
   offset = offset+1
else
   fit_rv = xtmp(offset+1:offset+mdata%nl0)
   offset = offset+mdata%nl0
end if

! Compute function
call define_fit(mdata%nc3,mdata%nl0r,mdata%nl0,mdata%l0rl0_to_l0,mdata%disth,mdata%distvr,fit_rh,fit_rv,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Cost
f = sum((mdata%obs-fit_pack)**2,mask=isnotmsr(mdata%obs))

! Bound penalty
do ix=1,mdata%nx
   if (mdata%bsup(ix)>mdata%binf(ix)) then
      xx = (xtmp(ix)-mdata%binf(ix))/(mdata%bsup(ix)-mdata%binf(ix))
      if (xx<0.0) then
         f = f+mdata%f_guess*xx**2
      elseif (xx>1.0) then
         f = f+mdata%f_guess*(xx-1.0)**2
      end if
   end if
end do

end subroutine func

!----------------------------------------------------------------------
! Subroutine: dummy
!> Purpose: dummy subroutine to avoid warnings
!----------------------------------------------------------------------
subroutine dummy(mdata)

implicit none

! Passed variables
type(mdatatype),intent(in) :: mdata !< Minimization data

! Local variables
real(kind_real) :: x(mdata%nx)
real(kind_real) :: f

if (.false.) call func(mdata,x,f)

end subroutine

end module hdiag_fit
