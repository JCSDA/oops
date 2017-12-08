!----------------------------------------------------------------------
! Module: module_fit.f90
!> Purpose: fit routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_fit

use omp_lib
use tools_const, only: gc99
use tools_display, only: msgerror,msgwarning,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_fit, only: fast_fit,ver_smooth
use tools_minim, only: minim
use tools_missing, only: msi,msr,isnotmsr
use type_curve, only: curvetype,curve_pack,curve_unpack
use type_geom, only: geomtype
use type_hdata, only: hdatatype
use type_min, only: mintype
use type_mpl, only: mpl,mpl_recv,mpl_send,mpl_bcast
use type_nam, only: namtype

implicit none

integer,parameter :: nsc = 50 !< Scaling optimization parameter

interface compute_fit
  module procedure compute_fit
  module procedure compute_fit_multi
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
integer :: il0r,il0,jl0,offset,isc
real(kind_real) :: distvr(hdata%nam%nl0r,hdata%geom%nl0),rawv(hdata%nam%nl0r)
real(kind_real) :: alpha,alpha_opt,mse,mse_opt
real(kind_real) :: fit_rh(hdata%geom%nl0),fit_rv(hdata%geom%nl0),fit(hdata%nam%nc,hdata%nam%nl0r,hdata%geom%nl0)
type(mintype) :: mindata

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
do jl0=1,geom%nl0
   do il0r=1,nam%nl0r
      il0 = bpar%il0rjl0ib_to_il0(il0r,jl0,ib)
      distvr(il0r,jl0) = geom%distv(il0,jl0)
   end do
end do

! Fast fit
do jl0=1,geom%nl0
   ! Get zero separation level
   il0r = bpar%il0rz(jl0,ib)

   ! Horizontal fast fit
   call fast_fit(nam%nc,1,geom%disth,curve%raw(:,il0r,jl0),curve%fit_rh(jl0))

   ! Vertical fast fit
   rawv = curve%raw(1,:,jl0)
   call fast_fit(nam%nl0r,il0r,distvr(:,jl0),rawv,curve%fit_rv(jl0))
end do
if (nam%lhomh) curve%fit_rh = sum(curve%fit_rh)/float(geom%nl0)
if (nam%lhomv) curve%fit_rv = sum(curve%fit_rv)/float(geom%nl0)

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
   call define_fit(nam%nc,nam%nl0r,geom%nl0,bpar%il0rjl0ib_to_il0(:,:,ib),geom%disth,distvr,fit_rh,fit_rv,fit)

   ! MSE
   mse = sum((fit-curve%raw)**2,mask=isnotmsr(curve%raw))
   if (mse<mse_opt) then
      mse_opt = mse
      alpha_opt = alpha
   end if
end do
curve%fit_rh = alpha_opt*curve%fit_rh
curve%fit_rv = alpha_opt*curve%fit_rv
write(mpl%unit,'(a7,a,f6.1,a)') '','Scaling optimization, cost function decrease:',abs(mse_opt-mse)/mse*100.0,'%'

select case (trim(nam%fit_type))
case ('nelder_mead','compass_search','praxis')
   ! Allocation
   mindata%nx = 0
   if (nam%lhomh) then
      mindata%nx = mindata%nx+1
   else
      mindata%nx = mindata%nx+geom%nl0
   end if
   if (nam%lhomv) then
      mindata%nx = mindata%nx+1
   else
      mindata%nx = mindata%nx+geom%nl0
   end if
   mindata%ny = nam%nc*nam%nl0r*geom%nl0
   allocate(mindata%x(mindata%nx))
   allocate(mindata%guess(mindata%nx))
   allocate(mindata%norm(mindata%nx))
   allocate(mindata%binf(mindata%nx))
   allocate(mindata%bsup(mindata%nx))
   allocate(mindata%obs(mindata%ny))
   allocate(mindata%il0rjl0_to_il0(nam%nl0r,nam%nl0r))
   allocate(mindata%disth(nam%nc))
   allocate(mindata%distvr(nam%nl0r,geom%nl0))

   ! Fill mindata
   offset = 0
   if (nam%lhomh) then
      mindata%guess(offset+1) = curve%fit_rh(1)
      offset = offset+1
   else
      mindata%guess(offset+1:offset+geom%nl0) = curve%fit_rh
      offset = offset+geom%nl0
   end if
   if (nam%lhomv) then
      mindata%guess(offset+1) = curve%fit_rv(1)
      offset = offset+1
   else
      mindata%guess(offset+1:offset+geom%nl0) = curve%fit_rv
      offset = offset+geom%nl0
   end if
   mindata%norm = mindata%guess
   mindata%binf = 0.75*mindata%guess
   mindata%bsup = 1.25*mindata%guess
   mindata%obs = pack(curve%raw,mask=.true.)
   mindata%fit_type = nam%fit_type
   mindata%nc = nam%nc
   mindata%nl0r = nam%nl0r
   mindata%nl0 = geom%nl0
   mindata%lhomh = nam%lhomh
   mindata%lhomv = nam%lhomv
   mindata%il0rjl0_to_il0 = bpar%il0rjl0ib_to_il0(:,:,ib)
   mindata%disth = geom%disth
   mindata%distvr = distvr

   ! Compute fit
   call minim(mindata,func,.true.)

   ! Copy parameters
   offset = 0
   if (nam%lhomh) then
      curve%fit_rh = mindata%x(offset+1)
      offset = offset+1
   else
      curve%fit_rh = mindata%x(offset+1:offset+geom%nl0)
      offset = offset+geom%nl0
   end if
   if (nam%lhomv) then
      curve%fit_rv = mindata%x(offset+1)
      offset = offset+1
   else
      curve%fit_rv = mindata%x(offset+1:offset+geom%nl0)
      offset = offset+geom%nl0
   end if

   ! Dummy call to avoid warnings
   call dummy(mindata)
end select

! Smooth vertically
call ver_smooth(geom%nl0,geom%vunit,nam%rvflt,curve%fit_rh)
call ver_smooth(geom%nl0,geom%vunit,nam%rvflt,curve%fit_rv)

! Rebuild fit
call define_fit(nam%nc,nam%nl0r,geom%nl0,bpar%il0rjl0ib_to_il0(:,:,ib),geom%disth,distvr,curve%fit_rh,curve%fit_rv,curve%fit)

! End associate
end associate

end subroutine compute_fit

!----------------------------------------------------------------------
! Subroutine: compute_fit_multi
!> Purpose: compute a semi-positive definite fit of a raw function, multiple curves
!----------------------------------------------------------------------
subroutine compute_fit_multi(hdata,ib,curve)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata               !< HDIAG data
integer,intent(in) :: ib               !< Block index
type(curvetype),intent(inout) :: curve(hdata%nc2) !< Curve

! Local variables
integer :: ic2,npack,progint
integer :: iproc,ic2_s(mpl%nproc),ic2_e(mpl%nproc),nc2_loc(mpl%nproc),ic2_loc
real(kind_real),allocatable :: rbuf(:),sbuf(:)
logical,allocatable :: done(:)

! MPI splitting
do iproc=1,mpl%nproc
   ic2_s(iproc) = (iproc-1)*(hdata%nc2/mpl%nproc+1)+1
   ic2_e(iproc) = min(iproc*(hdata%nc2/mpl%nproc+1),hdata%nc2)
   nc2_loc(iproc) = ic2_e(iproc)-ic2_s(iproc)+1
end do

! Allocation
npack = curve(1)%npack
allocate(rbuf(hdata%nc2*npack))
allocate(done(nc2_loc(mpl%myproc)))

! Loop over points
call prog_init(progint,done)
do ic2_loc=1,nc2_loc(mpl%myproc)
   ic2 = ic2_s(mpl%myproc)+ic2_loc-1

   ! Compute fit
   call compute_fit(hdata,ib,curve(ic2))

   ! Print progression
   done(ic2_loc) = .true.
   call prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Format data
         do ic2_loc=1,nc2_loc(iproc)
            ic2 = ic2_s(iproc)+ic2_loc-1
            call curve_pack(hdata,curve(ic2),rbuf((ic2-1)*npack+1:ic2*npack))
         end do
      else
         ! Receive data on ioproc
         call mpl_recv(nc2_loc(iproc)*npack, &
       & rbuf((ic2_s(iproc)-1)*npack+1:ic2_e(iproc)*npack),iproc,mpl%tag)
      end if
   end do
else
   ! Allocation
   allocate(sbuf(nc2_loc(mpl%myproc)*npack))

   ! Format data
   do ic2_loc=1,nc2_loc(mpl%myproc)
      ic2 = ic2_s(mpl%myproc)+ic2_loc-1
      call curve_pack(hdata,curve(ic2),sbuf((ic2_loc-1)*npack+1:ic2_loc*npack))
   end do

   ! Send data to ioproc
   call mpl_send(nc2_loc(mpl%myproc)*npack,sbuf,mpl%ioproc,mpl%tag)

   ! Release memory
   deallocate(sbuf)
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl_bcast(rbuf,mpl%ioproc)

! Format data
do ic2=1,hdata%nc2
   call curve_unpack(hdata,curve(ic2),rbuf((ic2-1)*npack+1:ic2*npack))
end do

! Release memory
deallocate(rbuf)

end subroutine compute_fit_multi

!----------------------------------------------------------------------
! Subroutine: define_fit
!> Purpose: define the fit
!----------------------------------------------------------------------
subroutine define_fit(nc,nl0r,nl0,il0rjl0_to_il0,disth,distvr,rh,rv,fit)

implicit none

! Passed variables
integer,intent(in) :: nc                        !< Number of classes
integer,intent(in) :: nl0r                      !< Reduced number of levels
integer,intent(in) :: nl0                       !< Number of levels
integer,intent(in) :: il0rjl0_to_il0(nl0r,nl0)  !< Reduced level to level
real(kind_real),intent(in) :: disth(nc)         !< Horizontal distance
real(kind_real),intent(in) :: distvr(nl0r,nl0)  !< Vertical distance
real(kind_real),intent(in) :: rh(nl0)           !< Horizontal support radius
real(kind_real),intent(in) :: rv(nl0)           !< Vertical support radius
real(kind_real),intent(out) :: fit(nc,nl0r,nl0) !< Fit

! Local variables
integer :: il0r,il0,jl0,kl0r,kl0,ic,kc,ip,jp,np,np_new
integer,allocatable :: plist(:,:),plist_new(:,:)
real(kind_real) :: rhsq,rvsq,distnorm,disttest
real(kind_real),allocatable :: dist(:,:)
logical :: add_to_front

! Initialization
fit = 0.0

do jl0=1,nl0
   ! Allocation
   allocate(plist(nc*nl0r,2))
   allocate(plist_new(nc*nl0r,2))
   allocate(dist(nc,nl0r))

   ! Initialize the front
   np = 1
   call msi(plist)
   plist(1,1) = 1
   do il0r=1,nl0r
      if (il0rjl0_to_il0(il0r,jl0)==jl0) plist(1,2) = il0r
   end do
   dist = 1.0
   dist(plist(1,1),plist(1,2)) = 0.0

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         ic = plist(ip,1)
         il0r = plist(ip,2)
         il0 = il0rjl0_to_il0(il0r,jl0)

         ! Loop over neighbors
         do kc=max(ic-1,1),min(ic+1,nc)
            do kl0r=max(il0r-1,1),min(il0r+1,nl0r)
               kl0 = il0rjl0_to_il0(kl0r,jl0)
               if (isnotmsr(rh(il0)).and.isnotmsr(rh(kl0))) then
                  rhsq = 0.5*(rh(il0)**2+rh(kl0)**2)
               else
                  rhsq = 0.0
               end if
               if (isnotmsr(rv(il0)).and.isnotmsr(rv(kl0))) then
                  rvsq = 0.5*(rv(il0)**2+rv(kl0)**2)
               else
                  rvsq = 0.0
               end if
               distnorm = 0.0
               if (rhsq>0.0) then
                  distnorm = distnorm+(disth(kc)-disth(ic))**2/rhsq
               else
                  distnorm = huge(1.0)
               end if
               if (rvsq>0.0) then
                  distnorm = distnorm+distvr(kl0r,il0)**2/rvsq
               elseif (kl0r/=il0r) then
                  distnorm = huge(1.0)
               end if
               distnorm = sqrt(distnorm)
               disttest = dist(ic,il0r)+distnorm
               if (disttest<1.0) then
                  ! Point is inside the support
                  if (disttest<dist(kc,kl0r)) then
                     ! Update distance
                     dist(kc,kl0r) = disttest

                     ! Check if the point should be added to the front (avoid duplicates)
                     add_to_front = .true.
                     do jp=1,np_new
                        if ((plist_new(jp,1)==kc).and.(plist_new(jp,1)==kl0r)) then
                           add_to_front = .false.
                           exit
                        end if
                     end do

                     if (add_to_front) then
                        ! Add point to the front
                        np_new = np_new+1
                        plist_new(np_new,1) = kc
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

   do il0r=1,nl0r
      do ic=1,nc
         ! Gaspari-Cohn (1999) function
         distnorm = dist(ic,il0r)
         if (distnorm<1.0) fit(ic,il0r,jl0) = gc99(distnorm)
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
subroutine func(mindata,x,f)

implicit none

! Passed variables
type(mintype),intent(in) :: mindata         !< Minimization data
real(kind_real),intent(in) :: x(mindata%nx) !< Control vector
real(kind_real),intent(out) :: f            !< Cost function value

! Local variables
integer :: offset,ix
real(kind_real) :: fit_rh(mindata%nl0),fit_rv(mindata%nl0)
real(kind_real) :: fit(mindata%nc,mindata%nl0r,mindata%nl0)
real(kind_real) :: xtmp(mindata%nx),fit_pack(mindata%ny),xx

! Renormalize
xtmp = x*mindata%norm

! Get data
offset = 0
if (mindata%lhomh) then
   fit_rh = xtmp(offset+1)
   offset = offset+1
else
   fit_rh = xtmp(offset+1:offset+mindata%nl0)
   offset = offset+mindata%nl0
end if
if (mindata%lhomv) then
   fit_rv = xtmp(offset+1)
   offset = offset+1
else
   fit_rv = xtmp(offset+1:offset+mindata%nl0)
   offset = offset+mindata%nl0
end if

! Compute function
call define_fit(mindata%nc,mindata%nl0r,mindata%nl0,mindata%il0rjl0_to_il0,mindata%disth,mindata%distvr,fit_rh,fit_rv,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Cost
f = sum((mindata%obs-fit_pack)**2,mask=isnotmsr(mindata%obs))

! Bound penalty
do ix=1,mindata%nx
   if (mindata%bsup(ix)>mindata%binf(ix)) then
      xx = (xtmp(ix)-mindata%binf(ix))/(mindata%bsup(ix)-mindata%binf(ix))
      if (xx<0.0) then
         f = f+mindata%f_guess*xx**2
      elseif (xx>1.0) then
         f = f+mindata%f_guess*(xx-1.0)**2
      end if
   end if
end do

end subroutine func

!----------------------------------------------------------------------
! Subroutine: dummy
!> Purpose: dummy subroutine to avoid warnings
!----------------------------------------------------------------------
subroutine dummy(mindata)

implicit none

! Passed variables
type(mintype),intent(in) :: mindata !< Minimization data

! Local variables
real(kind_real) :: x(mindata%nx)
real(kind_real) :: f

if (.false.) call func(mindata,x,f)

end subroutine

end module module_fit
