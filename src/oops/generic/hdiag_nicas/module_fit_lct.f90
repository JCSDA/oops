!----------------------------------------------------------------------
! Module: module_fit_lct.f90
!> Purpose: LCT fit routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_fit_lct

use omp_lib
use tools_const, only: lonmod
use tools_diffusion, only: matern
use tools_display, only: msgerror,msgwarning,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_minim, only: minim
use tools_missing, only: msr,isnotmsr
use type_hdata, only: hdatatype
use type_lct, only: lcttype,lct_alloc,lct_pack,lct_unpack
use type_min, only: mintype
use type_mpl, only: mpl,mpl_recv,mpl_send,mpl_bcast

implicit none

real(kind_real),parameter :: Hmin = 1.0e-12 !< Minimum tensor diagonal value
real(kind_real),parameter :: Hscale = 10.0  !< Typical factor between LCT scales
integer,parameter :: M = 0                  !< Number of implicit itteration for the Matern function (Gaussian function if M = 0)

interface compute_fit_lct
  module procedure compute_fit_lct
  module procedure compute_fit_lct_multi
end interface

private
public :: compute_fit_lct

contains

!----------------------------------------------------------------------
! Subroutine: compute_fit_lct
!> Purpose: compute a semi-positive definite fit of a raw function
!----------------------------------------------------------------------
subroutine compute_fit_lct(hdata,ib,dx,dy,dz,dmask,lct)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                               !< HDIAG data
integer,intent(in) :: ib                                          !< Block index
real(kind_real),intent(in) :: dx(hdata%nam%nc,hdata%bpar%nl0(ib)) !< Zonal separation
real(kind_real),intent(in) :: dy(hdata%nam%nc,hdata%bpar%nl0(ib)) !< Meridian separation
real(kind_real),intent(in) :: dz(hdata%bpar%nl0(ib))              !< Vertical separation
logical,intent(in) :: dmask(hdata%nam%nc,hdata%bpar%nl0(ib))      !< Mask
type(lcttype),intent(inout) :: lct                                !< LCT

! Local variables
integer :: il0,ic,iscales,offset
real(kind_real) :: distsq,Hh(hdata%nam%nc),Hv(hdata%bpar%nl0(ib)),Hhbar,Hvbar,det
logical :: spd
type(lcttype) :: lct_guess,lct_norm,lct_binf,lct_bsup
type(mintype) :: mindata

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Approximate homogeneous horizontal length-scale
call msr(Hh)
do il0=1,bpar%nl0(ib)
   if (.not.(abs(dz(il0))>0.0)) then
      do ic=1,nam%nc
         if (dmask(ic,il0)) then
            distsq = dx(ic,il0)**2+dy(ic,il0)**2
            if ((lct%raw(ic,il0)>0.0).and.(distsq>0.0)) Hh(ic) = -2.0*log(lct%raw(ic,il0))/distsq
         end if
      end do
   end if
end do
if (count(isnotmsr(Hh))>0) then
   Hhbar = sum(Hh,mask=isnotmsr(Hh))/float(count(isnotmsr(Hh)))
else
   return
end if
if (lct%nscales>1) Hhbar = Hhbar*Hscale

! Approximate homogeneous vertical length-scale
call msr(Hv)
ic = 1
do il0=1,bpar%nl0(ib)
   distsq = dz(il0)**2
   if ((lct%raw(ic,il0)>0.0).and.(distsq>0.0)) Hv(il0) = -2.0*log(lct%raw(ic,il0))/distsq
end do
if (bpar%nl0(ib)>0) then
   Hvbar = 1.0
else
   if (count(isnotmsr(Hv))>0) then
      Hvbar = sum(Hv,mask=isnotmsr(Hv))/float(count(isnotmsr(Hv)))
   else
     return
   end if
end if
if (lct%nscales>1) Hvbar = Hvbar*Hscale

! Allocation
mindata%nx = lct%nscales*lct%ncomp+lct%nscales
mindata%ny = nam%nc*bpar%nl0(ib)
allocate(mindata%x(mindata%nx))
allocate(mindata%guess(mindata%nx))
allocate(mindata%norm(mindata%nx))
allocate(mindata%binf(mindata%nx))
allocate(mindata%bsup(mindata%nx))
allocate(mindata%obs(mindata%ny))
allocate(mindata%dx(nam%nc,bpar%nl0(ib)))
allocate(mindata%dy(nam%nc,bpar%nl0(ib)))
allocate(mindata%dz(bpar%nl0(ib)))
allocate(mindata%dmask(nam%nc,bpar%nl0(ib)))
call lct_alloc(hdata,lct_guess)
call lct_alloc(hdata,lct_norm)
call lct_alloc(hdata,lct_binf)
call lct_alloc(hdata,lct_bsup)

! Define norm and bounds
lct_guess%H(1,1:3) = (/Hhbar,Hhbar,Hvbar/)
lct_norm%H(1,1:3) = (/Hhbar,Hhbar,Hvbar/)
lct_binf%H(1,1:3) = (/1.0/sqrt(Hscale),1.0/sqrt(Hscale),1.0/sqrt(Hscale)/)*lct_guess%H(1,1:3)
lct_bsup%H(1,1:3) = (/sqrt(Hscale),sqrt(Hscale),sqrt(Hscale)/)*lct_guess%H(1,1:3)
do iscales=2,lct%nscales
   lct_guess%H(iscales,1:3) = lct_guess%H(1,1:3)/Hscale**(iscales-1)
   lct_norm%H(iscales,1:3) = lct_norm%H(1,1:3)/Hscale**(iscales-1)
   lct_binf%H(iscales,1:3) = lct_binf%H(1,1:3)/Hscale**(iscales-1)
   lct_bsup%H(iscales,1:3) = lct_bsup%H(1,1:3)/Hscale**(iscales-1)
end do
if (lct%ncomp==4) then
   lct_guess%H(:,4) = 0.0
   lct_norm%H(:,4) = 1.0
   lct_binf%H(:,4) = -1.0
   lct_bsup%H(:,4) = 1.0
end if
do iscales=1,lct%nscales
   lct_guess%coef(iscales) = 1.0/float(lct%nscales)
   lct_norm%coef(iscales) = 1.0/float(lct%nscales)
   lct_binf%coef(iscales) = 0.0
   lct_bsup%coef(iscales) = 1.0
end do

! Fill mindata
offset = 0
do iscales=1,lct%nscales
   mindata%guess(offset+1:offset+lct%ncomp) = lct_guess%H(iscales,:)
   mindata%norm(offset+1:offset+lct%ncomp) = lct_norm%H(iscales,:)
   mindata%binf(offset+1:offset+lct%ncomp) = lct_binf%H(iscales,:)
   mindata%bsup(offset+1:offset+lct%ncomp) = lct_bsup%H(iscales,:)
   offset = offset+lct%ncomp
end do
mindata%guess(offset+1:offset+lct%nscales) = lct_guess%coef
mindata%norm(offset+1:offset+lct%nscales) = lct_norm%coef
mindata%binf(offset+1:offset+lct%nscales) = lct_binf%coef
mindata%bsup(offset+1:offset+lct%nscales) = lct_bsup%coef
mindata%obs = pack(lct%raw,.true.)
mindata%fit_type = trim(nam%fit_type)
mindata%nc = nam%nc
mindata%nl0 = bpar%nl0(ib)
mindata%dx = dx
mindata%dy = dy
mindata%dz = dz
mindata%dmask = dmask
mindata%nscales = lct%nscales
mindata%ncomp = lct%ncomp

! Compute fit
call minim(mindata,func,.false.)

! Copy parameters
offset = 0
do iscales=1,mindata%nscales
   lct%H(iscales,:) = mindata%x(offset+1:offset+mindata%ncomp)
   offset = offset+mindata%ncomp
end do
lct%coef = mindata%x(offset+1:offset+mindata%nscales)

! Dummy call to avoid warnings
call dummy(mindata)

! Fixed positive value for the 2D case
if (bpar%nl0(ib)==1) lct%H(:,3) = 1.0

! Check positive-definiteness
spd = .true.
do iscales=1,lct%nscales
   if (lct%ncomp==3) then
      det = lct%H(iscales,1)*lct%H(iscales,2)
   else
      det = lct%H(iscales,1)*lct%H(iscales,2)-lct%H(iscales,4)**2
   end if
   det = det*lct%H(iscales,3)
   spd = spd.and.(det>0.0)
   
   if (lct%coef(iscales)<0.0) lct%coef(iscales) = 0.0
end do
if (lct%nscales==1) then
   lct%coef(1) = 1.0
else
   lct%coef(lct%nscales) = 1.0-sum(lct%coef(1:lct%nscales-1))
end if
if (spd) then
   ! Rebuild fit
   call define_fit(nam%nc,bpar%nl0(ib),dx,dy,dz,dmask,lct%nscales,lct%ncomp,lct%H,lct%coef,lct%fit)
else
   ! Set as missing
   call msr(lct%H)
   call msr(lct%fit)
end if

! End associate
end associate

end subroutine compute_fit_lct

!----------------------------------------------------------------------
! Subroutine: compute_fit_lct
!> Purpose: compute a semi-positive definite fit of a raw function
!----------------------------------------------------------------------
subroutine compute_fit_lct_multi(hdata,ib,lct)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                              !< HDIAG data
integer,intent(in) :: ib                                         !< Block index
type(lcttype),intent(inout) :: lct(hdata%nam%nc1,hdata%geom%nl0) !< LCT

! Local variables
integer :: jl0,ic1,il0,il0r,ic,npack,progint
integer :: iproc,ic1_s(mpl%nproc),ic1_e(mpl%nproc),nc1_loc(mpl%nproc),ic1_loc
real(kind_real),allocatable :: dx(:,:),dy(:,:),dz(:)
real(kind_real),allocatable :: rbuf(:),sbuf(:)
logical,allocatable :: dmask(:,:),done(:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! MPI splitting
do iproc=1,mpl%nproc
   ic1_s(iproc) = (iproc-1)*(nam%nc1/mpl%nproc+1)+1
   ic1_e(iproc) = min(iproc*(nam%nc1/mpl%nproc+1),nam%nc1)
   nc1_loc(iproc) = ic1_e(iproc)-ic1_s(iproc)+1
end do

! Allocation
npack = lct(1,1)%npack
allocate(dx(nam%nc,bpar%nl0(ib)))
allocate(dy(nam%nc,bpar%nl0(ib)))
allocate(dz(bpar%nl0(ib)))
allocate(dmask(nam%nc,bpar%nl0(ib)))
allocate(rbuf(nam%nc1*npack))
allocate(done(nc1_loc(mpl%myproc)))

! Loop over levels
do jl0=1,geom%nl0
   write(mpl%unit,'(a13,a,i3,a)',advance='no') '','Level ',nam%levs(jl0),':'

   ! Loop over points
   call prog_init(progint,done)
   do ic1_loc=1,nc1_loc(mpl%myproc)
      ic1 = ic1_s(mpl%myproc)+ic1_loc-1

      ! Prepare vectors
      do il0r=1,bpar%nl0(ib)
         il0 = bpar%il0rjl0ib_to_il0(il0r,jl0,ib)
         do ic=1,nam%nc
            dmask(ic,il0r) = hdata%ic1il0_log(ic1,jl0).and.hdata%ic1icil0_log(ic1,ic,il0)
            if (dmask(ic,il0r)) then
               dx(ic,il0r) = lonmod(geom%lon(hdata%ic1icil0_to_ic0(ic1,ic,jl0))-geom%lon(hdata%ic1icil0_to_ic0(ic1,1,jl0))) &
                           & /cos(geom%lat(hdata%ic1icil0_to_ic0(ic1,1,jl0)))
               dy(ic,il0r) = geom%lat(hdata%ic1icil0_to_ic0(ic1,ic,jl0))-geom%lat(hdata%ic1icil0_to_ic0(ic1,1,jl0))
            end if
         end do
         dz(il0r) = float(nam%levs(il0)-nam%levs(jl0))
      end do

      ! Compute fit
      call compute_fit_lct(hdata,ib,dx,dy,dz,dmask,lct(ic1,jl0))

      ! Print progression
      done(ic1_loc) = .true.
      call prog_print(progint,done)
   end do
   write(mpl%unit,'(a)') '100%'

   ! Communication
   if (mpl%main) then
      do iproc=1,mpl%nproc
         if (iproc==mpl%ioproc) then
            ! Format data
            do ic1_loc=1,nc1_loc(iproc)
               ic1 = ic1_s(iproc)+ic1_loc-1
               call lct_pack(hdata,ib,lct(ic1,jl0),rbuf((ic1-1)*npack+1:ic1*npack))
            end do
         else
            ! Receive data on ioproc
            call mpl_recv(nc1_loc(iproc)*npack, &
          & rbuf((ic1_s(iproc)-1)*npack+1:ic1_e(iproc)*npack),iproc,mpl%tag)
         end if
      end do
   else
      ! Allocation
      allocate(sbuf(nc1_loc(mpl%myproc)*npack))

      ! Format data
      do ic1_loc=1,nc1_loc(mpl%myproc)
         ic1 = ic1_s(mpl%myproc)+ic1_loc-1
         call lct_pack(hdata,ib,lct(ic1,jl0),sbuf((ic1_loc-1)*npack+1:ic1_loc*npack))
      end do

      ! Send data to ioproc
      call mpl_send(nc1_loc(mpl%myproc)*npack,sbuf,mpl%ioproc,mpl%tag)

      ! Release memory
      deallocate(sbuf)
   end if
   mpl%tag = mpl%tag+1

   ! Broadcast data
   call mpl_bcast(rbuf,mpl%ioproc)

   ! Format data
   do ic1=1,nam%nc1
      call lct_unpack(hdata,ib,lct(ic1,jl0),rbuf((ic1-1)*npack+1:ic1*npack))
   end do
end do

! Release memory
deallocate(dx)
deallocate(dy)
deallocate(dz)
deallocate(dmask)
deallocate(rbuf)
deallocate(done)

! End associate
end associate

end subroutine compute_fit_lct_multi

!----------------------------------------------------------------------
! Subroutine: define_fit
!> Purpose: define the fit
!----------------------------------------------------------------------
subroutine define_fit(nc,nl0,dx,dy,dz,dmask,nscales,ncomp,H,coef,fit)

implicit none

! Passed variables
integer,intent(in) :: nc                       !< Number of classes
integer,intent(in) :: nl0                      !< Number of levels
real(kind_real),intent(in) :: dx(nc,nl0)       !< Zonal separation
real(kind_real),intent(in) :: dy(nc,nl0)       !< Meridian separation
real(kind_real),intent(in) :: dz(nl0)          !< Vertical separation
logical,intent(in) :: dmask(nc,nl0)            !< Mask
integer,intent(in) :: nscales                  !< Number of LCT scales
integer,intent(in) :: ncomp                    !< Number of LCT components
real(kind_real),intent(in) :: H(nscales,ncomp) !< LCT components
real(kind_real),intent(in) :: coef(nscales)  !< LCT coefficients
real(kind_real),intent(out) :: fit(nc,nl0)     !< Fit

! Local variables
integer :: il0,ic,iscales
real(kind_real) :: H11,H22,H33,Hc12,rsq

! Initialization
call msr(fit)

do iscales=1,nscales
   ! Force positive definiteness
   H11 = max(Hmin,H(iscales,1))
   H22 = max(Hmin,H(iscales,2))
   H33 = max(Hmin,H(iscales,3))
   call msr(Hc12)
   if (ncomp==4) Hc12 = max(-1.0_kind_real,min(H(iscales,4),1.0_kind_real))

   ! Homogeneous anisotropic approximation
   do il0=1,nl0
      do ic=1,nc
         if (dmask(ic,il0)) then
            ! Initialization
            if (iscales==1) fit(ic,il0) = 0.0
               
            ! Squared distance
            rsq = H11*dx(ic,il0)**2+H22*dy(ic,il0)**2+H33*dz(il0)**2
            if (ncomp==4) rsq = rsq+2.0*sqrt(H11*H22)*Hc12*dx(ic,il0)*dy(ic,il0)

            if (M==0) then
               ! Gaussian function
               fit(ic,il0) = fit(ic,il0)+coef(iscales)*exp(-0.5*rsq)
            else
               ! Matern function
               fit(ic,il0) = fit(ic,il0)+coef(iscales)*matern(M,sqrt(rsq))
            end if
         end if
      end do
   end do
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
integer :: ix,iscales,offset
real(kind_real) :: H(mindata%nscales,mindata%ncomp),coef(mindata%nscales)
real(kind_real) :: fit(mindata%nc,mindata%nl0)
real(kind_real) :: xtmp(mindata%nx),fit_pack(mindata%ny),xx

! Renormalize
xtmp = x*mindata%norm

! Get data
offset = 0
do iscales=1,mindata%nscales
   H(iscales,:) = xtmp(offset+1:offset+mindata%ncomp)
   offset = offset+mindata%ncomp
end do
coef = xtmp(offset+1:offset+mindata%nscales)

! Compute function
call define_fit(mindata%nc,mindata%nl0,mindata%dx,mindata%dy,mindata%dz,mindata%dmask,mindata%nscales,mindata%ncomp,H,coef,fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Cost
f = sum((mindata%obs-fit_pack)**2,mask=isnotmsr(fit_pack))

! Bound penalty
do ix=1,mindata%nx
   xx = (xtmp(ix)-mindata%binf(ix))/(mindata%bsup(ix)-mindata%binf(ix))
   if (xx<0.0) then
      f = f+mindata%f_guess*xx**2
   elseif (xx>1.0) then
      f = f+mindata%f_guess*(xx-1.0)**2
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

end module module_fit_lct
