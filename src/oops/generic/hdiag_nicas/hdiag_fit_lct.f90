!----------------------------------------------------------------------
! Module: hdiag_fit_lct.f90
!> Purpose: LCT fit routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_fit_lct

use omp_lib
use tools_const, only: lonlatmod
use tools_diffusion, only: matern
use tools_display, only: msgerror,msgwarning,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_minim, only: minim
use tools_missing, only: msr,isnotmsr
use type_hdata, only: hdatatype
use type_lct, only: lcttype,lct_alloc,lct_pack,lct_unpack
use type_mdata, only: mdatatype
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast,mpl_split

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
type(hdatatype),intent(in) :: hdata                                 !< HDIAG data
integer,intent(in) :: ib                                            !< Block index
real(kind_real),intent(in) :: dx(hdata%nam%nc3,hdata%bpar%nl0r(ib)) !< Zonal separation
real(kind_real),intent(in) :: dy(hdata%nam%nc3,hdata%bpar%nl0r(ib)) !< Meridian separation
real(kind_real),intent(in) :: dz(hdata%bpar%nl0r(ib))               !< Vertical separation
logical,intent(in) :: dmask(hdata%nam%nc3,hdata%bpar%nl0r(ib))      !< Mask
type(lcttype),intent(inout) :: lct                                  !< LCT

! Local variables
integer :: jl0r,jc3,iscales,offset
real(kind_real) :: distsq,Hh(hdata%nam%nc3),Hv(hdata%bpar%nl0r(ib)),Hhbar,Hvbar,det
logical :: spd
type(lcttype) :: lct_guess,lct_norm,lct_binf,lct_bsup
type(mdatatype) :: mdata

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Approximate homogeneous horizontal length-scale
call msr(Hh)
do jl0r=1,bpar%nl0r(ib)
   if (.not.(abs(dz(jl0r))>0.0)) then
      do jc3=1,nam%nc3
         if (dmask(jc3,jl0r)) then
            distsq = dx(jc3,jl0r)**2+dy(jc3,jl0r)**2
            if ((lct%raw(jc3,jl0r)>0.0).and.(distsq>0.0)) Hh(jc3) = -2.0*log(lct%raw(jc3,jl0r))/distsq
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
jc3 = 1
do jl0r=1,bpar%nl0r(ib)
   distsq = dz(jl0r)**2
   if ((lct%raw(jc3,jl0r)>0.0).and.(distsq>0.0)) Hv(jl0r) = -2.0*log(lct%raw(jc3,jl0r))/distsq
end do
if (bpar%nl0r(ib)>0) then
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
mdata%nx = sum(lct%ncomp)+lct%nscales
mdata%ny = nam%nc3*bpar%nl0r(ib)
allocate(mdata%x(mdata%nx))
allocate(mdata%guess(mdata%nx))
allocate(mdata%norm(mdata%nx))
allocate(mdata%binf(mdata%nx))
allocate(mdata%bsup(mdata%nx))
allocate(mdata%obs(mdata%ny))
allocate(mdata%dx(nam%nc3,bpar%nl0r(ib)))
allocate(mdata%dy(nam%nc3,bpar%nl0r(ib)))
allocate(mdata%dz(bpar%nl0r(ib)))
allocate(mdata%dmask(nam%nc3,bpar%nl0r(ib)))
allocate(mdata%ncomp(lct%nscales))
call lct_alloc(hdata,lct_guess)
call lct_alloc(hdata,lct_norm)
call lct_alloc(hdata,lct_binf)
call lct_alloc(hdata,lct_bsup)

! Define norm and bounds
offset = 0
do iscales=1,lct%nscales
   lct_guess%H(offset+1:offset+3) = (/Hhbar,Hhbar,Hvbar/)/Hscale**(iscales-1)
   lct_norm%H(offset+1:offset+3) = (/Hhbar,Hhbar,Hvbar/)/Hscale**(iscales-1)
   lct_binf%H(offset+1:offset+3) = (/1.0/sqrt(Hscale),1.0/sqrt(Hscale),1.0/sqrt(Hscale)/)*lct_guess%H(1:3)/Hscale**(iscales-1)
   lct_bsup%H(offset+1:offset+3) = (/sqrt(Hscale),sqrt(Hscale),sqrt(Hscale)/)*lct_guess%H(1:3)/Hscale**(iscales-1)
   offset = offset+3
   if (lct%ncomp(iscales)==4) then
      lct_guess%H(offset+1) = 0.0
      lct_norm%H(offset+1) = 1.0
      lct_binf%H(offset+1) = -1.0
      lct_bsup%H(offset+1) = 1.0
      offset = offset+1
   end if
   lct_guess%coef(iscales) = 1.0/float(lct%nscales)
   lct_norm%coef(iscales) = 1.0/float(lct%nscales)
   lct_binf%coef(iscales) = 0.0
   lct_bsup%coef(iscales) = 1.0
end do

! Fill mdata
mdata%guess(1:sum(lct%ncomp)) = lct_guess%H
mdata%norm(1:sum(lct%ncomp)) = lct_norm%H
mdata%binf(1:sum(lct%ncomp)) = lct_binf%H
mdata%bsup(1:sum(lct%ncomp)) = lct_bsup%H
mdata%guess(sum(lct%ncomp)+1:sum(lct%ncomp)+lct%nscales) = lct_guess%coef
mdata%norm(sum(lct%ncomp)+1:sum(lct%ncomp)+lct%nscales) = lct_norm%coef
mdata%binf(sum(lct%ncomp)+1:sum(lct%ncomp)+lct%nscales) = lct_binf%coef
mdata%bsup(sum(lct%ncomp)+1:sum(lct%ncomp)+lct%nscales) = lct_bsup%coef
mdata%obs = pack(lct%raw,.true.)
mdata%fit_type = trim(nam%fit_type)
mdata%nc3 = nam%nc3
mdata%nl0 = bpar%nl0r(ib)
mdata%dx = dx
mdata%dy = dy
mdata%dz = dz
mdata%dmask = dmask
mdata%nscales = lct%nscales
mdata%ncomp = lct%ncomp

! Compute fit
call minim(mdata,func,.false.)

! Copy parameters
lct%H = mdata%x(1:sum(lct%ncomp))
lct%coef = mdata%x(sum(lct%ncomp)+1:sum(lct%ncomp)+lct%nscales)

! Dummy call to avoid warnings
call dummy(mdata)

! Fixed positive value for the 2D case
if (bpar%nl0r(ib)==1) then
   offset = 0
   do iscales=1,lct%nscales
      lct%H(offset+3) = 1.0
      offset = offset+lct%ncomp(iscales)
   end do
end if

! Check positive-definiteness
spd = .true.
do iscales=1,lct%nscales
   offset = 0
   if (lct%ncomp(iscales)==3) then
      det = lct%H(offset+1)*lct%H(offset+2)
   else
      det = lct%H(offset+1)*lct%H(offset+2)-lct%H(offset+4)**2
   end if
   det = det*lct%H(offset+3)
   spd = spd.and.(det>0.0)
   if (lct%coef(iscales)<0.0) lct%coef(iscales) = 0.0
   offset = offset+lct%ncomp(iscales)
end do
if (lct%nscales==1) then
   lct%coef(1) = 1.0
else
   lct%coef(lct%nscales) = 1.0-sum(lct%coef(1:lct%nscales-1))
end if
if (spd) then
   ! Rebuild fit
   call define_fit(nam%nc3,bpar%nl0r(ib),dx,dy,dz,dmask,lct%nscales,lct%ncomp,lct%H,lct%coef,lct%fit)
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
integer :: il0,jl0r,jl0,ic1,jc3,npack,progint,iproc
integer :: ic1_s(mpl%nproc),ic1_e(mpl%nproc),nc1_loc(mpl%nproc),ic1_loc
real(kind_real),allocatable :: dx(:,:),dy(:,:),dz(:)
real(kind_real),allocatable :: rbuf(:),sbuf(:)
logical,allocatable :: dmask(:,:),done(:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! MPI splitting
call mpl_split(nam%nc1,ic1_s,ic1_e,nc1_loc)

! Allocation
npack = lct(1,1)%npack
allocate(dx(nam%nc3,bpar%nl0r(ib)))
allocate(dy(nam%nc3,bpar%nl0r(ib)))
allocate(dz(bpar%nl0r(ib)))
allocate(dmask(nam%nc3,bpar%nl0r(ib)))
allocate(sbuf(nc1_loc(mpl%myproc)*npack))
allocate(rbuf(nam%nc1*npack))
allocate(done(nc1_loc(mpl%myproc)))

! Loop over levels
do il0=1,geom%nl0
   write(mpl%unit,'(a13,a,i3,a)',advance='no') '','Level ',nam%levs(il0),':'

   ! Loop over points
   call prog_init(progint,done)
   do ic1_loc=1,nc1_loc(mpl%myproc)
      ic1 = ic1_s(mpl%myproc)+ic1_loc-1

      ! Prepare vectors
      do jl0r=1,bpar%nl0r(ib)
         jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
         do jc3=1,nam%nc3
            dmask(jc3,jl0r) = hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)
            if (dmask(jc3,jl0r)) then
               dx(jc3,jl0r) = geom%lon(hdata%c1c3_to_c0(ic1,jc3))-geom%lon(hdata%c1c3_to_c0(ic1,1))
               dy(jc3,jl0r) = geom%lat(hdata%c1c3_to_c0(ic1,jc3))-geom%lat(hdata%c1c3_to_c0(ic1,1))
               call lonlatmod(dx(jc3,jl0r),dy(jc3,jl0r))
               dx(jc3,jl0r) = dx(jc3,jl0r)/cos(geom%lat(hdata%c1c3_to_c0(ic1,1)))
            end if
         end do
         dz(jl0r) = float(nam%levs(jl0)-nam%levs(il0))
      end do

      ! Compute fit
      call compute_fit_lct(hdata,ib,dx,dy,dz,dmask,lct(ic1,il0))

      ! Print progression
      done(ic1_loc) = .true.
      call prog_print(progint,done)
   end do
   write(mpl%unit,'(a)') '100%'

   ! Prepare buffer
   do ic1_loc=1,nc1_loc(mpl%myproc)
      ic1 = ic1_s(mpl%myproc)+ic1_loc-1
      call lct_pack(hdata,ib,lct(ic1,il0),sbuf((ic1_loc-1)*npack+1:ic1_loc*npack))
   end do

   ! Communication
   if (mpl%main) then
      do iproc=1,mpl%nproc
         if (nc1_loc(iproc)*npack>0) then
            if (iproc==mpl%ioproc) then
               ! Copy data
               rbuf(ic1_s(iproc)*npack:ic1_e(iproc)*npack) = sbuf
            else
               ! Receive data on ioproc
               call mpl_recv(nc1_loc(iproc)*npack,rbuf(ic1_s(iproc)*npack:ic1_e(iproc)*npack),iproc,mpl%tag)
            end if
         end if
      end do
   else
      if (nc1_loc(mpl%myproc)*npack>0) then
         ! Send data to ioproc
         call mpl_send(nc1_loc(mpl%myproc)*npack,sbuf,mpl%ioproc,mpl%tag)
      end if
   end if
   mpl%tag = mpl%tag+1

   ! Broadcast data
   call mpl_bcast(rbuf,mpl%ioproc)

   ! Format data
   do ic1=1,nam%nc1
      call lct_unpack(hdata,ib,lct(ic1,il0),rbuf((ic1-1)*npack+1:ic1*npack))
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
integer,intent(in) :: nc                    !< Number of classes
integer,intent(in) :: nl0                   !< Number of levels
real(kind_real),intent(in) :: dx(nc,nl0)    !< Zonal separation
real(kind_real),intent(in) :: dy(nc,nl0)    !< Meridian separation
real(kind_real),intent(in) :: dz(nl0)       !< Vertical separation
logical,intent(in) :: dmask(nc,nl0)         !< Mask
integer,intent(in) :: nscales               !< Number of LCT scales
integer,intent(in) :: ncomp(nscales)        !< Number of LCT components
real(kind_real),intent(in) :: H(sum(ncomp)) !< LCT components
real(kind_real),intent(in) :: coef(nscales) !< LCT coefficients
real(kind_real),intent(out) :: fit(nc,nl0)  !< Fit

! Local variables
integer :: jl0,jc3,iscales,offset
real(kind_real) :: H11,H22,H33,Hc12,rsq

! Initialization
call msr(fit)

offset = 0
do iscales=1,nscales
   ! Force positive definiteness
   H11 = max(Hmin,H(offset+1))
   H22 = max(Hmin,H(offset+2))
   H33 = max(Hmin,H(offset+3))
   call msr(Hc12)
   if (ncomp(iscales)==4) Hc12 = max(-1.0_kind_real,min(H(offset+4),1.0_kind_real))

   ! Homogeneous anisotropic approximation
   do jl0=1,nl0
      do jc3=1,nc
         if (dmask(jc3,jl0)) then
            ! Initialization
            if (iscales==1) fit(jc3,jl0) = 0.0

            ! Squared distance
            rsq = H11*dx(jc3,jl0)**2+H22*dy(jc3,jl0)**2+H33*dz(jl0)**2
            if (ncomp(iscales)==4) rsq = rsq+2.0*sqrt(H11*H22)*Hc12*dx(jc3,jl0)*dy(jc3,jl0)

            if (M==0) then
               ! Gaussian function
               fit(jc3,jl0) = fit(jc3,jl0)+coef(iscales)*exp(-0.5*rsq)
            else
               ! Matern function
               fit(jc3,jl0) = fit(jc3,jl0)+coef(iscales)*matern(M,sqrt(rsq))
            end if
         end if
      end do
   end do
   offset = offset+ncomp(iscales)
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
integer :: ix
real(kind_real) :: fit(mdata%nc3,mdata%nl0)
real(kind_real) :: xtmp(mdata%nx),fit_pack(mdata%ny),xx

! Renormalize
xtmp = x*mdata%norm

! Compute function
call define_fit(mdata%nc3,mdata%nl0,mdata%dx,mdata%dy,mdata%dz,mdata%dmask,mdata%nscales,mdata%ncomp, &
 & xtmp(1:sum(mdata%ncomp)),xtmp(sum(mdata%ncomp)+1:sum(mdata%ncomp)+mdata%nscales),fit)

! Pack
fit_pack = pack(fit,mask=.true.)

! Cost
f = sum((mdata%obs-fit_pack)**2,mask=isnotmsr(fit_pack))

! Bound penalty
do ix=1,mdata%nx
   xx = (xtmp(ix)-mdata%binf(ix))/(mdata%bsup(ix)-mdata%binf(ix))
   if (xx<0.0) then
      f = f+mdata%f_guess*xx**2
   elseif (xx>1.0) then
      f = f+mdata%f_guess*(xx-1.0)**2
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

end module hdiag_fit_lct
