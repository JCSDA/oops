!----------------------------------------------------------------------
! Module: type_lct_blk
!> Purpose: LCT data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_lct_blk

use omp_lib
use tools_func, only: lonlatmod,fit_lct
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_minim, only: minim_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! LCT block data derived type
type lct_blk_type
   ! Attributes
   integer :: ic1a                          !< Local index
   integer :: il0                           !< Level index
   integer :: ib                            !< Block index
   integer :: nscales                       !< Number of LCT scales
   integer,allocatable :: ncomp(:)          !< Number of LCT components

   ! Data
   real(kind_real),allocatable :: H(:)      !< LCT components
   real(kind_real),allocatable :: coef(:)   !< LCT coefficients
   real(kind_real),allocatable :: raw(:,:)  !< Raw correlations
   real(kind_real),allocatable :: norm(:,:) !< Norm to take nsub into account
   real(kind_real),allocatable :: fit(:,:)  !< Fitted correlations
contains
   procedure :: lct_blk_alloc_base
   procedure :: lct_blk_alloc_block
   generic :: alloc => lct_blk_alloc_base,lct_blk_alloc_block
   procedure :: dealloc => lct_blk_dealloc
   procedure :: correlation => lct_blk_correlation
   procedure :: fitting => lct_blk_fitting
end type lct_blk_type

real(kind_real),parameter :: Hscale = 10.0  !< Typical factor between LCT scales
logical :: lprt = .false.                   !< Optimization print

private
public :: lct_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: lct_blk_alloc_base
!> Purpose: lct_blk object base allocation
!----------------------------------------------------------------------
subroutine lct_blk_alloc_base(lct_blk,nam)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk !< LCT block
type(nam_type),intent(in) :: nam             !< Namelist

! Local variables
integer :: iscales

! Attributes
lct_blk%nscales = nam%lct_nscales
allocate(lct_blk%ncomp(lct_blk%nscales))
do iscales=1,lct_blk%nscales
   if (nam%lct_diag(iscales)) then
      lct_blk%ncomp(iscales) = 3
   else
      lct_blk%ncomp(iscales) = 4
   end if
end do

! Allocation
allocate(lct_blk%H(sum(lct_blk%ncomp)))
allocate(lct_blk%coef(lct_blk%nscales))

! Initialization
call msr(lct_blk%H)
call msr(lct_blk%coef)

end subroutine lct_blk_alloc_base

!----------------------------------------------------------------------
! Subroutine: lct_blk_alloc_block
!> Purpose: lct_blk object block allocation
!----------------------------------------------------------------------
subroutine lct_blk_alloc_block(lct_blk,nam,bpar,ic1a,il0,ib)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk !< LCT block
type(nam_type),intent(in) :: nam             !< Namelist
type(bpar_type),intent(in) :: bpar           !< Block parameters
integer,intent(in) :: ic1a                   !< Local index
integer,intent(in) :: il0                    !< Level index
integer,intent(in) :: ib                     !< Block index

! Basic allocation
call lct_blk%alloc(nam)

! Attributes
lct_blk%ic1a = ic1a
lct_blk%il0 = il0
lct_blk%ib = ib

! Allocation
allocate(lct_blk%raw(nam%nc3,bpar%nl0r(ib)))
allocate(lct_blk%norm(nam%nc3,bpar%nl0r(ib)))
allocate(lct_blk%fit(nam%nc3,bpar%nl0r(ib)))

! Initialization
lct_blk%raw = 0.0
lct_blk%norm = 0.0
call msr(lct_blk%fit)

end subroutine lct_blk_alloc_block

!----------------------------------------------------------------------
! Subroutine: lct_blk_dealloc
!> Purpose: lct_blk object deallocation
!----------------------------------------------------------------------
subroutine lct_blk_dealloc(lct_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk !< LCT block

! Release memory
if (allocated(lct_blk%ncomp)) deallocate(lct_blk%ncomp)
if (allocated(lct_blk%H)) deallocate(lct_blk%H)
if (allocated(lct_blk%coef)) deallocate(lct_blk%coef)
if (allocated(lct_blk%raw)) deallocate(lct_blk%raw)
if (allocated(lct_blk%norm)) deallocate(lct_blk%norm)
if (allocated(lct_blk%fit)) deallocate(lct_blk%fit)

end subroutine lct_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_blk_correlation
!> Purpose: compute raw correlation
!----------------------------------------------------------------------
subroutine lct_blk_correlation(lct_blk,nam,bpar,mom_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk !< LCT block
type(nam_type),intent(in) :: nam             !< Namelist
type(bpar_type),intent(in) :: bpar           !< Block parameters
type(mom_blk_type),intent(in) :: mom_blk     !< Moments block

! Local variables
integer :: jsub,jl0,jc3
real(kind_real) :: den

! Associate
associate(ic1a=>lct_blk%ic1a,il0=>lct_blk%il0,ib=>lct_blk%ib)

! Sum
do jsub=1,mom_blk%nsub
   do jl0=1,bpar%nl0r(ib)
      do jc3=1,nam%nc3
         den = mom_blk%m2_1(ic1a,jc3,jl0,il0,jsub)*mom_blk%m2_2(ic1a,jc3,jl0,il0,jsub)
         if (den>0.0) then
            lct_blk%raw(jc3,jl0) = lct_blk%raw(jc3,jl0)+mom_blk%m11(ic1a,jc3,jl0,il0,jsub)/sqrt(den)
            lct_blk%norm(jc3,jl0) = lct_blk%norm(jc3,jl0)+1.0
         end if
      end do
   end do
end do

! Normalize
do jl0=1,bpar%nl0r(ib)
   do jc3=1,nam%nc3
      if (lct_blk%norm(jc3,jl0)>0.0) lct_blk%raw(jc3,jl0) = lct_blk%raw(jc3,jl0)/lct_blk%norm(jc3,jl0)
   end do
end do

! End associate
end associate

end subroutine lct_blk_correlation

!----------------------------------------------------------------------
! Subroutine: lct_blk_fitting
!> Purpose: fitting LCT
!----------------------------------------------------------------------
subroutine lct_blk_fitting(lct_blk,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk !< LCT block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
type(hdata_type),intent(in) :: hdata         !< HDIAG data

! Local variables
integer :: jl0r,jl0,ic1,jc3,iscales,offset
real(kind_real) :: distsq,Hhbar,Hvbar
real(kind_real),allocatable :: Hh(:),Hv(:),dx(:,:),dy(:,:),dz(:)
logical :: spd
logical,allocatable :: dmask(:,:),done(:)
type(lct_blk_type) :: lct_guess,lct_norm,lct_binf,lct_bsup
type(minim_type) :: minim

! Associate
associate(ic1a=>lct_blk%ic1a,il0=>lct_blk%il0,ib=>lct_blk%ib)

! Allocation
allocate(Hh(nam%nc3))
allocate(Hv(bpar%nl0r(ib)))
allocate(dx(nam%nc3,bpar%nl0r(ib)))
allocate(dy(nam%nc3,bpar%nl0r(ib)))
allocate(dz(bpar%nl0r(ib)))
allocate(dmask(nam%nc3,bpar%nl0r(ib)))
allocate(done(hdata%nc1a))

! Global index
ic1 = hdata%c1a_to_c1(ic1a)

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

! Approximate homogeneous horizontal length-scale
call msr(Hh)
do jl0r=1,bpar%nl0r(ib)
   if (.not.(abs(dz(jl0r))>0.0)) then
      do jc3=1,nam%nc3
         if (dmask(jc3,jl0r)) then
            distsq = dx(jc3,jl0r)**2+dy(jc3,jl0r)**2
            if ((lct_blk%raw(jc3,jl0r)>0.0).and.(distsq>0.0)) Hh(jc3) = -2.0*log(lct_blk%raw(jc3,jl0r))/distsq
         end if
      end do
   end if
end do
if (count(isnotmsr(Hh))>0) then
   Hhbar = sum(Hh,mask=isnotmsr(Hh))/float(count(isnotmsr(Hh)))
else
   return
end if
if (lct_blk%nscales>1) Hhbar = Hhbar*Hscale

! Approximate homogeneous vertical length-scale
call msr(Hv)
jc3 = 1
do jl0r=1,bpar%nl0r(ib)
   distsq = dz(jl0r)**2
   if ((lct_blk%raw(jc3,jl0r)>0.0).and.(distsq>0.0)) Hv(jl0r) = -2.0*log(lct_blk%raw(jc3,jl0r))/distsq
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
if (lct_blk%nscales>1) Hvbar = Hvbar*Hscale

! Allocation
minim%nx = sum(lct_blk%ncomp)+lct_blk%nscales
minim%ny = nam%nc3*bpar%nl0r(ib)
allocate(minim%x(minim%nx))
allocate(minim%guess(minim%nx))
allocate(minim%norm(minim%nx))
allocate(minim%binf(minim%nx))
allocate(minim%bsup(minim%nx))
allocate(minim%obs(minim%ny))
allocate(minim%dx(nam%nc3,bpar%nl0r(ib)))
allocate(minim%dy(nam%nc3,bpar%nl0r(ib)))
allocate(minim%dz(bpar%nl0r(ib)))
allocate(minim%dmask(nam%nc3,bpar%nl0r(ib)))
allocate(minim%ncomp(lct_blk%nscales))
call lct_guess%alloc(nam)
call lct_norm%alloc(nam)
call lct_binf%alloc(nam)
call lct_bsup%alloc(nam)

! Define norm and bounds
offset = 0
do iscales=1,lct_blk%nscales
   lct_guess%H(offset+1:offset+3) = (/Hhbar,Hhbar,Hvbar/)/Hscale**(iscales-1)
   lct_norm%H(offset+1:offset+3) = (/Hhbar,Hhbar,Hvbar/)/Hscale**(iscales-1)
   lct_binf%H(offset+1:offset+3) = (/1.0/sqrt(Hscale),1.0/sqrt(Hscale),1.0/sqrt(Hscale)/)*lct_guess%H(1:3)/Hscale**(iscales-1)
   lct_bsup%H(offset+1:offset+3) = (/sqrt(Hscale),sqrt(Hscale),sqrt(Hscale)/)*lct_guess%H(1:3)/Hscale**(iscales-1)
   offset = offset+3
   if (lct_blk%ncomp(iscales)==4) then
      lct_guess%H(offset+1) = 0.0
      lct_norm%H(offset+1) = 1.0
      lct_binf%H(offset+1) = -1.0
      lct_bsup%H(offset+1) = 1.0
      offset = offset+1
   end if
   lct_guess%coef(iscales) = 1.0/float(lct_blk%nscales)
   lct_norm%coef(iscales) = 1.0/float(lct_blk%nscales)
   lct_binf%coef(iscales) = 0.0
   lct_bsup%coef(iscales) = 1.0
end do

! Fill minim
minim%guess(1:sum(lct_blk%ncomp)) = lct_guess%H
minim%norm(1:sum(lct_blk%ncomp)) = lct_norm%H
minim%binf(1:sum(lct_blk%ncomp)) = lct_binf%H
minim%bsup(1:sum(lct_blk%ncomp)) = lct_bsup%H
minim%guess(sum(lct_blk%ncomp)+1:sum(lct_blk%ncomp)+lct_blk%nscales) = lct_guess%coef
minim%norm(sum(lct_blk%ncomp)+1:sum(lct_blk%ncomp)+lct_blk%nscales) = lct_norm%coef
minim%binf(sum(lct_blk%ncomp)+1:sum(lct_blk%ncomp)+lct_blk%nscales) = lct_binf%coef
minim%bsup(sum(lct_blk%ncomp)+1:sum(lct_blk%ncomp)+lct_blk%nscales) = lct_bsup%coef
minim%obs = pack(lct_blk%raw,.true.)
minim%cost_function = 'fit_lct'
minim%algo = trim(nam%minim_algo)
minim%nc3 = nam%nc3
minim%nl0 = bpar%nl0r(ib)
minim%dx = dx
minim%dy = dy
minim%dz = dz
minim%dmask = dmask
minim%nscales = lct_blk%nscales
minim%ncomp = lct_blk%ncomp

! Compute fit
minim%cost_function = 'fit_lct'
call minim%compute(lprt)

! Copy parameters
lct_blk%H = minim%x(1:sum(lct_blk%ncomp))
lct_blk%coef = minim%x(sum(lct_blk%ncomp)+1:sum(lct_blk%ncomp)+lct_blk%nscales)

! Fixed positive value for the 2D case
if (bpar%nl0r(ib)==1) then
   offset = 0
   do iscales=1,lct_blk%nscales
      lct_blk%H(offset+3) = 1.0
      offset = offset+lct_blk%ncomp(iscales)
   end do
end if

! Check positive-definiteness
spd = .true.
offset = 0
do iscales=1,lct_blk%nscales
   spd = spd.and.(lct_blk%H(offset+1)>0.0).and.(lct_blk%H(offset+2)>0.0).and.(lct_blk%H(offset+3)>0.0)
   if (lct_blk%ncomp(iscales)==4) spd = spd.and.(lct_blk%H(offset+4)>-1.0).and.(lct_blk%H(offset+4)<1.0)
   if (iscales<lct_blk%nscales) spd = spd.and.(lct_blk%coef(iscales)>0.0)
   offset = offset+lct_blk%ncomp(iscales)
end do
if (spd) then
   ! Rebuild fit
   call fit_lct(nam%nc3,bpar%nl0r(ib),dx,dy,dz,dmask,lct_blk%nscales,lct_blk%ncomp,lct_blk%H,lct_blk%coef,lct_blk%fit)

   ! Last coefficient
   if (lct_blk%nscales==1) then
       lct_blk%coef(1) = 1.0
   else
      lct_blk%coef(lct_blk%nscales) = 1.0-sum(lct_blk%coef(1:lct_blk%nscales-1))
   end if
else
   ! Missing values
   call msr(lct_blk%H)
   call msr(lct_blk%coef)
   call msr(lct_blk%fit)
end if

! End associate
end associate

end subroutine lct_blk_fitting

end module type_lct_blk
