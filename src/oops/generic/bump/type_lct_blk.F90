!----------------------------------------------------------------------
! Module: type_lct_blk
! Purpose: LCT data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_lct_blk

!$ use omp_lib
use tools_func, only: lonlatmod,fit_lct,check_cond
use tools_kinds, only: kind_real
use tools_repro, only: rth,inf,sup
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_minim, only: minim_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

real(kind_real),parameter :: cor_min = 0.5_kind_real ! Minimum relevant correlation for first guess
real(kind_real),parameter :: Dscale = 10.0_kind_real ! Typical factor between diffusion scales
logical,parameter :: lprt = .false.                  ! Optimization print

! LCT block data derived type
type lct_blk_type
   ! Attributes
   integer :: ib                                    ! Block index
   integer :: nscales                               ! Number of LCT scales

   ! Correlation/variances
   real(kind_real),allocatable :: raw(:,:,:,:)      ! Raw correlations
   real(kind_real),allocatable :: m2(:,:)           ! Variances
   real(kind_real),allocatable :: m2flt(:,:)        ! Filtered variances

   ! Diffusion data
   real(kind_real),allocatable :: D(:,:,:,:)        ! Diffusion components
   real(kind_real),allocatable :: coef(:,:,:)       ! Multi-scale coefficients
   real(kind_real),allocatable :: fit(:,:,:,:)      ! Fitted correlations

   ! Filtered diffusion data
   real(kind_real),allocatable :: D_filt(:,:,:,:)   ! Diffusion components
   real(kind_real),allocatable :: coef_filt(:,:,:)  ! Multi-scale coefficients
   real(kind_real),allocatable :: fit_filt(:,:,:,:) ! Fitted correlations

   ! Output data
   real(kind_real),allocatable :: D11(:,:,:)        ! Daley tensor, component 11
   real(kind_real),allocatable :: D22(:,:,:)        ! Daley tensor, component 22
   real(kind_real),allocatable :: D33(:,:,:)        ! Daley tensor, component 33
   real(kind_real),allocatable :: D12(:,:,:)        ! Daley tensor, component 12
   real(kind_real),allocatable :: H11(:,:,:)        ! Local correlation tensor, component 11
   real(kind_real),allocatable :: H22(:,:,:)        ! Local correlation tensor, component 22
   real(kind_real),allocatable :: H33(:,:,:)        ! Local correlation tensor, component 33
   real(kind_real),allocatable :: H12(:,:,:)        ! Local correlation tensor, component 12
   real(kind_real),allocatable :: Dcoef(:,:,:)      ! Tensor coefficient
   real(kind_real),allocatable :: DLh(:,:,:)        ! Tensor length-scale
   real(kind_real),allocatable :: coef_ens(:,:)     ! Variances
contains
   procedure :: alloc => lct_blk_alloc
   procedure :: dealloc => lct_blk_dealloc
   procedure :: correlation => lct_blk_correlation
   procedure :: fitting => lct_blk_fitting
end type lct_blk_type

private
public :: lct_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: lct_blk_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine lct_blk_alloc(lct_blk,nam,geom,bpar,samp,ib)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling
integer,intent(in) :: ib                     ! Block index

! Attributes
lct_blk%ib = ib
lct_blk%nscales = nam%lct_nscales

! Allocation
allocate(lct_blk%raw(nam%nc3,bpar%nl0r(ib),samp%nc1a,geom%nl0))
allocate(lct_blk%m2(samp%nc1a,geom%nl0))
allocate(lct_blk%m2flt(samp%nc1a,geom%nl0))
allocate(lct_blk%D(4,lct_blk%nscales,samp%nc1a,geom%nl0))
allocate(lct_blk%coef(lct_blk%nscales,samp%nc1a,geom%nl0))
allocate(lct_blk%fit(nam%nc3,bpar%nl0r(ib),samp%nc1a,geom%nl0))
if (nam%diag_rhflt>0.0) then
   allocate(lct_blk%D_filt(4,lct_blk%nscales,samp%nc1a,geom%nl0))
   allocate(lct_blk%coef_filt(lct_blk%nscales,samp%nc1a,geom%nl0))
   allocate(lct_blk%fit_filt(nam%nc3,bpar%nl0r(ib),samp%nc1a,geom%nl0))
end if
allocate(lct_blk%D11(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D22(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D33(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D12(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H11(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H22(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H33(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H12(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%Dcoef(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%DLh(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%coef_ens(geom%nc0a,geom%nl0))

end subroutine lct_blk_alloc

!----------------------------------------------------------------------
! Subroutine: lct_blk_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine lct_blk_dealloc(lct_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block

! Release memory
if (allocated(lct_blk%raw)) deallocate(lct_blk%raw)
if (allocated(lct_blk%m2)) deallocate(lct_blk%m2)
if (allocated(lct_blk%m2flt)) deallocate(lct_blk%m2flt)
if (allocated(lct_blk%D)) deallocate(lct_blk%D)
if (allocated(lct_blk%coef)) deallocate(lct_blk%coef)
if (allocated(lct_blk%fit)) deallocate(lct_blk%fit)
if (allocated(lct_blk%D_filt)) deallocate(lct_blk%D_filt)
if (allocated(lct_blk%coef_filt)) deallocate(lct_blk%coef_filt)
if (allocated(lct_blk%fit_filt)) deallocate(lct_blk%fit_filt)
if (allocated(lct_blk%D11)) deallocate(lct_blk%D11)
if (allocated(lct_blk%D22)) deallocate(lct_blk%D22)
if (allocated(lct_blk%D33)) deallocate(lct_blk%D33)
if (allocated(lct_blk%D12)) deallocate(lct_blk%D12)
if (allocated(lct_blk%H11)) deallocate(lct_blk%H11)
if (allocated(lct_blk%H22)) deallocate(lct_blk%H22)
if (allocated(lct_blk%H33)) deallocate(lct_blk%H33)
if (allocated(lct_blk%H12)) deallocate(lct_blk%H12)
if (allocated(lct_blk%Dcoef)) deallocate(lct_blk%Dcoef)
if (allocated(lct_blk%DLh)) deallocate(lct_blk%DLh)
if (allocated(lct_blk%coef_ens)) deallocate(lct_blk%coef_ens)

end subroutine lct_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_blk_correlation
! Purpose: compute raw correlation
!----------------------------------------------------------------------
subroutine lct_blk_correlation(lct_blk,nam,geom,bpar,samp,mom_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling
type(mom_blk_type),intent(in) :: mom_blk     ! Moments block

! Local variables
integer :: jsub,il0,jl0r,jl0,jc3,ic1a,ic1
real(kind_real) :: den
real(kind_real),allocatable :: norm_raw(:,:,:,:),norm_m2(:,:)

! Associate
associate(ib=>lct_blk%ib)

! Allocation
allocate(norm_raw(nam%nc3,bpar%nl0r(ib),samp%nc1a,geom%nl0))
allocate(norm_m2(samp%nc1a,geom%nl0))

! Initialize
lct_blk%raw = 0.0
lct_blk%m2 = 0.0
norm_raw = 0.0
norm_m2 = 0.0

! Sum over jsub
do jsub=1,mom_blk%nsub
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
         do jc3=1,nam%nc3
            do ic1a=1,samp%nc1a
               ic1 = samp%c1a_to_c1(ic1a)
               if (samp%c1l0_log(ic1,il0).and.samp%c1c3l0_log(ic1,jc3,jl0)) then
                  den = mom_blk%m2_1(ic1a,il0,jsub)*mom_blk%m2_2(ic1a,jc3,jl0,jsub)
                  if (den>0.0) then
                     lct_blk%raw(jc3,jl0r,ic1a,il0) = lct_blk%raw(jc3,jl0r,ic1a,il0)+mom_blk%m11(ic1a,jc3,jl0r,il0,jsub)/sqrt(den)
                     norm_raw(jc3,jl0r,ic1a,il0) = norm_raw(jc3,jl0r,ic1a,il0)+1.0
                  end if
               end if
            end do
         end do
      end do
   end do
   do il0=1,geom%nl0
      do ic1a=1,samp%nc1a
         ic1 = samp%c1a_to_c1(ic1a)
         if (samp%c1l0_log(ic1,il0)) then
            lct_blk%m2(ic1a,il0) = lct_blk%m2(ic1a,il0)+mom_blk%m2_1(ic1a,il0,jsub)
            norm_m2(ic1a,il0) = norm_m2(ic1a,il0)+1.0
         end if
      end do
   end do
end do

! Normalize
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,nam%nc3
         do ic1a=1,samp%nc1a
            ic1 = samp%c1a_to_c1(ic1a)
            if (samp%c1l0_log(ic1,il0)) then
               if (norm_raw(jc3,jl0r,ic1a,il0)>0.0) lct_blk%raw(jc3,jl0r,ic1a,il0) = lct_blk%raw(jc3,jl0r,ic1a,il0) &
             & /norm_raw(jc3,jl0r,ic1a,il0)
            end if
         end do
      end do
   end do
end do
do il0=1,geom%nl0
   do ic1a=1,samp%nc1a
      ic1 = samp%c1a_to_c1(ic1a)
      if (samp%c1l0_log(ic1,il0)) then
         if (norm_m2(ic1a,il0)>0.0) lct_blk%m2(ic1a,il0) = lct_blk%m2(ic1a,il0)/norm_m2(ic1a,il0)
      end if
   end do
end do

! Release memory
deallocate(norm_raw)
deallocate(norm_m2)

! End associate
end associate

end subroutine lct_blk_correlation

!----------------------------------------------------------------------
! Subroutine: lct_blk_fitting
! Purpose: fitting LCT
!----------------------------------------------------------------------
subroutine lct_blk_fitting(lct_blk,mpl,nam,geom,bpar,samp)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling

! Local variables
integer :: il0,jl0r,jl0,ic1a,ic1,ic0,jc3,jc0,iscales,icomp
real(kind_real) :: distsq,Dhbar,Dvbar,diag_rescale
real(kind_real),allocatable :: Dh(:),Dv(:),dx(:,:),dy(:,:),dz(:,:)
logical :: valid
logical,allocatable :: dmask(:,:),Dv_valid(:)
type(minim_type) :: minim

! Associate
associate(ib=>lct_blk%ib)

! Allocation
allocate(Dh(nam%nc3))
allocate(Dv(bpar%nl0r(ib)))
allocate(dx(nam%nc3,bpar%nl0r(ib)))
allocate(dy(nam%nc3,bpar%nl0r(ib)))
allocate(dz(nam%nc3,bpar%nl0r(ib)))
allocate(dmask(nam%nc3,bpar%nl0r(ib)))
allocate(Dv_valid(bpar%nl0r(ib)))
minim%nx = lct_blk%nscales*4
if (lct_blk%nscales>1) minim%nx = minim%nx+lct_blk%nscales-1
minim%ny = nam%nc3*bpar%nl0r(ib)
allocate(minim%x(minim%nx))
allocate(minim%guess(minim%nx))
allocate(minim%binf(minim%nx))
allocate(minim%bsup(minim%nx))
allocate(minim%obs(minim%ny))
allocate(minim%dx(nam%nc3,bpar%nl0r(ib)))
allocate(minim%dy(nam%nc3,bpar%nl0r(ib)))
allocate(minim%dz(nam%nc3,bpar%nl0r(ib)))
allocate(minim%dmask(nam%nc3,bpar%nl0r(ib)))

do il0=1,geom%nl0
   ! Initialization
   write(mpl%info,'(a13,a,i3,a)') '','Level ',nam%levs(il0),':'
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc1a)

   do ic1a=1,samp%nc1a
      ! Global index
      ic1 = samp%c1a_to_c1(ic1a)

      if (samp%c1l0_log(ic1,il0)) then
         ! Compute deltas
         ic0 = samp%c1_to_c0(ic1)
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               dmask(jc3,jl0r) = samp%c1l0_log(ic1,il0).and.samp%c1c3l0_log(ic1,jc3,jl0)
               if (dmask(jc3,jl0r)) then
                  jc0 = samp%c1c3_to_c0(ic1,jc3)
                  call geom%compute_deltas(ic0,il0,jc0,jl0,dx(jc3,jl0r),dy(jc3,jl0r),dz(jc3,jl0r))
               end if
            end do
         end do

         ! Approximate homogeneous horizontal length-scale
         Dh = mpl%msv%valr
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            if (il0==jl0) then
               do jc3=1,nam%nc3
                  if (dmask(jc3,jl0r)) then
                     distsq = dx(jc3,jl0r)**2+dy(jc3,jl0r)**2
                     if (sup(lct_blk%raw(jc3,jl0r,ic1a,il0),cor_min).and.inf(lct_blk%raw(jc3,jl0r,ic1a,il0),1.0_kind_real) &
                   & .and.(distsq>0.0)) Dh(jc3) = -distsq/(2.0*log(lct_blk%raw(jc3,jl0r,ic1a,il0)))
                  end if
               end do
            end if
         end do
         Dhbar = mpl%msv%valr
         if (count(mpl%msv%isnotr(Dh))>0) Dhbar = sum(Dh,mask=mpl%msv%isnotr(Dh))/real(count(mpl%msv%isnotr(Dh)),kind_real)

         ! Approximate homogeneous vertical length-scale
         Dv = mpl%msv%valr
         Dv_valid = .false.
         jc3 = 1
         do jl0r=1,bpar%nl0r(ib)
            if (dmask(jc3,jl0r)) then
               distsq = dz(jc3,jl0r)**2
               if (inf(abs(lct_blk%raw(jc3,jl0r,ic1a,il0)),1.0_kind_real).and.(distsq>0.0)) then
                  if (lct_blk%raw(jc3,jl0r,ic1a,il0)>0.0) Dv(jl0r) = -distsq/(2.0*log(lct_blk%raw(jc3,jl0r,ic1a,il0)))
                  Dv_valid(jl0r) = .true.
               end if
            end if
         end do
         Dvbar = mpl%msv%valr
         if (bpar%nl0r(ib)>1) then
            if (count(Dv_valid)>0) then
               if (count(mpl%msv%isnotr(Dv))>0) then
                  Dvbar = sum(Dv,mask=mpl%msv%isnotr(Dv))/real(count(mpl%msv%isnotr(Dv)),kind_real)
               else
                  Dvbar = 1.0/rth
               end if
            end if
         else
             Dvbar = 0.0
         end if

         if (mpl%msv%isnotr(Dhbar).and.(mpl%msv%isnotr(Dvbar))) then
            ! Define norm and bounds
            do iscales=1,lct_blk%nscales
               minim%guess((iscales-1)*4+1:(iscales-1)*4+3) = (/Dhbar,Dhbar,Dvbar/)*Dscale**(iscales-1)
               if (lct_blk%nscales==1) then
                  minim%binf((iscales-1)*4+1:(iscales-1)*4+3) = (/1.0/Dscale,1.0/Dscale,1.0/Dscale/)*minim%guess(1:3)
                  minim%bsup((iscales-1)*4+1:(iscales-1)*4+3) = (/Dscale,Dscale,Dscale/)*minim%guess(1:3)
               else
                  minim%binf((iscales-1)*4+1:(iscales-1)*4+3) = (/1.0/sqrt(Dscale),1.0/sqrt(Dscale),1.0/sqrt(Dscale)/) &
                                                              & *minim%guess(1:3)*Dscale**(iscales-1)
                  minim%bsup((iscales-1)*4+1:(iscales-1)*4+3) = (/sqrt(Dscale),sqrt(Dscale),sqrt(Dscale)/)*minim%guess(1:3) &
                                                              & *Dscale**(iscales-1)
               end if
               minim%guess((iscales-1)*4+4) = 0.0
               minim%binf((iscales-1)*4+4) = -1.0
               minim%bsup((iscales-1)*4+4) = 1.0
            end do
            do iscales=1,lct_blk%nscales-1
               minim%guess(lct_blk%nscales*4+1) = 1.0/real(lct_blk%nscales,kind_real)
               minim%binf(lct_blk%nscales*4+1) = 0.1
               minim%bsup(lct_blk%nscales*4+1) = 1.0
            end do

            ! Fill minim
            minim%obs = pack(lct_blk%raw(:,:,ic1a,il0),.true.)
            minim%cost_function = 'fit_lct'
            minim%algo = trim(nam%minim_algo)
            minim%nc3 = nam%nc3
            minim%nl0 = bpar%nl0r(ib)
            minim%dx = dx
            minim%dy = dy
            minim%dz = dz
            minim%dmask = dmask
            minim%nscales = lct_blk%nscales

            ! Compute fit
            call minim%compute(mpl,lprt)

            ! Copy parameters
            do iscales=1,lct_blk%nscales
               do icomp=1,4
                   lct_blk%D(icomp,iscales,ic1a,il0) = minim%x((iscales-1)*4+icomp)
               end do
            end do
            if (lct_blk%nscales>1) then
               lct_blk%coef(1:lct_blk%nscales-1,ic1a,il0) = minim%x(lct_blk%nscales*4+1:lct_blk%nscales*4+lct_blk%nscales-1)
               lct_blk%coef(lct_blk%nscales,ic1a,il0) = 1.0-sum(lct_blk%coef(1:lct_blk%nscales-1,ic1a,il0))
            else
               lct_blk%coef(1,ic1a,il0) = 1.0
            end if

            ! Set vertical value at zero for the 2D case
            if (bpar%nl0r(ib)==1) lct_blk%D(3,:,ic1a,il0) = 0.0

            ! Check tensor validity
            valid = .true.
            do iscales=1,lct_blk%nscales
               if (valid) then
                  call check_cond(lct_blk%D(1,iscales,ic1a,il0),lct_blk%D(2,iscales,ic1a,il0),lct_blk%D(4,iscales,ic1a,il0),valid)
                  if (bpar%nl0r(ib)>1) valid = valid.and.(lct_blk%D(3,iscales,ic1a,il0)>0.0)
                  valid = valid.and.(lct_blk%coef(iscales,ic1a,il0)>0.0)
                  if (lct_blk%nscales>1) valid = valid.and.(inf(lct_blk%coef(iscales,ic1a,il0),1.0_kind_real))
               end if
            end do
            if (valid) then
               do iscales=1,lct_blk%nscales
                  if (nam%lct_diag(iscales)) then
                     ! Rescale diagonal tensor
                     diag_rescale = sqrt(1.0-lct_blk%D(4,iscales,ic1a,il0)**2)
                     lct_blk%D(1,iscales,ic1a,il0) = lct_blk%D(1,iscales,ic1a,il0)*diag_rescale
                     lct_blk%D(2,iscales,ic1a,il0) = lct_blk%D(2,iscales,ic1a,il0)*diag_rescale
                     lct_blk%D(4,iscales,ic1a,il0) = 0.0
                  end if
               end do

               ! Rebuild fit
               call fit_lct(mpl,nam%nc3,bpar%nl0r(ib),dx,dy,dz,dmask,lct_blk%nscales, &
             & lct_blk%D(:,:,ic1a,il0),lct_blk%coef(:,ic1a,il0),lct_blk%fit(:,:,ic1a,il0))
            else
               ! Missing values
               lct_blk%D(:,:,ic1a,il0) = mpl%msv%valr
               lct_blk%coef(:,ic1a,il0) = mpl%msv%valr
               lct_blk%fit(:,:,ic1a,il0) = mpl%msv%valr
            end if
         else
            ! Missing values
            lct_blk%D(:,:,ic1a,il0) = mpl%msv%valr
            lct_blk%coef(:,ic1a,il0) = mpl%msv%valr
            lct_blk%fit(:,:,ic1a,il0) = mpl%msv%valr
         end if
      else
         ! Missing values
         lct_blk%D(:,:,ic1a,il0) = mpl%msv%valr
         lct_blk%coef(:,ic1a,il0) = mpl%msv%valr
         lct_blk%fit(:,:,ic1a,il0) = mpl%msv%valr
      end if

      ! Update
      call mpl%prog_print(ic1a)
   end do
   call mpl%prog_final
end do

! Release memory
deallocate(Dh)
deallocate(Dv)
deallocate(dx)
deallocate(dy)
deallocate(dz)
deallocate(dmask)
deallocate(minim%x)
deallocate(minim%guess)
deallocate(minim%binf)
deallocate(minim%bsup)
deallocate(minim%obs)
deallocate(minim%dx)
deallocate(minim%dy)
deallocate(minim%dz)
deallocate(minim%dmask)

! End associate
end associate

end subroutine lct_blk_fitting

end module type_lct_blk
