!----------------------------------------------------------------------
! Module: type_lct_blk
!> Purpose: LCT data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_lct_blk

!$ use omp_lib
use tools_display, only: prog_init,prog_print
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
   integer :: ib                                !< Block index
   integer :: nscales                           !< Number of LCT scales
   integer,allocatable :: ncomp(:)              !< Number of LCT components

   ! LCT data
   real(kind_real),allocatable :: H(:,:,:)      !< LCT components
   real(kind_real),allocatable :: coef(:,:,:)   !< LCT coefficients
   real(kind_real),allocatable :: raw(:,:,:,:)  !< Raw correlations
   real(kind_real),allocatable :: norm(:,:,:,:) !< Norm to take nsub into account
   real(kind_real),allocatable :: fit(:,:,:,:)  !< Fitted correlations

   ! DT data
   real(kind_real),allocatable :: Dcoef(:,:,:)
   real(kind_real),allocatable :: D11(:,:,:)
   real(kind_real),allocatable :: D22(:,:,:)
   real(kind_real),allocatable :: D33(:,:,:)
   real(kind_real),allocatable :: D12(:,:,:)
contains
   procedure :: alloc => lct_blk_alloc
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
! Subroutine: lct_blk_alloc
!> Purpose: lct_blk object block allocation
!----------------------------------------------------------------------
subroutine lct_blk_alloc(lct_blk,nam,geom,bpar,hdata,ib)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk !< LCT block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
type(hdata_type),intent(in) :: hdata         !< HDIAG data
integer,intent(in) :: ib                     !< Block index

! Local variables
integer :: iscales

! Attributes
lct_blk%ib = ib
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
allocate(lct_blk%H(sum(lct_blk%ncomp),hdata%nc1a,geom%nl0))
allocate(lct_blk%coef(lct_blk%nscales,hdata%nc1a,geom%nl0))
allocate(lct_blk%raw(nam%nc3,bpar%nl0r(ib),hdata%nc1a,geom%nl0))
allocate(lct_blk%norm(nam%nc3,bpar%nl0r(ib),hdata%nc1a,geom%nl0))
allocate(lct_blk%fit(nam%nc3,bpar%nl0r(ib),hdata%nc1a,geom%nl0))
allocate(lct_blk%Dcoef(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D11(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D22(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D33(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D12(geom%nc0a,geom%nl0,lct_blk%nscales))

! Initialization
call msr(lct_blk%H)
call msr(lct_blk%coef)
call msr(lct_blk%raw)
call msr(lct_blk%norm)
call msr(lct_blk%fit)
call msr(lct_blk%Dcoef)
call msr(lct_blk%D11)
call msr(lct_blk%D22)
call msr(lct_blk%D33)
call msr(lct_blk%D12)

end subroutine lct_blk_alloc

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
if (allocated(lct_blk%Dcoef)) deallocate(lct_blk%Dcoef)
if (allocated(lct_blk%D11)) deallocate(lct_blk%D11)
if (allocated(lct_blk%D22)) deallocate(lct_blk%D22)
if (allocated(lct_blk%D33)) deallocate(lct_blk%D33)
if (allocated(lct_blk%D12)) deallocate(lct_blk%D12)

end subroutine lct_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_blk_correlation
!> Purpose: compute raw correlation
!----------------------------------------------------------------------
subroutine lct_blk_correlation(lct_blk,nam,geom,bpar,hdata,mom_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk !< LCT block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
type(hdata_type),intent(in) :: hdata         !< HDIAG data
type(mom_blk_type),intent(in) :: mom_blk     !< Moments block

! Local variables
integer :: jsub,il0,jl0r,jc3,ic1a,ic1
real(kind_real) :: den

! Associate
associate(ib=>lct_blk%ib)

! Initialize
lct_blk%raw = 0.0
lct_blk%norm = 0.0

! Sum over jsub
do jsub=1,mom_blk%nsub
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,nam%nc3
            do ic1a=1,hdata%nc1a
               ic1 = hdata%c1a_to_c1(ic1a)
               if (hdata%c1l0_log(ic1,il0)) then
                  den = mom_blk%m2_1(ic1a,jc3,jl0r,il0,jsub)*mom_blk%m2_2(ic1a,jc3,jl0r,il0,jsub)
                  if (den>0.0) then
                     lct_blk%raw(jc3,jl0r,ic1a,il0) = lct_blk%raw(jc3,jl0r,ic1a,il0)+mom_blk%m11(ic1a,jc3,jl0r,il0,jsub)/sqrt(den)
                     lct_blk%norm(jc3,jl0r,ic1a,il0) = lct_blk%norm(jc3,jl0r,ic1a,il0)+1.0
                  end if
               end if
            end do
         end do
      end do
   end do
end do

! Normalize
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,nam%nc3
         do ic1a=1,hdata%nc1a
            ic1 = hdata%c1a_to_c1(ic1a)
            if (hdata%c1l0_log(ic1,il0)) then
               if (lct_blk%norm(jc3,jl0r,ic1a,il0)>0.0) lct_blk%raw(jc3,jl0r,ic1a,il0) = lct_blk%raw(jc3,jl0r,ic1a,il0) &
             & /lct_blk%norm(jc3,jl0r,ic1a,il0)
            end if
         end do
      end do
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
integer :: il0,jl0r,jl0,ic1a,ic1,jc3,iscales,offset,progint
real(kind_real) :: distsq,Hhbar,Hvbar
real(kind_real),allocatable :: Hh(:),Hv(:),dx(:,:),dy(:,:),dz(:)
logical :: spd
logical,allocatable :: dmask(:,:),done(:)
type(minim_type) :: minim

! Associate
associate(ib=>lct_blk%ib)

! Allocation
allocate(Hh(nam%nc3))
allocate(Hv(bpar%nl0r(ib)))
allocate(dx(nam%nc3,bpar%nl0r(ib)))
allocate(dy(nam%nc3,bpar%nl0r(ib)))
allocate(dz(bpar%nl0r(ib)))
allocate(dmask(nam%nc3,bpar%nl0r(ib)))
allocate(done(hdata%nc1a))
minim%nx = sum(lct_blk%ncomp)
if (lct_blk%nscales>1) minim%nx = minim%nx+lct_blk%nscales-1
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

do il0=1,geom%nl0
   write(mpl%unit,'(a13,a,i3,a)',advance='no') '','Level ',nam%levs(il0),':'
   call flush(mpl%unit)

   ! Initialization
   call prog_init(progint,done)

   do ic1a=1,hdata%nc1a
      ! Global index
      ic1 = hdata%c1a_to_c1(ic1a)

      if (hdata%c1l0_log(ic1,il0)) then
         ! Prepare vectors
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               dmask(jc3,jl0r) = hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)
               if (dmask(jc3,jl0r)) then
                  dx(jc3,jl0r) = geom%lon(hdata%c1c3_to_c0(ic1,jc3))-geom%lon(hdata%c1c3_to_c0(ic1,1))
                  dy(jc3,jl0r) = geom%lat(hdata%c1c3_to_c0(ic1,jc3))-geom%lat(hdata%c1c3_to_c0(ic1,1))
                  call lonlatmod(dx(jc3,jl0r),dy(jc3,jl0r))
                  dx(jc3,jl0r) = dx(jc3,jl0r)*cos(geom%lat(hdata%c1c3_to_c0(ic1,1)))
               end if
            end do
            dz(jl0r) = real(nam%levs(jl0)-nam%levs(il0),kind_real)
         end do

         ! Approximate homogeneous horizontal length-scale
         call msr(Hh)
         do jl0r=1,bpar%nl0r(ib)
            if (.not.(abs(dz(jl0r))>0.0)) then
               do jc3=1,nam%nc3
                  if (dmask(jc3,jl0r)) then
                     distsq = dx(jc3,jl0r)**2+dy(jc3,jl0r)**2
                     if ((lct_blk%raw(jc3,jl0r,ic1a,il0)>0.0).and.(distsq>0.0)) Hh(jc3) = -2.0*log(lct_blk%raw(jc3,jl0r,ic1a,il0)) &
                   & /distsq
                  end if
               end do
            end if
         end do
         if (count(isnotmsr(Hh))>0) then
            Hhbar = sum(Hh,mask=isnotmsr(Hh))/real(count(isnotmsr(Hh)),kind_real)
         else
            return
         end if
         if (lct_blk%nscales>1) Hhbar = Hhbar*Hscale

         ! Approximate homogeneous vertical length-scale
         call msr(Hv)
         jc3 = 1
         do jl0r=1,bpar%nl0r(ib)
            distsq = dz(jl0r)**2
            if ((lct_blk%raw(jc3,jl0r,ic1a,il0)>0.0).and.(distsq>0.0)) Hv(jl0r) = -2.0*log(lct_blk%raw(jc3,jl0r,ic1a,il0))/distsq
         end do
         if (bpar%nl0r(ib)>0) then
            Hvbar = 1.0
         else
            if (count(isnotmsr(Hv))>0) then
               Hvbar = sum(Hv,mask=isnotmsr(Hv))/real(count(isnotmsr(Hv)),kind_real)
            else
              return
            end if
         end if
         if (lct_blk%nscales>1) Hvbar = Hvbar*Hscale

         ! Define norm and bounds
         offset = 0
         do iscales=1,lct_blk%nscales
            minim%guess(offset+1:offset+3) = (/Hhbar,Hhbar,Hvbar/)/Hscale**(iscales-1)
            minim%norm(offset+1:offset+3) = (/Hhbar,Hhbar,Hvbar/)/Hscale**(iscales-1)
            minim%binf(offset+1:offset+3) = (/1.0/sqrt(Hscale),1.0/sqrt(Hscale),1.0/sqrt(Hscale)/)*minim%guess(1:3) &
                                          & /Hscale**(iscales-1)
            minim%bsup(offset+1:offset+3) = (/sqrt(Hscale),sqrt(Hscale),sqrt(Hscale)/)*minim%guess(1:3)/Hscale**(iscales-1)
            offset = offset+3
            if (lct_blk%ncomp(iscales)==4) then
               minim%guess(offset+1) = 0.0
               minim%norm(offset+1) = 1.0
               minim%binf(offset+1) = -1.0
               minim%bsup(offset+1) = 1.0
               offset = offset+1
            end if
            if (lct_blk%nscales>1) then
               minim%guess(offset+1) = 1.0/real(lct_blk%nscales,kind_real)
               minim%norm(offset+1) = 1.0/real(lct_blk%nscales,kind_real)
               minim%binf(offset+1) = 0.0
               minim%bsup(offset+1) = 1.0
               offset = offset+1
            end if
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
         minim%ncomp = lct_blk%ncomp

         ! Compute fit
         minim%cost_function = 'fit_lct'
         call minim%compute(lprt)

         ! Copy parameters
         lct_blk%H(:,ic1a,il0) = minim%x(1:sum(lct_blk%ncomp))
         if (lct_blk%nscales>1) then
            lct_blk%coef(1:lct_blk%nscales-1,ic1a,il0) = minim%x(sum(lct_blk%ncomp)+1:sum(lct_blk%ncomp)+lct_blk%nscales-1)
            lct_blk%coef(lct_blk%nscales,ic1a,il0) = 1.0-sum(lct_blk%coef(1:lct_blk%nscales-1,ic1a,il0))
         else
            lct_blk%coef(1,ic1a,il0) = 1.0
         end if

         ! Fixed positive value for the 2D case
         if (bpar%nl0r(ib)==1) then
            offset = 0
            do iscales=1,lct_blk%nscales
               lct_blk%H(offset+3,ic1a,il0) = 1.0
               offset = offset+lct_blk%ncomp(iscales)
            end do
         end if

         ! Check positive-definiteness
         spd = .true.
         offset = 0
         do iscales=1,lct_blk%nscales
            spd = spd.and.(lct_blk%H(offset+1,ic1a,il0)>0.0).and.(lct_blk%H(offset+2,ic1a,il0)>0.0) &
                & .and.(lct_blk%H(offset+3,ic1a,il0)>0.0)
            if (lct_blk%ncomp(iscales)==4) spd = spd.and.(lct_blk%H(offset+4,ic1a,il0)>-1.0).and.(lct_blk%H(offset+4,ic1a,il0)<1.0)
            if (iscales<lct_blk%nscales) spd = spd.and.(lct_blk%coef(iscales,ic1a,il0)>0.0)
            offset = offset+lct_blk%ncomp(iscales)
         end do
         if (spd) then
            ! Rebuild fit
            call fit_lct(nam%nc3,bpar%nl0r(ib),dx,dy,dz,dmask,lct_blk%nscales,lct_blk%ncomp,lct_blk%H(:,ic1a,il0), &
          & lct_blk%coef(:,ic1a,il0),lct_blk%fit(:,:,ic1a,il0))
         else
            ! Missing values
            call msr(lct_blk%H(:,ic1a,il0))
            call msr(lct_blk%coef(:,ic1a,il0))
            call msr(lct_blk%fit(:,:,ic1a,il0))
         end if
      end if

      ! Update
      done(ic1a) = .true.
      call prog_print(progint,done)
   end do
   write(mpl%unit,'(a)') '100%'
   call flush(mpl%unit)
end do

! End associate
end associate

end subroutine lct_blk_fitting

end module type_lct_blk
