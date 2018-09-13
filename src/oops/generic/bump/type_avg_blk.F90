!----------------------------------------------------------------------
! Module: type_avg_blk
!> Purpose: averaged statistics block derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_avg_blk

use tools_kinds, only: kind_real
use tools_missing, only: msr,isanynotmsr,isnotmsr,ismsr
use tools_repro, only: sup,inf
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
implicit none

real(kind_real),parameter :: var_min = 1.0e-24_kind_real !< Minimum variance for correlation computation
real(kind_real),parameter :: nc1a_cor_th = 0.5           !< Threshold on the effective sampling size for asymptotic statistics

! Averaged statistics block derived type
type avg_blk_type
   integer :: ic2                                        !< Global index
   integer :: ib                                         !< Block index
   integer :: ne                                         !< Ensemble size
   integer :: nsub                                       !< Sub-ensembles number
   real(kind_real),allocatable :: m2(:,:)                !< Variance
   real(kind_real),allocatable :: m4(:,:)                !< Fourth-order centered moment
   real(kind_real),allocatable :: m2flt(:)               !< Filtered variance
   real(kind_real),allocatable :: nc1a(:,:,:)            !< Number of points in subset Sc1 on halo A
   real(kind_real),allocatable :: m11(:,:,:)             !< Covariance average
   real(kind_real),allocatable :: m11m11(:,:,:,:,:)      !< Product of covariances average
   real(kind_real),allocatable :: m2m2(:,:,:,:,:)        !< Product of variances average
   real(kind_real),allocatable :: m22(:,:,:,:)           !< Fourth-order centered moment average
   real(kind_real),allocatable :: nc1a_cor(:,:,:)        !< Number of points in subset Sc1 on halo A with valid correlations
   real(kind_real),allocatable :: cor(:,:,:)             !< Correlation average
   real(kind_real),allocatable :: m11asysq(:,:,:)        !< Squared asymptotic covariance average
   real(kind_real),allocatable :: m2m2asy(:,:,:)         !< Product of asymptotic variances average
   real(kind_real),allocatable :: m22asy(:,:,:)          !< Asymptotic fourth-order centered moment average
   real(kind_real),allocatable :: m11sq(:,:,:)           !< Squared covariance average for several ensemble sizes
   real(kind_real),allocatable :: m11sta(:,:,:)          !< Ensemble covariance/static covariance product
   real(kind_real),allocatable :: stasq(:,:,:)           !< Squared static covariance
   real(kind_real),allocatable :: m11lrm11sub(:,:,:,:,:) !< LR covariance/HR covariance product average
   real(kind_real),allocatable :: m11lrm11(:,:,:)        !< LR covariance/HR covariance product average, averaged over sub-ensembles
   real(kind_real),allocatable :: m11lrm11asy(:,:,:)     !< LR covariance/HR asymptotic covariance product average
contains
   procedure :: alloc => avg_blk_alloc
   procedure :: dealloc => avg_blk_dealloc
   procedure :: copy => avg_blk_copy
   procedure :: compute => avg_blk_compute
   procedure :: compute_asy => avg_blk_compute_asy
   procedure :: compute_lr => avg_blk_compute_lr
   procedure :: compute_asy_lr => avg_blk_compute_asy_lr
end type avg_blk_type

private
public :: avg_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: avg_blk_alloc
!> Purpose: averaged statistics block data allocation
!----------------------------------------------------------------------
subroutine avg_blk_alloc(avg_blk,nam,geom,bpar,ic2,ib,ne,nsub)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk !< Averaged statistics block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
integer,intent(in) :: ne                     !< Ensemble size
integer,intent(in) :: nsub                   !< Sub-ensembles number
integer,intent(in) :: ic2                    !< Global index
integer,intent(in) :: ib                     !< Block index

! Set attributes
avg_blk%ic2 = ic2
avg_blk%ib = ib
avg_blk%ne = ne
avg_blk%nsub = nsub

! Allocation
if (.not.allocated(avg_blk%nc1a)) then
   if ((ic2==0).or.(nam%var_diag)) then
      allocate(avg_blk%m2(geom%nl0,avg_blk%nsub))
      if (nam%var_filter) then
         if (.not.nam%gau_approx) allocate(avg_blk%m4(geom%nl0,avg_blk%nsub))
         allocate(avg_blk%m2flt(geom%nl0))
      end if
   end if
   if ((ic2==0).or.(nam%local_diag)) then
      allocate(avg_blk%nc1a(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11m11(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk%nsub,avg_blk%nsub))
      allocate(avg_blk%m2m2(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk%nsub,avg_blk%nsub))
      if (.not.nam%gau_approx) allocate(avg_blk%m22(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk%nsub))
      allocate(avg_blk%nc1a_cor(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%cor(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11asysq(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m2m2asy(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      if (.not.nam%gau_approx) allocate(avg_blk%m22asy(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11sq(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd')
         allocate(avg_blk%m11sta(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
         allocate(avg_blk%stasq(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      case ('dual-ens')
         allocate(avg_blk%m11lrm11sub(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk%nsub,avg_blk%nsub))
         allocate(avg_blk%m11lrm11(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
         allocate(avg_blk%m11lrm11asy(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      end select
   end if
end if

! Initialization
if ((ic2==0).or.(nam%var_diag)) then
   call msr(avg_blk%m2)
   if (nam%var_filter) then
      if (.not.nam%gau_approx) call msr(avg_blk%m4)
      call msr(avg_blk%m2flt)
   end if
end if
if ((ic2==0).or.(nam%local_diag)) then
   call msr(avg_blk%nc1a)
   call msr(avg_blk%m11)
   call msr(avg_blk%m11m11)
   call msr(avg_blk%m2m2)
   if (.not.nam%gau_approx) call msr(avg_blk%m22)
   call msr(avg_blk%nc1a_cor)
   call msr(avg_blk%cor)
   call msr(avg_blk%m11asysq)
   call msr(avg_blk%m2m2asy)
   if (.not.nam%gau_approx) call msr(avg_blk%m22asy)
   call msr(avg_blk%m11sq)
   select case (trim(nam%method))
   case ('hyb-avg_blk','hyb-rnd')
      call msr(avg_blk%m11sta)
      call msr(avg_blk%stasq)
   case ('dual-ens')
      call msr(avg_blk%m11lrm11sub)
      call msr(avg_blk%m11lrm11)
      call msr(avg_blk%m11lrm11asy)
   end select
end if

end subroutine avg_blk_alloc

!----------------------------------------------------------------------
! Subroutine: avg_blk_dealloc
!> Purpose: averaged statistics block data deallocation
!----------------------------------------------------------------------
subroutine avg_blk_dealloc(avg_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk !< Averaged statistics block

! Allocation
if (allocated(avg_blk%m2)) deallocate(avg_blk%m2)
if (allocated(avg_blk%m4)) deallocate(avg_blk%m4)
if (allocated(avg_blk%m2flt)) deallocate(avg_blk%m2flt)
if (allocated(avg_blk%nc1a)) deallocate(avg_blk%nc1a)
if (allocated(avg_blk%m11)) deallocate(avg_blk%m11)
if (allocated(avg_blk%m11m11)) deallocate(avg_blk%m11m11)
if (allocated(avg_blk%m2m2)) deallocate(avg_blk%m2m2)
if (allocated(avg_blk%nc1a_cor)) deallocate(avg_blk%nc1a_cor)
if (allocated(avg_blk%cor)) deallocate(avg_blk%cor)
if (allocated(avg_blk%m11asysq)) deallocate(avg_blk%m11asysq)
if (allocated(avg_blk%m2m2asy)) deallocate(avg_blk%m2m2asy)
if (allocated(avg_blk%m11sq)) deallocate(avg_blk%m11sq)
if (allocated(avg_blk%m11sta)) deallocate(avg_blk%m11sta)
if (allocated(avg_blk%stasq)) deallocate(avg_blk%stasq)
if (allocated(avg_blk%m11lrm11sub)) deallocate(avg_blk%m11lrm11sub)
if (allocated(avg_blk%m11lrm11)) deallocate(avg_blk%m11lrm11)
if (allocated(avg_blk%m11lrm11asy)) deallocate(avg_blk%m11lrm11asy)

end subroutine avg_blk_dealloc

!----------------------------------------------------------------------
! Function: avg_blk_copy
!> Purpose: averaged statistics block copy
!----------------------------------------------------------------------
type(avg_blk_type) function avg_blk_copy(avg_blk,nam,geom,bpar)

implicit none

! Passed variables
class(avg_blk_type),intent(in) :: avg_blk !< Averaged statistics block
type(nam_type),intent(in) :: nam          !< Namelist
type(geom_type),intent(in) :: geom        !< Geometry
type(bpar_type),intent(in) :: bpar        !< Block parameters

! Deallocation
call avg_blk_copy%dealloc

! Allocation
call avg_blk_copy%alloc(nam,geom,bpar,avg_blk%ic2,avg_blk%ib,avg_blk%ne,avg_blk%nsub)

! Copy data
if ((avg_blk%ic2==0).or.(nam%var_diag)) then
   avg_blk_copy%m2 = avg_blk%m2
   if (nam%var_filter) then
      if (.not.nam%gau_approx) avg_blk_copy%m4 = avg_blk%m4
      avg_blk_copy%m2flt = avg_blk%m2flt
   end if
end if
if ((avg_blk%ic2==0).or.(nam%local_diag)) then
   avg_blk_copy%nc1a = avg_blk%nc1a
   avg_blk_copy%m11 = avg_blk%m11
   avg_blk_copy%m11m11 = avg_blk%m11m11
   avg_blk_copy%m2m2 = avg_blk%m2m2
   if (.not.nam%gau_approx) avg_blk_copy%m22 = avg_blk%m22
   avg_blk_copy%nc1a_cor = avg_blk%nc1a_cor
   avg_blk_copy%cor = avg_blk%cor
   avg_blk_copy%m11asysq = avg_blk%m11asysq
   avg_blk_copy%m2m2asy = avg_blk%m2m2asy
   if (.not.nam%gau_approx) avg_blk_copy%m22asy = avg_blk%m22asy
   avg_blk_copy%m11sq = avg_blk%m11sq
   select case (trim(nam%method))
   case ('hyb-avg_blk','hyb-rnd')
      avg_blk_copy%m11sta = avg_blk%m11sta
      avg_blk_copy%stasq = avg_blk%stasq
   case ('dual-ens')
      avg_blk_copy%m11lrm11sub = avg_blk%m11lrm11sub
      avg_blk_copy%m11lrm11sub = avg_blk%m11lrm11sub
      avg_blk_copy%m11lrm11asy = avg_blk%m11lrm11asy
   end select
end if

end function avg_blk_copy

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute
!> Purpose: compute averaged statistics via spatial-angular erogodicity assumption
!----------------------------------------------------------------------
subroutine avg_blk_compute(avg_blk,nam,geom,bpar,hdata,mom_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk !< Averaged statistics block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
type(hdata_type),intent(in) :: hdata         !< HDIAG data
type(mom_blk_type),intent(in) :: mom_blk     !< Moments

! Local variables
integer :: il0,jl0,jl0r,jc3,isub,jsub,ic1a,ic1,nc1amax,nc1a
real(kind_real) :: m2_1,m2_2
real(kind_real),allocatable :: list_m11(:),list_m11m11(:,:,:),list_m2m2(:,:,:),list_m22(:,:),list_cor(:)
logical :: involved,valid

! Associate
associate(ic2=>avg_blk%ic2,ib=>avg_blk%ib)

if ((ic2==0).or.(nam%var_diag)) then
   ! Copy variance
   if (ic2==0) then
      do isub=1,avg_blk%nsub
         do il0=1,geom%nl0
            jl0r = bpar%il0rz(il0,ib)
            avg_blk%m2(il0,isub) = sum(mom_blk%m2_1(:,1,il0,isub))/real(nam%nc1,kind_real)
            if (nam%var_filter.and.(.not.nam%gau_approx)) avg_blk%m4(il0,isub) = sum(mom_blk%m22(:,1,jl0r,il0,isub)) &
                                                                               & /real(nam%nc1,kind_real)
         end do
      end do
   else
      avg_blk%m2 = 0.0
      if (nam%var_filter.and.(.not.nam%gau_approx)) avg_blk%m4 = 0.0
      do ic1a=1,hdata%nc1a
         ic1 = hdata%c1a_to_c1(ic1a)
         if (ic1==hdata%c2_to_c1(ic2)) then
            do isub=1,avg_blk%nsub
               do il0=1,geom%nl0
                  jl0r = bpar%il0rz(il0,ib)
                  avg_blk%m2(il0,isub) = mom_blk%m2_1(ic1a,1,il0,isub)
                  if (nam%var_filter.and.(.not.nam%gau_approx)) avg_blk%m4(il0,isub) = mom_blk%m22(ic1a,1,jl0r,il0,isub)
               end do
            end do
         end if
      end do
   end if
end if

if ((ic2==0).or.(nam%local_diag)) then
   ! Check whether this task is involved
   if (ic2>0) then
      involved = any(hdata%local_mask(hdata%c1a_to_c1,ic2))
   else
      involved = .true.
   end if

   if (involved) then
      ! Average
      !$omp parallel do schedule(static) private(il0,jl0r,jl0,nc1amax,jc3,nc1a,ic1a,ic1,valid,m2_1,m2_2,isub,jsub), &
      !$omp&                             firstprivate(list_m11,list_m11m11,list_m2m2,list_m22,list_cor)
      do il0=1,geom%nl0
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)

            ! Allocation
            if (ic2>0) then
               nc1amax = count(hdata%local_mask(hdata%c1a_to_c1,ic2))
            else
               nc1amax = hdata%nc1a
            end if
            allocate(list_m11(nc1amax))
            allocate(list_m11m11(nc1amax,avg_blk%nsub,avg_blk%nsub))
            allocate(list_m2m2(nc1amax,avg_blk%nsub,avg_blk%nsub))
            allocate(list_m22(nc1amax,avg_blk%nsub))
            allocate(list_cor(nc1amax))

            do jc3=1,bpar%nc3(ib)
               ! Fill lists
               nc1a = 0
               do ic1a=1,hdata%nc1a
                  ! Index
                  ic1 = hdata%c1a_to_c1(ic1a)

                  ! Check validity
                  valid = hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)
                  if (ic2>0) valid = valid.and.hdata%local_mask(ic1,ic2)

                  if (valid) then
                     ! Update
                     nc1a = nc1a+1

                     ! Averages for diagnostics
                     list_m11(nc1a) = sum(mom_blk%m11(ic1a,jc3,jl0r,il0,:))/real(avg_blk%nsub,kind_real)
                     do isub=1,avg_blk%nsub
                        do jsub=1,avg_blk%nsub
                           list_m11m11(nc1a,jsub,isub) = mom_blk%m11(ic1a,jc3,jl0r,il0,isub)*mom_blk%m11(ic1a,jc3,jl0r,il0,jsub)
                           list_m2m2(nc1a,jsub,isub) = mom_blk%m2_1(ic1a,jc3,il0,isub)*mom_blk%m2_2(ic1a,jc3,jl0,jsub)
                        end do
                        if (.not.nam%gau_approx) list_m22(nc1a,isub) = mom_blk%m22(ic1a,jc3,jl0r,il0,isub)
                     end do

                     ! Correlation
                     m2_1 = sum(mom_blk%m2_1(ic1a,jc3,il0,:))/real(avg_blk%nsub,kind_real)
                     m2_2 = sum(mom_blk%m2_2(ic1a,jc3,jl0,:))/real(avg_blk%nsub,kind_real)
                     if (sup(m2_1,var_min).and.sup(m2_2,var_min)) then
                        list_cor(nc1a) = list_m11(nc1a)/sqrt(m2_1*m2_2)
                        if (sup(abs(list_cor(nc1a)),1.0_kind_real)) call msr(list_cor(nc1a))
                     else
                        call msr(list_cor(nc1a))
                     end if
                  end if
               end do

               ! Average
               avg_blk%nc1a(jc3,jl0r,il0) = real(nc1a,kind_real)
               if (nc1a>0) then
                  avg_blk%m11(jc3,jl0r,il0) = sum(list_m11(1:nc1a))
                  do isub=1,avg_blk%nsub
                     do jsub=1,avg_blk%nsub
                        avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) = sum(list_m11m11(1:nc1a,jsub,isub))
                        avg_blk%m2m2(jc3,jl0r,il0,jsub,isub) = sum(list_m2m2(1:nc1a,jsub,isub))
                     end do
                     if (.not.nam%gau_approx) avg_blk%m22(jc3,jl0r,il0,isub) = sum(list_m22(1:nc1a,isub))
                  end do
                  avg_blk%nc1a_cor(jc3,jl0r,il0) = real(count(isnotmsr(list_cor(1:nc1a))),kind_real)
                  if (avg_blk%nc1a_cor(jc3,jl0r,il0)>0.0) then
                     avg_blk%cor(jc3,jl0r,il0) = sum(list_cor(1:nc1a),mask=isnotmsr(list_cor(1:nc1a)))
                  else
                     call msr(avg_blk%cor(jc3,jl0r,il0))
                  end if
               else
                  avg_blk%m11(jc3,jl0r,il0) = 0.0
                  do isub=1,avg_blk%nsub
                     do jsub=1,avg_blk%nsub
                        avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) = 0.0
                        avg_blk%m2m2(jc3,jl0r,il0,jsub,isub) = 0.0
                     end do
                     if (.not.nam%gau_approx) avg_blk%m22(jc3,jl0r,il0,isub) = 0.0
                  end do
                  avg_blk%nc1a_cor(jc3,jl0r,il0) = 0.0
                  avg_blk%cor(jc3,jl0r,il0) = 0.0
               end if
            end do

            ! Release memory
            deallocate(list_m11)
            deallocate(list_m11m11)
            deallocate(list_m2m2)
            deallocate(list_m22)
            deallocate(list_cor)
         end do
      end do
      !$omp end parallel do
   else
      ! Set to zero
      avg_blk%nc1a = 0
      avg_blk%m11 = 0.0
      avg_blk%m11m11 = 0.0
      avg_blk%m2m2 = 0.0
      if (.not.nam%gau_approx) avg_blk%m22 = 0.0
      avg_blk%nc1a_cor = 0.0
      avg_blk%cor = 0.0
   end if
end if

! End associate
end associate

end subroutine avg_blk_compute

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_asy
!> Purpose: compute asymptotic statistics
!----------------------------------------------------------------------
subroutine avg_blk_compute_asy(avg_blk,nam,geom,bpar,ne)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk !< Averaged statistics block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
integer,intent(in) :: ne                     !< Ensemble size

! Local variables
integer :: il0,jl0r,jc3,isub,jsub,n
real(kind_real) :: P1,P3,P4,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17
real(kind_real),allocatable :: m11asysq(:,:),m2m2asy(:,:),m22asy(:)

! Associate
associate(ic2=>avg_blk%ic2,ib=>avg_blk%ib)

if ((ic2==0).or.(nam%local_diag)) then
   ! Ensemble size-dependent coefficients
   n = ne
   P1 = 1.0/real(n,kind_real)
   P3 = 1.0/real(n*(n-1),kind_real)
   P4 = 1.0/real(n-1,kind_real)
   P14 = real(n**2-2*n+2,kind_real)/real(n*(n-1),kind_real)
   P16 = real(n,kind_real)/real(n-1,kind_real)

   ! Ensemble/sub-ensemble size-dependent coefficients
   n = avg_blk%ne/avg_blk%nsub
   P7 = real((n-1)*(n**2-3*n+1),kind_real)/real(n*(n-2)*(n-3),kind_real)
   P8 = real(n-1,kind_real)/real(n*(n-2)*(n-3),kind_real)
   P9 = -real(n,kind_real)/real((n-2)*(n-3),kind_real)
   P10 = -real((n-1)*(2*n-3),kind_real)/real(n*(n-2)*(n-3),kind_real)
   P11 = real(n*(n**2-2*n+3),kind_real)/real((n-1)*(n-2)*(n-3),kind_real)
   P12 = real(n*(n-1),kind_real)/real((n-2)*(n+1),kind_real)
   P13 = -real(n-1,kind_real)/real((n-2)*(n+1),kind_real)
   P15 = real((n-1)**2,kind_real)/real(n*(n-3),kind_real)
   P17 = real((n-1)**2,kind_real)/real((n-2)*(n+1),kind_real)

   ! Asymptotic statistics
   !$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub) firstprivate(m11asysq,m2m2asy,m22asy)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,bpar%nc3(ib)
            if (avg_blk%nc1a_cor(jc3,jl0r,il0)>nc1a_cor_th*avg_blk%nc1a(jc3,jl0r,il0)) then
               ! Allocation
               allocate(m11asysq(avg_blk%nsub,avg_blk%nsub))
               allocate(m2m2asy(avg_blk%nsub,avg_blk%nsub))
               allocate(m22asy(avg_blk%nsub))

               ! Asymptotic statistics
               do isub=1,avg_blk%nsub
                  do jsub=1,avg_blk%nsub
                     if (isub==jsub) then
                        ! Diagonal terms
                        if (nam%gau_approx) then
                           ! Gaussian approximation
                           m11asysq(jsub,isub) = P17*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) &
                                               & +P13*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)
                           m2m2asy(jsub,isub) = 2.0*P13*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) &
                                              & +P12*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)
                        else
                           ! General case
                           m11asysq(jsub,isub) = P15*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) &
                                               & +P8*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)+P9*avg_blk%m22(jc3,jl0r,il0,jsub)
                           m2m2asy(jsub,isub) = 2.0*P8*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) &
                                              & +P7*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)+P9*avg_blk%m22(jc3,jl0r,il0,jsub)
                           m22asy(jsub) = P10*(2.0*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub)+avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)) &
                                        & +P11*avg_blk%m22(jc3,jl0r,il0,jsub)
                        end if
                     else
                        ! Off-diagonal terms
                        m11asysq(jsub,isub) = avg_blk%m11m11(jc3,jl0r,il0,jsub,isub)
                        m2m2asy(jsub,isub) = avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)
                     end if
                  end do
               end do

               ! Sum
               avg_blk%m11asysq(jc3,jl0r,il0) = sum(m11asysq)/real(avg_blk%nsub**2,kind_real)
               avg_blk%m2m2asy(jc3,jl0r,il0) = sum(m2m2asy)/real(avg_blk%nsub**2,kind_real)
               if (.not.nam%gau_approx) avg_blk%m22asy(jc3,jl0r,il0) = sum(m22asy)/real(avg_blk%nsub,kind_real)

               ! Check positivity
               if (avg_blk%m11asysq(jc3,jl0r,il0)<0.0) call msr(avg_blk%m11asysq(jc3,jl0r,il0))
               if (avg_blk%m2m2asy(jc3,jl0r,il0)<0.0) call msr(avg_blk%m2m2asy(jc3,jl0r,il0))
               if (.not.nam%gau_approx) then
                  if (avg_blk%m22asy(jc3,jl0r,il0)<0.0) call msr(avg_blk%m22asy(jc3,jl0r,il0))
               end if

               ! Squared covariance average
               if (nam%gau_approx) then
                  ! Gaussian approximation
                  if (isnotmsr(avg_blk%m11asysq(jc3,jl0r,il0)).and.isnotmsr(avg_blk%m2m2asy(jc3,jl0r,il0))) &
                & avg_blk%m11sq(jc3,jl0r,il0) = P16*avg_blk%m11asysq(jc3,jl0r,il0)+P4*avg_blk%m2m2asy(jc3,jl0r,il0)
               else
                  ! General case
                  if (isnotmsr(avg_blk%m22asy(jc3,jl0r,il0)).and.isnotmsr(avg_blk%m11asysq(jc3,jl0r,il0)) &
                & .and.isnotmsr(avg_blk%m2m2asy(jc3,jl0r,il0))) &
                & avg_blk%m11sq(jc3,jl0r,il0) = P1*avg_blk%m22asy(jc3,jl0r,il0)+P14*avg_blk%m11asysq(jc3,jl0r,il0) &
                                         & +P3*avg_blk%m2m2asy(jc3,jl0r,il0)
               end if

               ! Check value
               if (ismsr(avg_blk%m11sq(jc3,jl0r,il0))) then
                  if (inf(avg_blk%m11sq(jc3,jl0r,il0),avg_blk%m11asysq(jc3,jl0r,il0))) call msr(avg_blk%m11sq(jc3,jl0r,il0))
                  if (inf(avg_blk%m11sq(jc3,jl0r,il0),avg_blk%m11(jc3,jl0r,il0)**2)) call msr(avg_blk%m11sq(jc3,jl0r,il0))
               end if

               ! Allocation
               deallocate(m11asysq)
               deallocate(m2m2asy)
               deallocate(m22asy)
            end if
         end do
      end do
   end do
   !$omp end parallel do
end if

! End associate
end associate

end subroutine avg_blk_compute_asy

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_lr
!> Purpose: compute averaged statistics via spatial-angular erogodicity assumption, for LR covariance/HR covariance and LR covariance/HR asymptotic covariance products
!----------------------------------------------------------------------
subroutine avg_blk_compute_lr(avg_blk_lr,mpl,nam,geom,bpar,hdata,mom_blk,mom_lr_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk_lr !< Low-resolution averaged statistics block
type(mpl_type),intent(in) :: mpl                !< MPI data
type(nam_type),intent(in) :: nam                !< Namelist
type(geom_type),intent(in) :: geom              !< Geometry
type(bpar_type),intent(in) :: bpar              !< Block parameters
type(hdata_type),intent(in) :: hdata            !< HDIAG data
type(mom_blk_type),intent(in) :: mom_blk        !< Moments block
type(mom_blk_type),intent(in) :: mom_lr_blk     !< Low-resolution moments block

! Local variables
integer :: il0,jl0,jl0r,jc3,isub,jsub,ic1a,ic1,nc1amax,nc1a
real(kind_real),allocatable :: list_m11lrm11(:,:,:)
logical :: valid

! Associate
associate(ic2=>avg_blk_lr%ic2,ib=>avg_blk_lr%ib)

if ((ic2==0).or.(nam%local_diag)) then
   ! Check number of sub-ensembles
   if (avg_blk_lr%nsub/=avg_blk_lr%nsub) call mpl%abort('different number of sub-ensembles')

   ! Average
   !$omp parallel do schedule(static) private(il0,jl0r,jl0,nc1amax,list_m11lrm11,jc3,nc1a,ic1a,ic1,valid,isub,jsub)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)

         ! Allocation
         if (ic2>0) then
            nc1amax = count(hdata%local_mask(hdata%c1a_to_c1,ic2))
         else
            nc1amax = hdata%nc1a
         end if
         allocate(list_m11lrm11(nc1amax,avg_blk_lr%nsub,avg_blk_lr%nsub))

         do jc3=1,bpar%nc3(ib)
            ! Fill lists
            nc1a = 0
            do ic1a=1,hdata%nc1a
               ! Index
               ic1 = hdata%c1a_to_c1(ic1a)

               ! Check validity
               valid = hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)
               if (ic2>0) valid = valid.and.hdata%local_mask(ic1,ic2)

               if (valid) then
                  ! Update
                  nc1a = nc1a+1

                  ! Averages for diagnostics
                  do isub=1,avg_blk_lr%nsub
                     do jsub=1,avg_blk_lr%nsub
                        list_m11lrm11(nc1a,jsub,isub) = mom_blk%m11(ic1a,jc3,jl0r,il0,isub)*mom_lr_blk%m11(ic1a,jc3,jl0r,il0,jsub)
                     end do
                  end do
               end if
            end do

            ! Average
            avg_blk_lr%nc1a(jc3,jl0r,il0) = real(nc1a,kind_real)
            do isub=1,avg_blk_lr%nsub
               do jsub=1,avg_blk_lr%nsub
                  avg_blk_lr%m11lrm11sub(jc3,jl0r,il0,jsub,isub) = sum(list_m11lrm11(1:nc1a,jsub,isub))
               end do
            end do
         end do

         ! Release memory
         deallocate(list_m11lrm11)
      end do
   end do
   !$omp end parallel do
end if

! End associate
end associate

end subroutine avg_blk_compute_lr

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_asy_lr
!> Purpose: compute LR covariance/HR asymptotic covariance products
!----------------------------------------------------------------------
subroutine avg_blk_compute_asy_lr(avg_blk_lr,nam,geom,bpar)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk_lr !< Low-resolution averaged statistics block
type(nam_type),intent(in) :: nam                !< Namelist
type(geom_type),intent(in) :: geom              !< Geometry
type(bpar_type),intent(in) :: bpar              !< Block parameters

! Local variables
integer :: il0,jl0r,jc3,isub,jsub

! Associate
associate(ic2=>avg_blk_lr%ic2,ib=>avg_blk_lr%ib)

if ((ic2==0).or.(nam%local_diag)) then
   ! Normalize
   !$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,bpar%nc3(ib)
            if (avg_blk_lr%nc1a(jc3,jl0r,il0)>0.0) then
               do isub=1,avg_blk_lr%nsub
                  do jsub=1,avg_blk_lr%nsub
                     avg_blk_lr%m11lrm11sub(jc3,jl0r,il0,jsub,isub) = avg_blk_lr%m11lrm11sub(jc3,jl0r,il0,jsub,isub) &
                                                                    & /avg_blk_lr%nc1a(jc3,jl0r,il0)
                  end do
               end do
            end if
         end do
      end do
   end do
   !$omp end parallel do

   ! Average over sub-ensembles
   !$omp parallel do schedule(static) private(il0,jl0r,jc3)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,bpar%nc3(ib)
            if (isanynotmsr(avg_blk_lr%m11lrm11sub(jc3,jl0r,il0,:,:))) then
               avg_blk_lr%m11lrm11(jc3,jl0r,il0) = sum(avg_blk_lr%m11lrm11sub(jc3,jl0r,il0,:,:), &
                                                 & mask=isnotmsr(avg_blk_lr%m11lrm11sub(jc3,jl0r,il0,:,:))) &
                                                 & /real(count(isnotmsr(avg_blk_lr%m11lrm11sub(jc3,jl0r,il0,:,:))),kind_real)
            end if
         end do
      end do
   end do
   !$omp end parallel do

   ! Define asymptotic covariance
   avg_blk_lr%m11lrm11asy = avg_blk_lr%m11lrm11
end if

! End associate
end associate

end subroutine avg_blk_compute_asy_lr

end module type_avg_blk
