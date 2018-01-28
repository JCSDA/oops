!----------------------------------------------------------------------
! Module: hdiag_average.f90
!> Purpose: average routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_average

use omp_lib
use tools_const, only: add,divide
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_qsort, only: qsort
use type_avg, only: avgtype,avg_alloc,avg_copy,avg_pack,avg_unpack
use type_mom, only: momtype
use type_mpl, only: mpl,mpl_allreduce_sum
use type_hdata, only: hdatatype
implicit none

private
public :: compute_avg,compute_avg_lr,compute_avg_asy,compute_bwavg

contains

!----------------------------------------------------------------------
! Subroutine: compute_avg
!> Purpose: compute averaged statistics via spatial-angular erogodicity assumption
!----------------------------------------------------------------------
subroutine compute_avg(hdata,ib,mom,ic2a,avg)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata !< HDIAG data
integer,intent(in) :: ib            !< Block index
type(momtype),intent(in) :: mom     !< Moments
integer,intent(in) :: ic2a          !< Subgrid index
type(avgtype),intent(inout) :: avg  !< Averaged statistics

! Local variables
integer :: il0,jl0,jl0r,ic2,jc3,isub,jsub,ic1a,ic1,nc1amax,nc1a
real(kind_real) :: m2m2
real(kind_real),allocatable :: list_m11(:),list_m11m11(:,:,:),list_m2m2(:,:,:),list_m22(:,:),list_cor(:),sbuf(:),rbuf(:)
logical :: valid

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Copy ensemble size
avg%ne = mom%ne
avg%nsub = mom%nsub

! Allocation
call avg_alloc(hdata,ib,avg)

! Global index
if (ic2a==0) then
   ic2 = 0
else
   ic2 = hdata%c2a_to_c2(ic2a)
end if

! Average
!$omp parallel do schedule(static) private(il0,jl0r,jl0,nc1amax,list_m11,list_m11m11,list_m2m2,list_m22,list_cor,jc3), &
!$omp&                             private(nc1a,ic1a,ic1,valid,m2m2,isub,jsub)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)

      ! Allocation
      if (ic2>0) then
         nc1amax = count(hdata%local_mask(hdata%c1a_to_c1,ic2,min(jl0,geom%nl0i)))
      else
         nc1amax = hdata%nc1a
      end if
      allocate(list_m11(nc1amax))
      allocate(list_m11m11(nc1amax,avg%nsub,avg%nsub))
      allocate(list_m2m2(nc1amax,avg%nsub,avg%nsub))
      allocate(list_m22(nc1amax,avg%nsub))
      allocate(list_cor(nc1amax))

      do jc3=1,bpar%nc3(ib)
         ! Fill lists
         nc1a = 0
         do ic1a=1,hdata%nc1a
            ! Index
            ic1 = hdata%c1a_to_c1(ic1a)

            ! Check validity
            valid = hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)
            if (ic2>0) valid = valid.and.hdata%local_mask(ic1,ic2,min(jl0,geom%nl0i))

            if (valid) then
               ! Update
               nc1a = nc1a+1

               ! Averages for diagnostics
               list_m11(nc1a) = sum(mom%m11(ic1a,jc3,jl0r,il0,:))/float(avg%nsub)
               do isub=1,avg%nsub
                  do jsub=1,avg%nsub
                     list_m11m11(nc1a,jsub,isub) = mom%m11(ic1a,jc3,jl0r,il0,isub)*mom%m11(ic1a,jc3,jl0r,il0,jsub)
                     list_m2m2(nc1a,jsub,isub) = mom%m2_1(ic1a,jc3,jl0r,il0,isub)*mom%m2_2(ic1a,jc3,jl0r,il0,jsub)
                  end do
                  if (.not.nam%gau_approx) list_m22(nc1a,isub) = mom%m22(ic1a,jc3,jl0r,il0,isub)
               end do

               ! Correlation
               m2m2 = sum(mom%m2_1(ic1a,jc3,jl0r,il0,:))*sum(mom%m2_2(ic1a,jc3,jl0r,il0,:))/float(avg%nsub**2)
               if (m2m2>0.0) then
                  list_cor(nc1a) = list_m11(nc1a)/sqrt(m2m2)
               else
                  call msr(list_cor(nc1a))
               end if
            end if
         end do

         ! Average
         avg%nc1a(jc3,jl0r,il0) = float(nc1a)
         avg%m11(jc3,jl0r,il0) = sum(list_m11(1:nc1a))
         do isub=1,avg%nsub
            do jsub=1,avg%nsub
               avg%m11m11(jc3,jl0r,il0,jsub,isub) = sum(list_m11m11(1:nc1a,jsub,isub))
               avg%m2m2(jc3,jl0r,il0,jsub,isub) = sum(list_m2m2(1:nc1a,jsub,isub))
            end do
            if (.not.nam%gau_approx) avg%m22(jc3,jl0r,il0,isub) = sum(list_m22(1:nc1a,isub))
         end do
         avg%cor(jc3,jl0r,il0) = sum(list_cor(1:nc1a))
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

if (ic2a==0) then
   ! Allocation
   allocate(sbuf(avg%npack))
   allocate(rbuf(avg%npack))

   ! Pack data
   call avg_pack(hdata,ib,avg,sbuf)

   ! Allreduce
   call mpl_allreduce_sum(sbuf,rbuf)

   ! Unpack data
   call avg_unpack(hdata,ib,avg,rbuf)
end if

! Normalize
!$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (avg%nc1a(jc3,jl0r,il0)>0.0) then
            avg%m11(jc3,jl0r,il0) = avg%m11(jc3,jl0r,il0)/avg%nc1a(jc3,jl0r,il0)
            do isub=1,avg%nsub
               do jsub=1,avg%nsub
                  avg%m11m11(jc3,jl0r,il0,jsub,isub) = avg%m11m11(jc3,jl0r,il0,jsub,isub)/avg%nc1a(jc3,jl0r,il0)
                  avg%m2m2(jc3,jl0r,il0,jsub,isub) = avg%m2m2(jc3,jl0r,il0,jsub,isub)/avg%nc1a(jc3,jl0r,il0)
               end do
               if (.not.nam%gau_approx) avg%m22(jc3,jl0r,il0,isub) = avg%m22(jc3,jl0r,il0,isub)/avg%nc1a(jc3,jl0r,il0)
            end do
            avg%cor(jc3,jl0r,il0) = avg%cor(jc3,jl0r,il0)/avg%nc1a(jc3,jl0r,il0)
         end if
      end do
   end do
end do
!$omp end parallel do

! End associate
end associate

end subroutine compute_avg

!----------------------------------------------------------------------
! Subroutine: compute_avg_lr
!> Purpose: compute averaged statistics via spatial-angular erogodicity assumption, for LR covariance/HR covariance and LR covariance/HR asymptotic covariance products
!----------------------------------------------------------------------
subroutine compute_avg_lr(hdata,ib,mom,mom_lr,ic2a,avg,avg_lr)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata   !< HDIAG data
integer,intent(in) :: ib              !< Block index
type(momtype),intent(in) :: mom       !< Moments
type(momtype),intent(in) :: mom_lr    !< Low-resolution moments
integer,intent(in) :: ic2a            !< Subgrid index
type(avgtype),intent(inout) :: avg    !< Averaged statistics
type(avgtype),intent(inout) :: avg_lr !< Low-resolution averaged statistics

! Local variables
integer :: il0,jl0,jl0r,ic2,jc3,isub,jsub,ic1a,ic1,nc1amax,nc1a
real(kind_real),allocatable :: m11lrm11(:,:,:,:,:),list_m11lrm11(:,:,:),sbuf(:),rbuf(:)
logical :: valid
logical,allocatable :: mask_unpack(:,:,:,:,:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Global index
if (ic2a==0) then
   ic2 = 0
else
   ic2 = hdata%c2a_to_c2(ic2a)
end if

! Allocation
allocate(m11lrm11(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_lr%nsub,avg%nsub))

! Average
!$omp parallel do schedule(static) private(il0,jl0r,jl0,nc1amax,list_m11lrm11,jc3,nc1a,ic1a,ic1,valid,isub,jsub)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)

      ! Allocation
      if (ic2>0) then
         nc1amax = count(hdata%local_mask(hdata%c1a_to_c1,ic2,min(jl0,geom%nl0i)))
      else
         nc1amax = hdata%nc1a
      end if
      allocate(list_m11lrm11(nc1amax,avg_lr%nsub,avg%nsub))

      do jc3=1,bpar%nc3(ib)
         ! Fill lists
         nc1a = 0
         do ic1a=1,hdata%nc1a
            ! Index
            ic1 = hdata%c1a_to_c1(ic1a)

            ! Check validity
            valid = hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)
            if (ic2>0) valid = valid.and.hdata%local_mask(ic1,ic2,min(jl0,geom%nl0i))

            if (valid) then
               ! Update
               nc1a = nc1a+1

               ! Averages for diagnostics
               do isub=1,avg%nsub
                  do jsub=1,avg_lr%nsub
                     list_m11lrm11(nc1a,jsub,isub) = mom%m11(ic1a,jc3,jl0r,il0,isub)*mom_lr%m11(ic1a,jc3,jl0r,il0,jsub)
                  end do
               end do
            end if
         end do

         ! Average
         avg%nc1a(jc3,jl0r,il0) = float(nc1a)
         do isub=1,avg%nsub
            do jsub=1,avg_lr%nsub
               m11lrm11(jc3,jl0r,il0,jsub,isub) = sum(list_m11lrm11(1:nc1a,jsub,isub))
            end do
         end do
      end do

      ! Release memory
      deallocate(list_m11lrm11)
   end do
end do
!$omp end parallel do

if (ic2a==0) then
   ! Allocation
   allocate(sbuf(bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_lr%nsub*avg%nsub))
   allocate(rbuf(bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_lr%nsub*avg%nsub))
   allocate(mask_unpack(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_lr%nsub,avg%nsub))

   ! Pack data
   sbuf = pack(m11lrm11,mask=.true.)

   ! Allreduce
   call mpl_allreduce_sum(sbuf,rbuf)

   ! Unpack data
   mask_unpack = .true.
   m11lrm11 = unpack(rbuf,mask_unpack,m11lrm11)
end if

! Normalize
!$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (avg%nc1a(jc3,jl0r,il0)>0.0) then
            do isub=1,avg%nsub
               do jsub=1,avg_lr%nsub
                  m11lrm11(jc3,jl0r,il0,jsub,isub) = m11lrm11(jc3,jl0r,il0,jsub,isub)/avg%nc1a(jc3,jl0r,il0)
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
         if (isanynotmsr(m11lrm11(jc3,jl0r,il0,:,:))) then
            avg_lr%m11lrm11(jc3,jl0r,il0) = sum(m11lrm11(jc3,jl0r,il0,:,:),mask=isnotmsr(m11lrm11(jc3,jl0r,il0,:,:))) &
                                          & /float(count(isnotmsr(m11lrm11(jc3,jl0r,il0,:,:))))
         end if
      end do
   end do
end do
!$omp end parallel do

! End associate
end associate

end subroutine compute_avg_lr

!----------------------------------------------------------------------
! Subroutine: compute_avg_asy
!> Purpose: compute averaged asymptotic statistics
!----------------------------------------------------------------------
subroutine compute_avg_asy(hdata,ib,ne,avg)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata !< HDIAG data
integer,intent(in) :: ib            !< Block index
integer,intent(in) :: ne            !< Ensemble size
type(avgtype),intent(inout) :: avg  !< Averaged statistics

! Local variables
integer :: n,il0,jl0r,jc3,isub,jsub
real(kind_real) :: P1,P3,P4,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17
real(kind_real),allocatable :: m11asysq(:,:),m2m2asy(:,:),m22asy(:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Ensemble size-dependent coefficients
n = ne
P1 = 1.0/float(n)
P3 = 1.0/float(n*(n-1))
P4 = 1.0/float(n-1)
P14 = float(n**2-2*n+2)/float(n*(n-1))
P16 = float(n)/float(n-1)

! Ensemble/sub-ensemble size-dependent coefficients
n = avg%ne/avg%nsub
P7 = float((n-1)*(n**2-3*n+1))/float(n*(n-2)*(n-3))
P8 = float(n-1)/float(n*(n-2)*(n-3))
P9 = -float(n)/float((n-2)*(n-3))
P10 = -float((n-1)*(2*n-3))/float(n*(n-2)*(n-3))
P11 = float(n*(n**2-2*n+3))/float((n-1)*(n-2)*(n-3))
P12 = float(n*(n-1))/float((n-2)*(n+1))
P13 = -float(n-1)/float((n-2)*(n+1))
P15 = float((n-1)**2)/float(n*(n-3))
P17 = float((n-1)**2)/float((n-2)*(n+1))

! Asymptotic statistics
!$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub,m11asysq,m2m2asy,m22asy)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (avg%nc1a(jc3,jl0r,il0)>0.0) then
            ! Allocation
            allocate(m11asysq(avg%nsub,avg%nsub))
            allocate(m2m2asy(avg%nsub,avg%nsub))
            allocate(m22asy(avg%nsub))

            ! Asymptotic statistics
            do isub=1,avg%nsub
               do jsub=1,avg%nsub
                  if (isub==jsub) then
                     ! Diagonal terms
                     if (nam%gau_approx) then
                        ! Gaussian approximation
                        m11asysq(jsub,isub) = P17*avg%m11m11(jc3,jl0r,il0,jsub,isub)+P13*avg%m2m2(jc3,jl0r,il0,jsub,isub)
                        m2m2asy(jsub,isub) = 2.0*P13*avg%m11m11(jc3,jl0r,il0,jsub,isub)+P12*avg%m2m2(jc3,jl0r,il0,jsub,isub)
                     else
                        ! General case
                        m11asysq(jsub,isub) = P15*avg%m11m11(jc3,jl0r,il0,jsub,isub)+P8*avg%m2m2(jc3,jl0r,il0,jsub,isub) &
                                            & +P9*avg%m22(jc3,jl0r,il0,jsub)
                        m2m2asy(jsub,isub) = 2.0*P8*avg%m11m11(jc3,jl0r,il0,jsub,isub)+P7*avg%m2m2(jc3,jl0r,il0,jsub,isub) &
                                           & +P9*avg%m22(jc3,jl0r,il0,jsub)
                        m22asy(jsub) = P10*(2.0*avg%m11m11(jc3,jl0r,il0,jsub,isub)+avg%m2m2(jc3,jl0r,il0,jsub,isub)) &
                                     & +P11*avg%m22(jc3,jl0r,il0,jsub)
                     end if
                  else
                     ! Off-diagonal terms
                     m11asysq(jsub,isub) = avg%m11m11(jc3,jl0r,il0,jsub,isub)
                     m2m2asy(jsub,isub) = avg%m2m2(jc3,jl0r,il0,jsub,isub)
                  end if
               end do
            end do

            ! Sum
            avg%m11asysq(jc3,jl0r,il0) = sum(m11asysq)/float(avg%nsub**2)
            avg%m2m2asy(jc3,jl0r,il0) = sum(m2m2asy)/float(avg%nsub**2)
            if (.not.nam%gau_approx) avg%m22asy(jc3,jl0r,il0) = sum(m22asy)/float(avg%nsub)

            ! Check positivity
            if (.not.(avg%m11asysq(jc3,jl0r,il0)>0.0)) call msr(avg%m11asysq(jc3,jl0r,il0))
            if (.not.(avg%m2m2asy(jc3,jl0r,il0)>0.0)) call msr(avg%m2m2asy(jc3,jl0r,il0))
            if (.not.nam%gau_approx) then
               if (.not.(avg%m22asy(jc3,jl0r,il0)>0.0)) call msr(avg%m22asy(jc3,jl0r,il0))
            end if

            ! Squared covariance average
            if (nam%gau_approx) then
               ! Gaussian approximation
               if (isnotmsr(avg%m11asysq(jc3,jl0r,il0)).and.isnotmsr(avg%m2m2asy(jc3,jl0r,il0))) &
             & avg%m11sq(jc3,jl0r,il0) = P16*avg%m11asysq(jc3,jl0r,il0)+P4*avg%m2m2asy(jc3,jl0r,il0)
            else
               ! General case
               if (isnotmsr(avg%m22asy(jc3,jl0r,il0)).and.isnotmsr(avg%m11asysq(jc3,jl0r,il0)) &
             & .and.isnotmsr(avg%m2m2asy(jc3,jl0r,il0))) &
             & avg%m11sq(jc3,jl0r,il0) = P1*avg%m22asy(jc3,jl0r,il0)+P14*avg%m11asysq(jc3,jl0r,il0) &
                                      & +P3*avg%m2m2asy(jc3,jl0r,il0)
            end if

            ! Check value
            if (.not.isnotmsr(avg%m11sq(jc3,jl0r,il0))) then
               if (avg%m11sq(jc3,jl0r,il0)<avg%m11asysq(jc3,jl0r,il0)) call msr(avg%m11sq(jc3,jl0r,il0))
               if (avg%m11sq(jc3,jl0r,il0)<avg%m11(jc3,jl0r,il0)**2) call msr(avg%m11sq(jc3,jl0r,il0))
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

! End associate
end associate

end subroutine compute_avg_asy

!----------------------------------------------------------------------
! Subroutine: compute_bwavg
!> Purpose: compute block-averaged statistics
!----------------------------------------------------------------------
subroutine compute_bwavg(hdata,avg_wgt,avg)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                   !< HDIAG data
type(avgtype),intent(in) :: avg_wgt(hdata%bpar%nb)    !< Averaged statistics for weights
type(avgtype),intent(inout) :: avg(hdata%bpar%nb+1)   !< Averaged statistics

! Local variables
integer :: ib,il0,jl0r,jc3
real(kind_real) :: bwgtsq
real(kind_real),allocatable :: cor(:,:,:),m11asysq(:,:,:),m11sq(:,:,:)
real(kind_real),allocatable :: m11sta(:,:,:),stasq(:,:,:)
real(kind_real),allocatable :: m11lrm11(:,:,:),m11lrm11asy(:,:,:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Copy ensemble size
avg(bpar%nb+1)%ne = avg(1)%ne
avg(bpar%nb+1)%nsub = avg(1)%nsub

! Allocation
call avg_alloc(hdata,bpar%nb+1,avg(bpar%nb+1))

! Allocation
allocate(cor(nam%nc3,nam%nl0r,geom%nl0))
allocate(m11asysq(nam%nc3,nam%nl0r,geom%nl0))
allocate(m11sq(nam%nc3,nam%nl0r,geom%nl0))
select case (trim(nam%method))
case ('hyb-avg','hyb-rnd')
   allocate(m11sta(nam%nc3,nam%nl0r,geom%nl0))
   allocate(stasq(nam%nc3,nam%nl0r,geom%nl0))
case ('dual-ens')
   allocate(m11lrm11(nam%nc3,nam%nl0r,geom%nl0))
   allocate(m11lrm11asy(nam%nc3,nam%nl0r,geom%nl0))
end select

! Initialization
avg(bpar%nb+1)%cor = 0.0
cor = 0.0
avg(bpar%nb+1)%m11asysq = 0.0
m11asysq = 0.0
avg(bpar%nb+1)%m11sq = 0.0
m11sq = 0.0
select case (trim(nam%method))
case ('hyb-avg','hyb-rnd')
   avg(bpar%nb+1)%m11sta = 0.0
   m11sta = 0.0
   avg(bpar%nb+1)%stasq = 0.0
   stasq = 0.0
case ('dual-ens')
   avg(bpar%nb+1)%m11lrm11 = 0.0
   m11lrm11 = 0.0
   avg(bpar%nb+1)%m11lrm11asy = 0.0
   m11lrm11asy = 0.0
end select

! Block averages
do ib=1,bpar%nb
   if (bpar%avg_block(ib)) then
      !$omp parallel do schedule(static) private(il0,jl0r,bwgtsq,jc3)
      do il0=1,geom%nl0
         do jl0r=1,nam%nl0r
            ! Weight
            if (avg(ib)%m2m2asy(1,jl0r,il0)>0.0) then
               bwgtsq = 1.0/avg_wgt(ib)%m2m2asy(1,jl0r,il0)
            else
               bwgtsq = 0.0
            end if

            ! Compute sum
            do jc3=1,nam%nc3
               call add(avg(ib)%cor(jc3,jl0r,il0),avg(bpar%nb+1)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
               call add(avg(ib)%m11asysq(jc3,jl0r,il0),avg(bpar%nb+1)%m11asysq(jc3,jl0r,il0),m11asysq(jc3,jl0r,il0),bwgtsq)
               call add(avg(ib)%m11sq(jc3,jl0r,il0),avg(bpar%nb+1)%m11sq(jc3,jl0r,il0),m11sq(jc3,jl0r,il0),bwgtsq)
               select case (trim(nam%method))
               case ('hyb-avg','hyb-rnd')
                  call add(avg(ib)%m11sta(jc3,jl0r,il0),avg(bpar%nb+1)%m11sta(jc3,jl0r,il0),m11sta(jc3,jl0r,il0),bwgtsq)
                  call add(avg(ib)%stasq(jc3,jl0r,il0),avg(bpar%nb+1)%stasq(jc3,jl0r,il0),stasq(jc3,jl0r,il0),bwgtsq)
               case ('dual-ens')
                  call add(avg(ib)%m11lrm11(jc3,jl0r,il0),avg(bpar%nb+1)%m11lrm11(jc3,jl0r,il0),m11lrm11(jc3,jl0r,il0),bwgtsq)
                  call add(avg(ib)%m11lrm11asy(jc3,jl0r,il0),avg(bpar%nb+1)%m11lrm11asy(jc3,jl0r,il0),m11lrm11asy(jc3,jl0r,il0), &
                & bwgtsq)
               end select
            end do
         end do
      end do
      !$omp end parallel do
   end if
end do

! Normalization
!$omp parallel do schedule(static) private(il0,jl0r,jc3)
do il0=1,geom%nl0
   do jl0r=1,nam%nl0r
      do jc3=1,nam%nc3
         call divide(avg(bpar%nb+1)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
         call divide(avg(bpar%nb+1)%m11asysq(jc3,jl0r,il0),m11asysq(jc3,jl0r,il0))
         call divide(avg(bpar%nb+1)%m11sq(jc3,jl0r,il0),m11sq(jc3,jl0r,il0))
         select case (trim(nam%method))
         case ('hyb-avg','hyb-rnd')
            call divide(avg(bpar%nb+1)%m11sta(jc3,jl0r,il0),m11sta(jc3,jl0r,il0))
            call divide(avg(bpar%nb+1)%stasq(jc3,jl0r,il0),stasq(jc3,jl0r,il0))
         case ('dual-ens')
            call divide(avg(bpar%nb+1)%m11lrm11(jc3,jl0r,il0),m11lrm11(jc3,jl0r,il0))
            call divide(avg(bpar%nb+1)%m11lrm11asy(jc3,jl0r,il0),m11lrm11asy(jc3,jl0r,il0))
         end select
      end do
   end do
end do
!$omp end parallel do

! End associate
end associate

end subroutine compute_bwavg

end module hdiag_average
