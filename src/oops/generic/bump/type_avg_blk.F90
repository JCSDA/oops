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
use tools_missing, only: msr,isanynotmsr,isnotmsr
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl
use type_nam, only: nam_type
implicit none

! Averaged statistics block derived type
type avg_blk_type
   integer :: ic2                                    !< Global index
   integer :: ib                                     !< Block index
   integer :: ne                                     !< Ensemble size
   integer :: nsub                                   !< Sub-ensembles number
   integer :: npack                                  !< Pack format size
   real(kind_real),allocatable :: nc1a(:,:,:)        !< Number of points in subset Sc1 on halo A
   real(kind_real),allocatable :: m11(:,:,:)         !< Covariance average
   real(kind_real),allocatable :: m11m11(:,:,:,:,:)  !< Product of covariances average
   real(kind_real),allocatable :: m2m2(:,:,:,:,:)    !< Product of variances average
   real(kind_real),allocatable :: m22(:,:,:,:)       !< Fourth-order centered moment average
   real(kind_real),allocatable :: nc1a_cor(:,:,:)    !< Number of points in subset Sc1 on halo A with valid correlations
   real(kind_real),allocatable :: cor(:,:,:)         !< Correlation average
   real(kind_real),allocatable :: m11asysq(:,:,:)    !< Squared asymptotic covariance average
   real(kind_real),allocatable :: m2m2asy(:,:,:)     !< Product of asymptotic variances average
   real(kind_real),allocatable :: m22asy(:,:,:)      !< Asymptotic fourth-order centered moment average
   real(kind_real),allocatable :: m11sq(:,:,:)       !< Squared covariance average for several ensemble sizes
   real(kind_real),allocatable :: m11sta(:,:,:)      !< Ensemble covariance/static covariance product
   real(kind_real),allocatable :: stasq(:,:,:)       !< Squared static covariance
   real(kind_real),allocatable :: m11lrm11(:,:,:)    !< LR covariance/HR covariance product average
   real(kind_real),allocatable :: m11lrm11asy(:,:,:) !< LR covariance/HR asymptotic covariance product average
contains
   procedure :: alloc => avg_blk_alloc
   procedure :: dealloc => avg_blk_dealloc
   procedure :: pack => avg_blk_pack
   procedure :: unpack => avg_blk_unpack
   procedure :: compute => avg_blk_compute
   procedure :: compute_lr => avg_blk_compute_lr
end type avg_blk_type

private
public :: avg_blk_type

real(kind_real),parameter :: var_min = 1.0e-24_kind_real

contains

!----------------------------------------------------------------------
! Subroutine: avg_blk_alloc
!> Purpose: averaged statistics object allocation
!----------------------------------------------------------------------
subroutine avg_blk_alloc(avg_blk,nam,geom,bpar,hdata,ic2a,ib,ne,nsub)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk !< Averaged statistics block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
type(hdata_type),intent(in) :: hdata         !< HDIAG data
integer,intent(in) :: ne                     !< Ensemble size
integer,intent(in) :: nsub                   !< Sub-ensembles number
integer,intent(in) :: ic2a                   !< Local index
integer,intent(in) :: ib                     !< Block index

! Set attributes
if (ic2a==0) then
   avg_blk%ic2 = 0
else
   avg_blk%ic2 = hdata%c2a_to_c2(ic2a)
end if
avg_blk%ib = ib
avg_blk%ne = ne
avg_blk%nsub = nsub
avg_blk%npack = (4+2*avg_blk%nsub**2)*bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
if (.not.nam%gau_approx) avg_blk%npack = avg_blk%npack+avg_blk%nsub*bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0

! Allocation
if (.not.allocated(avg_blk%nc1a)) then
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
   case ('hyb-avg_blk','hyb-rnd')
      allocate(avg_blk%m11sta(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%stasq(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
   case ('dual-ens')
      allocate(avg_blk%m11lrm11(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11lrm11asy(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
   end select
end if

! Initialization
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
   call msr(avg_blk%m11lrm11)
   call msr(avg_blk%m11lrm11asy)
end select

end subroutine avg_blk_alloc

!----------------------------------------------------------------------
! Subroutine: avg_blk_dealloc
!> Purpose: averaged statistics object deallocation
!----------------------------------------------------------------------
subroutine avg_blk_dealloc(avg_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk !< Averaged statistics block

! Allocation
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
if (allocated(avg_blk%m11lrm11)) deallocate(avg_blk%m11lrm11)
if (allocated(avg_blk%m11lrm11asy)) deallocate(avg_blk%m11lrm11asy)

end subroutine avg_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: avg_blk_pack
!> Purpose: averaged statistics object packing
!----------------------------------------------------------------------
subroutine avg_blk_pack(avg_blk,nam,geom,bpar,buf)

implicit none

! Passed variables
class(avg_blk_type),intent(in) :: avg_blk         !< Averaged statistics block
type(nam_type),intent(in) :: nam                  !< Namelist
type(geom_type),intent(in) :: geom                !< Geometry
type(bpar_type),intent(in) :: bpar                !< Block parameters
real(kind_real),intent(out) :: buf(avg_blk%npack) !< Buffer

! Local variables
integer :: offset

! Associate
associate(ib=>avg_blk%ib)

! Pack
offset = 0
buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg_blk%nc1a,.true.)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg_blk%m11,.true.)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2) = pack(avg_blk%m11m11,.true.)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2
buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2) = pack(avg_blk%m2m2,.true.)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2
if (.not.nam%gau_approx) then
   buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub) = pack(avg_blk%m22,.true.)
   offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub
end if
buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg_blk%nc1a_cor,.true.)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg_blk%cor,.true.)

! End associate
end associate

end subroutine avg_blk_pack

!----------------------------------------------------------------------
! Subroutine: avg_blk_unpack
!> Purpose: averaged statistics object unpacking
!----------------------------------------------------------------------
subroutine avg_blk_unpack(avg_blk,nam,geom,bpar,buf)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk     !< Averaged statistics block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(bpar_type),intent(in) :: bpar               !< Block parameters
real(kind_real),intent(in) :: buf(avg_blk%npack) !< Buffer

! Local variables
integer :: offset
logical,allocatable :: mask_0(:,:,:),mask_1(:,:,:,:),mask_2(:,:,:,:,:)

! Associate
associate(ib=>avg_blk%ib)

! Allocation
allocate(mask_0(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
if (.not.nam%gau_approx) allocate(mask_1(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk%nsub))
allocate(mask_2(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk%nsub,avg_blk%nsub))
mask_0 = .true.
if (.not.nam%gau_approx) mask_1 = .true.
mask_2 = .true.

! Unpack
offset = 0
avg_blk%nc1a = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg_blk%m11)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
avg_blk%m11 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg_blk%m11)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
avg_blk%m11m11 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2),mask_2,avg_blk%m11m11)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2
avg_blk%m2m2 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2),mask_2,avg_blk%m2m2)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub**2
if (.not.nam%gau_approx) then
   avg_blk%m22 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub),mask_1,avg_blk%m22)
   offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk%nsub
end if
avg_blk%nc1a_cor = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg_blk%nc1a_cor)
offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
avg_blk%cor = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg_blk%cor)

! End associate
end associate

end subroutine avg_blk_unpack

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute
!> Purpose: compute averaged statistics via spatial-angular erogodicity assumption
!----------------------------------------------------------------------
subroutine avg_blk_compute(avg_blk,nam,geom,bpar,hdata,mom_blk,ne)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk !< Averaged statistics block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
type(hdata_type),intent(in) :: hdata         !< HDIAG data
type(mom_blk_type),intent(in) :: mom_blk     !< Moments
integer,intent(in) :: ne                     !< Ensemble size

! Local variables
integer :: il0,jl0,jl0r,jc3,isub,jsub,ic1a,ic1,nc1amax,nc1a,n
real(kind_real) :: m2_1,m2_2
real(kind_real) :: P1,P3,P4,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17
real(kind_real),allocatable :: list_m11(:),list_m11m11(:,:,:),list_m2m2(:,:,:),list_m22(:,:),list_cor(:),sbuf(:),rbuf(:)
real(kind_real),allocatable :: m11asysq(:,:),m2m2asy(:,:),m22asy(:)
logical :: valid

! Associate
associate(ic2=>avg_blk%ic2,ib=>avg_blk%ib)

! Average
!$omp parallel do schedule(static) private(il0,jl0r,jl0,nc1amax,jc3,nc1a,ic1a,ic1,valid,m2_1,m2_2,isub,jsub), &
!$omp&                             firstprivate(list_m11,list_m11m11,list_m2m2,list_m22,list_cor)
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
            if (ic2>0) valid = valid.and.hdata%local_mask(ic1,ic2,min(jl0,geom%nl0i))

            if (valid) then
               ! Update
               nc1a = nc1a+1

               ! Averages for diagnostics
               list_m11(nc1a) = sum(mom_blk%m11(ic1a,jc3,jl0r,il0,:))/real(avg_blk%nsub,kind_real)
               do isub=1,avg_blk%nsub
                  do jsub=1,avg_blk%nsub
                     list_m11m11(nc1a,jsub,isub) = mom_blk%m11(ic1a,jc3,jl0r,il0,isub)*mom_blk%m11(ic1a,jc3,jl0r,il0,jsub)
                     list_m2m2(nc1a,jsub,isub) = mom_blk%m2_1(ic1a,jc3,jl0r,il0,isub)*mom_blk%m2_2(ic1a,jc3,jl0r,il0,jsub)
                  end do
                  if (.not.nam%gau_approx) list_m22(nc1a,isub) = mom_blk%m22(ic1a,jc3,jl0r,il0,isub)
               end do

               ! Correlation
               m2_1 = sum(mom_blk%m2_1(ic1a,jc3,jl0r,il0,:))/real(avg_blk%nsub,kind_real)
               m2_2 = sum(mom_blk%m2_2(ic1a,jc3,jl0r,il0,:))/real(avg_blk%nsub,kind_real)
               if ((m2_1>var_min).and.(m2_2>var_min)) then
                  list_cor(nc1a) = list_m11(nc1a)/sqrt(m2_1*m2_2)
                  if (abs(list_cor(nc1a))>1.0) call msr(list_cor(nc1a))
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
            if (avg_blk%nc1a_cor(jc3,jl0r,il0)>0.0) then
               avg_blk%cor(jc3,jl0r,il0) = sum(list_cor(1:nc1a),mask=isnotmsr(list_cor(1:nc1a)))
            else
               call msr(avg_blk%cor(jc3,jl0r,il0))
            end if
         else
            ! Average
            call msr(avg_blk%m11(jc3,jl0r,il0))
            do isub=1,avg_blk%nsub
               do jsub=1,avg_blk%nsub
                  call msr(avg_blk%m11m11(jc3,jl0r,il0,jsub,isub))
                  call msr(avg_blk%m2m2(jc3,jl0r,il0,jsub,isub))
               end do
               if (.not.nam%gau_approx) call msr(avg_blk%m22(jc3,jl0r,il0,isub))
            end do
            avg_blk%nc1a_cor(jc3,jl0r,il0) = 0.0
            call msr(avg_blk%cor(jc3,jl0r,il0))
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

if (ic2==0) then
   ! Allocation
   allocate(sbuf(avg_blk%npack))
   allocate(rbuf(avg_blk%npack))

   ! Pack data
   call avg_blk%pack(nam,geom,bpar,sbuf)

   ! Allreduce
   call mpl%allreduce_sum(sbuf,rbuf)

   ! Unpack data
   call avg_blk%unpack(nam,geom,bpar,rbuf)
end if

! Normalize
!$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (avg_blk%nc1a(jc3,jl0r,il0)>0.0) then
            avg_blk%m11(jc3,jl0r,il0) = avg_blk%m11(jc3,jl0r,il0)/avg_blk%nc1a(jc3,jl0r,il0)
            do isub=1,avg_blk%nsub
               do jsub=1,avg_blk%nsub
                  avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) = avg_blk%m11m11(jc3,jl0r,il0,jsub,isub)/avg_blk%nc1a(jc3,jl0r,il0)
                  avg_blk%m2m2(jc3,jl0r,il0,jsub,isub) = avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)/avg_blk%nc1a(jc3,jl0r,il0)
               end do
               if (.not.nam%gau_approx) avg_blk%m22(jc3,jl0r,il0,isub) = avg_blk%m22(jc3,jl0r,il0,isub)/avg_blk%nc1a(jc3,jl0r,il0)
            end do
         end if
         if (avg_blk%nc1a_cor(jc3,jl0r,il0)>0.0) then
            avg_blk%cor(jc3,jl0r,il0) = avg_blk%cor(jc3,jl0r,il0)/avg_blk%nc1a_cor(jc3,jl0r,il0)
         else
            call msr(avg_blk%cor(jc3,jl0r,il0))
         end if
      end do
   end do
end do
!$omp end parallel do

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
         if (avg_blk%nc1a(jc3,jl0r,il0)>0.0) then
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
                        m11asysq(jsub,isub) = P17*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub)+P13*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)
                        m2m2asy(jsub,isub) = 2.0*P13*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub)+P12*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub)
                     else
                        ! General case
                        m11asysq(jsub,isub) = P15*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub)+P8*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub) &
                                            & +P9*avg_blk%m22(jc3,jl0r,il0,jsub)
                        m2m2asy(jsub,isub) = 2.0*P8*avg_blk%m11m11(jc3,jl0r,il0,jsub,isub)+P7*avg_blk%m2m2(jc3,jl0r,il0,jsub,isub) &
                                           & +P9*avg_blk%m22(jc3,jl0r,il0,jsub)
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
            if (.not.(avg_blk%m11asysq(jc3,jl0r,il0)>0.0)) call msr(avg_blk%m11asysq(jc3,jl0r,il0))
            if (.not.(avg_blk%m2m2asy(jc3,jl0r,il0)>0.0)) call msr(avg_blk%m2m2asy(jc3,jl0r,il0))
            if (.not.nam%gau_approx) then
               if (.not.(avg_blk%m22asy(jc3,jl0r,il0)>0.0)) call msr(avg_blk%m22asy(jc3,jl0r,il0))
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
            if (.not.isnotmsr(avg_blk%m11sq(jc3,jl0r,il0))) then
               if (avg_blk%m11sq(jc3,jl0r,il0)<avg_blk%m11asysq(jc3,jl0r,il0)) call msr(avg_blk%m11sq(jc3,jl0r,il0))
               if (avg_blk%m11sq(jc3,jl0r,il0)<avg_blk%m11(jc3,jl0r,il0)**2) call msr(avg_blk%m11sq(jc3,jl0r,il0))
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

end subroutine avg_blk_compute

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_lr
!> Purpose: compute averaged statistics via spatial-angular erogodicity assumption, for LR covariance/HR covariance and LR covariance/HR asymptotic covariance products
!----------------------------------------------------------------------
subroutine avg_blk_compute_lr(avg_blk_lr,geom,bpar,hdata,mom_blk,mom_lr_blk,avg_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk_lr !< Low-resolution averaged statistics block
type(geom_type),intent(in) :: geom              !< Geometry
type(bpar_type),intent(in) :: bpar              !< Block parameters
type(hdata_type),intent(in) :: hdata            !< HDIAG data
type(mom_blk_type),intent(in) :: mom_blk        !< Moments block
type(mom_blk_type),intent(in) :: mom_lr_blk     !< Low-resolution moments block
type(avg_blk_type),intent(inout) :: avg_blk     !< Averaged statistics block

! Local variables
integer :: il0,jl0,jl0r,jc3,isub,jsub,ic1a,ic1,nc1amax,nc1a
real(kind_real),allocatable :: m11lrm11(:,:,:,:,:),list_m11lrm11(:,:,:),sbuf(:),rbuf(:)
logical :: valid
logical,allocatable :: mask_unpack(:,:,:,:,:)

! Associate
associate(ic2=>avg_blk_lr%ic2,ib=>avg_blk_lr%ib)

! Allocation
allocate(m11lrm11(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk_lr%nsub,avg_blk%nsub))

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
      allocate(list_m11lrm11(nc1amax,avg_blk_lr%nsub,avg_blk%nsub))

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
               do isub=1,avg_blk%nsub
                  do jsub=1,avg_blk_lr%nsub
                     list_m11lrm11(nc1a,jsub,isub) = mom_blk%m11(ic1a,jc3,jl0r,il0,isub)*mom_lr_blk%m11(ic1a,jc3,jl0r,il0,jsub)
                  end do
               end do
            end if
         end do

         ! Average
         avg_blk%nc1a(jc3,jl0r,il0) = real(nc1a,kind_real)
         do isub=1,avg_blk%nsub
            do jsub=1,avg_blk_lr%nsub
               m11lrm11(jc3,jl0r,il0,jsub,isub) = sum(list_m11lrm11(1:nc1a,jsub,isub))
            end do
         end do
      end do

      ! Release memory
      deallocate(list_m11lrm11)
   end do
end do
!$omp end parallel do

if (ic2==0) then
   ! Allocation
   allocate(sbuf(bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk_lr%nsub*avg_blk%nsub))
   allocate(rbuf(bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_blk_lr%nsub*avg_blk%nsub))
   allocate(mask_unpack(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_blk_lr%nsub,avg_blk%nsub))

   ! Pack data
   sbuf = pack(m11lrm11,mask=.true.)

   ! Allreduce
   call mpl%allreduce_sum(sbuf,rbuf)

   ! Unpack data
   mask_unpack = .true.
   m11lrm11 = unpack(rbuf,mask_unpack,m11lrm11)
end if

! Normalize
!$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (avg_blk%nc1a(jc3,jl0r,il0)>0.0) then
            do isub=1,avg_blk%nsub
               do jsub=1,avg_blk_lr%nsub
                  m11lrm11(jc3,jl0r,il0,jsub,isub) = m11lrm11(jc3,jl0r,il0,jsub,isub)/avg_blk%nc1a(jc3,jl0r,il0)
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
            avg_blk_lr%m11lrm11(jc3,jl0r,il0) = sum(m11lrm11(jc3,jl0r,il0,:,:),mask=isnotmsr(m11lrm11(jc3,jl0r,il0,:,:))) &
                                          & /real(count(isnotmsr(m11lrm11(jc3,jl0r,il0,:,:))),kind_real)
         end if
      end do
   end do
end do
!$omp end parallel do

! End associate
end associate

end subroutine avg_blk_compute_lr

end module type_avg_blk
