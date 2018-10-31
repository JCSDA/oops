!----------------------------------------------------------------------
! Module: type_avg
! Purpose: average routines
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_avg

use fckit_mpi_module, only: fckit_mpi_sum
!$ use omp_lib
use tools_const,only: reqkm
use tools_func, only: add,divide
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_qsort, only: qsort
use type_avg_blk, only: avg_blk_type
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_mom, only: mom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

! Averaged statistics derived type
type avg_type
   integer :: ne                              ! Ensemble size
   integer :: nsub                            ! Number of sub-ensembles
   type(avg_blk_type),allocatable :: blk(:,:) ! Averaged statistics blocks
contains
   procedure :: alloc => avg_alloc
   procedure :: dealloc => avg_dealloc
   procedure :: copy => avg_copy
   procedure :: gather => avg_gather
   procedure :: normalize => avg_normalize
   procedure :: gather_lr => avg_gather_lr
   procedure :: normalize_lr => avg_normalize_lr
   procedure :: var_filter => avg_var_filter
   procedure :: compute => avg_compute
   procedure :: compute_hyb => avg_compute_hyb
   procedure :: copy_wgt => avg_copy_wgt
   procedure :: compute_bwavg => avg_compute_bwavg
end type avg_type

private
public :: avg_type

contains

!----------------------------------------------------------------------
! Subroutine: avg_alloc
! Purpose: averaged statistics allocation
!----------------------------------------------------------------------
subroutine avg_alloc(avg,nam,geom,bpar,ne,nsub)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
integer,intent(in) :: ne             ! Ensemble size
integer,intent(in) :: nsub           ! Number of sub-ensembles

! Local variables
integer :: ib,ic2

! Set attributes
avg%ne = ne
avg%nsub = nsub

! Allocation
allocate(avg%blk(0:nam%nc2,bpar%nbe))
do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         call avg%blk(ic2,ib)%alloc(nam,geom,bpar,ic2,ib,ne,nsub)
      end do
   end if
end do

end subroutine avg_alloc

!----------------------------------------------------------------------
! Subroutine: avg_dealloc
! Purpose: averaged statistics deallocation
!----------------------------------------------------------------------
subroutine avg_dealloc(avg,nam,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(nam_type),intent(in) :: nam     ! Namelist
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib,ic2

! Allocation
if (allocated(avg%blk)) then
   do ib=1,bpar%nbe
      if (bpar%diag_block(ib)) then
         do ic2=0,nam%nc2
            call avg%blk(ic2,ib)%dealloc
         end do
      end if
   end do
   deallocate(avg%blk)
end if

end subroutine avg_dealloc

!----------------------------------------------------------------------
! Function: avg_copy
! Purpose: averaged statistics copy
!----------------------------------------------------------------------
type(avg_type) function avg_copy(avg,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(in) :: avg ! Averaged statistics
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib,ic2

! Deallocation
call avg_copy%dealloc(nam,bpar)

! Allocation
call avg_copy%alloc(nam,geom,bpar,avg%ne,avg%nsub)

! Copy
do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         avg_copy%blk(ic2,ib) = avg%blk(ic2,ib)%copy(nam,geom,bpar)
      end do
   end if
end do

end function avg_copy

!----------------------------------------------------------------------
! Subroutine: avg_gather
! Purpose: gather averaged statistics data
!----------------------------------------------------------------------
subroutine avg_gather(avg,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(in) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: npack,offset,ib,ic2
real(kind_real),allocatable :: sbuf(:),rbuf(:)
logical,allocatable :: mask_0(:,:,:),mask_1(:,:,:,:),mask_2(:,:,:,:,:)

! Packing size
npack = 0
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         if ((ic2==0).or.(nam%var_diag)) then
            npack = npack+geom%nl0
            if (nam%var_filter.and.(.not.nam%gau_approx)) npack = npack+geom%nl0
         end if
         if ((ic2==0).or.(nam%local_diag)) then
            npack = npack+(4+2*avg%blk(ic2,ib)%nsub**2)*bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
            if (.not.nam%gau_approx) npack = npack+avg%blk(ic2,ib)%nsub*bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         end if
      end do
   end if
end do

! Allocation
allocate(sbuf(npack))
allocate(rbuf(npack))

! Pack data
offset = 0
sbuf = 0.0
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         if ((ic2==0).or.(nam%var_diag)) then
            sbuf(offset+1:offset+geom%nl0) = pack(avg%blk(ic2,ib)%m2,.true.)
            offset = offset+geom%nl0
            if (nam%var_filter.and.(.not.nam%gau_approx)) then
               sbuf(offset+1:offset+geom%nl0) = pack(avg%blk(ic2,ib)%m4,.true.)
               offset = offset+geom%nl0
            end if
         end if
         if ((ic2==0).or.(nam%local_diag)) then
            sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%nc1a,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
            sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%m11,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
            sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2) = pack(avg%blk(ic2,ib)%m11m11,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
            sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2) = pack(avg%blk(ic2,ib)%m2m2,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
            if (.not.nam%gau_approx) then
               sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub) = pack(avg%blk(ic2,ib)%m22,.true.)
               offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub
            end if
            sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%nc1a_cor,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
            sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%cor,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         end if
      end do
   end if
end do

! Reduce data
call mpl%f_comm%allreduce(sbuf,rbuf,fckit_mpi_sum())

! Unpack data
offset = 0
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      ! Allocation
      allocate(mask_0(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(mask_1(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg%blk(0,ib)%nsub))
      allocate(mask_2(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg%blk(0,ib)%nsub,avg%blk(0,ib)%nsub))

      ! Initialization
      mask_0 = .true.
      mask_1 = .true.
      mask_2 = .true.

      do ic2=0,nam%nc2
         ! Unpack
         if ((ic2==0).or.(nam%var_diag)) then
            avg%blk(ic2,ib)%m2 = unpack(rbuf(offset+1:offset+geom%nl0),mask_1(1,1,:,:),avg%blk(ic2,ib)%m2)
            offset = offset+geom%nl0
            if (nam%var_filter.and.(.not.nam%gau_approx)) then
               avg%blk(ic2,ib)%m4 = unpack(rbuf(offset+1:offset+geom%nl0),mask_1(1,1,:,:),avg%blk(ic2,ib)%m4)
               offset = offset+geom%nl0
            end if
         end if
         if ((ic2==0).or.(nam%local_diag)) then
            avg%blk(ic2,ib)%nc1a = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%m11)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
            avg%blk(ic2,ib)%m11 = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%m11)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
            avg%blk(ic2,ib)%m11m11 = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2), &
          & mask_2,avg%blk(ic2,ib)%m11m11)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
            avg%blk(ic2,ib)%m2m2 = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2), &
          & mask_2,avg%blk(ic2,ib)%m2m2)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
            if (.not.nam%gau_approx) then
               avg%blk(ic2,ib)%m22 = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub), &
             & mask_1,avg%blk(ic2,ib)%m22)
               offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub
            end if
            avg%blk(ic2,ib)%nc1a_cor = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0), &
                                     & mask_0,avg%blk(ic2,ib)%nc1a_cor)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
            avg%blk(ic2,ib)%cor = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%cor)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         end if
      end do

      ! Release memory
      deallocate(mask_0)
      deallocate(mask_1)
      deallocate(mask_2)
   end if
end do

end subroutine avg_gather

!----------------------------------------------------------------------
! Subroutine: avg_normalize
! Purpose: normalize averaged statistics data
!----------------------------------------------------------------------
subroutine avg_normalize(avg,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib,ic2,il0,jl0r,jc3,isub,jsub
real(kind_real) :: norm

! Normalize
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         if ((ic2==0).or.(nam%local_diag)) then
            !$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub,norm)
            do il0=1,geom%nl0
               do jl0r=1,bpar%nl0r(ib)
                  do jc3=1,bpar%nc3(ib)
                     if (avg%blk(ic2,ib)%nc1a(jc3,jl0r,il0)>0.0) then
                        norm = 1.0/avg%blk(ic2,ib)%nc1a(jc3,jl0r,il0)
                        avg%blk(ic2,ib)%m11(jc3,jl0r,il0) = avg%blk(ic2,ib)%m11(jc3,jl0r,il0)*norm
                        do isub=1,avg%blk(ic2,ib)%nsub
                           do jsub=1,avg%blk(ic2,ib)%nsub
                              avg%blk(ic2,ib)%m11m11(jc3,jl0r,il0,jsub,isub) = avg%blk(ic2,ib)%m11m11(jc3,jl0r,il0,jsub,isub)*norm
                              avg%blk(ic2,ib)%m2m2(jc3,jl0r,il0,jsub,isub) = avg%blk(ic2,ib)%m2m2(jc3,jl0r,il0,jsub,isub)*norm
                            end do
                           if (.not.nam%gau_approx) avg%blk(ic2,ib)%m22(jc3,jl0r,il0,isub) &
                         & = avg%blk(ic2,ib)%m22(jc3,jl0r,il0,isub)*norm
                        end do
                     else
                        call msr(avg%blk(ic2,ib)%m11(jc3,jl0r,il0))
                        do isub=1,avg%blk(ic2,ib)%nsub
                           do jsub=1,avg%blk(ic2,ib)%nsub
                              call msr(avg%blk(ic2,ib)%m11m11(jc3,jl0r,il0,jsub,isub))
                              call msr(avg%blk(ic2,ib)%m2m2(jc3,jl0r,il0,jsub,isub))
                           end do
                           if (.not.nam%gau_approx) call msr(avg%blk(ic2,ib)%m22(jc3,jl0r,il0,isub))
                        end do
                     end if
                     if (avg%blk(ic2,ib)%nc1a_cor(jc3,jl0r,il0)>0.0) then
                        avg%blk(ic2,ib)%cor(jc3,jl0r,il0) = avg%blk(ic2,ib)%cor(jc3,jl0r,il0) &
                                                          & /avg%blk(ic2,ib)%nc1a_cor(jc3,jl0r,il0)
                     else
                        call msr(avg%blk(ic2,ib)%cor(jc3,jl0r,il0))
                     end if
                  end do
               end do
             end do
             !$omp end parallel do
         end if
      end do
   end if
end do

end subroutine avg_normalize

!----------------------------------------------------------------------
! Subroutine: avg_gather_lr
! Purpose: gather low-resolution averaged statistics data
!----------------------------------------------------------------------
subroutine avg_gather_lr(avg_lr,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_lr ! Averaged statistics, low resolution
type(mpl_type),intent(in) :: mpl        ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: npack,offset,ib,ic2
real(kind_real),allocatable :: sbuf(:),rbuf(:)
logical,allocatable :: mask_unpack(:,:,:,:,:)

! Packing size
npack = 0
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         if ((ic2==0).or.(nam%local_diag)) npack = npack+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_lr%blk(ic2,ib)%nsub**2
      end do
   end if
end do

! Allocation
allocate(sbuf(npack))
allocate(rbuf(npack))

! Pack data
offset = 0
sbuf = 0.0
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         if ((ic2==0).or.(nam%local_diag)) then
            sbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_lr%blk(ic2,ib)%nsub**2) = &
          & pack(avg_lr%blk(ic2,ib)%m11lrm11sub,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_lr%blk(ic2,ib)%nsub**2
         end if
      end do
   end if
end do

! Reduce data
call mpl%f_comm%allreduce(sbuf,rbuf,fckit_mpi_sum())

! Unpack data
offset = 0
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      ! Allocation
      allocate(mask_unpack(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg_lr%nsub,avg_lr%nsub))

      ! Initialization
      mask_unpack = .true.

      do ic2=0,nam%nc2
         ! Unpack
         if ((ic2==0).or.(nam%local_diag)) then
            avg_lr%blk(ic2,ib)%m11lrm11sub = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)* &
                                           & geom%nl0*avg_lr%blk(ic2,ib)%nsub**2),mask_unpack,avg_lr%blk(ic2,ib)%m11lrm11sub)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg_lr%blk(ic2,ib)%nsub**2
         end if
      end do

      ! Release memory
      deallocate(mask_unpack)
   end if
end do

end subroutine avg_gather_lr

!----------------------------------------------------------------------
! Subroutine: avg_normalize_lr
! Purpose: normalize low-resolution averaged statistics data
!----------------------------------------------------------------------
subroutine avg_normalize_lr(avg_lr,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_lr ! Averaged statistics, low resolution
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: ib,ic2,il0,jl0r,jc3,isub,jsub
real(kind_real) :: norm

! Normalize
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,nam%nc2
         if ((ic2==0).or.(nam%local_diag)) then
            !$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub,norm)
            do il0=1,geom%nl0
               do jl0r=1,bpar%nl0r(ib)
                  do jc3=1,bpar%nc3(ib)
                     if (avg_lr%blk(ic2,ib)%nc1a(jc3,jl0r,il0)>0.0) then
                        norm = 1.0/avg_lr%blk(ic2,ib)%nc1a(jc3,jl0r,il0)
                        do isub=1,avg_lr%blk(ic2,ib)%nsub
                           do jsub=1,avg_lr%blk(ic2,ib)%nsub
                              avg_lr%blk(ic2,ib)%m11lrm11sub(jc3,jl0r,il0,jsub,isub) &
                            & = avg_lr%blk(ic2,ib)%m11lrm11sub(jc3,jl0r,il0,jsub,isub)*norm
                            end do
                        end do
                     else
                        do isub=1,avg_lr%blk(ic2,ib)%nsub
                           do jsub=1,avg_lr%blk(ic2,ib)%nsub
                              call msr(avg_lr%blk(ic2,ib)%m11lrm11sub(jc3,jl0r,il0,jsub,isub))
                           end do
                        end do
                     end if
                  end do
               end do
             end do
             !$omp end parallel do
         end if
      end do
   end if
end do

end subroutine avg_normalize_lr

!----------------------------------------------------------------------
! Subroutine: avg_var_filter
! Purpose: filter variance
!----------------------------------------------------------------------
subroutine avg_var_filter(avg,mpl,nam,geom,bpar,samp)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(samp_type),intent(in) :: samp   ! Sampling

! Local variables
integer :: n,ib,il0,ic2a,ic2,iter
real(kind_real) :: P9,P20,P21
real(kind_real) :: m2sq,m4,m2sqasy,rhflt,drhflt
real(kind_real) :: m2_ini(samp%nc2a),m2(samp%nc2a),m2prod,m2prod_tot
logical :: dichotomy,convergence

! Ensemble/sub-ensemble size-dependent coefficients
n = avg%ne/avg%nsub
P9 = -real(n,kind_real)/real((n-2)*(n-3),kind_real)
P20 = real((n-1)*(n**2-3*n+3),kind_real)/real(n*(n-2)*(n-3),kind_real)
P21 = real(n-1,kind_real)/real(n+1,kind_real)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%info)

      do il0=1,geom%nl0
         write(mpl%info,'(a16,a,i3,a)') '','Level ',nam%levs(il0),':'
         call flush(mpl%info)

         ! Global averages
         m2sq = 0.0
         m4 = 0.0
         do ic2=1,nam%nc2
            m2sq = m2sq+sum(avg%blk(ic2,ib)%m2(il0,:)**2)/real(avg%nsub,kind_real)
            if (.not.nam%gau_approx) m4 = m4+sum(avg%blk(ic2,ib)%m4(il0,:))/real(avg%nsub,kind_real)
         end do

         ! Asymptotic statistics
         if (nam%gau_approx) then
            ! Gaussian approximation
            m2sqasy = P21*m2sq
         else
            ! General case
            m2sqasy = P20*m2sq+P9*m4
         end if

         ! Dichotomy initialization
         do ic2a=1,samp%nc2a
            ic2 = samp%c2a_to_c2(ic2a)
            m2_ini(ic2a) = sum(avg%blk(ic2,ib)%m2(il0,:))/real(avg%nsub,kind_real)
         end do
         convergence = .true.
         dichotomy = .false.
         rhflt = nam%var_rhflt
         drhflt = rhflt

         do iter=1,nam%var_niter
            ! Copy initial value
            m2 = m2_ini

            ! Median filter to remove extreme values
            call samp%diag_filter(mpl,nam,geom,il0,'median',rhflt,m2)

            ! Average filter to smooth displacement
            call samp%diag_filter(mpl,nam,geom,il0,'gc99',rhflt,m2)

            ! Compute product
            m2prod = sum(m2*m2_ini)

            ! Reduce product
            call mpl%f_comm%allreduce(m2prod,m2prod_tot,fckit_mpi_sum())

            ! Print result
            write(mpl%info,'(a19,a,i2,a,f10.2,a,f10.2)') '','Iteration ',iter,': rhflt = ', &
          & rhflt*reqkm,' km, difference = ',m2prod_tot-m2sqasy
            call flush(mpl%info)

            ! Update support radius
            if (m2prod_tot>m2sqasy) then
               ! Increase filtering support radius
               if (dichotomy) then
                  drhflt = 0.5*drhflt
                  rhflt = rhflt+drhflt
               else
                  convergence = .false.
                  rhflt = rhflt+drhflt
                  drhflt = 2.0*drhflt
               end if
            else
               ! Convergence
               convergence = .true.

               ! Change dichotomy status
               if (.not.dichotomy) then
                  dichotomy = .true.
                  drhflt = 0.5*drhflt
               end if

               ! Decrease filtering support radius
               drhflt = 0.5*drhflt
               rhflt = rhflt-drhflt
            end if
         end do

         ! Copy final result
         avg%blk(0,ib)%m2flt(il0) = sum(avg%blk(0,ib)%m2(il0,:))/real(avg%nsub,kind_real)
         do ic2a=1,samp%nc2a
            ic2 = samp%c2a_to_c2(ic2a)
            avg%blk(ic2,ib)%m2flt(il0) = m2(ic2a)
         end do
      end do
   end if
end do

end subroutine avg_var_filter

!----------------------------------------------------------------------
! Subroutine: avg_compute
! Purpose: compute averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute(avg,mpl,nam,geom,bpar,samp,mom,ne)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(samp_type),intent(in) :: samp   ! Sampling
type(mom_type),intent(in) :: mom     ! Moments
integer,intent(in) :: ne             ! Ensemble size

! Local variables
integer :: ib,ic2

! Allocation
call avg%alloc(nam,geom,bpar,mom%ne,mom%nsub)

! Compute averaged statistics
write(mpl%info,'(a10,a)') '','Compute averaged statistics'
call flush(mpl%info)
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%info)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) call avg%blk(ic2,ib)%compute(nam,geom,bpar,samp,mom%blk(ib))
         call mpl%prog_print(ic2+1)
      end do
      write(mpl%info,'(a)') '100%'
      call flush(mpl%info)
   end if
end do

if (mpl%nproc>1) then
   ! Gather averaged statistics
   write(mpl%info,'(a10,a)') '','Gather averaged statistics'
   call flush(mpl%info)
   call avg%gather(mpl,nam,geom,bpar)
end if

! Normalize averaged statistics
write(mpl%info,'(a10,a)') '','Normalize averaged statistics'
call flush(mpl%info)
call avg%normalize(nam,geom,bpar)

if (nam%var_filter) then
   ! Filter variance
   write(mpl%info,'(a10,a)') '','Filter variance'
   call flush(mpl%info)
   call avg%var_filter(mpl,nam,geom,bpar,samp)
end if

! Compute asymptotic statistics
write(mpl%info,'(a10,a)') '','Compute asymptotic statistics:'
call flush(mpl%info)
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%info)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) call avg%blk(ic2,ib)%compute_asy(nam,geom,bpar,ne)
         call mpl%prog_print(ic2+1)
      end do
      write(mpl%info,'(a)') '100%'
      call flush(mpl%info)
   end if
end do

end subroutine avg_compute

!----------------------------------------------------------------------
! Subroutine: avg_compute_hyb
! Purpose: compute hybrid averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_hyb(avg_2,mpl,nam,geom,bpar,samp,mom_1,mom_2,avg_1)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_2 ! Ensemble 2 averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(samp_type),intent(in) :: samp     ! Sampling
type(mom_type),intent(in) :: mom_1     ! Ensemble 2 moments
type(mom_type),intent(in) :: mom_2     ! Ensemble 1 moments
class(avg_type),intent(inout) :: avg_1 ! Ensemble 1 averaged statistics

! Local variables
integer :: ib,ic2

! Allocation
if (.not.allocated(avg_2%blk)) call avg_2%alloc(nam,geom,bpar,avg_2%ne,avg_2%nsub)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%info)

      ! Compute averaged statistics
      write(mpl%info,'(a13,a)',advance='no') '','Compute averaged statistics:'
      call flush(mpl%info)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) then
            select case (trim(nam%method))
            case ('hyb-avg')
               ! Static covariance = ensemble covariance
               avg_2%blk(ic2,ib)%m11sta = avg_1%blk(ic2,ib)%m11**2
               avg_2%blk(ic2,ib)%stasq = avg_1%blk(ic2,ib)%m11**2
            case ('hyb-rnd')
               ! Static covariance = randomized covariance
               avg_2%blk(ic2,ib)%m11sta = avg_1%blk(ic2,ib)%m11*avg_2%blk(ic2,ib)%m11
               avg_2%blk(ic2,ib)%stasq = avg_2%blk(ic2,ib)%m11**2
            case ('dual-ens')
               ! LR covariance/HR covariance product average
               call avg_2%blk(ic2,ib)%compute_lr(mpl,nam,geom,bpar,samp,mom_1%blk(ib),mom_2%blk(ib))
            end select
         end if
         call mpl%prog_print(ic2+1)
      end do
      write(mpl%info,'(a)') '100%'
      call flush(mpl%info)
   end if
end do

if (trim(nam%method)=='dual-ens') then
   if (mpl%nproc>1) then
      ! Gather averaged statistics
      write(mpl%info,'(a13,a)') '','Gather averaged statistics'
      call flush(mpl%info)
      call avg_2%gather_lr(mpl,nam,geom,bpar)
    end if

   ! Normalize averaged statistics
   write(mpl%info,'(a10,a)') '','Normalize averaged statistics'
   call flush(mpl%info)
   call avg_2%normalize_lr(nam,geom,bpar)

   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         ! Compute asymptotic statistics
         write(mpl%info,'(a13,a)',advance='no') '','Compute asymptotic statistics:'
         call flush(mpl%info)
         call mpl%prog_init(nam%nc2+1)
         do ic2=0,nam%nc2
            if ((ic2==0).or.nam%local_diag) call avg_2%blk(ic2,ib)%compute_asy_lr(nam,geom,bpar)
            call mpl%prog_print(ic2+1)
         end do
         write(mpl%info,'(a)') '100%'
         call flush(mpl%info)
      end if
   end do
end if

end subroutine avg_compute_hyb

!----------------------------------------------------------------------
! Function: avg_copy_wgt
! Purpose: averaged statistics data copy for weight definition
!----------------------------------------------------------------------
type(avg_type) function avg_copy_wgt(avg,geom,bpar)

implicit none

! Passed variables
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
class(avg_type),intent(inout) :: avg ! Averaged statistics

! Local variables
integer :: ib

if (bpar%diag_block(bpar%nbe)) then
   ! Allocation
   allocate(avg_copy_wgt%blk(0:0,bpar%nb))

   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         ! Restricted allocation
         allocate(avg_copy_wgt%blk(0,ib)%m2m2asy(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))

         ! Restricted copy
         avg_copy_wgt%blk(0,ib)%m2m2asy = avg%blk(0,ib)%m2m2asy
      end if
   end do
end if

end function avg_copy_wgt

!----------------------------------------------------------------------
! Subroutine: avg_compute_bwavg
! Purpose: compute block-averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_bwavg(avg,mpl,nam,geom,bpar,avg_wgt)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(avg_type),intent(in) :: avg_wgt ! Averaged statistics for weights

! Local variables
integer :: ib,ic2,il0,jl0r,jc3
real(kind_real) :: bwgtsq
real(kind_real),allocatable :: cor(:,:,:),m11asysq(:,:,:),m11sq(:,:,:)
real(kind_real),allocatable :: m11sta(:,:,:),stasq(:,:,:)
real(kind_real),allocatable :: m11lrm11(:,:,:),m11lrm11asy(:,:,:)

! Allocation
allocate(cor(nam%nc3,bpar%nl0rmax,geom%nl0))
allocate(m11asysq(nam%nc3,bpar%nl0rmax,geom%nl0))
allocate(m11sq(nam%nc3,bpar%nl0rmax,geom%nl0))
select case (trim(nam%method))
case ('hyb-avg','hyb-rnd')
   allocate(m11sta(nam%nc3,bpar%nl0rmax,geom%nl0))
   allocate(stasq(nam%nc3,bpar%nl0rmax,geom%nl0))
case ('dual-ens')
   allocate(m11lrm11(nam%nc3,bpar%nl0rmax,geom%nl0))
   allocate(m11lrm11asy(nam%nc3,bpar%nl0rmax,geom%nl0))
end select

write(mpl%info,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call flush(mpl%info)

! Initialization
call mpl%prog_init(nam%nc2+1)

do ic2=0,nam%nc2
   ! Copy ensemble size
   avg%blk(ic2,bpar%nbe)%ne = avg%blk(ic2,1)%ne
   avg%blk(ic2,bpar%nbe)%nsub = avg%blk(ic2,1)%nsub

   if ((ic2==0).or.(nam%local_diag)) then
      ! Initialization
      avg%blk(ic2,bpar%nbe)%cor = 0.0
      cor = 0.0
      avg%blk(ic2,bpar%nbe)%m11asysq = 0.0
      m11asysq = 0.0
      avg%blk(ic2,bpar%nbe)%m11sq = 0.0
      m11sq = 0.0
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd')
         avg%blk(ic2,bpar%nbe)%m11sta = 0.0
         m11sta = 0.0
         avg%blk(ic2,bpar%nbe)%stasq = 0.0
         stasq = 0.0
      case ('dual-ens')
         avg%blk(ic2,bpar%nbe)%m11lrm11 = 0.0
         m11lrm11 = 0.0
         avg%blk(ic2,bpar%nbe)%m11lrm11asy = 0.0
         m11lrm11asy = 0.0
      end select

      ! Block averages
      do ib=1,bpar%nb
         if (bpar%avg_block(ib)) then
            !$omp parallel do schedule(static) private(il0,jl0r,bwgtsq,jc3)
            do il0=1,geom%nl0
               do jl0r=1,bpar%nl0r(ib)
                  ! Weight
                  if (avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)>0.0) then
                     bwgtsq = 1.0/avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)
                  else
                     bwgtsq = 0.0
                  end if

                  ! Compute sum
                  do jc3=1,nam%nc3
                     call add(avg%blk(ic2,ib)%cor(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
                     call add(avg%blk(ic2,ib)%m11asysq(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11asysq(jc3,jl0r,il0), &
                   & m11asysq(jc3,jl0r,il0),bwgtsq)
                     call add(avg%blk(ic2,ib)%m11sq(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11sq(jc3,jl0r,il0), &
                   & m11sq(jc3,jl0r,il0),bwgtsq)
                     select case (trim(nam%method))
                     case ('hyb-avg','hyb-rnd')
                        call add(avg%blk(ic2,ib)%m11sta(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11sta(jc3,jl0r,il0), &
                      & m11sta(jc3,jl0r,il0),bwgtsq)
                        call add(avg%blk(ic2,ib)%stasq(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%stasq(jc3,jl0r,il0), &
                      & stasq(jc3,jl0r,il0),bwgtsq)
                     case ('dual-ens')
                        call add(avg%blk(ic2,ib)%m11lrm11(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11lrm11(jc3,jl0r,il0), &
                      & m11lrm11(jc3,jl0r,il0),bwgtsq)
                        call add(avg%blk(ic2,ib)%m11lrm11asy(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11lrm11asy(jc3,jl0r,il0), &
                      & m11lrm11asy(jc3,jl0r,il0),bwgtsq)
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
         do jl0r=1,bpar%nl0r(ib)
            do jc3=1,nam%nc3
               call divide(avg%blk(ic2,bpar%nbe)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
               call divide(avg%blk(ic2,bpar%nbe)%m11asysq(jc3,jl0r,il0),m11asysq(jc3,jl0r,il0))
               call divide(avg%blk(ic2,bpar%nbe)%m11sq(jc3,jl0r,il0),m11sq(jc3,jl0r,il0))
               select case (trim(nam%method))
               case ('hyb-avg','hyb-rnd')
                  call divide(avg%blk(ic2,bpar%nbe)%m11sta(jc3,jl0r,il0),m11sta(jc3,jl0r,il0))
                  call divide(avg%blk(ic2,bpar%nbe)%stasq(jc3,jl0r,il0),stasq(jc3,jl0r,il0))
               case ('dual-ens')
                  call divide(avg%blk(ic2,bpar%nbe)%m11lrm11(jc3,jl0r,il0),m11lrm11(jc3,jl0r,il0))
                  call divide(avg%blk(ic2,bpar%nbe)%m11lrm11asy(jc3,jl0r,il0),m11lrm11asy(jc3,jl0r,il0))
               end select
            end do
         end do
      end do
      !$omp end parallel do
   end if

   ! Update
   call mpl%prog_print(ic2+1)
end do
write(mpl%info,'(a)') '100%'
call flush(mpl%info)

end subroutine avg_compute_bwavg

end module type_avg
