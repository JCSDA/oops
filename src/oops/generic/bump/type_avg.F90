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
   character(len=1024) :: prefix              ! Prefix
   integer :: ne                              ! Ensemble size
   integer :: nsub                            ! Number of sub-ensembles
   type(avg_blk_type),allocatable :: blk(:,:) ! Averaged statistics blocks
contains
   procedure :: alloc => avg_alloc
   procedure :: dealloc => avg_dealloc
   procedure :: copy => avg_copy
   procedure :: write => avg_write
   procedure :: gather => avg_gather
   procedure :: normalize => avg_normalize
   procedure :: var_filter => avg_var_filter
   procedure :: compute => avg_compute
   procedure :: compute_hyb => avg_compute_hyb
   procedure :: compute_deh => avg_compute_deh
   procedure :: copy_wgt => avg_copy_wgt
   procedure :: compute_bwavg => avg_compute_bwavg
   procedure :: compute_bwavg_hyb => avg_compute_bwavg_hyb
   procedure :: compute_bwavg_deh => avg_compute_bwavg_deh
end type avg_type

private
public :: avg_type

contains

!----------------------------------------------------------------------
! Subroutine: avg_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine avg_alloc(avg,nam,geom,bpar,ne,nsub,prefix)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg  ! Averaged statistics
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Block parameters
integer,intent(in) :: ne              ! Ensemble size
integer,intent(in) :: nsub            ! Number of sub-ensembles
character(len=*),intent(in) :: prefix ! Prefix

! Local variables
integer :: ib,ic2

! Set attributes
avg%prefix = trim(prefix)
avg%ne = ne
avg%nsub = nsub

! Allocation
allocate(avg%blk(0:nam%nc2,bpar%nbe))
do ib=1,bpar%nbe
   do ic2=0,nam%nc2
      call avg%blk(ic2,ib)%alloc(nam,geom,bpar,ic2,ib,ne,nsub,prefix)
   end do
end do

end subroutine avg_alloc

!----------------------------------------------------------------------
! Subroutine: avg_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine avg_dealloc(avg)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics

! Local variables
integer :: ib,ic2

! Allocation
if (allocated(avg%blk)) then
   do ib=1,size(avg%blk,2)
      do ic2=0,size(avg%blk,1)-1
         call avg%blk(ic2,ib)%dealloc
      end do
    end do
   deallocate(avg%blk)
end if

end subroutine avg_dealloc

!----------------------------------------------------------------------
! Subroutine: avg_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine avg_copy(avg_out,avg_in)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_out ! Output averaged statistics
class(avg_type),intent(in) :: avg_in     ! Input averaged statistics

! Local variables
integer :: ib,ic2

! Copy
do ib=1,size(avg_in%blk,2)
   do ic2=0,size(avg_in%blk,1)-1
      call avg_out%blk(ic2,ib)%copy(avg_in%blk(ic2,ib))
   end do
end do

end subroutine avg_copy

!----------------------------------------------------------------------
! Subroutine: avg_write
! Purpose: write
!----------------------------------------------------------------------
subroutine avg_write(avg,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Diagnostic
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib
character(len=1024) :: filename

if (mpl%main) then
   filename = trim(nam%prefix)//'_avg'
   do ib=1,bpar%nbe
      if (bpar%diag_block(ib)) call avg%blk(0,ib)%write(mpl,nam,geom,bpar,filename)
   end do
end if

end subroutine avg_write

!----------------------------------------------------------------------
! Subroutine: avg_gather
! Purpose: gather averaged statistics data
!----------------------------------------------------------------------
subroutine avg_gather(avg,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
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
         if ((ic2==0).or.nam%local_diag) then
            npack = npack+geom%nl0*avg%blk(ic2,ib)%nsub
            if (nam%var_filter.and.(.not.nam%gau_approx)) npack = npack+geom%nl0*avg%blk(ic2,ib)%nsub
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
         if ((ic2==0).or.nam%local_diag) then
            sbuf(offset+1:offset+geom%nl0*avg%blk(ic2,ib)%nsub) = pack(avg%blk(ic2,ib)%m2,.true.)
            offset = offset+geom%nl0*avg%blk(ic2,ib)%nsub
            if (nam%var_filter.and.(.not.nam%gau_approx)) then
              sbuf(offset+1:offset+geom%nl0*avg%blk(ic2,ib)%nsub) = pack(avg%blk(ic2,ib)%m4,.true.)
               offset = offset+geom%nl0*avg%blk(ic2,ib)%nsub
            end if
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
         if ((ic2==0).or.nam%local_diag) then
            avg%blk(ic2,ib)%m2 = unpack(rbuf(offset+1:offset+geom%nl0*avg%blk(ic2,ib)%nsub),mask_1(1,1,:,:),avg%blk(ic2,ib)%m2)
            offset = offset+geom%nl0*avg%blk(ic2,ib)%nsub
            if (nam%var_filter.and.(.not.nam%gau_approx)) then
               avg%blk(ic2,ib)%m4 = unpack(rbuf(offset+1:offset+geom%nl0*avg%blk(ic2,ib)%nsub),mask_1(1,1,:,:),avg%blk(ic2,ib)%m4)
               offset = offset+geom%nl0*avg%blk(ic2,ib)%nsub
            end if
            avg%blk(ic2,ib)%nc1a = unpack(rbuf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%nc1a)
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
subroutine avg_normalize(avg,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
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
         if ((ic2==0).or.nam%local_diag) then
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
                        avg%blk(ic2,ib)%m11(jc3,jl0r,il0) = mpl%msv%valr
                        do isub=1,avg%blk(ic2,ib)%nsub
                           do jsub=1,avg%blk(ic2,ib)%nsub
                              avg%blk(ic2,ib)%m11m11(jc3,jl0r,il0,jsub,isub) = mpl%msv%valr
                              avg%blk(ic2,ib)%m2m2(jc3,jl0r,il0,jsub,isub) = mpl%msv%valr
                           end do
                           if (.not.nam%gau_approx) avg%blk(ic2,ib)%m22(jc3,jl0r,il0,isub) = mpl%msv%valr
                        end do
                     end if
                     if (avg%blk(ic2,ib)%nc1a_cor(jc3,jl0r,il0)>0.0) then
                        avg%blk(ic2,ib)%cor(jc3,jl0r,il0) = avg%blk(ic2,ib)%cor(jc3,jl0r,il0) &
                                                          & /avg%blk(ic2,ib)%nc1a_cor(jc3,jl0r,il0)
                     else
                        avg%blk(ic2,ib)%cor(jc3,jl0r,il0) = mpl%msv%valr
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
integer :: n,ib,il0,isub,ic2a,ic2,iter
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
      call mpl%flush

      do il0=1,geom%nl0
         write(mpl%info,'(a16,a,i3,a)') '','Level ',nam%levs(il0),':'
         call mpl%flush

         do isub=1,avg%nsub
            ! Global sum
            m2sq = 0.0
            m4 = 0.0
            do ic2=1,nam%nc2
               if (mpl%msv%isnotr(avg%blk(ic2,ib)%m2(il0,isub))) then
                  m2sq = m2sq+avg%blk(ic2,ib)%m2(il0,isub)**2
                  if (.not.nam%gau_approx) m4 = m4+avg%blk(ic2,ib)%m4(il0,isub)
               end if
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
               m2_ini(ic2a) = avg%blk(ic2,ib)%m2(il0,isub)
            end do
            convergence = .true.
            dichotomy = .false.
            rhflt = nam%var_rhflt
            drhflt = rhflt

            do iter=1,nam%var_niter
               ! Copy initial value
               m2 = m2_ini

               ! Median filter to remove extreme values
               call samp%diag_filter(mpl,nam,'median',rhflt,m2)

               ! Average filter to smooth values
               call samp%diag_filter(mpl,nam,'gc99',rhflt,m2)

               ! Compute product
               m2prod = 0.0
               do ic2a=1,samp%nc2a
                  if (mpl%msv%isnotr(m2_ini(ic2a))) m2prod = m2prod+m2(ic2a)*m2_ini(ic2a)
               end do

               ! Reduce product
               call mpl%f_comm%allreduce(m2prod,m2prod_tot,fckit_mpi_sum())

               ! Print result
               write(mpl%info,'(a19,a,i2,a,f10.2,a,e12.5)') '','Iteration ',iter,': rhflt = ', &
             & rhflt*reqkm,' km, rel. diff. = ',(m2prod_tot-m2sqasy)/m2sqasy
               call mpl%flush

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
            avg%blk(0,ib)%m2flt(il0,isub) = avg%blk(0,ib)%m2(il0,isub)
            do ic2a=1,samp%nc2a
               ic2 = samp%c2a_to_c2(ic2a)
               avg%blk(ic2,ib)%m2flt(il0,isub) = m2(ic2a)
            end do
         end do
      end do
   end if
end do

end subroutine avg_var_filter

!----------------------------------------------------------------------
! Subroutine: avg_compute
! Purpose: compute averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute(avg,mpl,nam,geom,bpar,samp,mom,ne,prefix)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg  ! Averaged statistics
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Block parameters
type(samp_type),intent(in) :: samp    ! Sampling
type(mom_type),intent(in) :: mom      ! Moments
integer,intent(in) :: ne              ! Ensemble size
character(len=*),intent(in) :: prefix ! Prefix

! Local variables
integer :: ib,ic2

! Allocation
call avg%alloc(nam,geom,bpar,mom%ne,mom%nsub,prefix)

! Compute averaged statistics
write(mpl%info,'(a10,a)') '','Compute averaged statistics'
call mpl%flush
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) call avg%blk(ic2,ib)%compute(mpl,nam,geom,bpar,samp,mom%blk(ib))
         call mpl%prog_print(ic2+1)
      end do
      call mpl%prog_final
   end if
end do

if (mpl%main.and.(nam%avg_nbins>0)) then
   ! Write histograms
   write(mpl%info,'(a10,a)') '','Write histograms'
   call mpl%flush
   call avg%write(mpl,nam,geom,bpar)
end if

if (mpl%nproc>1) then
   ! Gather averaged statistics
   write(mpl%info,'(a10,a)') '','Gather averaged statistics'
   call mpl%flush
   call avg%gather(mpl,nam,geom,bpar)
end if

! Normalize averaged statistics
write(mpl%info,'(a10,a)') '','Normalize averaged statistics'
call mpl%flush
call avg%normalize(mpl,nam,geom,bpar)

if (nam%var_filter) then
   ! Filter variance
   write(mpl%info,'(a10,a)') '','Filter variance'
   call mpl%flush
   call avg%var_filter(mpl,nam,geom,bpar,samp)
end if

! Compute asymptotic statistics
write(mpl%info,'(a10,a)') '','Compute asymptotic statistics:'
call mpl%flush
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) call avg%blk(ic2,ib)%compute_asy(mpl,nam,geom,bpar,ne)
         call mpl%prog_print(ic2+1)
      end do
      call mpl%prog_final
   end if
end do

end subroutine avg_compute

!----------------------------------------------------------------------
! Subroutine: avg_compute_hyb
! Purpose: compute hybrid averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_hyb(avg_1,mpl,nam,geom,bpar,avg_2)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_1 ! Ensemble 1 averaged statistics
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(avg_type),intent(in) :: avg_2     ! Ensemble 2 averaged statistics

! Local variables
integer :: ib,ic2

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush

      ! Compute averaged statistics
      write(mpl%info,'(a13,a)') '','Compute averaged statistics:'
      call mpl%flush(.false.)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) then
            select case (trim(nam%method))
            case ('hyb-avg')
               ! Static covariance = ensemble covariance
               call avg_1%blk(ic2,ib)%compute_hyb(mpl,geom,bpar,avg_1%blk(ic2,ib))
            case ('hyb-rnd')
               ! Static covariance = randomized covariance
               call avg_1%blk(ic2,ib)%compute_hyb(mpl,geom,bpar,avg_2%blk(ic2,ib))
            end select
         end if
         call mpl%prog_print(ic2+1)
      end do
      call mpl%prog_final
   end if
end do

end subroutine avg_compute_hyb

!----------------------------------------------------------------------
! Subroutine: avg_compute_deh
! Purpose: compute dual-ensemble hybrid averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_deh(avg_1,mpl,nam,geom,bpar,samp,mom_1,mom_2)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_1 ! Ensemble 1 averaged statistics
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(samp_type),intent(in) :: samp     ! Sampling
type(mom_type),intent(in) :: mom_1     ! Ensemble 1 moments
type(mom_type),intent(in) :: mom_2     ! Ensemble 2 moments

! Local variables
integer :: ib,ic2

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush

      ! Compute averaged statistics
      write(mpl%info,'(a13,a)') '','Compute averaged statistics:'
      call mpl%flush(.false.)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) call avg_1%blk(ic2,ib)%compute_deh(mpl,nam,geom,bpar,samp,mom_1%blk(ib),mom_2%blk(ib))
         call mpl%prog_print(ic2+1)
      end do
      call mpl%prog_final
   end if
end do

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      ! Compute asymptotic statistics
      write(mpl%info,'(a13,a)') '','Compute asymptotic statistics at low resolution:'
      call mpl%flush(.false.)
      call mpl%prog_init(nam%nc2+1)
      do ic2=0,nam%nc2
         if ((ic2==0).or.nam%local_diag) call avg_1%blk(ic2,ib)%compute_asy_deh(mpl,nam,geom,bpar)
         call mpl%prog_print(ic2+1)
      end do
      call mpl%prog_final
   end if
end do

end subroutine avg_compute_deh

!----------------------------------------------------------------------
! Subroutine: avg_copy_wgt
! Purpose: averaged statistics data copy for weight definition
!----------------------------------------------------------------------
subroutine avg_copy_wgt(avg_out,geom,bpar,avg_in)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_out ! Output averaged statistics
type(geom_type),intent(in) :: geom       ! Geometry
type(bpar_type),intent(in) :: bpar       ! Block parameters
type(avg_type),intent(in) :: avg_in      ! Input averaged statistics

! Local variables
integer :: ib

if (bpar%diag_block(bpar%nbe)) then
   ! Allocation
   allocate(avg_out%blk(0:0,bpar%nb))

   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         ! Restricted allocation
         allocate(avg_out%blk(0,ib)%m2m2asy(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))

         ! Restricted copy
         avg_out%blk(0,ib)%m2m2asy = avg_in%blk(0,ib)%m2m2asy
      end if
   end do
end if

end subroutine avg_copy_wgt

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
real(kind_real) :: cor(nam%nc3,bpar%nl0rmax,geom%nl0),nc1a_cor(nam%nc3,bpar%nl0rmax,geom%nl0)
real(kind_real) :: m11asysq(nam%nc3,bpar%nl0rmax,geom%nl0),m11sq(nam%nc3,bpar%nl0rmax,geom%nl0),nc1a(nam%nc3,bpar%nl0rmax,geom%nl0)

! Initialization
write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call mpl%flush(.false.)
call mpl%prog_init(nam%nc2+1)

do ic2=0,nam%nc2
   ! Copy ensemble size
   avg%blk(ic2,bpar%nbe)%ne = avg%blk(ic2,1)%ne
   avg%blk(ic2,bpar%nbe)%nsub = avg%blk(ic2,1)%nsub

   if ((ic2==0).or.nam%local_diag) then
      ! Initialization
      if (nam%var_filter) then
         avg%blk(ic2,bpar%nbe)%m2flt = 1.0
      else
         avg%blk(ic2,bpar%nbe)%m2 = 1.0
      end if
      avg%blk(ic2,bpar%nbe)%cor = 0.0
      cor = 0.0
      avg%blk(ic2,bpar%nbe)%nc1a_cor = 0.0
      nc1a_cor = 0.0
      avg%blk(ic2,bpar%nbe)%m11asysq = 0.0
      m11asysq = 0.0
      avg%blk(ic2,bpar%nbe)%m11sq = 0.0
      m11sq = 0.0
      avg%blk(ic2,bpar%nbe)%nc1a = 0.0
      nc1a = 0.0

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
                     call add(mpl,avg%blk(ic2,ib)%cor(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
                     call add(mpl,avg%blk(ic2,ib)%nc1a_cor(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%nc1a_cor(jc3,jl0r,il0), &
                   & nc1a_cor(jc3,jl0r,il0))
                     call add(mpl,avg%blk(ic2,ib)%m11asysq(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11asysq(jc3,jl0r,il0), &
                   & m11asysq(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2,ib)%m11sq(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11sq(jc3,jl0r,il0), &
                   & m11sq(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2,ib)%nc1a(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%nc1a(jc3,jl0r,il0), &
                   & nc1a(jc3,jl0r,il0),bwgtsq)
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
               call divide(mpl,avg%blk(ic2,bpar%nbe)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2,bpar%nbe)%nc1a_cor(jc3,jl0r,il0),nc1a_cor(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2,bpar%nbe)%m11asysq(jc3,jl0r,il0),m11asysq(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2,bpar%nbe)%m11sq(jc3,jl0r,il0),m11sq(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2,bpar%nbe)%nc1a(jc3,jl0r,il0),nc1a(jc3,jl0r,il0))
            end do
         end do
      end do
      !$omp end parallel do
   end if

   ! Update
   call mpl%prog_print(ic2+1)
end do
call mpl%prog_final

end subroutine avg_compute_bwavg

!----------------------------------------------------------------------
! Subroutine: avg_compute_bwavg_hyb
! Purpose: compute hybrid block-averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_bwavg_hyb(avg,mpl,nam,geom,bpar,avg_wgt)

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
real(kind_real) :: m11sta(nam%nc3,bpar%nl0rmax,geom%nl0),stasq(nam%nc3,bpar%nl0rmax,geom%nl0)

! Initialization
write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call mpl%flush(.false.)
call mpl%prog_init(nam%nc2+1)

do ic2=0,nam%nc2
   if ((ic2==0).or.nam%local_diag) then
      ! Initialization
      avg%blk(ic2,bpar%nbe)%m11sta = 0.0
      m11sta = 0.0
      avg%blk(ic2,bpar%nbe)%stasq = 0.0
      stasq = 0.0

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
                     call add(mpl,avg%blk(ic2,ib)%m11sta(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11sta(jc3,jl0r,il0), &
                   & m11sta(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2,ib)%stasq(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%stasq(jc3,jl0r,il0), &
                   & stasq(jc3,jl0r,il0),bwgtsq)
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
               call divide(mpl,avg%blk(ic2,bpar%nbe)%m11sta(jc3,jl0r,il0),m11sta(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2,bpar%nbe)%stasq(jc3,jl0r,il0),stasq(jc3,jl0r,il0))
            end do
         end do
      end do
      !$omp end parallel do
   end if

   ! Update
   call mpl%prog_print(ic2+1)
end do
call mpl%prog_final

end subroutine avg_compute_bwavg_hyb

!----------------------------------------------------------------------
! Subroutine: avg_compute_bwavg_deh
! Purpose: compute dual-ensemble hybrid block-averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_bwavg_deh(avg,mpl,nam,geom,bpar,avg_wgt)

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
real(kind_real) :: m11lrm11(nam%nc3,bpar%nl0rmax,geom%nl0),m11lrm11asy(nam%nc3,bpar%nl0rmax,geom%nl0)

! Initialization
write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call mpl%flush(.false.)
call mpl%prog_init(nam%nc2+1)

do ic2=0,nam%nc2
   if ((ic2==0).or.nam%local_diag) then
      ! Initialization
      avg%blk(ic2,bpar%nbe)%m11lrm11 = 0.0
      m11lrm11 = 0.0
      avg%blk(ic2,bpar%nbe)%m11lrm11asy = 0.0
      m11lrm11asy = 0.0

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
                     call add(mpl,avg%blk(ic2,ib)%m11lrm11(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11lrm11(jc3,jl0r,il0), &
                   & m11lrm11(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2,ib)%m11lrm11asy(jc3,jl0r,il0),avg%blk(ic2,bpar%nbe)%m11lrm11asy(jc3,jl0r,il0), &
                   & m11lrm11asy(jc3,jl0r,il0),bwgtsq)
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
               call divide(mpl,avg%blk(ic2,bpar%nbe)%m11lrm11(jc3,jl0r,il0),m11lrm11(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2,bpar%nbe)%m11lrm11asy(jc3,jl0r,il0),m11lrm11asy(jc3,jl0r,il0))
            end do
         end do
      end do
      !$omp end parallel do
   end if

   ! Update
   call mpl%prog_print(ic2+1)
end do
call mpl%prog_final

end subroutine avg_compute_bwavg_deh

end module type_avg
