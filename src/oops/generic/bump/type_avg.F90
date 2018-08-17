!----------------------------------------------------------------------
! Module: type_avg
!> Purpose: average routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_avg

!$ use omp_lib
use tools_const,only: reqkm
use tools_func, only: add,divide
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_qsort, only: qsort
use type_avg_blk, only: avg_blk_type
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_mom, only: mom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

! Averaged statistics derived type
type avg_type
   integer :: ne                              !< Ensemble size
   integer :: nsub                            !< Number of sub-ensembles
   type(avg_blk_type),allocatable :: blk(:,:) !< Averaged statistics blocks
contains
   procedure :: alloc => avg_alloc
   procedure :: dealloc => avg_dealloc
   procedure :: copy => avg_copy
   procedure :: gather => avg_gather
   procedure :: gather_lr => avg_gather_lr
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
!> Purpose: averaged statistics allocation
!----------------------------------------------------------------------
subroutine avg_alloc(avg,nam,geom,bpar,ne,nsub)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
integer,intent(in) :: ne             !< Ensemble size
integer,intent(in) :: nsub           !< Number of sub-ensembles

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
!> Purpose: averaged statistics deallocation
!----------------------------------------------------------------------
subroutine avg_dealloc(avg,nam,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(nam_type),intent(in) :: nam     !< Namelist
type(bpar_type),intent(in) :: bpar   !< Block parameters

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
!> Purpose: averaged statistics copy
!----------------------------------------------------------------------
type(avg_type) function avg_copy(avg,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(in) :: avg !< Averaged statistics
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters

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
!> Purpose: gather averaged statistics data
!----------------------------------------------------------------------
subroutine avg_gather(avg,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(mpl_type),intent(in) :: mpl     !< MPI data
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters

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
call mpl%allreduce_sum(sbuf,rbuf)

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

! Normalize for the average
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      avg%blk(0,ib)%m2 = avg%blk(0,ib)%m2/real(mpl%nproc,kind_real)
      if (nam%var_filter.and.(.not.nam%gau_approx)) avg%blk(0,ib)%m4 = avg%blk(0,ib)%m4/real(mpl%nproc,kind_real)
      avg%blk(0,ib)%nc1a = avg%blk(0,ib)%nc1a/real(mpl%nproc,kind_real)
      avg%blk(0,ib)%m11 = avg%blk(0,ib)%m11/real(mpl%nproc,kind_real)
      avg%blk(0,ib)%m11m11 = avg%blk(0,ib)%m11m11/real(mpl%nproc,kind_real)
      avg%blk(0,ib)%m2m2 = avg%blk(0,ib)%m2m2/real(mpl%nproc,kind_real)
      if (.not.nam%gau_approx) avg%blk(0,ib)%m22 = avg%blk(0,ib)%m22/real(mpl%nproc,kind_real)
      avg%blk(0,ib)%nc1a_cor = avg%blk(0,ib)%nc1a_cor/real(mpl%nproc,kind_real)
      avg%blk(0,ib)%cor = avg%blk(0,ib)%cor/real(mpl%nproc,kind_real)
   end if
end do

end subroutine avg_gather

!----------------------------------------------------------------------
! Subroutine: avg_gather_lr
!> Purpose: gather low-resolutin averaged statistics data
!----------------------------------------------------------------------
subroutine avg_gather_lr(avg_lr,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_lr !< Averaged statistics, low resolution
type(mpl_type),intent(in) :: mpl        !< MPI data
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters

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
call mpl%allreduce_sum(sbuf,rbuf)

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

! Normalize for the average
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      avg_lr%blk(0,ib)%m11lrm11sub = avg_lr%blk(0,ib)%m11lrm11sub/real(mpl%nproc,kind_real)
   end if
end do

end subroutine avg_gather_lr

!----------------------------------------------------------------------
! Subroutine: avg_var_filter
!> Purpose: filter variance
!----------------------------------------------------------------------
subroutine avg_var_filter(avg,mpl,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(mpl_type),intent(in) :: mpl     !< MPI data
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data

! Local variables
integer :: n,ib,il0,ic2a,ic2,iter
real(kind_real) :: P9,P20,P21
real(kind_real) :: m2sq,m4,m2sqasy,rhflt,drhflt
real(kind_real) :: m2_ini(hdata%nc2a),m2(hdata%nc2a),m2prod,m2prod_tot
logical :: dichotomy,convergence
 
! Ensemble/sub-ensemble size-dependent coefficients
n = avg%ne/avg%nsub
P9 = -real(n,kind_real)/real((n-2)*(n-3),kind_real)
P20 = real((n-1)*(n**2-3*n+3),kind_real)/real(n*(n-2)*(n-3),kind_real)
P21 = real(n-1,kind_real)/real(n+1,kind_real)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%unit)

      do il0=1,geom%nl0
         write(mpl%unit,'(a16,a,i3,a)') '','Level ',nam%levs(il0),':'
         call flush(mpl%unit)

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
         do ic2a=1,hdata%nc2a
            ic2 = hdata%c2a_to_c2(ic2a)
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
            call hdata%diag_filter(mpl,nam,geom,il0,'median',rhflt,m2)
 
            ! Average filter to smooth displacement
            call hdata%diag_filter(mpl,nam,geom,il0,'gc99',rhflt,m2)

            ! Compute product
            m2prod = sum(m2*m2_ini)

            ! Reduce product
            call mpl%allreduce_sum(m2prod,m2prod_tot)

            ! Print result
            write(mpl%unit,'(a19,a,i2,a,f10.2,a,f10.2)') '','Iteration ',iter,': rhflt = ', &
          & rhflt*reqkm,' km, difference = ',m2prod_tot-m2sqasy
            call flush(mpl%unit)

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
         do ic2a=1,hdata%nc2a
            ic2 = hdata%c2a_to_c2(ic2a)
            avg%blk(ic2,ib)%m2flt(il0) = m2(ic2a)
         end do
      end do
   end if
end do

end subroutine avg_var_filter

!----------------------------------------------------------------------
! Subroutine: avg_compute
!> Purpose: compute averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute(avg,mpl,nam,geom,bpar,hdata,mom,ne)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(mpl_type),intent(in) :: mpl     !< MPI data
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data
type(mom_type),intent(in) :: mom     !< Moments
integer,intent(in) :: ne             !< Ensemble size

! Local variables
integer :: ib,ic2,progint
logical,allocatable :: done(:)

! Allocation
call avg%alloc(nam,geom,bpar,mom%ne,mom%nsub)
allocate(done(0:nam%nc2))

! Compute averaged statistics
write(mpl%unit,'(a10,a)') '','Compute averaged statistics'
call flush(mpl%unit)
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a13,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%unit)
      call mpl%prog_init(progint,done)
      do ic2=0,nam%nc2
         call avg%blk(ic2,ib)%compute(nam,geom,bpar,hdata,mom%blk(ib))
         done(ic2) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)
   end if
end do

if (mpl%nproc>1) then
   ! Gather statistics
   write(mpl%unit,'(a10,a)') '','Gather data'
   call flush(mpl%unit)
   call avg%gather(mpl,nam,geom,bpar)
end if

if (nam%var_filter) then
   ! Filter variance
   write(mpl%unit,'(a10,a)') '','Filter variance'
   call flush(mpl%unit)
   call avg%var_filter(mpl,nam,geom,bpar,hdata)
end if

! Compute asymptotic statistics
write(mpl%unit,'(a10,a)') '','Compute asymptotic statistics:'
call flush(mpl%unit)
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a13,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%unit)
      call mpl%prog_init(progint,done)
      do ic2=0,nam%nc2
         call avg%blk(ic2,ib)%compute_asy(nam,geom,bpar,ne)
         done(ic2) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)
   end if
end do

end subroutine avg_compute

!----------------------------------------------------------------------
! Subroutine: avg_compute_hyb
!> Purpose: compute hybrid averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_hyb(avg_2,mpl,nam,geom,bpar,hdata,mom_1,mom_2,avg_1)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_2 !< Ensemble 2 averaged statistics
type(mpl_type),intent(in) :: mpl       !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(mom_type),intent(in) :: mom_1     !< Ensemble 2 moments
type(mom_type),intent(in) :: mom_2     !< Ensemble 1 moments
class(avg_type),intent(inout) :: avg_1 !< Ensemble 1 averaged statistics

! Local variables
integer :: ib,ic2,progint
logical,allocatable :: done(:)

! Allocation
if (.not.allocated(avg_2%blk)) call avg_2%alloc(nam,geom,bpar,avg_2%ne,avg_2%nsub)
allocate(done(0:nam%nc2))

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%unit)

      ! Compute averaged statistics
      write(mpl%unit,'(a13,a)',advance='no') '','Compute averaged statistics:'
      call flush(mpl%unit)
      call mpl%prog_init(progint,done)
      do ic2=0,nam%nc2
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
            call avg_2%blk(ic2,ib)%compute_lr(mpl,nam,geom,bpar,hdata,mom_1%blk(ib),mom_2%blk(ib),avg_1%blk(ic2,ib))
         end select
         done(ic2) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)
   end if
end do

if (trim(nam%method)=='dual-ens') then
   if (mpl%nproc>1) then
      ! Gather statistics
      write(mpl%unit,'(a13,a)') '','Gather data'
      call flush(mpl%unit)
      call avg_2%gather_lr(mpl,nam,geom,bpar)
    end if

   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         ! Compute asymptotic statistics
         write(mpl%unit,'(a13,a)',advance='no') '','Compute asymptotic statistics:'
         call flush(mpl%unit)
         call mpl%prog_init(progint,done)
         do ic2=0,nam%nc2
            call avg_2%blk(ic2,ib)%compute_asy_lr(nam,geom,bpar)
            done(ic2) = .true.
            call mpl%prog_print(progint,done)
         end do
         write(mpl%unit,'(a)') '100%'
         call flush(mpl%unit)
      end if
   end do
end if

end subroutine avg_compute_hyb

!----------------------------------------------------------------------
! Function: avg_copy_wgt
!> Purpose: averaged statistics data copy for weight definition
!----------------------------------------------------------------------
type(avg_type) function avg_copy_wgt(avg,geom,bpar)

implicit none

! Passed variables
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
class(avg_type),intent(inout) :: avg !< Averaged statistics

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
!> Purpose: compute block-averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_bwavg(avg,mpl,nam,geom,bpar,avg_wgt)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(mpl_type),intent(in) :: mpl     !< MPI data
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(avg_type),intent(in) :: avg_wgt !< Averaged statistics for weights

! Local variables
integer :: ib,ic2,il0,jl0r,jc3,progint
real(kind_real) :: bwgtsq
real(kind_real),allocatable :: cor(:,:,:),m11asysq(:,:,:),m11sq(:,:,:)
real(kind_real),allocatable :: m11sta(:,:,:),stasq(:,:,:)
real(kind_real),allocatable :: m11lrm11(:,:,:),m11lrm11asy(:,:,:)
logical,allocatable :: done(:)

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
allocate(done(0:nam%nc2))

write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call flush(mpl%unit)

! Initialization
call mpl%prog_init(progint,done)

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
               do jl0r=1,nam%nl0r
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
         do jl0r=1,nam%nl0r
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
   done(ic2) = .true.
   call mpl%prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'
call flush(mpl%unit)

end subroutine avg_compute_bwavg

end module type_avg
