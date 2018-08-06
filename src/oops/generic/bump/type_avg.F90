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
   integer :: nc2                             !< Number of local points
   integer :: ne                              !< Ensemble size
   integer :: nsub                            !< Number of sub-ensembles
   integer :: npack                           !< Pack format size
   type(avg_blk_type),allocatable :: blk(:,:) !< Averaged statistics blocks
contains
   procedure :: alloc => avg_alloc
   procedure :: pack => avg_pack
   procedure :: unpack => avg_unpack
   procedure :: gather => avg_gather
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
!> Purpose: allocation
!----------------------------------------------------------------------
subroutine avg_alloc(avg,nam,geom,bpar,hdata,ne,nsub)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data
integer,intent(in) :: ne             !< Ensemble size
integer,intent(in) :: nsub           !< Number of sub-ensembles

! Local variables
integer :: ib,ic2

! Set attributes
if (nam%local_diag) then
   avg%nc2 = hdata%nc2
else
   avg%nc2 = 0
end if
avg%ne = ne
avg%nsub = nsub

! Allocation
allocate(avg%blk(0:avg%nc2,bpar%nb+1))
do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      do ic2=0,avg%nc2
         call avg%blk(ic2,ib)%alloc(nam,geom,bpar,ic2,ib,ne,nsub)
      end do
   end if
end do

! Packing size
avg%npack = 0
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,avg%nc2
         avg%npack = avg%npack+(4+2*avg%blk(ic2,ib)%nsub**2)*bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         if (.not.nam%gau_approx) avg%npack = avg%npack+avg%blk(ic2,ib)%nsub*bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
      end do
   end if
end do

end subroutine avg_alloc

!----------------------------------------------------------------------
! Subroutine: avg_pack
!> Purpose: averaged statistics data packing
!----------------------------------------------------------------------
subroutine avg_pack(avg,nam,geom,bpar,buf)

implicit none

! Passed variables
class(avg_type),intent(in) :: avg             !< Averaged statistics block
type(nam_type),intent(in) :: nam              !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry
type(bpar_type),intent(in) :: bpar            !< Block parameters
real(kind_real),intent(out) :: buf(avg%npack) !< Buffer

! Local variables
integer :: offset,ib,ic2

! Initialization
offset = 0

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      do ic2=0,avg%nc2
         ! Pack
         buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%nc1a,.true.)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%m11,.true.)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2) = pack(avg%blk(ic2,ib)%m11m11,.true.)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
         buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2) = pack(avg%blk(ic2,ib)%m2m2,.true.)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
         if (.not.nam%gau_approx) then
            buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub) = pack(avg%blk(ic2,ib)%m22,.true.)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub
         end if
         buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%nc1a_cor,.true.)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0) = pack(avg%blk(ic2,ib)%cor,.true.)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
      end do
   end if
end do

end subroutine avg_pack

!----------------------------------------------------------------------
! Subroutine: avg_unpack
!> Purpose: averaged statistics data unpacking
!----------------------------------------------------------------------
subroutine avg_unpack(avg,nam,geom,bpar,buf)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg         !< Averaged statistics block
type(nam_type),intent(in) :: nam             !< Namelist
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
real(kind_real),intent(in) :: buf(avg%npack) !< Buffer

! Local variables
integer :: offset,ib,ic2
logical,allocatable :: mask_0(:,:,:),mask_1(:,:,:,:),mask_2(:,:,:,:,:)

! Initialization
offset = 0

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      ! Allocation
      allocate(mask_0(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      if (.not.nam%gau_approx) allocate(mask_1(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg%blk(0,ib)%nsub))
      allocate(mask_2(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,avg%blk(0,ib)%nsub,avg%blk(0,ib)%nsub))

      ! Initialization
      mask_0 = .true.
      if (.not.nam%gau_approx) mask_1 = .true.
      mask_2 = .true.

      do ic2=0,avg%nc2
         ! Unpack
         avg%blk(ic2,ib)%nc1a = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%m11)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         avg%blk(ic2,ib)%m11 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%m11)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         avg%blk(ic2,ib)%m11m11 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2), &
       & mask_2,avg%blk(ic2,ib)%m11m11)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
         avg%blk(ic2,ib)%m2m2 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2), &
       & mask_2,avg%blk(ic2,ib)%m2m2)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub**2
         if (.not.nam%gau_approx) then
            avg%blk(ic2,ib)%m22 = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub), &
          & mask_1,avg%blk(ic2,ib)%m22)
            offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0*avg%blk(ic2,ib)%nsub
         end if
         avg%blk(ic2,ib)%nc1a_cor = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%nc1a_cor)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
         avg%blk(ic2,ib)%cor = unpack(buf(offset+1:offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0),mask_0,avg%blk(ic2,ib)%cor)
         offset = offset+bpar%nc3(ib)*bpar%nl0r(ib)*geom%nl0
      end do

      ! Release memory
      deallocate(mask_0)
      if (.not.nam%gau_approx) deallocate(mask_1)
      deallocate(mask_2)
   end if
end do

end subroutine avg_unpack

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
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Allocation
allocate(sbuf(avg%npack))
allocate(rbuf(avg%npack))

! Pack data
call avg%pack(nam,geom,bpar,sbuf)

! Reduce data
call mpl%allreduce_sum(sbuf,rbuf)

! Unpack data
call avg%unpack(nam,geom,bpar,rbuf)

end subroutine avg_gather

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
call avg%alloc(nam,geom,bpar,hdata,mom%ne,mom%nsub)
allocate(done(0:avg%nc2))

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%unit)

      ! Compute averaged statistics
      write(mpl%unit,'(a13,a)',advance='no') '','Compute averaged statistics:'
      call mpl%prog_init(progint,done)
      do ic2=0,avg%nc2
         call avg%blk(ic2,ib)%compute(mpl,nam,geom,bpar,hdata,mom%blk(ib))
         done(ic2) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)

      if (mpl%nproc>1) then
         ! Gather statistics
         write(mpl%unit,'(a13,a)') '','Gather data'
        call avg%gather(mpl,nam,geom,bpar)
      end if

      ! Compute asymptotic statistics
      write(mpl%unit,'(a13,a)',advance='no') '','Compute asymptotic statistics:'
      call mpl%prog_init(progint,done)
      do ic2=0,avg%nc2
         call avg%blk(ic2,ib)%compute_asy(mpl,nam,geom,bpar,ne)
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
if (.not.allocated(avg_2%blk)) call avg_2%alloc(nam,geom,bpar,hdata,avg_1%ne,avg_1%nsub)
allocate(done(0:avg_2%nc2))

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%unit)

      ! Initialization
      call mpl%prog_init(progint,done)

      do ic2=0,avg_2%nc2
         select case (trim(nam%method))
         case ('hyb-avg')
            ! Static covariance = ensemble covariance
            avg_2%blk(ic2,ib)%m11sta = avg_1%blk(ic2,ib)%m11*avg_1%blk(ic2,ib)%m11
            avg_2%blk(ic2,ib)%stasq = avg_1%blk(ic2,ib)%m11**2
         case ('hyb-rnd')
            ! Static covariance = randomized covariance
            avg_2%blk(ic2,ib)%m11sta = avg_1%blk(ic2,ib)%m11*avg_2%blk(ic2,ib)%m11
            avg_2%blk(ic2,ib)%stasq = avg_2%blk(ic2,ib)%m11**2
         case ('dual-ens')
            ! LR covariance/HR covariance product average
            call avg_2%blk(ic2,ib)%compute_lr(mpl,geom,bpar,hdata,mom_1%blk(ib),mom_2%blk(ib),avg_1%blk(ic2,ib))
         end select

         ! Update
         done(ic2) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)
   end if
end do

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

if (bpar%diag_block(bpar%nb+1)) then
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
allocate(done(0:avg%nc2))

write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(bpar%nb+1)),':'
call flush(mpl%unit)

! Initialization
call mpl%prog_init(progint,done)

do ic2=0,avg%nc2
   ! Copy ensemble size
   avg%blk(ic2,bpar%nb+1)%ne = avg%blk(ic2,1)%ne
   avg%blk(ic2,bpar%nb+1)%nsub = avg%blk(ic2,1)%nsub

   ! Initialization
   avg%blk(ic2,bpar%nb+1)%cor = 0.0
   cor = 0.0
   avg%blk(ic2,bpar%nb+1)%m11asysq = 0.0
   m11asysq = 0.0
   avg%blk(ic2,bpar%nb+1)%m11sq = 0.0
   m11sq = 0.0
   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      avg%blk(ic2,bpar%nb+1)%m11sta = 0.0
      m11sta = 0.0
      avg%blk(ic2,bpar%nb+1)%stasq = 0.0
      stasq = 0.0
   case ('dual-ens')
      avg%blk(ic2,bpar%nb+1)%m11lrm11 = 0.0
      m11lrm11 = 0.0
      avg%blk(ic2,bpar%nb+1)%m11lrm11asy = 0.0
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
                  call add(avg%blk(ic2,ib)%cor(jc3,jl0r,il0),avg%blk(ic2,bpar%nb+1)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
                  call add(avg%blk(ic2,ib)%m11asysq(jc3,jl0r,il0),avg%blk(ic2,bpar%nb+1)%m11asysq(jc3,jl0r,il0), &
                & m11asysq(jc3,jl0r,il0),bwgtsq)
                  call add(avg%blk(ic2,ib)%m11sq(jc3,jl0r,il0),avg%blk(ic2,bpar%nb+1)%m11sq(jc3,jl0r,il0), &
                & m11sq(jc3,jl0r,il0),bwgtsq)
                  select case (trim(nam%method))
                  case ('hyb-avg','hyb-rnd')
                     call add(avg%blk(ic2,ib)%m11sta(jc3,jl0r,il0),avg%blk(ic2,bpar%nb+1)%m11sta(jc3,jl0r,il0), &
                   & m11sta(jc3,jl0r,il0),bwgtsq)
                     call add(avg%blk(ic2,ib)%stasq(jc3,jl0r,il0),avg%blk(ic2,bpar%nb+1)%stasq(jc3,jl0r,il0), &
                   & stasq(jc3,jl0r,il0),bwgtsq)
                  case ('dual-ens')
                     call add(avg%blk(ic2,ib)%m11lrm11(jc3,jl0r,il0),avg%blk(ic2,bpar%nb+1)%m11lrm11(jc3,jl0r,il0), &
                   & m11lrm11(jc3,jl0r,il0),bwgtsq)
                     call add(avg%blk(ic2,ib)%m11lrm11asy(jc3,jl0r,il0),avg%blk(ic2,bpar%nb+1)%m11lrm11asy(jc3,jl0r,il0), &
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
            call divide(avg%blk(ic2,bpar%nb+1)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
            call divide(avg%blk(ic2,bpar%nb+1)%m11asysq(jc3,jl0r,il0),m11asysq(jc3,jl0r,il0))
            call divide(avg%blk(ic2,bpar%nb+1)%m11sq(jc3,jl0r,il0),m11sq(jc3,jl0r,il0))
            select case (trim(nam%method))
            case ('hyb-avg','hyb-rnd')
               call divide(avg%blk(ic2,bpar%nb+1)%m11sta(jc3,jl0r,il0),m11sta(jc3,jl0r,il0))
               call divide(avg%blk(ic2,bpar%nb+1)%stasq(jc3,jl0r,il0),stasq(jc3,jl0r,il0))
            case ('dual-ens')
               call divide(avg%blk(ic2,bpar%nb+1)%m11lrm11(jc3,jl0r,il0),m11lrm11(jc3,jl0r,il0))
               call divide(avg%blk(ic2,bpar%nb+1)%m11lrm11asy(jc3,jl0r,il0),m11lrm11asy(jc3,jl0r,il0))
            end select
         end do
      end do
   end do
   !$omp end parallel do

   ! Update
   done(ic2) = .true.
   call mpl%prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'
call flush(mpl%unit)

end subroutine avg_compute_bwavg

end module type_avg
