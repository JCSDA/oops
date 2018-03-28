!----------------------------------------------------------------------
! Module: type_avg.f90
!> Purpose: average routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_avg

use omp_lib
use tools_display, only: prog_init,prog_print,peach,black
use tools_func, only: add,divide
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_qsort, only: qsort
use type_avg_blk, only: avg_blk_type
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_mom, only: mom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! Averaged statistics derived type
type avg_type
   integer :: nc2a                            !< Number of local points
   integer :: ne                              !< Ensemble size
   integer :: nsub                            !< Number of sub-ensembles
   type(avg_blk_type),allocatable :: blk(:,:) !< Averaged statistics blocks
contains
   procedure :: alloc => avg_alloc
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
integer,intent(in) :: ne !<
integer,intent(in) :: nsub !<

! Local variables
integer :: ib,ic2a

! Set attributes
if (nam%local_diag) then
   avg%nc2a = hdata%nc2a
else
   avg%nc2a = 0
end if
avg%ne = ne
avg%nsub = nsub

! Allocation
allocate(avg%blk(0:avg%nc2a,bpar%nb+1))
do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      do ic2a=0,avg%nc2a
         call avg%blk(ic2a,ib)%alloc(nam,geom,bpar,hdata,ic2a,ib,ne,nsub)
      end do
   end if
end do

end subroutine avg_alloc

!----------------------------------------------------------------------
! Subroutine: avg_compute
!> Purpose: compute averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute(avg,nam,geom,bpar,hdata,mom,ne)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data
type(mom_type),intent(in) :: mom     !< Moments
integer,intent(in) :: ne             !< Ensemble size

! Local variables
integer :: ib,ic2a,progint
logical,allocatable :: done(:)

! Allocation
call avg%alloc(nam,geom,bpar,hdata,mom%ne,mom%nsub)
allocate(done(0:avg%nc2a))

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'

      ! Initialization
      call prog_init(progint,done)

      do ic2a=0,avg%nc2a
         ! Compute averaged statistics
         call avg%blk(ic2a,ib)%compute(nam,geom,bpar,hdata,mom%blk(ib),ne)

         ! Update
         done(ic2a) = .true.
         call prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
   end if
end do

end subroutine avg_compute

!----------------------------------------------------------------------
! Subroutine: avg_compute_hyb
!> Purpose: compute hybrid averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_hyb(avg_2,nam,geom,bpar,hdata,mom_1,mom_2,avg_1)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_2 !< Ensemble 2 averaged statistics
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(mom_type),intent(in) :: mom_1     !< Ensemble 2 moments
type(mom_type),intent(in) :: mom_2     !< Ensemble 1 moments
class(avg_type),intent(inout) :: avg_1 !< Ensemble 1 averaged statistics

! Local variables
integer :: ib,ic2a,progint
logical,allocatable :: done(:)

! Allocation
if (.not.allocated(avg_2%blk)) call avg_2%alloc(nam,geom,bpar,hdata,avg_1%ne,avg_1%nsub)
allocate(done(0:avg_2%nc2a))

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'

      ! Initialization
      call prog_init(progint,done)

      do ic2a=0,avg_2%nc2a
         select case (trim(nam%method))
         case ('hyb-avg')
            ! Static covariance = ensemble covariance
            avg_2%blk(ic2a,ib)%m11sta = avg_1%blk(ic2a,ib)%m11*avg_1%blk(ic2a,ib)%m11
            avg_2%blk(ic2a,ib)%stasq = avg_1%blk(ic2a,ib)%m11**2
         case ('hyb-rnd')
            ! Static covariance = randomized covariance
            avg_2%blk(ic2a,ib)%m11sta = avg_1%blk(ic2a,ib)%m11*avg_2%blk(ic2a,ib)%m11
            avg_2%blk(ic2a,ib)%stasq = avg_2%blk(ic2a,ib)%m11**2
         case ('dual-ens')
            ! LR covariance/HR covariance product average
            call avg_2%blk(ic2a,ib)%compute_lr(geom,bpar,hdata,mom_1%blk(ib),mom_2%blk(ib),avg_1%blk(ic2a,ib))
         end select

         ! Update
         done(ic2a) = .true.
         call prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
   end if
end do

end subroutine avg_compute_hyb

!----------------------------------------------------------------------
! Function: avg_copy_wgt
!> Purpose: averaged statistics object copy for weight definition
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
subroutine avg_compute_bwavg(avg,nam,geom,bpar,avg_wgt)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg !< Averaged statistics
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(avg_type),intent(in) :: avg_wgt !< Averaged statistics for weights

! Local variables
integer :: ib,ic2a,il0,jl0r,jc3,progint
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
allocate(done(0:avg%nc2a))

write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(bpar%nb+1)),':'

! Initialization
call prog_init(progint,done)

do ic2a=0,avg%nc2a
   ! Copy ensemble size
   avg%blk(ic2a,bpar%nb+1)%ne = avg%blk(ic2a,1)%ne
   avg%blk(ic2a,bpar%nb+1)%nsub = avg%blk(ic2a,1)%nsub

   ! Initialization
   avg%blk(ic2a,bpar%nb+1)%cor = 0.0
   cor = 0.0
   avg%blk(ic2a,bpar%nb+1)%m11asysq = 0.0
   m11asysq = 0.0
   avg%blk(ic2a,bpar%nb+1)%m11sq = 0.0
   m11sq = 0.0
   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      avg%blk(ic2a,bpar%nb+1)%m11sta = 0.0
      m11sta = 0.0
      avg%blk(ic2a,bpar%nb+1)%stasq = 0.0
      stasq = 0.0
   case ('dual-ens')
      avg%blk(ic2a,bpar%nb+1)%m11lrm11 = 0.0
      m11lrm11 = 0.0
      avg%blk(ic2a,bpar%nb+1)%m11lrm11asy = 0.0
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
                  call add(avg%blk(ic2a,ib)%cor(jc3,jl0r,il0),avg%blk(ic2a,bpar%nb+1)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
                  call add(avg%blk(ic2a,ib)%m11asysq(jc3,jl0r,il0),avg%blk(ic2a,bpar%nb+1)%m11asysq(jc3,jl0r,il0), &
                & m11asysq(jc3,jl0r,il0),bwgtsq)
                  call add(avg%blk(ic2a,ib)%m11sq(jc3,jl0r,il0),avg%blk(ic2a,bpar%nb+1)%m11sq(jc3,jl0r,il0), &
                & m11sq(jc3,jl0r,il0),bwgtsq)
                  select case (trim(nam%method))
                  case ('hyb-avg','hyb-rnd')
                     call add(avg%blk(ic2a,ib)%m11sta(jc3,jl0r,il0),avg%blk(ic2a,bpar%nb+1)%m11sta(jc3,jl0r,il0), &
                   & m11sta(jc3,jl0r,il0),bwgtsq)
                     call add(avg%blk(ic2a,ib)%stasq(jc3,jl0r,il0),avg%blk(ic2a,bpar%nb+1)%stasq(jc3,jl0r,il0), &
                   & stasq(jc3,jl0r,il0),bwgtsq)
                  case ('dual-ens')
                     call add(avg%blk(ic2a,ib)%m11lrm11(jc3,jl0r,il0),avg%blk(ic2a,bpar%nb+1)%m11lrm11(jc3,jl0r,il0), &
                   & m11lrm11(jc3,jl0r,il0),bwgtsq)
                     call add(avg%blk(ic2a,ib)%m11lrm11asy(jc3,jl0r,il0),avg%blk(ic2a,bpar%nb+1)%m11lrm11asy(jc3,jl0r,il0), &
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
            call divide(avg%blk(ic2a,bpar%nb+1)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
            call divide(avg%blk(ic2a,bpar%nb+1)%m11asysq(jc3,jl0r,il0),m11asysq(jc3,jl0r,il0))
            call divide(avg%blk(ic2a,bpar%nb+1)%m11sq(jc3,jl0r,il0),m11sq(jc3,jl0r,il0))
            select case (trim(nam%method))
            case ('hyb-avg','hyb-rnd')
               call divide(avg%blk(ic2a,bpar%nb+1)%m11sta(jc3,jl0r,il0),m11sta(jc3,jl0r,il0))
               call divide(avg%blk(ic2a,bpar%nb+1)%stasq(jc3,jl0r,il0),stasq(jc3,jl0r,il0))
            case ('dual-ens')
               call divide(avg%blk(ic2a,bpar%nb+1)%m11lrm11(jc3,jl0r,il0),m11lrm11(jc3,jl0r,il0))
               call divide(avg%blk(ic2a,bpar%nb+1)%m11lrm11asy(jc3,jl0r,il0),m11lrm11asy(jc3,jl0r,il0))
            end select
         end do
      end do
   end do
   !$omp end parallel do

   ! Update
   done(ic2a) = .true.
   call prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'

end subroutine avg_compute_bwavg

end module type_avg
