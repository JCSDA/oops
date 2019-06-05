!----------------------------------------------------------------------
! Module: type_hdiag
! Purpose: hybrid diagnostics derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_hdiag

use tools_kinds, only: kind_real
use type_adv, only: adv_type
use type_avg, only: avg_type
use type_bpar, only: bpar_type
use type_diag, only: diag_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_mom, only: mom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type
use type_samp, only: samp_type

implicit none

! Hybrid diagnostics derived type
type hdiag_type
   type(avg_type) :: avg_1   ! Averaged statistics, first ensemble
   type(avg_type) :: avg_2   ! Averaged statistics, second ensemble
   type(avg_type) :: avg_wgt ! Averaged statistics weights
   type(samp_type) :: samp   ! Sampling
   type(adv_type) :: adv     ! Advection
   type(mom_type) :: mom_1   ! Moments, first ensemble
   type(mom_type) :: mom_2   ! Moments, second ensemble
   type(diag_type) :: cov_1  ! Covariance, first ensemble
   type(diag_type) :: cov_2  ! Covariance, second ensemble
   type(diag_type) :: cor_1  ! Correlation, first ensemble
   type(diag_type) :: cor_2  ! Correlation, second ensemble
   type(diag_type) :: loc_1  ! Localization, first ensemble
   type(diag_type) :: loc_2  ! Localization, second ensemble
   type(diag_type) :: loc_3  ! Localization, low-resolution ensemble
contains
   procedure :: dealloc => hdiag_dealloc
   procedure :: run_hdiag => hdiag_run_hdiag
end type hdiag_type

private
public :: hdiag_type

contains

!----------------------------------------------------------------------
! Subroutine: hdiag_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine hdiag_dealloc(hdiag)

implicit none

! Passed variables
class(hdiag_type),intent(inout) :: hdiag ! Hybrid diagnostics

! Release memory
call hdiag%avg_1%dealloc
call hdiag%avg_2%dealloc
call hdiag%avg_wgt%dealloc
call hdiag%samp%dealloc
call hdiag%adv%dealloc
call hdiag%mom_1%dealloc
call hdiag%mom_2%dealloc
call hdiag%cov_1%dealloc
call hdiag%cov_2%dealloc
call hdiag%cor_1%dealloc
call hdiag%cor_2%dealloc
call hdiag%loc_1%dealloc
call hdiag%loc_2%dealloc
call hdiag%loc_3%dealloc

end subroutine hdiag_dealloc

!----------------------------------------------------------------------
! Subroutine: hdiag_run_hdiag
! Purpose: HDIAG driver
!----------------------------------------------------------------------
subroutine hdiag_run_hdiag(hdiag,mpl,rng,nam,geom,bpar,io,ens1,ens2)

implicit none

! Passed variables
class(hdiag_type),intent(inout) :: hdiag   ! Hybrid diagnostics
type(mpl_type),intent(inout) :: mpl        ! MPI data
type(rng_type),intent(inout) :: rng        ! Random number generator
type(nam_type),intent(inout) :: nam        ! Namelist
type(geom_type),intent(in) :: geom         ! Geometry
type(bpar_type),intent(in) :: bpar         ! Block parameters
type(io_type),intent(in) :: io             ! I/O
type(ens_type),intent(in) :: ens1          ! Ensemble 1
type(ens_type),intent(in),optional :: ens2 ! Ensemble 2

! Setup sampling
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'
call mpl%flush
call hdiag%samp%setup_sampling(mpl,rng,nam,geom,bpar,io,ens1)

! Compute MPI distribution, halo A
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute MPI distribution, halos A'
call mpl%flush
call hdiag%samp%compute_mpi_a(mpl,nam,geom)

if (nam%new_lct.or.nam%local_diag.or.nam%adv_diag) then
   ! Compute MPI distribution, halos A-B
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute MPI distribution, halos A-B'
   call mpl%flush
   call hdiag%samp%compute_mpi_ab(mpl,nam,geom)
end if

if (nam%adv_diag) then
   ! Compute advection diagnostic
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute advection diagnostic'
   call mpl%flush
   call hdiag%adv%compute(mpl,rng,nam,geom,hdiag%samp,ens1)
end if

! Compute MPI distribution, halo C
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute MPI distribution, halo C'
call mpl%flush
call hdiag%samp%compute_mpi_c(mpl,nam,geom)

if ((nam%local_diag.or.nam%adv_diag).and.(nam%diag_rhflt>0.0)) then
   ! Compute MPI distribution, halo F
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute MPI distribution, halo F'
   call mpl%flush
   call hdiag%samp%compute_mpi_f(mpl,nam)
end if

if (nam%new_mom) then
   ! Compute sample moments
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute sample moments'
   call mpl%flush

   ! Compute ensemble 1 sample moments
   write(mpl%info,'(a7,a)') '','Ensemble 1:'
   call mpl%flush
   call hdiag%mom_1%compute(mpl,nam,geom,bpar,hdiag%samp,ens1,'mom_1')

   select case(trim(nam%method))
   case ('hyb-rnd','dual-ens')
      ! Compute ensemble 2 sample moments
      write(mpl%info,'(a7,a)') '','Ensemble 2:'
      call mpl%flush
      call hdiag%mom_2%compute(mpl,nam,geom,bpar,hdiag%samp,ens2,'mom_2')
   end select
elseif (nam%load_mom) then
   ! Load sample moments
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Load sample moments'
   call mpl%flush

   ! Load ensemble 1 sample moments
   write(mpl%info,'(a7,a)') '','Ensemble 1'
   call mpl%flush
   call hdiag%mom_1%read(mpl,nam,geom,bpar,hdiag%samp,ens1,'mom_1')

   select case(trim(nam%method))
   case ('hyb-rnd','dual-ens')
      ! Load ensemble 2 sample moments
      write(mpl%info,'(a7,a)') '','Ensemble 2'
      call mpl%flush
      call hdiag%mom_2%read(mpl,nam,geom,bpar,hdiag%samp,ens2,'mom_2')
   end select
end if

! Compute statistics
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute statistics'
call mpl%flush

! Compute ensemble 1 statistics
write(mpl%info,'(a7,a)') '','Ensemble 1:'
call mpl%flush
call hdiag%avg_1%compute(mpl,nam,geom,bpar,hdiag%samp,hdiag%mom_1,nam%ne,'avg_1')

select case(trim(nam%method))
case ('hyb-rnd','dual-ens')
   ! Compute ensemble 2 statistics
   write(mpl%info,'(a7,a)') '','Ensemble 2:'
   call mpl%flush
   call hdiag%avg_2%compute(mpl,nam,geom,bpar,hdiag%samp,hdiag%mom_2,nam%ens2_ne,'avg_2')
case ('hyb-avg')
   ! Copy ensemble 1 statistics
   call hdiag%avg_2%alloc(nam,geom,bpar,nam%ens2_ne,nam%ens2_nsub,'avg_2')
   call hdiag%avg_2%copy(hdiag%avg_1)
end select

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd','dual-ens')
   ! Compute cross-ensembles statistics
   write(mpl%info,'(a7,a)') '','Cross-ensembles:'
   call mpl%flush
   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      call hdiag%avg_1%compute_hyb(mpl,nam,geom,bpar,hdiag%avg_2)
   case ('dual-ens')
      call hdiag%avg_1%compute_deh(mpl,nam,geom,bpar,hdiag%samp,hdiag%mom_1,hdiag%mom_2)
   end select
end select

if ((bpar%nbe>bpar%nb).and.bpar%diag_block(bpar%nbe)) then
   ! Compute block-averaged statistics
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute block-averaged statistics'
   call mpl%flush

   ! Copy weight
   write(mpl%info,'(a7,a)') '','Copy weight'
   call mpl%flush
   call hdiag%avg_wgt%copy_wgt(geom,bpar,hdiag%avg_1)

   ! Compute ensemble 1 block-averaged statistics
   write(mpl%info,'(a7,a)') '','Ensemble 1:'
   call mpl%flush
   call hdiag%avg_1%compute_bwavg(mpl,nam,geom,bpar,hdiag%avg_wgt)

   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd','dual-ens')
      ! Compute ensemble 2 block-averaged statistics
      write(mpl%info,'(a7,a)') '','Ensemble 2:'
      call mpl%flush
      call hdiag%avg_2%compute_bwavg(mpl,nam,geom,bpar,hdiag%avg_wgt)

      ! Compute cross-ensembles block-averaged statistics
      write(mpl%info,'(a7,a)') '','Cross-ensembles:'
      call mpl%flush
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd')
         call hdiag%avg_1%compute_bwavg_hyb(mpl,nam,geom,bpar,hdiag%avg_wgt)
      case ('dual-ens')
         call hdiag%avg_1%compute_bwavg_deh(mpl,nam,geom,bpar,hdiag%avg_wgt)
      end select
   end select
end if

write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute covariance'
call mpl%flush

! Compute ensemble 1 covariance
write(mpl%info,'(a7,a)') '','Ensemble 1:'
call mpl%flush
call hdiag%cov_1%covariance(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_1,'cov')

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd','dual-ens')
   ! Compute ensemble 2 covariance
   write(mpl%info,'(a7,a)') '','Ensemble 2:'
   call mpl%flush
   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      call hdiag%cov_2%covariance(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_2,'cov_sta')
   case ('dual-ens')
      call hdiag%cov_2%covariance(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_2,'cov_lr')
   end select
end select

write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute correlation'
call mpl%flush

! Compute ensemble 1 correlation
write(mpl%info,'(a7,a)') '','Ensemble 1:'
call mpl%flush
call hdiag%cor_1%correlation(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_1,'cor')

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd','dual-ens')
   ! Compute ensemble 2 correlation
   write(mpl%info,'(a7,a)') '','Ensemble 2:'
   call mpl%flush
   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      call hdiag%cor_2%correlation(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_2,'cor_sta')
   case ('dual-ens')
      call hdiag%cor_2%correlation(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_2,'cor_lr')
   end select
end select

select case (trim(nam%method))
case ('loc','hyb-avg','hyb-rnd','dual-ens')
   ! Compute localization
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute localization'
   call mpl%flush
   write(mpl%info,'(a7,a)') '','Ensemble 1:'
   call mpl%flush
   call hdiag%loc_1%localization(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_1,'loc')
end select

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd')
   ! Compute static hybridization
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute static hybridization'
   call mpl%flush
   write(mpl%info,'(a7,a)') '','Ensemble 1 and 2:'
   call mpl%flush
   call hdiag%loc_2%hybridization(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_1,'loc_hyb')
end select

if (trim(nam%method)=='dual-ens') then
   ! Compute dual-ensemble hybridization diagnostic and fit
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute dual-ensemble hybridization'
   call mpl%flush
   write(mpl%info,'(a7,a)') '','Ensembles 1 and 2:'
   call mpl%flush
   call hdiag%loc_2%dualens(mpl,nam,geom,bpar,io,hdiag%samp,hdiag%avg_1,hdiag%avg_2,hdiag%loc_3,'loc_deh','loc_deh_lr')
end if

if (nam%write_hdiag) then
   ! Write data
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Write data'
   call mpl%flush

   ! Advection
   if (nam%adv_diag) call hdiag%adv%write(mpl,nam,geom,bpar,io,hdiag%samp)
end if

end subroutine hdiag_run_hdiag

end module type_hdiag
