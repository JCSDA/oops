!----------------------------------------------------------------------
! Module: type_avg_blk
! Purpose: averaged statistics block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_avg_blk

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max,fckit_mpi_status
use netcdf
use tools_func, only: histogram
use tools_kinds, only: kind_real,nc_kind_real
use tools_repro, only: sup,inf
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

! Averaged statistics block derived type
type avg_blk_type
   integer :: ic2                                        ! Global index
   integer :: ib                                         ! Block index
   character(len=1024) :: name                           ! Name
   integer :: ne                                         ! Ensemble size
   integer :: nsub                                       ! Sub-ensembles number
   real(kind_real),allocatable :: m2(:,:)                ! Variance
   real(kind_real),allocatable :: m4(:,:)                ! Fourth-order centered moment
   real(kind_real),allocatable :: m2flt(:,:)             ! Filtered variance
   real(kind_real),allocatable :: nc1a(:,:,:)            ! Number of points in subset Sc1 on halo A
   real(kind_real),allocatable :: m11(:,:,:)             ! Covariance average
   real(kind_real),allocatable :: m11m11(:,:,:,:,:)      ! Product of covariances average
   real(kind_real),allocatable :: m2m2(:,:,:,:,:)        ! Product of variances average
   real(kind_real),allocatable :: m22(:,:,:,:)           ! Fourth-order centered moment average
   real(kind_real),allocatable :: nc1a_cor(:,:,:)        ! Number of points in subset Sc1 on halo A with valid correlations
   real(kind_real),allocatable :: cor(:,:,:)             ! Correlation average
   real(kind_real),allocatable :: m11asysq(:,:,:)        ! Squared asymptotic covariance average
   real(kind_real),allocatable :: m2m2asy(:,:,:)         ! Product of asymptotic variances average
   real(kind_real),allocatable :: m22asy(:,:,:)          ! Asymptotic fourth-order centered moment average
   real(kind_real),allocatable :: m11sq(:,:,:)           ! Squared covariance average for several ensemble sizes
   real(kind_real),allocatable :: m11sta(:,:,:)          ! Ensemble covariance/static covariance product
   real(kind_real),allocatable :: stasq(:,:,:)           ! Squared static covariance
   real(kind_real),allocatable :: m11lrm11sub(:,:,:,:,:) ! LR covariance/HR covariance product average
   real(kind_real),allocatable :: m11lrm11(:,:,:)        ! LR covariance/HR covariance product average, averaged over sub-ensembles
   real(kind_real),allocatable :: m11lrm11asy(:,:,:)     ! LR covariance/HR asymptotic covariance product average
   real(kind_real),allocatable :: m11_bins(:,:,:,:)      ! Covariance histrogram bins
   real(kind_real),allocatable :: m11_hist(:,:,:,:)      ! Covariance histrogram values
   real(kind_real),allocatable :: m11m11_bins(:,:,:,:)   ! Product of covariances  histrogram bins
   real(kind_real),allocatable :: m11m11_hist(:,:,:,:)   ! Product of covariances  histrogram values
   real(kind_real),allocatable :: m2m2_bins(:,:,:,:)     ! Product of variances  histrogram bins
   real(kind_real),allocatable :: m2m2_hist(:,:,:,:)     ! Product of variances  histrogram values
   real(kind_real),allocatable :: m22_bins(:,:,:,:)      ! Fourth-order centered moment  histrogram bins
   real(kind_real),allocatable :: m22_hist(:,:,:,:)      ! Fourth-order centered moment  histrogram values
   real(kind_real),allocatable :: cor_bins(:,:,:,:)      ! Correlation histrogram bins
   real(kind_real),allocatable :: cor_hist(:,:,:,:)      ! Correlation histrogram values
contains
   procedure :: alloc => avg_blk_alloc
   procedure :: dealloc => avg_blk_dealloc
   procedure :: copy => avg_blk_copy
   procedure :: write => avg_blk_write
   procedure :: compute => avg_blk_compute
   procedure :: compute_asy => avg_blk_compute_asy
   procedure :: compute_hyb => avg_blk_compute_hyb
   procedure :: compute_deh => avg_blk_compute_deh
   procedure :: compute_asy_deh => avg_blk_compute_asy_deh
end type avg_blk_type

private
public :: avg_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: avg_blk_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine avg_blk_alloc(avg_blk,nam,geom,bpar,ic2,ib,ne,nsub,prefix)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk ! Averaged statistics block

type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
integer,intent(in) :: ic2                    ! Global index
integer,intent(in) :: ib                     ! Block index
integer,intent(in) :: ne                     ! Ensemble size
integer,intent(in) :: nsub                   ! Sub-ensembles number
character(len=*),intent(in) :: prefix        ! Prefix

! Set attributes
avg_blk%ic2 = ic2
avg_blk%ib = ib
avg_blk%ne = ne
avg_blk%nsub = nsub
avg_blk%name = trim(prefix)//'_'//trim(bpar%blockname(ib))

! Allocation
if (bpar%diag_block(ib).and.(.not.allocated(avg_blk%nc1a))) then
   if ((ic2==0).or.nam%local_diag) then
      allocate(avg_blk%nc1a(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m2(geom%nl0,avg_blk%nsub))
      if (nam%var_filter) then
         if (.not.nam%gau_approx) allocate(avg_blk%m4(geom%nl0,avg_blk%nsub))
         allocate(avg_blk%m2flt(geom%nl0,avg_blk%nsub))
      end if
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
   if ((nam%avg_nbins>0).and.(ic2==0)) then
      allocate(avg_blk%m11_bins(nam%avg_nbins+1,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11_hist(nam%avg_nbins,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11m11_bins(nam%avg_nbins+1,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m11m11_hist(nam%avg_nbins,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m2m2_bins(nam%avg_nbins+1,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m2m2_hist(nam%avg_nbins,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m22_bins(nam%avg_nbins+1,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%m22_hist(nam%avg_nbins,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%cor_bins(nam%avg_nbins+1,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
      allocate(avg_blk%cor_hist(nam%avg_nbins,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
   end if
end if

end subroutine avg_blk_alloc

!----------------------------------------------------------------------
! Subroutine: avg_blk_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine avg_blk_dealloc(avg_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk ! Averaged statistics block

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
if (allocated(avg_blk%m11_bins)) deallocate(avg_blk%m11_bins)
if (allocated(avg_blk%m11_hist)) deallocate(avg_blk%m11_hist)
if (allocated(avg_blk%m11m11_bins)) deallocate(avg_blk%m11m11_bins)
if (allocated(avg_blk%m11m11_hist)) deallocate(avg_blk%m11m11_hist)
if (allocated(avg_blk%m2m2_bins)) deallocate(avg_blk%m2m2_bins)
if (allocated(avg_blk%m2m2_hist)) deallocate(avg_blk%m2m2_hist)
if (allocated(avg_blk%m22_bins)) deallocate(avg_blk%m22_bins)
if (allocated(avg_blk%m22_hist)) deallocate(avg_blk%m22_hist)
if (allocated(avg_blk%cor_bins)) deallocate(avg_blk%cor_bins)
if (allocated(avg_blk%cor_hist)) deallocate(avg_blk%cor_hist)

end subroutine avg_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: avg_blk_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine avg_blk_copy(avg_blk_out,avg_blk_in)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk_out ! Output averaged statistics block
type(avg_blk_type),intent(in) :: avg_blk_in      ! Input averaged statistics block

! Copy data
if (allocated(avg_blk_in%m2)) avg_blk_out%m2 = avg_blk_in%m2
if (allocated(avg_blk_in%m4)) avg_blk_out%m4 = avg_blk_in%m4
if (allocated(avg_blk_in%m2flt)) avg_blk_out%m2flt = avg_blk_in%m2flt
if (allocated(avg_blk_in%nc1a)) avg_blk_out%nc1a = avg_blk_in%nc1a
if (allocated(avg_blk_in%m11)) avg_blk_out%m11 = avg_blk_in%m11
if (allocated(avg_blk_in%m11m11)) avg_blk_out%m11m11 = avg_blk_in%m11m11
if (allocated(avg_blk_in%m2m2)) avg_blk_out%m2m2 = avg_blk_in%m2m2
if (allocated(avg_blk_in%m22)) avg_blk_out%m22 = avg_blk_in%m22
if (allocated(avg_blk_in%nc1a_cor)) avg_blk_out%nc1a_cor = avg_blk_in%nc1a_cor
if (allocated(avg_blk_in%cor)) avg_blk_out%cor = avg_blk_in%cor
if (allocated(avg_blk_in%m11asysq)) avg_blk_out%m11asysq = avg_blk_in%m11asysq
if (allocated(avg_blk_in%m2m2asy)) avg_blk_out%m2m2asy = avg_blk_in%m2m2asy
if (allocated(avg_blk_in%m22asy)) avg_blk_out%m22asy = avg_blk_in%m22asy
if (allocated(avg_blk_in%m11sq)) avg_blk_out%m11sq = avg_blk_in%m11sq
if (allocated(avg_blk_in%m11sta)) avg_blk_out%m11sta = avg_blk_in%m11sta
if (allocated(avg_blk_in%stasq)) avg_blk_out%stasq = avg_blk_in%stasq
if (allocated(avg_blk_in%m11lrm11sub)) avg_blk_out%m11lrm11sub = avg_blk_in%m11lrm11sub
if (allocated(avg_blk_in%m11lrm11sub)) avg_blk_out%m11lrm11sub = avg_blk_in%m11lrm11sub
if (allocated(avg_blk_in%m11lrm11asy)) avg_blk_out%m11lrm11asy = avg_blk_in%m11lrm11asy
if (allocated(avg_blk_in%m11_bins)) avg_blk_out%m11_bins = avg_blk_in%m11_bins
if (allocated(avg_blk_in%m11_hist)) avg_blk_out%m11_hist = avg_blk_in%m11_hist
if (allocated(avg_blk_in%m11m11_bins)) avg_blk_out%m11m11_bins = avg_blk_in%m11m11_bins
if (allocated(avg_blk_in%m11m11_hist)) avg_blk_out%m11m11_hist = avg_blk_in%m11m11_hist
if (allocated(avg_blk_in%m2m2_bins)) avg_blk_out%m2m2_bins = avg_blk_in%m2m2_bins
if (allocated(avg_blk_in%m2m2_hist)) avg_blk_out%m2m2_hist = avg_blk_in%m2m2_hist
if (allocated(avg_blk_in%m22_bins)) avg_blk_out%m22_bins = avg_blk_in%m22_bins
if (allocated(avg_blk_in%m22_hist)) avg_blk_out%m22_hist = avg_blk_in%m22_hist
if (allocated(avg_blk_in%cor_bins)) avg_blk_out%cor_bins = avg_blk_in%cor_bins
if (allocated(avg_blk_in%cor_hist)) avg_blk_out%cor_hist = avg_blk_in%cor_hist

end subroutine avg_blk_copy

!----------------------------------------------------------------------
! Subroutine: avg_blk_write
! Purpose: write
!----------------------------------------------------------------------
subroutine avg_blk_write(avg_blk,mpl,nam,geom,bpar,filename)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
character(len=*),intent(in) :: filename      ! File name

! Local variables
integer :: info,info_coord,ncid,nc3_id,nl0r_id,nl0_id,nbinsp1_id,nbins_id,disth_id,vunit_id
integer :: m11_bins_id,m11_hist_id,m11m11_bins_id,m11m11_hist_id,m2m2_bins_id,m2m2_hist_id,m22_bins_id,m22_hist_id
integer :: cor_bins_id,cor_hist_id
character(len=1024),parameter :: subr = 'avg_blk_write'

! Associate
associate(ib=>avg_blk%ib)

! Check if the file exists
info = nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_noclobber,nf90_64bit_offset),ncid)
if (info==nf90_noerr) then
   ! Write namelist parameters
   call nam%write(mpl,ncid)
else
   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_write,ncid))

   ! Redef mode
   call mpl%ncerr(subr,nf90_redef(ncid))
end if

! Define dimensions and coordinates if necessary
nc3_id = mpl%ncdimcheck(subr,ncid,'nc3',nam%nc3,.true.)
nl0r_id = mpl%ncdimcheck(subr,ncid,'nl0r',bpar%nl0rmax,.true.)
nl0_id = mpl%ncdimcheck(subr,ncid,'nl0',geom%nl0,.true.)
nbinsp1_id = mpl%ncdimcheck(subr,ncid,'nbinsp1',nam%avg_nbins+1,.true.)
nbins_id = mpl%ncdimcheck(subr,ncid,'nbins',nam%avg_nbins,.true.)
info_coord = nf90_inq_varid(ncid,'disth',disth_id)
if (info_coord/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'disth',nc_kind_real,(/nc3_id/),disth_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'vunit',nc_kind_real,(/nl0_id/),vunit_id))
end if

! Define variables if necessary
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m11_bins',m11_bins_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m11_bins',nc_kind_real, &
 & (/nbinsp1_id,nc3_id,nl0r_id,nl0_id/),m11_bins_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m11_bins_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m11_hist',m11_hist_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m11_hist',nc_kind_real, &
 & (/nbins_id,nc3_id,nl0r_id,nl0_id/),m11_hist_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m11_hist_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m11m11_bins',m11m11_bins_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m11m11_bins',nc_kind_real, &
 & (/nbinsp1_id,nc3_id,nl0r_id,nl0_id/),m11m11_bins_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m11m11_bins_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m11m11_hist',m11m11_hist_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m11m11_hist',nc_kind_real, &
 & (/nbins_id,nc3_id,nl0r_id,nl0_id/),m11m11_hist_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m11m11_hist_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m2m2_bins',m2m2_bins_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m2m2_bins',nc_kind_real, &
 & (/nbinsp1_id,nc3_id,nl0r_id,nl0_id/),m2m2_bins_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m2m2_bins_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m2m2_hist',m2m2_hist_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m2m2_hist',nc_kind_real, &
 & (/nbins_id,nc3_id,nl0r_id,nl0_id/),m2m2_hist_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m2m2_hist_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m22_bins',m22_bins_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m22_bins',nc_kind_real, &
 & (/nbinsp1_id,nc3_id,nl0r_id,nl0_id/),m22_bins_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m22_bins_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_m22_hist',m22_hist_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_m22_hist',nc_kind_real, &
 & (/nbins_id,nc3_id,nl0r_id,nl0_id/),m22_hist_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,m22_hist_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_cor_bins',cor_bins_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_cor_bins',nc_kind_real, &
 & (/nbinsp1_id,nc3_id,nl0r_id,nl0_id/),cor_bins_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,cor_bins_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,trim(avg_blk%name)//'_cor_hist',cor_hist_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(avg_blk%name)//'_cor_hist',nc_kind_real, &
 & (/nbins_id,nc3_id,nl0r_id,nl0_id/),cor_hist_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,cor_hist_id,'_FillValue',mpl%msv%valr))
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write coordinates if necessary
if (info_coord/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_put_var(ncid,disth_id,geom%disth(1:nam%nc3)))
   call mpl%ncerr(subr,nf90_put_var(ncid,vunit_id,geom%vunitavg))
end if

! Write variables
call mpl%ncerr(subr,nf90_put_var(ncid,m11_bins_id,avg_blk%m11_bins))
call mpl%ncerr(subr,nf90_put_var(ncid,m11_hist_id,avg_blk%m11_hist))
call mpl%ncerr(subr,nf90_put_var(ncid,m11m11_bins_id,avg_blk%m11m11_bins))
call mpl%ncerr(subr,nf90_put_var(ncid,m11m11_hist_id,avg_blk%m11m11_hist))
call mpl%ncerr(subr,nf90_put_var(ncid,m2m2_bins_id,avg_blk%m2m2_bins))
call mpl%ncerr(subr,nf90_put_var(ncid,m2m2_hist_id,avg_blk%m2m2_hist))
call mpl%ncerr(subr,nf90_put_var(ncid,m22_bins_id,avg_blk%m22_bins))
call mpl%ncerr(subr,nf90_put_var(ncid,m22_hist_id,avg_blk%m22_hist))
call mpl%ncerr(subr,nf90_put_var(ncid,cor_bins_id,avg_blk%cor_bins))
call mpl%ncerr(subr,nf90_put_var(ncid,cor_hist_id,avg_blk%cor_hist))

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine avg_blk_write

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute
! Purpose: compute averaged statistics via spatial-angular erogodicity assumption
!----------------------------------------------------------------------
subroutine avg_blk_compute(avg_blk,mpl,nam,geom,bpar,samp,mom_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling
type(mom_blk_type),intent(in) :: mom_blk     ! Moments

! Local variables
integer :: il0,jl0,jl0r,jc3,isub,jsub,ic1a,ic1,nc1a,nc1a_cor
real(kind_real) :: m2_1,m2_2
real(kind_real) :: min_m11,max_m11,min_m11m11,max_m11m11,min_m2m2,max_m2m2,min_m22,max_m22,min_cor,max_cor
real(kind_real) :: min_m11_tot,max_m11_tot,min_m11m11_tot,max_m11m11_tot,min_m2m2_tot,max_m2m2_tot,min_m22_tot,max_m22_tot
real(kind_real) :: min_cor_tot,max_cor_tot
real(kind_real),allocatable :: list_m11(:),list_m11m11(:,:,:),list_m2m2(:,:,:),list_m22(:,:),list_cor(:),hist(:,:,:,:)
logical :: involved,valid

! Associate
associate(ic2=>avg_blk%ic2,ib=>avg_blk%ib)

if ((ic2==0).or.nam%local_diag) then
   ! Copy variance
   if (ic2==0) then
      do isub=1,avg_blk%nsub
         do il0=1,geom%nl0
            jl0r = bpar%il0rz(il0,ib)
            avg_blk%m2(il0,isub) = sum(mom_blk%m2_1(:,il0,isub),mask=mpl%msv%isnotr(mom_blk%m2_1(:,il0,isub))) &
                                 & /real(count(samp%c1l0_log(:,il0)),kind_real)
            if (nam%var_filter.and.(.not.nam%gau_approx)) avg_blk%m4(il0,isub) = sum(mom_blk%m22(:,1,jl0r,il0,isub), &
          & mask=mpl%msv%isnotr(mom_blk%m22(:,1,jl0r,il0,isub)))/real(count(samp%c1l0_log(:,il0)),kind_real)
         end do
      end do
   else
      avg_blk%m2 = 0.0
      if (nam%var_filter.and.(.not.nam%gau_approx)) avg_blk%m4 = 0.0
      do ic1a=1,samp%nc1a
         ic1 = samp%c1a_to_c1(ic1a)
         if (ic1==samp%c2_to_c1(ic2)) then
            do isub=1,avg_blk%nsub
               do il0=1,geom%nl0
                  jl0r = bpar%il0rz(il0,ib)
                  avg_blk%m2(il0,isub) = mom_blk%m2_1(ic1a,il0,isub)
                  if (nam%var_filter.and.(.not.nam%gau_approx)) avg_blk%m4(il0,isub) = mom_blk%m22(ic1a,1,jl0r,il0,isub)
               end do
            end do
         end if
      end do
   end if

   ! Check whether this task is involved
   if (ic2>0) then
      involved = any(samp%local_mask(:,ic2))
   else
      involved = .true.
   end if

   if (involved) then
      ! Allocation
      allocate(list_m11(samp%nc1a))
      allocate(list_m11m11(samp%nc1a,avg_blk%nsub,avg_blk%nsub))
      allocate(list_m2m2(samp%nc1a,avg_blk%nsub,avg_blk%nsub))
      if (.not.nam%gau_approx) allocate(list_m22(samp%nc1a,avg_blk%nsub))
      allocate(list_cor(samp%nc1a))

      ! Average
      do il0=1,geom%nl0
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)

            do jc3=1,bpar%nc3(ib)
               ! Fill lists
               !$omp parallel do schedule(static) private(ic1a,ic1,valid,m2_1,m2_2,isub,jsub)
               do ic1a=1,samp%nc1a
                  ! Index
                  ic1 = samp%c1a_to_c1(ic1a)

                  ! Check validity
                  valid = samp%c1l0_log(ic1,il0).and.samp%c1c3l0_log(ic1,jc3,jl0)
                  if (ic2>0) valid = valid.and.samp%local_mask(ic1a,ic2)
                  if (trim(nam%mask_type)=='stddev') then
                     m2_1 = sum(mom_blk%m2_1(ic1a,il0,:))/real(avg_blk%nsub,kind_real)
                     m2_2 = sum(mom_blk%m2_2(ic1a,jc3,jl0,:))/real(avg_blk%nsub,kind_real)
                     if (trim(nam%mask_lu)=='lower') then
                        valid = valid.and.(m2_1>nam%mask_th**2).and.(m2_2>nam%mask_th**2)
                     elseif (trim(nam%mask_lu)=='upper') then
                        valid = valid.and.(m2_1<nam%mask_th**2).and.(m2_2<nam%mask_th**2)
                     end if
                  end if

                  if (valid) then
                     ! Averages for diagnostics
                     list_m11(ic1a) = sum(mom_blk%m11(ic1a,jc3,jl0r,il0,:))/real(avg_blk%nsub,kind_real)
                     do isub=1,avg_blk%nsub
                        do jsub=1,avg_blk%nsub
                           list_m11m11(ic1a,jsub,isub) = mom_blk%m11(ic1a,jc3,jl0r,il0,isub)*mom_blk%m11(ic1a,jc3,jl0r,il0,jsub)
                           list_m2m2(ic1a,jsub,isub) = mom_blk%m2_1(ic1a,il0,isub)*mom_blk%m2_2(ic1a,jc3,jl0,jsub)
                        end do
                        if (.not.nam%gau_approx) list_m22(ic1a,isub) = mom_blk%m22(ic1a,jc3,jl0r,il0,isub)
                     end do

                     ! Correlation
                     m2_1 = sum(mom_blk%m2_1(ic1a,il0,:))/real(avg_blk%nsub,kind_real)
                     m2_2 = sum(mom_blk%m2_2(ic1a,jc3,jl0,:))/real(avg_blk%nsub,kind_real)
                     if ((m2_1>0.0).and.(m2_2>0.0)) then
                        list_cor(ic1a) = list_m11(ic1a)/sqrt(m2_1*m2_2)
                        if (sup(abs(list_cor(ic1a)),1.0_kind_real)) list_cor(ic1a) = mpl%msv%valr
                     else
                        list_cor(ic1a) = mpl%msv%valr
                     end if
                  else
                     ! Missing value
                     list_m11(ic1a) = mpl%msv%valr
                     do isub=1,avg_blk%nsub
                        do jsub=1,avg_blk%nsub
                           list_m11m11(ic1a,jsub,isub) = mpl%msv%valr
                           list_m2m2(ic1a,jsub,isub) = mpl%msv%valr
                        end do
                        if (.not.nam%gau_approx) list_m22(ic1a,isub) = mpl%msv%valr
                     end do
                     list_cor(ic1a) = mpl%msv%valr
                  end if
               end do
               !$omp end parallel do

               ! Number of valid points
               nc1a = count(mpl%msv%isnotr(list_m11))
               nc1a_cor = count(mpl%msv%isnotr(list_cor))
               avg_blk%nc1a(jc3,jl0r,il0) = real(nc1a,kind_real)
               avg_blk%nc1a_cor(jc3,jl0r,il0) = real(nc1a_cor,kind_real)

               ! Average
               if (nc1a>0) then
                  avg_blk%m11(jc3,jl0r,il0) = sum(list_m11,mask=mpl%msv%isnotr(list_m11))
                  do isub=1,avg_blk%nsub
                     do jsub=1,avg_blk%nsub
                        avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) = sum(list_m11m11(:,jsub,isub),mask=mpl%msv%isnotr(list_m11))
                        avg_blk%m2m2(jc3,jl0r,il0,jsub,isub) = sum(list_m2m2(:,jsub,isub),mask=mpl%msv%isnotr(list_m11))
                     end do
                     if (.not.nam%gau_approx) avg_blk%m22(jc3,jl0r,il0,isub) = sum(list_m22(:,isub),mask=mpl%msv%isnotr(list_m11))
                  end do
               else
                  avg_blk%m11(jc3,jl0r,il0) = 0.0
                  do isub=1,avg_blk%nsub
                     do jsub=1,avg_blk%nsub
                        avg_blk%m11m11(jc3,jl0r,il0,jsub,isub) = 0.0
                        avg_blk%m2m2(jc3,jl0r,il0,jsub,isub) = 0.0
                     end do
                     if (.not.nam%gau_approx) avg_blk%m22(jc3,jl0r,il0,isub) = 0.0
                  end do
               end if
               if (nc1a_cor>0) then
                  avg_blk%cor(jc3,jl0r,il0) = sum(list_cor,mask=mpl%msv%isnotr(list_cor))
               else
                  avg_blk%cor(jc3,jl0r,il0) = 0.0
               end if

               if ((nam%avg_nbins>0).and.(ic2==0)) then
                  ! Get min/max values for each task
                  if (nc1a>0) then
                     min_m11 = minval(list_m11,mask=mpl%msv%isnotr(list_m11))
                     max_m11 = maxval(list_m11,mask=mpl%msv%isnotr(list_m11))
                     min_m11m11 = minval(list_m11m11,mask=mpl%msv%isnotr(list_m11m11))
                     max_m11m11 = maxval(list_m11m11,mask=mpl%msv%isnotr(list_m11m11))
                     min_m2m2 = minval(list_m2m2,mask=mpl%msv%isnotr(list_m2m2))
                     max_m2m2 = maxval(list_m2m2,mask=mpl%msv%isnotr(list_m2m2))
                     if (.not.nam%gau_approx) then
                        min_m22 = minval(list_m22,mask=mpl%msv%isnotr(list_m22))
                        max_m22 = maxval(list_m22,mask=mpl%msv%isnotr(list_m22))
                     end if
                  else
                     min_m11 = huge(1.0)
                     max_m11 = -huge(1.0)
                     min_m11m11 = huge(1.0)
                     max_m11m11 = -huge(1.0)
                     min_m2m2 = huge(1.0)
                     max_m2m2 = -huge(1.0)
                     if (.not.nam%gau_approx) then
                        min_m22 = huge(1.0)
                        max_m22 = -huge(1.0)
                     end if
                  end if
                  if (nc1a>0) then
                     min_cor = minval(list_cor,mask=mpl%msv%isnotr(list_cor))
                     max_cor = maxval(list_cor,mask=mpl%msv%isnotr(list_cor))
                  else
                     min_m11 = huge(1.0)
                     max_cor = -huge(1.0)
                  end if

                  ! Gather min/max values
                  call mpl%f_comm%allreduce(min_m11,min_m11_tot,fckit_mpi_min())
                  call mpl%f_comm%allreduce(max_m11,max_m11_tot,fckit_mpi_max())
                  call mpl%f_comm%allreduce(min_m11m11,min_m11m11_tot,fckit_mpi_min())
                  call mpl%f_comm%allreduce(max_m11m11,max_m11m11_tot,fckit_mpi_max())
                  call mpl%f_comm%allreduce(min_m2m2,min_m2m2_tot,fckit_mpi_min())
                  call mpl%f_comm%allreduce(max_m2m2,max_m2m2_tot,fckit_mpi_max())
                  if (.not.nam%gau_approx) then
                     call mpl%f_comm%allreduce(min_m22,min_m22_tot,fckit_mpi_min())
                     call mpl%f_comm%allreduce(max_m22,max_m22_tot,fckit_mpi_max())
                  end if
                  call mpl%f_comm%allreduce(min_cor,min_cor_tot,fckit_mpi_min())
                  call mpl%f_comm%allreduce(max_cor,max_cor_tot,fckit_mpi_max())

                  ! Compute histograms on each task
                  call histogram(mpl,samp%nc1a,list_m11,nam%avg_nbins,min_m11_tot,max_m11_tot,avg_blk%m11_bins(:,jc3,jl0r,il0), &
                & avg_blk%m11_hist(:,jc3,jl0r,il0))
                  call histogram(mpl,samp%nc1a*avg_blk%nsub**2,pack(list_m11m11,mask=.true.),nam%avg_nbins,min_m11m11_tot, &
                & max_m11m11_tot,avg_blk%m11m11_bins(:,jc3,jl0r,il0),avg_blk%m11m11_hist(:,jc3,jl0r,il0))
                  call histogram(mpl,samp%nc1a*avg_blk%nsub**2,pack(list_m2m2,mask=.true.),nam%avg_nbins,min_m2m2_tot, &
                & max_m2m2_tot,avg_blk%m2m2_bins(:,jc3,jl0r,il0),avg_blk%m2m2_hist(:,jc3,jl0r,il0))
                  if (.not.nam%gau_approx) call histogram(mpl,samp%nc1a*avg_blk%nsub,pack(list_m22,mask=.true.), &
                   & nam%avg_nbins,min_m22_tot,max_m22_tot,avg_blk%m22_bins(:,jc3,jl0r,il0),avg_blk%m22_hist(:,jc3,jl0r,il0))
                  call histogram(mpl,samp%nc1a,list_cor,nam%avg_nbins,min_cor_tot,max_cor_tot,avg_blk%cor_bins(:,jc3,jl0r,il0), &
                   & avg_blk%cor_hist(:,jc3,jl0r,il0))
               end if
            end do
         end do
      end do

      ! Release memory
      deallocate(list_m11)
      deallocate(list_m11m11)
      deallocate(list_m2m2)
      if (.not.nam%gau_approx) deallocate(list_m22)
      deallocate(list_cor)
   else
      ! Set to zero for this task
      avg_blk%nc1a = 0
      avg_blk%m11 = 0.0
      avg_blk%m11m11 = 0.0
      avg_blk%m2m2 = 0.0
      if (.not.nam%gau_approx) avg_blk%m22 = 0.0
      avg_blk%nc1a_cor = 0.0
      avg_blk%cor = 0.0
   end if
end if


! Gather histograms
if ((nam%avg_nbins>0).and.(ic2==0)) then
   ! Allocation
   allocate(hist(nam%avg_nbins,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))

   ! Gather data
   hist = avg_blk%m11_hist
   call mpl%f_comm%allreduce(hist,avg_blk%m11_hist,fckit_mpi_sum())
   hist = avg_blk%m11m11_hist
   call mpl%f_comm%allreduce(hist,avg_blk%m11m11_hist,fckit_mpi_sum())
   hist = avg_blk%m2m2_hist
   call mpl%f_comm%allreduce(hist,avg_blk%m2m2_hist,fckit_mpi_sum())
   if (.not.nam%gau_approx) then
      hist = avg_blk%m22_hist
      call mpl%f_comm%allreduce(hist,avg_blk%m22_hist,fckit_mpi_sum())
   end if
   hist = avg_blk%cor_hist
   call mpl%f_comm%allreduce(hist,avg_blk%cor_hist,fckit_mpi_sum())

   ! Normalization
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,bpar%nc3(ib)
            if (sum(avg_blk%m11_hist(:,jc3,jl0r,il0))>0.0) avg_blk%m11_hist(:,jc3,jl0r,il0) = &
          & avg_blk%m11_hist(:,jc3,jl0r,il0)/sum(avg_blk%m11_hist(:,jc3,jl0r,il0))
            if (sum(avg_blk%m11m11_hist(:,jc3,jl0r,il0))>0.0) avg_blk%m11m11_hist(:,jc3,jl0r,il0) = &
          & avg_blk%m11m11_hist(:,jc3,jl0r,il0)/sum(avg_blk%m11m11_hist(:,jc3,jl0r,il0))
            if (sum(avg_blk%m2m2_hist(:,jc3,jl0r,il0))>0.0) avg_blk%m2m2_hist(:,jc3,jl0r,il0) = &
          & avg_blk%m2m2_hist(:,jc3,jl0r,il0)/sum(avg_blk%m2m2_hist(:,jc3,jl0r,il0))
            if (.not.nam%gau_approx) then
               if (sum(avg_blk%m22_hist(:,jc3,jl0r,il0))>0.0) avg_blk%m22_hist(:,jc3,jl0r,il0) = &
             & avg_blk%m22_hist(:,jc3,jl0r,il0)/sum(avg_blk%m22_hist(:,jc3,jl0r,il0))
            end if
            if (sum(avg_blk%cor_hist(:,jc3,jl0r,il0))>0.0) avg_blk%cor_hist(:,jc3,jl0r,il0) = &
          & avg_blk%cor_hist(:,jc3,jl0r,il0)/sum(avg_blk%cor_hist(:,jc3,jl0r,il0))
         end do
      end do
   end do

   ! Release memory
   deallocate(hist)
end if

! End associate
end associate

end subroutine avg_blk_compute

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_asy
! Purpose: compute asymptotic statistics
!----------------------------------------------------------------------
subroutine avg_blk_compute_asy(avg_blk,mpl,nam,geom,bpar,ne)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
integer,intent(in) :: ne                     ! Ensemble size

! Local variables
integer :: il0,jl0r,jc3,isub,jsub,n
real(kind_real) :: P1,P3,P4,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17
real(kind_real),allocatable :: m11asysq(:,:),m2m2asy(:,:),m22asy(:)

! Associate
associate(ic2=>avg_blk%ic2,ib=>avg_blk%ib)

if ((ic2==0).or.nam%local_diag) then
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
            if (avg_blk%m11asysq(jc3,jl0r,il0)<0.0) avg_blk%m11asysq(jc3,jl0r,il0) = mpl%msv%valr
            if (avg_blk%m2m2asy(jc3,jl0r,il0)<0.0) avg_blk%m2m2asy(jc3,jl0r,il0) = mpl%msv%valr
            if (.not.nam%gau_approx) then
               if (avg_blk%m22asy(jc3,jl0r,il0)<0.0) avg_blk%m22asy(jc3,jl0r,il0) = mpl%msv%valr
            end if

            ! Squared covariance average
            if (nam%gau_approx) then
               ! Gaussian approximation
               if (mpl%msv%isnotr(avg_blk%m11asysq(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m2m2asy(jc3,jl0r,il0))) then
                  avg_blk%m11sq(jc3,jl0r,il0) = P16*avg_blk%m11asysq(jc3,jl0r,il0)+P4*avg_blk%m2m2asy(jc3,jl0r,il0)
               else
                  avg_blk%m11sq(jc3,jl0r,il0) = mpl%msv%valr
               end if
            else
               ! General case
               if (mpl%msv%isnotr(avg_blk%m22asy(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m11asysq(jc3,jl0r,il0)) &
             & .and.mpl%msv%isnotr(avg_blk%m2m2asy(jc3,jl0r,il0))) then
                  avg_blk%m11sq(jc3,jl0r,il0) = P1*avg_blk%m22asy(jc3,jl0r,il0)+P14*avg_blk%m11asysq(jc3,jl0r,il0) &
                                      & +P3*avg_blk%m2m2asy(jc3,jl0r,il0)
               else
                  avg_blk%m11sq(jc3,jl0r,il0) = mpl%msv%valr
               end if
            end if

            ! Check value
            if (mpl%msv%isr(avg_blk%m11sq(jc3,jl0r,il0))) then
               if (inf(avg_blk%m11sq(jc3,jl0r,il0),avg_blk%m11asysq(jc3,jl0r,il0))) avg_blk%m11sq(jc3,jl0r,il0) = mpl%msv%valr
               if (inf(avg_blk%m11sq(jc3,jl0r,il0),avg_blk%m11(jc3,jl0r,il0)**2)) avg_blk%m11sq(jc3,jl0r,il0) = mpl%msv%valr
            end if

            ! Release memory
            deallocate(m11asysq)
            deallocate(m2m2asy)
            deallocate(m22asy)
         end do
      end do
   end do
   !$omp end parallel do
end if

! End associate
end associate

end subroutine avg_blk_compute_asy

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_hyb
! Purpose: compute averaged statistics via spatial-angular erogodicity assumption, for hybrid covariance products
!----------------------------------------------------------------------
subroutine avg_blk_compute_hyb(avg_blk_hyb,mpl,geom,bpar,avg_blk_sta)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk_hyb ! Hybrid averaged statistics block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(geom_type),intent(in) :: geom               ! Geometry
type(bpar_type),intent(in) :: bpar               ! Block parameters
type(avg_blk_type),intent(in) :: avg_blk_sta     ! Static averaged statistics block

! Local variables
integer :: il0,jl0r,jc3

! Associate
associate(ib=>avg_blk_hyb%ib)

do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (mpl%msv%isnotr(avg_blk_hyb%m11(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk_sta%m11(jc3,jl0r,il0))) then
            ! Compute product
            avg_blk_hyb%m11sta(jc3,jl0r,il0) = avg_blk_hyb%m11(jc3,jl0r,il0)*avg_blk_sta%m11(jc3,jl0r,il0)
            avg_blk_hyb%stasq(jc3,jl0r,il0) = avg_blk_sta%m11(jc3,jl0r,il0)**2
         else
            ! Missing values
            avg_blk_hyb%m11sta(jc3,jl0r,il0) = mpl%msv%valr
            avg_blk_hyb%stasq(jc3,jl0r,il0) = mpl%msv%valr
         end if
      end do
   end do
end do

! End associate
end associate

end subroutine avg_blk_compute_hyb

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_deh
! Purpose: compute averaged statistics via spatial-angular erogodicity assumption, for LR covariance/HR covariance and LR covariance/HR asymptotic covariance products
!----------------------------------------------------------------------
subroutine avg_blk_compute_deh(avg_blk,mpl,nam,geom,bpar,samp,mom_blk,mom_lr_blk)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling
type(mom_blk_type),intent(in) :: mom_blk     ! Moments block
type(mom_blk_type),intent(in) :: mom_lr_blk  ! Low-resolution moments block

! Local variables
integer :: il0,jl0,jl0r,jc3,isub,jsub,ic1a,ic1,nc1amax,nc1a
real(kind_real),allocatable :: list_m11lrm11(:,:,:)
logical :: valid
character(len=1024),parameter :: subr = 'avg_blk_compute_deh'

! Associate
associate(ic2=>avg_blk%ic2,ib=>avg_blk%ib)

if ((ic2==0).or.nam%local_diag) then
   ! Check number of sub-ensembles
   if (avg_blk%nsub/=avg_blk%nsub) call mpl%abort(subr,'different number of sub-ensembles')

   ! Average
   !$omp parallel do schedule(static) private(il0,jl0r,jl0,nc1amax,list_m11lrm11,jc3,nc1a,ic1a,ic1,valid,isub,jsub)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)

         ! Allocation
         if (ic2>0) then
            nc1amax = count(samp%local_mask(:,ic2))
         else
            nc1amax = samp%nc1a
         end if
         allocate(list_m11lrm11(nc1amax,avg_blk%nsub,avg_blk%nsub))

         do jc3=1,bpar%nc3(ib)
            ! Fill lists
            nc1a = 0
            do ic1a=1,samp%nc1a
               ! Index
               ic1 = samp%c1a_to_c1(ic1a)

               ! Check validity
               valid = samp%c1l0_log(ic1,il0).and.samp%c1c3l0_log(ic1,jc3,jl0)
               if (ic2>0) valid = valid.and.samp%local_mask(ic1a,ic2)

               if (valid) then
                  ! Update
                  nc1a = nc1a+1

                  ! Averages for diagnostics
                  do isub=1,avg_blk%nsub
                     do jsub=1,avg_blk%nsub
                        list_m11lrm11(nc1a,jsub,isub) = mom_blk%m11(ic1a,jc3,jl0r,il0,isub)*mom_lr_blk%m11(ic1a,jc3,jl0r,il0,jsub)
                     end do
                  end do
               end if
            end do

            ! Average
            avg_blk%nc1a(jc3,jl0r,il0) = real(nc1a,kind_real)
            do isub=1,avg_blk%nsub
               do jsub=1,avg_blk%nsub
                  if (avg_blk%nc1a(jc3,jl0r,il0)>0) then
                     avg_blk%m11lrm11sub(jc3,jl0r,il0,jsub,isub) = sum(list_m11lrm11(1:nc1a,jsub,isub))
                  else
                     avg_blk%m11lrm11sub(jc3,jl0r,il0,jsub,isub) = mpl%msv%valr
                  end if
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

end subroutine avg_blk_compute_deh

!----------------------------------------------------------------------
! Subroutine: avg_blk_compute_asy_deh
! Purpose: compute LR covariance/HR asymptotic covariance products
!----------------------------------------------------------------------
subroutine avg_blk_compute_asy_deh(avg_blk,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_blk_type),intent(inout) :: avg_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters

! Local variables
integer :: il0,jl0r,jc3,isub,jsub

! Associate
associate(ic2=>avg_blk%ic2,ib=>avg_blk%ib)

if ((ic2==0).or.nam%local_diag) then
   ! Normalize
   !$omp parallel do schedule(static) private(il0,jl0r,jc3,isub,jsub)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,bpar%nc3(ib)
            if (avg_blk%nc1a(jc3,jl0r,il0)>0.0) then
               do isub=1,avg_blk%nsub
                  do jsub=1,avg_blk%nsub
                     avg_blk%m11lrm11sub(jc3,jl0r,il0,jsub,isub) = avg_blk%m11lrm11sub(jc3,jl0r,il0,jsub,isub) &
                                                                 & /avg_blk%nc1a(jc3,jl0r,il0)
                  end do
               end do
            else
               avg_blk%m11lrm11sub(jc3,jl0r,il0,:,:) = mpl%msv%valr
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
            if (mpl%msv%isanynotr(avg_blk%m11lrm11sub(jc3,jl0r,il0,:,:))) then
               avg_blk%m11lrm11(jc3,jl0r,il0) = sum(avg_blk%m11lrm11sub(jc3,jl0r,il0,:,:), &
                                              & mask=mpl%msv%isnotr(avg_blk%m11lrm11sub(jc3,jl0r,il0,:,:))) &
                                              & /real(count(mpl%msv%isnotr(avg_blk%m11lrm11sub(jc3,jl0r,il0,:,:))),kind_real)
            else
               avg_blk%m11lrm11(jc3,jl0r,il0) = mpl%msv%valr
            end if
         end do
      end do
   end do
   !$omp end parallel do

   ! Define asymptotic covariance
   avg_blk%m11lrm11asy = avg_blk%m11lrm11
end if

! End associate
end associate

end subroutine avg_blk_compute_asy_deh

end module type_avg_blk
