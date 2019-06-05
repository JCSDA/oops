!----------------------------------------------------------------------
! Module: type_diag_blk
! Purpose: diagnostic block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_diag_blk

use netcdf
!$ use omp_lib
use tools_fit, only: fast_fit,ver_fill
use tools_func, only: fit_diag,fit_diag_dble
use tools_kinds, only: kind_real,nc_kind_real
use tools_repro, only: sup
use type_avg_blk, only: avg_blk_type
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_minim, only: minim_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

integer,parameter :: nsc = 50                          ! Number of iterations for the scaling optimization
logical :: lprt = .false.                              ! Optimization print
real(kind_real),parameter :: maxfactor = 2.0_kind_real ! Maximum factor for diagnostics with respect to the origin

! Diagnostic block derived type
type diag_blk_type
   integer :: ic2a                                ! Local index
   integer :: ib                                  ! Block index
   character(len=1024) :: name                    ! Name
   logical :: double_fit                          ! Double fit flag

   real(kind_real),allocatable :: raw(:,:,:)      ! Raw diagnostic
   real(kind_real),allocatable :: valid(:,:,:)    ! Number of valid couples
   real(kind_real),allocatable :: raw_coef_ens(:) ! Raw ensemble coefficient
   real(kind_real) :: raw_coef_sta                ! Raw static coefficient
   real(kind_real),allocatable :: fit(:,:,:)      ! Fit
   real(kind_real),allocatable :: fit_rh(:)       ! Horizontal fit support radius
   real(kind_real),allocatable :: fit_rv(:)       ! Vertical fit support radius
   real(kind_real),allocatable :: fit_rv_rfac(:)  ! Vertical fit support radius ratio for the positive component (for double-radius fit)
   real(kind_real),allocatable :: fit_rv_coef(:)  ! Vertical fit coefficient (for double-radius fit)
   real(kind_real),allocatable :: distv(:,:)      ! Reduced vertical distance
contains
   procedure :: alloc => diag_blk_alloc
   procedure :: dealloc => diag_blk_dealloc
   procedure :: write => diag_blk_write
   procedure :: normalization => diag_blk_normalization
   procedure :: fitting => diag_blk_fitting
   procedure :: localization => diag_blk_localization
   procedure :: hybridization => diag_blk_hybridization
   procedure :: dualens => diag_blk_dualens
end type diag_blk_type

private
public :: diag_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: diag_blk_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine diag_blk_alloc(diag_blk,mpl,nam,geom,bpar,samp,ic2a,ib,prefix,double_fit)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk ! Diagnostic block
type(mpl_type),intent(inout) :: mpl            ! MPI data
type(nam_type),intent(in) :: nam               ! Namelist
type(geom_type),intent(in) :: geom             ! Geometry
type(bpar_type),intent(in) :: bpar             ! Block parameters
type(samp_type),intent(in) :: samp             ! Sampling
integer,intent(in) :: ic2a                     ! Local index
integer,intent(in) :: ib                       ! Block index
character(len=*),intent(in) :: prefix          ! Block prefix
logical,intent(in) :: double_fit               ! Double fit

! Local variables
integer :: ic0,ic2,il0,jl0
real(kind_real) :: vunit(geom%nl0)

! Set attributes
diag_blk%ic2a = ic2a
diag_blk%ib = ib
diag_blk%name = trim(prefix)//'_'//trim(bpar%blockname(ib))
diag_blk%double_fit = double_fit

! Allocation
allocate(diag_blk%raw_coef_ens(geom%nl0))
if ((ic2a==0).or.nam%local_diag) then
   allocate(diag_blk%raw(nam%nc3,bpar%nl0r(ib),geom%nl0))
   allocate(diag_blk%valid(nam%nc3,bpar%nl0r(ib),geom%nl0))
   if (trim(nam%minim_algo)/='none') then
      allocate(diag_blk%fit(nam%nc3,bpar%nl0r(ib),geom%nl0))
      allocate(diag_blk%fit_rh(geom%nl0))
      allocate(diag_blk%fit_rv(geom%nl0))
      if (diag_blk%double_fit) then
         allocate(diag_blk%fit_rv_rfac(geom%nl0))
         allocate(diag_blk%fit_rv_coef(geom%nl0))
      end if
      allocate(diag_blk%distv(geom%nl0,geom%nl0))
   end if
end if

! Initialization
diag_blk%raw_coef_ens = mpl%msv%valr
diag_blk%raw_coef_sta = mpl%msv%valr
if ((ic2a==0).or.nam%local_diag) then
   diag_blk%raw = mpl%msv%valr
   diag_blk%valid = mpl%msv%valr
   if (trim(nam%minim_algo)/='none') then
      diag_blk%fit = mpl%msv%valr
      diag_blk%fit_rh = mpl%msv%valr
      diag_blk%fit_rv = mpl%msv%valr
      if (diag_blk%double_fit) then
         diag_blk%fit_rv_rfac = mpl%msv%valr
         diag_blk%fit_rv_coef = mpl%msv%valr
      end if
   end if
end if

! Vertical distance
if (((ic2a==0).or.nam%local_diag).and.(trim(nam%minim_algo)/='none')) then
   if (ic2a==0) then
      vunit = geom%vunitavg
   else
      ic2 = samp%c2a_to_c2(ic2a)
      ic0 = samp%c2_to_c0(ic2)
      vunit = geom%vunit_c0(ic0,:)
   end if
   do il0=1,geom%nl0
      do jl0=1,geom%nl0
         diag_blk%distv(jl0,il0) = abs(vunit(il0)-vunit(jl0))
      end do
   end do
end if

end subroutine diag_blk_alloc

!----------------------------------------------------------------------
! Subroutine: diag_blk_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine diag_blk_dealloc(diag_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk ! Diagnostic block

! Release memory
if (allocated(diag_blk%raw)) deallocate(diag_blk%raw)
if (allocated(diag_blk%valid)) deallocate(diag_blk%valid)
if (allocated(diag_blk%raw_coef_ens)) deallocate(diag_blk%raw_coef_ens)
if (allocated(diag_blk%fit)) deallocate(diag_blk%fit)
if (allocated(diag_blk%fit_rh)) deallocate(diag_blk%fit_rh)
if (allocated(diag_blk%fit_rv)) deallocate(diag_blk%fit_rv)
if (allocated(diag_blk%fit_rv_rfac)) deallocate(diag_blk%fit_rv_rfac)
if (allocated(diag_blk%fit_rv_coef)) deallocate(diag_blk%fit_rv_coef)
if (allocated(diag_blk%distv)) deallocate(diag_blk%distv)

end subroutine diag_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: diag_blk_write
! Purpose: write
!----------------------------------------------------------------------
subroutine diag_blk_write(diag_blk,mpl,nam,geom,bpar,filename)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk ! Diagnostic block
type(mpl_type),intent(inout) :: mpl            ! MPI data
type(nam_type),intent(in) :: nam               ! Namelist
type(geom_type),intent(in) :: geom             ! Geometry
type(bpar_type),intent(in) :: bpar             ! Block parameters
character(len=*),intent(in) :: filename        ! File name

! Local variables
integer :: info,info_coord,ncid,one_id,nc3_id,nl0r_id,nl0_1_id,nl0_2_id,disth_id,vunit_id
integer :: raw_id,valid_id,raw_zs_id,raw_coef_ens_id,raw_coef_sta_id,l0rl0_to_l0_id
integer :: fit_id,fit_zs_id,fit_rh_id,fit_rv_id,fit_rv_rfac_id,fit_rv_coef_id
integer :: il0,jl0r,jl0
character(len=1024),parameter :: subr = 'diag_blk_write'

! Associate
associate(ib=>diag_blk%ib,ic2a=>diag_blk%ic2a)

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
one_id = mpl%ncdimcheck(subr,ncid,'one',1,.true.)
nc3_id = mpl%ncdimcheck(subr,ncid,'nc3',nam%nc3,.true.)
nl0r_id = mpl%ncdimcheck(subr,ncid,'nl0r',bpar%nl0rmax,.true.)
nl0_1_id = mpl%ncdimcheck(subr,ncid,'nl0_1',geom%nl0,.true.)
nl0_2_id = mpl%ncdimcheck(subr,ncid,'nl0_2',geom%nl0,.true.)
info_coord = nf90_inq_varid(ncid,'disth',disth_id)
if (info_coord/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'disth',nc_kind_real,(/nc3_id/),disth_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'vunit',nc_kind_real,(/nl0_1_id/),vunit_id))
end if

! Define variables if necessary
if (mpl%msv%isanynotr(diag_blk%raw_coef_ens)) then
  info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_raw_coef_ens',raw_coef_ens_id)
  if (info/=nf90_noerr) then
     call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_raw_coef_ens',nc_kind_real,(/nl0_1_id/),raw_coef_ens_id))
     call mpl%ncerr(subr,nf90_put_att(ncid,raw_coef_ens_id,'_FillValue',mpl%msv%valr))
  end if
end if
if ((ic2a==0).or.nam%local_diag) then
   if (mpl%msv%isanynotr(diag_blk%raw)) then
      info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_raw',raw_id)
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_raw',nc_kind_real,(/nc3_id,nl0r_id,nl0_1_id/),raw_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,raw_id,'_FillValue',mpl%msv%valr))
      end if
      info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_valid',valid_id)
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_valid',nc_kind_real,(/nc3_id,nl0r_id,nl0_1_id/),valid_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,valid_id,'_FillValue',mpl%msv%valr))
      end if
      if (bpar%nl0rmax/=geom%nl0) then
         info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_raw_zs',raw_zs_id)
         if (info/=nf90_noerr) then
            call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_raw_zs',nc_kind_real,(/nl0_2_id,nl0_1_id/),raw_zs_id))
            call mpl%ncerr(subr,nf90_put_att(ncid,raw_zs_id,'_FillValue',mpl%msv%valr))
         end if
      end if
   end if
   if (mpl%msv%isnotr(diag_blk%raw_coef_sta)) then
      info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_raw_coef_sta',raw_coef_sta_id)
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_raw_coef_sta',nc_kind_real,(/one_id/),raw_coef_sta_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,raw_coef_sta_id,'_FillValue',mpl%msv%valr))
      end if
   end if
   if ((trim(nam%minim_algo)/='none').and.(mpl%msv%isanynotr(diag_blk%fit))) then
      info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_fit',fit_id)
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit',nc_kind_real,(/nc3_id,nl0r_id,nl0_1_id/),fit_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,fit_id,'_FillValue',mpl%msv%valr))
      end if
      if (bpar%nl0rmax/=geom%nl0) then
         info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_fit_zs',fit_zs_id)
         if (info/=nf90_noerr) then
            call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit_zs',nc_kind_real,(/nl0_2_id,nl0_1_id/),fit_zs_id))
            call mpl%ncerr(subr,nf90_put_att(ncid,fit_zs_id,'_FillValue',mpl%msv%valr))
         end if
      end if
      info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_fit_rh',fit_rh_id)
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit_rh',nc_kind_real,(/nl0_1_id/),fit_rh_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,fit_rh_id,'_FillValue',mpl%msv%valr))
      end if
      info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_fit_rv',fit_rv_id)
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit_rv',nc_kind_real,(/nl0_1_id/),fit_rv_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,fit_rv_id,'_FillValue',mpl%msv%valr))
      end if
      if (diag_blk%double_fit) then
         info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_fit_rv_rfac',fit_rv_rfac_id)
         if (info/=nf90_noerr) then
            call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit_rv_rfac',nc_kind_real,(/nl0_1_id/),fit_rv_rfac_id))
            call mpl%ncerr(subr,nf90_put_att(ncid,fit_rv_rfac_id,'_FillValue',mpl%msv%valr))
         end if
         info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_fit_rv_coef',fit_rv_coef_id)
         if (info/=nf90_noerr) then
            call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit_rv_coef',nc_kind_real,(/nl0_1_id/),fit_rv_coef_id))
            call mpl%ncerr(subr,nf90_put_att(ncid,fit_rv_coef_id,'_FillValue',mpl%msv%valr))
         end if
      end if
   end if
   info = nf90_inq_varid(ncid,trim(diag_blk%name)//'_l0rl0_to_l0',l0rl0_to_l0_id)
   if (info/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_l0rl0_to_l0',nf90_int,(/nl0r_id,nl0_1_id/),l0rl0_to_l0_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,l0rl0_to_l0_id,'_FillValue',mpl%msv%vali))
   end if
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write coordinates if necessary
if (info_coord/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_put_var(ncid,disth_id,geom%disth(1:nam%nc3)))
   call mpl%ncerr(subr,nf90_put_var(ncid,vunit_id,geom%vunitavg))
end if

! Write variables
if (mpl%msv%isanynotr(diag_blk%raw_coef_ens)) call mpl%ncerr(subr,nf90_put_var(ncid,raw_coef_ens_id,diag_blk%raw_coef_ens))
if ((ic2a==0).or.nam%local_diag) then
   if (mpl%msv%isanynotr(diag_blk%raw)) then
      call mpl%ncerr(subr,nf90_put_var(ncid,raw_id,diag_blk%raw))
      call mpl%ncerr(subr,nf90_put_var(ncid,valid_id,diag_blk%valid))
      if (bpar%nl0rmax/=geom%nl0) then
         do il0=1,geom%nl0
            do jl0r=1,bpar%nl0rmax
               jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
               call mpl%ncerr(subr,nf90_put_var(ncid,raw_zs_id,diag_blk%raw(1,jl0r,il0),(/jl0,il0/)))
            end do
         end do
      end if
   end if
   if (mpl%msv%isnotr(diag_blk%raw_coef_sta)) call mpl%ncerr(subr,nf90_put_var(ncid,raw_coef_sta_id,diag_blk%raw_coef_sta))
   if ((trim(nam%minim_algo)/='none').and.(mpl%msv%isanynotr(diag_blk%fit))) then
      call mpl%ncerr(subr,nf90_put_var(ncid,fit_id,diag_blk%fit))
      if (bpar%nl0rmax/=geom%nl0) then
         do il0=1,geom%nl0
            do jl0r=1,bpar%nl0rmax
               jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
               call mpl%ncerr(subr,nf90_put_var(ncid,fit_zs_id,diag_blk%fit(1,jl0r,il0),(/jl0,il0/)))
            end do
          end do
      end if
      call mpl%ncerr(subr,nf90_put_var(ncid,fit_rh_id,diag_blk%fit_rh))
      call mpl%ncerr(subr,nf90_put_var(ncid,fit_rv_id,diag_blk%fit_rv))
      if (diag_blk%double_fit) then
         call mpl%ncerr(subr,nf90_put_var(ncid,fit_rv_rfac_id,diag_blk%fit_rv_rfac))
         call mpl%ncerr(subr,nf90_put_var(ncid,fit_rv_coef_id,diag_blk%fit_rv_coef))
      end if
   end if
   call mpl%ncerr(subr,nf90_put_var(ncid,l0rl0_to_l0_id,bpar%l0rl0b_to_l0(:,:,ib)))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine diag_blk_write

!----------------------------------------------------------------------
! Subroutine: diag_blk_normalization
! Purpose: compute diagnostic block normalization
!----------------------------------------------------------------------
subroutine diag_blk_normalization(diag_blk,mpl,geom,bpar,remove_max)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk ! Diagnostic block
type(mpl_type),intent(inout) :: mpl            ! MPI data
type(geom_type),intent(in) :: geom             ! Geometry
type(bpar_type),intent(in) :: bpar             ! Block parameters
logical,intent(in),optional :: remove_max      ! Remove excessive values

! Local variables
integer :: il0,jl0r,jl0,jc3
logical :: lremove_max

! Associate
associate(ib=>diag_blk%ib)

! Get diagonal values
do il0=1,geom%nl0
   jl0r = bpar%il0rz(il0,ib)
   if (mpl%msv%isnotr(diag_blk%raw(1,jl0r,il0))) then
      diag_blk%raw_coef_ens(il0) = diag_blk%raw(1,jl0r,il0)
   else
      diag_blk%raw_coef_ens(il0) = mpl%msv%valr
   end if
end do

! Normalize
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
      do jc3=1,bpar%nc3(ib)
         if (mpl%msv%isnotr(diag_blk%raw(jc3,jl0r,il0)).and.mpl%msv%isnotr(diag_blk%raw_coef_ens(il0)) &
       & .and.mpl%msv%isnotr(diag_blk%raw_coef_ens(jl0))) &
       & diag_blk%raw(jc3,jl0r,il0) = diag_blk%raw(jc3,jl0r,il0)/sqrt(diag_blk%raw_coef_ens(il0)*diag_blk%raw_coef_ens(jl0))
      end do
   end do
end do

! Remove excessive values compared to the origin point
if (present(remove_max)) then
   lremove_max = remove_max
else
   lremove_max = .false.
end if
if (lremove_max) then
   !$omp parallel do schedule(static) private(il0,jl0r,jc3)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,bpar%nc3(ib)
            if (sup(diag_blk%raw(jc3,jl0r,il0),maxfactor)) diag_blk%raw(jc3,jl0r,il0) = mpl%msv%valr
         end do
      end do
   end do
   !$omp end parallel do
end if

! End associate
end associate

end subroutine diag_blk_normalization

!----------------------------------------------------------------------
! Subroutine: diag_blk_fitting
! Purpose: compute a semi-positive definite fit of a raw function
!----------------------------------------------------------------------
subroutine diag_blk_fitting(diag_blk,mpl,nam,geom,bpar,samp)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk ! Diagnostic block
type(mpl_type),intent(inout) :: mpl            ! MPI data
type(nam_type),intent(in) :: nam               ! Namelist
type(geom_type),intent(in) :: geom             ! Geometry
type(bpar_type),intent(in) :: bpar             ! Block parameters
type(samp_type),intent(in) :: samp             ! Sampling

! Local variables
integer :: ic2,ic0,il0,jl0r,offset,isc
real(kind_real) :: alpha,alpha_opt,mse,mse_opt
real(kind_real) :: vunit(geom%nl0),fit_rh(geom%nl0),fit_rv(geom%nl0)
real(kind_real),allocatable :: rawv(:),distv(:),fit(:,:,:)
type(minim_type) :: minim
character(len=1024),parameter :: subr = 'diag_blk_fitting'

! Associate
associate(ic2a=>diag_blk%ic2a,ib=>diag_blk%ib)

! Check
if (trim(nam%minim_algo)=='none') call mpl%abort(subr,'cannot compute fit if minim_algo = none')

! Allocation
allocate(rawv(bpar%nl0r(ib)))
allocate(distv(bpar%nl0r(ib)))
allocate(fit(nam%nc3,bpar%nl0r(ib),geom%nl0))

! Initialization
diag_blk%fit_rh = mpl%msv%valr
diag_blk%fit_rv = mpl%msv%valr
diag_blk%fit = mpl%msv%valr

! Vertical unit
if (ic2a==0) then
   vunit = geom%vunitavg
else
   ic2 = samp%c2a_to_c2(ic2a)
   ic0 = samp%c2_to_c0(ic2)
   vunit = geom%vunit_c0(ic0,:)
end if

! Fast fit
do il0=1,geom%nl0
   ! Get zero separation level
    jl0r = bpar%il0rz(il0,ib)

   ! Horizontal fast fit
   call fast_fit(mpl,nam%nc3,1,geom%disth,diag_blk%raw(:,jl0r,il0),diag_blk%fit_rh(il0))

   ! Vertical fast fit
   rawv = diag_blk%raw(1,:,il0)
   distv = diag_blk%distv(bpar%l0rl0b_to_l0(:,il0,ib),il0)
   call fast_fit(mpl,bpar%nl0r(ib),jl0r,distv,rawv,diag_blk%fit_rv(il0))
end do

if (any(mpl%msv%isnotr(diag_blk%fit_rh)).and.any(mpl%msv%isnotr(diag_blk%fit_rv))) then
   ! Fill missing values
   call ver_fill(mpl,geom%nl0,vunit,diag_blk%fit_rh)
   call ver_fill(mpl,geom%nl0,vunit,diag_blk%fit_rv)

   ! Vertically homogeneous case
   if (nam%lhomh) diag_blk%fit_rh = sum(diag_blk%fit_rh,mask=mpl%msv%isnotr(diag_blk%fit_rh)) &
 & /real(count(mpl%msv%isnotr(diag_blk%fit_rh)),kind_real)
   if (nam%lhomv) diag_blk%fit_rv = sum(diag_blk%fit_rv,mask=mpl%msv%isnotr(diag_blk%fit_rv)) &
 & /real(count(mpl%msv%isnotr(diag_blk%fit_rv)),kind_real)

   ! Double-fit default parameters
   if (diag_blk%double_fit) then
      diag_blk%fit_rv_rfac = 0.3
      diag_blk%fit_rv_coef = 0.5
   end if

   ! Scaling optimization (brute-force)
   if (all(mpl%msv%isnotr(diag_blk%fit_rh)).and.all(mpl%msv%isnotr(diag_blk%fit_rv))) then
      mse_opt = huge(1.0)
      alpha_opt = 1.0
      do isc=1,nsc
         ! Scaling factor
         alpha = 0.5+real(isc-1,kind_real)/real(nsc-1,kind_real)*(2.0-0.5)

         ! Scaled radii
         fit_rh = alpha*diag_blk%fit_rh
         fit_rv = alpha*diag_blk%fit_rv

         ! Define fit
         if (diag_blk%double_fit) then
            call fit_diag_dble(mpl,nam%nc3,bpar%nl0r(ib),geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth,diag_blk%distv, &
          & fit_rh,fit_rv,diag_blk%fit_rv_rfac,diag_blk%fit_rv_coef,fit)
         else
            call fit_diag(mpl,nam%nc3,bpar%nl0r(ib),geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth,diag_blk%distv,fit_rh,fit_rv,fit)
         end if

         ! MSE
         mse = sum((fit-diag_blk%raw)**2,mask=mpl%msv%isnotr(diag_blk%raw))
         if (mse<mse_opt) then
            mse_opt = mse
            alpha_opt = alpha
         end if
      end do
      diag_blk%fit_rh = alpha_opt*diag_blk%fit_rh
      diag_blk%fit_rv = alpha_opt*diag_blk%fit_rv

      if (lprt) then
         write(mpl%info,'(a)') ''
         call mpl%flush
         write(mpl%info,'(a13,a,f6.1,a)') '','Scaling optimization, cost function decrease:',abs(mse_opt-mse)/mse*100.0,'%'
         call mpl%flush
      end if
   end if

   select case (trim(nam%minim_algo))
   case ('hooke')
      ! Allocation
      minim%nx = 0
      if (nam%lhomh) then
         minim%nx = minim%nx+1
      else
         minim%nx = minim%nx+geom%nl0
      end if
      if (nam%lhomv) then
         minim%nx = minim%nx+1
         if (diag_blk%double_fit) minim%nx = minim%nx+2
      else
         minim%nx = minim%nx+geom%nl0
         if (diag_blk%double_fit) minim%nx = minim%nx+2*geom%nl0
      end if
      minim%ny = nam%nc3*bpar%nl0r(ib)*geom%nl0
      allocate(minim%x(minim%nx))
      allocate(minim%guess(minim%nx))
      allocate(minim%binf(minim%nx))
      allocate(minim%bsup(minim%nx))
      allocate(minim%obs(minim%ny))
      allocate(minim%l0rl0_to_l0(bpar%nl0r(ib),geom%nl0))
      allocate(minim%disth(nam%nc3))
      allocate(minim%distv(geom%nl0,geom%nl0))

      ! Fill minim
      offset = 0
      if (nam%lhomh) then
         minim%guess(offset+1) = diag_blk%fit_rh(1)
         minim%binf(offset+1) = 0.5*minim%guess(offset+1)
         minim%bsup(offset+1) = 1.5*minim%guess(offset+1)
         offset = offset+1
      else
         minim%guess(offset+1:offset+geom%nl0) = diag_blk%fit_rh
         minim%binf(offset+1:offset+geom%nl0) = 0.5*minim%guess(offset+1:offset+geom%nl0)
         minim%bsup(offset+1:offset+geom%nl0) = 1.5*minim%guess(offset+1:offset+geom%nl0)
         offset = offset+geom%nl0
      end if
      if (nam%lhomv) then
         minim%guess(offset+1) = diag_blk%fit_rv(1)
         minim%binf(offset+1) = 0.5*minim%guess(offset+1)
         minim%bsup(offset+1) = 1.5*minim%guess(offset+1)
         offset = offset+1
         if (diag_blk%double_fit) then
            minim%guess(offset+1) = diag_blk%fit_rv_rfac(1)
            minim%binf(offset+1) = 0.0
            minim%bsup(offset+1) = 1.0
            offset = offset+1
            minim%guess(offset+1) = diag_blk%fit_rv_coef(1)
            minim%binf(offset+1) = 0.0
            minim%bsup(offset+1) = 20.0
            offset = offset+1
         end if
      else
         minim%guess(offset+1:offset+geom%nl0) = diag_blk%fit_rv
         minim%binf(offset+1:offset+geom%nl0) = 0.5*minim%guess(offset+1:offset+geom%nl0)
         minim%bsup(offset+1:offset+geom%nl0) = 1.5*minim%guess(offset+1:offset+geom%nl0)
         offset = offset+geom%nl0
         if (diag_blk%double_fit) then
            minim%guess(offset+1:offset+geom%nl0) = diag_blk%fit_rv_rfac
            minim%binf(offset+1:offset+geom%nl0) = 0.0
            minim%bsup(offset+1:offset+geom%nl0) = 1.0
            offset = offset+geom%nl0
            minim%guess(offset+1:offset+geom%nl0) = diag_blk%fit_rv_coef
            minim%binf(offset+1:offset+geom%nl0) = 0.0
            minim%bsup(offset+1:offset+geom%nl0) = 20.0
            offset = offset+geom%nl0
         end if
      end if
      minim%obs = pack(diag_blk%raw,mask=.true.)
      if (diag_blk%double_fit) then
         minim%cost_function = 'fit_diag_dble'
      else
         minim%cost_function = 'fit_diag'
      end if
      minim%algo = nam%minim_algo
      minim%nc3 = nam%nc3
      minim%nl0r = bpar%nl0r(ib)
      minim%nl0 = geom%nl0
      minim%lhomh = nam%lhomh
      minim%lhomv = nam%lhomv
      minim%l0rl0_to_l0 = bpar%l0rl0b_to_l0(:,:,ib)
      minim%disth = geom%disth
      minim%distv = diag_blk%distv

      ! Compute fit
      call minim%compute(mpl,lprt)

      ! Apply bounds
      minim%x = max(minim%binf,min(minim%x,minim%bsup))

      ! Copy parameters
      offset = 0
      if (nam%lhomh) then
         diag_blk%fit_rh = minim%x(offset+1)
         offset = offset+1
      else
         diag_blk%fit_rh = minim%x(offset+1:offset+geom%nl0)
         offset = offset+geom%nl0
      end if
      if (nam%lhomv) then
         diag_blk%fit_rv = minim%x(offset+1)
         offset = offset+1
         if (diag_blk%double_fit) then
            diag_blk%fit_rv_rfac = minim%x(offset+1)
            offset = offset+1
            diag_blk%fit_rv_coef = minim%x(offset+1)
            offset = offset+1
         end if
      else
         diag_blk%fit_rv = minim%x(offset+1:offset+geom%nl0)
         offset = offset+geom%nl0
         if (diag_blk%double_fit) then
            diag_blk%fit_rv_rfac = minim%x(offset+1:offset+geom%nl0)
            offset = offset+geom%nl0
            diag_blk%fit_rv_coef = minim%x(offset+1:offset+geom%nl0)
            offset = offset+geom%nl0
         end if
      end if

      ! Release memory
      deallocate(minim%x)
      deallocate(minim%guess)
      deallocate(minim%binf)
      deallocate(minim%bsup)
      deallocate(minim%obs)
      deallocate(minim%l0rl0_to_l0)
      deallocate(minim%disth)
      deallocate(minim%distv)
   end select
end if

! Set to missing values if no point available
do il0=1,geom%nl0
   if (mpl%msv%isallr(diag_blk%raw(:,:,il0))) then
      diag_blk%fit_rh(il0) = mpl%msv%valr
      diag_blk%fit_rv(il0) = mpl%msv%valr
      if (diag_blk%double_fit) then
         diag_blk%fit_rh(il0) = mpl%msv%valr
         diag_blk%fit_rv(il0) = mpl%msv%valr
      end if
   end if
end do

! Release memory
deallocate(rawv)
deallocate(distv)
deallocate(fit)

! End associate
end associate

end subroutine diag_blk_fitting

!----------------------------------------------------------------------
! Subroutine: diag_blk_localization
! Purpose: diag_blk localization
!----------------------------------------------------------------------
subroutine diag_blk_localization(diag_blk,mpl,geom,bpar,avg_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk ! Diagnostic block (localization)
type(mpl_type),intent(inout) :: mpl            ! MPI data
type(geom_type),intent(in) :: geom             ! Geometry
type(bpar_type),intent(in) :: bpar             ! Block parameters
type(avg_blk_type),intent(in) :: avg_blk       ! Averaged statistics block

! Local variables
integer :: il0,jl0r,jc3

! Associate
associate(ib=>diag_blk%ib)

!$omp parallel do schedule(static) private(il0,jl0r,jc3)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (mpl%msv%isnotr(avg_blk%m11asysq(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m11sq(jc3,jl0r,il0))) then
            ! Compute localization
            diag_blk%raw(jc3,jl0r,il0) = avg_blk%m11asysq(jc3,jl0r,il0)/avg_blk%m11sq(jc3,jl0r,il0)
            diag_blk%valid(jc3,jl0r,il0) = avg_blk%nc1a(jc3,jl0r,il0)
         else
            ! Missing value
            diag_blk%raw(jc3,jl0r,il0) = mpl%msv%valr
            diag_blk%valid(jc3,jl0r,il0) = mpl%msv%valr
         end if
      end do
   end do
end do
!$omp end parallel do

! Hybrid weight
diag_blk%raw_coef_sta = mpl%msv%valr

! End associate
end associate

end subroutine diag_blk_localization

!----------------------------------------------------------------------
! Subroutine: diag_blk_hybridization
! Purpose: diag_blk hybridization
!----------------------------------------------------------------------
subroutine diag_blk_hybridization(diag_blk,mpl,geom,bpar,avg_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk ! Diagnostic block (localization)
type(mpl_type),intent(inout) :: mpl            ! MPI data
type(geom_type),intent(in) :: geom             ! Geometry
type(bpar_type),intent(in) :: bpar             ! Block parameters
type(avg_blk_type),intent(in) :: avg_blk       ! Averaged statistics block

! Local variables
integer :: il0,jl0r,jl0,jc3
real(kind_real) :: wgt,num,den

! Associate
associate(ib=>diag_blk%ib)

! Compute raw hybridization
num = 0.0
den = 0.0
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
      do jc3=1,bpar%nc3(ib)
         if (mpl%msv%isnotr(avg_blk%m11asysq(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m11sq(jc3,jl0r,il0)) &
       & .and.mpl%msv%isnotr(avg_blk%m11sta(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%stasq(jc3,jl0r,il0))) then
            wgt = geom%disth(jc3)*diag_blk%distv(jl0,il0)/real(bpar%nl0r(ib)+bpar%nc3(ib),kind_real)
            num = num+wgt*(1.0-avg_blk%m11asysq(jc3,jl0r,il0)/avg_blk%m11sq(jc3,jl0r,il0))*avg_blk%m11sta(jc3,jl0r,il0)
            den = den+wgt*(avg_blk%stasq(jc3,jl0r,il0)-avg_blk%m11sta(jc3,jl0r,il0)**2/avg_blk%m11sq(jc3,jl0r,il0))
         end if
      end do
   end do
end do
if ((num>0.0).and.(den>0.0)) then
   ! Valid numerator and denominator
   diag_blk%raw_coef_sta = num/den

   !$omp parallel do schedule(static) private(il0,jl0r,jc3)
   do il0=1,geom%nl0
      do jl0r=1,bpar%nl0r(ib)
         do jc3=1,bpar%nc3(ib)
            if (mpl%msv%isnotr(avg_blk%m11asysq(jc3,jl0r,il0)).and.mpl%msv%isnotr(diag_blk%raw_coef_sta) &
          & .and.mpl%msv%isnotr(avg_blk%m11sta(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m11sq(jc3,jl0r,il0))) then
               ! Compute localization
               diag_blk%raw(jc3,jl0r,il0) = (avg_blk%m11asysq(jc3,jl0r,il0)-diag_blk%raw_coef_sta &
                                          & *avg_blk%m11sta(jc3,jl0r,il0))/avg_blk%m11sq(jc3,jl0r,il0)
               diag_blk%valid(jc3,jl0r,il0) = avg_blk%nc1a(jc3,jl0r,il0)

               ! Lower bound
               if (diag_blk%raw(jc3,jl0r,il0)<0.0) then
                  diag_blk%raw(jc3,jl0r,il0) = mpl%msv%valr
                  diag_blk%valid(jc3,jl0r,il0) = mpl%msv%valr
               end if
            else
               ! Missing value
               diag_blk%raw(jc3,jl0r,il0) = mpl%msv%valr
               diag_blk%valid(jc3,jl0r,il0) = mpl%msv%valr
            end if
         end do
      end do
   end do
   !$omp end parallel do
else
   ! Missing values
   diag_blk%raw_coef_sta = mpl%msv%valr
   diag_blk%raw = mpl%msv%valr
end if

! End associate
end associate

end subroutine diag_blk_hybridization

!----------------------------------------------------------------------
! Subroutine: diag_blk_dualens
! Purpose: diag_blk dualens
!----------------------------------------------------------------------
subroutine diag_blk_dualens(diag_blk,mpl,geom,bpar,avg_blk,avg_lr_blk,diag_lr_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk   ! Diagnostic block (localization)
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(geom_type),intent(in) :: geom               ! Geometry
type(bpar_type),intent(in) :: bpar               ! Block parameters
type(avg_blk_type),intent(in) :: avg_blk         ! Averaged statistics block
type(avg_blk_type),intent(in) :: avg_lr_blk      ! LR averaged statistics block
type(diag_blk_type),intent(inout) :: diag_lr_blk ! Diagnostic block (LR localization)

! Local variables
integer :: il0,jl0r,jc3
real(kind_real),allocatable :: num(:),num_lr(:),den(:)

! Associate
associate(ib=>diag_blk%ib)

! Allocation
allocate(num(bpar%nc3(ib)))
allocate(num_lr(bpar%nc3(ib)))
allocate(den(bpar%nc3(ib)))

! Compute raw dual-ensemble hybridization
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (mpl%msv%isnotr(avg_blk%m11asysq(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m11sq(jc3,jl0r,il0)) &
       & .and.mpl%msv%isnotr(avg_blk%m11lrm11asy(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m11lrm11(jc3,jl0r,il0)) &
       & .and.mpl%msv%isnotr(avg_lr_blk%m11sq(jc3,jl0r,il0)).and.mpl%msv%isnotr(avg_blk%m11lrm11asy(jc3,jl0r,il0))) then
            num(jc3) = avg_blk%m11asysq(jc3,jl0r,il0)*avg_lr_blk%m11sq(jc3,jl0r,il0) &
                    & -avg_blk%m11lrm11asy(jc3,jl0r,il0)*avg_blk%m11lrm11(jc3,jl0r,il0)
            num_lr(jc3) = avg_blk%m11lrm11asy(jc3,jl0r,il0)*avg_blk%m11sq(jc3,jl0r,il0) &
                       & -avg_blk%m11asysq(jc3,jl0r,il0)*avg_blk%m11lrm11(jc3,jl0r,il0)
            den(jc3) = avg_blk%m11sq(jc3,jl0r,il0)*avg_lr_blk%m11sq(jc3,jl0r,il0)-avg_blk%m11lrm11(jc3,jl0r,il0)**2
            if ((num(jc3)>0.0).and.(den(jc3)>0.0)) then
               ! Compute localization
               diag_blk%raw(jc3,jl0r,il0) = num(jc3)/den(jc3)
               diag_lr_blk%raw(jc3,jl0r,il0) = num_lr(jc3)/den(jc3)
               diag_blk%valid(jc3,jl0r,il0) = avg_blk%nc1a(jc3,jl0r,il0)
               diag_lr_blk%valid(jc3,jl0r,il0) = avg_blk%nc1a(jc3,jl0r,il0)
            else
               ! Missing value
               diag_blk%raw(jc3,jl0r,il0) = mpl%msv%valr
               diag_lr_blk%raw(jc3,jl0r,il0) = mpl%msv%valr
               diag_blk%valid(jc3,jl0r,il0) = mpl%msv%valr
               diag_lr_blk%valid(jc3,jl0r,il0) = mpl%msv%valr
            end if
         end if
      end do
   end do
end do

! Hybrid weight
diag_blk%raw_coef_sta = mpl%msv%valr

! Release memory
deallocate(num)
deallocate(num_lr)
deallocate(den)

! End associate
end associate

end subroutine diag_blk_dualens

end module type_diag_blk
