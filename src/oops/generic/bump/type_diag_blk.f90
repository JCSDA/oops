!----------------------------------------------------------------------
! Module: type_diag_blk
!> Purpose: diagnostic block derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_diag_blk

use netcdf
!$ use omp_lib
use tools_const, only: msvali,msvalr
use tools_display, only: msgerror
use tools_fit, only: fast_fit,ver_smooth
use tools_func, only: fit_diag
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_avg_blk, only: avg_blk_type
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_minim, only: minim_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! Diagnostic block derived type
type diag_blk_type
   integer :: ic2a                                !< Local index
   integer :: ib                                  !< Block index
   character(len=1024) :: name                    !< Name
   integer :: npack                               !< Pack buffer size

   real(kind_real),allocatable :: raw(:,:,:)      !< Raw diagnostic
   real(kind_real),allocatable :: raw_coef_ens(:) !< Raw ensemble coefficient
   real(kind_real) :: raw_coef_sta                !< Raw static coefficient
   real(kind_real),allocatable :: fit(:,:,:)      !< Fit
   real(kind_real),allocatable :: fit_rh(:)       !< Fit support radius
   real(kind_real),allocatable :: fit_rv(:)       !< Fit support radius
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

integer,parameter :: nsc = 50 !< Scaling optimization parameter
logical :: lprt = .false.     !< Optimization print

private
public :: diag_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: diag_blk_alloc
!> Purpose: diagnostic block object allocation
!----------------------------------------------------------------------
subroutine diag_blk_alloc(diag_blk,nam,geom,bpar,ic2a,ib,prefix)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk !< Diagnostic block
type(nam_type),intent(in) :: nam               !< Namelist
type(geom_type),intent(in) :: geom             !< Geometry
type(bpar_type),intent(in) :: bpar             !< Block parameters
integer,intent(in) :: ic2a                     !< Local index
integer,intent(in) :: ib                       !< Block index
character(len=*),intent(in) :: prefix          !< Block prefix

! Set attributes
diag_blk%ic2a = ic2a
diag_blk%ib = ib
diag_blk%name = trim(prefix)//'_'//trim(bpar%blockname(ib))
diag_blk%npack = nam%nc3*geom%nl0*bpar%nl0r(ib)+2*geom%nl0

! Allocation
allocate(diag_blk%raw(nam%nc3,bpar%nl0r(ib),geom%nl0))
allocate(diag_blk%raw_coef_ens(geom%nl0))

! Initialization
call msr(diag_blk%raw)
call msr(diag_blk%raw_coef_ens)
call msr(diag_blk%raw_coef_sta)

if (trim(nam%minim_algo)/='none') then
   ! Allocation
   allocate(diag_blk%fit(nam%nc3,bpar%nl0r(ib),geom%nl0))
   allocate(diag_blk%fit_rh(geom%nl0))
   allocate(diag_blk%fit_rv(geom%nl0))

   ! Initialization
   call msr(diag_blk%fit)
   call msr(diag_blk%fit_rh)
   call msr(diag_blk%fit_rv)
end if

end subroutine diag_blk_alloc

!----------------------------------------------------------------------
! Subroutine: diag_blk_dealloc
!> Purpose: diag block object deallocation
!----------------------------------------------------------------------
subroutine diag_blk_dealloc(diag_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk !< Diagnostic block

! Deallocation
if (allocated(diag_blk%raw)) deallocate(diag_blk%raw)
if (allocated(diag_blk%raw_coef_ens)) deallocate(diag_blk%raw_coef_ens)
if (allocated(diag_blk%fit)) deallocate(diag_blk%fit)
if (allocated(diag_blk%fit_rh)) deallocate(diag_blk%fit_rh)
if (allocated(diag_blk%fit_rv)) deallocate(diag_blk%fit_rv)

end subroutine diag_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: diag_blk_write
!> Purpose: write a diagnostic
!----------------------------------------------------------------------
subroutine diag_blk_write(diag_blk,nam,geom,bpar,filename)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk !< Diagnostic block
type(nam_type),intent(in) :: nam               !< Namelist
type(geom_type),intent(in) :: geom             !< Geometry
type(bpar_type),intent(in) :: bpar             !< Block parameters
character(len=*),intent(in) :: filename        !< File name

! Local variables
integer :: info,ncid,one_id,nc3_id,nl0r_id,nl0_id,disth_id,vunit_id
integer :: raw_id,raw_coef_ens_id,raw_coef_sta_id,l0rl0_to_l0_id
integer :: fit_id,fit_rh_id,fit_rv_id
character(len=1024) :: subr = 'diag_blk_write'

! Associate
associate(ib=>diag_blk%ib)

! Check if the file exists
info = nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_noclobber,nf90_64bit_offset),ncid)
if (info==nf90_noerr) then
   ! Write namelist parameters
   call nam%ncwrite(ncid)

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,'one',1,one_id))
   call ncerr(subr,nf90_def_dim(ncid,'nc3',nam%nc3,nc3_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl0r',nam%nl0r,nl0r_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,'disth',ncfloat,(/nc3_id/),disth_id))
   call ncerr(subr,nf90_def_var(ncid,'vunit',ncfloat,(/nl0_id/),vunit_id))
else
   ! Open file
   call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_write,ncid))

   ! Get dimensions ID
   call ncerr(subr,nf90_inq_dimid(ncid,'one',one_id))
   call ncerr(subr,nf90_inq_dimid(ncid,'nc3',nc3_id))
   call ncerr(subr,nf90_inq_dimid(ncid,'nl0r',nl0r_id))
   call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))

   ! Redef mode
   call ncerr(subr,nf90_redef(ncid))
end if

! Define variables
call ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_raw',ncfloat,(/nc3_id,nl0r_id,nl0_id/),raw_id))
call ncerr(subr,nf90_put_att(ncid,raw_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_raw_coef_ens',ncfloat,(/nl0_id/),raw_coef_ens_id))
call ncerr(subr,nf90_put_att(ncid,raw_coef_ens_id,'_FillValue',msvalr))
if (isnotmsr(diag_blk%raw_coef_sta)) then
   call ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_raw_coef_sta',ncfloat,(/one_id/),raw_coef_sta_id))
   call ncerr(subr,nf90_put_att(ncid,raw_coef_sta_id,'_FillValue',msvalr))
end if
if ((trim(nam%minim_algo)/='none').and.(isanynotmsr(diag_blk%fit))) then
   call ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit',ncfloat,(/nc3_id,nl0r_id,nl0_id/),fit_id))
   call ncerr(subr,nf90_put_att(ncid,fit_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit_rh',ncfloat,(/nl0_id/),fit_rh_id))
   call ncerr(subr,nf90_put_att(ncid,fit_rh_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_fit_rv',ncfloat,(/nl0_id/),fit_rv_id))
   call ncerr(subr,nf90_put_att(ncid,fit_rv_id,'_FillValue',msvalr))
end if
call ncerr(subr,nf90_def_var(ncid,trim(diag_blk%name)//'_l0rl0_to_l0',nf90_int,(/nl0r_id,nl0_id/),l0rl0_to_l0_id))
call ncerr(subr,nf90_put_att(ncid,l0rl0_to_l0_id,'_FillValue',msvali))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write variables
if (info==nf90_noerr) then
   call ncerr(subr,nf90_put_var(ncid,disth_id,geom%disth(1:nam%nc3)))
   call ncerr(subr,nf90_put_var(ncid,vunit_id,sum(geom%vunit,mask=geom%mask,dim=1)/real(count(geom%mask,dim=1),kind_real)))
end if
call ncerr(subr,nf90_put_var(ncid,raw_id,diag_blk%raw))
call ncerr(subr,nf90_put_var(ncid,raw_coef_ens_id,diag_blk%raw_coef_ens))
if (isnotmsr(diag_blk%raw_coef_sta)) call ncerr(subr,nf90_put_var(ncid,raw_coef_sta_id,diag_blk%raw_coef_sta))
if ((trim(nam%minim_algo)/='none').and.(isanynotmsr(diag_blk%fit))) then
   call ncerr(subr,nf90_put_var(ncid,fit_id,diag_blk%fit))
   call ncerr(subr,nf90_put_var(ncid,fit_rh_id,diag_blk%fit_rh))
   call ncerr(subr,nf90_put_var(ncid,fit_rv_id,diag_blk%fit_rv))
end if
call ncerr(subr,nf90_put_var(ncid,l0rl0_to_l0_id,bpar%l0rl0b_to_l0(:,:,ib)))

! Close file
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine diag_blk_write

!----------------------------------------------------------------------
! Subroutine: diag_blk_normalization
!> Purpose: compute diagnostic block normalization
!----------------------------------------------------------------------
subroutine diag_blk_normalization(diag_blk,geom,bpar)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk !< Diagnostic block
type(geom_type),intent(in) :: geom             !< Geometry
type(bpar_type),intent(in) :: bpar             !< Block parameters

! Local variables
integer :: il0r,il0,jl0,jc3

! Associate
associate(ib=>diag_blk%ib)

! Get diagonal values
do jl0=1,geom%nl0
   il0r = bpar%il0rz(jl0,ib)
   if (isnotmsr(diag_blk%raw(1,il0r,jl0))) diag_blk%raw_coef_ens(jl0) = diag_blk%raw(1,il0r,jl0)
end do

! Normalize
do jl0=1,geom%nl0
   do il0r=1,bpar%nl0r(ib)
      il0 = bpar%l0rl0b_to_l0(il0r,jl0,ib)
      do jc3=1,bpar%nc3(ib)
         if (isnotmsr(diag_blk%raw(jc3,il0r,jl0)).and.isnotmsr(diag_blk%raw_coef_ens(il0)) &
       & .and.isnotmsr(diag_blk%raw_coef_ens(jl0))) &
       & diag_blk%raw(jc3,il0r,jl0) = diag_blk%raw(jc3,il0r,jl0)/sqrt(diag_blk%raw_coef_ens(il0)*diag_blk%raw_coef_ens(jl0))
      end do
   end do
end do

! End associate
end associate

end subroutine diag_blk_normalization

!----------------------------------------------------------------------
! Subroutine: diag_blk_fitting
!> Purpose: compute a semi-positive definite fit of a raw function
!----------------------------------------------------------------------
subroutine diag_blk_fitting(diag_blk,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk !< Diagnostic block
type(nam_type),intent(in) :: nam               !< Namelist
type(geom_type),intent(in) :: geom             !< Geometry
type(bpar_type),intent(in) :: bpar             !< Block parameters
type(hdata_type),intent(in) :: hdata           !< HDIAG data
! Local variables
integer :: ic2,ic0,il0,jl0r,jl0,offset,isc
real(kind_real) :: distvr(nam%nl0r,geom%nl0),rawv(nam%nl0r)
real(kind_real) :: alpha,alpha_opt,mse,mse_opt
real(kind_real) :: fit_rh(geom%nl0),fit_rv(geom%nl0),fit(nam%nc3,nam%nl0r,geom%nl0)
type(minim_type) :: minim

! Associate
associate(ic2a=>diag_blk%ic2a,ib=>diag_blk%ib)

! Check
if (trim(nam%minim_algo)=='none') call msgerror('cannot compute fit if minim_algo = none')

! Initialization
call msr(diag_blk%fit_rh)
call msr(diag_blk%fit_rv)
call msr(diag_blk%fit)

! Reduced vertical distance
call msr(distvr)
do il0=1,geom%nl0
   do jl0r=1,nam%nl0r
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
      if (ic2a==0) then
         distvr(jl0,il0) = abs(geom%vunitavg(il0)-geom%vunitavg(jl0))
      else
         ic2 = hdata%c2a_to_c2(ic2a)
         ic0 = hdata%c2_to_c0(ic2)
         distvr(jl0,il0) = abs(geom%vunit(ic0,il0)-geom%vunit(ic0,jl0))
      end if
   end do
end do

! Fast fit
do il0=1,geom%nl0
   ! Get zero separation level
   jl0r = bpar%il0rz(il0,ib)

   ! Horizontal fast fit
   call fast_fit(nam%nc3,1,geom%disth,diag_blk%raw(:,jl0r,il0),diag_blk%fit_rh(il0))

   ! Vertical fast fit
   rawv = diag_blk%raw(1,:,il0)
   call fast_fit(nam%nl0r,jl0r,distvr(:,il0),rawv,diag_blk%fit_rv(il0))
end do
if (nam%lhomh) diag_blk%fit_rh = sum(diag_blk%fit_rh,mask=isnotmsr(diag_blk%fit_rh)) &
 & /real(count(isnotmsr(diag_blk%fit_rh)),kind_real)
if (nam%lhomv) diag_blk%fit_rv = sum(diag_blk%fit_rv,mask=isnotmsr(diag_blk%fit_rv)) &
 & /real(count(isnotmsr(diag_blk%fit_rv)),kind_real)

! Scaling optimization (brute-force)
mse_opt = huge(1.0)
alpha_opt = 1.0
do isc=1,nsc
   ! Scaling factor
   alpha = 0.5+real(isc-1,kind_real)/real(nsc-1,kind_real)*(2.0-0.5)

   ! Scaled radii
   fit_rh = alpha*diag_blk%fit_rh
   fit_rv = alpha*diag_blk%fit_rv

   ! Define fit
   call fit_diag(nam%nc3,nam%nl0r,geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth,distvr,fit_rh,fit_rv,fit)

   ! MSE
   mse = sum((fit-diag_blk%raw)**2,mask=isnotmsr(diag_blk%raw))
   if (mse<mse_opt) then
      mse_opt = mse
      alpha_opt = alpha
   end if
end do
do il0=1,geom%nl0
   if (isnotmsr(diag_blk%fit_rh(il0))) diag_blk%fit_rh(il0) = alpha_opt*diag_blk%fit_rh(il0)
   if (isnotmsr(diag_blk%fit_rv(il0))) diag_blk%fit_rv(il0) = alpha_opt*diag_blk%fit_rv(il0)
end do
if (lprt) then
   write(mpl%unit,'(a)') ''
   write(mpl%unit,'(a13,a,f6.1,a)') '','Scaling optimization, cost function decrease:',abs(mse_opt-mse)/mse*100.0,'%'
   call flush(mpl%unit)
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
   else
      minim%nx = minim%nx+geom%nl0
   end if
   minim%ny = nam%nc3*nam%nl0r*geom%nl0
   allocate(minim%x(minim%nx))
   allocate(minim%guess(minim%nx))
   allocate(minim%norm(minim%nx))
   allocate(minim%binf(minim%nx))
   allocate(minim%bsup(minim%nx))
   allocate(minim%obs(minim%ny))
   allocate(minim%l0rl0_to_l0(nam%nl0r,geom%nl0))
   allocate(minim%disth(nam%nc3))
   allocate(minim%distvr(nam%nl0r,geom%nl0))

   ! Fill minim
   offset = 0
   if (nam%lhomh) then
      minim%guess(offset+1) = diag_blk%fit_rh(1)
      offset = offset+1
   else
      minim%guess(offset+1:offset+geom%nl0) = diag_blk%fit_rh
      offset = offset+geom%nl0
   end if
   if (nam%lhomv) then
      minim%guess(offset+1) = diag_blk%fit_rv(1)
      offset = offset+1
   else
      minim%guess(offset+1:offset+geom%nl0) = diag_blk%fit_rv
      offset = offset+geom%nl0
   end if
   minim%norm = minim%guess
   minim%binf = 0.5*minim%guess
   minim%bsup = 1.5*minim%guess
   minim%obs = pack(diag_blk%raw,mask=.true.)
   minim%cost_function = 'fit_diag'
   minim%algo = nam%minim_algo
   minim%nc3 = nam%nc3
   minim%nl0r = nam%nl0r
   minim%nl0 = geom%nl0
   minim%lhomh = nam%lhomh
   minim%lhomv = nam%lhomv
   minim%l0rl0_to_l0 = bpar%l0rl0b_to_l0(:,:,ib)
   minim%disth = geom%disth
   minim%distvr = distvr

   ! Compute fit
   minim%cost_function = 'fit_diag'
   call minim%compute(lprt)

   ! Apply bounds
   minim%x = max(minim%binf,min(minim%x,minim%bsup))

   ! Copy parameters
   offset = 0
   if (nam%lhomh) then
      diag_blk%fit_rh = minim%x(offset+1)
      offset = offset+1
   else
      do il0=1,geom%nl0
         if (isnotmsr(diag_blk%fit_rh(il0))) diag_blk%fit_rh(il0) = minim%x(offset+il0)
      end do
      offset = offset+geom%nl0
   end if
   if (nam%lhomv) then
      diag_blk%fit_rv = minim%x(offset+1)
      offset = offset+1
   else
      do il0=1,geom%nl0
         if (isnotmsr(diag_blk%fit_rv(il0))) diag_blk%fit_rv(il0) = minim%x(offset+il0)
      end do
      offset = offset+geom%nl0
   end if
end select

! Smooth vertically
call ver_smooth(geom%nl0,geom%vunitavg,nam%rvflt,diag_blk%fit_rh)
call ver_smooth(geom%nl0,geom%vunitavg,nam%rvflt,diag_blk%fit_rv)

! Rebuild fit
call fit_diag(nam%nc3,nam%nl0r,geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth,distvr,diag_blk%fit_rh,diag_blk%fit_rv,diag_blk%fit)

! End associate
end associate

end subroutine diag_blk_fitting

!----------------------------------------------------------------------
! Subroutine: diag_blk_localization
!> Purpose: diag_blk localization
!----------------------------------------------------------------------
subroutine diag_blk_localization(diag_blk,geom,bpar,avg_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk !< Diagnostic block (localization)
type(geom_type),intent(in) :: geom             !< Geometry
type(bpar_type),intent(in) :: bpar             !< Block parameters
type(avg_blk_type),intent(in) :: avg_blk       !< Averaged statistics block

! Local variables
integer :: il0,jl0r,jc3

! Associate
associate(ib=>diag_blk%ib)

! Compute raw localization
!$omp parallel do schedule(static) private(il0,jl0r,jc3)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (isnotmsr(avg_blk%m11asysq(jc3,jl0r,il0)).and.isnotmsr(avg_blk%m11sq(jc3,jl0r,il0))) &
       & diag_blk%raw(jc3,jl0r,il0) = avg_blk%m11asysq(jc3,jl0r,il0)/avg_blk%m11sq(jc3,jl0r,il0)
      end do
   end do
end do
!$omp end parallel do

! End associate
end associate

end subroutine diag_blk_localization

!----------------------------------------------------------------------
! Subroutine: diag_blk_hybridization
!> Purpose: diag_blk hybridization
!----------------------------------------------------------------------
subroutine diag_blk_hybridization(diag_blk,geom,bpar,avg_blk,avg_sta_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk !< Diagnostic block (localization)
type(geom_type),intent(in) :: geom             !< Geometry
type(bpar_type),intent(in) :: bpar             !< Block parameters
type(avg_blk_type),intent(in) :: avg_blk       !< Averaged statistics block
type(avg_blk_type),intent(in) :: avg_sta_blk   !< Static averaged statistics block

! Local variables
integer :: il0,jl0r,jc3
real(kind_real) :: num,den

! Associate
associate(ib=>diag_blk%ib)

! Compute raw hybridization
num = 0.0
den = 0.0
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (isnotmsr(avg_blk%m11asysq(jc3,jl0r,il0)).and.isnotmsr(avg_blk%m11sq(jc3,jl0r,il0)) &
       & .and.isnotmsr(avg_sta_blk%m11sta(jc3,jl0r,il0)).and.isnotmsr(avg_sta_blk%stasq(jc3,jl0r,il0))) then
            num = num+geom%disth(jc3)*(1.0-avg_blk%m11asysq(jc3,jl0r,il0)/avg_blk%m11sq(jc3,jl0r,il0)) &
                & *avg_sta_blk%m11sta(jc3,jl0r,il0)
            den = den+geom%disth(jc3)*(avg_sta_blk%stasq(jc3,jl0r,il0)-avg_sta_blk%m11sta(jc3,jl0r,il0)**2 &
                & /avg_blk%m11sq(jc3,jl0r,il0))
         end if
      end do
   end do
end do
if ((num>0.0).and.(den>0.0)) diag_blk%raw_coef_sta = num/den
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (isnotmsr(avg_blk%m11asysq(jc3,jl0r,il0)).and.isnotmsr(diag_blk%raw_coef_sta) &
       & .and.isnotmsr(avg_sta_blk%m11sta(jc3,jl0r,il0)).and.isnotmsr(avg_sta_blk%stasq(jc3,jl0r,il0)) &
       & .and.isnotmsr(avg_blk%m11sq(jc3,jl0r,il0))) then
            ! Compute localization
            diag_blk%raw(jc3,jl0r,il0) = (avg_blk%m11asysq(jc3,jl0r,il0)-diag_blk%raw_coef_sta &
                                    & *avg_sta_blk%m11sta(jc3,jl0r,il0))/avg_blk%m11sq(jc3,jl0r,il0)

            ! Lower bound
            if (diag_blk%raw(jc3,jl0r,il0)<0.0) call msr(diag_blk%raw(jc3,jl0r,il0))
         end if
      end do
   end do
end do

! End associate
end associate

end subroutine diag_blk_hybridization

!----------------------------------------------------------------------
! Subroutine: diag_blk_dualens
!> Purpose: diag_blk dualens
!----------------------------------------------------------------------
subroutine diag_blk_dualens(diag_blk,geom,bpar,avg_blk,avg_lr_blk,diag_lr_blk)

implicit none

! Passed variables
class(diag_blk_type),intent(inout) :: diag_blk   !< Diagnostic block (localization)
type(geom_type),intent(in) :: geom               !< Geometry
type(bpar_type),intent(in) :: bpar               !< Block parameters
type(avg_blk_type),intent(in) :: avg_blk         !< Averaged statistics block
type(avg_blk_type),intent(in) :: avg_lr_blk      !< LR averaged statistics block
type(diag_blk_type),intent(inout) :: diag_lr_blk !< Diagnostic block (LR localization)

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
         if (isnotmsr(avg_blk%m11asysq(jc3,jl0r,il0)).and.isnotmsr(avg_blk%m11sq(jc3,jl0r,il0)) &
       & .and.isnotmsr(avg_lr_blk%m11lrm11asy(jc3,jl0r,il0)).and.isnotmsr(avg_lr_blk%m11lrm11(jc3,jl0r,il0)) &
       & .and.isnotmsr(avg_lr_blk%m11sq(jc3,jl0r,il0)).and.isnotmsr(avg_lr_blk%m11lrm11asy(jc3,jl0r,il0))) then
            num(jc3) = avg_blk%m11asysq(jc3,jl0r,il0)*avg_lr_blk%m11sq(jc3,jl0r,il0) &
                    & -avg_lr_blk%m11lrm11asy(jc3,jl0r,il0)*avg_lr_blk%m11lrm11(jc3,jl0r,il0)
            num_lr(jc3) = avg_lr_blk%m11lrm11asy(jc3,jl0r,il0)*avg_blk%m11sq(jc3,jl0r,il0) &
                       & -avg_blk%m11asysq(jc3,jl0r,il0)*avg_lr_blk%m11lrm11(jc3,jl0r,il0)
            den(jc3) = avg_blk%m11sq(jc3,jl0r,il0)*avg_lr_blk%m11sq(jc3,jl0r,il0)-avg_lr_blk%m11lrm11(jc3,jl0r,il0)**2
            if ((num(jc3)>0.0).and.(den(jc3)>0.0)) then
               diag_blk%raw(jc3,jl0r,il0) = num(jc3)/den(jc3)
               diag_lr_blk%raw(jc3,jl0r,il0) = num_lr(jc3)/den(jc3)
            end if
         end if
      end do
   end do
end do

! End associate
end associate

end subroutine diag_blk_dualens

end module type_diag_blk
