!----------------------------------------------------------------------
! Module: nicas_apply_localization.f90
!> Purpose: apply localization
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_apply_localization

use nicas_apply_nicas, only: apply_nicas,apply_nicas_sqrt,apply_nicas_sqrt_ad
use nicas_apply_adv, only: apply_adv,apply_adv_ad
use tools_asa007, only: cholesky
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr
use type_bpar, only: bpartype
use type_cv, only: cvtype,cv_alloc,cv_random
use type_geom, only: geomtype
use type_mpl, only: mpl
use type_nam, only: namtype
use type_ndata, only: ndatatype
use yomhook, only: lhook,dr_hook

implicit none

logical :: lcoef_ens = .false. !< Apply ensemble coefficient (will reduce variance)

private
public :: apply_localization,apply_localization_sqrt,apply_localization_sqrt_ad,apply_localization_from_sqrt
public :: randomize_localization

contains

!----------------------------------------------------------------------
! Subroutine: apply_localization
!> Purpose: apply localization
!----------------------------------------------------------------------
subroutine apply_localization(nam,geom,bpar,ndata,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                !< Geometry
type(bpartype),target,intent(in) :: bpar                                !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                          !< NICAS data
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ib,its,iv,jv,il0,ic0a
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:)

if (lhook) call dr_hook('apply_localization',0,zhook_handle)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Adjoint displacement
   if (nam%displ_diag) call apply_adv_ad(nam,geom,ndata(bpar%nb+1),fld)

   ! Sum product over variables and timeslots
   fld_3d = 0.0
   !$omp parallel do schedule(static) private(il0,ic0a,its,iv)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         do its=1,nam%nts
            do iv=1,nam%nv
               fld_3d(ic0a,il0) = fld_3d(ic0a,il0)+fld(ic0a,il0,iv,its)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   if (lcoef_ens) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(ndata(bpar%nb+1)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
   end if

   ! Apply common localization
   call apply_nicas(geom,ndata(bpar%nb+1),fld_3d)

   if (lcoef_ens) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(ndata(bpar%nb+1)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
   end if

   ! Build final vector
   !$omp parallel do schedule(static) private(il0,ic0a,its,iv)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         do its=1,nam%nts
            do iv=1,nam%nv
               fld(ic0a,il0,iv,its) = fld_3d(ic0a,il0)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   ! Displacement
   if (nam%displ_diag) call apply_adv(nam,geom,ndata(bpar%nb+1),fld)
case ('specific_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas(geom,ndata(ib),fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('specific_multivariate')
   call msgerror('specific multivariate strategy should not be called from apply_localization')
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))

   ! Copy weights
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         if (isnotmsr(ndata(ib)%wgt)) then
            wgt(iv,jv) = ndata(ib)%wgt
         else
            wgt(iv,jv) = 0.0
         end if
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do iv=1,nam%nv
      ! Apply common ensemble coefficient square-root
      if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(bpar%nb+1)%coef_ens)

      ! Apply common localization
      call apply_nicas(geom,ndata(bpar%nb+1),fld_4d(:,:,iv))

      ! Apply common ensemble coefficient square-root
      if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(bpar%nb+1)%coef_ens)
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=1,nam%nv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d_tmp
   end do
end select

if (lhook) call dr_hook('apply_localization',1,zhook_handle)

end subroutine apply_localization

!----------------------------------------------------------------------
! Subroutine: apply_localization_sqrt
!> Purpose: apply localization square-root
!----------------------------------------------------------------------
subroutine apply_localization_sqrt(nam,geom,bpar,ndata,cv,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                !< Namelist
type(geomtype),target,intent(in) :: geom                              !< Geometry
type(bpartype),target,intent(in) :: bpar                              !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                        !< NICAS data
type(cvtype),intent(in) :: cv(bpar%nb+1)                              !< Control variable
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ib,its,iv,jv,i,nullty,info
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),a(:),u(:)

if (lhook) call dr_hook('apply_localization_sqrt',0,zhook_handle)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Apply common localization
   call apply_nicas_sqrt(geom,ndata(bpar%nb+1),cv(bpar%nb+1)%alpha,fld_3d)

   ! Apply common ensemble coefficient square-root
   if (lcoef_ens) fld_3d = fld_3d*sqrt(ndata(bpar%nb+1)%coef_ens)

   ! Build final vector
   do its=1,nam%nts
      do iv=1,nam%nv
         fld(:,:,iv,its) = fld_3d
      end do
   end do
case ('specific_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt(geom,ndata(ib),cv(ib)%alpha,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('specific_multivariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt(geom,ndata(ib),cv(bpar%nb+1)%alpha,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))
   allocate(a((nam%nv*(nam%nv+1))/2))
   allocate(u((nam%nv*(nam%nv+1))/2))

   ! Copy weights
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         if (isnotmsr(ndata(ib)%wgt)) then
            wgt(iv,jv) = ndata(ib)%wgt
         else
            wgt(iv,jv) = 0.0
         end if
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Cholesky decomposition
   i = 0
   do iv=1,nam%nv
      do jv=1,iv
         i = i+1
         a(i) = wgt(iv,jv)
      end do
   end do
   call cholesky(a,nam%nv,(nam%nv*(nam%nv+1))/2,u,nullty,info)
   if (info==2) call msgerror('weights are not positive definite in apply_localization')
   i = 0
   wgt = 0.0
   do iv=1,nam%nv
      do jv=1,iv
         i = i+1
         wgt(iv,jv) = u(i)
      end do
   end do

   do ib=1,bpar%nb
      if (bpar%cv_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt(geom,ndata(bpar%nb+1),cv(ib)%alpha,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(bpar%nb+1)%coef_ens)
      end if
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=1,iv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d_tmp
   end do
end select

if (lhook) call dr_hook('apply_localization_sqrt',1,zhook_handle)

end subroutine apply_localization_sqrt

!----------------------------------------------------------------------
! Subroutine: apply_localization_sqrt_ad
!> Purpose: apply localization square-root, adjoint
!----------------------------------------------------------------------
subroutine apply_localization_sqrt_ad(nam,geom,bpar,ndata,fld,cv)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                               !< Namelist
type(geomtype),target,intent(in) :: geom                             !< Geometry
type(bpartype),target,intent(in) :: bpar                             !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                       !< NICAS data
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field
type(cvtype),intent(out) :: cv(bpar%nb+1)                            !< Control variable

! Local variable
integer :: ib,its,iv,jv,i,nullty,info
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),a(:),u(:)
type(cvtype) :: cv_tmp(bpar%nb+1)

if (lhook) call dr_hook('apply_localization_sqrt_ad',0,zhook_handle)

! Allocation
call cv_alloc(bpar,ndata,cv)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Sum product over variables and timeslots
   fld_3d = 0.0
   do its=1,nam%nts
      do iv=1,nam%nv
         fld_3d = fld_3d+fld(:,:,iv,its)
      end do
   end do

   ! Apply common ensemble coefficient square-root
   if (lcoef_ens) fld_3d = fld_3d*sqrt(ndata(bpar%nb+1)%coef_ens)

   ! Apply common localization
   call apply_nicas_sqrt_ad(geom,ndata(bpar%nb+1),fld_3d,cv(bpar%nb+1)%alpha)
case ('specific_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt_ad(geom,ndata(ib),fld_4d(:,:,iv),cv(ib)%alpha)
      end if
   end do
case ('specific_multivariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   call cv_alloc(bpar,ndata,cv_tmp)

   ! Initialization
   cv(bpar%nb+1)%alpha = 0.0

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndata(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt_ad(geom,ndata(ib),fld_4d(:,:,iv),cv_tmp(bpar%nb+1)%alpha)

         ! Sum control variable
         cv(bpar%nb+1)%alpha = cv(bpar%nb+1)%alpha+cv_tmp(bpar%nb+1)%alpha
      end if
   end do
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))
   allocate(a((nam%nv*(nam%nv+1))/2))
   allocate(u((nam%nv*(nam%nv+1))/2))

   ! Copy weights
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         if (isnotmsr(ndata(ib)%wgt)) then
            wgt(iv,jv) = ndata(ib)%wgt
         else
            wgt(iv,jv) = 0.0
         end if
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Cholesky decomposition
   i = 0
   do iv=1,nam%nv
      do jv=1,iv
         i = i+1
         a(i) = wgt(iv,jv)
      end do
   end do
   call cholesky(a,nam%nv,(nam%nv*(nam%nv+1))/2,u,nullty,info)
   if (info==2) call msgerror('weights are not positive definite in apply_localization')
   i = 0
   wgt = 0.0
   do jv=1,nam%nv
      do iv=1,jv
         i = i+1
         wgt(iv,jv) = u(i)
      end do
   end do

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=iv,nam%nv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   do ib=1,bpar%nb
      if (bpar%cv_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)*sqrt(ndata(bpar%nb+1)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt_ad(geom,ndata(bpar%nb+1),fld_4d_tmp(:,:,iv),cv(ib)%alpha)
      end if
   end do
end select

if (lhook) call dr_hook('apply_localization_sqrt_ad',1,zhook_handle)

end subroutine apply_localization_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: apply_localization_from_sqrt
!> Purpose: apply localization from square-root
!----------------------------------------------------------------------
subroutine apply_localization_from_sqrt(nam,geom,bpar,ndata,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                !< Geometry
type(bpartype),target,intent(in) :: bpar                                !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                          !< NICAS data
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
real(kind_real) :: zhook_handle
type(cvtype) :: cv(bpar%nb+1)

if (lhook) call dr_hook('apply_localization_from_sqrt',0,zhook_handle)

! Apply square-root adjoint
call apply_localization_sqrt_ad(nam,geom,bpar,ndata,fld,cv)

! Apply square-root
call apply_localization_sqrt(nam,geom,bpar,ndata,cv,fld)

if (lhook) call dr_hook('apply_localization_from_sqrt',1,zhook_handle)

end subroutine apply_localization_from_sqrt

!----------------------------------------------------------------------
! Subroutine: randomize_localization
!> Purpose: randomize localization from square-root
!----------------------------------------------------------------------
subroutine randomize_localization(nam,geom,bpar,ndata,ne,ens)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                   !< Namelist
type(geomtype),target,intent(in) :: geom                                 !< Geometry
type(bpartype),target,intent(in) :: bpar                                 !< Blocal parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                           !< NICAS data
integer,intent(in) :: ne                                                 !< Number of members
real(kind_real),intent(out) :: ens(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne) !< Ensemble

! Local variable
integer :: ie
real(kind_real) :: mean(geom%nc0a,geom%nl0,nam%nv,nam%nts),std(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: zhook_handle
type(cvtype) :: cv(bpar%nb+1,ne)

if (lhook) call dr_hook('randomize_localization',0,zhook_handle)

! Random control vector
call cv_random(bpar,ndata,ne,cv)

do ie=1,ne
   ! Apply square-root
   call apply_localization_sqrt(nam,geom,bpar,ndata,cv(:,ie),ens(:,:,:,:,ie))
end do

! Normalize ensemble
mean = sum(ens,dim=5)/float(ne)
do ie=1,ne
   ens(:,:,:,:,ie) = ens(:,:,:,:,ie)-mean
end do
std = sqrt(sum(ens**2,dim=5)/float(ne-1))
do ie=1,ne
   ens(:,:,:,:,ie) = ens(:,:,:,:,ie)/std
end do

if (lhook) call dr_hook('randomize_localization',1,zhook_handle)

end subroutine randomize_localization

end module nicas_apply_localization
