!----------------------------------------------------------------------
! Module: module_apply_localization.f90
!> Purpose: apply 4D localization
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_localization

use module_apply_nicas, only: apply_nicas,apply_nicas_sqrt,apply_nicas_sqrt_ad
use tools_asa007, only: cholesky
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr
use type_bpar, only: bpartype
use type_cv, only: cvtype,cv_alloc,cv_random
use type_geom, only: geomtype
use type_mpl, only: mpl
use type_nam, only: namtype
use type_ndata, only: ndataloctype

implicit none

private
public :: apply_localization,apply_localization_sqrt,apply_localization_sqrt_ad,apply_localization_from_sqrt
public :: randomize_localization

contains

!----------------------------------------------------------------------
! Subroutine: apply_localization
!> Purpose: apply localization
!----------------------------------------------------------------------
subroutine apply_localization(nam,geom,bpar,ndataloc,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                !< Geometry
type(bpartype),target,intent(in) :: bpar                                !< Blocal parameters
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1)                    !< NICAS data, local
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field, local

! Local variable
integer :: ib,its,iv,jv,i,nullty,info
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),a(:),u(:)

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
   fld_3d = fld_3d*sqrt(ndataloc(bpar%nb+1)%coef_ens)

   ! Apply common localization
   call apply_nicas(nam,geom,ndataloc(bpar%nb+1),fld_3d)

   ! Apply common ensemble coefficient square-root
   fld_3d = fld_3d*sqrt(ndataloc(bpar%nb+1)%coef_ens)

   ! Build final vector
   do its=1,nam%nts
      do iv=1,nam%nv
         fld(:,:,iv,its) = fld_3d
      end do
   end do
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
         iv = bpar%ib_to_iv(ib)

         ! Apply common ensemble coefficient square-root
         fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndataloc(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas(nam,geom,ndataloc(ib),fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndataloc(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('specific_multivariate')
   call msgerror('not implemented yet in apply_localization')
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))
   allocate(a(nam%nv*(nam%nv+1)/2))
   allocate(u(nam%nv*(nam%nv+1)/2))

   ! Copy weights
   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         ! Variable indices
         iv = bpar%ib_to_iv(ib)
         jv = bpar%ib_to_jv(ib)
         if (isnotmsr(ndataloc(ib)%wgt)) then
            wgt(iv,jv) = ndataloc(ib)%wgt
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

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=1,nam%nv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do
   fld_4d = fld_4d_tmp

   do iv=1,nam%nv
      ! Apply common ensemble coefficient square-root
      fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndataloc(bpar%nb+1)%coef_ens)

      ! Apply common localization
      call apply_nicas(nam,geom,ndataloc(bpar%nb+1),fld_4d(:,:,iv))

      ! Apply common ensemble coefficient square-root
      fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndataloc(bpar%nb+1)%coef_ens)
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
end select

end subroutine apply_localization

!----------------------------------------------------------------------
! Subroutine: apply_localization_sqrt
!> Purpose: apply localization square-root
!----------------------------------------------------------------------
subroutine apply_localization_sqrt(nam,geom,bpar,ndataloc,cv,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                !< Namelist
type(geomtype),target,intent(in) :: geom                              !< Geometry
type(bpartype),target,intent(in) :: bpar                              !< Blocal parameters
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1)                  !< NICAS data, local
type(cvtype),intent(in) :: cv(bpar%nb+1)                              !< Control variable
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field, local

! Local variable
integer :: ib,its,iv
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Apply common localization
   call apply_nicas_sqrt(nam,geom,ndataloc(bpar%nb+1),cv(bpar%nb+1)%alpha,fld_3d)

   ! Apply common ensemble coefficient square-root
   fld_3d = fld_3d*sqrt(ndataloc(bpar%nb+1)%coef_ens)

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
         iv = bpar%ib_to_iv(ib)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt(nam,geom,ndataloc(ib),cv(ib)%alpha,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndataloc(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('specific_multivariate')
   call msgerror('not implemented yet in apply_localization')
case ('common_weighted')
   call msgerror('not implemented yet in apply_localization')
end select

end subroutine apply_localization_sqrt

!----------------------------------------------------------------------
! Subroutine: apply_localization_sqrt_ad
!> Purpose: apply localization square-root, adjoint
!----------------------------------------------------------------------
subroutine apply_localization_sqrt_ad(nam,geom,bpar,ndataloc,fld,cv)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                               !< Namelist
type(geomtype),target,intent(in) :: geom                             !< Geometry
type(bpartype),target,intent(in) :: bpar                             !< Block parameters
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1)                 !< NICAS data, local
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field, local
type(cvtype),intent(out) :: cv(bpar%nb+1)                            !< Control variable

! Local variable
integer :: ib,its,iv
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:)

! Allocation
call cv_alloc(bpar,ndataloc,cv)

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
   fld_3d = fld_3d*sqrt(ndataloc(bpar%nb+1)%coef_ens)

   ! Apply common localization
   call apply_nicas_sqrt_ad(nam,geom,ndataloc(bpar%nb+1),fld_3d,cv(bpar%nb+1)%alpha)
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
         iv = bpar%ib_to_iv(ib)

         ! Apply common ensemble coefficient square-root
         fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(ndataloc(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call apply_nicas_sqrt_ad(nam,geom,ndataloc(ib),fld_4d(:,:,iv),cv(ib)%alpha)
      end if
   end do
case ('specific_multivariate')
   call msgerror('not implemented yet in apply_localization')
case ('common_weighted')
   call msgerror('not implemented yet in apply_localization')
end select

end subroutine apply_localization_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: apply_localization_from_sqrt
!> Purpose: apply localization from square-root
!----------------------------------------------------------------------
subroutine apply_localization_from_sqrt(nam,geom,bpar,ndataloc,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                !< Geometry
type(bpartype),target,intent(in) :: bpar                                !< Blocal parameters
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1)                    !< NICAS data, local
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field, local

! Local variable
type(cvtype) :: cv(bpar%nb+1)

! Apply square-root adjoint
call apply_localization_sqrt_ad(nam,geom,bpar,ndataloc,fld,cv)

! Apply square-root
call apply_localization_sqrt(nam,geom,bpar,ndataloc,cv,fld)

end subroutine apply_localization_from_sqrt

!----------------------------------------------------------------------
! Subroutine: randomize_localization
!> Purpose: randomize localization from square-root
!----------------------------------------------------------------------
subroutine randomize_localization(nam,geom,bpar,ndataloc,ne,ens)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                   !< Namelist
type(geomtype),target,intent(in) :: geom                                 !< Geometry
type(bpartype),target,intent(in) :: bpar                                 !< Blocal parameters
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1)                     !< NICAS data, local
integer,intent(in) :: ne                                                 !< Number of members
real(kind_real),intent(out) :: ens(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne) !< Ensemble, local

! Local variable
integer :: ie
type(cvtype) :: cv(bpar%nb+1,ne)

! Random control vector
call cv_random(bpar,ndataloc,ne,cv)

do ie=1,ne
   ! Apply square-root
   call apply_localization_sqrt(nam,geom,bpar,ndataloc,cv(:,ie),ens(:,:,:,:,ie))
end do

end subroutine randomize_localization

end module module_apply_localization
