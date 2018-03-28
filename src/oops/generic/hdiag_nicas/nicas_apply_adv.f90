!----------------------------------------------------------------------
! Module: nicas_apply_adv.f90
!> Purpose: apply advection
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_apply_adv

use omp_lib
use tools_kinds, only: kind_real
use type_com, only: com_ext,com_red
use type_geom, only: geomtype
use type_linop, only: apply_linop,apply_linop_ad
use type_mpl, only: mpl,mpl_barrier
use type_nam, only: namtype
use type_ndata, only: ndatatype
use yomhook, only: lhook,dr_hook

implicit none

private
public :: apply_adv,apply_adv_ad

contains

!----------------------------------------------------------------------
! Subroutine: apply_adv
!> Purpose: apply advection
!----------------------------------------------------------------------
subroutine apply_adv(nam,geom,ndata,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                !< Geometry
type(ndatatype),intent(in) :: ndata                                     !< NICAS data
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real),allocatable :: fld_tmp(:,:)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('apply_adv',0,zhook_handle)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Allocation
      allocate(fld_tmp(geom%nc0a,geom%nl0))

      ! Copy
      fld_tmp = fld(:,:,iv,its)

      ! Halo extension from zone A to zone D
      call com_ext(ndata%AD,fld_tmp)

      ! Interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call apply_linop(ndata%d(il0,its-1),fld_tmp,fld(:,:,iv,its))
      end do
      !$omp end parallel do

      ! Release memory
      deallocate(fld_tmp)
   end do
end do

if (lhook) call dr_hook('apply_adv',1,zhook_handle)

end subroutine apply_adv

!----------------------------------------------------------------------
! Subroutine: apply_adv_ad
!> Purpose: apply advection
!----------------------------------------------------------------------
subroutine apply_adv_ad(nam,geom,ndata,fld)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                !< Geometry
type(ndatatype),intent(in) :: ndata                                     !< NICAS data
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real),allocatable :: fld_tmp(:,:)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('apply_adv_ad',0,zhook_handle)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Allocation
      allocate(fld_tmp(ndata%nc0d,geom%nl0))

      ! Adjoint interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call apply_linop_ad(ndata%d(il0,its-1),fld(:,:,iv,its),fld_tmp)
      end do
      !$omp end parallel do

      ! Halo reduction from zone D to zone A
      call com_red(ndata%AD,fld_tmp)

      ! Copy
      fld(:,:,iv,its) = fld_tmp

      ! Release memory
      deallocate(fld_tmp)
   end do
end do

if (lhook) call dr_hook('apply_adv_ad',1,zhook_handle)

end subroutine apply_adv_ad

end module nicas_apply_adv
