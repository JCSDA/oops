!----------------------------------------------------------------------
! Module: module_apply_obsop.f90
!> Purpose: apply observation operator
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_obsop

use module_namelist, only: namtype
use omp_lib
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_com, only: com_ext,com_red
use type_geom, only: geomtype
use type_linop, only: apply_linop,apply_linop_ad
use type_mpl, only: mpl
use type_odata, only: odatatype,odataloctype

implicit none

interface apply_obsop
   module procedure apply_obsop_global
   module procedure apply_obsop_local
end interface

interface apply_obsop_ad
   module procedure apply_obsop_ad_global
   module procedure apply_obsop_ad_local
end interface

private
public :: apply_obsop,apply_obsop_ad

contains

!----------------------------------------------------------------------
! Subroutine: apply_obsop_global
!> Purpose: observation operator interpolation, global
!----------------------------------------------------------------------
subroutine apply_obsop_global(odata,fld,obs)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata !< Observation operator data
real(kind_real),intent(in) :: fld(odata%nc0,odata%nl0) !< Field
real(kind_real),intent(out) :: obs(odata%nobs,odata%nl0)  !< Observations columns

! Local variables
integer :: il0

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,odata%nl0
   call apply_linop(odata%interp,fld(:,il0),obs(:,il0))
end do
!$omp end parallel do

end subroutine apply_obsop_global

!----------------------------------------------------------------------
! Subroutine: apply_obsop_local
!> Purpose: observation operator interpolation, local
!----------------------------------------------------------------------
subroutine apply_obsop_local(odataloc,fld,obs)

implicit none

! Passed variables
type(odataloctype),intent(in) :: odataloc !< Observation operator data
real(kind_real),intent(in) :: fld(odataloc%nc0a,odataloc%nl0) !< Field
real(kind_real),intent(out) :: obs(odataloc%nobsa,odataloc%nl0)  !< Observations columns

! Local variables
integer :: il0
real(kind_real) :: fld_ext(odataloc%nc0b,odataloc%nl0)
real(kind_real),allocatable :: slab(:)

! Allocation
allocate(slab(odataloc%nc0a))

! Halo extension
do il0=1,odataloc%nl0
   slab = fld(:,il0)
   call com_ext(odataloc%com,slab)
   fld_ext(:,il0) = slab
end do

! Horizontal odataloc%interpolation
!$omp parallel do private(il0)
do il0=1,odataloc%nl0
   call apply_linop(odataloc%interp,fld_ext(:,il0),obs(:,il0))
end do
!$omp end parallel do

end subroutine apply_obsop_local

!----------------------------------------------------------------------
! Subroutine: apply_obsop_ad_global
!> Purpose: observation operator interpolation adjoint, global
!----------------------------------------------------------------------
subroutine apply_obsop_ad_global(odata,obs,fld)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata !< Observation operator data
real(kind_real),intent(in) :: obs(odata%nobs,odata%nl0)  !< Observations columns
real(kind_real),intent(out) :: fld(odata%nc0,odata%nl0) !< Field

! Local variables
integer :: il0

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,odata%nl0
   call apply_linop_ad(odata%interp,obs(:,il0),fld(:,il0))
end do
!$omp end parallel do

end subroutine apply_obsop_ad_global

!----------------------------------------------------------------------
! Subroutine: apply_obsop_ad_local
!> Purpose: observation operator interpolation adjoint, local
!----------------------------------------------------------------------
subroutine apply_obsop_ad_local(odataloc,obs,fld)

implicit none

! Passed variables
type(odataloctype),intent(in) :: odataloc !< Observation operator data
real(kind_real),intent(in) :: obs(odataloc%nobsa,odataloc%nl0)  !< Observations columns
real(kind_real),intent(out) :: fld(odataloc%nc0a,odataloc%nl0) !< Field

! Local variables
integer :: il0
real(kind_real) :: fld_ext(odataloc%nc0b,odataloc%nl0)
real(kind_real),allocatable :: slab(:)

! Allocation
allocate(slab(odataloc%nc0a))

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,odataloc%nl0
   call apply_linop_ad(odataloc%interp,obs(:,il0),fld_ext(:,il0))
end do
!$omp end parallel do

! Halo reduction
do il0=1,odataloc%nl0
   slab = fld_ext(:,il0)
   call com_red(odataloc%com,slab)
   fld(:,il0) = slab
end do

end subroutine apply_obsop_ad_local

end module module_apply_obsop
