!----------------------------------------------------------------------
! Module: module_apply_nicas_sqrt.f90
!> Purpose: apply NICAS method, square-root method
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_nicas_sqrt

use module_apply_com, only: alpha_com_ab,alpha_com_ba
use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_namelist, only: nam
use type_fields, only: fldtype,alphatype
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype

implicit none

interface apply_nicas_sqrt
   module procedure apply_nicas_sqrt_global
   module procedure apply_nicas_sqrt_local
end interface

interface apply_nicas_sqrt_ad
   module procedure apply_nicas_sqrt_ad_global
   module procedure apply_nicas_sqrt_ad_local
end interface

private
public :: apply_nicas_sqrt,apply_nicas_sqrt_ad

contains

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_global
!> Purpose: apply NICAS method square-root, global
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_global(ndata,alpha,fld)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< Sampling data
type(alphatype),intent(inout) :: alpha !< Subgrid variable
type(fldtype),intent(inout) :: fld  !< Field

! Convolution
call convol(ndata,alpha)

! Interpolation
call interp(ndata,alpha,fld)

end subroutine apply_nicas_sqrt_global

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_local
!> Purpose: apply NICAS method square-root, local
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_local(ndataloc,alpha,fld)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(alphatype),intent(inout) :: alpha       !< Subgrid variable
type(fldtype),intent(inout) :: fld        !< Field

! Halo reduction from zone A to zone B
call alpha_com_AB(ndataloc,alpha)

! Convolution
call convol(ndataloc,alpha)

! Interpolation
call interp(ndataloc,alpha,fld)

end subroutine apply_nicas_sqrt_local

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_ad, global
!> Purpose: apply NICAS method square-root adjoint, global
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_ad_global(ndata,fld,alpha)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata    !< Sampling data
type(fldtype),intent(in) :: fld        !< Field
type(alphatype),intent(inout) :: alpha !< Subgrid variable

! Adjoint interpolation
call interp_ad(ndata,fld,alpha)

! Convolution
call convol(ndata,alpha)

end subroutine apply_nicas_sqrt_ad_global

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_ad_local
!> Purpose: apply NICAS method square-root adjoint, local
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_ad_local(ndataloc,fld,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(fldtype),intent(in) :: fld           !< Field
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Adjoint interpolation
call interp_ad(ndataloc,fld,alpha)

! Convolution
call convol(ndataloc,alpha)

! Halo reduction from zone B to zone A
call alpha_com_BA(ndataloc,alpha)

end subroutine apply_nicas_sqrt_ad_local

end module module_apply_nicas_sqrt
