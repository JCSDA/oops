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

use module_apply_com, only: alpha_com_AB,alpha_com_BA,alpha_com_AC,alpha_com_CA
use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_namelist, only: nam
use tools_kinds, only: kind_real
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
real(kind_real),intent(in) :: alpha(ndata%ns) !< Subgrid variable
real(kind_real),intent(out) :: fld(ndata%nc0,ndata%nl0)  !< Field

! Local variable
real(kind_real) :: alpha_tmp(ndata%ns)

! Copy
alpha_tmp = alpha

! Convolution
call convol(ndata,alpha_tmp)

! Interpolation
call interp(ndata,alpha_tmp,fld)

end subroutine apply_nicas_sqrt_global

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_local
!> Purpose: apply NICAS method square-root, local
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_local(ndataloc,alpha,fld)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
real(kind_real),intent(in) :: alpha(ndataloc%nsa) !< Subgrid variable
real(kind_real),intent(out) :: fld(ndataloc%nc0a,ndataloc%nl0)  !< Field

! Local variable
real(kind_real),allocatable :: alpha_tmp(:)

! Allocation
allocate(alpha_tmp(ndataloc%nsa))

! Copy
alpha_tmp = alpha

! Halo extension from zone A to zone C
call alpha_com_AC(ndataloc,alpha_tmp)

! Convolution
call convol(ndataloc,alpha_tmp)

! Halo reduction from zone C to zone A
call alpha_com_CA(ndataloc,alpha_tmp)

! Halo extension from zone A to zone B
call alpha_com_AB(ndataloc,alpha_tmp)

! Interpolation
call interp(ndataloc,alpha_tmp,fld)

! Release memory
deallocate(alpha_tmp)

end subroutine apply_nicas_sqrt_local

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_ad, global
!> Purpose: apply NICAS method square-root adjoint, global
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_ad_global(ndata,fld,alpha)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata    !< Sampling data
real(kind_real),intent(in) :: fld(ndata%nc0,ndata%nl0)  !< Field
real(kind_real),intent(out) :: alpha(ndata%ns) !< Subgrid variable

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
real(kind_real),intent(in) :: fld(ndataloc%nc0a,ndataloc%nl0)  !< Field
real(kind_real),intent(out) :: alpha(ndataloc%nsa) !< Subgrid variable

! Local variable
real(kind_real),allocatable :: alpha_tmp(:)

! Allocation
allocate(alpha_tmp(ndataloc%nsb))

! Adjoint interpolation
call interp_ad(ndataloc,fld,alpha_tmp)

! Halo reduction from zone B to zone A
call alpha_com_BA(ndataloc,alpha_tmp)

! Halo extension from zone A to zone C
call alpha_com_AC(ndataloc,alpha_tmp)

! Convolution
call convol(ndataloc,alpha_tmp)

! Halo reduction from zone C to zone A
call alpha_com_CA(ndataloc,alpha_tmp)

! Copy
alpha = alpha_tmp

! Release memory
deallocate(alpha_tmp)

end subroutine apply_nicas_sqrt_ad_local

end module module_apply_nicas_sqrt
