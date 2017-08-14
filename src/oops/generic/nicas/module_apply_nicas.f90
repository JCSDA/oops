!----------------------------------------------------------------------
! Module: module_apply_nicas.f90
!> Purpose: apply NICAS method
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_nicas

use module_apply_com, only: alpha_com_ab,alpha_com_ba,alpha_com_ca,alpha_copy_bc,alpha_copy_ac
use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_apply_nicas_sqrt, only: apply_nicas_sqrt,apply_nicas_sqrt_ad
use module_namelist, only: nam
use type_fields, only: fldtype,alphatype
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype

implicit none

interface apply_nicas
   module procedure apply_nicas_global
   module procedure apply_nicas_local
end interface

interface apply_nicas_from_sqrt
   module procedure apply_nicas_from_sqrt_global
   module procedure apply_nicas_from_sqrt_local
end interface

private
public :: apply_nicas,apply_nicas_from_sqrt

contains

!----------------------------------------------------------------------
! Subroutine: apply_nicas_global
!> Purpose: apply NICAS method, global
!----------------------------------------------------------------------
subroutine apply_nicas_global(ndata,fld)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< Sampling data
type(fldtype),intent(inout) :: fld  !< Field

! Local variables
type(alphatype) :: alpha

! Allocation
allocate(alpha%val(ndata%ns))

! Adjoint interpolation
call interp_ad(ndata,fld,alpha)

! Convolution
call convol(ndata,alpha)

! Interpolation
call interp(ndata,alpha,fld)

! Release memory
deallocate(alpha%val)

end subroutine apply_nicas_global

!----------------------------------------------------------------------
! Subroutine: apply_nicas_local
!> Purpose: apply NICAS method, local
!----------------------------------------------------------------------
subroutine apply_nicas_local(ndataloc,fld)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(fldtype),intent(inout) :: fld     !< Field

! Local variables
type(alphatype) :: alpha

! Allocation
allocate(alpha%valb(ndataloc%nsb))

! Adjoint interpolation
call interp_ad(ndataloc,fld,alpha)

! Communication
if (ndataloc%mpicom==1) then
   ! Copy zone B into zone C
   call alpha_copy_BC(ndataloc,alpha)
elseif (ndataloc%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call alpha_com_BA(ndataloc,alpha)

   ! Copy zone A into zone C
   call alpha_copy_AC(ndataloc,alpha)
end if

! Convolution
call convol(ndataloc,alpha)

! Halo reduction from zone C to zone A
call alpha_com_CA(ndataloc,alpha)

! Halo extension from zone A to zone B
call alpha_com_AB(ndataloc,alpha)

! Interpolation
call interp(ndataloc,alpha,fld)

! Release memory
deallocate(alpha%valb)

end subroutine apply_nicas_local

!----------------------------------------------------------------------
! Subroutine: apply_nicas_from_sqrt_global
!> Purpose: apply NICAS method from its square-root formulation, global
!----------------------------------------------------------------------
subroutine apply_nicas_from_sqrt_global(ndata,fld)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< Sampling data
type(fldtype),intent(inout) :: fld  !< Field

! Local variables
type(alphatype) :: alpha

! Apply square-root adjoint
call apply_nicas_sqrt_ad(ndata,fld,alpha)

! Apply square-root
call apply_nicas_sqrt(ndata,alpha,fld)

end subroutine apply_nicas_from_sqrt_global

!----------------------------------------------------------------------
! Subroutine: apply_nicas_from_sqrt_local
!> Purpose: apply NICAS method from its square-root formulation, local
!----------------------------------------------------------------------
subroutine apply_nicas_from_sqrt_local(ndataloc,fld)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(fldtype),intent(inout) :: fld        !< Field

! Local variables
type(alphatype) :: alpha

! Apply square-root adjoint
call apply_nicas_sqrt_ad(ndataloc,fld,alpha)

! Apply square-root
call apply_nicas_sqrt(ndataloc,alpha,fld)

end subroutine apply_nicas_from_sqrt_local

end module module_apply_nicas
