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

use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_apply_nicas_sqrt, only: apply_nicas_sqrt,apply_nicas_sqrt_ad
use module_namelist, only: namtype
use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_com, only: com_ext,com_red
use type_mpl, only: mpl,mpl_barrier
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
subroutine apply_nicas_global(nam,ndata,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(in) :: ndata !< Sampling data
real(kind_real),intent(inout) :: fld(ndata%nc0,ndata%nl0)  !< Field

! Local variables
real(kind_real) :: alpha(ndata%ns)

! Adjoint interpolation
call interp_ad(nam,ndata,fld,alpha)

! Convolution
call convol(ndata,alpha)

! Interpolation
call interp(nam,ndata,alpha,fld)

end subroutine apply_nicas_global

!----------------------------------------------------------------------
! Subroutine: apply_nicas_local
!> Purpose: apply NICAS method, local
!----------------------------------------------------------------------
subroutine apply_nicas_local(nam,ndataloc,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
real(kind_real),intent(inout) :: fld(ndataloc%nc0a,ndataloc%nl0)  !< Field

! Local variables
real(kind_real),allocatable :: alpha(:),alpha_tmp(:)

! Allocation
allocate(alpha(ndataloc%nsb))

! Adjoint interpolation
call interp_ad(nam,ndataloc,fld,alpha)

! Communication
if (nam%mpicom==1) then
   ! Allocation 
   allocate(alpha_tmp(ndataloc%nsb))

   ! Copy zone B
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndataloc%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone B into zone C
   alpha(ndataloc%isb_to_isc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
elseif (nam%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call com_red(ndataloc%AB,alpha)

   ! Allocation 
   allocate(alpha_tmp(ndataloc%nsb))

   ! Copy zone A
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndataloc%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone A into zone C
   alpha(ndataloc%isa_to_isc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
end if

! Convolution
call convol(ndataloc,alpha)

! Halo reduction from zone C to zone A
call com_red(ndataloc%AC,alpha)

! Halo extension from zone A to zone B
call com_ext(ndataloc%AB,alpha)

! Interpolation
call interp(nam,ndataloc,alpha,fld)

! Release memory
deallocate(alpha)

end subroutine apply_nicas_local

!----------------------------------------------------------------------
! Subroutine: apply_nicas_from_sqrt_global
!> Purpose: apply NICAS method from its square-root formulation, global
!----------------------------------------------------------------------
subroutine apply_nicas_from_sqrt_global(nam,ndata,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(in) :: ndata !< Sampling data
real(kind_real),intent(inout) :: fld(ndata%nc0,ndata%nl0)  !< Field

! Local variables
real(kind_real) :: alpha(ndata%ns)

! Apply square-root adjoint
call apply_nicas_sqrt_ad(nam,ndata,fld,alpha)

! Apply square-root
call apply_nicas_sqrt(nam,ndata,alpha,fld)

end subroutine apply_nicas_from_sqrt_global

!----------------------------------------------------------------------
! Subroutine: apply_nicas_from_sqrt_local
!> Purpose: apply NICAS method from its square-root formulation, local
!----------------------------------------------------------------------
subroutine apply_nicas_from_sqrt_local(nam,ndataloc,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
real(kind_real),intent(inout) :: fld(ndataloc%nc0a,ndataloc%nl0)  !< Field

! Local variables
real(kind_real) :: alpha(ndataloc%nsa)

! Apply square-root adjoint
call apply_nicas_sqrt_ad(nam,ndataloc,fld,alpha)

! Apply square-root
call apply_nicas_sqrt(nam,ndataloc,alpha,fld)

end subroutine apply_nicas_from_sqrt_local

end module module_apply_nicas
