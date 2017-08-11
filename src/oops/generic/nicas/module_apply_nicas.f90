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

use module_apply_com, only: fld_com_gl,fld_com_lg,alpha_copy_AB
use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_namelist, only: nam
use type_fields, only: fldtype,alphatype
use type_mpl, only: mpl
use type_sdata, only: sdatatype

implicit none

private
public :: apply_nicas,apply_nicas_sqrt,apply_nicas_sqrt_ad,apply_nicas_from_sqrt

contains

!----------------------------------------------------------------------
! Subroutine: apply_nicas
!> Purpose: apply NICAS method
!----------------------------------------------------------------------
subroutine apply_nicas(sdata,fld)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata !< Sampling data
type(fldtype),intent(inout) :: fld  !< Field

! Local variables
integer :: iproc
type(fldtype) :: fld_tmp(sdata%nproc)
type(alphatype) :: alpha(sdata%nproc)

! Allocation
if (nam%nproc==0) then
   allocate(fld_tmp(1)%val(sdata%nc0,sdata%nl0))
   allocate(alpha(1)%val(sdata%ns))
elseif (nam%nproc>0) then
   do iproc=1,sdata%nproc
      allocate(alpha(iproc)%valb(sdata%mpi(iproc)%nsb))
   end do
end if

if (nam%nproc==0) then
   ! Copy data
   fld_tmp(1)%val = fld%val
elseif (nam%nproc>0) then
   ! Global to local
   call fld_com_gl(sdata,fld,fld_tmp)
end if

! Adjoint interpolation
call interp_ad(sdata,fld_tmp,alpha)

! Convolution
call convol(sdata,alpha)

! Interpolation
call interp(sdata,alpha,fld_tmp)

if (nam%nproc==0) then
   ! Copy data
   fld%val = fld_tmp(1)%val
elseif (nam%nproc>0) then
   ! Local to global
   call fld_com_lg(sdata,fld_tmp,fld)
end if

! Release memory
if (nam%nproc==0) then
   deallocate(fld_tmp(1)%val)
   deallocate(alpha(1)%val)
elseif (nam%nproc>0) then
   do iproc=1,sdata%nproc
      deallocate(alpha(iproc)%valb)
   end do
end if
end subroutine apply_nicas

!----------------------------------------------------------------------
! Subroutine: apply_nicas_from_sqrt
!> Purpose: apply NICAS method from its square-root formulation
!----------------------------------------------------------------------
subroutine apply_nicas_from_sqrt(sdata,fld)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata !< Sampling data
type(fldtype),intent(inout) :: fld  !< Field

! Local variables
type(alphatype) :: alpha(sdata%nproc)

! Apply square-root adjoint
call apply_nicas_sqrt_ad(sdata,fld,alpha)

! Apply square-root
call apply_nicas_sqrt(sdata,alpha,fld)

end subroutine apply_nicas_from_sqrt

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt
!> Purpose: apply NICAS method square-root
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt(sdata,alpha,fld)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata              !< Sampling data
type(alphatype),intent(in) :: alpha(sdata%nproc) !< Subgrid variable
type(fldtype),intent(inout) :: fld               !< Field

! Local variables
integer :: iproc
type(fldtype) :: fld_tmp(sdata%nproc)
type(alphatype) :: alpha_tmp(sdata%nproc)

! Allocation
if (nam%nproc==0) then
   allocate(fld_tmp(1)%val(sdata%nc0,sdata%nl0))
   allocate(alpha_tmp(1)%val(sdata%ns))
elseif (nam%nproc>0) then
   do iproc=1,sdata%nproc
      allocate(fld_tmp(iproc)%vala(sdata%mpi(iproc)%nc0a,sdata%mpi(iproc)%nl0))
      allocate(alpha_tmp(iproc)%vala(sdata%mpi(iproc)%nsa))
   end do
end if
if (.not.allocated(fld%val)) allocate(fld%val(sdata%nc0,sdata%nl0))

if (nam%nproc==0) then
   ! Copy data
   alpha_tmp(1)%val = alpha(1)%val
elseif (nam%nproc>0) then
   ! Copy from zone A to zone B
   do iproc=1,sdata%nproc
      alpha_tmp(iproc)%vala = alpha(iproc)%vala
      call alpha_copy_AB(sdata%mpi(iproc),alpha_tmp(iproc))
   end do
end if

! Convolution
call convol(sdata,alpha_tmp)

! Interpolation
call interp(sdata,alpha_tmp,fld_tmp)

if (nam%nproc==0) then
   ! Copy data
   fld%val = fld_tmp(1)%val
elseif (nam%nproc>0) then
   ! Local to global
   call fld_com_lg(sdata,fld_tmp,fld)
end if

! Release memory
if (nam%nproc==0) then
   deallocate(fld_tmp(1)%val)
   deallocate(alpha_tmp(1)%val)
elseif (nam%nproc>0) then
   do iproc=1,sdata%nproc
      deallocate(fld_tmp(iproc)%vala)
      deallocate(alpha_tmp(iproc)%vala)
   end do
end if

end subroutine apply_nicas_sqrt

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_ad
!> Purpose: apply NICAS method square-root adjoint
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_ad(sdata,fld,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                 !< Sampling data
type(fldtype),intent(in) :: fld                     !< Field
type(alphatype),intent(inout) :: alpha(sdata%nproc) !< Subgrid variable

! Local variables
integer :: iproc
type(fldtype) :: fld_copy,fld_tmp(sdata%nproc)

! Allocation
if (nam%nproc==0) then
   allocate(fld_tmp(1)%val(sdata%nc0,sdata%nl0))
   if (.not.allocated(alpha(1)%val)) allocate(alpha(1)%val(sdata%ns))
elseif (nam%nproc>0) then
   allocate(fld_copy%val(sdata%nc0,sdata%nl0))
   do iproc=1,sdata%nproc
      if (.not.allocated(alpha(iproc)%valb)) allocate(alpha(iproc)%valb(sdata%mpi(iproc)%nsb))
   end do
end if

if (nam%nproc==0) then
   ! Copy data
   fld_tmp(1)%val = fld%val
elseif (nam%nproc>0) then
   ! Global to local
   fld_copy%val = fld%val
   call fld_com_gl(sdata,fld_copy,fld_tmp)
end if

! Adjoint interpolation
call interp_ad(sdata,fld_tmp,alpha)

! Convolution
call convol(sdata,alpha)

! Release memory
if (nam%nproc==0) then
   deallocate(fld_tmp(1)%val)
elseif (nam%nproc>0) then
   deallocate(fld_copy%val)
   do iproc=1,sdata%nproc
      deallocate(fld_tmp(iproc)%vala)
   end do
end if

end subroutine apply_nicas_sqrt_ad

end module module_apply_nicas
