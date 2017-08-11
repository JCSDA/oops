!----------------------------------------------------------------------
! Module: module_apply_interp
!> Purpose: interpolation routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_interp

use module_apply_com, only: alpha_com_AB
use module_namelist, only: nam
use omp_lib
use tools_kinds,only: kind_real
use tools_missing, only: msr
use type_fields, only: fldtype,alphatype
use type_linop, only: apply_linop,apply_linop_ad
use type_sdata, only: sdatatype,sdatampitype
implicit none

private
public :: interp,interp_ad

contains

!----------------------------------------------------------------------
! Subroutine: interp
!> Purpose: interpolation
!----------------------------------------------------------------------
subroutine interp(sdata,alpha,fld)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata              !< Sampling data
type(alphatype),intent(in) :: alpha(sdata%nproc) !< Subgrid variable
type(fldtype),intent(inout) :: fld(sdata%nproc)  !< Field

! Local variables
integer :: iproc

if (nam%nproc==0) then
   ! Apply global interpolation
   call interp_global(sdata,alpha(1),fld(1))
elseif (nam%nproc>0) then
   ! Apply local interpolation
   do iproc=1,sdata%nproc
      call interp_local(sdata%mpi(iproc),alpha(iproc),fld(iproc))
   end do
end if

end subroutine interp

!----------------------------------------------------------------------
! Subroutine: interp_global
!> Purpose: interpolation, global
!----------------------------------------------------------------------
subroutine interp_global(sdata,alpha,fld)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata !< Sampling data
type(alphatype),intent(in) :: alpha !< Subgrid variable
type(fldtype),intent(inout) :: fld  !< Field

! Local variables
integer :: is,il1,ic1,il0
real(kind_real) :: beta(sdata%nc1,sdata%nl1),gamma(sdata%nc1,sdata%nl1),delta(sdata%nc1,sdata%nl0)

!$omp parallel do private(is)
do is=1,sdata%ns
   if (nam%lsqrt) then
      ! Internal normalization
      beta(sdata%is_to_ic2(is),sdata%is_to_il1(is)) = alpha%val(is)*sdata%norm_sqrt%val(is)
   else
      ! Copy
      beta(sdata%is_to_ic2(is),sdata%is_to_il1(is)) = alpha%val(is)
   end if
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,sdata%nl1
   call apply_linop(sdata%s(il1),beta(1:sdata%nc2(il1),il1),gamma(:,il1))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1)
do ic1=1,sdata%nc1
   call apply_linop(sdata%v(sdata%vbot(ic1)),gamma(ic1,:),delta(ic1,:))
end do
!$omp end parallel do

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,sdata%nl0
   call apply_linop(sdata%h(min(il0,sdata%nl0i)),delta(:,il0),fld%val(:,il0))
end do
!$omp end parallel do

! Normalization
fld%val = fld%val*sdata%norm%val

end subroutine interp_global

!----------------------------------------------------------------------
! Subroutine: interp_local
!> Purpose: interpolation, local
!----------------------------------------------------------------------
subroutine interp_local(sdatampi,alpha,fld)

implicit none

! Passed variables
type(sdatampitype),intent(in) :: sdatampi !< Sampling data
type(alphatype),intent(in) :: alpha       !< Subgrid variable
type(fldtype),intent(inout) :: fld        !< Field

! Local variables
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(sdatampi%nc1b,sdatampi%nl1),gamma(sdatampi%nc1b,sdatampi%nl1),delta(sdatampi%nc1b,sdatampi%nl0)

!$omp parallel do private(isb)
do isb=1,sdatampi%nsb
   if (nam%lsqrt) then
      ! Internal normalization
      beta(sdatampi%isb_to_ic2b(isb),sdatampi%isb_to_il1(isb)) = alpha%valb(isb)*sdatampi%norm_sqrt%valb(isb)
   else
      ! Copy
      beta(sdatampi%isb_to_ic2b(isb),sdatampi%isb_to_il1(isb)) = alpha%valb(isb)
   end if
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,sdatampi%nl1
   call apply_linop(sdatampi%s(il1),beta(1:sdatampi%nc2b(il1),il1),gamma(:,il1))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1b)
do ic1b=1,sdatampi%nc1b
   call apply_linop(sdatampi%v(sdatampi%vbot(ic1b)),gamma(ic1b,:),delta(ic1b,:))
end do
!$omp end parallel do

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,sdatampi%nl0
   call apply_linop(sdatampi%h(min(il0,sdatampi%nl0i)),delta(:,il0),fld%vala(:,il0))
end do
!$omp end parallel do

! Normalization
fld%vala = fld%vala*sdatampi%norm%vala

end subroutine interp_local

!----------------------------------------------------------------------
! Subroutine: interp_ad
!> Purpose: interpolation adjoint
!----------------------------------------------------------------------
subroutine interp_ad(sdata,fld,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                 !< Sampling data
type(fldtype),intent(in) :: fld(sdata%nproc)        !< Field
type(alphatype),intent(inout) :: alpha(sdata%nproc) !< Subgrid variable

! Local variables
integer :: iproc

if (nam%nproc==0) then
   ! Apply global adjoint interpolation
   call interp_ad_global(sdata,fld(1),alpha(1))
elseif (nam%nproc>0) then
   ! Apply local adjoint interpolation
   do iproc=1,sdata%nproc
      call interp_ad_local(sdata%mpi(iproc),fld(iproc),alpha(iproc))
   end do
end if

end subroutine interp_ad

!----------------------------------------------------------------------
! Subroutine: interp_ad_global
!> Purpose: interpolation adjoint, global
!----------------------------------------------------------------------
subroutine interp_ad_global(sdata,fld,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata    !< Sampling data
type(fldtype),intent(in) :: fld        !< Field
type(alphatype),intent(inout) :: alpha !< Subgrid variable

! Local variables
integer :: is,il1,ic1,il0
real(kind_real) :: beta(sdata%nc1,sdata%nl1),gamma(sdata%nc1,sdata%nl1),delta(sdata%nc1,sdata%nl0)
type(fldtype) :: fld_tmp

! Allocation
allocate(fld_tmp%val(sdata%nc0,sdata%nl0))

! Normalization
fld_tmp%val = fld%val*sdata%norm%val

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,sdata%nl0
   call apply_linop_ad(sdata%h(min(il0,sdata%nl0i)),fld_tmp%val(:,il0),delta(:,il0))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1)
do ic1=1,sdata%nc1
   call apply_linop_ad(sdata%v(sdata%vbot(ic1)),delta(ic1,:),gamma(ic1,:))
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,sdata%nl1
   call apply_linop_ad(sdata%s(il1),gamma(:,il1),beta(1:sdata%nc2(il1),il1))
end do
!$omp end parallel do

!$omp parallel do private(is)
do is=1,sdata%ns
   if (nam%lsqrt) then
      ! Internal normalization
      alpha%val(is) = beta(sdata%is_to_ic2(is),sdata%is_to_il1(is))*sdata%norm_sqrt%val(is)
   else
      ! Copy
      alpha%val(is) = beta(sdata%is_to_ic2(is),sdata%is_to_il1(is))
   end if
end do
!$omp end parallel do

end subroutine interp_ad_global

!----------------------------------------------------------------------
! Subroutine: interp_ad_local
!> Purpose: interpolation adjoint, local
!----------------------------------------------------------------------
subroutine interp_ad_local(sdatampi,fld,alpha)

implicit none

! Passed variables
type(sdatampitype),intent(in) :: sdatampi  !< Sampling data
type(fldtype),intent(in) :: fld            !< Field
type(alphatype),intent(inout) :: alpha     !< Subgrid variable

! Local variables
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(sdatampi%nc1b,sdatampi%nl1),gamma(sdatampi%nc1b,sdatampi%nl1),delta(sdatampi%nc1b,sdatampi%nl0)
type(fldtype) :: fld_tmp

! Allocation
allocate(fld_tmp%vala(sdatampi%nc0a,sdatampi%nl0))

! Normalization
fld_tmp%vala = fld%vala*sdatampi%norm%vala

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,sdatampi%nl0
   call apply_linop_ad(sdatampi%h(min(il0,sdatampi%nl0i)),fld_tmp%vala(:,il0),delta(:,il0))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1b)
do ic1b=1,sdatampi%nc1b
   call apply_linop_ad(sdatampi%v(sdatampi%vbot(ic1b)),delta(ic1b,:),gamma(ic1b,:))
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,sdatampi%nl1
   call apply_linop_ad(sdatampi%s(il1),gamma(:,il1),beta(1:sdatampi%nc2b(il1),il1))
end do
!$omp end parallel do

!$omp parallel do private(isb)
do isb=1,sdatampi%nsb
   if (nam%lsqrt) then
      ! Internal normalization
      alpha%valb(isb) = beta(sdatampi%isb_to_ic2b(isb),sdatampi%isb_to_il1(isb))*sdatampi%norm_sqrt%valb(isb)
   else
      ! Copy
      alpha%valb(isb) = beta(sdatampi%isb_to_ic2b(isb),sdatampi%isb_to_il1(isb))
   end if
end do
!$omp end parallel do

end subroutine interp_ad_local

end module module_apply_interp
