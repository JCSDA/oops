!----------------------------------------------------------------------
! Module: nicas_apply_nicas.f90
!> Purpose: apply NICAS method
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_apply_nicas

use nicas_apply_convol, only: apply_convol
use nicas_apply_interp, only: apply_interp,apply_interp_ad
use omp_lib
use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_com, only: com_ext,com_red
use type_geom, only: geomtype
use type_mpl, only: mpl,mpl_barrier
use type_nam, only: namtype
use type_ndata, only: ndatatype
use yomhook, only: lhook,dr_hook

implicit none

private
public :: apply_nicas,apply_nicas_sqrt,apply_nicas_sqrt_ad,apply_nicas_from_sqrt

contains

!----------------------------------------------------------------------
! Subroutine: apply_nicas
!> Purpose: apply NICAS method
!----------------------------------------------------------------------
subroutine apply_nicas(geom,ndata,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                        !< Geometry
type(ndatatype),intent(in) :: ndata                      !< NICAS data
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: alpha(:),alpha_tmp(:)

if (lhook) call dr_hook('apply_nicas',0,zhook_handle)

! Allocation
allocate(alpha(ndata%nsb))

! Adjoint interpolation
call apply_interp_ad(geom,ndata,fld,alpha)

! Communication
if (ndata%mpicom==1) then
   ! Allocation
   allocate(alpha_tmp(ndata%nsb))

   ! Copy zone B
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndata%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone B into zone C
   alpha(ndata%sb_to_sc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
elseif (ndata%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call com_red(ndata%AB,alpha)

   ! Allocation
   allocate(alpha_tmp(ndata%nsa))

   ! Copy zone A
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndata%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone A into zone C
   alpha(ndata%sa_to_sc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
end if

! Convolution
call apply_convol(ndata,alpha)

! Halo reduction from zone C to zone A
call com_red(ndata%AC,alpha)

! Halo extension from zone A to zone B
call com_ext(ndata%AB,alpha)

! Interpolation
call apply_interp(geom,ndata,alpha,fld)

! Release memory
deallocate(alpha)

if (lhook) call dr_hook('apply_nicas',1,zhook_handle)

end subroutine apply_nicas

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt
!> Purpose: apply NICAS method square-root
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt(geom,ndata,alpha,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                      !< Geometry
type(ndatatype),intent(in) :: ndata                    !< NICAS data
real(kind_real),intent(in) :: alpha(ndata%nsa)         !< Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variable
integer :: isa
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: alpha_copy(:),alpha_tmp(:)

if (lhook) call dr_hook('apply_nicas_sqrt',0,zhook_handle)

! Allocation
allocate(alpha_copy(ndata%nsa))
allocate(alpha_tmp(ndata%nsa))

! Copy
!$omp parallel do schedule(static) private(isa)
do isa=1,ndata%nsa
   alpha_copy(isa) = alpha(isa)
end do
!$omp end parallel do

! Copy zone A
!$omp parallel do schedule(static) private(isa)
do isa=1,ndata%nsa
   alpha_tmp(isa) = alpha(isa)
end do
!$omp end parallel do!$omp end parallel do

! Reallocation
deallocate(alpha_copy)
allocate(alpha_copy(ndata%nsc))

! Initialize
alpha_copy = 0.0

! Copy zone A into zone C
!$omp parallel do schedule(static) private(isa)
do isa=1,ndata%nsa
   alpha_copy(ndata%sa_to_sc(isa)) = alpha_tmp(isa)
end do
!$omp end parallel do

! Release memory
deallocate(alpha_tmp)

! Convolution
call apply_convol(ndata,alpha_copy)

! Halo reduction from zone C to zone A
call com_red(ndata%AC,alpha_copy)

! Halo extension from zone A to zone B
call com_ext(ndata%AB,alpha_copy)

! Interpolation
call apply_interp(geom,ndata,alpha_copy,fld)

! Release memory
deallocate(alpha_copy)

if (lhook) call dr_hook('apply_nicas_sqrt',1,zhook_handle)

end subroutine apply_nicas_sqrt

!----------------------------------------------------------------------
! Subroutine: apply_nicas_sqrt_ad
!> Purpose: apply NICAS method square-root adjoint
!----------------------------------------------------------------------
subroutine apply_nicas_sqrt_ad(geom,ndata,fld,alpha)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                     !< Geometry
type(ndatatype),intent(in) :: ndata                   !< NICAS data
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field
real(kind_real),intent(out) :: alpha(ndata%nsa)       !< Subgrid field

! Local variable
integer :: isa
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: alpha_copy(:),alpha_tmp(:)

if (lhook) call dr_hook('apply_nicas_sqrt_ad',0,zhook_handle)

! Allocation
allocate(alpha_copy(ndata%nsb))
allocate(alpha_tmp(ndata%nsa))

! Adjoint interpolation
call apply_interp_ad(geom,ndata,fld,alpha_copy)

! Halo reduction from zone B to zone A
call com_red(ndata%AB,alpha_copy)

! Copy zone A
!$omp parallel do schedule(static) private(isa)
do isa=1,ndata%nsa
   alpha_tmp(isa) = alpha_copy(isa)
end do
!$omp end parallel do

! Reallocation
deallocate(alpha_copy)
allocate(alpha_copy(ndata%nsc))

! Initialize
alpha_copy = 0.0

! Copy zone A into zone C
!$omp parallel do schedule(static) private(isa)
do isa=1,ndata%nsa
   alpha_copy(ndata%sa_to_sc(isa)) = alpha_tmp(isa)
end do
!$omp end parallel do

! Release memory
deallocate(alpha_tmp)

! Convolution
call apply_convol(ndata,alpha_copy)

! Halo reduction from zone C to zone A
call com_red(ndata%AC,alpha_copy)

! Copy
!$omp parallel do schedule(static) private(isa)
do isa=1,ndata%nsa
   alpha(isa) = alpha_copy(isa)
end do
!$omp end parallel do

! Release memory
deallocate(alpha_copy)

if (lhook) call dr_hook('apply_nicas_sqrt_ad',1,zhook_handle)

end subroutine apply_nicas_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: apply_nicas_from_sqrt
!> Purpose: apply NICAS method from jts square-root formulation
!----------------------------------------------------------------------
subroutine apply_nicas_from_sqrt(geom,ndata,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                        !< Geometry
type(ndatatype),intent(in) :: ndata                      !< NICAS data
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: zhook_handle
real(kind_real) :: alpha(ndata%nsa)

if (lhook) call dr_hook('apply_nicas_from_sqrt',0,zhook_handle)

! Apply square-root adjoint
call apply_nicas_sqrt_ad(geom,ndata,fld,alpha)

! Apply square-root
call apply_nicas_sqrt(geom,ndata,alpha,fld)

if (lhook) call dr_hook('apply_nicas_from_sqrt',1,zhook_handle)

end subroutine apply_nicas_from_sqrt

end module nicas_apply_nicas
