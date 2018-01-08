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
use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_com, only: com_ext,com_red
use type_geom, only: geomtype
use type_mpl, only: mpl,mpl_barrier
use type_nam, only: namtype
use type_ndata, only: ndatatype

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
integer :: isa,isb
real(kind_real),allocatable :: alpha(:),alpha_tmp(:)

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
   do isb=1,ndata%nsb
      alpha(ndata%sb_to_sc(isb)) = alpha_tmp(isb)
   end do

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
   do isa=1,ndata%nsa
      alpha(ndata%sa_to_sc(isa)) = alpha_tmp(isa)
   end do

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
real(kind_real),allocatable :: alpha_tmp(:),alpha_tmp2(:)

! Allocation
allocate(alpha_tmp(ndata%nsa))
allocate(alpha_tmp2(ndata%nsa))

! Copy
alpha_tmp = alpha

! Copy zone A
alpha_tmp2 = alpha_tmp

! Reallocation
deallocate(alpha_tmp)
allocate(alpha_tmp(ndata%nsc))

! Initialize
alpha_tmp = 0.0

! Copy zone A into zone C
alpha_tmp(ndata%sa_to_sc) = alpha_tmp2

! Release memory
deallocate(alpha_tmp2)

! Convolution
call apply_convol(ndata,alpha_tmp)

! Halo reduction from zone C to zone A
call com_red(ndata%AC,alpha_tmp)

! Halo extension from zone A to zone B
call com_ext(ndata%AB,alpha_tmp)

! Interpolation
call apply_interp(geom,ndata,alpha_tmp,fld)

! Release memory
deallocate(alpha_tmp)

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
real(kind_real),allocatable :: alpha_tmp(:),alpha_tmp2(:)

! Allocation
allocate(alpha_tmp(ndata%nsb))
allocate(alpha_tmp2(ndata%nsa))

! Adjoint interpolation
call apply_interp_ad(geom,ndata,fld,alpha_tmp)

! Halo reduction from zone B to zone A
call com_red(ndata%AB,alpha_tmp)

! Copy zone A
alpha_tmp2 = alpha_tmp

! Reallocation
deallocate(alpha_tmp)
allocate(alpha_tmp(ndata%nsc))

! Initialize
alpha_tmp = 0.0

! Copy zone A into zone C
alpha_tmp(ndata%sa_to_sc) = alpha_tmp2

! Release memory
deallocate(alpha_tmp2)

! Convolution
call apply_convol(ndata,alpha_tmp)

! Halo reduction from zone C to zone A
call com_red(ndata%AC,alpha_tmp)

! Copy
alpha = alpha_tmp

! Release memory
deallocate(alpha_tmp)

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
real(kind_real) :: alpha(ndata%nsa)

! Apply square-root adjoint
call apply_nicas_sqrt_ad(geom,ndata,fld,alpha)

! Apply square-root
call apply_nicas_sqrt(geom,ndata,alpha,fld)

end subroutine apply_nicas_from_sqrt

end module nicas_apply_nicas
