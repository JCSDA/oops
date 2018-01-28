!----------------------------------------------------------------------
! Module: obsop_apply.f90
!> Purpose: apply observation operator
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module obsop_apply

use omp_lib
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: isnotmsi
use type_com, only: com_ext,com_red
use type_geom, only: geomtype
use type_linop, only: apply_linop,apply_linop_ad
use type_mpl, only: mpl
use type_odata, only: odatatype

implicit none

private
public :: apply_obsop,apply_obsop_ad

contains

!----------------------------------------------------------------------
! Subroutine: apply_obsop
!> Purpose: observation operator interpolation
!----------------------------------------------------------------------
subroutine apply_obsop(geom,odata,fld,obs)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                        !< Geometry
type(odatatype),intent(in) :: odata                      !< Observation operator data
real(kind_real),intent(in) :: fld(odata%nc0a,geom%nl0)   !< Field
real(kind_real),intent(out) :: obs(odata%nobsa,geom%nl0) !< Observations columns

! Local variables
integer :: il0
real(kind_real) :: fld_ext(odata%nc0b,geom%nl0)
real(kind_real),allocatable :: slab(:)

! Halo extension
do il0=1,geom%nl0
   ! Allocation
   allocate(slab(odata%nc0a))

   ! Copy reduced slab
   slab = fld(:,il0)

   ! Extend halo
   call com_ext(odata%com,slab)

   ! Copy extended slab
   fld_ext(:,il0) = slab

   ! Release memory
   deallocate(slab)
end do

if (odata%nobsa>0) then
   ! Horizontal interpolation
   !$omp parallel do schedule(static) private(il0)
   do il0=1,geom%nl0
      call apply_linop(odata%h,fld_ext(:,il0),obs(:,il0))
   end do
   !$omp end parallel do
end if

end subroutine apply_obsop

!----------------------------------------------------------------------
! Subroutine: apply_obsop_ad
!> Purpose: observation operator interpolation adjoint
!----------------------------------------------------------------------
subroutine apply_obsop_ad(geom,odata,obs,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                       !< Geometry
type(odatatype),intent(in) :: odata                     !< Observation operator data
real(kind_real),intent(in) :: obs(odata%nobsa,geom%nl0) !< Observations columns
real(kind_real),intent(out) :: fld(odata%nc0a,geom%nl0) !< Field

! Local variables
integer :: il0
real(kind_real) :: fld_ext(odata%nc0b,geom%nl0)
real(kind_real),allocatable :: slab(:)

if (odata%nobsa>0) then
   ! Horizontal interpolation
   !$omp parallel do schedule(static) private(il0)
   do il0=1,geom%nl0
      call apply_linop_ad(odata%h,obs(:,il0),fld_ext(:,il0))
   end do
   !$omp end parallel do
end if

! Halo reduction
do il0=1,geom%nl0
   ! Allocation
   allocate(slab(odata%nc0a))

   ! Copy extended slab
   slab = fld_ext(:,il0)

   ! Reduce halo
   call com_red(odata%com,slab)

   ! Copy reduced slab
   fld(:,il0) = slab

   ! Release memory
   deallocate(slab)
end do

end subroutine apply_obsop_ad

end module obsop_apply
