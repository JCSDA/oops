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
use type_com, only: com_ext,com_red
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
subroutine apply_obsop(odata,fld,obs)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata                            !< Observation operator data
real(kind_real),intent(in) :: fld(odata%nc0a,odata%geom%nl0)   !< Field
real(kind_real),intent(out) :: obs(odata%nobsa,odata%geom%nl0) !< Observations columns

! Local variables
integer :: il0
real(kind_real) :: fld_ext(odata%nc0b,odata%geom%nl0)
real(kind_real),allocatable :: slab(:)

! Associate
associate(geom=>odata%geom)

! Allocation
allocate(slab(odata%nc0a))

! Halo extension
do il0=1,geom%nl0
   slab = fld(:,il0)
   call com_ext(odata%com,slab)
   fld_ext(:,il0) = slab
end do

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop(odata%h,fld_ext(:,il0),obs(:,il0))
end do
!$omp end parallel do

! End associate
end associate

end subroutine apply_obsop

!----------------------------------------------------------------------
! Subroutine: apply_obsop_ad
!> Purpose: observation operator interpolation adjoint
!----------------------------------------------------------------------
subroutine apply_obsop_ad(odata,obs,fld)

implicit none

! Passed variables
type(odatatype),intent(in) :: odata                           !< Observation operator data
real(kind_real),intent(in) :: obs(odata%nobsa,odata%geom%nl0) !< Observations columns
real(kind_real),intent(out) :: fld(odata%nc0a,odata%geom%nl0) !< Field

! Local variables
integer :: il0
real(kind_real) :: fld_ext(odata%nc0b,odata%geom%nl0)
real(kind_real),allocatable :: slab(:)

! Associate
associate(geom=>odata%geom)

! Allocation
allocate(slab(odata%nc0a))

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop_ad(odata%h,obs(:,il0),fld_ext(:,il0))
end do
!$omp end parallel do

! Halo reduction
do il0=1,geom%nl0
   slab = fld_ext(:,il0)
   call com_red(odata%com,slab)
   fld(:,il0) = slab
end do

! End associate
end associate

end subroutine apply_obsop_ad

end module obsop_apply
