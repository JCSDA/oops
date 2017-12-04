!----------------------------------------------------------------------
! Module: module_parameters_obsop.f90
!> Purpose: compute observation operator interpolation
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_parameters_obsop

use tools_display, only: msgerror
use tools_interp, only: compute_interp_bilin
use tools_kinds, only: kind_real
use type_ctree, only: ctreetype,create_ctree
use type_linop, only: linop_reorder
use type_mesh, only: create_mesh
use type_mpl, only: mpl
use type_odata, only: odatatype

implicit none

private
public :: compute_parameters_obsop

contains

!----------------------------------------------------------------------
! Subroutine: compute_parameters_obsop
!> Purpose: compute observation operator interpolation
!----------------------------------------------------------------------
subroutine compute_parameters_obsop(odata)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata !< Observation operator data

! Local variables
logical,allocatable :: mask_ctree(:),maskobs(:)
type(ctreetype) :: ctree

! Associate
associate(geom=>odata%geom)

! Allocation
allocate(maskobs(odata%nobs))

! Generate random observation network
maskobs = .true.

! Create mesh
call create_mesh(geom%nc0,geom%lon,geom%lat,.false.,geom%mesh)

! Compute cover tree
allocate(mask_ctree(geom%mesh%nnr))
mask_ctree = .true.
ctree = create_ctree(geom%mesh%nnr,dble(geom%lon(geom%mesh%order)),dble(geom%lat(geom%mesh%order)),mask_ctree)
deallocate(mask_ctree)

! Compute interpolation
odata%interp%prefix = 'o'
write(mpl%unit,'(a7,a)') '','Single level:'
call compute_interp_bilin(geom%mesh,ctree,geom%nc0,any(geom%mask,dim=2),odata%nobs,odata%lonobs,odata%latobs,maskobs,odata%interp)

! Reorder interpolation
call linop_reorder(odata%interp)

! End associate
end associate

end subroutine compute_parameters_obsop

end module module_parameters_obsop
