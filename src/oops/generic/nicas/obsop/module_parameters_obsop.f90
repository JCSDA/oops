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
use tools_interp, only: interp_horiz
use tools_kinds, only: kind_real
use type_ctree, only: ctreetype,create_ctree
use type_geom, only: geomtype
use type_linop, only: linop_reorder
use type_mesh, only: create_mesh
use type_mpl, only: mpl
use type_odata, only: odatatype
use type_randgen, only: rng,rand_real

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
type(odatatype),intent(inout) :: odata

! Local variables
integer,allocatable :: mask_ctree(:)
real(kind_real),allocatable :: lonobs(:),latobs(:)
logical,allocatable :: maskobs(:)
type(ctreetype) :: ctree

! Define number of observations
odata%nobs = int(1.0e-2*float(odata%nc0))

! Allocation
allocate(lonobs(odata%nobs))
allocate(latobs(odata%nobs))
allocate(maskobs(odata%nobs))

! Generate random observation network
call rand_real(rng,-180.0_kind_real,180.0_kind_real,.true.,lonobs) 
call rand_real(rng,-90.0_kind_real,90.0_kind_real,.true.,latobs) 
maskobs = .true.

! Create mesh
call create_mesh(rng,odata%nc0,odata%geom%lon,odata%geom%lat,.false.,odata%geom%mesh)

! Compute cover tree
allocate(mask_ctree(odata%geom%mesh%nnr))
mask_ctree = 1
ctree = create_ctree(odata%geom%mesh%nnr,dble(odata%geom%lon(odata%geom%mesh%order)), &
 & dble(odata%geom%lat(odata%geom%mesh%order)),mask_ctree)
deallocate(mask_ctree)

! Compute interpolation
odata%interp%prefix = 'o'
write(mpl%unit,'(a7,a)') '','Single level:'
call interp_horiz(odata%geom%mesh,ctree,odata%nc0,any(odata%geom%mask,dim=2),odata%nobs,lonobs,latobs,maskobs,odata%interp)

! Reorder interpolation
call linop_reorder(odata%interp)

! Print results
write(mpl%unit,'(a7,a,i8)') '','Number of observations: ',odata%nobs

end subroutine compute_parameters_obsop

end module module_parameters_obsop
