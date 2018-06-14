!----------------------------------------------------------------------
! Module: type_cmat_blk
!> Purpose: correlation matrix derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_cmat_blk

use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_nam, only: nam_type

implicit none

! C matrix block data derived type
type cmat_blk_type
   ! Block index and name
   integer :: ib                                !< Block index
   character(len=1024) :: name                  !< Name

   ! Data
   real(kind_real),allocatable :: coef_ens(:,:)    !< Ensemble coefficient
   real(kind_real),allocatable :: coef_sta(:,:)    !< Static coefficient
   real(kind_real),allocatable :: rh0(:,:)         !< Fit support radius
   real(kind_real),allocatable :: rv0(:,:)         !< Fit support radius
   real(kind_real),allocatable :: rh0s(:,:)        !< Fit support radius  for sampling
   real(kind_real),allocatable :: rv0s(:,:)        !< Fit support radius, for sampling
   real(kind_real) :: wgt                          !< Block weight
   real(kind_real),allocatable :: displ_lon(:,:,:) !< Displaced longitude
   real(kind_real),allocatable :: displ_lat(:,:,:) !< Displaced latitude
contains
   procedure :: alloc => cmat_blk_alloc
   procedure :: dealloc => cmat_blk_dealloc
end type cmat_blk_type

private
public :: cmat_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: cmat_blk_alloc
!> Purpose: cmat block object allocation
!----------------------------------------------------------------------
subroutine cmat_blk_alloc(cmat_blk,nam,geom,bpar,prefix)

implicit none

! Passed variables
class(cmat_blk_type),intent(inout) :: cmat_blk !< C matrix data block
type(nam_type),target,intent(in) :: nam        !< Namelist
type(geom_type),target,intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar             !< Block parameters
character(len=*),intent(in) :: prefix          !< Prefix

! Associate
associate(ib=>cmat_blk%ib)

! Set block name
cmat_blk%name = trim(prefix)//'_'//trim(bpar%blockname(ib))

if (bpar%diag_block(ib)) then
   ! Allocation
   allocate(cmat_blk%coef_ens(geom%nc0a,geom%nl0))
   allocate(cmat_blk%coef_sta(geom%nc0a,geom%nl0))
   allocate(cmat_blk%rh0(geom%nc0a,geom%nl0))
   allocate(cmat_blk%rv0(geom%nc0a,geom%nl0))
   allocate(cmat_blk%rh0s(geom%nc0a,geom%nl0))
   allocate(cmat_blk%rv0s(geom%nc0a,geom%nl0))

   ! Initialization
   call msr(cmat_blk%coef_ens)
   call msr(cmat_blk%coef_sta)
   call msr(cmat_blk%rh0)
   call msr(cmat_blk%rv0)
   call msr(cmat_blk%rh0s)
   call msr(cmat_blk%rv0s)
   call msr(cmat_blk%wgt)
end if

if ((ib==bpar%nb+1).and.nam%displ_diag) then
   ! Allocation
   allocate(cmat_blk%displ_lon(geom%nc0a,geom%nl0,2:nam%nts))
   allocate(cmat_blk%displ_lat(geom%nc0a,geom%nl0,2:nam%nts))

   ! Initialization
   if (nam%displ_diag) then
      call msr(cmat_blk%displ_lon)
      call msr(cmat_blk%displ_lat)
   end if
end if

! End associate
end associate

end subroutine cmat_blk_alloc

!----------------------------------------------------------------------
! Subroutine: cmat_blk_dealloc
!> Purpose: cmat block object deallocation
!----------------------------------------------------------------------
subroutine cmat_blk_dealloc(cmat_blk)

implicit none

! Passed variables
class(cmat_blk_type),intent(inout) :: cmat_blk !< C matrix data block

! Release memory
if (allocated(cmat_blk%coef_ens)) deallocate(cmat_blk%coef_ens)
if (allocated(cmat_blk%coef_sta)) deallocate(cmat_blk%coef_sta)
if (allocated(cmat_blk%rh0)) deallocate(cmat_blk%rh0)
if (allocated(cmat_blk%rv0)) deallocate(cmat_blk%rv0)
if (allocated(cmat_blk%rh0s)) deallocate(cmat_blk%rh0s)
if (allocated(cmat_blk%rv0s)) deallocate(cmat_blk%rv0s)
if (allocated(cmat_blk%displ_lon)) deallocate(cmat_blk%displ_lon)
if (allocated(cmat_blk%displ_lat)) deallocate(cmat_blk%displ_lat)

end subroutine cmat_blk_dealloc

end module type_cmat_blk
