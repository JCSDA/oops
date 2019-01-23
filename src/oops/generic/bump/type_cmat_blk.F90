!----------------------------------------------------------------------
! Module: type_cmat_blk
! Purpose: correlation matrix derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_cmat_blk

use tools_kinds, only: kind_real
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

! C matrix block data derived type
type cmat_blk_type
   ! Block index and name
   integer :: ib                                     ! Block index
   character(len=1024) :: name                       ! Name
   logical :: double_fit                             ! Double fit
   logical :: anisotropic                            ! Anisoptropic tensor

   ! Read data
   real(kind_real),allocatable :: oops_coef_ens(:,:) ! OOPS ensemble coefficient
   real(kind_real),allocatable :: oops_coef_sta(:,:) ! OOPS static coefficient
   real(kind_real),allocatable :: oops_rh(:,:)       ! OOPS horizontal fit support radius
   real(kind_real),allocatable :: oops_rv(:,:)       ! OOPS vertical fit support radius
   real(kind_real),allocatable :: oops_rv_rfac(:,:)  ! OOPS vertical fit support radius factor
   real(kind_real),allocatable :: oops_rv_coef(:,:)  ! OOPS vertical fit coefficient

   ! Data
   real(kind_real),allocatable :: coef_ens(:,:)      ! Ensemble coefficient
   real(kind_real),allocatable :: coef_sta(:,:)      ! Static coefficient
   real(kind_real),allocatable :: rh(:,:)            ! Horizontal fit support radius
   real(kind_real),allocatable :: rv(:,:)            ! Vertical fit support radius
   real(kind_real),allocatable :: rv_rfac(:,:)       ! Vertical fit support radius factor
   real(kind_real),allocatable :: rv_coef(:,:)       ! Vertical fit coefficient
   real(kind_real),allocatable :: rhs(:,:)           ! Fit support radius  for sampling
   real(kind_real),allocatable :: rvs(:,:)           ! Fit support radius, for sampling
   real(kind_real) :: wgt                            ! Block weight
   real(kind_real),allocatable :: H11(:,:)           ! LCT component 11
   real(kind_real),allocatable :: H22(:,:)           ! LCT component 22
   real(kind_real),allocatable :: H33(:,:)           ! LCT component 33
   real(kind_real),allocatable :: H12(:,:)           ! LCT component 12
   real(kind_real),allocatable :: Hcoef(:,:)         ! LCT scales coefficients
   real(kind_real),allocatable :: displ_lon(:,:,:)   ! Displaced longitude
   real(kind_real),allocatable :: displ_lat(:,:,:)   ! Displaced latitude
contains
   procedure :: alloc => cmat_blk_alloc
   procedure :: dealloc => cmat_blk_dealloc
   procedure :: copy => cmat_blk_copy
end type cmat_blk_type

private
public :: cmat_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: cmat_blk_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine cmat_blk_alloc(cmat_blk,nam,geom,bpar)

implicit none

! Passed variables
class(cmat_blk_type),intent(inout) :: cmat_blk ! C matrix data block
type(nam_type),intent(in) :: nam               ! Namelist
type(geom_type),intent(in) :: geom             ! Geometry
type(bpar_type),intent(in) :: bpar             ! Block parameters

! Associate
associate(ib=>cmat_blk%ib)

if (bpar%diag_block(ib)) then
   ! Allocation
   allocate(cmat_blk%coef_ens(geom%nc0a,geom%nl0))
   allocate(cmat_blk%coef_sta(geom%nc0a,geom%nl0))
   allocate(cmat_blk%rh(geom%nc0a,geom%nl0))
   allocate(cmat_blk%rv(geom%nc0a,geom%nl0))
   if (cmat_blk%double_fit) then
      allocate(cmat_blk%rv_rfac(geom%nc0a,geom%nl0))
      allocate(cmat_blk%rv_coef(geom%nc0a,geom%nl0))
   end if
   allocate(cmat_blk%rhs(geom%nc0a,geom%nl0))
   allocate(cmat_blk%rvs(geom%nc0a,geom%nl0))
   if (cmat_blk%anisotropic) then
      allocate(cmat_blk%H11(geom%nc0a,geom%nl0))
      allocate(cmat_blk%H22(geom%nc0a,geom%nl0))
      allocate(cmat_blk%H33(geom%nc0a,geom%nl0))
      allocate(cmat_blk%H12(geom%nc0a,geom%nl0))
      allocate(cmat_blk%Hcoef(geom%nc0a,geom%nl0))
   end if
end if

if ((ib==bpar%nbe).and.nam%displ_diag) then
   ! Allocation
   allocate(cmat_blk%displ_lon(geom%nc0a,geom%nl0,2:nam%nts))
   allocate(cmat_blk%displ_lat(geom%nc0a,geom%nl0,2:nam%nts))
end if

! End associate
end associate

end subroutine cmat_blk_alloc

!----------------------------------------------------------------------
! Subroutine: cmat_blk_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine cmat_blk_dealloc(cmat_blk)

implicit none

! Passed variables
class(cmat_blk_type),intent(inout) :: cmat_blk ! C matrix data block

! Release memory
if (allocated(cmat_blk%oops_coef_ens)) deallocate(cmat_blk%oops_coef_ens)
if (allocated(cmat_blk%oops_coef_sta)) deallocate(cmat_blk%oops_coef_sta)
if (allocated(cmat_blk%oops_rh)) deallocate(cmat_blk%oops_rh)
if (allocated(cmat_blk%oops_rv)) deallocate(cmat_blk%oops_rv)
if (allocated(cmat_blk%oops_rv_rfac)) deallocate(cmat_blk%oops_rv_rfac)
if (allocated(cmat_blk%oops_rv_coef)) deallocate(cmat_blk%oops_rv_coef)
if (allocated(cmat_blk%coef_ens)) deallocate(cmat_blk%coef_ens)
if (allocated(cmat_blk%coef_sta)) deallocate(cmat_blk%coef_sta)
if (allocated(cmat_blk%rh)) deallocate(cmat_blk%rh)
if (allocated(cmat_blk%rv)) deallocate(cmat_blk%rv)
if (allocated(cmat_blk%rv_rfac)) deallocate(cmat_blk%rv_rfac)
if (allocated(cmat_blk%rv_coef)) deallocate(cmat_blk%rv_coef)
if (allocated(cmat_blk%rhs)) deallocate(cmat_blk%rhs)
if (allocated(cmat_blk%rvs)) deallocate(cmat_blk%rvs)
if (allocated(cmat_blk%H11)) deallocate(cmat_blk%H11)
if (allocated(cmat_blk%H22)) deallocate(cmat_blk%H22)
if (allocated(cmat_blk%H33)) deallocate(cmat_blk%H33)
if (allocated(cmat_blk%H12)) deallocate(cmat_blk%H12)
if (allocated(cmat_blk%Hcoef)) deallocate(cmat_blk%Hcoef)
if (allocated(cmat_blk%displ_lon)) deallocate(cmat_blk%displ_lon)
if (allocated(cmat_blk%displ_lat)) deallocate(cmat_blk%displ_lat)

end subroutine cmat_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: cmat_blk_copy
! Purpose: copy
!----------------------------------------------------------------------
type(cmat_blk_type) function cmat_blk_copy(cmat_blk)

implicit none

! Passed variables
class(cmat_blk_type),intent(in) :: cmat_blk ! C matrix data block

! Copy data
cmat_blk_copy%ib = cmat_blk%ib
cmat_blk_copy%name = cmat_blk%name
cmat_blk_copy%double_fit = cmat_blk%double_fit
cmat_blk_copy%anisotropic = cmat_blk%anisotropic
if (allocated(cmat_blk%oops_coef_ens)) cmat_blk_copy%oops_coef_ens = cmat_blk%oops_coef_ens
if (allocated(cmat_blk%oops_coef_sta)) cmat_blk_copy%oops_coef_sta = cmat_blk%oops_coef_sta
if (allocated(cmat_blk%oops_rh)) cmat_blk_copy%oops_rh = cmat_blk%oops_rh
if (allocated(cmat_blk%oops_rv)) cmat_blk_copy%oops_rv = cmat_blk%oops_rv
if (allocated(cmat_blk%oops_rv_rfac)) cmat_blk_copy%oops_rv_rfac = cmat_blk%oops_rv_rfac
if (allocated(cmat_blk%oops_rv_coef)) cmat_blk_copy%oops_rv_coef = cmat_blk%oops_rv_coef
if (allocated(cmat_blk%coef_ens)) cmat_blk_copy%coef_ens = cmat_blk%coef_ens
if (allocated(cmat_blk%coef_sta)) cmat_blk_copy%coef_sta = cmat_blk%coef_sta
if (allocated(cmat_blk%rh)) cmat_blk_copy%rh = cmat_blk%rh
if (allocated(cmat_blk%rv)) cmat_blk_copy%rv = cmat_blk%rv
if (allocated(cmat_blk%rv_rfac)) cmat_blk_copy%rv_rfac = cmat_blk%rv_rfac
if (allocated(cmat_blk%rv_coef)) cmat_blk_copy%rv_coef = cmat_blk%rv_coef
if (allocated(cmat_blk%rhs)) cmat_blk_copy%rhs = cmat_blk%rhs
if (allocated(cmat_blk%rvs)) cmat_blk_copy%rvs = cmat_blk%rvs
if (allocated(cmat_blk%H11)) cmat_blk_copy%H11 = cmat_blk%H11
if (allocated(cmat_blk%H22)) cmat_blk_copy%H22 = cmat_blk%H22
if (allocated(cmat_blk%H33)) cmat_blk_copy%H33 = cmat_blk%H33
if (allocated(cmat_blk%H12)) cmat_blk_copy%H12 = cmat_blk%H12
if (allocated(cmat_blk%Hcoef)) cmat_blk_copy%Hcoef = cmat_blk%Hcoef
if (allocated(cmat_blk%displ_lon)) cmat_blk_copy%displ_lon = cmat_blk%displ_lon
if (allocated(cmat_blk%displ_lat)) cmat_blk_copy%displ_lat = cmat_blk%displ_lat

end function cmat_blk_copy

end module type_cmat_blk
