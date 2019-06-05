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
   real(kind_real),allocatable :: bump_coef_ens(:,:) ! BUMP ensemble coefficient
   real(kind_real),allocatable :: bump_coef_sta(:,:) ! BUMP static coefficient
   real(kind_real),allocatable :: bump_rh(:,:)       ! BUMP horizontal fit support radius
   real(kind_real),allocatable :: bump_rv(:,:)       ! BUMP vertical fit support radius
   real(kind_real),allocatable :: bump_rv_rfac(:,:)  ! BUMP vertical fit support radius factor
   real(kind_real),allocatable :: bump_rv_coef(:,:)  ! BUMP vertical fit coefficient
   real(kind_real),allocatable :: bump_D11(:,:)      ! BUMP Daley tensor component 11
   real(kind_real),allocatable :: bump_D22(:,:)      ! BUMP Daley tensor component 22
   real(kind_real),allocatable :: bump_D33(:,:)      ! BUMP Daley tensor component 33
   real(kind_real),allocatable :: bump_D12(:,:)      ! BUMP Daley tensor component 12
   real(kind_real),allocatable :: bump_Dcoef(:,:)    ! BUMP Daley tensor scales coefficients

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
   real(kind_real),allocatable :: adv_lon(:,:,:)     ! Advected longitude
   real(kind_real),allocatable :: adv_lat(:,:,:)     ! Advected latitude
contains
   procedure :: alloc => cmat_blk_alloc
   procedure :: init => cmat_blk_init
   procedure :: dealloc => cmat_blk_dealloc
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

if ((ib==bpar%nbe).and.nam%adv_diag) then
   ! Allocation
   allocate(cmat_blk%adv_lon(geom%nc0a,geom%nl0,nam%nts))
   allocate(cmat_blk%adv_lat(geom%nc0a,geom%nl0,nam%nts))
end if

! End associate
end associate

end subroutine cmat_blk_alloc

!----------------------------------------------------------------------
! Subroutine: cmat_blk_init
! Purpose: initialization
!----------------------------------------------------------------------
subroutine cmat_blk_init(cmat_blk,mpl,nam,bpar)

implicit none

! Passed variables
class(cmat_blk_type),intent(inout) :: cmat_blk ! C matrix data block
type(mpl_type),intent(in) :: mpl               ! MPI data
type(nam_type),intent(in) :: nam               ! Namelist
type(bpar_type),intent(in) :: bpar             ! Block parameters

! Associate
associate(ib=>cmat_blk%ib)

if (bpar%diag_block(ib)) then
   ! Initialization
   cmat_blk%coef_ens = mpl%msv%valr
   cmat_blk%coef_sta = mpl%msv%valr
   cmat_blk%rh = mpl%msv%valr
   cmat_blk%rv = mpl%msv%valr
   if (cmat_blk%double_fit) then
      cmat_blk%rv_rfac = mpl%msv%valr
      cmat_blk%rv_coef = mpl%msv%valr
   end if
   cmat_blk%rhs = mpl%msv%valr
   cmat_blk%rvs = mpl%msv%valr
   if (cmat_blk%anisotropic) then
      cmat_blk%H11 = mpl%msv%valr
      cmat_blk%H22 = mpl%msv%valr
      cmat_blk%H33 = mpl%msv%valr
      cmat_blk%H12 = mpl%msv%valr
      cmat_blk%Hcoef = mpl%msv%valr
   end if
end if

if ((ib==bpar%nbe).and.nam%adv_diag) then
   ! Initialization
   cmat_blk%adv_lon = mpl%msv%valr
   cmat_blk%adv_lat = mpl%msv%valr
end if

! End associate
end associate

end subroutine cmat_blk_init

!----------------------------------------------------------------------
! Subroutine: cmat_blk_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine cmat_blk_dealloc(cmat_blk)

implicit none

! Passed variables
class(cmat_blk_type),intent(inout) :: cmat_blk ! C matrix data block

! Release memory
if (allocated(cmat_blk%bump_coef_ens)) deallocate(cmat_blk%bump_coef_ens)
if (allocated(cmat_blk%bump_coef_sta)) deallocate(cmat_blk%bump_coef_sta)
if (allocated(cmat_blk%bump_rh)) deallocate(cmat_blk%bump_rh)
if (allocated(cmat_blk%bump_rv)) deallocate(cmat_blk%bump_rv)
if (allocated(cmat_blk%bump_rv_rfac)) deallocate(cmat_blk%bump_rv_rfac)
if (allocated(cmat_blk%bump_rv_coef)) deallocate(cmat_blk%bump_rv_coef)
if (allocated(cmat_blk%bump_D11)) deallocate(cmat_blk%bump_D11)
if (allocated(cmat_blk%bump_D22)) deallocate(cmat_blk%bump_D22)
if (allocated(cmat_blk%bump_D33)) deallocate(cmat_blk%bump_D33)
if (allocated(cmat_blk%bump_D12)) deallocate(cmat_blk%bump_D12)
if (allocated(cmat_blk%bump_Dcoef)) deallocate(cmat_blk%bump_Dcoef)
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
if (allocated(cmat_blk%adv_lon)) deallocate(cmat_blk%adv_lon)
if (allocated(cmat_blk%adv_lat)) deallocate(cmat_blk%adv_lat)

end subroutine cmat_blk_dealloc

end module type_cmat_blk
