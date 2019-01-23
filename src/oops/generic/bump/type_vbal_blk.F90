!----------------------------------------------------------------------
! Module: type_vbal_blk
! Purpose: vertical balance block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_vbal_blk

use tools_kinds, only: kind_real
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

! Vertical balance block derived type
type vbal_blk_type
   integer :: iv                                  ! First variable index
   integer :: jv                                  ! Second variable index
   character(len=1024) :: name                    ! Name
   real(kind_real),allocatable :: auto(:,:,:)     ! Auto-covariance
   real(kind_real),allocatable :: cross(:,:,:)    ! Cross-covariance
   real(kind_real),allocatable :: auto_inv(:,:,:) ! Inverse auto-covariance
   real(kind_real),allocatable :: reg(:,:,:)      ! Regression
contains
   procedure :: alloc => vbal_blk_alloc
   procedure :: partial_dealloc => vbal_blk_partial_dealloc
   procedure :: dealloc => vbal_blk_dealloc
   procedure :: apply => vbal_blk_apply
   procedure :: apply_ad => vbal_blk_apply_ad
end type vbal_blk_type

private
public :: vbal_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: vbal_blk_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine vbal_blk_alloc(vbal_blk,nam,geom,nc2b,iv,jv)

implicit none

! Passed variables
class(vbal_blk_type),intent(inout) :: vbal_blk ! Vertical balance block
type(nam_type),intent(in) :: nam               ! Namelist
type(geom_type),intent(in) :: geom             ! Geometry
integer,intent(in) :: nc2b                     ! Subset Sc2 size, halo B
integer,intent(in) :: iv                       ! First variable index
integer,intent(in) :: jv                       ! Second variable index

! Set attributes
vbal_blk%iv = iv
vbal_blk%jv = jv
vbal_blk%name = trim(nam%varname(jv))//'-'//trim(nam%varname(iv))

! Allocation
allocate(vbal_blk%auto(nc2b,geom%nl0,geom%nl0))
allocate(vbal_blk%cross(nc2b,geom%nl0,geom%nl0))
allocate(vbal_blk%auto_inv(nc2b,geom%nl0,geom%nl0))
allocate(vbal_blk%reg(nc2b,geom%nl0,geom%nl0))

end subroutine vbal_blk_alloc

!----------------------------------------------------------------------
! Subroutine: vbal_blk_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine vbal_blk_partial_dealloc(vbal_blk)

implicit none

! Passed variables
class(vbal_blk_type),intent(inout) :: vbal_blk ! Vertical balance block

! Release memory
if (allocated(vbal_blk%auto)) deallocate(vbal_blk%auto)
if (allocated(vbal_blk%cross)) deallocate(vbal_blk%cross)
if (allocated(vbal_blk%auto_inv)) deallocate(vbal_blk%auto_inv)

end subroutine vbal_blk_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: vbal_blk_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine vbal_blk_dealloc(vbal_blk)

implicit none

! Passed variables
class(vbal_blk_type),intent(inout) :: vbal_blk ! Vertical balance block

! Release memory
call vbal_blk%partial_dealloc
if (allocated(vbal_blk%reg)) deallocate(vbal_blk%reg)

end subroutine vbal_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: vbal_blk_apply
! Purpose: apply vertical balance block
!----------------------------------------------------------------------
subroutine vbal_blk_apply(vbal_blk,geom,np,h_n_s,h_c2b,h_S,fld)

implicit none

! Passed variables
class(vbal_blk_type),intent(in) :: vbal_blk               ! Vertical balance block
type(geom_type),intent(in) :: geom                        ! Geometry
integer,intent(in) :: np                                  ! Maximum number of neighbors
integer,intent(in) :: h_n_s(geom%nc0a,geom%nl0i)          ! Number of neighbors for the horizontal interpolation
integer,intent(in) :: h_c2b(np,geom%nc0a,geom%nl0i)       ! Index of neighbors for the horizontal interpolation
real(kind_real),intent(in) :: h_S(np,geom%nc0a,geom%nl0i) ! Weight of neighbors for the horizontal interpolation
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0)  ! Source/destination vector

! Local variables
integer :: ic0a,il0,jl0,i_s,ic2b
real(kind_real) :: S,fld_tmp(geom%nc0a,geom%nl0)

! Initialization
fld_tmp = 0.0

do ic0a=1,geom%nc0a
   do il0=1,geom%nl0
      do jl0=1,geom%nl0
         ! Loop over neighbors
         do i_s=1,h_n_s(ic0a,min(il0,geom%nl0i))
            ! Find neighbor index and weight
            ic2b = h_c2b(i_s,ic0a,min(il0,geom%nl0i))
            S = h_S(i_s,ic0a,min(il0,geom%nl0i))

            ! Apply regression coefficient weighted by the neighbor weight
            fld_tmp(ic0a,il0) = fld_tmp(ic0a,il0)+S*vbal_blk%reg(ic2b,il0,jl0)*fld(ic0a,jl0)
         end do
      end do
   end do
end do

! Final copy
fld = fld_tmp

end subroutine vbal_blk_apply

!----------------------------------------------------------------------
! Subroutine: vbal_blk_apply_ad
! Purpose: apply adjoint vertical balance block
!----------------------------------------------------------------------
subroutine vbal_blk_apply_ad(vbal_blk,geom,np,h_n_s,h_c2b,h_S,fld)

implicit none

! Passed variables
class(vbal_blk_type),intent(in) :: vbal_blk               ! Vertical balance block
type(geom_type),intent(in) :: geom                        ! Geometry
integer,intent(in) :: np                                  ! Maximum number of neighbors
integer,intent(in) :: h_n_s(geom%nc0a,geom%nl0i)          ! Number of neighbors for the horizontal interpolation
integer,intent(in) :: h_c2b(np,geom%nc0a,geom%nl0i)       ! Index of neighbors for the horizontal interpolation
real(kind_real),intent(in) :: h_S(np,geom%nc0a,geom%nl0i) ! Weight of neighbors for the horizontal interpolation
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0)  ! Source/destination vector

! Local variables
integer :: ic0a,il0,jl0,i_s,ic2b
real(kind_real) :: S,fld_tmp(geom%nc0a,geom%nl0)

! Initialization
fld_tmp = 0.0

do ic0a=1,geom%nc0a
   do il0=1,geom%nl0
      do jl0=1,geom%nl0
         ! Loop over neighbors
         do i_s=1,h_n_s(ic0a,min(il0,geom%nl0i))
            ! Find neighbor index and weight
            ic2b = h_c2b(i_s,ic0a,min(il0,geom%nl0i))
            S = h_S(i_s,ic0a,min(il0,geom%nl0i))

            ! Apply regression coefficient weighted by the neighbor weight
            fld_tmp(ic0a,jl0) = fld_tmp(ic0a,jl0)+S*vbal_blk%reg(ic2b,il0,jl0)*fld(ic0a,il0)
         end do
      end do
   end do
end do

! Final copy
fld = fld_tmp

end subroutine vbal_blk_apply_ad

end module type_vbal_blk
