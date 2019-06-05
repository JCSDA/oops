!----------------------------------------------------------------------
! Module: type_mom_blk
! Purpose: moments block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mom_blk

use tools_kinds, only: kind_real

implicit none

! Moments block derived type
type mom_blk_type
   integer :: ne                                  ! Ensemble size
   integer :: nsub                                ! Number of sub-ensembles
   real(kind_real),allocatable :: m2_1(:,:,:)     ! Variance for variable 1
   real(kind_real),allocatable :: m2_2(:,:,:,:)   ! Variance for variable 2
   real(kind_real),allocatable :: m11(:,:,:,:,:)  ! Covariance
   real(kind_real),allocatable :: m22(:,:,:,:,:)  ! Fourth-order centered moment
end type mom_blk_type

private
public :: mom_blk_type

end module type_mom_blk
