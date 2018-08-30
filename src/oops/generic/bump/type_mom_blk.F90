!----------------------------------------------------------------------
! Module: type_mom_blk
!> Purpose: moments block derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_mom_blk

use tools_kinds, only: kind_real

implicit none

! Moments block derived type
type mom_blk_type
   integer :: ne                                  !< Ensemble size
   integer :: nsub                                !< Number of sub-ensembles
   real(kind_real),allocatable :: m1_1(:,:,:,:,:) !< Mean
   real(kind_real),allocatable :: m2_1(:,:,:,:,:) !< Variance
   real(kind_real),allocatable :: m1_2(:,:,:,:,:) !< Mean
   real(kind_real),allocatable :: m2_2(:,:,:,:,:) !< Variance
   real(kind_real),allocatable :: m11(:,:,:,:,:)  !< Covariance
   real(kind_real),allocatable :: m12(:,:,:,:,:)  !< Third-order centered moment
   real(kind_real),allocatable :: m21(:,:,:,:,:)  !< Third-order centered moment
   real(kind_real),allocatable :: m22(:,:,:,:,:)  !< Fourth-order centered moment
   real(kind_real),allocatable :: m1full(:,:,:)   !< Full mean
   real(kind_real),allocatable :: m2full(:,:,:)   !< Full variance
end type mom_blk_type

private
public :: mom_blk_type

end module type_mom_blk
