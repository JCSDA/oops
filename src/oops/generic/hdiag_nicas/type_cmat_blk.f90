!----------------------------------------------------------------------
! Module: type_cmat_blk
!> Purpose: correlation matrix derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_cmat_blk

use tools_kinds, only: kind_real

implicit none

! B matrix block data derived type
type cmat_blk_type
   ! Block name
   character(len=1024) :: name                     !< Block name

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
end type cmat_blk_type

private
public :: cmat_blk_type

end module type_cmat_blk
