!----------------------------------------------------------------------
! Module: type_min
!> Purpose: minimization data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_min

use tools_kinds, only: kind_real

implicit none

! Minimization data derived type
type mintype
   ! Generic data
   integer :: nx                              !< Control vector size
   integer :: ny                              !< Function output size
   real(kind_real),allocatable :: x(:)        !< Control vector
   real(kind_real),allocatable :: guess(:)    !< Control vector guess
   real(kind_real),allocatable :: norm(:)     !< Control vector norm
   real(kind_real),allocatable :: binf(:)     !< Control vector lower bound
   real(kind_real),allocatable :: bsup(:)     !< Control vector upper bound
   real(kind_real),allocatable :: obs(:)      !< Observation
   real(kind_real) :: f_guess                 !< Guess cost
   character(len=1024) :: fit_type            !< Fit type

   ! Common data
   integer :: nl0                             !< Number of levels
   integer :: nc                              !< Number of classes

   ! Specific data (fit)
   integer :: nl0r                            !< Reduced number of levels
   logical :: lhomh                           !< Vertically homogenous horizontal support radius key
   logical :: lhomv                           !< Vertically homogenous vertical support radius key
   integer,allocatable :: il0rjl0_to_il0(:,:) !< Reduced level to level
   real(kind_real),allocatable :: distvr(:,:) !< Vertical distance
   real(kind_real),allocatable :: disth(:)    !< Horizontal distance

   ! Specific data (LCT)
   integer :: nscales                         !< Number of LCT scales
   integer :: ncomp                           !< Number of LCT components
   real(kind_real),allocatable :: dx(:,:)     !< Zonal separation
   real(kind_real),allocatable :: dy(:,:)     !< Meridian separation
   real(kind_real),allocatable :: dz(:)       !< Vertical separation
   logical,allocatable :: dmask(:,:)          !< Mask
   logical :: lct_diag                        !< Diagonal tensor key
end type mintype

private
public :: mintype

end module type_min
