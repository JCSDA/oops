!----------------------------------------------------------------------
! Module: type_cv
!> Purpose: control vector derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_cv

use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_bpar, only: bpartype
use type_ndata, only: ndataloctype
use type_randgen, only: rand_gau

implicit none

! Control vector derived type
type cvtype
   real(kind_real),allocatable :: alpha(:) !< Control vector field
end type cvtype

private
public :: cvtype
public :: cv_alloc,cv_random

contains

!----------------------------------------------------------------------
! Subroutine: cv_alloc
!> Purpose: control vector object allocation
!----------------------------------------------------------------------
subroutine cv_alloc(bpar,ndataloc,cv)

implicit none

! Passed variables
type(bpartype),target,intent(in) :: bpar             !< Block parameters
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1) !< NICAS data, local
type(cvtype),intent(inout) :: cv(bpar%nb+1)          !< Control vector

! Local variables
integer :: ib

! Allocation
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) allocate(cv(ib)%alpha(ndataloc(ib)%nsa))
end do

end subroutine cv_alloc

!----------------------------------------------------------------------
! Subroutine: cv_random
!> Purpose: generate a random control vector
!----------------------------------------------------------------------
subroutine cv_random(bpar,ndataloc,ne,cv)

implicit none

! Passed variables
type(bpartype),target,intent(in) :: bpar             !< Block parameters
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1) !< NICAS data, local
integer,intent(in) :: ne                             !< Ensemble size
type(cvtype),intent(inout) :: cv(bpar%nb+1,ne)       !< Control vector

! Local variables
integer :: ib,ie
type(cvtype) :: mean(bpar%nb+1)

! Allocation
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      allocate(mean(ib)%alpha(ndataloc(ib)%nsa))
      do ie=1,ne
         allocate(cv(ib,ie)%alpha(ndataloc(ib)%nsa))
      end do
   end if
end do

! Random initialization
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      mean(ib)%alpha = 0.0
      do ie=1,ne
         call rand_gau(cv(ib,ie)%alpha)
         mean(ib)%alpha = mean(ib)%alpha+cv(ib,ie)%alpha
      end do
      mean(ib)%alpha = mean(ib)%alpha/float(ne)
   end if
end do

! Remove mean
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      do ie=1,ne
         cv(ib,ie)%alpha = cv(ib,ie)%alpha-mean(ib)%alpha
      end do
   end if
end do

end subroutine cv_random

end module type_cv
