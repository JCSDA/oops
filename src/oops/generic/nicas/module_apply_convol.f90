!----------------------------------------------------------------------
! Module: module_apply_convol.f90
!> Purpose: convolution routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_convol

use omp_lib
use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_linop, only: apply_linop_sym
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype

implicit none

interface convol
   module procedure convol_global
   module procedure convol_local
end interface

private
public :: convol

contains

!----------------------------------------------------------------------
! Subroutine: convol_global
!> Purpose: convolution, global
!----------------------------------------------------------------------
subroutine convol_global(ndata,alpha)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata    !< Sampling data
real(kind_real),intent(inout) :: alpha(ndata%ns) !< Subgrid variable

! Apply linear operator, symmetric
call apply_linop_sym(ndata%c,alpha)

end subroutine convol_global

!----------------------------------------------------------------------
! Subroutine: convol_local
!> Purpose: convolution, local
!----------------------------------------------------------------------
subroutine convol_local(ndataloc,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
real(kind_real),intent(inout) :: alpha(ndataloc%nsc) !< Subgrid variable

! Apply linear operator, symmetric
call apply_linop_sym(ndataloc%c,alpha)

end subroutine convol_local

end module module_apply_convol
