!----------------------------------------------------------------------
! Module: nicas_apply_convol.f90
!> Purpose: convolution routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_apply_convol

use omp_lib
use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_linop, only: apply_linop_sym
use type_mpl, only: mpl
use type_ndata, only: ndatatype

implicit none

private
public :: apply_convol

contains

!----------------------------------------------------------------------
! Subroutine: apply_convol
!> Purpose: apply convolution
!----------------------------------------------------------------------
subroutine apply_convol(ndata,alpha)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata               !< NICAS data
real(kind_real),intent(inout) :: alpha(ndata%nsc) !< Subgrid field

! Apply linear operator, symmetric
call apply_linop_sym(ndata%c,alpha)

end subroutine apply_convol

end module nicas_apply_convol
