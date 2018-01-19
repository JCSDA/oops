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

use tools_kinds, only: kind_real
use tools_missing, only: msr
use type_linop, only: apply_linop_sym
use type_mpl, only: mpl
use type_ndata, only: ndatatype
use yomhook, only: lhook,dr_hook

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

! Local variables
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('apply_convol',0,zhook_handle)

! Apply linear operator, symmetric
call apply_linop_sym(ndata%c,alpha)

if (lhook) call dr_hook('apply_convol',1,zhook_handle)

end subroutine apply_convol

end module nicas_apply_convol
