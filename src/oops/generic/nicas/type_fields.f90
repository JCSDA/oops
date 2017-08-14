!----------------------------------------------------------------------
! Module: type_fields
!> Purpose: fields derived types
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_fields

use tools_kinds, only: kind_real

implicit none

! Full field derived type
type fldtype
   real(kind_real),allocatable :: val(:,:)  !< Global field data
   real(kind_real),allocatable :: vala(:,:) !< Local field data, halo A
end type fldtype

! Reduced field derived type
type alphatype
   real(kind_real),allocatable :: val(:)  !< Global subgrid variable data
   real(kind_real),allocatable :: vala(:) !< Local subgrid variable data, halo A
   real(kind_real),allocatable :: valb(:) !< Local subgrid variable data, halo B
   real(kind_real),allocatable :: valc(:) !< Local subgrid variable data, halo C
end type alphatype

private
public :: fldtype,alphatype

contains

end module type_fields
