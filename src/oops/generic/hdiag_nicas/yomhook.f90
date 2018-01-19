!----------------------------------------------------------------------
! Module: yomhook
!> Purpose: dummy module for compatibility with IFS
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module yomhook

use tools_kinds, only: kind_real

implicit none

logical :: lhook = .false. !< Dr Hook activation

private
public :: lhook
public :: dr_hook

contains

!----------------------------------------------------------------------
! Subroutine: dr_hook
!> Purpose: dr_hook
!----------------------------------------------------------------------
subroutine dr_hook(cdname,kswitch,pkey)

implicit none

! Passed variables
character(len=*),intent(in) :: cdname !< Subroutine name
integer,intent(in) :: kswitch         !< Switch
real(kind_real),intent(inout) :: pkey !< Handle

! Local variables
character(len=len(cdname)) :: cdname_loc
integer :: kswitch_loc
real(kind_real) :: pkey_loc

! To avoid warnings
cdname_loc = cdname
kswitch_loc = kswitch
pkey_loc = pkey

end subroutine dr_hook

end module yomhook
