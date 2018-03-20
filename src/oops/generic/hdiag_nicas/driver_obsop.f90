!----------------------------------------------------------------------
! Module: driver_obsop
!> Purpose: observation operator driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_obsop

use tools_display, only: msgerror
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type
use type_obsop, only: obsop_type

implicit none

private
public :: run_obsop

contains

!----------------------------------------------------------------------
! Subroutine: run_obsop
!> Purpose: observation operator
!----------------------------------------------------------------------
subroutine run_obsop(nam,geom,obsop)

implicit none

! Passed variables
type(nam_type),target,intent(in) :: nam   !< Namelist
type(geom_type),target,intent(in) :: geom !< Geometry
type(obsop_type),intent(inout) :: obsop   !< Observation operator data

if (nam%new_obsop) then
   ! Compute observation operator parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute observation operator parameters'
   call obsop%parameters(nam,geom)
   call flush(mpl%unit)

   if (nam%check_adjoints) then
      ! Test adjoints
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Test observation operator adjoint'
      call obsop%test_adjoint(geom)
      call flush(mpl%unit)
   end if

   ! Test precision
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test observation operator precision'
   call obsop%test_accuracy(geom)
   call flush(mpl%unit)
end if

end subroutine run_obsop

end module driver_obsop
