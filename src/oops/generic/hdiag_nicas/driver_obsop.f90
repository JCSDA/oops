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

use obsop_parameters, only: compute_parameters
use obsop_test, only: test_adjoint,test_accuracy
use tools_display, only: msgerror
use type_geom, only: geomtype
use type_mpl, only: mpl
use type_nam, only: namtype
use type_odata, only: odatatype

implicit none

private
public :: run_obsop

contains

!----------------------------------------------------------------------
! Subroutine: run_obsop
!> Purpose: observation operator
!----------------------------------------------------------------------
subroutine run_obsop(nam,geom,odata)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam       !< Namelist
type(geomtype),target,intent(in) :: geom     !< Geometry
type(odatatype),intent(inout) :: odata       !< Observation operator data

if (nam%new_obsop) then
   ! Set namelist and geometry
   odata%nam => nam
   odata%geom => geom

   ! Compute observation operator parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute observation operator parameters'
   call compute_parameters(odata)
   call flush(mpl%unit)

   if (nam%check_adjoints) then
      ! Test adjoints
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Test observation operator adjoint'
      call test_adjoint(odata)
      call flush(mpl%unit)
   end if

   ! Test precision
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test observation operator precision'
   call test_accuracy(odata)
   call flush(mpl%unit)
end if

end subroutine run_obsop

end module driver_obsop
