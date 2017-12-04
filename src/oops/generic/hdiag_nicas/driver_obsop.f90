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

use module_mpi_obsop, only: compute_mpi_obsop
use module_parameters_obsop, only: compute_parameters_obsop
use module_test_obsop, only: test_adjoint_obsop,test_mpi_obsop,test_mpi_obsop_ad,test_obsop
use tools_display, only: msgerror
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_geom, only: geomtype
use type_mpl, only: mpl
use type_nam, only: namtype
use type_odata, only: odatatype,odataloctype

implicit none

private
public :: run_obsop

contains

!----------------------------------------------------------------------
! Subroutine: run_obsop
!> Purpose: observation operator
!----------------------------------------------------------------------
subroutine run_obsop(nam,geom,odata,odataloc)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam       !< Namelist
type(geomtype),target,intent(in) :: geom     !< Geometry
type(odatatype),intent(inout) :: odata       !< Observation operator data
type(odataloctype),intent(inout) :: odataloc !< Observation operator data, local

if (nam%new_obsop) then
   ! Set namelist
   odata%nam => nam
   odataloc%nam => nam

   ! Set geometry
   odata%geom => geom
   odataloc%geom => geom

   ! Compute observation operator parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute observation operator parameters'
   call compute_parameters_obsop(odata)
   call flush(mpl%unit)

   ! Compute observation operator MPI distribution
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute observation operator MPI distribution'
   call compute_mpi_obsop(odata,odataloc)
   call flush(mpl%unit)

   if (nam%check_adjoints) then
      ! Test adjoints
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Test observation operator adjoint'
      call test_adjoint_obsop(odata)
      call flush(mpl%unit)
   end if

   if (nam%check_mpi.and.(mpl%nproc>0)) then
      ! Test single/multi-procs equivalence
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Test observation operator single/multi-procs equivalence'
      call test_mpi_obsop(odata,odataloc)
      call test_mpi_obsop_ad(odata,odataloc)
      call flush(mpl%unit)
   end if

   ! Test precision
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test observation operator precision'   
   call test_obsop(odata,odataloc)
   call flush(mpl%unit)
end if

end subroutine run_obsop

end module driver_obsop
