!----------------------------------------------------------------------
! Module: module_nicas
!> Purpose: nicas driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_nicas

use module_mpi, only: compute_mpi
use module_namelist, only: nam
use module_normalization, only: compute_normalization
use module_parameters, only: compute_parameters
use module_test, only: test_adjoints,test_pos_def,test_mpi,test_dirac,test_perf
use tools_const, only: pi
use type_mpl, only: mpl
use type_sdata, only: sdatatype,sdata_read_param,sdata_read_local,sdata_read_mpi,sdata_write_param,sdata_write_mpi
implicit none

private
public :: nicas_driver

contains

!----------------------------------------------------------------------
! Subroutine: nicas_driver
!> Purpose: NICAS driver (separated from main program for OOPS interfacing)
!----------------------------------------------------------------------
subroutine nicas_driver(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0
logical :: same_mask

! Determine whether the model is regional
sdata%regional = all(sdata%area<4.0*pi)

! Check whether the mask is the same for all levels
same_mask = .true.
do il0=2,sdata%nl0
   same_mask = same_mask.and.(all((sdata%mask(:,il0).and.sdata%mask(:,1)).or.(.not.sdata%mask(:,il0).and..not.sdata%mask(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   sdata%nl0i = 1
else
   sdata%nl0i = sdata%nl0
end if
write(mpl%unit,'(a7,a,i3)') '','Number of independent levels: ',sdata%nl0i

! Initialize sdata parameters from namelist
sdata%nproc = max(nam%nproc,1)
sdata%mpicom = nam%mpicom 

if (nam%new_param) then
   !----------------------------------------------------------------------
   ! Compute NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS parameters'

   call compute_parameters(sdata)

   !----------------------------------------------------------------------
   ! Compute NICAS normalization
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS normalization'

   call compute_normalization(sdata)

   !----------------------------------------------------------------------
   ! Write NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write NICAS parameters'

   call sdata_write_param(sdata)
else
   !----------------------------------------------------------------------
   ! Read NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS parameters'

   call sdata_read_param(sdata)
end if

if (nam%new_mpi) then
   !----------------------------------------------------------------------
   ! Read NICAS local distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS local distribution'

   call sdata_read_local(sdata)

   !----------------------------------------------------------------------
   ! Compute NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS MPI distribution'

   sdata%mpicom = nam%mpicom
   call compute_mpi(sdata)

   !----------------------------------------------------------------------
   ! Write NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write NICAS MPI distribution'

   call sdata_write_mpi(sdata)
else
   !----------------------------------------------------------------------
   ! Read NICAS local distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS local distribution'

   call sdata_read_local(sdata)

   !----------------------------------------------------------------------
   ! Read NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS MPI distribution'

   call sdata_read_mpi(sdata)
end if

if (nam%check_adjoints) then
   !----------------------------------------------------------------------
   ! Test adjoints
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test adjoints'

   call test_adjoints(sdata)
end if

if (nam%check_pos_def) then
   !----------------------------------------------------------------------
   ! Test NICAS positive definiteness
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS positive definiteness'

   call test_pos_def(sdata)
end if

if (nam%check_mpi.and.(nam%nproc>0)) then
   !----------------------------------------------------------------------
   ! Test single/multi-procs equivalence
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test single/multi-procs equivalence'

  call test_mpi(sdata)
end if

if (nam%check_dirac) then
   !----------------------------------------------------------------------
   ! Apply NICAS to diracs
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'

   call test_dirac(sdata)
end if

if (.true.) then
   !----------------------------------------------------------------------
   ! Test NICAS performance
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS performance'

   call test_perf(sdata)
end if

end subroutine nicas_driver

end module module_nicas
