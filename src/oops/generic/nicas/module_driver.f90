!----------------------------------------------------------------------
! Module: module_driver
!> Purpose: nicas driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_driver

use module_mpi, only: compute_mpi
use module_namelist, only: nam
use module_normalization, only: compute_normalization
use module_parameters, only: compute_parameters
use module_test, only: test_adjoints,test_pos_def,test_mpi,test_dirac,test_perf
use tools_const, only: pi
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype,ndataloc_dealloc,ndataloc_copy, &
& ndata_read_param,ndata_read_local,ndata_read_mpi, &
& ndata_write_param,ndata_write_mpi,ndata_write_mpi_summary

implicit none

private
public :: nicas_driver

contains

!----------------------------------------------------------------------
! Subroutine: nicas_driver
!> Purpose: NICAS driver (separated from main program for OOPS interfacing)
!----------------------------------------------------------------------
subroutine nicas_driver(ndata,ndataloc)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata             !< Sampling data
type(ndataloctype),intent(inout) :: ndataloc !< Sampling data,local

! Local variables
integer :: il0,iproc
logical :: same_mask
type(ndataloctype),allocatable :: ndataloc_arr(:)

! Check whether the mask is the same for all levels
same_mask = .true.
do il0=2,ndata%nl0
   same_mask = same_mask.and.(all((ndata%mask(:,il0).and.ndata%mask(:,1)).or.(.not.ndata%mask(:,il0).and..not.ndata%mask(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   ndata%nl0i = 1
else
   ndata%nl0i = ndata%nl0
end if
write(mpl%unit,'(a7,a,i3)') '','Number of independent levels: ',ndata%nl0i

if (nam%new_param) then
   !----------------------------------------------------------------------
   ! Compute NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS parameters'

   call compute_parameters(ndata)

   !----------------------------------------------------------------------
   ! Compute NICAS normalization
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS normalization'

   call compute_normalization(ndata)

   if (mpl%main) then
      !----------------------------------------------------------------------
      ! Write NICAS parameters
      !----------------------------------------------------------------------

      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Write NICAS parameters'

      call ndata_write_param(ndata)
   end if
else
   !----------------------------------------------------------------------
   ! Read NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS parameters'

   call ndata_read_param(ndata)
end if
call flush(mpl%unit)

!----------------------------------------------------------------------
! Read NICAS local distribution
!----------------------------------------------------------------------

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Read NICAS local distribution'

call ndata_read_local(ndata)
call flush(mpl%unit)

if (nam%new_mpi) then
   ! Allocation
   allocate(ndataloc_arr(nam%nproc))

   !----------------------------------------------------------------------
   ! Compute NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS MPI distribution'

   call compute_mpi(ndata,ndataloc_arr)

   if (mpl%main) then
      !----------------------------------------------------------------------
      ! Write NICAS MPI distribution
      !----------------------------------------------------------------------

      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Write NICAS MPI distribution'

      call ndata_write_mpi(ndataloc_arr)

      !----------------------------------------------------------------------
      ! Write NICAS MPI summary
      !----------------------------------------------------------------------
  
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Write NICAS MPI summary'

      call ndata_write_mpi_summary(ndata,ndataloc_arr)
   end if

   do iproc=1,nam%nproc
      ! Copy ndataloc for the concerned processor
      if (iproc==mpl%myproc) call ndataloc_copy(ndataloc_arr(iproc),ndataloc)

      ! Release memory
      call ndataloc_dealloc(ndataloc_arr(iproc))     
   end do

   ! Release memory
   deallocate(ndataloc_arr)
else
   !----------------------------------------------------------------------
   ! Read NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS MPI distribution'

   call ndata_read_mpi(ndataloc)
end if
call flush(mpl%unit)

if (nam%check_adjoints) then
   !----------------------------------------------------------------------
   ! Test adjoints
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test adjoints'

   if (mpl%main) call test_adjoints(ndata)
end if
call flush(mpl%unit)

if (nam%check_pos_def) then
   !----------------------------------------------------------------------
   ! Test NICAS positive definiteness
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS positive definiteness'

   if (mpl%main) call test_pos_def(ndata)
end if
call flush(mpl%unit)

if (nam%check_mpi.and.(nam%nproc>0)) then
   !----------------------------------------------------------------------
   ! Test single/multi-procs equivalence
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test single/multi-procs equivalence'

  call test_mpi(ndata,ndataloc)
end if
call flush(mpl%unit)

if (nam%check_dirac) then
   !----------------------------------------------------------------------
   ! Apply NICAS to diracs
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'

   call test_dirac(ndata,ndataloc)
end if
call flush(mpl%unit)

if (nam%check_perf) then
   !----------------------------------------------------------------------
   ! Test NICAS performance
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS performance'

   call test_perf(ndataloc)
end if
call flush(mpl%unit)

end subroutine nicas_driver

end module module_driver
