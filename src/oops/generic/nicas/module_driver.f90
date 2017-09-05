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
use module_mpi_obsop, only: compute_mpi_obsop
use module_namelist, only: namtype
use module_normalization, only: compute_normalization
use module_parameters, only: compute_parameters
use module_parameters_obsop, only: compute_parameters_obsop
use module_test, only: test_adjoints,test_pos_def,test_mpi,test_dirac,test_perf
use module_test_obsop, only: test_adjoint_obsop,test_mpi_obsop,test_mpi_obsop_ad
use tools_const, only: pi
use type_geom, only: geomtype,geom_read_local,compute_grid_mesh
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype,ndata_read_param,ndata_read_mpi, &
  & ndata_write_param,ndata_write_mpi,ndata_write_mpi_summary
use type_odata, only: odatatype,odataloctype

implicit none

private
public :: nicas_driver,obsop_driver

contains

!----------------------------------------------------------------------
! Subroutine: nicas_driver
!> Purpose: NICAS driver
!----------------------------------------------------------------------
subroutine nicas_driver(nam,geom,ndataloc)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(geomtype),target,intent(in) :: geom    !< Sampling data
type(ndataloctype),intent(inout) :: ndataloc !< Sampling data,local

! Local variables
integer :: il0
logical :: same_mask
type(ndatatype) :: ndata

! Set ndata geometry
ndata%geom => geom
ndata%nc0 = ndata%geom%nc0
ndata%nl0 = ndata%geom%nl0

! Compute grid mesh
call compute_grid_mesh(nam,ndata%geom)

! Check whether the mask is the same for all levels
same_mask = .true.
do il0=2,ndata%nl0
   same_mask = same_mask.and.(all((ndata%geom%mask(:,il0).and.ndata%geom%mask(:,1)) &
             & .or.(.not.ndata%geom%mask(:,il0).and..not.ndata%geom%mask(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   ndata%nl0i = 1
else
   ndata%nl0i = ndata%nl0
end if
write(mpl%unit,'(a7,a,i3)') '','Number of independent levels: ',ndata%nl0i

if (nam%new_param) then
   ! Compute NICAS parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS parameters'
   call compute_parameters(nam,ndata)

   ! Compute NICAS normalization
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS normalization'
   call compute_normalization(nam,ndata)

   if (mpl%main) then
      ! Write NICAS parameters
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Write NICAS parameters'
      call ndata_write_param(nam,ndata)
   end if
else
   ! Read NICAS parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS parameters'
   call ndata_read_param(nam,ndata)
end if
call flush(mpl%unit)

! Read local distribution
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Read local distribution'
call geom_read_local(nam,ndata%geom)
call flush(mpl%unit)

if (nam%new_mpi) then
   ! Compute NICAS MPI distribution
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS MPI distribution'
   call compute_mpi(nam,ndata,ndataloc)

   ! Write NICAS MPI distribution
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write NICAS MPI distribution'
   call ndata_write_mpi(nam,ndataloc)

   if (mpl%main.and.(nam%nproc>1)) then
      ! Write NICAS MPI summary
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Write NICAS MPI summary'
      call ndata_write_mpi_summary(nam,ndata)
   end if
else
   ! Read NICAS MPI distribution
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS MPI distribution'
   call ndata_read_mpi(nam,ndataloc)
end if
call flush(mpl%unit)

if (nam%check_adjoints) then
   ! Test adjoints
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS adjoints'
   if (mpl%main) call test_adjoints(nam,ndata)
   call flush(mpl%unit)
end if

if (nam%check_pos_def) then
   ! Test NICAS positive definiteness
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS positive definiteness'
   if (mpl%main) call test_pos_def(nam,ndata)
   call flush(mpl%unit)
end if

if (nam%check_mpi.and.(nam%nproc>0)) then
   ! Test single/multi-procs equivalence
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS single/multi-procs equivalence'
   call test_mpi(nam,ndata,ndataloc)
   call flush(mpl%unit)
end if

if (nam%check_dirac) then
   ! Apply NICAS to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'
   call test_dirac(nam,ndata,ndataloc)
   call flush(mpl%unit)
end if

if (nam%check_perf) then
   ! Test NICAS performance
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS performance'
   call test_perf(nam,ndataloc)
   call flush(mpl%unit)
end if

end subroutine nicas_driver

!----------------------------------------------------------------------
! Subroutine: obsop_driver
!> Purpose: observation operator driver
!----------------------------------------------------------------------
subroutine obsop_driver(nam,geom,odataloc)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(geomtype),target,intent(inout) :: geom         !< Sampling data
type(odataloctype),intent(inout) :: odataloc !< Sampling data

! Local variables
type(odatatype) :: odata

! Set odata geometry
odata%geom => geom
odata%nc0 = odata%geom%nc0
odata%nl0 = odata%geom%nl0

!if (nam%new_param) then
   ! Compute observation operator parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute observation operator parameters'
   call compute_parameters_obsop(odata)

!   if (mpl%main) then
!      ! Write observation operator parameters
!      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
!      write(mpl%unit,'(a)') '--- Write observation operator parameters'
!      call odata_write_param(odata)
!   end if
!else
!   ! Read observation operator parameters
!   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
!   write(mpl%unit,'(a)') '--- Read observation operator parameters'
!   call odata_read_param(odata)
!end if
call flush(mpl%unit)

! Read local distribution
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Read local distribution'
call geom_read_local(nam,geom)
call flush(mpl%unit)

!if (nam%new_mpi) then
   ! Compute observation operator MPI distribution
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute observation operator MPI distribution'
   call compute_mpi_obsop(nam,odata,odataloc)

!   if (mpl%main) then
!      ! Write observation operator MPI distribution
!      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
!      write(mpl%unit,'(a)') '--- Write observation operator MPI distribution'
!      call odata_write_mpi(odataloc_arr)
!   end if
!else
!   ! Read observation operator MPI distribution
!   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
!   write(mpl%unit,'(a)') '--- Read observation operator MPI distribution'
!   call odata_read_mpi(odataloc)
!end if
call flush(mpl%unit)

if (nam%check_adjoints) then
   ! Test adjoints
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test observation operator adjoint'
   call test_adjoint_obsop(odata)
   call flush(mpl%unit)
end if

if (nam%check_mpi.and.(nam%nproc>0)) then
   ! Test single/multi-procs equivalence
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test observation operator single/multi-procs equivalence'
   call test_mpi_obsop(odata,odataloc)
   call test_mpi_obsop_ad(odata,odataloc)
   call flush(mpl%unit)
end if

end subroutine obsop_driver

end module module_driver
