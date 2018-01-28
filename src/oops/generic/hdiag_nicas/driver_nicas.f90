!---------------------------------------------------------------------
! Module: driver_nicas
!> Purpose: nicas driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_nicas

use nicas_parameters, only: compute_parameters
use nicas_test, only: test_nicas_adjoint,test_nicas_pos_def,test_nicas_sqrt,test_loc_adjoint, &
 & test_loc_sqrt,test_nicas_dirac,test_loc_dirac,test_randomization,test_consistency,test_optimality
use tools_const, only: pi
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_bdata, only: bdatatype
use type_bpar, only: bpartype
use type_geom, only: geomtype
use type_mpl, only: mpl
use type_nam, only: namtype
use type_ndata, only: ndatatype,ndata_read,ndata_write,ndata_write_mpi_summary

implicit none

private
public :: run_nicas

contains

!----------------------------------------------------------------------
! Subroutine: run_nicas
!> Purpose: NICAS
!----------------------------------------------------------------------
subroutine run_nicas(nam,geom,bpar,bdata,ndata,ens1)

implicit none

! Passed variables
type(namtype),target,intent(inout) :: nam             !< Namelist
type(geomtype),target,intent(inout) :: geom           !< Geometry
type(bpartype),target,intent(in) :: bpar              !< Block parameters
type(bdatatype),intent(in) :: bdata(bpar%nb+1)        !< B data
type(ndatatype),allocatable,intent(inout) :: ndata(:) !< NICAS data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ib,ic0,ic0a

! Allocate ndata
allocate(ndata(bpar%nb+1))
do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      ! Set name, namelist and geometry
      write(ndata(ib)%cname,'(a,i1,a,i4.4,a,i4.4,a,a)') 'ndata_',nam%mpicom,'_',mpl%nproc,'-',mpl%myproc, &
    & '_',trim(bpar%blockname(ib))
      ndata(ib)%nam => nam
      ndata(ib)%geom => geom
   end if
end do

do ib=1,bpar%nb+1
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))

   if (nam%new_param) then
      if (bpar%nicas_block(ib)) then
         ! Compute NICAS parameters
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Compute NICAS parameters'
         call compute_parameters(bdata(ib),ndata(ib))
      end if

      if (bpar%B_block(ib)) then
         ! Copy weights
         ndata(ib)%wgt = bdata(ib)%wgt
         if (bpar%nicas_block(ib)) then
            allocate(ndata(ib)%coef_ens(geom%nc0a,geom%nl0))
            do ic0=1,geom%nc0
               if (geom%c0_to_proc(ic0)==mpl%myproc) then
                  ic0a = geom%c0_to_c0a(ic0)
                  ndata(ib)%coef_ens(ic0a,:) = bdata(ib)%coef_ens(ic0,:)
               end if
            end do
         end if

         ! Write NICAS parameters
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Write NICAS parameters'
         call ndata_write(nam,geom,ndata(ib),bpar%nicas_block(ib),bpar%auto_block(ib))

         ! Write NICAS MPI summary
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Write NICAS MPI summary'
         call ndata_write_mpi_summary(ndata(ib))
      end if
   elseif (nam%new_param.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_sqrt.or.nam%check_dirac.or. &
 & nam%check_consistency.or.nam%check_optimality) then
      if (bpar%B_block(ib)) then
         ! Read NICAS parameters
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Read NICAS parameters'
         call ndata_read(nam,geom,ndata(ib),bpar%nicas_block(ib),bpar%auto_block(ib))
      end if
      call flush(mpl%unit)
   end if

   if (bpar%B_block(ib).and.nam%transform.and.bpar%auto_block(ib)) then
      ! Allocation
      allocate(ndata(ib)%trans(geom%nl0,geom%nl0))
      allocate(ndata(ib)%transinv(geom%nl0,geom%nl0))

      ! Copy
      ndata(ib)%trans = bdata(ib)%trans
      ndata(ib)%transinv = bdata(ib)%transinv
   end if

   if (bpar%nicas_block(ib)) then
      if (nam%check_adjoints) then
         ! Test adjoint
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Test NICAS adjoint'
         call test_nicas_adjoint(ndata(ib))
         call flush(mpl%unit)
      end if

      if (nam%check_pos_def) then
         ! Test NICAS positive definiteness
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Test NICAS positive definiteness'
         call test_nicas_pos_def(ndata(ib))
         call flush(mpl%unit)
      end if

      if (nam%check_sqrt) then
         ! Test NICAS full/square-root equivalence
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Test NICAS full/square-root equivalence'
         call test_nicas_sqrt(trim(bpar%blockname(ib)),bdata(ib),ndata(ib))
         call flush(mpl%unit)
      end if

      if (nam%check_dirac) then
         ! Apply NICAS to diracs
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Apply NICAS to diracs'
         write(mpl%unit,'(a7,a)') '','Dirac test for block: '//trim(bpar%blockname(ib))
         call test_nicas_dirac(trim(bpar%blockname(ib)),ndata(ib))
         call flush(mpl%unit)
      end if
   end if
end do

if (nam%check_adjoints) then
   ! Test localization adjoint
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test localization adjoint'
   if (present(ens1)) then
      call test_loc_adjoint(nam,geom,bpar,ndata,ens1)
   else
      call test_loc_adjoint(nam,geom,bpar,ndata)
   end if
   call flush(mpl%unit)
end if

if (nam%check_sqrt) then
   ! Test localization full/square-root equivalence
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test localization full/square-root equivalence'
   if (present(ens1)) then
      call test_loc_sqrt(nam,geom,bpar,bdata,ndata,ens1)
   else
      call test_loc_sqrt(nam,geom,bpar,bdata,ndata)
   end if
   call flush(mpl%unit)
end if

if (nam%check_dirac) then
   ! Apply localization to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply localization to diracs'
   if (present(ens1)) then
      call test_loc_dirac(nam,geom,bpar,ndata,ens1)
   else
      call test_loc_dirac(nam,geom,bpar,ndata)
   end if
   call flush(mpl%unit)
end if

if (nam%check_consistency.or.nam%check_optimality) then
   ! Test randomization
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test randomization'
   call test_randomization(nam,geom,bpar,ndata)
   call flush(mpl%unit)
end if

if (nam%check_consistency) then
   ! Test HDIAG_NICAS consistency
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test HDIAG_NICAS consistency'
   call test_consistency(nam,geom,bpar,bdata,ndata)
   call flush(mpl%unit)
end if

if (nam%check_optimality) then
   ! Test HDIAG optimality
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test HDIAG optimality'
   call test_optimality(nam,geom,bpar,ndata)
   call flush(mpl%unit)
end if

end subroutine run_nicas

end module driver_nicas
