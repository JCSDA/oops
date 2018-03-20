!----------------------------------------------------------------------
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

use tools_const, only: pi
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_cmat, only: cmat_type
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_mpl, only: mpl
use type_nam, only: nam_type
use type_nicas, only: nicas_type

implicit none

private
public :: run_nicas

contains

!----------------------------------------------------------------------
! Subroutine: run_nicas
!> Purpose: NICAS
!----------------------------------------------------------------------
subroutine run_nicas(nam,geom,bpar,cmat,nicas,ens1)

implicit none

! Passed variables
type(nam_type),intent(inout) :: nam                                                        !< Namelist
type(geom_type),intent(inout) :: geom                                                      !< Geometry
type(bpar_type),intent(in) :: bpar                                                         !< Block parameters
type(cmat_type),intent(in) :: cmat                                                         !< C matrix data
type(nicas_type),intent(out) :: nicas                                                      !< NICAS data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ib,ic0,ic0a

! Allocation
call nicas%alloc(nam,bpar,'nicas')

if (nam%new_param) then
   ! Compute NICAS parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS parameters'

   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib).or.((ib==bpar%nb+1).and.nam%displ_diag)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))
      end if

      ! NICAS parameters
      if (bpar%nicas_block(ib)) call nicas%blk(ib)%compute_parameters(nam,geom,cmat%blk(ib))

      ! Advection
      if ((ib==bpar%nb+1).and.nam%displ_diag) call nicas%blk(ib)%compute_adv(nam,geom,cmat%blk(ib))

      if (bpar%B_block(ib)) then
         ! Copy weights
         nicas%blk(ib)%wgt = cmat%blk(ib)%wgt
         if (bpar%nicas_block(ib)) then
            allocate(nicas%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
            do ic0=1,geom%nc0
               if (geom%c0_to_proc(ic0)==mpl%myproc) then
                  ic0a = geom%c0_to_c0a(ic0)
                  nicas%blk(ib)%coef_ens(ic0a,:) = cmat%blk(ib)%coef_ens(ic0,:)
               end if
            end do
         end if
      end if
   end do

   ! Write NICAS parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write NICAS parameters'
   call nicas%write(nam,geom,bpar)

   ! Write NICAS MPI summary
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write NICAS MPI summary'
   call nicas%write_mpi_summary(nam,geom,bpar)
elseif (nam%new_param.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_sqrt.or.nam%check_dirac.or. &
 & nam%check_randomization.or.nam%check_consistency.or.nam%check_optimality) then
   ! Read NICAS parameters
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS parameters'
   call nicas%read(nam,geom,bpar)
   call flush(mpl%unit)
end if

if (nam%check_adjoints) then
   ! Test adjoint
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS adjoint'

   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))

         call nicas%blk(ib)%test_adjoint(nam,geom)
         call flush(mpl%unit)
      end if
   end do

   ! Test localization adjoint
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test localization adjoint'
   if (present(ens1)) then
      call nicas%test_adjoint(nam,geom,bpar,ens1)
   else
      call nicas%test_adjoint(nam,geom,bpar)
   end if
   call flush(mpl%unit)
end if

if (nam%check_pos_def) then
   ! Test NICAS positive definiteness
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS positive definiteness'

   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))

         call nicas%blk(ib)%test_pos_def(nam,geom)
         call flush(mpl%unit)
      end if
   end do
end if

if (nam%check_sqrt) then
   ! Test NICAS full/square-root equivalence
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS full/square-root equivalence'

   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))

         call nicas%blk(ib)%test_sqrt(nam,geom,bpar,cmat%blk(ib))
         call flush(mpl%unit)
      end if
   end do

   ! Test localization full/square-root equivalence
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test localization full/square-root equivalence'
   if (present(ens1)) then
      call nicas%test_sqrt(nam,geom,bpar,cmat,ens1)
   else
      call nicas%test_sqrt(nam,geom,bpar,cmat)
   end if
   call flush(mpl%unit)
end if

if (nam%check_dirac) then
   ! Apply NICAS to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'

   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))

         call nicas%blk(ib)%test_dirac(nam,geom,bpar)
         call flush(mpl%unit)
      end if
   end do

   ! Apply localization to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply localization to diracs'
   if (present(ens1)) then
      call nicas%test_dirac(nam,geom,bpar,ens1)
   else
      call nicas%test_dirac(nam,geom,bpar)
   end if
   call flush(mpl%unit)
end if

if (nam%check_randomization) then
   ! Test NICAS randomization
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS randomization'

   call nicas%test_randomization(nam,geom,bpar)
   call flush(mpl%unit)
end if

if (nam%check_consistency) then
   ! Test HDIAG_NICAS consistency
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test HDIAG_NICAS consistency'

   call nicas%test_consistency(nam,geom,bpar,cmat)
   call flush(mpl%unit)
end if

if (nam%check_optimality) then
   ! Test HDIAG optimality
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test HDIAG optimality'

   call nicas%test_optimality(nam,geom,bpar)
   call flush(mpl%unit)
end if

end subroutine run_nicas

end module driver_nicas
