!----------------------------------------------------------------------
! Module: driver_test
!> Purpose: test driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_test

use module_test, only: test_loc_adjoint,test_nicas_dirac,test_loc_dirac,test_loc_ens_dirac,test_nicas_perf,test_hdiag
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use type_bdata, only: bdatatype
use type_bpar, only: bpartype
use type_geom, only: geomtype
use type_mpl, only: mpl
use type_ndata, only: ndataloctype
use type_nam, only: namtype

implicit none

private
public :: run_test

contains

!----------------------------------------------------------------------
! Subroutine: run_test
!> Purpose: test NICAS method
!----------------------------------------------------------------------
subroutine run_test(nam,geom,bpar,bdata,ndataloc,ens1)

implicit none

! Passed variables
type(namtype),intent(inout) :: nam                                                         !< Namelist
type(geomtype),intent(in) :: geom                                                          !< Geometry
type(bpartype),intent(in) :: bpar                                                          !< Block parameters
type(bdatatype),intent(in) :: bdata(bpar%nb+1)                                             !< B data
type(ndataloctype),intent(in) :: ndataloc(bpar%nb+1)                                       !< NICAS data,local
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ib

if (nam%check_dirac) then
   ! Apply NICAS to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'
   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a7,a)') '','Dirac test for block: '//trim(bpar%blockname(ib))
         call test_nicas_dirac(nam,geom,trim(bpar%blockname(ib)),ndataloc(ib))
      end if
   end do
   call flush(mpl%unit)
end if

if (nam%check_perf) then
   ! Test NICAS performance
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS performance'
   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a7,a)') '','Performance results (elapsed time) for block: '//trim(bpar%blockname(ib))
         call test_nicas_perf(nam,geom,ndataloc(ib))
      end if
   end do
   call flush(mpl%unit)
end if

if (nam%check_adjoints) then
   ! Test adjoints
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test localization adjoint'
   call test_loc_adjoint(nam,geom,bpar,ndataloc)
   call flush(mpl%unit)
end if

if (nam%check_dirac) then
   ! Apply localization to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply localization to diracs'
   call test_loc_dirac(nam,geom,bpar,ndataloc)
   call flush(mpl%unit)

   if (present(ens1)) then
      ! Apply localized ensemble covariance to diracs
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Apply localized ensemble covariance  to diracs'
      call test_loc_ens_dirac(nam,geom,bpar,ndataloc,ens1)
      call flush(mpl%unit)
   end if
end if

if (nam%check_hdiag) then
   ! Test hdiag consistency
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test hdiag consistency'
   call test_hdiag(nam,geom,bpar,bdata,ndataloc)
   call flush(mpl%unit)
end if

end subroutine run_test

end module driver_test
