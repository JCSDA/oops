!----------------------------------------------------------------------
! Module: driver_lct
!> Purpose: hybrid_diag driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_lct

use hdiag_moments, only: compute_moments
use hdiag_sampling, only: setup_sampling
use hdiag_lct, only: compute_lct
use tools_kinds, only: kind_real
use type_bpar, only: bpartype
use type_geom, only: geomtype
use type_hdata, only: hdatatype
use type_lct, only: lcttype,lct_write
use type_mom, only: momtype
use type_mpl, only: mpl
use type_nam, only: namtype

implicit none

private
public :: run_lct

contains

!----------------------------------------------------------------------
! Subroutine: run_lct
!> Purpose: LCT diagnostics
!----------------------------------------------------------------------
subroutine run_lct(nam,geom,bpar,ens1)

implicit none

! Passed variables
type(namtype),target,intent(inout) :: nam                                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                                   !< Geometry
type(bpartype),target,intent(in) :: bpar                                                   !< Block parameters
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
type(hdatatype) :: hdata
type(momtype) :: mom(bpar%nb)
type(lcttype) :: lct(nam%nc1,geom%nl0,bpar%nb)

if (nam%new_lct) then
   ! Set pointers
   hdata%nam => nam
   hdata%geom => geom
   hdata%bpar => bpar

   ! Setup sampling
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'

   ! Set artificially small local radius
   nam%local_rad = 1.0e-12

   ! Setup sampling
   call setup_sampling(hdata)

   ! Compute sample moments
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute sample moments'

   if (present(ens1)) then
      call compute_moments(hdata,'ens1',mom,ens1)
   else
      call compute_moments(hdata,'ens1',mom)
   end if

   ! Compute LCT
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute LCT'

   call compute_lct(hdata,mom,lct)

   ! Write LCT
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write LCT'

   call lct_write(hdata,lct)
end if

end subroutine run_lct

end module driver_lct
