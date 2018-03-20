!----------------------------------------------------------------------
! Module: driver_lct
!> Purpose: LCT driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_lct

use tools_kinds, only: kind_real
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_lct, only: lct_type
use type_mom, only: mom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

logical,parameter :: write_cor = .true. !< Write raw and fitted correlations

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
type(nam_type),target,intent(inout) :: nam                                                 !< Namelist
type(geom_type),target,intent(in) :: geom                                                  !< Geometry
type(bpar_type),target,intent(in) :: bpar                                                  !< Block parameters
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
type(hdata_type) :: hdata
type(mom_type) :: mom
type(lct_type) :: lct

if (nam%new_lct) then
   ! Setup sampling
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'

   ! Set artificially small local radius
   nam%local_rad = 1.0e-12

   ! Setup sampling
   call hdata%setup_sampling(nam,geom)
   call flush(mpl%unit)

   ! Compute MPI distribution, halo A
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A'

   call hdata%compute_mpi_a(nam,geom)
   call flush(mpl%unit)

   ! Compute MPI distribution, halos A-B
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A-B'

   call hdata%compute_mpi_ab(geom)
   call flush(mpl%unit)

   ! Compute MPI distribution, halo C
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute MPI distribution, halo C'

   call hdata%compute_mpi_c(nam,geom)
   call flush(mpl%unit)

   ! Compute sample moments
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute sample moments'

   if (present(ens1)) then
      call mom%compute(nam,geom,bpar,hdata,'ens1',ens1)
   else
      call mom%compute(nam,geom,bpar,hdata,'ens1')
   end if

   ! Compute LCT
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute LCT'

   call lct%compute(nam,geom,bpar,hdata,mom)

   ! Filter LCT
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Filter LCT'

   call lct%filter(nam,geom,bpar,hdata)

   ! LCT RMSE
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- LCT RMSE'

   call lct%rmse(nam,geom,bpar,hdata)

   ! Write LCT
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write LCT'

   call lct%write(nam,geom,bpar,hdata)

   if (write_cor) then
      ! Write correlation and LCT fit
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Write correlation and LCT fit'

      call lct%write_cor(nam,geom,bpar,hdata)
   end if
end if

end subroutine run_lct

end module driver_lct
