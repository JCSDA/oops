!----------------------------------------------------------------------
! Module: driver_hdiag
!> Purpose: HDIAG driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_hdiag

use model_interface, only: model_write
use netcdf
use tools_const, only: reqkm
use tools_display, only: msgerror,msgwarning
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msr,isnotmsi,isnotmsr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_avg, only: avg_type
use type_cmat, only: cmat_type
use type_bpar, only: bpar_type
use type_diag, only: diag_type
use type_displ, only: displ_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_mom, only: mom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

private
public :: run_hdiag

contains

!----------------------------------------------------------------------
! Subroutine: run_hdiag
!> Purpose: HDIAG
!----------------------------------------------------------------------
subroutine run_hdiag(nam,geom,bpar,cmat,ens1)

implicit none

! Passed variables
type(nam_type),target,intent(inout) :: nam                                                 !< Namelist
type(geom_type),target,intent(in) :: geom                                                  !< Geometry
type(bpar_type),target,intent(in) :: bpar                                                  !< Block parameters
type(cmat_type),intent(out) :: cmat                                                        !< C matrix data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ib
character(len=1024) :: filename
type(avg_type) :: avg_1,avg_2,avg_wgt
type(diag_type) :: cov_1,cov_2,cor_1,cor_2,loc_1,loc_2,loc_3
type(displ_type) :: displ
type(hdata_type) :: hdata
type(mom_type) :: mom_1,mom_2

if (nam%new_hdiag) then
   ! Setup sampling
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'
   call hdata%setup_sampling(nam,geom)
   call flush(mpl%unit)

   ! Compute MPI distribution, halo A
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A'

   call hdata%compute_mpi_a(nam,geom)
   call flush(mpl%unit)

   if (nam%local_diag.or.nam%displ_diag) then
      ! Compute MPI distribution, halos A-B
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A-B'

      call hdata%compute_mpi_ab(geom)
      call flush(mpl%unit)
   end if

   if (nam%displ_diag) then
      ! Compute MPI distribution, halo D
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute MPI distribution, halo D'

      call hdata%compute_mpi_d(nam,geom)
      call flush(mpl%unit)

      ! Compute displacement diagnostic
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute displacement diagnostic'

      if (present(ens1)) then
         call displ%compute(nam,geom,hdata,ens1)
      else
        call displ%compute(nam,geom,hdata)
      end if
   end if
   call flush(mpl%unit)

   ! Compute MPI distribution, halo C
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute MPI distribution, halo C'

   call hdata%compute_mpi_c(nam,geom)
   call flush(mpl%unit)

   ! Compute sample moments
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute sample moments'

   ! Compute ensemble 1 sample moments
   write(mpl%unit,'(a7,a)') '','Ensemble 1:'
   if (present(ens1)) then
      call mom_1%compute(nam,geom,bpar,hdata,'ens1',ens1)
   else
      call mom_1%compute(nam,geom,bpar,hdata,'ens1')
   end if

   if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) then
      ! Compute randomized sample moments
      write(mpl%unit,'(a7,a)') '','Ensemble 2:'
      call mom_2%compute(nam,geom,bpar,hdata,'ens2')
   end if
   call flush(mpl%unit)

   ! Compute statistics
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute statistics'

   ! Compute ensemble 1 statistics
   write(mpl%unit,'(a7,a)') '','Ensemble 1:'
   call avg_1%compute(nam,geom,bpar,hdata,mom_1,nam%ne)

   if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) then
      ! Compute randomized sample moments
      write(mpl%unit,'(a7,a)') '','Ensemble 2:'
      call avg_2%compute(nam,geom,bpar,hdata,mom_2,nam%ens2_ne)
   end if
   call flush(mpl%unit)

   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd','dual-ens')
      ! Compute hybrid statistics
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute hybrid statistics'

      call avg_2%compute_hyb(nam,geom,bpar,hdata,mom_1,mom_2,avg_1)
   end select

   if (bpar%diag_block(bpar%nb+1)) then
      ! Compute block-averaged statistics
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute block-averaged statistics'

      ! Copy block-averaging weight
      avg_wgt = avg_1%copy_wgt(geom,bpar)

      call avg_1%compute_bwavg(nam,geom,bpar,avg_wgt)
      if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) call avg_2%compute_bwavg(nam,geom,bpar,avg_wgt)
   end if

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute covariance'

   ! Compute ensemble 1 covariance
   write(mpl%unit,'(a7,a)') '','Ensemble 1:'
   call cov_1%covariance(nam,geom,bpar,hdata,avg_1,'cov')

   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd','dual-ens')
      ! Compute ensemble 2 covariance
      write(mpl%unit,'(a7,a)') '','Ensemble 2:'
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd')
         call cov_2%covariance(nam,geom,bpar,hdata,avg_2,'cov_sta')
      case ('dual-ens')
         call cov_2%covariance(nam,geom,bpar,hdata,avg_2,'cov_lr')
      end select
      call flush(mpl%unit)
   end select

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute correlation'

   ! Compute ensemble 1 correlation
   write(mpl%unit,'(a7,a)') '','Ensemble 1:'
   call cor_1%correlation(nam,geom,bpar,hdata,avg_1,'cor')

   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd','dual-ens')
      ! Compute ensemble 2 correlation
      write(mpl%unit,'(a7,a)') '','Ensemble 2:'
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd')
         call cor_2%correlation(nam,geom,bpar,hdata,avg_2,'cor_sta')
      case ('dual-ens')
         call cor_2%correlation(nam,geom,bpar,hdata,avg_2,'cor_lr')
      end select
      call flush(mpl%unit)
   end select

   select case (trim(nam%method))
   case ('loc','hyb-avg','hyb-rnd','dual-ens')
      ! Compute localization
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute localization'
      write(mpl%unit,'(a7,a)') '','Ensemble 1:'
      call loc_1%localization(nam,geom,bpar,hdata,avg_1,'loc')
   end select
   call flush(mpl%unit)

   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      ! Compute static hybridization
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute static hybridization'
      write(mpl%unit,'(a7,a)') '','Ensemble 1 and 2:'
      call loc_2%hybridization(nam,geom,bpar,hdata,avg_1,avg_2,'loc_hyb')
   end select
   call flush(mpl%unit)

   if (trim(nam%method)=='dual-ens') then
      ! Compute dual-ensemble hybridization diagnostic and fit
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute dual-ensemble hybridization'
      write(mpl%unit,'(a7,a)') '','Ensembles 1 and 2:'
      call loc_2%dualens(nam,geom,bpar,hdata,avg_1,avg_2,loc_3,'loc_deh','loc_deh_lr')
   end if
   call flush(mpl%unit)

   if (trim(nam%minim_algo)/='none') then
      ! Copy diagnostics into C matrix data
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Copy diagnostics into C matrix data'

      select case (trim(nam%method))
      case ('cor')
         call cmat%from_diag(nam,geom,bpar,hdata,cor_1)
      case ('loc')
         call cmat%from_diag(nam,geom,bpar,hdata,loc_1)
      case default
         call msgerror('cmat not implemented yet for this method')
      end select

      if (mpl%main) then
         ! Write C matrix data
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a,i5,a)') '--- Write C matrix data'

         call cmat%write(nam,geom,bpar)
      end if
   end if
   call flush(mpl%unit)

   ! Write data
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write data'

   ! Displacement
   if (nam%displ_diag) call displ%write(nam,geom,hdata,trim(nam%prefix)//'_displ_diag.nc')

   ! Full variances
   if (nam%full_var) then
      filename = trim(nam%prefix)//'_full_var_gridded.nc'
      do ib=1,bpar%nb
         if (bpar%diag_block(ib)) call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_var', &
       & sum(mom_1%blk(ib)%m2full,dim=3)/float(mom_1%blk(ib)%nsub))
      end do
   end if
   call flush(mpl%unit)
elseif (nam%new_param) then
   ! Read C matrix data
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,i5,a)') '--- Read C matrix data'

   call cmat%read(nam,geom,bpar)
   call flush(mpl%unit)
end if

end subroutine run_hdiag

end module driver_hdiag
