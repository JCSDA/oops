!----------------------------------------------------------------------
! Module: type_cmat
!> Purpose: C matrix derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_cmat

use netcdf
use tools_const, only: rad2deg,reqkm
use tools_display, only: msgwarning,msgerror,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msr,isnotmsi,isnotmsr,isallnotmsr,isanynotmsr
use type_avg, only: avg_type
use type_bpar, only: bpar_type
use type_cmat_blk, only: cmat_blk_type
use type_diag, only: diag_type
use type_displ, only: displ_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_io, only: io_type
use type_mom, only: mom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! C matrix data derived type
type cmat_type
   character(len=1024) :: prefix             !< Prefix
   type(cmat_blk_type),allocatable :: blk(:) !< C matrix blocks
contains
   procedure :: alloc => cmat_alloc
   procedure :: copy => cmat_copy
   procedure :: read => cmat_read
   procedure :: write => cmat_write
   procedure :: run_hdiag => cmat_run_hdiag
   procedure :: cmat_from_fields
   procedure :: cmat_from_diag
   generic :: from => cmat_from_fields,cmat_from_diag
end type cmat_type

private
public :: cmat_blk_type,cmat_type

contains

!----------------------------------------------------------------------
! Subroutine: cmat_alloc
!> Purpose: cmat object allocation
!----------------------------------------------------------------------
subroutine cmat_alloc(cmat,nam,geom,bpar,prefix)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat    !< C matrix data
type(nam_type),target,intent(in) :: nam   !< Namelist
type(geom_type),target,intent(in) :: geom !< Geometry
type(bpar_type),intent(in) :: bpar        !< Block parameters
character(len=*),intent(in) :: prefix     !< Prefix

! Local variables
integer :: ib

! Copy prefix
cmat%prefix = prefix

! Allocation
allocate(cmat%blk(bpar%nb+1))

do ib=1,bpar%nb+1
   ! Set block name
   cmat%blk(ib)%name = trim(prefix)//'_'//trim(bpar%blockname(ib))

   if (bpar%diag_block(ib)) then
      ! Allocation
      allocate(cmat%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
      allocate(cmat%blk(ib)%coef_sta(geom%nc0a,geom%nl0))
      allocate(cmat%blk(ib)%rh0(geom%nc0a,geom%nl0))
      allocate(cmat%blk(ib)%rv0(geom%nc0a,geom%nl0))
      allocate(cmat%blk(ib)%rh0s(geom%nc0a,geom%nl0))
      allocate(cmat%blk(ib)%rv0s(geom%nc0a,geom%nl0))

      ! Initialization
      call msr(cmat%blk(ib)%coef_ens)
      call msr(cmat%blk(ib)%coef_sta)
      call msr(cmat%blk(ib)%rh0)
      call msr(cmat%blk(ib)%rv0)
      call msr(cmat%blk(ib)%rh0s)
      call msr(cmat%blk(ib)%rv0s)
      call msr(cmat%blk(ib)%wgt)
   end if

   if ((ib==bpar%nb+1).and.nam%displ_diag) then
      ! Allocation
      allocate(cmat%blk(ib)%displ_lon(geom%nc0a,geom%nl0,2:nam%nts))
      allocate(cmat%blk(ib)%displ_lat(geom%nc0a,geom%nl0,2:nam%nts))

      ! Initialization
      if (nam%displ_diag) then
         call msr(cmat%blk(ib)%displ_lon)
         call msr(cmat%blk(ib)%displ_lat)
      end if
   end if
end do

end subroutine cmat_alloc

!----------------------------------------------------------------------
! Subroutine: cmat_copy
!> Purpose: cmat object copy
!----------------------------------------------------------------------
type(cmat_type) function cmat_copy(cmat,nam,geom,bpar)

implicit none

! Passed variables
class(cmat_type),intent(in) :: cmat       !< C matrix data
type(nam_type),target,intent(in) :: nam   !< Namelist
type(geom_type),target,intent(in) :: geom !< Geometry
type(bpar_type),intent(in) :: bpar        !< Block parameters

! Local variables
integer :: ib

! Allocation
call cmat_copy%alloc(nam,geom,bpar,trim(cmat%prefix))

! Copy
do ib=1,bpar%nb+1
   if (allocated(cmat%blk(ib)%coef_ens)) cmat_copy%blk(ib)%coef_ens = cmat%blk(ib)%coef_ens
   if (allocated(cmat%blk(ib)%coef_sta)) cmat_copy%blk(ib)%coef_sta = cmat%blk(ib)%coef_sta
   if (allocated(cmat%blk(ib)%rh0)) cmat_copy%blk(ib)%rh0 = cmat%blk(ib)%rh0
   if (allocated(cmat%blk(ib)%rv0)) cmat_copy%blk(ib)%rv0 = cmat%blk(ib)%rv0
   if (allocated(cmat%blk(ib)%rh0s)) cmat_copy%blk(ib)%rh0s = cmat%blk(ib)%rh0s
   if (allocated(cmat%blk(ib)%rv0s)) cmat_copy%blk(ib)%rv0s = cmat%blk(ib)%rv0s
   if (allocated(cmat%blk(ib)%displ_lon)) cmat_copy%blk(ib)%displ_lon = cmat%blk(ib)%displ_lon
   if (allocated(cmat%blk(ib)%displ_lat)) cmat_copy%blk(ib)%displ_lat = cmat%blk(ib)%displ_lat
end do

end function cmat_copy

!----------------------------------------------------------------------
! Subroutine: cmat_read
!> Purpose: read cmat object
!----------------------------------------------------------------------
subroutine cmat_read(cmat,nam,geom,bpar,io)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat !< C matrix data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(io_type),intent(in) :: io         !< I/O

! Local variables
integer :: ib
character(len=1024) :: filename

! Allocation
call cmat%alloc(nam,geom,bpar,'cmat')

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      filename = trim(nam%prefix)//'_'//trim(cmat%blk(ib)%name)//'.nc'
      if (bpar%nicas_block(ib)) then
         call io%fld_read(nam,geom,filename,'coef_ens',cmat%blk(ib)%coef_ens)
         call io%fld_read(nam,geom,filename,'coef_sta',cmat%blk(ib)%coef_sta)
         call io%fld_read(nam,geom,filename,'rh0',cmat%blk(ib)%rh0)
         call io%fld_read(nam,geom,filename,'rv0',cmat%blk(ib)%rv0)
         call io%fld_read(nam,geom,filename,'rh0s',cmat%blk(ib)%rh0s)
         call io%fld_read(nam,geom,filename,'rv0s',cmat%blk(ib)%rv0s)
      end if
      if ((ib==bpar%nb+1).and.nam%displ_diag) then
         call io%fld_read(nam,geom,filename,'displ_lon',cmat%blk(ib)%displ_lon)
         call io%fld_read(nam,geom,filename,'displ_lat',cmat%blk(ib)%displ_lat)
      end if

      ! Check
      if (bpar%nicas_block(ib)) then
         if (any((cmat%blk(ib)%rh0<0.0).and.isnotmsr(cmat%blk(ib)%rh0))) call msgerror('rh0 should be positive')
         if (any((cmat%blk(ib)%rv0<0.0).and.isnotmsr(cmat%blk(ib)%rv0))) call msgerror('rv0 should be positive')
         if (any((cmat%blk(ib)%rh0s<0.0).and.isnotmsr(cmat%blk(ib)%rh0s))) call msgerror('rh0s should be positive')
         if (any((cmat%blk(ib)%rv0s<0.0).and.isnotmsr(cmat%blk(ib)%rv0s))) call msgerror('rv0s should be positive')
      end if
   end if
end do

end subroutine cmat_read

!----------------------------------------------------------------------
! Subroutine: cmat_write
!> Purpose: write cmat object
!----------------------------------------------------------------------
subroutine cmat_write(cmat,nam,geom,bpar,io)

implicit none

! Passed variables
class(cmat_type),intent(in) :: cmat !< C matrix data
type(nam_type),intent(in) :: nam    !< Namelist
type(geom_type),intent(in) :: geom  !< Geometry
type(bpar_type),intent(in) :: bpar  !< Block parameters
type(io_type),intent(in) :: io      !< I/O

! Local variables
integer :: ib
character(len=1024) :: filename

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      filename = trim(nam%prefix)//'_'//trim(cmat%blk(ib)%name)//'.nc'
      if (bpar%nicas_block(ib)) then
         call io%fld_write(nam,geom,filename,'coef_ens',cmat%blk(ib)%coef_ens)
         call io%fld_write(nam,geom,filename,'coef_sta',cmat%blk(ib)%coef_sta)
         call io%fld_write(nam,geom,filename,'rh0',cmat%blk(ib)%rh0)
         call io%fld_write(nam,geom,filename,'rv0',cmat%blk(ib)%rv0)
         call io%fld_write(nam,geom,filename,'rh0s',cmat%blk(ib)%rh0s)
         call io%fld_write(nam,geom,filename,'rv0s',cmat%blk(ib)%rv0s)
      end if
      if ((ib==bpar%nb+1).and.nam%displ_diag) then
         call io%fld_write(nam,geom,filename,'displ_lon',cmat%blk(ib)%displ_lon)
         call io%fld_write(nam,geom,filename,'displ_lat',cmat%blk(ib)%displ_lat)
      end if
   end if
end do

end subroutine cmat_write

!----------------------------------------------------------------------
! Subroutine: cmat_run_hdiag
!> Purpose: HDIAG driver
!----------------------------------------------------------------------
subroutine cmat_run_hdiag(cmat,nam,geom,bpar,io,ens1,ens2)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat     !< C matrix data
type(nam_type),intent(inout) :: nam        !< Namelist
type(geom_type),intent(in) :: geom         !< Geometry
type(bpar_type),intent(in) :: bpar         !< Block parameters
type(io_type),intent(in) :: io             !< I/O
type(ens_type),intent(in) :: ens1          !< Ensemble 1
type(ens_type),intent(in),optional :: ens2 !< Ensemble 2

! Local variables
integer :: ib
character(len=1024) :: filename
type(avg_type) :: avg_1,avg_2,avg_wgt
type(diag_type) :: cov_1,cov_2,cor_1,cor_2,loc_1,loc_2,loc_3
type(displ_type) :: displ
type(hdata_type) :: hdata
type(mom_type) :: mom_1,mom_2

! Setup sampling
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'
call flush(mpl%unit)
call hdata%setup_sampling(nam,geom,io)

! Compute MPI distribution, halo A
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A'
call flush(mpl%unit)
call hdata%compute_mpi_a(nam,geom)

if (nam%local_diag.or.nam%displ_diag) then
   ! Compute MPI distribution, halos A-B
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A-B'
   call flush(mpl%unit)
   call hdata%compute_mpi_ab(geom)
end if

if (nam%displ_diag) then
   ! Compute MPI distribution, halo D
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute MPI distribution, halo D'
   call flush(mpl%unit)
   call hdata%compute_mpi_d(nam,geom)

   ! Compute displacement diagnostic
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute displacement diagnostic'
   call flush(mpl%unit)
   call displ%compute(nam,geom,hdata,ens1)
end if

! Compute MPI distribution, halo C
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute MPI distribution, halo C'
call flush(mpl%unit)
call hdata%compute_mpi_c(nam,geom)

! Compute sample moments
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute sample moments'
call flush(mpl%unit)

! Compute ensemble 1 sample moments
write(mpl%unit,'(a7,a)') '','Ensemble 1:'
call flush(mpl%unit)
call mom_1%compute(nam,geom,bpar,hdata,ens1)

if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) then
   ! Compute randomized sample moments
   write(mpl%unit,'(a7,a)') '','Ensemble 2:'
   call flush(mpl%unit)
   call mom_2%compute(nam,geom,bpar,hdata,ens2)
end if

! Compute statistics
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute statistics'
call flush(mpl%unit)

! Compute ensemble 1 statistics
write(mpl%unit,'(a7,a)') '','Ensemble 1:'
call flush(mpl%unit)
call avg_1%compute(nam,geom,bpar,hdata,mom_1,nam%ne)

if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) then
   ! Compute randomized sample moments
   write(mpl%unit,'(a7,a)') '','Ensemble 2:'
   call flush(mpl%unit)
   call avg_2%compute(nam,geom,bpar,hdata,mom_2,nam%ens2_ne)
end if

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd','dual-ens')
   ! Compute hybrid statistics
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute hybrid statistics'
   call flush(mpl%unit)
   call avg_2%compute_hyb(nam,geom,bpar,hdata,mom_1,mom_2,avg_1)
end select

if (bpar%diag_block(bpar%nb+1)) then
   ! Compute block-averaged statistics
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute block-averaged statistics'
   call flush(mpl%unit)
   avg_wgt = avg_1%copy_wgt(geom,bpar)
   call avg_1%compute_bwavg(nam,geom,bpar,avg_wgt)
   if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) call avg_2%compute_bwavg(nam,geom,bpar,avg_wgt)
end if

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute covariance'
call flush(mpl%unit)

! Compute ensemble 1 covariance
write(mpl%unit,'(a7,a)') '','Ensemble 1:'
call flush(mpl%unit)
call cov_1%covariance(nam,geom,bpar,io,hdata,avg_1,'cov')

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd','dual-ens')
   ! Compute ensemble 2 covariance
   write(mpl%unit,'(a7,a)') '','Ensemble 2:'
   call flush(mpl%unit)
   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      call cov_2%covariance(nam,geom,bpar,io,hdata,avg_2,'cov_sta')
   case ('dual-ens')
      call cov_2%covariance(nam,geom,bpar,io,hdata,avg_2,'cov_lr')
   end select
end select

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute correlation'
call flush(mpl%unit)

! Compute ensemble 1 correlation
write(mpl%unit,'(a7,a)') '','Ensemble 1:'
call flush(mpl%unit)
call cor_1%correlation(nam,geom,bpar,io,hdata,avg_1,'cor')

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd','dual-ens')
   ! Compute ensemble 2 correlation
   write(mpl%unit,'(a7,a)') '','Ensemble 2:'
   call flush(mpl%unit)
   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      call cor_2%correlation(nam,geom,bpar,io,hdata,avg_2,'cor_sta')
   case ('dual-ens')
      call cor_2%correlation(nam,geom,bpar,io,hdata,avg_2,'cor_lr')
   end select
end select

select case (trim(nam%method))
case ('loc','hyb-avg','hyb-rnd','dual-ens')
   ! Compute localization
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute localization'
   write(mpl%unit,'(a7,a)') '','Ensemble 1:'
   call flush(mpl%unit)
   call loc_1%localization(nam,geom,bpar,io,hdata,avg_1,'loc')
end select

select case (trim(nam%method))
case ('hyb-avg','hyb-rnd')
   ! Compute static hybridization
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute static hybridization'
   write(mpl%unit,'(a7,a)') '','Ensemble 1 and 2:'
   call flush(mpl%unit)
   call loc_2%hybridization(nam,geom,bpar,io,hdata,avg_1,avg_2,'loc_hyb')
end select

if (trim(nam%method)=='dual-ens') then
   ! Compute dual-ensemble hybridization diagnostic and fit
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute dual-ensemble hybridization'
   write(mpl%unit,'(a7,a)') '','Ensembles 1 and 2:'
   call flush(mpl%unit)
   call loc_2%dualens(nam,geom,bpar,io,hdata,avg_1,avg_2,loc_3,'loc_deh','loc_deh_lr')
end if

if (trim(nam%minim_algo)/='none') then
   ! Copy diagnostics into C matrix data
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Copy diagnostics into C matrix data'
   call flush(mpl%unit)
   select case (trim(nam%method))
   case ('cor')
      call cmat%from(nam,geom,bpar,hdata,cor_1)
   case ('loc')
      call cmat%from(nam,geom,bpar,hdata,loc_1)
   case default
      call msgerror('cmat not implemented yet for this method')
   end select

   ! Write C matrix data
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write C matrix data'
   call flush(mpl%unit)
   call cmat%write(nam,geom,bpar,io)
end if

! Write data
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Write data'
call flush(mpl%unit)

! Displacement
if (nam%displ_diag) call displ%write(nam,geom,hdata,trim(nam%prefix)//'_displ_diag.nc')

! Full variances
if (nam%full_var) then
   filename = trim(nam%prefix)//'_full_var.nc'
   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) call io%fld_write(nam,geom,filename,trim(bpar%blockname(ib))//'_var', &
    & sum(mom_1%blk(ib)%m2full,dim=3)/real(mom_1%blk(ib)%nsub,kind_real))
   end do
end if

end subroutine cmat_run_hdiag

!----------------------------------------------------------------------
! Subroutine: cmat_from_fields
!> Purpose: copy radii into cmat object
!----------------------------------------------------------------------
subroutine cmat_from_fields(cmat,nam,geom,bpar,rh0,rv0)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat                               !< C matrix data
type(nam_type),intent(in) :: nam                                     !< Namelist
type(geom_type),intent(in) :: geom                                   !< Geometry
type(bpar_type),intent(in) :: bpar                                   !< Block parameters
real(kind_real),intent(in) :: rh0(geom%nmga,geom%nl0,nam%nv,nam%nts) !< Horizontal support radius on model grid, halo A
real(kind_real),intent(in) :: rv0(geom%nmga,geom%nl0,nam%nv,nam%nts) !< Vertical support radius on model grid, halo A

! Local variables
integer :: ib,iv,jv,its,jts,il0,ic0

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Copy radii into C matrix'
call flush(mpl%unit)

! Allocation
call cmat%alloc(nam,geom,bpar,'cmat')

! Convolution parameters
do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      if (bpar%nicas_block(ib)) then
         ! Indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         its = bpar%b_to_ts1(ib)
         jts = bpar%b_to_ts2(ib)
         if ((iv/=jv).or.(its/=jts)) call msgerror('only diagonal blocks for cmat_from_radii')

         ! Copy support radii
         do il0=1,geom%nl0
            cmat%blk(ib)%rh0(:,il0) = rh0(geom%c0a_to_mga,il0,iv,its)
            cmat%blk(ib)%rv0(:,il0) = rv0(geom%c0a_to_mga,il0,iv,its)
         end do

         ! Set coefficients
         cmat%blk(ib)%coef_ens = 1.0
         cmat%blk(ib)%coef_sta = 0.0
         cmat%blk(ib)%wgt = 1.0
      end if
   end if
end do

! Sampling parameters
if (trim(nam%strategy)=='specific_multivariate') then
   ! Initialization
   cmat%blk(ib)%rh0s = huge(1.0)
   cmat%blk(ib)%rv0s = huge(1.0)

   ! Get minimum
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         do il0=1,geom%nl0
            do ic0=1,geom%nc0
               cmat%blk(ib)%rh0s(ic0,il0) = min(cmat%blk(ib)%rh0s(ic0,il0),cmat%blk(ib)%rh0(ic0,il0))
               cmat%blk(ib)%rv0s(ic0,il0) = min(cmat%blk(ib)%rv0s(ic0,il0),cmat%blk(ib)%rv0(ic0,il0))
            end do
         end do
      end if
   end do
else
   ! Copy
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         cmat%blk(ib)%rh0s = cmat%blk(ib)%rh0
         cmat%blk(ib)%rv0s = cmat%blk(ib)%rv0
      end if
   end do
end if

end subroutine cmat_from_fields

!----------------------------------------------------------------------
! Subroutine: cmat_from_diag
!> Purpose: transform diagnostics into cmat object
!----------------------------------------------------------------------
subroutine cmat_from_diag(cmat,nam,geom,bpar,hdata,diag)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat !< C matrix data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(diag_type),intent(in) :: diag     !< Diagnostics

! Local variables
integer :: ib,i,ic0a,il0,il0i,ic2a,its
real(kind_real) :: fld_c2a(hdata%nc2a,geom%nl0),fld_c2b(hdata%nc2b,geom%nl0),fld_c0a(geom%nc0a,geom%nl0)

! Allocation
call cmat%alloc(nam,geom,bpar,'cmat')

! Convolution parameters
do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      if (bpar%nicas_block(ib)) then
         if (nam%local_diag) then
            do i=1,4
               ! Copy data
               do ic2a=1,hdata%nc2a
                  if (i==1) then
                     fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%raw_coef_ens
                  elseif (i==2) then
                     select case (trim(nam%method))
                     case ('cor','loc')
                        fld_c2a(ic2a,:) = 0.0
                     case ('hyb-avg','hyb-rnd')
                        fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%raw_coef_sta
                     case ('dual-ens')
                        call msgerror('dual-ens not ready yet for C matrix data')
                     end select
                  elseif (i==3) then
                     fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rh
                  elseif (i==4) then
                     fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rv
                  end if
               end do

               ! Median filter
               do il0=1,geom%nl0
                  call hdata%diag_filter(geom,il0,'median',nam%diag_rhflt,fld_c2a(:,il0))
               end do

               ! Interpolate
               call hdata%com_AB%ext(geom%nl0,fld_c2a,fld_c2b)
               do il0=1,geom%nl0
                  il0i = min(il0,geom%nl0i)
                  call hdata%h(il0i)%apply(fld_c2b(:,il0),fld_c0a(:,il0))
               end do

               if (i==1) then
                  cmat%blk(ib)%coef_ens = fld_c0a
                  call mpl%allreduce_sum(sum(cmat%blk(ib)%coef_ens,mask=geom%mask(geom%c0a_to_c0,:)),cmat%blk(ib)%wgt)
                  cmat%blk(ib)%wgt = cmat%blk(ib)%wgt/real(count(geom%mask),kind_real)
               elseif (i==2) then
                  cmat%blk(ib)%coef_sta = fld_c0a
               elseif (i==3) then
                  cmat%blk(ib)%rh0 = fld_c0a
               elseif (i==4) then
                  cmat%blk(ib)%rv0 = fld_c0a
               end if
            end do
         else
            do il0=1,geom%nl0
               cmat%blk(ib)%coef_ens(:,il0) = diag%blk(0,ib)%raw_coef_ens(il0)
               cmat%blk(ib)%rh0(:,il0) = diag%blk(0,ib)%fit_rh(il0)
               cmat%blk(ib)%rv0(:,il0) = diag%blk(0,ib)%fit_rv(il0)
               select case (trim(nam%method))
               case ('cor','loc')
                  cmat%blk(ib)%coef_sta(:,il0) = 0.0
               case ('hyb-avg','hyb-rnd')
                  cmat%blk(ib)%coef_sta(:,il0) = diag%blk(0,ib)%raw_coef_sta
               case ('dual-ens')
                  call msgerror('dual-ens not ready yet for C matrix data')
               end select
            end do
            cmat%blk(ib)%wgt = sum(diag%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
         end if
      else
         cmat%blk(ib)%wgt = sum(diag%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
      end if
   end if
end do

! Sampling parameters
if (trim(nam%strategy)=='specific_multivariate') then
   ! Initialization
   cmat%blk(ib)%rh0s = huge(1.0)
   cmat%blk(ib)%rv0s = huge(1.0)

   ! Get minimum
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               cmat%blk(ib)%rh0s(ic0a,il0) = min(cmat%blk(ib)%rh0s(ic0a,il0),cmat%blk(ib)%rh0(ic0a,il0))
               cmat%blk(ib)%rv0s(ic0a,il0) = min(cmat%blk(ib)%rv0s(ic0a,il0),cmat%blk(ib)%rv0(ic0a,il0))
            end do
         end do
      end if
   end do
else
   ! Copy
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         cmat%blk(ib)%rh0s = cmat%blk(ib)%rh0
         cmat%blk(ib)%rv0s = cmat%blk(ib)%rv0
      end if
   end do
end if

! Displacement
if (nam%displ_diag) then
   do its=2,nam%nts
      cmat%blk(bpar%nb+1)%displ_lon(:,:,its) = hdata%displ_lon(:,:,its)
      cmat%blk(bpar%nb+1)%displ_lat(:,:,its) = hdata%displ_lat(:,:,its)
   end do
end if

end subroutine cmat_from_diag

end module type_cmat
