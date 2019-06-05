!----------------------------------------------------------------------
! Module: type_samp
! Purpose: sampling derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_samp

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg
use tools_func, only: gc99,sphere_dist
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: eq,inf
use type_bpar, only: bpar_type
use type_com, only: com_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_tree, only: tree_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Sampling derived type
type samp_type
   ! Sampling
   real(kind_real),allocatable :: rh_c0(:)          ! Sampling radius on subset Sc0
   logical,allocatable :: mask_c0(:,:)              ! Mask on subset Sc0
   logical,allocatable :: mask_hor_c0(:)            ! Union of horizontal masks on subset Sc0, global
   integer,allocatable :: nc0_mask(:)               ! Horizontal mask size on subset Sc0
   integer,allocatable :: c1_to_c0(:)               ! First sampling index
   logical,allocatable :: c1l0_log(:,:)             ! Log for the first sampling index
   integer,allocatable :: c1c3_to_c0(:,:)           ! Second horizontal sampling index
   logical,allocatable :: c1c3l0_log(:,:,:)         ! Log for the second horizontal sampling index
   integer,allocatable :: c2_to_c1(:)               ! Subgrid to diagnostic points
   integer,allocatable :: c2_to_c0(:)               ! Subgrid to grid
   logical,allocatable :: mask_c2(:,:)              ! Mask on subset Sc2

   ! Local data
   logical,allocatable ::  vbal_mask(:,:)           ! Vertical balance mask
   logical,allocatable ::  local_mask(:,:)          ! Local mask
   integer,allocatable :: nn_c2a_index(:,:)         ! Nearest diagnostic neighbors from diagnostic points
   real(kind_real),allocatable :: nn_c2a_dist(:,:)  ! Nearest diagnostic neighbors distance from diagnostic points

   ! Forced points
   integer :: nfor                                  ! Number of forced points
   integer,allocatable :: ldwv_to_c0(:)             ! Local diagnostics profiles to subset Sc0
   integer,allocatable :: ldwv_to_c2(:)             ! Local diagnostics profiles to subset Sc2

   ! Sampling mesh
   type(mesh_type) :: mesh                          ! Sampling mesh

   ! Advection
   integer,allocatable :: adv_nn(:)                 ! Number of nearest neighbors inside search radius
   integer :: adv_nnmax                             ! Maximum number of nearest neighbors inside search radius
   integer,allocatable :: adv_nn_index(:,:)         ! Index of nearest neighbors inside search radius
   real(kind_real),allocatable :: adv_lon(:,:,:)    ! Interpolated advected longitude
   real(kind_real),allocatable :: adv_lat(:,:,:)    ! Interpolated advected latitude

   ! Interpolations
   type(linop_type),allocatable :: hfull(:)         ! Horizontal interpolation from Sc2 to Sc0 (global)
   type(linop_type),allocatable :: h(:)             ! Horizontal interpolation from Sc2 to Sc0 (local)
   type(linop_type),allocatable :: d(:,:)           ! Advection interpolation

   ! MPI distribution
   integer :: nc0c                                  ! Number of points in subset Sc0, halo C
   integer :: nc1a                                  ! Number of points in subset Sc1, halo A
   integer :: nc2a                                  ! Number of points in subset Sc2, halo A
   integer :: nc2b                                  ! Number of points in subset Sc2, halo B
   integer :: nc2f                                  ! Number of points in subset Sc2, halo F
   logical,allocatable :: lcheck_c0a(:)             ! Detection of halo A on subset Sc0
   logical,allocatable :: lcheck_c0c(:)             ! Detection of halo C on subset Sc0
   logical,allocatable :: lcheck_c1a(:)             ! Detection of halo A on subset Sc1
   logical,allocatable :: lcheck_c2a(:)             ! Detection of halo A on subset Sc2
   logical,allocatable :: lcheck_c2b(:)             ! Detection of halo B on subset Sc2
   logical,allocatable :: lcheck_c2f(:)             ! Detection of halo F on subset Sc2
   logical,allocatable :: lcheck_h(:,:)             ! Detection of horizontal interpolation coefficients
   logical,allocatable :: lcheck_d(:,:,:)           ! Detection of advection interpolation coefficients
   integer,allocatable :: c0c_to_c0(:)              ! Subset Sc0, halo C to global
   integer,allocatable :: c0_to_c0c(:)              ! Subset Sc0, global to halo C
   integer,allocatable :: c0a_to_c0c(:)             ! Subset Sc0, halo A to halo C
   integer,allocatable :: c1a_to_c1(:)              ! Subset Sc1, halo A to global
   integer,allocatable :: c1_to_c1a(:)              ! Subset Sc1, global to halo A
   integer,allocatable :: c2a_to_c2(:)              ! Subset Sc2, halo A to global
   integer,allocatable :: c2_to_c2a(:)              ! Subset Sc2, global to halo A
   logical,allocatable :: mask_c2a(:,:)             ! Mask on subset Sc2, halo A
   integer,allocatable :: c2b_to_c2(:)              ! Subset Sc2, halo B to global
   integer,allocatable :: c2_to_c2b(:)              ! Subset Sc2, global to halo B
   integer,allocatable :: c2a_to_c2b(:)             ! Subset Sc2, halo A to halo B
   integer,allocatable :: c2f_to_c2(:)              ! Subset Sc2, halo B to global
   integer,allocatable :: c2_to_c2f(:)              ! Subset Sc2, global to halo B
   integer,allocatable :: c2a_to_c2f(:)             ! Subset Sc2, halo A to halo B
   integer,allocatable :: c2_to_proc(:)             ! Subset Sc2, global to processor
   integer,allocatable :: proc_to_nc2a(:)           ! Number of points in subset Sc2, halo A, for each processor
   type(com_type) :: com_AC                         ! Communication between halos A and C
   type(com_type) :: com_AB                         ! Communication between halos A and B
   type(com_type) :: com_AF                         ! Communication between halos A and F (filtering)
contains
   procedure :: samp_alloc_mask
   procedure :: samp_alloc_other
   generic :: alloc => samp_alloc_mask,samp_alloc_other
   procedure :: dealloc => samp_dealloc
   procedure :: read => samp_read
   procedure :: write => samp_write
   procedure :: setup_sampling => samp_setup_sampling
   procedure :: compute_mask => samp_compute_mask
   procedure :: compute_sampling_zs => samp_compute_sampling_zs
   procedure :: compute_sampling_ps => samp_compute_sampling_ps
   procedure :: compute_sampling_lct => samp_compute_sampling_lct
   procedure :: check_mask => samp_check_mask
   procedure :: compute_mpi_a => samp_compute_mpi_a
   procedure :: compute_mpi_ab => samp_compute_mpi_ab
   procedure :: compute_mpi_c => samp_compute_mpi_c
   procedure :: compute_mpi_f => samp_compute_mpi_f
   procedure :: diag_filter => samp_diag_filter
   procedure :: diag_fill => samp_diag_fill
end type samp_type

private
public :: samp_type

contains

!----------------------------------------------------------------------
! Subroutine: samp_alloc_mask
! Purpose: allocation for mask
!----------------------------------------------------------------------
subroutine samp_alloc_mask(samp,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(geom_type),intent(in) :: geom     ! Geometry

! Allocation
allocate(samp%mask_c0(geom%nc0,geom%nl0))
allocate(samp%mask_hor_c0(geom%nc0))
allocate(samp%nc0_mask(geom%nl0))

end subroutine samp_alloc_mask

!----------------------------------------------------------------------
! Subroutine: samp_alloc_other
! Purpose: allocation for other variables
!----------------------------------------------------------------------
subroutine samp_alloc_other(samp,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Allocation
allocate(samp%rh_c0(geom%nc0))
allocate(samp%c1_to_c0(nam%nc1))
allocate(samp%c1l0_log(nam%nc1,geom%nl0))
allocate(samp%c1c3_to_c0(nam%nc1,nam%nc3))
allocate(samp%c1c3l0_log(nam%nc1,nam%nc3,geom%nl0))
if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%local_diag.or.nam%adv_diag))) then
   allocate(samp%c2_to_c1(nam%nc2))
   allocate(samp%c2_to_c0(nam%nc2))
   allocate(samp%mask_c2(nam%nc2,geom%nl0))
   allocate(samp%hfull(geom%nl0i))
   call samp%mesh%alloc(nam%nc2)
end if
if (nam%adv_diag) then
   allocate(samp%adv_lon(geom%nc0a,geom%nl0,nam%nts))
   allocate(samp%adv_lat(geom%nc0a,geom%nl0,nam%nts))
end if

end subroutine samp_alloc_other

!----------------------------------------------------------------------
! Subroutine: samp_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine samp_dealloc(samp)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling

! Local variables
integer :: il0

! Release memory
if (allocated(samp%rh_c0)) deallocate(samp%rh_c0)
if (allocated(samp%mask_c0)) deallocate(samp%mask_c0)
if (allocated(samp%mask_hor_c0)) deallocate(samp%mask_hor_c0)
if (allocated(samp%nc0_mask)) deallocate(samp%nc0_mask)
if (allocated(samp%c1_to_c0)) deallocate(samp%c1_to_c0)
if (allocated(samp%c1l0_log)) deallocate(samp%c1l0_log)
if (allocated(samp%c1c3_to_c0)) deallocate(samp%c1c3_to_c0)
if (allocated(samp%c1c3l0_log)) deallocate(samp%c1c3l0_log)
if (allocated(samp%c2_to_c1)) deallocate(samp%c2_to_c1)
if (allocated(samp%c2_to_c0)) deallocate(samp%c2_to_c0)
if (allocated(samp%mask_c2)) deallocate(samp%mask_c2)
if (allocated(samp%c2a_to_c2)) deallocate(samp%c2a_to_c2)
if (allocated(samp%c2_to_c2a)) deallocate(samp%c2_to_c2a)
if (allocated(samp%mask_c2a)) deallocate(samp%mask_c2a)
if (allocated(samp%vbal_mask)) deallocate(samp%vbal_mask)
if (allocated(samp%local_mask)) deallocate(samp%local_mask)
if (allocated(samp%nn_c2a_index)) deallocate(samp%nn_c2a_index)
if (allocated(samp%nn_c2a_dist)) deallocate(samp%nn_c2a_dist)
if (allocated(samp%adv_nn)) deallocate(samp%adv_nn)
if (allocated(samp%adv_nn_index)) deallocate(samp%adv_nn_index)
if (allocated(samp%hfull)) then
   do il0=1,size(samp%hfull)
      call samp%hfull(il0)%dealloc
   end do
   deallocate(samp%hfull)
end if
if (allocated(samp%h)) then
   do il0=1,size(samp%h)
      call samp%h(il0)%dealloc
   end do
   deallocate(samp%h)
end if
if (allocated(samp%ldwv_to_c0)) deallocate(samp%ldwv_to_c0)
if (allocated(samp%ldwv_to_c2)) deallocate(samp%ldwv_to_c2)
if (allocated(samp%adv_lon)) deallocate(samp%adv_lon)
if (allocated(samp%adv_lat)) deallocate(samp%adv_lat)

end subroutine samp_dealloc

!----------------------------------------------------------------------
! Subroutine: samp_read
! Purpose: read
!----------------------------------------------------------------------
subroutine samp_read(samp,mpl,nam,geom,bpar,new_sampling)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
logical,intent(out) :: new_sampling    ! Status flag

! Local variables
integer :: il0,il0i,ic1,jc3
integer :: info,ncid,nl0_id,nl0r_id,nc3_id,nc1_id,nc2_id
integer :: c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id
integer :: c1l0_logint(nam%nc1,geom%nl0),c1c3l0_logint(nam%nc1,nam%nc3,geom%nl0)
character(len=3) :: il0ichar
character(len=1024),parameter :: subr = 'samp_read'

! Initialization
new_sampling = .false.

! Open file
info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',nf90_nowrite,ncid)
if (info/=nf90_noerr) then
   call mpl%warning(subr,'cannot find sampling to read, recomputing sampling')
   new_sampling = .true.
   return
end if

! Check dimensions
nl0_id = mpl%ncdimcheck(subr,ncid,'nl0',geom%nl0,.false.)
nl0r_id = mpl%ncdimcheck(subr,ncid,'nl0r',bpar%nl0rmax,.false.)
nc3_id = mpl%ncdimcheck(subr,ncid,'nc3',nam%nc3,.false.)
nc1_id = mpl%ncdimcheck(subr,ncid,'nc1',nam%nc1,.false.)
if (nam%new_lct.or.nam%local_diag.or.nam%adv_diag) then
   nc2_id = mpl%ncdimcheck(subr,ncid,'nc2',nam%nc2,.false.)
else
   nc2_id = 0
end if
if (mpl%msv%isi(nl0_id).or.mpl%msv%isi(nl0r_id).or.mpl%msv%isi(nc3_id).or.mpl%msv%isi(nc1_id) &
 & .or.mpl%msv%isi(nc2_id)) then
   call mpl%warning(subr,'wrong dimension when reading sampling, recomputing sampling')
   call mpl%ncerr(subr,nf90_close(ncid))
   new_sampling = .true.
   return
end if
write(mpl%info,'(a7,a)') '','Read sampling'
call mpl%flush

! Get arrays ID
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1_to_c0',c1_to_c0_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1l0_log',c1l0_log_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1c3_to_c0',c1c3_to_c0_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1c3l0_log',c1c3l0_log_id))
if (nam%new_lct.or.nam%local_diag.or.nam%adv_diag) then
   call mpl%ncerr(subr,nf90_inq_varid(ncid,'c2_to_c1',c2_to_c1_id))
   call mpl%ncerr(subr,nf90_inq_varid(ncid,'c2_to_c0',c2_to_c0_id))
end if

! Read arrays
call mpl%ncerr(subr,nf90_get_var(ncid,c1_to_c0_id,samp%c1_to_c0))
call mpl%ncerr(subr,nf90_get_var(ncid,c1l0_log_id,c1l0_logint))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (c1l0_logint(ic1,il0)==0) then
         samp%c1l0_log(ic1,il0) = .false.
      else if (c1l0_logint(ic1,il0)==1) then
         samp%c1l0_log(ic1,il0) = .true.
      else
         call mpl%abort(subr,'wrong c1l0_log')
      end if
   end do
end do
call mpl%ncerr(subr,nf90_get_var(ncid,c1c3_to_c0_id,samp%c1c3_to_c0))
call mpl%ncerr(subr,nf90_get_var(ncid,c1c3l0_log_id,c1c3l0_logint))
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (c1c3l0_logint(ic1,jc3,il0)==0) then
            samp%c1c3l0_log(ic1,jc3,il0) = .false.
         else if (c1c3l0_logint(ic1,jc3,il0)==1) then
            samp%c1c3l0_log(ic1,jc3,il0) = .true.
         else
            call mpl%abort(subr,'wrong c1c3l0_log')
         end if
      end do
   end do
end do
if (nam%new_lct.or.nam%local_diag.or.nam%adv_diag) then
   call mpl%ncerr(subr,nf90_get_var(ncid,c2_to_c1_id,samp%c2_to_c1))
   call mpl%ncerr(subr,nf90_get_var(ncid,c2_to_c0_id,samp%c2_to_c0))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%local_diag.or.nam%adv_diag))) then
   ! Read interpolation
   do il0i=1,geom%nl0i
      ! Open file
      write(il0ichar,'(i3.3)') il0i
      info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc',nf90_nowrite,ncid)
      if (info/=nf90_noerr) then
         call mpl%warning(subr,'cannot find interpolation data to read, recomputing sampling')
         new_sampling = .true.
         return
      end if

      ! Interpolation
      write(samp%hfull(il0i)%prefix,'(a,i3.3)') 'hfull_',il0i
      call samp%hfull(il0i)%read(mpl,ncid)

      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   end do
end if

end subroutine samp_read

!----------------------------------------------------------------------
! Subroutine: samp_write
! Purpose: write
!----------------------------------------------------------------------
subroutine samp_write(samp,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(bpar_type),intent(in) :: bpar  ! Block parameters

! Local variables
integer :: il0,il0i,ic1,jc3
integer :: ncid,nl0_id,nl0r_id,nc1_id,nc2_id,nc3_id
integer :: lat_id,lon_id,smax_id,c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id
integer :: c1l0_logint(nam%nc1,geom%nl0),c1c3l0_logint(nam%nc1,nam%nc3,geom%nl0)
real(kind_real) :: lon(nam%nc1,nam%nc3,geom%nl0),lat(nam%nc1,nam%nc3,geom%nl0)
character(len=3) :: il0ichar
character(len=1024),parameter :: subr = 'samp_write'

! Processor verification
if (.not.mpl%main) call mpl%abort(subr,'only I/O proc should enter '//trim(subr))

! Create file
write(mpl%info,'(a7,a)') '','Write sampling'
call mpl%flush
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%write(mpl,ncid)

! Define dimensions
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0r',bpar%nl0rmax,nl0r_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc3',nam%nc3,nc3_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%local_diag.or.nam%adv_diag))) &
 & call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2',nam%nc2,nc2_id))

! Define variables
call mpl%ncerr(subr,nf90_def_var(ncid,'lat',nc_kind_real,(/nc1_id,nc3_id,nl0_id/),lat_id))
call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',mpl%msv%valr))
call mpl%ncerr(subr,nf90_def_var(ncid,'lon',nc_kind_real,(/nc1_id,nc3_id,nl0_id/),lon_id))
call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',mpl%msv%valr))
call mpl%ncerr(subr,nf90_def_var(ncid,'smax',nc_kind_real,(/nc3_id,nl0_id/),smax_id))
call mpl%ncerr(subr,nf90_put_att(ncid,smax_id,'_FillValue',mpl%msv%valr))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1_to_c0',nf90_int,(/nc1_id/),c1_to_c0_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1_to_c0_id,'_FillValue',mpl%msv%vali))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1l0_log',nf90_int,(/nc1_id,nl0_id/),c1l0_log_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1l0_log_id,'_FillValue',mpl%msv%vali))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1c3_to_c0',nf90_int,(/nc1_id,nc3_id/),c1c3_to_c0_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1c3_to_c0_id,'_FillValue',mpl%msv%vali))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1c3l0_log',nf90_int,(/nc1_id,nc3_id,nl0_id/),c1c3l0_log_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1c3l0_log_id,'_FillValue',mpl%msv%vali))
if (nam%new_lct.or.nam%local_diag.or.nam%adv_diag) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'c2_to_c1',nf90_int,(/nc2_id/),c2_to_c1_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,c2_to_c1_id,'_FillValue',mpl%msv%vali))
   call mpl%ncerr(subr,nf90_def_var(ncid,'c2_to_c0',nf90_int,(/nc2_id/),c2_to_c0_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,c2_to_c0_id,'_FillValue',mpl%msv%vali))
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Convert data
lon = mpl%msv%valr
lat = mpl%msv%valr
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (samp%c1c3l0_log(ic1,jc3,il0)) then
            lon(ic1,jc3,il0) = geom%lon(samp%c1c3_to_c0(ic1,jc3))*rad2deg
            lat(ic1,jc3,il0) = geom%lat(samp%c1c3_to_c0(ic1,jc3))*rad2deg
            c1c3l0_logint(ic1,jc3,il0) = 1
         else
            c1c3l0_logint(ic1,jc3,il0) = 0
         end if
      end do
   end do
   do ic1=1,nam%nc1
      if (samp%c1l0_log(ic1,il0)) then
         c1l0_logint(ic1,il0) = 1
      else
         c1l0_logint(ic1,il0) = 0
      end if
   end do
end do

! Write variables
call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_put_var(ncid,smax_id,real(count(samp%c1c3l0_log,dim=1),kind_real)))
call mpl%ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,samp%c1_to_c0))
call mpl%ncerr(subr,nf90_put_var(ncid,c1l0_log_id,c1l0_logint))
call mpl%ncerr(subr,nf90_put_var(ncid,c1c3_to_c0_id,samp%c1c3_to_c0))
call mpl%ncerr(subr,nf90_put_var(ncid,c1c3l0_log_id,c1c3l0_logint))
if (nam%new_lct.or.nam%local_diag.or.nam%adv_diag) then
   call mpl%ncerr(subr,nf90_put_var(ncid,c2_to_c1_id,samp%c2_to_c1))
   call mpl%ncerr(subr,nf90_put_var(ncid,c2_to_c0_id,samp%c2_to_c0))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%local_diag.or.nam%adv_diag))) then
   ! Write interpolation
   do il0i=1,geom%nl0i
      ! Create file
      write(il0ichar,'(i3.3)') il0i
      call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc', &
    & or(nf90_clobber,nf90_64bit_offset),ncid))

      ! Write namelist parameters
      call nam%write(mpl,ncid)

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      call samp%hfull(il0i)%write(mpl,ncid)

      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   end do
end if

end subroutine samp_write

!----------------------------------------------------------------------
! Subroutine: samp_setup_sampling
! Purpose: setup sampling
!----------------------------------------------------------------------
subroutine samp_setup_sampling(samp,mpl,rng,nam,geom,bpar,io,ens)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O
type(ens_type),intent(in) :: ens       ! Ensemble

! Local variables
integer :: ic0,il0,ic1,ic2,ildw,jc3,il0i,ifor,ival,nn_index(geom%nc0)
integer,allocatable :: vbot(:),vtop(:),for(:)
real(kind_real) :: rh_c0a(geom%nc0a,geom%nl0)
real(kind_real),allocatable :: lon_c1(:),lat_c1(:)
real(kind_real),allocatable :: lon_c2(:),lat_c2(:)
real(kind_real),allocatable :: rh_c1(:)
logical :: new_sampling
logical,allocatable :: mask_c1(:)
character(len=8) :: ivalformat
character(len=1024) :: filename,color
character(len=1024),parameter :: subr = 'samp_setup_sampling'
type(linop_type) :: hbase

! Allocation
call samp%alloc(geom)

! Compute nearest neighbors for local diagnostics output
if (nam%nldwv>0) then
   write(mpl%info,'(a7,a)') '','Compute local diagnostics locations:'
   call mpl%flush

   ! Allocation
   allocate(samp%ldwv_to_c0(nam%nldwv))

   ! Find nearest neighbors
   samp%ldwv_to_c0 = mpl%msv%vali
   do ildw=1,nam%nldwv
      ic0 = 0
      do while (mpl%msv%isi(samp%ldwv_to_c0(ildw)))
         ic0 = ic0+1
         call geom%tree%find_nearest_neighbors(nam%lon_ldwv(ildw),nam%lat_ldwv(ildw),ic0,nn_index(1:ic0))
         if (geom%mask_hor_c0(nn_index(ic0))) samp%ldwv_to_c0(ildw) = nn_index(ic0)
      end do
      write(mpl%info,'(a10,a,e15.8,a,e15.8)') '','Profile '//trim(nam%name_ldwv(ildw))//' computed at lon/lat: ', &
    & geom%lon(samp%ldwv_to_c0(ildw))*rad2deg,' / ',geom%lat(samp%ldwv_to_c0(ildw))*rad2deg
      call mpl%flush
   end do
end if

! Compute sampling mask
call samp%compute_mask(mpl,nam,geom,ens)

! Check subsampling size
if (nam%nc1>count(samp%mask_hor_c0)) then
   ! Not enough points remaining in the sampling mask
   call mpl%warning(subr,'not enough points remaining in sampling mask, resetting nc1 to the largest possible value')
   nam%nc1 = count(samp%mask_hor_c0)
end if
if (nam%nc2>nam%nc1) then
   ! Subsampling should have less points
   call mpl%warning(subr,'subsampling should have less points, resetting nc2 to nc1')
   nam%nc2 = nam%nc1
end if

! Allocation
call samp%alloc(nam,geom)
allocate(lon_c1(nam%nc1))
allocate(lat_c1(nam%nc1))
allocate(mask_c1(nam%nc1))

! Initialization
samp%c1_to_c0 = mpl%msv%vali
samp%c1c3_to_c0 = mpl%msv%vali

! Read or compute sampling data
new_sampling = .true.
if (nam%sam_read) then
   call samp%read(mpl,nam,geom,bpar,new_sampling)
   if (new_sampling) nam%sam_write = .true.
end if
if (new_sampling) then
   ! Compute zero-separation sampling
   call samp%compute_sampling_zs(mpl,rng,nam,geom)

   if (nam%new_lct) then
      ! Compute LCT sampling
      call samp%compute_sampling_lct(mpl,nam,geom)
   elseif (nam%new_vbal.or.nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
      ! Compute positive separation sampling
      call samp%compute_sampling_ps(mpl,rng,nam,geom)
   end if

   ! Check sampling mask
   call samp%check_mask(mpl,nam,geom)
end if

if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%local_diag.or.nam%adv_diag))) then
   if (new_sampling) then
      ! Allocation
      allocate(rh_c1(nam%nc1))
      allocate(for(samp%nfor))

      ! Define subsampling
      write(mpl%info,'(a7,a)') '','Define subsampling:'
      call mpl%flush(.false.)
      lon_c1 = geom%lon(samp%c1_to_c0)
      lat_c1 = geom%lat(samp%c1_to_c0)
      mask_c1 = any(samp%c1l0_log(:,:),dim=2)
      do ifor=1,samp%nfor
         for(ifor) = ifor
      end do
      rh_c1 = 1.0
      call rng%initialize_sampling(mpl,nam%nc1,lon_c1,lat_c1,mask_c1,samp%nfor,for,rh_c1, &
    & nam%ntry,nam%nrep,nam%nc2,samp%c2_to_c1)
      samp%c2_to_c0 = samp%c1_to_c0(samp%c2_to_c1)

      ! Release memory
      deallocate(rh_c1)
      deallocate(for)
   end if

   ! Define mask on subset Sc2
   do il0=1,geom%nl0
      do ic2=1,nam%nc2
         ic1 = samp%c2_to_c1(ic2)
         samp%mask_c2(ic2,il0) = samp%c1l0_log(ic1,il0)
      end do
   end do

   ! Allocation
   allocate(lon_c2(nam%nc2))
   allocate(lat_c2(nam%nc2))

   ! Initialization
   lon_c2 = geom%lon(samp%c2_to_c0)
   lat_c2 = geom%lat(samp%c2_to_c0)
   call samp%mesh%init(mpl,rng,lon_c2,lat_c2)

   ! Compute triangles list
   write(mpl%info,'(a7,a)') '','Compute triangles list '
   call mpl%flush
   call samp%mesh%trlist(mpl)

   ! Find boundary nodes
   write(mpl%info,'(a7,a)') '','Find boundary nodes'
   call mpl%flush
   call samp%mesh%bnodes(mpl)

   if (new_sampling) then
      ! Allocation
      allocate(vbot(nam%nc2))
      allocate(vtop(nam%nc2))

      ! Initialize vbot and vtop
      vbot = 1
      vtop = geom%nl0

      ! Compute grid interpolation
      write(mpl%info,'(a7,a)') '','Compute grid interpolation'
      call mpl%flush
      do il0i=1,geom%nl0i
         ! Compute grid interpolation
         write(samp%hfull(il0i)%prefix,'(a,i3.3)') 'hfull_',il0i
         call samp%hfull(il0i)%interp(mpl,rng,geom,il0i,nam%nc2,samp%c2_to_c0,nam%mask_check,vbot,vtop,nam%diag_interp,hbase)
      end do

      ! Release memory
      deallocate(vbot)
      deallocate(vtop)
   end if

   ! Release memory
   deallocate(lon_c2)
   deallocate(lat_c2)
end if

! Write sampling data
if (nam%sam_write) then
   if (mpl%main) call samp%write(mpl,nam,geom,bpar)

   ! Write rh_c0
   if (trim(nam%draw_type)=='random_coast') then
      rh_c0a = mpl%msv%valr
      call mpl%glb_to_loc(geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,samp%rh_c0,geom%nc0a,rh_c0a(:,1))
      filename = trim(nam%prefix)//'_sampling_rh_c0'
      call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)
      call io%fld_write(mpl,nam,geom,filename,'rh_c0',rh_c0a)
   end if
end if

! Print results
write(mpl%info,'(a7,a)') '','Sampling efficiency (%):'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a3)') '','Level ',nam%levs(il0),' ~>'
   call mpl%flush(.false.)
   do jc3=1,nam%nc3
      ival = int(100.0*real(count(samp%c1c3l0_log(:,jc3,il0)),kind_real)/real(nam%nc1,kind_real))
      if (ival==100) then
         ivalformat = '(a,i3,a)'
      else
         ivalformat = '(a,i2,a)'
      end if
      if (count(samp%c1c3l0_log(:,jc3,il0))>=nam%nc1/2) then
         ! Sucessful sampling
         color = mpl%green
      else
         ! Insufficient sampling
         color = mpl%peach
      end if
      if (jc3==1) color = ' '//trim(color)
      write(mpl%info,ivalformat) trim(color),ival,trim(mpl%black)
      call mpl%flush(.false.)
      if (jc3<nam%nc3) then
         write(mpl%info,'(a)') '-'
         call mpl%flush(.false.)
      end if
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
end do

end subroutine samp_setup_sampling

!----------------------------------------------------------------------
! Subroutine: samp_compute_mask
! Purpose: compute mask
!----------------------------------------------------------------------
subroutine samp_compute_mask(samp,mpl,nam,geom,ens)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(ens_type),intent(in) :: ens       ! Ensemble

! Local variables
integer :: ic0a,ic0,il0,ildw,iv,its,ie,ncontig,ncontigmax
integer :: latmin,latmax
real(kind_real) :: dist
real(kind_real),allocatable :: var(:,:,:,:)
logical :: valid,mask_c0a(geom%nc0a,geom%nl0)
character(len=1024),parameter :: subr = 'samp_compute_mask'

! Compute sampling mask
write(mpl%info,'(a7,a)') '','Compute sampling mask'
call mpl%flush

! Copy geometry mask
mask_c0a = geom%mask_c0a
if (allocated(geom%smask_c0a)) mask_c0a = mask_c0a.and.geom%smask_c0a

! Mask restriction
if (nam%mask_type(1:3)=='lat') then
   ! Latitude band
   read(nam%mask_type(4:6),'(i3)') latmin
   read(nam%mask_type(7:9),'(i3)') latmax
   write(mpl%info,'(a10,a,i3,a,i3)') '','Latitude band between ',latmin,' and ',latmax
   call mpl%flush
   if (latmin>=latmax) call mpl%abort(subr,'latmin should be lower than latmax')
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      valid = (geom%lat(ic0)>=real(latmin,kind_real)*deg2rad).and.(geom%lat(ic0)<=real(latmax,kind_real)*deg2rad)
      do il0=1,geom%nl0
         mask_c0a(ic0a,il0) = mask_c0a(ic0a,il0).and.valid
      end do
   end do
elseif (trim(nam%mask_type)=='ldwv') then
   ! Disk around vertical diagnostic points
   write(mpl%info,'(a10,a,e10.3,a)') '','Disk of ',1.1*nam%local_rad*reqkm,' km around vertical diagonstic points'
   call mpl%flush
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      valid = .false.
      do ildw=1,nam%nldwv
         call sphere_dist(geom%lon(samp%ldwv_to_c0(ildw)),geom%lat(samp%ldwv_to_c0(ildw)),geom%lon(ic0),geom%lat(ic0),dist)
         valid = valid.or.(dist<1.1*nam%local_rad)
      end do
      do il0=1,geom%nl0
         mask_c0a(ic0a,il0) = mask_c0a(ic0a,il0).and.valid
      end do
   end do
elseif (trim(nam%mask_type)=='stddev') then
   ! Standard-deviation threshold
   write(mpl%info,'(a10,a,e10.3,a)') '','Threshold ',nam%mask_th,' used as a '//trim(nam%mask_lu)//' bound for standard-deviation'
   call mpl%flush

   ! Allocation
   allocate(var(geom%nc0a,geom%nl0,nam%nv,nam%nts))

   ! Compute variances
   var = 0.0
   do ie=1,ens%ne
      var = var+ens%fld(:,:,:,:,ie)**2
   end do
   var = var/real(ens%ne-ens%nsub,kind_real)

   ! Check standard-deviation value
   do its=1,nam%nts
      do iv=1,nam%nv
         if (trim(nam%mask_lu)=='lower') then
            mask_c0a = mask_c0a.and.(var(:,:,iv,its)>nam%mask_th**2)
         elseif (trim(nam%mask_lu)=='upper') then
            mask_c0a = mask_c0a.and.(var(:,:,iv,its)<nam%mask_th**2)
         end if
      end do
   end do

   ! Release memory
   deallocate(var)
elseif (trim(nam%mask_type)=='none') then
   ! Nothing to do
else
   if (.not.allocated(geom%smask_c0a)) call mpl%abort(subr,'mask_type not recognized')
end if

! Check vertically contiguous points
if (nam%ncontig_th>0) then
   write(mpl%info,'(a10,a,i3,a)') '','Mask restricted with at least ',min(nam%ncontig_th,geom%nl0),' vertically contiguous points'
   call mpl%flush
   do ic0a=1,geom%nc0a
      ncontig = 0
      ncontigmax = 0
      do il0=1,geom%nl0
         if (mask_c0a(ic0a,il0)) then
            ncontig = ncontig+1
         else
            ncontig = 0
         end if
         if (ncontig>ncontigmax) ncontigmax = ncontig
      end do
      mask_c0a(ic0a,:) = mask_c0a(ic0a,:).and.(ncontigmax>=min(nam%ncontig_th,geom%nl0))
   end do
end if

! Local to global
call mpl%loc_to_glb(geom%nl0,geom%nc0a,mask_c0a,geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,samp%mask_c0)

! Check mask size
if (count(samp%mask_c0)==0) call mpl%abort(subr,'no more points in the sampling mask')

! Other masks
samp%mask_hor_c0 = any(samp%mask_c0,dim=2)
samp%nc0_mask = count(samp%mask_c0,dim=1)

! Print results
write(mpl%info,'(a7,a)') '','HDIAG sampling valid points (% of domain mask):'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a,f5.1,a)') '','Level',nam%levs(il0),' ~> ',100.0*real(count(samp%mask_c0(:,il0)),kind_real) &
 & /real(count(geom%mask_c0(:,il0)),kind_real),'%'
   call mpl%flush
end do

end subroutine samp_compute_mask

!----------------------------------------------------------------------
! Subroutine: samp_compute_sampling_zs
! Purpose: compute zero-separation sampling
!----------------------------------------------------------------------
subroutine samp_compute_sampling_zs(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0,jc0,ic1,il0i,jb,ildwv,jldwv,ifor
integer,allocatable :: for(:)
real(kind_real) :: norm
logical :: valid
character(len=1024),parameter :: subr = 'samp_compute_sampling_zs'

! Compute subset
if (count(samp%mask_hor_c0)>nam%nc1) then
   write(mpl%info,'(a7,a)') '','Compute horizontal subset C1: '
   call mpl%flush(.false.)
   select case (trim(nam%draw_type))
   case ('random_uniform')
      ! Random draw
      do ic0=1,geom%nc0
         if (samp%mask_hor_c0(ic0)) samp%rh_c0(ic0) = 1.0
      end do
   case ('random_coast')
      ! More points around coasts
      if (all(samp%mask_c0)) call mpl%abort(subr,'random_coast is not relevant if there is no coast')
      samp%rh_c0 = 0.0
      norm = 0.0
      do il0i=1,geom%nl0i
         if (any(.not.samp%mask_c0(:,il0i))) then
            do ic0=1,geom%nc0
               if (samp%mask_c0(ic0,il0i)) then
                  samp%rh_c0(ic0) = samp%rh_c0(ic0)+exp(-geom%mdist(ic0,il0i)/nam%Lcoast)
               else
                  samp%rh_c0(ic0) = samp%rh_c0(ic0)+1.0
               end if
            end do
            norm = norm+1.0
         end if
      end do
      samp%rh_c0 = nam%rcoast+(1.0-nam%rcoast)*(1.0-samp%rh_c0/norm)
   end select

   ! Count forced points
   samp%nfor = 0
   do ildwv=1,nam%nldwv
      ic0 = samp%ldwv_to_c0(ildwv)
      valid = .true.
      do jb=1,geom%mesh%nb
         jc0 = geom%mesh%order(geom%mesh%bnd(jb))
         if (eq(geom%lon(ic0),geom%lon(jc0)).and.eq(geom%lat(ic0),geom%lat(jc0))) valid = .false.
      end do
      do jldwv=1,ildwv-1
         jc0 = samp%ldwv_to_c0(jldwv)
         if (eq(geom%lon(ic0),geom%lon(jc0)).and.eq(geom%lat(ic0),geom%lat(jc0))) valid = .false.
      end do
      if (valid) then
         if (samp%mask_hor_c0(ic0)) samp%nfor = samp%nfor+1
      end if
   end do

   ! Allocation
   allocate(for(samp%nfor))

   if (samp%nfor>0) then
      ! Initialization
      ifor = 0

      ! Add local diagnostic profiles
      do ildwv=1,nam%nldwv
         ic0 = samp%ldwv_to_c0(ildwv)
         valid = .true.
         do jb=1,geom%mesh%nb
            jc0 = geom%mesh%order(geom%mesh%bnd(jb))
            if (eq(geom%lon(ic0),geom%lon(jc0)).and.eq(geom%lat(ic0),geom%lat(jc0))) valid = .false.
         end do
         do jldwv=1,ildwv-1
            jc0 = samp%ldwv_to_c0(jldwv)
            if (eq(geom%lon(ic0),geom%lon(jc0)).and.eq(geom%lat(ic0),geom%lat(jc0))) valid = .false.
         end do
         if (valid) then
            if (samp%mask_hor_c0(ic0)) then
               ifor = ifor+1
               for(ifor) = ic0
            end if
         end if
      end do
   end if

   ! Initialize sampling
   call rng%initialize_sampling(mpl,geom%nc0,geom%lon,geom%lat,samp%mask_hor_c0,samp%nfor,for,samp%rh_c0, &
 & nam%ntry,nam%nrep,nam%nc1,samp%c1_to_c0)
elseif (count(samp%mask_hor_c0)==nam%nc1) then
   ! Keep all remaining points
   ic1 = 0
   do ic0=1,geom%nc0
      if (samp%mask_hor_c0(ic0)) then
         ic1 = ic1+1
         samp%c1_to_c0(ic1) = ic0
      end if
   end do
else
   ! Not enough points remaining in the sampling mask
   call mpl%abort(subr,'not enough points remaining in sampling mask')
end if

end subroutine samp_compute_sampling_zs

!----------------------------------------------------------------------
! Subroutine: samp_compute_sampling_ps
! Purpose: compute positive separation sampling
!----------------------------------------------------------------------
subroutine samp_compute_sampling_ps(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: irmaxloc,jc3,ic1,ir,ic0,jc0,i,nvc0,ivc0,icinf,icsup,ictest
integer,allocatable :: vic0(:)
real(kind_real) :: d
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical :: found

! First class
samp%c1c3_to_c0(:,1) = samp%c1_to_c0

if (nam%nc3>1) then
   ! Initialize
   do jc3=1,nam%nc3
      if (jc3/=1) samp%c1c3_to_c0(:,jc3) = mpl%msv%vali
   end do

   ! Define valid nodes vector
   nvc0 = count(samp%mask_hor_c0)
   allocate(vic0(nvc0))
   ivc0 = 0
   do ic0=1,geom%nc0
      if (samp%mask_hor_c0(ic0)) then
         ivc0 = ivc0+1
         vic0(ivc0) = ic0
      end if
   end do

   ! Sample classes of positive separation
   write(mpl%info,'(a7,a)') '','Compute positive separation sampling: '
   call mpl%flush(.false.)
   call mpl%prog_init(nam%nc3*nam%nc1)
   ir = 0
   irmaxloc = nam%irmax
   do while ((.not.all(mpl%msv%isnoti(samp%c1c3_to_c0))).and.(nvc0>1).and.(ir<=irmaxloc))
      ! Try a random point
      if (mpl%main) call rng%rand_integer(1,nvc0,i)
      call mpl%f_comm%broadcast(i,mpl%ioproc-1)
      ir = ir+1
      jc0 = vic0(i)

      !$omp parallel do schedule(static) private(ic1,ic0,d,jc3,icinf,icsup,found,ictest) firstprivate(x,y,z,v1,v2,va,vp,t)
      do ic1=1,nam%nc1
         ! Allocation
         allocate(x(2))
         allocate(y(2))
         allocate(z(2))
         allocate(v1(3))
         allocate(v2(3))
         allocate(va(3))
         allocate(vp(3))
         allocate(t(4))

         ! Check if there is a valid first point
         if (mpl%msv%isnoti(samp%c1_to_c0(ic1))) then
            ! Compute the distance
            ic0 = samp%c1_to_c0(ic1)
            call sphere_dist(geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),d)

            ! Find the class (dichotomy method)
            if ((d>0.0).and.(d<(real(nam%nc3,kind_real)-0.5)*nam%dc)) then
               jc3 = 1
               icinf = 1
               icsup = nam%nc3
               found = .false.
               do while (.not.found)
                  ! New value
                  ictest = (icsup+icinf)/2

                  ! Update
                  if (d<(real(ictest,kind_real)-0.5)*nam%dc) icsup = ictest
                  if (d>(real(ictest,kind_real)-0.5)*nam%dc) icinf = ictest

                  ! Exit test
                  if (icsup==icinf+1) then
                     if (abs((real(icinf,kind_real)-0.5)*nam%dc-d)<abs((real(icsup,kind_real)-0.5)*nam%dc-d)) then
                        jc3 = icinf
                     else
                        jc3 = icsup
                     end if
                     found = .true.
                  end if
               end do

               ! Find if this class has not been aready filled
               if ((jc3/=1).and.(mpl%msv%isi(samp%c1c3_to_c0(ic1,jc3)))) samp%c1c3_to_c0(ic1,jc3) = jc0
            end if
         end if

         ! Release memory
         deallocate(x)
         deallocate(y)
         deallocate(z)
         deallocate(v1)
         deallocate(v2)
         deallocate(va)
         deallocate(vp)
         deallocate(t)
      end do
      !$omp end parallel do

      ! Update valid nodes vector
      vic0(i) = vic0(nvc0)
      nvc0 = nvc0-1

      ! Update
      mpl%done = pack(mpl%msv%isnoti(samp%c1c3_to_c0),mask=.true.)
      call mpl%prog_print
   end do
   call mpl%prog_final

   ! Release memory
   deallocate(vic0)
end if

end subroutine samp_compute_sampling_ps

!----------------------------------------------------------------------
! Subroutine: samp_compute_sampling_lct
! Purpose: compute LCT sampling
!----------------------------------------------------------------------
subroutine samp_compute_sampling_lct(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: il0,ic1,ic0,jc0,jc3
integer :: nn_index(nam%nc3)
integer :: nc1_loc(0:mpl%nproc),ic1_loc
real(kind_real) :: nn_dist(nam%nc3)

! MPI splitting
call mpl%split(nam%nc1,nc1_loc)

! Initialization
write(mpl%info,'(a7,a)') '','Compute LCT sampling: '
call mpl%flush(.false.)
call mpl%prog_init(nc1_loc(mpl%myproc))

do ic1_loc=1,nc1_loc(mpl%myproc)
   ! MPI offset
   ic1 = sum(nc1_loc(0:mpl%myproc-1))+ic1_loc

   ! Check location validity
   if (mpl%msv%isnoti(samp%c1_to_c0(ic1))) then
      ! Find neighbors
      call geom%tree%find_nearest_neighbors(geom%lon(samp%c1_to_c0(ic1)),geom%lat(samp%c1_to_c0(ic1)), &
    & nam%nc3,nn_index,nn_dist)

      ! Copy neighbor index
      do jc3=1,nam%nc3
         jc0 = nn_index(jc3)
         samp%c1c3_to_c0(ic1,jc3) = nn_index(jc3)
         do il0=1,geom%nl0
            samp%c1c3l0_log(ic1,jc3,il0) = samp%mask_c0(jc0,il0)
         end do
      end do

      if (nam%mask_check) then
         ! Check that great circle to neighbors is not crossing mask boundaries
         do il0=1,geom%nl0
            !$omp parallel do schedule(static) private(jc3,ic0,jc0)
            do jc3=1,nam%nc3
               ! Indices
               ic0 = samp%c1_to_c0(ic1)
               jc0 = samp%c1c3_to_c0(ic1,jc3)

               ! Check arc validity
               call geom%check_arc(mpl,il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),samp%c1c3l0_log(ic1,jc3,il0))
            end do
            !$omp end parallel do
         end do
      end if
   end if

   ! Update
   call mpl%prog_print(ic1_loc)
end do
call mpl%prog_final

! MPI sharing
call mpl%share(nam%nc1,nam%nc3,nc1_loc,samp%c1c3_to_c0)
call mpl%share(nam%nc1,nam%nc3,geom%nl0,nc1_loc,samp%c1c3l0_log)

end subroutine samp_compute_sampling_lct

!----------------------------------------------------------------------
! Subroutine: samp_check_mask
! Purpose: check sampling mask
!----------------------------------------------------------------------
subroutine samp_check_mask(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: il0,jc3,ic1,ic0,jc0,i
integer :: nc1_loc(0:mpl%nproc),ic1_loc
logical :: valid

! First point
do il0=1,geom%nl0
   samp%c1l0_log(:,il0) = samp%mask_c0(samp%c1_to_c0,il0)
end do

! MPL split
call mpl%split(nam%nc1,nc1_loc)

! Second point
if (nam%mask_check) then
   write(mpl%info,'(a7,a)') '','Check mask boundaries: '
   call mpl%flush(.false.)
   call mpl%prog_init(nc1_loc(mpl%myproc)*nam%nc3*geom%nl0)
   i = 0
end if

! Check mask
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1_loc=1,nc1_loc(mpl%myproc)
         ! MPI offset
         ic1 = sum(nc1_loc(0:mpl%myproc-1))+ic1_loc

         ! Indices
         ic0 = samp%c1c3_to_c0(ic1,jc3)
         jc0 = samp%c1c3_to_c0(ic1,1)

         ! Check point index
         valid = mpl%msv%isnoti(ic0).and.mpl%msv%isnoti(jc0)

         if (valid) then
            ! Check mask
            valid = samp%mask_c0(ic0,il0).and.samp%mask_c0(jc0,il0)

            ! Check mask bounds
            if (nam%mask_check) then
               if (valid) call geom%check_arc(mpl,il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),valid)
            end if
         end if
         samp%c1c3l0_log(ic1,jc3,il0) = valid

         ! Update
         if (nam%mask_check) then
            i = i+1
            call mpl%prog_print(i)
         end if
      end do
   end do
end do
if (nam%mask_check) call mpl%prog_final

! MPI sharing
call mpl%share(nam%nc1,nam%nc3,geom%nl0,nc1_loc,samp%c1c3l0_log)

end subroutine samp_check_mask

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_a
! Purpose: compute sampling MPI distribution, halo A
!----------------------------------------------------------------------
subroutine samp_compute_mpi_a(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0a,ic0,ic1a,ic1

! Allocation
allocate(samp%lcheck_c0a(geom%nc0))
allocate(samp%lcheck_c1a(nam%nc1))

! Halo definitions

! Halo A
samp%lcheck_c0a = .false.
samp%lcheck_c1a = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   if (geom%c0_to_proc(ic0)==mpl%myproc) samp%lcheck_c0a(ic0) = .true.
end do
do ic1=1,nam%nc1
   ic0 = samp%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) samp%lcheck_c1a(ic1) = .true.
end do

! Halo sizes
samp%nc1a = count(samp%lcheck_c1a)

! Global <-> local conversions for fields

! Halo A
allocate(samp%c1a_to_c1(samp%nc1a))
allocate(samp%c1_to_c1a(nam%nc1))
ic1a = 0
do ic1=1,nam%nc1
   if (samp%lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      samp%c1a_to_c1(ic1a) = ic1
      samp%c1_to_c1a(ic1) = ic1a
   end if
end do

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0 =        ',geom%nc0
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0a =       ',geom%nc0a
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nl0 =        ',geom%nl0
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1 =        ',nam%nc1
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1a =       ',samp%nc1a
call mpl%flush

end subroutine samp_compute_mpi_a

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_ab
! Purpose: compute sampling MPI distribution, halos A-B
!----------------------------------------------------------------------
subroutine samp_compute_mpi_ab(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: iproc,ic0,ic2a,ic2b,ic2,jc2,i_s,i_s_loc,h_n_s_max,il0i,il0,h_n_s_max_loc,jc1,kc1a,kc1,kc0
integer,allocatable :: interph_lg(:,:),nn_c1_index(:)
real(kind_real) :: lon_c2(nam%nc2),lat_c2(nam%nc2)
real(kind_real),allocatable :: lon_c1(:),lat_c1(:),nn_c1_dist(:)
logical,allocatable :: mask_c1(:)
type(tree_type) :: tree

! Allocation
h_n_s_max = 0
do il0i=1,geom%nl0i
   h_n_s_max = max(h_n_s_max,samp%hfull(il0i)%n_s)
end do
allocate(samp%c2_to_proc(nam%nc2))
allocate(samp%proc_to_nc2a(mpl%nproc))
allocate(samp%h(geom%nl0i))
allocate(samp%lcheck_c2a(nam%nc2))
allocate(samp%lcheck_c2b(nam%nc2))
allocate(samp%lcheck_h(h_n_s_max,geom%nl0i))

! Halo definitions

! Halo A
samp%lcheck_c2a = .false.
do ic2=1,nam%nc2
   ic0 = samp%c2_to_c0(ic2)
   if (geom%c0_to_proc(ic0)==mpl%myproc) samp%lcheck_c2a(ic2) = .true.
end do

! Halo B
samp%lcheck_h = .false.
samp%lcheck_c2b = .false.
do il0i=1,geom%nl0i
   do i_s=1,samp%hfull(il0i)%n_s
      ic0 = samp%hfull(il0i)%row(i_s)
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%myproc) then
         jc2 = samp%hfull(il0i)%col(i_s)
         samp%lcheck_h(i_s,il0i) = .true.
         samp%lcheck_c2b(jc2) = .true.
      end if
   end do
end do

! Halo sizes
samp%nc2a = count(samp%lcheck_c2a)
do il0i=1,geom%nl0i
   samp%h(il0i)%n_s = count(samp%lcheck_h(:,il0i))
end do
samp%nc2b = count(samp%lcheck_c2b)

! Global <-> local conversions for fields

! Halo A
allocate(samp%c2a_to_c2(samp%nc2a))
allocate(samp%c2_to_c2a(nam%nc2))
ic2a = 0
do ic2=1,nam%nc2
   if (samp%lcheck_c2a(ic2)) then
      ic2a = ic2a+1
      samp%c2a_to_c2(ic2a) = ic2
   end if
end do
call mpl%glb_to_loc_index(samp%nc2a,samp%c2a_to_c2,nam%nc2,samp%c2_to_c2a)

! Halo B
allocate(samp%c2b_to_c2(samp%nc2b))
allocate(samp%c2_to_c2b(nam%nc2))
ic2b = 0
do ic2=1,nam%nc2
   if (samp%lcheck_c2b(ic2)) then
      ic2b = ic2b+1
      samp%c2b_to_c2(ic2b) = ic2
      samp%c2_to_c2b(ic2) = ic2b
   end if
end do

! Inter-halo conversions
allocate(samp%c2a_to_c2b(samp%nc2a))
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic2b = samp%c2_to_c2b(ic2)
   samp%c2a_to_c2b(ic2a) = ic2b
end do

! Global <-> local conversions for data
h_n_s_max_loc = 0
do il0i=1,geom%nl0i
   h_n_s_max_loc = max(h_n_s_max_loc,samp%h(il0i)%n_s)
end do
allocate(interph_lg(h_n_s_max_loc,geom%nl0i))
do il0i=1,geom%nl0i
   i_s_loc = 0
   do i_s=1,samp%hfull(il0i)%n_s
      if (samp%lcheck_h(i_s,il0i)) then
         i_s_loc = i_s_loc+1
         interph_lg(i_s_loc,il0i) = i_s
      end if
   end do
end do

! Local data

! Horizontal interpolation
do il0i=1,geom%nl0i
   write(samp%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
   samp%h(il0i)%n_src = samp%nc2b
   samp%h(il0i)%n_dst = geom%nc0a
   call samp%h(il0i)%alloc
   do i_s_loc=1,samp%h(il0i)%n_s
      i_s = interph_lg(i_s_loc,il0i)
      samp%h(il0i)%row(i_s_loc) = geom%c0_to_c0a(samp%hfull(il0i)%row(i_s))
      samp%h(il0i)%col(i_s_loc) = samp%c2_to_c2b(samp%hfull(il0i)%col(i_s))
      samp%h(il0i)%S(i_s_loc) = samp%hfull(il0i)%S(i_s)
   end do
   call samp%h(il0i)%reorder(mpl)
end do

! Release memory
deallocate(interph_lg)

! MPI splitting
do ic2=1,nam%nc2
   ic0 = samp%c2_to_c0(ic2)
   samp%c2_to_proc(ic2) = geom%c0_to_proc(ic0)
end do
do iproc=1,mpl%nproc
   samp%proc_to_nc2a(iproc) = count(samp%c2_to_proc==iproc)
end do

! Setup communications
call samp%com_AB%setup(mpl,'com_AB',nam%nc2,samp%nc2a,samp%nc2b,samp%c2b_to_c2,samp%c2a_to_c2b,samp%c2_to_proc, &
 & samp%c2_to_c2a)

! Mask on subset Sc2, halo A
allocate(samp%mask_c2a(samp%nc2a,geom%nl0))
do il0=1,geom%nl0
   do ic2a=1,samp%nc2a
      ic2 = samp%c2a_to_c2(ic2a)
      samp%mask_c2a(ic2a,il0) = samp%mask_c2(ic2,il0)
   end do
end do


! Find nearest neighbors

! Allocation
allocate(samp%nn_c2a_index(nam%nc2,samp%nc2a))
allocate(samp%nn_c2a_dist(nam%nc2,samp%nc2a))
call tree%alloc(mpl,nam%nc2)

! Initialization
lon_c2 = geom%lon(samp%c2_to_c0)
lat_c2 = geom%lat(samp%c2_to_c0)
call tree%init(lon_c2,lat_c2)

! Find nearest neighbors
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic0 = samp%c2_to_c0(ic2)
   call tree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nam%nc2,samp%nn_c2a_index(:,ic2a),samp%nn_c2a_dist(:,ic2a))
end do

! Release memory
call tree%dealloc

! Compute local masks
if (nam%new_vbal.or.nam%local_diag) then
   ! Allocation
   if (nam%new_vbal) allocate(samp%vbal_mask(samp%nc1a,nam%nc2))
   if (nam%local_diag) allocate(samp%local_mask(samp%nc1a,nam%nc2))

   ! Initialization
   if (nam%new_vbal) samp%vbal_mask = .false.
   if (nam%local_diag) samp%local_mask = .false.

   ! Allocation
   allocate(nn_c1_index(nam%nc1))
   allocate(nn_c1_dist(nam%nc1))
   allocate(lon_c1(nam%nc1))
   allocate(lat_c1(nam%nc1))
   allocate(mask_c1(nam%nc1))

   ! Initialization
   mask_c1 = any(samp%c1l0_log,dim=2)
   lon_c1 = geom%lon(samp%c1_to_c0)
   lat_c1 = geom%lat(samp%c1_to_c0)

   ! Allocation
   call tree%alloc(mpl,nam%nc1,mask=mask_c1)

   ! Initialization
   call tree%init(lon_c1,lat_c1)

   do ic2=1,nam%nc2
      ! Find nearest neighbors
      ic0 = samp%c2_to_c0(ic2)
      call tree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nam%nc1,nn_c1_index,nn_c1_dist)
      do jc1=1,nam%nc1
         kc1 = nn_c1_index(jc1)
         kc0 = samp%c1_to_c0(kc1)
         if (geom%c0_to_proc(kc0)==mpl%myproc) then
            kc1a = samp%c1_to_c1a(kc1)
            if (nam%new_vbal) samp%vbal_mask(kc1a,ic2) = (jc1==1).or.(nn_c1_dist(jc1)<nam%vbal_rad)
            if (nam%local_diag) samp%local_mask(kc1a,ic2) = (jc1==1).or.(nn_c1_dist(jc1)<nam%local_rad)
         end if
      end do
   end do

   ! Release memory
   deallocate(nn_c1_index)
   deallocate(nn_c1_dist)
   deallocate(lon_c1)
   deallocate(lat_c1)
   deallocate(mask_c1)
   call tree%dealloc
end if

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2 =        ',nam%nc2
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2a =       ',samp%nc2a
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2b =       ',samp%nc2b
do il0i=1,geom%nl0i
   write(mpl%info,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',samp%h(il0i)%n_s
   call mpl%flush
end do

end subroutine samp_compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_c
! Purpose: compute sampling MPI distribution, halo C
!----------------------------------------------------------------------
subroutine samp_compute_mpi_c(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp        ! Sampling
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry

! Local variables
integer :: jc3,ic0,ic0a,ic0c,ic1,ic1a,its,il0,d_n_s_max,d_n_s_max_loc,i_s,i_s_loc
integer,allocatable :: interpd_lg(:,:,:)
real(kind_real),allocatable :: adv_lon(:,:),adv_lat(:,:),lon_c1(:),lat_c1(:)
type(linop_type),allocatable :: dfull(:,:)

if (nam%adv_diag) then
   write(mpl%info,'(a7,a)') '','Compute advection interpolation'
   call mpl%flush

   ! Allocation
   allocate(dfull(geom%nl0,nam%nts))
   allocate(adv_lon(geom%nc0,geom%nl0))
   allocate(adv_lat(geom%nc0,geom%nl0))
   allocate(lon_c1(nam%nc1))
   allocate(lat_c1(nam%nc1))

   ! Prepare advection interpolation
   do its=1,nam%nts
      ! Local to global
      call mpl%loc_to_glb(geom%nl0,geom%nc0a,samp%adv_lon(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,adv_lon)
      call mpl%loc_to_glb(geom%nl0,geom%nc0a,samp%adv_lat(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,adv_lat)

      do il0=1,geom%nl0
         ! Copy Sc1 points
         do ic1=1,nam%nc1
            ic0 = samp%c1_to_c0(ic1)
            lon_c1(ic1) = adv_lon(ic0,il0)
            lat_c1(ic1) = adv_lat(ic0,il0)
         end do

         ! Compute interpolation
         call dfull(il0,its)%interp(mpl,geom%mesh,geom%tree,geom%nc0,geom%mask_c0(:,il0),nam%nc1,lon_c1,lat_c1, &
       & samp%c1l0_log(:,il0),nam%diag_interp)
      end do
   end do

   ! Allocation
   d_n_s_max = 0
   do its=1,nam%nts
      do il0=1,geom%nl0
         d_n_s_max = max(d_n_s_max,dfull(il0,its)%n_s)
      end do
   end do
   allocate(samp%lcheck_d(d_n_s_max,geom%nl0,nam%nts))
   allocate(samp%d(geom%nl0,nam%nts))

   ! Release memory
   deallocate(adv_lon)
   deallocate(adv_lat)
   deallocate(lon_c1)
   deallocate(lat_c1)
end if
allocate(samp%lcheck_c0c(geom%nc0))

! Halo C
samp%lcheck_c0c = samp%lcheck_c0a
do jc3=1,nam%nc3
   do ic1a=1,samp%nc1a
      ic1 = samp%c1a_to_c1(ic1a)
      if (any(samp%c1c3l0_log(ic1,jc3,:))) then
         ic0 = samp%c1c3_to_c0(ic1,jc3)
         samp%lcheck_c0c(ic0) = .true.
      end if
   end do
end do
if (nam%adv_diag) then
   samp%lcheck_d = .false.
   do its=1,nam%nts
      do il0=1,geom%nl0
         do i_s=1,dfull(il0,its)%n_s
            ic0 = dfull(il0,its)%col(i_s)
            ic1 = dfull(il0,its)%row(i_s)
            if (samp%lcheck_c1a(ic1)) then
               samp%lcheck_c0c(ic0) = .true.
               samp%lcheck_d(i_s,il0,its) = .true.
            end if
         end do
      end do
   end do
end if
samp%nc0c = count(samp%lcheck_c0c)
if (nam%adv_diag) then
   do its=1,nam%nts
      do il0=1,geom%nl0
         samp%d(il0,its)%n_s = count(samp%lcheck_d(:,il0,its))
      end do
   end do
end if

! Global <-> local conversions for fields

! Halo C
allocate(samp%c0c_to_c0(samp%nc0c))
allocate(samp%c0_to_c0c(geom%nc0))
samp%c0_to_c0c = mpl%msv%vali
ic0c = 0
do ic0=1,geom%nc0
   if (samp%lcheck_c0c(ic0)) then
      ic0c = ic0c+1
      samp%c0c_to_c0(ic0c) = ic0
      samp%c0_to_c0c(ic0) = ic0c
   end if
end do

! Inter-halo conversions
allocate(samp%c0a_to_c0c(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0c = samp%c0_to_c0c(ic0)
   samp%c0a_to_c0c(ic0a) = ic0c
end do

if (nam%adv_diag) then
   ! Global <-> local conversions for data
   d_n_s_max_loc = 0
   do its=1,nam%nts
      do il0=1,geom%nl0
         d_n_s_max_loc = max(d_n_s_max_loc,samp%d(il0,its)%n_s)
      end do
   end do
   allocate(interpd_lg(d_n_s_max_loc,geom%nl0,nam%nts))
   do its=1,nam%nts
      do il0=1,geom%nl0
         i_s_loc = 0
         do i_s=1,dfull(il0,its)%n_s
            if (samp%lcheck_d(i_s,il0,its)) then
               i_s_loc = i_s_loc+1
               interpd_lg(i_s_loc,il0,its) = i_s
            end if
         end do
      end do
   end do

   ! Local data
   do its=1,nam%nts
      do il0=1,geom%nl0
         write(samp%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
         samp%d(il0,its)%n_src = samp%nc0c
         samp%d(il0,its)%n_dst = samp%nc1a
         call samp%d(il0,its)%alloc
         do i_s_loc=1,samp%d(il0,its)%n_s
            i_s = interpd_lg(i_s_loc,il0,its)
            samp%d(il0,its)%row(i_s_loc) = samp%c1_to_c1a(dfull(il0,its)%row(i_s))
            samp%d(il0,its)%col(i_s_loc) = samp%c0_to_c0c(dfull(il0,its)%col(i_s))
            samp%d(il0,its)%S(i_s_loc) = dfull(il0,its)%S(i_s)
         end do
         call samp%d(il0,its)%reorder(mpl)
      end do
   end do

   ! Release memory
   deallocate(interpd_lg)
end if

! Setup communications
call samp%com_AC%setup(mpl,'com_AC',geom%nc0,geom%nc0a,samp%nc0c,samp%c0c_to_c0,samp%c0a_to_c0c,geom%c0_to_proc,geom%c0_to_c0a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0c =      ',samp%nc0c
call mpl%flush

end subroutine samp_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_f
! Purpose: compute sampling MPI distribution, halo F
!----------------------------------------------------------------------
subroutine samp_compute_mpi_f(samp,mpl,nam)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: ic2a,ic2f,ic2,jc2,kc2

! Allocation
allocate(samp%lcheck_c2f(nam%nc2))

! Halo definitions

! Halo F
samp%lcheck_c2f = samp%lcheck_c2a
do ic2a=1,samp%nc2a
   jc2 = 1
   do while (samp%nn_c2a_dist(jc2,ic2a)<nam%diag_rhflt)
      kc2 = samp%nn_c2a_index(jc2,ic2a)
      samp%lcheck_c2f(kc2) = .true.
      jc2 = jc2+1
      if (jc2>nam%nc2) exit
   end do
end do
samp%nc2f = count(samp%lcheck_c2f)

! Global <-> local conversions for fields

! Halo F
allocate(samp%c2f_to_c2(samp%nc2f))
allocate(samp%c2_to_c2f(nam%nc2))
samp%c2_to_c2f = mpl%msv%vali
ic2f = 0
do ic2=1,nam%nc2
   if (samp%lcheck_c2f(ic2)) then
      ic2f = ic2f+1
      samp%c2f_to_c2(ic2f) = ic2
      samp%c2_to_c2f(ic2) = ic2f
   end if
end do

! Inter-halo conversions
allocate(samp%c2a_to_c2f(samp%nc2a))
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic2f = samp%c2_to_c2f(ic2)
   samp%c2a_to_c2f(ic2a) = ic2f
end do

! Setup communications
call samp%com_AF%setup(mpl,'com_AF',nam%nc2,samp%nc2a,samp%nc2f,samp%c2f_to_c2,samp%c2a_to_c2f,samp%c2_to_proc, &
 & samp%c2_to_c2a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2f =       ',samp%nc2f
call mpl%flush

end subroutine samp_compute_mpi_f

!----------------------------------------------------------------------
! Subroutine: samp_diag_filter
! Purpose: filter diagnostics
!----------------------------------------------------------------------
subroutine samp_diag_filter(samp,mpl,nam,filter_type,rflt,diag,val)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp                   ! Sampling
type(mpl_type),intent(inout) :: mpl                   ! MPI data
type(nam_type),intent(in) :: nam                      ! Namelist
character(len=*),intent(in) :: filter_type            ! Filter type
real(kind_real),intent(in) :: rflt                    ! Filter support radius
real(kind_real),intent(inout) :: diag(samp%nc2a)      ! Filtered diagnostic
real(kind_real),intent(in),optional :: val(samp%nc2a) ! Useful value for filtering

! Local variables
integer :: ic2a,jc2,nc2eff,kc2,kc2_glb
integer,allocatable :: order(:)
real(kind_real) :: distnorm,norm,wgt
real(kind_real),allocatable :: diag_glb(:),diag_eff(:),diag_eff_dist(:)
real(kind_real),allocatable :: val_glb(:),val_eff(:)
logical :: nam_rad
character(len=1024),parameter :: subr = 'samp_diag_filter'

! Check radius
nam_rad = .not.(abs(rflt-nam%diag_rhflt)>0.0)

if (rflt>0.0) then
   if (nam_rad) then
      ! Allocation
      allocate(diag_glb(samp%nc2f))
      if (present(val)) allocate(val_glb(samp%nc2f))

      ! Communication
      call samp%com_AF%ext(mpl,diag,diag_glb)
      if (present(val)) call samp%com_AF%ext(mpl,diag,diag_glb)
   else
      ! Allocation
      allocate(diag_glb(nam%nc2))
      if (present(val)) allocate(val_glb(nam%nc2))

      ! Local to global
      call mpl%loc_to_glb(samp%nc2a,diag,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.true.,diag_glb)
      if (present(val)) call mpl%loc_to_glb(samp%nc2a,val,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.true.,val_glb)
   end if

   !$omp parallel do schedule(static) private(ic2a,nc2eff,jc2,kc2,kc2_glb,distnorm,norm,wgt), &
   !$omp&                             firstprivate(diag_eff,diag_eff_dist,val_eff,order)
   do ic2a=1,samp%nc2a
      ! Allocation
      allocate(diag_eff(nam%nc2))
      allocate(diag_eff_dist(nam%nc2))
      if (present(val)) allocate(val_eff(nam%nc2))

      ! Build diag_eff of valid points
      nc2eff = 0
      jc2 = 1
      do while (inf(samp%nn_c2a_dist(jc2,ic2a),rflt))
         ! Check the point validity
         kc2 = samp%nn_c2a_index(jc2,ic2a)
         if (nam_rad) then
            kc2_glb = samp%c2_to_c2f(kc2)
         else
            kc2_glb = kc2
         end if
         if (mpl%msv%isnotr(diag_glb(kc2_glb))) then
            nc2eff = nc2eff+1
            diag_eff(nc2eff) = diag_glb(kc2_glb)
            diag_eff_dist(nc2eff) = samp%nn_c2a_dist(jc2,ic2a)
            if (present(val)) val_eff(nc2eff) = val_glb(kc2_glb)
         end if
         jc2 = jc2+1
         if (jc2>nam%nc2) exit
      end do

      ! Apply filter
      if (nc2eff>0) then
         select case (trim(filter_type))
         case ('average')
            ! Compute average
            diag(ic2a) = sum(diag_eff(1:nc2eff))/real(nc2eff,kind_real)
         case ('gc99')
            ! Gaspari-Cohn (1999) kernel
            diag(ic2a) = 0.0
            norm = 0.0
            do jc2=1,nc2eff
               distnorm = diag_eff_dist(jc2)/rflt
               wgt = gc99(mpl,distnorm)
               diag(ic2a) = diag(ic2a)+wgt*diag_eff(jc2)
               norm = norm+wgt
            end do
            if (norm>0.0) diag(ic2a) = diag(ic2a)/norm
         case ('median')
            ! Compute median
            allocate(order(nc2eff))
            if (present(val)) then
               ! Use external value
               call qsort(nc2eff,val_eff(1:nc2eff),order)
               diag_eff = diag_eff(order)
            else
               ! Use diagnostic value
               call qsort(nc2eff,diag_eff(1:nc2eff),order)
            end if
            if (mod(nc2eff,2)==0) then
               diag(ic2a) = 0.5*(diag_eff(nc2eff/2)+diag_eff(nc2eff/2+1))
            else
               diag(ic2a) = diag_eff((nc2eff+1)/2)
            end if
            deallocate(order)
         case default
            ! Wrong filter
            call mpl%abort(subr,'wrong filter type')
         end select
      else
         diag(ic2a) = mpl%msv%valr
      end if

      ! Release memory
      deallocate(diag_eff)
      deallocate(diag_eff_dist)
      if (present(val)) deallocate(val_eff)
   end do
   !$omp end parallel do

   ! Release memory
   deallocate(diag_glb)
   if (present(val)) deallocate(val_glb)
end if

end subroutine samp_diag_filter

!----------------------------------------------------------------------
! Subroutine: samp_diag_fill
! Purpose: fill diagnostics missing values
!----------------------------------------------------------------------
subroutine samp_diag_fill(samp,mpl,nam,diag)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp              ! Sampling
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
real(kind_real),intent(inout) :: diag(samp%nc2a) ! Filtered diagnostic

! Local variables
integer :: nmsr,nmsr_tot,ic2a,ic2,jc2,kc2
real(kind_real),allocatable :: diag_glb(:)

! Count missing points
if (samp%nc2a>0) then
   nmsr = count(mpl%msv%isr(diag))
else
   nmsr = 0
end if
call mpl%f_comm%allreduce(nmsr,nmsr_tot,fckit_mpi_sum())


if (nmsr_tot>0) then
   ! Allocation
   allocate(diag_glb(nam%nc2))

   ! Local to global
   call mpl%loc_to_glb(samp%nc2a,diag,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.true.,diag_glb)

   do ic2a=1,samp%nc2a
      ic2 = samp%c2a_to_c2(ic2a)
      jc2 = 1
      do while (mpl%msv%isr(diag(ic2a)))
         kc2 = samp%nn_c2a_index(jc2,ic2a)
         if (mpl%msv%isnotr(diag_glb(kc2))) diag(ic2a) = diag_glb(kc2)
         jc2 = jc2+1
         if (jc2>nam%nc2) exit
      end do
   end do

   ! Release memory
   deallocate(diag_glb)
end if

end subroutine samp_diag_fill

end module type_samp
