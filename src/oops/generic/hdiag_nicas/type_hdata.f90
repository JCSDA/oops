!----------------------------------------------------------------------
! Module: type_hdata
!> Purpose: sample data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_hdata

use netcdf
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg
use tools_display, only: prog_init,prog_print,msgerror,msgwarning,black,green,peach
use tools_func, only: gc99,sphere_dist,vector_product,vector_triple_product
use tools_icos, only: closest_icos,build_icos
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsi,isnotmsr,isanynotmsr,isallnotmsr
use tools_nc, only: ncerr,ncfloat
use tools_qsort, only: qsort
use tools_stripack, only: trans
use type_com, only: com_type
use type_ctree, only: ctree_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl
use type_nam, only: nam_type
use type_rng, only: rng

implicit none

! HDIAG data derived type
type hdata_type
   ! Sampling
   integer,allocatable :: c1_to_c0(:)               !< First sampling index
   logical,allocatable :: c1l0_log(:,:)             !< Log for the first sampling index
   integer,allocatable :: c1c3_to_c0(:,:)           !< Second horizontal sampling index
   logical,allocatable :: c1c3l0_log(:,:,:)         !< Log for the second horizontal sampling index
   integer :: nc2                                   !< Subgrid size
   integer,allocatable :: c2_to_c1(:)               !< Subgrid to diagnostic points
   integer,allocatable :: c2_to_c0(:)               !< Subgrid to grid

   ! Cover tree and nearest neighbors
   logical,allocatable ::  local_mask(:,:,:)        !< Local mask
   logical,allocatable ::  displ_mask(:,:,:)        !< Displacement mask
   integer,allocatable :: nn_c2_index(:,:,:)        !< Nearest diagnostic neighbors from diagnostic points
   real(kind_real),allocatable :: nn_c2_dist(:,:,:) !< Nearest diagnostic neighbors distance from diagnostic points
   integer,allocatable :: nn_ldwv_index(:)          !< Nearest diagnostic neighbors for local diagnostics profiles

   ! Sampling mesh
   type(mesh_type) :: mesh                          !< Sampling mesh

   ! Displacement
   real(kind_real),allocatable :: displ_lon(:,:,:)  !< Interpolated displaced longitude
   real(kind_real),allocatable :: displ_lat(:,:,:)  !< Interpolated displaced latitude

   ! Interpolations
   type(linop_type),allocatable :: hfull(:)         !< Horizontal interpolation from Sc2 to Sc0 (global)
   type(linop_type),allocatable :: h(:)             !< Horizontal interpolation from Sc2 to Sc0 (local)
   type(linop_type),allocatable :: d(:,:)           !< Displacement interpolation

   ! MPI distribution
   integer :: nc0c                                  !< Number of points in subset Sc0, halo C
   integer :: nc0d                                  !< Number of points in subset Sc0, halo D
   integer :: nc1a                                  !< Number of points in subset Sc1, halo A
   integer :: nc2a                                  !< Number of points in subset Sc2, halo A
   integer :: nc2b                                  !< Number of points in subset Sc2, halo B
   logical,allocatable :: lcheck_c0a(:)             !< Detection of halo A on subset Sc0
   logical,allocatable :: lcheck_c0c(:)             !< Detection of halo C on subset Sc0
   logical,allocatable :: lcheck_c0d(:)             !< Detection of halo D on subset Sc0
   logical,allocatable :: lcheck_c1a(:)             !< Detection of halo A on subset Sc1
   logical,allocatable :: lcheck_c2a(:)             !< Detection of halo A on subset Sc2
   logical,allocatable :: lcheck_c2b(:)             !< Detection of halo B on subset Sc2
   logical,allocatable :: lcheck_h(:,:)             !< Detection of horizontal interpolation coefficients
   logical,allocatable :: lcheck_d(:,:,:)           !< Detection of displacement interpolation coefficients
   integer,allocatable :: c0c_to_c0(:)              !< Subset Sc0, halo C to global
   integer,allocatable :: c0_to_c0c(:)              !< Subset Sc0, global to halo C
   integer,allocatable :: c0a_to_c0c(:)             !< Subset Sc0, halo A to halo C
   integer,allocatable :: c0d_to_c0(:)              !< Subset Sc0, halo D to global
   integer,allocatable :: c0_to_c0d(:)              !< Subset Sc0, global to halo D
   integer,allocatable :: c0a_to_c0d(:)             !< Subset Sc0, halo A to halo D
   integer,allocatable :: c1a_to_c1(:)              !< Subset Sc1, halo A to global
   integer,allocatable :: c1_to_c1a(:)              !< Subset Sc1, global to halo A
   integer,allocatable :: c2a_to_c2(:)              !< Subset Sc2, halo A to global
   integer,allocatable :: c2_to_c2a(:)              !< Subset Sc2, global to halo A
   integer,allocatable :: c2b_to_c2(:)              !< Subset Sc2, halo B to global
   integer,allocatable :: c2_to_c2b(:)              !< Subset Sc2, global to halo B
   integer,allocatable :: c2a_to_c2b(:)             !< Subset Sc2, halo A to halo B
   integer,allocatable :: c2_to_proc(:)             !< Subset Sc2, global to processor
   integer,allocatable :: proc_to_nc2a(:)           !< Number of points in subset Sc2, halo A, for each processor
   type(com_type) :: com_AC                         !< Communication between halos A and C
   type(com_type) :: com_AB                         !< Communication between halos A and B
   type(com_type) :: com_AD                         !< Communication between halos A and D
contains
   procedure :: alloc => hdata_alloc
   procedure :: dealloc => hdata_dealloc
   procedure :: read => hdata_read
   procedure :: write => hdata_write
   procedure :: setup_sampling => hdata_setup_sampling
   procedure :: compute_sampling_zs => hdata_compute_sampling_zs
   procedure :: compute_sampling_ps => hdata_compute_sampling_ps
   procedure :: compute_sampling_lct => hdata_compute_sampling_lct
   procedure :: compute_sampling_mask => hdata_compute_sampling_mask
   procedure :: compute_mpi_a => hdata_compute_mpi_a
   procedure :: compute_mpi_ab => hdata_compute_mpi_ab
   procedure :: compute_mpi_d => hdata_compute_mpi_d
   procedure :: compute_mpi_c => hdata_compute_mpi_c
   procedure :: diag_filter => hdata_diag_filter
   procedure :: diag_com_lg => hdata_diag_com_lg
end type hdata_type

integer,parameter :: irmax = 10000 !< Maximum number of random number draws

private
public :: hdata_type

contains

!----------------------------------------------------------------------
! Subroutine: hdata_alloc
!> Purpose: hdata object allocation
!----------------------------------------------------------------------
subroutine hdata_alloc(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Allocation
allocate(hdata%c1_to_c0(nam%nc1))
allocate(hdata%c1l0_log(nam%nc1,geom%nl0))
allocate(hdata%c1c3_to_c0(nam%nc1,nam%nc3))
allocate(hdata%c1c3l0_log(nam%nc1,nam%nc3,geom%nl0))
if (nam%local_diag.or.nam%displ_diag) then
   allocate(hdata%c2_to_c1(hdata%nc2))
   allocate(hdata%c2_to_c0(hdata%nc2))
   allocate(hdata%local_mask(nam%nc1,hdata%nc2,geom%nl0i))
   allocate(hdata%displ_mask(nam%nc1,hdata%nc2,geom%nl0i))
   allocate(hdata%nn_c2_index(hdata%nc2,hdata%nc2,geom%nl0i))
   allocate(hdata%nn_c2_dist(hdata%nc2,hdata%nc2,geom%nl0i))
   allocate(hdata%hfull(geom%nl0i))
end if
if (nam%displ_diag) then
   allocate(hdata%displ_lon(geom%nc0,geom%nl0,nam%nts))
   allocate(hdata%displ_lat(geom%nc0,geom%nl0,nam%nts))
end if

! Initialization
call msi(hdata%c1_to_c0)
hdata%c1l0_log = .false.
call msi(hdata%c1c3_to_c0)
hdata%c1c3l0_log = .false.
if (nam%local_diag.or.nam%displ_diag) then
   call msi(hdata%c2_to_c1)
   call msi(hdata%c2_to_c0)
   hdata%local_mask = .false.
   hdata%displ_mask = .false.
   call msi(hdata%nn_c2_index)
   call msr(hdata%nn_c2_dist)
end if
if (nam%displ_diag) then
   call msr(hdata%displ_lon)
   call msr(hdata%displ_lat)
end if

end subroutine hdata_alloc

!----------------------------------------------------------------------
! Subroutine: hdata_dealloc
!> Purpose: hdata object deallocation
!----------------------------------------------------------------------
subroutine hdata_dealloc(hdata,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: il0

! Release memory
if (allocated(hdata%c1_to_c0)) deallocate(hdata%c1_to_c0)
if (allocated(hdata%c1l0_log)) deallocate(hdata%c1l0_log)
if (allocated(hdata%c1c3_to_c0)) deallocate(hdata%c1c3_to_c0)
if (allocated(hdata%c1c3l0_log)) deallocate(hdata%c1c3l0_log)
if (allocated(hdata%c2_to_c1)) deallocate(hdata%c2_to_c1)
if (allocated(hdata%c2_to_c0)) deallocate(hdata%c2_to_c0)
if (allocated(hdata%c2a_to_c2)) deallocate(hdata%c2a_to_c2)
if (allocated(hdata%c2_to_c2a)) deallocate(hdata%c2_to_c2a)
if (allocated(hdata%local_mask)) deallocate(hdata%local_mask)
if (allocated(hdata%displ_mask)) deallocate(hdata%displ_mask)
if (allocated(hdata%nn_c2_index)) deallocate(hdata%nn_c2_index)
if (allocated(hdata%nn_c2_dist)) deallocate(hdata%nn_c2_dist)
if (allocated(hdata%hfull)) then
   do il0=1,geom%nl0
      call hdata%hfull(il0)%dealloc
   end do
   deallocate(hdata%hfull)
end if
if (allocated(hdata%h)) then
   do il0=1,geom%nl0
      call hdata%h(il0)%dealloc
   end do
   deallocate(hdata%h)
end if
if (allocated(hdata%nn_ldwv_index)) deallocate(hdata%nn_ldwv_index)
if (allocated(hdata%displ_lon)) deallocate(hdata%displ_lon)
if (allocated(hdata%displ_lat)) deallocate(hdata%displ_lat)

end subroutine hdata_dealloc

!----------------------------------------------------------------------
! Subroutine: hdata_read
!> Purpose: read hdata object
!----------------------------------------------------------------------
subroutine hdata_read(hdata,nam,geom,ios)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(inout) :: nam      !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
integer,intent(out) :: ios               !< Status flag

! Local variables
integer :: il0,il0i,ic1,jc3,ic2
integer :: nl0_test,nl0r_test,nc_test,nc1_test,nc2_test
integer :: info,ncid,nl0_id,nc_id,nc1_id,nc2_id
integer :: c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id,local_mask_id,displ_mask_id,nn_c2_index_id,nn_c2_dist_id
integer :: c1l0_logint(nam%nc1,geom%nl0),c1c3l0_logint(nam%nc1,nam%nc3,geom%nl0)
integer,allocatable :: local_maskint(:,:),displ_maskint(:,:)
character(len=3) :: il0ichar
character(len=1024) :: subr = 'hdata_read'

! Initialization
ios = 0

! Open file
info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',nf90_nowrite,ncid)
if (info/=nf90_noerr) then
   call msgwarning('cannot find HDIAG data to read, recomputing HDIAG sampling')
   nam%sam_write = .true.
   ios = 1
   return
end if

! Check dimensions
call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=nl0_test))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'nl0r',nl0r_test))
call ncerr(subr,nf90_inq_dimid(ncid,'nc',nc_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc_id,len=nc_test))
call ncerr(subr,nf90_inq_dimid(ncid,'nc1',nc1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc1_id,len=nc1_test))
if (nam%local_diag.or.nam%displ_diag) then
   info = nf90_inq_dimid(ncid,'nc2_1',nc2_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nc2_id,len=nc2_test))
   else
      call msgwarning('cannot find nc2 when reading HDIAG sampling, recomputing HDIAG sampling')
      nam%sam_write = .true.
      ios = 2
   end if
end if
if ((geom%nl0/=nl0_test).or.(nam%nl0r/=nl0r_test).or.(nam%nc3/=nc_test).or.(nam%nc1/=nc1_test)) then
   call msgwarning('wrong dimension when reading HDIAG sampling, recomputing HDIAG sampling')
   nam%sam_write = .true.
   call ncerr(subr,nf90_close(ncid))
   ios = 1
   return
end if
if (nam%local_diag.or.nam%displ_diag) then
   if (hdata%nc2/=nc2_test) then
      call msgwarning('wrong dimension when reading HDIAG sampling, recomputing HDIAG sampling')
      nam%sam_write = .true.
      ios = 2
   end if
end if

write(mpl%unit,'(a7,a)') '','Read HDIAG sampling'

! Get arrays ID
call ncerr(subr,nf90_inq_varid(ncid,'c1_to_c0',c1_to_c0_id))
call ncerr(subr,nf90_inq_varid(ncid,'c1l0_log',c1l0_log_id))
call ncerr(subr,nf90_inq_varid(ncid,'c1c3_to_c0',c1c3_to_c0_id))
call ncerr(subr,nf90_inq_varid(ncid,'c1c3l0_log',c1c3l0_log_id))
if ((ios==0).and.(nam%local_diag.or.nam%displ_diag)) then
   call ncerr(subr,nf90_inq_varid(ncid,'c2_to_c1',c2_to_c1_id))
   call ncerr(subr,nf90_inq_varid(ncid,'c2_to_c0',c2_to_c0_id))
end if

! Read arrays
call ncerr(subr,nf90_get_var(ncid,c1_to_c0_id,hdata%c1_to_c0))
call ncerr(subr,nf90_get_var(ncid,c1l0_log_id,c1l0_logint))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (c1l0_logint(ic1,il0)==0) then
         hdata%c1l0_log(ic1,il0) = .false.
      else if (c1l0_logint(ic1,il0)==1) then
         hdata%c1l0_log(ic1,il0) = .true.
      end if
   end do
end do
call ncerr(subr,nf90_get_var(ncid,c1c3_to_c0_id,hdata%c1c3_to_c0))
call ncerr(subr,nf90_get_var(ncid,c1c3l0_log_id,c1c3l0_logint))
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (c1c3l0_logint(ic1,jc3,il0)==0) then
            hdata%c1c3l0_log(ic1,jc3,il0) = .false.
         else if (c1c3l0_logint(ic1,jc3,il0)==1) then
            hdata%c1c3l0_log(ic1,jc3,il0) = .true.
         end if
      end do
   end do
end do
if ((ios==0).and.(nam%local_diag.or.nam%displ_diag)) then
   call ncerr(subr,nf90_get_var(ncid,c2_to_c1_id,hdata%c2_to_c1))
   call ncerr(subr,nf90_get_var(ncid,c2_to_c0_id,hdata%c2_to_c0))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

! Read nearest neighbors and interpolation
if ((ios==0).and.(nam%local_diag.or.nam%displ_diag)) then
   ! Allocation
   allocate(local_maskint(nam%nc1,hdata%nc2))
   allocate(displ_maskint(nam%nc1,hdata%nc2))

   do il0i=1,geom%nl0i
      write(il0ichar,'(i3.3)') il0i
      info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc',nf90_nowrite,ncid)
      if (info/=nf90_noerr) then
         call msgwarning('cannot find nearest neighbors and interpolation data to read, recomputing HDIAG sampling')
         nam%sam_write = .true.
         ios = 3
         return
      end if
      call ncerr(subr,nf90_inq_varid(ncid,'local_mask',local_mask_id))
      call ncerr(subr,nf90_inq_varid(ncid,'displ_mask',displ_mask_id))
      call ncerr(subr,nf90_get_var(ncid,local_mask_id,local_maskint))
      call ncerr(subr,nf90_get_var(ncid,displ_mask_id,displ_maskint))
      do ic2=1,hdata%nc2
         do ic1=1,nam%nc1
            if (local_maskint(ic1,ic2)==1) then
               hdata%local_mask(ic1,ic2,il0i) = .true.
            elseif (local_maskint(ic1,ic2)==0) then
               hdata%local_mask(ic1,ic2,il0i) = .false.
            else
               call msgerror('wrong local_mask')
            end if
            if (displ_maskint(ic1,ic2)==1) then
               hdata%displ_mask(ic1,ic2,il0i) = .true.
            elseif (displ_maskint(ic1,ic2)==0) then
               hdata%displ_mask(ic1,ic2,il0i) = .false.
            else
               call msgerror('wrong displ_mask')
            end if
         end do
      end do
      info = nf90_inq_varid(ncid,'nn_c2_index',nn_c2_index_id)
      if (info==nf90_noerr) then
         call ncerr(subr,nf90_inq_varid(ncid,'nn_c2_dist',nn_c2_dist_id))
         call ncerr(subr,nf90_get_var(ncid,nn_c2_index_id,hdata%nn_c2_index(:,:,il0i)))
         call ncerr(subr,nf90_get_var(ncid,nn_c2_dist_id,hdata%nn_c2_dist(:,:,il0i)))
      else
         call msgwarning('cannot find nc2 nearest neighbors data to read, recomputing HDIAG sampling')
         nam%sam_write = .true.
         ios = 4
      end if
      write(hdata%hfull(il0i)%prefix,'(a,i3.3)') 'hfull_',il0i
      call hdata%hfull(il0i)%read(ncid)
      call ncerr(subr,nf90_close(ncid))
   end do

   ! Release memory
   deallocate(local_maskint)
   deallocate(displ_maskint)
end if

end subroutine hdata_read

!----------------------------------------------------------------------
! Subroutine: hdata_write
!> Purpose: write hdata object
!----------------------------------------------------------------------
subroutine hdata_write(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(in) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(in) :: geom    !< Geometry

! Local variables
integer :: il0,il0i,ic1,jc3,ic2
integer :: ncid,nl0_id,nc1_id,nc2_1_id,nc2_2_id,nc_id
integer :: lat_id,lon_id,smax_id,c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id,local_mask_id,displ_mask_id,nn_c2_index_id,nn_c2_dist_id
integer :: c1l0_logint(nam%nc1,geom%nl0),c1c3l0_logint(nam%nc1,nam%nc3,geom%nl0)
integer,allocatable :: local_maskint(:,:),displ_maskint(:,:)
real(kind_real) :: lon(nam%nc1,nam%nc3,geom%nl0),lat(nam%nc1,nam%nc3,geom%nl0)
character(len=3) :: il0ichar
character(len=1024) :: subr = 'hdata_write'

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Create file
write(mpl%unit,'(a7,a)') '','Write HDIAG sampling'
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0r',nam%nl0r))
call ncerr(subr,nf90_def_dim(ncid,'nc',nam%nc3,nc_id))
call ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
if (nam%local_diag.or.nam%displ_diag) call ncerr(subr,nf90_def_dim(ncid,'nc2_1',hdata%nc2,nc2_1_id))

! Define arrays
call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc1_id,nc_id,nl0_id/),lat_id))
call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc1_id,nc_id,nl0_id/),lon_id))
call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'smax',ncfloat,(/nc_id,nl0_id/),smax_id))
call ncerr(subr,nf90_put_att(ncid,smax_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'c1_to_c0',nf90_int,(/nc1_id/),c1_to_c0_id))
call ncerr(subr,nf90_put_att(ncid,c1_to_c0_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'c1l0_log',nf90_int,(/nc1_id,nl0_id/),c1l0_log_id))
call ncerr(subr,nf90_put_att(ncid,c1l0_log_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'c1c3_to_c0',nf90_int,(/nc1_id,nc_id/),c1c3_to_c0_id))
call ncerr(subr,nf90_put_att(ncid,c1c3_to_c0_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'c1c3l0_log',nf90_int,(/nc1_id,nc_id,nl0_id/),c1c3l0_log_id))
call ncerr(subr,nf90_put_att(ncid,c1c3l0_log_id,'_FillValue',msvali))
if (nam%local_diag.or.nam%displ_diag) then
   call ncerr(subr,nf90_def_var(ncid,'c2_to_c1',nf90_int,(/nc2_1_id/),c2_to_c1_id))
   call ncerr(subr,nf90_put_att(ncid,c2_to_c1_id,'_FillValue',msvali))
   call ncerr(subr,nf90_def_var(ncid,'c2_to_c0',nf90_int,(/nc2_1_id/),c2_to_c0_id))
   call ncerr(subr,nf90_put_att(ncid,c2_to_c0_id,'_FillValue',msvali))
end if

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write arrays
call msr(lon)
call msr(lat)
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (hdata%c1c3l0_log(ic1,jc3,il0)) then
            lon(ic1,jc3,il0) = geom%lon(hdata%c1c3_to_c0(ic1,jc3))*rad2deg
            lat(ic1,jc3,il0) = geom%lat(hdata%c1c3_to_c0(ic1,jc3))*rad2deg
         end if
      end do
   end do
end do
call ncerr(subr,nf90_put_var(ncid,lon_id,lon))
call ncerr(subr,nf90_put_var(ncid,lat_id,lat))
call ncerr(subr,nf90_put_var(ncid,smax_id,float(count(hdata%c1c3l0_log,dim=1))))
call ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,hdata%c1_to_c0))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (hdata%c1l0_log(ic1,il0)) then
         c1l0_logint(ic1,il0) = 1
      else
         c1l0_logint(ic1,il0) = 0
      end if
   end do
end do
call ncerr(subr,nf90_put_var(ncid,c1l0_log_id,c1l0_logint))
call ncerr(subr,nf90_put_var(ncid,c1c3_to_c0_id,hdata%c1c3_to_c0))
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (hdata%c1c3l0_log(ic1,jc3,il0)) then
            c1c3l0_logint(ic1,jc3,il0) = 1
         else
            c1c3l0_logint(ic1,jc3,il0) = 0
         end if
      end do
   end do
end do
call ncerr(subr,nf90_put_var(ncid,c1c3l0_log_id,c1c3l0_logint))
if (nam%local_diag.or.nam%displ_diag) then
   call ncerr(subr,nf90_put_var(ncid,c2_to_c1_id,hdata%c2_to_c1))
   call ncerr(subr,nf90_put_var(ncid,c2_to_c0_id,hdata%c2_to_c0))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

! Write nearest neighbors and interpolation
if (nam%local_diag.or.nam%displ_diag) then
   ! Allocation
   allocate(local_maskint(nam%nc1,hdata%nc2))
   allocate(displ_maskint(nam%nc1,hdata%nc2))

   do il0i=1,geom%nl0i
      write(il0ichar,'(i3.3)') il0i
      call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc', &
    & or(nf90_clobber,nf90_64bit_offset),ncid))
      call ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
      call ncerr(subr,nf90_def_dim(ncid,'nc2_1',hdata%nc2,nc2_1_id))
      call ncerr(subr,nf90_def_var(ncid,'local_mask',nf90_int,(/nc1_id,nc2_1_id/),local_mask_id))
      call ncerr(subr,nf90_def_var(ncid,'displ_mask',nf90_int,(/nc1_id,nc2_1_id/),displ_mask_id))
      call ncerr(subr,nf90_put_att(ncid,local_mask_id,'_FillValue',msvali))
      call ncerr(subr,nf90_put_att(ncid,displ_mask_id,'_FillValue',msvali))
      call ncerr(subr,nf90_def_dim(ncid,'nc2_2',hdata%nc2,nc2_2_id))
      call ncerr(subr,nf90_def_var(ncid,'nn_c2_index',nf90_int,(/nc2_1_id,nc2_2_id/),nn_c2_index_id))
      call ncerr(subr,nf90_put_att(ncid,nn_c2_index_id,'_FillValue',msvali))
      call ncerr(subr,nf90_def_var(ncid,'nn_c2_dist',ncfloat,(/nc2_1_id,nc2_2_id/),nn_c2_dist_id))
      call ncerr(subr,nf90_put_att(ncid,nn_c2_dist_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
      do ic2=1,hdata%nc2
         do ic1=1,nam%nc1
            if (hdata%local_mask(ic1,ic2,il0i)) then
               local_maskint(ic1,ic2) = 1
            else
               local_maskint(ic1,ic2) = 0
            end if
            if (hdata%displ_mask(ic1,ic2,il0i)) then
               displ_maskint(ic1,ic2) = 1
            else
               displ_maskint(ic1,ic2) = 0
            end if
         end do
      end do
      call ncerr(subr,nf90_put_var(ncid,local_mask_id,local_maskint))
      call ncerr(subr,nf90_put_var(ncid,displ_mask_id,displ_maskint))
      call ncerr(subr,nf90_put_var(ncid,nn_c2_index_id,hdata%nn_c2_index(:,:,il0i)))
      call ncerr(subr,nf90_put_var(ncid,nn_c2_dist_id,hdata%nn_c2_dist(:,:,il0i)))
      call hdata%hfull(il0i)%write(ncid)
      call ncerr(subr,nf90_close(ncid))
   end do

   ! Release memory
   deallocate(local_maskint)
   deallocate(displ_maskint)
end if

end subroutine hdata_write

!----------------------------------------------------------------------
! Subroutine: hdata_setup_sampling
!> Purpose: setup sampling
!----------------------------------------------------------------------
subroutine hdata_setup_sampling(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(inout) :: nam      !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: info,ic0,il0,ic1,ic2,ildw,jc3,il0i,jc1,kc1
integer :: mask_ind(nam%nc1)
integer,allocatable :: vbot(:),vtop(:),nn_c1_index(:)
real(kind_real) :: rh0(geom%nc0,geom%nl0),dum(1)
real(kind_real),allocatable :: nn_c1_dist(:)
type(ctree_type) :: ctree_diag
type(linop_type) :: hbase

! Check subsampling size
if (nam%nc1>maxval(count(geom%mask,dim=1))) then
   call msgwarning('nc1 is too large for then mask, reset nc1 to the largest possible value')
   nam%nc1 = maxval(count(geom%mask,dim=1))
end if

! Define nc2
if (nam%new_lct) then
   hdata%nc2 = nam%nc1
elseif (nam%local_diag) then
   hdata%nc2 = int(2.0*maxval(geom%area)/(sqrt(3.0)*(nam%local_rad)**2))
   write(mpl%unit,'(a7,a,i8)') '','Estimated nc2 from local diagnostic radius: ',hdata%nc2
   hdata%nc2 = min(hdata%nc2,nam%nc1)
   write(mpl%unit,'(a7,a,i8)') '','Final nc2: ',hdata%nc2
elseif (nam%displ_diag) then
   hdata%nc2 = nam%nc1
end if

! Allocation
call hdata%alloc(nam,geom)

! Read or compute sampling data
info = 1
if (nam%sam_read) call hdata%read(nam,geom,info)
if (info==1) then
   ! Compute zero-separation sampling
   call hdata%compute_sampling_zs(nam,geom)

   if (nam%new_lct) then
      ! Compute LCT sampling
      call hdata%compute_sampling_lct(nam,geom)
   else
      ! Compute positive separation sampling
      call hdata%compute_sampling_ps(nam,geom)
   end if

   ! Compute sampling mask
   call hdata%compute_sampling_mask(nam,geom)
end if

if (nam%local_diag.or.nam%displ_diag) then
   if ((info==1).or.(info==2)) then
      ! Define subsampling
      write(mpl%unit,'(a7,a)') '','Define subsampling'
      mask_ind = 1
      rh0 = 1.0
      if (mpl%main) call rng%initialize_sampling(nam%nc1,dble(geom%lon(hdata%c1_to_c0)),dble(geom%lat(hdata%c1_to_c0)),mask_ind, &
    & rh0,nam%ntry,nam%nrep,hdata%nc2,hdata%c2_to_c1)
      call mpl%bcast(hdata%c2_to_c1,mpl%ioproc)
      hdata%c2_to_c0 = hdata%c1_to_c0(hdata%c2_to_c1)
   end if

   if ((info==1).or.(info==2).or.(info==3).or.(info==4)) then
      ! Create cover trees
      write(mpl%unit,'(a7,a)') '','Create cover trees'
      do il0=1,geom%nl0
         if ((il0==1).or.(geom%nl0i>1)) then
            write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(il0)
            call ctree_diag%create(hdata%nc2,geom%lon(hdata%c2_to_c0),geom%lat(hdata%c2_to_c0),hdata%c1l0_log(hdata%c2_to_c1,il0))
            do ic2=1,hdata%nc2
               ic1 = hdata%c2_to_c1(ic2)
               ic0 = hdata%c2_to_c0(ic2)
               if (hdata%c1l0_log(ic1,il0)) call ctree_diag%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0), &
                & hdata%nc2,hdata%nn_c2_index(:,ic2,il0),hdata%nn_c2_dist(:,ic2,il0))
            end do
            call ctree_diag%delete
         end if
      end do
   end if

   ! Compute sampling mesh
   write(mpl%unit,'(a7,a)') '','Compute sampling mesh'
   call hdata%mesh%create(hdata%nc2,geom%lon(hdata%c2_to_c0),geom%lat(hdata%c2_to_c0),.false.)

   ! Compute triangles list
   write(mpl%unit,'(a7,a)') '','Compute triangles list '
   call hdata%mesh%trlist

   ! Find boundary nodes
   write(mpl%unit,'(a7,a)') '','Find boundary nodes'
   call hdata%mesh%bnodes

   if ((info==1).or.(info==2).or.(info==3)) then
      ! Allocation
      allocate(nn_c1_index(nam%nc1))
      allocate(nn_c1_dist(nam%nc1))
      allocate(vbot(hdata%nc2))
      allocate(vtop(hdata%nc2))

      ! Compute nearest neighbors
      write(mpl%unit,'(a7,a)') '','Compute nearest neighbors'
      do il0i=1,geom%nl0i
         write(mpl%unit,'(a10,a,i3)') '','Independent level ',il0i
         call ctree_diag%create(nam%nc1,geom%lon(hdata%c1_to_c0),geom%lat(hdata%c1_to_c0),hdata%c1l0_log(:,il0i))
         do ic2=1,hdata%nc2
            ic1 = hdata%c2_to_c1(ic2)
            ic0 = hdata%c2_to_c0(ic2)
            if (hdata%c1l0_log(ic1,il0i)) then
               ! Find nearest neighbors
               call ctree_diag%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nam%nc1,nn_c1_index,nn_c1_dist)

               do jc1=1,nam%nc1
                  kc1 = nn_c1_index(jc1)
                  hdata%local_mask(kc1,ic2,il0i) = (jc1==1).or.(nn_c1_dist(jc1)<min(nam%local_rad,hdata%mesh%bdist(ic2)))
                  hdata%displ_mask(kc1,ic2,il0i) = (jc1==1).or.(nn_c1_dist(jc1)<min(nam%displ_rad,hdata%mesh%bdist(ic2)))
               end do
            end if
         end do
         call ctree_diag%delete
      end do

      ! Initialize vbot and vtop
      vbot = 1
      vtop = geom%nl0

      do il0i=1,geom%nl0i
         ! Compute grid interpolation
         write(hdata%hfull(il0i)%prefix,'(a,i3.3)') 'hfull_',il0i
         call hdata%hfull(il0i)%interp(geom,il0i,hdata%nc2,hdata%c2_to_c0,nam%mask_check,vbot,vtop,nam%diag_interp,hbase)
      end do

      ! Release memory
      deallocate(nn_c1_index)
      deallocate(nn_c1_dist)
      deallocate(vbot)
      deallocate(vtop)
   end if
end if

! Write sampling data
if (nam%sam_write.and.mpl%main) call hdata%write(nam,geom)

! Compute nearest neighbors for local diagnostics output
if (nam%local_diag.and.(nam%nldwv>0)) then
   write(mpl%unit,'(a7,a)') '','Compute nearest neighbors for local diagnostics output'
   allocate(hdata%nn_ldwv_index(nam%nldwv))
   call ctree_diag%create(hdata%nc2,geom%lon(hdata%c2_to_c0), &
                geom%lat(hdata%c2_to_c0),hdata%c1l0_log(hdata%c2_to_c1,1))
   do ildw=1,nam%nldwv
      call ctree_diag%find_nearest_neighbors(nam%lon_ldwv(ildw)*deg2rad,nam%lat_ldwv(ildw)*deg2rad, &
    & 1,hdata%nn_ldwv_index(ildw:ildw),dum)
   end do
   call ctree_diag%delete
end if

! Print results
write(mpl%unit,'(a7,a)') '','Sampling efficiency (%):'
do il0=1,geom%nl0
   write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),' ~> '
   do jc3=1,nam%nc3
      if (count(hdata%c1c3l0_log(:,jc3,il0))>=nam%nc1/2) then
         ! Sucessful sampling
         write(mpl%unit,'(a,i3,a)',advance='no') trim(green), &
       & int(100.0*float(count(hdata%c1c3l0_log(:,jc3,il0)))/float(nam%nc1)),trim(black)
      else
         ! Insufficient sampling
         write(mpl%unit,'(a,i3,a)',advance='no') trim(peach), &
       & int(100.0*float(count(hdata%c1c3l0_log(:,jc3,il0)))/float(nam%nc1)),trim(black)
      end if
      if (jc3<nam%nc3) write(mpl%unit,'(a)',advance='no') '-'
   end do
   write(mpl%unit,'(a)') ' '
end do

end subroutine hdata_setup_sampling

!----------------------------------------------------------------------
! Subroutine: hdata_compute_sampling_zs
!> Purpose: compute zero-separation sampling
!----------------------------------------------------------------------
subroutine hdata_compute_sampling_zs(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: ic0,ic1,fac,np,ip
integer :: mask_ind_col(geom%nc0),nn_index(1)
real(kind_real) :: rh0(geom%nc0),dum(1)
real(kind_real),allocatable :: lon(:),lat(:)
character(len=5) :: ic1char

! Initialize mask
mask_ind_col = 0
do ic0=1,geom%nc0
   if (any(geom%mask(ic0,:))) mask_ind_col(ic0) = 1
end do

! Initialize support radius to 1.0
rh0 = 1.0

! Compute subset
write(mpl%unit,'(a7,a)') '','Compute horizontal subset C1'
if (nam%nc1<maxval(count(geom%mask,dim=1))) then
   if (mpl%main) then
      if (.true.) then
         ! Random draw
         call rng%initialize_sampling(geom%nc0,dble(geom%lon),dble(geom%lat),mask_ind_col,rh0,nam%ntry,nam%nrep, &
       & nam%nc1,hdata%c1_to_c0)
      else
         ! Compute icosahedron size
         call closest_icos(nam%nc1,fac,np)

         ! Allocation
         allocate(lon(np))
         allocate(lat(np))

         ! Compute icosahedron
         call build_icos(fac,np,lon,lat)

         ! Fill c1_to_c0
         ic1 = 0
         do ip=1,np
            ! Find nearest neighbor
            call geom%ctree%find_nearest_neighbors(lon(ip),lat(ip),1,nn_index,dum)
            ic0 = nn_index(1)

            ! Check mask
            if (mask_ind_col(ic0)==1) then
               ic1 = ic1+1
               if (ic1<=nam%nc1) hdata%c1_to_c0(ic1) = ic0
            end if
         end do

         ! Check size
         if (ic1<nam%nc1) then
            write(ic1char,'(i5)') ic1
            call msgerror('nc1 should be decreased to '//ic1char)
         end if
         if (ic1>nam%nc1) then
            write(ic1char,'(i5)') ic1
            call msgwarning('nc1 could be increased to '//ic1char)
         end if

         ! Release memory
         deallocate(lon)
         deallocate(lat)
      end if
   end if
   call mpl%bcast(hdata%c1_to_c0,mpl%ioproc)
else
   ic1 = 0
   do ic0=1,geom%nc0
      if (any(geom%mask(ic0,:))) then
         ic1 = ic1+1
         hdata%c1_to_c0(ic1) = ic0
      end if
   end do
end if

end subroutine hdata_compute_sampling_zs

!----------------------------------------------------------------------
! Subroutine: hdata_compute_sampling_ps
!> Purpose: compute positive separation sampling
!----------------------------------------------------------------------
subroutine hdata_compute_sampling_ps(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: irmaxloc,progint,jc3,ic1,ir,ic0,jc0,i,nvc0,ivc0,icinf,icsup,ictest
integer,allocatable :: vic0(:)
real(kind_real) :: d
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical :: found,done(nam%nc3*nam%nc1)

! First class
hdata%c1c3_to_c0(:,1) = hdata%c1_to_c0

if (nam%nc3>1) then
   write(mpl%unit,'(a7,a)',advance='no') '','Compute positive separation sampling: '

   ! Initialize
   do jc3=1,nam%nc3
      if (jc3/=1) call msi(hdata%c1c3_to_c0(:,jc3))
   end do

   ! Define valid nodes vector
   nvc0 = count(any(geom%mask,dim=2))
   allocate(vic0(nvc0))
   ivc0 = 0
   do ic0=1,geom%nc0
      if (any(geom%mask(ic0,:))) then
         ivc0 = ivc0+1
         vic0(ivc0) = ic0
      end if
   end do

   ! Sample classes of positive separation
   call prog_init(progint)
   ir = 0
   irmaxloc = irmax
   do while ((.not.all(isnotmsi(hdata%c1c3_to_c0))).and.(nvc0>1).and.(ir<=irmaxloc))
      ! Try a random point
      if (mpl%main) call rng%rand_integer(1,nvc0,i)
      call mpl%bcast(i,mpl%ioproc)
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
         if (isnotmsi(hdata%c1_to_c0(ic1))) then
            ! Compute the distance
            ic0 = hdata%c1_to_c0(ic1)
            call sphere_dist(geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),d)

            ! Find the class (dichotomy method)
            if ((d>0.0).and.(d<(float(nam%nc3)-0.5)*nam%dc)) then
               jc3 = 1
               icinf = 1
               icsup = nam%nc3
               found = .false.
               do while (.not.found)
                  ! New value
                  ictest = (icsup+icinf)/2

                  ! Update
                  if (d<(float(ictest)-0.5)*nam%dc) icsup = ictest
                  if (d>(float(ictest)-0.5)*nam%dc) icinf = ictest

                  ! Exit test
                  if (icsup==icinf+1) then
                     if (abs((float(icinf)-0.5)*nam%dc-d)<abs((float(icsup)-0.5)*nam%dc-d)) then
                        jc3 = icinf
                     else
                        jc3 = icsup
                     end if
                     found = .true.
                  end if
               end do

               ! Find if this class has not been aready filled
               if ((jc3/=1).and.(.not.isnotmsi(hdata%c1c3_to_c0(ic1,jc3)))) hdata%c1c3_to_c0(ic1,jc3) = jc0
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

      ! Print progression
      done = pack(isnotmsi(hdata%c1c3_to_c0),mask=.true.)
      call prog_print(progint,done)
   end do
   write(mpl%unit,'(a)') '100%'

   ! Release memory
   deallocate(vic0)
end if

end subroutine hdata_compute_sampling_ps

!----------------------------------------------------------------------
! Subroutine: hdata_compute_sampling_lct
!> Purpose: compute LCT sampling
!----------------------------------------------------------------------
subroutine hdata_compute_sampling_lct(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: i,il0,ic1,ic0,jc0,ibnd,ic3,progint
integer :: nn(nam%nc3)
integer :: iproc,ic1_s(mpl%nproc),ic1_e(mpl%nproc),nc1_loc(mpl%nproc),ic1_loc
integer,allocatable :: sbufi(:),rbufi(:)
real(kind_real) :: dum(nam%nc3)
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: sbufl(:),rbufl(:),done(:)

write(mpl%unit,'(a7,a)',advance='no') '','Compute LCT sampling: '

! MPI splitting
call mpl%split(nam%nc1,ic1_s,ic1_e,nc1_loc)

! Allocation
allocate(done(nc1_loc(mpl%myproc)))

! Initialization
call prog_init(progint)

do ic1_loc=1,nc1_loc(mpl%myproc)
   ! MPI offset
   ic1 = ic1_s(mpl%myproc)+ic1_loc-1

   ! Check location validity
   if (isnotmsi(hdata%c1_to_c0(ic1))) then
      ! Find neighbors
      call geom%ctree%find_nearest_neighbors(dble(geom%lon(hdata%c1_to_c0(ic1))),dble(geom%lat(hdata%c1_to_c0(ic1))), &
    & nam%nc3,nn,dum)

      ! Copy neighbor index
      do ic3=1,nam%nc3
         jc0 = nn(ic3)
         hdata%c1c3_to_c0(ic1,ic3) = nn(ic3)
         do il0=1,geom%nl0
            hdata%c1c3l0_log(ic1,ic3,il0) = geom%mask(jc0,il0)
         end do
      end do

      if (nam%mask_check) then
         ! Check that great circle to neighbors is not crossing mask boundaries
         do il0=1,geom%nl0
            !$omp parallel do schedule(static) private(ic3,ic0,jc0) firstprivate(x,y,z,v1,v2,va,vp,t)
            do ic3=1,nam%nc3
               ! Allocation
               allocate(x(2))
               allocate(y(2))
               allocate(z(2))
               allocate(v1(3))
               allocate(v2(3))
               allocate(va(3))
               allocate(vp(3))
               allocate(t(4))

               ! Indices
               ic0 = hdata%c1_to_c0(ic1)
               jc0 = hdata%c1c3_to_c0(ic1,ic3)

               ! Transform to cartesian coordinates
               call trans(2,geom%lat((/ic0,jc0/)),geom%lon((/ic0,jc0/)),x,y,z)

               ! Compute arc orthogonal vector
               v1 = (/x(1),y(1),z(1)/)
               v2 = (/x(2),y(2),z(2)/)
               call vector_product(v1,v2,va)

               ! Check if arc is crossing boundary arcs
               do ibnd=1,geom%nbnd(il0)
                  call vector_product(va,geom%vbnd(:,ibnd,il0),vp)
                  v1 = (/x(1),y(1),z(1)/)
                  call vector_triple_product(v1,va,vp,t(1))
                  v1 = (/x(2),y(2),z(2)/)
                  call vector_triple_product(v1,va,vp,t(2))
                  v1 = (/geom%xbnd(1,ibnd,il0),geom%ybnd(1,ibnd,il0),geom%zbnd(1,ibnd,il0)/)
                  call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(3))
                  v1 = (/geom%xbnd(2,ibnd,il0),geom%ybnd(2,ibnd,il0),geom%zbnd(2,ibnd,il0)/)
                  call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(4))
                  t(1) = -t(1)
                  t(3) = -t(3)
                  if (all(t>0).or.(all(t<0))) then
                     hdata%c1c3l0_log(ic1,ic3,il0) = .false.
                     exit
                  end if
               end do

               ! Memory release
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
         end do
      end if
   end if

   ! Print progression
   done(ic1_loc) = .true.
   call prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc/=mpl%ioproc) then
         ! Allocation
         allocate(rbufi(nc1_loc(iproc)*nam%nc3))
         allocate(rbufl(nc1_loc(iproc)*nam%nc3*geom%nl0))

         ! Receive data on ioproc
         call mpl%recv(nc1_loc(iproc)*nam%nc3,rbufi,iproc,mpl%tag)
         call mpl%recv(nc1_loc(iproc)*nam%nc3*geom%nl0,rbufl,iproc,mpl%tag+1)

         ! Format data
         i = 0
         do ic3=1,nam%nc3
            do ic1_loc=1,nc1_loc(iproc)
               i = i+1
               ic1 = ic1_s(iproc)+ic1_loc-1
               hdata%c1c3_to_c0(ic1,ic3) = rbufi(i)
            end do
         end do
         i = 0
         do il0=1,geom%nl0
            do ic3=1,nam%nc3
               do ic1_loc=1,nc1_loc(iproc)
                  i = i+1
                  ic1 = ic1_s(iproc)+ic1_loc-1
                  hdata%c1c3l0_log(ic1,ic3,il0) = rbufl(i)
               end do
            end do
         end do

         ! Release memory
         deallocate(rbufi)
         deallocate(rbufl)
      end if
   end do
else
   ! Allocation
   allocate(sbufi(nc1_loc(mpl%myproc)*nam%nc3))
   allocate(sbufl(nc1_loc(mpl%myproc)*nam%nc3*geom%nl0))

   ! Prepare buffers
   i = 0
   do ic3=1,nam%nc3
      do ic1_loc=1,nc1_loc(mpl%myproc)
         i = i+1
         ic1 = ic1_s(mpl%myproc)+ic1_loc-1
         sbufi(i) = hdata%c1c3_to_c0(ic1,ic3)
      end do
   end do
   i = 0
   do il0=1,geom%nl0
      do ic3=1,nam%nc3
         do ic1_loc=1,nc1_loc(mpl%myproc)
            i = i+1
            ic1 = ic1_s(mpl%myproc)+ic1_loc-1
            sbufl(i) = hdata%c1c3l0_log(ic1,ic3,il0)
         end do
      end do
   end do

   ! Send data to ioproc
   call mpl%send(nc1_loc(mpl%myproc)*nam%nc3,sbufi,mpl%ioproc,mpl%tag)
   call mpl%send(nc1_loc(mpl%myproc)*nam%nc3*geom%nl0,sbufl,mpl%ioproc,mpl%tag+1)

   ! Release memory
   deallocate(sbufi)
   deallocate(sbufl)
end if
mpl%tag = mpl%tag+2

! Broadcast data
call mpl%bcast(hdata%c1c3_to_c0,mpl%ioproc)
call mpl%bcast(hdata%c1c3l0_log,mpl%ioproc)

end subroutine hdata_compute_sampling_lct

!----------------------------------------------------------------------
! Subroutine: hdata_compute_sampling_mask
!> Purpose: compute sampling mask
!----------------------------------------------------------------------
subroutine hdata_compute_sampling_mask(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: jc3,ic1,ic0,jc0,il0
logical :: valid

! First point
do il0=1,geom%nl0
   hdata%c1l0_log(:,il0) = geom%mask(hdata%c1_to_c0,il0)
end do

! Second point
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         ! Indices
         ic0 = hdata%c1c3_to_c0(ic1,jc3)
         jc0 = hdata%c1c3_to_c0(ic1,1)

         ! Check point index
         valid = isnotmsi(ic0).and.isnotmsi(jc0)

         if (valid) then
            ! Check mask
            valid = geom%mask(ic0,il0).and.geom%mask(jc0,il0)

            ! Check mask bounds
            if (nam%mask_check.and.valid) call geom%check_arc(il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),valid)
         end if
         hdata%c1c3l0_log(ic1,jc3,il0) = valid
      end do
   end do
end do

end subroutine hdata_compute_sampling_mask

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_a
!> Purpose: compute HDIAG MPI distribution, halo A
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_a(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: ic0a,ic0,ic1a,ic1

! Allocation
allocate(hdata%lcheck_c0a(geom%nc0))
allocate(hdata%lcheck_c1a(nam%nc1))

! Halo definitions

! Halo A
hdata%lcheck_c0a = .false.
hdata%lcheck_c1a = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   if (geom%c0_to_proc(ic0)==mpl%myproc) hdata%lcheck_c0a(ic0) = .true.
end do
do ic1=1,nam%nc1
   ic0 = hdata%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) hdata%lcheck_c1a(ic1) = .true.
end do

! Halo sizes
hdata%nc1a = count(hdata%lcheck_c1a)

! Global <-> local conversions for fields

! Halo A
allocate(hdata%c1a_to_c1(hdata%nc1a))
allocate(hdata%c1_to_c1a(nam%nc1))
ic1a = 0
do ic1=1,nam%nc1
   if (hdata%lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      hdata%c1a_to_c1(ic1a) = ic1
      hdata%c1_to_c1a(ic1) = ic1a
   end if
end do

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0 =        ',geom%nc0
write(mpl%unit,'(a10,a,i8)') '','nc0a =       ',geom%nc0a
write(mpl%unit,'(a10,a,i8)') '','nl0 =        ',geom%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =        ',nam%nc1
write(mpl%unit,'(a10,a,i8)') '','nc1a =       ',hdata%nc1a

end subroutine hdata_compute_mpi_a

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_ab
!> Purpose: compute HDIAG MPI distribution, halos A-B
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_ab(hdata,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: iproc,ic0,ic2a,ic2b,ic2,jc2,i_s,i_s_loc,h_n_s_max,il0i,h_n_s_max_loc,nc2a,nc2b
integer,allocatable :: interph_lg(:,:)
integer,allocatable :: c2a_to_c2(:),c2b_to_c2(:),c2a_to_c2b(:)
type(com_type) :: com_AB(mpl%nproc)

! Allocation
h_n_s_max = 0
do il0i=1,geom%nl0i
   h_n_s_max = max(h_n_s_max,hdata%hfull(il0i)%n_s)
end do
allocate(hdata%c2_to_proc(hdata%nc2))
allocate(hdata%proc_to_nc2a(mpl%nproc))
allocate(hdata%h(geom%nl0i))
allocate(hdata%lcheck_c2a(hdata%nc2))
allocate(hdata%lcheck_c2b(hdata%nc2))
allocate(hdata%lcheck_h(h_n_s_max,geom%nl0i))

! Halo definitions

! Halo A
hdata%lcheck_c2a = .false.
do ic2=1,hdata%nc2
   ic0 = hdata%c2_to_c0(ic2)
   if (geom%c0_to_proc(ic0)==mpl%myproc) hdata%lcheck_c2a(ic2) = .true.
end do

! Halo B
hdata%lcheck_h = .false.
hdata%lcheck_c2b = .false.
do il0i=1,geom%nl0i
   do i_s=1,hdata%hfull(il0i)%n_s
      ic0 = hdata%hfull(il0i)%row(i_s)
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%myproc) then
         jc2 = hdata%hfull(il0i)%col(i_s)
         hdata%lcheck_h(i_s,il0i) = .true.
         hdata%lcheck_c2b(jc2) = .true.
      end if
   end do
end do

! Halo sizes
hdata%nc2a = count(hdata%lcheck_c2a)
do il0i=1,geom%nl0i
   hdata%h(il0i)%n_s = count(hdata%lcheck_h(:,il0i))
end do
hdata%nc2b = count(hdata%lcheck_c2b)

! Global <-> local conversions for fields

allocate(hdata%c2a_to_c2(hdata%nc2a))
allocate(hdata%c2_to_c2a(hdata%nc2))
ic2a = 0
do ic2=1,hdata%nc2
   if (hdata%lcheck_c2a(ic2)) then
      ic2a = ic2a+1
      hdata%c2a_to_c2(ic2a) = ic2
   end if
end do

! Halo B
allocate(hdata%c2b_to_c2(hdata%nc2b))
allocate(hdata%c2_to_c2b(hdata%nc2))
ic2b = 0
do ic2=1,hdata%nc2
   if (hdata%lcheck_c2b(ic2)) then
      ic2b = ic2b+1
      hdata%c2b_to_c2(ic2b) = ic2
      hdata%c2_to_c2b(ic2) = ic2b
   end if
end do

! Inter-halo conversions
allocate(hdata%c2a_to_c2b(hdata%nc2a))
do ic2a=1,hdata%nc2a
   ic2 = hdata%c2a_to_c2(ic2a)
   ic2b = hdata%c2_to_c2b(ic2)
   hdata%c2a_to_c2b(ic2a) = ic2b
end do

! Global <-> local conversions for data
h_n_s_max_loc = 0
do il0i=1,geom%nl0i
   h_n_s_max_loc = max(h_n_s_max_loc,hdata%h(il0i)%n_s)
end do
allocate(interph_lg(h_n_s_max_loc,geom%nl0i))
do il0i=1,geom%nl0i
   i_s_loc = 0
   do i_s=1,hdata%hfull(il0i)%n_s
      if (hdata%lcheck_h(i_s,il0i)) then
         i_s_loc = i_s_loc+1
         interph_lg(i_s_loc,il0i) = i_s
      end if
   end do
end do

! Local data

! Horizontal interpolation
do il0i=1,geom%nl0i
   write(hdata%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
   hdata%h(il0i)%n_src = hdata%nc2b
   hdata%h(il0i)%n_dst = geom%nc0a
   call hdata%h(il0i)%alloc
   do i_s_loc=1,hdata%h(il0i)%n_s
      i_s = interph_lg(i_s_loc,il0i)
      hdata%h(il0i)%row(i_s_loc) = geom%c0_to_c0a(hdata%hfull(il0i)%row(i_s))
      hdata%h(il0i)%col(i_s_loc) = hdata%c2_to_c2b(hdata%hfull(il0i)%col(i_s))
      hdata%h(il0i)%S(i_s_loc) = hdata%hfull(il0i)%S(i_s)
   end do
   call hdata%h(il0i)%reorder
end do

! Get global distribution of the subgrid on ioproc
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         nc2a = hdata%nc2a
      else
         ! Receive dimension on ioproc
         call mpl%recv(nc2a,iproc,mpl%tag)
      end if

      ! Allocation
      allocate(c2a_to_c2(nc2a))

      if (iproc==mpl%ioproc) then
         ! Copy data
         c2a_to_c2 = hdata%c2a_to_c2
      else
         ! Receive data on ioproc
         call mpl%recv(nc2a,c2a_to_c2,iproc,mpl%tag+1)
      end if

      ! Fill c2_to_c2a
      do ic2a=1,nc2a
         hdata%c2_to_c2a(c2a_to_c2(ic2a)) = ic2a
      end do

      ! Release memory
      deallocate(c2a_to_c2)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(hdata%nc2a,mpl%ioproc,mpl%tag)

   ! Send data to ioproc
   call mpl%send(hdata%nc2a,hdata%c2a_to_c2,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc2a = hdata%nc2a
         nc2b = hdata%nc2b
      else
         ! Receive dimensions on ioproc
         call mpl%recv(nc2a,iproc,mpl%tag)
         call mpl%recv(nc2b,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(c2b_to_c2(nc2b))
      allocate(c2a_to_c2b(nc2a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c2b_to_c2 = hdata%c2b_to_c2
         c2a_to_c2b = hdata%c2a_to_c2b
      else
         ! Receive data on ioproc
         call mpl%recv(nc2b,c2b_to_c2,iproc,mpl%tag+2)
         call mpl%recv(nc2a,c2a_to_c2b,iproc,mpl%tag+3)
      end if

      ! Allocation
      com_AB(iproc)%nred = nc2a
      com_AB(iproc)%next = nc2b
      allocate(com_AB(iproc)%ext_to_proc(com_AB(iproc)%next))
      allocate(com_AB(iproc)%ext_to_red(com_AB(iproc)%next))
      allocate(com_AB(iproc)%red_to_ext(com_AB(iproc)%nred))

      ! AB communication
      do ic2b=1,nc2b
         ic2 = c2b_to_c2(ic2b)
         ic0 = hdata%c2_to_c0(ic2)
         com_AB(iproc)%ext_to_proc(ic2b) = geom%c0_to_proc(ic0)
         ic2a = hdata%c2_to_c2a(ic2)
         com_AB(iproc)%ext_to_red(ic2b) = ic2a
      end do
      com_AB(iproc)%red_to_ext = c2a_to_c2b

      ! Release memory
      deallocate(c2b_to_c2)
      deallocate(c2a_to_c2b)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(hdata%nc2a,mpl%ioproc,mpl%tag)
   call mpl%send(hdata%nc2b,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl%send(hdata%nc2b,hdata%c2b_to_c2,mpl%ioproc,mpl%tag+2)
   call mpl%send(hdata%nc2a,hdata%c2a_to_c2b,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4
call hdata%com_AB%setup(com_AB,'com_AB')

! MPI splitting
do ic2=1,hdata%nc2
   ic0 = hdata%c2_to_c0(ic2)
   hdata%c2_to_proc(ic2) = geom%c0_to_proc(ic0)
end do
do iproc=1,mpl%nproc
   hdata%proc_to_nc2a(iproc) = count(hdata%c2_to_proc==iproc)
end do

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc2 =        ',hdata%nc2
write(mpl%unit,'(a10,a,i8)') '','nc2a =       ',hdata%nc2a
write(mpl%unit,'(a10,a,i8)') '','nc2b =       ',hdata%nc2b
do il0i=1,geom%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',hdata%h(il0i)%n_s
end do

end subroutine hdata_compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_d
!> Purpose: compute HDIAG MPI distribution, halo D
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_d(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: iproc,ic0,ic0a,ic0d,ic1,ic2,ic2a,jc0,jc1
integer :: nc0a,nc0d
integer,allocatable :: c0d_to_c0(:),c0a_to_c0d(:)
type(com_type) :: com_AD(mpl%nproc)

! Allocation
allocate(hdata%lcheck_c0d(geom%nc0))

! Halo definitions

! Halo D
hdata%lcheck_c0d = hdata%lcheck_c0a
do ic2a=1,hdata%nc2a
   ic2 = hdata%c2a_to_c2(ic2a)
   ic1 = hdata%c2_to_c1(ic2)
   if (any(hdata%c1l0_log(ic1,:))) then
      do jc1=1,nam%nc1
         if (any(hdata%displ_mask(jc1,ic2,:))) then
            jc0 = hdata%c1_to_c0(jc1)
            hdata%lcheck_c0d(jc0) = .true.
         end if
      end do
   end if
end do
hdata%nc0d = count(hdata%lcheck_c0d)

! Global <-> local conversions for fields

! Halo D
allocate(hdata%c0d_to_c0(hdata%nc0d))
allocate(hdata%c0_to_c0d(geom%nc0))
call msi(hdata%c0_to_c0d)
ic0d = 0
do ic0=1,geom%nc0
   if (hdata%lcheck_c0d(ic0)) then
      ic0d = ic0d+1
      hdata%c0d_to_c0(ic0d) = ic0
      hdata%c0_to_c0d(ic0) = ic0d
   end if
end do

! Inter-halo conversions
allocate(hdata%c0a_to_c0d(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0d = hdata%c0_to_c0d(ic0)
   hdata%c0a_to_c0d(ic0a) = ic0d
end do

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = geom%nc0a
         nc0d = hdata%nc0d
      else
         ! Receive dimensions on ioproc
         call mpl%recv(nc0a,iproc,mpl%tag)
         call mpl%recv(nc0d,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(c0d_to_c0(nc0d))
      allocate(c0a_to_c0d(nc0a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0d_to_c0 = hdata%c0d_to_c0
         c0a_to_c0d = hdata%c0a_to_c0d
      else
         ! Receive data on ioproc
         call mpl%recv(nc0d,c0d_to_c0,iproc,mpl%tag+2)
         call mpl%recv(nc0a,c0a_to_c0d,iproc,mpl%tag+3)
      end if

      ! Allocation
      com_AD(iproc)%nred = nc0a
      com_AD(iproc)%next = nc0d
      allocate(com_AD(iproc)%ext_to_proc(com_AD(iproc)%next))
      allocate(com_AD(iproc)%ext_to_red(com_AD(iproc)%next))
      allocate(com_AD(iproc)%red_to_ext(com_AD(iproc)%nred))

      ! AD communication
      do ic0d=1,nc0d
         ic0 = c0d_to_c0(ic0d)
         com_AD(iproc)%ext_to_proc(ic0d) = geom%c0_to_proc(ic0)
         ic0a = geom%c0_to_c0a(ic0)
         com_AD(iproc)%ext_to_red(ic0d) = ic0a
      end do
      com_AD(iproc)%red_to_ext = c0a_to_c0d

      ! Release memory
      deallocate(c0d_to_c0)
      deallocate(c0a_to_c0d)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(geom%nc0a,mpl%ioproc,mpl%tag)
   call mpl%send(hdata%nc0d,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl%send(hdata%nc0d,hdata%c0d_to_c0,mpl%ioproc,mpl%tag+2)
   call mpl%send(geom%nc0a,hdata%c0a_to_c0d,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4
call hdata%com_AD%setup(com_AD,'com_AD')

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0d =       ',hdata%nc0d

end subroutine hdata_compute_mpi_d

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_c
!> Purpose: compute HDIAG MPI distribution, halo C
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_c(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: iproc,jc3,ic0,ic0a,ic0c,ic1,ic1a,its,il0,d_n_s_max,d_n_s_max_loc,i_s,i_s_loc
integer :: nc0a,nc0c
integer,allocatable :: interpd_lg(:,:,:),c0c_to_c0(:),c0a_to_c0c(:)
real(kind_real),allocatable :: lon_c1(:),lat_c1(:)
type(com_type) :: com_AC(mpl%nproc)
type(linop_type),allocatable :: dfull(:,:)

if (nam%displ_diag) then
   ! Allocation
   allocate(dfull(geom%nl0,nam%nts))
   allocate(lon_c1(nam%nc1))
   allocate(lat_c1(nam%nc1))

   ! Prepare displacement interpolation
   do its=1,nam%nts
      do il0=1,geom%nl0
         ! Copy Sc1 points
         do ic1=1,nam%nc1
            ic0 = hdata%c1_to_c0(ic1)
            lon_c1(ic1) = hdata%displ_lon(ic0,il0,its)
            lat_c1(ic1) = hdata%displ_lat(ic0,il0,its)
         end do

         ! Compute interpolation
         call dfull(il0,its)%interp(geom%mesh,geom%ctree,geom%nc0,geom%mask(:,il0),nam%nc1,lon_c1,lat_c1,hdata%c1l0_log(:,il0), &
       & nam%diag_interp)
      end do
   end do
end if

! Allocation
if (nam%displ_diag) then
   d_n_s_max = 0
   do its=1,nam%nts
      do il0=1,geom%nl0
         d_n_s_max = max(d_n_s_max,dfull(il0,its)%n_s)
      end do
   end do
   allocate(hdata%lcheck_d(d_n_s_max,geom%nl0,nam%nts))
   allocate(hdata%d(geom%nl0,nam%nts))
end if
allocate(hdata%lcheck_c0c(geom%nc0))

! Halo C
hdata%lcheck_c0c = hdata%lcheck_c0a
do jc3=1,nam%nc3
   do ic1a=1,hdata%nc1a
      ic1 = hdata%c1a_to_c1(ic1a)
      if (any(hdata%c1c3l0_log(ic1,jc3,:))) then
         ic0 = hdata%c1c3_to_c0(ic1,jc3)
         hdata%lcheck_c0c(ic0) = .true.
      end if
   end do
end do
if (nam%displ_diag) then
   hdata%lcheck_d = .false.
   do its=1,nam%nts
      do il0=1,geom%nl0
         do i_s=1,dfull(il0,its)%n_s
            ic0 = dfull(il0,its)%col(i_s)
            ic1 = dfull(il0,its)%row(i_s)
            if (hdata%lcheck_c1a(ic1)) then
               hdata%lcheck_c0c(ic0) = .true.
               hdata%lcheck_d(i_s,il0,its) = .true.
            end if
         end do
      end do
   end do
end if
hdata%nc0c = count(hdata%lcheck_c0c)
if (nam%displ_diag) then
   do its=1,nam%nts
      do il0=1,geom%nl0
         hdata%d(il0,its)%n_s = count(hdata%lcheck_d(:,il0,its))
      end do
   end do
end if

! Global <-> local conversions for fields

! Halo C
allocate(hdata%c0c_to_c0(hdata%nc0c))
allocate(hdata%c0_to_c0c(geom%nc0))
call msi(hdata%c0_to_c0c)
ic0c = 0
do ic0=1,geom%nc0
   if (hdata%lcheck_c0c(ic0)) then
      ic0c = ic0c+1
      hdata%c0c_to_c0(ic0c) = ic0
      hdata%c0_to_c0c(ic0) = ic0c
   end if
end do

! Inter-halo conversions
allocate(hdata%c0a_to_c0c(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0c = hdata%c0_to_c0c(ic0)
   hdata%c0a_to_c0c(ic0a) = ic0c
end do

if (nam%displ_diag) then
   ! Global <-> local conversions for data
   d_n_s_max_loc = 0
   do its=1,nam%nts
      do il0=1,geom%nl0
         d_n_s_max_loc = max(d_n_s_max_loc,hdata%d(il0,its)%n_s)
      end do
   end do
   allocate(interpd_lg(d_n_s_max_loc,geom%nl0,nam%nts))
   do its=1,nam%nts
      do il0=1,geom%nl0
         i_s_loc = 0
         do i_s=1,dfull(il0,its)%n_s
            if (hdata%lcheck_d(i_s,il0,its)) then
               i_s_loc = i_s_loc+1
               interpd_lg(i_s_loc,il0,its) = i_s
            end if
         end do
      end do
   end do

   ! Local data
   do its=1,nam%nts
      do il0=1,geom%nl0
         write(hdata%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
         hdata%d(il0,its)%n_src = hdata%nc0c
         hdata%d(il0,its)%n_dst = hdata%nc1a
         call hdata%d(il0,its)%alloc
         do i_s_loc=1,hdata%d(il0,its)%n_s
            i_s = interpd_lg(i_s_loc,il0,its)
            hdata%d(il0,its)%row(i_s_loc) = hdata%c1_to_c1a(dfull(il0,its)%row(i_s))
            hdata%d(il0,its)%col(i_s_loc) = hdata%c0_to_c0c(dfull(il0,its)%col(i_s))
            hdata%d(il0,its)%S(i_s_loc) = dfull(il0,its)%S(i_s)
         end do
         call hdata%d(il0,its)%reorder
      end do
   end do
end if

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = geom%nc0a
         nc0c = hdata%nc0c
      else
         ! Receive dimensions on ioproc
         call mpl%recv(nc0a,iproc,mpl%tag)
         call mpl%recv(nc0c,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(c0c_to_c0(nc0c))
      allocate(c0a_to_c0c(nc0a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0c_to_c0 = hdata%c0c_to_c0
         c0a_to_c0c = hdata%c0a_to_c0c
      else
         ! Receive data on ioproc
         call mpl%recv(nc0c,c0c_to_c0,iproc,mpl%tag+2)
         call mpl%recv(nc0a,c0a_to_c0c,iproc,mpl%tag+3)
      end if

      ! Allocation
      com_AC(iproc)%nred = nc0a
      com_AC(iproc)%next = nc0c
      allocate(com_AC(iproc)%ext_to_proc(com_AC(iproc)%next))
      allocate(com_AC(iproc)%ext_to_red(com_AC(iproc)%next))
      allocate(com_AC(iproc)%red_to_ext(com_AC(iproc)%nred))

      ! AC communication
      do ic0c=1,nc0c
         ic0 = c0c_to_c0(ic0c)
         com_AC(iproc)%ext_to_proc(ic0c) = geom%c0_to_proc(ic0)
         ic0a = geom%c0_to_c0a(ic0)
         com_AC(iproc)%ext_to_red(ic0c) = ic0a
      end do
      com_AC(iproc)%red_to_ext = c0a_to_c0c

      ! Release memory
      deallocate(c0c_to_c0)
      deallocate(c0a_to_c0c)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(geom%nc0a,mpl%ioproc,mpl%tag)
   call mpl%send(hdata%nc0c,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl%send(hdata%nc0c,hdata%c0c_to_c0,mpl%ioproc,mpl%tag+2)
   call mpl%send(geom%nc0a,hdata%c0a_to_c0c,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4
call hdata%com_AC%setup(com_AC,'com_AC')

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0c =      ',hdata%nc0c

end subroutine hdata_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: hdata_diag_filter
!> Purpose: filter diagnostics
!----------------------------------------------------------------------
subroutine hdata_diag_filter(hdata,geom,il0,filter_type,r,diag)

implicit none

! Passed variables
class(hdata_type),intent(in) :: hdata             !< HDIAG data
type(geom_type),intent(in) :: geom                !< Geometry
integer,intent(in) :: il0                         !< Level
character(len=*),intent(in) :: filter_type        !< Filter type
real(kind_real),intent(in) :: r                   !< Filter support radius
real(kind_real),intent(inout) :: diag(hdata%nc2a) !< Filtered diagnostics

! Local variables
integer :: ic2a,ic2,ic1,jc2,nc2eff
integer,allocatable :: order(:)
real(kind_real) :: diag_tmp(hdata%nc2),distnorm,norm,wgt
real(kind_real),allocatable :: diag_eff(:),diag_eff_dist(:)

! Local to global
call hdata%diag_com_lg(1,diag,diag_tmp)

! Broadcast
call mpl%bcast(diag_tmp,mpl%ioproc)

!$omp parallel do schedule(static) private(ic2a,ic2,ic1,nc2eff,jc2,distnorm,norm,wgt) firstprivate(diag_eff,diag_eff_dist,order)
do ic2a=1,hdata%nc2a
   ic2 = hdata%c2a_to_c2(ic2a)
   ic1 = hdata%c2_to_c1(ic2)
   if (hdata%c1l0_log(ic1,il0)) then
      ! Allocation
      allocate(diag_eff(hdata%nc2))
      allocate(diag_eff_dist(hdata%nc2))

      ! Build diag_eff of valid points
      nc2eff = 0
      jc2 = 1
      do while (hdata%nn_c2_dist(jc2,ic2,min(il0,geom%nl0i))<r)
         ! Check the point validity
         if (isnotmsr(diag_tmp(hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i))))) then
            nc2eff = nc2eff+1
            diag_eff(nc2eff) = diag_tmp(hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i)))
            diag_eff_dist(nc2eff) = hdata%nn_c2_dist(jc2,ic2,min(il0,geom%nl0i))
         end if
         jc2 = jc2+1
         if (jc2>hdata%nc2) exit
      end do

      ! Apply filter
      if (nc2eff>0) then
         select case (trim(filter_type))
         case ('average')
            ! Compute average
            diag(ic2a) = sum(diag_eff(1:nc2eff))/float(nc2eff)
         case ('gc99')
            ! Gaspari-Cohn (1999) kernel
            diag(ic2a) = 0.0
            norm = 0.0
            do jc2=1,nc2eff
               distnorm = diag_eff_dist(jc2)/r
               wgt = gc99(distnorm)
               diag(ic2a) = diag(ic2a)+wgt*diag_eff(jc2)
               norm = norm+wgt
            end do
            diag(ic2a) = diag(ic2a)/norm
         case ('median')
            ! Compute median
            allocate(order(nc2eff))
            call qsort(nc2eff,diag_eff(1:nc2eff),order)
            if (mod(nc2eff,2)==0) then
               diag(ic2a) = 0.5*(diag_eff(nc2eff/2)+diag_eff(nc2eff/2+1))
            else
               diag(ic2a) = diag_eff((nc2eff+1)/2)
            end if
            deallocate(order)
         case default
            ! Wrong filter
            call msgerror('wrong filter type')
         end select
      else
         call msr(diag(ic2a))
      end if

      ! Release memory
      deallocate(diag_eff)
      deallocate(diag_eff_dist)
   end if
end do
!$omp end parallel do

end subroutine hdata_diag_filter

!----------------------------------------------------------------------
! Subroutine: hdata_diag_com_lg
!> Purpose: communicate diagnostic from local to global distribution
!----------------------------------------------------------------------
subroutine hdata_diag_com_lg(hdata,nl,diag_loc,diag_glb)

implicit none

! Passed variables
class(hdata_type),intent(in) :: hdata                 !< HDIAG data
integer,intent(in) :: nl                              !< Number of levels
real(kind_real),intent(in) :: diag_loc(hdata%nc2a,nl) !< Diagnostic (local)
real(kind_real),intent(out) :: diag_glb(hdata%nc2,nl) !< Diagnostic (global)

! Local variables
integer :: ic2,ic2a,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:),diag_tmp(:,:)
logical,allocatable :: mask_unpack(:,:)

! Allocation
allocate(sbuf(hdata%nc2a*nl))

! Prepare buffer
sbuf = pack(diag_loc,mask=.true.)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      allocate(diag_tmp(hdata%proc_to_nc2a(iproc),nl))
      allocate(mask_unpack(hdata%proc_to_nc2a(iproc),nl))
      allocate(rbuf(hdata%proc_to_nc2a(iproc)*nl))

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Receive data from iproc
         call mpl%recv(hdata%proc_to_nc2a(iproc)*nl,rbuf,iproc,mpl%tag)
      end if

      ! Copy from buffer
      mask_unpack = .true.
      diag_tmp = unpack(rbuf,mask_unpack,diag_tmp)
      do ic2=1,hdata%nc2
         jproc = hdata%c2_to_proc(ic2)
         if (jproc==iproc) then
            ic2a = hdata%c2_to_c2a(ic2)
            diag_glb(ic2,:) = diag_tmp(ic2a,:)
         end if
      end do

      ! Release memory
      deallocate(diag_tmp)
      deallocate(mask_unpack)
      deallocate(rbuf)
   end do
else
   ! Sending data to iproc
   call mpl%send(hdata%nc2a*nl,sbuf,mpl%ioproc,mpl%tag)

   ! Setting at missing value
   call msr(diag_glb)
end if
mpl%tag = mpl%tag+1

end subroutine hdata_diag_com_lg

end module type_hdata
