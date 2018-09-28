!----------------------------------------------------------------------
! Module: type_hdata
!> Purpose: sample data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_hdata

use netcdf
!$ use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,msvali,msvalr
use tools_func, only: gc99,sphere_dist,vector_product,vector_triple_product
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr,ismsi,ismsr
use tools_nc, only: ncfloat
use tools_qsort, only: qsort
use tools_stripack, only: trans
use type_com, only: com_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_kdtree, only: kdtree_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

integer,parameter :: irmax = 10000                           !< Maximum number of random number draws
real(kind_real),parameter :: Lcoast = 1000.0e3_kind_real/req !< Length-scale to increase sampling density along coasts
real(kind_real),parameter :: rcoast = 0.2_kind_real          !< Minimum value to increase sampling density along coasts

! HDIAG data derived type
type hdata_type
   ! Sampling
   real(kind_real),allocatable :: rh_c0(:,:)        !< Sampling radius
   integer,allocatable :: c1_to_c0(:)               !< First sampling index
   logical,allocatable :: c1l0_log(:,:)             !< Log for the first sampling index
   integer,allocatable :: c1c3_to_c0(:,:)           !< Second horizontal sampling index
   logical,allocatable :: c1c3l0_log(:,:,:)         !< Log for the second horizontal sampling index
   integer,allocatable :: c2_to_c1(:)               !< Subgrid to diagnostic points
   integer,allocatable :: c2_to_c0(:)               !< Subgrid to grid

   ! Local data
   logical,allocatable ::  vbal_mask(:,:)           !< Vertical balance mask
   logical,allocatable ::  local_mask(:,:)          !< Local mask
   logical,allocatable ::  displ_mask(:,:)          !< Displacement mask
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
   integer :: nc2f                                  !< Number of points in subset Sc2, halo F
   logical,allocatable :: lcheck_c0a(:)             !< Detection of halo A on subset Sc0
   logical,allocatable :: lcheck_c0c(:)             !< Detection of halo C on subset Sc0
   logical,allocatable :: lcheck_c0d(:)             !< Detection of halo D on subset Sc0
   logical,allocatable :: lcheck_c1a(:)             !< Detection of halo A on subset Sc1
   logical,allocatable :: lcheck_c2a(:)             !< Detection of halo A on subset Sc2
   logical,allocatable :: lcheck_c2b(:)             !< Detection of halo B on subset Sc2
   logical,allocatable :: lcheck_c2f(:)             !< Detection of halo F on subset Sc2
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
   integer,allocatable :: c2f_to_c2(:)              !< Subset Sc2, halo B to global
   integer,allocatable :: c2_to_c2f(:)              !< Subset Sc2, global to halo B
   integer,allocatable :: c2a_to_c2f(:)             !< Subset Sc2, halo A to halo B
   integer,allocatable :: c2_to_proc(:)             !< Subset Sc2, global to processor
   integer,allocatable :: proc_to_nc2a(:)           !< Number of points in subset Sc2, halo A, for each processor
   type(com_type) :: com_AC                         !< Communication between halos A and C
   type(com_type) :: com_AB                         !< Communication between halos A and B
   type(com_type) :: com_AD                         !< Communication between halos A and D
   type(com_type) :: com_AF                         !< Communication between halos A and F (filtering)
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
   procedure :: compute_mpi_f => hdata_compute_mpi_f
   procedure :: diag_filter => hdata_diag_filter
   procedure :: diag_fill => hdata_diag_fill
end type hdata_type

private
public :: hdata_type

contains

!----------------------------------------------------------------------
! Subroutine: hdata_alloc
!> Purpose: HDIAG data allocation
!----------------------------------------------------------------------
subroutine hdata_alloc(hdata,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Allocation
allocate(hdata%rh_c0(geom%nc0,geom%nl0))
allocate(hdata%c1_to_c0(nam%nc1))
allocate(hdata%c1l0_log(nam%nc1,geom%nl0))
allocate(hdata%c1c3_to_c0(nam%nc1,nam%nc3))
allocate(hdata%c1c3l0_log(nam%nc1,nam%nc3,geom%nl0))
if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%var_diag.or.nam%local_diag.or.nam%displ_diag))) then
   allocate(hdata%c2_to_c1(nam%nc2))
   allocate(hdata%c2_to_c0(nam%nc2))
   allocate(hdata%nn_c2_index(nam%nc2,nam%nc2,geom%nl0i))
   allocate(hdata%nn_c2_dist(nam%nc2,nam%nc2,geom%nl0i))
   allocate(hdata%hfull(geom%nl0i))
end if
if (nam%new_vbal) allocate(hdata%vbal_mask(nam%nc1,nam%nc2))
if (nam%local_diag) allocate(hdata%local_mask(nam%nc1,nam%nc2))
if (nam%displ_diag) then
   allocate(hdata%displ_mask(nam%nc1,nam%nc2))
   allocate(hdata%displ_lon(geom%nc0a,geom%nl0,nam%nts))
   allocate(hdata%displ_lat(geom%nc0a,geom%nl0,nam%nts))
end if

! Initialization
call msr(hdata%rh_c0)
call msi(hdata%c1_to_c0)
hdata%c1l0_log = .false.
call msi(hdata%c1c3_to_c0)
hdata%c1c3l0_log = .false.
if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%var_diag.or.nam%local_diag.or.nam%displ_diag))) then
   call msi(hdata%c2_to_c1)
   call msi(hdata%c2_to_c0)
   call msi(hdata%nn_c2_index)
   call msr(hdata%nn_c2_dist)
end if
if (nam%new_vbal) hdata%vbal_mask = .false.
if (nam%local_diag) hdata%local_mask = .false.
if (nam%displ_diag) then
   hdata%displ_mask = .false.
   call msr(hdata%displ_lon)
   call msr(hdata%displ_lat)
end if

end subroutine hdata_alloc

!----------------------------------------------------------------------
! Subroutine: hdata_dealloc
!> Purpose: HDIAG data deallocation
!----------------------------------------------------------------------
subroutine hdata_dealloc(hdata,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: il0

! Release memory
if (allocated(hdata%rh_c0)) deallocate(hdata%rh_c0)
if (allocated(hdata%c1_to_c0)) deallocate(hdata%c1_to_c0)
if (allocated(hdata%c1l0_log)) deallocate(hdata%c1l0_log)
if (allocated(hdata%c1c3_to_c0)) deallocate(hdata%c1c3_to_c0)
if (allocated(hdata%c1c3l0_log)) deallocate(hdata%c1c3l0_log)
if (allocated(hdata%c2_to_c1)) deallocate(hdata%c2_to_c1)
if (allocated(hdata%c2_to_c0)) deallocate(hdata%c2_to_c0)
if (allocated(hdata%c2a_to_c2)) deallocate(hdata%c2a_to_c2)
if (allocated(hdata%c2_to_c2a)) deallocate(hdata%c2_to_c2a)
if (allocated(hdata%vbal_mask)) deallocate(hdata%vbal_mask)
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
!> Purpose: read HDIAG data
!----------------------------------------------------------------------
subroutine hdata_read(hdata,mpl,nam,geom,ios)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(in) :: mpl         !< MPI data
type(nam_type),intent(inout) :: nam      !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
integer,intent(out) :: ios               !< Status flag

! Local variables
integer :: il0,il0i,ic1,jc3,ic2
integer :: nl0_test,nl0r_test,nc_test,nc1_test,nc2_test,nc2_1_test,nc2_2_test
integer :: info,ncid,nl0_id,nc3_id,nc1_id,nc2_id,nc2_1_id,nc2_2_id
integer :: c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id,vbal_mask_id,local_mask_id,displ_mask_id,nn_c2_index_id,nn_c2_dist_id
integer :: c1l0_logint(nam%nc1,geom%nl0),c1c3l0_logint(nam%nc1,nam%nc3,geom%nl0)
integer,allocatable :: vbal_maskint(:,:),local_maskint(:,:),displ_maskint(:,:)
character(len=3) :: il0ichar
character(len=1024) :: subr = 'hdata_read'

! Initialization
ios = 0

! Open file
info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',nf90_nowrite,ncid)
if (info/=nf90_noerr) then
   call mpl%warning('cannot find HDIAG data to read, recomputing HDIAG sampling')
   nam%sam_write = .true.
   ios = 1
   return
end if

! Check dimensions
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=nl0_test))
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'nl0r',nl0r_test))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc3',nc3_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc3_id,len=nc_test))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc1',nc1_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc1_id,len=nc1_test))
if (nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag) then
   info = nf90_inq_dimid(ncid,'nc2',nc2_id)
   if (info==nf90_noerr) then
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc2_id,len=nc2_test))
   else
      call mpl%warning('cannot find nc2 when reading HDIAG sampling, recomputing HDIAG sampling')
      nam%sam_write = .true.
      ios = 2
   end if
end if
if ((geom%nl0/=nl0_test).or.(nam%nl0r/=nl0r_test).or.(nam%nc3/=nc_test).or.(nam%nc1/=nc1_test)) then
   call mpl%warning('wrong dimension when reading HDIAG sampling, recomputing HDIAG sampling')
   nam%sam_write = .true.
   call mpl%ncerr(subr,nf90_close(ncid))
   ios = 1
   return
end if
if (nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag) then
   if (nam%nc2/=nc2_test) then
      call mpl%warning('wrong dimension when reading HDIAG sampling, recomputing HDIAG sampling')
      nam%sam_write = .true.
      ios = 2
   end if
end if

write(mpl%info,'(a7,a)') '','Read HDIAG sampling'
call flush(mpl%info)

! Get arrays ID
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1_to_c0',c1_to_c0_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1l0_log',c1l0_log_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1c3_to_c0',c1c3_to_c0_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'c1c3l0_log',c1c3l0_log_id))
if ((ios==0).and.(nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag)) then
   call mpl%ncerr(subr,nf90_inq_varid(ncid,'c2_to_c1',c2_to_c1_id))
   call mpl%ncerr(subr,nf90_inq_varid(ncid,'c2_to_c0',c2_to_c0_id))
end if

! Read arrays
call mpl%ncerr(subr,nf90_get_var(ncid,c1_to_c0_id,hdata%c1_to_c0))
call mpl%ncerr(subr,nf90_get_var(ncid,c1l0_log_id,c1l0_logint))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (c1l0_logint(ic1,il0)==0) then
         hdata%c1l0_log(ic1,il0) = .false.
      else if (c1l0_logint(ic1,il0)==1) then
         hdata%c1l0_log(ic1,il0) = .true.
      end if
   end do
end do
call mpl%ncerr(subr,nf90_get_var(ncid,c1c3_to_c0_id,hdata%c1c3_to_c0))
call mpl%ncerr(subr,nf90_get_var(ncid,c1c3l0_log_id,c1c3l0_logint))
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
if ((ios==0).and.(nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag)) then
   call mpl%ncerr(subr,nf90_get_var(ncid,c2_to_c1_id,hdata%c2_to_c1))
   call mpl%ncerr(subr,nf90_get_var(ncid,c2_to_c0_id,hdata%c2_to_c0))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Read nearest neighbors and interpolation
do il0i=1,geom%nl0i
   if ((ios==0).and.(nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%var_diag.or.nam%local_diag.or.nam%displ_diag)))) then
      ! Open file
      write(il0ichar,'(i3.3)') il0i
      info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc',nf90_nowrite,ncid)
      if (info/=nf90_noerr) then
         call mpl%warning('cannot find nearest neighbors and interpolation data to read, recomputing HDIAG sampling')
         nam%sam_write = .true.
         ios = 3
         return
      end if

      ! Check dimensions
      call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc1',nc1_id))
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc1_id,len=nc1_test))
      call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc2_1',nc2_1_id))
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc2_1_id,len=nc2_1_test))
      call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc2_2',nc2_2_id))
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc2_2_id,len=nc2_2_test))
      if ((nam%nc1/=nc1_test).or.(nam%nc2/=nc2_1_test).or.(nam%nc2/=nc2_2_test)) then
         call mpl%warning('wrong dimension when reading HDIAG sampling, recomputing HDIAG sampling')
         nam%sam_write = .true.
         ios = 2
      end if
      info = nf90_inq_varid(ncid,'nn_c2_index',nn_c2_index_id)
      if (info==nf90_noerr) then
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'nn_c2_dist',nn_c2_dist_id))
         call mpl%ncerr(subr,nf90_get_var(ncid,nn_c2_index_id,hdata%nn_c2_index(:,:,il0i)))
         call mpl%ncerr(subr,nf90_get_var(ncid,nn_c2_dist_id,hdata%nn_c2_dist(:,:,il0i)))
      else
         call mpl%warning('cannot find nc2 nearest neighbors data to read, recomputing HDIAG sampling')
         nam%sam_write = .true.
         ios = 4
      end if
      write(hdata%hfull(il0i)%prefix,'(a,i3.3)') 'hfull_',il0i
      call hdata%hfull(il0i)%read(mpl,ncid)
   end if

   if ((ios==0).and.nam%new_vbal) then
      ! Allocation
      allocate(vbal_maskint(nam%nc1,nam%nc2))

      call mpl%ncerr(subr,nf90_inq_varid(ncid,'vbal_mask',vbal_mask_id))
      call mpl%ncerr(subr,nf90_get_var(ncid,vbal_mask_id,vbal_maskint))
      do ic2=1,nam%nc2
         do ic1=1,nam%nc1
            if (vbal_maskint(ic1,ic2)==1) then
               hdata%vbal_mask(ic1,ic2) = .true.
            elseif (vbal_maskint(ic1,ic2)==0) then
               hdata%vbal_mask(ic1,ic2) = .false.
            else
               call mpl%abort('wrong vbal_mask')
           end if
         end do
      end do

      ! Release memory
      deallocate(vbal_maskint)
   end if

   if ((ios==0).and.nam%local_diag) then
      ! Allocation
      allocate(local_maskint(nam%nc1,nam%nc2))

      call mpl%ncerr(subr,nf90_inq_varid(ncid,'local_mask',local_mask_id))
      call mpl%ncerr(subr,nf90_get_var(ncid,local_mask_id,local_maskint))
      do ic2=1,nam%nc2
         do ic1=1,nam%nc1
            if (local_maskint(ic1,ic2)==1) then
               hdata%local_mask(ic1,ic2) = .true.
            elseif (local_maskint(ic1,ic2)==0) then
               hdata%local_mask(ic1,ic2) = .false.
            else
               call mpl%abort('wrong local_mask')
           end if
         end do
      end do

      ! Release memory
      deallocate(local_maskint)
   end if

   if ((ios==0).and.nam%displ_diag) then
      ! Allocation
      allocate(displ_maskint(nam%nc1,nam%nc2))

      call mpl%ncerr(subr,nf90_inq_varid(ncid,'displ_mask',displ_mask_id))
      call mpl%ncerr(subr,nf90_get_var(ncid,displ_mask_id,displ_maskint))
      do ic2=1,nam%nc2
         do ic1=1,nam%nc1
            if (displ_maskint(ic1,ic2)==1) then
               hdata%displ_mask(ic1,ic2) = .true.
            elseif (displ_maskint(ic1,ic2)==0) then
               hdata%displ_mask(ic1,ic2) = .false.
            else
                call mpl%abort('wrong displ_mask')
            end if
         end do
      end do

      ! Release memory
      deallocate(displ_maskint)
   end if

   if ((ios==0).and.(nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag)) then
      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   end if
end do

end subroutine hdata_read

!----------------------------------------------------------------------
! Subroutine: hdata_write
!> Purpose: write HDIAG data
!----------------------------------------------------------------------
subroutine hdata_write(hdata,mpl,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(in) :: hdata !< HDIAG data
type(mpl_type),intent(in) :: mpl         !< MPI data
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(in) :: geom    !< Geometry

! Local variables
integer :: il0,il0i,ic1,jc3,ic2
integer :: ncid,nl0_id,nc1_id,nc2_id,nc2_1_id,nc2_2_id,nc3_id
integer :: lat_id,lon_id,smax_id,c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id,vbal_mask_id,local_mask_id,displ_mask_id,nn_c2_index_id,nn_c2_dist_id
integer :: c1l0_logint(nam%nc1,geom%nl0),c1c3l0_logint(nam%nc1,nam%nc3,geom%nl0)
integer,allocatable :: vbal_maskint(:,:),local_maskint(:,:),displ_maskint(:,:)
real(kind_real) :: lon(nam%nc1,nam%nc3,geom%nl0),lat(nam%nc1,nam%nc3,geom%nl0)
character(len=3) :: il0ichar
character(len=1024) :: subr = 'hdata_write'

! Processor verification
if (.not.mpl%main) call mpl%abort('only I/O proc should enter '//trim(subr))

! Create file
write(mpl%info,'(a7,a)') '','Write HDIAG sampling'
call flush(mpl%info)
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%ncwrite(mpl,ncid)

! Define dimensions
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0r',nam%nl0r))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc3',nam%nc3,nc3_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%var_diag.or.nam%local_diag.or.nam%displ_diag))) &
 & call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2',nam%nc2,nc2_id))

! Define variables
call mpl%ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc1_id,nc3_id,nl0_id/),lat_id))
call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
call mpl%ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc1_id,nc3_id,nl0_id/),lon_id))
call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
call mpl%ncerr(subr,nf90_def_var(ncid,'smax',ncfloat,(/nc3_id,nl0_id/),smax_id))
call mpl%ncerr(subr,nf90_put_att(ncid,smax_id,'_FillValue',msvalr))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1_to_c0',nf90_int,(/nc1_id/),c1_to_c0_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1_to_c0_id,'_FillValue',msvali))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1l0_log',nf90_int,(/nc1_id,nl0_id/),c1l0_log_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1l0_log_id,'_FillValue',msvali))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1c3_to_c0',nf90_int,(/nc1_id,nc3_id/),c1c3_to_c0_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1c3_to_c0_id,'_FillValue',msvali))
call mpl%ncerr(subr,nf90_def_var(ncid,'c1c3l0_log',nf90_int,(/nc1_id,nc3_id,nl0_id/),c1c3l0_log_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1c3l0_log_id,'_FillValue',msvali))
if (nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'c2_to_c1',nf90_int,(/nc2_id/),c2_to_c1_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,c2_to_c1_id,'_FillValue',msvali))
   call mpl%ncerr(subr,nf90_def_var(ncid,'c2_to_c0',nf90_int,(/nc2_id/),c2_to_c0_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,c2_to_c0_id,'_FillValue',msvali))
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Convert data
call msr(lon)
call msr(lat)
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (hdata%c1c3l0_log(ic1,jc3,il0)) then
            lon(ic1,jc3,il0) = geom%lon(hdata%c1c3_to_c0(ic1,jc3))*rad2deg
            lat(ic1,jc3,il0) = geom%lat(hdata%c1c3_to_c0(ic1,jc3))*rad2deg
            c1c3l0_logint(ic1,jc3,il0) = 1
         else
            c1c3l0_logint(ic1,jc3,il0) = 0
         end if
      end do
   end do
   do ic1=1,nam%nc1
      if (hdata%c1l0_log(ic1,il0)) then
         c1l0_logint(ic1,il0) = 1
      else
         c1l0_logint(ic1,il0) = 0
      end if
   end do
end do

! Write variables
call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_put_var(ncid,smax_id,real(count(hdata%c1c3l0_log,dim=1),kind_real)))
call mpl%ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,hdata%c1_to_c0))
call mpl%ncerr(subr,nf90_put_var(ncid,c1l0_log_id,c1l0_logint))
call mpl%ncerr(subr,nf90_put_var(ncid,c1c3_to_c0_id,hdata%c1c3_to_c0))
call mpl%ncerr(subr,nf90_put_var(ncid,c1c3l0_log_id,c1c3l0_logint))
if (nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag) then
   call mpl%ncerr(subr,nf90_put_var(ncid,c2_to_c1_id,hdata%c2_to_c1))
   call mpl%ncerr(subr,nf90_put_var(ncid,c2_to_c0_id,hdata%c2_to_c0))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Write nearest neighbors and interpolation
do il0i=1,geom%nl0i
   if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%var_diag.or.nam%local_diag.or.nam%displ_diag))) then
      ! Create file
      write(il0ichar,'(i3.3)') il0i
      call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc', &
    & or(nf90_clobber,nf90_64bit_offset),ncid))

      ! Write namelist parameters
      call nam%ncwrite(mpl,ncid)

      ! Define dimensions
      call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
      call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2_1',nam%nc2,nc2_1_id))
      call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2_2',nam%nc2,nc2_2_id))

      ! Define variables
      call mpl%ncerr(subr,nf90_def_var(ncid,'nn_c2_index',nf90_int,(/nc2_1_id,nc2_2_id/),nn_c2_index_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,nn_c2_index_id,'_FillValue',msvali))
      call mpl%ncerr(subr,nf90_def_var(ncid,'nn_c2_dist',ncfloat,(/nc2_1_id,nc2_2_id/),nn_c2_dist_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,nn_c2_dist_id,'_FillValue',msvalr))

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      call mpl%ncerr(subr,nf90_put_var(ncid,nn_c2_index_id,hdata%nn_c2_index(:,:,il0i)))
      call mpl%ncerr(subr,nf90_put_var(ncid,nn_c2_dist_id,hdata%nn_c2_dist(:,:,il0i)))
      call hdata%hfull(il0i)%write(mpl,ncid)
   end if

   if (nam%new_vbal) then
      ! Allocation
      allocate(vbal_maskint(nam%nc1,nam%nc2))

      ! Definition mode
      call mpl%ncerr(subr,nf90_redef(ncid))

      ! Define variables
      call mpl%ncerr(subr,nf90_def_var(ncid,'vbal_mask',nf90_int,(/nc1_id,nc2_1_id/),vbal_mask_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,vbal_mask_id,'_FillValue',msvali))

      ! Convert data
      do ic2=1,nam%nc2
         do ic1=1,nam%nc1
            if (hdata%vbal_mask(ic1,ic2)) then
               vbal_maskint(ic1,ic2) = 1
            else
               vbal_maskint(ic1,ic2) = 0
            end if
         end do
      end do

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      call mpl%ncerr(subr,nf90_put_var(ncid,vbal_mask_id,vbal_maskint))

      ! Release memory
      deallocate(vbal_maskint)
   end if

   if (nam%local_diag) then
      ! Allocation
      allocate(local_maskint(nam%nc1,nam%nc2))

      ! Definition mode
      call mpl%ncerr(subr,nf90_redef(ncid))

      ! Define variables
      call mpl%ncerr(subr,nf90_def_var(ncid,'local_mask',nf90_int,(/nc1_id,nc2_1_id/),local_mask_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,local_mask_id,'_FillValue',msvali))

      ! Convert data
      do ic2=1,nam%nc2
         do ic1=1,nam%nc1
            if (hdata%local_mask(ic1,ic2)) then
               local_maskint(ic1,ic2) = 1
            else
               local_maskint(ic1,ic2) = 0
            end if
         end do
      end do

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      call mpl%ncerr(subr,nf90_put_var(ncid,local_mask_id,local_maskint))

      ! Release memory
      deallocate(local_maskint)
   end if

   if (nam%displ_diag) then
      ! Allocation
      allocate(displ_maskint(nam%nc1,nam%nc2))

      ! Definition mode
      call mpl%ncerr(subr,nf90_redef(ncid))

      ! Define variables
      call mpl%ncerr(subr,nf90_def_var(ncid,'displ_mask',nf90_int,(/nc1_id,nc2_1_id/),displ_mask_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,displ_mask_id,'_FillValue',msvali))

      ! Convert data
      do ic2=1,nam%nc2
         do ic1=1,nam%nc1
            if (hdata%displ_mask(ic1,ic2)) then
               displ_maskint(ic1,ic2) = 1
            else
               displ_maskint(ic1,ic2) = 0
            end if
         end do
      end do

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      call mpl%ncerr(subr,nf90_put_var(ncid,displ_mask_id,displ_maskint))

      ! Release memory
      deallocate(displ_maskint)
   end if

   if (nam%new_lct.or.nam%var_diag.or.nam%local_diag.or.nam%displ_diag) then
      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   end if
end do

end subroutine hdata_write

!----------------------------------------------------------------------
! Subroutine: hdata_setup_sampling
!> Purpose: setup sampling
!----------------------------------------------------------------------
subroutine hdata_setup_sampling(hdata,mpl,rng,nam,geom,io)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(rng_type),intent(inout) :: rng      !< Random number generator
type(nam_type),intent(inout) :: nam      !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
type(io_type),intent(in) :: io           !< I/O

! Local variables
integer :: ios,ic0,il0,ic1,ic2,ildw,jc3,il0i,jc1,kc1,nc2_eff
integer,allocatable :: vbot(:),vtop(:),nn_c1_index(:)
real(kind_real) :: lon_c1(nam%nc1),lat_c1(nam%nc1)
real(kind_real) :: lon_c2(nam%nc2),lat_c2(nam%nc2)
real(kind_real) :: rh_c0a(geom%nc0a,geom%nl0),nn_dist(1)
real(kind_real),allocatable :: rh_c1(:),nn_c1_dist(:)
logical :: mask_c1(nam%nc1),mask_c2(nam%nc2)
character(len=1024) :: filename
type(kdtree_type) :: kdtree
type(linop_type) :: hbase

! Check subsampling size
if (nam%nc1>maxval(geom%nc0_mask)) then
   call mpl%warning('nc1 is too large for then mask, reset nc1 to the largest possible value')
   nam%nc1 = maxval(geom%nc0_mask)
end if

! Allocation
call hdata%alloc(nam,geom)

! Read or compute sampling data
ios = 1
if (nam%sam_read) call hdata%read(mpl,nam,geom,ios)
if (ios==1) then
   ! Compute zero-separation sampling
   call hdata%compute_sampling_zs(mpl,rng,nam,geom)

   if (nam%new_lct) then
      ! Compute LCT sampling
      call hdata%compute_sampling_lct(mpl,nam,geom)
   elseif (nam%new_vbal.or.nam%new_hdiag) then
      ! Compute positive separation sampling
      call hdata%compute_sampling_ps(mpl,rng,nam,geom)
   end if

   ! Compute sampling mask
   call hdata%compute_sampling_mask(nam,geom)
end if

if (nam%new_vbal.or.nam%new_lct.or.(nam%new_hdiag.and.(nam%var_diag.or.nam%local_diag.or.nam%displ_diag))) then
   if ((ios==1).or.(ios==2)) then
      ! Define subsampling
      write(mpl%info,'(a7,a)',advance='no') '','Define subsampling:'
      call flush(mpl%info)
      allocate(rh_c1(nam%nc1))
      lon_c1 = geom%lon(hdata%c1_to_c0)
      lat_c1 = geom%lat(hdata%c1_to_c0)
      mask_c1 = .true.
      rh_c1 = 1.0
      call rng%initialize_sampling(mpl,nam%nc1,lon_c1,lat_c1,mask_c1,rh_c1,nam%ntry,nam%nrep,nam%nc2,hdata%c2_to_c1)
      hdata%c2_to_c0 = hdata%c1_to_c0(hdata%c2_to_c1)
      deallocate(rh_c1)
   end if

   if ((ios==1).or.(ios==2).or.(ios==3).or.(ios==4)) then
      ! Find nearest neighbors
      write(mpl%info,'(a7,a)') '','Find nearest neighbors'
      call flush(mpl%info)
      do il0=1,geom%nl0
         if ((il0==1).or.(geom%nl0i>1)) then
            write(mpl%info,'(a10,a,i3)') '','Level ',nam%levs(il0)
            call flush(mpl%info)
            mask_c2 = hdata%c1l0_log(hdata%c2_to_c1,il0)
            if (any(mask_c2)) then
               lon_c2 = geom%lon(hdata%c2_to_c0)
               lat_c2 = geom%lat(hdata%c2_to_c0)
               call kdtree%create(mpl,nam%nc2,lon_c2,lat_c2,mask=mask_c2)
               do ic2=1,nam%nc2
                  ic1 = hdata%c2_to_c1(ic2)
                  ic0 = hdata%c2_to_c0(ic2)
                  nc2_eff = min(nam%nc2,count(hdata%c1l0_log(hdata%c2_to_c1,il0)))
                  if (hdata%c1l0_log(ic1,il0)) call kdtree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0), &
                   & nc2_eff,hdata%nn_c2_index(1:nc2_eff,ic2,il0),hdata%nn_c2_dist(1:nc2_eff,ic2,il0))
               end do
               call kdtree%dealloc
            end if
         end if
      end do
   end if

   ! Compute sampling mesh
   write(mpl%info,'(a7,a)') '','Compute sampling mesh'
   call flush(mpl%info)
   lon_c2 = geom%lon(hdata%c2_to_c0)
   lat_c2 = geom%lat(hdata%c2_to_c0)
   call hdata%mesh%create(mpl,rng,nam%nc2,lon_c2,lat_c2)

   ! Compute triangles list
   write(mpl%info,'(a7,a)') '','Compute triangles list '
   call flush(mpl%info)
   call hdata%mesh%trlist

   ! Find boundary nodes
   write(mpl%info,'(a7,a)') '','Find boundary nodes'
   call flush(mpl%info)
   call hdata%mesh%bnodes

   ! Find boundary arcs
   write(mpl%info,'(a7,a)') '','Find boundary arcs'
   call flush(mpl%info)
   call hdata%mesh%barcs

   if ((ios==1).or.(ios==2).or.(ios==3)) then
      if (nam%new_vbal.or.nam%local_diag.or.nam%displ_diag) then
         ! Allocation
         allocate(nn_c1_index(nam%nc1))
         allocate(nn_c1_dist(nam%nc1))

         ! Compute local masks
         write(mpl%info,'(a7,a)') '','Compute local masks'
         call flush(mpl%info)
         lon_c1 = geom%lon(hdata%c1_to_c0)
         lat_c1 = geom%lat(hdata%c1_to_c0)
         mask_c1 = any(hdata%c1l0_log,dim=2)
         call kdtree%create(mpl,nam%nc1,lon_c1,lat_c1,mask=mask_c1)
         do ic2=1,nam%nc2
            ! Inidices
            ic1 = hdata%c2_to_c1(ic2)
            ic0 = hdata%c2_to_c0(ic2)

            ! Find nearest neighbors
            call kdtree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nam%nc1,nn_c1_index,nn_c1_dist)
            do jc1=1,nam%nc1
               kc1 = nn_c1_index(jc1)
               if (nam%new_vbal) hdata%vbal_mask(kc1,ic2) = (jc1==1) &
             & .or.(nn_c1_dist(jc1)<nam%vbal_rad)
               if (nam%local_diag) hdata%local_mask(kc1,ic2) = (jc1==1) &
             & .or.(nn_c1_dist(jc1)<nam%local_rad)
               if (nam%displ_diag) hdata%displ_mask(kc1,ic2) = (jc1==1) &
             & .or.(nn_c1_dist(jc1)<min(nam%displ_rad,hdata%mesh%bdist(ic2)))
            end do
         end do
         call kdtree%dealloc

         ! Release memory
         deallocate(nn_c1_index)
         deallocate(nn_c1_dist)
      end if

      ! Allocation
      allocate(vbot(nam%nc2))
      allocate(vtop(nam%nc2))

      ! Initialize vbot and vtop
      vbot = 1
      vtop = geom%nl0

      ! Compute grid interpolation
      write(mpl%info,'(a7,a)') '','Compute grid interpolation'
      do il0i=1,geom%nl0i
         ! Compute grid interpolation
         write(hdata%hfull(il0i)%prefix,'(a,i3.3)') 'hfull_',il0i
         call hdata%hfull(il0i)%interp(mpl,rng,geom,il0i,nam%nc2,hdata%c2_to_c0,nam%mask_check,vbot,vtop,nam%diag_interp,hbase)
      end do

      ! Release memory
      deallocate(vbot)
      deallocate(vtop)
   end if
end if

! Write sampling data
if (nam%sam_write) then
   if (mpl%main) call hdata%write(mpl,nam,geom)

   ! Write rh_c0
   if (trim(nam%draw_type)=='random_coast') then
      call mpl%glb_to_loc(geom%nl0,geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,hdata%rh_c0,geom%nc0a,rh_c0a)
      filename = trim(nam%prefix)//'_sampling_rh_c0'
      call io%fld_write(mpl,nam,geom,filename,'rh_c0',rh_c0a)
   end if
end if

! Compute nearest neighbors for local diagnostics output
if (nam%local_diag.and.(nam%nldwv>0)) then
   write(mpl%info,'(a7,a)') '','Compute nearest neighbors for local diagnostics output'
   call flush(mpl%info)
   allocate(hdata%nn_ldwv_index(nam%nldwv))
   call kdtree%create(mpl,nam%nc2,geom%lon(hdata%c2_to_c0), &
                geom%lat(hdata%c2_to_c0),mask=hdata%c1l0_log(hdata%c2_to_c1,1))
   do ildw=1,nam%nldwv
      call kdtree%find_nearest_neighbors(nam%lon_ldwv(ildw),nam%lat_ldwv(ildw),1,hdata%nn_ldwv_index(ildw:ildw),nn_dist)
   end do
   call kdtree%dealloc
end if

! Print results
write(mpl%info,'(a7,a)') '','Sampling efficiency (%):'
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),' ~> '
   do jc3=1,nam%nc3
      if (count(hdata%c1c3l0_log(:,jc3,il0))>=nam%nc1/2) then
         ! Sucessful sampling
         write(mpl%info,'(a,i3,a)',advance='no') trim(mpl%green), &
       & int(100.0*real(count(hdata%c1c3l0_log(:,jc3,il0)),kind_real)/real(nam%nc1,kind_real)),trim(mpl%black)
      else
         ! Insufficient sampling
         write(mpl%info,'(a,i3,a)',advance='no') trim(mpl%peach), &
       & int(100.0*real(count(hdata%c1c3l0_log(:,jc3,il0)),kind_real)/real(nam%nc1,kind_real)),trim(mpl%black)
      end if
      if (jc3<nam%nc3) write(mpl%info,'(a)',advance='no') '-'
   end do
   write(mpl%info,'(a)') ' '
end do
call flush(mpl%info)

end subroutine hdata_setup_sampling

!----------------------------------------------------------------------
! Subroutine: hdata_compute_sampling_zs
!> Purpose: compute zero-separation sampling
!----------------------------------------------------------------------
subroutine hdata_compute_sampling_zs(hdata,mpl,rng,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(rng_type),intent(inout) :: rng      !< Random number generator
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: ic0,ic1,il0
integer :: nn_index(1)
real(kind_real) :: nn_dist(1)
type(kdtree_type) :: kdtree

! Compute subset
if (nam%nc1<maxval(geom%nc0_mask)) then
   write(mpl%info,'(a7,a)',advance='no') '','Compute horizontal subset C1: '
   call flush(mpl%info)
   select case (trim(nam%draw_type))
   case ('random_uniform')
      ! Random draw
      do ic0=1,geom%nc0
         if (geom%mask_hor_c0(ic0)) hdata%rh_c0(ic0,1) = 1.0
      end do
   case ('random_coast')
      ! More points around coasts
      do ic0=1,geom%nc0
         if (geom%mask_hor_c0(ic0)) hdata%rh_c0(ic0,1) = 0.0
      end do
      do il0=1,geom%nl0
         call kdtree%create(mpl,geom%nc0,geom%lon,geom%lat,mask=.not.geom%mask_c0(:,il0))
         do ic0=1,geom%nc0
            if (geom%mask_c0(ic0,il0)) then
               call kdtree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),1,nn_index,nn_dist)
               hdata%rh_c0(ic0,1) = hdata%rh_c0(ic0,1)+exp(-nn_dist(1)/Lcoast)
            else
               hdata%rh_c0(ic0,1) = hdata%rh_c0(ic0,1)+1.0
            end if
         end do
         call kdtree%dealloc
      end do
      hdata%rh_c0(:,1) = rcoast+(1.0-rcoast)*(1.0-hdata%rh_c0(:,1)/real(geom%nl0,kind_real))
   end select

   ! Initialize sampling
   call rng%initialize_sampling(mpl,geom%nc0,geom%lon,geom%lat,geom%mask_hor_c0,hdata%rh_c0(:,1),nam%ntry,nam%nrep, &
 & nam%nc1,hdata%c1_to_c0)
else
   ic1 = 0
   do ic0=1,geom%nc0
      if (geom%mask_hor_c0(ic0)) then
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
subroutine hdata_compute_sampling_ps(hdata,mpl,rng,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(rng_type),intent(inout) :: rng      !< Random number generator
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: irmaxloc,jc3,ic1,ir,ic0,jc0,i,nvc0,ivc0,icinf,icsup,ictest
integer,allocatable :: vic0(:)
real(kind_real) :: d
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical :: found

! First class
hdata%c1c3_to_c0(:,1) = hdata%c1_to_c0

if (nam%nc3>1) then
   write(mpl%info,'(a7,a)',advance='no') '','Compute positive separation sampling: '
   call flush(mpl%info)

   ! Initialize
   do jc3=1,nam%nc3
      if (jc3/=1) call msi(hdata%c1c3_to_c0(:,jc3))
   end do

   ! Define valid nodes vector
   nvc0 = count(geom%mask_hor_c0)
   allocate(vic0(nvc0))
   ivc0 = 0
   do ic0=1,geom%nc0
      if (geom%mask_hor_c0(ic0)) then
         ivc0 = ivc0+1
         vic0(ivc0) = ic0
      end if
   end do

   ! Sample classes of positive separation
   call mpl%prog_init(nam%nc3*nam%nc1)
   ir = 0
   irmaxloc = irmax
   do while ((.not.all(isnotmsi(hdata%c1c3_to_c0))).and.(nvc0>1).and.(ir<=irmaxloc))
      ! Try a random point
      if (mpl%main) call rng%rand_integer(1,nvc0,i)
      call mpl%bcast(i)
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
               if ((jc3/=1).and.(ismsi(hdata%c1c3_to_c0(ic1,jc3)))) hdata%c1c3_to_c0(ic1,jc3) = jc0
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
      mpl%done = pack(isnotmsi(hdata%c1c3_to_c0),mask=.true.)
      call mpl%prog_print
   end do
   write(mpl%info,'(a)') '100%'
   call flush(mpl%info)

   ! Release memory
   deallocate(vic0)
end if

end subroutine hdata_compute_sampling_ps

!----------------------------------------------------------------------
! Subroutine: hdata_compute_sampling_lct
!> Purpose: compute LCT sampling
!----------------------------------------------------------------------
subroutine hdata_compute_sampling_lct(hdata,mpl,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: i,il0,ic1,ic0,jc0,ibnd,ic3
integer :: nn_index(nam%nc3)
integer :: iproc,ic1_s(mpl%nproc),ic1_e(mpl%nproc),nc1_loc(mpl%nproc),ic1_loc
integer,allocatable :: sbufi(:),rbufi(:)
real(kind_real) :: nn_dist(nam%nc3)
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: sbufl(:),rbufl(:)

write(mpl%info,'(a7,a)',advance='no') '','Compute LCT sampling: '
call flush(mpl%info)

! MPI splitting
call mpl%split(nam%nc1,ic1_s,ic1_e,nc1_loc)

! Initialization
call mpl%prog_init(nc1_loc(mpl%myproc))

do ic1_loc=1,nc1_loc(mpl%myproc)
   ! MPI offset
   ic1 = ic1_s(mpl%myproc)+ic1_loc-1

   ! Check location validity
   if (isnotmsi(hdata%c1_to_c0(ic1))) then
      ! Find neighbors
      call geom%kdtree%find_nearest_neighbors(geom%lon(hdata%c1_to_c0(ic1)),geom%lat(hdata%c1_to_c0(ic1)), &
    & nam%nc3,nn_index,nn_dist)

      ! Copy neighbor index
      do ic3=1,nam%nc3
         jc0 = nn_index(ic3)
         hdata%c1c3_to_c0(ic1,ic3) = nn_index(ic3)
         do il0=1,geom%nl0
            hdata%c1c3l0_log(ic1,ic3,il0) = geom%mask_c0(jc0,il0)
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

   ! Update
   call mpl%prog_print(ic1_loc)
end do
write(mpl%info,'(a)') '100%'
call flush(mpl%info)

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
call mpl%update_tag(2)

! Broadcast data
call mpl%bcast(hdata%c1c3_to_c0)
call mpl%bcast(hdata%c1c3l0_log)

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
   hdata%c1l0_log(:,il0) = geom%mask_c0(hdata%c1_to_c0,il0)
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
            valid = geom%mask_c0(ic0,il0).and.geom%mask_c0(jc0,il0)

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
subroutine hdata_compute_mpi_a(hdata,mpl,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(in) :: mpl         !< MPI data
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
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%info,'(a10,a,i8)') '','nc0 =        ',geom%nc0
write(mpl%info,'(a10,a,i8)') '','nc0a =       ',geom%nc0a
write(mpl%info,'(a10,a,i8)') '','nl0 =        ',geom%nl0
write(mpl%info,'(a10,a,i8)') '','nc1 =        ',nam%nc1
write(mpl%info,'(a10,a,i8)') '','nc1a =       ',hdata%nc1a
call flush(mpl%info)

end subroutine hdata_compute_mpi_a

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_ab
!> Purpose: compute HDIAG MPI distribution, halos A-B
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_ab(hdata,mpl,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: iproc,ic0,ic2a,ic2b,ic2,jc2,i_s,i_s_loc,h_n_s_max,il0i,h_n_s_max_loc
integer,allocatable :: interph_lg(:,:)

! Allocation
h_n_s_max = 0
do il0i=1,geom%nl0i
   h_n_s_max = max(h_n_s_max,hdata%hfull(il0i)%n_s)
end do
allocate(hdata%c2_to_proc(nam%nc2))
allocate(hdata%proc_to_nc2a(mpl%nproc))
allocate(hdata%h(geom%nl0i))
allocate(hdata%lcheck_c2a(nam%nc2))
allocate(hdata%lcheck_c2b(nam%nc2))
allocate(hdata%lcheck_h(h_n_s_max,geom%nl0i))

! Halo definitions

! Halo A
hdata%lcheck_c2a = .false.
do ic2=1,nam%nc2
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

! Halo A
allocate(hdata%c2a_to_c2(hdata%nc2a))
allocate(hdata%c2_to_c2a(nam%nc2))
ic2a = 0
do ic2=1,nam%nc2
   if (hdata%lcheck_c2a(ic2)) then
      ic2a = ic2a+1
      hdata%c2a_to_c2(ic2a) = ic2
   end if
end do
call mpl%glb_to_loc_index(hdata%nc2a,hdata%c2a_to_c2,nam%nc2,hdata%c2_to_c2a)

! Halo B
allocate(hdata%c2b_to_c2(hdata%nc2b))
allocate(hdata%c2_to_c2b(nam%nc2))
ic2b = 0
do ic2=1,nam%nc2
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
   call hdata%h(il0i)%reorder(mpl)
end do

! MPI splitting
do ic2=1,nam%nc2
   ic0 = hdata%c2_to_c0(ic2)
   hdata%c2_to_proc(ic2) = geom%c0_to_proc(ic0)
end do
do iproc=1,mpl%nproc
   hdata%proc_to_nc2a(iproc) = count(hdata%c2_to_proc==iproc)
end do

! Setup communications
call hdata%com_AB%setup(mpl,'com_AB',nam%nc2,hdata%nc2a,hdata%nc2b,hdata%c2b_to_c2,hdata%c2a_to_c2b,hdata%c2_to_proc, &
 & hdata%c2_to_c2a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%info,'(a10,a,i8)') '','nc2 =        ',nam%nc2
write(mpl%info,'(a10,a,i8)') '','nc2a =       ',hdata%nc2a
write(mpl%info,'(a10,a,i8)') '','nc2b =       ',hdata%nc2b
do il0i=1,geom%nl0i
   write(mpl%info,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',hdata%h(il0i)%n_s
end do
call flush(mpl%info)

end subroutine hdata_compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_d
!> Purpose: compute HDIAG MPI distribution, halo D
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_d(hdata,mpl,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: ic0,ic0a,ic0d,ic1,ic2,ic2a,jc0,jc1

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
         if (hdata%displ_mask(jc1,ic2)) then
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
call hdata%com_AD%setup(mpl,'com_AD',geom%nc0,geom%nc0a,hdata%nc0d,hdata%c0d_to_c0,hdata%c0a_to_c0d,geom%c0_to_proc,geom%c0_to_c0a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%info,'(a10,a,i8)') '','nc0d =       ',hdata%nc0d
call flush(mpl%info)

end subroutine hdata_compute_mpi_d

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_c
!> Purpose: compute HDIAG MPI distribution, halo C
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_c(hdata,mpl,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: jc3,ic0,ic0a,ic0c,ic1,ic1a,its,il0,d_n_s_max,d_n_s_max_loc,i_s,i_s_loc
integer,allocatable :: interpd_lg(:,:,:)
real(kind_real),allocatable :: displ_lon(:,:),displ_lat(:,:),lon_c1(:),lat_c1(:)
type(linop_type),allocatable :: dfull(:,:)

if (nam%displ_diag) then
   ! Allocation
   allocate(dfull(geom%nl0,nam%nts))
   allocate(displ_lon(geom%nc0,geom%nl0))
   allocate(displ_lat(geom%nc0,geom%nl0))
   allocate(lon_c1(nam%nc1))
   allocate(lat_c1(nam%nc1))

   ! Prepare displacement interpolation
   do its=1,nam%nts
      ! Local to global
      call mpl%loc_to_glb(geom%nl0,geom%nc0a,hdata%displ_lon(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,displ_lon)
      call mpl%loc_to_glb(geom%nl0,geom%nc0a,hdata%displ_lat(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,displ_lat)

      do il0=1,geom%nl0
         ! Copy Sc1 points
         do ic1=1,nam%nc1
            ic0 = hdata%c1_to_c0(ic1)
            lon_c1(ic1) = displ_lon(ic0,il0)
            lat_c1(ic1) = displ_lat(ic0,il0)
         end do

         ! Compute interpolation
         call dfull(il0,its)%interp(mpl,geom%mesh,geom%kdtree,geom%nc0,geom%mask_c0(:,il0),nam%nc1,lon_c1,lat_c1, &
       & hdata%c1l0_log(:,il0),nam%diag_interp)
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
         call hdata%d(il0,its)%reorder(mpl)
      end do
   end do
end if

! Setup communications
call hdata%com_AC%setup(mpl,'com_AC',geom%nc0,geom%nc0a,hdata%nc0c,hdata%c0c_to_c0,hdata%c0a_to_c0c,geom%c0_to_proc,geom%c0_to_c0a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%info,'(a10,a,i8)') '','nc0c =      ',hdata%nc0c
call flush(mpl%info)

end subroutine hdata_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: hdata_compute_mpi_f
!> Purpose: compute HDIAG MPI distribution, halo F
!----------------------------------------------------------------------
subroutine hdata_compute_mpi_f(hdata,mpl,nam,geom)

implicit none

! Passed variables
class(hdata_type),intent(inout) :: hdata !< HDIAG data
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry

! Local variables
integer :: ic2,ic1,jc2,ic2a,ic2f,il0,kc2

! Allocation
allocate(hdata%lcheck_c2f(nam%nc2))

! Halo definitions

! Halo F
do il0=1,geom%nl0
   hdata%lcheck_c2f = hdata%lcheck_c2a
   do ic2a=1,hdata%nc2a
      ic2 = hdata%c2a_to_c2(ic2a)
      ic1 = hdata%c2_to_c1(ic2)
      jc2 = 1
      do while (isnotmsi(hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i))) &
    & .and.(hdata%nn_c2_dist(jc2,ic2,min(il0,geom%nl0i))<nam%diag_rhflt))
          kc2 = hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i))
          hdata%lcheck_c2f(kc2) = .true.
          jc2 = jc2+1
          if (jc2>nam%nc2) exit
      end do
   end do
end do
hdata%nc2f = count(hdata%lcheck_c2f)

! Global <-> local conversions for fields

! Halo F
allocate(hdata%c2f_to_c2(hdata%nc2f))
allocate(hdata%c2_to_c2f(nam%nc2))
call msi(hdata%c2_to_c2f)
ic2f = 0
do ic2=1,nam%nc2
   if (hdata%lcheck_c2f(ic2)) then
      ic2f = ic2f+1
      hdata%c2f_to_c2(ic2f) = ic2
      hdata%c2_to_c2f(ic2) = ic2f
   end if
end do

! Inter-halo conversions
allocate(hdata%c2a_to_c2f(hdata%nc2a))
do ic2a=1,hdata%nc2a
   ic2 = hdata%c2a_to_c2(ic2a)
   ic2f = hdata%c2_to_c2f(ic2)
   hdata%c2a_to_c2f(ic2a) = ic2f
end do

! Setup communications
call hdata%com_AF%setup(mpl,'com_AF',nam%nc2,hdata%nc2a,hdata%nc2f,hdata%c2f_to_c2,hdata%c2a_to_c2f,hdata%c2_to_proc, &
 & hdata%c2_to_c2a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%info,'(a10,a,i8)') '','nc2f =       ',hdata%nc2f
call flush(mpl%info)

end subroutine hdata_compute_mpi_f

!----------------------------------------------------------------------
! Subroutine: hdata_diag_filter
!> Purpose: filter diagnostics
!----------------------------------------------------------------------
subroutine hdata_diag_filter(hdata,mpl,nam,geom,il0,filter_type,rflt,diag)

implicit none

! Passed variables
class(hdata_type),intent(in) :: hdata             !< HDIAG data
type(mpl_type),intent(inout) :: mpl               !< MPI data
type(nam_type),intent(in) :: nam                  !< Namelist
type(geom_type),intent(in) :: geom                !< Geometry
integer,intent(in) :: il0                         !< Level
character(len=*),intent(in) :: filter_type        !< Filter type
real(kind_real),intent(in) :: rflt                !< Filter support radius
real(kind_real),intent(inout) :: diag(hdata%nc2a) !< Filtered diagnostic

! Local variables
integer :: ic2a,ic2,ic1,jc2,nc2eff,kc2,kc2_glb
integer,allocatable :: order(:)
real(kind_real) :: distnorm,norm,wgt
real(kind_real),allocatable :: diag_glb(:),diag_eff(:),diag_eff_dist(:)
logical :: nam_rad

! Check radius
nam_rad = .not.(abs(rflt-nam%diag_rhflt)>0.0)

if (rflt>0.0) then
   if (nam_rad) then
      ! Allocation
      allocate(diag_glb(hdata%nc2f))

      ! Communication
      call hdata%com_AF%ext(mpl,diag,diag_glb)
   else
      ! Allocation
      allocate(diag_glb(nam%nc2))

      ! Local to global
     call mpl%loc_to_glb(hdata%nc2a,diag,nam%nc2,hdata%c2_to_proc,hdata%c2_to_c2a,.true.,diag_glb)
   end if

   !$omp parallel do schedule(static) private(ic2a,ic2,ic1,nc2eff,jc2,distnorm,norm,wgt) firstprivate(diag_eff,diag_eff_dist,order)
   do ic2a=1,hdata%nc2a
      ic2 = hdata%c2a_to_c2(ic2a)
      ic1 = hdata%c2_to_c1(ic2)
      if (hdata%c1l0_log(ic1,il0)) then
         ! Allocation
         allocate(diag_eff(nam%nc2))
         allocate(diag_eff_dist(nam%nc2))

         ! Build diag_eff of valid points
         nc2eff = 0
         jc2 = 1
         do while (isnotmsi(hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i))).and.(hdata%nn_c2_dist(jc2,ic2,min(il0,geom%nl0i))<rflt))
            ! Check the point validity
            kc2 = hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i))
            if (nam_rad) then
               kc2_glb = hdata%c2_to_c2f(kc2)
            else
               kc2_glb = kc2
            end if
            if (isnotmsr(diag_glb(kc2_glb))) then
               nc2eff = nc2eff+1
               diag_eff(nc2eff) = diag_glb(kc2_glb)
               diag_eff_dist(nc2eff) = hdata%nn_c2_dist(jc2,ic2,min(il0,geom%nl0i))
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
               call mpl%abort('wrong filter type')
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
end if

end subroutine hdata_diag_filter

!----------------------------------------------------------------------
! Subroutine: hdata_diag_fill
!> Purpose: fill diagnostics missing values
!----------------------------------------------------------------------
subroutine hdata_diag_fill(hdata,mpl,nam,geom,il0,diag)

implicit none

! Passed variables
class(hdata_type),intent(in) :: hdata             !< HDIAG data
type(mpl_type),intent(inout) :: mpl               !< MPI data
type(nam_type),intent(in) :: nam                  !< Namelist
type(geom_type),intent(in) :: geom                !< Geometry
integer,intent(in) :: il0                         !< Level
real(kind_real),intent(inout) :: diag(hdata%nc2a) !< Filtered diagnostic

! Local variables
integer :: nmsr,nmsr_tot,ic2,jc2,kc2
real(kind_real),allocatable :: diag_glb(:)

! Count missing points
if (hdata%nc2a>0) then
   nmsr = count(ismsr(diag))
else
   nmsr = 0
end if
call mpl%allreduce_sum(nmsr,nmsr_tot)

if (nmsr_tot>0) then
   ! Allocation
   if (mpl%main) allocate(diag_glb(nam%nc2))

   ! Local to global
   call mpl%loc_to_glb(hdata%nc2a,diag,nam%nc2,hdata%c2_to_proc,hdata%c2_to_c2a,.false.,diag_glb)

   if (mpl%main) then
      do ic2=1,nam%nc2
         jc2 = 1
         do while (ismsr(diag_glb(ic2)))
            kc2 = hdata%nn_c2_index(jc2,ic2,min(il0,geom%nl0i))
            if (isnotmsr(diag_glb(kc2))) diag_glb(ic2) = diag_glb(kc2)
            jc2 = jc2+1
            if (jc2>nam%nc2) exit
         end do
      end do
   end if

   ! Global to local
   call mpl%glb_to_loc(nam%nc2,hdata%c2_to_proc,hdata%c2_to_c2a,diag_glb,hdata%nc2a,diag)
end if

end subroutine hdata_diag_fill

end module type_hdata
