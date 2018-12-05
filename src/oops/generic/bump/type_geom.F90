!----------------------------------------------------------------------
! Module: type_geom
! Purpose: geometry derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_geom

use netcdf
use tools_const, only: pi,req,deg2rad,rad2deg,reqkm
use tools_func, only: lonlatmod,sphere_dist,vector_product,vector_triple_product
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,ismsi,ismsr,isanynotmsr
use tools_nc, only: ncfloat
use tools_qsort, only: qsort
use tools_stripack, only: areas,trans
use type_com, only: com_type
use type_kdtree, only: kdtree_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type
use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max,fckit_mpi_status

implicit none

integer,parameter :: nredmax = 10                 ! Maximum number of similar redundant points
real(kind_real),parameter :: distminred = 1.0e-12 ! Minimum distance (in rad) to consider points as redundant
logical,parameter :: test_no_point = .false.      ! Test BUMP with no grid point on the last MPI task

! Geometry derived type
type geom_type
   ! Offline geometry data
   integer :: nlon                            ! Longitude size
   integer :: nlat                            ! Latitude size
   integer :: nlev                            ! Number of levels
   integer,allocatable :: c0_to_lon(:)        ! Subset Sc0 to longitude index
   integer,allocatable :: c0_to_lat(:)        ! Subset Sc0 to latgitude index
   integer,allocatable :: c0_to_tile(:)       ! Subset Sc0 to tile index

   ! Number of points and levels
   integer :: nmg                             ! Number of model grid points
   integer :: nc0                             ! Number of points in subset Sc0
   integer :: nl0                             ! Number of levels in subset Sl0f
   integer :: nl0i                            ! Number of independent levels in subset Sl0

   ! Basic geometry data
   real(kind_real),allocatable :: lon(:)      ! Longitudes
   real(kind_real),allocatable :: lat(:)      ! Latitudes
   real(kind_real),allocatable :: area(:)     ! Domain area
   real(kind_real),allocatable :: vunit(:,:)  ! Vertical unit
   real(kind_real),allocatable :: vunitavg(:) ! Averaged vertical unit
   real(kind_real),allocatable :: disth(:)    ! Horizontal distance

   ! Masks
   logical,allocatable :: mask_hor_mg(:)      ! Horizontal mask on model grid, global
   logical,allocatable :: mask_c0(:,:)        ! Mask on subset Sc0, global
   logical,allocatable :: mask_c0a(:,:)       ! Mask on subset Sc0, halo A
   logical,allocatable :: mask_hor_c0(:)      ! Union of horizontal masks on subset Sc0, global
   logical,allocatable :: mask_hor_c0a(:)     ! Union of horizontal masks on subset Sc0, halo A
   logical,allocatable :: mask_ver_c0(:)      ! Union of vertical masks
   integer,allocatable :: nc0_mask(:)         ! Horizontal mask size on subset Sc0
   logical :: mask_del                        ! Remove subset Sc0 points that are masked for all level

   ! Mesh
   type(mesh_type) :: mesh                    ! Mesh

   ! KD-tree
   type(kdtree_type) :: kdtree                ! KD-tree

   ! Boundary nodes
   integer,allocatable :: nbnd(:)             ! Number of boundary nodes
   real(kind_real),allocatable :: xbnd(:,:,:) ! Boundary nodes, x-coordinate
   real(kind_real),allocatable :: ybnd(:,:,:) ! Boundary nodes, y-coordinate
   real(kind_real),allocatable :: zbnd(:,:,:) ! Boundary nodes, z-coordinate
   real(kind_real),allocatable :: vbnd(:,:,:) ! Boundary nodes, orthogonal vector

   ! Gripoints and subset Sc0
   integer,allocatable :: redundant(:)        ! Redundant points array
   integer,allocatable :: c0_to_mg(:)         ! Subset Sc0 to model grid
   integer,allocatable :: mg_to_c0(:)         ! Model grid to subset Sc0

   ! Dirac information
   integer :: ndir                            ! Number of valid Dirac points
   real(kind_real),allocatable :: londir(:)   ! Dirac longitude
   real(kind_real),allocatable :: latdir(:)   ! Dirac latitude
   integer,allocatable :: iprocdir(:)         ! Dirac processor
   integer,allocatable :: ic0adir(:)          ! Dirac gridpoint
   integer,allocatable :: il0dir(:)           ! Dirac level
   integer,allocatable :: ivdir(:)            ! Dirac variable
   integer,allocatable :: itsdir(:)           ! Dirac timeslot

   ! MPI distribution
   integer :: nmga                            ! Halo A size for model grid
   integer :: nc0a                            ! Halo A size for subset Sc0
   integer,allocatable :: mg_to_proc(:)       ! Model grid to local task
   integer,allocatable :: mg_to_mga(:)        ! Model grid, global to halo A
   integer,allocatable :: mga_to_mg(:)        ! Model grid, halo A to global
   integer,allocatable :: proc_to_nmga(:)     ! Halo A size for each proc
   integer,allocatable :: c0_to_proc(:)       ! Subset Sc0 to local task
   integer,allocatable :: c0_to_c0a(:)        ! Subset Sc0, global to halo A
   integer,allocatable :: c0a_to_c0(:)        ! Subset Sc0, halo A to global
   integer,allocatable :: proc_to_nc0a(:)     ! Halo A size for each proc
   integer,allocatable :: mga_to_c0(:)        ! Model grid, halo A to subset Sc0, global
   integer,allocatable :: c0a_to_mga(:)       ! Subset Sc0 to model grid, halo A
   type(com_type) :: com_mg                   ! Communication between subset Sc0 and model grid
contains
   procedure :: alloc => geom_alloc
   procedure :: dealloc => geom_dealloc
   procedure :: setup_online => geom_setup_online
   procedure :: find_sc0 => geom_find_sc0
   procedure :: init => geom_init
   procedure :: compute_area => geom_compute_area
   procedure :: define_dirac => geom_define_dirac
   procedure :: define_distribution => geom_define_distribution
   procedure :: check_arc => geom_check_arc
   procedure :: copy_c0a_to_mga => geom_copy_c0a_to_mga
   procedure :: copy_mga_to_c0a => geom_copy_mga_to_c0a
   procedure :: compute_deltas => geom_compute_deltas
end type geom_type

private
public :: geom_type

contains

!----------------------------------------------------------------------
! Subroutine: geom_alloc
! Purpose: geometry allocation
!----------------------------------------------------------------------
subroutine geom_alloc(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry

! Allocation
allocate(geom%c0_to_proc(geom%nc0))
allocate(geom%c0_to_c0a(geom%nc0))
allocate(geom%c0_to_lon(geom%nc0))
allocate(geom%c0_to_lat(geom%nc0))
allocate(geom%c0_to_tile(geom%nc0))
allocate(geom%lon(geom%nc0))
allocate(geom%lat(geom%nc0))
allocate(geom%area(geom%nl0))
allocate(geom%vunit(geom%nc0,geom%nl0))
allocate(geom%vunitavg(geom%nl0))
allocate(geom%mask_c0(geom%nc0,geom%nl0))
allocate(geom%mask_hor_c0(geom%nc0))
allocate(geom%mask_ver_c0(geom%nl0))
allocate(geom%nc0_mask(geom%nl0))

! Initialization
call msi(geom%c0_to_proc)
call msi(geom%c0_to_c0a)
call msi(geom%c0_to_lon)
call msi(geom%c0_to_lat)
call msi(geom%c0_to_tile)
call msr(geom%lon)
call msr(geom%lat)
call msr(geom%area)
call msr(geom%vunit)
call msr(geom%vunitavg)
geom%mask_c0 = .false.
geom%mask_hor_c0 = .false.
geom%mask_ver_c0 = .false.
call msi(geom%nc0_mask)

end subroutine geom_alloc

!----------------------------------------------------------------------
! Subroutine: geom_dealloc
! Purpose: geometry deallocation
!----------------------------------------------------------------------
subroutine geom_dealloc(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry

! Release memory
if (allocated(geom%c0_to_lon)) deallocate(geom%c0_to_lon)
if (allocated(geom%c0_to_lat)) deallocate(geom%c0_to_lat)
if (allocated(geom%c0_to_tile)) deallocate(geom%c0_to_tile)
if (allocated(geom%lon)) deallocate(geom%lon)
if (allocated(geom%lat)) deallocate(geom%lat)
if (allocated(geom%area)) deallocate(geom%area)
if (allocated(geom%vunit)) deallocate(geom%vunit)
if (allocated(geom%vunitavg)) deallocate(geom%vunitavg)
if (allocated(geom%disth)) deallocate(geom%disth)
if (allocated(geom%mask_c0)) deallocate(geom%mask_c0)
if (allocated(geom%mask_c0a)) deallocate(geom%mask_c0a)
if (allocated(geom%mask_hor_c0)) deallocate(geom%mask_hor_c0)
if (allocated(geom%mask_hor_c0a)) deallocate(geom%mask_hor_c0a)
if (allocated(geom%mask_ver_c0)) deallocate(geom%mask_ver_c0)
if (allocated(geom%nc0_mask)) deallocate(geom%nc0_mask)
call geom%mesh%dealloc
call geom%kdtree%dealloc
if (allocated(geom%redundant)) deallocate(geom%redundant)
if (allocated(geom%nbnd)) deallocate(geom%nbnd)
if (allocated(geom%xbnd)) deallocate(geom%xbnd)
if (allocated(geom%ybnd)) deallocate(geom%ybnd)
if (allocated(geom%zbnd)) deallocate(geom%zbnd)
if (allocated(geom%vbnd)) deallocate(geom%vbnd)
if (allocated(geom%c0_to_mg)) deallocate(geom%c0_to_mg)
if (allocated(geom%mg_to_c0)) deallocate(geom%mg_to_c0)
if (allocated(geom%mg_to_proc)) deallocate(geom%mg_to_proc)
if (allocated(geom%mg_to_mga)) deallocate(geom%mg_to_mga)
if (allocated(geom%mga_to_mg)) deallocate(geom%mga_to_mg)
if (allocated(geom%proc_to_nmga)) deallocate(geom%proc_to_nmga)
if (allocated(geom%c0_to_proc)) deallocate(geom%c0_to_proc)
if (allocated(geom%c0_to_c0a)) deallocate(geom%c0_to_c0a)
if (allocated(geom%c0a_to_c0)) deallocate(geom%c0a_to_c0)
if (allocated(geom%proc_to_nc0a)) deallocate(geom%proc_to_nc0a)
if (allocated(geom%mga_to_c0)) deallocate(geom%mga_to_c0)
if (allocated(geom%c0a_to_mga)) deallocate(geom%c0a_to_mga)
call geom%com_mg%dealloc

end subroutine geom_dealloc

!----------------------------------------------------------------------
! Subroutine: geom_setup_online
! Purpose: setup online geometry
!----------------------------------------------------------------------
subroutine geom_setup_online(geom,mpl,rng,nam,nmga,nl0,lon,lat,area,vunit,lmask)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom        ! Geometry
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(rng_type),intent(inout) :: rng           ! Random number generator
type(nam_type),intent(in) :: nam              ! Namelist
integer,intent(in) :: nmga                    ! Halo A size
integer,intent(in) :: nl0                     ! Number of levels in subset Sl0
real(kind_real),intent(in) :: lon(nmga)       ! Longitudes
real(kind_real),intent(in) :: lat(nmga)       ! Latitudes
real(kind_real),intent(in) :: area(nmga)      ! Area
real(kind_real),intent(in) :: vunit(nmga,nl0) ! Vertical unit
logical,intent(in) :: lmask(nmga,nl0)         ! Mask

! Local variables
integer :: ic0,ic0a,il0,offset,iproc,img,imga
integer,allocatable :: order(:),order_inv(:)
real(kind_real),allocatable :: lon_mg(:),lat_mg(:),area_mg(:),vunit_mg(:,:),list(:)
logical,allocatable :: lmask_mg(:,:)
type(fckit_mpi_status) :: status

! Copy geometry variables
geom%nmga = nmga
geom%nl0 = nl0
geom%nlev = nl0

! Allocation
allocate(geom%proc_to_nmga(mpl%nproc))

! Communication
call mpl%f_comm%allgather(geom%nmga,geom%proc_to_nmga)

! Global number of model grid points
geom%nmg = sum(geom%proc_to_nmga)

! Allocation
allocate(lon_mg(geom%nmg))
allocate(lat_mg(geom%nmg))
allocate(area_mg(geom%nmg))
allocate(vunit_mg(geom%nmg,geom%nl0))
allocate(lmask_mg(geom%nmg,geom%nl0))
allocate(geom%mg_to_proc(geom%nmg))
allocate(geom%mg_to_mga(geom%nmg))
allocate(geom%mga_to_mg(geom%nmga))

! Communication of model grid points
if (mpl%main) then
   ! Allocation
   offset = 0
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         lon_mg(offset+1:offset+geom%proc_to_nmga(iproc)) = lon
         lat_mg(offset+1:offset+geom%proc_to_nmga(iproc)) = lat
         area_mg(offset+1:offset+geom%proc_to_nmga(iproc)) = area
         do il0=1,geom%nl0
            vunit_mg(offset+1:offset+geom%proc_to_nmga(iproc),il0) = vunit(:,il0)
            lmask_mg(offset+1:offset+geom%proc_to_nmga(iproc),il0) = lmask(:,il0)
         end do
      else
         ! Receive data on ioproc
         call mpl%f_comm%receive(lon_mg(offset+1:offset+geom%proc_to_nmga(iproc)),iproc-1,mpl%tag,status)
         call mpl%f_comm%receive(lat_mg(offset+1:offset+geom%proc_to_nmga(iproc)),iproc-1,mpl%tag+1,status)
         call mpl%f_comm%receive(area_mg(offset+1:offset+geom%proc_to_nmga(iproc)),iproc-1,mpl%tag+2,status)
         do il0=1,geom%nl0
            call mpl%f_comm%receive(vunit_mg(offset+1:offset+geom%proc_to_nmga(iproc),il0),iproc-1,mpl%tag+2+il0,status)
            call mpl%f_comm%receive(lmask_mg(offset+1:offset+geom%proc_to_nmga(iproc),il0),iproc-1, &
          & mpl%tag+2+geom%nl0+il0,status)
         end do
      end if

      !  Update offset
      offset = offset+geom%proc_to_nmga(iproc)
   end do
else
   ! Send data to ioproc
   call mpl%f_comm%send(lon,mpl%ioproc-1,mpl%tag)
   call mpl%f_comm%send(lat,mpl%ioproc-1,mpl%tag+1)
   call mpl%f_comm%send(area,mpl%ioproc-1,mpl%tag+2)
   do il0=1,geom%nl0
      call mpl%f_comm%send(vunit(:,il0),mpl%ioproc-1,mpl%tag+2+il0)
      call mpl%f_comm%send(lmask(:,il0),mpl%ioproc-1,mpl%tag+2+geom%nl0+il0)
   end do
end if
call mpl%update_tag(3+2*geom%nl0)

if (mpl%main) then
   ! Convert to radians
   lon_mg = lon_mg*deg2rad
   lat_mg = lat_mg*deg2rad
end if

! Broadcast data
call mpl%f_comm%broadcast(lon_mg,mpl%ioproc-1)
call mpl%f_comm%broadcast(lat_mg,mpl%ioproc-1)
call mpl%f_comm%broadcast(area_mg,mpl%ioproc-1)
call mpl%f_comm%broadcast(vunit_mg,mpl%ioproc-1)
call mpl%f_comm%broadcast(lmask_mg,mpl%ioproc-1)

! Find subset Sc0 points
call geom%find_sc0(mpl,rng,lon_mg,lat_mg,lmask_mg,.true.,nam%mask_check,.false.)

! Allocation
call geom%alloc
allocate(geom%proc_to_nc0a(mpl%nproc))
allocate(list(geom%nc0))
allocate(order(geom%nc0))
allocate(order_inv(geom%nc0))

! Model grid conversions and Sc0 size on halo A
img = 0
geom%proc_to_nc0a = 0
do iproc=1,mpl%nproc
   do imga=1,geom%proc_to_nmga(iproc)
      img = img+1
      geom%mg_to_proc(img) = iproc
      geom%mg_to_mga(img) = imga
      if (iproc==mpl%myproc) geom%mga_to_mg(imga) = img
      if (geom%mask_hor_mg(img)) geom%proc_to_nc0a(iproc) = geom%proc_to_nc0a(iproc)+1
   end do
end do
geom%nc0a = geom%proc_to_nc0a(mpl%myproc)

! Subset Sc0 conversions
allocate(geom%c0a_to_c0(geom%nc0a))
ic0 = 0
do iproc=1,mpl%nproc
   do ic0a=1,geom%proc_to_nc0a(iproc)
      ic0 = ic0+1
      geom%c0_to_proc(ic0) = iproc
      if (iproc==mpl%myproc) geom%c0a_to_c0(ic0a) = ic0
   end do
end do
call mpl%glb_to_loc_index(geom%nc0a,geom%c0a_to_c0,geom%nc0,geom%c0_to_c0a)

! Inter-halo conversions
allocate(geom%mga_to_c0(geom%nmga))
allocate(geom%c0a_to_mga(geom%nc0a))
do imga=1,geom%nmga
   img = geom%mga_to_mg(imga)
   ic0 = geom%mg_to_c0(img)
   geom%mga_to_c0(imga) = ic0
end do
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   img = geom%c0_to_mg(ic0)
   imga = geom%mg_to_mga(img)
   geom%c0a_to_mga(ic0a) = imga
end do

! Deal with mask on redundant points
do il0=1,geom%nl0
   do img=1,geom%nmg
      if (isnotmsi(geom%redundant(img))) lmask_mg(img,il0) = lmask_mg(img,il0).or.lmask_mg(geom%redundant(img),il0)
   end do
end do

! Remove redundant points
geom%lon = lon_mg(geom%c0_to_mg)
geom%lat = lat_mg(geom%c0_to_mg)
do il0=1,geom%nl0
   geom%area(il0) = sum(area_mg(geom%c0_to_mg),lmask_mg(geom%c0_to_mg,il0))/req**2
   geom%vunit(:,il0) = vunit_mg(geom%c0_to_mg,il0)
   geom%mask_c0(:,il0) = lmask_mg(geom%c0_to_mg,il0)
end do

! Reorder points based on lon/lat
do ic0=1,geom%nc0
   list(ic0) = aint(abs(geom%lon(ic0))*1.0e6)+abs(geom%lat(ic0))*1.0e-1
   if (geom%lon(ic0)<0.0) list(ic0) = list(ic0)+2.0e7
   if (geom%lat(ic0)<0.0) list(ic0) = list(ic0)+1.0e7
end do
call qsort(geom%nc0,list,order)
do ic0=1,geom%nc0
   order_inv(order(ic0)) = ic0
end do
geom%c0_to_proc = geom%c0_to_proc(order)
geom%c0_to_c0a = geom%c0_to_c0a(order)
geom%c0a_to_c0 = order_inv(geom%c0a_to_c0)
geom%lon = geom%lon(order)
geom%lat = geom%lat(order)
do il0=1,geom%nl0
   geom%vunit(:,il0) = geom%vunit(order,il0)
   geom%mask_c0(:,il0) = geom%mask_c0(order,il0)
end do
geom%c0_to_mg = geom%c0_to_mg(order)
geom%mg_to_c0 = order_inv(geom%mg_to_c0)
geom%mga_to_c0 = order_inv(geom%mga_to_c0)

! Other masks
allocate(geom%mask_c0a(geom%nc0a,geom%nl0))
allocate(geom%mask_hor_c0a(geom%nc0a))
geom%mask_c0a = geom%mask_c0(geom%c0a_to_c0,:)
geom%mask_hor_c0 = any(geom%mask_c0,dim=2)
geom%mask_hor_c0a = geom%mask_hor_c0(geom%c0a_to_c0)
geom%mask_ver_c0 = any(geom%mask_c0,dim=1)
geom%nc0_mask = count(geom%mask_c0,dim=1)

! Setup redundant points communication
call geom%com_mg%setup(mpl,'com_mg',geom%nc0,geom%nc0a,geom%nmga,geom%mga_to_c0,geom%c0a_to_mga,geom%c0_to_proc,geom%c0_to_c0a)

end subroutine geom_setup_online

!----------------------------------------------------------------------
! Subroutine: geom_find_sc0
! Purpose: find subset Sc0 points
!----------------------------------------------------------------------
subroutine geom_find_sc0(geom,mpl,rng,lon,lat,lmask,red_check,mask_check,mask_del)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom         ! Geometry
type(mpl_type),intent(in) :: mpl               ! MPI data
type(rng_type),intent(inout) :: rng            ! Random number generator
real(kind_real),intent(in) :: lon(geom%nmg)    ! Longitudes
real(kind_real),intent(in) :: lat(geom%nmg)    ! Latitudes
logical,intent(in) :: lmask(geom%nmg,geom%nl0) ! Mask
logical,intent(in) :: red_check                ! Check redundant points
logical,intent(in) :: mask_check               ! Check mask boundaries
logical,intent(in) :: mask_del                 ! Delete masked points

! Local variables
integer :: img,ic0,ired,nn_index(nredmax),nc0full,ic0full,jc0full,kc0full,i,j,k,iend,ibnd,il0
integer,allocatable :: c0full_to_mg(:),mg_to_c0full(:),ic0full_bnd(:,:,:)
real(kind_real) :: nn_dist(nredmax),latbnd(2),lonbnd(2),v1(3),v2(3)
real(kind_real),allocatable :: lon_c0full(:),lat_c0full(:)
logical :: init
logical,allocatable :: lmask_c0full(:,:)
type(kdtree_type) :: kdtree
type(mesh_type) :: mesh

! Allocation
allocate(geom%mask_hor_mg(geom%nmg))
allocate(geom%redundant(geom%nmg))
allocate(geom%mg_to_c0(geom%nmg))

! Initialization
call msi(geom%redundant)

! Look for redundant points
if (red_check) then
   write(mpl%info,'(a7,a)') '','Look for redundant points in the model grid'

   ! Create KD-tree
   call kdtree%create(mpl,geom%nmg,lon,lat)

   ! Find redundant points
   do img=1,geom%nmg
      ! Find nearest neighbors
      call kdtree%find_nearest_neighbors(lon(img),lat(img),nredmax,nn_index,nn_dist)

      ! Count redundant points
      do ired=1,nredmax
         if ((nn_dist(ired)>distminred).or.(nn_index(ired)>=img)) nn_index(ired) = geom%nmg+1
      end do

      if (any(nn_index<=geom%nmg)) then
         ! Redundant point
         geom%redundant(img) = minval(nn_index)
      end if
   end do

   ! Check for successive redundant points
   do img=1,geom%nmg
      if (isnotmsi(geom%redundant(img))) then
         do while (isnotmsi(geom%redundant(geom%redundant(img))))
            geom%redundant(img) = geom%redundant(geom%redundant(img))
         end do
      end if
   end do

   ! Release memory
   call kdtree%dealloc
end if

! Horizontal model grid mask
geom%mask_hor_mg = ismsi(geom%redundant)

! Allocation
nc0full = count(geom%mask_hor_mg)
allocate(c0full_to_mg(nc0full))
allocate(mg_to_c0full(geom%nmg))
allocate(lon_c0full(nc0full))
allocate(lat_c0full(nc0full))
allocate(lmask_c0full(nc0full,geom%nl0))

! Conversion
ic0full = 0
do img=1,geom%nmg
   if (ismsi(geom%redundant(img))) then
      ic0full = ic0full+1
      c0full_to_mg(ic0full) = img
      mg_to_c0full(img) = ic0full
   end if
end do
lon_c0full = lon(c0full_to_mg)
lat_c0full = lat(c0full_to_mg)
lmask_c0full = lmask(c0full_to_mg,:)

! Delete points flag
geom%mask_del = mask_del.and.(any(.not.lmask_c0full))

if (geom%mask_del.or.mask_check) then
   ! Create mesh
   call mesh%create(mpl,rng,nc0full,lon_c0full,lat_c0full)

   ! Allocation
   allocate(geom%nbnd(geom%nl0))
   allocate(ic0full_bnd(2,mesh%n,geom%nl0))

   ! Find border points
   do il0=1,geom%nl0
      geom%nbnd(il0) = 0
      do i=1,mesh%n
         ! Check mask points only
         ic0full = mesh%order(i)
         if (.not.lmask_c0full(ic0full,il0)) then
            iend = mesh%lend(i)
            init = .true.
            do while ((iend/=mesh%lend(i)).or.init)
               j = abs(mesh%list(iend))
               k = abs(mesh%list(mesh%lptr(iend)))
               jc0full = mesh%order(j)
               kc0full = mesh%order(k)
               if (.not.lmask_c0full(jc0full,il0).and.lmask_c0full(kc0full,il0)) then
                  ! Create a new boundary arc
                  geom%nbnd(il0) = geom%nbnd(il0)+1
                  if (geom%nbnd(il0)>mesh%n) call mpl%abort('too many boundary arcs')
                  ic0full_bnd(1,geom%nbnd(il0),il0) = ic0full
                  ic0full_bnd(2,geom%nbnd(il0),il0) = jc0full
               end if
               iend = mesh%lptr(iend)
               init = .false.
            end do
         end if
      end do
   end do

   ! Allocation
   allocate(geom%xbnd(2,maxval(geom%nbnd),geom%nl0))
   allocate(geom%ybnd(2,maxval(geom%nbnd),geom%nl0))
   allocate(geom%zbnd(2,maxval(geom%nbnd),geom%nl0))
   allocate(geom%vbnd(3,maxval(geom%nbnd),geom%nl0))

   do il0=1,geom%nl0
      ! Compute boundary arcs
      do ibnd=1,geom%nbnd(il0)
         latbnd = lat_c0full(ic0full_bnd(:,ibnd,il0))
         lonbnd = lon_c0full(ic0full_bnd(:,ibnd,il0))
         call trans(2,latbnd,lonbnd,geom%xbnd(:,ibnd,il0),geom%ybnd(:,ibnd,il0),geom%zbnd(:,ibnd,il0))
      end do
      do ibnd=1,geom%nbnd(il0)
         v1 = (/geom%xbnd(1,ibnd,il0),geom%ybnd(1,ibnd,il0),geom%zbnd(1,ibnd,il0)/)
         v2 = (/geom%xbnd(2,ibnd,il0),geom%ybnd(2,ibnd,il0),geom%zbnd(2,ibnd,il0)/)
         call vector_product(v1,v2,geom%vbnd(:,ibnd,il0))
      end do
   end do
end if

if (geom%mask_del) then
   ! Remove subset Sc0 points that are masked for all levels
   do img=1,geom%nmg
      if (geom%mask_hor_mg(img)) then
         ic0full = mg_to_c0full(img)
         geom%mask_hor_mg(img) = any(lmask_c0full(ic0full,:))
      end if
   end do
end if

! Allocation
geom%nc0 = count(geom%mask_hor_mg)
allocate(geom%c0_to_mg(geom%nc0))

! Initialization
call msi(geom%mg_to_c0)

! Conversion
ic0 = 0
do img=1,geom%nmg
   if (geom%mask_hor_mg(img)) then
      ic0 = ic0+1
      geom%c0_to_mg(ic0) = img
      geom%mg_to_c0(img) = ic0
   end if
end do

 ! Deal with redundant points
do img=1,geom%nmg
   if (isnotmsi(geom%redundant(img))) geom%mg_to_c0(img) = geom%mg_to_c0(geom%redundant(img))
end do

! Print results
if (red_check) write(mpl%info,'(a7,a,i6,a,f6.2,a)') '','Number of redundant points:    ',(geom%nmg-nc0full), &
 & ' (',real(geom%nmg-nc0full,kind_real)/real(geom%nmg,kind_real)*100.0,'%)'
if (mask_del) write(mpl%info,'(a7,a,i6,a,f6.2,a)') '','Number of masked points       :',(nc0full-geom%nc0), &
 & ' (',real(nc0full-geom%nc0,kind_real)/real(geom%nmg,kind_real)*100.0,'%)'
if (red_check.and.mask_del) write(mpl%info,'(a7,a,i6,a,f6.2,a)') '','Total number of removed points:',(geom%nmg-geom%nc0), &
 & ' (',real(geom%nmg-geom%nc0,kind_real)/real(geom%nmg,kind_real)*100.0,'%)'
call flush(mpl%info)
write(mpl%info,'(a7,a,i8)') '','Model grid size:         ',geom%nmg
if (red_check) write(mpl%info,'(a7,a,i8)') '','Non-redundant grid size: ',nc0full
write(mpl%info,'(a7,a,i8)') '','Subset Sc0 size:         ',geom%nc0

end subroutine geom_find_sc0

!----------------------------------------------------------------------
! Subroutine: geom_init
! Purpose: initialize geometry
!----------------------------------------------------------------------
subroutine geom_init(geom,mpl,rng,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(in) :: mpl       ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: ic0,il0,jc3,iproc
logical :: same_mask

! Set longitude and latitude bounds
do ic0=1,geom%nc0
   call lonlatmod(geom%lon(ic0),geom%lat(ic0))
end do

! Averaged vertical unit
do il0=1,geom%nl0
   if (geom%mask_ver_c0(il0)) then
      geom%vunitavg(il0) = sum(geom%vunit(:,il0),geom%mask_c0(:,il0))/real(geom%nc0_mask(il0),kind_real)
   else
      geom%vunitavg(il0) = 0.0
   end if
end do

! Create mesh
call geom%mesh%create(mpl,rng,geom%nc0,geom%lon,geom%lat)
call geom%mesh%bnodes

! Compute area
if (.not.isanynotmsr(geom%area)) call geom%compute_area(mpl)

! Check whether the mask is the same for all levels
same_mask = .true.
do il0=2,geom%nl0
   same_mask = same_mask.and.(all((geom%mask_c0(:,il0).and.geom%mask_c0(:,1)) &
             & .or.(.not.geom%mask_c0(:,il0).and..not.geom%mask_c0(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   geom%nl0i = 1
else
   geom%nl0i = geom%nl0
end if
write(mpl%info,'(a7,a,i3)') '','Number of independent levels: ',geom%nl0i
call flush(mpl%info)

! Create KD-tree
call geom%kdtree%create(mpl,geom%nc0,geom%lon,geom%lat)

! Horizontal distance
allocate(geom%disth(nam%nc3))
do jc3=1,nam%nc3
   geom%disth(jc3) = real(jc3-1,kind_real)*nam%dc
end do

! Define dirac points
if (nam%check_dirac.and.(nam%ndir>0)) call geom%define_dirac(mpl,nam)

! Print summary
write(mpl%info,'(a10,a,f7.1,a,f7.1)') '','Min. / max. longitudes:',minval(geom%lon)*rad2deg,' / ',maxval(geom%lon)*rad2deg
write(mpl%info,'(a10,a,f7.1,a,f7.1)') '','Min. / max. latitudes: ',minval(geom%lat)*rad2deg,' / ',maxval(geom%lat)*rad2deg
write(mpl%info,'(a10,a)') '','Unmasked area (% of Earth area) / masked points / vertical unit:'
do il0=1,geom%nl0
   write(mpl%info,'(a13,a,i3,a,f5.1,a,f5.1,a,f12.1,a)') '','Level ',nam%levs(il0),' ~> ',geom%area(il0)/(4.0*pi)*100.0,'% / ', &
 & real(count(.not.geom%mask_c0(:,il0)),kind_real)/real(geom%nc0,kind_real)*100.0,'% / ',geom%vunitavg(il0),' '//trim(mpl%vunitchar)
end do
write(mpl%info,'(a7,a)') '','Distribution summary:'
do iproc=1,mpl%nproc
   write(mpl%info,'(a10,a,i3,a,i8,a)') '','Proc #',iproc,': ',geom%proc_to_nc0a(iproc),' grid-points'
end do
call flush(mpl%info)

end subroutine geom_init

!----------------------------------------------------------------------
! Subroutine: geom_compute_area
! Purpose: compute domain area
!----------------------------------------------------------------------
subroutine geom_compute_area(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(in) :: mpl       ! MPI data

! Local variables
integer :: il0,it
real(kind_real) :: area,frac

write(mpl%info,'(a7,a)') '','Compute area (might take a long time)'

! Create triangles list
if (.not.allocated(geom%mesh%ltri)) call geom%mesh%trlist

! Compute area
geom%area = 0.0
do it=1,geom%mesh%nt
   area = areas((/geom%mesh%x(geom%mesh%ltri(1,it)),geom%mesh%y(geom%mesh%ltri(1,it)),geom%mesh%z(geom%mesh%ltri(1,it))/), &
              & (/geom%mesh%x(geom%mesh%ltri(2,it)),geom%mesh%y(geom%mesh%ltri(2,it)),geom%mesh%z(geom%mesh%ltri(2,it))/), &
              & (/geom%mesh%x(geom%mesh%ltri(3,it)),geom%mesh%y(geom%mesh%ltri(3,it)),geom%mesh%z(geom%mesh%ltri(3,it))/))
   do il0=1,geom%nl0
      frac = real(count(geom%mask_c0(geom%mesh%order(geom%mesh%ltri(1:3,it)),il0)),kind_real)/3.0
      geom%area(il0) = geom%area(il0)+frac*area
   end do
end do

end subroutine geom_compute_area

!----------------------------------------------------------------------
! Subroutine: geom_define_dirac
! Purpose: define dirac indices
!----------------------------------------------------------------------
subroutine geom_define_dirac(geom,mpl,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(in) :: mpl       ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: idir,il0,nn_index(1),ic0dir,il0dir
real(kind_real) :: nn_dist(1)
logical :: valid

! Allocation
allocate(geom%londir(nam%ndir))
allocate(geom%latdir(nam%ndir))
allocate(geom%iprocdir(nam%ndir))
allocate(geom%ic0adir(nam%ndir))
allocate(geom%il0dir(nam%ndir))
allocate(geom%ivdir(nam%ndir))
allocate(geom%itsdir(nam%ndir))

! Initialization
geom%ndir = 0
do idir=1,nam%ndir
   ! Find level
   call msi(il0dir)
   do il0=1,geom%nl0
      if (nam%levs(il0)==nam%levdir(idir)) il0dir = il0
   end do
   if (ismsi(il0dir)) call mpl%abort('impossible to find the Dirac level')

   ! Find nearest neighbor
   call geom%kdtree%find_nearest_neighbors(nam%londir(idir),nam%latdir(idir),1,nn_index,nn_dist)
   ic0dir = nn_index(1)

   ! Check arc
   if (geom%mask_del) then
      call geom%check_arc(il0dir,nam%londir(idir),nam%latdir(idir),geom%lon(ic0dir),geom%lat(ic0dir),valid)
   else
      valid = geom%mask_c0(ic0dir,il0dir)
   end if

   if (valid) then
      ! Add valid dirac point
      geom%ndir = geom%ndir+1
      geom%londir(geom%ndir) = nam%londir(idir)
      geom%latdir(geom%ndir) = nam%latdir(idir)
      geom%iprocdir(geom%ndir) = geom%c0_to_proc(ic0dir)
      geom%ic0adir(geom%ndir) = geom%c0_to_c0a(ic0dir)
      geom%il0dir(geom%ndir) = il0dir
      geom%ivdir(geom%ndir) = nam%ivdir(idir)
      geom%itsdir(geom%ndir) = nam%itsdir(idir)
   end if
end do

end subroutine geom_define_dirac

!----------------------------------------------------------------------
! Subroutine: geom_define_distribution
! Purpose: define local distribution
!----------------------------------------------------------------------
subroutine geom_define_distribution(geom,mpl,nam,rng)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(rng_type),intent(inout) :: rng    ! Random number generator

! Local variables
integer :: ic0,il0,info,iproc,ic0a,nc0a,ny,nres,iy,delta,ix
integer :: ncid,nc0_id,c0_to_proc_id,c0_to_c0a_id,lon_id,lat_id
integer :: c0_reorder(geom%nc0),nn_index(1)
integer,allocatable :: center_to_c0(:),nx(:),ic0a_arr(:)
real(kind_real) :: nn_dist(1),dlat,dlon
real(kind_real),allocatable :: rh_c0(:),lon_center(:),lat_center(:)
logical,allocatable :: mask_hor_c0(:)
character(len=4) :: nprocchar
character(len=1024) :: filename_nc
character(len=1024) :: subr = 'geom_define_distribution'
type(kdtree_type) :: kdtree

if (mpl%nproc==1) then
   ! All points on a single processor
   geom%c0_to_proc = 1
   do ic0=1,geom%nc0
      geom%c0_to_c0a(ic0) = ic0
   end do
elseif (mpl%nproc>1) then
   if (mpl%main) then
      ! Open file
      write(nprocchar,'(i4.4)') mpl%nproc
      filename_nc = trim(nam%prefix)//'_distribution_'//nprocchar//'.nc'
      info = nf90_open(trim(nam%datadir)//'/'//trim(filename_nc),nf90_nowrite,ncid)
   end if
   call mpl%f_comm%broadcast(info,mpl%ioproc-1)

   if (info==nf90_noerr) then
      ! Read local distribution
      write(mpl%info,'(a7,a,i4,a)') '','Read local distribution for: ',mpl%nproc,' MPI tasks'
      call flush(mpl%info)

      if (mpl%main) then
         ! Get variables ID
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'c0_to_proc',c0_to_proc_id))
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'c0_to_c0a',c0_to_c0a_id))

         ! Read varaibles
         call mpl%ncerr(subr,nf90_get_var(ncid,c0_to_proc_id,geom%c0_to_proc))
         call mpl%ncerr(subr,nf90_get_var(ncid,c0_to_c0a_id,geom%c0_to_c0a))

         ! Close file
         call mpl%ncerr(subr,nf90_close(ncid))
      end if

      ! Broadcast distribution
      call mpl%f_comm%broadcast(geom%c0_to_proc,mpl%ioproc-1)
      call mpl%f_comm%broadcast(geom%c0_to_c0a,mpl%ioproc-1)

      ! Check
      if (maxval(geom%c0_to_proc)>mpl%nproc) call mpl%abort('wrong distribution')
   else
      ! Generate a distribution

      ! Allocation
      allocate(lon_center(mpl%nproc))
      allocate(lat_center(mpl%nproc))
      allocate(ic0a_arr(mpl%nproc))

      ! Define distribution centers
      if (.false.) then
         ! Using a random sampling

         ! Allocation
         allocate(mask_hor_c0(geom%nc0))
         allocate(rh_c0(geom%nc0))
         allocate(center_to_c0(mpl%nproc))

         ! Initialization
         mask_hor_c0 = any(geom%mask_c0,dim=2)
         rh_c0 = 1.0

         ! Compute sampling
         write(mpl%info,'(a7,a)',advance='no') '','Define distribution centers:'
         call rng%initialize_sampling(mpl,geom%nc0,geom%lon,geom%lat,mask_hor_c0,rh_c0,nam%ntry,nam%nrep,mpl%nproc,center_to_c0)

         ! Define centers coordinates
         lon_center = geom%lon(center_to_c0)
         lat_center = geom%lat(center_to_c0)
      else
         ! Using a regular splitting

         ! Allocation
         ny = nint(sqrt(real(mpl%nproc,kind_real)))
         if (ny**2<mpl%nproc) ny = ny+1
         allocate(nx(ny))
         nres = mpl%nproc
         do iy=1,ny
            delta = mpl%nproc/ny
            if (nres>(ny-iy+1)*delta) delta = delta+1
            nx(iy) = delta
            nres = nres-delta
         end do
         if (sum(nx)/=mpl%nproc) call mpl%abort('wrong number of tiles in define_distribution')
         dlat = (maxval(geom%lat)-minval(geom%lat))/ny
         iproc = 0
         do iy=1,ny
            dlon = (maxval(geom%lon)-minval(geom%lon))/nx(iy)
            do ix=1,nx(iy)
               iproc = iproc+1
               lat_center(iproc) = minval(geom%lat)+(real(iy,kind_real)-0.5)*dlat
               lon_center(iproc) = minval(geom%lon)+(real(ix,kind_real)-0.5)*dlon
            end do
         end do
      end if

      if (mpl%main) then
         ! Define kdtree
         call kdtree%create(mpl,mpl%nproc,lon_center,lat_center)

         ! Local processor
         do ic0=1,geom%nc0
            call kdtree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),1,nn_index,nn_dist)
            geom%c0_to_proc(ic0) = nn_index(1)
         end do

         ! Local index
         ic0a_arr = 0
         do ic0=1,geom%nc0
            iproc = geom%c0_to_proc(ic0)
            ic0a_arr(iproc) = ic0a_arr(iproc)+1
            geom%c0_to_c0a(ic0) = ic0a_arr(iproc)
         end do
      end if

      ! Broadcast distribution
      call mpl%f_comm%broadcast(geom%c0_to_proc,mpl%ioproc-1)
      call mpl%f_comm%broadcast(geom%c0_to_c0a,mpl%ioproc-1)


      if (test_no_point) then
         ! Count points on the penultimate processor
         nc0a = count(geom%c0_to_proc==mpl%nproc-1)

         ! Move all point from the last to the penultimate processor
         do ic0=1,geom%nc0
            if (geom%c0_to_proc(ic0)==mpl%nproc) then
               nc0a = nc0a+1
               geom%c0_to_proc(ic0) = mpl%nproc-1
               geom%c0_to_c0a(ic0) = nc0a
            end if
         end do
      end if

      ! Write distribution
      if (mpl%main) then
         ! Create file
         call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename_nc),or(nf90_clobber,nf90_64bit_offset),ncid))

         ! Write namelist parameters
         call nam%ncwrite(mpl,ncid)

         ! Define dimension
         call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))

         ! Define variables
         call mpl%ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'c0_to_proc',nf90_int,(/nc0_id/),c0_to_proc_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'c0_to_c0a',nf90_int,(/nc0_id/),c0_to_c0a_id))

         ! End definition mode
         call mpl%ncerr(subr,nf90_enddef(ncid))

         ! Write variables
         call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
         call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
         call mpl%ncerr(subr,nf90_put_var(ncid,c0_to_proc_id,geom%c0_to_proc))
         call mpl%ncerr(subr,nf90_put_var(ncid,c0_to_c0a_id,geom%c0_to_c0a))

         ! Close file
         call mpl%ncerr(subr,nf90_close(ncid))
      end if
   end if
end if

! Size of tiles
allocate(geom%proc_to_nc0a(mpl%nproc))
do iproc=1,mpl%nproc
   geom%proc_to_nc0a(iproc) = count(geom%c0_to_proc==iproc)
end do
geom%nc0a = geom%proc_to_nc0a(mpl%myproc)

! Conversion
allocate(geom%c0a_to_c0(geom%nc0a))
ic0a = 0
do ic0=1,geom%nc0
   if (geom%c0_to_proc(ic0)==mpl%myproc) then
      ic0a = ic0a+1
      geom%c0a_to_c0(ic0a) = ic0
   end if
end do

! Reorder Sc0 points to improve communication efficiency
do ic0=1,geom%nc0
   iproc = geom%c0_to_proc(ic0)
   ic0a = geom%c0_to_c0a(ic0)
   if (iproc==1) then
      c0_reorder(ic0) = ic0a
   else
      c0_reorder(ic0) = sum(geom%proc_to_nc0a(1:iproc-1))+ic0a
   end if
end do
geom%c0_to_lon(c0_reorder) = geom%c0_to_lon
geom%c0_to_lat(c0_reorder) = geom%c0_to_lat
geom%c0_to_tile(c0_reorder) = geom%c0_to_tile
geom%lon(c0_reorder) = geom%lon
geom%lat(c0_reorder) = geom%lat
do il0=1,geom%nl0
   geom%vunit(c0_reorder,il0) = geom%vunit(:,il0)
   geom%mask_c0(c0_reorder,il0) = geom%mask_c0(:,il0)
end do
geom%c0_to_proc(c0_reorder) = geom%c0_to_proc
geom%c0_to_c0a(c0_reorder) = geom%c0_to_c0a
do ic0a=1,geom%nc0a
   geom%c0a_to_c0(ic0a) = c0_reorder(geom%c0a_to_c0(ic0a))
end do

! Other masks
allocate(geom%mask_c0a(geom%nc0a,geom%nl0))
allocate(geom%mask_hor_c0a(geom%nc0a))
geom%mask_c0a = geom%mask_c0(geom%c0a_to_c0,:)
geom%mask_hor_c0 = any(geom%mask_c0,dim=2)
geom%mask_hor_c0a = geom%mask_hor_c0(geom%c0a_to_c0)
geom%mask_ver_c0 = any(geom%mask_c0,dim=1)
geom%nc0_mask = count(geom%mask_c0,dim=1)

end subroutine geom_define_distribution

!----------------------------------------------------------------------
! Subroutine: geom_check_arc
! Purpose: check if an arc is crossing boundaries
!----------------------------------------------------------------------
subroutine geom_check_arc(geom,il0,lon_s,lat_s,lon_e,lat_e,valid)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: il0           ! Level
real(kind_real),intent(in) :: lon_s ! First point longitude
real(kind_real),intent(in) :: lat_s ! First point latitude
real(kind_real),intent(in) :: lon_e ! Second point longitude
real(kind_real),intent(in) :: lat_e ! Second point latitude
logical,intent(out) :: valid        ! True for valid arcs

! Local variables
integer :: ibnd
real(kind_real) :: x(2),y(2),z(2),v1(3),v2(3),va(3),vp(3),t(4)

! Transform to cartesian coordinates
call trans(2,(/lat_s,lat_e/),(/lon_s,lon_e/),x,y,z)

! Compute arc orthogonal vector
v1 = (/x(1),y(1),z(1)/)
v2 = (/x(2),y(2),z(2)/)
call vector_product(v1,v2,va)

! Check if arc is crossing boundary arcs
valid = .true.
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
      valid = .false.
      exit
   end if
end do

end subroutine geom_check_arc

!----------------------------------------------------------------------
! Subroutine: geom_copy_c0a_to_mga
! Purpose: copy from subset Sc0 to model grid, halo A
!----------------------------------------------------------------------
subroutine geom_copy_c0a_to_mga(geom,mpl,fld_c0a,fld_mga)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                        ! Geometry
type(mpl_type),intent(in) :: mpl                           ! MPI data
real(kind_real),intent(in) :: fld_c0a(geom%nc0a,geom%nl0)  ! Field on subset Sc0, halo A
real(kind_real),intent(out) :: fld_mga(geom%nmga,geom%nl0) ! Field on model grid, halo A

if (geom%nc0==geom%nmg) then
   ! Model grid and subset Sc0 are identical
   fld_mga = fld_c0a
else
   ! Extend subset Sc0 to model grid
   call geom%com_mg%ext(mpl,geom%nl0,fld_c0a,fld_mga)
end if

end subroutine geom_copy_c0a_to_mga

!----------------------------------------------------------------------
! Subroutine: geom_copy_mga_to_c0a
! Purpose: copy from model grid to subset Sc0, halo A
!----------------------------------------------------------------------
subroutine geom_copy_mga_to_c0a(geom,mpl,fld_mga,fld_c0a)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                        ! Geometry
type(mpl_type),intent(in) :: mpl                           ! MPI data
real(kind_real),intent(in) :: fld_mga(geom%nmga,geom%nl0)  ! Field on model grid, halo A
real(kind_real),intent(out) :: fld_c0a(geom%nc0a,geom%nl0) ! Field on subset Sc0, halo A

! Local variables
integer :: ic0a,imga,il0
real(kind_real) :: fld_mga_zero(geom%nmga,geom%nl0)

if (geom%nc0==geom%nmg) then
   ! Model grid and subset Sc0 are identical
   fld_c0a = fld_mga
else
   ! Initialization
   fld_mga_zero = fld_mga
   do ic0a=1,geom%nc0a
      imga = geom%c0a_to_mga(ic0a)
      fld_mga_zero(imga,:) = 0.0
   end do

   ! Reduce model grid to subset Sc0
   call geom%com_mg%red(mpl,geom%nl0,fld_mga_zero,fld_c0a)

   ! Copy non-redundant points
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         imga = geom%c0a_to_mga(ic0a)
         if (abs(fld_mga(imga,il0)-fld_c0a(ic0a,il0))>0.0) then
            ! Values are different
            if (ismsr(fld_c0a(ic0a,il0))) then
               ! Subset Sc0 value is missing
               fld_c0a(ic0a,il0) = fld_mga(imga,il0)
            elseif (ismsr(fld_mga(imga,il0))) then
               ! Nothing to do
            else
               ! Both values are not missing, check for zero value
               if (.not.(abs(fld_c0a(ic0a,il0))>0.0)) then
                  ! Subset Sc0 value is zero
                  fld_c0a(ic0a,il0) = fld_mga(imga,il0)
               elseif (.not.(abs(fld_mga(imga,il0))>0.0)) then
                  ! Nothing to do
               else
                  call mpl%abort('both redundant values are different, not missing and nonzero')
               end if
            end if
         end if
      end do
   end do
end if

end subroutine geom_copy_mga_to_c0a

!----------------------------------------------------------------------
! Subroutine: geom_compute_deltas
! Purpose: compute deltas for LCT definition
!----------------------------------------------------------------------
subroutine geom_compute_deltas(geom,ic0,il0,jc0,jl0,dx,dy,dz)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: ic0           ! First horizontal index
integer,intent(in) :: il0           ! First vertical index
integer,intent(in) :: jc0           ! Second horizontal index
integer,intent(in) :: jl0           ! Second vertical index
real(kind_real),intent(out) :: dx   ! Longitude delta
real(kind_real),intent(out) :: dy   ! Latitude delta
real(kind_real),intent(out) :: dz   ! Altitude delta

! Compute deltas
dx = geom%lon(jc0)-geom%lon(ic0)
dy = geom%lat(jc0)-geom%lat(ic0)
call lonlatmod(dx,dy)
dx = dx*cos(geom%lat(ic0))
dz = real(geom%vunit(ic0,jl0)-geom%vunit(ic0,il0),kind_real)

end subroutine geom_compute_deltas

end module type_geom
