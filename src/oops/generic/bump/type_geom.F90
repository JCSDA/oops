!----------------------------------------------------------------------
! Module: type_geom
!> Purpose: geometry derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_geom

use netcdf
use tools_const, only: pi,req,deg2rad,rad2deg,reqkm,rth
use tools_func, only: lonlatmod,sphere_dist,vector_product,vector_triple_product
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,ismsi
use tools_nc, only: ncfloat
use tools_qsort, only: qsort
use tools_stripack, only: areas,trans
use type_com, only: com_type
use type_kdtree, only: kdtree_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

integer :: nredmax = 10 !< Maximum number of similar redundant points

! Geometry derived type
type geom_type
   ! Offline geometry data
   integer :: nlon                            !< Longitude size
   integer :: nlat                            !< Latitude size
   integer :: nlev                            !< Number of levels
   integer,allocatable :: c0_to_lon(:)        !< Subset Sc0 to longitude index
   integer,allocatable :: c0_to_lat(:)        !< Subset Sc0 to latgitude index

   ! Number of points and levels
   integer :: nmg                             !< Number of model grid points
   integer :: nc0                             !< Number of points in subset Sc0
   integer :: nl0                             !< Number of levels in subset Sl0
   integer :: nl0i                            !< Number of independent levels in subset Sl0

   ! Basic geometry data
   real(kind_real),allocatable :: lon(:)      !< Longitudes
   real(kind_real),allocatable :: lat(:)      !< Latitudes
   logical,allocatable :: mask(:,:)           !< Mask
   real(kind_real),allocatable :: area(:)     !< Domain area
   real(kind_real),allocatable :: vunit(:,:)  !< Vertical unit
   real(kind_real),allocatable :: vunitavg(:) !< Averaged vertical unit
   real(kind_real),allocatable :: disth(:)    !< Horizontal distance

   ! Mesh
   type(mesh_type) :: mesh                    !< Mesh

   ! KD-tree
   type(kdtree_type) :: kdtree                !< KD-tree

   ! Boundary nodes
   integer,allocatable :: nbnd(:)             !< Number of boundary nodes
   real(kind_real),allocatable :: xbnd(:,:,:) !< Boundary nodes, x-coordinate
   real(kind_real),allocatable :: ybnd(:,:,:) !< Boundary nodes, y-coordinate
   real(kind_real),allocatable :: zbnd(:,:,:) !< Boundary nodes, z-coordinate
   real(kind_real),allocatable :: vbnd(:,:,:) !< Boundary nodes, orthogonal vector

   ! Gripoints and subset Sc0
   integer,allocatable :: redundant(:)        !< Redundant points array
   integer,allocatable :: c0_to_mg(:)         !< Subset Sc0 to model grid
   integer,allocatable :: mg_to_c0(:)         !< Model grid to subset Sc0

   ! MPI distribution
   integer :: nmga                            !< Halo A size for model grid
   integer :: nc0a                            !< Halo A size for subset Sc0
   integer,allocatable :: mg_to_proc(:)       !< Model grid to local task
   integer,allocatable :: mg_to_mga(:)        !< Model grid, global to halo A
   integer,allocatable :: mga_to_mg(:)        !< Model grid, halo A to global
   integer,allocatable :: proc_to_nmga(:)     !< Halo A size for each proc
   integer,allocatable :: c0_to_proc(:)       !< Subset Sc0 to local task
   integer,allocatable :: c0_to_c0a(:)        !< Subset Sc0, global to halo A
   integer,allocatable :: c0a_to_c0(:)        !< Subset Sc0, halo A to global
   integer,allocatable :: proc_to_nc0a(:)     !< Halo A size for each proc
   integer,allocatable :: c0a_to_mga(:)       !< Subset Sc0 to model grid, halo A
   type(com_type) :: com_mg                   !< Communication between subset Sc0 and model grid
contains
   procedure :: alloc => geom_alloc
   procedure :: dealloc => geom_dealloc
   procedure :: setup_online => geom_setup_online
   procedure :: find_redundant => geom_find_redundant
   procedure :: init => geom_init
   procedure :: define_mask => geom_define_mask
   procedure :: compute_area => geom_compute_area
   procedure :: compute_mask_boundaries => geom_compute_mask_boundaries
   procedure :: define_distribution => geom_define_distribution
   procedure :: check_arc => geom_check_arc
end type geom_type

private
public :: geom_type

contains

!----------------------------------------------------------------------
! Subroutine: geom_alloc
!> Purpose: geometry allocation
!----------------------------------------------------------------------
subroutine geom_alloc(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry

! Allocation
allocate(geom%c0_to_lon(geom%nc0))
allocate(geom%c0_to_lat(geom%nc0))
allocate(geom%lon(geom%nc0))
allocate(geom%lat(geom%nc0))
allocate(geom%area(geom%nl0))
allocate(geom%vunit(geom%nc0,geom%nl0))
allocate(geom%vunitavg(geom%nl0))
allocate(geom%mask(geom%nc0,geom%nl0))

! Initialization
call msi(geom%c0_to_lon)
call msi(geom%c0_to_lat)
call msr(geom%lon)
call msr(geom%lat)
call msr(geom%area)
call msr(geom%vunit)
call msr(geom%vunitavg)
geom%mask = .false.

end subroutine geom_alloc

!----------------------------------------------------------------------
! Subroutine: geom_dealloc
!> Purpose: geometry deallocation
!----------------------------------------------------------------------
subroutine geom_dealloc(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry

! Release memory
if (allocated(geom%c0_to_lon)) deallocate(geom%c0_to_lon)
if (allocated(geom%c0_to_lat)) deallocate(geom%c0_to_lat)
if (allocated(geom%lon)) deallocate(geom%lon)
if (allocated(geom%lat)) deallocate(geom%lat)
if (allocated(geom%mask)) deallocate(geom%mask)
if (allocated(geom%area)) deallocate(geom%area)
if (allocated(geom%vunit)) deallocate(geom%vunit)
if (allocated(geom%vunitavg)) deallocate(geom%vunitavg)
if (allocated(geom%disth)) deallocate(geom%disth)
call geom%mesh%dealloc
call geom%kdtree%dealloc
if (allocated(geom%nbnd)) deallocate(geom%nbnd)
if (allocated(geom%xbnd)) deallocate(geom%xbnd)
if (allocated(geom%ybnd)) deallocate(geom%ybnd)
if (allocated(geom%zbnd)) deallocate(geom%zbnd)
if (allocated(geom%vbnd)) deallocate(geom%vbnd)
if (allocated(geom%redundant)) deallocate(geom%redundant)
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
if (allocated(geom%c0a_to_mga)) deallocate(geom%c0a_to_mga)
call geom%com_mg%dealloc

end subroutine geom_dealloc

!----------------------------------------------------------------------
! Subroutine: geom_setup_online
!> Purpose: setup online geometry
!----------------------------------------------------------------------
subroutine geom_setup_online(geom,mpl,nmga,nl0,lon,lat,area,vunit,lmask)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom        !< Geometry
type(mpl_type),intent(inout) :: mpl           !< MPI data
integer,intent(in) :: nmga                    !< Halo A size
integer,intent(in) :: nl0                     !< Number of levels in subset Sl0
real(kind_real),intent(in) :: lon(nmga)       !< Longitudes
real(kind_real),intent(in) :: lat(nmga)       !< Latitudes
real(kind_real),intent(in) :: area(nmga)      !< Area
real(kind_real),intent(in) :: vunit(nmga,nl0) !< Vertical unit
logical,intent(in) :: lmask(nmga,nl0)         !< Mask

! Local variables
integer :: ic0,ic0a,il0,offset,iproc,img,imga
real(kind_real),allocatable :: lon_mg(:),lat_mg(:),area_mg(:),vunit_mg(:,:)
logical,allocatable :: lmask_mg(:,:)

! Copy geometry variables
geom%nmga = nmga
geom%nl0 = nl0
geom%nlev = nl0

! Allocation
allocate(geom%proc_to_nmga(mpl%nproc))

! Communication
call mpl%allgather(1,(/geom%nmga/),geom%proc_to_nmga)

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
         call mpl%recv(geom%proc_to_nmga(iproc),lon_mg(offset+1:offset+geom%proc_to_nmga(iproc)),iproc,mpl%tag)
         call mpl%recv(geom%proc_to_nmga(iproc),lat_mg(offset+1:offset+geom%proc_to_nmga(iproc)),iproc,mpl%tag+1)
         call mpl%recv(geom%proc_to_nmga(iproc),area_mg(offset+1:offset+geom%proc_to_nmga(iproc)),iproc,mpl%tag+2)
         do il0=1,geom%nl0
            call mpl%recv(geom%proc_to_nmga(iproc),vunit_mg(offset+1:offset+geom%proc_to_nmga(iproc),il0),iproc,mpl%tag+2+il0)
            call mpl%recv(geom%proc_to_nmga(iproc),lmask_mg(offset+1:offset+geom%proc_to_nmga(iproc),il0),iproc, &
          & mpl%tag+2+geom%nl0+il0)
         end do
      end if

      !  Update offset
      offset = offset+geom%proc_to_nmga(iproc)
   end do
else
   ! Send data to ioproc
   call mpl%send(geom%nmga,lon,mpl%ioproc,mpl%tag)
   call mpl%send(geom%nmga,lat,mpl%ioproc,mpl%tag+1)
   call mpl%send(geom%nmga,area,mpl%ioproc,mpl%tag+2)
   do il0=1,geom%nl0
      call mpl%send(geom%nmga,vunit(:,il0),mpl%ioproc,mpl%tag+2+il0)
      call mpl%send(geom%nmga,lmask(:,il0),mpl%ioproc,mpl%tag+2+geom%nl0+il0)
   end do
end if
mpl%tag = mpl%tag+3+2*geom%nl0

if (mpl%main) then
   ! Convert to radians
   lon_mg = lon_mg*deg2rad
   lat_mg = lat_mg*deg2rad
end if

! Broadcast data
call mpl%bcast(lon_mg)
call mpl%bcast(lat_mg)
call mpl%bcast(area_mg)
call mpl%bcast(vunit_mg)
call mpl%bcast(lmask_mg)

! Find redundant points
call geom%find_redundant(mpl,lon_mg,lat_mg)

! Allocation
call geom%alloc
allocate(geom%proc_to_nc0a(mpl%nproc))
allocate(geom%c0_to_proc(geom%nc0))

! Model grid conversions and Sc0 size on halo A
img = 0
geom%proc_to_nc0a = 0
do iproc=1,mpl%nproc
   do imga=1,geom%proc_to_nmga(iproc)
      img = img+1
      geom%mg_to_proc(img) = iproc
      geom%mg_to_mga(img) = imga
      if (iproc==mpl%myproc) geom%mga_to_mg(imga) = img
      if (ismsi(geom%redundant(img))) geom%proc_to_nc0a(iproc) = geom%proc_to_nc0a(iproc)+1
   end do
end do
geom%nc0a = geom%proc_to_nc0a(mpl%myproc)

! Subset Sc0 conversions
allocate(geom%c0a_to_c0(geom%nc0a))
allocate(geom%c0_to_c0a(geom%nc0))
ic0 = 0
do iproc=1,mpl%nproc
   do ic0a=1,geom%proc_to_nc0a(iproc)
      ic0 = ic0+1
      geom%c0_to_proc(ic0) = iproc
      if (iproc==mpl%myproc) geom%c0a_to_c0(ic0a) = ic0
   end do
end do
call mpl%glb_to_loc(geom%nc0a,geom%c0a_to_c0,geom%nc0,geom%c0_to_c0a)

! Inter-halo conversions
allocate(geom%c0a_to_mga(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   img = geom%c0_to_mg(ic0)
   imga = geom%mg_to_mga(img)
   geom%c0a_to_mga(ic0a) = imga
end do

! Setup communications
call geom%com_mg%setup(mpl,'com_mg',geom%nmg,geom%nc0a,geom%nmga,geom%mga_to_mg,geom%c0a_to_mga, &
 & geom%c0_to_proc(geom%mg_to_c0),geom%c0_to_c0a)

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
   geom%mask(:,il0) = lmask_mg(geom%c0_to_mg,il0)
end do

end subroutine geom_setup_online

!----------------------------------------------------------------------
! Subroutine: geom_find_redundant
!> Purpose: find redundant model grid points
!----------------------------------------------------------------------
subroutine geom_find_redundant(geom,mpl,lon,lat)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom               !< Geometry
type(mpl_type),intent(in) :: mpl                     !< MPI data
real(kind_real),intent(in),optional :: lon(geom%nmg) !< Longitudes
real(kind_real),intent(in),optional :: lat(geom%nmg) !< Latitudes

! Local variables
integer :: img,ic0,ired,nn_index(nredmax)
real(kind_real) :: nn_dist(nredmax)
type(kdtree_type) :: kdtree

! Allocation
allocate(geom%redundant(geom%nmg))
call msi(geom%redundant)

! Look for redundant points
if (present(lon).and.present(lat)) then
   write(mpl%unit,'(a7,a)') '','Look for redundant points in the model grid'

   ! Create KD-tree
   call kdtree%create(mpl,geom%nmg,lon,lat)

   ! Find redundant points
   do img=1,geom%nmg
      ! Find nearest neighbors
      call kdtree%find_nearest_neighbors(lon(img),lat(img),nredmax,nn_index,nn_dist)

      ! Count redundant points
      do ired=1,nredmax
         if ((nn_dist(ired)>rth).or.(nn_index(ired)>=img)) nn_index(ired) = geom%nmg+1
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
end if
geom%nc0 = count(ismsi(geom%redundant))
write(mpl%unit,'(a7,a,i8)') '','Model grid size:         ',geom%nmg
write(mpl%unit,'(a7,a,i8)') '','Subset Sc0 size:         ',geom%nc0
write(mpl%unit,'(a7,a,i6,a,f6.2,a)') '','Number of redundant points:',(geom%nmg-geom%nc0), &
 & ' (',real(geom%nmg-geom%nc0,kind_real)/real(geom%nmg,kind_real)*100.0,'%)'
call flush(mpl%unit)

! Conversion
allocate(geom%c0_to_mg(geom%nc0))
allocate(geom%mg_to_c0(geom%nmg))
ic0 = 0
do img=1,geom%nmg
   if (ismsi(geom%redundant(img))) then
      ic0 = ic0+1
      geom%c0_to_mg(ic0) = img
      geom%mg_to_c0(img) = ic0
   end if
end do
do img=1,geom%nmg
   if (isnotmsi(geom%redundant(img))) geom%mg_to_c0(img) = geom%mg_to_c0(geom%redundant(img))
end do

end subroutine geom_find_redundant

!----------------------------------------------------------------------
! Subroutine: geom_init
!> Purpose: initialize geometry
!----------------------------------------------------------------------
subroutine geom_init(geom,mpl,rng,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(mpl_type),intent(in) :: mpl       !< MPI data
type(rng_type),intent(inout) :: rng    !< Random number generator
type(nam_type),intent(in) :: nam       !< Namelist

! Local variables
integer :: ic0,il0,jc3,iproc
logical :: same_mask

! Set longitude and latitude bounds
do ic0=1,geom%nc0
   call lonlatmod(geom%lon(ic0),geom%lat(ic0))
end do

! Define mask
call geom%define_mask(mpl,nam)

! Averaged vertical unit
do il0=1,geom%nl0
   if (any(geom%mask(:,il0))) then
      geom%vunitavg(il0) = sum(geom%vunit(:,il0),geom%mask(:,il0))/real(count(geom%mask(:,il0)),kind_real)
   else
      geom%vunitavg(il0) = 0.0
   end if
end do

! Create mesh
call geom%mesh%create(mpl,rng,geom%nc0,geom%lon,geom%lat)
call geom%mesh%bnodes

! Compute area
if ((.not.any(geom%area>0.0))) call geom%compute_area

! Compute mask boundaries
if ((nam%new_param.or.nam%new_lct).and.nam%mask_check) call geom%compute_mask_boundaries(mpl)

! Check whether the mask is the same for all levels
same_mask = .true.
do il0=2,geom%nl0
   same_mask = same_mask.and.(all((geom%mask(:,il0).and.geom%mask(:,1)) &
             & .or.(.not.geom%mask(:,il0).and..not.geom%mask(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   geom%nl0i = 1
else
   geom%nl0i = geom%nl0
end if
write(mpl%unit,'(a7,a,i3)') '','Number of independent levels: ',geom%nl0i
call flush(mpl%unit)

! Create KD-tree
call geom%kdtree%create(mpl,geom%nc0,geom%lon,geom%lat)

! Horizontal distance
allocate(geom%disth(nam%nc3))
do jc3=1,nam%nc3
   geom%disth(jc3) = real(jc3-1,kind_real)*nam%dc
end do

! Print summary
write(mpl%unit,'(a10,a,f7.1,a,f7.1)') '','Min. / max. longitudes:',minval(geom%lon)*rad2deg,' / ',maxval(geom%lon)*rad2deg
write(mpl%unit,'(a10,a,f7.1,a,f7.1)') '','Min. / max. latitudes: ',minval(geom%lat)*rad2deg,' / ',maxval(geom%lat)*rad2deg
write(mpl%unit,'(a10,a)') '','Unmasked area (% of Earth area) / masked points / vertical unit:'
do il0=1,geom%nl0
   write(mpl%unit,'(a13,a,i3,a,f5.1,a,f5.1,a,f12.1,a)') '','Level ',nam%levs(il0),' ~> ',geom%area(il0)/(4.0*pi)*100.0,'% / ', &
 & real(count(.not.geom%mask(:,il0)),kind_real)/real(geom%nc0,kind_real)*100.0,'% / ',geom%vunitavg(il0),' '//trim(mpl%vunitchar)
end do
write(mpl%unit,'(a7,a)') '','Distribution summary:'
do iproc=1,mpl%nproc
   write(mpl%unit,'(a10,a,i3,a,i8,a)') '','Proc #',iproc,': ',geom%proc_to_nc0a(iproc),' grid-points'
end do
call flush(mpl%unit)

end subroutine geom_init

!----------------------------------------------------------------------
! Subroutine: geom_define_mask
!> Purpose: define mask
!----------------------------------------------------------------------
subroutine geom_define_mask(geom,mpl,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(mpl_type),intent(in) :: mpl       !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist

! Local variables
integer :: latmin,latmax,il0,ic0,ildw
integer :: ncid,nlon_id,nlon_test,nlat_id,nlat_test,mask_id
real(kind_real) :: dist
real(kind_real),allocatable :: hydmask(:,:)
logical :: mask_test
character(len=3) :: il0char
character(len=1024) :: subr = 'geom_define_mask'

! Mask restriction
if (nam%mask_type(1:3)=='lat') then
   ! Latitude mask
   read(nam%mask_type(4:6),'(i3)') latmin
   read(nam%mask_type(7:9),'(i3)') latmax
   if (latmin>=latmax) call mpl%abort('latmin should be lower than latmax')
   do il0=1,geom%nl0
      geom%mask(:,il0) = geom%mask(:,il0).and.(geom%lat>=real(latmin,kind_real)*deg2rad) &
                       & .and.(geom%lat<=real(latmax,kind_real)*deg2rad)
   end do
elseif (trim(nam%mask_type)=='hyd') then
   ! Read from hydrometeors mask file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_hyd.nc',nf90_nowrite,ncid))
   if (trim(nam%model)=='aro') then
      call mpl%ncerr(subr,nf90_inq_dimid(ncid,'X',nlon_id))
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=nlon_test))
      call mpl%ncerr(subr,nf90_inq_dimid(ncid,'Y',nlat_id))
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=nlat_test))
      if ((nlon_test/=geom%nlon).or.(nlat_test/=geom%nlat)) call mpl%abort('wrong dimensions in the mask')
      allocate(hydmask(geom%nlon,geom%nlat))
      do il0=1,geom%nl0
         write(il0char,'(i3.3)') nam%levs(il0)
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'S'//il0char//'MASK',mask_id))
         call mpl%ncerr(subr,nf90_get_var(ncid,mask_id,hydmask,(/1,1/),(/geom%nlon,geom%nlat/)))
         geom%mask(:,il0) = geom%mask(:,il0).and.pack(real(hydmask,kind(1.0))>nam%mask_th,mask=.true.)
      end do
      deallocate(hydmask)
      call mpl%ncerr(subr,nf90_close(ncid))
   end if
elseif (trim(nam%mask_type)=='ldwv') then
   ! Compute distance to the vertical diagnostic points
   do ic0=1,geom%nc0
      if (any(geom%mask(ic0,:))) then
         mask_test = .false.
         do ildw=1,nam%nldwv
            call sphere_dist(nam%lon_ldwv(ildw),nam%lat_ldwv(ildw),geom%lon(ic0),geom%lat(ic0),dist)
            mask_test = mask_test.or.(dist<1.1*nam%local_rad)
         end do
         do il0=1,geom%nl0
            if (geom%mask(ic0,il0)) geom%mask(ic0,:) = mask_test
         end do
      end if
   end do
end if

end subroutine geom_define_mask

!----------------------------------------------------------------------
! Subroutine: geom_compute_area
!> Purpose: compute domain area
!----------------------------------------------------------------------
subroutine geom_compute_area(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: il0,it
real(kind_real) :: area,frac

! Create triangles list
if (.not.allocated(geom%mesh%ltri)) call geom%mesh%trlist

! Compute area
geom%area = 0.0
do it=1,geom%mesh%nt
   area = areas((/geom%mesh%x(geom%mesh%ltri(1,it)),geom%mesh%y(geom%mesh%ltri(1,it)),geom%mesh%z(geom%mesh%ltri(1,it))/), &
              & (/geom%mesh%x(geom%mesh%ltri(2,it)),geom%mesh%y(geom%mesh%ltri(2,it)),geom%mesh%z(geom%mesh%ltri(2,it))/), &
              & (/geom%mesh%x(geom%mesh%ltri(3,it)),geom%mesh%y(geom%mesh%ltri(3,it)),geom%mesh%z(geom%mesh%ltri(3,it))/))
   do il0=1,geom%nl0
      frac = real(count(geom%mask(geom%mesh%order(geom%mesh%ltri(1:3,it)),il0)),kind_real)/3.0
      geom%area(il0) = geom%area(il0)+frac*area
   end do
end do

end subroutine geom_compute_area

!----------------------------------------------------------------------
! Subroutine: geom_compute_mask_boundaries
!> Purpose: compute domain area
!----------------------------------------------------------------------
subroutine geom_compute_mask_boundaries(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(mpl_type),intent(in) :: mpl       !< MPI data

! Local variables
integer :: i,j,k,iend,ic0,jc0,kc0,ibnd,il0
integer,allocatable :: ic0_bnd(:,:,:)
real(kind_real) :: latbnd(2),lonbnd(2),v1(3),v2(3)
logical :: init

! Allocation
allocate(geom%nbnd(geom%nl0))
allocate(ic0_bnd(2,geom%mesh%n,geom%nl0))

! Find border points
do il0=1,geom%nl0
   geom%nbnd(il0) = 0
   do i=1,geom%mesh%n
      ! Check mask points only
      ic0 = geom%mesh%order(i)
      if (.not.geom%mask(ic0,il0)) then
         iend = geom%mesh%lend(i)
         init = .true.
         do while ((iend/=geom%mesh%lend(i)).or.init)
            j = abs(geom%mesh%list(iend))
            k = abs(geom%mesh%list(geom%mesh%lptr(iend)))
            jc0 = geom%mesh%order(j)
            kc0 = geom%mesh%order(k)
            if (.not.geom%mask(jc0,il0).and.geom%mask(kc0,il0)) then
               ! Create a new boundary arc
               geom%nbnd(il0) = geom%nbnd(il0)+1
               if (geom%nbnd(il0)>geom%mesh%n) call mpl%abort('too many boundary arcs')
               ic0_bnd(1,geom%nbnd(il0),il0) = ic0
               ic0_bnd(2,geom%nbnd(il0),il0) = jc0
            end if
            iend = geom%mesh%lptr(iend)
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
      latbnd = geom%lat(ic0_bnd(:,ibnd,il0))
      lonbnd = geom%lon(ic0_bnd(:,ibnd,il0))
      call trans(2,latbnd,lonbnd,geom%xbnd(:,ibnd,il0),geom%ybnd(:,ibnd,il0),geom%zbnd(:,ibnd,il0))
   end do
   do ibnd=1,geom%nbnd(il0)
      v1 = (/geom%xbnd(1,ibnd,il0),geom%ybnd(1,ibnd,il0),geom%zbnd(1,ibnd,il0)/)
      v2 = (/geom%xbnd(2,ibnd,il0),geom%ybnd(2,ibnd,il0),geom%zbnd(2,ibnd,il0)/)
      call vector_product(v1,v2,geom%vbnd(:,ibnd,il0))
   end do
end do

end subroutine geom_compute_mask_boundaries

!----------------------------------------------------------------------
! Subroutine: geom_define_distribution
!> Purpose: define local distribution
!----------------------------------------------------------------------
subroutine geom_define_distribution(geom,mpl,nam,rng)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(mpl_type),intent(in) :: mpl       !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(rng_type),intent(inout) :: rng    !< Random number generator

! Local variables
integer :: ic0,il0,i,j,iend,info,iproc,ic0a,nc0amax,lunit
integer :: ncid,nc0_id,c0_to_proc_id,c0_to_c0a_id,lon_id,lat_id
integer :: c0_reorder(geom%nc0)
integer,allocatable :: nr_to_proc(:),ic0a_arr(:)
logical :: init,ismetis
character(len=4) :: nprocchar
character(len=1024) :: filename_nc,filename_metis
character(len=1024) :: subr = 'geom_define_distribution'
type(mesh_type) :: mesh

! Allocation
allocate(geom%c0_to_proc(geom%nc0))
allocate(geom%c0_to_c0a(geom%nc0))

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
   call mpl%bcast(info)

   if (info==nf90_noerr) then
      ! Read local distribution
      write(mpl%unit,'(a7,a,i4,a)') '','Read local distribution for: ',mpl%nproc,' MPI tasks'
      call flush(mpl%unit)

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
      call mpl%bcast(geom%c0_to_proc)
      call mpl%bcast(geom%c0_to_c0a)

      ! Check
      if (maxval(geom%c0_to_proc)>mpl%nproc) call mpl%abort('wrong distribution')
   else
      ! Generate a distribution
      if (nam%use_metis) then
         write(mpl%unit,'(a7,a,i4,a)') '','Try to use METIS for ',mpl%nproc,' MPI tasks'
         call flush(mpl%unit)

         ! Compute graph
         call mesh%create(mpl,rng,geom%nc0,geom%lon,geom%lat)
         call mesh%bnodes

         if (mpl%main) then
            ! Open file
            filename_metis = trim(nam%prefix)//'_metis'
            call mpl%newunit(lunit)
            open(unit=lunit,file=trim(nam%datadir)//'/'//trim(filename_metis),status='replace')

            ! Write header
            write(lunit,*) mesh%n,mesh%na-mesh%nb/2

            ! Write connectivity
            do i=1,mesh%n
               iend = mesh%lend(i)
               init = .true.
               do while ((iend/=mesh%lend(i)).or.init)
                  j = mesh%list(iend)
                  if (j>0) write(lunit,'(i7)',advance='no') j
                  iend = mesh%lptr(iend)
                  init = .false.
               end do
               write(lunit,*) ''
            end do

            ! Close file
            close(unit=lunit)

            ! Call METIS
            write(nprocchar,'(i4)') mpl%nproc
            call system('gpmetis '//trim(nam%datadir)//'/'//trim(filename_metis)//' '//adjustl(nprocchar)//' > '// &
          & trim(nam%datadir)//'/'//trim(filename_metis)//'.out')

            ! Check for METIS output
            inquire(file=trim(nam%datadir)//'/'//trim(filename_metis)//'.part.'//adjustl(nprocchar),exist=ismetis)
            if (.not.ismetis) call mpl%warning('METIS not available to generate the local distribution')
         end if
         call mpl%bcast(ismetis)
      else
         ! No METIS
         ismetis = .false.
      end if

      if (ismetis) then
         write(mpl%unit,'(a7,a)') '','Use METIS to generate the local distribution'
         call flush(mpl%unit)

         if (mpl%main) then
            ! Allocation
            allocate(nr_to_proc(mesh%n))
            allocate(ic0a_arr(mpl%nproc))

            ! Read METIS file
            call mpl%newunit(lunit)
            open(unit=lunit,file=trim(nam%datadir)//'/'//trim(filename_metis)//'.part.'//adjustl(nprocchar),status='old')
            do i=1,mesh%n
               read(lunit,*) nr_to_proc(i)
            end do
            close(unit=lunit)

            ! Reorder and offset
            do ic0=1,geom%nc0
               i = mesh%order_inv(ic0)
               geom%c0_to_proc(ic0) = nr_to_proc(i)+1
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
         call mpl%bcast(geom%c0_to_proc)
         call mpl%bcast(geom%c0_to_c0a)
      else
         write(mpl%unit,'(a7,a)') '','Define a basic local distribution'
         call flush(mpl%unit)

         ! Basic distribution
         nc0amax = geom%nc0/mpl%nproc
         if (nc0amax*mpl%nproc<geom%nc0) nc0amax = nc0amax+1
         iproc = 1
         ic0a = 1
         do ic0=1,geom%nc0
            geom%c0_to_proc(ic0) = iproc
            geom%c0_to_c0a(ic0) = ic0a
            ic0a = ic0a+1
            if (ic0a>nc0amax) then
               ! Change proc
               iproc = iproc+1
               ic0a = 1
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
geom%lon(c0_reorder) = geom%lon
geom%lat(c0_reorder) = geom%lat
do il0=1,geom%nl0
   geom%vunit(c0_reorder,il0) = geom%vunit(:,il0)
   geom%mask(c0_reorder,il0) = geom%mask(:,il0)
end do
geom%c0_to_proc(c0_reorder) = geom%c0_to_proc
geom%c0_to_c0a(c0_reorder) = geom%c0_to_c0a
do ic0a=1,geom%nc0a
   geom%c0a_to_c0(ic0a) = c0_reorder(geom%c0a_to_c0(ic0a))
end do

end subroutine geom_define_distribution

!----------------------------------------------------------------------
! Subroutine: geom_check_arc
!> Purpose: check if an arc is crossing boundaries
!----------------------------------------------------------------------
subroutine geom_check_arc(geom,il0,lon_s,lat_s,lon_e,lat_e,valid)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom !< Geometry
integer,intent(in) :: il0           !< Level
real(kind_real),intent(in) :: lon_s !< First point longitude
real(kind_real),intent(in) :: lat_s !< First point latitude
real(kind_real),intent(in) :: lon_e !< Second point longitude
real(kind_real),intent(in) :: lat_e !< Second point latitude
logical,intent(out) :: valid        !< True for valid arcs

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

end module type_geom
