!----------------------------------------------------------------------
! Module: type_bump
!> Purpose: BUMP derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_bump

use netcdf
use model_interface, only: model_coord
use tools_const, only: req,deg2rad
use tools_func, only: sphere_dist
use tools_kinds,only: kind_real
use type_bpar, only: bpar_type
use type_cmat, only: cmat_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_lct, only: lct_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_nicas, only: nicas_type
use type_obsop, only: obsop_type
use type_rng, only: rng_type
use type_timer, only: timer_type

implicit none

logical,parameter :: write_online = .false. !< Write online data for tests

! BUMP derived type
type bump_type
  type(bpar_type) :: bpar
  type(cmat_type) :: cmat
  type(ens_type) :: ens1
  type(ens_type) :: ens2
  type(geom_type) :: geom
  type(io_type) :: io
  type(lct_type) :: lct
  type(mpl_type) :: mpl
  type(nam_type) :: nam
  type(nicas_type) :: nicas
  type(obsop_type) :: obsop
  type(rng_type) :: rng
  real(kind_real),allocatable :: rh(:,:,:,:)
  real(kind_real),allocatable :: rv(:,:,:,:)
  integer :: nx
  integer :: ny
contains
   procedure :: setup_offline => bump_setup_offline
   procedure :: bump_setup_online
   procedure :: bump_setup_online_oops
   procedure :: bump_setup_online_nemovar
   generic :: setup_online => bump_setup_online,bump_setup_online_oops,bump_setup_online_nemovar
   procedure :: setup_generic => bump_setup_generic
   procedure :: run_drivers => bump_run_drivers
   procedure :: apply_nicas => bump_apply_nicas
   procedure :: apply_obsop => bump_apply_obsop
   procedure :: apply_obsop_ad => bump_apply_obsop_ad
   procedure :: back_to => bump_back_to
   procedure :: dealloc => bump_dealloc
end type bump_type

private
public :: bump_type

contains

!----------------------------------------------------------------------
! Subroutine: bump_setup_offline
!> Purpose: offline setup
!----------------------------------------------------------------------
subroutine bump_setup_offline(bump,mpi_comm,namelname)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump   !< BUMP
integer,intent(in) :: mpi_comm           !< MPI communicator
character(len=*),intent(in) :: namelname !< Namelist name

! Local variables
type(timer_type) :: timer

! Initialize MPL
call bump%mpl%init(mpi_comm)

! Initialize timer
if (bump%mpl%main) call timer%start

! Initialize, read and broadcast namelist
call bump%nam%init
call bump%nam%read(bump%mpl,namelname)
call bump%nam%bcast(bump%mpl)

! Initialize listing
call bump%mpl%init_listing(bump%nam%prefix,bump%nam%model,bump%nam%colorlog,bump%nam%logpres)

! Generic setup, first step
call bump%setup_generic

! Initialize geometry
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize geometry'
call flush(bump%mpl%unit)
call model_coord(bump%mpl,bump%rng,bump%nam,bump%geom)
call bump%geom%init(bump%mpl,bump%rng,bump%nam)

if (bump%nam%grid_output) then
   ! Initialize fields regridding
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize fields regridding'
   call flush(bump%mpl%unit)
   call bump%io%grid_init(bump%mpl,bump%rng,bump%nam,bump%geom)
end if

! Initialize block parameters
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize block parameters'
call bump%bpar%alloc(bump%nam,bump%geom)

if (bump%nam%new_hdiag.or.bump%nam%new_lct.or.bump%nam%load_ensemble) then
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Load ensemble 1'
   call flush(bump%mpl%unit)
   call bump%ens1%load(bump%mpl,bump%nam,bump%geom,'ens1')
end if

if (bump%nam%new_hdiag.and.((trim(bump%nam%method)=='hyb-rnd').or.(trim(bump%nam%method)=='dual-ens'))) then
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Load ensemble 2'
   call flush(bump%mpl%unit)
   call bump%ens2%load(bump%mpl,bump%nam,bump%geom,'ens2')
end if

if (bump%nam%new_obsop) then
   ! Generate observations locations
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Generate observations locations'
   call flush(bump%mpl%unit)
   call bump%obsop%generate(bump%mpl,bump%rng,bump%nam,bump%geom)
end if

! Run drivers
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Run drivers'
call flush(bump%mpl%unit)
call bump%run_drivers

! Execution stats
if (bump%mpl%main) then
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Execution stats'
   call timer%display(bump%mpl)
end if
call flush(bump%mpl%unit)

! Close listings
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Close listings'
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
call flush(bump%mpl%unit)
close(unit=bump%mpl%unit)

end subroutine bump_setup_offline

!----------------------------------------------------------------------
! Subroutine: bump_setup_online
!> Purpose: online setup
!----------------------------------------------------------------------
subroutine bump_setup_online(bump,mpi_comm,nmga,nl0,nv,nts,lon,lat,area,vunit,lmask,ens1_ne,ens1,ens2_ne,ens2,rh,rv, &
                                & nobs,lonobs,latobs)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                 !< BUMP
integer,intent(in) :: mpi_comm                         !< MPI communicator
integer,intent(in) :: nmga                             !< Halo A size
integer,intent(in) :: nl0                              !< Number of levels in subset Sl0
integer,intent(in) :: nv                               !< Number of variables
integer,intent(in) :: nts                              !< Number of time slots
real(kind_real),intent(in) :: lon(nmga)                !< Longitude (in degrees: -180 to 180)
real(kind_real),intent(in) :: lat(nmga)                !< Latitude (in degrees: -90 to 90)
real(kind_real),intent(in) :: area(nmga)               !< Area (in m^2)
real(kind_real),intent(in) :: vunit(nmga,nl0)          !< Vertical unit
logical,intent(in) :: lmask(nmga,nl0)                  !< Mask
integer,intent(in),optional :: ens1_ne                 !< Ensemble 1 size
real(kind_real),intent(in),optional :: ens1(:,:,:,:,:) !< Ensemble 1
integer,intent(in),optional :: ens2_ne                 !< Ensemble 2 size
real(kind_real),intent(in),optional :: ens2(:,:,:,:,:) !< Ensemble 2
real(kind_real),intent(in),optional :: rh(:,:,:,:)     !< Horizontal support radius for covariance (in m)
real(kind_real),intent(in),optional :: rv(:,:,:,:)     !< Vertical support radius for covariance
integer,intent(in),optional :: nobs                    !< Number of observations
real(kind_real),intent(in),optional :: lonobs(:)       !< Observations longitude (in degrees: -180 to 180)
real(kind_real),intent(in),optional :: latobs(:)       !< Observations latitude (in degrees: -90 to 90)

! Local variables
logical :: test(3)

! Initialize MPL
call bump%mpl%init(mpi_comm)

! Set internal namelist parameters
if (present(ens1_ne).and.present(ens2_ne)) then
   call bump%nam%setup_internal(nl0,nv,nts,ens1_ne,ens2_ne)
elseif (present(ens1_ne)) then
   call bump%nam%setup_internal(nl0,nv,nts,ens1_ne)
else
   call bump%nam%setup_internal(nl0,nv,nts)
end if

! Initialize listing
call bump%mpl%init_listing(bump%nam%prefix,bump%nam%model,bump%nam%colorlog,bump%nam%logpres)

! Generic setup
call bump%setup_generic

! Check arguments consistency
test(1:2) = (/present(ens1_ne),present(ens1)/)
if (.not.(all(test(1:2)).or.all(.not.test(1:2)))) call bump%mpl%abort('ens1_ne and ens1 should be present together')
test(1:2) = (/present(ens2_ne),present(ens2)/)
if (.not.(all(test(1:2)).or.all(.not.test(1:2)))) call bump%mpl%abort('ens2_ne and ens2 should be present together')
test(1:2) = (/present(rh),present(rv)/)
if (.not.(all(test(1:2)).or.all(.not.test(1:2)))) call bump%mpl%abort('rh and rv should be present together')
test(1:3) = (/present(nobs),present(lonobs),present(latobs)/)
if (.not.(all(test(1:3)).or.all(.not.test(1:3)))) call bump%mpl%abort('nobs, lonobs and latobs should be present together')

! Check sizes consistency
if (present(ens1_ne).and.present(ens1)) then
   if (size(ens1)/=nmga*nl0*nv*nts*ens1_ne) call bump%mpl%abort('wrong size for ens1')
end if
if (present(ens2_ne).and.present(ens2)) then
   if (size(ens2)/=nmga*nl0*nv*nts*ens2_ne) call bump%mpl%abort('wrong size for ens2')
end if
if (present(rh).and.present(rv)) then
   if (size(rh)/=nmga*nl0*nv*nts) call bump%mpl%abort('wrong size for rh')
   if (size(rv)/=nmga*nl0*nv*nts) call bump%mpl%abort('wrong size for rv')
end if
if (present(nobs).and.present(lonobs).and.present(latobs)) then
   if (size(lonobs)/=nobs) call bump%mpl%abort('wrong size for lonobs')
   if (size(latobs)/=nobs) call bump%mpl%abort('wrong size for latobs')
end if

! Initialize geometry
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize geometry'
call flush(bump%mpl%unit)
call bump%geom%setup_online(bump%mpl,nmga,nl0,lon,lat,area,vunit,lmask)
call bump%geom%init(bump%mpl,bump%rng,bump%nam)

if (bump%nam%grid_output) then
   ! Initialize fields regridding
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize fields regridding'
   call flush(bump%mpl%unit)
   call bump%io%grid_init(bump%mpl,bump%rng,bump%nam,bump%geom)
end if

! Initialize block parameters
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize block parameters'
call bump%bpar%alloc(bump%nam,bump%geom)

if (present(ens1_ne).and.present(ens1)) then
   ! Initialize ensemble 1
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize ensemble 1'
   call flush(bump%mpl%unit)
   call bump%ens1%from(bump%nam,bump%geom,ens1_ne,ens1)
end if

if (present(ens2_ne).and.present(ens2)) then
   ! Initialize ensemble 2
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize ensemble 2'
   call flush(bump%mpl%unit)
   call bump%ens2%from(bump%nam,bump%geom,ens2_ne,ens2)
end if

if (present(rh).and.present(rv)) then
   ! Initialize C matrix from support radii
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize C matrix from support radii'
   call flush(bump%mpl%unit)
   call bump%cmat%from(bump%mpl,bump%nam,bump%geom,bump%bpar,rh,rv)
end if

if (present(nobs).and.present(lonobs).and.present(latobs)) then
   ! Initialize observations locations
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize observations locations'
   call flush(bump%mpl%unit)
   call bump%obsop%from(nobs,lonobs,latobs)
end if

! Run drivers
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Run drivers'
call flush(bump%mpl%unit)
call bump%run_drivers

! Close listings
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Close listings'
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
call flush(bump%mpl%unit)
close(unit=bump%mpl%unit)

end subroutine bump_setup_online

!----------------------------------------------------------------------
! Subroutine: bump_setup_online_oops
!> Purpose: online setup for OOPS
!----------------------------------------------------------------------
subroutine bump_setup_online_oops(bump,mpi_comm,nmga,nl0,nv,nts,lon,lat,area,vunit_vec,imask_vec,ens1_ne,ens1_vec)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump             !< BUMP
integer,intent(in) :: mpi_comm                     !< MPI communicator
integer,intent(in) :: nmga                         !< Halo A size
integer,intent(in) :: nl0                          !< Number of levels in subset Sl0
integer,intent(in) :: nv                           !< Number of variables
integer,intent(in) :: nts                          !< Number of time slots
real(kind_real),intent(in) :: lon(nmga)            !< Longitude (in degrees: -180 to 180)
real(kind_real),intent(in) :: lat(nmga)            !< Latitude (in degrees: -90 to 90)
real(kind_real),intent(in) :: area(nmga)           !< Area (in m^2)
real(kind_real),intent(in) :: vunit_vec(nmga*nl0)  !< Vertical unit
integer,intent(in) :: imask_vec(nmga*nl0)          !< Mask
integer,intent(in),optional :: ens1_ne             !< Ensemble 1 size
real(kind_real),intent(in),optional :: ens1_vec(:) !< Ensemble 1

! Local variables
integer :: il0,imga,offset
real(kind_real) :: vunit(nmga,nl0)
logical :: test(3),lmask(nmga,nl0)

! Initialize MPL
call bump%mpl%init(mpi_comm)

! Set internal namelist parameters
if (present(ens1_ne)) then
   call bump%nam%setup_internal(nl0,nv,nts,ens1_ne)
else
   call bump%nam%setup_internal(nl0,nv,nts)
end if

! Initialize listing
call bump%mpl%init_listing(bump%nam%prefix,bump%nam%model,bump%nam%colorlog,bump%nam%logpres)

! Generic setup
call bump%setup_generic

! Check arguments consistency
test(1:2) = (/present(ens1_ne),present(ens1_vec)/)
if (.not.(all(test(1:2)).or.all(.not.test(1:2)))) call bump%mpl%abort('ens1_ne and ens1 should be present together')

! Check sizes consistency
if (present(ens1_ne).and.present(ens1_vec)) then
   if (size(ens1_vec)/=nmga*nl0*nv*nts*ens1_ne) call bump%mpl%abort('wrong size for ens1')
end if

! Convert vunit and mask
do il0=1,nl0
   offset = (il0-1)*nmga
   do imga=1,nmga
      vunit(imga,il0) = vunit_vec(offset+imga)
      if (imask_vec(offset+imga)==0) then
         lmask(imga,il0) = .false.
      elseif (imask_vec(offset+imga)==1) then
         lmask(imga,il0) = .true.
      end if
   end do
end do

! Initialize geometry
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize geometry'
call flush(bump%mpl%unit)
call bump%geom%setup_online(bump%mpl,nmga,nl0,lon,lat,area,vunit,lmask)
call bump%geom%init(bump%mpl,bump%rng,bump%nam)

if (bump%nam%grid_output) then
   ! Initialize fields regridding
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize fields regridding'
   call flush(bump%mpl%unit)
   call bump%io%grid_init(bump%mpl,bump%rng,bump%nam,bump%geom)
end if

! Initialize block parameters
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize block parameters'
call bump%bpar%alloc(bump%nam,bump%geom)

if (present(ens1_ne).and.present(ens1_vec)) then
   ! Initialize ensemble 1
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize ensemble 1'
   call flush(bump%mpl%unit)
   call bump%ens1%from(bump%nam,bump%geom,ens1_ne,ens1_vec)
end if

! Run drivers
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Run drivers'
call flush(bump%mpl%unit)
call bump%run_drivers

! Close listings
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Close listings'
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
call flush(bump%mpl%unit)
close(unit=bump%mpl%unit)

end subroutine bump_setup_online_oops

!----------------------------------------------------------------------
! Subroutine: bump_setup_online_nemovar
!> Purpose: online setup for NEMOVAR
!----------------------------------------------------------------------
subroutine bump_setup_online_nemovar(bump,mpi_comm,lunit,nx,ny,nl0,lon,lat,area,vunit,lmask,nens,ncyc,ens1_2d,ens1_3d)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                    !< BUMP
integer,intent(in) :: mpi_comm                            !< MPI communicator
integer,intent(in) :: lunit                               !< Main listing unit
integer,intent(in) :: nx                                  !< X-axis size
integer,intent(in) :: ny                                  !< Y-axis size
integer,intent(in) :: nl0                                 !< Number of levels
real(kind_real),intent(in) :: lon(nx,ny)                  !< Longitude (in degrees: -180 to 180)
real(kind_real),intent(in) :: lat(nx,ny)                  !< Latitude (in degrees: -90 to 90)
real(kind_real),intent(in) :: area(nx,ny)                 !< Area (in m^2)
real(kind_real),intent(in) :: vunit(nx,ny,nl0)            !< Vertical unit
logical,intent(in) :: lmask(nx,ny,nl0)                    !< Mask
integer,intent(in) :: nens                                !< Number of members
integer,intent(in) :: ncyc                                !< Number of cycles
real(kind_real),intent(in),optional :: ens1_2d(:,:,:,:)   !< Ensemble 1, 2d
real(kind_real),intent(in),optional :: ens1_3d(:,:,:,:,:) !< Ensemble 1, 3d

! Local variables
integer :: nmga,nts,nv,ens1_ne,il0
real(kind_real),allocatable :: lon_mga(:),lat_mga(:),area_mga(:),vunit_mga(:,:)
logical,allocatable :: lmask_mga(:,:)

! Initialize MPL
call bump%mpl%init(mpi_comm)

! Copy sizes
bump%nx = nx
bump%ny = ny

! Sizes
nmga = nx*ny
nv = 1
nts = 1
ens1_ne = nens*ncyc

! Set internal namelist parameters
call bump%nam%setup_internal(nl0,nv,nts,ens1_ne)

! Initialize listing
call bump%mpl%init_listing(bump%nam%prefix,bump%nam%model,bump%nam%colorlog,bump%nam%logpres,lunit=lunit)

! Generic setup
call bump%setup_generic

! Allocation
allocate(lon_mga(nmga))
allocate(lat_mga(nmga))
allocate(area_mga(nmga))
allocate(vunit_mga(nmga,nl0))
allocate(lmask_mga(nmga,nl0))

! Pack
lon_mga = pack(lon,.true.)
lat_mga = pack(lat,.true.)
area_mga = pack(area,.true.)
do il0=1,nl0
   vunit_mga(:,il0) = pack(vunit(:,:,il0),.true.)
   lmask_mga(:,il0) = pack(lmask(:,:,il0),.true.)
end do

! Initialize geometry
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize geometry'
call flush(bump%mpl%unit)
call bump%geom%setup_online(bump%mpl,nmga,nl0,lon_mga,lat_mga,area_mga,vunit_mga,lmask_mga)
call bump%geom%init(bump%mpl,bump%rng,bump%nam)

if (bump%nam%grid_output) then
   ! Initialize fields regridding
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Initialize fields regridding'
   call flush(bump%mpl%unit)
   call bump%io%grid_init(bump%mpl,bump%rng,bump%nam,bump%geom)
end if

! Initialize block parameters
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize block parameters'
call bump%bpar%alloc(bump%nam,bump%geom)

! Initialize ensemble 1
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize ensemble 1'
call flush(bump%mpl%unit)
if (present(ens1_2d)) then
   call bump%ens1%from(bump%mpl,bump%nam,bump%geom,nx,ny,nens,ncyc,ens_2d=ens1_2d)
elseif (present(ens1_3d)) then
   call bump%ens1%from(bump%mpl,bump%nam,bump%geom,nx,ny,nens,ncyc,ens_3d=ens1_3d)
end if

! Run drivers
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Run drivers'
call flush(bump%mpl%unit)
call bump%run_drivers

if (bump%mpl%main) then
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- BUMP done'
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
else
   ! Close listings
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Close listings'
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   call flush(bump%mpl%unit)
   close(unit=bump%mpl%unit)
end if

end subroutine bump_setup_online_nemovar

!----------------------------------------------------------------------
! Subroutine: bump_setup_generic
!> Purpose: generic setup
!----------------------------------------------------------------------
subroutine bump_setup_generic(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump !< BUMP

! Header
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- You are running bump ------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Author: Benjamin Menetrier ------------------------------------'
write(bump%mpl%unit,'(a)') '--- Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE -----------'
call flush(bump%mpl%unit)

! Check namelist parameters
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Check namelist parameters'
call flush(bump%mpl%unit)
call bump%nam%check(bump%mpl)

! Write parallel setup
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a,i3,a,i2,a)') '--- Parallelization with ',bump%mpl%nproc,' MPI tasks and ', &
 & bump%mpl%nthread,' OpenMP threads'
call flush(bump%mpl%unit)

! Initialize random number generator
write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%unit,'(a)') '--- Initialize random number generator'
call flush(bump%mpl%unit)
call bump%rng%init(bump%mpl,bump%nam)

! Initialize allocation flags
bump%cmat%allocated = .false.
bump%lct%allocated = .false.
bump%nicas%allocated = .false.
bump%obsop%allocated = .false.

end subroutine bump_setup_generic

!----------------------------------------------------------------------
! Subroutine: bump_run_drivers
!> Purpose: run drivers
!----------------------------------------------------------------------
subroutine bump_run_drivers(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump !< BUMP

! Reset seed
if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

! Check inconsistencies
if (bump%nam%new_hdiag.and.(.not.allocated(bump%ens1%fld))) call bump%mpl%abort('new_hdiag requires ensemble 1')
if (bump%nam%new_hdiag.and.(trim(bump%nam%method)=='hyb-rnd').or.(trim(bump%nam%method)=='dual-ens') &
 & .and.(.not.allocated(bump%ens2%fld))) call bump%mpl%abort('new_hdiag requires ensemble 1')
if (bump%nam%new_hdiag.and.(allocated(bump%rh).and.allocated(bump%rv))) &
 & call bump%mpl%abort('rh and rv should not be provided if new_hdiag is active')

if (.not.bump%cmat%allocated) then
   if (bump%nam%new_hdiag) then
      ! Call HDIAG driver
      write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(bump%mpl%unit,'(a)') '--- Run HDIAG driver'
      call flush(bump%mpl%unit)
      if ((trim(bump%nam%method)=='hyb-rnd').or.(trim(bump%nam%method)=='dual-ens')) then
         call bump%cmat%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,bump%ens2)
      else
         call bump%cmat%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
      end if
   elseif (bump%nam%new_param.and..not.(allocated(bump%rh).and.allocated(bump%rv))) then
      if (bump%nam%forced_radii) then
         ! Copy namelist support radii into C matrix
         call bump%cmat%from(bump%mpl,bump%nam,bump%geom,bump%bpar)
      else
         ! Read C matrix data
         write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(bump%mpl%unit,'(a)') '--- Read C matrix data'
         call flush(bump%mpl%unit)
         call bump%cmat%read(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io)
      end if
   end if
end if

if (.not.bump%nicas%allocated) then
   if (bump%nam%new_param) then
      ! Call NICAS driver
      write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(bump%mpl%unit,'(a)') '--- Call NICAS driver'
      call flush(bump%mpl%unit)
      call bump%nicas%run_nicas(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%cmat)
   elseif (bump%nam%check_adjoints.or.bump%nam%check_pos_def.or.bump%nam%check_sqrt.or.bump%nam%check_dirac.or. &
    & bump%nam%check_randomization.or.bump%nam%check_consistency.or.bump%nam%check_optimality) then
      ! Read NICAS parameters
      write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(bump%mpl%unit,'(a)') '--- Read NICAS parameters'
      call flush(bump%mpl%unit)
      call bump%nicas%read(bump%mpl,bump%nam,bump%geom,bump%bpar)
   end if
end if

if (bump%nam%check_adjoints.or.bump%nam%check_pos_def.or.bump%nam%check_sqrt.or.bump%nam%check_dirac.or. &
 & bump%nam%check_randomization.or.bump%nam%check_consistency.or.bump%nam%check_optimality) then
   ! Call NICAS tests driver
   write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%unit,'(a)') '--- Call NICAS tests driver'
   call flush(bump%mpl%unit)
   if (allocated(bump%ens1%fld)) then
      call bump%nicas%run_nicas_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%cmat,bump%ens1)
   else
      call bump%nicas%run_nicas_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%cmat)
   end if
end if

if (.not.bump%lct%allocated) then
   if (bump%nam%new_lct) then
      ! Call LCT driver
      write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(bump%mpl%unit,'(a)') '--- Call LCT driver'
      call flush(bump%mpl%unit)
      call bump%lct%run_lct(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
   end if
end if

if (.not.bump%obsop%allocated) then
   if (bump%nam%new_obsop) then
      ! Call observation operator driver
      write(bump%mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(bump%mpl%unit,'(a)') '--- Call observation operator driver'
      call flush(bump%mpl%unit)
      call bump%obsop%run_obsop(bump%mpl,bump%rng,bump%nam,bump%geom)
   end if
end if

end subroutine bump_run_drivers

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas
!> Purpose: NICAS application
!----------------------------------------------------------------------
subroutine bump_apply_nicas(bump,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                         !< BUMP
real(kind_real),intent(inout) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

! Apply NICAS
if (bump%nam%lsqrt) then
   call bump%nicas%apply_from_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,fld)
else
   call bump%nicas%apply(bump%mpl,bump%nam,bump%geom,bump%bpar,fld)
end if

end subroutine bump_apply_nicas

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop
!> Purpose: observation operator application
!----------------------------------------------------------------------
subroutine bump_apply_obsop(bump,fld,obs)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                !< BUMP
real(kind_real),intent(in) :: fld(bump%geom%nc0a,bump%geom%nl0)    !< Field
real(kind_real),intent(out) :: obs(bump%obsop%nobsa,bump%geom%nl0) !< Observations columns

! Apply observation operator
if (bump%obsop%nobsa>0) call bump%obsop%apply(bump%mpl,bump%geom,fld,obs)

end subroutine bump_apply_obsop

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop_ad
!> Purpose: observation operator adjoint application
!----------------------------------------------------------------------
subroutine bump_apply_obsop_ad(bump,obs,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                               !< BUMP
real(kind_real),intent(in) :: obs(bump%obsop%nobsa,bump%geom%nl0) !< Observations columns
real(kind_real),intent(out) :: fld(bump%geom%nc0a,bump%geom%nl0)  !< Field

! Apply observation operator adjoint
if (bump%obsop%nobsa>0) then
   call bump%obsop%apply_ad(bump%mpl,bump%geom,obs,fld)
else
   fld = 0.0
end if

end subroutine bump_apply_obsop_ad

!----------------------------------------------------------------------
! Subroutine: bump_back_to
!> Purpose: conversion from subset Sc0 to XY model grid, halo A
!----------------------------------------------------------------------
subroutine bump_back_to(bump,fld_c0a,fld_xya)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                     !< BUMP
real(kind_real),intent(in) :: fld_c0a(bump%geom%nc0a)   !< Field on subset Sc0, halo A
real(kind_real),intent(out) :: fld_xya(bump%nx,bump%ny) !< Field on XY model grid, halo A

! Local variables
real(kind_real) :: fld_mga(bump%geom%nmga)
logical :: mask_unpack(bump%nx,bump%ny)

! Halo extension from subset Sc0 to model grid, halo A
call bump%geom%com_mg%ext(bump%mpl,fld_c0a,fld_mga)

! Unpack to XY grid
mask_unpack = .true.
fld_xya = unpack(fld_mga,mask_unpack,fld_xya)

end subroutine bump_back_to

!----------------------------------------------------------------------
! Subroutine: bump_dealloc
!> Purpose: deallocation of BUMP fields
!----------------------------------------------------------------------
subroutine bump_dealloc(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump !< BUMP

! Release memory
call bump%bpar%dealloc
call bump%cmat%dealloc(bump%bpar)
call bump%ens1%dealloc
call bump%ens2%dealloc
call bump%geom%dealloc
call bump%io%dealloc
call bump%lct%dealloc(bump%bpar)
call bump%nicas%dealloc(bump%nam,bump%geom,bump%bpar)
call bump%obsop%dealloc
if (allocated(bump%rh)) deallocate(bump%rh)
if (allocated(bump%rv)) deallocate(bump%rv)

end subroutine bump_dealloc

end module type_bump
