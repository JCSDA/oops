!----------------------------------------------------------------------
! Module: type_bump
! Purpose: BUMP derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_bump

use netcdf
use tools_const, only: req,deg2rad
use tools_func, only: sphere_dist,lct_r2d
use tools_kinds,only: kind_real
use type_bpar, only: bpar_type
use type_cmat, only: cmat_type
use type_cv, only: cv_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_hdiag, only: hdiag_type
use type_io, only: io_type
use type_lct, only: lct_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_nicas, only: nicas_type
use type_obsop, only: obsop_type
use type_rng, only: rng_type
use type_vbal, only: vbal_type

implicit none

! BUMP derived type
type bump_type
   type(bpar_type) :: bpar
   type(cmat_type) :: cmat
   type(ens_type) :: ens1
   type(ens_type) :: ens1u
   type(ens_type) :: ens2
   type(geom_type) :: geom
   type(hdiag_type) :: hdiag
   type(io_type) :: io
   type(lct_type) :: lct
   type(mpl_type) :: mpl
   type(nam_type) :: nam
   type(nicas_type) :: nicas
   type(obsop_type) :: obsop
   type(rng_type) :: rng
   type(vbal_type) :: vbal
contains
   procedure :: setup_online => bump_setup_online
   procedure :: run_drivers => bump_run_drivers
   procedure :: add_member => bump_add_member
   procedure :: apply_vbal => bump_apply_vbal
   procedure :: apply_vbal_inv => bump_apply_vbal_inv
   procedure :: apply_vbal_ad => bump_apply_vbal_ad
   procedure :: apply_vbal_inv_ad => bump_apply_vbal_inv_ad
   procedure :: apply_nicas => bump_apply_nicas
   procedure :: get_cv_size => bump_get_cv_size
   procedure :: apply_nicas_sqrt => bump_apply_nicas_sqrt
   procedure :: apply_nicas_sqrt_ad => bump_apply_nicas_sqrt_ad
   procedure :: randomize => bump_randomize
   procedure :: apply_obsop => bump_apply_obsop
   procedure :: apply_obsop_ad => bump_apply_obsop_ad
   procedure :: get_parameter => bump_get_parameter
   procedure :: copy_to_field => bump_copy_to_field
   procedure :: set_parameter => bump_set_parameter
   procedure :: copy_from_field => bump_copy_from_field
   procedure :: bump_crtm_neighbors_3d
   procedure :: bump_crtm_neighbors_2d
   generic :: crtm_neighbors => bump_crtm_neighbors_3d,bump_crtm_neighbors_2d
   procedure :: dealloc => bump_dealloc
   procedure :: partial_dealloc => bump_partial_dealloc
end type bump_type

private
public :: bump_type

contains

!----------------------------------------------------------------------
! Subroutine: bump_setup_online
! Purpose: online setup
!----------------------------------------------------------------------
subroutine bump_setup_online(bump,nmga,nl0,nv,nts,lon,lat,area,vunit,gmask,smask,ens1_ne,ens1_nsub,ens2_ne,ens2_nsub, &
                           & nobs,lonobs,latobs,namelname,lunit,msvali,msvalr)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump            ! BUMP
integer,intent(in) :: nmga                        ! Halo A size
integer,intent(in) :: nl0                         ! Number of levels in subset Sl0
integer,intent(in) :: nv                          ! Number of variables
integer,intent(in) :: nts                         ! Number of time slots
real(kind_real),intent(in) :: lon(nmga)           ! Longitude (in degrees: -180 to 180)
real(kind_real),intent(in) :: lat(nmga)           ! Latitude (in degrees: -90 to 90)
real(kind_real),intent(in) :: area(nmga)          ! Area (in m^2)
real(kind_real),intent(in) :: vunit(nmga,nl0)     ! Vertical unit
logical,intent(in) :: gmask(nmga,nl0)             ! Geometry mask
logical,intent(in),optional :: smask(nmga,nl0)    ! Sampling mask
integer,intent(in),optional :: ens1_ne            ! Ensemble 1 size
integer,intent(in),optional :: ens1_nsub          ! Ensemble 1 number of sub-ensembles
integer,intent(in),optional :: ens2_ne            ! Ensemble 2 size
integer,intent(in),optional :: ens2_nsub          ! Ensemble 2 size of sub-ensembles
integer,intent(in),optional :: nobs               ! Number of observations
real(kind_real),intent(in),optional :: lonobs(:)  ! Observations longitude (in degrees: -180 to 180)
real(kind_real),intent(in),optional :: latobs(:)  ! Observations latitude (in degrees: -90 to 90)
character(len=*),intent(in),optional :: namelname ! Namelist name
integer,intent(in),optional :: lunit              ! Listing unit
integer,intent(in),optional :: msvali             ! Missing value for integers
real(kind_real),intent(in),optional :: msvalr     ! Missing value for reals

! Local variables
integer :: lmsvali,lens1_ne,lens1_nsub,lens2_ne,lens2_nsub
real(kind_real) :: lmsvalr
logical :: lgmask(nmga,nl0)
character(len=1024),parameter :: subr = 'bump_setup_online'

! Set missing values
lmsvali = -999
if (present(msvali)) lmsvali = msvali
lmsvalr = -999.0
if (present(msvalr)) lmsvalr = msvalr
call bump%mpl%msv%init(lmsvali,lmsvalr)

! Initialize MPL
call bump%mpl%init

if (present(namelname)) then
   ! Read and broadcast namelist
   call bump%nam%read(bump%mpl,namelname)
   call bump%nam%bcast(bump%mpl)
end if

! Set internal namelist parameters
lens1_ne = 0
lens1_nsub = 1
lens2_ne = 0
lens2_nsub = 1
if (present(ens1_ne)) lens1_ne = ens1_ne
if (present(ens1_nsub)) lens1_nsub = ens1_nsub
if (present(ens2_ne)) lens2_ne = ens2_ne
if (present(ens2_ne)) lens2_nsub = ens2_nsub
call bump%nam%setup_internal(nl0,nv,nts,lens1_ne,lens1_nsub,lens2_ne,lens2_nsub)

! Initialize listing
bump%mpl%lunit = bump%mpl%msv%vali
if (present(lunit)) bump%mpl%lunit = lunit
bump%mpl%verbosity = bump%nam%verbosity
if (bump%nam%colorlog) then
   bump%mpl%black = char(27)//'[0;0m'
   bump%mpl%green = char(27)//'[0;32m'
   bump%mpl%peach = char(27)//'[1;91m'
   bump%mpl%aqua = char(27)//'[1;36m'
   bump%mpl%purple = char(27)//'[1;35m'
   bump%mpl%err = char(27)//'[0;37;41;1m'
   bump%mpl%wng = char(27)//'[0;37;42;1m'
else
   bump%mpl%black = ' '
   bump%mpl%green = ' '
   bump%mpl%peach = ' '
   bump%mpl%aqua = ' '
   bump%mpl%purple = ' '
   bump%mpl%err = ' '
   bump%mpl%wng = ' '
end if

! Header
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- You are running the BUMP library ------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Author: Benjamin Menetrier ------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT -----'
call bump%mpl%flush

! Check namelist parameters
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Check namelist parameters'
call bump%mpl%flush
call bump%nam%check(bump%mpl)
call bump%nam%write(bump%mpl)

! Write parallel setup
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a,i3,a,i2,a)') '--- Parallelization with ',bump%mpl%nproc,' MPI tasks and ', &
 & bump%mpl%nthread,' OpenMP threads'
call bump%mpl%flush

! Initialize random number generator
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Initialize random number generator'
call bump%mpl%flush
call bump%rng%init(bump%mpl,bump%nam)

! Initialize allocation flags
bump%cmat%allocated = .false.
bump%lct%allocated = .false.
bump%nicas%allocated = .false.

! Initialize geometry
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Initialize geometry'
call bump%mpl%flush
lgmask = gmask.or.bump%nam%nomask
call bump%geom%setup(bump%mpl,bump%rng,bump%nam,nmga,nl0,lon,lat,area,vunit,lgmask)
if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

if (bump%nam%grid_output) then
   ! Initialize fields regridding
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize fields regridding'
   call bump%mpl%flush
   call bump%io%grid_init(bump%mpl,bump%rng,bump%nam,bump%geom)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

! Initialize block parameters
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Initialize block parameters'
call bump%mpl%flush
call bump%bpar%alloc(bump%nam,bump%geom)
call bump%bpar%init(bump%mpl,bump%nam,bump%geom)

if (bump%nam%ens1_ne>0) then
   ! Initialize ensemble 1
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize ensemble 1'
   call bump%mpl%flush
   call bump%ens1%alloc(bump%nam,bump%geom,bump%nam%ens1_ne,bump%nam%ens1_nsub)
else
   call bump%ens1%set_att(bump%nam%ens1_ne,bump%nam%ens1_nsub)
end if

if (bump%nam%ens2_ne>0) then
   ! Initialize ensemble 2
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize ensemble 2'
   call bump%mpl%flush
   call bump%ens2%alloc(bump%nam,bump%geom,bump%nam%ens2_ne,bump%nam%ens2_nsub)
else
   call bump%ens2%set_att(bump%nam%ens2_ne,bump%nam%ens2_nsub)
end if

if (present(nobs)) then
   ! Check arguments consistency
   if ((.not.present(lonobs)).or.(.not.present(latobs))) call bump%mpl%abort(subr,'lonobs and latobs are missing')

   ! Check sizes consistency
   if (size(lonobs)/=nobs) call bump%mpl%abort(subr,'wrong size for lonobs')
   if (size(latobs)/=nobs) call bump%mpl%abort(subr,'wrong size for latobs')

   ! Initialize observations locations
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize observations locations'
   call bump%mpl%flush
   call bump%obsop%from(nobs,lonobs,latobs)
end if

! Copy sampling mask
if (present(smask)) then
   allocate(bump%geom%smask_c0a(bump%geom%nc0a,bump%geom%nl0))
   bump%geom%smask_c0a = smask(bump%geom%c0a_to_mga,:)
end if

end subroutine bump_setup_online

!----------------------------------------------------------------------
! Subroutine: bump_run_drivers
! Purpose: run drivers
!----------------------------------------------------------------------
subroutine bump_run_drivers(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

if (bump%nam%ens1_ne>0) then
   ! Finalize ensemble 1
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Finalize ensemble 1'
   call bump%mpl%flush
   call bump%ens1%remove_mean
end if

if (bump%nam%ens2_ne>0) then
   ! Finalize ensemble 2
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Finalize ensemble 2'
   call bump%mpl%flush
   call bump%ens2%remove_mean
end if

if (bump%nam%new_cortrack) then
   ! Run correlation tracker
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run correlation tracker'
   call bump%mpl%flush
   call bump%ens1%cortrack(bump%mpl,bump%rng,bump%nam,bump%geom,bump%io)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

if (bump%nam%new_vbal) then
   ! Run vertical balance driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run vertical balance driver'
   call bump%mpl%flush
   call bump%vbal%run_vbal(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,bump%ens1u)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
elseif (bump%nam%load_vbal) then
   ! Read vertical balance
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read vertical balance'
   call bump%mpl%flush
   call bump%vbal%read(bump%mpl,bump%nam,bump%geom,bump%bpar)
end if

if (bump%nam%check_vbal) then
   ! Run vertical balance tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run vertical balance tests driver'
   call bump%mpl%flush
   call bump%vbal%run_vbal_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

if (bump%nam%new_hdiag) then
   ! Run HDIAG driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run HDIAG driver'
   call bump%mpl%flush
   if ((trim(bump%nam%method)=='hyb-rnd').or.(trim(bump%nam%method)=='dual-ens')) then
      call bump%hdiag%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,bump%ens2)
   else
      call bump%hdiag%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
   end if
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Copy HDIAG into C matrix
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Copy HDIAG into C matrix'
   call bump%mpl%flush
   call bump%cmat%from_hdiag(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%hdiag)
end if

if (bump%nam%new_lct) then
   ! Run LCT driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run LCT driver'
   call bump%mpl%flush
   call bump%lct%run_lct(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Copy LCT into C matrix
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Copy LCT into C matrix'
   call bump%mpl%flush
   call bump%cmat%from_lct(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%lct)
end if

if (bump%nam%load_cmat) then
   ! Read C matrix
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read C matrix'
   call bump%mpl%flush
   call bump%cmat%read(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io)
else
   if (bump%nam%forced_radii) then
      ! Copy namelist support radii into C matrix
      write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
      call bump%mpl%flush
      write(bump%mpl%info,'(a)') '--- Copy namelist support radii into C matrix'
      call bump%mpl%flush
      call bump%cmat%from_nam(bump%mpl,bump%nam,bump%geom,bump%bpar)
   end if
end if

if (bump%cmat%allocated.or.bump%nam%new_nicas) then
   ! Get C matrix from BUMP interface
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Get C matrix from BUMP interface'
   call bump%mpl%flush
   call bump%cmat%from_bump(bump%mpl,bump%nam,bump%geom,bump%bpar)

   ! Setup C matrix sampling
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Setup C matrix sampling'
   call bump%mpl%flush
   call bump%cmat%setup_sampling(bump%nam,bump%geom,bump%bpar)

   if (bump%nam%write_cmat) then
      ! Write C matrix
      write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
      call bump%mpl%flush
      write(bump%mpl%info,'(a)') '--- Write C matrix'
      call bump%mpl%flush
      call bump%cmat%write(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io)
   end if
end if

if (bump%nam%new_nicas) then
   ! Run NICAS driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run NICAS driver'
   call bump%mpl%flush
   call bump%nicas%run_nicas(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%cmat)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
elseif (bump%nam%load_nicas) then
   ! Read NICAS parameters
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read NICAS parameters'
   call bump%mpl%flush
   call bump%nicas%read(bump%mpl,bump%nam,bump%geom,bump%bpar)
end if

if (bump%nam%check_adjoints.or.bump%nam%check_pos_def.or.bump%nam%check_dirac.or.bump%nam%check_randomization.or. &
 & bump%nam%check_consistency.or.bump%nam%check_optimality) then
   ! Run NICAS tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run NICAS tests driver'
   call bump%mpl%flush
   call bump%nicas%run_nicas_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

if (bump%nam%new_obsop) then
   ! Run observation operator driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run observation operator driver'
   call bump%mpl%flush
   call bump%obsop%run_obsop(bump%mpl,bump%rng,bump%nam,bump%geom)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
elseif (bump%nam%load_obsop) then
   ! Read observation operator
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read observation operator'
   call bump%mpl%flush
   call bump%obsop%read(bump%mpl,bump%nam)
end if

if (bump%nam%check_obsop) then
   ! Run observation operator tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run observation operator tests driver'
   call bump%mpl%flush
   call bump%obsop%run_obsop_tests(bump%mpl,bump%rng,bump%geom)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

end subroutine bump_run_drivers

!----------------------------------------------------------------------
! Subroutine: bump_add_member
! Purpose: add member into bump%ens[1,2]
!----------------------------------------------------------------------
subroutine bump_add_member(bump,fld_mga,ie,iens)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field
integer,intent(in) :: ie                                                                        ! Member index
integer,intent(in) :: iens                                                                      ! Ensemble number

! Local variables
integer :: its,iv,nnonzero,nzero,nmask
real(kind_real) :: norm,fld_c0a(bump%geom%nc0a,bump%geom%nl0)
character(len=1024),parameter :: subr = 'bump_add_member'

! Add member
write(bump%mpl%info,'(a7,a,i3,a,i1)') '','Member ',ie,' added to ensemble ',iens
do its=1,bump%nam%nts
   do iv=1,bump%nam%nv
      ! Model grid to subset Sc0
      call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a)

      ! Copy to ensemble structure
      if (iens==1) then
         bump%ens1%fld(:,:,iv,its,ie) = fld_c0a
      elseif (iens==2) then
         bump%ens2%fld(:,:,iv,its,ie) = fld_c0a
      else
         call bump%mpl%abort(subr,'wrong ensemble number')
      end if

      ! Print norm
      norm = sum(fld_c0a**2,mask=bump%geom%mask_c0a)
      write(bump%mpl%info,'(a10,a,i2,a,i2,a,e9.2)') '','Local norm for variable ',iv,' and timeslot ',its,': ',norm
      nnonzero = count((abs(fld_c0a)>0.0).and.bump%geom%mask_c0a)
      nzero = count((.not.(abs(fld_c0a)>0.0)).and.bump%geom%mask_c0a)
      nmask = count(.not.bump%geom%mask_c0a)
      write(bump%mpl%info,'(a10,a,i8,a,i8,a,i8,a,i8)') '','Total / non-zero / zero / masked points: ',bump%geom%nc0a,' / ', &
    & nnonzero,' / ',nzero,' / ',nmask
   end do
end do

end subroutine bump_add_member

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal
! Purpose: vertical balance application
!----------------------------------------------------------------------
subroutine bump_apply_vbal(bump,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)

do its=1,bump%nam%nts
   if (bump%geom%nc0==bump%geom%nmg) then
       ! Apply vertical balance
      call bump%vbal%apply(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance
      call bump%vbal%apply(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

end subroutine bump_apply_vbal

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv
! Purpose: vertical balance application, inverse
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv(bump,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)

do its=1,bump%nam%nts
   if (bump%geom%nc0==bump%geom%nmg) then
      ! Apply vertical balance, inverse
      call bump%vbal%apply_inv(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance, inverse
      call bump%vbal%apply_inv(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

end subroutine bump_apply_vbal_inv

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_ad
! Purpose: vertical balance application, adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_ad(bump,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)

do its=1,bump%nam%nts
   if (bump%geom%nc0==bump%geom%nmg) then
      ! Apply vertical balance, adjoint
      call bump%vbal%apply_ad(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance, adjoint
      call bump%vbal%apply_ad(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

end subroutine bump_apply_vbal_ad

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv_ad
! Purpose: vertical balance application, inverse adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv_ad(bump,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)

do its=1,bump%nam%nts
   if (bump%geom%nc0==bump%geom%nmg) then
      ! Apply vertical balance, inverse adjoint
      call bump%vbal%apply_inv_ad(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance, inverse adjoint
      call bump%vbal%apply_inv_ad(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

end subroutine bump_apply_vbal_inv_ad

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas
! Purpose: NICAS application
!----------------------------------------------------------------------
subroutine bump_apply_nicas(bump,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

if (bump%geom%nc0==bump%geom%nmg) then
   ! Apply NICAS
   if (bump%nam%lsqrt) then
      call bump%nicas%apply_from_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga)
   else
      call bump%nicas%apply(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga)
   end if
else
   ! Model grid to subset Sc0
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv,its))
      end do
   end do

   ! Apply NICAS
   if (bump%nam%lsqrt) then
      call bump%nicas%apply_from_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a)
   else
      call bump%nicas%apply(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a)
   end if

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

end subroutine bump_apply_nicas

!----------------------------------------------------------------------
! Subroutine: bump_get_cv_size
! Purpose: get control variable size
!----------------------------------------------------------------------
subroutine bump_get_cv_size(bump,n)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP
integer,intent(out) :: n               ! Control variable size

! Local variables
type(cv_type) :: cv

! Allocate control variable
call bump%nicas%alloc_cv(bump%mpl,bump%bpar,cv,getsizeonly=.true.)

! Copy size
n = cv%n

end subroutine bump_get_cv_size

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt
! Purpose: NICAS square-root application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt(bump,pcv,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(in) :: pcv(:)                                                            ! Packed control variable
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024),parameter :: subr = 'bump_apply_nicas_sqrt'
type(cv_type) :: cv

! Allocation
call bump%nicas%alloc_cv(bump%mpl,bump%bpar,cv)

! Check dimension
if (size(pcv)==cv%n) then
   ! Unpack control variable
   call cv%unpack(pcv)
else
   call bump%mpl%abort(subr,'wrong control variable size in bump_apply_nicas_sqrt')
end if

if (bump%geom%nc0==bump%geom%nmg) then
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_mga)
else
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_c0a)

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

end subroutine bump_apply_nicas_sqrt

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt_ad
! Purpose: NICAS square-root adjoint application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt_ad(bump,fld_mga,pcv)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                   ! BUMP
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field
real(kind_real),intent(inout) :: pcv(:)                                                  ! Packed control variable

! Local variables
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024),parameter :: subr = 'bump_apply_nicas_sqrt_ad'
type(cv_type) :: cv

if (bump%geom%nc0==bump%geom%nmg) then
   ! Apply NICAS square-root adjoint
   call bump%nicas%apply_sqrt_ad(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga,cv)
else
   ! Model grid to subset Sc0
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv,its))
      end do
   end do

   ! Apply NICAS square-root adjoint
   call bump%nicas%apply_sqrt_ad(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a,cv)
end if

! Check dimension
if (size(pcv)==cv%n) then
   ! Pack control variable
   call cv%pack(pcv)
else
   call bump%mpl%abort(subr,'wrong control variable size in bump_apply_nicas_sqrt_ad')
end if

end subroutine bump_apply_nicas_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: bump_randomize
! Purpose: NICAS randomization
!----------------------------------------------------------------------
subroutine bump_randomize(bump,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                        ! BUMP
real(kind_real),intent(out) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
type(cv_type) :: cv

! Generate random control vector
call bump%nicas%random_cv(bump%mpl,bump%rng,bump%bpar,cv)

if (bump%geom%nc0==bump%geom%nmg) then
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_mga)
else
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_c0a)

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

end subroutine bump_randomize

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop
! Purpose: observation operator application
!----------------------------------------------------------------------
subroutine bump_apply_obsop(bump,fld_mga,obs)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                              ! BUMP
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field
real(kind_real),intent(out) :: obs(bump%obsop%nobsa,bump%geom%nl0)  ! Observations columns

! Local variables
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0)

if (bump%geom%nc0==bump%geom%nmg) then
   ! Apply observation operator
   call bump%obsop%apply(bump%mpl,bump%geom,fld_mga,obs)
else
   ! Model grid to subset Sc0
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,fld_c0a)

   ! Apply observation operator
   call bump%obsop%apply(bump%mpl,bump%geom,fld_c0a,obs)
end if

end subroutine bump_apply_obsop

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop_ad
! Purpose: observation operator adjoint application
!----------------------------------------------------------------------
subroutine bump_apply_obsop_ad(bump,obs,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                               ! BUMP
real(kind_real),intent(in) :: obs(bump%obsop%nobsa,bump%geom%nl0)    ! Observations columns
real(kind_real),intent(out) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field

! Local variables
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0)

if (bump%geom%nc0==bump%geom%nmg) then
   ! Apply observation operator adjoint
   call bump%obsop%apply_ad(bump%mpl,bump%geom,obs,fld_mga)
else
   ! Apply observation operator adjoint
   call bump%obsop%apply_ad(bump%mpl,bump%geom,obs,fld_c0a)

   ! Subset Sc0 to model grid
   call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a,fld_mga)
end if

end subroutine bump_apply_obsop_ad

!----------------------------------------------------------------------
! Subroutine: bump_get_parameter
! Purpose: get a parameter
!----------------------------------------------------------------------
subroutine bump_get_parameter(bump,param,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                        ! BUMP
character(len=*),intent(in) :: param                                                          ! Parameter
real(kind_real),intent(out) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variables
integer :: ib,iv,jv,its,jts

select case (trim(param))
case ('var','cor_rh','cor_rv','cor_rv_rfac','cor_rv_coef','loc_coef','loc_rh','loc_rv','hyb_coef')
   select case (trim(bump%nam%strategy))
   case ('specific_univariate','specific_multivariate')
      do ib=1,bump%bpar%nb
         ! Get indices
         iv = bump%bpar%b_to_v1(ib)
         jv = bump%bpar%b_to_v2(ib)
         its = bump%bpar%b_to_ts1(ib)
         jts = bump%bpar%b_to_ts2(ib)

         ! Copy to field
         if ((iv==jv).and.(its==jts)) call bump%copy_to_field(param,ib,fld_mga(:,:,iv,its))
      end do
   case ('common','common_univariate','common_weighted')
      ! Set common index
      ib = bump%bpar%nbe

      do its=1,bump%nam%nts
         do iv=1,bump%nam%nv
            ! Copy to field
            call bump%copy_to_field(param,ib,fld_mga(:,:,iv,its))
         end do
      end do
   end select
case default
   do ib=1,bump%bpar%nb
      ! Get indices
      iv = bump%bpar%b_to_v1(ib)
      jv = bump%bpar%b_to_v2(ib)
      its = bump%bpar%b_to_ts1(ib)
      jts = bump%bpar%b_to_ts2(ib)

      ! Copy to field
      if ((iv==jv).and.(its==jts)) call bump%copy_to_field(param,ib,fld_mga(:,:,iv,its))
   end do
end select

end subroutine bump_get_parameter

!----------------------------------------------------------------------
! Subroutine: bump_copy_to_field
! Purpose: copy to field
!----------------------------------------------------------------------
subroutine bump_copy_to_field(bump,param,ib,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                               ! BUMP
character(len=*),intent(in) :: param                                 ! Parameter
integer,intent(in) :: ib                                             ! Block index
real(kind_real),intent(out) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field

! Local variables
integer :: iscales,ie,imga,il0,iv,its
real(kind_real) :: tmp
character(len=1024),parameter :: subr = 'bump_copy_to_field'

! Select parameter
select case (trim(param))
case ('var')
   if (.not.allocated(bump%cmat%blk(ib)%coef_ens)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_ens,fld_mga)
case ('cor_rh')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
      end do
   end do
case ('cor_rv')
   if (.not.allocated(bump%cmat%blk(ib)%rv)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld_mga)
case ('cor_rv_rfac')
   if (.not.allocated(bump%cmat%blk(ib)%rv_rfac)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv_rfac,fld_mga)
case ('cor_rv_coef')
   if (.not.allocated(bump%cmat%blk(ib)%rv_coef)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv_coef,fld_mga)
case ('loc_coef')
   if (.not.allocated(bump%cmat%blk(ib)%coef_ens)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_ens,fld_mga)
case ('loc_rh')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
      end do
   end do
case ('loc_rv')
   if (.not.allocated(bump%cmat%blk(ib)%rv)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld_mga)
case ('hyb_coef')
   if (.not.allocated(bump%cmat%blk(ib)%coef_sta)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_sta,fld_mga)
case ('loc_D11','loc_D22')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
      end do
   end do
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         tmp = fld_mga(imga,il0)
         call lct_r2d(tmp,fld_mga(imga,il0))
      end do
   end do
case ('loc_D12')
   fld_mga = 0.0
case ('loc_D33')
   if (.not.allocated(bump%cmat%blk(ib)%rv)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         tmp = fld_mga(imga,il0)
         call lct_r2d(tmp,fld_mga(imga,il0))
      end do
   end do
case default
   select case (param(1:4))
   case ('D11_')
      if (.not.allocated(bump%lct%blk(ib)%D11)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
      read(param(5:5),'(i1)') iscales
      if (iscales>size(bump%lct%blk(ib)%D11,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D11(:,:,iscales),fld_mga)
      do il0=1,bump%geom%nl0
         do imga=1,bump%geom%nmga
            if (bump%mpl%msv%isnotr(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req**2
         end do
      end do
   case ('D22_')
      if (.not.allocated(bump%lct%blk(ib)%D22)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
      read(param(5:5),'(i1)') iscales
      if (iscales>size(bump%lct%blk(ib)%D22,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D22(:,:,iscales),fld_mga)
      do il0=1,bump%geom%nl0
         do imga=1,bump%geom%nmga
            if (bump%mpl%msv%isnotr(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req**2
         end do
      end do
   case ('D33_')
      if (.not.allocated(bump%lct%blk(ib)%D33)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
      read(param(5:5),'(i1)') iscales
      if (iscales>size(bump%lct%blk(ib)%D33,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D33(:,:,iscales),fld_mga)
   case ('D12_')
      if (.not.allocated(bump%lct%blk(ib)%D12)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
      read(param(5:5),'(i1)') iscales
      if (iscales>size(bump%lct%blk(ib)%D12,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D12(:,:,iscales),fld_mga)
      do il0=1,bump%geom%nl0
         do imga=1,bump%geom%nmga
            if (bump%mpl%msv%isnotr(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req**2
         end do
      end do
   case ('Dcoe')
      if (.not.allocated(bump%lct%blk(ib)%Dcoef)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
      read(param(7:7),'(i1)') iscales
      if (iscales>size(bump%lct%blk(ib)%Dcoef,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%Dcoef(:,:,iscales),fld_mga)
   case ('DLh_')
      if (.not.allocated(bump%lct%blk(ib)%DLh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
      read(param(5:5),'(i1)') iscales
      if (iscales>size(bump%lct%blk(ib)%DLh,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%DLh(:,:,iscales),fld_mga)
      do il0=1,bump%geom%nl0
         do imga=1,bump%geom%nmga
            if (bump%mpl%msv%isnotr(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
         end do
      end do
   case default
      if (param(1:6)=='ens1u_') then
         if (.not.allocated(bump%ens1u%fld)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
         read(param(7:10),'(i4.4)') ie
         if (ie>size(bump%ens1u%fld,5)) call bump%mpl%abort(subr,trim(param)//' has fewer members in bump%copy_to_field')
         iv = bump%bpar%b_to_v1(ib)
         its = bump%bpar%b_to_ts1(ib)
         call bump%geom%copy_c0a_to_mga(bump%mpl,bump%ens1u%fld(:,:,iv,its,ie),fld_mga)
      else
         call bump%mpl%abort(subr,'parameter '//trim(param)//' not yet implemented in get_parameter')
      end if
   end select
end select

end subroutine bump_copy_to_field

!----------------------------------------------------------------------
! Subroutine: bump_set_parameter
! Purpose: set a parameter
!----------------------------------------------------------------------
subroutine bump_set_parameter(bump,param,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                       ! BUMP
character(len=*),intent(in) :: param                                                         ! Parameter
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variables
integer :: ib,iv,jv,its,jts

select case (trim(param))
case ('var','cor_rh','cor_rv','cor_rv_rfac','cor_rv_coef','loc_coef','loc_rh','loc_rv','hyb_coef')
   select case (trim(bump%nam%strategy))
   case ('specific_univariate','specific_multivariate')
      do ib=1,bump%bpar%nb
         ! Get indices
         iv = bump%bpar%b_to_v1(ib)
         jv = bump%bpar%b_to_v2(ib)
         its = bump%bpar%b_to_ts1(ib)
         jts = bump%bpar%b_to_ts2(ib)

         ! Copy to field
         if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
      end do
   case ('common','common_univariate','common_weighted')
      ! Set common index
      ib = bump%bpar%nbe

      do its=1,bump%nam%nts
         do iv=1,bump%nam%nv
            ! Copy to field
            call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
         end do
      end do
   end select
case default
   do ib=1,bump%bpar%nb
      ! Get indices
      iv = bump%bpar%b_to_v1(ib)
      jv = bump%bpar%b_to_v2(ib)
      its = bump%bpar%b_to_ts1(ib)
      jts = bump%bpar%b_to_ts2(ib)

      ! Copy to field
      if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
   end do
end select

end subroutine bump_set_parameter

!----------------------------------------------------------------------
! Subroutine: bump_copy_from_field
! Purpose: copy from field
!----------------------------------------------------------------------
subroutine bump_copy_from_field(bump,param,ib,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                              ! BUMP
character(len=*),intent(in) :: param                                ! Parameter
integer,intent(in) :: ib                                            ! Block index
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field

! Local variables
integer :: imga,il0
character(len=1024),parameter :: subr = 'bump_copy_from_field'

! Allocation
if (.not.allocated(bump%cmat%blk)) allocate(bump%cmat%blk(bump%bpar%nbe))

! Select parameter
select case (trim(param))
case ('var')
   if (.not.allocated(bump%cmat%blk(ib)%bump_coef_ens)) allocate(bump%cmat%blk(ib)%bump_coef_ens(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_coef_ens)
case ('cor_rh')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rh)) allocate(bump%cmat%blk(ib)%bump_rh(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rh)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(bump%cmat%blk(ib)%bump_rh(imga,il0))) &
       & bump%cmat%blk(ib)%bump_rh(imga,il0) = bump%cmat%blk(ib)%bump_rh(imga,il0)/req
      end do
   end do
case ('cor_rv')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rv)) allocate(bump%cmat%blk(ib)%bump_rv(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rv)
case ('cor_rv_rfac')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rv_rfac)) allocate(bump%cmat%blk(ib)%bump_rv_rfac(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rv_rfac)
case ('cor_rv_coef')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rv_coef)) allocate(bump%cmat%blk(ib)%bump_rv_coef(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rv_coef)
case ('loc_coef')
   if (.not.allocated(bump%cmat%blk(ib)%bump_coef_ens)) allocate(bump%cmat%blk(ib)%bump_coef_ens(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_coef_ens)
case ('loc_rh')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rh)) allocate(bump%cmat%blk(ib)%bump_rh(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rh)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(bump%cmat%blk(ib)%bump_rh(imga,il0))) &
       & bump%cmat%blk(ib)%bump_rh(imga,il0) = bump%cmat%blk(ib)%bump_rh(imga,il0)/req
      end do
   end do
case ('loc_rv')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rv)) allocate(bump%cmat%blk(ib)%bump_rv(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rv)
case ('hyb_coef')
   if (.not.allocated(bump%cmat%blk(ib)%bump_coef_sta)) allocate(bump%cmat%blk(ib)%bump_coef_sta(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_coef_sta)
case ('D11')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D11)) allocate(bump%cmat%blk(ib)%bump_D11(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D11)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(bump%cmat%blk(ib)%bump_D11(imga,il0))) &
       & bump%cmat%blk(ib)%bump_D11(imga,il0) = bump%cmat%blk(ib)%bump_D11(imga,il0)/req**2
      end do
   end do
case ('D22')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D22)) allocate(bump%cmat%blk(ib)%bump_D22(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D22)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(bump%cmat%blk(ib)%bump_D22(imga,il0))) &
       & bump%cmat%blk(ib)%bump_D22(imga,il0) = bump%cmat%blk(ib)%bump_D22(imga,il0)/req**2
      end do
   end do
case ('D33')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D33)) allocate(bump%cmat%blk(ib)%bump_D33(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D33)
case ('D12')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D12)) allocate(bump%cmat%blk(ib)%bump_D12(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D12)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnotr(bump%cmat%blk(ib)%bump_D12(imga,il0))) &
       & bump%cmat%blk(ib)%bump_D12(imga,il0) = bump%cmat%blk(ib)%bump_D12(imga,il0)/req**2
      end do
   end do
case ('Dcoef')
   if (.not.allocated(bump%cmat%blk(ib)%bump_Dcoef)) allocate(bump%cmat%blk(ib)%bump_Dcoef(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_Dcoef)
case default
   call bump%mpl%abort(subr,'parameter '//trim(param)//' not yet implemented in set_parameter')
end select

end subroutine bump_copy_from_field

!----------------------------------------------------------------------
! Subroutine: bump_crtm_neighbors_3d
! Purpose: find nearest neighbors for CRTM, 3D
!----------------------------------------------------------------------
subroutine bump_crtm_neighbors_3d(bump,fld_mga,nobs,lon,lat,nn,nn_val,nn_dist)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                              ! BUMP
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field
integer,intent(in) :: nobs                                          ! Number of observations
real(kind_real),intent(in) :: lon(nobs)                             ! Observation longiutde
real(kind_real),intent(in) :: lat(nobs)                             ! Observation latitude
integer,intent(in) :: nn                                            ! Number of nearest neighbors
real(kind_real),intent(out) :: nn_val(bump%geom%nl0,nn,nobs)        ! Nearest neighbors values
real(kind_real),intent(out) :: nn_dist(nn,nobs)                     ! Nearest neighbors distances

! Local variables
integer :: iobs,nn_index_c0a(nn),i,nn_index_mga,il0

do iobs=1,nobs
   ! Get neighbors indices and distances
   call bump%geom%tree%find_nearest_neighbors(lon(iobs),lat(iobs),nn,nn_index_c0a,nn_dist(:,iobs))

   do i=1,nn
      ! Convert indices
      nn_index_mga = bump%geom%c0a_to_mga(nn_index_c0a(i))

      ! Get neighbors
      do il0=1,bump%geom%nl0
         nn_val(il0,i,iobs) = fld_mga(nn_index_mga,il0)
      end do
   end do
end do

end subroutine bump_crtm_neighbors_3d

!----------------------------------------------------------------------
! Subroutine: bump_crtm_neighbors_2d
! Purpose: find nearest neighbors for CRTM, 2D
!----------------------------------------------------------------------
subroutine bump_crtm_neighbors_2d(bump,fld_mga,nobs,lon,lat,nn,nn_val,nn_dist)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                ! BUMP
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga) ! Field
integer,intent(in) :: nobs                            ! Number of observations
real(kind_real),intent(in) :: lon(nobs)               ! Observation longiutde
real(kind_real),intent(in) :: lat(nobs)               ! Observation latitude
integer,intent(in) :: nn                              ! Number of nearest neighbors
real(kind_real),intent(out) :: nn_val(nn,nobs)        ! Nearest neighbors values
real(kind_real),intent(out) :: nn_dist(nn,nobs)       ! Nearest neighbors distances

! Local variables
integer :: iobs,nn_index_c0a(nn),i,nn_index_mga

do iobs=1,nobs
   ! Get neighbors indices and distances
   call bump%geom%tree%find_nearest_neighbors(lon(iobs),lat(iobs),nn,nn_index_c0a,nn_dist(:,iobs))

   do i=1,nn
      ! Convert indices
      nn_index_mga = bump%geom%c0a_to_mga(nn_index_c0a(i))

      ! Get neighbors
      nn_val(i,iobs) = fld_mga(nn_index_mga)
   end do
end do

end subroutine bump_crtm_neighbors_2d

!----------------------------------------------------------------------
! Subroutine: bump_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine bump_partial_dealloc(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

! Release memory
call bump%bpar%dealloc
call bump%cmat%dealloc
call bump%ens1%dealloc
call bump%ens1u%dealloc
call bump%ens2%dealloc
call bump%geom%dealloc
call bump%hdiag%dealloc
call bump%io%dealloc
call bump%lct%dealloc
call bump%nicas%partial_dealloc
call bump%obsop%partial_dealloc
call bump%vbal%partial_dealloc

end subroutine bump_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: bump_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine bump_dealloc(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

! Release memory
call bump%bpar%dealloc
call bump%cmat%dealloc
call bump%ens1%dealloc
call bump%ens1u%dealloc
call bump%ens2%dealloc
call bump%geom%dealloc
call bump%hdiag%dealloc
call bump%io%dealloc
call bump%lct%dealloc
call bump%nicas%dealloc
call bump%obsop%dealloc
call bump%vbal%dealloc

end subroutine bump_dealloc

end module type_bump
