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
use tools_const, only: req,deg2rad
use tools_func, only: sphere_dist
use tools_kinds,only: kind_real
use type_bpar, only: bpar_type
use type_cmat, only: cmat_type
use type_cv, only: cv_type
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
   type(io_type) :: io
   type(lct_type) :: lct
   type(mpl_type) :: mpl
   type(nam_type) :: nam
   type(nicas_type) :: nicas
   type(obsop_type) :: obsop
   type(rng_type) :: rng
   type(vbal_type) :: vbal
   logical :: close_listing
contains
   procedure :: setup_online => bump_setup_online
   procedure :: setup_generic => bump_setup_generic
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
   procedure :: apply_obsop => bump_apply_obsop
   procedure :: apply_obsop_ad => bump_apply_obsop_ad
   procedure :: get_parameter => bump_get_parameter
   procedure :: copy_to_field => bump_copy_to_field
   procedure :: set_parameter => bump_set_parameter
   procedure :: copy_from_field => bump_copy_from_field
   procedure :: dealloc => bump_dealloc
end type bump_type

private
public :: bump_type

contains

!----------------------------------------------------------------------
! Subroutine: bump_setup_online
!> Purpose: online setup
!----------------------------------------------------------------------
subroutine bump_setup_online(bump,mpi_comm,nmga,nl0,nv,nts,lon,lat,area,vunit,lmask,ens1_ne,ens1_nsub,ens2_ne,ens2_nsub, &
                           & nobs,lonobs,latobs,namelname,lunit)

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
integer,intent(in),optional :: ens1_nsub               !< Ensemble 1 number of sub-ensembles
integer,intent(in),optional :: ens2_ne                 !< Ensemble 2 size
integer,intent(in),optional :: ens2_nsub               !< Ensemble 2 size of sub-ensembles
integer,intent(in),optional :: nobs                    !< Number of observations
real(kind_real),intent(in),optional :: lonobs(:)       !< Observations longitude (in degrees: -180 to 180)
real(kind_real),intent(in),optional :: latobs(:)       !< Observations latitude (in degrees: -90 to 90)
character(len=*),intent(in),optional :: namelname      !< Namelist name
integer,intent(in),optional :: lunit                   !< Listing unit

! Local variables
integer :: lens1_ne,lens1_nsub,lens2_ne,lens2_nsub

! Initialize MPL
call bump%mpl%init(mpi_comm)

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
if (present(lunit)) then
   call bump%mpl%init_listing(bump%nam%prefix,bump%nam%model,bump%nam%colorlog,bump%nam%logpres,lunit)
   bump%close_listing = .false.
else
   call bump%mpl%init_listing(bump%nam%prefix,bump%nam%model,bump%nam%colorlog,bump%nam%logpres)
   bump%close_listing = (trim(bump%nam%model)=='online')
end if

! Generic setup
call bump%setup_generic

! Initialize geometry
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Initialize geometry'
call flush(bump%mpl%info)
call bump%geom%setup_online(bump%mpl,nmga,nl0,lon,lat,area,vunit,lmask)
call bump%geom%init(bump%mpl,bump%rng,bump%nam)

if (bump%nam%grid_output) then
   ! Initialize fields regridding
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Initialize fields regridding'
   call flush(bump%mpl%info)
   call bump%io%grid_init(bump%mpl,bump%rng,bump%nam,bump%geom)
end if

! Initialize block parameters
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Initialize block parameters'
call bump%bpar%alloc(bump%nam,bump%geom)

! Initialize ensemble 1 setup
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Initialize ensemble 1'
call bump%ens1%alloc(bump%nam,bump%geom,bump%nam%ens1_ne,bump%nam%ens1_nsub)

! Initialize ensemble 2 setup
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Initialize ensemble 2'
call bump%ens2%alloc(bump%nam,bump%geom,bump%nam%ens2_ne,bump%nam%ens2_nsub)

if (present(nobs)) then
   ! Check arguments consistency
   if ((.not.present(lonobs)).or.(.not.present(latobs))) call bump%mpl%abort('lonobs and latobs are missing')

   ! Check sizes consistency
   if (size(lonobs)/=nobs) call bump%mpl%abort('wrong size for lonobs')
   if (size(latobs)/=nobs) call bump%mpl%abort('wrong size for latobs')

   ! Initialize observations locations
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Initialize observations locations'
   call flush(bump%mpl%info)
   call bump%obsop%from(nobs,lonobs,latobs)

   ! Run drivers
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run drivers'
   call flush(bump%mpl%info)
   call bump%run_drivers

   ! Close listings
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Close listings'
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call flush(bump%mpl%info)
   close(unit=bump%mpl%info)
   call flush(bump%mpl%test)
   close(unit=bump%mpl%test)
   call bump%mpl%delete_empty_test(bump%nam%prefix)
end if

if ((bump%nam%ens1_ne>0).or.(bump%nam%ens2_ne>0)) then
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Add members to BUMP ensembles'
end if

end subroutine bump_setup_online

!----------------------------------------------------------------------
! Subroutine: bump_setup_generic
!> Purpose: generic setup
!----------------------------------------------------------------------
subroutine bump_setup_generic(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump !< BUMP

! Header
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- You are running bump ------------------------------------------'
write(bump%mpl%info,'(a)') '--- Author: Benjamin Menetrier ------------------------------------'
write(bump%mpl%info,'(a)') '--- Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE -----------'
call flush(bump%mpl%info)

! Check namelist parameters
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Check namelist parameters'
call flush(bump%mpl%info)
call bump%nam%check(bump%mpl)

! Write parallel setup
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a,i3,a,i2,a)') '--- Parallelization with ',bump%mpl%nproc,' MPI tasks and ', &
 & bump%mpl%nthread,' OpenMP threads'
call flush(bump%mpl%info)

! Initialize random number generator
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Initialize random number generator'
call flush(bump%mpl%info)
call bump%rng%init(bump%mpl,bump%nam)

! Initialize allocation flags
bump%cmat%allocated = .false.
bump%lct%allocated = .false.
bump%nicas%allocated = .false.

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

! Finalize ensemble 1
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Finalize ensemble 1'
call flush(bump%mpl%info)
call bump%ens1%remove_mean

! Finalize ensemble 2
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
write(bump%mpl%info,'(a)') '--- Finalize ensemble 2'
call flush(bump%mpl%info)
call bump%ens2%remove_mean

if (bump%nam%new_vbal) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run vertical balance driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run vertical balance driver'
   call flush(bump%mpl%info)
   call bump%vbal%run_vbal(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,bump%ens1u)
elseif (bump%nam%load_vbal) then
   ! Read vertical balance data
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Read vertical balance data'
   call flush(bump%mpl%info)
   call bump%vbal%read(bump%mpl,bump%nam,bump%geom,bump%bpar)
end if

if (bump%nam%check_vbal) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run vertical balance tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run vertical balance tests driver'
   call flush(bump%mpl%info)
   call bump%vbal%run_vbal_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar)
end if

if (bump%nam%new_hdiag) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run HDIAG driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run HDIAG driver'
   call flush(bump%mpl%info)
   if ((trim(bump%nam%method)=='hyb-rnd').or.(trim(bump%nam%method)=='dual-ens')) then
      call bump%cmat%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,bump%ens2)
   else
      call bump%cmat%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
   end if
end if

if (bump%nam%new_lct) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run LCT driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run LCT driver'
   call flush(bump%mpl%info)
   call bump%lct%run_lct(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)

   ! Copy LCT into C matrix
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Copy LCT into C matrix'
   call flush(bump%mpl%info)
   call bump%cmat%from_lct(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io,bump%lct)
end if

if (bump%nam%load_cmat) then
   ! Read C matrix data
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Read C matrix data'
   call flush(bump%mpl%info)
   call bump%cmat%read(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io)
else
   if (bump%nam%forced_radii) then
      ! Copy namelist support radii into C matrix
      write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
      write(bump%mpl%info,'(a)') '--- Copy namelist support radii into C matrix'
      call flush(bump%mpl%info)
      call bump%cmat%from_nam(bump%mpl,bump%nam,bump%geom,bump%bpar)
   end if
end if

if (allocated(bump%cmat%blk)) then
   ! Get C matrix data from OOPS
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Get C matrix data from OOPS'
   call bump%cmat%from_oops(bump%mpl,bump%geom,bump%bpar)

   ! Setup C matrix sampling
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Setup C matrix sampling'
   call bump%cmat%setup_sampling(bump%nam,bump%geom,bump%bpar)

   ! Write C matrix data
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Write C matrix data'
   call flush(bump%mpl%info)
   call bump%cmat%write(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io)
end if

if (bump%nam%new_nicas) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run NICAS driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run NICAS driver'
   call flush(bump%mpl%info)
   call bump%nicas%run_nicas(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%cmat)
elseif (bump%nam%load_nicas) then
   ! Read NICAS parameters
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Read NICAS parameters'
   call flush(bump%mpl%info)
   call bump%nicas%read(bump%mpl,bump%nam,bump%geom,bump%bpar)
end if

if (bump%nam%check_adjoints.or.bump%nam%check_pos_def.or.bump%nam%check_sqrt.or.bump%nam%check_dirac.or. &
 & bump%nam%check_randomization.or.bump%nam%check_consistency.or.bump%nam%check_optimality) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run NICAS tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run NICAS tests driver'
   call flush(bump%mpl%info)
   call bump%nicas%run_nicas_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%cmat,bump%ens1)
end if

if (bump%nam%new_obsop) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run observation operator driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run observation operator driver'
   call flush(bump%mpl%info)
   call bump%obsop%run_obsop(bump%mpl,bump%rng,bump%nam,bump%geom)
elseif (bump%nam%load_obsop) then
   ! Read observation operator
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Read observation operator'
   call flush(bump%mpl%info)
   call bump%obsop%read(bump%mpl,bump%nam)
end if

if (bump%nam%check_obsop) then
   ! Reseed random number generator
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Run observation operator tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Run observation operator tests driver'
   call flush(bump%mpl%info)
   call bump%obsop%run_obsop_tests(bump%mpl,bump%rng,bump%geom)
end if

if (bump%close_listing) then
   ! Close listings
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   write(bump%mpl%info,'(a)') '--- Close listings'
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call flush(bump%mpl%info)
   close(unit=bump%mpl%info)
   call flush(bump%mpl%test)
   close(unit=bump%mpl%test)
   call bump%mpl%delete_empty_test(bump%nam%prefix)
end if

end subroutine bump_run_drivers

!----------------------------------------------------------------------
! Subroutine: bump_add_member
!> Purpose: add member into bump%ens[1,2]
!----------------------------------------------------------------------
subroutine bump_add_member(bump,fld,ie,iens)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                      !< BUMP
real(kind_real),intent(inout) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field
integer,intent(in) :: ie                                                                    !< Member index
integer,intent(in) :: iens                                                                  !< Ensemble number

! Local variables
integer :: its,iv
real(kind_real) :: norm

! Add member
if (iens==1) then
   bump%ens1%fld(:,:,:,:,ie) = fld
elseif (iens==2) then
   bump%ens2%fld(:,:,:,:,ie) = fld
else
   call bump%mpl%abort('wrong ensemble number')
end if

! Print norm
write(bump%mpl%info,'(a7,a,i3,a,i1)') '','Member ',ie,' added to ensemble ',iens
do its=1,bump%nam%nts
   do iv=1,bump%nam%nv
      norm = sum(fld(:,:,iv,its)**2)
      write(bump%mpl%info,'(a10,a,i2,a,i2,a,e9.2)') '','Local norm for variable ',iv,' and timeslot ',its,': ',norm
   end do
end do

end subroutine bump_add_member

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal
!> Purpose: vertical balance application
!----------------------------------------------------------------------
subroutine bump_apply_vbal(bump,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                         !< BUMP
real(kind_real),intent(inout) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

! Local variable
integer :: its

! Apply vertical balance
do its=1,bump%nam%nts
   call bump%vbal%apply(bump%nam,bump%geom,bump%bpar,fld(:,:,:,its))
end do

end subroutine bump_apply_vbal

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv
!> Purpose: vertical balance application, inverse
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv(bump,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                         !< BUMP
real(kind_real),intent(inout) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

! Local variable
integer :: its

! Apply vertical balance, inverse
do its=1,bump%nam%nts
   call bump%vbal%apply_inv(bump%nam,bump%geom,bump%bpar,fld(:,:,:,its))
end do

end subroutine bump_apply_vbal_inv

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_ad
!> Purpose: vertical balance application, adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_ad(bump,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                         !< BUMP
real(kind_real),intent(inout) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

! Local variable
integer :: its

! Apply vertical balance, adjoint
do its=1,bump%nam%nts
   call bump%vbal%apply_ad(bump%nam,bump%geom,bump%bpar,fld(:,:,:,its))
end do

end subroutine bump_apply_vbal_ad

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv_ad
!> Purpose: vertical balance application, inverse adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv_ad(bump,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                         !< BUMP
real(kind_real),intent(inout) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

! Local variable
integer :: its

! Apply vertical balance, inverse adjoint
do its=1,bump%nam%nts
   call bump%vbal%apply_inv_ad(bump%nam,bump%geom,bump%bpar,fld(:,:,:,its))
end do

end subroutine bump_apply_vbal_inv_ad

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
! Subroutine: bump_get_cv_size
!> Purpose: get control variable size
!----------------------------------------------------------------------
subroutine bump_get_cv_size(bump,n)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump !< BUMP
integer,intent(out) :: n            !< Control variable size

! Local variables
type(cv_type) :: cv

! Allocate control variable
call bump%nicas%alloc_cv(bump%bpar,cv,getsizeonly=.true.)

! Copy size
n = cv%n

end subroutine bump_get_cv_size

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt
!> Purpose: NICAS square-root application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt(bump,pcv,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                         !< BUMP
real(kind_real),intent(in) :: pcv(:)                                                        !< Packed control variable
real(kind_real),intent(inout) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

! Local variables
type(cv_type) :: cv

! Allocation
call bump%nicas%alloc_cv(bump%bpar,cv)

! Check dimension
if (size(pcv)==cv%n) then
   ! Unpack control variable
   call cv%unpack(pcv)
else
   call bump%mpl%abort('wrong control variable size in bump_apply_nicas_sqrt')
end if

! Apply NICAS square-root
call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld)

end subroutine bump_apply_nicas_sqrt

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt_ad
!> Purpose: NICAS square-root adjoint application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt_ad(bump,fld,pcv)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                      !< BUMP
real(kind_real),intent(in) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field
real(kind_real),intent(inout) :: pcv(:)                                                  !< Packed control variable

! Local variables
type(cv_type) :: cv

! Apply NICAS square-root adjoint
call bump%nicas%apply_sqrt_ad(bump%mpl,bump%nam,bump%geom,bump%bpar,fld,cv)

! Check dimension
if (size(pcv)==cv%n) then
   ! Pack control variable
   call cv%pack(pcv)
else
   call bump%mpl%abort('wrong control variable size in bump_apply_nicas_sqrt_ad')
end if

end subroutine bump_apply_nicas_sqrt_ad

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
! Subroutine: bump_get_parameter
!> Purpose: get a parameter
!----------------------------------------------------------------------
subroutine bump_get_parameter(bump,param,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                                                       !< BUMP
character(len=*),intent(in) :: param                                                      !< Parameter
real(kind_real),intent(out) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

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
         if ((iv==jv).and.(its==jts)) call bump%copy_to_field(param,ib,fld(:,:,iv,its))
      end do
   case ('common','common_univariate','common_weighted')
      ! Set common index
      ib = bump%bpar%nbe

      do its=1,bump%nam%nts
         do iv=1,bump%nam%nv
            ! Copy to field
            call bump%copy_to_field(param,ib,fld(:,:,iv,its))
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
      if ((iv==jv).and.(its==jts)) call bump%copy_to_field(param,ib,fld(:,:,iv,its))
   end do
end select

end subroutine bump_get_parameter

!----------------------------------------------------------------------
! Subroutine: bump_copy_to_field
!> Purpose: copy to field
!----------------------------------------------------------------------
subroutine bump_copy_to_field(bump,param,ib,fld)

implicit none

! Passed variables
class(bump_type),intent(in) :: bump                              !< BUMP
character(len=*),intent(in) :: param                             !< Parameter
integer,intent(in) :: ib                                         !< Block index
real(kind_real),intent(out) :: fld(bump%geom%nmga,bump%geom%nl0) !< Field

! Local variables
integer :: iscales,ie,iv,its

! Select parameter
select case (trim(param))
case ('var')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_ens,fld)
case ('cor_rh')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld)
   fld = fld*req
case ('cor_rv')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld)
case ('cor_rv_rfac')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv_rfac,fld)
case ('cor_rv_coef')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv_coef,fld)
case ('loc_coef')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_ens,fld)
case ('loc_rh')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld)
   fld = fld*req
case ('loc_rv')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld)
case ('hyb_coef')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_sta,fld)
case default
   select case (param(1:4))
   case ('D11_')
      read(param(5:5),'(i1)') iscales
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D11(:,:,iscales),fld)
      fld = fld*req**2
   case ('D22_')
      read(param(5:5),'(i1)') iscales
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D22(:,:,iscales),fld)
      fld = fld*req**2
   case ('D33_')
      read(param(5:5),'(i1)') iscales
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D33(:,:,iscales),fld)
   case ('D12_')
      read(param(5:5),'(i1)') iscales
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D12(:,:,iscales),fld)
      fld = fld*req**2
   case ('Dcoe')
      read(param(7:7),'(i1)') iscales
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%Dcoef(:,:,iscales),fld)
   case ('DLh_')
      read(param(5:5),'(i1)') iscales
      call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%DLh(:,:,iscales),fld)
      fld = fld*req
   case default
      if (param(1:6)=='ens1u_') then
         read(param(7:10),'(i4.4)') ie
         iv = bump%bpar%b_to_v1(ib)
         its = bump%bpar%b_to_ts1(ib)
         call bump%geom%copy_c0a_to_mga(bump%mpl,bump%ens1u%fld(:,:,iv,its,ie),fld)
      else
         call bump%mpl%abort('parameter '//trim(param)//' not yet implemented in get_parameter')
      end if
   end select
end select

end subroutine bump_copy_to_field

!----------------------------------------------------------------------
! Subroutine: bump_set_parameter
!> Purpose: set a parameter
!----------------------------------------------------------------------
subroutine bump_set_parameter(bump,param,fld)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                   !< BUMP
character(len=*),intent(in) :: param                                                     !< Parameter
real(kind_real),intent(in) :: fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts) !< Field

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
         if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld(:,:,iv,its))
      end do
   case ('common','common_univariate','common_weighted')
      ! Set common index
      ib = bump%bpar%nbe
   
      do its=1,bump%nam%nts
         do iv=1,bump%nam%nv
            ! Copy to field
            call bump%copy_from_field(param,ib,fld(:,:,iv,its))
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
      if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld(:,:,iv,its))
   end do
end select

end subroutine bump_set_parameter

!----------------------------------------------------------------------
! Subroutine: bump_copy_from_field
!> Purpose: copy from field
!----------------------------------------------------------------------
subroutine bump_copy_from_field(bump,param,ib,fld)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                          !< BUMP
character(len=*),intent(in) :: param                            !< Parameter
integer,intent(in) :: ib                                        !< Block index
real(kind_real),intent(in) :: fld(bump%geom%nmga,bump%geom%nl0) !< Field

! Allocation
if (.not.allocated(bump%cmat%blk)) allocate(bump%cmat%blk(bump%bpar%nbe))

! Select parameter
select case (trim(param))
case ('var')
   if (.not.allocated(bump%cmat%blk(ib)%oops_coef_ens)) allocate(bump%cmat%blk(ib)%oops_coef_ens(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_coef_ens)
case ('cor_rh')
   if (.not.allocated(bump%cmat%blk(ib)%oops_rh)) allocate(bump%cmat%blk(ib)%oops_rh(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_rh)
   bump%cmat%blk(ib)%oops_rh = bump%cmat%blk(ib)%oops_rh/req
case ('cor_rv')
   if (.not.allocated(bump%cmat%blk(ib)%oops_rv)) allocate(bump%cmat%blk(ib)%oops_rv(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_rv)
case ('cor_rv_rfac')
   if (.not.allocated(bump%cmat%blk(ib)%oops_rv_rfac)) allocate(bump%cmat%blk(ib)%oops_rv_rfac(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_rv_rfac)
case ('cor_rv_coef')
   if (.not.allocated(bump%cmat%blk(ib)%oops_rv_coef)) allocate(bump%cmat%blk(ib)%oops_rv_coef(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_rv_coef)
case ('loc_coef')
   if (.not.allocated(bump%cmat%blk(ib)%oops_coef_ens)) allocate(bump%cmat%blk(ib)%oops_coef_ens(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_coef_ens)
case ('loc_rh')
   if (.not.allocated(bump%cmat%blk(ib)%oops_rh)) allocate(bump%cmat%blk(ib)%oops_rh(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_rh)
   bump%cmat%blk(ib)%oops_rh = bump%cmat%blk(ib)%oops_rh/req
case ('loc_rv')
   if (.not.allocated(bump%cmat%blk(ib)%oops_rv)) allocate(bump%cmat%blk(ib)%oops_rv(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_rv)
case ('hyb_coef')
   if (.not.allocated(bump%cmat%blk(ib)%oops_coef_sta)) allocate(bump%cmat%blk(ib)%oops_coef_sta(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld,bump%cmat%blk(ib)%oops_coef_sta)
case default
   call bump%mpl%abort('parameter '//trim(param)//' not yet implemented in set_parameter')
end select

end subroutine bump_copy_from_field

!----------------------------------------------------------------------
! Subroutine: bump_dealloc
!> Purpose: deallocation of BUMP fields
!----------------------------------------------------------------------
subroutine bump_dealloc(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump !< BUMP

! Release memory
call bump%cmat%dealloc(bump%bpar)
call bump%ens1%dealloc
call bump%ens1u%dealloc
call bump%ens2%dealloc
call bump%io%dealloc
call bump%lct%dealloc(bump%bpar)
call bump%nicas%dealloc(bump%nam,bump%geom,bump%bpar)
call bump%obsop%dealloc
call bump%vbal%dealloc(bump%nam)

! Final memory release (because objects required for previous memory releases)
call bump%bpar%dealloc
call bump%geom%dealloc

end subroutine bump_dealloc

end module type_bump
