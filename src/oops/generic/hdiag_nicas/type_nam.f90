!----------------------------------------------------------------------
! Module: type_nam
!> Purpose: namelist derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_nam

use iso_c_binding
use netcdf, only: nf90_put_att,nf90_global
use omp_lib, only: omp_get_num_procs
use tools_const, only: req
use tools_display, only: msgerror,msgwarning
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr
use tools_nc, only: ncerr,put_att
use type_mpl, only: mpl

implicit none

! Namelist parameters maximum sizes
integer,parameter :: nvmax = 20                     !< Maximum number of variables
integer,parameter :: ntsmax = 20                    !< Maximum number of time slots
integer,parameter :: nlmax = 200                    !< Maximum number of levels
integer,parameter :: nc3max = 1000                  !< Maximum number of classes
integer,parameter :: nscalesmax = 5                 !< Maximum number of variables
integer,parameter :: nldwvmax = 100                 !< Maximum number of local diagnostic profiles
integer,parameter :: ndirmax = 100                  !< Maximum number of diracs

type nam_type
   ! general_param
   character(len=1024) :: datadir                   !< Data directory
   character(len=1024) :: prefix                    !< Files prefix
   character(len=1024) :: model                     !< Model name ('aro', 'arp', 'gem', 'geos', 'gfs', 'ifs', 'mpas', 'nemo' or 'wrf')
   logical :: colorlog                              !< Add colors to the log (for display on terminal)
   logical :: default_seed                          !< Default seed for random numbers
   logical :: load_ensemble                         !< Load ensemble before computations
   logical :: use_metis                             !< Use METIS to split the domain between processors

   ! driver_param
   character(len=1024) :: method                    !< Localization/hybridization to compute ('cor', 'loc', 'hyb-avg', 'hyb-rnd' or 'dual-ens')
   character(len=1024) :: strategy                  !< Localization strategy ('common', 'specific_univariate', 'specific_multivariate' or 'common_weighted')
   logical :: new_hdiag                             !< Compute new HDIAG diagnostics (if false, read file)
   logical :: new_param                             !< Compute new NICAS parameters (if false, read file)
   logical :: check_adjoints                        !< Test adjoints
   logical :: check_pos_def                         !< Test positive definiteness
   logical :: check_sqrt                            !< Test full/square-root equivalence
   logical :: check_dirac                           !< Test NICAS application on diracs
   logical :: check_randomization                   !< Test NICAS randomization
   logical :: check_consistency                     !< Test HDIAG_NICAS consistency
   logical :: check_optimality                      !< Test HDIAG optimality
   logical :: new_lct                               !< Compute new LCT
   logical :: new_obsop                             !< Compute observation operator

   ! model_param
   integer :: nl                                    !< Number of levels
   integer :: levs(nlmax)                           !< Levels
   logical :: logpres                               !< Use pressure logarithm as vertical coordinate (model level if .false.)
   integer :: nv                                    !< Number of variables
   character(len=1024),dimension(nvmax) :: varname  !< Variables names
   character(len=1024),dimension(nvmax) :: addvar2d !< Additionnal 2d variables names
   integer :: nts                                   !< Number of time slots
   integer,dimension(ntsmax) :: timeslot            !< Timeslots

   ! ens1_param
   integer :: ens1_ne                               !< Ensemble 1 size
   integer :: ens1_ne_offset                        !< Ensemble 1 index offset
   integer :: ens1_nsub                             !< Ensemble 1 sub-ensembles number

   ! ens2_param
   integer :: ens2_ne                               !< Ensemble 2 size
   integer :: ens2_ne_offset                        !< Ensemble 2 index offset
   integer :: ens2_nsub                             !< Ensemble 2 sub-ensembles number

   ! sampling_param
   logical :: sam_write                             !< Write sampling
   logical :: sam_read                              !< Read sampling
   character(len=1024) :: mask_type                 !< Mask restriction type
   real(kind_real) ::  mask_th                      !< Mask threshold
   logical :: mask_check                            !< Check that sampling couples and interpolations do not cross mask boundaries
   integer :: nc1                                   !< Number of sampling points
   integer :: ntry                                  !< Number of tries to get the most separated point for the zero-separation sampling
   integer :: nrep                                  !< Number of replacement to improve homogeneity of the zero-separation sampling
   integer :: nc3                                   !< Number of classes
   real(kind_real) ::  dc                           !< Class size (for sam_type='hor'), should be larger than the typical grid cell size
   integer :: nl0r                                  !< Reduced number of levels for diagnostics

   ! diag_param
   integer :: ne                                    !< Ensemble sizes
   logical :: gau_approx                            !< Gaussian approximation for asymptotic quantities
   logical :: full_var                              !< Compute full variances
   logical :: local_diag                            !< Activate local diagnostics
   real(kind_real) ::  local_rad                    !< Local diagnostics calculation radius
   logical :: displ_diag                            !< Activate displacement diagnostics
   real(kind_real) ::  displ_rad                    !< Displacement diagnostics calculation radius
   integer :: displ_niter                           !< Number of iteration for the displacement filtering (for displ_diag = .true.)
   real(kind_real) ::  displ_rhflt                  !< Displacement initial filtering support radius (for displ_diag = .true.)
   real(kind_real) ::  displ_tol                    !< Displacement tolerance for mesh check (for displ_diag = .true.)

   ! fit_param
   character(len=1024) :: minim_algo                !< Minimization algorithm ('none', 'fast' or 'hooke')
   logical :: lhomh                                 !< Vertically homogenous horizontal support radius
   logical :: lhomv                                 !< Vertically homogenous vertical support radius
   real(kind_real) ::  rvflt                        !< Vertical smoother support radius
   integer :: lct_nscales                           !< Number of LCT scales
   logical :: lct_diag(nscalesmax)                  !< Diagnostic of diagonal LCT components only

   ! output_param
   integer :: nldwh                                 !< Number of local diagnostics fields to write (for local_diag = .true.)
   integer :: il_ldwh(nlmax*nc3max)                 !< Levels of local diagnostics fields to write (for local_diag = .true.)
   integer :: ic_ldwh(nlmax*nc3max)                 !< Classes of local diagnostics fields to write (for local_diag = .true.)
   integer :: nldwv                                 !< Number of local diagnostics profiles to write (for local_diag = .true.)
   real(kind_real) ::  lon_ldwv(nldwvmax)           !< Longitudes (in degrees) local diagnostics profiles to write (for local_diag = .true.)
   real(kind_real) ::  lat_ldwv(nldwvmax)           !< Latitudes (in degrees) local diagnostics profiles to write (for local_diag = .true.)
   real(kind_real) ::  diag_rhflt                   !< Diagnostics filtering radius
   character(len=1024) :: diag_interp               !< Diagnostics interpolation type

   ! nicas_param
   logical :: lsqrt                                 !< Square-root formulation
   real(kind_real) :: rh(nlmax)                     !< Default horizontal support radius
   real(kind_real) :: rv(nlmax)                     !< Default vertical support radius
   real(kind_real) :: resol                         !< Resolution
   character(len=1024) :: nicas_interp              !< NICAS interpolation type
   logical :: network                               !< Network-base convolution calculation (distance-based if false)
   integer :: mpicom                                !< Number of communication steps
   integer :: advmode                               !< Advection mode (1: direct, -1: direct and inverse)
   integer :: ndir                                  !< Number of Diracs
   real(kind_real) :: londir(ndirmax)               !< Diracs longitudes
   real(kind_real) :: latdir(ndirmax)               !< Diracs latitudes
   integer :: levdir(ndirmax)                       !< Diracs level
   integer :: ivdir(ndirmax)                        !< Diracs variable
   integer :: itsdir(ndirmax)                       !< Diracs timeslot

   ! obsop_param
   integer :: nobs                                  !< Number of observations
   real(kind_real) :: obsdis                        !< Observation distribution parameter
   character(len=1024) :: obsop_interp              !< Observation operator interpolation type
contains
   procedure :: read => nam_read
   procedure :: check => nam_check
   procedure :: ncwrite => nam_ncwrite
end type nam_type

private
public :: nam_type

contains

!----------------------------------------------------------------------
! Subroutine: nam_read
!> Purpose: read and check namelist parameters
!----------------------------------------------------------------------
subroutine nam_read(nam)

implicit none

! Passed variable
class(nam_type),intent(out) :: nam !< Namelist

! Local variables
integer :: iv

! Namelist variables
integer :: nl,levs(nlmax),nv,nts,timeslot(ntsmax),ens1_ne,ens1_ne_offset,ens1_nsub,ens2_ne,ens2_ne_offset,ens2_nsub
integer :: nc1,ntry,nrep,nc3,nl0r,ne,displ_niter,lct_nscales,nldwh,il_ldwh(nlmax*nc3max),ic_ldwh(nlmax*nc3max),nldwv
integer :: mpicom,advmode,ndir,levdir(ndirmax),ivdir(ndirmax),itsdir(ndirmax),nobs
logical :: colorlog,default_seed,load_ensemble,use_metis
logical :: new_hdiag,new_param,check_adjoints,check_pos_def,check_sqrt,check_dirac,check_randomization,check_consistency
logical :: check_optimality,new_lct,new_obsop,logpres,sam_write,sam_read,mask_check,gau_approx,full_var,local_diag
logical :: displ_diag,lhomh,lhomv,lct_diag(nscalesmax),lsqrt,network
real(kind_real) :: mask_th,dc,local_rad,displ_rad,displ_rhflt,displ_tol,rvflt,lon_ldwv(nldwvmax),lat_ldwv(nldwvmax),diag_rhflt
real(kind_real) :: rh(nlmax),rv(nlmax),resol,londir(ndirmax),latdir(ndirmax),obsdis
character(len=1024) :: datadir,prefix,model,strategy,method,mask_type,minim_algo,diag_interp,nicas_interp,obsop_interp
character(len=1024),dimension(nvmax) :: varname,addvar2d

! Namelist blocks
namelist/general_param/datadir,prefix,model,colorlog,default_seed,load_ensemble,use_metis
namelist/driver_param/method,strategy,new_hdiag,new_param,check_adjoints,check_pos_def,check_sqrt,check_dirac, &
                    & check_randomization,check_consistency,check_optimality,new_lct,new_obsop
namelist/model_param/nl,levs,logpres,nv,varname,addvar2d,nts,timeslot
namelist/ens1_param/ens1_ne,ens1_ne_offset,ens1_nsub
namelist/ens2_param/ens2_ne,ens2_ne_offset,ens2_nsub
namelist/sampling_param/sam_write,sam_read,mask_type,mask_th,mask_check,nc1,ntry,nrep,nc3,dc,nl0r
namelist/diag_param/ne,gau_approx,full_var,local_diag,local_rad,displ_diag,displ_rad,displ_niter,displ_rhflt,displ_tol
namelist/fit_param/minim_algo,lhomh,lhomv,rvflt,lct_nscales,lct_diag
namelist/output_param/nldwh,il_ldwh,ic_ldwh,nldwv,lon_ldwv,lat_ldwv,diag_rhflt,diag_interp
namelist/nicas_param/lsqrt,rh,rv,resol,nicas_interp,network,mpicom,advmode,ndir,londir,latdir,levdir,ivdir,itsdir
namelist/obsop_param/nobs,obsdis,obsop_interp

! Default initialization

! general_param default
datadir = ''
prefix = ''
model = ''
colorlog = .false.
default_seed = .false.
load_ensemble = .false.
use_metis = .false.

! driver_param default
method = ''
strategy = ''
new_hdiag = .false.
new_param = .false.
check_adjoints = .false.
check_pos_def = .false.
check_sqrt = .false.
check_dirac = .false.
check_randomization = .false.
check_consistency = .false.
check_optimality = .false.
new_lct = .false.
new_obsop = .false.

! model_param default
call msi(nl)
call msi(levs)
logpres = .false.
call msi(nv)
do iv=1,nvmax
   varname = ''
   addvar2d = ''
end do
call msi(nts)
call msi(timeslot)

! ens1_param default
call msi(ens1_ne)
call msi(ens1_ne_offset)
call msi(ens1_nsub)

! ens2_param default
call msi(ens2_ne)
call msi(ens2_ne_offset)
call msi(ens2_nsub)

! sampling_param default
sam_write = .false.
sam_read = .false.
mask_type = ''
call msr(mask_th)
mask_check = .false.
call msi(nc1)
call msi(ntry)
call msi(nrep)
call msi(nc3)
call msr(dc)
call msi(nl0r)

! diag_param default
call msi(ne)
gau_approx = .false.
full_var = .false.
local_diag = .false.
call msr(local_rad)
displ_diag = .false.
call msr(displ_rad)
call msi(displ_niter)
call msr(displ_rhflt)
call msr(displ_tol)

! fit_param default
minim_algo = ''
lhomh = .false.
lhomv = .false.
call msr(rvflt)
call msi(lct_nscales)
lct_diag = .false.

! output_param default
call msi(nldwh)
call msi(il_ldwh)
call msi(ic_ldwh)
call msi(nldwv)
call msr(lon_ldwv)
call msr(lat_ldwv)
call msr(diag_rhflt)
diag_interp = ''

! nicas_param default
lsqrt = .false.
call msr(rh)
call msr(rv)
call msr(resol)
nicas_interp = ''
network = .false.
call msi(mpicom)
call msi(advmode)
call msi(ndir)
call msr(londir)
call msr(latdir)
call msi(levdir)
call msi(ivdir)
call msi(itsdir)

! obsop_param default
call msi(nobs)
call msr(obsdis)
obsop_interp = ''

if (mpl%main) then
   ! Read namelist and copy into derived type

   ! general_param
   read(*,nml=general_param)
   nam%datadir = datadir
   nam%prefix = prefix
   nam%model = model
   nam%colorlog = colorlog
   nam%default_seed = default_seed
   nam%load_ensemble = load_ensemble
   nam%use_metis = use_metis

   ! driver_param
   read(*,nml=driver_param)
   nam%method = method
   nam%strategy = strategy
   nam%new_hdiag = new_hdiag
   nam%new_param = new_param
   nam%check_adjoints = check_adjoints
   nam%check_pos_def = check_pos_def
   nam%check_sqrt = check_sqrt
   nam%check_dirac = check_dirac
   nam%check_randomization = check_randomization
   nam%check_consistency = check_consistency
   nam%check_optimality = check_optimality
   nam%new_lct = new_lct
   nam%new_obsop = new_obsop

   ! model_param
   read(*,nml=model_param)
   nam%nl = nl
   nam%levs = levs
   nam%logpres = logpres
   nam%nv = nv
   nam%varname = varname
   nam%addvar2d = addvar2d
   nam%nts = nts
   nam%timeslot = timeslot

   ! ens1_param
   read(*,nml=ens1_param)
   nam%ens1_ne = ens1_ne
   nam%ens1_ne_offset = ens1_ne_offset
   nam%ens1_nsub = ens1_nsub

   ! ens2_param
   read(*,nml=ens2_param)
   nam%ens2_ne = ens2_ne
   nam%ens2_ne_offset = ens2_ne_offset
   nam%ens2_nsub = ens2_nsub

   ! sampling_param
   read(*,nml=sampling_param)
   nam%sam_write = sam_write
   nam%sam_read = sam_read
   nam%mask_type = mask_type
   nam%mask_th = mask_th
   nam%mask_check = mask_check
   nam%nc1 = nc1
   nam%ntry = ntry
   nam%nrep = nrep
   nam%nc3 = nc3
   nam%dc = dc/req
   nam%nl0r = nl0r

   ! diag_param
   read(*,nml=diag_param)
   nam%ne = ne
   nam%gau_approx = gau_approx
   nam%full_var = full_var
   nam%local_diag = local_diag
   nam%local_rad = local_rad/req
   nam%displ_diag = displ_diag
   nam%displ_rad = displ_rad/req
   nam%displ_niter = displ_niter
   nam%displ_rhflt = displ_rhflt/req
   nam%displ_tol = displ_tol

   ! fit_param
   read(*,nml=fit_param)
   nam%minim_algo = minim_algo
   nam%lhomh = lhomh
   nam%lhomv = lhomv
   nam%rvflt = rvflt
   nam%lct_nscales = lct_nscales
   nam%lct_diag = lct_diag

   ! output_param
   read(*,nml=output_param)
   nam%nldwh = nldwh
   nam%il_ldwh = il_ldwh
   nam%ic_ldwh = ic_ldwh
   nam%nldwv = nldwv
   nam%lon_ldwv = lon_ldwv
   nam%lat_ldwv = lat_ldwv
   nam%diag_rhflt = diag_rhflt/req
   nam%diag_interp = diag_interp

   ! nicas_param
   read(*,nml=nicas_param)
   nam%lsqrt = lsqrt
   nam%rh = rh/req
   nam%rv = rv
   nam%resol = resol
   nam%nicas_interp = nicas_interp
   nam%network = network
   nam%mpicom = mpicom
   nam%advmode = advmode
   nam%ndir = ndir
   nam%londir = londir
   nam%latdir = latdir
   nam%levdir = levdir
   nam%ivdir = ivdir
   nam%itsdir = itsdir

   ! obsop_param
   read(*,nml=obsop_param)
   nam%nobs = nobs
   nam%obsdis = obsdis
   nam%obsop_interp = obsop_interp
end if

! Broadcast parameters

! general_param
call mpl%bcast(nam%datadir,mpl%ioproc)
call mpl%bcast(nam%prefix,mpl%ioproc)
call mpl%bcast(nam%model,mpl%ioproc)
call mpl%bcast(nam%colorlog,mpl%ioproc)
call mpl%bcast(nam%default_seed,mpl%ioproc)
call mpl%bcast(nam%load_ensemble,mpl%ioproc)
call mpl%bcast(nam%use_metis,mpl%ioproc)

! driver_param
call mpl%bcast(nam%method,mpl%ioproc)
call mpl%bcast(nam%strategy,mpl%ioproc)
call mpl%bcast(nam%new_hdiag,mpl%ioproc)
call mpl%bcast(nam%new_param,mpl%ioproc)
call mpl%bcast(nam%check_adjoints,mpl%ioproc)
call mpl%bcast(nam%check_pos_def,mpl%ioproc)
call mpl%bcast(nam%check_sqrt,mpl%ioproc)
call mpl%bcast(nam%check_dirac,mpl%ioproc)
call mpl%bcast(nam%check_randomization,mpl%ioproc)
call mpl%bcast(nam%check_consistency,mpl%ioproc)
call mpl%bcast(nam%check_optimality,mpl%ioproc)
call mpl%bcast(nam%new_lct,mpl%ioproc)
call mpl%bcast(nam%new_obsop,mpl%ioproc)

! model_param
call mpl%bcast(nam%nl,mpl%ioproc)
call mpl%bcast(nam%levs,mpl%ioproc)
call mpl%bcast(nam%logpres,mpl%ioproc)
call mpl%bcast(nam%nv,mpl%ioproc)
call mpl%bcast(nam%varname,mpl%ioproc)
call mpl%bcast(nam%addvar2d,mpl%ioproc)
call mpl%bcast(nam%nts,mpl%ioproc)
call mpl%bcast(nam%timeslot,mpl%ioproc)

! ens1_param
call mpl%bcast(nam%ens1_ne,mpl%ioproc)
call mpl%bcast(nam%ens1_ne_offset,mpl%ioproc)
call mpl%bcast(nam%ens1_nsub,mpl%ioproc)

! ens2_param
call mpl%bcast(nam%ens2_ne,mpl%ioproc)
call mpl%bcast(nam%ens2_ne_offset,mpl%ioproc)
call mpl%bcast(nam%ens2_nsub,mpl%ioproc)

! sampling_param
call mpl%bcast(nam%sam_write,mpl%ioproc)
call mpl%bcast(nam%sam_read,mpl%ioproc)
call mpl%bcast(nam%mask_type,mpl%ioproc)
call mpl%bcast(nam%mask_th,mpl%ioproc)
call mpl%bcast(nam%mask_check,mpl%ioproc)
call mpl%bcast(nam%nc1,mpl%ioproc)
call mpl%bcast(nam%ntry,mpl%ioproc)
call mpl%bcast(nam%nrep,mpl%ioproc)
call mpl%bcast(nam%nc3,mpl%ioproc)
call mpl%bcast(nam%dc,mpl%ioproc)
call mpl%bcast(nam%nl0r,mpl%ioproc)

! diag_param
call mpl%bcast(nam%ne,mpl%ioproc)
call mpl%bcast(nam%gau_approx,mpl%ioproc)
call mpl%bcast(nam%full_var,mpl%ioproc)
call mpl%bcast(nam%local_diag,mpl%ioproc)
call mpl%bcast(nam%local_rad,mpl%ioproc)
call mpl%bcast(nam%displ_diag,mpl%ioproc)
call mpl%bcast(nam%displ_rad,mpl%ioproc)
call mpl%bcast(nam%displ_niter,mpl%ioproc)
call mpl%bcast(nam%displ_rhflt,mpl%ioproc)
call mpl%bcast(nam%displ_tol,mpl%ioproc)

! fit_param
call mpl%bcast(nam%minim_algo,mpl%ioproc)
call mpl%bcast(nam%lhomh,mpl%ioproc)
call mpl%bcast(nam%lhomv,mpl%ioproc)
call mpl%bcast(nam%rvflt,mpl%ioproc)
call mpl%bcast(nam%lct_nscales,mpl%ioproc)
call mpl%bcast(nam%lct_diag,mpl%ioproc)

! output_param
call mpl%bcast(nam%nldwh,mpl%ioproc)
call mpl%bcast(nam%il_ldwh,mpl%ioproc)
call mpl%bcast(nam%ic_ldwh,mpl%ioproc)
call mpl%bcast(nam%nldwv,mpl%ioproc)
call mpl%bcast(nam%lon_ldwv,mpl%ioproc)
call mpl%bcast(nam%lat_ldwv,mpl%ioproc)
call mpl%bcast(nam%diag_rhflt,mpl%ioproc)
call mpl%bcast(nam%diag_interp,mpl%ioproc)

! nicas_param
call mpl%bcast(nam%lsqrt,mpl%ioproc)
call mpl%bcast(nam%rh,mpl%ioproc)
call mpl%bcast(nam%rv,mpl%ioproc)
call mpl%bcast(nam%resol,mpl%ioproc)
call mpl%bcast(nam%nicas_interp,mpl%ioproc)
call mpl%bcast(nam%network,mpl%ioproc)
call mpl%bcast(nam%mpicom,mpl%ioproc)
call mpl%bcast(nam%advmode,mpl%ioproc)
call mpl%bcast(nam%ndir,mpl%ioproc)
call mpl%bcast(nam%londir,mpl%ioproc)
call mpl%bcast(nam%latdir,mpl%ioproc)
call mpl%bcast(nam%levdir,mpl%ioproc)
call mpl%bcast(nam%ivdir,mpl%ioproc)
call mpl%bcast(nam%itsdir,mpl%ioproc)

! obsop_param
call mpl%bcast(nam%nobs,mpl%ioproc)
call mpl%bcast(nam%obsdis,mpl%ioproc)
call mpl%bcast(nam%obsop_interp,mpl%ioproc)

end subroutine nam_read

!----------------------------------------------------------------------
! Subroutine: nam_check
!> Purpose: check namelist parameters
!----------------------------------------------------------------------
subroutine nam_check(nam)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam !< Namelist

! Local variables
integer :: iv,its,il,idir,ildw,itest
character(len=2) :: ivchar
character(len=4) :: itestchar
character(len=7) :: lonchar,latchar
character(len=1024) :: filename

! Check general_param
if (trim(nam%datadir)=='') call msgerror('datadir not specified')
if (trim(nam%prefix)=='') call msgerror('prefix not specified')
select case (trim(nam%model))
case ('aro','arp','gem','geos','gfs','ifs','mpas','nemo','oops','wrf')
case default
   call msgerror('wrong model')
end select

! Check driver_param
select case (trim(nam%method))
case ('cor','loc','hyb-avg','hyb-rnd','dual-ens')
case default
   call msgerror('wrong method')
end select
select case (trim(nam%strategy))
case ('common','specific_univariate','common_weighted')
case ('specific_multivariate')
   if (.not.nam%lsqrt) call msgerror('specific multivariate strategy requires a square-root formulation')
case default
   call msgerror('wrong strategy')
end select
if (nam%check_sqrt.and.(.not.nam%new_param)) call msgerror('square-root check requires new parameters calculation')
if (nam%check_randomization) then
   if (.not.nam%lsqrt) call msgerror('lsqrt required for check_randomization')
end if
if (nam%check_consistency) then
   if (.not.nam%new_hdiag) call msgerror('new_hdiag required for check_consistency')
   if (.not.nam%new_param) call msgerror('new_param required for check_consistency')
   if (.not.nam%lsqrt) call msgerror('lsqrt required for check_consistency')
end if
if (nam%check_optimality) then
   if (.not.nam%new_hdiag) call msgerror('new_hdiag required for check_optimality')
   if (.not.nam%new_param) call msgerror('new_param required for check_optimality')
   if (.not.nam%lsqrt) call msgerror('lsqrt required for check_optimality')
end if
if (nam%new_lct) then
   if (nam%new_hdiag.or.nam%new_param.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_sqrt.or.nam%check_dirac.or. &
 & nam%check_randomization.or.nam%check_consistency.or.nam%check_optimality) call msgerror('new_lct should be executed alone')
   if (.not.nam%local_diag) then
      call msgwarning('new_lct requires local_diag, resetting local_diag to .true.')
      nam%local_diag = .true.
   end if
   if (nam%displ_diag) call msgerror('new_lct requires displ_diag deactivated')
end if

! Check model_param
if (nam%nl<=0) call msgerror('nl should be positive')
do il=1,nam%nl
   if (nam%levs(il)<=0) call msgerror('levs should be positive')
   if (count(nam%levs(1:nam%nl)==nam%levs(il))>1) call msgerror('redundant levels')
end do
if (nam%logpres) then
   select case (trim(nam%model))
   case ('nemo')
      call msgwarning('pressure logarithm vertical coordinate is not available for this model, resetting to model level index')
      nam%logpres = .false.
   end select
end if
if (nam%nv<=0) call msgerror('nv should be positive')
do iv=1,nam%nv
   write(ivchar,'(i2.2)') iv
   if (trim(nam%varname(iv))=='') call msgerror('varname not specified for variable '//ivchar)
end do
do its=1,nam%nts
   if (nam%timeslot(its)<0) call msgerror('timeslot should be non-negative')
end do

! Surface level
do iv=1,nam%nv
   if (trim(nam%addvar2d(iv))/='') nam%levs(nam%nl+1) = maxval(nam%levs(1:nam%nl))+1
end do

if (nam%new_hdiag.or.nam%new_lct) then
   ! Check ens1_param
   if (nam%ens1_ne_offset<0) call msgerror('ens1_ne_offset should be non-negative')
   if (nam%ens1_nsub<1) call msgerror('ens1_nsub should be positive')
   if (mod(nam%ens1_ne,nam%ens1_nsub)/=0) call msgerror('ens1_nsub should be a divider of ens1_ne')
   if (nam%ens1_ne/nam%ens1_nsub<=3) call msgerror('ens1_ne/ens1_nsub should be larger than 3')

   ! Check ens2_param
   select case (trim(nam%method))
   case ('hyb-rnd','dual-ens')
      if (nam%ens2_ne_offset<0) call msgerror('ens2_ne_offset should be non-negative')
      if (nam%ens2_nsub<1) call msgerror('ens2_nsub should be non-negative')
      if (mod(nam%ens2_ne,nam%ens2_nsub)/=0) call msgerror('ens2_nsub should be a divider of ens2_ne')
      if (nam%ens2_ne/nam%ens2_nsub<=3) call msgerror('ens2_ne/ens2_nsub should be larger than 3')
   end select

   ! Check sampling_param
   if (nam%sam_write.and.nam%sam_read) call msgerror('sam_write and sam_read are both true')
   if (nam%nc1<=0) call msgerror('nc1 should be positive')
   if (nam%ntry<=0) call msgerror('ntry should be positive')
   if (nam%nrep<0) call msgerror('nrep should be non-negative')
   if (nam%nc3<=0) call msgerror('nc3 should be positive')
   if (nam%dc<0.0) call msgerror('dc should be positive')
   if (nam%nl0r<1) call msgerror ('nl0r should be positive')
   if (any(nam%addvar2d(1:nam%nv)/='')) then
      if (nam%nl0r>nam%nl+1) then
         call msgwarning('nl0r should be lower that nl+1, resetting nl0r to nl+1 or the lower odd number')
         nam%nl0r = nam%nl+1
         if (mod(nam%nl0r,2)<1) nam%nl0r = nam%nl0r-1
      end if
   else
      if (nam%nl0r>nam%nl) then
         call msgwarning('nl0r should be lower that nl, resetting nl0r to nl or the lower odd number')
         nam%nl0r = nam%nl
         if (mod(nam%nl0r,2)<1) nam%nl0r = nam%nl0r-1
      end if
   end if
   if (mod(nam%nl0r,2)<1) call msgerror ('nl0r should be odd')

   ! Check diag_param
   if (nam%ne<=3) call msgerror('ne should be larger than 3')
   if (nam%local_diag.or.nam%displ_diag) then
      if (nam%displ_rad<0.0) call msgerror('displ_rad should be non-negative')
   end if
   if (nam%displ_diag) then
      if (nam%local_rad<0.0) call msgerror('local_rad should be non-negative')
      if (nam%displ_niter<0) call msgerror('displ_niter should be positive')
      if (nam%displ_rhflt<0.0) call msgerror('displ_rhflt should be non-negative')
      if (nam%displ_tol<0.0) call msgerror('displ_tol should be non-negative')
   end if

   ! Check fit_param
   select case (trim(nam%minim_algo))
   case ('none','fast','hooke')
   case default
      call msgerror('wrong minim_algo')
   end select
   if (nam%new_lct.and.((trim(nam%minim_algo)=='none').or.(trim(nam%minim_algo)=='fast'))) call msgerror('wrong minim_algo for LCT')
   if (nam%rvflt<0) call msgerror('rvflt should be non-negative')
   if (nam%new_lct) then
      if (nam%lct_nscales<0) call msgerror('lct_nscales should be non-negative')
   end if

   ! Check output_param
   if (nam%local_diag) then
      if (nam%nldwh<0) call msgerror('nldwh should be non-negative')
      if (any(nam%il_ldwh(1:nam%nldwh)<0)) call msgerror('il_ldwh should be non-negative')
      if (any(nam%il_ldwh(1:nam%nldwh)>nam%nl)) call msgerror('il_ldwh should be lower than nl')
      if (any(nam%ic_ldwh(1:nam%nldwh)<0)) call msgerror('ic_ldwh should be non-negative')
      if (any(nam%ic_ldwh(1:nam%nldwh)>nam%nc3)) call msgerror('ic_ldwh should be lower than nc3')
      if (nam%nldwv<0) call msgerror('nldwv should be non-negative')
      if (any(nam%lon_ldwv(1:nam%nldwv)<-180.0).or.any(nam%lon_ldwv(1:nam%nldwv)>180.0)) call msgerror('wrong lon_ldwv')
      if (any(nam%lat_ldwv(1:nam%nldwv)<-90.0).or.any(nam%lat_ldwv(1:nam%nldwv)>90.0)) call msgerror('wrong lat_ldwv')
   end if
   if (nam%local_diag.or.nam%displ_diag) then
      if (nam%diag_rhflt<0.0) call msgerror('diag_rhflt should be non-negative')
   end if
   select case (trim(nam%diag_interp))
   case ('bilin','natural')
   case default
      call msgerror('wrong interpolation for diagnostics')
   end select

   ! Check ensemble sizes
   if (trim(nam%method)/='cor') then
      if (nam%ne>nam%ens1_ne) call msgwarning('ensemble size larger than ens1_ne (might enhance sampling noise)')
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd','dual-ens')
         if (nam%ne>nam%ens2_ne) call msgwarning('ensemble size larger than ens2_ne (might enhance sampling noise)')
      end select
   end if
end if

! Check nicas_param
if (nam%lsqrt) then
   if (nam%mpicom==1) call msgerror('mpicom should be 2 for square-root application')
end if
if (nam%new_param) then
   if (.not.(nam%resol>0.0)) call msgerror('resol should be positive')
end if
if (nam%new_param.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_sqrt.or.nam%check_dirac.or. &
 & nam%check_randomization.or.nam%check_consistency.or.nam%check_optimality) then
   if ((nam%mpicom/=1).and.(nam%mpicom/=2)) call msgerror('mpicom should be 1 or 2')
end if
if (abs(nam%advmode)>1) call msgerror('nam%advmode should be -1, 0 or 1')
if (nam%check_dirac) then
   if (nam%ndir<1) call msgerror('ndir should be positive')
   do idir=1,nam%ndir
      if ((nam%londir(idir)<-180.0).or.(nam%londir(idir)>180.0)) call msgerror('Dirac longitude should lie between -180 and 180')
      if ((nam%latdir(idir)<-90.0).or.(nam%latdir(idir)>90.0)) call msgerror('Dirac latitude should lie between -90 and 90')
      if (.not.(any(nam%levs(1:nam%nl)==nam%levdir(idir)).or.(any(nam%addvar2d(1:nam%nv)/='') &
    & .and.(nam%levs(nam%nl+1)==nam%levdir(idir))))) call msgerror('wrong level for a Dirac')
      if ((nam%ivdir(idir)<1).or.(nam%ivdir(idir)>nam%nv)) call msgerror('wrong variable for a Dirac')
      if ((nam%itsdir(idir)<1).or.(nam%itsdir(idir)>nam%nts)) call msgerror('wrong timeslot for a Dirac')
   end do
end if
select case (trim(nam%diag_interp))
case ('bilin','natural')
case default
   call msgerror('wrong interpolation for NICAS')
end select

! Check obsop_param
if (nam%new_obsop) then
   if (nam%nobs<1) call msgerror('nobs should be positive')
   if (nam%obsdis>1.0) call msgerror('obsdis should be lower than 1')
   select case (trim(nam%diag_interp))
   case ('bilin','natural')
   case default
      call msgerror('wrong interpolation for observation operator')
   end select
end if

! Clean files
if (mpl%main) then
   ! Diagnostics
   if (nam%new_hdiag) then
      filename = trim(nam%prefix)//'_diag.nc'
      call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))

      if (nam%local_diag) then
         filename = trim(nam%prefix)//'_local_diag_*.nc'
         call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      end if

      do ildw=1,nam%nldwv
         write(lonchar,'(f7.2)') nam%lon_ldwv(ildw)
         write(latchar,'(f7.2)') nam%lat_ldwv(ildw)
         filename = trim(nam%prefix)//'_diag_'//trim(adjustl(lonchar))//'-'//trim(adjustl(latchar))//'.nc'
         call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      end do
   end if

   ! Dirac test
   if (nam%check_dirac) then
      filename = trim(nam%prefix)//'_dirac.nc'
      call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      filename = trim(nam%prefix)//'_dirac_gridded.nc'
      call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
   end if

   ! Randomization test
   if (nam%check_randomization) then
      do itest=1,10
         write(itestchar,'(i4.4)') itest
         filename = trim(nam%prefix)//'_randomize_'//itestchar//'_gridded.nc'
         call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      end do
   end if

   ! LCT
   if (nam%new_lct) then
      filename = trim(nam%prefix)//'_lct_gridded.nc'
      call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
   end if
end if

end subroutine nam_check

!----------------------------------------------------------------------
! Subroutine: nam_ncwrite
!> Purpose: write namelist parameters as NetCDF attributes
!----------------------------------------------------------------------
subroutine nam_ncwrite(nam,ncid)

implicit none

! Passed variable
class(nam_type),intent(in) :: nam !< Namelist
integer,intent(in) :: ncid        !< NetCDF file ID

! general_param
call put_att(ncid,'datadir',trim(nam%datadir))
call put_att(ncid,'prefix',trim(nam%prefix))
call put_att(ncid,'model',trim(nam%model))
call put_att(ncid,'colorlog',nam%colorlog)
call put_att(ncid,'default_seed',nam%default_seed)
call put_att(ncid,'load_ensemble',nam%load_ensemble)
call put_att(ncid,'use_metis',nam%use_metis)

! driver_param
call put_att(ncid,'method',trim(nam%method))
call put_att(ncid,'strategy',trim(nam%strategy))
call put_att(ncid,'new_hdiag',nam%new_hdiag)
call put_att(ncid,'new_param',nam%new_param)
call put_att(ncid,'check_adjoints',nam%check_adjoints)
call put_att(ncid,'check_pos_def',nam%check_pos_def)
call put_att(ncid,'check_sqrt',nam%check_sqrt)
call put_att(ncid,'check_dirac',nam%check_dirac)
call put_att(ncid,'check_randomization',nam%check_randomization)
call put_att(ncid,'check_consistency',nam%check_consistency)
call put_att(ncid,'check_optimality',nam%check_optimality)
call put_att(ncid,'new_lct',nam%new_lct)

! model_param
call put_att(ncid,'nl',nam%nl)
call put_att(ncid,'levs',nam%nl,nam%levs(1:nam%nl))
call put_att(ncid,'logpres',nam%logpres)
call put_att(ncid,'nv',nam%nv)
call put_att(ncid,'varname',nam%nv,nam%varname(1:nam%nv))
call put_att(ncid,'addvar2d',nam%nv,nam%addvar2d(1:nam%nv))
call put_att(ncid,'nts',nam%nts)
call put_att(ncid,'timeslot',nam%nts,nam%timeslot(1:nam%nts))

! ens1_param
call put_att(ncid,'ens1_ne',nam%ens1_ne)
call put_att(ncid,'ens1_ne_offset',nam%ens1_ne_offset)
call put_att(ncid,'ens1_nsub',nam%ens1_nsub)

! ens2_param
call put_att(ncid,'ens2_ne',nam%ens2_ne)
call put_att(ncid,'ens2_ne_offset',nam%ens2_ne_offset)
call put_att(ncid,'ens2_nsub',nam%ens2_nsub)

! sampling_param
call put_att(ncid,'sam_write',nam%sam_write)
call put_att(ncid,'sam_read',nam%sam_read)
call put_att(ncid,'mask_type',nam%mask_type)
call put_att(ncid,'mask_th',nam%mask_th)
call put_att(ncid,'mask_check',nam%mask_check)
call put_att(ncid,'nc1',nam%nc1)
call put_att(ncid,'ntry',nam%ntry)
call put_att(ncid,'nrep',nam%nrep)
call put_att(ncid,'nc3',nam%nc3)
call put_att(ncid,'dc',nam%dc)
call put_att(ncid,'nl0r',nam%nl0r)

! diag_param
call put_att(ncid,'ne',nam%ne)
call put_att(ncid,'gau_approx',nam%gau_approx)
call put_att(ncid,'full_var',nam%full_var)
call put_att(ncid,'local_diag',nam%local_diag)
call put_att(ncid,'local_rad',nam%local_rad)
call put_att(ncid,'displ_diag',nam%displ_diag)
call put_att(ncid,'displ_rad',nam%displ_rad)
call put_att(ncid,'displ_niter',nam%displ_niter)
call put_att(ncid,'displ_rhflt',nam%displ_rhflt)
call put_att(ncid,'displ_tol',nam%displ_tol)

! fit_param
call put_att(ncid,'minim_algo',nam%minim_algo)
call put_att(ncid,'lhomh',nam%lhomh)
call put_att(ncid,'lhomv',nam%lhomv)
call put_att(ncid,'rvflt',nam%rvflt)
call put_att(ncid,'lct_nscales',nam%lct_nscales)
call put_att(ncid,'lct_diag',nam%lct_nscales,nam%lct_diag)

! output_param
call put_att(ncid,'nldwh',nam%nldwh)
call put_att(ncid,'il_ldwh',nam%nldwh,nam%il_ldwh(1:nam%nldwh))
call put_att(ncid,'ic_ldwh',nam%nldwh,nam%ic_ldwh(1:nam%nldwh))
call put_att(ncid,'nldwv',nam%nldwv)
call put_att(ncid,'lon_ldwv',nam%nldwv,nam%lon_ldwv(1:nam%nldwv))
call put_att(ncid,'lat_ldwv',nam%nldwv,nam%lat_ldwv(1:nam%nldwv))
call put_att(ncid,'diag_rhflt',nam%diag_rhflt)
call put_att(ncid,'diag_interp',nam%diag_interp)

! nicas_param
call put_att(ncid,'lsqrt',nam%lsqrt)
call put_att(ncid,'rh',nam%nl,nam%rh(1:nam%nl))
call put_att(ncid,'rv',nam%nl,nam%rv(1:nam%nl))
call put_att(ncid,'resol',nam%resol)
call put_att(ncid,'nicas_interp',nam%nicas_interp)
call put_att(ncid,'network',nam%network)
call put_att(ncid,'mpicom',nam%mpicom)
call put_att(ncid,'advmode',nam%advmode)
call put_att(ncid,'ndir',nam%ndir)
call put_att(ncid,'londir',nam%ndir,nam%londir(1:nam%ndir))
call put_att(ncid,'latdir',nam%ndir,nam%latdir(1:nam%ndir))
call put_att(ncid,'levdir',nam%ndir,nam%levdir(1:nam%ndir))
call put_att(ncid,'ivdir',nam%ndir,nam%ivdir(1:nam%ndir))
call put_att(ncid,'itsdir',nam%ndir,nam%itsdir(1:nam%ndir))

! obsop_param
call put_att(ncid,'nobs',nam%nobs)
call put_att(ncid,'obsdis',nam%obsdis)
call put_att(ncid,'obsop_interp',nam%obsop_interp)

end subroutine nam_ncwrite

end module type_nam
