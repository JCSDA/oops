!----------------------------------------------------------------------
! Module: type_nam
! Purpose: namelist derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_nam

use iso_c_binding
use netcdf, only: nf90_put_att,nf90_global
!$ use omp_lib, only: omp_get_num_procs
use tools_const, only: pi,req,deg2rad,rad2deg
use tools_kinds,only: kind_real
use type_mpl, only: mpl_type

implicit none

integer,parameter :: nvmax = 20     ! Maximum number of variables
integer,parameter :: ntsmax = 99    ! Maximum number of time slots
integer,parameter :: nlmax = 200    ! Maximum number of levels
integer,parameter :: nc3max = 1000  ! Maximum number of classes
integer,parameter :: nscalesmax = 5 ! Maximum number of variables
integer,parameter :: ndirmax = 300  ! Maximum number of diracs
integer,parameter :: nldwvmax = 99  ! Maximum number of local diagnostic profiles

type nam_type
   ! general_param
   character(len=1024) :: datadir                       ! Data directory
   character(len=1024) :: prefix                        ! Files prefix
   character(len=1024) :: model                         ! Model name ('aro', 'arp', 'fv3', 'gem', 'geos', 'gfs', 'ifs', 'mpas', 'nemo', 'qg, 'res' or 'wrf')
   character(len=1024) :: verbosity                     ! Verbosity level ('all', 'main' or 'none')
   logical :: colorlog                                  ! Add colors to the log (for display on terminal)
   logical :: default_seed                              ! Default seed for random numbers

   ! driver_param
   character(len=1024) :: method                        ! Localization/hybridization to compute ('cor', 'loc', 'hyb-avg', 'hyb-rnd' or 'dual-ens')
   character(len=1024) :: strategy                      ! Localization strategy ('diag_all', 'common', 'common_univariate', 'common_weighted', 'specific_univariate' or 'specific_multivariate')
   logical :: new_cortrack                              ! New correlation tracker
   logical :: new_vbal                                  ! Compute new vertical balance operator
   logical :: load_vbal                                 ! Load existing vertical balance operator
   logical :: write_vbal                                ! Write vertical balance operator
   logical :: new_mom                                   ! Compute new sample moments
   logical :: load_mom                                  ! Load sample moments
   logical :: write_mom                                 ! Write sample moments
   logical :: new_hdiag                                 ! Compute new HDIAG diagnostics
   logical :: write_hdiag                               ! Write HDIAG diagnostics
   logical :: new_lct                                   ! Compute new LCT
   logical :: write_lct                                 ! Write LCT
   logical :: load_cmat                                 ! Load existing C matrix
   logical :: write_cmat                                ! Write existing C matrix
   logical :: new_nicas                                 ! Compute new NICAS parameters
   logical :: load_nicas                                ! Load existing NICAS parameters
   logical :: write_nicas                               ! Write NICAS parameters
   logical :: new_obsop                                 ! Compute new observation operator
   logical :: load_obsop                                ! Load existing observation operator
   logical :: write_obsop                               ! Write observation operator
   logical :: check_vbal                                ! Test vertical balance inverse and adjoint
   logical :: check_adjoints                            ! Test NICAS adjoints
   logical :: check_pos_def                             ! Test NICAS positive definiteness
   logical :: check_dirac                               ! Test NICAS application on diracs
   logical :: check_randomization                       ! Test NICAS randomization
   logical :: check_consistency                         ! Test HDIAG-NICAS consistency
   logical :: check_optimality                          ! Test HDIAG optimality
   logical :: check_obsop                               ! Test observation operator

   ! model_param
   integer :: nl                                        ! Number of levels
   integer :: levs(nlmax)                               ! Levels
   logical :: logpres                                   ! Use pressure logarithm as vertical coordinate (model level if .false.)
   integer :: nv                                        ! Number of variables
   character(len=1024),dimension(nvmax) :: varname      ! Variables names
   character(len=1024),dimension(nvmax) :: addvar2d     ! Additionnal 2d variables names
   integer :: nts                                       ! Number of time slots
   integer,dimension(ntsmax) :: timeslot                ! Timeslots
   logical :: nomask                                    ! Do not use geometry mask

   ! ens1_param
   integer :: ens1_ne                                   ! Ensemble 1 size
   integer :: ens1_nsub                                 ! Ensemble 1 sub-ensembles number

   ! ens2_param
   integer :: ens2_ne                                   ! Ensemble 2 size
   integer :: ens2_nsub                                 ! Ensemble 2 sub-ensembles number

   ! sampling_param
   logical :: sam_write                                 ! Write sampling
   logical :: sam_read                                  ! Read sampling
   character(len=1024) :: mask_type                     ! Mask restriction type
   character(len=1024) :: mask_lu                       ! Mask threshold side ("lower" if mask_th is the lower bound, "upper" if mask_th is the upper bound)
   real(kind_real) ::  mask_th                          ! Mask threshold
   integer :: ncontig_th                                ! Threshold on vertically contiguous points for sampling mask (0 to skip the test)
   logical :: mask_check                                ! Check that sampling couples and interpolations do not cross mask boundaries
   character(len=1024) :: draw_type                     ! Sampling draw type ('random_uniform','random_coast' or 'icosahedron')
   real(kind_real) ::  Lcoast                           ! Length-scale to increase sampling density along coasts
   real(kind_real) ::  rcoast                           ! Minimum value to increase sampling density along coasts
   integer :: nc1                                       ! Number of sampling points
   integer :: nc2                                       ! Number of diagnostic points
   integer :: ntry                                      ! Number of tries to get the most separated point for the zero-separation sampling
   integer :: nrep                                      ! Number of replacement to improve homogeneity of the zero-separation sampling
   integer :: nc3                                       ! Number of classes
   real(kind_real) ::  dc                               ! Class size (for sam_type='hor'), should be larger than the typical grid cell size
   integer :: nl0r                                      ! Reduced number of levels for diagnostics
   integer :: irmax                                     ! Maximum number of random number draws

   ! diag_param
   integer :: ne                                        ! Ensemble size
   logical :: gau_approx                                ! Gaussian approximation for asymptotic quantities
   integer :: avg_nbins                                 ! Number of bins for averaged statistics histograms
   logical :: vbal_block(nvmax*(nvmax-1)/2)             ! Activation of vertical balance (ordered line by line in the lower triangular formulation)
   real(kind_real) :: vbal_rad                          ! Vertical balance diagnostic radius
   logical :: var_filter                                ! Filter variances
   integer :: var_niter                                 ! Number of iteration for the variances filtering (for var_filter = .true.)
   real(kind_real) ::  var_rhflt                        ! Variances initial filtering support radius (for var_filter = .true.)
   logical :: local_diag                                ! Activate local diagnostics
   real(kind_real) ::  local_rad                        ! Local diagnostics calculation radius (for local_rad = .true.)
   logical :: adv_diag                                  ! Activate advection diagnostic
   real(kind_real) ::  adv_rad                          ! Advection diagnostic calculation radius
   integer :: adv_niter                                 ! Number of iteration for the advection filtering (for adv_diag = .true.)
   real(kind_real) ::  adv_rhflt                        ! Advection initial filtering support radius (for adv_diag = .true.)

   ! fit_param
   character(len=1024) :: minim_algo                    ! Minimization algorithm ('none', 'fast' or 'hooke')
   logical :: double_fit(0:nvmax)                       ! Double fit to introduce negative lobes on the vertical
   logical :: lhomh                                     ! Vertically homogenous horizontal support radius
   logical :: lhomv                                     ! Vertically homogenous vertical support radius
   real(kind_real) ::  rvflt                            ! Vertical smoother support radius
   integer :: lct_nscales                               ! Number of LCT scales
   logical :: lct_diag(nscalesmax)                      ! Diagnostic of diagonal LCT components only

   ! nicas_param
   logical :: nonunit_diag                              ! Non-unit diagonal for the NICAS application
   logical :: lsqrt                                     ! Square-root formulation
   real(kind_real) :: resol                             ! Resolution
   logical :: fast_sampling                             ! Fast sampling flag
   character(len=1024) :: subsamp                       ! Subsampling structure ('h', 'hv', 'vh' or 'hvh')
   character(len=1024) :: nicas_interp                  ! NICAS interpolation type
   logical :: network                                   ! Network-base convolution calculation (distance-based if false)
   integer :: mpicom                                    ! Number of communication steps
   integer :: adv_mode                                  ! Advection mode (1: direct, -1: direct and inverse)
   logical :: forced_radii                              ! Force specific support radii
   real(kind_real) :: rh                                ! Forced horizontal support radius
   real(kind_real) :: rv                                ! Forced vertical support radius
   logical :: write_grids                               ! Write NICAS grids
   integer :: ndir                                      ! Number of Diracs
   real(kind_real) :: londir(ndirmax)                   ! Diracs longitudes (in degrees)
   real(kind_real) :: latdir(ndirmax)                   ! Diracs latitudes (in degrees)
   integer :: levdir(ndirmax)                           ! Diracs level
   integer :: ivdir(ndirmax)                            ! Diracs variable
   integer :: itsdir(ndirmax)                           ! Diracs timeslot

   ! obsop_param
   integer :: nobs                                      ! Number of observations
   character(len=1024) :: obsdis                        ! Observation distribution parameter
   character(len=1024) :: obsop_interp                  ! Observation operator interpolation type

   ! output_param
   integer :: nldwv                                     ! Number of local diagnostics profiles to write (for local_diag = .true.)
   integer ::  img_ldwv(nldwvmax)                       ! Index on model grid of the local diagnostics profiles to write
   real(kind_real) ::  lon_ldwv(nldwvmax)               ! Longitudes (in degrees) of the local diagnostics profiles to write
   real(kind_real) ::  lat_ldwv(nldwvmax)               ! Latitudes (in degrees) of the local diagnostics profiles to write
   character(len=1024),dimension(nldwvmax) :: name_ldwv ! Name of the local diagnostics profiles to write
   real(kind_real) ::  diag_rhflt                       ! Diagnostics filtering radius
   character(len=1024) :: diag_interp                   ! Diagnostics interpolation type
   logical :: grid_output                               ! Write regridded fields
   real(kind_real) :: grid_resol                        ! Regridded fields resolution
   character(len=1024) :: grid_interp                   ! Regridding interpolation type
contains
   procedure :: init => nam_init
   procedure :: read => nam_read
   procedure :: bcast => nam_bcast
   procedure :: setup_internal => nam_setup_internal
   procedure :: check => nam_check
   procedure :: write => nam_write
end type nam_type

private
public :: nvmax,ntsmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax
public :: nam_type

contains

!----------------------------------------------------------------------
! Subroutine: nam_init
! Purpose: intialize
!----------------------------------------------------------------------
subroutine nam_init(nam)

implicit none

! Passed variable
class(nam_type),intent(out) :: nam ! Namelist

! Local variable
integer :: iv,ildwv

! general_param default
nam%datadir = '.'
nam%prefix = ''
nam%model = ''
nam%verbosity = 'all'
nam%colorlog = .false.
nam%default_seed = .false.

! driver_param default
nam%method = ''
nam%strategy = ''
nam%new_cortrack = .false.
nam%new_vbal = .false.
nam%load_vbal = .false.
nam%write_vbal = .true.
nam%new_mom = .true.
nam%load_mom = .false.
nam%write_mom = .false.
nam%new_hdiag = .false.
nam%write_hdiag = .true.
nam%new_lct = .false.
nam%write_lct = .true.
nam%load_cmat = .false.
nam%write_cmat = .true.
nam%new_nicas = .false.
nam%load_nicas = .false.
nam%write_nicas = .true.
nam%new_obsop = .false.
nam%load_obsop = .false.
nam%write_obsop = .true.
nam%check_vbal = .false.
nam%check_adjoints = .false.
nam%check_pos_def = .false.
nam%check_dirac = .false.
nam%check_randomization = .false.
nam%check_consistency = .false.
nam%check_optimality = .false.
nam%check_obsop = .false.

! model_param default
nam%nl = 0
nam%levs = 0
nam%logpres = .false.
nam%nv = 0
do iv=1,nvmax
   nam%varname(iv) = ''
   nam%addvar2d(iv) = ''
end do
nam%nts = 0
nam%timeslot = 0
nam%nomask = .false.

! ens1_param default
nam%ens1_ne = 0
nam%ens1_nsub = 1

! ens2_param default
nam%ens2_ne = 0
nam%ens2_nsub = 1

! sampling_param default
nam%sam_write = .false.
nam%sam_read = .false.
nam%mask_type = 'none'
nam%mask_lu = 'lower'
nam%mask_th = 0.0
nam%ncontig_th = 0
nam%mask_check = .false.
nam%draw_type = 'random_uniform'
nam%Lcoast = 0.0
nam%rcoast = 0.0
nam%nc1 = 0
nam%nc2 = 0
nam%ntry = 0
nam%nrep = 0
nam%nc3 = 0
nam%dc = 0.0
nam%nl0r = 0
nam%irmax = 10000

! diag_param default
nam%ne = 0
nam%gau_approx = .false.
nam%avg_nbins = 0
do iv=1,nvmax*(nvmax-1)/2
   nam%vbal_block(iv) = .false.
end do
nam%vbal_rad = 0.0
nam%var_filter = .false.
nam%var_niter = 0
nam%var_rhflt = 0.0
nam%local_diag = .false.
nam%local_rad = 0.0
nam%adv_diag = .false.
nam%adv_rad = 0.0
nam%adv_niter = 0
nam%adv_rhflt = 0.0

! fit_param default
nam%minim_algo = 'hooke'
do iv=0,nvmax
   nam%double_fit(iv) = .false.
end do
nam%lhomh = .false.
nam%lhomv = .false.
nam%rvflt = 0.0
nam%lct_nscales = 0
nam%lct_diag = .false.

! nicas_param default
nam%nonunit_diag = .false.
nam%lsqrt = .true.
nam%resol = 0.0
nam%fast_sampling = .false.
nam%subsamp = 'hv'
nam%nicas_interp = 'bilin'
nam%network = .false.
nam%mpicom = 0
nam%adv_mode = 0
nam%forced_radii = .false.
nam%rh = 0.0
nam%rv = 0.0
nam%write_grids = .false.
nam%ndir = 0
nam%londir = 0.0
nam%latdir = 0.0
nam%levdir = 0
nam%ivdir = 0
nam%itsdir = 0

! obsop_param default
nam%nobs = 0
nam%obsdis = ''
nam%obsop_interp = 'bilin'

! output_param default
nam%nldwv = 0
nam%img_ldwv = 0
nam%lon_ldwv = 0.0
nam%lat_ldwv = 0.0
do ildwv=1,nldwvmax
   nam%name_ldwv(ildwv) = ''
end do
nam%diag_rhflt = 0.0
nam%diag_interp = 'bilin'
nam%grid_output = .false.
nam%grid_resol = 0.0
nam%grid_interp = 'bilin'

end subroutine nam_init

!----------------------------------------------------------------------
! Subroutine: nam_read
! Purpose: read
!----------------------------------------------------------------------
subroutine nam_read(nam,mpl,namelname)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam     ! Namelist
type(mpl_type),intent(inout) :: mpl      ! MPI data
character(len=*),intent(in) :: namelname ! Namelist name

! Local variables
integer :: iv
character(len=1024),parameter :: subr = 'nam_read'

! Namelist variables
integer :: lunit
integer :: nl,levs(nlmax),nv,nts,timeslot(ntsmax),ens1_ne,ens1_nsub,ens2_ne,ens2_nsub
integer :: ncontig_th,nc1,nc2,ntry,nrep,nc3,nl0r,irmax,ne,avg_nbins,var_niter,adv_niter,lct_nscales,mpicom,adv_mode,ndir
integer :: levdir(ndirmax),ivdir(ndirmax),itsdir(ndirmax),nobs,nldwv,img_ldwv(nldwvmax),ildwv
logical :: colorlog,default_seed
logical :: new_cortrack,new_vbal,load_vbal,write_vbal,new_mom,load_mom,write_mom,new_hdiag,write_hdiag,new_lct,write_lct,load_cmat
logical :: write_cmat,new_nicas,load_nicas,write_nicas,new_obsop,load_obsop,write_obsop,check_vbal,check_adjoints,check_pos_def
logical :: check_dirac,check_randomization,check_consistency,check_optimality,check_obsop,logpres,nomask
logical :: sam_write,sam_read,mask_check,vbal_block(nvmax*(nvmax-1)/2),var_filter,gau_approx,local_diag,adv_diag
logical :: double_fit(0:nvmax),lhomh,lhomv,lct_diag(nscalesmax),nonunit_diag,lsqrt,fast_sampling,network,forced_radii,write_grids
logical :: grid_output
real(kind_real) :: mask_th,Lcoast,rcoast,dc,vbal_rad,var_rhflt,local_rad,adv_rad,adv_rhflt,rvflt,lon_ldwv(nldwvmax)
real(kind_real) :: lat_ldwv(nldwvmax),diag_rhflt,resol,rh,rv,londir(ndirmax),latdir(ndirmax),grid_resol
character(len=1024) :: datadir,prefix,model,verbosity,strategy,method,mask_type,mask_lu,draw_type,minim_algo,subsamp,nicas_interp
character(len=1024) :: obsdis,obsop_interp,diag_interp,grid_interp
character(len=1024),dimension(nvmax) :: varname,addvar2d
character(len=1024),dimension(nldwvmax) :: name_ldwv

! Namelist blocks
namelist/general_param/datadir,prefix,model,verbosity,colorlog,default_seed
namelist/driver_param/method,strategy,new_cortrack,new_vbal,load_vbal,new_mom,load_mom,write_mom,write_vbal,new_hdiag, &
                    & write_hdiag,new_lct,write_lct,load_cmat,write_cmat,new_nicas,load_nicas,write_nicas,new_obsop,load_obsop, &
                    & write_obsop,check_vbal,check_adjoints,check_pos_def,check_dirac,check_randomization,check_consistency, &
                    & check_optimality,check_obsop
namelist/model_param/nl,levs,logpres,nv,varname,addvar2d,nts,timeslot,nomask
namelist/ens1_param/ens1_ne,ens1_nsub
namelist/ens2_param/ens2_ne,ens2_nsub
namelist/sampling_param/sam_write,sam_read,mask_type,mask_lu,mask_th,ncontig_th,mask_check,draw_type,Lcoast,rcoast,nc1,nc2,ntry, &
                      & nrep,nc3,dc,nl0r,irmax
namelist/diag_param/ne,gau_approx,avg_nbins,vbal_block,vbal_rad,var_filter,var_niter,var_rhflt,local_diag,local_rad, &
                  & adv_diag,adv_rad,adv_niter,adv_rhflt
namelist/fit_param/minim_algo,double_fit,lhomh,lhomv,rvflt,lct_nscales,lct_diag
namelist/nicas_param/nonunit_diag,lsqrt,resol,fast_sampling,subsamp,nicas_interp,network,mpicom,adv_mode,forced_radii,rh,rv, &
                   & write_grids,ndir,londir,latdir,levdir,ivdir,itsdir
namelist/obsop_param/nobs,obsdis,obsop_interp
namelist/output_param/nldwv,img_ldwv,lon_ldwv,lat_ldwv,name_ldwv,diag_rhflt,diag_interp,grid_output,grid_resol,grid_interp

if (mpl%main) then
   ! general_param default
   datadir = '.'
   prefix = ''
   model = ''
   verbosity = 'all'
   colorlog = .false.
   default_seed = .false.

   ! driver_param default
   method = ''
   strategy = ''
   new_cortrack = .false.
   new_vbal = .false.
   load_vbal = .false.
   write_vbal = .true.
   new_mom = .true.
   load_mom = .false.
   write_mom = .false.
   new_hdiag = .false.
   write_hdiag = .true.
   new_lct = .false.
   write_lct = .true.
   load_cmat = .false.
   write_cmat = .true.
   new_nicas = .false.
   load_nicas = .false.
   write_nicas = .true.
   new_obsop = .false.
   load_obsop = .false.
   write_obsop = .true.
   check_vbal = .false.
   check_adjoints = .false.
   check_pos_def = .false.
   check_dirac = .false.
   check_randomization = .false.
   check_consistency = .false.
   check_optimality = .false.
   check_obsop = .false.

   ! model_param default
   nl = 0
   levs = 0
   logpres = .false.
   nv = 0
   do iv=1,nvmax
      varname(iv) = ''
      addvar2d(iv) = ''
   end do
   nts = 0
   timeslot = 0
   nomask = .false.

   ! ens1_param default
   ens1_ne = 0
   ens1_nsub = 1

   ! ens2_param default
   ens2_ne = 0
   ens2_nsub = 1

   ! sampling_param default
   sam_write = .false.
   sam_read = .false.
   mask_type = 'none'
   mask_lu = 'lower'
   mask_th = 0.0
   ncontig_th = 0
   mask_check = .false.
   draw_type = 'random_uniform'
   Lcoast = 0.0
   rcoast = 0.0
   nc1 = 0
   nc2 = 0
   ntry = 0
   nrep = 0
   nc3 = 0
   dc = 0.0
   nl0r = 0
   irmax = 10000

   ! diag_param default
   ne = 0
   gau_approx = .false.
   avg_nbins = 0
   do iv=1,nvmax*(nvmax-1)/2
      vbal_block(iv) = .false.
   end do
   vbal_rad = 0.0
   var_filter = .false.
   var_niter = 0
   var_rhflt = 0.0
   local_diag = .false.
   local_rad = 0.0
   adv_diag = .false.
   adv_rad = 0.0
   adv_niter = 0
   adv_rhflt = 0.0

   ! fit_param default
   minim_algo = 'hooke'
   do iv=0,nvmax
      double_fit(iv) = .false.
   end do
   lhomh = .false.
   lhomv = .false.
   rvflt = 0.0
   lct_nscales = 0
   lct_diag = .false.

   ! nicas_param default
   nonunit_diag = .false.
   lsqrt = .true.
   resol = 0.0
   fast_sampling = .false.
   subsamp = 'hv'
   nicas_interp = 'bilin'
   network = .false.
   mpicom = 0
   adv_mode = 0
   forced_radii = .false.
   rh = 0.0
   rv = 0.0
   write_grids = .false.
   ndir = 0
   londir = 0.0
   latdir = 0.0
   levdir = 0
   ivdir = 0
   itsdir = 0

   ! obsop_param default
   nobs = 0
   obsdis = ''
   obsop_interp = 'bilin'

   ! output_param default
   nldwv = 0
   img_ldwv = 0
   lon_ldwv = 0.0
   lat_ldwv = 0.0
   do ildwv=1,nldwvmax
      name_ldwv(ildwv) = ''
   end do
   diag_rhflt = 0.0
   diag_interp = 'bilin'
   grid_output = .false.
   grid_resol = 0.0
   grid_interp = 'bilin'

   ! Open namelist
   call mpl%newunit(lunit)
   open(unit=lunit,file=trim(namelname),status='old',action='read')

   ! general_param
   read(lunit,nml=general_param)
   nam%datadir = datadir
   nam%prefix = prefix
   nam%model = model
   nam%verbosity = verbosity
   nam%colorlog = colorlog
   nam%default_seed = default_seed

   ! driver_param
   read(lunit,nml=driver_param)
   nam%method = method
   nam%strategy = strategy
   nam%new_cortrack = new_cortrack
   nam%new_vbal = new_vbal
   nam%load_vbal = load_vbal
   nam%write_vbal = write_vbal
   nam%new_mom = new_mom
   nam%load_mom = load_mom
   nam%write_mom = write_mom
   nam%new_hdiag = new_hdiag
   nam%write_hdiag = write_hdiag
   nam%new_lct = new_lct
   nam%write_lct = write_lct
   nam%load_cmat = load_cmat
   nam%write_cmat = write_cmat
   nam%new_nicas = new_nicas
   nam%load_nicas = load_nicas
   nam%write_nicas = write_nicas
   nam%new_obsop = new_obsop
   nam%load_obsop = load_obsop
   nam%write_obsop = write_obsop
   nam%check_vbal = check_vbal
   nam%check_adjoints = check_adjoints
   nam%check_pos_def = check_pos_def
   nam%check_dirac = check_dirac
   nam%check_randomization = check_randomization
   nam%check_consistency = check_consistency
   nam%check_optimality = check_optimality
   nam%check_obsop = check_obsop

   ! model_param
   read(lunit,nml=model_param)
   if (nl>nlmax) call mpl%abort(subr,'nl is too large')
   if (nv>nvmax) call mpl%abort(subr,'nv is too large')
   if (nts>ntsmax) call mpl%abort(subr,'nts is too large')
   nam%nl = nl
   if (nl>0) nam%levs(1:nl) = levs(1:nl)
   nam%logpres = logpres
   nam%nv = nv
   if (nv>0) nam%varname(1:nv) = varname(1:nv)
   if (nv>0) nam%addvar2d(1:nv) = addvar2d(1:nv)
   nam%nts = nts
   if (nts>0) nam%timeslot(1:nts) = timeslot(1:nts)
   nam%nomask = nomask

   ! ens1_param
   read(lunit,nml=ens1_param)
   nam%ens1_ne = ens1_ne
   nam%ens1_nsub = ens1_nsub

   ! ens2_param
   read(lunit,nml=ens2_param)
   nam%ens2_ne = ens2_ne
   nam%ens2_nsub = ens2_nsub

   ! sampling_param
   read(lunit,nml=sampling_param)
   if (nc3>nc3max) call mpl%abort(subr,'nc3 is too large')
   nam%sam_write = sam_write
   nam%sam_read = sam_read
   nam%mask_type = mask_type
   nam%mask_lu = mask_lu
   nam%mask_th = mask_th
   nam%ncontig_th = ncontig_th
   nam%mask_check = mask_check
   nam%draw_type = draw_type
   nam%Lcoast = Lcoast
   nam%rcoast = rcoast
   nam%nc1 = nc1
   nam%nc2 = nc2
   nam%ntry = ntry
   nam%nrep = nrep
   nam%nc3 = nc3
   nam%dc = dc
   nam%nl0r = nl0r
   nam%irmax = irmax

   ! diag_param
   read(lunit,nml=diag_param)
   nam%ne = ne
   nam%gau_approx = gau_approx
   nam%avg_nbins = avg_nbins
   if (nv>1) nam%vbal_block(1:nam%nv*(nam%nv-1)/2) = vbal_block(1:nam%nv*(nam%nv-1)/2)
   nam%vbal_rad = vbal_rad
   nam%var_filter = var_filter
   nam%var_niter = var_niter
   nam%var_rhflt = var_rhflt
   nam%local_diag = local_diag
   nam%local_rad = local_rad
   nam%adv_diag = adv_diag
   nam%adv_rad = adv_rad
   nam%adv_niter = adv_niter
   nam%adv_rhflt = adv_rhflt

   ! fit_param
   read(lunit,nml=fit_param)
   if (lct_nscales>nscalesmax) call mpl%abort(subr,'lct_nscales is too large')
   nam%minim_algo = minim_algo
   if (nv>0) nam%double_fit(1:nv) = double_fit(1:nv)
   nam%lhomh = lhomh
   nam%lhomv = lhomv
   nam%rvflt = rvflt
   nam%lct_nscales = lct_nscales
   if (lct_nscales>0) nam%lct_diag(1:lct_nscales) = lct_diag(1:lct_nscales)

   ! nicas_param
   read(lunit,nml=nicas_param)
   if (ndir>ndirmax) call mpl%abort(subr,'ndir is too large')
   nam%nonunit_diag = nonunit_diag
   nam%lsqrt = lsqrt
   nam%resol = resol
   nam%fast_sampling = fast_sampling
   nam%subsamp = subsamp
   nam%nicas_interp = nicas_interp
   nam%network = network
   nam%mpicom = mpicom
   nam%adv_mode = adv_mode
   nam%forced_radii = forced_radii
   nam%rh = rh
   nam%rv = rv
   nam%write_grids = write_grids
   nam%ndir = ndir
   if (ndir>0) nam%londir(1:ndir) = londir(1:ndir)
   if (ndir>0) nam%latdir(1:ndir) = latdir(1:ndir)
   if (ndir>0) nam%levdir(1:ndir) = levdir(1:ndir)
   if (ndir>0) nam%ivdir(1:ndir) = ivdir(1:ndir)
   if (ndir>0) nam%itsdir(1:ndir) = itsdir(1:ndir)

   ! obsop_param
   read(lunit,nml=obsop_param)
   nam%nobs = nobs
   nam%obsdis = obsdis
   nam%obsop_interp = obsop_interp

   ! output_param
   read(lunit,nml=output_param)
   nam%nldwv = nldwv
   if (nldwv>0) then
      nam%img_ldwv(1:nldwv) = img_ldwv(1:nldwv)
      nam%lon_ldwv(1:nldwv) = lon_ldwv(1:nldwv)
      nam%lat_ldwv(1:nldwv) = lat_ldwv(1:nldwv)
      nam%name_ldwv(1:nldwv) = name_ldwv(1:nldwv)
   end if
   nam%diag_rhflt = diag_rhflt
   nam%diag_interp = diag_interp
   nam%grid_output = grid_output
   nam%grid_resol = grid_resol
   nam%grid_interp = grid_interp

   ! Close namelist
   close(unit=lunit)
end if

end subroutine nam_read

!----------------------------------------------------------------------
! Subroutine: nam_bcast
! Purpose: broadcast
!----------------------------------------------------------------------
subroutine nam_bcast(nam,mpl)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam ! Namelist
type(mpl_type),intent(inout) :: mpl  ! MPI data

! general_param
call mpl%f_comm%broadcast(nam%datadir,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%prefix,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%model,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%verbosity,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%colorlog,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%default_seed,mpl%ioproc-1)

! driver_param
call mpl%f_comm%broadcast(nam%method,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%strategy,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%new_cortrack,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%new_vbal,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%load_vbal,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_vbal,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%new_mom,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%load_mom,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_mom,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%new_hdiag,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_hdiag,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%new_lct,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_lct,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%load_cmat,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_cmat,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%new_nicas,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%load_nicas,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_nicas,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%new_obsop,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%load_obsop,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_obsop,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_vbal,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_adjoints,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_pos_def,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_dirac,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_randomization,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_consistency,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_optimality,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%check_obsop,mpl%ioproc-1)

! model_param
call mpl%f_comm%broadcast(nam%nl,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%levs,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%logpres,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nv,mpl%ioproc-1)
call mpl%bcast(nam%varname,mpl%ioproc-1)
call mpl%bcast(nam%addvar2d,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nts,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%timeslot,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nomask,mpl%ioproc-1)

! ens1_param
call mpl%f_comm%broadcast(nam%ens1_ne,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%ens1_nsub,mpl%ioproc-1)

! ens2_param
call mpl%f_comm%broadcast(nam%ens2_ne,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%ens2_nsub,mpl%ioproc-1)

! sampling_param
call mpl%f_comm%broadcast(nam%sam_write,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%sam_read,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%mask_type,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%mask_lu,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%mask_th,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%ncontig_th,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%mask_check,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%draw_type,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%Lcoast,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%rcoast,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nc1,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nc2,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%ntry,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nrep,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nc3,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%dc,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nl0r,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%irmax,mpl%ioproc-1)

! diag_param
call mpl%f_comm%broadcast(nam%ne,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%gau_approx,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%avg_nbins,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%vbal_block,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%vbal_rad,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%var_filter,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%var_niter,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%var_rhflt,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%local_diag,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%local_rad,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%adv_diag,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%adv_rad,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%adv_niter,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%adv_rhflt,mpl%ioproc-1)

! fit_param
call mpl%f_comm%broadcast(nam%minim_algo,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%double_fit,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%lhomh,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%lhomv,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%rvflt,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%lct_nscales,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%lct_diag,mpl%ioproc-1)

! nicas_param
call mpl%f_comm%broadcast(nam%nonunit_diag,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%lsqrt,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%resol,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%fast_sampling,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%subsamp,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%nicas_interp,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%network,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%mpicom,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%adv_mode,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%forced_radii,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%rh,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%rv,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%write_grids,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%ndir,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%londir,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%latdir,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%levdir,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%ivdir,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%itsdir,mpl%ioproc-1)

! obsop_param
call mpl%f_comm%broadcast(nam%nobs,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%obsdis,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%obsop_interp,mpl%ioproc-1)

! output_param
call mpl%f_comm%broadcast(nam%nldwv,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%img_ldwv,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%lon_ldwv,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%lat_ldwv,mpl%ioproc-1)
call mpl%bcast(nam%name_ldwv,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%diag_rhflt,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%diag_interp,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%grid_output,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%grid_resol,mpl%ioproc-1)
call mpl%f_comm%broadcast(nam%grid_interp,mpl%ioproc-1)

end subroutine nam_bcast

!----------------------------------------------------------------------
! Subroutine: nam_setup_internal
! Purpose: setup namelist parameters internally (model 'online')
!----------------------------------------------------------------------
subroutine nam_setup_internal(nam,nl0,nv,nts,ens1_ne,ens1_nsub,ens2_ne,ens2_nsub)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam      ! Namelist
integer,intent(in) :: nl0                 ! Number of levels
integer,intent(in) :: nv                  ! Number of variables
integer,intent(in) :: nts                 ! Number of time-slots
integer,intent(in) :: ens1_ne             ! Ensemble 1 size
integer,intent(in) :: ens1_nsub           ! Ensemble 1 number of sub-ensembles
integer,intent(in) :: ens2_ne             ! Ensemble 2 size
integer,intent(in) :: ens2_nsub           ! Ensemble 2 size of sub-ensembles

! Local variables
integer :: il,iv,its

if (trim(nam%model)=='') then
   nam%model = 'online'
   nam%nl = nl0
   do il=1,nam%nl
      nam%levs(il) = il
   end do
   nam%nv = nv
   do iv=1,nam%nv
      write(nam%varname(iv),'(a,i2.2)') 'var_',iv
      nam%addvar2d(iv) = ''
   end do
   nam%nts = nts
   do its=1,nts
      nam%timeslot(its) = its
   end do
   nam%ens1_ne = ens1_ne
   nam%ens1_nsub = ens1_nsub
   nam%ens2_ne = ens2_ne
   nam%ens2_nsub = ens2_nsub
end if

end subroutine nam_setup_internal

!----------------------------------------------------------------------
! Subroutine: nam_check
! Purpose: check namelist parameters
!----------------------------------------------------------------------
subroutine nam_check(nam,mpl)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam ! Namelist
type(mpl_type),intent(inout) :: mpl  ! MPI data

! Local variables
integer :: iv,its,il,idir,ildwv
character(len=2) :: ivchar,ildwvchar
character(len=1024),parameter :: subr = 'nam_check'

! Check maximum sizes
if (nam%nl>nlmax) call mpl%abort(subr,'nl is too large')
if (nam%nv>nvmax) call mpl%abort(subr,'nv is too large')
if (nam%nts>ntsmax) call mpl%abort(subr,'nts is too large')
if (nam%nc3>nc3max) call mpl%abort(subr,'nc3 is too large')
if (nam%lct_nscales>nscalesmax) call mpl%abort(subr,'lct_nscales is too large')
if (nam%ndir>ndirmax) call mpl%abort(subr,'ndir is too large')
if (nam%nldwv>nldwvmax) call mpl%abort(subr,'nldwv is too large')

! Namelist parameters normalization (meters to radians and degrees to radians)
nam%Lcoast = nam%Lcoast/req
nam%dc = nam%dc/req
nam%vbal_rad = nam%vbal_rad/req
nam%var_rhflt = nam%var_rhflt/req
nam%local_rad = nam%local_rad/req
nam%adv_rad = nam%adv_rad/req
nam%adv_rhflt = nam%adv_rhflt/req
nam%rh = nam%rh/req
if (nam%ndir>0) nam%londir(1:nam%ndir) = nam%londir(1:nam%ndir)*deg2rad
if (nam%ndir>0) nam%latdir(1:nam%ndir) = nam%latdir(1:nam%ndir)*deg2rad
if (nam%nldwv>0) nam%lon_ldwv(1:nam%nldwv) = nam%lon_ldwv(1:nam%nldwv)*deg2rad
if (nam%nldwv>0) nam%lat_ldwv(1:nam%nldwv) = nam%lat_ldwv(1:nam%nldwv)*deg2rad
nam%diag_rhflt = nam%diag_rhflt/req
nam%grid_resol = nam%grid_resol/req

! Check general_param
if (trim(nam%datadir)=='') call mpl%abort(subr,'datadir not specified')
if (trim(nam%prefix)=='') call mpl%abort(subr,'prefix not specified')
select case (trim(nam%model))
case ('aro','arp','fv3','gem','geos','gfs','ifs','mpas','nemo','online','qg','res','wrf')
case default
   call mpl%abort(subr,'wrong model')
end select
select case (trim(nam%verbosity))
case ('all','main','none')
case default
   call mpl%abort(subr,'wrong verbosity level')
end select

! Check driver_param
if (nam%new_hdiag.or.nam%check_optimality) then
   select case (trim(nam%method))
   case ('cor','loc','hyb-avg','hyb-rnd','dual-ens')
   case default
      call mpl%abort(subr,'wrong method')
   end select
end if
if (nam%new_lct) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for new_lct')
end if
if (nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat.or.nam%new_nicas.or.nam%load_nicas) then
   select case (trim(nam%strategy))
   case ('diag_all','common','common_univariate','common_weighted','specific_univariate','specific_multivariate')
   case default
      call mpl%abort(subr,'wrong strategy')
   end select
end if
if (nam%new_vbal.and.nam%load_vbal) call mpl%abort(subr,'new_vbal and load_vbal are exclusive')
if (nam%new_mom.and.nam%load_mom) call mpl%abort(subr,'new_mom and load_mom are exclusive')
if (nam%new_hdiag.and.nam%new_lct) call mpl%abort(subr,'new_hdiag and new_lct are exclusive')
if ((nam%new_hdiag.or.nam%new_lct).and.nam%load_cmat) call mpl%abort(subr,'new_hdiag or new_lct and load_cmat are exclusive')
if (nam%new_nicas.and.nam%load_nicas) call mpl%abort(subr,'new_nicas and load_nicas are exclusive')
if (nam%new_obsop.and.nam%load_obsop) call mpl%abort(subr,'new_obsop and load_obsop are exclusive')
if (nam%check_vbal.and..not.(nam%new_vbal.or.nam%load_vbal)) call mpl%abort(subr,'new_vbal or load_vbal required for check_vbal')
if ((nam%new_hdiag.or.nam%new_lct).and.(.not.(nam%new_mom.or.nam%load_mom))) &
 & call mpl%abort(subr,'new_mom or load_mom required for new_hdiag and new_lct')
if (nam%check_adjoints.and..not.(nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort(subr,'new_nicas or load_nicas required for check_adjoints')
if (nam%check_pos_def.and..not.(nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort(subr,'new_nicas or load_nicas required for check_pos_def')
if (nam%check_dirac.and..not.(nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort(subr,'new_nicas or load_nicas required for check_dirac')
if (nam%check_randomization) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for check_randomization')
   if (.not.nam%new_nicas) call mpl%abort(subr,'new_nicas required for check_randomization')
end if
if (nam%check_consistency) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for check_consistency')
   if (.not.nam%new_nicas) call mpl%abort(subr,'new_nicas required for check_consistency')
end if
if (nam%check_optimality) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for check_optimality')
   if (.not.nam%new_nicas) call mpl%abort(subr,'new_nicas required for check_optimality')
   if (.not.nam%write_hdiag) call mpl%abort(subr,'write_hdiag required for check_optimality')
end if
if (nam%check_obsop.and..not.(nam%new_obsop.or.nam%load_obsop)) &
 & call mpl%abort(subr,'new_obsop or load_obsop required for check_obsop')

! Check model_param
if (nam%nl<=0) call mpl%abort(subr,'nl should be positive')
do il=1,nam%nl
   if (nam%levs(il)<=0) call mpl%abort(subr,'levs should be positive')
   if (count(nam%levs(1:nam%nl)==nam%levs(il))>1) call mpl%abort(subr,'redundant levels')
end do
if (nam%new_vbal.or.nam%load_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat.or.nam%new_nicas.or.nam%load_nicas) then
   if (nam%nv<=0) call mpl%abort(subr,'nv should be positive')
   do iv=1,nam%nv
      write(ivchar,'(i2.2)') iv
      if (trim(nam%varname(iv))=='') call mpl%abort(subr,'varname not specified for variable '//ivchar)
   end do
   if (nam%nts<=0) call mpl%abort(subr,'nts should be positive')
   do its=1,nam%nts
      if (nam%timeslot(its)<0) call mpl%abort(subr,'timeslot should be non-negative')
   end do
   do iv=1,nam%nv
      if (trim(nam%addvar2d(iv))/='') nam%levs(nam%nl+1) = maxval(nam%levs(1:nam%nl))+1
   end do
end if

! Check ens1_param
if (nam%new_cortrack.or.nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%check_randomization &
 & .or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%ens1_nsub<1) call mpl%abort(subr,'ens1_nsub should be positive')
   if (mod(nam%ens1_ne,nam%ens1_nsub)/=0) call mpl%abort(subr,'ens1_nsub should be a divider of ens1_ne')
   if (nam%ens1_ne/nam%ens1_nsub<=3) call mpl%abort(subr,'ens1_ne/ens1_nsub should be larger than 3')
end if

! Check ens2_param
if (nam%new_hdiag.and.((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens'))) then
   if (nam%ens2_nsub<1) call mpl%abort(subr,'ens2_nsub should be non-negative')
   if (mod(nam%ens2_ne,nam%ens2_nsub)/=0) call mpl%abort(subr,'ens2_nsub should be a divider of ens2_ne')
   if (nam%ens2_ne/nam%ens2_nsub<=3) call mpl%abort(subr,'ens2_ne/ens2_nsub should be larger than 3')
end if

! Check sampling_param
if (nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%sam_write.and.nam%sam_read) call mpl%abort(subr,'sam_write and sam_read are both true')
   select case (trim(nam%draw_type))
   case ('random_uniform')
   case ('random_coast')
      if (.not.(nam%Lcoast>0.0)) call mpl%abort (subr,'Lcoast should be positive')
      if (.not.(nam%rcoast>0.0)) call mpl%abort (subr,'rcoast should be positive')
   case default
      call mpl%abort(subr,'wrong draw_type')
   end select
   select case (trim(nam%mask_type))
   case ('ldwv')
      if (nam%nldwv<=0) call mpl%abort(subr,'nldwv should not be negative for mask_type = ldwv')
   case ('stddev')
      select case (trim(nam%mask_lu))
      case ('lower','upper')
      case default
         call mpl%abort(subr,'wrong mask_lu')
      end select
   end select
   if (nam%nc1<3) call mpl%abort(subr,'nc1 should be larger than 2')
   if (nam%new_vbal.or.(nam%new_hdiag.and.(nam%local_diag.or.nam%adv_diag))) then
      if (nam%nc2<3) call mpl%abort(subr,'nc2 should be larger than 2')
   else
      if (nam%nc2<0) then
          call mpl%warning(subr,'nc2 should be set non-negative, resetting nc2 to zero')
          nam%nc2 = 0
      end if
   end if
   if (nam%new_lct) then
      if (nam%nc2/=nam%nc1) then
         call mpl%warning(subr,'nc2 should be equal to nc1 for new_lct, resetting nc2 to nc1')
         nam%nc2 = nam%nc1
      end if
   end if
end if
if (nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%new_nicas) then
   if (nam%ntry<=0) call mpl%abort(subr,'ntry should be positive')
   if (nam%nrep<0) call mpl%abort(subr,'nrep should be non-negative')
end if
if (nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%nc3<=0) call mpl%abort(subr,'nc3 should be positive')
end if
if (nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%dc<0.0) call mpl%abort(subr,'dc should be positive')
end if
if (nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%nl0r<1) call mpl%abort (subr,'nl0r should be positive')
end if
if (nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%irmax<1) call mpl%abort (subr,'irmax should be positive')
end if

! Check diag_param
if (nam%new_vbal) then
   if (nam%nv<2) call mpl%abort(subr,'at least two variables required to diagnose vertical balance')
   if (.not.(any(nam%vbal_block(1:nam%nv*(nam%nv-1)/2)))) &
 & call mpl%abort(subr,'no block selected for the vertical balance diagnostics')

   if (nam%vbal_rad<0.0) call mpl%abort(subr,'vbal_rad should be non-negative')
end if
if (nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
   select case (trim(nam%method))
   case ('loc','hyb-avg','hyb-rnd','dual-ens')
      if (nam%ne<=3) call mpl%abort(subr,'ne should be larger than 3')
   end select
   if (nam%var_filter.and.(.not.nam%local_diag)) call mpl%abort(subr,'local_diag required for var_filter')
   if (nam%var_filter) then
      if (nam%var_niter<=0) call mpl%abort(subr,'var_niter should be positive')
      if (nam%var_rhflt<0.0) call mpl%abort(subr,'var_rhflt should be non-negative')
   end if
   if (nam%local_diag) then
      if (nam%local_rad<0.0) call mpl%abort(subr,'adv_rad should be non-negative')
   end if
   if (nam%adv_diag) then
      if (nam%adv_rad<0.0) call mpl%abort(subr,'local_rad should be non-negative')
      if (nam%adv_niter<=0) call mpl%abort(subr,'adv_niter should be positive')
      if (nam%adv_rhflt<0.0) call mpl%abort(subr,'adv_rhflt should be non-negative')
   end if
end if

! Check fit_param
if (nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   select case (trim(nam%minim_algo))
   case ('none','fast','hooke')
   case default
      call mpl%abort(subr,'wrong minim_algo')
   end select
   if (nam%new_lct.and.((trim(nam%minim_algo)=='none').or.(trim(nam%minim_algo)=='fast'))) &
 & call mpl%abort(subr,'wrong minim_algo for LCT')
   if (nam%rvflt<0) call mpl%abort(subr,'rvflt should be non-negative')
end if
if (nam%new_lct) then
   if (nam%lct_nscales<=0) call mpl%abort(subr,'lct_nscales should be postive')
end if

! Check ensemble sizes
if (nam%new_hdiag) then
   if (trim(nam%method)/='cor') then
      if (nam%ne>nam%ens1_ne) call mpl%warning(subr,'ensemble size larger than ens1_ne (might enhance sampling noise)')
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd','dual-ens')
         if (nam%ne>nam%ens2_ne) call mpl%warning(subr,'ensemble size larger than ens2_ne (might enhance sampling noise)')
      end select
   end if
end if

! Check nicas_param
if (nam%new_nicas.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_dirac.or.nam%check_randomization) then
   if (nam%lsqrt) then
      if (nam%mpicom==1) call mpl%abort(subr,'mpicom should be 2 for square-root application')
   end if
   if (trim(nam%method)=='specific_multivariate') then
      if (.not.nam%lsqrt) call mpl%abort(subr,'square-root formulation required for specific multivariate strategy')
   end if
   if (nam%check_randomization) then
      if (.not.nam%lsqrt) call mpl%abort(subr,'lsqrt required for check_randomization')
      if (.not.nam%forced_radii) call mpl%abort(subr,'forced_radii required for check_randomization')
   end if
   if (nam%check_consistency) then
      if (.not.nam%lsqrt) call mpl%abort(subr,'lsqrt required for check_consistency')
      if (.not.nam%forced_radii) call mpl%abort(subr,'forced_radii required for check_consistency')
   end if
   if (nam%check_optimality) then
      if (.not.nam%lsqrt) call mpl%abort(subr,'lsqrt required for check_optimality')
      if (.not.nam%forced_radii) call mpl%abort(subr,'forced_radii required for check_optimality')
   end if
   if (nam%new_nicas) then
      if (.not.(nam%resol>0.0)) call mpl%abort(subr,'resol should be positive')
   end if
   if (nam%new_nicas.or.nam%load_nicas) then
      if ((nam%mpicom/=1).and.(nam%mpicom/=2)) call mpl%abort(subr,'mpicom should be 1 or 2')
   end if
   if (nam%forced_radii) then
      if (nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat) &
    & call mpl%abort(subr,'new_hdiag, new_lct and load_cmat forbidden for forced_radii')
      if (nam%rh<0.0) call mpl%abort(subr,'rh should be non-negative')
      if (nam%rv<0.0) call mpl%abort(subr,'rv should be non-negative')
   end if
   if (abs(nam%adv_mode)>1) call mpl%abort(subr,'nam%adv_mode should be -1, 0 or 1')
   if (nam%check_dirac) then
      if (nam%ndir<1) call mpl%abort(subr,'ndir should be positive')
      do idir=1,nam%ndir
         if ((nam%londir(idir)<-pi).or.(nam%londir(idir)>pi)) &
       & call mpl%abort(subr,'londir should lie between -180 and 180')
         if ((nam%latdir(idir)<-0.5*pi).or.(nam%latdir(idir)>0.5*pi)) &
       & call mpl%abort(subr,'latdir should lie between -90 and 90')
         if (.not.(any(nam%levs(1:nam%nl)==nam%levdir(idir)).or.(any(nam%addvar2d(1:nam%nv)/='') &
       & .and.(nam%levs(nam%nl+1)==nam%levdir(idir))))) call mpl%abort(subr,'wrong level for a Dirac')
         if ((nam%ivdir(idir)<1).or.(nam%ivdir(idir)>nam%nv)) call mpl%abort(subr,'wrong variable for a Dirac')
         if ((nam%itsdir(idir)<1).or.(nam%itsdir(idir)>nam%nts)) call mpl%abort(subr,'wrong timeslot for a Dirac')
      end do
   end if
   select case (trim(nam%subsamp))
   case ('h','hv','vh','hvh')
   case default
      call mpl%abort(subr,'wrong subsampling structure for NICAS')
   end select
   select case (trim(nam%nicas_interp))
   case ('bilin','natural')
   case default
      call mpl%abort(subr,'wrong interpolation for NICAS')
   end select
end if

! Check obsop_param
if (nam%new_obsop) then
   select case (trim(nam%obsdis))
   case('','random','local','adjusted')
   case default
      call mpl%abort(subr,'wrong observation distribution')
   end select
   select case (trim(nam%obsop_interp))
   case ('bilin','natural')
   case default
      call mpl%abort(subr,'wrong interpolation for observation operator')
   end select
end if

! Check output_param
if (nam%new_hdiag) then
   if (nam%nldwv<0) call mpl%abort(subr,'nldwv should be non-negative')
   if (nam%nldwv>0) then
      if (.not.nam%local_diag) call mpl%abort(subr,'local_diag required for nldwv>0')
         if (.not.all(nam%img_ldwv(1:nam%nldwv)>0)) then
         if (any(nam%lon_ldwv(1:nam%nldwv)<-pi).or.any(nam%lon_ldwv(1:nam%nldwv)>pi)) &
       & call mpl%abort(subr,'lon_ldwv should lie between -180 and 180')
         if (any(nam%lat_ldwv(1:nam%nldwv)<-0.5*pi).or.any(nam%lat_ldwv(1:nam%nldwv)>0.5*pi)) &
       & call mpl%abort(subr,'lat_ldwv should lie between -90 and 90')
         do ildwv=1,nam%nldwv
            write(ildwvchar,'(i2.2)') ildwv
            if (trim(nam%name_ldwv(ildwv))=='') call mpl%abort(subr,'name_ldwv not specified for profile '//ildwvchar)
         end do
      end if
   end if
   if (nam%local_diag.or.nam%adv_diag) then
      if (nam%diag_rhflt<0.0) call mpl%abort(subr,'diag_rhflt should be non-negative')
   end if
end if
if (nam%new_hdiag.or.nam%new_lct) then
   select case (trim(nam%diag_interp))
   case ('bilin','natural')
   case default
      call mpl%abort(subr,'wrong interpolation for diagnostics')
   end select
end if
if (nam%new_hdiag.or.nam%new_nicas.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_dirac &
 & .or.nam%check_randomization.or.nam%new_lct) then
   if (nam%grid_output) then
      if (.not.(nam%grid_resol>0.0)) call mpl%abort(subr,'grid_resol should be positive')
      select case (trim(nam%grid_interp))
      case ('bilin','natural')
      case default
         call mpl%abort(subr,'wrong interpolation for fields regridding')
      end select
   end if
end if

end subroutine nam_check

!----------------------------------------------------------------------
! Subroutine: nam_write
! Purpose: write namelist parameters into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine nam_write(nam,mpl,ncid)

implicit none

! Passed variable
class(nam_type),intent(in) :: nam   ! Namelist
type(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in),optional :: ncid ! NetCDF file ID

! Local variables
integer :: lncid
real(kind_real),allocatable :: londir(:),latdir(:),lon_ldwv(:),lat_ldwv(:)

! Set ncid
lncid = mpl%msv%vali
if (present(ncid)) lncid = ncid

! general_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','General parameters'
   call mpl%flush
end if
call mpl%write(lncid,'datadir',nam%datadir)
call mpl%write(lncid,'prefix',nam%prefix)
call mpl%write(lncid,'model',nam%model)
call mpl%write(lncid,'verbosity',nam%verbosity)
call mpl%write(lncid,'colorlog',nam%colorlog)
call mpl%write(lncid,'default_seed',nam%default_seed)

! driver_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Driver parameters'
   call mpl%flush
end if
call mpl%write(lncid,'method',nam%method)
call mpl%write(lncid,'strategy',nam%strategy)
call mpl%write(lncid,'new_cortrack',nam%new_cortrack)
call mpl%write(lncid,'new_vbal',nam%new_vbal)
call mpl%write(lncid,'load_vbal',nam%load_vbal)
call mpl%write(lncid,'write_vbal',nam%write_vbal)
call mpl%write(lncid,'new_mom',nam%new_mom)
call mpl%write(lncid,'load_mom',nam%load_mom)
call mpl%write(lncid,'write_mom',nam%write_mom)
call mpl%write(lncid,'new_hdiag',nam%new_hdiag)
call mpl%write(lncid,'write_hdiag',nam%write_hdiag)
call mpl%write(lncid,'new_lct',nam%new_lct)
call mpl%write(lncid,'write_lct',nam%write_lct)
call mpl%write(lncid,'load_cmat',nam%load_cmat)
call mpl%write(lncid,'write_cmat',nam%write_cmat)
call mpl%write(lncid,'new_nicas',nam%new_nicas)
call mpl%write(lncid,'load_nicas',nam%load_nicas)
call mpl%write(lncid,'write_nicas',nam%write_nicas)
call mpl%write(lncid,'new_obsop',nam%new_obsop)
call mpl%write(lncid,'load_obsop',nam%load_obsop)
call mpl%write(lncid,'write_obsop',nam%write_obsop)
call mpl%write(lncid,'check_vbal',nam%check_vbal)
call mpl%write(lncid,'check_adjoints',nam%check_adjoints)
call mpl%write(lncid,'check_pos_def',nam%check_pos_def)
call mpl%write(lncid,'check_dirac',nam%check_dirac)
call mpl%write(lncid,'check_randomization',nam%check_randomization)
call mpl%write(lncid,'check_consistency',nam%check_consistency)
call mpl%write(lncid,'check_optimality',nam%check_optimality)
call mpl%write(lncid,'check_obsop',nam%check_obsop)

! model_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Model parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nl',nam%nl)
call mpl%write(lncid,'levs',nam%nl,nam%levs(1:nam%nl))
call mpl%write(lncid,'logpres',nam%logpres)
call mpl%write(lncid,'nv',nam%nv)
call mpl%write(lncid,'varname',nam%nv,nam%varname(1:nam%nv))
call mpl%write(lncid,'addvar2d',nam%nv,nam%addvar2d(1:nam%nv))
call mpl%write(lncid,'nts',nam%nts)
call mpl%write(lncid,'timeslot',nam%nts,nam%timeslot(1:nam%nts))
call mpl%write(lncid,'nomask',nam%nomask)

! ens1_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Ensemble 1 parameters'
   call mpl%flush
end if
call mpl%write(lncid,'ens1_ne',nam%ens1_ne)
call mpl%write(lncid,'ens1_nsub',nam%ens1_nsub)

! ens2_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Ensemble 2 parameters'
   call mpl%flush
end if
call mpl%write(lncid,'ens2_ne',nam%ens2_ne)
call mpl%write(lncid,'ens2_nsub',nam%ens2_nsub)

! sampling_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Sampling parameters'
   call mpl%flush
end if
call mpl%write(lncid,'sam_write',nam%sam_write)
call mpl%write(lncid,'sam_read',nam%sam_read)
call mpl%write(lncid,'mask_type',nam%mask_type)
call mpl%write(lncid,'mask_lu',nam%mask_lu)
call mpl%write(lncid,'mask_th',nam%mask_th)
call mpl%write(lncid,'ncontig_th',nam%ncontig_th)
call mpl%write(lncid,'mask_check',nam%mask_check)
call mpl%write(lncid,'draw_type',nam%draw_type)
call mpl%write(lncid,'Lcoast',nam%Lcoast*req)
call mpl%write(lncid,'rcoast',nam%rcoast)
call mpl%write(lncid,'nc1',nam%nc1)
call mpl%write(lncid,'ntry',nam%ntry)
call mpl%write(lncid,'nrep',nam%nrep)
call mpl%write(lncid,'nc3',nam%nc3)
call mpl%write(lncid,'dc',nam%dc*req)
call mpl%write(lncid,'nl0r',nam%nl0r)
call mpl%write(lncid,'irmax',nam%irmax)

! diag_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Diagnostics parameters'
   call mpl%flush
end if
call mpl%write(lncid,'ne',nam%ne)
call mpl%write(lncid,'gau_approx',nam%gau_approx)
call mpl%write(lncid,'avg_nbins',nam%avg_nbins)
call mpl%write(lncid,'vbal_block',nam%nv*(nam%nv-1)/2,nam%vbal_block(1:nam%nv*(nam%nv-1)/2))
call mpl%write(lncid,'var_filter',nam%var_filter)
call mpl%write(lncid,'var_niter',nam%var_niter)
call mpl%write(lncid,'var_rhflt',nam%var_rhflt*req)
call mpl%write(lncid,'local_diag',nam%local_diag)
call mpl%write(lncid,'local_rad',nam%local_rad*req)
call mpl%write(lncid,'adv_diag',nam%adv_diag)
call mpl%write(lncid,'adv_rad',nam%adv_rad*req)
call mpl%write(lncid,'adv_niter',nam%adv_niter)
call mpl%write(lncid,'adv_rhflt',nam%adv_rhflt*req)

! fit_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Fit parameters'
   call mpl%flush
end if
call mpl%write(lncid,'minim_algo',nam%minim_algo)
call mpl%write(lncid,'double_fit',nam%nv+1,nam%double_fit(0:nam%nv))
call mpl%write(lncid,'lhomh',nam%lhomh)
call mpl%write(lncid,'lhomv',nam%lhomv)
call mpl%write(lncid,'rvflt',nam%rvflt)
call mpl%write(lncid,'lct_nscales',nam%lct_nscales)
call mpl%write(lncid,'lct_diag',nam%lct_nscales,nam%lct_diag(1:nam%lct_nscales))

! nicas_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','NICAS parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nonunit_diag',nam%nonunit_diag)
call mpl%write(lncid,'lsqrt',nam%lsqrt)
call mpl%write(lncid,'resol',nam%resol)
call mpl%write(lncid,'fast_sampling',nam%fast_sampling)
call mpl%write(lncid,'subsamp',nam%subsamp)
call mpl%write(lncid,'nicas_interp',nam%nicas_interp)
call mpl%write(lncid,'network',nam%network)
call mpl%write(lncid,'mpicom',nam%mpicom)
call mpl%write(lncid,'adv_mode',nam%adv_mode)
call mpl%write(lncid,'forced_radii',nam%forced_radii)
call mpl%write(lncid,'rh',nam%rh)
call mpl%write(lncid,'rv',nam%rv)
call mpl%write(lncid,'write_grids',nam%write_grids)
call mpl%write(lncid,'ndir',nam%ndir)
allocate(londir(nam%ndir))
allocate(latdir(nam%ndir))
if (nam%ndir>0) then
   londir = nam%londir(1:nam%ndir)*rad2deg
   latdir = nam%latdir(1:nam%ndir)*rad2deg
end if
call mpl%write(lncid,'londir',nam%ndir,londir)
call mpl%write(lncid,'latdir',nam%ndir,latdir)
call mpl%write(lncid,'levdir',nam%ndir,nam%levdir(1:nam%ndir))
call mpl%write(lncid,'ivdir',nam%ndir,nam%ivdir(1:nam%ndir))
call mpl%write(lncid,'itsdir',nam%ndir,nam%itsdir(1:nam%ndir))

! obsop_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Observation operator parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nobs',nam%nobs)
call mpl%write(lncid,'obsdis',nam%obsdis)
call mpl%write(lncid,'obsop_interp',nam%obsop_interp)

! output_param
if (mpl%msv%isi(lncid)) then
   write(mpl%info,'(a7,a)') '','Output parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nldwv',nam%nldwv)
call mpl%write(lncid,'img_ldwv',nam%nldwv,nam%img_ldwv(1:nam%nldwv))
allocate(lon_ldwv(nam%nldwv))
allocate(lat_ldwv(nam%nldwv))
if (nam%nldwv>0) then
   lon_ldwv = nam%lon_ldwv(1:nam%nldwv)*rad2deg
   lat_ldwv = nam%lat_ldwv(1:nam%nldwv)*rad2deg
end if
call mpl%write(lncid,'lon_ldwv',nam%nldwv,lon_ldwv)
call mpl%write(lncid,'lat_ldwv',nam%nldwv,lat_ldwv)
call mpl%write(lncid,'name_ldwv',nam%nldwv,nam%name_ldwv(1:nam%nldwv))
call mpl%write(lncid,'diag_rhflt',nam%diag_rhflt*req)
call mpl%write(lncid,'diag_interp',nam%diag_interp)
call mpl%write(lncid,'grid_output',nam%grid_output)
call mpl%write(lncid,'grid_resol',nam%grid_resol*req)
call mpl%write(lncid,'grid_interp',nam%grid_interp)

! Release memory
deallocate(londir)
deallocate(latdir)
deallocate(lon_ldwv)
deallocate(lat_ldwv)

end subroutine nam_write

end module type_nam
