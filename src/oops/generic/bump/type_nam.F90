!----------------------------------------------------------------------
! Module: type_nam
!> Purpose: namelist derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_nam

use iso_c_binding
use netcdf, only: nf90_put_att,nf90_global
!$ use omp_lib, only: omp_get_num_procs
use tools_const, only: req,deg2rad,rad2deg
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr
use tools_nc, only: put_att
use type_mpl, only: mpl_type

implicit none

integer,parameter :: nvmax = 20     !< Maximum number of variables
integer,parameter :: ntsmax = 20    !< Maximum number of time slots
integer,parameter :: nlmax = 200    !< Maximum number of levels
integer,parameter :: nc3max = 1000  !< Maximum number of classes
integer,parameter :: nscalesmax = 5 !< Maximum number of variables
integer,parameter :: ndirmax = 100  !< Maximum number of diracs
integer,parameter :: nldwvmax = 100 !< Maximum number of local diagnostic profiles

type nam_type
   ! general_param
   character(len=1024) :: datadir                   !< Data directory
   character(len=1024) :: prefix                    !< Files prefix
   character(len=1024) :: model                     !< Model name ('aro', 'arp', 'fv3', 'gem', 'geos', 'gfs', 'ifs', 'mpas', 'nemo' or 'wrf')
   logical :: colorlog                              !< Add colors to the log (for display on terminal)
   logical :: default_seed                          !< Default seed for random numbers
   logical :: use_metis                             !< Use METIS to split the domain between processors

   ! driver_param
   character(len=1024) :: method                    !< Localization/hybridization to compute ('cor', 'loc_norm', 'loc', 'hyb-avg', 'hyb-rnd' or 'dual-ens')
   character(len=1024) :: strategy                  !< Localization strategy ('diag_all', 'common', 'common_univariate', 'common_weighted', 'specific_univariate' or 'specific_multivariate')
   logical :: new_vbal                              !< Compute new vertical balance operator
   logical :: load_vbal                             !< Load existing vertical balance operator
   logical :: new_hdiag                             !< Compute new HDIAG diagnostics
   logical :: new_lct                               !< Compute new LCT
   logical :: load_cmat                             !< Load existing C matrix
   logical :: new_nicas                             !< Compute new NICAS parameters
   logical :: load_nicas                            !< Load existing NICAS parameters
   logical :: new_obsop                             !< Compute new observation operator
   logical :: load_obsop                            !< Load existing observation operator
   logical :: check_vbal                            !< Test vertical balance inverse and adjoint
   logical :: check_adjoints                        !< Test NICAS adjoints
   logical :: check_pos_def                         !< Test NICAS positive definiteness
   logical :: check_sqrt                            !< Test NICAS full/square-root equivalence
   logical :: check_dirac                           !< Test NICAS application on diracs
   logical :: check_randomization                   !< Test NICAS randomization
   logical :: check_consistency                     !< Test HDIAG-NICAS consistency
   logical :: check_optimality                      !< Test HDIAG optimality
   logical :: check_obsop                           !< Test observation operator

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
   character(len=1024) :: draw_type                 !< Sampling draw type ('random_uniform','random_coast' or 'icosahedron')
   integer :: nc1                                   !< Number of sampling points
   integer :: nc2                                   !< Number of diagnostic points
   integer :: ntry                                  !< Number of tries to get the most separated point for the zero-separation sampling
   integer :: nrep                                  !< Number of replacement to improve homogeneity of the zero-separation sampling
   integer :: nc3                                   !< Number of classes
   real(kind_real) ::  dc                           !< Class size (for sam_type='hor'), should be larger than the typical grid cell size
   integer :: nl0r                                  !< Reduced number of levels for diagnostics

   ! diag_param
   integer :: ne                                    !< Ensemble size
   logical :: gau_approx                            !< Gaussian approximation for asymptotic quantities
   logical :: vbal_block(nvmax*(nvmax-1)/2)         !< Activation of vertical balance (ordered line by line in the lower triangular formulation)
   real(kind_real) :: vbal_rad                      !< Vertical balance diagnostic radius
   logical :: var_diag                              !< Compute variances
   logical :: var_filter                            !< Filter variances
   integer :: var_niter                             !< Number of iteration for the variances filtering (for var_filter = .true.)
   real(kind_real) ::  var_rhflt                    !< Variances initial filtering support radius (for var_filter = .true.)
   logical :: var_full                              !< Compute variances on full grid
   logical :: local_diag                            !< Activate local diagnostics
   real(kind_real) ::  local_rad                    !< Local diagnostics calculation radius (for local_rad = .true.)
   logical :: displ_diag                            !< Activate displacement diagnostics
   real(kind_real) ::  displ_rad                    !< Displacement diagnostics calculation radius
   integer :: displ_niter                           !< Number of iteration for the displacement filtering (for displ_diag = .true.)
   real(kind_real) ::  displ_rhflt                  !< Displacement initial filtering support radius (for displ_diag = .true.)
   real(kind_real) ::  displ_tol                    !< Displacement tolerance for mesh check (for displ_diag = .true.)

   ! fit_param
   character(len=1024) :: minim_algo                !< Minimization algorithm ('none', 'fast' or 'hooke')
   logical :: double_fit(0:nvmax)                   !< Double fit to introduce negative lobes on the vertical
   logical :: lhomh                                 !< Vertically homogenous horizontal support radius
   logical :: lhomv                                 !< Vertically homogenous vertical support radius
   real(kind_real) ::  rvflt                        !< Vertical smoother support radius
   integer :: lct_nscales                           !< Number of LCT scales
   logical :: lct_diag(nscalesmax)                  !< Diagnostic of diagonal LCT components only

   ! nicas_param
   logical :: lsqrt                                 !< Square-root formulation
   real(kind_real) :: resol                         !< Resolution
   character(len=1024) :: nicas_interp              !< NICAS interpolation type
   logical :: network                               !< Network-base convolution calculation (distance-based if false)
   integer :: mpicom                                !< Number of communication steps
   integer :: advmode                               !< Advection mode (1: direct, -1: direct and inverse)
   logical :: forced_radii                          !< Force specific support radii
   real(kind_real) :: rh                            !< Forced horizontal support radius
   real(kind_real) :: rv                            !< Forced vertical support radius
   integer :: ndir                                  !< Number of Diracs
   real(kind_real) :: londir(ndirmax)               !< Diracs longitudes (in degrees)
   real(kind_real) :: latdir(ndirmax)               !< Diracs latitudes (in degrees)
   integer :: levdir(ndirmax)                       !< Diracs level
   integer :: ivdir(ndirmax)                        !< Diracs variable
   integer :: itsdir(ndirmax)                       !< Diracs timeslot

   ! obsop_param
   integer :: nobs                                  !< Number of observations
   character(len=1024) :: obsdis                    !< Observation distribution parameter
   character(len=1024) :: obsop_interp              !< Observation operator interpolation type

   ! output_param
   integer :: nldwh                                 !< Number of local diagnostics fields to write (for local_diag = .true.)
   integer :: il_ldwh(nlmax*nc3max)                 !< Levels of local diagnostics fields to write (for local_diag = .true.)
   integer :: ic_ldwh(nlmax*nc3max)                 !< Classes of local diagnostics fields to write (for local_diag = .true.)
   integer :: nldwv                                 !< Number of local diagnostics profiles to write (for local_diag = .true.)
   real(kind_real) ::  lon_ldwv(nldwvmax)           !< Longitudes (in degrees) local diagnostics profiles to write (for local_diag = .true.)
   real(kind_real) ::  lat_ldwv(nldwvmax)           !< Latitudes (in degrees) local diagnostics profiles to write (for local_diag = .true.)
   real(kind_real) ::  diag_rhflt                   !< Diagnostics filtering radius
   character(len=1024) :: diag_interp               !< Diagnostics interpolation type
   logical :: field_io                              !< Field I/O
   logical :: split_io                              !< Split I/O (each task read and write its own file)
   logical :: grid_output                           !< Write regridded fields
   real(kind_real) :: grid_resol                    !< Regridded fields resolution
   character(len=1024) :: grid_interp               !< Regridding interpolation type
contains
   procedure :: init => nam_init
   procedure :: read => nam_read
   procedure :: bcast => nam_bcast
   procedure :: setup_internal => nam_setup_internal
   procedure :: check => nam_check
   procedure :: ncwrite => nam_ncwrite
end type nam_type

private
public :: nam_type

contains

!----------------------------------------------------------------------
! Subroutine: nam_init
!> Purpose: intialize namelist parameters
!----------------------------------------------------------------------
subroutine nam_init(nam)

implicit none

! Passed variable
class(nam_type),intent(out) :: nam !< Namelist

! Local variable
integer :: iv

! general_param default
nam%datadir = ''
nam%prefix = ''
nam%model = ''
nam%colorlog = .false.
nam%default_seed = .false.
nam%use_metis = .false.

! driver_param default
nam%method = ''
nam%strategy = ''
nam%new_vbal = .false.
nam%load_vbal = .false.
nam%new_hdiag = .false.
nam%new_lct = .false.
nam%load_cmat = .false.
nam%new_nicas = .false.
nam%load_nicas = .false.
nam%new_obsop = .false.
nam%load_obsop = .false.
nam%check_vbal = .false.
nam%check_adjoints = .false.
nam%check_pos_def = .false.
nam%check_sqrt = .false.
nam%check_dirac = .false.
nam%check_randomization = .false.
nam%check_consistency = .false.
nam%check_optimality = .false.
nam%check_obsop = .false.

! model_param default
call msi(nam%nl)
call msi(nam%levs)
nam%logpres = .false.
call msi(nam%nv)
do iv=1,nvmax
   nam%varname(iv) = ''
   nam%addvar2d(iv) = ''
end do
call msi(nam%nts)
call msi(nam%timeslot)

! ens1_param default
call msi(nam%ens1_ne)
call msi(nam%ens1_ne_offset)
call msi(nam%ens1_nsub)

! ens2_param default
call msi(nam%ens2_ne)
call msi(nam%ens2_ne_offset)
call msi(nam%ens2_nsub)

! sampling_param default
nam%sam_write = .false.
nam%sam_read = .false.
nam%mask_type = ''
call msr(nam%mask_th)
nam%mask_check = .false.
nam%draw_type = ''
call msi(nam%nc1)
call msi(nam%nc2)
call msi(nam%ntry)
call msi(nam%nrep)
call msi(nam%nc3)
call msr(nam%dc)
call msi(nam%nl0r)

! diag_param default
call msi(nam%ne)
nam%gau_approx = .false.
do iv=1,nvmax*(nvmax-1)/2
   nam%vbal_block(iv) = .false.
end do
call msr(nam%vbal_rad)
nam%var_diag = .false.
nam%var_filter = .false.
nam%var_full = .false.
call msi(nam%var_niter)
call msr(nam%var_rhflt)
nam%local_diag = .false.
call msr(nam%local_rad)
nam%displ_diag = .false.
call msr(nam%displ_rad)
call msi(nam%displ_niter)
call msr(nam%displ_rhflt)
call msr(nam%displ_tol)

! fit_param default
nam%minim_algo = ''
do iv=0,nvmax
   nam%double_fit(iv) = .false.
end do
nam%lhomh = .false.
nam%lhomv = .false.
call msr(nam%rvflt)
call msi(nam%lct_nscales)
nam%lct_diag = .false.

! nicas_param default
nam%lsqrt = .false.
call msr(nam%resol)
nam%nicas_interp = ''
nam%network = .false.
call msi(nam%mpicom)
call msi(nam%advmode)
nam%forced_radii = .false.
call msr(nam%rh)
call msr(nam%rv)
call msi(nam%ndir)
call msr(nam%londir)
call msr(nam%latdir)
call msi(nam%levdir)
call msi(nam%ivdir)
call msi(nam%itsdir)

! obsop_param default
call msi(nam%nobs)
nam%obsdis = ''
nam%obsop_interp = ''

! output_param default
call msi(nam%nldwh)
call msi(nam%il_ldwh)
call msi(nam%ic_ldwh)
call msi(nam%nldwv)
call msr(nam%lon_ldwv)
call msr(nam%lat_ldwv)
call msr(nam%diag_rhflt)
nam%diag_interp = ''
nam%field_io = .true.
nam%split_io = .false.
nam%grid_output = .false.
call msr(nam%grid_resol)
nam%grid_interp = ''

end subroutine nam_init

!----------------------------------------------------------------------
! Subroutine: nam_read
!> Purpose: read namelist parameters
!----------------------------------------------------------------------
subroutine nam_read(nam,mpl,namelname)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam     !< Namelist
type(mpl_type),intent(in) :: mpl         !< MPI data
character(len=*),intent(in) :: namelname !< Namelist name

! Local variables
integer :: iv

! Namelist variables
integer :: lunit
integer :: nl,levs(nlmax),nv,nts,timeslot(ntsmax),ens1_ne,ens1_ne_offset,ens1_nsub,ens2_ne,ens2_ne_offset,ens2_nsub
integer :: nc1,nc2,ntry,nrep,nc3,nl0r,ne,var_niter,displ_niter,lct_nscales,mpicom,advmode,ndir,levdir(ndirmax),ivdir(ndirmax)
integer :: itsdir(ndirmax),nobs,nldwh,il_ldwh(nlmax*nc3max),ic_ldwh(nlmax*nc3max),nldwv
logical :: colorlog,default_seed,use_metis
logical :: new_vbal,load_vbal,new_hdiag,new_lct,load_cmat,new_nicas,load_nicas,new_obsop,load_obsop
logical :: check_vbal,check_adjoints,check_pos_def,check_sqrt,check_dirac,check_randomization,check_consistency,check_optimality
logical :: check_obsop,logpres,sam_write,sam_read,mask_check
logical :: vbal_block(nvmax*(nvmax-1)/2),var_diag,var_filter,var_full,gau_approx,local_diag,displ_diag,double_fit(0:nvmax)
logical :: lhomh,lhomv,lct_diag(nscalesmax),lsqrt,network,forced_radii,field_io,split_io,grid_output
real(kind_real) :: mask_th,dc,vbal_rad,var_rhflt,local_rad,displ_rad,displ_rhflt,displ_tol,rvflt,lon_ldwv(nldwvmax)
real(kind_real) :: lat_ldwv(nldwvmax),diag_rhflt,resol,rh,rv,londir(ndirmax),latdir(ndirmax),grid_resol
character(len=1024) :: datadir,prefix,model,strategy,method,mask_type,draw_type,minim_algo,nicas_interp
character(len=1024) :: obsdis,obsop_interp,diag_interp,grid_interp
character(len=1024),dimension(nvmax) :: varname,addvar2d

! Namelist blocks
namelist/general_param/datadir,prefix,model,colorlog,default_seed,use_metis
namelist/driver_param/method,strategy,new_vbal,load_vbal,new_hdiag,new_lct,load_cmat,new_nicas,load_nicas,new_obsop,load_obsop, &
                    & check_vbal,check_adjoints,check_pos_def,check_sqrt,check_dirac,check_randomization,check_consistency, &
                    & check_optimality,check_obsop
namelist/model_param/nl,levs,logpres,nv,varname,addvar2d,nts,timeslot
namelist/ens1_param/ens1_ne,ens1_ne_offset,ens1_nsub
namelist/ens2_param/ens2_ne,ens2_ne_offset,ens2_nsub
namelist/sampling_param/sam_write,sam_read,mask_type,mask_th,mask_check,draw_type,nc1,nc2,ntry,nrep,nc3,dc,nl0r
namelist/diag_param/ne,gau_approx,vbal_block,vbal_rad,var_diag,var_filter,var_full,var_niter,var_rhflt,local_diag,local_rad, &
                  & displ_diag,displ_rad,displ_niter,displ_rhflt,displ_tol
namelist/fit_param/minim_algo,double_fit,lhomh,lhomv,rvflt,lct_nscales,lct_diag
namelist/nicas_param/lsqrt,resol,nicas_interp,network,mpicom,advmode,forced_radii,rh,rv,ndir,londir,latdir,levdir,ivdir,itsdir
namelist/obsop_param/nobs,obsdis,obsop_interp
namelist/output_param/nldwh,il_ldwh,ic_ldwh,nldwv,lon_ldwv,lat_ldwv,diag_rhflt,diag_interp,field_io,split_io, &
                    & grid_output,grid_resol,grid_interp

if (mpl%main) then
   ! general_param default
   datadir = ''
   prefix = ''
   model = ''
   colorlog = .false.
   default_seed = .false.
   use_metis = .false.

   ! driver_param default
   method = ''
   strategy = ''
   new_vbal = .false.
   load_vbal = .false.
   new_hdiag = .false.
   new_lct = .false.
   load_cmat = .false.
   new_nicas = .false.
   load_nicas = .false.
   new_obsop = .false.
   load_obsop = .false.
   check_vbal = .false.
   check_adjoints = .false.
   check_pos_def = .false.
   check_sqrt = .false.
   check_dirac = .false.
   check_randomization = .false.
   check_consistency = .false.
   check_optimality = .false.
   check_obsop = .false.

   ! model_param default
   call msi(nl)
   call msi(levs)
   logpres = .false.
   call msi(nv)
   do iv=1,nvmax
      varname(iv) = ''
      addvar2d(iv) = ''
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
   draw_type = ''
   call msi(nc1)
   call msi(nc2)
   call msi(ntry)
   call msi(nrep)
   call msi(nc3)
   call msr(dc)
   call msi(nl0r)

   ! diag_param default
   call msi(ne)
   gau_approx = .false.
   do iv=1,nvmax*(nvmax-1)/2
      vbal_block(iv) = .false.
   end do
   call msr(vbal_rad)
   var_diag = .false.
   var_filter = .false.
   var_full = .false.
   call msi(var_niter)
   call msr(var_rhflt)
   local_diag = .false.
   call msr(local_rad)
   displ_diag = .false.
   call msr(displ_rad)
   call msi(displ_niter)
   call msr(displ_rhflt)
   call msr(displ_tol)

   ! fit_param default
   minim_algo = ''
   do iv=0,nvmax
      double_fit(iv) = .false.
   end do
   lhomh = .false.
   lhomv = .false.
   call msr(rvflt)
   call msi(lct_nscales)
   lct_diag = .false.

   ! nicas_param default
   lsqrt = .false.
   call msr(resol)
   nicas_interp = ''
   network = .false.
   call msi(mpicom)
   call msi(advmode)
   forced_radii = .false.
   call msr(rh)
   call msr(rv)
   call msi(ndir)
   call msr(londir)
   call msr(latdir)
   call msi(levdir)
   call msi(ivdir)
   call msi(itsdir)

   ! obsop_param default
   call msi(nobs)
   obsdis = ''
   obsop_interp = ''

   ! output_param default
   call msi(nldwh)
   call msi(il_ldwh)
   call msi(ic_ldwh)
   call msi(nldwv)
   call msr(lon_ldwv)
   call msr(lat_ldwv)
   call msr(diag_rhflt)
   diag_interp = ''
   field_io = .true.
   split_io = .false.
   grid_output = .false.
   call msr(grid_resol)
   grid_interp = ''

   ! Open namelist
   call mpl%newunit(lunit)
   open(unit=lunit,file=trim(namelname),status='old',action='read')

   ! general_param
   read(lunit,nml=general_param)
   nam%datadir = datadir
   nam%prefix = prefix
   nam%model = model
   nam%colorlog = colorlog
   nam%default_seed = default_seed
   nam%use_metis = use_metis

   ! driver_param
   read(lunit,nml=driver_param)
   nam%method = method
   nam%strategy = strategy
   nam%new_vbal = new_vbal
   nam%load_vbal = load_vbal
   nam%new_hdiag = new_hdiag
   nam%new_lct = new_lct
   nam%load_cmat = load_cmat
   nam%new_nicas = new_nicas
   nam%load_nicas = load_nicas
   nam%new_obsop = new_obsop
   nam%load_obsop = load_obsop
   nam%check_vbal = check_vbal
   nam%check_adjoints = check_adjoints
   nam%check_pos_def = check_pos_def
   nam%check_sqrt = check_sqrt
   nam%check_dirac = check_dirac
   nam%check_randomization = check_randomization
   nam%check_consistency = check_consistency
   nam%check_optimality = check_optimality
   nam%check_obsop = check_obsop

   ! model_param
   read(lunit,nml=model_param)
   if (nl>nlmax) call mpl%abort('nl is too large')
   if (nv>nvmax) call mpl%abort('nv is too large')
   if (nts>ntsmax) call mpl%abort('nts is too large')
   nam%nl = nl
   if (nl>0) nam%levs(1:nl) = levs(1:nl)
   nam%logpres = logpres
   nam%nv = nv
   if (nv>0) nam%varname(1:nv) = varname(1:nv)
   if (nv>0) nam%addvar2d(1:nv) = addvar2d(1:nv)
   nam%nts = nts
   if (nts>0) nam%timeslot(1:nts) = timeslot(1:nts)

   ! ens1_param
   read(lunit,nml=ens1_param)
   nam%ens1_ne = ens1_ne
   nam%ens1_ne_offset = ens1_ne_offset
   nam%ens1_nsub = ens1_nsub

   ! ens2_param
   read(lunit,nml=ens2_param)
   nam%ens2_ne = ens2_ne
   nam%ens2_ne_offset = ens2_ne_offset
   nam%ens2_nsub = ens2_nsub

   ! sampling_param
   read(lunit,nml=sampling_param)
   if (nc3>nc3max) call mpl%abort('nc3 is too large')
   nam%sam_write = sam_write
   nam%sam_read = sam_read
   nam%mask_type = mask_type
   nam%mask_th = mask_th
   nam%mask_check = mask_check
   nam%draw_type = draw_type
   nam%nc1 = nc1
   nam%nc2 = nc2
   nam%ntry = ntry
   nam%nrep = nrep
   nam%nc3 = nc3
   nam%dc = dc
   nam%nl0r = nl0r

   ! diag_param
   read(lunit,nml=diag_param)
   nam%ne = ne
   nam%gau_approx = gau_approx
   if (nv>1) nam%vbal_block(1:nam%nv*(nam%nv-1)/2) = vbal_block(1:nam%nv*(nam%nv-1)/2)
   nam%vbal_rad = vbal_rad
   nam%var_diag = var_diag
   nam%var_filter = var_filter
   nam%var_full = var_full
   nam%var_niter = var_niter
   nam%var_rhflt = var_rhflt
   nam%local_diag = local_diag
   nam%local_rad = local_rad
   nam%displ_diag = displ_diag
   nam%displ_rad = displ_rad
   nam%displ_niter = displ_niter
   nam%displ_rhflt = displ_rhflt
   nam%displ_tol = displ_tol

   ! fit_param
   read(lunit,nml=fit_param)
   if (lct_nscales>nscalesmax) call mpl%abort('lct_nscales is too large')
   nam%minim_algo = minim_algo
   if (nv>0) nam%double_fit(1:nv) = double_fit(1:nv)
   nam%lhomh = lhomh
   nam%lhomv = lhomv
   nam%rvflt = rvflt
   nam%lct_nscales = lct_nscales
   if (lct_nscales>0) nam%lct_diag(1:lct_nscales) = lct_diag(1:lct_nscales)

   ! nicas_param
   read(lunit,nml=nicas_param)
   if (ndir>ndirmax) call mpl%abort('ndir is too large')
   nam%lsqrt = lsqrt
   nam%resol = resol
   nam%nicas_interp = nicas_interp
   nam%network = network
   nam%mpicom = mpicom
   nam%advmode = advmode
   nam%forced_radii = forced_radii
   nam%rh = rh
   nam%rv = rv
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
   if (nldwh>nlmax*nc3max) call mpl%abort('nldwh is too large')
   if (nldwv>nldwvmax) call mpl%abort('nldwv is too large')
   nam%nldwh = nldwh
   if (nldwh>0) nam%il_ldwh(1:nldwh) = il_ldwh(1:nldwh)
   if (nldwh>0) nam%ic_ldwh(1:nldwh) = ic_ldwh(1:nldwh)
   nam%nldwv = nldwv
   if (nldwv>0) nam%lon_ldwv(1:nldwv) = lon_ldwv(1:nldwv)
   if (nldwv>0) nam%lat_ldwv(1:nldwv) = lat_ldwv(1:nldwv)
   nam%diag_rhflt = diag_rhflt
   nam%diag_interp = diag_interp
   nam%field_io = field_io
   nam%split_io = split_io
   nam%grid_output = grid_output
   nam%grid_resol = grid_resol
   nam%grid_interp = grid_interp

   ! Close namelist
   close(unit=lunit)
end if

end subroutine nam_read

!----------------------------------------------------------------------
! Subroutine: nam_bcast
!> Purpose: broadcast namelist parameters
!----------------------------------------------------------------------
subroutine nam_bcast(nam,mpl)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam !< Namelist
type(mpl_type),intent(in) :: mpl     !< MPI data

! general_param
call mpl%bcast(nam%datadir)
call mpl%bcast(nam%prefix)
call mpl%bcast(nam%model)
call mpl%bcast(nam%colorlog)
call mpl%bcast(nam%default_seed)
call mpl%bcast(nam%use_metis)

! driver_param
call mpl%bcast(nam%method)
call mpl%bcast(nam%strategy)
call mpl%bcast(nam%new_vbal)
call mpl%bcast(nam%load_vbal)
call mpl%bcast(nam%new_hdiag)
call mpl%bcast(nam%new_lct)
call mpl%bcast(nam%load_cmat)
call mpl%bcast(nam%new_nicas)
call mpl%bcast(nam%load_nicas)
call mpl%bcast(nam%new_obsop)
call mpl%bcast(nam%load_obsop)
call mpl%bcast(nam%check_vbal)
call mpl%bcast(nam%check_adjoints)
call mpl%bcast(nam%check_pos_def)
call mpl%bcast(nam%check_sqrt)
call mpl%bcast(nam%check_dirac)
call mpl%bcast(nam%check_randomization)
call mpl%bcast(nam%check_consistency)
call mpl%bcast(nam%check_optimality)
call mpl%bcast(nam%check_obsop)

! model_param
call mpl%bcast(nam%nl)
call mpl%bcast(nam%levs)
call mpl%bcast(nam%logpres)
call mpl%bcast(nam%nv)
call mpl%bcast(nam%varname)
call mpl%bcast(nam%addvar2d)
call mpl%bcast(nam%nts)
call mpl%bcast(nam%timeslot)

! ens1_param
call mpl%bcast(nam%ens1_ne)
call mpl%bcast(nam%ens1_ne_offset)
call mpl%bcast(nam%ens1_nsub)

! ens2_param
call mpl%bcast(nam%ens2_ne)
call mpl%bcast(nam%ens2_ne_offset)
call mpl%bcast(nam%ens2_nsub)

! sampling_param
call mpl%bcast(nam%sam_write)
call mpl%bcast(nam%sam_read)
call mpl%bcast(nam%mask_type)
call mpl%bcast(nam%mask_th)
call mpl%bcast(nam%mask_check)
call mpl%bcast(nam%draw_type)
call mpl%bcast(nam%nc1)
call mpl%bcast(nam%nc2)
call mpl%bcast(nam%ntry)
call mpl%bcast(nam%nrep)
call mpl%bcast(nam%nc3)
call mpl%bcast(nam%dc)
call mpl%bcast(nam%nl0r)

! diag_param
call mpl%bcast(nam%ne)
call mpl%bcast(nam%gau_approx)
call mpl%bcast(nam%vbal_block)
call mpl%bcast(nam%vbal_rad)
call mpl%bcast(nam%var_diag)
call mpl%bcast(nam%var_filter)
call mpl%bcast(nam%var_full)
call mpl%bcast(nam%var_niter)
call mpl%bcast(nam%var_rhflt)
call mpl%bcast(nam%local_diag)
call mpl%bcast(nam%local_rad)
call mpl%bcast(nam%displ_diag)
call mpl%bcast(nam%displ_rad)
call mpl%bcast(nam%displ_niter)
call mpl%bcast(nam%displ_rhflt)
call mpl%bcast(nam%displ_tol)

! fit_param
call mpl%bcast(nam%minim_algo)
call mpl%bcast(nam%double_fit)
call mpl%bcast(nam%lhomh)
call mpl%bcast(nam%lhomv)
call mpl%bcast(nam%rvflt)
call mpl%bcast(nam%lct_nscales)
call mpl%bcast(nam%lct_diag)

! nicas_param
call mpl%bcast(nam%lsqrt)
call mpl%bcast(nam%resol)
call mpl%bcast(nam%nicas_interp)
call mpl%bcast(nam%network)
call mpl%bcast(nam%mpicom)
call mpl%bcast(nam%advmode)
call mpl%bcast(nam%forced_radii)
call mpl%bcast(nam%rh)
call mpl%bcast(nam%rv)
call mpl%bcast(nam%ndir)
call mpl%bcast(nam%londir)
call mpl%bcast(nam%latdir)
call mpl%bcast(nam%levdir)
call mpl%bcast(nam%ivdir)
call mpl%bcast(nam%itsdir)

! obsop_param
call mpl%bcast(nam%nobs)
call mpl%bcast(nam%obsdis)
call mpl%bcast(nam%obsop_interp)

! output_param
call mpl%bcast(nam%nldwh)
call mpl%bcast(nam%il_ldwh)
call mpl%bcast(nam%ic_ldwh)
call mpl%bcast(nam%nldwv)
call mpl%bcast(nam%lon_ldwv)
call mpl%bcast(nam%lat_ldwv)
call mpl%bcast(nam%diag_rhflt)
call mpl%bcast(nam%diag_interp)
call mpl%bcast(nam%field_io)
call mpl%bcast(nam%split_io)
call mpl%bcast(nam%grid_output)
call mpl%bcast(nam%grid_resol)
call mpl%bcast(nam%grid_interp)

end subroutine nam_bcast

!----------------------------------------------------------------------
! Subroutine: nam_setup_internal
!> Purpose: setup namelist parameters internally (model 'online')
!----------------------------------------------------------------------
subroutine nam_setup_internal(nam,nl0,nv,nts,ens1_ne,ens1_nsub,ens2_ne,ens2_nsub)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam      !< Namelist
integer,intent(in) :: nl0                 !< Number of levels
integer,intent(in) :: nv                  !< Number of variables
integer,intent(in) :: nts                 !< Number of time-slots
integer,intent(in) :: ens1_ne             !< Ensemble 1 size
integer,intent(in) :: ens1_nsub           !< Ensemble 1 number of sub-ensembles
integer,intent(in) :: ens2_ne             !< Ensemble 2 size
integer,intent(in) :: ens2_nsub           !< Ensemble 2 size of sub-ensembles

! Local variables
integer :: il,iv

if (trim(nam%datadir)=='') nam%datadir = '.'
nam%model = 'online'
nam%colorlog = .false.
nam%use_metis = .false.
nam%nl = nl0
do il=1,nam%nl
   nam%levs(il) = il
end do
nam%logpres = .false.
nam%nv = nv
do iv=1,nam%nv
   write(nam%varname(iv),'(a,i2.2)') 'var_',iv
   nam%addvar2d(iv) = ''
end do
nam%nts = nts
nam%timeslot = 0
nam%ens1_ne = ens1_ne
nam%ens1_ne_offset = 0
nam%ens1_nsub = ens1_nsub
nam%ens2_ne = ens2_ne
nam%ens2_ne_offset = 0
nam%ens2_nsub = ens2_nsub

end subroutine nam_setup_internal

!----------------------------------------------------------------------
! Subroutine: nam_check
!> Purpose: check namelist parameters
!----------------------------------------------------------------------
subroutine nam_check(nam,mpl)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam !< Namelist
type(mpl_type),intent(in) :: mpl     !< MPI data

! Local variables
integer :: iv,its,il,idir,ildw,itest
character(len=2) :: ivchar
character(len=4) :: itestchar
character(len=7) :: lonchar,latchar
character(len=1024) :: filename

! Check maximum sizes
if (nam%nl>nlmax) call mpl%abort('nl is too large')
if (nam%nv>nvmax) call mpl%abort('nv is too large')
if (nam%nts>ntsmax) call mpl%abort('nts is too large')
if (nam%nc3>nc3max) call mpl%abort('nc3 is too large')
if (nam%lct_nscales>nscalesmax) call mpl%abort('lct_nscales is too large')
if (nam%ndir>ndirmax) call mpl%abort('ndir is too large')
if (nam%nldwh>nlmax*nc3max) call mpl%abort('nldwh is too large')
if (nam%nldwv>nldwvmax) call mpl%abort('nldwv is too large')

! Namelist parameters normalization (meters to radians and degrees to radians)
nam%dc = nam%dc/req
nam%vbal_rad = nam%vbal_rad/req
nam%var_rhflt = nam%var_rhflt/req
nam%local_rad = nam%local_rad/req
nam%displ_rad = nam%displ_rad/req
nam%displ_rhflt = nam%displ_rhflt/req
nam%rh = nam%rh/req
if (nam%ndir>0) nam%londir(1:nam%ndir) = nam%londir(1:nam%ndir)*deg2rad
if (nam%ndir>0) nam%latdir(1:nam%ndir) = nam%latdir(1:nam%ndir)*deg2rad
if (nam%nldwv>0) nam%lon_ldwv(1:nam%nldwv) = nam%lon_ldwv(1:nam%nldwv)*deg2rad
if (nam%nldwv>0) nam%lat_ldwv(1:nam%nldwv) = nam%lat_ldwv(1:nam%nldwv)*deg2rad
nam%diag_rhflt = nam%diag_rhflt/req
nam%grid_resol = nam%grid_resol/req

! Check general_param
if (trim(nam%datadir)=='') call mpl%abort('datadir not specified')
if (trim(nam%prefix)=='') call mpl%abort('prefix not specified')
select case (trim(nam%model))
case ('aro','arp','fv3','gem','geos','gfs','ifs','mpas','nemo','online','wrf')
case default
   call mpl%abort('wrong model')
end select

! Check driver_param
if (nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
   select case (trim(nam%method))
   case ('cor','loc_norm','loc','hyb-avg','hyb-rnd','dual-ens')
   case default
      call mpl%abort('wrong method')
   end select
end if
if (nam%new_lct) then
   if (trim(nam%method)/='cor') call mpl%abort('new_lct requires cor method')
end if
if (nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat.or.nam%new_nicas.or.nam%load_nicas) then
   select case (trim(nam%strategy))
   case ('diag_all','common','common_univariate','common_weighted','specific_univariate')
   case ('specific_multivariate')
      if (.not.nam%lsqrt) call mpl%abort('specific multivariate strategy requires a square-root formulation')
   case default
      call mpl%abort('wrong strategy')
   end select
end if
if (nam%new_vbal.and.nam%load_vbal) call mpl%abort('new_vbal and load_vbal are exclusive')
if (nam%new_hdiag.and.nam%new_lct) call mpl%abort('new_hdiag and new_lct are exclusive')
if ((nam%new_hdiag.or.nam%new_lct).and.nam%load_cmat) call mpl%abort('new_hdiag or new_lct and load_cmat are exclusive')
if (nam%new_nicas.and.nam%load_nicas) call mpl%abort('new_nicas and load_nicas are exclusive')
if (nam%new_obsop.and.nam%load_obsop) call mpl%abort('new_obsop and load_obsop are exclusive')
if (nam%check_vbal.and..not.(nam%new_vbal.or.nam%load_vbal)) call mpl%abort('check_vbal requires new_vbal or load_vbal')
if (nam%check_adjoints.and..not.(nam%new_nicas.or.nam%load_nicas)) call mpl%abort('check_adjoint requires new_nicas or load_nicas')
if (nam%check_pos_def.and..not.(nam%new_nicas.or.nam%load_nicas)) call mpl%abort('check_pos_def requires new_nicas or load_nicas')
if (nam%check_sqrt.and..not.(nam%new_nicas.or.nam%load_nicas)) call mpl%abort('check_sqrt requires new_nicas or load_nicas')
if (nam%check_dirac.and..not.(nam%new_nicas.or.nam%load_nicas)) call mpl%abort('check_dirac requires new_nicas or load_nicas')
if (nam%check_randomization.and..not.(nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort('check_randomization requires new_nicas or load_nicas')
if (nam%check_consistency.and..not.((nam%new_hdiag.or.nam%load_cmat).and.nam%new_nicas)) &
 & call mpl%abort('check_adjoint requires new_nicas or load_nicas and new_nicas')
if (nam%check_optimality.and..not.(nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort('check_optimality requires new_nicas or load_nicas')
if (nam%check_obsop.and..not.(nam%new_obsop.or.nam%load_obsop)) call mpl%abort('check_obsop requires new_obsop or load_obsop')

! Check model_param
if (nam%nl<=0) call mpl%abort('nl should be positive')
do il=1,nam%nl
   if (nam%levs(il)<=0) call mpl%abort('levs should be positive')
   if (count(nam%levs(1:nam%nl)==nam%levs(il))>1) call mpl%abort('redundant levels')
end do
if (nam%logpres) then
   select case (trim(nam%model))
   case ('nemo')
      call mpl%warning('pressure logarithm vertical coordinate is not available for this model, resetting to model level index')
      nam%logpres = .false.
   end select
end if
if (nam%new_vbal.or.nam%load_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat.or.nam%new_nicas.or.nam%load_nicas) then
   if (nam%nv<=0) call mpl%abort('nv should be positive')
   do iv=1,nam%nv
      write(ivchar,'(i2.2)') iv
      if (trim(nam%varname(iv))=='') call mpl%abort('varname not specified for variable '//ivchar)
   end do
   if (nam%nts<=0) call mpl%abort('nts should be positive')
   do its=1,nam%nts
      if (nam%timeslot(its)<0) call mpl%abort('timeslot should be non-negative')
   end do
   do iv=1,nam%nv
      if (trim(nam%addvar2d(iv))/='') nam%levs(nam%nl+1) = maxval(nam%levs(1:nam%nl))+1
   end do
end if

! Check ens1_param
if (nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct) then
   if (nam%ens1_ne_offset<0) call mpl%abort('ens1_ne_offset should be non-negative')
   if (nam%ens1_nsub<1) call mpl%abort('ens1_nsub should be positive')
   if (mod(nam%ens1_ne,nam%ens1_nsub)/=0) call mpl%abort('ens1_nsub should be a divider of ens1_ne')
   if (nam%ens1_ne/nam%ens1_nsub<=3) call mpl%abort('ens1_ne/ens1_nsub should be larger than 3')
end if

! Check ens2_param
if (nam%new_hdiag) then
   select case (trim(nam%method))
   case ('hyb-rnd','dual-ens')
      if (nam%ens2_ne_offset<0) call mpl%abort('ens2_ne_offset should be non-negative')
      if (nam%ens2_nsub<1) call mpl%abort('ens2_nsub should be non-negative')
      if (mod(nam%ens2_ne,nam%ens2_nsub)/=0) call mpl%abort('ens2_nsub should be a divider of ens2_ne')
      if (nam%ens2_ne/nam%ens2_nsub<=3) call mpl%abort('ens2_ne/ens2_nsub should be larger than 3')
   end select
end if

! Check sampling_param
if (nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct) then
   if (nam%sam_write.and.nam%sam_read) call mpl%abort('sam_write and sam_read are both true')
   select case (trim(nam%draw_type))
   case ('random_uniform','random_coast','icosahedron')
   case default
      call mpl%abort('wrong draw_type')
   end select
   if (nam%nc1<3) call mpl%abort('nc1 should be larger than 2')
   if (nam%new_vbal.or.(nam%new_hdiag.and.(nam%var_diag.or.nam%local_diag.or.nam%displ_diag))) then
      if (nam%nc2<3) call mpl%abort('nc2 should be larger than 2')
   else
      if (nam%nc2<0) then
          call mpl%warning('nc2 should be set non-negative, resetting nc2 to zero')
          nam%nc2 = 0
      end if
   end if
   if (nam%new_lct) then
      call mpl%warning('nc2 should be equal to nc2 for new_lct, resetting nc2 to nc1')
      nam%nc2 = nam%nc1
   end if
end if
if (nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%new_nicas) then
   if (nam%ntry<=0) call mpl%abort('ntry should be positive')
   if (nam%nrep<0) call mpl%abort('nrep should be non-negative')
end if
if (nam%new_hdiag.or.nam%new_lct) then
   if (nam%nc3<=0) call mpl%abort('nc3 should be positive')
end if
if (nam%new_vbal.or.nam%load_vbal.or.nam%new_hdiag.or.nam%new_lct) then
   if (nam%nl0r<1) call mpl%abort ('nl0r should be positive')
end if
if (nam%new_hdiag) then
   if (nam%dc<0.0) call mpl%abort('dc should be positive')
end if

! Check diag_param
if (nam%new_vbal) then
   if (nam%nv<2) call mpl%abort('at least two variables required to diagnose vertical balance')
   if (.not.(any(nam%vbal_block(1:nam%nv*(nam%nv-1)/2)))) &
 & call mpl%abort('no block selected for the vertical balance diagnostics')

   if (nam%vbal_rad<0.0) call mpl%abort('vbal_rad should be non-negative')
end if
if (nam%new_hdiag) then
   select case (trim(nam%method))
   case ('loc_norm','loc','hyb-avg','hyb-rnd','dual-ens')
      if (nam%ne<=3) call mpl%abort('ne should be larger than 3')
   end select
   if (nam%var_diag.and.(.not.trim(nam%method)=='cor')) call mpl%abort('var_diag requires method = cor')
   if (nam%var_filter.and.(.not.nam%var_diag)) call mpl%abort('var_filter requires var_diag')
   if (nam%var_filter) then
      if (nam%var_niter<=0) call mpl%abort('var_niter should be positive')
      if (nam%var_rhflt<0.0) call mpl%abort('var_rhflt should be non-negative')
   end if
   if (nam%local_diag) then
      if (nam%local_rad<0.0) call mpl%abort('displ_rad should be non-negative')
   end if
   if (nam%displ_diag) then
      if (nam%displ_rad<0.0) call mpl%abort('local_rad should be non-negative')
      if (nam%displ_niter<=0) call mpl%abort('displ_niter should be positive')
      if (nam%displ_rhflt<0.0) call mpl%abort('displ_rhflt should be non-negative')
      if (nam%displ_tol<0.0) call mpl%abort('displ_tol should be non-negative')
   end if
end if

! Check fit_param
if (nam%new_hdiag.or.nam%new_lct) then
   select case (trim(nam%minim_algo))
   case ('none','fast','hooke')
   case default
      call mpl%abort('wrong minim_algo')
   end select
   if (nam%new_lct.and.((trim(nam%minim_algo)=='none').or.(trim(nam%minim_algo)=='fast'))) &
 & call mpl%abort('wrong minim_algo for LCT')
   if (nam%rvflt<0) call mpl%abort('rvflt should be non-negative')
end if
if (nam%new_lct) then
   if (nam%lct_nscales<=0) call mpl%abort('lct_nscales should be postive')
end if

! Check ensemble sizes
if (nam%new_hdiag) then
   if (trim(nam%method)/='cor') then
      if (nam%ne>nam%ens1_ne) call mpl%warning('ensemble size larger than ens1_ne (might enhance sampling noise)')
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd','dual-ens')
         if (nam%ne>nam%ens2_ne) call mpl%warning('ensemble size larger than ens2_ne (might enhance sampling noise)')
      end select
   end if
end if

! Check nicas_param
if (nam%new_nicas.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_sqrt.or.nam%check_dirac.or.nam%check_randomization &
 & .or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%lsqrt) then
      if (nam%mpicom==1) call mpl%abort('mpicom should be 2 for square-root application')
   end if
   if (nam%check_randomization) then
      if (.not.nam%lsqrt) call mpl%abort('lsqrt required for check_randomization')
   end if
   if (nam%check_consistency) then
      if (.not.nam%lsqrt) call mpl%abort('lsqrt required for check_consistency')
   end if
   if (nam%check_optimality) then
      if (.not.nam%lsqrt) call mpl%abort('lsqrt required for check_optimality')
   end if
   if (nam%new_nicas) then
      if (.not.(nam%resol>0.0)) call mpl%abort('resol should be positive')
   end if
   if (nam%new_nicas.or.nam%load_nicas) then
      if ((nam%mpicom/=1).and.(nam%mpicom/=2)) call mpl%abort('mpicom should be 1 or 2')
   end if
   if (nam%forced_radii) then
      if (nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat) &
    & call mpl%abort('forced_radii requires new_hdiag, new_lct and load_cmat to be inactive')
      if (nam%rh<0.0) call mpl%abort('rh should be non-negative')
      if (nam%rv<0.0) call mpl%abort('rv should be non-negative')
   end if
   if (abs(nam%advmode)>1) call mpl%abort('nam%advmode should be -1, 0 or 1')
   if (nam%check_dirac) then
      if (nam%ndir<1) call mpl%abort('ndir should be positive')
      do idir=1,nam%ndir
         if ((nam%londir(idir)<-180.0).or.(nam%londir(idir)>180.0)) &
       & call mpl%abort('Dirac longitude should lie between -180 and 180')
         if ((nam%latdir(idir)<-90.0).or.(nam%latdir(idir)>90.0)) call mpl%abort('Dirac latitude should lie between -90 and 90')
         if (.not.(any(nam%levs(1:nam%nl)==nam%levdir(idir)).or.(any(nam%addvar2d(1:nam%nv)/='') &
       & .and.(nam%levs(nam%nl+1)==nam%levdir(idir))))) call mpl%abort('wrong level for a Dirac')
         if ((nam%ivdir(idir)<1).or.(nam%ivdir(idir)>nam%nv)) call mpl%abort('wrong variable for a Dirac')
         if ((nam%itsdir(idir)<1).or.(nam%itsdir(idir)>nam%nts)) call mpl%abort('wrong timeslot for a Dirac')
      end do
   end if
   select case (trim(nam%nicas_interp))
   case ('bilin','natural')
   case default
      call mpl%abort('wrong interpolation for NICAS')
   end select
end if

! Check obsop_param
if (nam%new_obsop) then
   if (nam%nobs<1) call mpl%abort('nobs should be positive')
   select case (trim(nam%obsop_interp))
   case ('bilin','natural')
   case default
      call mpl%abort('wrong interpolation for observation operator')
   end select
end if

! Check output_param
if (nam%new_hdiag) then
   if (nam%local_diag) then
      if (nam%nldwh<0) call mpl%abort('nldwh should be non-negative')
      if (nam%nldwh>0) then
         if (any(nam%il_ldwh(1:nam%nldwh)<0)) call mpl%abort('il_ldwh should be non-negative')
         if (any(nam%il_ldwh(1:nam%nldwh)>nam%nl)) call mpl%abort('il_ldwh should be lower than nl')
         if (any(nam%ic_ldwh(1:nam%nldwh)<0)) call mpl%abort('ic_ldwh should be non-negative')
         if (any(nam%ic_ldwh(1:nam%nldwh)>nam%nc3)) call mpl%abort('ic_ldwh should be lower than nc3')
      end if
      if (nam%nldwv<0) call mpl%abort('nldwv should be non-negative')
      if (nam%nldwv>0) then
         if (any(nam%lon_ldwv(1:nam%nldwv)<-180.0).or.any(nam%lon_ldwv(1:nam%nldwv)>180.0)) call mpl%abort('wrong lon_ldwv')
         if (any(nam%lat_ldwv(1:nam%nldwv)<-90.0).or.any(nam%lat_ldwv(1:nam%nldwv)>90.0)) call mpl%abort('wrong lat_ldwv')
      end if
   end if
   if (nam%local_diag.or.nam%displ_diag) then
      if (nam%diag_rhflt<0.0) call mpl%abort('diag_rhflt should be non-negative')
   end if
end if
if (nam%new_hdiag.or.nam%new_lct) then
   select case (trim(nam%diag_interp))
   case ('bilin','natural')
   case default
      call mpl%abort('wrong interpolation for diagnostics')
   end select
end if
if (nam%new_hdiag.or.nam%new_nicas.or.nam%check_adjoints.or.nam%check_pos_def.or.nam%check_sqrt.or.nam%check_dirac &
 & .or.nam%check_randomization.or.nam%check_consistency.or.nam%check_optimality.or.nam%new_lct) then
   if (nam%grid_output) then
      if (.not.(nam%grid_resol>0.0)) call mpl%abort('grid_resol should be positive')
      select case (trim(nam%grid_interp))
      case ('bilin','natural')
      case default
         call mpl%abort('wrong interpolation for fields regridding')
      end select
   end if
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

      ! rh0
      if (nam%sam_write.and.(trim(nam%draw_type)=='random_coast')) then
         filename = trim(nam%prefix)//'_sampling_rh0.nc'
         call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      end if
   end if

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

      ! rh0
      if (nam%sam_write.and.(trim(nam%draw_type)=='random_coast')) then
         filename = trim(nam%prefix)//'_sampling_rh0.nc'
         call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      end if

      ! C matrix
      if (trim(nam%minim_algo)/='none') then
         filename = trim(nam%prefix)//'_cmat_*.nc'
         call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      end if
   end if

   ! Dirac test
   if (nam%check_dirac) then
      filename = trim(nam%prefix)//'_dirac.nc'
      call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
   end if

   ! Randomization test
   if (nam%check_randomization) then
      do itest=1,10
         write(itestchar,'(i4.4)') itest
         filename = trim(nam%prefix)//'_randomize_'//itestchar//'.nc'
         call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
      end do
   end if

   ! LCT
   if (nam%new_lct) then
      filename = trim(nam%prefix)//'_lct.nc'
      call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
   end if
end if

end subroutine nam_check

!----------------------------------------------------------------------
! Subroutine: nam_ncwrite
!> Purpose: write namelist parameters as NetCDF attributes
!----------------------------------------------------------------------
subroutine nam_ncwrite(nam,mpl,ncid)

implicit none

! Passed variable
class(nam_type),intent(in) :: nam !< Namelist
type(mpl_type),intent(in) :: mpl  !< MPI data
integer,intent(in) :: ncid        !< NetCDF file ID

! general_param
call put_att(mpl,ncid,'datadir',trim(nam%datadir))
call put_att(mpl,ncid,'prefix',trim(nam%prefix))
call put_att(mpl,ncid,'model',trim(nam%model))
call put_att(mpl,ncid,'colorlog',nam%colorlog)
call put_att(mpl,ncid,'default_seed',nam%default_seed)
call put_att(mpl,ncid,'use_metis',nam%use_metis)

! driver_param
call put_att(mpl,ncid,'method',trim(nam%method))
call put_att(mpl,ncid,'strategy',trim(nam%strategy))
call put_att(mpl,ncid,'new_vbal',nam%new_vbal)
call put_att(mpl,ncid,'load_vbal',nam%load_vbal)
call put_att(mpl,ncid,'new_hdiag',nam%new_hdiag)
call put_att(mpl,ncid,'new_lct',nam%new_lct)
call put_att(mpl,ncid,'load_cmat',nam%load_cmat)
call put_att(mpl,ncid,'new_nicas',nam%new_nicas)
call put_att(mpl,ncid,'load_nicas',nam%load_nicas)
call put_att(mpl,ncid,'new_obsop',nam%new_obsop)
call put_att(mpl,ncid,'load_obsop',nam%load_obsop)
call put_att(mpl,ncid,'check_vbal',nam%check_vbal)
call put_att(mpl,ncid,'check_adjoints',nam%check_adjoints)
call put_att(mpl,ncid,'check_pos_def',nam%check_pos_def)
call put_att(mpl,ncid,'check_sqrt',nam%check_sqrt)
call put_att(mpl,ncid,'check_dirac',nam%check_dirac)
call put_att(mpl,ncid,'check_randomization',nam%check_randomization)
call put_att(mpl,ncid,'check_consistency',nam%check_consistency)
call put_att(mpl,ncid,'check_optimality',nam%check_optimality)
call put_att(mpl,ncid,'check_obsop',nam%check_obsop)

! model_param
call put_att(mpl,ncid,'nl',nam%nl)
call put_att(mpl,ncid,'levs',nam%nl,nam%levs(1:nam%nl))
call put_att(mpl,ncid,'logpres',nam%logpres)
call put_att(mpl,ncid,'nv',nam%nv)
call put_att(mpl,ncid,'varname',nam%nv,nam%varname(1:nam%nv))
call put_att(mpl,ncid,'addvar2d',nam%nv,nam%addvar2d(1:nam%nv))
call put_att(mpl,ncid,'nts',nam%nts)
call put_att(mpl,ncid,'timeslot',nam%nts,nam%timeslot(1:nam%nts))

! ens1_param
call put_att(mpl,ncid,'ens1_ne',nam%ens1_ne)
call put_att(mpl,ncid,'ens1_ne_offset',nam%ens1_ne_offset)
call put_att(mpl,ncid,'ens1_nsub',nam%ens1_nsub)

! ens2_param
call put_att(mpl,ncid,'ens2_ne',nam%ens2_ne)
call put_att(mpl,ncid,'ens2_ne_offset',nam%ens2_ne_offset)
call put_att(mpl,ncid,'ens2_nsub',nam%ens2_nsub)

! sampling_param
call put_att(mpl,ncid,'sam_write',nam%sam_write)
call put_att(mpl,ncid,'sam_read',nam%sam_read)
call put_att(mpl,ncid,'mask_type',nam%mask_type)
call put_att(mpl,ncid,'mask_th',nam%mask_th)
call put_att(mpl,ncid,'mask_check',nam%mask_check)
call put_att(mpl,ncid,'draw_type',nam%draw_type)
call put_att(mpl,ncid,'nc1',nam%nc1)
call put_att(mpl,ncid,'ntry',nam%ntry)
call put_att(mpl,ncid,'nrep',nam%nrep)
call put_att(mpl,ncid,'nc3',nam%nc3)
call put_att(mpl,ncid,'dc',nam%dc*req)
call put_att(mpl,ncid,'nl0r',nam%nl0r)

! vbal_param
call put_att(mpl,ncid,'vbal_block',nam%nv*(nam%nv-1)/2,nam%vbal_block(1:nam%nv*(nam%nv-1)/2))

! diag_param
call put_att(mpl,ncid,'ne',nam%ne)
call put_att(mpl,ncid,'gau_approx',nam%gau_approx)
call put_att(mpl,ncid,'var_diag',nam%var_diag)
call put_att(mpl,ncid,'var_filter',nam%var_filter)
call put_att(mpl,ncid,'var_full',nam%var_full)
call put_att(mpl,ncid,'var_niter',nam%var_niter)
call put_att(mpl,ncid,'var_rhflt',nam%var_rhflt*req)
call put_att(mpl,ncid,'local_diag',nam%local_diag)
call put_att(mpl,ncid,'local_rad',nam%local_rad*req)
call put_att(mpl,ncid,'displ_diag',nam%displ_diag)
call put_att(mpl,ncid,'displ_rad',nam%displ_rad*req)
call put_att(mpl,ncid,'displ_niter',nam%displ_niter)
call put_att(mpl,ncid,'displ_rhflt',nam%displ_rhflt*req)
call put_att(mpl,ncid,'displ_tol',nam%displ_tol)

! fit_param
call put_att(mpl,ncid,'minim_algo',nam%minim_algo)
call put_att(mpl,ncid,'double_fit',nam%nv+1,nam%double_fit(0:nam%nv))
call put_att(mpl,ncid,'lhomh',nam%lhomh)
call put_att(mpl,ncid,'lhomv',nam%lhomv)
call put_att(mpl,ncid,'rvflt',nam%rvflt)
call put_att(mpl,ncid,'lct_nscales',nam%lct_nscales)
call put_att(mpl,ncid,'lct_diag',nam%lct_nscales,nam%lct_diag)

! nicas_param
call put_att(mpl,ncid,'lsqrt',nam%lsqrt)
call put_att(mpl,ncid,'resol',nam%resol)
call put_att(mpl,ncid,'nicas_interp',nam%nicas_interp)
call put_att(mpl,ncid,'network',nam%network)
call put_att(mpl,ncid,'mpicom',nam%mpicom)
call put_att(mpl,ncid,'advmode',nam%advmode)
call put_att(mpl,ncid,'ndir',nam%ndir)
call put_att(mpl,ncid,'londir',nam%ndir,nam%londir(1:nam%ndir)*rad2deg)
call put_att(mpl,ncid,'latdir',nam%ndir,nam%latdir(1:nam%ndir)*rad2deg)
call put_att(mpl,ncid,'levdir',nam%ndir,nam%levdir(1:nam%ndir))
call put_att(mpl,ncid,'ivdir',nam%ndir,nam%ivdir(1:nam%ndir))
call put_att(mpl,ncid,'itsdir',nam%ndir,nam%itsdir(1:nam%ndir))

! obsop_param
call put_att(mpl,ncid,'nobs',nam%nobs)
call put_att(mpl,ncid,'obsdis',nam%obsdis)
call put_att(mpl,ncid,'obsop_interp',nam%obsop_interp)

! output_param
call put_att(mpl,ncid,'nldwh',nam%nldwh)
call put_att(mpl,ncid,'il_ldwh',nam%nldwh,nam%il_ldwh(1:nam%nldwh))
call put_att(mpl,ncid,'ic_ldwh',nam%nldwh,nam%ic_ldwh(1:nam%nldwh))
call put_att(mpl,ncid,'nldwv',nam%nldwv)
call put_att(mpl,ncid,'lon_ldwv',nam%nldwv,nam%lon_ldwv(1:nam%nldwv)*rad2deg)
call put_att(mpl,ncid,'lat_ldwv',nam%nldwv,nam%lat_ldwv(1:nam%nldwv)*rad2deg)
call put_att(mpl,ncid,'diag_rhflt',nam%diag_rhflt*req)
call put_att(mpl,ncid,'diag_interp',nam%diag_interp)
call put_att(mpl,ncid,'field_io',nam%field_io)
call put_att(mpl,ncid,'split_io',nam%split_io)
call put_att(mpl,ncid,'grid_output',nam%grid_output)
call put_att(mpl,ncid,'grid_resol',nam%grid_resol*req)
call put_att(mpl,ncid,'grid_interp',nam%grid_interp)

end subroutine nam_ncwrite

end module type_nam
