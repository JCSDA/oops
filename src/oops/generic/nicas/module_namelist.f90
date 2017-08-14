!----------------------------------------------------------------------
! Module: module_namelist
!> Purpose: namelist parameters management
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_namelist

use netcdf, only: nf90_put_att,nf90_global
use omp_lib, only: omp_get_num_procs
use tools_display, only: msgerror,msgwarning
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr
use tools_nc, only: ncerr
use type_mpl, only: mpl,mpl_bcast

implicit none

! Namelist parameters maximum sizes
integer,parameter :: nlmax = 200   !< Maximum number of levels
integer,parameter :: ndirmax = 100 !< Maximum number of diracs

type namtype
! general_param
character(len=1024) :: datadir     !< Data directory
character(len=1024) :: prefix      !< Files prefix
character(len=1024) :: model       !< Model name ('aro','arp', 'gfs', 'ifs','mpas', 'nemo' or 'wrf')
logical :: colorlog                !< Add colors to the log (for display on terminal)
integer :: nl                      !< Number of levels
integer :: levs(nlmax)             !< Levels
logical :: new_param               !< Compute new parameters (if false, read file)
logical :: new_mpi                 !< Compute new mpi splitting (if false, read file)
logical :: check_adjoints          !< Test adjoints
logical :: check_pos_def           !< Test positive definiteness
logical :: check_mpi               !< Test single proc/multi-procs equivalence
logical :: check_dirac             !< Test NICAS application on diracs
logical :: check_perf              !< Test NICAS performance
integer :: ndir                    !< Number of Diracs
integer :: dirlev                  !< Diracs level
real(kind_real) :: dirlon(ndirmax) !< Diracs longitudes
real(kind_real) :: dirlat(ndirmax) !< Diracs latitudes

! sampling_param
logical :: sam_default_seed        !< Default seed for random numbers
logical :: mask_check              !< Check that interpolations do not cross mask boundaries
integer :: ntry                    !< Number of tries to get the most separated point for the zero-separation sampling
integer :: nrep                    !< Number of replacement to improve homogeneity of the zero-separation sampling
logical :: logpres                 !< Use pressure logarithm as vertical coordinate (model level if .false.)

! nicas_param
logical :: lsqrt                   !< Square-root formulation
character(len=1024) :: Lbh_file    !< Horizontal length-scale file
real(kind_real) :: Lbh(nlmax)      !< Horizontal length-scale
character(len=1024) :: Lbv_file    !< Vertical length-scale file
real(kind_real) :: Lbv(nlmax)      !< Vertical length-scale
real(kind_real) :: resol           !< Resolution
logical :: network                 !< Network-base convolution calculation (distance-based if false)
integer :: nproc                   !< Number of tasks
integer :: mpicom                  !< Number of communication steps
end type namtype

type(namtype) :: nam

interface namncwrite_param
  module procedure namncwrite_integer
  module procedure namncwrite_integer_array
  module procedure namncwrite_real
  module procedure namncwrite_real_array
  module procedure namncwrite_logical
  module procedure namncwrite_string
end interface

private
public :: nam
public :: namread,namcheck,namncwrite

contains

!----------------------------------------------------------------------
! Subroutine: namread
!> Purpose: read and check namelist parameters
!----------------------------------------------------------------------
subroutine namread

implicit none

! Local variables
character(len=1024) :: datadir
character(len=1024) :: prefix
character(len=1024) :: model
logical :: colorlog
integer :: nl
integer :: levs(nlmax)
logical :: new_param
logical :: new_mpi
logical :: check_adjoints
logical :: check_pos_def
logical :: check_mpi
logical :: check_dirac
logical :: check_perf
integer :: ndir
integer :: dirlev
real(kind_real) :: dirlon(ndirmax)
real(kind_real) :: dirlat(ndirmax)
logical :: sam_default_seed
logical :: mask_check
integer :: ntry
integer :: nrep
logical :: logpres
logical :: lsqrt
character(len=1024) :: Lbh_file
real(kind_real) :: Lbh(nlmax)
character(len=1024) :: Lbv_file
real(kind_real) :: Lbv(nlmax)
real(kind_real) :: resol
logical :: network
integer :: nproc
integer :: mpicom 

! Namelist blocks
namelist/general_param/datadir,prefix,colorlog,model,nl,levs,new_param,new_mpi, &
 & check_adjoints,check_pos_def,check_mpi,check_dirac,check_perf,ndir,dirlev,dirlon,dirlat
namelist/sampling_param/sam_default_seed,mask_check,ntry,nrep,logpres
namelist/nicas_param/lsqrt,Lbh_file,Lbh,Lbv_file,Lbv,resol,network,nproc,mpicom

! Default initialization

! general_param
datadir = ''
prefix = ''
colorlog = .false.
model = ''
call msi(nl)
call msi(levs)
new_param = .false.
new_mpi = .false.
check_adjoints = .false.
check_pos_def = .false.
check_mpi = .false.
check_dirac = .false.
check_dirac = .false.
call msi(ndir)
call msi(dirlev)
call msr(dirlon)
call msr(dirlat)

! sampling_param
sam_default_seed = .false.
mask_check = .false.
call msi(ntry)
call msi(nrep)
logpres = .false.

! nicas_param
lsqrt = .false.
Lbh_file = ''
call msr(Lbh)
Lbv_file = ''
call msr(Lbv)
call msr(resol)
network = .false.
call msi(nproc)
call msi(mpicom)

if (mpl%main) then
   ! Read namelist
   read(*,nml=general_param)
   read(*,nml=sampling_param)
   read(*,nml=nicas_param)

   ! Copy into derived type
   nam%datadir = datadir
   nam%prefix = prefix
   nam%colorlog = colorlog
   nam%model = model
   nam%nl = nl
   nam%levs = levs
   nam%new_param = new_param
   nam%new_mpi = new_mpi
   nam%check_adjoints = check_adjoints
   nam%check_pos_def = check_pos_def
   nam%check_mpi = check_mpi
   nam%check_dirac = check_dirac
   nam%check_perf = check_perf
   nam%ndir = ndir
   nam%dirlev = dirlev
   nam%dirlon = dirlon
   nam%dirlat = -dirlat
   nam%sam_default_seed = sam_default_seed
   nam%mask_check = mask_check
   nam%ntry = ntry
   nam%nrep = nrep
   nam%logpres = logpres
   nam%lsqrt = lsqrt
   nam%Lbh_file = Lbh_file
   nam%Lbh = Lbh
   nam%Lbv_file = Lbv_file
   nam%Lbv = Lbv
   nam%resol = resol
   nam%network = network
   nam%nproc = nproc
   nam%mpicom = mpicom
end if

! Broadcast parameters
call mpl_bcast(nam%datadir,mpl%ioproc)
call mpl_bcast(nam%prefix,mpl%ioproc)
call mpl_bcast(nam%colorlog,mpl%ioproc)
call mpl_bcast(nam%model,mpl%ioproc)
call mpl_bcast(nam%nl,mpl%ioproc)
call mpl_bcast(nam%levs,mpl%ioproc)
call mpl_bcast(nam%new_param,mpl%ioproc)
call mpl_bcast(nam%new_mpi,mpl%ioproc)
call mpl_bcast(nam%check_adjoints,mpl%ioproc)
call mpl_bcast(nam%check_pos_def,mpl%ioproc)
call mpl_bcast(nam%check_mpi,mpl%ioproc)
call mpl_bcast(nam%check_dirac,mpl%ioproc)
call mpl_bcast(nam%check_perf,mpl%ioproc)
call mpl_bcast(nam%ndir,mpl%ioproc)
call mpl_bcast(nam%dirlev,mpl%ioproc)
call mpl_bcast(nam%dirlon,mpl%ioproc)
call mpl_bcast(nam%dirlat,mpl%ioproc)
call mpl_bcast(nam%sam_default_seed,mpl%ioproc)
call mpl_bcast(nam%mask_check,mpl%ioproc)
call mpl_bcast(nam%ntry,mpl%ioproc)
call mpl_bcast(nam%nrep,mpl%ioproc)
call mpl_bcast(nam%logpres,mpl%ioproc)
call mpl_bcast(nam%lsqrt,mpl%ioproc)
call mpl_bcast(nam%Lbh_file,mpl%ioproc)
call mpl_bcast(nam%Lbh,mpl%ioproc)
call mpl_bcast(nam%Lbv_file,mpl%ioproc)
call mpl_bcast(nam%Lbv,mpl%ioproc)
call mpl_bcast(nam%resol,mpl%ioproc)
call mpl_bcast(nam%network,mpl%ioproc)
call mpl_bcast(nam%nproc,mpl%ioproc)
call mpl_bcast(nam%mpicom,mpl%ioproc)

end subroutine namread

!----------------------------------------------------------------------
! Subroutine: namcheck
!> Purpose: check namelist parameters
!----------------------------------------------------------------------
subroutine namcheck

implicit none

! Local variables
integer :: il,idir

! Check general_param
if (trim(nam%datadir)=='') call msgerror('datadir not specified')
if (trim(nam%prefix)=='') call msgerror('prefix not specified')
select case (trim(nam%model))
case ('aro','arp','gem','geos','gfs','ifs','mpas','nemo','oops','wrf')
case default
   call msgerror('wrong model')
end select
if (nam%nl<=0) call msgerror('nl should be positive')
do il=1,nam%nl
   if (nam%levs(il)<=0) call msgerror('levs should be positive')
   if (count(nam%levs(1:nam%nl)==nam%levs(il))>1) call msgerror('redundant levels')
end do
if (nam%new_param.and.(.not.nam%new_mpi)) then
   call msgwarning('new parameters calculation implies new MPI splitting, resetting new_mpi to .true.')
   nam%new_mpi = .true.
end if
if (nam%check_dirac) then
   if (nam%ndir<1) call msgerror('ndir should be positive')
   if (.not.any(nam%dirlev==nam%levs(1:nam%nl))) call msgerror('wrong level for a Dirac')
   do idir=1,nam%ndir
      if ((nam%dirlon(idir)<-180.0).or.(nam%dirlon(idir)>180.0)) call msgerror('Dirac longitude should lie between -180 and 180')
      if ((nam%dirlat(idir)<-90.0).or.(nam%dirlat(idir)>90.0)) call msgerror('Dirac latitude should lie between -90 and 90')
   end do
end if

! Check sampling_param
if (nam%ntry<=0) call msgerror('ntry should be positive')
if (nam%nrep<0) call msgerror('nrep should be non-negative')
if (nam%logpres) then
   select case (trim(nam%model))
   case ('aro','arp','gem','geos','gfs','mpas','oops','wrf')
   case default
      call msgwarning('pressure logarithm vertical coordinate is not available for this model, resetting to model level index')
      nam%logpres = .false.
   end select
end if

! Check nicas_param
if (trim(nam%Lbh_file)=='none') then
   nam%Lbh(1:nam%nl) = nam%Lbh(nam%levs(1:nam%nl))
   do il=1,nam%nl
      if (nam%Lbh(il)<tiny(1.0)) call msgerror('Lbh should be positive')
   end do
end if
if (trim(nam%Lbv_file)=='none') then
   nam%Lbv(1:nam%nl) = nam%Lbv(nam%levs(1:nam%nl))
   do il=1,nam%nl
      if (nam%Lbv(il)<tiny(1.0)) call msgerror('Lbv should be positive')
   end do
end if
if (.not.(nam%resol>0.0)) call msgerror('resol should be positive')
if (nam%nproc<1) call msgerror('nproc should be positive')
if ((nam%mpicom/=1).and.(nam%mpicom/=2)) call msgerror('mpicom should be 1 or 2')

! Cross-check
if (nam%check_mpi.or.nam%check_dirac) then
   if (nam%nproc/=mpl%nproc) call msgerror('nam%nproc should be equal to mpl%nproc for parallel tests')
end if

end subroutine namcheck

!----------------------------------------------------------------------
! Subroutine: namncwrite
!> Purpose: write namelist parameters as NetCDF attributes
!----------------------------------------------------------------------
subroutine namncwrite(ncid)

implicit none

! Passed variables
integer,intent(in) :: ncid !< NetCDF file id

! Local variables
character(len=1024) :: subr = 'namncwrite'

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! general_param
call namncwrite_param(ncid,'general_param_datadir',trim(nam%datadir))
call namncwrite_param(ncid,'general_param_prefix',trim(nam%prefix))
call namncwrite_param(ncid,'general_param_colorlog',nam%colorlog)
call namncwrite_param(ncid,'general_param_model',trim(nam%model))
call namncwrite_param(ncid,'general_param_nl',nam%nl)
call namncwrite_param(ncid,'general_param_levs',nam%nl,nam%levs)
call namncwrite_param(ncid,'general_param_new_param',nam%new_param)
call namncwrite_param(ncid,'general_param_new_mpi',nam%new_mpi)
call namncwrite_param(ncid,'general_param_check_adjoints',nam%check_adjoints)
call namncwrite_param(ncid,'general_param_check_pos_def',nam%check_pos_def)
call namncwrite_param(ncid,'general_param_check_mpi',nam%check_mpi)
call namncwrite_param(ncid,'general_param_check_dirac',nam%check_dirac)
call namncwrite_param(ncid,'general_param_ndir',nam%ndir)
call namncwrite_param(ncid,'general_param_dirlon',nam%ndir,nam%dirlon)
call namncwrite_param(ncid,'general_param_dirlat',nam%ndir,nam%dirlat)

! sampling_param
call namncwrite_param(ncid,'sampling_param_sam_default_seed',nam%sam_default_seed)
call namncwrite_param(ncid,'sampling_param_mask_check',nam%mask_check)
call namncwrite_param(ncid,'sampling_param_ntry',nam%ntry)
call namncwrite_param(ncid,'sampling_param_nrep',nam%nrep)
call namncwrite_param(ncid,'sampling_param_logpres',nam%logpres)

! nicas_param
call namncwrite_param(ncid,'nicas_param_lsqrt',nam%lsqrt)
call namncwrite_param(ncid,'nicas_param_Lbh_file',trim(nam%Lbh_file))
call namncwrite_param(ncid,'nicas_param_Lbh',nam%nl,nam%Lbh)
call namncwrite_param(ncid,'nicas_param_Lbv_file',trim(nam%Lbv_file))
call namncwrite_param(ncid,'nicas_param_Lbv',nam%nl,nam%Lbv)
call namncwrite_param(ncid,'nicas_param_resol',nam%resol)
call namncwrite_param(ncid,'nicas_param_network',nam%network)
call namncwrite_param(ncid,'nicas_param_nproc',nam%nproc)
call namncwrite_param(ncid,'nicas_param_mpicom',nam%mpicom)

end subroutine namncwrite

!----------------------------------------------------------------------
! Subroutine: namncwrite_integer
!> Purpose: write namelist integer as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_integer(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: var              !< Integer

! Local variables
character(len=1024) :: subr='namncwrite_integer'

! Write integer
call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))

end subroutine namncwrite_integer

!----------------------------------------------------------------------
! Subroutine: namncwrite_integer_array
!> Purpose: write namelist integer array as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_integer_array(ncid,varname,n,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< Integer array size
integer,intent(in) :: var(n)           !< Integer array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024) :: subr='namncwrite_integer_array'

! Write integer array as a string
if (n>0) then
   write(fullstr,'(i3.3)') var(1)
   do i=2,n
      write(str,'(i3.3)') var(i)
      fullstr = trim(fullstr)//':'//trim(str)
   end do
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine namncwrite_integer_array

!----------------------------------------------------------------------
! Subroutine: namncwrite_real
!> Purpose: write namelist real as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_real(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
real(kind_real),intent(in) :: var      !< Real

! Local variables
character(len=1024) :: subr='namncwrite_real'

! Write real
call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))

end subroutine namncwrite_real

!----------------------------------------------------------------------
! Subroutine: namncwrite_real_array
!> Purpose: write namelist real array as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_real_array(ncid,varname,n,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< Real array size
real(kind_real),intent(in) :: var(n)   !< Real array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024) :: subr='namncwrite_real_array'

! Write real array as a string
if (n>0) then
   write(fullstr,'(e10.3)') var(1)
   do i=2,n
      write(str,'(e10.3)') var(i)
      fullstr = trim(fullstr)//':'//trim(str)
   end do
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine namncwrite_real_array

!----------------------------------------------------------------------
! Subroutine: namncwrite_logical
!> Purpose: write namelist logical as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_logical(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
logical,intent(in) :: var              !< Logical

! Local variables
character(len=1024) :: subr='namncwrite_logical'

! Write logical as a string
if (var) then
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.true.'))
else
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.false.'))
end if

end subroutine namncwrite_logical

!----------------------------------------------------------------------
! Subroutine: namncwrite_string
!> Purpose: write namelist string as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_string(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
character(len=*),intent(in) :: var     !< String

! Local variables
character(len=1024) :: subr='namncwrite_string'

! Write string
call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(var)))

end subroutine namncwrite_string

end module module_namelist
