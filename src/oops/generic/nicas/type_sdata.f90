!----------------------------------------------------------------------
! Module: type_sdata
!> Purpose: sampling data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_sdata

use esmf
use module_namelist, only: nam,namncwrite
use netcdf
use tools_const, only: rad2deg
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr
use tools_nc, only: ncerr,ncfloat
use type_com, only: comtype,com_read,com_write
use type_fields, only: fldtype,alphatype
use type_linop, only: linoptype,linop_alloc,linop_read,linop_write
use type_mpl, only: mpl
use type_randgen, only: randgentype,create_randgen,delete_randgen
implicit none

! Local sampling data derived type
type sdatampitype
   ! Number of points
   integer :: nc0a                       !< Number of points in subset Sc0 on halo A
   integer :: nl0                        !< Number of levels in subset Sl0
   integer :: nc1b                       !< Number of points in subset Sc1 on halo B
   integer :: nl1                        !< Number of levels in subset Sl1
   integer,allocatable :: nc2b(:)        !< Number of points in subset Sc2 on halo B
   integer :: nsa                        !< Number of subgrid nodes on halo A
   integer :: nsb                        !< Number of subgrid nodes on halo B
   integer :: nsc                        !< Number of subgrid nodes on halo C

   ! Inter-halo conversions
   integer,allocatable :: isa_to_isb(:)  !< Subgrid, halo A to halo B
   integer,allocatable :: isa_to_isc(:)  !< Subgrid, halo A to halo C
   integer,allocatable :: isb_to_isc(:)  !< Subgrid, halo B to halo B

   ! Level-related
   integer :: nl0i                       !< Number of independent levels
   integer,allocatable :: vbot(:)        !< Bottom level

   ! Linear operators
   type(linoptype) :: c                  !< Convolution
   type(linoptype),allocatable :: h(:)   !< Horizontal interpolation
   type(linoptype),allocatable :: v(:)   !< Vertical interpolation
   type(linoptype),allocatable :: s(:)   !< Subsample interpolation

   ! Copy conversions
   integer,allocatable :: isb_to_ic2b(:) !< Subgrid to subset Sc2 on halo B
   integer,allocatable :: isb_to_il1(:)  !< Subgrid to subset Sl1 on halo B

   ! Normalization
   type(fldtype) :: norm                 !< Normalization factor
   type(alphatype) :: norm_sqrt          !< Internal normalization factor for the square-root formulation

   ! Communications
   type(comtype) :: AB                   !< Communication between halos A and B
   type(comtype) :: AC                   !< Communication between halos A and C

   ! MPI-related conversions
   integer,allocatable :: isa_to_is(:)   !< Subgrid, halo A to global
   integer,allocatable :: is_to_isa(:)   !< Subgrid, global to halo A
   integer,allocatable :: isb_to_is(:)   !< Subgrid, halo B to global
   integer,allocatable :: is_to_isb(:)   !< Subgrid, global to halo B
   integer,allocatable :: isc_to_is(:)   !< Subgrid, halo C to global
   integer,allocatable :: is_to_isc(:)   !< Subgrid, global to halo C
end type sdatampitype

! Sampling data derived type
type sdatatype
   ! Vector sizes
   integer :: nlon                         !< Longitude size
   integer :: nlat                         !< Latitude size
   integer :: nlev                         !< Number of levels
   logical,allocatable :: rgmask(:,:)      !< Reduced Gaussian grid mask
   real(kind_real),allocatable :: area(:)             !< Domain area
   logical :: regional                     !< LAM (if .true.) or global

   ! Vector coordinates
   real(kind_real),allocatable :: lon(:)              !< Cells longitude
   real(kind_real),allocatable :: lat(:)              !< Cells latitude
   logical,allocatable :: mask(:,:)        !< Cells mask
   real(kind_real),allocatable :: vunit(:)            !< Vertical unit
   integer :: nfor                         !< Number of forced sampling points
   integer,allocatable :: ifor(:)          !< Forced sampling points indices

   ! Sampling properties
   logical,allocatable :: llev(:)          !< Vertical interpolation key
   integer,allocatable :: na(:)            !< Number of arcs
   integer,allocatable :: larc(:,:,:)      !< Arcs indices
   integer,allocatable :: grid_nnb(:)      !< Number of neighbors on the full grid
   integer,allocatable :: grid_inb(:,:)    !< Neighbors indices on the full grid
   real(kind_real),allocatable :: grid_dnb(:,:)       !< Neighbors distances on the full grid
   integer,allocatable :: subgrid_nnb(:)   !< Number of neighbors on the subgrid
   integer,allocatable :: subgrid_inb(:,:) !< Neighbors indices on the subgrid
   real(kind_real),allocatable :: subgrid_dnb(:,:)    !< Neighbors distances on the subgrid

   ! ESMF fields
   type(esmf_field) :: c0field
   type(esmf_field) :: c1field
   type(esmf_field),allocatable :: c2field(:)

   ! Boundary nodes
   integer,allocatable :: nbnd(:)          !< Number of boundary nodes
   real(kind_real),allocatable :: xbnd(:,:,:)         !< Boundary nodes, x-coordinate
   real(kind_real),allocatable :: ybnd(:,:,:)         !< Boundary nodes, y-coordinate
   real(kind_real),allocatable :: zbnd(:,:,:)         !< Boundary nodes, z-coordinate
   real(kind_real),allocatable :: vbnd(:,:,:)         !< Boundary nodes, orthogonal vector

   ! Random number generator
   type(randgentype) :: rng                !< Random number generator

   ! NICAS global data

   ! Number of points
   integer :: nc0                          !< Number of points in subset Sc0
   integer :: nl0                          !< Number of levels in subset Sl0
   integer :: nc1                          !< Number of points in subset Sc1
   integer :: nl1                          !< Number of levels in subset Sl1
   integer,allocatable :: nc2(:)           !< Number of points in subset Sc2
   integer :: ns                           !< Number of subgrid nodes

   ! Level-related
   integer :: nl0i                         !< Number of independent levels
   integer,allocatable :: vbot(:)          !< Bottom level

   ! Linear operators
   type(linoptype) :: c                    !< Convolution
   type(linoptype),allocatable :: h(:)     !< Horizontal interpolation
   type(linoptype),allocatable :: v(:)     !< Vertical interpolation
   type(linoptype),allocatable :: s(:)     !< Subsample interpolation

   ! Normalization
   type(fldtype) :: norm                   !< Normalization factor
   type(alphatype) :: norm_sqrt            !< Internal normalization factor for the square-root formulation

   ! Other data

   ! Parameters/normalization conversion
   integer,allocatable :: is_to_ic1(:)       !< Subgrid to subset Sc1
   integer,allocatable :: is_to_il1(:)       !< Subgrid to subset Sl1
   integer,allocatable :: is_to_ic2(:)       !< Subgrid to subset sc2
   integer,allocatable :: ic1_to_ic0(:)      !< Subset Sc1 to subset Sc0
   integer,allocatable :: il1_to_il0(:)      !< Subset Sl1 to subset Sl0
   integer,allocatable :: ic2il1_to_ic0(:,:) !< Grid Gs to subset Sc0
   integer,allocatable :: ic0_to_ic1(:)      !< Subset Sc0 to subset Sc1
   integer,allocatable :: il0_to_il1(:)      !< Subset Sl0 to subset Sl1
   integer,allocatable :: ic0il0_to_is(:,:)  !< Grid Gf to subgrid
   integer,allocatable :: ic2il1_to_is(:,:)  !< Grid Gs to subgrid
   integer,allocatable :: ic2il1_to_ic1(:,:) !< Grid Gs to subset Sc1
   integer,allocatable :: ic1il1_to_is(:,:)  !< Grid Gv to subgrid

   ! NICAS local data
   integer :: nproc                          !< Number of tasks
   integer,allocatable :: ic0_to_iproc(:)    !< Subset Sc0 to local task
   integer,allocatable :: ic0_to_ic0a(:)     !< Subset Sc0, global to halo A
   integer :: mpicom                         !< Number of communication steps
   type(sdatampitype),allocatable :: mpi(:)  !< Local sampling data

   ! Illustration
   integer,allocatable :: halo(:)            !< Halo points
end type sdatatype

private
public :: sdatampitype,sdatatype
public :: sdata_alloc,sdata_read_param,sdata_read_local,sdata_read_mpi,sdata_write_param,sdata_write_mpi

contains

!----------------------------------------------------------------------
! Subroutine: sdata_alloc
!> Purpose: sdata object allocation for grid parameters
!----------------------------------------------------------------------
subroutine sdata_alloc(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Allocation
allocate(sdata%lon(sdata%nc0))
allocate(sdata%lat(sdata%nc0))
allocate(sdata%mask(sdata%nc0,sdata%nl0))
allocate(sdata%vunit(sdata%nl0))

! Initialization
call msr(sdata%lon)
call msr(sdata%lat)
sdata%mask = .false.
call msr(sdata%vunit)
sdata%rng = create_randgen()

end subroutine sdata_alloc

!----------------------------------------------------------------------
! Subroutine: sdata_read_param
!> Purpose: read sdata object
!----------------------------------------------------------------------
subroutine sdata_read_param(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ncid
integer :: nc0_id,nl0_id,nc1_id,nl1_id,ns_id
integer :: vbot_id,nc2_id,is_to_ic1_id,is_to_il1_id,is_to_ic2_id,ic0il0_to_is_id,ic2il1_to_ic0_id
integer :: ic2il1_to_is_id,ic1_to_ic0_id,il1_to_il0_id,ic0_to_ic1_id,il0_to_il1_id,ic2il1_to_ic1_id
integer :: ic1il1_to_is_id,norm_id,norm_sqrt_id
character(len=1024) :: filename
character(len=1024) :: subr = 'sdata_read_param'

! Open file and get dimensions
filename = trim(nam%prefix)//'_param.nc'
call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'nc0',nc0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc0_id,len=sdata%nc0))
call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=sdata%nl0))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'nl0i',sdata%nl0i))
call ncerr(subr,nf90_inq_dimid(ncid,'nc1',nc1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc1_id,len=sdata%nc1))
call ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=sdata%nl1))
call ncerr(subr,nf90_inq_dimid(ncid,'ns',ns_id))
call ncerr(subr,nf90_inquire_dimension(ncid,ns_id,len=sdata%ns))

! Allocation
allocate(sdata%vbot(sdata%nc1))
allocate(sdata%nc2(sdata%nl1))
allocate(sdata%is_to_ic1(sdata%ns))
allocate(sdata%is_to_il1(sdata%ns))
allocate(sdata%is_to_ic2(sdata%ns))
allocate(sdata%ic0il0_to_is(sdata%nc0,sdata%nl0))
allocate(sdata%ic2il1_to_ic0(sdata%nc1,sdata%nl1))
allocate(sdata%ic2il1_to_is(sdata%nc1,sdata%nl1))
allocate(sdata%ic1_to_ic0(sdata%nc1))
allocate(sdata%il1_to_il0(sdata%nl1))
allocate(sdata%ic0_to_ic1(sdata%nc0))
allocate(sdata%il0_to_il1(sdata%nl0))
allocate(sdata%ic2il1_to_ic1(sdata%nc1,sdata%nl1))
allocate(sdata%ic1il1_to_is(sdata%nc1,sdata%nl1))
allocate(sdata%norm%val(sdata%nc0,sdata%nl0))
if (nam%lsqrt) allocate(sdata%norm_sqrt%val(sdata%ns))

! Read data
call ncerr(subr,nf90_inq_varid(ncid,'vbot',vbot_id))
call ncerr(subr,nf90_inq_varid(ncid,'nc2',nc2_id))
call ncerr(subr,nf90_inq_varid(ncid,'is_to_ic1',is_to_ic1_id))
call ncerr(subr,nf90_inq_varid(ncid,'is_to_il1',is_to_il1_id))
call ncerr(subr,nf90_inq_varid(ncid,'is_to_ic2',is_to_ic2_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic0il0_to_is',ic0il0_to_is_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic2il1_to_ic0',ic2il1_to_ic0_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic2il1_to_is',ic2il1_to_is_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic1_to_ic0',ic1_to_ic0_id))
call ncerr(subr,nf90_inq_varid(ncid,'il1_to_il0',il1_to_il0_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic0_to_ic1',ic0_to_ic1_id))
call ncerr(subr,nf90_inq_varid(ncid,'il0_to_il1',il0_to_il1_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic2il1_to_ic1',ic2il1_to_ic1_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic1il1_to_is',ic1il1_to_is_id))
call ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
if (nam%lsqrt) call ncerr(subr,nf90_inq_varid(ncid,'norm_sqrt',norm_sqrt_id))
call ncerr(subr,nf90_get_var(ncid,vbot_id,sdata%vbot))
call ncerr(subr,nf90_get_var(ncid,nc2_id,sdata%nc2))
call ncerr(subr,nf90_get_var(ncid,is_to_ic1_id,sdata%is_to_ic1))
call ncerr(subr,nf90_get_var(ncid,is_to_il1_id,sdata%is_to_il1))
call ncerr(subr,nf90_get_var(ncid,is_to_ic2_id,sdata%is_to_ic2))
call ncerr(subr,nf90_get_var(ncid,ic0il0_to_is_id,sdata%ic0il0_to_is))
call ncerr(subr,nf90_get_var(ncid,ic2il1_to_ic0_id,sdata%ic2il1_to_ic0))
call ncerr(subr,nf90_get_var(ncid,ic2il1_to_is_id,sdata%ic2il1_to_is))
call ncerr(subr,nf90_get_var(ncid,ic1_to_ic0_id,sdata%ic1_to_ic0))
call ncerr(subr,nf90_get_var(ncid,il1_to_il0_id,sdata%il1_to_il0))
call ncerr(subr,nf90_get_var(ncid,ic0_to_ic1_id,sdata%ic0_to_ic1))
call ncerr(subr,nf90_get_var(ncid,il0_to_il1_id,sdata%il0_to_il1))
call ncerr(subr,nf90_get_var(ncid,ic2il1_to_ic1_id,sdata%ic2il1_to_ic1))
call ncerr(subr,nf90_get_var(ncid,ic1il1_to_is_id,sdata%ic1il1_to_is))
call ncerr(subr,nf90_get_var(ncid,norm_id,sdata%norm%val))
if (nam%lsqrt) call ncerr(subr,nf90_get_var(ncid,norm_sqrt_id,sdata%norm_sqrt%val))

! Read linear operators
call linop_read(ncid,'c',sdata%c)
call linop_read(ncid,'h',sdata%h)
call linop_read(ncid,'v',sdata%v)
call linop_read(ncid,'s',sdata%s)

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine sdata_read_param

!----------------------------------------------------------------------
! Subroutine: sdata_read_local
!> Purpose: read sdata object
!----------------------------------------------------------------------
subroutine sdata_read_local(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ic0,info
integer :: ncid,ic0_to_iproc_id,ic0_to_ic0a_id
character(len=4) :: nprocchar
character(len=1024) :: filename
character(len=1024) :: subr = 'sdata_read_local'

! Allocation
allocate(sdata%ic0_to_iproc(sdata%nc0))
allocate(sdata%ic0_to_ic0a(sdata%nc0))

if (sdata%nproc==1) then
   ! All points on a single processor
   sdata%ic0_to_iproc = 1
   do ic0=1,sdata%nc0
      sdata%ic0_to_ic0a(ic0) = ic0
   end do
elseif (sdata%nproc>1) then
   ! Open file
   write(nprocchar,'(i4.4)') sdata%nproc
   filename = trim(nam%prefix)//'_distribution_'//nprocchar//'.nc'
   info = nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid)

   if (info==nf90_noerr) then
      ! Read data and close file
      call ncerr(subr,nf90_inq_varid(ncid,'ic0_to_iproc',ic0_to_iproc_id))
      call ncerr(subr,nf90_inq_varid(ncid,'ic0_to_ic0a',ic0_to_ic0a_id))
      call ncerr(subr,nf90_get_var(ncid,ic0_to_iproc_id,sdata%ic0_to_iproc))
      call ncerr(subr,nf90_get_var(ncid,ic0_to_ic0a_id,sdata%ic0_to_ic0a))
      call ncerr(subr,nf90_close(ncid))
   else
      ! Generate a distribution
      call msgerror('distribution not implemented yet')
   end if
end if

! Check
if (maxval(sdata%ic0_to_iproc)>sdata%nproc) call msgerror('wrong distribution')

end subroutine sdata_read_local
subroutine sdata_read_mpi(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: iproc
integer :: ncid,info
integer :: nc0a_id,nl0_id,nc1b_id,nl1_id,nsa_id,nsb_id
integer :: vbot_id,nc2b_id,isb_to_ic2b_id,isb_to_il1_id
integer :: isa_to_isb_id,isa_to_isc_id,isb_to_isc_id
integer :: norm_id
character(len=1) :: mpicomchar
character(len=4) :: nprocchar,iprocchar
character(len=1024) :: filename
character(len=1024) :: subr = 'sdata_read_mpi'

! Allocation
allocate(sdata%mpi(sdata%nproc))
write(mpicomchar,'(i1)') sdata%mpicom
write(nprocchar,'(i4.4)') sdata%nproc

do iproc=1,sdata%nproc
   ! Open file and get dimensions
   write(iprocchar,'(i4.4)') iproc
   filename = trim(nam%prefix)//'_mpi-'//mpicomchar//'_'//nprocchar//'-'//iprocchar//'.nc'
   call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
   call ncerr(subr,nf90_inq_dimid(ncid,'nc0a',nc0a_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,nc0a_id,len=sdata%mpi(iproc)%nc0a))
   call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=sdata%mpi(iproc)%nl0))
   call ncerr(subr,nf90_get_att(ncid,nf90_global,'nl0i',sdata%mpi(iproc)%nl0i))
   info = nf90_inq_dimid(ncid,'nc1b',nc1b_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nc1b_id,len=sdata%mpi(iproc)%nc1b))
   else
      sdata%mpi(iproc)%nc1b = 0
   end if
   call ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=sdata%mpi(iproc)%nl1))
   info = nf90_inq_dimid(ncid,'nsa',nsa_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nsa_id,len=sdata%mpi(iproc)%nsa))
   else
      sdata%mpi(iproc)%nsa = 0
   end if
   info = nf90_inq_dimid(ncid,'nsb',nsb_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nsb_id,len=sdata%mpi(iproc)%nsb))
   else
      sdata%mpi(iproc)%nsb = 0
   end if
   call ncerr(subr,nf90_get_att(ncid,nf90_global,'nsc',sdata%mpi(iproc)%nsc))

   ! Allocation
   if (sdata%mpi(iproc)%nc1b>0) allocate(sdata%mpi(iproc)%vbot(sdata%mpi(iproc)%nc1b))
   allocate(sdata%mpi(iproc)%nc2b(sdata%mpi(iproc)%nl1))
   if (sdata%mpi(iproc)%nsb>0) allocate(sdata%mpi(iproc)%isb_to_ic2b(sdata%mpi(iproc)%nsb))
   if (sdata%mpi(iproc)%nsb>0) allocate(sdata%mpi(iproc)%isb_to_il1(sdata%mpi(iproc)%nsb))
   if (sdata%mpi(iproc)%nsa>0) allocate(sdata%mpi(iproc)%isa_to_isb(sdata%mpi(iproc)%nsa))
   if (sdata%mpi(iproc)%nsa>0) allocate(sdata%mpi(iproc)%isa_to_isc(sdata%mpi(iproc)%nsa))
   if (sdata%mpi(iproc)%nsb>0) allocate(sdata%mpi(iproc)%isb_to_isc(sdata%mpi(iproc)%nsb))
   allocate(sdata%mpi(iproc)%norm%vala(sdata%mpi(iproc)%nc0a,sdata%mpi(iproc)%nl0))

   ! Read data
   if (sdata%mpi(iproc)%nc1b>0) call ncerr(subr,nf90_inq_varid(ncid,'vbot',vbot_id))
   call ncerr(subr,nf90_inq_varid(ncid,'nc2b',nc2b_id))
   if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'isb_to_ic2b',isb_to_ic2b_id))
   if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'isb_to_il1',isb_to_il1_id))
   if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'isa_to_isb',isa_to_isb_id))
   if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'isa_to_isc',isa_to_isc_id))
   if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'isb_to_isc',isb_to_isc_id))
   call ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
   if (sdata%mpi(iproc)%nc1b>0) call ncerr(subr,nf90_get_var(ncid,vbot_id,sdata%mpi(iproc)%vbot))
   call ncerr(subr,nf90_get_var(ncid,nc2b_id,sdata%mpi(iproc)%nc2b))
   if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_get_var(ncid,isb_to_ic2b_id,sdata%mpi(iproc)%isb_to_ic2b))
   if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_get_var(ncid,isb_to_il1_id,sdata%mpi(iproc)%isb_to_il1))
   if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_get_var(ncid,isa_to_isb_id,sdata%mpi(iproc)%isa_to_isb))
   if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_get_var(ncid,isa_to_isc_id,sdata%mpi(iproc)%isa_to_isc))
   if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_get_var(ncid,isb_to_isc_id,sdata%mpi(iproc)%isb_to_isc))
   call ncerr(subr,nf90_get_var(ncid,norm_id,sdata%mpi(iproc)%norm%vala))

   ! Read communications
   if (sdata%nproc>1) then
      call com_read(ncid,'AB',sdata%mpi(iproc)%AB)
      call com_read(ncid,'AC',sdata%mpi(iproc)%AC)
   end if

   ! Read linear operators
   call linop_read(ncid,'c',sdata%mpi(iproc)%c)
   call linop_read(ncid,'h',sdata%mpi(iproc)%h)
   call linop_read(ncid,'v',sdata%mpi(iproc)%v)
   call linop_read(ncid,'s',sdata%mpi(iproc)%s)

   ! Close file
   call ncerr(subr,nf90_close(ncid))
end do

end subroutine sdata_read_mpi

!----------------------------------------------------------------------
! Subroutine: sdata_write_param
!> Purpose: write sdata object
!----------------------------------------------------------------------
subroutine sdata_write_param(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ncid
integer :: nc0_id,nc1_id,nl0_id,nl1_id,ns_id
integer :: vbot_id,nc2_id,is_to_ic1_id,is_to_il1_id,is_to_ic2_id,ic0il0_to_is_id,ic2il1_to_ic0_id
integer :: ic2il1_to_is_id,ic1_to_ic0_id,il1_to_il0_id,ic0_to_ic1_id,il0_to_il1_id,ic2il1_to_ic1_id
integer :: ic1il1_to_is_id,norm_id,norm_sqrt_id
character(len=1024) :: filename
character(len=1024) :: subr = 'sdata_write_param'

if (mpl%main) then
   ! Create file
   filename = trim(nam%prefix)//'_param.nc'
   call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call namncwrite(ncid)

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,'nc0',sdata%nc0,nc0_id))
   call ncerr(subr,nf90_def_dim(ncid,'nc1',sdata%nc1,nc1_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl0',sdata%nl0,nl0_id))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',sdata%nl0i))
   call ncerr(subr,nf90_def_dim(ncid,'nl1',sdata%nl1,nl1_id))
   call ncerr(subr,nf90_def_dim(ncid,'ns',sdata%ns,ns_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,'vbot',nf90_int,(/nc1_id/),vbot_id))
   call ncerr(subr,nf90_def_var(ncid,'nc2',nf90_int,(/nl1_id/),nc2_id))
   call ncerr(subr,nf90_def_var(ncid,'is_to_ic1',nf90_int,(/ns_id/),is_to_ic1_id))
   call ncerr(subr,nf90_def_var(ncid,'is_to_il1',nf90_int,(/ns_id/),is_to_il1_id))
   call ncerr(subr,nf90_def_var(ncid,'is_to_ic2',nf90_int,(/ns_id/),is_to_ic2_id))
   call ncerr(subr,nf90_def_var(ncid,'ic0il0_to_is',nf90_int,(/nc0_id,nl0_id/),ic0il0_to_is_id))
   call ncerr(subr,nf90_def_var(ncid,'ic2il1_to_ic0',nf90_int,(/nc1_id,nl1_id/),ic2il1_to_ic0_id))
   call ncerr(subr,nf90_def_var(ncid,'ic2il1_to_is',nf90_int,(/nc1_id,nl1_id/),ic2il1_to_is_id))
   call ncerr(subr,nf90_def_var(ncid,'ic1_to_ic0',nf90_int,(/nc1_id/),ic1_to_ic0_id))
   call ncerr(subr,nf90_def_var(ncid,'il1_to_il0',nf90_int,(/nl1_id/),il1_to_il0_id))
   call ncerr(subr,nf90_def_var(ncid,'ic0_to_ic1',nf90_int,(/nc0_id/),ic0_to_ic1_id))
   call ncerr(subr,nf90_def_var(ncid,'il0_to_il1',nf90_int,(/nl0_id/),il0_to_il1_id))
   call ncerr(subr,nf90_def_var(ncid,'ic2il1_to_ic1',nf90_int,(/nc1_id,nl1_id/),ic2il1_to_ic1_id))
   call ncerr(subr,nf90_def_var(ncid,'ic1il1_to_is',nf90_int,(/nc1_id,nl1_id/),ic1il1_to_is_id))
   call ncerr(subr,nf90_def_var(ncid,'norm',ncfloat,(/nc0_id,nl0_id/),norm_id))
   if (nam%lsqrt) call ncerr(subr,nf90_def_var(ncid,'norm_sqrt',ncfloat,(/ns_id/),norm_sqrt_id))
   call ncerr(subr,nf90_enddef(ncid))

   ! Write variables
   call ncerr(subr,nf90_put_var(ncid,vbot_id,sdata%vbot))
   call ncerr(subr,nf90_put_var(ncid,nc2_id,sdata%nc2))
   call ncerr(subr,nf90_put_var(ncid,is_to_ic1_id,sdata%is_to_ic1))
   call ncerr(subr,nf90_put_var(ncid,is_to_il1_id,sdata%is_to_il1))
   call ncerr(subr,nf90_put_var(ncid,is_to_ic2_id,sdata%is_to_ic2))
   call ncerr(subr,nf90_put_var(ncid,ic0il0_to_is_id,sdata%ic0il0_to_is))
   call ncerr(subr,nf90_put_var(ncid,ic2il1_to_ic0_id,sdata%ic2il1_to_ic0))
   call ncerr(subr,nf90_put_var(ncid,ic2il1_to_is_id,sdata%ic2il1_to_is))
   call ncerr(subr,nf90_put_var(ncid,ic1_to_ic0_id,sdata%ic1_to_ic0))
   call ncerr(subr,nf90_put_var(ncid,il1_to_il0_id,sdata%il1_to_il0))
   call ncerr(subr,nf90_put_var(ncid,ic0_to_ic1_id,sdata%ic0_to_ic1))
   call ncerr(subr,nf90_put_var(ncid,il0_to_il1_id,sdata%il0_to_il1))
   call ncerr(subr,nf90_put_var(ncid,ic2il1_to_ic1_id,sdata%ic2il1_to_ic1))
   call ncerr(subr,nf90_put_var(ncid,ic1il1_to_is_id,sdata%ic1il1_to_is))
   call ncerr(subr,nf90_put_var(ncid,norm_id,sdata%norm%val))
   if (nam%lsqrt) call ncerr(subr,nf90_put_var(ncid,norm_sqrt_id,sdata%norm_sqrt%val))

   ! Write linear operators
   call linop_write(ncid,sdata%c)
   call linop_write(ncid,sdata%h)
   call linop_write(ncid,sdata%v)
   call linop_write(ncid,sdata%s)

   ! Close file
   call ncerr(subr,nf90_close(ncid))
end if

end subroutine sdata_write_param

!----------------------------------------------------------------------
! Subroutine: sdata_write_mpi
!> Purpose: write sdata object
!----------------------------------------------------------------------
subroutine sdata_write_mpi(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: iproc,il0i,il1,h_sum,s_sum
integer :: ncid
integer :: nc0a_id,nl0_id,nc1b_id,nl1_id,nsa_id,nsb_id
integer :: vbot_id,nc2b_id,isb_to_ic2b_id,isb_to_il1_id
integer :: c_n_s_id,h_n_s_id,s_n_s_id
integer :: isa_to_isb_id,isa_to_isc_id,isb_to_isc_id
integer :: norm_id
integer :: AB_jhalocounts_id,AB_jexclcounts_id
integer :: AC_jhalocounts_id,AC_jexclcounts_id
integer :: nc0_id,nproc1_id,nproc2_id
integer :: lon_id,lat_id,ic0_to_iproc_id,halo_id,nc0a_per_proc_id,nsa_per_proc_id
character(len=1) :: mpicomchar
character(len=4) :: nprocchar,iprocchar
character(len=1024) :: filename
character(len=1024) :: subr = 'sdata_write_mpi'

if (mpl%main) then
   write(mpicomchar,'(i1)') sdata%mpicom
   write(nprocchar,'(i4.4)') sdata%nproc

   do iproc=1,sdata%nproc
      ! Create file
      write(iprocchar,'(i4.4)') iproc
      filename = trim(nam%prefix)//'_mpi-'//mpicomchar//'_'//nprocchar//'-'//iprocchar//'.nc'
      call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

      ! Write namelist parameters
      call namncwrite(ncid)

      ! Define dimensions
      call ncerr(subr,nf90_def_dim(ncid,'nc0a',sdata%mpi(iproc)%nc0a,nc0a_id))
      call ncerr(subr,nf90_def_dim(ncid,'nl0',sdata%mpi(iproc)%nl0,nl0_id))
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',sdata%mpi(iproc)%nl0i))
      if (sdata%mpi(iproc)%nc1b>0) call ncerr(subr,nf90_def_dim(ncid,'nc1b',sdata%mpi(iproc)%nc1b,nc1b_id))
      call ncerr(subr,nf90_def_dim(ncid,'nl1',sdata%mpi(iproc)%nl1,nl1_id))
      if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_def_dim(ncid,'nsa',sdata%mpi(iproc)%nsa,nsa_id))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_def_dim(ncid,'nsb',sdata%mpi(iproc)%nsb,nsb_id))
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'nsc',sdata%mpi(iproc)%nsc))

      ! Define variables
      if (sdata%mpi(iproc)%nc1b>0) call ncerr(subr,nf90_def_var(ncid,'vbot',nf90_int,(/nc1b_id/),vbot_id))
      call ncerr(subr,nf90_def_var(ncid,'nc2b',nf90_int,(/nl1_id/),nc2b_id))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'isb_to_ic2b',nf90_int,(/nsb_id/),isb_to_ic2b_id))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'isb_to_il1',nf90_int,(/nsb_id/),isb_to_il1_id))
      if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_def_var(ncid,'isa_to_isb',nf90_int,(/nsa_id/),isa_to_isb_id))
      if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_def_var(ncid,'isa_to_isc',nf90_int,(/nsa_id/),isa_to_isc_id))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'isb_to_isc',nf90_int,(/nsb_id/),isb_to_isc_id))
      call ncerr(subr,nf90_def_var(ncid,'norm',ncfloat,(/nc0a_id,nl0_id/),norm_id))
      if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_put_att(ncid,isa_to_isb_id,'_FillValue',msvali))
      if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_put_att(ncid,isa_to_isc_id,'_FillValue',msvali))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_put_att(ncid,isb_to_isc_id,'_FillValue',msvali))
      call ncerr(subr,nf90_put_att(ncid,norm_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      if (sdata%mpi(iproc)%nc1b>0) call ncerr(subr,nf90_put_var(ncid,vbot_id,sdata%mpi(iproc)%vbot))
      call ncerr(subr,nf90_put_var(ncid,nc2b_id,sdata%mpi(iproc)%nc2b))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_put_var(ncid,isb_to_ic2b_id,sdata%mpi(iproc)%isb_to_ic2b))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_put_var(ncid,isb_to_il1_id,sdata%mpi(iproc)%isb_to_il1))
      if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_put_var(ncid,isa_to_isb_id,sdata%mpi(iproc)%isa_to_isb))
      if (sdata%mpi(iproc)%nsa>0) call ncerr(subr,nf90_put_var(ncid,isa_to_isc_id,sdata%mpi(iproc)%isa_to_isc))
      if (sdata%mpi(iproc)%nsb>0) call ncerr(subr,nf90_put_var(ncid,isb_to_isc_id,sdata%mpi(iproc)%isb_to_isc))
      call ncerr(subr,nf90_put_var(ncid,norm_id,sdata%mpi(iproc)%norm%vala))

      ! Write communications
      if (sdata%nproc>1) then
         call com_write(ncid,sdata%mpi(iproc)%AB)
         call com_write(ncid,sdata%mpi(iproc)%AC)
      end if

      ! Write linear operators
      call linop_write(ncid,sdata%mpi(iproc)%c)
      call linop_write(ncid,sdata%mpi(iproc)%h)
      call linop_write(ncid,sdata%mpi(iproc)%v)
      call linop_write(ncid,sdata%mpi(iproc)%s)

      ! Close file
      call ncerr(subr,nf90_close(ncid))
   end do

   ! Create summary file
   write(mpicomchar,'(i1)') sdata%mpicom
   write(nprocchar,'(i4.4)') sdata%nproc
   filename = trim(nam%prefix)//'_mpi-'//mpicomchar//'_'//nprocchar//'_summary.nc'
   call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call namncwrite(ncid)

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,'nc0',sdata%nc0,nc0_id))
   call ncerr(subr,nf90_def_dim(ncid,'nproc1',sdata%nproc,nproc1_id))
   call ncerr(subr,nf90_def_dim(ncid,'nproc2',sdata%nproc,nproc2_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
   call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
   call ncerr(subr,nf90_def_var(ncid,'ic0_to_iproc',nf90_int,(/nc0_id/),ic0_to_iproc_id))
   call ncerr(subr,nf90_def_var(ncid,'halo',nf90_int,(/nc0_id/),halo_id))
   call ncerr(subr,nf90_def_var(ncid,'nc0a_per_proc',nf90_int,(/nproc1_id/),nc0a_per_proc_id))
   call ncerr(subr,nf90_def_var(ncid,'nsa_per_proc',nf90_int,(/nproc1_id/),nsa_per_proc_id))
   call ncerr(subr,nf90_def_var(ncid,'c_n_s',nf90_int,(/nproc1_id/),c_n_s_id))
   call ncerr(subr,nf90_def_var(ncid,'h_n_s',nf90_int,(/nproc1_id/),h_n_s_id))
   call ncerr(subr,nf90_def_var(ncid,'s_n_s',nf90_int,(/nproc1_id/),s_n_s_id))
   call ncerr(subr,nf90_def_var(ncid,'AB_jhalocounts',nf90_int,(/nproc1_id,nproc2_id/),AB_jhalocounts_id))
   call ncerr(subr,nf90_def_var(ncid,'AB_jexclcounts',nf90_int,(/nproc1_id,nproc2_id/),AB_jexclcounts_id))
   call ncerr(subr,nf90_def_var(ncid,'AC_jhalocounts',nf90_int,(/nproc1_id,nproc2_id/),AC_jhalocounts_id))
   call ncerr(subr,nf90_def_var(ncid,'AC_jexclcounts',nf90_int,(/nproc1_id,nproc2_id/),AC_jexclcounts_id))
   call ncerr(subr,nf90_enddef(ncid))

   ! Write variables
   call ncerr(subr,nf90_put_var(ncid,lon_id,sdata%lon*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,lat_id,sdata%lat*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,ic0_to_iproc_id,sdata%ic0_to_iproc))
   if (sdata%nproc>1) call ncerr(subr,nf90_put_var(ncid,halo_id,sdata%halo))
   do iproc=1,sdata%nproc
      call ncerr(subr,nf90_put_var(ncid,nc0a_per_proc_id,sdata%mpi(iproc)%nc0a,(/iproc/)))
      call ncerr(subr,nf90_put_var(ncid,nsa_per_proc_id,sdata%mpi(iproc)%nsa,(/iproc/)))
      call ncerr(subr,nf90_put_var(ncid,c_n_s_id,sdata%mpi(iproc)%c%n_s,(/iproc/)))
      h_sum = 0
      do il0i=1,sdata%mpi(iproc)%nl0i
         h_sum = h_sum+sdata%mpi(iproc)%h(il0i)%n_s
      end do
      call ncerr(subr,nf90_put_var(ncid,h_n_s_id,h_sum,(/iproc/)))
      s_sum = 0
      do il1=1,sdata%mpi(iproc)%nl1
         s_sum = s_sum+sdata%mpi(iproc)%s(il1)%n_s
      end do
      call ncerr(subr,nf90_put_var(ncid,s_n_s_id,s_sum,(/iproc/)))
      call ncerr(subr,nf90_put_var(ncid,AB_jhalocounts_id,sdata%mpi(iproc)%AB%jhalocounts,(/iproc,1/),(/1,sdata%nproc/)))
      call ncerr(subr,nf90_put_var(ncid,AB_jexclcounts_id,sdata%mpi(iproc)%AB%jexclcounts,(/iproc,1/),(/1,sdata%nproc/)))
      call ncerr(subr,nf90_put_var(ncid,AC_jhalocounts_id,sdata%mpi(iproc)%AC%jhalocounts,(/iproc,1/),(/1,sdata%nproc/)))
      call ncerr(subr,nf90_put_var(ncid,AC_jexclcounts_id,sdata%mpi(iproc)%AC%jexclcounts,(/iproc,1/),(/1,sdata%nproc/)))
   end do

   ! Close summary file
   call ncerr(subr,nf90_close(ncid))
end if

end subroutine sdata_write_mpi

end module type_sdata
