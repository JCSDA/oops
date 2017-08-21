!----------------------------------------------------------------------
! Module: type_ndata
!> Purpose: sampling data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_ndata

use module_namelist, only: nam,namncwrite
use netcdf
use tools_const, only: rad2deg
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr
use tools_nc, only: ncerr,ncfloat
use type_com, only: comtype,com_dealloc,com_read,com_write
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_read,linop_write
use type_mpl, only: mpl
use type_randgen, only: randgentype,create_randgen,delete_randgen

implicit none

! Sampling data derived type
type ndatatype
   ! Vector sizes
   integer :: nlon                         !< Longitude size
   integer :: nlat                         !< Latitude size
   integer :: nlev                         !< Number of levels
   logical,allocatable :: rgmask(:,:)      !< Reduced Gaussian grid mask
   real(kind_real),allocatable :: area(:)  !< Domain area

   ! Vector coordinates
   real(kind_real),allocatable :: lon(:)              !< Cells longitude
   real(kind_real),allocatable :: lat(:)              !< Cells latitude
   logical,allocatable :: mask(:,:)        !< Cells mask
   real(kind_real),allocatable :: vunit(:)            !< Vertical unit
   integer :: nfor                         !< Number of forced sampling points
   integer,allocatable :: ifor(:)          !< Forced sampling points indices

   ! Sampling properties
   logical,allocatable :: llev(:)          !< Vertical interpolation key
   integer,allocatable :: net_nnb(:)      !< Number of neighbors on the full grid
   integer,allocatable :: net_inb(:,:)    !< Neighbors indices on the full grid
   real(kind_real),allocatable :: net_dnb(:,:)       !< Neighbors distances on the full grid

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
   real(kind_real),allocatable :: norm(:,:) !< Normalization factor
   real(kind_real),allocatable :: norm_sqrt(:) !< Internal normalization factor for the square-root formulation

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
   integer,allocatable :: ic0_to_iproc(:)    !< Subset Sc0 to local task
   integer,allocatable :: ic0_to_ic0a(:)     !< Subset Sc0, global to halo A
   integer :: nc0amax                        !< Maximum size of halo A

   ! Illustration
   integer,allocatable :: halo(:)            !< Halo points for illustration
end type ndatatype

! Local sampling data derived type
type ndataloctype
   ! Number of points
   integer :: nc0a                       !< Number of points in subset Sc0 on halo A
   integer :: nl0                        !< Number of levels in subset Sl0
   integer :: nc1b                       !< Number of points in subset Sc1 on halo B
   integer :: nl1                        !< Number of levels in subset Sl1
   integer :: nl0i                       !< Number of independent levels
   integer,allocatable :: vbot(:)        !< Bottom level
   integer,allocatable :: nc2b(:)        !< Number of points in subset Sc2 on halo B
   integer :: nsa                        !< Number of subgrid nodes on halo A
   integer :: nsb                        !< Number of subgrid nodes on halo B
   integer :: nsc                        !< Number of subgrid nodes on halo C

   ! Inter-halo conversions
   integer,allocatable :: isa_to_isb(:)  !< Subgrid, halo A to halo B
   integer,allocatable :: isa_to_isc(:)  !< Subgrid, halo A to halo C
   integer,allocatable :: isb_to_isc(:)  !< Subgrid, halo B to halo B

   ! Linear operators
   type(linoptype) :: c                  !< Convolution
   type(linoptype),allocatable :: h(:)   !< Horizontal interpolation
   type(linoptype),allocatable :: v(:)   !< Vertical interpolation
   type(linoptype),allocatable :: s(:)   !< Subsample interpolation

   ! Copy conversions
   integer,allocatable :: isb_to_ic2b(:) !< Subgrid to subset Sc2 on halo B
   integer,allocatable :: isb_to_il1(:)  !< Subgrid to subset Sl1 on halo B

   ! Normalization
   real(kind_real),allocatable :: norm(:,:) !< Normalization factor
   real(kind_real),allocatable :: norm_sqrt(:) !< Internal normalization factor for the square-root formulation

   ! Communications
   type(comtype) :: AB                   !< Communication between halos A and B
   type(comtype) :: AC                   !< Communication between halos A and C
end type ndataloctype

private
public :: ndatatype,ndataloctype
public :: ndata_alloc,ndataloc_dealloc, ndataloc_copy, &
 & ndata_read_param,ndata_read_local,ndata_read_mpi, &
 & ndata_write_param,ndata_write_mpi,ndata_write_mpi_summary

contains

!----------------------------------------------------------------------
! Subroutine: ndata_alloc
!> Purpose: ndata object allocation for grid parameters
!----------------------------------------------------------------------
subroutine ndata_alloc(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Allocation
allocate(ndata%lon(ndata%nc0))
allocate(ndata%lat(ndata%nc0))
allocate(ndata%area(ndata%nl0))
allocate(ndata%mask(ndata%nc0,ndata%nl0))
allocate(ndata%vunit(ndata%nl0))

! Initialization
call msr(ndata%lon)
call msr(ndata%lat)
ndata%mask = .false.
call msr(ndata%vunit)
ndata%rng = create_randgen()

end subroutine ndata_alloc

!----------------------------------------------------------------------
! Subroutine: ndataloc_dealloc
!> Purpose: ndataloc object deallocation
!----------------------------------------------------------------------
subroutine ndataloc_dealloc(ndataloc)

implicit none

! Passed variables
type(ndataloctype),intent(inout) :: ndataloc !< Sampling data, local

! Local variables
integer :: il0i,il1

! Release memory
if (allocated(ndataloc%nc2b)) deallocate(ndataloc%nc2b)
if (allocated(ndataloc%isa_to_isb)) deallocate(ndataloc%isa_to_isb)
if (allocated(ndataloc%isa_to_isc)) deallocate(ndataloc%isa_to_isc)
if (allocated(ndataloc%isb_to_isc)) deallocate(ndataloc%isb_to_isc)
if (allocated(ndataloc%vbot)) deallocate(ndataloc%vbot)
call linop_dealloc(ndataloc%c)
if (allocated(ndataloc%h)) then
   do il0i=1,ndataloc%nl0i
      call linop_dealloc(ndataloc%h(il0i))
   end do
   deallocate(ndataloc%h)
end if
if (allocated(ndataloc%v)) then
   do il0i=1,ndataloc%nl0i
      call linop_dealloc(ndataloc%v(il0i))
   end do
   deallocate(ndataloc%v)
end if
if (allocated(ndataloc%s)) then
   do il1=1,ndataloc%nl1
      call linop_dealloc(ndataloc%s(il1))
   end do
   deallocate(ndataloc%s)
end if
if (allocated(ndataloc%isb_to_ic2b)) deallocate(ndataloc%isb_to_ic2b)
if (allocated(ndataloc%isb_to_il1)) deallocate(ndataloc%isb_to_il1)
if (allocated(ndataloc%norm)) deallocate(ndataloc%norm)
if (allocated(ndataloc%norm_sqrt)) deallocate(ndataloc%norm_sqrt)
call com_dealloc(ndataloc%AB)
call com_dealloc(ndataloc%AC)

end subroutine ndataloc_dealloc

!----------------------------------------------------------------------
! Subroutine: ndataloc_copy
!> Purpose: linear operator copy
!----------------------------------------------------------------------
subroutine ndataloc_copy(ndataloc_in,ndataloc_out)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc_in     !< Input linear operator
type(ndataloctype),intent(inout) :: ndataloc_out !< Output linear operator

! Local variables
integer :: il0i,il1

! Copy attributes
ndataloc_out%nc0a = ndataloc_in%nc0a
ndataloc_out%nl0 = ndataloc_in%nl0
ndataloc_out%nc1b = ndataloc_in%nc1b
ndataloc_out%nl1 = ndataloc_in%nl1
ndataloc_out%nl0i = ndataloc_in%nl0i
ndataloc_out%nsa = ndataloc_in%nsa
ndataloc_out%nsb = ndataloc_in%nsb
ndataloc_out%nsc = ndataloc_in%nsc

! Deallocation
call ndataloc_dealloc(ndataloc_out)

! Allocation
allocate(ndataloc_out%vbot(ndataloc_out%nc1b))
allocate(ndataloc_out%nc2b(ndataloc_out%nl1))
allocate(ndataloc_out%isa_to_isb(ndataloc_out%nsa))
allocate(ndataloc_out%isa_to_isc(ndataloc_out%nsa))
allocate(ndataloc_out%isb_to_isc(ndataloc_out%nsa))
allocate(ndataloc_out%h(ndataloc_out%nl0i))
allocate(ndataloc_out%v(ndataloc_out%nl0i))
allocate(ndataloc_out%s(ndataloc_out%nl1))
allocate(ndataloc_out%isb_to_ic2b(ndataloc_out%nsb))
allocate(ndataloc_out%isb_to_il1(ndataloc_out%nsb))
allocate(ndataloc_out%norm(ndataloc_out%nc0a,ndataloc_out%nl0))
if (nam%lsqrt) allocate(ndataloc_out%norm_sqrt(ndataloc_out%nsb))

! Copy
ndataloc_out%vbot = ndataloc_in%vbot
ndataloc_out%nc2b = ndataloc_in%nc2b
ndataloc_out%isa_to_isb = ndataloc_in%isa_to_isb
ndataloc_out%isa_to_isc = ndataloc_in%isa_to_isc
ndataloc_out%isb_to_isc = ndataloc_in%isb_to_isc
ndataloc_out%c = ndataloc_in%c
do il0i=1,ndataloc_out%nl0i
   call linop_copy(ndataloc_in%h(il0i),ndataloc_out%h(il0i))
end do
do il0i=1,ndataloc_out%nl0i
   call linop_copy(ndataloc_in%v(il0i),ndataloc_out%v(il0i))
end do
do il1=1,ndataloc_out%nl1
   call linop_copy(ndataloc_in%s(il1),ndataloc_out%s(il1))
end do
ndataloc_out%isb_to_ic2b = ndataloc_in%isb_to_ic2b
ndataloc_out%isb_to_il1 = ndataloc_in%isb_to_il1
ndataloc_out%norm = ndataloc_in%norm
if (nam%lsqrt) ndataloc_out%norm_sqrt = ndataloc_in%norm_sqrt
ndataloc_out%AB = ndataloc_in%AB
ndataloc_out%AC = ndataloc_in%AC

end subroutine ndataloc_copy

!----------------------------------------------------------------------
! Subroutine: ndata_read_param
!> Purpose: read ndata object
!----------------------------------------------------------------------
subroutine ndata_read_param(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ncid
integer :: nc0_id,nl0_id,nc1_id,nl1_id,ns_id
integer :: vbot_id,nc2_id,is_to_ic1_id,is_to_il1_id,is_to_ic2_id,ic0il0_to_is_id,ic2il1_to_ic0_id
integer :: ic2il1_to_is_id,ic1_to_ic0_id,il1_to_il0_id,ic0_to_ic1_id,il0_to_il1_id,ic2il1_to_ic1_id
integer :: ic1il1_to_is_id,norm_id,norm_sqrt_id
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_read_param'

! Open file and get dimensions
filename = trim(nam%prefix)//'_param.nc'
call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'nc0',nc0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc0_id,len=ndata%nc0))
call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=ndata%nl0))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'nl0i',ndata%nl0i))
call ncerr(subr,nf90_inq_dimid(ncid,'nc1',nc1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc1_id,len=ndata%nc1))
call ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=ndata%nl1))
call ncerr(subr,nf90_inq_dimid(ncid,'ns',ns_id))
call ncerr(subr,nf90_inquire_dimension(ncid,ns_id,len=ndata%ns))

! Allocation
allocate(ndata%vbot(ndata%nc1))
allocate(ndata%nc2(ndata%nl1))
allocate(ndata%is_to_ic1(ndata%ns))
allocate(ndata%is_to_il1(ndata%ns))
allocate(ndata%is_to_ic2(ndata%ns))
allocate(ndata%ic0il0_to_is(ndata%nc0,ndata%nl0))
allocate(ndata%ic2il1_to_ic0(ndata%nc1,ndata%nl1))
allocate(ndata%ic2il1_to_is(ndata%nc1,ndata%nl1))
allocate(ndata%ic1_to_ic0(ndata%nc1))
allocate(ndata%il1_to_il0(ndata%nl1))
allocate(ndata%ic0_to_ic1(ndata%nc0))
allocate(ndata%il0_to_il1(ndata%nl0))
allocate(ndata%ic2il1_to_ic1(ndata%nc1,ndata%nl1))
allocate(ndata%ic1il1_to_is(ndata%nc1,ndata%nl1))
allocate(ndata%norm(ndata%nc0,ndata%nl0))
if (nam%lsqrt) allocate(ndata%norm_sqrt(ndata%ns))

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
call ncerr(subr,nf90_get_var(ncid,vbot_id,ndata%vbot))
call ncerr(subr,nf90_get_var(ncid,nc2_id,ndata%nc2))
call ncerr(subr,nf90_get_var(ncid,is_to_ic1_id,ndata%is_to_ic1))
call ncerr(subr,nf90_get_var(ncid,is_to_il1_id,ndata%is_to_il1))
call ncerr(subr,nf90_get_var(ncid,is_to_ic2_id,ndata%is_to_ic2))
call ncerr(subr,nf90_get_var(ncid,ic0il0_to_is_id,ndata%ic0il0_to_is))
call ncerr(subr,nf90_get_var(ncid,ic2il1_to_ic0_id,ndata%ic2il1_to_ic0))
call ncerr(subr,nf90_get_var(ncid,ic2il1_to_is_id,ndata%ic2il1_to_is))
call ncerr(subr,nf90_get_var(ncid,ic1_to_ic0_id,ndata%ic1_to_ic0))
call ncerr(subr,nf90_get_var(ncid,il1_to_il0_id,ndata%il1_to_il0))
call ncerr(subr,nf90_get_var(ncid,ic0_to_ic1_id,ndata%ic0_to_ic1))
call ncerr(subr,nf90_get_var(ncid,il0_to_il1_id,ndata%il0_to_il1))
call ncerr(subr,nf90_get_var(ncid,ic2il1_to_ic1_id,ndata%ic2il1_to_ic1))
call ncerr(subr,nf90_get_var(ncid,ic1il1_to_is_id,ndata%ic1il1_to_is))
call ncerr(subr,nf90_get_var(ncid,norm_id,ndata%norm))
if (nam%lsqrt) call ncerr(subr,nf90_get_var(ncid,norm_sqrt_id,ndata%norm_sqrt))

! Read linear operators
call linop_read(ncid,'c',ndata%c)
call linop_read(ncid,'h',ndata%h)
call linop_read(ncid,'v',ndata%v)
call linop_read(ncid,'s',ndata%s)

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine ndata_read_param

!----------------------------------------------------------------------
! Subroutine: ndata_read_local
!> Purpose: read ndata object
!----------------------------------------------------------------------
subroutine ndata_read_local(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ic0,info,iproc,ic0a
integer :: ncid,ic0_to_iproc_id,ic0_to_ic0a_id
character(len=4) :: nprocchar
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_read_local'

if (.not.allocated(ndata%ic0_to_iproc)) then
   ! Allocation
   allocate(ndata%ic0_to_iproc(ndata%nc0))
   allocate(ndata%ic0_to_ic0a(ndata%nc0))
   
   if (nam%nproc==1) then
      ! All points on a single processor
      ndata%ic0_to_iproc = 1
      do ic0=1,ndata%nc0
         ndata%ic0_to_ic0a(ic0) = ic0
      end do
      ndata%nc0amax = ndata%nc0
   elseif (nam%nproc>1) then
      ! Open file
      write(nprocchar,'(i4.4)') nam%nproc
      filename = trim(nam%prefix)//'_distribution_'//nprocchar//'.nc'
      info = nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid)
   
      if (info==nf90_noerr) then
         ! Read data and close file
         call ncerr(subr,nf90_inq_varid(ncid,'ic0_to_iproc',ic0_to_iproc_id))
         call ncerr(subr,nf90_inq_varid(ncid,'ic0_to_ic0a',ic0_to_ic0a_id))
         call ncerr(subr,nf90_get_var(ncid,ic0_to_iproc_id,ndata%ic0_to_iproc))
         call ncerr(subr,nf90_get_var(ncid,ic0_to_ic0a_id,ndata%ic0_to_ic0a))
         call ncerr(subr,nf90_close(ncid))
         ndata%nc0amax = maxval(ndata%ic0_to_ic0a)
      else
         ! Generate a distribution (use METIS one day?)
         ndata%nc0amax = ndata%nc0/nam%nproc
         if (ndata%nc0amax*nam%nproc<ndata%nc0) ndata%nc0amax = ndata%nc0amax+1
         iproc = 1
         ic0a = 1
         do ic0=1,ndata%nc0
            ndata%ic0_to_iproc(ic0) = iproc
            ndata%ic0_to_ic0a(ic0) = ic0a
            ic0a = ic0a+1
            if (ic0a>ndata%nc0amax) then
               ! Change proc
               iproc = iproc+1
               ic0a = 1
            end if
         end do
      end if
   end if
end if

! Check
if (maxval(ndata%ic0_to_iproc)>nam%nproc) call msgerror('wrong distribution')

end subroutine ndata_read_local

!----------------------------------------------------------------------
! Subroutine: ndata_read_mpi
!> Purpose: read ndata object
!----------------------------------------------------------------------
subroutine ndata_read_mpi(ndataloc)

implicit none

! Passed variables
type(ndataloctype),intent(inout) :: ndataloc !< Sampling data, local

! Local variables
integer :: ncid,info
integer :: nc0a_id,nl0_id,nc1b_id,nl1_id,nsa_id,nsb_id
integer :: vbot_id,nc2b_id,isb_to_ic2b_id,isb_to_il1_id
integer :: isa_to_isb_id,isa_to_isc_id,isb_to_isc_id
integer :: norm_id
character(len=1) :: mpicomchar
character(len=4) :: nprocchar,myprocchar
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_read_mpi'

! Open file and get dimensions
write(mpicomchar,'(i1)') nam%mpicom
write(nprocchar,'(i4.4)') nam%nproc
write(myprocchar,'(i4.4)') mpl%myproc
filename = trim(nam%prefix)//'_mpi-'//mpicomchar//'_'//nprocchar//'-'//myprocchar//'.nc'
call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'nc0a',nc0a_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc0a_id,len=ndataloc%nc0a))
call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=ndataloc%nl0))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'nl0i',ndataloc%nl0i))
info = nf90_inq_dimid(ncid,'nc1b',nc1b_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,nc1b_id,len=ndataloc%nc1b))
else
   ndataloc%nc1b = 0
end if
call ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=ndataloc%nl1))
info = nf90_inq_dimid(ncid,'nsa',nsa_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,nsa_id,len=ndataloc%nsa))
else
   ndataloc%nsa = 0
end if
info = nf90_inq_dimid(ncid,'nsb',nsb_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,nsb_id,len=ndataloc%nsb))
else
   ndataloc%nsb = 0
end if
call ncerr(subr,nf90_get_att(ncid,nf90_global,'nsc',ndataloc%nsc))

! Allocation
if (ndataloc%nc1b>0) allocate(ndataloc%vbot(ndataloc%nc1b))
allocate(ndataloc%nc2b(ndataloc%nl1))
if (ndataloc%nsb>0) allocate(ndataloc%isb_to_ic2b(ndataloc%nsb))
if (ndataloc%nsb>0) allocate(ndataloc%isb_to_il1(ndataloc%nsb))
if (ndataloc%nsa>0) allocate(ndataloc%isa_to_isb(ndataloc%nsa))
if (ndataloc%nsa>0) allocate(ndataloc%isa_to_isc(ndataloc%nsa))
if (ndataloc%nsb>0) allocate(ndataloc%isb_to_isc(ndataloc%nsb))
allocate(ndataloc%norm(ndataloc%nc0a,ndataloc%nl0))

! Read data
if (ndataloc%nc1b>0) call ncerr(subr,nf90_inq_varid(ncid,'vbot',vbot_id))
call ncerr(subr,nf90_inq_varid(ncid,'nc2b',nc2b_id))
if (ndataloc%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'isb_to_ic2b',isb_to_ic2b_id))
if (ndataloc%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'isb_to_il1',isb_to_il1_id))
if (ndataloc%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'isa_to_isb',isa_to_isb_id))
if (ndataloc%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'isa_to_isc',isa_to_isc_id))
if (ndataloc%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'isb_to_isc',isb_to_isc_id))
call ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
if (ndataloc%nc1b>0) call ncerr(subr,nf90_get_var(ncid,vbot_id,ndataloc%vbot))
call ncerr(subr,nf90_get_var(ncid,nc2b_id,ndataloc%nc2b))
if (ndataloc%nsb>0) call ncerr(subr,nf90_get_var(ncid,isb_to_ic2b_id,ndataloc%isb_to_ic2b))
if (ndataloc%nsb>0) call ncerr(subr,nf90_get_var(ncid,isb_to_il1_id,ndataloc%isb_to_il1))
if (ndataloc%nsa>0) call ncerr(subr,nf90_get_var(ncid,isa_to_isb_id,ndataloc%isa_to_isb))
if (ndataloc%nsa>0) call ncerr(subr,nf90_get_var(ncid,isa_to_isc_id,ndataloc%isa_to_isc))
if (ndataloc%nsb>0) call ncerr(subr,nf90_get_var(ncid,isb_to_isc_id,ndataloc%isb_to_isc))
call ncerr(subr,nf90_get_var(ncid,norm_id,ndataloc%norm))

! Read communications
call com_read(ncid,'AB',ndataloc%AB)
call com_read(ncid,'AC',ndataloc%AC)

! Read linear operators
call linop_read(ncid,'c',ndataloc%c)
call linop_read(ncid,'h',ndataloc%h)
call linop_read(ncid,'v',ndataloc%v)
call linop_read(ncid,'s',ndataloc%s)

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine ndata_read_mpi

!----------------------------------------------------------------------
! Subroutine: ndata_write_param
!> Purpose: write ndata object
!----------------------------------------------------------------------
subroutine ndata_write_param(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ncid
integer :: nc0_id,nc1_id,nl0_id,nl1_id,ns_id
integer :: vbot_id,nc2_id,is_to_ic1_id,is_to_il1_id,is_to_ic2_id,ic0il0_to_is_id,ic2il1_to_ic0_id
integer :: ic2il1_to_is_id,ic1_to_ic0_id,il1_to_il0_id,ic0_to_ic1_id,il0_to_il1_id,ic2il1_to_ic1_id
integer :: ic1il1_to_is_id,norm_id,norm_sqrt_id
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_write_param'

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Create file
filename = trim(nam%prefix)//'_param.nc'
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call namncwrite(ncid)

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nc0',ndata%nc0,nc0_id))
call ncerr(subr,nf90_def_dim(ncid,'nc1',ndata%nc1,nc1_id))
call ncerr(subr,nf90_def_dim(ncid,'nl0',ndata%nl0,nl0_id))
call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',ndata%nl0i))
call ncerr(subr,nf90_def_dim(ncid,'nl1',ndata%nl1,nl1_id))
call ncerr(subr,nf90_def_dim(ncid,'ns',ndata%ns,ns_id))

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
call ncerr(subr,nf90_put_var(ncid,vbot_id,ndata%vbot))
call ncerr(subr,nf90_put_var(ncid,nc2_id,ndata%nc2))
call ncerr(subr,nf90_put_var(ncid,is_to_ic1_id,ndata%is_to_ic1))
call ncerr(subr,nf90_put_var(ncid,is_to_il1_id,ndata%is_to_il1))
call ncerr(subr,nf90_put_var(ncid,is_to_ic2_id,ndata%is_to_ic2))
call ncerr(subr,nf90_put_var(ncid,ic0il0_to_is_id,ndata%ic0il0_to_is))
call ncerr(subr,nf90_put_var(ncid,ic2il1_to_ic0_id,ndata%ic2il1_to_ic0))
call ncerr(subr,nf90_put_var(ncid,ic2il1_to_is_id,ndata%ic2il1_to_is))
call ncerr(subr,nf90_put_var(ncid,ic1_to_ic0_id,ndata%ic1_to_ic0))
call ncerr(subr,nf90_put_var(ncid,il1_to_il0_id,ndata%il1_to_il0))
call ncerr(subr,nf90_put_var(ncid,ic0_to_ic1_id,ndata%ic0_to_ic1))
call ncerr(subr,nf90_put_var(ncid,il0_to_il1_id,ndata%il0_to_il1))
call ncerr(subr,nf90_put_var(ncid,ic2il1_to_ic1_id,ndata%ic2il1_to_ic1))
call ncerr(subr,nf90_put_var(ncid,ic1il1_to_is_id,ndata%ic1il1_to_is))
call ncerr(subr,nf90_put_var(ncid,norm_id,ndata%norm))
if (nam%lsqrt) call ncerr(subr,nf90_put_var(ncid,norm_sqrt_id,ndata%norm_sqrt))

! Write linear operators
call linop_write(ncid,ndata%c)
call linop_write(ncid,ndata%h)
call linop_write(ncid,ndata%v)
call linop_write(ncid,ndata%s)

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine ndata_write_param

!----------------------------------------------------------------------
! Subroutine: ndata_write_mpi
!> Purpose: write ndata object
!----------------------------------------------------------------------
subroutine ndata_write_mpi(ndataloc)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc(nam%nproc) !< Sampling data, local

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
character(len=1024) :: subr = 'ndata_write_mpi'

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Filename suffix
write(mpicomchar,'(i1)') nam%mpicom
write(nprocchar,'(i4.4)') nam%nproc

do iproc=1,nam%nproc
   ! Create file
   write(iprocchar,'(i4.4)') iproc
   filename = trim(nam%prefix)//'_mpi-'//mpicomchar//'_'//nprocchar//'-'//iprocchar//'.nc'
   call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call namncwrite(ncid)

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,'nc0a',ndataloc(iproc)%nc0a,nc0a_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl0',ndataloc(iproc)%nl0,nl0_id))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',ndataloc(iproc)%nl0i))
   if (ndataloc(iproc)%nc1b>0) call ncerr(subr,nf90_def_dim(ncid,'nc1b',ndataloc(iproc)%nc1b,nc1b_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl1',ndataloc(iproc)%nl1,nl1_id))
   if (ndataloc(iproc)%nsa>0) call ncerr(subr,nf90_def_dim(ncid,'nsa',ndataloc(iproc)%nsa,nsa_id))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_def_dim(ncid,'nsb',ndataloc(iproc)%nsb,nsb_id))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'nsc',ndataloc(iproc)%nsc))

   ! Define variables
   if (ndataloc(iproc)%nc1b>0) call ncerr(subr,nf90_def_var(ncid,'vbot',nf90_int,(/nc1b_id/),vbot_id))
   call ncerr(subr,nf90_def_var(ncid,'nc2b',nf90_int,(/nl1_id/),nc2b_id))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'isb_to_ic2b',nf90_int,(/nsb_id/),isb_to_ic2b_id))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'isb_to_il1',nf90_int,(/nsb_id/),isb_to_il1_id))
   if (ndataloc(iproc)%nsa>0) call ncerr(subr,nf90_def_var(ncid,'isa_to_isb',nf90_int,(/nsa_id/),isa_to_isb_id))
   if (ndataloc(iproc)%nsa>0) call ncerr(subr,nf90_def_var(ncid,'isa_to_isc',nf90_int,(/nsa_id/),isa_to_isc_id))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'isb_to_isc',nf90_int,(/nsb_id/),isb_to_isc_id))
   call ncerr(subr,nf90_def_var(ncid,'norm',ncfloat,(/nc0a_id,nl0_id/),norm_id))
   if (ndataloc(iproc)%nsa>0) call ncerr(subr,nf90_put_att(ncid,isa_to_isb_id,'_FillValue',msvali))
   if (ndataloc(iproc)%nsa>0) call ncerr(subr,nf90_put_att(ncid,isa_to_isc_id,'_FillValue',msvali))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_put_att(ncid,isb_to_isc_id,'_FillValue',msvali))
   call ncerr(subr,nf90_put_att(ncid,norm_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))

   ! Write variables
   if (ndataloc(iproc)%nc1b>0) call ncerr(subr,nf90_put_var(ncid,vbot_id,ndataloc(iproc)%vbot))
   call ncerr(subr,nf90_put_var(ncid,nc2b_id,ndataloc(iproc)%nc2b))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_put_var(ncid,isb_to_ic2b_id,ndataloc(iproc)%isb_to_ic2b))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_put_var(ncid,isb_to_il1_id,ndataloc(iproc)%isb_to_il1))
   if (ndataloc(iproc)%nsa>0) call ncerr(subr,nf90_put_var(ncid,isa_to_isb_id,ndataloc(iproc)%isa_to_isb))
   if (ndataloc(iproc)%nsa>0) call ncerr(subr,nf90_put_var(ncid,isa_to_isc_id,ndataloc(iproc)%isa_to_isc))
   if (ndataloc(iproc)%nsb>0) call ncerr(subr,nf90_put_var(ncid,isb_to_isc_id,ndataloc(iproc)%isb_to_isc))
   call ncerr(subr,nf90_put_var(ncid,norm_id,ndataloc(iproc)%norm))

   ! Write communications
   call com_write(ncid,ndataloc(iproc)%AB)
   call com_write(ncid,ndataloc(iproc)%AC)

   ! Write linear operators
   call linop_write(ncid,ndataloc(iproc)%c)
   call linop_write(ncid,ndataloc(iproc)%h)
   call linop_write(ncid,ndataloc(iproc)%v)
   call linop_write(ncid,ndataloc(iproc)%s)

   ! Close file
   call ncerr(subr,nf90_close(ncid))
end do

end subroutine ndata_write_mpi

!----------------------------------------------------------------------
! Subroutine: ndata_write_mpi_summary
!> Purpose: write ndata object
!----------------------------------------------------------------------
subroutine ndata_write_mpi_summary(ndata,ndataloc)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata                  !< Sampling data
type(ndataloctype),intent(in) :: ndataloc(nam%nproc) !< Sampling data, local

! Local variables
integer :: iproc,il0i,il1,h_sum,s_sum
integer :: ncid
integer :: nc0a_id,nl0_id,nc1b_id,nl1_id,nsa_id,nsb_id
integer :: vbot_id,nc2b_id,isb_to_ic2b_id,isb_to_il1_id
integer :: c_n_s_id,h_n_s_id,s_n_s_id
integer :: isa_to_isb_id,isa_to_isc_id,isb_to_isc_id
integer :: AB_jhalocounts_id,AB_jexclcounts_id
integer :: AC_jhalocounts_id,AC_jexclcounts_id
integer :: nc0_id,nproc1_id,nproc2_id
integer :: lon_id,lat_id,ic0_to_iproc_id,halo_id,nc0a_per_proc_id,nsa_per_proc_id
character(len=1) :: mpicomchar
character(len=4) :: nprocchar,iprocchar
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_write_mpi_summary'

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Filename suffix
write(mpicomchar,'(i1)') nam%mpicom
write(nprocchar,'(i4.4)') nam%nproc

! Create summary file
write(mpicomchar,'(i1)') nam%mpicom
write(nprocchar,'(i4.4)') nam%nproc
filename = trim(nam%prefix)//'_mpi-'//mpicomchar//'_'//nprocchar//'_summary.nc'
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call namncwrite(ncid)

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nc0',ndata%nc0,nc0_id))
call ncerr(subr,nf90_def_dim(ncid,'nproc1',nam%nproc,nproc1_id))
call ncerr(subr,nf90_def_dim(ncid,'nproc2',nam%nproc,nproc2_id))

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
call ncerr(subr,nf90_put_var(ncid,lon_id,ndata%lon*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_id,ndata%lat*rad2deg))
call ncerr(subr,nf90_put_var(ncid,ic0_to_iproc_id,ndata%ic0_to_iproc))
if (nam%nproc>1) call ncerr(subr,nf90_put_var(ncid,halo_id,ndata%halo))
do iproc=1,nam%nproc
   call ncerr(subr,nf90_put_var(ncid,nc0a_per_proc_id,ndataloc(iproc)%nc0a,(/iproc/)))
   call ncerr(subr,nf90_put_var(ncid,nsa_per_proc_id,ndataloc(iproc)%nsa,(/iproc/)))
   call ncerr(subr,nf90_put_var(ncid,c_n_s_id,ndataloc(iproc)%c%n_s,(/iproc/)))
   h_sum = 0
   do il0i=1,ndataloc(iproc)%nl0i
      h_sum = h_sum+ndataloc(iproc)%h(il0i)%n_s
   end do
   call ncerr(subr,nf90_put_var(ncid,h_n_s_id,h_sum,(/iproc/)))
   s_sum = 0
   do il1=1,ndataloc(iproc)%nl1
      s_sum = s_sum+ndataloc(iproc)%s(il1)%n_s
   end do
   call ncerr(subr,nf90_put_var(ncid,s_n_s_id,s_sum,(/iproc/)))
   call ncerr(subr,nf90_put_var(ncid,AB_jhalocounts_id,ndataloc(iproc)%AB%jhalocounts,(/iproc,1/),(/1,nam%nproc/)))
   call ncerr(subr,nf90_put_var(ncid,AB_jexclcounts_id,ndataloc(iproc)%AB%jexclcounts,(/iproc,1/),(/1,nam%nproc/)))
   call ncerr(subr,nf90_put_var(ncid,AC_jhalocounts_id,ndataloc(iproc)%AC%jhalocounts,(/iproc,1/),(/1,nam%nproc/)))
   call ncerr(subr,nf90_put_var(ncid,AC_jexclcounts_id,ndataloc(iproc)%AC%jexclcounts,(/iproc,1/),(/1,nam%nproc/)))
end do

! Close summary file
call ncerr(subr,nf90_close(ncid))

end subroutine ndata_write_mpi_summary

end module type_ndata
