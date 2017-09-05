!----------------------------------------------------------------------
! Module: type_geom
!> Purpose: geometry derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_geom

use module_namelist, only: namtype
use netcdf
use tools_const, only: req,sphere_dist,vector_product
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsi
use tools_nc, only: ncerr
use tools_stripack, only: areas,trans,trlist
use type_mesh, only: meshtype,create_mesh
use type_mpl, only: mpl,mpl_recv,mpl_send
use type_randgen, only: rng

implicit none

! Sampling data derived type
type geomtype
   ! Vector sizes
   integer :: nlon                         !< Longitude size
   integer :: nlat                         !< Latitude size
   integer :: nlev                         !< Number of levels
   logical,allocatable :: rgmask(:,:)      !< Reduced Gaussian grid mask
   real(kind_real),allocatable :: area(:)  !< Domain area

   ! Number of points and levels
   integer :: nc0                          !< Number of points in subset Sc0
   integer :: nl0                          !< Number of levels in subset Sl0

   ! Vector coordinates
   real(kind_real),allocatable :: lon(:)              !< Cells longitude
   real(kind_real),allocatable :: lat(:)              !< Cells latitude
   logical,allocatable :: mask(:,:)                !< Cells mask
   real(kind_real),allocatable :: vunit(:)            !< Vertical unit

   ! Mesh
   type(meshtype) :: mesh

   ! Boundary nodes
   integer,allocatable :: nbnd(:)          !< Number of boundary nodes
   real(kind_real),allocatable :: xbnd(:,:,:)         !< Boundary nodes, x-coordinate
   real(kind_real),allocatable :: ybnd(:,:,:)         !< Boundary nodes, y-coordinate
   real(kind_real),allocatable :: zbnd(:,:,:)         !< Boundary nodes, z-coordinate
   real(kind_real),allocatable :: vbnd(:,:,:)         !< Boundary nodes, orthogonal vector

   ! Neighbors
   integer,allocatable :: net_nnb(:)      !< Number of neighbors on the full grid
   integer,allocatable :: net_inb(:,:)    !< Neighbors indices on the full grid
   real(kind_real),allocatable :: net_dnb(:,:)       !< Neighbors distances on the full grid

   ! MPI distribution
   integer,allocatable :: ic0_to_iproc(:)    !< Subset Sc0 to local task
   integer,allocatable :: ic0_to_ic0a(:)     !< Subset Sc0, global to halo A
   integer,allocatable :: iproc_to_nc0a(:)     !< Halo A size for each proc
end type geomtype

private
public :: geomtype
public :: geom_alloc,compute_grid_mesh,geom_read_local,fld_com_gl,fld_com_lg

contains

!----------------------------------------------------------------------
! Subroutine: geom_alloc
!> Purpose: geom object allocation
!----------------------------------------------------------------------
subroutine geom_alloc(geom)

implicit none

! Passed variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Allocation
allocate(geom%lon(geom%nc0))
allocate(geom%lat(geom%nc0))
allocate(geom%area(geom%nl0))
allocate(geom%mask(geom%nc0,geom%nl0))
allocate(geom%vunit(geom%nl0))

! Initialization
call msr(geom%lon)
call msr(geom%lat)
geom%mask = .false.
call msr(geom%vunit)

end subroutine geom_alloc

!----------------------------------------------------------------------
! Subroutine: compute_grid_mesh
!> Purpose: compute grid mesh
!----------------------------------------------------------------------
subroutine compute_grid_mesh(nam,geom)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Create mesh
if ((.not.all(geom%area>0.0)).or.(nam%new_param.and.(nam%mask_check.or.nam%network))) &
 & call create_mesh(rng,geom%nc0,geom%lon,geom%lat,.true.,geom%mesh)

! Compute area
if ((.not.all(geom%area>0.0))) call compute_area(geom)

! Compute mask boundaries
if (nam%new_param.and.nam%mask_check) call compute_mask_boundaries(geom)

! Find grid neighbors
if (nam%new_param.and.nam%network.and..not.allocated(geom%net_nnb)) call find_grid_neighbors(geom)
   
! Compute distances between neighbors
if (nam%new_param.and.nam%network) call compute_grid_neighbors_distances(geom)

end subroutine compute_grid_mesh

!----------------------------------------------------------------------
! Subroutine: compute_area
!> Purpose: compute domain area
!----------------------------------------------------------------------
subroutine compute_area(geom)

implicit none

! Passed variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Local variables
integer :: info,il0,nt,it
integer,allocatable :: ltri(:,:)
real(kind_real) :: area,frac

! Allocation
allocate(ltri(6,2*(geom%mesh%nnr-2)))

! Create triangles list
call trlist(geom%mesh%nnr,geom%mesh%list,geom%mesh%lptr,geom%mesh%lend,6,nt,ltri,info)

! Compute area
geom%area = 0.0
do it=1,nt
   area = areas((/geom%mesh%x(ltri(1,it)),geom%mesh%y(ltri(1,it)),geom%mesh%z(ltri(1,it))/), &
              & (/geom%mesh%x(ltri(2,it)),geom%mesh%y(ltri(2,it)),geom%mesh%z(ltri(2,it))/), &
              & (/geom%mesh%x(ltri(3,it)),geom%mesh%y(ltri(3,it)),geom%mesh%z(ltri(3,it))/))
   do il0=1,geom%nl0
      frac = float(count(geom%mask(geom%mesh%order(ltri(1:3,it)),il0)))/3.0
      geom%area(il0) = geom%area(il0)+frac*area
   end do  
end do

! Release memory
deallocate(ltri)

end subroutine compute_area

!----------------------------------------------------------------------
! Subroutine: compute_mask_boundaries
!> Purpose: compute domain area
!----------------------------------------------------------------------
subroutine compute_mask_boundaries(geom)

implicit none

! Passed variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Local variables
integer :: ic0,jc0,kc0,i,ibnd,il0
integer,allocatable :: ic0_bnd(:,:,:)
real(kind_real) :: latbnd(2),lonbnd(2),v1(3),v2(3)
logical :: init

! Allocation
allocate(geom%nbnd(geom%nl0))
allocate(ic0_bnd(2,geom%mesh%nnr,geom%nl0))

! Find border points
do il0=1,geom%nl0
   geom%nbnd(il0) = 0
   do ic0=1,geom%mesh%nnr
      ! Check mask points only
      if (.not.geom%mask(geom%mesh%order(ic0),il0)) then
         i = geom%mesh%lend(ic0)
         init = .true.
         do while ((i/=geom%mesh%lend(ic0)).or.init)
            jc0 = abs(geom%mesh%list(i))
            kc0 = abs(geom%mesh%list(geom%mesh%lptr(i)))
            if (.not.geom%mask(geom%mesh%order(jc0),il0).and.geom%mask(geom%mesh%order(kc0),il0)) then
               ! Create a new boundary arc
               geom%nbnd(il0) = geom%nbnd(il0)+1
               if (geom%nbnd(il0)>geom%mesh%nnr) call msgerror('too many boundary arcs')
               ic0_bnd(1,geom%nbnd(il0),il0) = geom%mesh%order(ic0)
               ic0_bnd(2,geom%nbnd(il0),il0) = geom%mesh%order(jc0)
            end if
            i = geom%mesh%lptr(i)
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

end subroutine compute_mask_boundaries

!----------------------------------------------------------------------
! Subroutine: find_grid_neighbors
!> Purpose: find full grid neighbors
!----------------------------------------------------------------------
subroutine find_grid_neighbors(geom)

implicit none

! Passed variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Local variables
integer :: ic0,i
logical :: init

! Allocation
allocate(geom%net_nnb(geom%nc0))
   
! Count neighbors
geom%net_nnb = 0
do ic0=1,geom%mesh%nnr
   i = geom%mesh%lend(ic0)
   init = .true.
   do while ((i/=geom%mesh%lend(ic0)).or.init)
      geom%net_nnb(geom%mesh%order(ic0)) = geom%net_nnb(geom%mesh%order(ic0))+1
      i = geom%mesh%lptr(i)
      init = .false.
   end do
end do
   
! Allocation
allocate(geom%net_inb(maxval(geom%net_nnb),geom%nc0))
   
! Find neighbors
geom%net_nnb = 0
do ic0=1,geom%mesh%nnr
   i = geom%mesh%lend(ic0)
   init = .true.
   do while ((i/=geom%mesh%lend(ic0)).or.init)
      geom%net_nnb(geom%mesh%order(ic0)) = geom%net_nnb(geom%mesh%order(ic0))+1
      geom%net_inb(geom%net_nnb(geom%mesh%order(ic0)),geom%mesh%order(ic0)) = geom%mesh%order(abs(geom%mesh%list(i)))
      i = geom%mesh%lptr(i)
      init = .false.
   end do
end do
   
! Copy neighbors for redudant points
do ic0=1,geom%nc0
   if (isnotmsi(geom%mesh%redundant(ic0))) then
      geom%net_nnb(ic0) = geom%net_nnb(geom%mesh%redundant(ic0))
      geom%net_inb(:,ic0) = geom%net_inb(:,geom%mesh%redundant(ic0))
   end if
end do

end subroutine find_grid_neighbors

!----------------------------------------------------------------------
! Subroutine: compute_grid_neighbors_distances
!> Purpose: compute distances between full grid neighbors
!----------------------------------------------------------------------
subroutine compute_grid_neighbors_distances(geom)

implicit none

! Passed variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Local variables
integer :: ic0,i


! Allocation
allocate(geom%net_dnb(maxval(geom%net_nnb),geom%nc0))

! Compute distances 
do ic0=1,geom%nc0
   do i=1,geom%net_nnb(ic0)
      call sphere_dist(geom%lon(ic0),geom%lat(ic0),geom%lon(geom%net_inb(i,ic0)), &
    & geom%lat(geom%net_inb(i,ic0)),geom%net_dnb(i,ic0))
      geom%net_dnb(i,ic0) = (geom%net_dnb(i,ic0)/req)**2
   end do
end do

end subroutine compute_grid_neighbors_distances

!----------------------------------------------------------------------
! Subroutine: geom_read_local
!> Purpose: read local distribution
!----------------------------------------------------------------------
subroutine geom_read_local(nam,geom)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(geomtype),intent(inout) :: geom !< Sampling data

! Local variables
integer :: ic0,info,iproc,ic0a,nc0amax
integer :: ncid,ic0_to_iproc_id,ic0_to_ic0a_id
character(len=4) :: nprocchar
character(len=1024) :: filename
character(len=1024) :: subr = 'geom_read_local'

if (.not.allocated(geom%ic0_to_iproc)) then
   ! Allocation
   allocate(geom%ic0_to_iproc(geom%nc0))
   allocate(geom%ic0_to_ic0a(geom%nc0))

   if (nam%nproc==1) then
      ! All points on a single processor
      geom%ic0_to_iproc = 1
      do ic0=1,geom%nc0
         geom%ic0_to_ic0a(ic0) = ic0
      end do
   elseif (nam%nproc>1) then
      ! Open file
      write(nprocchar,'(i4.4)') nam%nproc
      filename = trim(nam%prefix)//'_distribution_'//nprocchar//'.nc'
      info = nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid)
   
      if (info==nf90_noerr) then
         ! Read data and close file
         call ncerr(subr,nf90_inq_varid(ncid,'ic0_to_iproc',ic0_to_iproc_id))
         call ncerr(subr,nf90_inq_varid(ncid,'ic0_to_ic0a',ic0_to_ic0a_id))
         call ncerr(subr,nf90_get_var(ncid,ic0_to_iproc_id,geom%ic0_to_iproc))
         call ncerr(subr,nf90_get_var(ncid,ic0_to_ic0a_id,geom%ic0_to_ic0a))
         call ncerr(subr,nf90_close(ncid))

         ! Check
         if (maxval(geom%ic0_to_iproc)>nam%nproc) call msgerror('wrong distribution')
      else
         ! Generate a distribution (use METIS one day?)
         nc0amax = geom%nc0/nam%nproc
         if (nc0amax*nam%nproc<geom%nc0) nc0amax = nc0amax+1
         iproc = 1
         ic0a = 1
         do ic0=1,geom%nc0
            geom%ic0_to_iproc(ic0) = iproc
            geom%ic0_to_ic0a(ic0) = ic0a
            ic0a = ic0a+1
            if (ic0a>nc0amax) then
               ! Change proc
               iproc = iproc+1
               ic0a = 1
            end if
         end do
      end if
   end if
end if

! Allocation
allocate(geom%iproc_to_nc0a(geom%nc0))

! Size of tiles
do iproc=1,nam%nproc
   geom%iproc_to_nc0a(iproc) = count(geom%ic0_to_iproc==iproc)
end do

end subroutine geom_read_local

!----------------------------------------------------------------------
! Subroutine: fld_com_gl
!> Purpose: communicate full field from global to local distribution
!----------------------------------------------------------------------
subroutine fld_com_gl(geom,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom       !< Sampling data
real(kind_real),allocatable,intent(inout) :: fld(:,:)        !< Field

! Local variables
integer :: ic0,ic0a,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      allocate(sbuf(geom%iproc_to_nc0a(iproc)*geom%nl0))

      ! Initialization
      call msr(sbuf)

      ! Prepare buffer
      do ic0=1,geom%nc0
         jproc = geom%ic0_to_iproc(ic0)
         if (jproc==iproc) then
            ic0a = geom%ic0_to_ic0a(ic0)
            sbuf((ic0a-1)*geom%nl0+1:ic0a*geom%nl0) = fld(ic0,:)
         end if
      end do

      if (iproc==mpl%ioproc) then
         ! Allocation
         allocate(rbuf(geom%iproc_to_nc0a(iproc)*geom%nl0))

         ! Copy data
         rbuf = sbuf
      else
         ! Send data to iproc
         call mpl_send(geom%iproc_to_nc0a(iproc)*geom%nl0,sbuf,iproc,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Allocation
   allocate(rbuf(geom%iproc_to_nc0a(mpl%myproc)*geom%nl0))

   ! Receive data from ioproc
   call mpl_recv(geom%iproc_to_nc0a(mpl%myproc)*geom%nl0,rbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Rellocation
if (allocated(fld)) deallocate(fld)
allocate(fld(geom%iproc_to_nc0a(mpl%myproc),geom%nl0))

! Copy from buffer
do ic0a=1,geom%iproc_to_nc0a(mpl%myproc)
   fld(ic0a,:) = rbuf((ic0a-1)*geom%nl0+1:ic0a*geom%nl0)
end do

end subroutine fld_com_gl

!----------------------------------------------------------------------
! Subroutine: fld_com_lg
!> Purpose: communicate full field from local to global distribution
!----------------------------------------------------------------------
subroutine fld_com_lg(geom,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom       !< Sampling data
real(kind_real),allocatable,intent(inout) :: fld(:,:)        !< Field

! Local variables
integer :: ic0,ic0a,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Allocation
allocate(sbuf(geom%iproc_to_nc0a(mpl%myproc)*geom%nl0))

! Prepare buffer
do ic0a=1,geom%iproc_to_nc0a(mpl%myproc)
   sbuf((ic0a-1)*geom%nl0+1:ic0a*geom%nl0) = fld(ic0a,:)
end do

! Release memory
deallocate(fld)

! Communication
if (mpl%main) then
   ! Allocation
   allocate(fld(geom%nc0,geom%nl0))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(rbuf(geom%iproc_to_nc0a(iproc)*geom%nl0))

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Receive data from iproc
         call mpl_recv(geom%iproc_to_nc0a(iproc)*geom%nl0,rbuf,iproc,mpl%tag)
      end if

      ! Copy from buffer
      do ic0=1,geom%nc0
         jproc = geom%ic0_to_iproc(ic0)
         if (jproc==iproc) then
            ic0a = geom%ic0_to_ic0a(ic0)
            fld(ic0,:) = rbuf((ic0a-1)*geom%nl0+1:ic0a*geom%nl0)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Sending data to iproc
   call mpl_send(geom%iproc_to_nc0a(mpl%myproc)*geom%nl0,sbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

end subroutine fld_com_lg

end module type_geom
