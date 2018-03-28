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

use netcdf
use tools_const, only: req,deg2rad,rad2deg
use tools_display, only: msgerror,msgwarning,newunit
use tools_func, only: lonlatmod,sphere_dist,vector_product,vector_triple_product
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsi,msvalr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use tools_stripack, only: areas,trans
use type_ctree, only: ctree_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! Geometry derived type
type geom_type
   ! Vector sizes
   integer :: nlon                            !< Longitude size
   integer :: nlat                            !< Latitude size
   integer :: nlev                            !< Number of levels
   logical,allocatable :: rgmask(:,:)         !< Reduced Gaussian grid mask
   real(kind_real),allocatable :: area(:)     !< Domain area

   ! Number of points and levels
   integer :: nc0                             !< Number of points in subset Sc0
   integer :: nl0                             !< Number of levels in subset Sl0
   integer :: nl0i                            !< Number of independent levels in subset Sl0

   ! Vector coordinates
   real(kind_real),allocatable :: lon(:)      !< Longitudes
   real(kind_real),allocatable :: lat(:)      !< Latitudes
   logical,allocatable :: mask(:,:)           !< Mask
   real(kind_real),allocatable :: vunit(:)    !< Vertical unit
   real(kind_real),allocatable :: disth(:)    !< Horizontal distance
   real(kind_real),allocatable :: distv(:,:)  !< Vertical distance
   logical :: redgrid                         !< Redundant grid

   ! Mesh
   type(mesh_type) :: mesh                    !< Mesh

   ! Cover tree
   type(ctree_type) :: ctree                  !< Cover tree

   ! Boundary nodes
   integer,allocatable :: nbnd(:)             !< Number of boundary nodes
   real(kind_real),allocatable :: xbnd(:,:,:) !< Boundary nodes, x-coordinate
   real(kind_real),allocatable :: ybnd(:,:,:) !< Boundary nodes, y-coordinate
   real(kind_real),allocatable :: zbnd(:,:,:) !< Boundary nodes, z-coordinate
   real(kind_real),allocatable :: vbnd(:,:,:) !< Boundary nodes, orthogonal vector

   ! MPI distribution
   integer :: nc0a                            !< Halo A size
   integer,allocatable :: c0_to_proc(:)       !< Subset Sc0 to local task
   integer,allocatable :: c0_to_c0a(:)        !< Subset Sc0, global to halo A
   integer,allocatable :: proc_to_nc0a(:)     !< Halo A size for each proc
   integer,allocatable :: c0a_to_c0(:)        !< Subset Sc0, halo A to global
contains
   procedure :: alloc => geom_alloc
   procedure :: compute_grid_mesh
   procedure :: define_mask
   procedure :: compute_area
   procedure :: compute_mask_boundaries
   procedure :: compute_metis_graph
   procedure :: read_local
   procedure :: check_arc
   procedure :: fld_com_gl
   procedure :: fld_com_lg
   procedure :: fld_write
end type geom_type

private
public :: geom_type

contains

!----------------------------------------------------------------------
! Subroutine: geom_alloc
!> Purpose: geom object allocation
!----------------------------------------------------------------------
subroutine geom_alloc(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry

! Allocation
allocate(geom%lon(geom%nc0))
allocate(geom%lat(geom%nc0))
allocate(geom%area(geom%nl0))
allocate(geom%mask(geom%nc0,geom%nl0))
allocate(geom%vunit(geom%nl0))
allocate(geom%distv(geom%nl0,geom%nl0))

! Initialization
call msr(geom%lon)
call msr(geom%lat)
geom%mask = .false.
call msr(geom%vunit)
call msr(geom%distv)
geom%redgrid = .true.

end subroutine geom_alloc

!----------------------------------------------------------------------
! Subroutine: compute_grid_mesh
!> Purpose: compute grid mesh
!----------------------------------------------------------------------
subroutine compute_grid_mesh(geom,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(nam_type),intent(in) :: nam       !< Namelist

! Local variables
integer :: ic0,il0,jl0,jc3
logical :: same_mask,ctree_mask(geom%nc0)

! Set longitude and latitude bounds
do ic0=1,geom%nc0
   call lonlatmod(geom%lon(ic0),geom%lat(ic0))
end do

! Define mask
call geom%define_mask(nam)

! Create mesh
if ((.not.all(geom%area>0.0)).or.(nam%new_hdiag.and.nam%displ_diag).or.((nam%new_param.or.nam%new_lct) &
 & .and.nam%mask_check).or.(nam%new_param.and.nam%network).or.nam%new_obsop) &
 & call geom%mesh%create(geom%nc0,geom%lon,geom%lat,geom%redgrid)

! Compute area
if ((.not.all(geom%area>0.0))) call geom%compute_area

! Compute mask boundaries
if ((nam%new_param.or.nam%new_lct).and.nam%mask_check) call geom%compute_mask_boundaries

! Check whether the mask is the same for all levels
same_mask = .true.
do il0=2,geom%nl0
   same_mask = same_mask.and.(all((geom%mask(:,il0).and.geom%mask(:,1)) &
             & .or.(.not.geom%mask(:,il0).and..not.geom%mask(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   geom%nl0i = 1
else
   geom%nl0i = geom%nl0
end if
write(mpl%unit,'(a7,a,i3)') '','Number of independent levels: ',geom%nl0i

! Create cover tree
ctree_mask = .true.
if (nam%new_hdiag.or.nam%check_dirac.or.nam%new_lct.or.nam%new_obsop) &
 & call geom%ctree%create(geom%nc0,geom%lon,geom%lat,ctree_mask)

! Vertical distance
do jl0=1,geom%nl0
   do il0=1,geom%nl0
      geom%distv(il0,jl0) = abs(geom%vunit(il0)-geom%vunit(jl0))
   end do
end do

! Horizontal distance
allocate(geom%disth(nam%nc3))
do jc3=1,nam%nc3
   geom%disth(jc3) = float(jc3-1)*nam%dc
end do

! Read local distribution
call geom%read_local(nam)

end subroutine compute_grid_mesh

!----------------------------------------------------------------------
! Subroutine: define_mask
!> Purpose: define mask
!----------------------------------------------------------------------
subroutine define_mask(geom,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(nam_type),intent(in) :: nam       !< Namelist

! Local variables
integer :: latmin,latmax,il0,ic0,ildw
integer :: ncid,nlon_id,nlon_test,nlat_id,nlat_test,mask_id
real(kind_real) :: dist
real(kind_real),allocatable :: hydmask(:,:)
logical :: mask_test
character(len=3) :: il0char
character(len=1024) :: subr = 'define_mask'

! Mask restriction
if (nam%mask_type(1:3)=='lat') then
   ! Latitude mask
   read(nam%mask_type(4:6),'(i3)') latmin
   read(nam%mask_type(7:9),'(i3)') latmax
   if (latmin>=latmax) call msgerror('latmin should be lower than latmax')
   do il0=1,geom%nl0
      geom%mask(:,il0) = geom%mask(:,il0).and.(geom%lat>=float(latmin)*deg2rad).and.(geom%lat<=float(latmax)*deg2rad)
   end do
elseif (trim(nam%mask_type)=='hyd') then
   ! Read from hydrometeors mask file
   call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_hyd.nc',nf90_nowrite,ncid))
   if (trim(nam%model)=='aro') then
      call ncerr(subr,nf90_inq_dimid(ncid,'X',nlon_id))
      call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=nlon_test))
      call ncerr(subr,nf90_inq_dimid(ncid,'Y',nlat_id))
      call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=nlat_test))
      if ((nlon_test/=geom%nlon).or.(nlat_test/=geom%nlat)) call msgerror('wrong dimensions in the mask')
      allocate(hydmask(geom%nlon,geom%nlat))
      do il0=1,geom%nl0
         write(il0char,'(i3.3)') nam%levs(il0)
         call ncerr(subr,nf90_inq_varid(ncid,'S'//il0char//'MASK',mask_id))
         call ncerr(subr,nf90_get_var(ncid,mask_id,hydmask,(/1,1/),(/geom%nlon,geom%nlat/)))
         geom%mask(:,il0) = geom%mask(:,il0).and.pack(real(hydmask,kind(1.0))>nam%mask_th,mask=.true.)
      end do
      deallocate(hydmask)
      call ncerr(subr,nf90_close(ncid))
   end if
elseif (trim(nam%mask_type)=='ldwv') then
   ! Compute distance to the vertical diagnostic points
   do ic0=1,geom%nc0
      if (any(geom%mask(ic0,:))) then
         mask_test = .false.
         do ildw=1,nam%nldwv
            call sphere_dist(nam%lon_ldwv(ildw)*deg2rad,nam%lat_ldwv(ildw)*deg2rad,geom%lon(ic0),geom%lat(ic0),dist)
            mask_test = mask_test.or.(dist<1.1*nam%local_rad)
         end do
         do il0=1,geom%nl0
            if (geom%mask(ic0,il0)) geom%mask(ic0,:) = mask_test
         end do
      end if
   end do
end if

end subroutine define_mask

!----------------------------------------------------------------------
! Subroutine: compute_area
!> Purpose: compute domain area
!----------------------------------------------------------------------
subroutine compute_area(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: il0,it
real(kind_real) :: area,frac

! Create triangles list
if (.not.allocated(geom%mesh%ltri)) call geom%mesh%trlist

! Compute area
geom%area = 0.0
do it=1,geom%mesh%nt
   area = areas((/geom%mesh%x(geom%mesh%ltri(1,it)),geom%mesh%y(geom%mesh%ltri(1,it)),geom%mesh%z(geom%mesh%ltri(1,it))/), &
              & (/geom%mesh%x(geom%mesh%ltri(2,it)),geom%mesh%y(geom%mesh%ltri(2,it)),geom%mesh%z(geom%mesh%ltri(2,it))/), &
              & (/geom%mesh%x(geom%mesh%ltri(3,it)),geom%mesh%y(geom%mesh%ltri(3,it)),geom%mesh%z(geom%mesh%ltri(3,it))/))
   do il0=1,geom%nl0
      frac = float(count(geom%mask(geom%mesh%order(geom%mesh%ltri(1:3,it)),il0)))/3.0
      geom%area(il0) = geom%area(il0)+frac*area
   end do
end do

end subroutine compute_area

!----------------------------------------------------------------------
! Subroutine: compute_mask_boundaries
!> Purpose: compute domain area
!----------------------------------------------------------------------
subroutine compute_mask_boundaries(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry

! Local variables
integer :: inr,jnr,knr,ic0,jc0,kc0,i,ibnd,il0
integer,allocatable :: ic0_bnd(:,:,:)
real(kind_real) :: latbnd(2),lonbnd(2),v1(3),v2(3)
logical :: init

! Allocation
allocate(geom%nbnd(geom%nl0))
allocate(ic0_bnd(2,geom%mesh%nnr,geom%nl0))

! Find border points
do il0=1,geom%nl0
   geom%nbnd(il0) = 0
   do inr=1,geom%mesh%nnr
      ! Check mask points only
      ic0 = geom%mesh%order(inr)
      if (.not.geom%mask(ic0,il0)) then
         i = geom%mesh%lend(inr)
         init = .true.
         do while ((i/=geom%mesh%lend(inr)).or.init)
            jnr = abs(geom%mesh%list(i))
            knr = abs(geom%mesh%list(geom%mesh%lptr(i)))
            jc0 = geom%mesh%order(jnr)
            kc0 = geom%mesh%order(knr)
            if (.not.geom%mask(jc0,il0).and.geom%mask(kc0,il0)) then
               ! Create a new boundary arc
               geom%nbnd(il0) = geom%nbnd(il0)+1
               if (geom%nbnd(il0)>geom%mesh%nnr) call msgerror('too many boundary arcs')
               ic0_bnd(1,geom%nbnd(il0),il0) = ic0
               ic0_bnd(2,geom%nbnd(il0),il0) = jc0
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
! Subroutine: compute_metis_graph
!> Purpose: compute METIS graph
!----------------------------------------------------------------------
subroutine compute_metis_graph(geom,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(nam_type),intent(in) :: nam       !< Namelist


! Local variables
integer :: inr,jnr,i,lunit
logical :: init
character(len=1024) :: filename

if (mpl%main) then
   ! Compute triangles list
   write(mpl%unit,'(a7,a)') '','Compute METIS graph'
   call geom%mesh%bnodes

   ! Open file
   filename = trim(nam%prefix)//'_metis'
   lunit = newunit()
   open(unit=lunit,file=trim(nam%datadir)//'/'//trim(filename),status='replace')

   ! Write header
   write(lunit,*) geom%mesh%nnr,geom%mesh%na-geom%mesh%nb/2

   ! Write connectivity
   do inr=1,geom%mesh%nnr
      i = geom%mesh%lend(inr)
      init = .true.
      do while ((i/=geom%mesh%lend(inr)).or.init)
         jnr = geom%mesh%list(i)
         if (jnr>0) write(lunit,'(i7)',advance='no') jnr
         i = geom%mesh%lptr(i)
         init = .false.
      end do
      write(lunit,*) ''
   end do

   ! Close file
   close(unit=lunit)
end if

end subroutine compute_metis_graph

!----------------------------------------------------------------------
! Subroutine: read_local
!> Purpose: read local distribution
!----------------------------------------------------------------------
subroutine read_local(geom,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom !< Geometry
type(nam_type),intent(in) :: nam       !< Namelist

! Local variables
integer :: ic0,inr,info,iproc,ic0a,nc0amax,lunit
integer :: ncid,nc0_id,c0_to_proc_id,c0_to_c0a_id,lon_id,lat_id
integer,allocatable :: nr_to_proc(:),ic0a_arr(:)
logical :: ismetis
character(len=4) :: nprocchar
character(len=1024) :: filename_nc,filename_metis
character(len=1024) :: subr = 'read_local'

if (.not.allocated(geom%c0_to_proc)) then
   ! Allocation
   allocate(geom%c0_to_proc(geom%nc0))
   allocate(geom%c0_to_c0a(geom%nc0))

   if (mpl%nproc==1) then
      ! All points on a single processor
      geom%c0_to_proc = 1
      do ic0=1,geom%nc0
         geom%c0_to_c0a(ic0) = ic0
      end do
   elseif (mpl%nproc>1) then
      ! Open file
      write(nprocchar,'(i4.4)') mpl%nproc
      filename_nc = trim(nam%prefix)//'_distribution_'//nprocchar//'.nc'
      info = nf90_open(trim(nam%datadir)//'/'//trim(filename_nc),nf90_nowrite,ncid)

      if (info==nf90_noerr) then
         ! Read local distribution
         write(mpl%unit,'(a7,a,i4,a)') '','Read local distribution for ',mpl%nproc,' MPI tasks'

         ! Get variables ID
         call ncerr(subr,nf90_inq_varid(ncid,'c0_to_proc',c0_to_proc_id))
         call ncerr(subr,nf90_inq_varid(ncid,'c0_to_c0a',c0_to_c0a_id))

         ! Read varaibles
         call ncerr(subr,nf90_get_var(ncid,c0_to_proc_id,geom%c0_to_proc))
         call ncerr(subr,nf90_get_var(ncid,c0_to_c0a_id,geom%c0_to_c0a))

         ! Close file
         call ncerr(subr,nf90_close(ncid))

         ! Check
         if (maxval(geom%c0_to_proc)>mpl%nproc) call msgerror('wrong distribution')
      else
         ! Create a mesh
         if (.not.allocated(geom%mesh%redundant)) call geom%mesh%create(geom%nc0,geom%lon,geom%lat,.true.)

         ! Generate a distribution
         if (mpl%main) then
            if (nam%use_metis) then
               ! Try to use METIS
               write(mpl%unit,'(a7,a,i4,a)') '','Try to use METIS for ',mpl%nproc,' MPI tasks'

               ! Compute graph
               call geom%compute_metis_graph(nam)

               ! Write graph
               filename_metis = trim(nam%prefix)//'_metis'
               write(nprocchar,'(i4)') mpl%nproc

               ! Call METIS
               call system('gpmetis '//trim(nam%datadir)//'/'//trim(filename_metis)//' '//adjustl(nprocchar)//' > '// &
             & trim(nam%datadir)//'/'//trim(filename_metis)//'.out')

               ! Check for METIS output
               inquire(file=trim(nam%datadir)//'/'//trim(filename_metis)//'.part.'//adjustl(nprocchar),exist=ismetis)
               if (.not.ismetis) call msgwarning('METIS not available to generate the local distribution')
            else
               ! No METIS
               ismetis = .false.
            end if

            if (ismetis) then
               write(mpl%unit,'(a7,a)') '','Use METIS to generate the local distribution'

               ! Allocation
               allocate(nr_to_proc(geom%mesh%nnr))
               allocate(ic0a_arr(mpl%nproc))

               ! Read METIS file
               lunit = newunit()
               open(unit=lunit,file=trim(nam%datadir)//'/'//trim(filename_metis)//'.part.'//adjustl(nprocchar),status='old')
               do inr=1,geom%mesh%nnr
                  read(lunit,*) nr_to_proc(inr)
               end do
               close(unit=lunit)

               ! Reorder and offset
               do ic0=1,geom%nc0
                  inr = geom%mesh%order_inv(ic0)
                  geom%c0_to_proc(ic0) = nr_to_proc(inr)+1
               end do

               ! Local index
               ic0a_arr = 0
               do ic0=1,geom%nc0
                  iproc = geom%c0_to_proc(ic0)
                  ic0a_arr(iproc) = ic0a_arr(iproc)+1
                  geom%c0_to_c0a(ic0) = ic0a_arr(iproc)
               end do
            else
               write(mpl%unit,'(a7,a)') '','Define a basic local distribution'

               ! Basic distribution
               nc0amax = geom%nc0/mpl%nproc
               if (nc0amax*mpl%nproc<geom%nc0) nc0amax = nc0amax+1
               iproc = 1
               ic0a = 1
               do ic0=1,geom%nc0
                  geom%c0_to_proc(ic0) = iproc
                  geom%c0_to_c0a(ic0) = ic0a
                  ic0a = ic0a+1
                  if (ic0a>nc0amax) then
                     ! Change proc
                     iproc = iproc+1
                     ic0a = 1
                  end if
               end do
            end if
         end if

         ! Broadcast distribution
         call mpl%bcast(geom%c0_to_proc,mpl%ioproc)
         call mpl%bcast(geom%c0_to_c0a,mpl%ioproc)

         ! Write distribution
         if (mpl%main) then
            ! Create file
            call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename_nc),or(nf90_clobber,nf90_64bit_offset),ncid))

            ! Write namelist parameters
            call nam%ncwrite(ncid)

            ! Define dimension
            call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))

            ! Define variables
            call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
            call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
            call ncerr(subr,nf90_def_var(ncid,'c0_to_proc',nf90_int,(/nc0_id/),c0_to_proc_id))
            call ncerr(subr,nf90_def_var(ncid,'c0_to_c0a',nf90_int,(/nc0_id/),c0_to_c0a_id))

            ! End definition mode
            call ncerr(subr,nf90_enddef(ncid))

            ! Write variables
            call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
            call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
            call ncerr(subr,nf90_put_var(ncid,c0_to_proc_id,geom%c0_to_proc))
            call ncerr(subr,nf90_put_var(ncid,c0_to_c0a_id,geom%c0_to_c0a))

            ! Close file
            call ncerr(subr,nf90_close(ncid))
         end if
      end if
   end if
end if

! Size of tiles
allocate(geom%proc_to_nc0a(mpl%nproc))
do iproc=1,mpl%nproc
   geom%proc_to_nc0a(iproc) = count(geom%c0_to_proc==iproc)
end do
geom%nc0a = geom%proc_to_nc0a(mpl%myproc)

! Conversion
allocate(geom%c0a_to_c0(geom%nc0a))
ic0a = 0
do ic0=1,geom%nc0
   if (geom%c0_to_proc(ic0)==mpl%myproc) then
      ic0a = ic0a+1
      geom%c0a_to_c0(ic0a) = ic0
   end if
end do

end subroutine read_local

!----------------------------------------------------------------------
! Subroutine: check_arc
!> Purpose: check if an arc is crossing boundaries
!----------------------------------------------------------------------
subroutine check_arc(geom,il0,lon_s,lat_s,lon_e,lat_e,valid)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom !< Geometry
integer,intent(in) :: il0           !< Level
real(kind_real),intent(in) :: lon_s !< First point longitude
real(kind_real),intent(in) :: lat_s !< First point latitude
real(kind_real),intent(in) :: lon_e !< Second point longitude
real(kind_real),intent(in) :: lat_e !< Second point latitude
logical,intent(out) :: valid        !< True for valid arcs

! Local variables
integer :: ibnd
real(kind_real) :: x(2),y(2),z(2),v1(3),v2(3),va(3),vp(3),t(4)

! Transform to cartesian coordinates
call trans(2,(/lat_s,lat_e/),(/lon_s,lon_e/),x,y,z)

! Compute arc orthogonal vector
v1 = (/x(1),y(1),z(1)/)
v2 = (/x(2),y(2),z(2)/)
call vector_product(v1,v2,va)

! Check if arc is crossing boundary arcs
valid = .true.
do ibnd=1,geom%nbnd(il0)
   call vector_product(va,geom%vbnd(:,ibnd,il0),vp)
   v1 = (/x(1),y(1),z(1)/)
   call vector_triple_product(v1,va,vp,t(1))
   v1 = (/x(2),y(2),z(2)/)
   call vector_triple_product(v1,va,vp,t(2))
   v1 = (/geom%xbnd(1,ibnd,il0),geom%ybnd(1,ibnd,il0),geom%zbnd(1,ibnd,il0)/)
   call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(3))
   v1 = (/geom%xbnd(2,ibnd,il0),geom%ybnd(2,ibnd,il0),geom%zbnd(2,ibnd,il0)/)
   call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(4))
   t(1) = -t(1)
   t(3) = -t(3)
   if (all(t>0).or.(all(t<0))) then
      valid = .false.
      exit
   end if
end do

end subroutine check_arc

!----------------------------------------------------------------------
! Subroutine: fld_com_gl
!> Purpose: communicate full field from global to local distribution
!----------------------------------------------------------------------
subroutine fld_com_gl(geom,fld_glb,fld_loc)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                        !< Geometry
real(kind_real),intent(in) :: fld_glb(geom%nc0,geom%nl0)   !< Field (global)
real(kind_real),intent(out) :: fld_loc(geom%nc0a,geom%nl0) !< Field (local)

! Local variables
integer :: il0,ic0,ic0a,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:),fld_tmp(:,:)
logical :: mask_unpack(geom%nc0a,geom%nl0)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      allocate(fld_tmp(geom%proc_to_nc0a(iproc),geom%nl0))
      allocate(sbuf(geom%proc_to_nc0a(iproc)*geom%nl0))

      ! Initialization
      call msr(sbuf)

      ! Prepare buffer
      do il0=1,geom%nl0
         do ic0=1,geom%nc0
            jproc = geom%c0_to_proc(ic0)
            if (jproc==iproc) then
               ic0a = geom%c0_to_c0a(ic0)
               fld_tmp(ic0a,il0) = fld_glb(ic0,il0)
            end if
         end do
      end do
      sbuf = pack(fld_tmp,mask=.true.)

      if (iproc==mpl%ioproc) then
         ! Allocation
         allocate(rbuf(geom%proc_to_nc0a(iproc)*geom%nl0))

         ! Copy data
         rbuf = sbuf
      else
         ! Send data to iproc
         call mpl%send(geom%proc_to_nc0a(iproc)*geom%nl0,sbuf,iproc,mpl%tag)
      end if

      ! Release memory
      deallocate(fld_tmp)
      deallocate(sbuf)
   end do
else
   ! Allocation
   allocate(rbuf(geom%nc0a*geom%nl0))

   ! Receive data from ioproc
   call mpl%recv(geom%nc0a*geom%nl0,rbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Copy from buffer
mask_unpack = .true.
call msr(fld_loc)
fld_loc = unpack(rbuf,mask_unpack,fld_loc)

end subroutine fld_com_gl

!----------------------------------------------------------------------
! Subroutine: fld_com_lg
!> Purpose: communicate full field from local to global distribution
!----------------------------------------------------------------------
subroutine fld_com_lg(geom,fld_loc,fld_glb)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                       !< Geometry
real(kind_real),intent(in) :: fld_loc(geom%nc0a,geom%nl0) !< Field (local)
real(kind_real),intent(out) :: fld_glb(geom%nc0,geom%nl0) !< Field (global)

! Local variables
integer :: il0,ic0,ic0a,iproc,jproc
real(kind_real),allocatable :: sbuf(:),rbuf(:),fld_tmp(:,:)
logical,allocatable :: mask_unpack(:,:)

! Allocation
allocate(sbuf(geom%nc0a*geom%nl0))

! Prepare buffer
sbuf = pack(fld_loc,mask=.true.)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      allocate(fld_tmp(geom%proc_to_nc0a(iproc),geom%nl0))
      allocate(mask_unpack(geom%proc_to_nc0a(iproc),geom%nl0))
      allocate(rbuf(geom%proc_to_nc0a(iproc)*geom%nl0))

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Receive data from iproc
         call mpl%recv(geom%proc_to_nc0a(iproc)*geom%nl0,rbuf,iproc,mpl%tag)
      end if

      ! Copy from buffer
      mask_unpack = .true.
      fld_tmp = unpack(rbuf,mask_unpack,fld_tmp)
      do il0=1,geom%nl0
         do ic0=1,geom%nc0
            jproc = geom%c0_to_proc(ic0)
            if (jproc==iproc) then
               ic0a = geom%c0_to_c0a(ic0)
               fld_glb(ic0,il0) = fld_tmp(ic0a,il0)
            end if
         end do
      end do

      ! Release memory
      deallocate(fld_tmp)
      deallocate(mask_unpack)
      deallocate(rbuf)
   end do
else
   ! Sending data to iproc
   call mpl%send(geom%nc0a*geom%nl0,sbuf,mpl%ioproc,mpl%tag)

   ! Set global field to missing value
   call msr(fld_glb)
end if
mpl%tag = mpl%tag+1

end subroutine fld_com_lg

!----------------------------------------------------------------------
! Subroutine: fld_write
!> Purpose: write field
!----------------------------------------------------------------------
subroutine fld_write(geom,nam,filename,varname,fld)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                   !< Geometry
type(nam_type),intent(in) :: nam                      !< Namelist
character(len=*),intent(in) :: filename               !< File name
character(len=*),intent(in) :: varname                !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
integer :: ic0a,ic0,jc0a,jc0,il0,info
integer :: ncid,nc0_id,nlev_id,fld_id,lon_id,lat_id
real(kind_real) :: fld_loc(geom%nc0a,geom%nl0),fld_glb(geom%nc0,geom%nl0)
character(len=1024) :: subr = 'fld_write'

! Apply mask
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      if (geom%mask(ic0,il0)) then
         fld_loc(ic0a,il0) = fld(ic0a,il0)
      else
         call msr(fld_loc(ic0a,il0))
      end if
   end do
end do

if (allocated(geom%mesh%redundant)) then
   ! Copy redundant points
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      jc0 = geom%mesh%redundant(ic0)
      if (isnotmsi(jc0)) then
         jc0a = geom%c0_to_c0a(jc0)
         fld_loc(ic0a,:) = fld_loc(jc0a,:)
      end if
   end do
end if

! Local to global
call geom%fld_com_lg(fld_loc,fld_glb)

if (mpl%main) then
   ! Check if the file exists
   info = nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_noclobber,nf90_64bit_offset),ncid)
   if (info==nf90_noerr) then
      ! Write namelist parameters
      call nam%ncwrite(ncid)

      ! Define attribute
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))

      ! End definition mode
      call ncerr(subr,nf90_enddef(ncid))
   else
      ! Open file
      call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_write,ncid))
   end if

   ! Get variable id
   info = nf90_inq_varid(ncid,trim(varname),fld_id)

   ! Define dimensions and variable if necessary
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      info = nf90_inq_dimid(ncid,'nc0',nc0_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
      info = nf90_inq_dimid(ncid,'nlev',nlev_id)
      if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nlev',geom%nl0,nlev_id))
      call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlev_id,nc0_id/),fld_id))
      call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
   end if

   ! Write data
   do il0=1,geom%nl0
      if (isanynotmsr(fld_glb(:,il0))) then
         call ncerr(subr,nf90_put_var(ncid,fld_id,fld_glb(:,il0),(/il0,1/),(/1,geom%nc0/)))
      end if
   end do

   ! Write coordinates
   info = nf90_inq_varid(ncid,'lon',lon_id)
   if (info/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
      call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
      call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
      call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
      call ncerr(subr,nf90_enddef(ncid))
      call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
      call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
   end if

   ! Close file
   call ncerr(subr,nf90_close(ncid))
end if

end subroutine fld_write

end module type_geom
