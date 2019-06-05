!----------------------------------------------------------------------
! Module: type_mesh
! Purpose: mesh derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mesh

!$ use omp_lib
use tools_const, only: pi,req
use tools_func, only: sphere_dist,lonlat2xyz,xyz2lonlat,vector_product
use tools_kinds, only: kind_real
use tools_stripack, only: addnod,areas,bnodes,crlist,inside,trfind,trlist,trmesh
use type_mpl, only: mpl_type
use type_rng, only: rng_type

implicit none

logical,parameter :: shuffle = .true. ! Shuffle mesh order (more efficient to compute the Delaunay triangulation)

! Mesh derived type
type mesh_type
   ! Mesh structure
   integer :: n                                 ! Number of points
   integer,allocatable :: order(:)              ! Order of shuffled points
   integer,allocatable :: order_inv(:)          ! Inverse order of shuffled points
   real(kind_real),allocatable :: lon(:)        ! Points longitudes
   real(kind_real),allocatable :: lat(:)        ! Points latitudes
   real(kind_real),allocatable :: x(:)          ! x-coordinate
   real(kind_real),allocatable :: y(:)          ! y-coordinate
   real(kind_real),allocatable :: z(:)          ! z-coordinate
   integer,allocatable :: list(:)               ! Stripack list
   integer,allocatable :: lptr(:)               ! Stripack list pointer
   integer,allocatable :: lend(:)               ! Stripack list end
   integer :: lnew                              ! Stripack pointer to the first empty location in list
   integer :: nb                                ! Number of boundary nodes
   integer,allocatable :: bnd(:)                ! Boundary nodes
   integer,allocatable :: barc(:,:)             ! Boundary arcs
   real(kind_real),allocatable :: barc_lon(:,:) ! Boundary arcs longitudes
   real(kind_real),allocatable :: barc_lat(:,:) ! Boundary arcs latitudes
   real(kind_real),allocatable :: barc_dist(:)  ! Boundary arcs distance
   real(kind_real),allocatable :: barc_vp(:,:)  ! Boundary arcs normal vector
   real(kind_real),allocatable :: bdist(:)      ! Distance to the closest boundary arc

   ! Triangles data
   integer :: nt                                ! Number of triangles
   integer :: na                                ! Number of arcs
   integer,allocatable :: ltri(:,:)             ! Triangles indices
   integer,allocatable :: larc(:,:)             ! Arcs indices
   logical,allocatable :: valid(:)              ! Valid mesh nodes
contains
   procedure :: alloc => mesh_alloc
   procedure :: init => mesh_init
   procedure :: dealloc => mesh_dealloc
   procedure :: copy => mesh_copy
   procedure :: store => mesh_store
   procedure :: trlist => mesh_trlist
   procedure :: bnodes => mesh_bnodes
   procedure :: find_bdist => mesh_find_bdist
   procedure :: check => mesh_check
   procedure :: inside => mesh_inside
   procedure :: barycentric => mesh_barycentric
   procedure :: addnode => mesh_addnode
   procedure :: polygon => mesh_polygon
end type mesh_type

private
public :: mesh_type

contains

!----------------------------------------------------------------------
! Subroutine: mesh_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine mesh_alloc(mesh,n)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh
integer,intent(in) :: n                ! Mesh size

! Allocation
mesh%n = n

! Allocation
allocate(mesh%order(mesh%n))
allocate(mesh%order_inv(mesh%n))
allocate(mesh%lon(mesh%n))
allocate(mesh%lat(mesh%n))
allocate(mesh%x(mesh%n))
allocate(mesh%y(mesh%n))
allocate(mesh%z(mesh%n))
allocate(mesh%list(6*(mesh%n-2)))
allocate(mesh%lptr(6*(mesh%n-2)))
allocate(mesh%lend(mesh%n))

end subroutine mesh_alloc

!----------------------------------------------------------------------
! Subroutine: mesh_init
! Purpose: intialization
!----------------------------------------------------------------------
subroutine mesh_init(mesh,mpl,rng,lon,lat)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh    ! Mesh
type(mpl_type),intent(inout) :: mpl       ! MPI data
type(rng_type),intent(inout) :: rng       ! Random number generator
real(kind_real),intent(in) :: lon(mesh%n) ! Longitudes
real(kind_real),intent(in) :: lat(mesh%n) ! Latitudes

! Local variables
integer :: i,k,info
integer :: jtab(mesh%n),near(mesh%n),next(mesh%n)
real(kind_real) :: dist(mesh%n)

! Points order
do i=1,mesh%n
   mesh%order(i) = i
end do

if (shuffle) then
   ! Shuffle order (more efficient to compute the Delaunay triangulation)
   if (mpl%main) call rng%rand_integer(1,mesh%n,jtab)
   call mpl%f_comm%broadcast(jtab,mpl%ioproc-1)
   do i=mesh%n,2,-1
      k = mesh%order(jtab(i))
      mesh%order(jtab(i)) = mesh%order(i)
      mesh%order(i) = k
   end do
end if

! Restrictive inverse order
mesh%order_inv = mpl%msv%vali
do i=1,mesh%n
   mesh%order_inv(mesh%order(i)) = i
end do

! Store coordinates
call mesh%store(mpl,lon,lat)

! Create mesh
mesh%list = 0
call trmesh(mpl,mesh%n,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,mesh%lnew,near,next,dist,info)

! Boundaries not computed yet
mesh%nb = mpl%msv%vali

end subroutine mesh_init

!----------------------------------------------------------------------
! Subroutine: mesh_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine mesh_dealloc(mesh)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh

! Release memory
if (allocated(mesh%order)) deallocate(mesh%order)
if (allocated(mesh%order_inv)) deallocate(mesh%order_inv)
if (allocated(mesh%lon)) deallocate(mesh%lon)
if (allocated(mesh%lat)) deallocate(mesh%lat)
if (allocated(mesh%x)) deallocate(mesh%x)
if (allocated(mesh%y)) deallocate(mesh%y)
if (allocated(mesh%z)) deallocate(mesh%z)
if (allocated(mesh%list)) deallocate(mesh%list)
if (allocated(mesh%lptr)) deallocate(mesh%lptr)
if (allocated(mesh%lend)) deallocate(mesh%lend)
if (allocated(mesh%bnd)) deallocate(mesh%bnd)
if (allocated(mesh%barc)) deallocate(mesh%barc)
if (allocated(mesh%barc_lon)) deallocate(mesh%barc_lon)
if (allocated(mesh%barc_lat)) deallocate(mesh%barc_lat)
if (allocated(mesh%barc_dist)) deallocate(mesh%barc_dist)
if (allocated(mesh%barc_vp)) deallocate(mesh%barc_vp)
if (allocated(mesh%bdist)) deallocate(mesh%bdist)
if (allocated(mesh%ltri)) deallocate(mesh%ltri)
if (allocated(mesh%larc)) deallocate(mesh%larc)
if (allocated(mesh%valid)) deallocate(mesh%valid)

end subroutine mesh_dealloc

!----------------------------------------------------------------------
! Subroutine: mesh_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine mesh_copy(mesh_out,mesh_in)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh_out ! Output mesh
type(mesh_type),intent(in) :: mesh_in      ! Input mesh

! Release memory
call mesh_out%dealloc

! Allocation
call mesh_out%alloc(mesh_in%n)
if (allocated(mesh_in%bnd)) allocate(mesh_out%bnd(mesh_in%nb))
if (allocated(mesh_in%ltri)) allocate(mesh_out%ltri(3,mesh_in%nt))
if (allocated(mesh_in%larc)) allocate(mesh_out%larc(2,mesh_in%na))
if (allocated(mesh_in%bdist)) allocate(mesh_out%bdist(mesh_in%n))

! Copy data
mesh_out%order = mesh_in%order
mesh_out%order_inv = mesh_in%order_inv
mesh_out%lon = mesh_in%lon
mesh_out%lat = mesh_in%lat
mesh_out%x = mesh_in%x
mesh_out%y = mesh_in%y
mesh_out%z = mesh_in%z
mesh_out%list = mesh_in%list
mesh_out%lptr = mesh_in%lptr
mesh_out%lend = mesh_in%lend
mesh_out%lnew = mesh_in%lnew
mesh_out%nb = mesh_in%nb
if (allocated(mesh_in%barc)) mesh_out%barc = mesh_in%barc
if (allocated(mesh_in%barc_lon)) mesh_out%barc_lon = mesh_in%barc_lon
if (allocated(mesh_in%barc_lat)) mesh_out%barc_lat = mesh_in%barc_lat
if (allocated(mesh_in%barc_dist)) mesh_out%barc_dist = mesh_in%barc_dist
if (allocated(mesh_in%barc_vp)) mesh_out%barc_vp = mesh_in%barc_vp
mesh_out%nt = mesh_in%nt
mesh_out%na = mesh_in%na
if (allocated(mesh_in%bnd)) mesh_out%bnd = mesh_in%bnd
if (allocated(mesh_in%ltri)) mesh_out%ltri = mesh_in%ltri
if (allocated(mesh_in%larc)) mesh_out%larc = mesh_in%larc
if (allocated(mesh_in%valid)) mesh_out%valid = mesh_in%valid

end subroutine mesh_copy

!----------------------------------------------------------------------
! Subroutine: mesh_store
! Purpose: store mesh cartesian coordinates
!----------------------------------------------------------------------
subroutine mesh_store(mesh,mpl,lon,lat)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh    ! Mesh
type(mpl_type),intent(inout) :: mpl       ! MPI data
real(kind_real),intent(in) :: lon(mesh%n) ! Longitude
real(kind_real),intent(in) :: lat(mesh%n) ! Latitude

! Local variables
integer :: i

! Initialize
mesh%lon = mpl%msv%valr
mesh%lat = mpl%msv%valr

! Copy lon/lat
mesh%lon = lon(mesh%order)
mesh%lat = lat(mesh%order)

! Transform to cartesian coordinates
do i=1,mesh%n
   call lonlat2xyz(mpl,mesh%lon(i),mesh%lat(i),mesh%x(i),mesh%y(i),mesh%z(i))
end do

end subroutine mesh_store

!----------------------------------------------------------------------
! Subroutine: mesh_trlist
! Purpose: compute triangle list, arc list
!----------------------------------------------------------------------
subroutine mesh_trlist(mesh,mpl)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh
type(mpl_type),intent(inout) :: mpl    ! MPI data

! Local variables
integer :: info,ia,it,i,i1,i2
integer :: ltri(9,2*(mesh%n-2))
character(len=6) :: notvalidchar
character(len=1024),parameter :: subr = 'mesh_trlist'

! Create triangles list
call trlist(mesh%n,mesh%list,mesh%lptr,mesh%lend,9,mesh%nt,ltri,info)

! Allocation
mesh%na = maxval(ltri(7:9,1:mesh%nt))
allocate(mesh%ltri(3,mesh%nt))
allocate(mesh%larc(2,mesh%na))
allocate(mesh%valid(mesh%n))

! Copy triangle list
mesh%ltri = ltri(1:3,1:mesh%nt)

! Copy arcs list
do ia=1,mesh%na
   it = 1
   do while (it<=mesh%nt)
      if (any(ltri(7:9,it)==ia)) exit
      it = it+1
   end do
   i = 1
   do while (i<=3)
      if (ltri(6+i,it)==ia) exit
      i = i+1
   end do
   i1 = mod(i+1,3)
   if (i1==0) i1 = 3
   i2 = mod(i+2,3)
   if (i2==0) i2 = 3
   mesh%larc(1,ia) = ltri(i1,it)
   mesh%larc(2,ia) = ltri(i2,it)
end do

! Check mesh
call mesh%check(mpl,mesh%valid)
if (.not.all(mesh%valid)) then
   write(notvalidchar,'(i6)') count(.not.mesh%valid)
   call mpl%warning(subr,'unvalid mesh at creation ('//notvalidchar//' points)')
end if

end subroutine mesh_trlist

!----------------------------------------------------------------------
! Subroutine: mesh_bnodes
! Purpose: find boundary nodes
!----------------------------------------------------------------------
subroutine mesh_bnodes(mesh,mpl,bdist)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh
type(mpl_type),intent(inout) :: mpl    ! MPI data
logical,intent(in),optional :: bdist   ! Find minimum distance a boundary arc

! Local variables
integer :: i,bnd(mesh%n)
real(kind_real) :: v1(3),v2(3)
logical :: lbdist

! Initialization
lbdist = .false.
if (present(bdist)) lbdist = bdist

! Find boundary nodes
bnd = mpl%msv%vali
call bnodes(mesh%n,mesh%list,mesh%lptr,mesh%lend,bnd,mesh%nb,mesh%na,mesh%nt)

! Allocation
allocate(mesh%bnd(mesh%nb))

! Copy
mesh%bnd = bnd(1:mesh%nb)

! Allocation
if (mesh%nb>0) then
   allocate(mesh%barc(2,mesh%nb))
   allocate(mesh%barc_lon(2,mesh%nb))
   allocate(mesh%barc_lat(2,mesh%nb))
   allocate(mesh%barc_dist(mesh%nb))
   allocate(mesh%barc_vp(3,mesh%nb))
end if
allocate(mesh%bdist(mesh%n))

! Define boundary arcs
if (mesh%nb>0) then
   do i=1,mesh%nb-1
      mesh%barc(1,i) = mesh%bnd(i)
      mesh%barc(2,i) = mesh%bnd(i+1)
   end do
   mesh%barc(1,mesh%nb) = mesh%bnd(mesh%nb)
   mesh%barc(2,mesh%nb) = mesh%bnd(1)
end if

! Compute boundary arcs properties
do i=1,mesh%nb
   mesh%barc_lon(:,i) = (/mesh%lon(mesh%barc(1,i)),mesh%lon(mesh%barc(2,i))/)
   mesh%barc_lat(:,i) = (/mesh%lat(mesh%barc(1,i)),mesh%lat(mesh%barc(2,i))/)
   call sphere_dist(mesh%barc_lon(1,i),mesh%barc_lat(1,i),mesh%barc_lon(2,i),mesh%barc_lat(2,i),mesh%barc_dist(i))
   v1 = (/mesh%x(mesh%barc(1,i)),mesh%y(mesh%barc(1,i)),mesh%z(mesh%barc(1,i))/)
   v2 = (/mesh%x(mesh%barc(2,i)),mesh%y(mesh%barc(2,i)),mesh%z(mesh%barc(2,i))/)
   call vector_product(v1,v2,mesh%barc_vp(:,i))
end do

if (lbdist) then
   ! Find minimal distance to a boundary arc
   do i=1,mesh%n
      call mesh%find_bdist(mpl,mesh%lon(i),mesh%lat(i),mesh%bdist(i))
   end do
else
   ! Missing
   mesh%bdist = mpl%msv%valr
end if

end subroutine mesh_bnodes

!----------------------------------------------------------------------
! Subroutine: mesh_find_bdist
! Purpose: find shortest distance to boundary arcs
!----------------------------------------------------------------------
subroutine mesh_find_bdist(mesh,mpl,lon,lat,bdist)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh  ! Mesh
type(mpl_type),intent(inout) :: mpl  ! MPI data
real(kind_real),intent(in) :: lon    ! Longitude
real(kind_real),intent(in) :: lat    ! Latitude
real(kind_real),intent(out) :: bdist ! Distance to boundary

! Local variables
integer :: i
real(kind_real) :: v(3),vf(3),vt(3),tlat,tlon,dist_t1,dist_t2
character(len=1024),parameter :: subr = 'mesh_find_bdist'

! Check
if (mpl%msv%isi(mesh%nb)) call mpl%abort(subr,'boundary arcs have not been computed')

! Initialization
bdist = pi

if (mesh%nb>0) then
   ! Transform to cartesian coordinates
   call lonlat2xyz(mpl,lon,lat,v(1),v(2),v(3))
end if

! Compute the shortest distance from each boundary arc great-circle
do i=1,mesh%nb
   ! Vector products
   call vector_product(v,mesh%barc_vp(:,i),vf)
   call vector_product(mesh%barc_vp(:,i),vf,vt)

   ! Back to spherical coordinates
   call xyz2lonlat(mpl,vt(1),vt(2),vt(3),tlon,tlat)

   ! Check whether T is on the arc
   call sphere_dist(tlon,tlat,mesh%barc_lon(1,i),mesh%barc_lat(1,i),dist_t1)
   call sphere_dist(tlon,tlat,mesh%barc_lon(2,i),mesh%barc_lat(2,i),dist_t2)
   if ((dist_t1<mesh%barc_dist(i)).and.(dist_t2<mesh%barc_dist(i))) then
      ! T is on the arc
      call sphere_dist(lon,lat,tlon,tlat,dist_t1)
      bdist = min(bdist,dist_t1)
   else
      ! T is not on the arc
      call sphere_dist(lon,lat,mesh%barc_lon(1,i),mesh%barc_lat(1,i),dist_t1)
      call sphere_dist(lon,lat,mesh%barc_lon(2,i),mesh%barc_lat(2,i),dist_t2)
      bdist = min(bdist,min(dist_t1,dist_t2))
   end if
end do

end subroutine mesh_find_bdist

!----------------------------------------------------------------------
! Subroutine: mesh_check
! Purpose: check whether the mesh is made of counter-clockwise triangles
!----------------------------------------------------------------------
subroutine mesh_check(mesh,mpl,valid,fix)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh
type(mpl_type),intent(inout) :: mpl    ! MPI data
logical,intent(out) :: valid(mesh%n)   ! Validity flag
logical,intent(in),optional :: fix     ! Fix mesh flag

! Local variables
integer :: it,i,iend,navg,j,ind_avg(40),ii
real(kind_real),allocatable :: a(:),b(:),c(:),cd(:),cp(:)
logical :: validt(mesh%nt),lfix,init
character(len=1024),parameter :: subr = 'mesh_check'

!$omp parallel do schedule(static) private(it) firstprivate(a,b,c,cd,cp)
do it=1,mesh%nt
   ! Allocation
   allocate(a(3))
   allocate(b(3))
   allocate(c(3))
   allocate(cd(3))
   allocate(cp(3))

   ! Check vertices status
   if (mpl%msv%isallnotr(mesh%x(mesh%ltri(:,it))).and.mpl%msv%isallnotr(mesh%y(mesh%ltri(:,it))) &
 & .and.mpl%msv%isallnotr(mesh%z(mesh%ltri(:,it)))) then
      ! Vertices
      a = (/mesh%x(mesh%ltri(1,it)),mesh%y(mesh%ltri(1,it)),mesh%z(mesh%ltri(1,it))/)
      b = (/mesh%x(mesh%ltri(2,it)),mesh%y(mesh%ltri(2,it)),mesh%z(mesh%ltri(2,it))/)
      c = (/mesh%x(mesh%ltri(3,it)),mesh%y(mesh%ltri(3,it)),mesh%z(mesh%ltri(3,it))/)

      ! Cross-product (c-b)x(a-b)
      call vector_product(c-b,a-b,cp)

      ! Centroid
      cd = (a+b+c)/3.0

      ! Compare the directions
      validt(it) = sum(cp*cd)>0.0
   else
      ! At least one vertex ic1 is missing
      validt(it) = .false.
   end if

   ! Release memory
   deallocate(a)
   deallocate(b)
   deallocate(c)
   deallocate(cd)
   deallocate(cp)
end do
!$omp end parallel do

! Check vertices
valid = .true.
do it=1,mesh%nt
   if (.not.validt(it)) valid(mesh%ltri(:,it)) = .false.
end do

! Set fix parameter
lfix = .false.
if (present(fix)) lfix = fix
lfix = (lfix.and.any(.not.validt))

! Fix mesh
if (lfix) then
   do it=1,mesh%nt
      if (.not.validt(it)) then
         do ii=1,3
            i = mesh%ltri(ii,it)
            iend = mesh%lend(i)
            init = .true.
            navg = 0
            do while ((iend/=mesh%lend(i)).or.init)
               j = abs(mesh%list(iend))
               if (valid(j)) then
                  navg = navg+1
                  if (navg>40) call mpl%abort(subr,'max. number of neighbors should be increased')
                  ind_avg(navg) = j
               end if
               iend = mesh%lptr(iend)
               init = .false.
            end do

            if (navg>0) then
               ! Average in cartesian coordinates
               mesh%x(i) = sum(mesh%x(ind_avg(1:navg)))/real(navg,kind_real)
               mesh%y(i) = sum(mesh%y(ind_avg(1:navg)))/real(navg,kind_real)
               mesh%z(i) = sum(mesh%z(ind_avg(1:navg)))/real(navg,kind_real)

               ! Back to spherical coordinates
               call xyz2lonlat(mpl,mesh%x(i),mesh%y(i),mesh%z(i),mesh%lon(i),mesh%lat(i))
            end if
         end do

         ! Allocation
         allocate(a(3))
         allocate(b(3))
         allocate(c(3))
         allocate(cd(3))
         allocate(cp(3))

         ! Check vertices status
         if (mpl%msv%isallnotr(mesh%x(mesh%ltri(:,it))).and.mpl%msv%isallnotr(mesh%y(mesh%ltri(:,it))) &
       & .and.mpl%msv%isallnotr(mesh%z(mesh%ltri(:,it)))) then
            ! Vertices
            a = (/mesh%x(mesh%ltri(1,it)),mesh%y(mesh%ltri(1,it)),mesh%z(mesh%ltri(1,it))/)
            b = (/mesh%x(mesh%ltri(2,it)),mesh%y(mesh%ltri(2,it)),mesh%z(mesh%ltri(2,it))/)
            c = (/mesh%x(mesh%ltri(3,it)),mesh%y(mesh%ltri(3,it)),mesh%z(mesh%ltri(3,it))/)

            ! Cross-product (c-b)x(a-b)
            call vector_product(c-b,a-b,cp)

            ! Centroid
            cd = (a+b+c)/3.0

            ! Compare the directions
            validt(it) = sum(cp*cd)>0.0
         else
            ! At least one vertex ic1 ismissing
            validt(it) = .false.
         end if

         if (validt(it)) then
            write(mpl%info,'(a19,a,i6,a)') '','Triangle ',it,' of the mesh has been fixed sucessfully'
         else
            write(mpl%info,'(a19,a,i6,a)') '','Triangle ',it,' of the mesh failed to be fixed'
            call mpl%abort(subr,'mesh%fix failed')
         end if

         ! Release memory
         deallocate(a)
         deallocate(b)
         deallocate(c)
         deallocate(cd)
         deallocate(cp)
      end if
   end do

   ! Check vertices
   valid = .true.
   do it=1,mesh%nt
      if (.not.validt(it)) valid(mesh%ltri(:,it)) = .false.
   end do
end if

end subroutine mesh_check

!----------------------------------------------------------------------
! Subroutine: mesh_inside
! Purpose: find whether a point is inside the mesh
!----------------------------------------------------------------------
subroutine mesh_inside(mesh,mpl,lon,lat,inside_mesh)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh ! Mesh
type(mpl_type),intent(inout) :: mpl ! MPI data
real(kind_real),intent(in) :: lon   ! Longitude
real(kind_real),intent(in) :: lat   ! Latitude
logical,intent(out) :: inside_mesh  ! True if the point is inside the mesh

! Local variables
integer :: info
real(kind_real) :: p(3)

if (mesh%nb>0) then
   ! Transform to cartesian coordinates
   call lonlat2xyz(mpl,lon,lat,p(1),p(2),p(3))

   ! Find whether the point is inside the convex hull
   inside_mesh = inside(p,mesh%n,mesh%x,mesh%y,mesh%z,mesh%nb,mesh%bnd,info)
else
   ! No boundary
   inside_mesh = .true.
end if

end subroutine mesh_inside

!----------------------------------------------------------------------
! Subroutine: mesh_barycentric
! Purpose: compute barycentric coordinates
!----------------------------------------------------------------------
subroutine mesh_barycentric(mesh,mpl,lon,lat,istart,b,ib)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh ! Mesh
type(mpl_type),intent(inout) :: mpl ! MPI data
real(kind_real),intent(in) :: lon   ! Longitude
real(kind_real),intent(in) :: lat   ! Latitude
integer,intent(in) :: istart        ! Starting index
real(kind_real),intent(out) :: b(3) ! Barycentric weights
integer,intent(out) :: ib(3)        ! Barycentric indices

! Local variables
real(kind_real) :: p(3)

! Transform to cartesian coordinates
call lonlat2xyz(mpl,lon,lat,p(1),p(2),p(3))

! Compute barycentric coordinates
b = 0.0
ib = 0
call trfind(istart,p,mesh%n,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,b(1),b(2),b(3),ib(1),ib(2),ib(3))

end subroutine mesh_barycentric

!----------------------------------------------------------------------
! Subroutine: mesh_addnode
! Purpose: add node to a mesh
!----------------------------------------------------------------------
subroutine mesh_addnode(mesh,mpl,lonnew,latnew)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh
type(mpl_type),intent(inout) :: mpl    ! MPI data
real(kind_real),intent(in) :: lonnew   ! Longitude
real(kind_real),intent(in) :: latnew   ! Latitude

! Local variables
integer :: info
integer,allocatable :: order(:),order_inv(:),list(:),lptr(:),lend(:)
real(kind_real),allocatable :: lon(:),lat(:),x(:),y(:),z(:)

! Allocation
allocate(order(mesh%n))
allocate(order_inv(mesh%n))
allocate(lon(mesh%n))
allocate(lat(mesh%n))
allocate(x(mesh%n))
allocate(y(mesh%n))
allocate(z(mesh%n))
allocate(list(6*(mesh%n-2)))
allocate(lptr(6*(mesh%n-2)))
allocate(lend(mesh%n))

! Copy
order = mesh%order
order_inv = mesh%order_inv
lon = mesh%lon
lat = mesh%lat
x = mesh%x
y = mesh%y
z = mesh%z
list = mesh%list
lptr = mesh%lptr
lend = mesh%lend

! Release memory
call mesh%dealloc

! Reallocation
mesh%n = mesh%n+1
allocate(mesh%order(mesh%n))
allocate(mesh%order_inv(mesh%n))
allocate(mesh%lon(mesh%n))
allocate(mesh%lat(mesh%n))
allocate(mesh%x(mesh%n))
allocate(mesh%y(mesh%n))
allocate(mesh%z(mesh%n))
allocate(mesh%list(6*(mesh%n-2)))
allocate(mesh%lptr(6*(mesh%n-2)))
allocate(mesh%lend(mesh%n))

! Copy
mesh%order(1:mesh%n-1) = order
mesh%order_inv(1:mesh%n-1) = order_inv
mesh%lon(1:mesh%n-1) = lon
mesh%lat(1:mesh%n-1) = lat
mesh%x(1:mesh%n-1) = x
mesh%y(1:mesh%n-1) = y
mesh%z(1:mesh%n-1) = z
mesh%list(1:6*(mesh%n-3)) = list
mesh%lptr(1:6*(mesh%n-3)) = lptr
mesh%lend(1:mesh%n-1) = lend

! Compute new element coordinates
call lonlat2xyz(mpl,lonnew,latnew,mesh%x(mesh%n),mesh%y(mesh%n),mesh%z(mesh%n))

! Update mesh
mesh%order(mesh%n) = mesh%n
mesh%order_inv(mesh%n) = mesh%n
mesh%lon(mesh%n) = lonnew
mesh%lat(mesh%n) = latnew
call addnod(mpl,1,mesh%n,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,mesh%lnew,info)

! Release memory
deallocate(order)
deallocate(order_inv)
deallocate(lon)
deallocate(lat)
deallocate(x)
deallocate(y)
deallocate(z)
deallocate(list)
deallocate(lptr)
deallocate(lend)

end subroutine mesh_addnode

!----------------------------------------------------------------------
! Subroutine: mesh_polygon
! Purpose: compute polygon area
!----------------------------------------------------------------------
subroutine mesh_polygon(mesh,np,plist,area_polygon)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh             ! Mesh
integer,intent(in) :: np                        ! Number of points
integer,intent(in) :: plist(np)                 ! List of indices
real(kind_real),intent(out) :: area_polygon(np) ! Area

! Local variables
integer :: ip,i_src,index_triangle,i_src_last,i_src_new,i_src_stop,vertex_last,vertex_new,nb,info
integer :: lnew_tmp,lptr_tmp(6*(mesh%n-2)),lbtri(6,mesh%n),listc(6*(mesh%n-2))
real(kind_real) :: area_triangle,v1(3),v2(3),v3(3)
real(kind_real) :: xc(6*(mesh%n-2)),yc(6*(mesh%n-2)),zc(6*(mesh%n-2)),rc(6*(mesh%n-2))
logical :: loop

! Copy
lptr_tmp = mesh%lptr
lnew_tmp = mesh%lnew

! Compute Voronoi polygons
call crlist(mesh%n,mesh%n,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lend,lptr_tmp,lnew_tmp,lbtri,listc,nb,xc,yc,zc,rc,info)

do ip=1,np
   i_src = plist(ip)
   area_polygon(ip) = 0.0
   index_triangle = 0
   i_src_stop = mesh%lend(i_src)
   i_src_new = i_src_stop
   vertex_new = listc(i_src_new)
   loop = .true.
   do while (loop)
      index_triangle = index_triangle+1
      i_src_last = i_src_new
      i_src_new = lptr_tmp(i_src_last)
      vertex_last = vertex_new
      vertex_new = listc(i_src_new)
      v1 = (/mesh%x(i_src),mesh%y(i_src),mesh%z(i_src)/)
      v2 = (/xc(vertex_last),yc(vertex_last),zc(vertex_last)/)
      v3 = (/xc(vertex_new),yc(vertex_new),zc(vertex_new)/)
      area_triangle = areas(v1,v2,v3)
      area_polygon(ip) = area_polygon(ip)+area_triangle
      loop = (i_src_new/=i_src_stop)
   end do
end do

end subroutine mesh_polygon

end module type_mesh
