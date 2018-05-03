!----------------------------------------------------------------------
! Module: type_mesh
!> Purpose: mesh derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_mesh

!$ use omp_lib
use tools_const, only: req
use tools_display, only: msgerror,msgwarning
use tools_func, only: sphere_dist,vector_product
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr,isallnotmsr
use tools_stripack, only: addnod,areas,bnodes,crlist,inside,scoord,trans,trfind,trlist,trmesh
use type_mpl, only: mpl
use type_rng, only: rng

implicit none

! Mesh derived type
type mesh_type
   ! Mesh structure
   integer :: n                            !< Number of points
   integer,allocatable :: order(:)         !< Order of shuffled points
   integer,allocatable :: order_inv(:)     !< Inverse order of shuffled points
   real(kind_real),allocatable :: lon(:)   !< Points longitudes
   real(kind_real),allocatable :: lat(:)   !< Points latitudes
   real(kind_real),allocatable :: x(:)     !< x-coordinate
   real(kind_real),allocatable :: y(:)     !< y-coordinate
   real(kind_real),allocatable :: z(:)     !< z-coordinate
   integer,allocatable :: list(:)          !< Stripack list
   integer,allocatable :: lptr(:)          !< Stripack list pointer
   integer,allocatable :: lend(:)          !< Stripack list end
   integer :: lnew                         !< Stripack pointer to the first empty location in list
   integer :: nb                           !< Number of boundary nodes
   integer,allocatable :: bnd(:)           !< Boundary nodes

   ! Triangles data
   integer :: nt                           !< Number of triangles
   integer :: na                           !< Number of arcs
   integer,allocatable :: ltri(:,:)        !< Triangles indices
   integer,allocatable :: larc(:,:)        !< Arcs indices
   real(kind_real),allocatable :: bdist(:) !< Distance to the closest boundary arc
contains
   procedure :: create => mesh_create
   procedure :: dealloc => mesh_dealloc
   procedure :: copy => mesh_copy
   procedure :: trans => mesh_trans
   procedure :: trlist => mesh_trlist
   procedure :: bnodes => mesh_bnodes
   procedure :: barcs => mesh_barcs
   procedure :: check => mesh_check
   procedure :: inside => mesh_inside
   procedure :: barycentric
   procedure :: addnode
   procedure :: polygon
end type mesh_type

logical,parameter :: shuffle = .true. !< Shuffle mesh order (more efficient to compute the Delaunay triangulation)

private
public :: mesh_type

contains

!----------------------------------------------------------------------
! Subroutine: mesh_create
!> Purpose: create mesh
!----------------------------------------------------------------------
subroutine mesh_create(mesh,n,lon,lat)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh !< Mesh
integer,intent(in) :: n                !< Mesh size
real(kind_real),intent(in) :: lon(n)   !< Longitudes
real(kind_real),intent(in) :: lat(n)   !< Latitudes

! Local variables
integer :: i,k,info
integer,allocatable :: jtab(:),near(:),next(:)
real(kind_real),allocatable :: dist(:)

! Allocation
mesh%n = n

! Allocation
allocate(mesh%order(mesh%n))
allocate(mesh%order_inv(mesh%n))
allocate(mesh%lon(mesh%n))
allocate(mesh%lat(mesh%n))
allocate(mesh%list(6*(mesh%n-2)))
allocate(mesh%lptr(6*(mesh%n-2)))
allocate(mesh%lend(mesh%n))
allocate(mesh%x(mesh%n))
allocate(mesh%y(mesh%n))
allocate(mesh%z(mesh%n))
allocate(near(mesh%n))
allocate(next(mesh%n))
allocate(dist(mesh%n))

! Points order
do i=1,mesh%n
   mesh%order(i) = i
end do

if (shuffle) then
   ! Shuffle order (more efficient to compute the Delaunay triangulation)
   allocate(jtab(mesh%n))
   if (mpl%main) call rng%rand_integer(1,mesh%n,jtab)
   call mpl%bcast(jtab)
   do i=mesh%n,2,-1
      k = mesh%order(jtab(i))
      mesh%order(jtab(i)) = mesh%order(i)
      mesh%order(i) = k
   end do
end if

! Restrictive inverse order
call msi(mesh%order_inv)
do i=1,mesh%n
   mesh%order_inv(mesh%order(i)) = i
end do

! Transform to cartesian coordinates
call mesh%trans(lon,lat)

! Create mesh
mesh%list = 0
call trmesh(mesh%n,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,mesh%lnew,near,next,dist,info)

end subroutine mesh_create

!----------------------------------------------------------------------
! Subroutine: mesh_dealloc
!> Purpose: deallocate mesh
!----------------------------------------------------------------------
subroutine mesh_dealloc(mesh)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh !< Mesh

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
if (allocated(mesh%ltri)) deallocate(mesh%ltri)
if (allocated(mesh%larc)) deallocate(mesh%larc)
if (allocated(mesh%bdist)) deallocate(mesh%bdist)

end subroutine mesh_dealloc

!----------------------------------------------------------------------
! Function: mesh_copy
!> Purpose: copy mesh
!----------------------------------------------------------------------
type(mesh_type) function mesh_copy(mesh)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh !< Input mesh

! Copy sizes
mesh_copy%n = mesh%n
if (allocated(mesh%ltri)) then
   mesh_copy%nt = mesh%nt
   mesh_copy%na = mesh%na
end if

! Release memory
call mesh_copy%dealloc

! Allocation
allocate(mesh_copy%order(mesh_copy%n))
allocate(mesh_copy%order_inv(mesh_copy%n))
allocate(mesh_copy%lon(mesh_copy%n))
allocate(mesh_copy%lat(mesh_copy%n))
allocate(mesh_copy%x(mesh_copy%n))
allocate(mesh_copy%y(mesh_copy%n))
allocate(mesh_copy%z(mesh_copy%n))
allocate(mesh_copy%list(6*(mesh_copy%n-2)))
allocate(mesh_copy%lptr(6*(mesh_copy%n-2)))
allocate(mesh_copy%lend(mesh_copy%n))
if (allocated(mesh%ltri)) then
   allocate(mesh_copy%ltri(3,mesh_copy%nt))
   allocate(mesh_copy%larc(2,mesh_copy%na))
   allocate(mesh_copy%bdist(mesh_copy%n))
end if

! Copy data
mesh_copy%order = mesh%order
mesh_copy%order_inv = mesh%order_inv
mesh_copy%lon = mesh%lon
mesh_copy%lat = mesh%lat
mesh_copy%x = mesh%x
mesh_copy%y = mesh%y
mesh_copy%z = mesh%z
mesh_copy%list = mesh%list
mesh_copy%lptr = mesh%lptr
mesh_copy%lend = mesh%lend
mesh_copy%lnew = mesh%lnew
if (allocated(mesh%ltri)) then
   mesh_copy%nb = mesh%nb
   mesh_copy%ltri = mesh%ltri
   mesh_copy%larc = mesh%larc
   mesh_copy%bdist = mesh%bdist
end if

end function mesh_copy

!----------------------------------------------------------------------
! Subroutine: mesh_trans
!> Purpose: transform to cartesian coordinates
!----------------------------------------------------------------------
subroutine mesh_trans(mesh,lon,lat)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh    !< Mesh
real(kind_real),intent(in) :: lon(mesh%n) !< Longitude
real(kind_real),intent(in) :: lat(mesh%n) !< Latitude

! Copy lon/lat
mesh%lon = lon(mesh%order)
mesh%lat = lat(mesh%order)

! Transform to cartesian coordinates
call trans(mesh%n,mesh%lat,mesh%lon,mesh%x,mesh%y,mesh%z)

end subroutine mesh_trans

!----------------------------------------------------------------------
! Subroutine: mesh_trlist
!> Purpose: compute triangle list, arc list
!----------------------------------------------------------------------
subroutine mesh_trlist(mesh)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh !< Mesh

! Local variables
integer :: info
integer,allocatable :: ltri(:,:)
integer :: ia,it,i,i1,i2

! Allocation
allocate(ltri(9,2*(mesh%n-2)))

! Create triangles list
call trlist(mesh%n,mesh%list,mesh%lptr,mesh%lend,9,mesh%nt,ltri,info)

! Allocation
mesh%na = maxval(ltri(7:9,1:mesh%nt))
allocate(mesh%ltri(3,mesh%nt))
allocate(mesh%larc(2,mesh%na))

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

end subroutine mesh_trlist

!----------------------------------------------------------------------
! Subroutine: mesh_bnodes
!> Purpose: find boundary nodes
!----------------------------------------------------------------------
subroutine mesh_bnodes(mesh)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh !< Mesh

! Allocation
allocate(mesh%bnd(mesh%n))

! Find boundary nodes
call msi(mesh%bnd)
call bnodes(mesh%n,mesh%list,mesh%lptr,mesh%lend,mesh%bnd,mesh%nb,mesh%na,mesh%nt)

end subroutine mesh_bnodes

!----------------------------------------------------------------------
! Subroutine: mesh_barcs
!> Purpose: find boundary arcs
!----------------------------------------------------------------------
subroutine mesh_barcs(mesh)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh !< Mesh

! Local variables
integer :: i
integer,allocatable :: larcb(:,:)
integer :: ia,nab,iab
real(kind_real) :: dist_12,v1(3),v2(3),vp(3),v(3),vf(3),vt(3),tlat,tlon,trad,dist_t1,dist_t2

! Allocation
allocate(larcb(2,3*(mesh%n-2)))
allocate(mesh%bdist(mesh%n))

if (mesh%nb>0) then
   ! Find boundary arcs
   nab = 0
   do ia=1,mesh%na
      if (any(mesh%bnd(1:mesh%nb)==mesh%larc(1,ia)).and.any(mesh%bnd(1:mesh%nb)==mesh%larc(2,ia))) then
         nab = nab+1
         larcb(:,nab) = mesh%larc(:,ia)
      end if
   end do

   ! Find minimal distance to a boundary arc
   mesh%bdist = huge(1.0)
   do iab=1,nab
      ! Distance
      call sphere_dist(mesh%lon(larcb(1,iab)),mesh%lat(larcb(1,iab)), &
    & mesh%lon(larcb(2,iab)),mesh%lat(larcb(2,iab)),dist_12)

      ! Vectors
      v1 = (/mesh%x(larcb(1,iab)),mesh%y(larcb(1,iab)),mesh%z(larcb(1,iab))/)
      v2 = (/mesh%x(larcb(2,iab)),mesh%y(larcb(2,iab)),mesh%z(larcb(2,iab))/)

      ! Compute normal vector to the boundary arc plane
      call vector_product(v1,v2,vp)

      ! Compute the shortest distance from each point to the boundary arc great-circle
      do i=1,mesh%n
         ! Vector
         v = (/mesh%x(i),mesh%y(i),mesh%z(i)/)

         ! Vector products
         call vector_product(v,vp,vf)
         call vector_product(vp,vf,vt)

         ! Back to spherical coordinates
         call scoord(vt(1),vt(2),vt(3),tlat,tlon,trad)

         ! Check whether T is on the arc
         call sphere_dist(tlon,tlat,mesh%lon(larcb(1,iab)),mesh%lat(larcb(1,iab)),dist_t1)
         call sphere_dist(tlon,tlat,mesh%lon(larcb(2,iab)),mesh%lat(larcb(2,iab)),dist_t2)
         if ((dist_t1<dist_12).and.(dist_t2<dist_12)) then
            ! T is on the arc
            call sphere_dist(mesh%lon(i),mesh%lat(i),tlon,tlat,dist_t1)
            mesh%bdist(i) = min(mesh%bdist(i),dist_t1)
         else
            ! T is not on the arc
            call sphere_dist(mesh%lon(i),mesh%lat(i),mesh%lon(larcb(1,iab)),mesh%lat(larcb(1,iab)),dist_t1)
            call sphere_dist(mesh%lon(i),mesh%lat(i),mesh%lon(larcb(2,iab)),mesh%lat(larcb(2,iab)),dist_t2)
            mesh%bdist(i) = min(mesh%bdist(i),min(dist_t1,dist_t2))
         end if
      end do
   end do
   mesh%bdist = mesh%bdist/req
else
   mesh%bdist = huge(1.0)
end if

end subroutine mesh_barcs

!----------------------------------------------------------------------
! Subroutine: mesh_check
!> Purpose: check whether the mesh is made of counter-clockwise triangles
!----------------------------------------------------------------------
subroutine mesh_check(mesh,valid)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh       !< Mesh
real(kind_real),intent(out) :: valid(mesh%n) !< Validity flag (1.0 if the vertex is valid, else 0.0)

! Local variables
integer :: it
real(kind_real),allocatable :: a(:),b(:),c(:),cd(:),cp(:)
logical :: validt(mesh%nt)

!$omp parallel do schedule(static) private(it) firstprivate(a,b,c,cd,cp)
do it=1,mesh%nt
   ! Allocation
   allocate(a(3))
   allocate(b(3))
   allocate(c(3))
   allocate(cd(3))
   allocate(cp(3))

   ! Check vertices status
   if (isallnotmsr(mesh%x(mesh%ltri(:,it))).and.isallnotmsr(mesh%y(mesh%ltri(:,it))).and.isallnotmsr(mesh%z(mesh%ltri(:,it)))) then
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
      ! At least one vertex ic1 ms
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
valid = 1.0
do it=1,mesh%nt
   if (.not.validt(it)) valid(mesh%ltri(:,it)) = 0.0
end do

end subroutine mesh_check

!----------------------------------------------------------------------
! Subroutine: mesh_inside
!> Purpose: find whether a point is inside the mesh
!----------------------------------------------------------------------
subroutine mesh_inside(mesh,lon,lat,inside_mesh)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh  !< Mesh
real(kind_real),intent(in) :: lon(1) !< Longitude
real(kind_real),intent(in) :: lat(1) !< Latitude
logical,intent(out) :: inside_mesh   !< True if the point is inside the mesh

! Local variables
integer :: info
real(kind_real) :: p(3)

if (mesh%nb>0) then
   ! Transform to cartesian coordinates
   call trans(1,lat,lon,p(1),p(2),p(3))

   ! Find answer
   inside_mesh = inside(p,mesh%n,mesh%x,mesh%y,mesh%z,mesh%nb,mesh%bnd,info)
else
   ! No boundary
   inside_mesh = .true.
end if

end subroutine mesh_inside

!----------------------------------------------------------------------
! Subroutine: barycentric
!> Purpose: compute barycentric coordinates
!----------------------------------------------------------------------
subroutine barycentric(mesh,lon,lat,istart,b,ib)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh  !< Mesh
real(kind_real),intent(in) :: lon(1) !< Longitude
real(kind_real),intent(in) :: lat(1) !< Latitude
integer,intent(in) :: istart         !< Starting index
real(kind_real),intent(out) :: b(3)  !< Barycentric weights
integer,intent(out) :: ib(3)         !< Barycentric indices

! Local variables
real(kind_real) :: p(3)

! Transform to cartesian coordinates
call trans(1,lat,lon,p(1),p(2),p(3))

! Compute barycentric coordinates
b = 0.0
ib = 0
call trfind(istart,p,mesh%n,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,b(1),b(2),b(3),ib(1),ib(2),ib(3))

end subroutine barycentric

!----------------------------------------------------------------------
! Subroutine: addnode
!> Purpose: add node to a mesh
!----------------------------------------------------------------------
subroutine addnode(mesh,lonnew,latnew)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh  !< Mesh
real(kind_real),intent(in) :: lonnew(1) !< Longitude
real(kind_real),intent(in) :: latnew(1) !< Latitude

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

! Reallocate
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
call trans(1,latnew,lonnew,mesh%x(mesh%n:mesh%n),mesh%y(mesh%n:mesh%n),mesh%z(mesh%n:mesh%n))

! Update mesh
mesh%order(mesh%n) = mesh%n
mesh%order_inv(mesh%n) = mesh%n
mesh%lon(mesh%n) = lonnew(1)
mesh%lat(mesh%n) = latnew(1)
call addnod(1,mesh%n,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,mesh%lnew,info)

end subroutine addnode

!----------------------------------------------------------------------
! Subroutine: polygon
!> Purpose: compute polygon area
!----------------------------------------------------------------------
subroutine polygon(mesh,np,plist,area_polygon)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh             !< Mesh
integer,intent(in) :: np                        !< Number of points
integer,intent(in) :: plist(np)                 !< List of indices
real(kind_real),intent(out) :: area_polygon(np) !< Area

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

end subroutine polygon

end module type_mesh
