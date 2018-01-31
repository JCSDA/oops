!----------------------------------------------------------------------
! Module: type_mesh
!> Purpose: mesh derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_mesh

use omp_lib
use tools_const, only: vector_product
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr,isallnotmsr
use tools_stripack, only: addnod,areas,crlist,trans,trfind,trmesh
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast,mpl_split
use type_randgen, only: rand_integer

implicit none

! Linear operator derived type
type meshtype
   integer :: n                        !< Number of points
   integer,allocatable :: redundant(:) !< Redundant points
   integer :: nnr                      !< Number of non-redundant points
   integer,allocatable :: order(:)     !< Order of shuffled points
   integer,allocatable :: order_inv(:) !< Inverse order of shuffled points
   real(kind_real),allocatable :: x(:) !< x-coordinate
   real(kind_real),allocatable :: y(:) !< y-coordinate
   real(kind_real),allocatable :: z(:) !< z-coordinate
   integer,allocatable :: list(:)      !< Stripack list
   integer,allocatable :: lptr(:)      !< Stripack list pointer
   integer,allocatable :: lend(:)      !< Stripack list end
   integer :: lnew                     !< Stripack pointer to the first empty location in list
end type meshtype

private
public :: meshtype
public :: create_mesh,mesh_dealloc,copy_mesh,check_mesh,barycentric,addnode,polygon

contains

!----------------------------------------------------------------------
! Subroutine: create_mesh
!> Purpose: create mesh
!----------------------------------------------------------------------
subroutine create_mesh(n,lon,lat,lred,mesh)

implicit none

! Passed variables
integer,intent(in) :: n               !< Mesh size
real(kind_real),intent(in) :: lon(n)  !< Longitudes
real(kind_real),intent(in) :: lat(n)  !< Latitudes
logical,intent(in) :: lred            !< Redundant points check
type(meshtype),intent(inout) :: mesh  !< Mesh

! Local variables
integer :: i,j,k,info,iproc
integer :: i_s(mpl%nproc),i_e(mpl%nproc),n_loc(mpl%nproc),i_loc,progint
integer,allocatable :: jtab(:),near(:),next(:),sbuf(:)
real(kind_real),allocatable :: dist(:)
logical,allocatable :: done(:)

! Allocation
mesh%n = n
allocate(mesh%redundant(mesh%n))
call msi(mesh%redundant)

! Look for redundant points
if (lred.and..false.) then
   ! MPI splitting
   call mpl_split(n,i_s,i_e,n_loc)

   ! Allocation
   allocate(done(n_loc(mpl%myproc)))
   allocate(sbuf(n_loc(mpl%myproc)))

   ! Loop over points
   write(mpl%unit,'(a7,a)',advance='no') '','Look for redundant points: '
   call prog_init(progint,done)
   call msi(sbuf)
   !$omp parallel do schedule(static) private(i_loc,i,j)
   do i_loc=1,n_loc(mpl%myproc)
      i = i_s(mpl%myproc)+i_loc-1
      do j=1,i-1
         if (.not.(abs(lon(i)-lon(j))>0.0).and..not.(abs(lat(i)-lat(j))>0.0)) then
            sbuf(i_loc) = j
            exit
         end if
      end do
      done(i_loc) = .true.
      call prog_print(progint,done)
   end do
   !$omp end parallel do
   write(mpl%unit,'(a)') '100%'

   ! Communication
   if (mpl%main) then
      do iproc=1,mpl%nproc
         if (n_loc(iproc)>0) then
            if (iproc==mpl%ioproc) then
               ! Copy data
               mesh%redundant(i_s(iproc):i_e(iproc)) = sbuf
            else
               ! Receive data on ioproc
               call mpl_recv(n_loc(iproc),mesh%redundant(i_s(iproc):i_e(iproc)),iproc,mpl%tag)
            end if
         end if
      end do
   else
      if (n_loc(mpl%myproc)>0) then
         ! Send data to ioproc
         call mpl_send(n_loc(mpl%myproc),sbuf,mpl%ioproc,mpl%tag)
      end if
   end if
   mpl%tag = mpl%tag+1

   ! Broadcast data
   call mpl_bcast(mesh%redundant,mpl%ioproc)

   ! Check for successive redundant points
   do i=1,mesh%n
      if (isnotmsi(mesh%redundant(i))) then
         do while (isnotmsi(mesh%redundant(mesh%redundant(i))))
            mesh%redundant(i) = mesh%redundant(mesh%redundant(i))
         end do
      end if
   end do
end if
mesh%nnr = count(.not.isnotmsi(mesh%redundant))

! Allocation
allocate(mesh%order(mesh%nnr))
allocate(mesh%order_inv(mesh%n))
allocate(mesh%list(6*(mesh%nnr-2)))
allocate(mesh%lptr(6*(mesh%nnr-2)))
allocate(mesh%lend(mesh%nnr))
allocate(mesh%x(mesh%nnr))
allocate(mesh%y(mesh%nnr))
allocate(mesh%z(mesh%nnr))
allocate(near(mesh%nnr))
allocate(next(mesh%nnr))
allocate(jtab(mesh%nnr))
allocate(dist(mesh%nnr))

! Shuffle arrays (more efficient to compute the Delaunay triangulation)
i = 0
do j=1,mesh%n
   if (.not.isnotmsi(mesh%redundant(j))) then
      i = i+1
      mesh%order(i) = j
   end if
end do
if (mpl%main) call rand_integer(1,mesh%nnr,jtab)
call mpl_bcast(jtab,mpl%ioproc)
do i=mesh%nnr,2,-1
   k = mesh%order(jtab(i))
   mesh%order(jtab(i)) = mesh%order(i)
   mesh%order(i) = k
end do

! Inverse order
call msi(mesh%order_inv)
do i=1,mesh%nnr
   mesh%order_inv(mesh%order(i)) = i
end do

! Transform to cartesian coordinates
call trans(mesh%nnr,lat(mesh%order),lon(mesh%order),mesh%x,mesh%y,mesh%z)

! Create mesh
mesh%list = 0
call trmesh(mesh%nnr,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,mesh%lnew,near,next,dist,info)

end subroutine create_mesh

!----------------------------------------------------------------------
! Subroutine: mesh_dealloc
!> Purpose: deallocate mesh
!----------------------------------------------------------------------
subroutine mesh_dealloc(mesh)

implicit none

! Passed variables
type(meshtype),intent(inout) :: mesh !< Mesh

! Release memory
if (allocated(mesh%redundant)) deallocate(mesh%redundant)
if (allocated(mesh%x)) deallocate(mesh%x)
if (allocated(mesh%y)) deallocate(mesh%y)
if (allocated(mesh%z)) deallocate(mesh%z)
if (allocated(mesh%list)) deallocate(mesh%list)
if (allocated(mesh%lptr)) deallocate(mesh%lptr)
if (allocated(mesh%lend)) deallocate(mesh%lend)
if (allocated(mesh%order)) deallocate(mesh%order)
if (allocated(mesh%order_inv)) deallocate(mesh%order_inv)

end subroutine mesh_dealloc

!----------------------------------------------------------------------
! Subroutine: copy_mesh
!> Purpose: copy mesh
!----------------------------------------------------------------------
subroutine copy_mesh(mesh_in,mesh_out)

implicit none

! Passed variables
type(meshtype),intent(in) :: mesh_in   !< Input mesh
type(meshtype),intent(out) :: mesh_out !< Output mesh

! Copy sizes
mesh_out%n = mesh_in%n
mesh_out%nnr = mesh_in%nnr

! Release memory
if (allocated(mesh_out%redundant)) deallocate(mesh_out%redundant)
if (allocated(mesh_out%order)) deallocate(mesh_out%order)
if (allocated(mesh_out%order_inv)) deallocate(mesh_out%order_inv)
if (allocated(mesh_out%x)) deallocate(mesh_out%x)
if (allocated(mesh_out%y)) deallocate(mesh_out%y)
if (allocated(mesh_out%z)) deallocate(mesh_out%z)
if (allocated(mesh_out%list)) deallocate(mesh_out%list)
if (allocated(mesh_out%lptr)) deallocate(mesh_out%lptr)
if (allocated(mesh_out%lend)) deallocate(mesh_out%lend)

! Allocation
allocate(mesh_out%redundant(mesh_out%n))
allocate(mesh_out%order(mesh_out%nnr))
allocate(mesh_out%order_inv(mesh_out%n))
allocate(mesh_out%x(mesh_out%nnr))
allocate(mesh_out%y(mesh_out%nnr))
allocate(mesh_out%z(mesh_out%nnr))
allocate(mesh_out%list(6*(mesh_out%nnr-2)))
allocate(mesh_out%lptr(6*(mesh_out%nnr-2)))
allocate(mesh_out%lend(mesh_out%nnr))

! Copy data
mesh_out%redundant = mesh_in%redundant
mesh_out%order = mesh_in%order
mesh_out%order_inv = mesh_in%order_inv
mesh_out%x = mesh_in%x
mesh_out%y = mesh_in%y
mesh_out%z = mesh_in%z
mesh_out%list = mesh_in%list
mesh_out%lptr = mesh_in%lptr
mesh_out%lend = mesh_in%lend
mesh_out%lnew = mesh_in%lnew

end subroutine copy_mesh

!----------------------------------------------------------------------
! Subroutine: check_mesh
!> Purpose: check whether the mesh is made of counter-clockwise triangles
!----------------------------------------------------------------------
subroutine check_mesh(n,lon,lat,nt,ltri,valid)

implicit none

! Passed variables
integer,intent(in) :: n                 !< Mesh size
real(kind_real),intent(in) :: lon(n)    !< Longitudes
real(kind_real),intent(in) :: lat(n)    !< Latitudes
integer,intent(in) :: nt                !< Number of triangles
integer,intent(in) :: ltri(3,nt)        !< Stripack triangles list
real(kind_real),intent(out) :: valid(n) !< 1.0 if the vertex ic1 valid, else 0.0

! Local variables
integer :: it
real(kind_real) :: x(n),y(n),z(n)
real(kind_real),allocatable :: a(:),b(:),c(:),cd(:),v1(:),v2(:),cp(:)
logical :: validt(nt)

! Transform to cartesian coordinates
call trans(n,lat,lon,x,y,z)

!$omp parallel do schedule(static) private(it,a,b,c,cd,v1,v2,cp)
do it=1,nt
   ! Check vertices status
   if (isallnotmsr(x(ltri(:,it))).and.isallnotmsr(y(ltri(:,it))).and.isallnotmsr(z(ltri(:,it)))) then
      ! Allocation
      allocate(a(3))
      allocate(b(3))
      allocate(c(3))
      allocate(cd(3))
      allocate(v1(3))
      allocate(v2(3))
      allocate(cp(3))

      ! Vertices
      a = (/x(ltri(1,it)),y(ltri(1,it)),z(ltri(1,it))/)
      b = (/x(ltri(2,it)),y(ltri(2,it)),z(ltri(2,it))/)
      c = (/x(ltri(3,it)),y(ltri(3,it)),z(ltri(3,it))/)

      ! Cross-product (c-b)x(a-b)
      call vector_product(c-b,a-b,cp)

      ! Centroid
      cd = (a+b+c)/3.0

      ! Compare the directions
      validt(it) = dot_product(cp,cd)>0.0

      ! Release memory
      deallocate(a)
      deallocate(b)
      deallocate(c)
      deallocate(cd)
      deallocate(v1)
      deallocate(v2)
      deallocate(cp)
   else
      ! At least one vertex ic1 ms
      validt(it) = .false.
   end if
end do
!$omp end parallel do

! Check vertices
valid = 1.0
do it=1,nt
   if (.not.validt(it)) valid(ltri(:,it)) = 0.0
end do

end subroutine check_mesh

!----------------------------------------------------------------------
! Subroutine: barycentric
!> Purpose: compute barycentric coordinates
!----------------------------------------------------------------------
subroutine barycentric(lon,lat,mesh,istart,b,ib)

implicit none

! Passed variables
real(kind_real),intent(in) :: lon(1) !< Longitude
real(kind_real),intent(in) :: lat(1) !< Latitude
type(meshtype),intent(in) :: mesh    !< Mesh
integer,intent(in) :: istart         !< Starting index
real(kind_real),intent(out) :: b(3)  !< Barycentric weights
integer,intent(out) :: ib(3)         !< Barycentric indices

! Local variables
real(kind_real) :: p(3)

! Transform to cartesian coordinates
call trans(1,lat,lon,p(1),p(2),p(3))

! Compute barycentric coordinates
call trfind(istart,dble(p),mesh%nnr,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,b(1),b(2),b(3),ib(1),ib(2),ib(3))

end subroutine barycentric

!----------------------------------------------------------------------
! Subroutine: addnode
!> Purpose: add node to a mesh
!----------------------------------------------------------------------
subroutine addnode(lon,lat,mesh)

implicit none

! Passed variables
real(kind_real),intent(in) :: lon(1) !< Longitude
real(kind_real),intent(in) :: lat(1) !< Latitude
type(meshtype),intent(inout) :: mesh !< Mesh

! Local variables
integer :: info
integer,allocatable :: order(:),order_inv(:),list(:),lptr(:),lend(:)
real(kind_real),allocatable :: x(:),y(:),z(:)

! Allocation
allocate(order(mesh%nnr))
allocate(order_inv(mesh%n))
allocate(x(mesh%nnr))
allocate(y(mesh%nnr))
allocate(z(mesh%nnr))
allocate(list(6*(mesh%nnr-2)))
allocate(lptr(6*(mesh%nnr-2)))
allocate(lend(mesh%nnr))

! Copy
order = mesh%order
order_inv = mesh%order_inv
x = mesh%x
y = mesh%y
z = mesh%z
list = mesh%list
lptr = mesh%lptr
lend = mesh%lend

! Reallocation
deallocate(mesh%order)
deallocate(mesh%order_inv)
deallocate(mesh%x)
deallocate(mesh%y)
deallocate(mesh%z)
deallocate(mesh%list)
deallocate(mesh%lptr)
deallocate(mesh%lend)
mesh%n = mesh%n+1
mesh%nnr = mesh%nnr+1
allocate(mesh%order(mesh%nnr))
allocate(mesh%order_inv(mesh%n))
allocate(mesh%x(mesh%nnr))
allocate(mesh%y(mesh%nnr))
allocate(mesh%z(mesh%nnr))
allocate(mesh%list(6*(mesh%nnr-2)))
allocate(mesh%lptr(6*(mesh%nnr-2)))
allocate(mesh%lend(mesh%nnr))

! Copy
mesh%order(1:mesh%nnr-1) = order
mesh%order_inv(1:mesh%n-1) = order_inv
mesh%x(1:mesh%nnr-1) = x
mesh%y(1:mesh%nnr-1) = y
mesh%z(1:mesh%nnr-1) = z
mesh%list(1:6*(mesh%nnr-3)) = list
mesh%lptr(1:6*(mesh%nnr-3)) = lptr
mesh%lend(1:mesh%nnr-1) = lend

! Compute new element coordinates
call trans(1,lat,lon,mesh%x(mesh%nnr:mesh%nnr),mesh%y(mesh%nnr:mesh%nnr),mesh%z(mesh%nnr:mesh%nnr))

! Update mesh
call addnod(1,mesh%nnr,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,mesh%lnew,info)

! Update order
mesh%order(mesh%nnr) = mesh%n
mesh%order_inv(mesh%n) = mesh%nnr

end subroutine addnode

!----------------------------------------------------------------------
! Subroutine: polygon
!> Purpose: compute polygon area
!----------------------------------------------------------------------
subroutine polygon(mesh,np,plist,area_polygon)

implicit none

! Passed variables
type(meshtype),intent(in) :: mesh               !< Mesh
integer,intent(in) :: np                        !< Number of points
integer,intent(in) :: plist(np)                 !< List of indices
real(kind_real),intent(out) :: area_polygon(np) !< Area

! Local variables
integer :: ip,i_src,index_triangle,i_src_last,i_src_new,i_src_stop,vertex_last,vertex_new,nb,info
integer :: lnew_tmp,lptr_tmp(6*(mesh%nnr-2)),lbtri(6,mesh%nnr),listc(6*(mesh%nnr-2))
real(kind_real) :: area_triangle,v1(3),v2(3),v3(3)
real(kind_real) :: xc(6*(mesh%nnr-2)),yc(6*(mesh%nnr-2)),zc(6*(mesh%nnr-2)),rc(6*(mesh%nnr-2))
logical :: loop

! Copy
lptr_tmp = mesh%lptr
lnew_tmp = mesh%lnew

! Compute Voronoi polygons
call crlist(mesh%nnr,mesh%nnr,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lend,lptr_tmp,lnew_tmp,lbtri,listc,nb,xc,yc,zc,rc,info)

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
