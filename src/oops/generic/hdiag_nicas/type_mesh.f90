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
use tools_stripack, only: trans,trmesh
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast
use type_randgen, only: rand_integer

implicit none

! Linear operator derived type
type meshtype
   integer,allocatable :: redundant(:) !< Redundant points
   integer :: nnr                      !< Number of non-redundant points
   real(kind_real),allocatable :: x(:) !< x-coordinate
   real(kind_real),allocatable :: y(:) !< y-coordinate
   real(kind_real),allocatable :: z(:) !< z-coordinate
   integer,allocatable :: list(:)      !< Stripack list
   integer,allocatable :: lptr(:)      !< Stripack list pointer
   integer,allocatable :: lend(:)      !< Stripack list end
   integer,allocatable :: order(:)     !< Order of shuffled points
   integer,allocatable :: order_inv(:) !< Inverse order of shuffled points
end type meshtype

private
public :: meshtype
public :: create_mesh,mesh_dealloc,check_mesh

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
integer :: i,j,k,lnew,info
integer :: n_loc(mpl%nproc),i_loc,iproc,progint
integer,allocatable :: jtab(:),near(:),next(:),i_glb(:,:),rbuf(:),sbuf(:)
real(kind_real),allocatable :: dist(:)
logical,allocatable :: done(:)

! Allocation
allocate(mesh%redundant(n))
call msi(mesh%redundant)

! Look for redundant points
if (lred) then
   ! MPI splitting
   iproc = 1
   n_loc = 0
   allocate(i_glb(n/mpl%nproc+1,mpl%nproc))
   do i=1,n
      n_loc(iproc) = n_loc(iproc)+1
      i_glb(n_loc(iproc),iproc) = i
      iproc = iproc+1
      if (iproc>mpl%nproc) iproc = 1
   end do

   ! Allocation
   allocate(done(n_loc(mpl%myproc)))

   ! Loop over points
   write(mpl%unit,'(a7,a)',advance='no') '','Look for redundant points: '
   call prog_init(progint,done)
   !$omp parallel do schedule(static) private(i_loc,i,j)
   do i_loc=1,n_loc(mpl%myproc)
      i = i_glb(i_loc,mpl%myproc)
      do j=1,i-1
         if (.not.(abs(lon(i)-lon(j))>0.0).and..not.(abs(lat(i)-lat(j))>0.0)) then
            mesh%redundant(i) = j
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
         if (iproc/=mpl%ioproc) then
            ! Allocation
            allocate(rbuf(n_loc(iproc)))

            ! Receive data on ioproc
            call mpl_recv(n_loc(iproc),rbuf,iproc,mpl%tag)

            ! Format data
            mesh%redundant(i_glb(1:n_loc(iproc),iproc)) = rbuf

            ! Release memory
            deallocate(rbuf)
         end if
      end do
   else
      ! Allocation
      allocate(sbuf(n_loc(mpl%myproc)))

      ! Format data
      sbuf = mesh%redundant(i_glb(1:n_loc(mpl%myproc),mpl%myproc))

      ! Send data to ioproc
      call mpl_send(n_loc(mpl%myproc),sbuf,mpl%ioproc,mpl%tag)

      ! Release memory
      deallocate(sbuf)
   end if
   mpl%tag = mpl%tag+1

   ! Broadcast data
   call mpl_bcast(mesh%redundant,mpl%ioproc)

   ! Check for successive redundant points
   do i=1,n
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
allocate(mesh%order_inv(n))
allocate(mesh%list(6*(mesh%nnr-2)))
allocate(mesh%lptr(6*(mesh%nnr-2)))
allocate(mesh%lend(mesh%nnr))
allocate(near(mesh%nnr))
allocate(next(mesh%nnr))
allocate(jtab(mesh%nnr))
allocate(mesh%x(mesh%nnr))
allocate(mesh%y(mesh%nnr))
allocate(mesh%z(mesh%nnr))
allocate(dist(mesh%nnr))

! Shuffle arrays (more efficient to compute the Delaunay triangulation)
i = 0
do j=1,n
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
call trmesh(mesh%nnr,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend,lnew,near,next,dist,info)

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
deallocate(mesh%redundant)
deallocate(mesh%x)
deallocate(mesh%y)
deallocate(mesh%z)
deallocate(mesh%list)
deallocate(mesh%lptr)
deallocate(mesh%lend)
deallocate(mesh%order)
deallocate(mesh%order_inv)

end subroutine mesh_dealloc

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

end module type_mesh
