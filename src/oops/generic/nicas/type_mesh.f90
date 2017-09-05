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
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr
use tools_stripack, only: trans,trmesh
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast
use type_randgen, only: randgentype,rand_integer

implicit none

! Linear operator derived type
type meshtype
   integer,allocatable :: redundant(:)
   integer :: nnr                 !< 
   real(kind_real),allocatable :: x(:)
   real(kind_real),allocatable :: y(:)
   real(kind_real),allocatable :: z(:)
   integer,allocatable :: list(:) 
   integer,allocatable :: lptr(:) 
   integer,allocatable :: lend(:) 
   integer,allocatable :: order(:)
end type meshtype

private
public :: meshtype
public :: create_mesh,mesh_dealloc

contains

!----------------------------------------------------------------------
! Subroutine: create_mesh
!> Purpose: create mesh
!----------------------------------------------------------------------
subroutine create_mesh(randgen,n,lon,lat,lred,mesh)

implicit none

! Passed variables
type(randgentype),intent(in) :: randgen
integer,intent(in) :: n
real(kind_real),intent(in) :: lon(n)
real(kind_real),intent(in) :: lat(n)
logical,intent(in) :: lred
type(meshtype),intent(inout) :: mesh

! Local variables
integer :: i,j,k,lnew,info
integer :: n_loc(mpl%nproc),i_loc,iproc,progint
integer,allocatable :: near(:),next(:),i_glb(:,:),rbuf(:),sbuf(:)
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
   !$omp parallel do private(i_loc,i,j)
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
end if
mesh%nnr = count(.not.isnotmsi(mesh%redundant))

! Allocation
allocate(mesh%order(mesh%nnr))
allocate(mesh%list(6*(mesh%nnr-2)))
allocate(mesh%lptr(6*(mesh%nnr-2)))
allocate(mesh%lend(mesh%nnr))
allocate(near(mesh%nnr))
allocate(next(mesh%nnr))
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
do i=mesh%nnr,2,-1
   call rand_integer(randgen,1,mesh%nnr,.true.,j)
   k = mesh%order(j)
   mesh%order(j) = mesh%order(i)
   mesh%order(i) = k
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
type(meshtype),intent(inout) :: mesh

! Release memory
deallocate(mesh%redundant)
deallocate(mesh%x)
deallocate(mesh%y)
deallocate(mesh%z)
deallocate(mesh%list)
deallocate(mesh%lptr)
deallocate(mesh%lend)
deallocate(mesh%order)

end subroutine mesh_dealloc

end module type_mesh
