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
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr
use type_mpl, only: mpl
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
public :: create_mesh

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
integer,allocatable :: near(:),next(:)
real(kind_real),allocatable :: dist(:)

! Look for redundant points TODO : change that, make it parallel?
allocate(mesh%redundant(n))
call msi(mesh%redundant)
if (lred) then
   do i=1,n
      if (.not.isnotmsi(mesh%redundant(i))) then
         do j=i+1,n
            if ((abs(lon(i)-lon(j))<tiny(1.0)).and.(abs(lat(i)-lat(j))<tiny(1.0))) mesh%redundant(j) = i
         end do
      end if
   end do
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
   call rand_integer(randgen,1,mesh%nnr,j)
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

end module type_mesh
