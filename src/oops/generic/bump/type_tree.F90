!----------------------------------------------------------------------
! Module: type_tree
! Purpose: tree derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_tree

use iso_c_binding, only: c_ptr
use fckit_kdtree_module, only: kdtree,kdtree_create,kdtree_destroy,kdtree_k_nearest_neighbors,kdtree_find_in_sphere
use tools_const, only: pi,rad2deg
use tools_func, only: lonlat2xyz,sphere_dist
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use tools_repro, only: rth,sup
use type_mpl, only: mpl_type

implicit none

! Tree derived type
type tree_type
    integer :: n                          ! Data size
    integer :: neff                       ! Effective tree size
    logical,allocatable :: mask(:)        ! Mask
    integer,allocatable :: from_eff(:)    ! Effective index conversion
    real(kind_real),allocatable :: lon(:) ! Longitudes
    real(kind_real),allocatable :: lat(:) ! Latitudes
    type(kdtree) :: kd                    ! KDTree from fckit
contains
    procedure :: alloc => tree_alloc
    procedure :: init => tree_init
    procedure :: dealloc => tree_dealloc
    procedure :: find_nearest_neighbors => tree_find_nearest_neighbors
    procedure :: count_nearest_neighbors => tree_count_nearest_neighbors
end type tree_type

real(kind_real),parameter :: nn_inc = 1.5 ! Increase factor for the nearest neighbors numbers search

private
public :: tree_type

contains

!----------------------------------------------------------------------
! Subroutine: tree_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine tree_alloc(tree,mpl,n,mask)

implicit none

! Passed variables
class(tree_type),intent(inout) :: tree ! Tree
type(mpl_type),intent(inout) :: mpl    ! MPI data
integer,intent(in) :: n                ! Number of points
logical,intent(in),optional :: mask(n) ! Mask

! Local variables
character(len=1024),parameter :: subr = 'tree_alloc'

! Allocation
tree%n = n
allocate(tree%mask(n))

! Mask
if (present(mask)) then
   tree%mask = mask
else
   tree%mask = .true.
end if

! Effective tree size
tree%neff = count(tree%mask)

! Check size
if (tree%neff<1) call mpl%abort(subr,'mask should have at least one valid point to create a tree')

! Allocation
allocate(tree%lon(tree%neff))
allocate(tree%lat(tree%neff))
allocate(tree%from_eff(tree%neff))

end subroutine tree_alloc

!----------------------------------------------------------------------
! Subroutine: tree_init
! Purpose: initialization
!----------------------------------------------------------------------
subroutine tree_init(tree,lon,lat)

implicit none

! Passed variables
class(tree_type),intent(inout) :: tree    ! Tree
real(kind_real),intent(in) :: lon(tree%n) ! Points longitudes (in radians)
real(kind_real),intent(in) :: lat(tree%n) ! Points latitudes (in radians)

! Local variable
integer :: i,ieff

! Loop over points
ieff = 0
do i=1,tree%n
   if (tree%mask(i)) then
      ieff = ieff+1

      ! Conversion
      tree%from_eff(ieff) = i

      ! Copy lon/lat
      tree%lon(ieff) = lon(i)
      tree%lat(ieff) = lat(i)
   end if
end do

! Create KDTree
tree%kd = kdtree_create(tree%neff,tree%lon*rad2deg,tree%lat*rad2deg)

end subroutine tree_init

!----------------------------------------------------------------------
! Subroutine: tree_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine tree_dealloc(tree)

implicit none

! Passed variables
class(tree_type),intent(inout) :: tree ! Tree

if (allocated(tree%from_eff)) then
   ! Release memory
   deallocate(tree%mask)
   deallocate(tree%from_eff)
   deallocate(tree%lon)
   deallocate(tree%lat)

   ! Destroy KDTree
   call kdtree_destroy(tree%kd)
end if

end subroutine tree_dealloc

!----------------------------------------------------------------------
! Subroutine: tree_find_nearest_neighbors
! Purpose: find nearest neighbors using a KDTree
!----------------------------------------------------------------------
subroutine tree_find_nearest_neighbors(tree,lon,lat,nn,nn_index,nn_dist)

implicit none

! Passed variables
class(tree_type),intent(in) :: tree                 ! Tree
real(kind_real),intent(in) :: lon                   ! Point longitude (in radians)
real(kind_real),intent(in) :: lat                   ! Point latitude (in radians)
integer,intent(in) :: nn                            ! Number of nearest neighbors to find
integer,intent(out) :: nn_index(nn)                 ! Neareast neighbors index
real(kind_real),intent(out),optional :: nn_dist(nn) ! Neareast neighbors distance

! Local variables
integer :: nn_tmp,i,j,nid
integer,allocatable :: order(:),nn_index_tmp(:)
real(kind_real) :: dist_ref,dist_last
real(kind_real),allocatable :: nn_dist_tmp(:)
logical :: separate

! Initialization
separate = .false.
nn_tmp = min(nn+1,tree%neff)

do while (.not.separate)
   ! Allocation
   allocate(nn_index_tmp(nn_tmp))

   ! Find neighbors
   call kdtree_k_nearest_neighbors(tree%kd,lon*rad2deg,lat*rad2deg,nn_tmp,nn_index_tmp)

   ! Check distance between reference and last nearest neighbors
   call sphere_dist(lon,lat,tree%lon(nn_index_tmp(nn)),tree%lat(nn_index_tmp(nn)),dist_ref)
   call sphere_dist(lon,lat,tree%lon(nn_index_tmp(nn_tmp)),tree%lat(nn_index_tmp(nn_tmp)),dist_last)
   if (sup(dist_last,dist_ref).or.(nn_tmp==tree%neff)) then
      ! Last neighbor is significantly further away (or all points are used)
      separate = .true.
   else
      ! Last neighbor is at the same distance as the reference neighbor, increase number of neighbors
      nn_tmp = min(int(nn_inc*real(nn_tmp,kind_real)),tree%neff)

      ! Release memory
      deallocate(nn_index_tmp)
   end if
end do

! Allocation
allocate(nn_dist_tmp(nn_tmp))

! Compute distance
do i=1,nn_tmp
   call sphere_dist(lon,lat,tree%lon(nn_index_tmp(i)),tree%lat(nn_index_tmp(i)),nn_dist_tmp(i))
end do

! Transform indices
nn_index_tmp = tree%from_eff(nn_index_tmp)

! Reorder neighbors based on their index
i = 1
do while (i<nn_tmp)
   ! Count indistinguishable neighbors
   nid = 1
   do j=i+1,nn_tmp
      if (abs(nn_dist_tmp(i)-nn_dist_tmp(j))<rth*nn_dist_tmp(i)) nid = nid+1
   end do

   ! Reorder
   if (nid>1) then
      allocate(order(nid))
      call qsort(nid,nn_index_tmp(i:i+nid-1),order)
      do j=1,nid
         nn_dist_tmp(i+j-1) = nn_dist_tmp(i+order(j)-1)
      end do
      deallocate(order)
   end if

   ! Update
   i = i+nid
end do

! Copy nn_index
nn_index = nn_index_tmp(1:nn)

! Copy nn_dist if required
if (present(nn_dist)) nn_dist = nn_dist_tmp(1:nn)

! Release memory
deallocate(nn_index_tmp)
deallocate(nn_dist_tmp)

end subroutine tree_find_nearest_neighbors

!----------------------------------------------------------------------
! Subroutine: tree_count_nearest_neighbors
! Purpose: count nearest neighbors using a tree
!----------------------------------------------------------------------
subroutine tree_count_nearest_neighbors(tree,lon,lat,sr,nn)

implicit none

! Passed variables
class(tree_type),intent(in) :: tree ! Tree
real(kind_real),intent(in) :: lon   ! Point longitude (in radians)
real(kind_real),intent(in) :: lat   ! Point latitude (in radians)
real(kind_real),intent(in) :: sr    ! Spherical radius (in radians)
integer,intent(out) :: nn           ! Number of nearest neighbors found

! Count nearest neighbors
call kdtree_find_in_sphere(tree%kd,lon*rad2deg,lat*rad2deg,sr,nn)

end subroutine tree_count_nearest_neighbors

end module type_tree
