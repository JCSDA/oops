!----------------------------------------------------------------------
! Module: type_kdtree
! Purpose: KD-tree derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_kdtree

use tools_const, only: pi
use tools_kdtree2, only: kdtree2,kdtree2_result,kdtree2_create,kdtree2_destroy, &
                       & kdtree2_n_nearest,kdtree2_r_count
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use tools_repro, only: rth
use tools_stripack, only: trans,scoord
use type_mpl, only: mpl_type

implicit none

! KD-tree derived type
type kdtree_type
    integer :: n                       ! Data size
    integer :: neff                    ! Effective KD-tree size
    logical,allocatable :: mask(:)     ! Mask
    type(kdtree2),pointer :: tp        ! KD-tree pointer
    integer,allocatable :: from_eff(:) ! Effective index conversion
contains
    procedure :: alloc => kdtree_alloc
    procedure :: init => kdtree_init
    procedure :: dealloc => kdtree_dealloc
    procedure :: find_nearest_neighbors => kdtree_find_nearest_neighbors
    procedure :: count_nearest_neighbors => kdtree_count_nearest_neighbors
end type kdtree_type

private
public :: kdtree_type

contains

!----------------------------------------------------------------------
! Subroutine: kdtree_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine kdtree_alloc(kdtree,mpl,n,mask)

implicit none

! Passed variables
class(kdtree_type),intent(inout) :: kdtree ! KD-tree
type(mpl_type),intent(inout) :: mpl        ! MPI data
integer,intent(in) :: n                    ! Number of points
logical,intent(in),optional :: mask(n)     ! Mask

! Local variables
character(len=1024),parameter :: subr = 'kdtree_alloc'

! Allocation
kdtree%n = n
allocate(kdtree%mask(n))

! Mask
if (present(mask)) then
   kdtree%mask = mask
else
   kdtree%mask = .true.
end if

! Effective tree size
kdtree%neff = count(kdtree%mask)

! Check size
if (kdtree%neff<1) call mpl%abort(subr,'mask should have at least one valid point to create a kdtree')

! Allocation
allocate(kdtree%from_eff(kdtree%neff))

end subroutine kdtree_alloc

!----------------------------------------------------------------------
! Subroutine: kdtree_init
! Purpose: initialization
!----------------------------------------------------------------------
subroutine kdtree_init(kdtree,mpl,lon,lat,sort,rearrange)

implicit none

! Passed variables
class(kdtree_type),intent(inout) :: kdtree  ! KD-tree
type(mpl_type),intent(inout) :: mpl         ! MPI data
real(kind_real),intent(in) :: lon(kdtree%n) ! Points longitudes
real(kind_real),intent(in) :: lat(kdtree%n) ! Points latitudes
logical,intent(in),optional :: sort         ! Sorting flag
logical,intent(in),optional :: rearrange    ! Rearranging flag

! Local variable
integer :: i,ieff
real(kind_real) :: input_data(3,kdtree%neff)
logical :: lsort,lrearrange

! Sorting flag
if (present(sort)) then
   lsort = sort
else
   lsort = .true.
end if

! Rearranging flag
if (present(rearrange)) then
   lrearrange = rearrange
else
   lrearrange = .true.
end if

! Loop over points
ieff = 0
do i=1,kdtree%n
   if (kdtree%mask(i)) then
      ieff = ieff+1

      ! Transform to cartesian coordinates
      call trans(mpl,1,lat(i),lon(i),input_data(1,ieff),input_data(2,ieff),input_data(3,ieff))

      ! Conversion
      kdtree%from_eff(ieff) = i
   end if
end do

! Create KD-tree
kdtree%tp => kdtree2_create(input_data,sort=lsort,rearrange=lrearrange)

end subroutine kdtree_init

!----------------------------------------------------------------------
! Subroutine: kdtree_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine kdtree_dealloc(kdtree)

implicit none

! Passed variables
class(kdtree_type),intent(inout) :: kdtree ! KD-tree

if (allocated(kdtree%from_eff)) then
   ! Release memory
   deallocate(kdtree%mask)
   deallocate(kdtree%from_eff)

   ! Delete KD-tree
   call kdtree2_destroy(kdtree%tp)
end if

end subroutine kdtree_dealloc

!----------------------------------------------------------------------
! Subroutine: kdtree_find_nearest_neighbors
! Purpose: find nearest neighbors using a KD-tree
!----------------------------------------------------------------------
subroutine kdtree_find_nearest_neighbors(kdtree,mpl,lon,lat,nn,nn_index,nn_dist)

implicit none

! Passed variables
class(kdtree_type),intent(in) :: kdtree    ! KD-tree
type(mpl_type),intent(inout) :: mpl        ! MPI data
real(kind_real),intent(in) :: lon          ! Point longitude
real(kind_real),intent(in) :: lat          ! Point latitude
integer,intent(in) :: nn                   ! Number of nearest neighbors to find
integer,intent(out) :: nn_index(nn)        ! Neareast neighbors index
real(kind_real),intent(out) :: nn_dist(nn) ! Neareast neighbors distance

! Local variables
integer :: i,j,nid
integer,allocatable :: order(:)
real(kind_real) :: lontmp(1),lattmp(1),qv(3)
type(kdtree2_result) :: results(nn)

! Copy lon/lat
lontmp(1) = lon
lattmp(1) = lat

! Transform to cartesian coordinates
call trans(mpl,1,lattmp,lontmp,qv(1),qv(2),qv(3))

! Find nearest neighbors
call kdtree2_n_nearest(kdtree%tp,qv,nn,results)
do i=1,nn
   nn_index(i) = kdtree%from_eff(results(i)%idx)
   nn_dist(i) = results(i)%sdis
end do

! Indistinguishability threshold for cross-plateform reproducibility
i = 1
do while (i<nn)
   ! Count indistinguishable neighbors
   nid = 1
   do j=i+1,nn
      if (abs(nn_dist(i)-nn_dist(j))<rth*nn_dist(i)) nid = nid+1
   end do

   ! Reorder
   if (nid>1) then
      allocate(order(nid))
      call qsort(nid,nn_index(i:i+nid-1),order)
      do j=1,nid
         nn_dist(i+j-1) = nn_dist(i+order(j)-1)
      end do
      deallocate(order)
   end if

   ! Update
   i = i+nid
end do

end subroutine kdtree_find_nearest_neighbors

!----------------------------------------------------------------------
! Subroutine: kdtree_count_nearest_neighbors
! Purpose: count nearest neighbors using a KD-tree
!----------------------------------------------------------------------
subroutine kdtree_count_nearest_neighbors(kdtree,mpl,lon,lat,sr,nn)

implicit none

! Passed variables
class(kdtree_type),intent(in) :: kdtree ! KD-tree
type(mpl_type),intent(inout) :: mpl     ! MPI data
real(kind_real),intent(in) :: lon       ! Point longitude
real(kind_real),intent(in) :: lat       ! Point latitude
real(kind_real),intent(in) :: sr        ! Spherical radius
integer,intent(out) :: nn               ! Number of nearest neighbors found

! Local variables
real(kind_real) :: lontmp(1),lattmp(1)
real(kind_real) :: qv(3),chordsq

! Copy lon/lat
lontmp(1) = lon
lattmp(1) = lat

! Transform to cartesian coordinates
call trans(mpl,1,lattmp,lontmp,qv(1),qv(2),qv(3))

! Convert radius on sphere to chord squared
chordsq = 4.0*sin(0.5*min(sr,pi))**2

! Count nearest neighbors
nn = kdtree2_r_count(kdtree%tp,qv,chordsq)

end subroutine kdtree_count_nearest_neighbors

end module type_kdtree
