!----------------------------------------------------------------------
! Module: type_kdtree
!> Purpose: KD-tree derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_kdtree

use tools_kdtree2, only: kdtree2,kdtree2_result,kdtree2_create,kdtree2_destroy, &
                       & kdtree2_n_nearest,kdtree2_r_count
use tools_kinds, only: kind_real
use tools_missing, only: msi,isnotmsi
use tools_qsort, only: qsort
use tools_repro, only: rth
use tools_stripack, only: trans,scoord
use type_mpl, only: mpl_type

implicit none

! KD-tree derived type
type kdtree_type
    type(kdtree2),pointer :: tp        !< KD-tree pointer
    integer,allocatable :: from_eff(:) !< Effective index conversion
contains
    procedure :: create => kdtree_create
    procedure :: dealloc => kdtree_dealloc
    procedure :: find_nearest_neighbors => kdtree_find_nearest_neighbors
    procedure :: count_nearest_neighbors => kdtree_count_nearest_neighbors
end type kdtree_type

private
public :: kdtree_type

contains

!----------------------------------------------------------------------
! Subroutine: kdtree_create
!> Purpose: create a KD-tree
!----------------------------------------------------------------------
subroutine kdtree_create(kdtree,mpl,n,lon,lat,mask,sort,rearrange)

implicit none

! Passed variables
class(kdtree_type),intent(inout) :: kdtree !< KD-tree
type(mpl_type),intent(in) :: mpl           !< MPI data
integer,intent(in) :: n                    !< Number of points
real(kind_real),intent(in) :: lon(n)       !< Points longitudes
real(kind_real),intent(in) :: lat(n)       !< Points latitudes
logical,intent(in),optional :: mask(n)     !< Mask
logical,intent(in),optional :: sort        !< Sorting flag
logical,intent(in),optional :: rearrange   !< Rearranging flag

! Local variable
integer :: neff,i,ieff
real(kind_real),allocatable :: input_data(:,:)
logical :: lmask(n),lsort,lrearrange

! Mask
if (present(mask)) then
   lmask = mask
else
   lmask = .true.
end if

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

! Effective tree size
neff = count(lmask)

! Check size
if (neff<1) call mpl%abort('mask should have at least one valid point to create a kdtree')

! Allocation
allocate(input_data(3,neff))
allocate(kdtree%from_eff(neff))

! Loop over points
ieff = 0
do i=1,n
   if (lmask(i)) then
      ieff = ieff+1

      ! Transform to cartesian coordinates
      call trans(1,lat(i),lon(i),input_data(1,ieff),input_data(2,ieff),input_data(3,ieff))

      ! Conversion
      kdtree%from_eff(ieff) = i
   end if
end do

! Create KD-tree
kdtree%tp => kdtree2_create(input_data,sort=lsort,rearrange=lrearrange)

end subroutine kdtree_create

!----------------------------------------------------------------------
! Subroutine: kdtree_dealloc
!> Purpose: deallocate KD-tree
!----------------------------------------------------------------------
subroutine kdtree_dealloc(kdtree)

implicit none

! Passed variables
class(kdtree_type),intent(inout) :: kdtree !< KD-tree

if (allocated(kdtree%from_eff)) then
   ! Deallocation
   deallocate(kdtree%from_eff)

   ! Delete KD-tree
   call kdtree2_destroy(kdtree%tp)
end if

end subroutine kdtree_dealloc

!----------------------------------------------------------------------
! Subroutine: kdtree_find_nearest_neighbors
!> Purpose: find nearest neighbors using a KD-tree
!----------------------------------------------------------------------
subroutine kdtree_find_nearest_neighbors(kdtree,lon,lat,nn,nn_index,nn_dist)

implicit none

! Passed variables
class(kdtree_type),intent(in) :: kdtree    !< KD-tree
real(kind_real),intent(in) :: lon          !< Point longitude
real(kind_real),intent(in) :: lat          !< Point latitude
integer,intent(in) :: nn                   !< Number of nearest neighbors to find
integer,intent(out) :: nn_index(nn)        !< Neareast neighbors index
real(kind_real),intent(out) :: nn_dist(nn) !< Neareast neighbors distance

! Local variables
integer :: i,j,nid
integer,allocatable :: order(:)
real(kind_real) :: lontmp(1),lattmp(1)
real(kind_real) :: qv(3)
type(kdtree2_result) :: results(nn)


! Copy lon/lat
lontmp(1) = lon
lattmp(1) = lat

! Transform to cartesian coordinates
call trans(1,lattmp,lontmp,qv(1),qv(2),qv(3))

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
!> Purpose: count nearest neighbors using a KD-tree
!----------------------------------------------------------------------
subroutine kdtree_count_nearest_neighbors(kdtree,lon,lat,sr,nn)

implicit none

! Passed variables
class(kdtree_type),intent(in) :: kdtree !< KD-tree
real(kind_real),intent(in) :: lon       !< Point longitude
real(kind_real),intent(in) :: lat       !< Point latitude
real(kind_real),intent(in) :: sr        !< Spherical radius
integer,intent(out) :: nn               !< Number of nearest neighbors found

! Local variables
real(kind_real) :: lontmp(1),lattmp(1)
real(kind_real) :: qv(3),brsq

! Copy lon/lat
lontmp(1) = lon
lattmp(1) = lat

! Transform to cartesian coordinates
call trans(1,lattmp,lontmp,qv(1),qv(2),qv(3))

! Convert spherical radius to squared ball radius
brsq = (1.0-cos(sr))**2+sin(sr)**2

! Count nearest neighbors
nn = kdtree2_r_count(kdtree%tp,qv,brsq)

end subroutine kdtree_count_nearest_neighbors

end module type_kdtree
