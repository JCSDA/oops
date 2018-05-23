!----------------------------------------------------------------------
! Module: type_rng
!> Purpose: random numbers generator derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_rng

use iso_fortran_env, only: int64
use tools_display, only: msgerror
use tools_func, only: sphere_dist
use tools_kinds, only: kind_real
use tools_missing, only: msi,isnotmsi
use type_kdtree, only: kdtree_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

type rng_type
   integer(kind=int64) :: seed
contains
   procedure :: init => rng_init
   procedure :: reseed => rng_reseed
   procedure :: lcg => rng_lcg
   procedure :: rand_integer_0d
   procedure :: rand_integer_1d
   generic :: rand_integer => rand_integer_0d,rand_integer_1d
   procedure :: rand_real_0d
   procedure :: rand_real_1d
   procedure :: rand_real_2d
   procedure :: rand_real_3d
   procedure :: rand_real_4d
   procedure :: rand_real_5d
   generic :: rand_real => rand_real_0d,rand_real_1d,rand_real_2d,rand_real_3d,rand_real_4d,rand_real_5d
   procedure :: rand_gau_1d
   procedure :: rand_gau_5d
   generic :: rand_gau => rand_gau_1d,rand_gau_5d
   procedure :: initialize_sampling
end type rng_type

type(rng_type) :: rng !< Random number generator

integer,parameter :: default_seed = 140587            !< Default seed
integer(kind=int64),parameter :: a = 1103515245_int64 !< Linear congruential multiplier
integer(kind=int64),parameter :: c = 12345_int64      !< Linear congruential offset
integer(kind=int64),parameter :: m = 2147483648_int64 !< Linear congruential modulo

private
public :: rng

contains

!----------------------------------------------------------------------
! Subroutine: rng_init
!> Purpose: initialize the random number generator
!----------------------------------------------------------------------
subroutine rng_init(rng,nam)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
type(nam_type),intent(in) :: nam     !< Namelist variables

! Local variable
integer :: seed

! Set seed
if (nam%default_seed) then
   ! Default seed
   seed = default_seed
else
   ! Time-based seed
   call system_clock(count=seed)
end if

! Different seed for each task
seed = seed+mpl%myproc

! Long integer
rng%seed = int(seed,kind=int64)

! Print result
if (nam%default_seed) then
   write(mpl%unit,'(a7,a)') '','Linear congruential generator initialized with a default seed'
else
   write(mpl%unit,'(a7,a)') '','Linear congruential generator initialized'
end if
call flush(mpl%unit)

end subroutine rng_init

!----------------------------------------------------------------------
! Subroutine: rng_reseed
!> Purpose: re-seed the random number generator
!----------------------------------------------------------------------
subroutine rng_reseed(rng)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng !< Random number generator

! Local variable
integer :: seed

! Default seed
seed = default_seed

! Different seed for each task
seed = seed+mpl%myproc

! Long integer
rng%seed = int(seed,kind=int64)

end subroutine rng_reseed

!----------------------------------------------------------------------
! Subroutine: rng_lcg
!> Purpose: linear congruential generator
!----------------------------------------------------------------------
subroutine rng_lcg(rng,x)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(out) :: x             !< Random number between 0 and 1

! Update seed
rng%seed = mod(a*rng%seed+c,m)

! Random number
x = real(rng%seed,kind_real)/real(m-1,kind_real)

end subroutine rng_lcg

!----------------------------------------------------------------------
! Subroutine: rand_integer_0d
!> Purpose: generate a random integer, 0d
!----------------------------------------------------------------------
subroutine rand_integer_0d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
integer,intent(in) :: binf           !< Lower bound
integer,intent(in) :: bsup           !< Upper bound
integer,intent(out) :: ir            !< Random integer

! Local variables
real(kind_real) :: x

! Generate random number between 0 and 1
call rng%lcg(x)

! Adapt range
x = x*real(bsup-binf+1,kind_real)

! Add offset
ir = binf+int(x)

end subroutine rand_integer_0d

!----------------------------------------------------------------------
! Subroutine: rand_integer_1d
!> Purpose: generate a random integer, 1d
!----------------------------------------------------------------------
subroutine rand_integer_1d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
integer,intent(in) :: binf           !< Lower bound
integer,intent(in) :: bsup           !< Upper bound
integer,intent(out) :: ir(:)         !< Random integer

! Local variables
integer :: i

do i=1,size(ir)
   call rng%rand_integer(binf,bsup,ir(i))
end do

end subroutine rand_integer_1d

!----------------------------------------------------------------------
! Subroutine: rand_real_0d
!> Purpose: generate a random real, 0d
!----------------------------------------------------------------------
subroutine rand_real_0d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(in) :: binf   !< Lower bound
real(kind_real),intent(in) :: bsup   !< Upper bound
real(kind_real),intent(out) :: rr    !< Random integer

! Local variables
real(kind_real) :: x

! Generate random number between 0 and 1
call rng%lcg(x)

! Adapt range
x = x*(bsup-binf)

! Add offset
rr = binf+x

end subroutine rand_real_0d

!----------------------------------------------------------------------
! Subroutine: rand_real_1d
!> Purpose: generate a random real, 1d
!----------------------------------------------------------------------
subroutine rand_real_1d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(in) :: binf   !< Lower bound
real(kind_real),intent(in) :: bsup   !< Upper bound
real(kind_real),intent(out) :: rr(:) !< Random integer

! Local variables
integer :: i

do i=1,size(rr)
   call rng%rand_real(binf,bsup,rr(i))
end do

end subroutine rand_real_1d

!----------------------------------------------------------------------
! Subroutine: rand_real_2d
!> Purpose: generate a random real, 2d
!----------------------------------------------------------------------
subroutine rand_real_2d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng    !< Random number generator
real(kind_real),intent(in) :: binf      !< Lower bound
real(kind_real),intent(in) :: bsup      !< Upper bound
real(kind_real),intent(out) :: rr(:,:)  !< Random integer

! Local variables
integer :: i,j

do i=1,size(rr,2)
   do j=1,size(rr,1)
   call rng%rand_real(binf,bsup,rr(j,i))
   end do
end do

end subroutine rand_real_2d

!----------------------------------------------------------------------
! Subroutine: rand_real_3d
!> Purpose: generate a random real, 3d
!----------------------------------------------------------------------
subroutine rand_real_3d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng     !< Random number generator
real(kind_real),intent(in) :: binf       !< Lower bound
real(kind_real),intent(in) :: bsup       !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:) !< Random integer

! Local variables
integer :: i,j,k

do i=1,size(rr,3)
   do j=1,size(rr,2)
      do k=1,size(rr,1)
         call rng%rand_real(binf,bsup,rr(k,j,i))
      end do
   end do
end do

end subroutine rand_real_3d

!----------------------------------------------------------------------
! Subroutine: rand_real_4d
!> Purpose: generate a random real, 4d
!----------------------------------------------------------------------
subroutine rand_real_4d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng       !< Random number generator
real(kind_real),intent(in) :: binf         !< Lower bound
real(kind_real),intent(in) :: bsup         !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:) !< Random integer

! Local variables
integer :: i,j,k,l

do i=1,size(rr,4)
   do j=1,size(rr,3)
      do k=1,size(rr,2)
         do l=1,size(rr,1)
            call rng%rand_real(binf,bsup,rr(l,k,j,i))
         end do
      end do
   end do
end do

end subroutine rand_real_4d

!----------------------------------------------------------------------
! Subroutine: rand_real_5d
!> Purpose: generate a random real, 5d
!----------------------------------------------------------------------
subroutine rand_real_5d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         !< Random number generator
real(kind_real),intent(in) :: binf           !< Lower bound
real(kind_real),intent(in) :: bsup           !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:,:) !< Random integer

! Local variables
integer :: i,j,k,l,m

do i=1,size(rr,5)
   do j=1,size(rr,4)
      do k=1,size(rr,3)
         do l=1,size(rr,3)
            do m=1,size(rr,1)
               call rng%rand_real(binf,bsup,rr(m,l,k,j,i))
            end do
         end do
      end do
   end do
end do

end subroutine rand_real_5d

!----------------------------------------------------------------------
! Subroutine: rand_gau_1d
!> Purpose: generate random Gaussian deviates, 1d
!----------------------------------------------------------------------
subroutine rand_gau_1d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(out) :: rr(:) !< Random integer

! Local variables
integer :: i,iset
real(kind_real) :: gasdev,fac,gset,rsq,v1,v2

! Normal distribution
iset = 0
do i=1,size(rr,1)
   if( iset==0) then
      rsq = 0.0
      do while((rsq>=1.0).or.(rsq<=0.0))
         call rng%rand_real(0.0_kind_real,1.0_kind_real,v1)
         v1 = 2.0*v1-1.0
         call rng%rand_real(0.0_kind_real,1.0_kind_real,v2)
         v2 = 2.0*v2-1.0
         rsq = v1**2+v2**2
      end do
      fac = sqrt(-2.0*log(rsq)/rsq)
      gset = v1*fac
      gasdev = v2*fac
      iset = 1
   else
      gasdev = gset
      iset = 0
   end if
   rr(i) = gasdev
end do

end subroutine rand_gau_1d

!----------------------------------------------------------------------
! Subroutine: rand_gau_5d
!> Purpose: generate random Gaussian deviates, 5d
!----------------------------------------------------------------------
subroutine rand_gau_5d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         !< Random number generator
real(kind_real),intent(out) :: rr(:,:,:,:,:) !< Random integer

! Local variables
integer :: i,j,k,l

do i=1,size(rr,5)
   do j=1,size(rr,4)
      do k=1,size(rr,3)
         do l=1,size(rr,2)
            call rng%rand_gau(rr(:,l,k,j,i))
         end do
      end do
   end do
end do

end subroutine rand_gau_5d

!----------------------------------------------------------------------
! Subroutine: initialize_sampling
!> Purpose: intialize sampling
!----------------------------------------------------------------------
subroutine initialize_sampling(rng,n,lon,lat,mask,rh,ntry,nrep,ns,ihor)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
integer,intent(in) :: n              !< Number of points
real(kind_real),intent(in) :: lon(n) !< Longitudes
real(kind_real),intent(in) :: lat(n) !< Latitudes
logical,intent(in) :: mask(n)        !< Mask
real(kind_real),intent(in) :: rh(n)  !< Horizontal support radius
integer,intent(in) :: ntry           !< Number of tries
integer,intent(in) :: nrep           !< Number of replacements
integer,intent(in) :: ns             !< Number of samplings points
integer,intent(out) :: ihor(ns)      !< Horizontal sampling index

! Local variables
integer :: is,js,i,irep,irmax,itry,ir,nn_index(2),ismin(1)
real(kind_real) :: distmax,d,dist(ns),nn_dist(2)
logical :: lmask(n),smask(n)
type(kdtree_type) :: kdtree

if (ns>count(mask)) then
   call msgerror('ns greater that mask size in initialize_sampling')
elseif (ns==count(mask)) then
   is = 0
   do i=1,n
      if (mask(i)) then
         is = is+1
         ihor(is) = i
      end if
   end do
else
   ! Initialization
   call msi(ihor)
   lmask = mask
   smask = .false.
   is = 1
   irep = 1

   ! Define sampling
   do while (is<=ns)
      ! Create KD-tree (unsorted)
      if (is>2) call kdtree%create(n,lon,lat,mask=smask,sort=.false.)

      ! Initialization
      distmax = 0.0
      irmax = 0
      itry = 1

      ! Find a new point
      do while (itry<=ntry)
         ! Generate a random index
         call rng%rand_integer(1,n,ir)

         ! Check point validity
         if (lmask(ir)) then
            if (is==1) then
               ! Accept point
               irmax = ir
            else
               if (is==2) then
                  ! Compute distance
                  call sphere_dist(lon(ir),lat(ir),lon(ihor(1)),lat(ihor(1)),d)
               else
                  ! Find nearest neighbor distance
                  call kdtree%find_nearest_neighbors(lon(ir),lat(ir),1,nn_index(1:1),nn_dist(1:1))
                  d = nn_dist(1)**2/(rh(ir)**2+rh(nn_index(1))**2)
               end if

               ! Check distance
               if (d>distmax) then
                  distmax = d
                  irmax = ir
               end if
            end if
         end if

         ! Update itry
         itry = itry+1
      end do

      ! Delete kdtree
      if (is>2) call kdtree%delete

      ! Add point to sampling
      if (irmax>0) then
         ihor(is) = irmax
         lmask(irmax) = .false.
         smask(irmax) = .true.
         is = is+1
      end if

      if (is==ns+1) then
         ! Try replacement
         if (irep<=nrep) then
            ! Create KD-tree (unsorted)
            call kdtree%create(n,lon,lat,mask=smask,sort=.false.)

            ! Get minimum distance
            do js=1,ns
               ! Find nearest neighbor distance
               call kdtree%find_nearest_neighbors(lon(ihor(js)),lat(ihor(js)),2,nn_index,nn_dist)
               if (nn_index(1)==ihor(js)) then
                  dist(js) = nn_dist(2)
               elseif (nn_index(2)==ihor(js)) then
                  dist(js) = nn_dist(1)
               else
                  call msgerror('wrong index in replacement')
               end if
               dist(js) = dist(js)**2/(rh(nn_index(1))**2+rh(nn_index(2))**2)
            end do

            ! Delete kdtree
            call kdtree%delete

            ! Remove worst point
            ismin = minloc(dist)
            smask(ihor(ismin(1))) = .false.
            call msi(ihor(ismin(1)))

            ! Shift sampling
            if (ismin(1)<ns) ihor(ismin(1):ns-1) = ihor(ismin(1)+1:ns)

            ! Reset is to ns and try again!
            is = ns
         end if

         ! Update irep
         irep = irep+1
      end if
   end do
end if

end subroutine initialize_sampling

end module type_rng
