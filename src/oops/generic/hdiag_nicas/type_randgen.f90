!----------------------------------------------------------------------
! Module: type_randgen
!> Purpose: random numbers generator derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_randgen

use iso_c_binding, only: c_ptr,c_int,c_double
use tools_kinds, only: kind_real
use tools_missing, only: msi
use type_mpl, only: mpl,mpl_bcast
use type_nam, only: namtype

implicit none

type randgentype
   type(c_ptr) :: ptr !< Pointer to the C++ class
end type randgentype

type(randgentype) :: rng !< Random number generator

! C++ interface
interface
   function create_randgen_c(default_seed) bind(C,name="create_randgen")
   use iso_c_binding
   implicit none
   type(c_ptr) :: create_randgen_c
   integer(c_int),value :: default_seed
   end function create_randgen_c
end interface
interface
   subroutine delete_randgen_c(randgen) bind(C,name="delete_randgen")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: randgen
   end subroutine delete_randgen_c
end interface
interface
   subroutine rand_integer_c(randgen,binf,bsup,ir) bind(C,name="rand_integer")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: randgen
   integer(c_int),value :: binf
   integer(c_int),value :: bsup
   integer(c_int) :: ir
   end subroutine rand_integer_c
end interface
interface
   subroutine rand_real_c(randgen,binf,bsup,rr) bind(C,name="rand_real")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: randgen
   real(c_double),value :: binf
   real(c_double),value :: bsup
   real(c_double) :: rr
   end subroutine rand_real_c
end interface
interface
   subroutine initialize_sampling_c(randgen,n,lon,lat,mask,rh,ntry,nrep,ns,ihor) bind(C,name="initialize_sampling")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: randgen
   integer(c_int),value :: n
   real(c_double) :: lon(*)
   real(c_double) :: lat(*)
   integer(c_int) :: mask(*)
   real(c_double) :: rh(*)
   integer(c_int),value :: ntry
   integer(c_int),value :: nrep
   integer(c_int),value :: ns
   integer(c_int) :: ihor(*)
   end subroutine initialize_sampling_c
end interface

interface rand_integer
  module procedure rand_integer_0d
  module procedure rand_integer_1d
end interface

interface rand_real
  module procedure rand_real_0d
  module procedure rand_real_1d
  module procedure rand_real_2d
  module procedure rand_real_3d
  module procedure rand_real_4d
end interface

interface rand_gau
  module procedure rand_gau_1d
end interface

private
public :: create_randgen,delete_randgen,rand_integer,rand_real,rand_gau,initialize_sampling

contains

!----------------------------------------------------------------------
! Subroutine: create_randgen
!> Purpose: create a random number generator
!----------------------------------------------------------------------
subroutine create_randgen(nam)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables

! Local variable
integer :: default_seed

! Set default seed key to integer
if (nam%sam_default_seed) then
   default_seed = 1
else
   default_seed = 0
end if

! Call C++ function
rng%ptr = create_randgen_c(default_seed)

end subroutine create_randgen

!----------------------------------------------------------------------
! Subroutine: delete_randgen
!> Purpose: delete a random number generator
!----------------------------------------------------------------------
subroutine delete_randgen()

implicit none

! Call C++ function
call delete_randgen_c(rng%ptr)

end subroutine delete_randgen

!----------------------------------------------------------------------
! Subroutine: rand_integer_0d
!> Purpose: generate a random integer, 0d
!----------------------------------------------------------------------
subroutine rand_integer_0d(binf,bsup,ir)

implicit none

! Passed variables
integer,intent(in) :: binf !< Lower bound
integer,intent(in) :: bsup !< Upper bound
integer,intent(out) :: ir  !< Random integer

! Call C++ function
call rand_integer_c(rng%ptr,binf,bsup,ir)

end subroutine rand_integer_0d

!----------------------------------------------------------------------
! Subroutine: rand_integer_1d
!> Purpose: generate a random integer, 1d
!----------------------------------------------------------------------
subroutine rand_integer_1d(binf,bsup,ir)

implicit none

! Passed variables
integer,intent(in) :: binf   !< Lower bound
integer,intent(in) :: bsup   !< Upper bound
integer,intent(out) :: ir(:) !< Random integer

! Local variables
integer :: i

do i=1,size(ir)
   ! Call C++ function
   call rand_integer_c(rng%ptr,binf,bsup,ir(i))
end do

end subroutine rand_integer_1d

!----------------------------------------------------------------------
! Subroutine: rand_real_0d
!> Purpose: generate a random real, 0d
!----------------------------------------------------------------------
subroutine rand_real_0d(binf,bsup,rr)

implicit none

! Passed variables
real(kind_real),intent(in) :: binf !< Lower bound
real(kind_real),intent(in) :: bsup !< Upper bound
real(kind_real),intent(out) :: rr  !< Random integer

! Call C++ function
call rand_real_c(rng%ptr,binf,bsup,rr)

end subroutine rand_real_0d

!----------------------------------------------------------------------
! Subroutine: rand_real_1d
!> Purpose: generate a random real, 1d
!----------------------------------------------------------------------
subroutine rand_real_1d(binf,bsup,rr)

implicit none

! Passed variables
real(kind_real),intent(in) :: binf    !< Lower bound
real(kind_real),intent(in) :: bsup    !< Upper bound
real(kind_real),intent(out) :: rr(:)  !< Random integer

! Local variables
integer :: i

do i=1,size(rr)
   ! Call C++ function
   call rand_real_c(rng%ptr,binf,bsup,rr(i))
end do

end subroutine rand_real_1d

!----------------------------------------------------------------------
! Subroutine: rand_real_2d
!> Purpose: generate a random real, 2d
!----------------------------------------------------------------------
subroutine rand_real_2d(binf,bsup,rr)

implicit none

! Passed variables
real(kind_real),intent(in) :: binf     !< Lower bound
real(kind_real),intent(in) :: bsup     !< Upper bound
real(kind_real),intent(out) :: rr(:,:) !< Random integer

! Local variables
integer :: i,j

do i=1,size(rr,1)
   do j=1,size(rr,2)
      ! Call C++ function
      call rand_real_c(rng%ptr,binf,bsup,rr(i,j))
   end do
end do

end subroutine rand_real_2d

!----------------------------------------------------------------------
! Subroutine: rand_real_3d
!> Purpose: generate a random real, 3d
!----------------------------------------------------------------------
subroutine rand_real_3d(binf,bsup,rr)

implicit none

! Passed variables
real(kind_real),intent(in) :: binf       !< Lower bound
real(kind_real),intent(in) :: bsup       !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:) !< Random integer

! Local variables
integer :: i,j,k

do i=1,size(rr,1)
   do j=1,size(rr,2)
      do k=1,size(rr,3)
         ! Call C++ function
         call rand_real_c(rng%ptr,binf,bsup,rr(i,j,k))
      end do
   end do
end do

end subroutine rand_real_3d

!----------------------------------------------------------------------
! Subroutine: rand_real_4d
!> Purpose: generate a random real, 4d
!----------------------------------------------------------------------
subroutine rand_real_4d(binf,bsup,rr)

implicit none

! Passed variables
real(kind_real),intent(in) :: binf         !< Lower bound
real(kind_real),intent(in) :: bsup         !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:) !< Random integer

! Local variables
integer :: i,j,k,l

do i=1,size(rr,1)
   do j=1,size(rr,2)
      do k=1,size(rr,3)
         do l=1,size(rr,4)
            ! Call C++ function
            call rand_real_c(rng%ptr,binf,bsup,rr(i,j,k,l))
         end do
      end do
   end do
end do


end subroutine rand_real_4d

!----------------------------------------------------------------------
! Subroutine: rand_gau_1d
!> Purpose: generate random Gaussian deviates, 1d
!----------------------------------------------------------------------
subroutine rand_gau_1d(rr)

implicit none

! Passed variables
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
         call rand_real_c(rng%ptr,0.0_kind_real,1.0_kind_real,v1)
         v1 = 2.0*v1-1.0
         call rand_real_c(rng%ptr,0.0_kind_real,1.0_kind_real,v2)
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
! Subroutine: initialize_sampling
!> Purpose: intialize sampling
!----------------------------------------------------------------------
subroutine initialize_sampling(n,lon,lat,mask,rh,ntry,nrep,ns,ihor)

implicit none

! Passed variables
integer,intent(in) :: n              !< Number of points
real(kind_real),intent(in) :: lon(n) !< Longitudes
real(kind_real),intent(in) :: lat(n) !< Latitudes
integer,intent(in) :: mask(n)        !< Mask
real(kind_real),intent(in) :: rh(n)  !< Horizontal support radius
integer,intent(in) :: ntry           !< Number of tries
integer,intent(in) :: nrep           !< Number of replacements
integer,intent(in) :: ns             !< Number of samplings points
integer,intent(out) :: ihor(ns)      !< Horizontal sampling index

! Local variables
integer :: i,is

! Call C++ function
if (ns>=sum(mask)) then
   call msi(ihor)
   is = 0
   do i=1,n
      is = is+mask(i)
      ihor(is) = i
   end do
else
   if (mpl%main) call initialize_sampling_c(rng%ptr,n,lon,lat,mask,rh,ntry,nrep,ns,ihor)
   call mpl_bcast(ihor,mpl%ioproc)
end if

end subroutine initialize_sampling

end module type_randgen
