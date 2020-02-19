!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Fortran interface to oops random number generator

module random_mod

use, intrinsic :: iso_c_binding  
use kinds

implicit none

private
public uniform_distribution
public normal_distribution

#include "random_interface.f"

!-------------------------------------------------------------------------------

!> interface for generating uniformly-distributed random numbers
!!
!! \param[inout] x Array that will contain the result
!! \param[in] minv minimum value of range
!! \param[in] maxv minimum value of range
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] sort if true, the result will be sorted from low to high
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
interface uniform_distribution
   module procedure uniform_float_distribution
   module procedure uniform_double_distribution
   module procedure uniform_int_distribution
   module procedure uniform_long_distribution
end interface uniform_distribution

!> interface for generating normally-distributed random numbers
!!
!! \param[inout] x Array that will contain the result
!! \param[in] mean mean
!! \param[in] sdev standard deviation
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
interface normal_distribution
   module procedure normal_float_distribution
   module procedure normal_float_distribution_2D
   module procedure normal_float_distribution_3D
   module procedure normal_double_distribution
   module procedure normal_double_distribution_2D
   module procedure normal_double_distribution_3D
end interface normal_distribution

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> class for generating uniformly-distributed random numbers (float)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] minv minimum value of range
!! \param[in] maxv minimum value of range
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] sort if true, the result will be sorted from low to high
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine uniform_float_distribution(x, minv, maxv, seed, sort, reset)
  implicit none
  real(c_float), intent(inout)  :: x(:)
  real(c_float), intent(in)     :: minv
  real(c_float), intent(in)     :: maxv
  integer, intent(in), optional :: seed
  logical, intent(in), optional :: sort
  logical, intent(in), optional :: reset
  integer(c_size_t) :: length, csort, creset
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif

  if (present(sort)) then
     if (sort) then
        csort = 1
     else
        csort = 0
     endif
  else
     csort = 0
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_uniform_float(length,minv,maxv,cseed,csort,creset,x)
  
end subroutine uniform_float_distribution

!> class for generating uniformly-distributed random numbers (double)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] minv minimum value of range
!! \param[in] maxv minimum value of range
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] sort if true, the result will be sorted from low to high
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine uniform_double_distribution(x, minv, maxv, seed, sort, reset)
  implicit none
  real(c_double), intent(inout)  :: x(:)
  real(c_double), intent(in)     :: minv
  real(c_double), intent(in)     :: maxv
  integer, intent(in), optional  :: seed
  logical, intent(in), optional  :: sort
  logical, intent(in), optional  :: reset
  integer(c_size_t) :: length, csort, creset
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif

  if (present(sort)) then
     if (sort) then
        csort = 1
     else
        csort = 0
     endif
  else
     csort = 0
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_uniform_double(length,minv,maxv,cseed,csort,creset,x)
  
end subroutine uniform_double_distribution
  
!> class for generating uniformly-distributed random numbers (integer)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] minv minimum value of range
!! \param[in] maxv minimum value of range
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] sort if true, the result will be sorted from low to high
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine uniform_int_distribution(x, minv, maxv, seed, sort, reset)
  implicit none
  integer(c_int32_t), intent(inout) :: x(:)
  integer(c_int32_t), intent(in)    :: minv
  integer(c_int32_t), intent(in)    :: maxv
  integer, intent(in), optional     :: seed
  logical, intent(in), optional     :: sort
  logical, intent(in), optional     :: reset
  integer(c_size_t) :: length, csort, creset
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif

  if (present(sort)) then
     if (sort) then
        csort = 1
     else
        csort = 0
     endif
  else
     csort = 0
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_uniform_int(length,minv,maxv,cseed,csort,creset,x)
  
end subroutine uniform_int_distribution
  
!> class for generating uniformly-distributed random numbers (long integer)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] minv minimum value of range
!! \param[in] maxv minimum value of range
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] sort if true, the result will be sorted from low to high
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine uniform_long_distribution(x, minv, maxv, seed, sort, reset)
  implicit none
  integer(c_int64_t), intent(inout) :: x(:)
  integer(c_int64_t), intent(in)    :: minv
  integer(c_int64_t), intent(in)    :: maxv
  integer, intent(in), optional     :: seed
  logical, intent(in), optional     :: sort
  logical, intent(in), optional     :: reset
  integer(c_size_t) :: length, csort, creset
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif

  if (present(sort)) then
     if (sort) then
        csort = 1
     else
        csort = 0
     endif
  else
     csort = 0
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_uniform_long(length,minv,maxv,cseed,csort,creset,x)
  
end subroutine uniform_long_distribution
  
!> class for generating normally-distributed random numbers (float)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] mean mean
!! \param[in] sdev standard deviation
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine normal_float_distribution(x, mean, sdev, seed, reset)
  implicit none
  real(c_float), intent(inout)  :: x(:)
  real(c_float), intent(in)     :: mean
  real(c_float), intent(in)     :: sdev
  integer, intent(in), optional :: seed
  logical, intent(in), optional :: reset
  integer(c_size_t) :: length, creset
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_normal_float(length,mean,sdev,cseed,creset,x)
  
end subroutine normal_float_distribution

!> class for generating normally-distributed random numbers (double)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] mean mean
!! \param[in] sdev standard deviation
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine normal_double_distribution(x, mean, sdev, seed, reset)
  implicit none
  real(c_double), intent(inout) :: x(:)
  real(c_double), intent(in)    :: mean
  real(c_double), intent(in)    :: sdev
  integer, intent(in), optional :: seed
  logical, intent(in), optional :: reset
  integer(c_size_t) :: length, creset
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_normal_double(length,mean,sdev,cseed,creset,x)
  
end subroutine normal_double_distribution

!> class for generating normally-distributed random numbers (double 2D)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] mean mean
!! \param[in] sdev standard deviation
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine normal_double_distribution_2D(x, mean, sdev, seed, reset)
  implicit none
  real(c_double), intent(inout) :: x(:,:)
  real(c_double), intent(in)    :: mean
  real(c_double), intent(in)    :: sdev
  integer, intent(in), optional :: seed
  logical, intent(in), optional :: reset
  integer(c_size_t) :: length, creset
  real(c_double), allocatable :: x_1d(:)
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_normal_double(length,mean,sdev,cseed,creset,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_double_distribution_2D

!> class for generating normally-distributed random numbers (double 3D)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] mean mean
!! \param[in] sdev standard deviation
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine normal_double_distribution_3D(x, mean, sdev, seed, reset)
  implicit none
  real(c_double), intent(inout)  :: x(:,:,:)
  real(c_double), intent(in)     :: mean
  real(c_double), intent(in)     :: sdev
  integer, intent(in), optional  :: seed
  logical, intent(in), optional  :: reset
  integer(c_size_t) :: length, creset
  real(c_double), allocatable :: x_1d(:)
  integer(c_int32_t) :: cseed
  real(kind_real) :: rseed

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_normal_double(length,mean,sdev,cseed,creset,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_double_distribution_3D

!> class for generating normally-distributed random numbers (float 2D)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] mean mean
!! \param[in] sdev standard deviation
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine normal_float_distribution_2D(x, mean, sdev, seed, reset)
  implicit none
  real(c_float), intent(inout)  :: x(:,:)
  real(c_float), intent(in)     :: mean
  real(c_float), intent(in)     :: sdev
  integer, intent(in), optional :: seed
  logical, intent(in), optional :: reset
  integer(c_size_t) :: length, creset
  integer(c_int32_t) :: cseed
  real(c_float), allocatable :: x_1d(:)
  real(kind_real) :: rseed

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_normal_float(length,mean,sdev,cseed,creset,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_float_distribution_2D

!> class for generating normally-distributed random numbers (float 3D)
!!
!! \param[inout] x Array that will contain the result
!! \param[in] mean mean
!! \param[in] sdev standard deviation
!! \param[in] seed The seed to use for the random number generator upon first
!!            call to the subroutine (optional).  If omitted, a seed
!!            will be generated based on the calendar time.  This parameter is 
!!            only used on the first instantiation of this class for a
!!            particular data type and/or if the **reset** flag is set to true.
!! \param[in] reset If this is set to **true** then this forces the generator
!!            to re-initialize itself with the specified seed.  Otherwise, the
!!            seed is only used on first instantiation (this is the default
!!            behavior).
!!  
subroutine normal_float_distribution_3D(x, mean, sdev, seed, reset)
  implicit none
  real(c_float), intent(inout)  :: x(:,:,:)
  real(c_float), intent(in)     :: mean
  real(c_float), intent(in)     :: sdev
  integer, intent(in), optional :: seed
  logical, intent(in), optional :: reset
  integer(c_size_t) :: length, creset
  integer(c_int32_t) :: cseed
  real(c_float), allocatable :: x_1d(:)
  real(kind_real) :: rseed

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     call random_number(rseed)
     cseed = int(1000*rseed,c_int32_t)
  endif
  if (present(reset)) then
     if (reset) then
        creset = 1
     else
        creset = 0
     endif
  else
     creset = 0
  endif
  call c_random_normal_float(length,mean,sdev,cseed,creset,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_float_distribution_3D

!-------------------------------------------------------------------------------
end module random_mod
