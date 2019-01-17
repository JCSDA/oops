!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Fortran interface to ObsSpace.

module random_mod

use, intrinsic :: iso_c_binding  
use kinds

implicit none

private
public uniform_distribution
public normal_distribution

#include "random_interface.f"

!-------------------------------------------------------------------------------

interface uniform_distribution
   module procedure uniform_float_distribution
   module procedure uniform_double_distribution
   module procedure uniform_int_distribution
   module procedure uniform_long_distribution
end interface uniform_distribution

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

subroutine uniform_float_distribution(x, minv, maxv, seed)
  implicit none
  real(c_float), intent(inout)  :: x(:)  !< Array that will contain the result
  real(c_float), intent(in)     :: minv  !< minimum value of range
  real(c_float), intent(in)     :: maxv  !< minimum value of range
  integer, intent(in), optional :: seed  !< random seed
  integer(c_size_t) :: length
  integer(c_int32_t) :: cseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_uniform_float(length,minv,maxv,cseed,x)
  
end subroutine uniform_float_distribution

subroutine uniform_double_distribution(x, minv, maxv, seed)
  implicit none
  real(c_double), intent(inout)  :: x(:)  !< Array that will contain the result
  real(c_double), intent(in)     :: minv  !< minimum value of range
  real(c_double), intent(in)     :: maxv  !< minimum value of range
  integer, intent(in), optional  :: seed  !< random seed
  integer(c_size_t) :: length
  integer(c_int32_t) :: cseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_uniform_double(length,minv,maxv,cseed,x)
  
end subroutine uniform_double_distribution
  
subroutine uniform_int_distribution(x, minv, maxv, seed)
  implicit none
  integer(c_int32_t), intent(inout) :: x(:)  !< Array that will contain the result
  integer(c_int32_t), intent(in)    :: minv  !< minimum value of range
  integer(c_int32_t), intent(in)    :: maxv  !< minimum value of range
  integer, intent(in), optional     :: seed  !< random seed
  integer(c_size_t) :: length
  integer(c_int32_t) :: cseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_uniform_int(length,minv,maxv,cseed,x)
  
end subroutine uniform_int_distribution
  
subroutine uniform_long_distribution(x, minv, maxv, seed)
  implicit none
  integer(c_int64_t), intent(inout) :: x(:)  !< Array that will contain the result
  integer(c_int64_t), intent(in)    :: minv  !< minimum value of range
  integer(c_int64_t), intent(in)    :: maxv  !< minimum value of range
  integer, intent(in), optional     :: seed  !< random seed
  integer(c_size_t) :: length
  integer(c_int32_t) :: cseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_uniform_long(length,minv,maxv,cseed,x)
  
end subroutine uniform_long_distribution
  
subroutine normal_float_distribution(x, mean, sdev, seed)
  implicit none
  real(c_float), intent(inout)  :: x(:)  !< Array that will contain the result
  real(c_float), intent(in)     :: mean  !< mean
  real(c_float), intent(in)     :: sdev  !< standard deviation
  integer, intent(in), optional :: seed  !< random seed
  integer(c_size_t) :: length  
  integer(c_int32_t) :: cseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_normal_float(length,mean,sdev,cseed,x)
  
end subroutine normal_float_distribution

subroutine normal_double_distribution(x, mean, sdev, seed)
  implicit none
  real(c_double), intent(inout) :: x(:)  !< Array that will contain the result
  real(c_double), intent(in)    :: mean  !< mean
  real(c_double), intent(in)    :: sdev  !< standard deviation
  integer, intent(in), optional :: seed  !< random seed
  integer(c_size_t) :: length  
  integer(c_int32_t) :: cseed
  
  length = size(x)
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_normal_double(length,mean,sdev,cseed,x)
  
end subroutine normal_double_distribution

subroutine normal_double_distribution_2D(x, mean, sdev, seed)
  implicit none
  real(c_double), intent(inout) :: x(:,:)  !< Array that will contain the result
  real(c_double), intent(in)    :: mean    !< mean
  real(c_double), intent(in)    :: sdev    !< standard deviation
  integer, intent(in), optional :: seed    !< random seed
  integer(c_size_t) :: length
  real(c_double), allocatable :: x_1d(:)
  integer(c_int32_t) :: cseed

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_normal_double(length,mean,sdev,cseed,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_double_distribution_2D

subroutine normal_double_distribution_3D(x, mean, sdev, seed)
  implicit none
  real(c_double), intent(inout)  :: x(:,:,:)  !< Array that will contain the result
  real(c_double), intent(in)     :: mean      !< mean
  real(c_double), intent(in)     :: sdev      !< standard deviation
  integer, intent(in), optional  :: seed      !< random seed
  integer(c_size_t) :: length
  real(c_double), allocatable :: x_1d(:)
  integer(c_int32_t) :: cseed

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_normal_double(length,mean,sdev,cseed,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_double_distribution_3D

subroutine normal_float_distribution_2D(x, mean, sdev, seed)
  implicit none
  real(c_float), intent(inout)  :: x(:,:)  !< Array that will contain the result
  real(c_float), intent(in)     :: mean    !< mean
  real(c_float), intent(in)     :: sdev    !< standard deviation
  integer, intent(in), optional :: seed    !< random seed
  integer(c_size_t) :: length
  integer(c_int32_t) :: cseed
  real(c_float), allocatable :: x_1d(:)

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_normal_float(length,mean,sdev,cseed,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_float_distribution_2D

subroutine normal_float_distribution_3D(x, mean, sdev, seed)
  implicit none
  real(c_float), intent(inout)  :: x(:,:,:)  !< Array that will contain the result
  real(c_float), intent(in)     :: mean      !< mean
  real(c_float), intent(in)     :: sdev      !< standard deviation
  integer, intent(in), optional :: seed      !< random seed
  integer(c_size_t) :: length
  integer(c_int32_t) :: cseed
  real(c_float), allocatable :: x_1d(:)

  length = size(x)
  allocate(x_1d(length))
  if (present(seed)) then
     cseed = seed
  else
     cseed = time()
  endif
  call c_random_normal_float(length,mean,sdev,cseed,x_1d)
  x = reshape(x_1d, shape(x))
  deallocate(x_1d)
  
end subroutine normal_float_distribution_3D

!-------------------------------------------------------------------------------
end module random_mod
