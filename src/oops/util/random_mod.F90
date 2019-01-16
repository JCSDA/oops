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
!public normal_distribution

#include "random_interface.f"

!-------------------------------------------------------------------------------

interface uniform_distribution
   module procedure uniform_float_distribution
   module procedure uniform_double_distribution
!   module procedure uniform_int_distribution
!   module procedure uniform_long_distribution
end interface uniform_distribution

!interface normal_distribution
!   module procedure normal_float_distribution
!   module procedure normal_double_distribution
!end interface normal_distribution

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine uniform_float_distribution(x, minv, maxv, seed)
  implicit none
  real(c_float), intent(inout)     :: x(:)  !< Array that will contain the result
  real(c_float), intent(in)        :: minv  !< minimum value of range
  real(c_float), intent(in)        :: maxv  !< minimum value of range
  integer, intent(inout), optional :: seed  !< random seed
  integer(c_size_t) :: length  
  
  length = size(x)
  if (.not. present(seed)) seed = time()
  call c_random_uniform_float(length,minv,maxv,int(seed,c_int32_t),x)
  
end subroutine uniform_float_distribution

subroutine uniform_double_distribution(x, minv, maxv, seed)
  implicit none
  real(c_double), intent(inout)     :: x(:)  !< Array that will contain the result
  real(c_double), intent(in)        :: minv  !< minimum value of range
  real(c_double), intent(in)        :: maxv  !< minimum value of range
  integer, intent(inout), optional  :: seed  !< random seed
  integer(c_size_t) :: length
  
  length = size(x)
  if (.not. present(seed)) seed = time()
  call c_random_uniform_double(length,minv,maxv,int(seed,c_int32_t),x)
  
end subroutine uniform_double_distribution
  
!-------------------------------------------------------------------------------
end module random_mod
