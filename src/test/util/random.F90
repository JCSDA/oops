!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

!> Test interface for C++ random number generators called from Fortran

module test_random

use, intrinsic :: iso_c_binding
use config_mod
use kinds
use random_mod

implicit none
private

integer, parameter :: max_string = 800, max_number = 30

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Test uniform real distribution
!
  
integer(c_int32_t) function c_test_uniform_real(c_conf) bind(c,name='test_uniform_real_f')  
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: N, seed, i
character(len=max_number) :: range(2)
real(kind_real), allocatable :: x(:), x_check(:)
real(kind_real) :: minv, maxv, tol
character(len=*), parameter :: myname_="test_uniform_real"
character(max_string) :: err_msg

N = config_get_int(c_conf, "N")
seed = config_get_int(c_conf, "seed")

if (size(config_get_string_vector(c_conf, max_string, "uniform_real_range")) == 2) then
   range = config_get_string_vector(c_conf, max_number, "uniform_real_range")
   read(range(1),*) minv
   read(range(2),*) maxv
else
   write(err_msg,*) myname_ // "error reading range"
   call abor1_ftn(err_msg)
endif

!> Compute random vector
allocate(x(N))
call uniform_distribution(x, minv, maxv, seed)

if (size(config_get_string_vector(c_conf, max_string, "uniform_real_answer")) == N) then
   allocate(x_check(N))
   call config_get_double_vector(c_conf, "uniform_real_answer", x_check)
   do i=1,N 
      write(*,*) "MSM xcheck: ", x(i), x_check(i)
   enddo
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

! Test uniform real distribution.
! The tolerance is based on the precision of the data type but multiply by 10 to
! allow for roundoff error in the last digit during normalization

tol = 10.0_kind_real * epsilon(minv)
write(*,*) "MSM tol: ", tol

! if (passes)
c_test_uniform_real = 0
! else
!c_test_uniform_real = 1

! clean up
deallocate(x,x_check)

end function c_test_uniform_real

!-------------------------------------------------------------------------------

end module test_random
