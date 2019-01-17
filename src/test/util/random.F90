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
use fckit_log_module, only : fckit_log

implicit none
private

integer, parameter :: max_string = 800

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Test uniform real distribution
!
  
integer(c_int32_t) function c_test_uniform_real(c_conf) bind(c,name='test_uniform_real_f')  
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, seed, i
real(kind_real) :: range(2)
real(kind_real), allocatable :: x(:), x_check(:)
real(kind_real) :: minv, maxv, tol, err
character(len=*), parameter :: myname_="test_uniform_real"
character(max_string) :: err_msg

n = config_get_int(c_conf, "N")
seed = config_get_int(c_conf, "seed")

if (size(config_get_string_vector(c_conf, max_string, "uniform_real_range")) == 2) then
   call config_get_double_vector(c_conf, "uniform_real_range", range, [0.0_kind_real,1.0_kind_real])
   minv = range(1)
   maxv = range(2)
else
   write(err_msg,*) myname_ // "error reading range"
   call abor1_ftn(err_msg)
endif

!> Compute random vector
allocate(x(n))
call uniform_distribution(x, minv, maxv, seed)

if (size(config_get_string_vector(c_conf, max_string, "uniform_real_answer")) == n) then
   allocate(x_check(n))
   call config_get_double_vector(c_conf, "uniform_real_answer", x_check)
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

! The tolerance is based on the precision of the data type but multiply by 10 to
! allow for roundoff error in the last digit during normalization
tol = 10.0_kind_real * epsilon(minv)

c_test_uniform_real = 0

do i = 1, n
   err = abs(x(i)-x_check(i))
   if (err > tol*abs(x_check(i))) then
       c_test_uniform_real = 1
       write(err_msg,*) myname_ // " error exceeds tolerence: ", err, tol
       call fckit_log%error(err_msg)
   endif
enddo

! clean up
deallocate(x,x_check)

end function c_test_uniform_real

!-------------------------------------------------------------------------------
!> Test uniform int distribution
!
  
integer(c_int32_t) function c_test_uniform_int(c_conf) bind(c,name='test_uniform_int_f')  
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, seed, i
integer(c_int32_t) :: range(2)
integer(c_int32_t), allocatable :: x(:), x_check(:)
integer(c_int32_t) :: minv, maxv
character(len=*), parameter :: myname_="test_uniform_int"
character(max_string) :: err_msg

n = config_get_int(c_conf, "N")
seed = config_get_int(c_conf, "seed")

if (size(config_get_string_vector(c_conf, max_string, "uniform_int_range")) == 2) then
   call config_get_int_vector(c_conf, "uniform_int_range", range, [1,100])
   minv = range(1)
   maxv = range(2)
else
   write(err_msg,*) myname_ // "error reading range"
   call abor1_ftn(err_msg)
endif

!> Compute random vector
allocate(x(n))
call uniform_distribution(x, minv, maxv, seed)

if (size(config_get_string_vector(c_conf, max_string, "uniform_int_answer")) == n) then
   allocate(x_check(n))
   call config_get_int_vector(c_conf, "uniform_int_answer", x_check)
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

c_test_uniform_int = 0

do i = 1, n
   if (x(i) /= x_check(i)) then
      c_test_uniform_int = 1
      write(err_msg,*) myname_ // " mismatch: ", x(i), x_check(i)
      call fckit_log%error(err_msg)
   endif
enddo

! clean up
deallocate(x,x_check)

end function c_test_uniform_int

!-------------------------------------------------------------------------------
!> Test normal real distribution
!
  
integer(c_int32_t) function c_test_normal_real(c_conf) bind(c,name='test_normal_real_f')  
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, seed, i
real(kind_real) :: mean, sdev, tol, err
real(kind_real), allocatable :: x(:), x_check(:)
character(len=*), parameter :: myname_="test_normal_real"
character(max_string) :: err_msg

n = config_get_int(c_conf, "N")
seed = config_get_int(c_conf, "seed")

mean = config_get_real(c_conf, "normal_mean",0.0_kind_real)
sdev = config_get_real(c_conf, "normal_sdev",1.0_kind_real)

!> Compute random vector
allocate(x(n))
call normal_distribution(x, mean, sdev, seed)

if (size(config_get_string_vector(c_conf, max_string, "normal_answer")) == n) then
   allocate(x_check(n))
   call config_get_double_vector(c_conf, "normal_answer", x_check)
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

! The tolerance is based on the precision of the data type but multiply by 10 to
! allow for roundoff error in the last digit during normalization
tol = 10.0_kind_real * epsilon(mean)

c_test_normal_real = 0

do i = 1, n
   err = abs(x(i)-x_check(i))
   if (err > tol*abs(x_check(i))) then
       c_test_normal_real = 1
       write(err_msg,*) myname_ // " error exceeds tolerence: ", err, tol
       call fckit_log%error(err_msg)
   endif
enddo

! clean up
deallocate(x,x_check)

end function c_test_normal_real

!-------------------------------------------------------------------------------

end module test_random
