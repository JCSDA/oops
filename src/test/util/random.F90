!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!> Test interface for C++ random number generators called from Fortran

module test_random

use fckit_configuration_module, only: fckit_configuration
use, intrinsic :: iso_c_binding
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
real(c_float) :: range(2)
real(c_float), allocatable :: x(:), x_check(:), real_array(:)
real(c_float) :: minv, maxv, tol, err
character(len=*), parameter :: myname_="test_uniform_real"
character(max_string) :: err_msg, ans
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("N",n)
call f_conf%get_or_die("seed",seed)

if (f_conf%get_size("uniform_real_range") == 2) then
   call f_conf%get_or_die("uniform_real_range",real_array)
   range = real_array
   minv = range(1)
   maxv = range(2)
else
   write(err_msg,*) myname_ // "error reading range"
   call abor1_ftn(err_msg)
endif

!> Compute random vector
allocate(x(n))
call uniform_distribution(x, minv, maxv, seed, reset=.true.)

! write result for inclusion in config file
err_msg = achar(10) // myname_ // " Testing oops::util::random.F90 Uniform Real Distribution: "
call fckit_log%info(err_msg)
do i = 1, n
   write(ans,'(F30.16)') x(i)
   call fckit_log%info("   - " // adjustl(ans))
enddo

if (f_conf%get_size("uniform_real_answer") == n) then
   allocate(x_check(n))
   call f_conf%get_or_die("uniform_real_answer",real_array)
   x_check = real_array
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
!> Test uniform double distribution
!

integer(c_int32_t) function c_test_uniform_double(c_conf) bind(c,name='test_uniform_double_f')
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, seed, i
real(c_double) :: range(2)
real(c_double), allocatable :: x(:), x_check(:), real_array(:)
real(c_double) :: minv, maxv, tol, err
character(len=*), parameter :: myname_="test_uniform_double"
character(max_string) :: err_msg, ans
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("N",n)
call f_conf%get_or_die("seed",seed)

if (f_conf%get_size("uniform_real_range") == 2) then
   call f_conf%get_or_die("uniform_real_range",real_array)
   range = real_array
   minv = range(1)
   maxv = range(2)
else
   write(err_msg,*) myname_ // "error reading range"
   call abor1_ftn(err_msg)
endif

!> Compute random vector
allocate(x(n))
call uniform_distribution(x, minv, maxv, seed, reset=.true.)

! write result for inclusion in config file
err_msg = achar(10) // myname_ // " Testing oops::util::random.F90 Uniform Double Distribution: "
call fckit_log%info(err_msg)
do i = 1, n
   write(ans,'(F30.16)') x(i)
   call fckit_log%info("   - " // adjustl(ans))
enddo

if (f_conf%get_size("uniform_real_answer") == n) then
   allocate(x_check(n))
   call f_conf%get_or_die("uniform_real_answer",real_array)
   x_check = real_array
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

! The tolerance is based on the precision of the data type but multiply by 10 to
! allow for roundoff error in the last digit during normalization
tol = 10.0_kind_real * epsilon(minv)

c_test_uniform_double = 0

do i = 1, n
   err = abs(x(i)-x_check(i))
   if (err > tol*abs(x_check(i))) then
       c_test_uniform_double = 1
       write(err_msg,*) myname_ // " error exceeds tolerence: ", err, tol
       call fckit_log%error(err_msg)
   endif
enddo

! clean up
deallocate(x,x_check)

end function c_test_uniform_double

!-------------------------------------------------------------------------------
!> Test uniform int distribution
!

integer(c_int32_t) function c_test_uniform_int(c_conf) bind(c,name='test_uniform_int_f')
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, seed, i
integer(c_int32_t) :: range(2)
integer(c_int32_t), allocatable :: x(:), x_check(:), integer_array(:)
integer(c_int32_t) :: minv, maxv
character(len=*), parameter :: myname_="test_uniform_int"
character(max_string) :: err_msg, ans
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("N",n)
call f_conf%get_or_die("seed",seed)

if (f_conf%get_size("uniform_int_range") == 2) then
   call f_conf%get_or_die("uniform_int_range",integer_array)
   range = integer_array
   minv = range(1)
   maxv = range(2)
else
   write(err_msg,*) myname_ // "error reading range"
   call abor1_ftn(err_msg)
endif

!> Compute random vector
allocate(x(n))
call uniform_distribution(x, minv, maxv, seed, reset=.true.)

! write result for inclusion in config file
err_msg = achar(10) // myname_ // " Testing oops::util::random.F90 Uniform Int Distribution: "
call fckit_log%info(err_msg)
do i = 1, n
   write(ans,'(I12)') x(i)
   call fckit_log%info("   - " // adjustl(ans))
enddo

if (f_conf%get_size("uniform_int_answer") == n) then
   allocate(x_check(n))
   call f_conf%get_or_die("uniform_int_answer",integer_array)
   x_check = integer_array
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
!> Test uniform long int distribution
!

integer(c_int32_t) function c_test_uniform_long(c_conf) bind(c,name='test_uniform_long_f')
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, seed, i
integer(c_int64_t) :: range(2)
integer(c_int64_t), allocatable :: x(:), x_check(:), integer_array(:)
integer(c_int64_t) :: minv, maxv
character(len=*), parameter :: myname_="test_uniform_int"
character(max_string) :: err_msg, ans
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("N",n)
call f_conf%get_or_die("seed",seed)

if (f_conf%get_size("uniform_int_range") == 2) then
   call f_conf%get_or_die("uniform_int_range",integer_array)
   range = integer_array
   minv = range(1)
   maxv = range(2)
else
   write(err_msg,*) myname_ // "error reading range"
   call abor1_ftn(err_msg)
endif

!> Compute random vector
allocate(x(n))
call uniform_distribution(x, minv, maxv, seed, reset=.true.)

! write result for inclusion in config file
err_msg = achar(10) // myname_ // " Testing oops::util::random.F90 Uniform Long Distribution: "
call fckit_log%info(err_msg)
do i = 1, n
   write(ans,'(I12)') x(i)
   call fckit_log%info("   - " // adjustl(ans))
enddo

if (f_conf%get_size("uniform_int_answer") == n) then
   allocate(x_check(n))
   call f_conf%get_or_die("uniform_int_answer",integer_array)
   x_check = integer_array
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

c_test_uniform_long = 0

do i = 1, n
   if (x(i) /= x_check(i)) then
      c_test_uniform_long = 1
      write(err_msg,*) myname_ // " mismatch: ", x(i), x_check(i)
      call fckit_log%error(err_msg)
   endif
enddo

! clean up
deallocate(x,x_check)

end function c_test_uniform_long

!-------------------------------------------------------------------------------
!> Test normal real distribution
!

integer(c_int32_t) function c_test_normal_real(c_conf) bind(c,name='test_normal_real_f')
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, nx, ny, seed, i, j, k, idx
real(c_float) :: mean, sdev, tol, err
real(c_float), allocatable :: x(:), x_check(:), real_array(:)
real(c_float), allocatable :: x_2d(:,:), x_3d(:,:,:)
character(len=*), parameter :: myname_="test_normal_real"
character(max_string) :: err_msg, ans
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("N",n)
call f_conf%get_or_die("seed",seed)

mean = 0.0
if (f_conf%has("normal_mean")) call f_conf%get_or_die("normal_mean",mean)
sdev = 1.0
if (f_conf%has("normal_sdev")) call f_conf%get_or_die("normal_sdev",sdev)

!> Compute random vector
allocate(x(n))
call normal_distribution(x, mean, sdev, seed, reset = .true.)

! write result for inclusion in config file
err_msg = achar(10) // myname_ // " Testing oops::util::random.F90 Normal Real Distribution: "
call fckit_log%info(err_msg)
do i = 1, n
   write(ans,'(F30.16)') x(i)
   call fckit_log%info("   - " // adjustl(ans))
enddo

if (f_conf%get_size("normal_real_answer") == n) then
   allocate(x_check(n))
   call f_conf%get_or_die("normal_real_answer",real_array)
   x_check = real_array
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

! The tolerance is based on the precision of the data type but multiply by 10 to
! allow for roundoff error in the last digit during normalization
tol = 10.0_kind_real * epsilon(mean)

c_test_normal_real = 0


call fckit_log%info("Testing 1D normal distribution")
do i = 1, n
   err = abs(x(i)-x_check(i))
   if (err > tol*abs(x_check(i))) then
       c_test_normal_real = 1
       write(err_msg,*) myname_ // " 1D error exceeds tolerence: ", err, tol
       call fckit_log%error(err_msg)
   endif
enddo

! Now test multidimensional distributions
call fckit_log%info("Testing 2D normal distribution")
nx = n / 2
ny = n / nx
allocate(x_2d(nx,ny))
call normal_distribution(x_2d, mean, sdev, seed, reset = .true.)

idx = 0
do j = 1, ny
   do i = 1, nx
      idx = idx + 1
      err = abs(x_2d(i,j)-x_check(idx))
      if (err > tol*abs(x_check(i))) then
          c_test_normal_real = 1
          write(err_msg,*) myname_ // " 2D error exceeds tolerence: ", err, tol
          call fckit_log%error(err_msg)
      endif
   enddo
enddo

call fckit_log%info("Testing 3D normal distribution")
allocate(x_3d(nx,1,ny))
call normal_distribution(x_3d, mean, sdev, seed, reset = .true.)

idx = 0
do j = 1, ny
   do i = 1, nx
      idx = idx + 1
      err = abs(x_3d(i,1,j)-x_check(idx))
      if (err > tol*abs(x_check(i))) then
          c_test_normal_real = 1
          write(err_msg,*) myname_ // " 3D error exceeds tolerence: ", err, tol
          call fckit_log%error(err_msg)
      endif
   enddo
enddo

! clean up
deallocate(x,x_2d,x_3d,x_check)

end function c_test_normal_real

!-------------------------------------------------------------------------------
!> Test normal double distribution
!

integer(c_int32_t) function c_test_normal_double(c_conf) bind(c,name='test_normal_double_f')
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: n, nx, ny, seed, i, j, k, idx
real(c_double) :: mean, sdev, tol, err
real(c_double), allocatable :: x(:), x_check(:), real_array(:)
real(c_double), allocatable :: x_2d(:,:), x_3d(:,:,:)
character(len=*), parameter :: myname_="test_normal_double"
character(max_string) :: err_msg, ans
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("N",n)
call f_conf%get_or_die("seed",seed)

mean = 0.0
if (f_conf%has("normal_mean")) call f_conf%get_or_die("normal_mean",mean)
sdev = 1.0
if (f_conf%has("normal_sdev")) call f_conf%get_or_die("normal_sdev",sdev)

!> Compute random vector
allocate(x(n))
call normal_distribution(x, mean, sdev, seed, reset = .true.)

! write result for inclusion in config file
err_msg = achar(10) // myname_ // " Testing oops::util::random.F90 Normal Double Distribution: "
call fckit_log%info(err_msg)
do i = 1, n
   write(ans,'(F30.16)') x(i)
   call fckit_log%info("   - " // adjustl(ans))
enddo

if (f_conf%get_size("normal_double_answer") == n) then
   allocate(x_check(n))
   call f_conf%get_or_die("normal_double_answer",real_array)
   x_check = real_array
else
   write(err_msg,*) myname_ // "error reading answer"
   call abor1_ftn(err_msg)
endif

! The tolerance is based on the precision of the data type but multiply by 10 to
! allow for roundoff error in the last digit during normalization
tol = 10.0_kind_real * epsilon(mean)

c_test_normal_double = 0

call fckit_log%info("Testing 1D normal distribution")
do i = 1, n
   err = abs(x(i)-x_check(i))
   if (err > tol*abs(x_check(i))) then
       c_test_normal_double = 1
       write(err_msg,*) myname_ // " error exceeds tolerence: ", err, tol
       call fckit_log%error(err_msg)
   endif
enddo

! Now test multidimensional distributions
call fckit_log%info("Testing 2D normal distribution")
nx = n / 2
ny = n / nx
allocate(x_2d(nx,ny))
call normal_distribution(x_2d, mean, sdev, seed, reset = .true.)

idx = 0
do j = 1, ny
   do i = 1, nx
      idx = idx + 1
      err = abs(x_2d(i,j)-x_check(idx))
      if (err > tol*abs(x_check(i))) then
          c_test_normal_double = 1
          write(err_msg,*) myname_ // " 2D error exceeds tolerence: ", err, tol
          call fckit_log%error(err_msg)
      endif
   enddo
enddo

call fckit_log%info("Testing 3D normal distribution")
allocate(x_3d(nx,1,ny))
call normal_distribution(x_3d, mean, sdev, seed, reset = .true.)

idx = 0
do j = 1, ny
   do i = 1, nx
      idx = idx + 1
      err = abs(x_3d(i,1,j)-x_check(idx))
      if (err > tol*abs(x_check(i))) then
          c_test_normal_double = 1
          write(err_msg,*) myname_ // " 3D error exceeds tolerence: ", err, tol
          call fckit_log%error(err_msg)
      endif
   enddo
enddo

! clean up
deallocate(x,x_2d,x_3d,x_check)

end function c_test_normal_double

!-------------------------------------------------------------------------------

end module test_random
