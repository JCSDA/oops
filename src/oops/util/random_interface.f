!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Define Fortran interface for random number generator

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

subroutine c_random_uniform_float(length, minv, maxv, seed, csort, creset, vect) &
                                  bind(C,name='random_uniform_float_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t, c_float
   implicit none
   integer(c_size_t), intent(in) :: length, csort, creset
   real(c_float), intent(in) :: minv, maxv   
   integer(c_int32_t) :: seed
   real(c_float), intent(inout) :: vect(length)
      
end subroutine c_random_uniform_float

subroutine c_random_uniform_double(length, minv, maxv, seed, csort, creset, vect) &
                                   bind(C,name='random_uniform_double_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t, c_double
   implicit none
   integer(c_size_t), intent(in) :: length, csort, creset
   real(c_double), intent(in) :: minv, maxv   
   integer(c_int32_t) :: seed
   real(c_double), intent(inout) :: vect(length)
      
end subroutine c_random_uniform_double

subroutine c_random_uniform_int(length, minv, maxv, seed, csort, creset, vect) &
                                bind(C,name='random_uniform_int_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t
   implicit none
   integer(c_size_t), intent(in) :: length, csort, creset
   integer(c_int32_t), intent(in) :: minv, maxv   
   integer(c_int32_t) :: seed
   integer(c_int32_t), intent(inout) :: vect(length)
      
end subroutine c_random_uniform_int

subroutine c_random_uniform_long(length, minv, maxv, seed, csort, creset, vect) &
                                 bind(C,name='random_uniform_long_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t, c_int64_t
   implicit none
   integer(c_size_t), intent(in) :: length, csort, creset
   integer(c_int64_t), intent(in) :: minv, maxv   
   integer(c_int32_t) :: seed
   integer(c_int64_t), intent(inout) :: vect(length)
      
end subroutine c_random_uniform_long

subroutine c_random_normal_float(length, mean, sdev, seed, creset, vect) &
                                 bind(C,name='random_normal_float_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t, c_float
   implicit none
   integer(c_size_t), intent(in) :: length, creset
   real(c_float), intent(in) :: mean, sdev
   integer(c_int32_t) :: seed
   real(c_float), intent(inout) :: vect(length)
      
end subroutine c_random_normal_float

subroutine c_random_normal_double(length, mean, sdev, seed, creset, vect) &
                                  bind(C,name='random_normal_double_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t, c_double
   implicit none
   integer(c_size_t), intent(in) :: length, creset
   real(c_double), intent(in) :: mean, sdev
   integer(c_int32_t) :: seed
   real(c_double), intent(inout) :: vect(length)
      
end subroutine c_random_normal_double

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
      
