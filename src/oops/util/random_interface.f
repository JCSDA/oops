!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Define Fortran interface for random number generator

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

subroutine c_random_uniform_float(length, minv, maxv, seed, vect) &
                                  bind(C,name='random_uniform_float_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t, c_float
   implicit none
   integer(c_size_t), intent(in) :: length
   real(c_float), intent(in) :: minv, maxv   
   integer(c_int32_t) :: seed
   real(c_float), intent(inout) :: vect(length)
      
end subroutine c_random_uniform_float

subroutine c_random_uniform_double(length, minv, maxv, seed, vect) &
                                   bind(C,name='random_uniform_float_f')     
   use, intrinsic :: iso_c_binding, only : c_size_t, c_int32_t, c_double
   implicit none
   integer(c_size_t), intent(in) :: length
   real(c_double), intent(in) :: minv, maxv   
   integer(c_int32_t) :: seed
   real(c_double), intent(inout) :: vect(length)
      
end subroutine c_random_uniform_double

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
      
