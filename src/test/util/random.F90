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

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

integer(c_int32_t) function c_test_uniform_real(c_conf) bind(c,name='test_uniform_real_f')  
implicit none
type(c_ptr), intent(in) :: c_conf
integer :: N, seed

N = config_get_int(c_conf, "N")
seed = config_get_int(c_conf, "seed")

write(*,*) "MSM N: ", N
write(*,*) "MSM seed: ", seed

! if (passes)
c_test_uniform_real = 0
! else
!c_test_uniform_real = 1

end function c_test_uniform_real

!-------------------------------------------------------------------------------

end module test_random
