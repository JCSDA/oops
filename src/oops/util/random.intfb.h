! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Define interface for C++ config code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

subroutine random_normal_c(kk, pp) bind(C,name='random_normal_f')
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), intent(in) :: kk
  real(kind=c_double), intent(out) :: pp
end subroutine random_normal_c

!-------------------------------------------------------------------------------

subroutine random_uniform_c(kk, pp) bind(C,name='random_uniform_f')
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), intent(in) :: kk
  real(kind=c_double), intent(out) :: pp
end subroutine random_uniform_c

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
