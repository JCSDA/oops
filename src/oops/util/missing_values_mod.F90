!
! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module missing_values_mod

use, intrinsic :: iso_c_binding
use kinds

implicit none

private
public missing_value

!>  Define interface for C++ code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------
real(kind=c_double) function c_missing_value_flt() bind(C,name='missing_value_flt_f')
  use, intrinsic :: iso_c_binding
  implicit none
end function c_missing_value_flt
!-------------------------------------------------------------------------------
real(kind=c_double) function c_missing_value_dbl() bind(C,name='missing_value_dbl_f')
  use, intrinsic :: iso_c_binding
  implicit none
end function c_missing_value_dbl
!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

!>  Return the missing value indicator

interface missing_value
!  module procedure missing_value_int32
!  module procedure missing_value_int64
  module procedure missing_value_float
  module procedure missing_value_double
!  module procedure missing_value_kind_real
end interface
   
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

real(c_float) function missing_value_float(fmvi)
  implicit none
  real(c_float), intent(inout) :: fmvi
  missing_value_float = c_missing_value_flt()
end function missing_value_float

!-------------------------------------------------------------------------------

real(c_double) function missing_value_double(dmvi)
  implicit none
  real(c_double), intent(inout) :: dmvi
  missing_value_double = c_missing_value_dbl()
end function missing_value_double

!-------------------------------------------------------------------------------

!real(kind_real) function missing_value_kind_real(rmvi)
!  implicit none
!  real(kind_real), intent(inout) :: rmvi
!  missing_value_kind_real = c_missing_value_dbl()
!end function missing_value_kind_real

!-------------------------------------------------------------------------------

end module missing_values_mod
