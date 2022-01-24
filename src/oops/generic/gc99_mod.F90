!
! (C) Copyright 2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Gaspari-Cohn (1999) function interface

module gc99_mod

use, intrinsic :: iso_c_binding, only : c_double

implicit none
private
public gc99

interface
  function c_gc99(value_in) bind(c,name='gc99') result(value_out)
  use iso_c_binding, only: c_double
  real(c_double) :: value_in
  real(c_double) :: value_out
  end function c_gc99
end interface

contains

function gc99(value_in) result(value_out)
  real(c_double), intent(in) :: value_in
  real(c_double) :: value_out
  value_out = c_gc99(value_in)
end function gc99

end module gc99_mod
