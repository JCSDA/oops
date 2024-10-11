! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Fortran module handling the interface to utils::Duration

module duration_mod

use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t, c_char
use string_f_c_mod

implicit none
private
public duration, assignment(=), c_f_duration, &
     & duration_to_string, duration_seconds

!>  Derived type encapsulating a C++ Duration pointer.

type duration
  private
  integer(kind=c_int64_t) :: seconds = 0
end type duration

interface assignment(=)
  module procedure duration_assign, duration_from_string, &
                 & duration_from_int, duration_from_int64
end interface assignment(=)

#include "duration.intfb.h"

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!>  Create a Duration from an ISO8601 string.

subroutine duration_from_string(self, fstring)
implicit none
type(duration), intent(out)  :: self
character(len=*), intent(in) :: fstring
character(kind=c_char, len=1), allocatable :: cstring(:)

call f_c_string(fstring, cstring)
self%seconds = c_duration_str_int(cstring)

end subroutine duration_from_string

!-------------------------------------------------------------------------------

subroutine duration_from_int(self, ksecs)
implicit none
type(duration), intent(inout) :: self
integer, intent(in)           :: ksecs
self%seconds = ksecs
end subroutine duration_from_int

!-------------------------------------------------------------------------------

subroutine duration_from_int64(self, ksecs)
implicit none
type(duration), intent(inout)       :: self
integer(kind=c_int64_t), intent(in) :: ksecs
self%seconds = ksecs
end subroutine duration_from_int64

!-------------------------------------------------------------------------------

subroutine duration_assign(self, other)
implicit none
type(duration), intent(inout) :: self
type(duration), intent(in)    :: other
self%seconds = other%seconds
end subroutine duration_assign

!-------------------------------------------------------------------------------

!> Convert a C++ Duration to a Fortran duration

subroutine c_f_duration(cdt, fdt)
implicit none
type(c_ptr), intent(in)     :: cdt
type(duration), intent(out) :: fdt
fdt%seconds = c_duration_int(cdt)
end subroutine c_f_duration

!-------------------------------------------------------------------------------

integer(kind=c_int64_t) function duration_seconds(self)
implicit none
type(duration), intent(in) :: self
duration_seconds = self%seconds
end function duration_seconds

!-------------------------------------------------------------------------------

subroutine duration_to_string(self, fstring)
implicit none
type(duration), intent(in)      :: self
character(len=*), intent(inout) :: fstring
character(kind=c_char, len=1) :: cstring(21)

call c_duration_int_str(self%seconds, cstring)
call c_f_string(cstring, fstring)

end subroutine duration_to_string

!-------------------------------------------------------------------------------

end module duration_mod
