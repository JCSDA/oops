! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Fortran module handling the interface to utils::DateTime

module datetime_mod

use, intrinsic :: iso_c_binding
use string_f_c_mod  ! TODO: replace with fckit equivalent (in fckit_c_interop_module)
use duration_mod

implicit none
private
public datetime, datetime_create, datetime_set, datetime_delete, &
     & assignment(=), c_f_datetime, datetime_to_string, &
     & operator(<), operator(<=), operator(>=), operator(>), &
     & datetime_update, datetime_diff, &
     & datetime_to_ifs, datetime_from_ifs

!>  Derived type encapsulating a C++ DateTime pointer.

type datetime
  private
  type(c_ptr) :: ptr = c_null_ptr
end type datetime

interface assignment(=)
  module procedure datetime_assign
end interface assignment(=)

interface operator(<)
  module procedure datetime_lt
end interface operator(<)

interface operator(<=)
  module procedure datetime_le
end interface operator(<=)

interface operator(>=)
  module procedure datetime_ge
end interface operator(>=)

interface operator(>)
  module procedure datetime_gt
end interface operator(>)

#include "datetime.intfb.h"

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!>  Create a DateTime from an ISO8601 string.

subroutine datetime_create(fstring, self)
use fckit_log_module, only : fckit_log
implicit none
type(datetime), intent(out)  :: self
character(len=*), intent(in) :: fstring
character(kind=c_char,len=1), allocatable :: cstring(:)
character(len=20) :: cltest
character(len=200) :: record

if (c_associated(self%ptr)) then
  call datetime_to_string(self, cltest)
  write(record,*) 'datetime_create self already is: ',cltest
  call fckit_log%error(record)
  call abor1_ftn('datetime_create: ptr already associated')
endif

call f_c_string(fstring, cstring)
self%ptr = c_datetime_construct(cstring)

end subroutine datetime_create

!-------------------------------------------------------------------------------

!> Delete a DateTime

subroutine datetime_delete(self)
implicit none
type(datetime), intent(inout) :: self

if (c_associated(self%ptr)) then
  call c_datetime_destruct(self%ptr)
  self%ptr = c_null_ptr
else
  call abor1_ftn('datetime_delete: nothing to destruct')
endif

end subroutine datetime_delete

!-------------------------------------------------------------------------------

!> Copy a DateTime

subroutine datetime_assign(self, other)
implicit none
type(datetime), intent(out) :: self
type(datetime), intent(in)  :: other
character(kind=c_char,len=1) :: cstring(21)

call c_datetime_string(other%ptr, cstring)

if (c_associated(self%ptr)) then
  call c_datetime_set(cstring, self%ptr)
else
  self%ptr = c_datetime_construct(cstring)
endif

end subroutine datetime_assign

!-------------------------------------------------------------------------------

!> Convert a C++ DateTime to a Fortran datetime

subroutine c_f_datetime(cdt, fdt)
implicit none
type(C_PTR), intent(in)     :: cdt
type(datetime), intent(out) :: fdt
fdt%ptr = cdt
end subroutine c_f_datetime

!-------------------------------------------------------------------------------

!>  Create a DateTime from an ISO8601 string.

subroutine datetime_set(fstring, self)
implicit none
type(datetime), intent(inout) :: self
character(len=*), intent(in)  :: fstring
character(kind=c_char,len=1), allocatable :: cstring(:)

call f_c_string(fstring, cstring)
call c_datetime_set(cstring, self%ptr)

end subroutine datetime_set

!-------------------------------------------------------------------------------

!> Get DateTime as string

subroutine datetime_to_string(fdt, fstring)
implicit none
type(datetime), intent(in)      :: fdt
character(len=*), intent(inout) :: fstring
character(kind=c_char,len=1) :: cstring(21)

call c_datetime_string(fdt%ptr, cstring)
call c_f_string(cstring, fstring)

end subroutine datetime_to_string

!-------------------------------------------------------------------------------

!> Get Date and Time as integers for IFS :-(

subroutine datetime_to_ifs(fdt, kdate, ksecs)
implicit none
type(datetime), intent(in)       :: fdt
integer(kind=c_int), intent(out) :: kdate
integer(kind=c_int), intent(out) :: ksecs

integer(kind=c_int64_t) :: idate

call c_datetime_getints(fdt%ptr, idate, ksecs)
kdate = idate

end subroutine datetime_to_ifs

!-------------------------------------------------------------------------------

!> Set Date and Time as integers for IFS :-(

subroutine datetime_from_ifs(fdt, kdate, ksecs)
implicit none
type(datetime), intent(inout)   :: fdt
integer(kind=c_int), intent(in) :: kdate
integer(kind=c_int), intent(in) :: ksecs

integer(kind=c_int64_t) :: idate

idate = kdate
call c_datetime_setints(fdt%ptr, idate, ksecs)

end subroutine datetime_from_ifs

!-------------------------------------------------------------------------------

!> Advance a DateTime

subroutine datetime_update(self, dt)
implicit none
type(datetime), intent(inout) :: self
type(duration), intent(in)    :: dt
call c_datetime_update(self%ptr, duration_seconds(dt))
end subroutine datetime_update

!-------------------------------------------------------------------------------

subroutine datetime_diff(t1, t2, dt)
implicit none
type(datetime), intent(in)  :: t1, t2
type(duration), intent(out) :: dt
dt = c_datetime_diff(t1%ptr, t2%ptr)
end subroutine datetime_diff

!-------------------------------------------------------------------------------

logical function datetime_lt(dt1, dt2)
implicit none
type(datetime), intent(in) :: dt1, dt2
datetime_lt = (c_datetime_diff(dt1%ptr, dt2%ptr) < 0)
end function datetime_lt

!-------------------------------------------------------------------------------

logical function datetime_le(dt1, dt2)
implicit none
type(datetime), intent(in) :: dt1, dt2
datetime_le = (c_datetime_diff(dt1%ptr, dt2%ptr) <= 0)
end function datetime_le

!-------------------------------------------------------------------------------

logical function datetime_ge(dt1, dt2)
implicit none
type(datetime), intent(in) :: dt1, dt2
datetime_ge = (c_datetime_diff(dt1%ptr, dt2%ptr) >= 0)
end function datetime_ge

!-------------------------------------------------------------------------------

logical function datetime_gt(dt1, dt2)
implicit none
type(datetime), intent(in) :: dt1, dt2
datetime_gt = (c_datetime_diff(dt1%ptr, dt2%ptr) > 0)
end function datetime_gt

!-------------------------------------------------------------------------------

end module datetime_mod
