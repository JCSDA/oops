! (C) Copyright 2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran interface to Logger

module logger_mod

use, intrinsic :: iso_c_binding
use fckit_c_interop_module, only : c_str_right_trim

implicit none

public :: oops_log

type :: oops_log_type
contains
  procedure, nopass, public :: info    !< Log to info channel
  procedure, nopass, public :: error   !< Log to error channel
  procedure, nopass, public :: warning !< Log to warning channel
  procedure, nopass, public :: debug   !< Log to debug channel
  procedure, nopass, public :: trace   !< Log to trace channel
  procedure, nopass, public :: stats   !< Log to stats channel
  procedure, nopass, public :: test    !< Log to test channel
end type oops_log_type

type(oops_log_type) :: oops_log

#include "logger_interface.f"

private

contains

subroutine info(msg,newl,flush)

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_info(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine info

!-------------------------------------------------------------------------------

subroutine error(msg,newl,flush)

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_error(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine error

!-------------------------------------------------------------------------------

subroutine warning(msg,newl,flush)

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_warning(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine warning

!-------------------------------------------------------------------------------

subroutine debug(msg,newl,flush)

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_debug(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine debug

!-------------------------------------------------------------------------------

subroutine trace(msg,newl,flush)

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_trace(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine trace

!-------------------------------------------------------------------------------

subroutine stats(msg,newl,flush)

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_stats(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine stats

!-------------------------------------------------------------------------------

subroutine test(msg,newl,flush)

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_test(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine test

end module logger_mod
