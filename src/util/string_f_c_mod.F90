! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Convert Fortran strings to/from C++

module string_f_c_mod

use, intrinsic :: iso_c_binding, only : c_char, c_null_char

implicit none
private
public f_c_string, c_f_string

contains

! ------------------------------------------------------------------------------
!> Convert Fortran string to allocatable C++ string.

subroutine f_c_string(fstring, cstring)
character(len=*), intent(in)               :: fstring
character(kind=c_char, len=1), allocatable, intent(inout) :: cstring(:)
integer :: jj

if (allocated(cstring)) deallocate(cstring)
allocate(cstring(len(trim(fstring))+1))

do jj=1,len(fstring)
  cstring(jj) = fstring(jj:jj)
enddo
cstring(len(fstring)+1) = c_null_char

end subroutine f_c_string

! ------------------------------------------------------------------------------
!> Convert C++ string to Fortran string.

subroutine c_f_string(cstring, fstring)
character(kind=c_char, len=1), intent(in) :: cstring(:)
character(len=*), intent(inout)           :: fstring
integer :: jj

if (len(fstring)<size(cstring)-1) call abor1_ftn('c_f_string: fstring too short')

fstring=''
do jj=1,size(cstring)-1
  if (IACHAR(cstring(jj)) == 0) EXIT
  fstring(jj:jj) = cstring(jj)
enddo

end subroutine c_f_string

! ------------------------------------------------------------------------------

end module string_f_c_mod
