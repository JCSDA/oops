! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Convert Fortran strings to/from C++

module string_f_c_mod

use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_horizontal_tab

implicit none
private
public f_c_string, c_f_string, f_c_string_vector

contains

! ------------------------------------------------------------------------------
!> Convert Fortran string to allocatable C++ string.

subroutine f_c_string(fstring, cstring)
character(len=*), intent(in)               :: fstring
character(kind=c_char, len=1), allocatable, intent(out) :: cstring(:)
integer :: jj

allocate(cstring(len_trim(fstring)+1))

do jj=1,len_trim(fstring)
  cstring(jj) = fstring(jj:jj)
enddo
cstring(len_trim(fstring)+1) = c_null_char

end subroutine f_c_string

! ------------------------------------------------------------------------------
!> Convert Fortran vector of strings to allocatable C++ string buffer.

subroutine f_c_string_vector(fstring_vec, cstring_vec)
use fckit_log_module, only : fckit_log
character(len=*), intent(in)               :: fstring_vec(:)
character(kind=c_char, len=1), intent(inout) :: cstring_vec(:)
integer :: ii,jj,idx

!> store string in buffer delineated by tabs
idx=0
do ii=1,size(fstring_vec)
   do jj=1,len_trim(fstring_vec(ii))
      idx=idx+1
      if (idx > size(cstring_vec)) &
          call abor1_ftn("oops/util::f_c_string_vector: string buffer too small")
      cstring_vec(idx) = fstring_vec(ii)(jj:jj)
   enddo
   idx=idx+1
   if (idx > size(cstring_vec)) &
       call abor1_ftn("oops/util::f_c_string_vector: string buffer too small")
   if (ii < size(fstring_vec)) then
      cstring_vec(idx) = c_horizontal_tab
   else
      cstring_vec(idx) = c_null_char
   endif
enddo

end subroutine f_c_string_vector

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
