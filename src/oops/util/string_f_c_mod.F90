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
public f_c_string, c_f_string, f_c_push_string_vector

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

subroutine c_push_string_to_vector(c_vec, vname) bind(C,name='push_string_to_vector_f')     
   use, intrinsic :: iso_c_binding, only : c_ptr, c_char
   implicit none
   type(c_ptr), value :: c_vec !< pointer to C++ std::vector<std::string> Object
   character(kind=c_char, len=1), intent(in) :: vname(*)      
end subroutine c_push_string_to_vector

end interface

!-------------------------------------------------------------------------------
! Fortran utilities
!-------------------------------------------------------------------------------
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
!> Push a string vector from Fortran to a C++ std::vector<std::string> object
!
subroutine f_c_push_string_vector(c_vec, vnames)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_char
  implicit none
  type(c_ptr), value, intent(in) :: c_vec !< pointer to C++ std::vector<std::string> Object
  character(len=*), intent(in) :: vnames(:)  !< names to be added to the variable list

  character(kind=c_char,len=1), allocatable :: c_vname(:)
  integer :: iname

  do iname = 1, size(vnames)
  
     call f_c_string(trim(vnames(iname)), c_vname)

     call c_push_string_to_vector(c_vec, c_vname)

     deallocate(c_vname)
     
  end do

end subroutine f_c_push_string_vector
  
! ------------------------------------------------------------------------------

end module string_f_c_mod
