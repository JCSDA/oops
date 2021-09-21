! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Convert Fortran strings to/from C++

module string_utils

use fckit_configuration_module, only: fckit_configuration

implicit none
private
public swap_name_member, replace_string

!---------------------------------------------------------------------------------------------------
! Fortran utilities
!---------------------------------------------------------------------------------------------------
contains

! --------------------------------------------------------------------------------------------------
!> Changes the string '%{member}%' to '{iii}'

subroutine swap_name_member(f_conf, str, ndigits)
type(fckit_configuration), intent(in)        :: f_conf
character(len=:), allocatable, intent(inout) :: str
integer, optional, intent(in)                :: ndigits

character(len=:), allocatable :: str2, member_folder
character(len=6)  :: my_fmt
character(len=10) :: mymember_str
integer :: mymember, member_index, ldigits

ldigits = 3
if (present(ndigits)) ldigits = ndigits
if (ldigits < 1 .or. ldigits > 9) call abor1_ftn('ndigits must be >1 and <=9')

write(my_fmt, '(a, i0, a)') '(I0.', ldigits, ')'
if ( f_conf%has("member") ) then
  call f_conf%get_or_die("member", mymember)
  write(mymember_str,my_fmt) mymember

  member_folder = "%{member}%"
  member_index = index(str, member_folder)

  if ( member_index>0 ) then
    str2 = str(1:member_index-1) // trim(mymember_str) // &
      str(member_index+len(member_folder):len(str))
    deallocate(str)
    str=str2
  end if
end if

end subroutine swap_name_member

!---------------------------------------------------------------------------------------------------

function replace_string (inputstr, search, replace) result(outputstr)

implicit none
character(len=*), intent(in) :: inputstr
character(len=*), intent(in) :: search
character(len=*), intent(in) :: replace
character(len(inputstr)+100) :: outputstr

! Locals
integer :: i, nt, nr

outputstr = inputstr
nt = len_trim(search)
nr = len_trim(replace)

do
  i = index(outputstr,search(:nt)) ; if (i == 0) exit
  outputstr = outputstr(:i-1) // replace(:nr) // outputstr(i+nt:)
end do

end function replace_string

! --------------------------------------------------------------------------------------------------

end module string_utils
