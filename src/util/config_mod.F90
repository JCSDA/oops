! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Fortran interface to Config.

module config_mod

use, intrinsic :: iso_c_binding
use string_f_c_mod
use kinds

implicit none

private
public config_element_exists, config_get_int, config_get_real, config_get_string

#include "config.intfb.h"

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Check if an element exists

logical function config_element_exists(c_dom,query)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  character(len=*), intent(in) :: query

  character(kind=c_char,len=1), allocatable :: c_query(:)

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query
  config_element_exists = LOGICAL(c_config_element_exists(c_dom,c_query))
  deallocate(c_query)
end function config_element_exists

!-------------------------------------------------------------------------------

!>  Return data from a Config element as a integer.

integer function config_get_int(c_dom,query)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  character(len=*), intent(in) :: query

  character(kind=c_char,len=1), allocatable :: c_query(:)

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query
  config_get_int = c_config_get_data_as_int(c_dom,c_query)
  deallocate(c_query)
end function config_get_int

!-------------------------------------------------------------------------------

!>  Return data from a Config element as a real.

real(kind=kind_real) function config_get_real(c_dom,query)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  character(len=*), intent(in) :: query

  character(kind=c_char,len=1), allocatable :: c_query(:)

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query
  config_get_real = c_config_get_data_as_double(c_dom,c_query)
  deallocate(c_query)
end function config_get_real

!-------------------------------------------------------------------------------

!>  Return the contents of a Config element as a string

function config_get_string(c_dom,length,query)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  integer, intent(in) :: length
  character(len=*), intent(in) :: query
  character(len=length) ::  config_get_string

  character(kind=c_char,len=1), allocatable :: c_query(:), c_data(:)
  integer nchars

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query
  nchars = c_config_get_data_length(c_dom,c_query)
  if (nchars > length) &
    call abor1_ftn('config_get_string: return argument too short')

  if (nchars>0) then
    allocate(c_data(nchars+1))
    call c_config_get_data(c_dom,c_query,c_data)
    call c_f_string(c_data, config_get_string)
  else
    call abor1_ftn('config_get_string: element not found')
!    config_get_string = ""
  endif

  deallocate(c_query)
end function config_get_string

!-------------------------------------------------------------------------------

end module config_mod
