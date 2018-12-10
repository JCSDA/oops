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
integer, parameter :: max_string=800
public config_element_exists, config_get_int, config_get_real, config_get_string, &
       config_get_string_vector

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

integer function config_get_int(c_dom,query,idefault)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  character(len=*), intent(in) :: query
  integer, intent(in), optional :: idefault
  character(kind=c_char,len=1), allocatable :: c_query(:)
  character(max_string) :: err_msg

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query
  if(LOGICAL(c_config_element_exists(c_dom,c_query))) then
    config_get_int = c_config_get_data_as_int(c_dom,c_query)
  else
    if (present(idefault)) then
      config_get_int = idefault
    else
      write(err_msg,*) "config_get_int: ", trim(query), " does not exist in config and no default"
      call abor1_ftn(err_msg)
    endif
  endif
  deallocate(c_query)
end function config_get_int

!-------------------------------------------------------------------------------

!>  Return data from a Config element as a real.

real(kind=kind_real) function config_get_real(c_dom,query,rdefault)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  character(len=*), intent(in) :: query
  real(kind_real), intent(in), optional :: rdefault
  character(kind=c_char,len=1), allocatable :: c_query(:)
  character(max_string) :: err_msg

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query
  if(LOGICAL(c_config_element_exists(c_dom,c_query))) then
    config_get_real = c_config_get_data_as_double(c_dom,c_query)
  else
    if (present(rdefault)) then
      config_get_real = rdefault
    else
      write(err_msg,*) "config_get_real: ", trim(query), " does not exist in config and no default"
      call abor1_ftn(err_msg)
    endif
  endif
  deallocate(c_query)
end function config_get_real

!-------------------------------------------------------------------------------

!>  Return the contents of a Config element as a string

function config_get_string(c_dom,length,query,sdefault)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  integer,     intent(in) :: length
  character(len=*), intent(in) :: query
  character(len=length)                       ::  config_get_string
  character(len=*), intent(in), optional ::  sdefault
  character(kind=c_char,len=1), allocatable   :: c_query(:), c_data(:)
  integer nchars
  character(max_string) :: err_msg

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query
  if(LOGICAL(c_config_element_exists(c_dom,c_query))) then

     nchars = c_config_get_data_length(c_dom,c_query)
     if (nchars > length) &
       call abor1_ftn('config_get_string: return argument too short')

     allocate(c_data(nchars+1))
     call c_config_get_data(c_dom,c_query,c_data)
     call c_f_string(c_data, config_get_string)

  else

    if (present(sdefault)) then
      config_get_string = sdefault
    else
      write(err_msg,*) "config_get_string: ", trim(query), " does not exist in config and no default"
      call abor1_ftn(err_msg)
    endif

  endif

  deallocate(c_query)
  if (allocated(c_data)) deallocate(c_data)
end function config_get_string

!-------------------------------------------------------------------------------

!>  Return the contents of a Config element as a string vector

function config_get_string_vector(c_dom, length, query)
  implicit none
  type(c_ptr), intent(in) :: c_dom
  integer, intent(in) :: length
  character(len=*), intent(in) :: query
  character(len=length), allocatable ::  config_get_string_vector(:)

  character(kind=c_char,len=1), allocatable :: c_query(:), c_data(:)
  integer :: nchars, ndim, index

  !  Translate query from Fortran string to C++ char[].
  call f_c_string(query,c_query)

  ! Call C++ to process the query the dimension of string vector
  ndim = c_config_get_data_dimension(c_dom, c_query)
  if (ndim>0) then
    allocate(config_get_string_vector(ndim))
  else
    call abor1_ftn('config_get_string_vector: element not found')
  endif

  ! Call C++ to process the query the string vector element one by one
  do index = 1, ndim
    nchars = c_config_get_data_element_length(c_dom, c_query, index)
    if (nchars > length) &
      call abor1_ftn('config_get_string: return argument too short')

    if (nchars>0) then
      allocate(c_data(nchars+1))
      call c_config_get_data_element(c_dom, c_query, index, c_data)
      call c_f_string(c_data, config_get_string_vector(index))
      deallocate(c_data)
    else
      call abor1_ftn('config_get_string: element not found')
    endif
  enddo

  deallocate(c_query)
  if (allocated(c_data)) deallocate(c_data)
end function config_get_string_vector

!-------------------------------------------------------------------------------

end module config_mod
