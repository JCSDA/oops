!
! (C) Copyright 2019-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Fortran interface to Variables

module oops_variables_mod
use iso_c_binding, only: c_ptr
implicit none

public :: oops_variables
integer, parameter, private :: MAXVARLEN = 100

type oops_variables
private
  type(c_ptr) :: ptr
contains
  procedure, public :: destruct

  procedure, private :: push_back_string
  procedure, private :: push_back_vector

  generic, public :: push_back => push_back_string, push_back_vector

  procedure, public :: nvars
  procedure, public :: variable
  procedure, public :: varlist

  procedure, public :: has
end type

interface oops_variables
  module procedure ctor_from_ptr
  module procedure empty_ctor
  module procedure copy_ctor
end interface

private

#include "oops/base/variables_interface.f"

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

function ctor_from_ptr(ptr) result(this)
  use iso_c_binding, only: c_ptr
  implicit none
  type(oops_variables)    :: this
  type(c_ptr), intent(in) :: ptr

  this%ptr = ptr
end function ctor_from_ptr

!-------------------------------------------------------------------------------

function empty_ctor() result(this)
  type(oops_variables) :: this

  this%ptr = c_variables_empty_ctor()
end function empty_ctor

!-------------------------------------------------------------------------------

function copy_ctor(other) result(this)
  type(oops_variables) :: this
  type(oops_variables), intent(in) :: other
  integer :: jj

  this%ptr = c_variables_empty_ctor()
  do jj = 1, other%nvars()
     call this%push_back(other%variable(jj))
  enddo

end function copy_ctor

!-------------------------------------------------------------------------------

subroutine destruct(this)
  use iso_c_binding, only: c_null_ptr
  implicit none
  class(oops_variables), intent(inout) :: this

  call c_variables_destruct(this%ptr)
  this%ptr = c_null_ptr
end subroutine destruct

!-------------------------------------------------------------------------------

subroutine push_back_string(this, varname)
  use iso_c_binding, only: c_char
  use string_f_c_mod
  implicit none
  class(oops_variables), intent(in)  :: this
  character(*), intent(in) :: varname

  character(kind=c_char,len=1), allocatable :: c_vname(:)

  call f_c_string(trim(varname), c_vname)
  call c_variables_push_back(this%ptr, c_vname)
  deallocate(c_vname)

end subroutine push_back_string

!-------------------------------------------------------------------------------

subroutine push_back_vector(this, varnames)
  use iso_c_binding, only: c_char
  use string_f_c_mod
  implicit none
  class(oops_variables), intent(in)  :: this
  character(*), intent(in) :: varnames(:)

  character(kind=c_char,len=1), allocatable :: c_vname(:)
  integer :: iname

  do iname = 1, size(varnames)
     call f_c_string(trim(varnames(iname)), c_vname)
     call c_variables_push_back(this%ptr, c_vname)
     deallocate(c_vname)
  end do

end subroutine push_back_vector

!-------------------------------------------------------------------------------

integer function nvars(this)
  implicit none
  class(oops_variables), intent(in) :: this

  nvars = c_variables_size(this%ptr)
end function nvars

!-------------------------------------------------------------------------------

function variable(this, jj) result(varname)
  use iso_c_binding, only: c_char, c_size_t
  use string_f_c_mod
  implicit none

  class(oops_variables), intent(in) :: this
  integer, intent(in)  :: jj
  character(MAXVARLEN) :: varname

  integer(c_size_t) :: lcname
  character(kind=c_char,len=1), allocatable :: cname(:)

  ! Fortran indices start from 1, C++ indices start from 0
  call c_variables_getvariablelength(this%ptr, int(jj-1, c_size_t), lcname)  
  allocate(cname(lcname+1))
  call c_variables_getvariable(this%ptr, int(jj-1, c_size_t), lcname, &
                               int(size(cname), c_size_t), cname)
  call c_f_string(cname, varname)
  deallocate(cname)

end function variable

!-------------------------------------------------------------------------------

function varlist(this)
  implicit none

  class(oops_variables), intent(in) :: this
  character(MAXVARLEN), allocatable :: varlist(:)
  integer :: jj

  allocate(varlist(this%nvars()))
  
  do jj = 1, this%nvars()
     varlist(jj) = this%variable(jj)
  enddo  

end function varlist
!-------------------------------------------------------------------------------

logical function has(this, var)
  use iso_c_binding, only: c_char
  use string_f_c_mod
  implicit none

  class(oops_variables), intent(in) :: this
  character(*), intent(in) :: var

  character(kind=c_char,len=1), allocatable :: c_var(:)

  call f_c_string(trim(var), c_var)
  has = c_variables_has(this%ptr, c_var)
  deallocate(c_var)
end function has

end module oops_variables_mod
