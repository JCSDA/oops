! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_vars_mod

use config_mod
use iso_c_binding

implicit none

private
public :: nvmax,varnames
public :: qg_vars,qg_vars_create
! ------------------------------------------------------------------------------
type :: qg_vars
  logical :: lx  !< Streamfunction flag
  logical :: lq  !< Potential vorticity flag
  logical :: lu  !< Zonal wind flag
  logical :: lv  !< Meridional wind flag
end type qg_vars

integer,parameter :: nvmax = 4                                                !< Number of possible variables
character(len=1),dimension(nvmax),parameter :: varnames = (/'x','q','u','v'/) !< Possible variables names
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Create variables
subroutine qg_vars_create(self,kvars)

implicit none

! Passed variables
type(qg_vars),intent(inout) :: self       !< Variables
integer(c_int),intent(in) :: kvars(nvmax) !< Variable flags array

! Local variables
integer :: iv

! Check kvars
if (any(kvars<0).or.any(kvars>1)) call abor1_ftn ('qg_vars_create: wrong kvars values')

! Initialization
self%lx = .false.
self%lq = .false.
self%lu = .false.
self%lv = .false.

! 3d fields
do iv=1,nvmax
  if (kvars(iv)==1) then
    select case (varnames(iv))
    case ('x')
      self%lx = .true.
    case ('q')
      self%lq = .true.
    case ('u')
      self%lu = .true.
    case ('v')
      self%lv = .true.
    endselect
  endif
enddo

end subroutine qg_vars_create
! ------------------------------------------------------------------------------
end module qg_vars_mod
