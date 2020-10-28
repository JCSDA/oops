!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Fortran interface to LibOOPS

module liboops_mod

implicit none

private
public liboops_initialise
public liboops_finalise

#include "oops/util/liboops_interface.f"

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine liboops_initialise()
  implicit none

  call c_liboops_initialise()
end subroutine liboops_initialise

!-------------------------------------------------------------------------------

subroutine liboops_finalise()
  implicit none

  call c_liboops_finalise()
end subroutine liboops_finalise

!-------------------------------------------------------------------------------

end module liboops_mod
