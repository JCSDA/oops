!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Define interface for C++ LibOOPS code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

subroutine c_liboops_initialise() bind(C, name='liboops_initialise_f')
end subroutine c_liboops_initialise

!-------------------------------------------------------------------------------

subroutine c_liboops_finalise() bind(C, name='liboops_finalise_f')
end subroutine c_liboops_finalise

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------

