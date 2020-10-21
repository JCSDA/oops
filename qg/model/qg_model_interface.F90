! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_model_interface

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use qg_fields_mod
use qg_model_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup model
subroutine qg_model_setup_c(c_key_self,c_conf) bind (c,name='qg_model_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Model configuration
type(c_ptr),value,intent(in) :: c_conf        !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_model_config),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
call qg_model_registry%init()
call qg_model_registry%add(c_key_self)
call qg_model_registry%get(c_key_self,self)

! Call Fortran
call qg_model_setup(self,f_conf)
   
end subroutine qg_model_setup_c
! ------------------------------------------------------------------------------
!> Delete the QG model
subroutine qg_delete_c(c_key_conf) bind (c,name='qg_model_delete_f90')

implicit none
integer(c_int),intent(inout) :: c_key_conf !< Model configuration

call qg_model_registry%remove(c_key_conf)

end subroutine qg_delete_c
! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model
subroutine qg_model_propagate_c(c_key_conf,c_key_state) bind(c,name='qg_model_propagate_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_conf  !< Model configuration
integer(c_int),intent(in) :: c_key_state !< State fields

! Local variables
type(qg_model_config),pointer :: conf
type(qg_fields),pointer :: fld

! Interface
call qg_model_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_state,fld)

! Call Fortran
call qg_model_propagate(conf,fld)

end subroutine qg_model_propagate_c
! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model - tangent linear
subroutine qg_model_propagate_tl_c(c_key_conf,c_key_traj,c_key_incr) bind(c,name='qg_model_propagate_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_conf !< Model configuration
integer(c_int),intent(in) :: c_key_traj !< Trajectory fields
integer(c_int),intent(in) :: c_key_incr !< Increment fields

! Local variables
type(qg_model_config),pointer :: conf
type(qg_fields),pointer :: traj
type(qg_fields),pointer :: fld

! Interface
call qg_model_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_traj,traj)
call qg_fields_registry%get(c_key_incr,fld)

! Call Fortran
call qg_model_propagate_tl(conf,traj,fld)

end subroutine qg_model_propagate_tl_c
! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model - adjoint
subroutine qg_model_propagate_ad_c(c_key_conf,c_key_traj,c_key_incr) bind(c,name='qg_model_propagate_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_conf !< Model configuration
integer(c_int),intent(in) :: c_key_traj !< Trajectory fields
integer(c_int),intent(in) :: c_key_incr !< Increment fields

! Local variables
type(qg_model_config),pointer :: conf
type(qg_fields),pointer :: fld
type(qg_fields),pointer :: traj

! Interface
call qg_model_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_traj,traj)
call qg_fields_registry%get(c_key_incr,fld)

! Call Fortran
call qg_model_propagate_ad(conf,traj,fld)

end subroutine qg_model_propagate_ad_c
! ------------------------------------------------------------------------------
end module qg_model_interface
