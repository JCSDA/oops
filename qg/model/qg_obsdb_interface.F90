! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2017-2019 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_obsdb_interface

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use duration_mod
use fckit_log_module, only: fckit_log
use iso_c_binding
use string_f_c_mod
use qg_locs_mod
use qg_obsdb_mod
use qg_obsvec_mod
use kinds

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup observation data
subroutine qg_obsdb_setup_c(c_key_self,c_conf,c_winbgn,c_winend) bind(c,name='qg_obsdb_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Observation data
type(c_ptr),value,intent(in) :: c_conf     !< Configuration
type(c_ptr),value,intent(in) :: c_winbgn   !< Start of window
type(c_ptr),value,intent(in) :: c_winend   !< End of window

! Local variables
type(fckit_configuration) :: f_conf
type(qg_obsdb),pointer :: self
type(datetime) :: winbgn
type(datetime) :: winend

! Interface
f_conf = fckit_configuration(c_conf)
call qg_obsdb_registry%init()
call qg_obsdb_registry%add(c_key_self)
call qg_obsdb_registry%get(c_key_self,self)
call c_f_datetime(c_winbgn,winbgn)
call c_f_datetime(c_winend,winend)

! Call Fortran
call qg_obsdb_setup(self,f_conf,winbgn,winend)

end subroutine qg_obsdb_setup_c
! ------------------------------------------------------------------------------
!> Delete observation data
subroutine qg_obsdb_delete_c(c_key_self) bind(c,name='qg_obsdb_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Observation data

! Local variables
type(qg_obsdb),pointer :: self

! Interface
call qg_obsdb_registry%get(c_key_self,self)

! Call Fortran
call qg_obsdb_delete(self)

! Clear interface
call qg_obsdb_registry%remove(c_key_self)

end subroutine qg_obsdb_delete_c
! ------------------------------------------------------------------------------
!> Save observation data
subroutine qg_obsdb_save_c(c_key_self) bind(c,name='qg_obsdb_save_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Observation data

! Local variables
type(qg_obsdb),pointer :: self

! Interface
call qg_obsdb_registry%get(c_key_self,self)
call qg_obsdb_save(self)

end subroutine qg_obsdb_save_c
! ------------------------------------------------------------------------------
!> Get observation data
subroutine qg_obsdb_get_c(c_key_self,lgrp,c_grp,lcol,c_col,c_key_ovec) bind(c,name='qg_obsdb_get_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self                  !< Observation data
integer(c_int),intent(in) :: lgrp                        !< Group size
character(kind=c_char,len=1),intent(in) :: c_grp(lgrp+1) !< Group name
integer(c_int),intent(in) :: lcol                        !< Column size
character(kind=c_char,len=1),intent(in) :: c_col(lcol+1) !< Column name
integer(c_int),intent(in) :: c_key_ovec                  !< Observation vector

! Local variables
type(qg_obsdb),pointer :: self
type(qg_obsvec),pointer :: ovec
character(len=lgrp) :: grp
character(len=lcol) :: col

! Interface
call qg_obsdb_registry%get(c_key_self,self)
call c_f_string(c_grp,grp)
call c_f_string(c_col,col)
call qg_obsvec_registry%get(c_key_ovec,ovec)

! Call Fortran
call qg_obsdb_get(self,trim(grp),trim(col),ovec)

end subroutine qg_obsdb_get_c
! ------------------------------------------------------------------------------
!> Put observation data
subroutine qg_obsdb_put_c(c_key_self,lgrp,c_grp,lcol,c_col,c_key_ovec) bind(c,name='qg_obsdb_put_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self                  !< Observation data
integer(c_int),intent(in) :: lgrp                        !< Group size
character(kind=c_char,len=1),intent(in) :: c_grp(lgrp+1) !< Group name
integer(c_int),intent(in) :: lcol                        !< Column size
character(kind=c_char,len=1),intent(in) :: c_col(lcol+1) !< Column name
integer(c_int),intent(in) :: c_key_ovec                  !< Observation vector

! Local variables
type(qg_obsdb),pointer :: self
type(qg_obsvec),pointer :: ovec
character(len=lgrp) :: grp
character(len=lcol) :: col

! Interface
call qg_obsdb_registry%get(c_key_self,self)
call c_f_string(c_grp,grp)
call c_f_string(c_col,col)
call qg_obsvec_registry%get(c_key_ovec,ovec)

! Call Fortran
call qg_obsdb_put(self,trim(grp),trim(col),ovec)

end subroutine qg_obsdb_put_c
! ------------------------------------------------------------------------------
!> Get locations from observation data
subroutine qg_obsdb_locations_c(c_key_self,lgrp,c_grp,c_fields,c_times) bind(c,name='qg_obsdb_locations_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self                  !< Observation data
integer(c_int),intent(in) :: lgrp                        !< Group size
character(kind=c_char,len=1),intent(in) :: c_grp(lgrp+1) !< Group name
type(c_ptr), intent(in), value :: c_fields               !< Locations fieldset
type(c_ptr), intent(in), value :: c_times                !< times

! Local variables
type(qg_obsdb),pointer :: self
character(len=lgrp) :: grp
type(atlas_fieldset) :: fields

! Interface
call qg_obsdb_registry%get(c_key_self,self)
call c_f_string(c_grp,grp)
fields = atlas_fieldset(c_fields)

! Call Fortran
call qg_obsdb_locations(self,grp,fields,c_times)

call fields%final()

end subroutine qg_obsdb_locations_c
! ------------------------------------------------------------------------------
!> Generate observation data
subroutine qg_obsdb_generate_c(c_key_self,lgrp,c_grp,c_conf,c_bgn,c_step,ktimes,kobs) bind(c,name='qg_obsdb_generate_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self                  !< Observation data
integer(c_int),intent(in) :: lgrp                        !< Group size
character(kind=c_char,len=1),intent(in) :: c_grp(lgrp+1) !< Group name
type(c_ptr),value,intent(in) :: c_conf                   !< Configuration
type(c_ptr),value,intent(in) :: c_bgn                    !< Start time
type(c_ptr),value,intent(in) :: c_step                   !< Time-step
integer(c_int),intent(in) :: ktimes                      !< Number of time-slots
integer(c_int),intent(inout) :: kobs                     !< Number of observations

! Mocal variables
type(fckit_configuration) :: f_conf
type(qg_obsdb),pointer :: self
character(len=lgrp) :: grp
type(datetime) :: bgn
type(duration) :: step

! Interface
f_conf = fckit_configuration(c_conf)
call qg_obsdb_registry%get(c_key_self,self)
call c_f_string(c_grp,grp)
call c_f_datetime(c_bgn,bgn)
call c_f_duration(c_step,step)

! Call Fortran
call qg_obsdb_generate(self,grp,f_conf,bgn,step,ktimes,kobs)

end subroutine qg_obsdb_generate_c
! ------------------------------------------------------------------------------
!> Get observation data size
subroutine qg_obsdb_nobs_c(c_key_self,lgrp,c_grp,kobs) bind(c,name='qg_obsdb_nobs_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self                  !< Observation data
integer(c_int),intent(in) :: lgrp                        !< Group size
character(kind=c_char,len=1),intent(in) :: c_grp(lgrp+1) !< Group name
integer(c_int),intent(inout) :: kobs                     !< Number of observations

! Local variables
type(qg_obsdb),pointer :: self
character(len=lgrp) :: grp

! Interface
call qg_obsdb_registry%get(c_key_self,self)
call c_f_string(c_grp,grp)

! Call Fortran
call qg_obsdb_nobs(self,grp,kobs)

end subroutine qg_obsdb_nobs_c
! ------------------------------------------------------------------------------
end module qg_obsdb_interface
