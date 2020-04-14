! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_wspeed_interface

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use qg_gom_mod
use qg_obsoper_mod
use qg_obsvec_mod
use qg_vars_mod
use qg_wspeed_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup wind speed observation
subroutine qg_wspeed_setup_c(c_key_self,c_conf) bind(c,name='qg_wspeed_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Wind speed observation
type(c_ptr),value,intent(in) :: c_conf     !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_obsoper),pointer :: self
character(len=1) :: svars(2) = (/'u','v'/)

! Interface
f_conf = fckit_configuration(c_conf)
call qg_obsoper_registry%init()
call qg_obsoper_registry%add(c_key_self)
call qg_obsoper_registry%get(c_key_self,self)

call qg_oper_setup(self,f_conf,svars,1)

end subroutine qg_wspeed_setup_c
! ------------------------------------------------------------------------------
!> Delete wind speed observation
subroutine qg_wspeed_delete_c(c_key_self) bind(c,name='qg_wspeed_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Wind speed observation

! Clear interface
call qg_obsoper_registry%remove(c_key_self)

end subroutine qg_wspeed_delete_c
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed
subroutine qg_wspeed_equiv_c(c_key_gom,c_key_hofx,c_bias) bind(c,name='qg_wspeed_equiv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
real(c_double),intent(in) :: c_bias     !< Bias

! Local variables
type(qg_gom),pointer  :: gom
type(qg_obsvec),pointer :: hofx

! Interface
call qg_gom_registry%get(c_key_gom,gom)
call qg_obsvec_registry%get(c_key_hofx,hofx)

! Call Fortran
call qg_wspeed_equiv(gom,hofx,c_bias)

end subroutine qg_wspeed_equiv_c
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - tangent linear
subroutine qg_wspeed_equiv_tl_c(c_key_gom,c_key_hofx,c_key_traj,c_bias) bind(c,name='qg_wspeed_equiv_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
integer(c_int),intent(in) :: c_key_traj !< GOM trajectory
real(c_double),intent(in) :: c_bias     !< Bias

! Local variables
type(qg_gom),pointer :: gom,traj
type(qg_obsvec),pointer :: hofx

! Interface
call qg_gom_registry%get(c_key_gom,gom)
call qg_obsvec_registry%get(c_key_hofx,hofx)
call qg_gom_registry%get(c_key_traj,traj)

! Call Fortran
call qg_wspeed_equiv_tl(gom,hofx,traj,c_bias)

end subroutine qg_wspeed_equiv_tl_c
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - adjoint
subroutine qg_wspeed_equiv_ad_c(c_key_gom,c_key_hofx,c_key_traj,c_bias) bind(c,name='qg_wspeed_equiv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
integer(c_int),intent(in) :: c_key_traj !< GOM trajectory
real(c_double),intent(inout) :: c_bias  !< Bias

! Local variables
type(qg_gom),pointer :: gom,traj
type(qg_obsvec),pointer :: hofx

! Interface
call qg_gom_registry%get(c_key_gom,gom)
call qg_obsvec_registry%get(c_key_hofx,hofx)
call qg_gom_registry%get(c_key_traj,traj)

! Call Fortran
call qg_wspeed_equiv_ad(gom,hofx,traj,c_bias)

end subroutine qg_wspeed_equiv_ad_c
! ------------------------------------------------------------------------------
!> Get wind speed trajectory
subroutine qg_wspeed_gettraj_c(c_nobs,c_vars,c_key_traj) bind(c,name='qg_wspeed_gettraj_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_nobs              !< Number of observations
integer(c_int),dimension(*),intent(in) :: c_vars !< Variables
integer(c_int),intent(inout) :: c_key_traj       !< GOM trajectory

! Local variables
type(qg_vars) :: vars
type(qg_gom),pointer :: traj

! Interface
call qg_vars_create(vars,c_vars)
call qg_gom_registry%init()
call qg_gom_registry%add(c_key_traj)
call qg_gom_registry%get(c_key_traj,traj)

! Call Fortran
call qg_wspeed_gettraj(c_nobs,vars,traj)

end subroutine qg_wspeed_gettraj_c
! ------------------------------------------------------------------------------
!> Set wind speed trajectory
subroutine qg_wspeed_settraj_c(c_key_gom,c_key_traj) bind(c,name='qg_wspeed_settraj_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_traj !< GOM trajectory

! Local variables
type(qg_gom),pointer :: gom
type(qg_gom),pointer :: traj

! Interface
call qg_gom_registry%get(c_key_gom,gom)
call qg_gom_registry%get(c_key_traj,traj)

! Call Fortran
call qg_wspeed_settraj(gom,traj)

end subroutine qg_wspeed_settraj_c
! ------------------------------------------------------------------------------
end module qg_wspeed_interface
