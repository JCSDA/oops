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
use qg_obsdb_mod
use qg_obsvec_mod
use qg_wspeed_mod
use oops_variables_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed
subroutine qg_wspeed_equiv_c(c_obsdb,c_gom,c_hofx,c_bias) bind(c,name='qg_wspeed_equiv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_obsdb  !< Observation data
integer(c_int),intent(in) :: c_gom    !< GOM
integer(c_int),intent(in) :: c_hofx   !< Observation vector
real(c_double),intent(in) :: c_bias   !< Bias

! Local variables
type(qg_obsdb),pointer :: obsdb
type(qg_gom),pointer  :: gom
type(qg_obsvec),pointer :: hofx

! Interface
call qg_obsdb_registry%get(c_obsdb,obsdb)
call qg_gom_registry%get(c_gom,gom)
call qg_obsvec_registry%get(c_hofx,hofx)

! Call Fortran
call qg_wspeed_equiv(obsdb,gom,hofx,c_bias)

end subroutine qg_wspeed_equiv_c
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - tangent linear
subroutine qg_wspeed_equiv_tl_c(c_obsdb,c_gom,c_hofx,c_traj,c_bias) bind(c,name='qg_wspeed_equiv_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_obsdb  !< Observation data
integer(c_int),intent(in) :: c_gom    !< GOM
integer(c_int),intent(in) :: c_hofx   !< Observation vector
integer(c_int),intent(in) :: c_traj   !< GOM trajectory
real(c_double),intent(in) :: c_bias   !< Bias

! Local variables
type(qg_obsdb),pointer :: obsdb
type(qg_gom),pointer :: gom
type(qg_obsvec),pointer :: hofx
type(wspeed_traj),pointer :: traj

! Interface
call qg_obsdb_registry%get(c_obsdb,obsdb)
call qg_gom_registry%get(c_gom,gom)
call qg_obsvec_registry%get(c_hofx,hofx)
call qg_wspeed_registry%get(c_traj,traj)

! Call Fortran
call qg_wspeed_equiv_tl(obsdb,gom,hofx,traj,c_bias)

end subroutine qg_wspeed_equiv_tl_c
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - adjoint
subroutine qg_wspeed_equiv_ad_c(c_obsdb,c_gom,c_hofx,c_traj,c_bias) bind(c,name='qg_wspeed_equiv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_obsdb    !< Observation data
integer(c_int),intent(in) :: c_gom      !< GOM
integer(c_int),intent(in) :: c_hofx     !< Observation vector
integer(c_int),intent(in) :: c_traj     !< GOM trajectory
real(c_double),intent(inout) :: c_bias  !< Bias

! Local variables
type(qg_obsdb),pointer :: obsdb
type(qg_gom),pointer :: gom
type(qg_obsvec),pointer :: hofx
type(wspeed_traj),pointer :: traj

! Interface
call qg_obsdb_registry%get(c_obsdb,obsdb)
call qg_gom_registry%get(c_gom,gom)
call qg_obsvec_registry%get(c_hofx,hofx)
call qg_wspeed_registry%get(c_traj,traj)

! Call Fortran
call qg_wspeed_equiv_ad(obsdb,gom,hofx,traj,c_bias)

end subroutine qg_wspeed_equiv_ad_c
! ------------------------------------------------------------------------------
!> Get wind speed trajectory
subroutine qg_wspeed_alloctraj_c(c_nobs,c_traj) bind(c,name='qg_wspeed_alloctraj_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_nobs     !< Number of observations
integer(c_int),intent(inout) :: c_traj  !< GOM trajectory

! Local variables
type(wspeed_traj),pointer :: traj
integer :: nobs

! Interface
call qg_wspeed_registry%init()
call qg_wspeed_registry%add(c_traj)
call qg_wspeed_registry%get(c_traj,traj)
nobs = c_nobs

! Call Fortran
call qg_wspeed_alloctraj(traj,nobs)

end subroutine qg_wspeed_alloctraj_c
! ------------------------------------------------------------------------------
!> Set wind speed trajectory
subroutine qg_wspeed_settraj_c(c_obsdb,c_gom,c_traj) bind(c,name='qg_wspeed_settraj_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_obsdb    !< Observation data
integer(c_int),intent(in) :: c_gom    !< GOM
integer(c_int),intent(in) :: c_traj   !< GOM trajectory

! Local variables
type(qg_obsdb),pointer :: obsdb
type(qg_gom),pointer :: gom
type(wspeed_traj),pointer :: traj

! Interface
call qg_obsdb_registry%get(c_obsdb,obsdb)
call qg_gom_registry%get(c_gom,gom)
call qg_wspeed_registry%get(c_traj,traj)

! Call Fortran
call qg_wspeed_settraj(obsdb,gom,traj)

end subroutine qg_wspeed_settraj_c
! ------------------------------------------------------------------------------
end module qg_wspeed_interface
