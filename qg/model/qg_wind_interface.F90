! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_wind_interface

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use qg_gom_mod
use qg_obsvec_mod
use qg_wind_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Get equivalent for wind
subroutine qg_wind_equiv_c(c_key_gom,c_key_hofx,c_bias) bind(c,name='qg_wind_equiv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
real(c_double),intent(in) :: c_bias(2)  !< Bias

! Local variables
type(qg_gom),pointer  :: gom
type(qg_obsvec),pointer :: hofx

! Interface
call qg_gom_registry%get(c_key_gom,gom) 
call qg_obsvec_registry%get(c_key_hofx,hofx)

! Call Fortran
call qg_wind_equiv(gom,hofx,c_bias)

end subroutine qg_wind_equiv_c
! ------------------------------------------------------------------------------
!> Get equivalent for wind - tangent linear
subroutine qg_wind_equiv_tl_c(c_key_gom,c_key_hofx,c_bias) bind(c,name='qg_wind_equiv_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom  !< GOM
integer(c_int),intent(in) :: c_key_hofx !< Observation vector
real(c_double),intent(in) :: c_bias(2)  !< Bias

! Local variables
type(qg_gom),pointer  :: gom
type(qg_obsvec),pointer :: hofx

! Interface
call qg_gom_registry%get(c_key_gom,gom) 
call qg_obsvec_registry%get(c_key_hofx,hofx)

! Call Fortran
call qg_wind_equiv(gom,hofx,c_bias)

end subroutine qg_wind_equiv_tl_c
! ------------------------------------------------------------------------------
!> Get equivalent for wind - adjoint
subroutine qg_wind_equiv_ad_c(c_key_gom,c_key_hofx,c_bias) bind(c,name='qg_wind_equiv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom    !< GOM
integer(c_int),intent(in) :: c_key_hofx   !< Observation vector
real(c_double),intent(inout) :: c_bias(2) !< Bias

! Local variables
type(qg_gom),pointer  :: gom
type(qg_obsvec),pointer :: hofx

! Interface
call qg_gom_registry%get(c_key_gom,gom)
call qg_obsvec_registry%get(c_key_hofx,hofx)

! Call Fortran
call qg_wind_equiv_ad(gom,hofx,c_bias)

end subroutine qg_wind_equiv_ad_c
!------------------------------------------------------------------------------
end module qg_wind_interface
