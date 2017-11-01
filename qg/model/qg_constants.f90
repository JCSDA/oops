! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Constants for the QG model

module qg_constants

use kinds
implicit none

!--- Dimensional parameters

real(kind=kind_real),parameter :: domain_zonal=12e6_kind_real  !< model domain (m) in zonal direction
real(kind=kind_real),parameter :: domain_meridional=6.3e6_kind_real  !< meridional model domain (m) 
real(kind=kind_real),parameter :: scale_length = 1e6_kind_real !< horizontal length scale (m)
real(kind=kind_real),parameter :: ubar = 10.0_kind_real       !< typical verlocity (m/s)
real(kind=kind_real),parameter :: ubar1 = 40.0_kind_real      !< mean zonal wind in the top layer (m/s)
real(kind=kind_real),parameter :: ubar2 = 10.0_kind_real      !< mean zonal wind in the bottom layer (m/s)
real(kind=kind_real),parameter :: dlogtheta = 0.1_kind_real   !< difference in log(pot. T) across boundary
real(kind=kind_real),parameter :: g=10.0_kind_real            !< gravity (m^2 s^{-2})
real(kind=kind_real),parameter :: f0 = 1e-4_kind_real         !< Coriolis parameter at southern boundary
real(kind=kind_real),parameter :: bet0 = 1.5e-11_kind_real    !< Meridional gradient of f (s^{-1} m^{-1})
real(kind=kind_real),parameter :: horog = 2000.0_kind_real    !< height of orography (m)
real(kind=kind_real),parameter :: worog = 1000e3_kind_real    !< e-folding width of orography (m)

!--- Non-dimensional parameters

real(kind=kind_real),parameter :: u1 = ubar1/ubar
real(kind=kind_real),parameter :: u2 = ubar2/ubar
real(kind=kind_real),parameter :: bet = bet0*scale_length*scale_length/ubar
real(kind=kind_real),parameter :: rossby_number = ubar/(f0*scale_length)

end module qg_constants
