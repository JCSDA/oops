! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_constants_mod

use kinds

implicit none

private
public :: pi,deg_to_rad,rad_to_deg
public :: req,omega,g
public :: domain_zonal,domain_meridional,domain_depth,f0,bet0,dlogtheta
public :: heating_amplitude,heating_scale

! Mathematical parameters
real(kind_real),parameter :: pi = acos(-1.0)                      !< Pi
real(kind_real),parameter :: deg_to_rad = pi/180.0_kind_real      !< Degrees to radians
real(kind_real),parameter :: rad_to_deg = 180.0_kind_real/pi      !< Radians to degrees

! Geophysical parameters
real(kind_real),parameter :: req = 6371229.0_kind_real            !< Earth radius at equator (m)
real(kind_real),parameter :: omega = 7.2921e-5_kind_real          !< Rotation rate of the Earth (rad/s)
real(kind_real),parameter :: g = 10.0_kind_real                   !< Gravity (m^2 s^{-2})

! Domain parameters
real(kind_real),parameter :: domain_zonal = 2.0*pi*req            !< Model domain in zonal direction (m)
real(kind_real),parameter :: domain_meridional = 0.5*pi*req       !< Model domain in meridional direction (m)
real(kind_real),parameter :: domain_depth = 1.0e4_kind_real       !< Model depth (m)
real(kind_real),parameter :: f0 = 1.0e-4_kind_real                !< Coriolis parameter at the center of the domain
real(kind_real),parameter :: bet0 = 1.5e-11_kind_real             !< Meridional gradient of f (s^{-1} m^{-1})
real(kind_real),parameter :: dlogtheta = 0.1_kind_real            !< Difference in log(pot. T) across boundary

! Heating term parameters
real(kind_real),parameter :: heating_amplitude = 5.0e-5_kind_real !< Heating term amplitude
real(kind_real),parameter :: heating_scale = 1000e3_kind_real     !< Heating term scale (m)

end module qg_constants_mod
