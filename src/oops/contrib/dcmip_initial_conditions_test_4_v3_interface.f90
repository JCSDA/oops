! (C) Crown Copyright 2023 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

subroutine test4_baroclinic_wave_c(moist_c,X_c,lon_c,lat_c,p_c,z_c, &
    & zcoords_c,u_c,v_c,w_c,t_c,phis_c,ps_c,rho_c,q_c,q1_c,q2_c) &
    & bind(c,name="test4_baroclinic_wave_f90")

use, intrinsic :: iso_c_binding, only: c_int,c_double
use :: dcmip_initial_conditions_test_4, only: test4_baroclinic_wave

implicit none

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
integer(c_int), intent(in)  :: moist_c      ! Moist (1) or non-moist (0) test case

real(c_double), intent(in)  :: &
            lon_c,        & ! Longitude (radians)
            lat_c,        & ! Latitude (radians)
            z_c,          & ! Height (m)
            X_c             ! Scale factor, not used in this version since unscaled EPV is selected

real(c_double), intent(inout) :: p_c        ! Pressure  (Pa)

integer(c_int), intent(in) :: zcoords_c     ! 0 or 1 see below

real(c_double), intent(out) :: &
            u_c,          & ! Zonal wind (m s^-1)
            v_c,          & ! Meridional wind (m s^-1)
            w_c,          & ! Vertical Velocity (m s^-1)
            t_c,          & ! Temperature (K)
            phis_c,       & ! Surface Geopotential (m^2 s^-2)
            ps_c,         & ! Surface Pressure (Pa)
            rho_c,        & ! density (kg m^-3)
            q_c,          & ! Specific Humidity (kg/kg)
            q1_c,         & ! Tracer q1 - Potential temperature (kg/kg)
            q2_c            ! Tracer q2 - Ertel's potential vorticity (kg/kg)

call test4_baroclinic_wave(moist_c,X_c,lon_c,lat_c,p_c,z_c, &
    & zcoords_c,u_c,v_c,w_c,t_c,phis_c,ps_c,rho_c,q_c,q1_c,q2_c)

end subroutine test4_baroclinic_wave_c

