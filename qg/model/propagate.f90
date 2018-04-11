! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Perform a timestep of the QG model

!> This routine is called from C++ to propagate the state.
!!
!! 
!! The timestep starts by advecting the potential vorticity using
!! a semi-Lagrangian method. The potential vorticity is then inverted
!! to determine the streamfunction. Velocity components are then
!! calculated, which are used at the next timestep to advect the
!! potential vorticity.
!!
!! Note that the state visible to C++ contains potential vorticity and
!! wind components, in addition to streamfunction. This makes these
!! fields available for use in the observation operators.

subroutine propagate(flds,config)

use qg_fields
use qg_configs
use qg_constants, only: bet
use kinds

implicit none
type(qg_field),  intent(inout) :: flds
type(qg_config), intent(in)    :: config

real(kind=kind_real), allocatable :: qnew(:,:,:)

! ------------------------------------------------------------------------------

allocate(qnew(flds%geom%nx,flds%geom%ny,2))

!--- advect the potential vorticity

call advect_pv(qnew,flds%q,flds%q_north,flds%q_south,flds%u,flds%v, &
             & flds%geom%nx,flds%geom%ny,config%deltax,config%deltay,config%dt)

!--- invert the potential vorticity to determine streamfunction

call invert_pv(flds%x,qnew,flds%x_north,flds%x_south,config%rs, &
             & flds%geom%nx,flds%geom%ny,config%deltax,config%deltay, &
             & config%f1,config%f2,bet)

! -- calculate potential vorticity and wind components

flds%q(:,:,:) = qnew(:,:,:)
call zonal_wind(flds%u,flds%x,flds%x_north,flds%x_south, &
              & flds%geom%nx,flds%geom%ny,config%deltay)
call meridional_wind(flds%v,flds%x, flds%geom%nx,flds%geom%ny,config%deltax)

deallocate(qnew)
! ------------------------------------------------------------------------------
return
end subroutine propagate
