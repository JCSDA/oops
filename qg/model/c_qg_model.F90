! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------
!> Setup the QG model

!> The C++ code only knows about dimensional variables, whereas the QG model
!! is expressed in terms of non-dimensional variables.
!! Here, we query the configuration for the dimensional values, calculate the
!! equivalent non-dimensional values, and store them in the config structure.
!!
!! The convention used in the QG model is that a zero is appended to the
!! names of dimensional variables, whereas non-dimensional variable names do
!! not end in zero. (An exception to this convention are the layer depths:
!! D1 and D2.)
!!
!! The non-dimensionalisation is the usual one:
!! \f{eqnarray*}{
!!t &=& \tilde t \frac{\overline U}{L} \\\\
!!x &=& \frac{\tilde x}{L} \\\\
!!y &=& \frac{\tilde y}{L} \\\\
!!u &=& \frac{\tilde u}{\overline U} \\\\
!!v &=& \frac{\tilde v}{\overline U} \\\\
!!F_1 &=& \frac{f_0^2 L^2}{D_1 g \Delta\theta / {\overline\theta}} \\\\
!!F_2 &=& \frac{f_0^2 L^2}{D_2 g \Delta\theta / {\overline\theta}} \\\\
!!\beta &=& \beta_0 \frac{L^2}{\overline U}
!! \f}
!! where \f$ \overline U \f$, \f$ L \f$,
!! \f$ \Delta\theta / {\overline\theta} \f$, etc. are defined in
!! qg_constants.f90 .
!!
!! The Rossby number is \f$ \epsilon = {\overline U} / f_0 L \f$.

subroutine c_qg_setup(c_confspec, c_key_geom, c_key_confdata) bind (c,name='qg_setup_f90')

use qg_constants
use qg_configs
use qg_geom_mod
use iso_c_binding
use config_mod
use duration_mod
use kinds
use fckit_log_module, only : fckit_log

implicit none
type(c_ptr), intent(in)    :: c_confspec  !< pointer to object of class Config
integer(c_int), intent(in) :: c_key_geom    !< Geometry
integer(c_int), intent(inout) :: c_key_confdata  !< Key to configuration data

type(qg_config), pointer :: config
type(qg_geom), pointer :: geom

integer :: icentre, jcentre, ii, jj
real(kind=kind_real) :: distx, disty
type(duration) :: dtstep
character(len=20) :: ststep
character(len=160) :: record
character(len=30) :: otype

! ------------------------------------------------------------------------------

call qg_geom_registry%get(c_key_geom, geom)
call qg_config_registry%init()
call qg_config_registry%add(c_key_confdata)
call qg_config_registry%get(c_key_confdata, config)

config%nx  = geom%nx
config%ny  = geom%ny
write(record,*)'c_qg_setup: nx, ny=',config%nx,config%ny
call fckit_log%info(record)

config%d1  = config_get_real(c_confspec,"top_layer_depth")
config%d2  = config_get_real(c_confspec,"bottom_layer_depth")
write(record,*)'c_qg_setup: d1, d2=',config%d1,config%d2
call fckit_log%info(record)

ststep = config_get_string(c_confspec,len(ststep),"tstep")
dtstep = trim(ststep)
config%dt0 = duration_seconds(dtstep)
write(record,*)'c_qg_setup: dt0=',config%dt0
call fckit_log%info(record)

config%dt = config%dt0 * ubar/scale_length
config%f1 = f0*f0*scale_length*scale_length/(g*dlogtheta*config%d1)
config%f2 = f0*f0*scale_length*scale_length/(g*dlogtheta*config%d2)
config%rsmax = horog/(rossby_number*config%d2)
config%deltax0 = domain_zonal/real(config%nx,kind_real)
config%deltay0 = domain_meridional/real(config%ny+1,kind_real)
config%deltax = config%deltax0/scale_length
config%deltay = config%deltay0/scale_length

!--- Orography 

allocate(config%rs(config%nx,config%ny))

if (config_element_exists(c_confspec,"orography")) then
   otype = trim(config_get_string(c_confspec,len(otype),"orography"))
else
   ! This default value is for backward compatibility
   otype = "bump"
endif

write(record,*)"qg_geom_mod: orography = "//otype
call fckit_log%info(record)

if (otype == "flat") then

   config%rs=0.0_kind_real

else ! we could add other options in the future 

   !--- Gaussian hill centred on (icentre,jcentre)
   icentre=config%nx/4
   jcentre=3*config%ny/4
   do jj=1,config%ny
      do ii=1,config%nx
         distx = real(min(icentre-ii,config%nx-(icentre-ii)),kind_real) &
              & *config%deltax0
         disty = real(abs(jj-jcentre),kind_real) * config%deltay0
         config%rs(ii,jj) = config%rsmax &
              & *exp(-(distx*distx+disty*disty)/(worog*worog))
      enddo
   enddo

endif
   
return
end subroutine c_qg_setup

! ------------------------------------------------------------------------------
!> Delete the QG model

subroutine c_qg_delete(c_key_conf) bind (c,name='qg_delete_f90')

use qg_configs
use iso_c_binding

implicit none
integer(c_int), intent(inout) :: c_key_conf !< Key to configuration structure
type(qg_config), pointer :: conf

call qg_config_registry%get(c_key_conf, conf)
deallocate(conf%rs)
call qg_config_registry%remove(c_key_conf)

return
end subroutine c_qg_delete

! ------------------------------------------------------------------------------
!> Prepare for an integration of the QG model.

!> At the start of a timestep of the QG model, the state must contain
!! streamfunction, potential vorticity and wind components. The control
!! variable for the analysis, however, contains only streamfunction.
!! This routine calculates potential vorticity and wind from the
!! streamfunction, and is called before an integration of the QG model.

subroutine c_qg_prepare_integration(c_key_conf, c_key_state) &
         & bind(c,name='qg_prepare_integration_f90')

use iso_c_binding
use qg_fields
use qg_configs
use qg_constants, only: bet

implicit none
integer(c_int), intent(in) :: c_key_conf  !< Configuration structure
integer(c_int), intent(in) :: c_key_state !< Model fields

type(qg_config), pointer :: conf
type(qg_field), pointer  :: flds

call qg_field_registry%get(c_key_state,flds)
call qg_config_registry%get(c_key_conf, conf)

! -- calculate potential vorticity and wind components

call calc_pv(flds%nx,flds%ny,flds%q,flds%x,flds%x_north,flds%x_south, &
           & conf%f1,conf%f2,conf%deltax,conf%deltay,bet,conf%rs)
call zonal_wind(flds%u,flds%x,flds%x_north,flds%x_south,flds%nx,flds%ny, &
              & conf%deltay)
call meridional_wind(flds%v,flds%x,flds%nx,flds%ny,conf%deltax)

end subroutine c_qg_prepare_integration

! ------------------------------------------------------------------------------
!> Prepare for an integration of the QG model - Adjoint.

subroutine c_qg_prepare_integration_ad(c_key_conf, c_key_incr) &
           bind(c,name='qg_prepare_integration_ad_f90')

use iso_c_binding
use qg_fields
use qg_configs
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf !< Configuration structure
integer(c_int), intent(in) :: c_key_incr !< Model fields

type(qg_config), pointer :: conf
type(qg_field), pointer  :: flds

call qg_config_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_incr,flds)

! -- calculate potential vorticity and wind components

call meridional_wind_ad(flds%v,flds%x,flds%nx,flds%ny, conf%deltax)
call zonal_wind_ad(flds%u,flds%x,flds%nx,flds%ny,conf%deltay)
call calc_pv_ad(flds%q,flds%x,flds%nx,flds%ny, &
              & conf%f1,conf%f2,conf%deltax,conf%deltay)
flds%q=0.0_kind_real;
flds%u=0.0_kind_real;
flds%v=0.0_kind_real;

return
end subroutine c_qg_prepare_integration_ad

! ------------------------------------------------------------------------------
!> Prepare for an integration of the QG model - Tangent Linear

subroutine c_qg_prepare_integration_tl(c_key_conf, c_key_incr) &
           bind(c,name='qg_prepare_integration_tl_f90')

use iso_c_binding
use qg_fields
use qg_configs

implicit none
integer(c_int), intent(in) :: c_key_conf  !< Configuration structure
integer(c_int), intent(in) :: c_key_incr  !< Model fields

type(qg_config), pointer :: conf
type(qg_field), pointer  :: flds

call qg_config_registry%get(c_key_conf, conf)
call qg_field_registry%get(c_key_incr, flds)

! -- calculate potential vorticity and wind components

call calc_pv_tl(flds%q,flds%x,flds%nx,flds%ny, &
              & conf%f1,conf%f2,conf%deltax,conf%deltay)
call zonal_wind_tl(flds%u,flds%x,flds%nx,flds%ny,conf%deltay)
call meridional_wind_tl(flds%v,flds%x,flds%nx,flds%ny, conf%deltax)

end subroutine c_qg_prepare_integration_tl

! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model

subroutine c_qg_propagate(c_key_conf, c_key_state) bind(c,name='qg_propagate_f90')

use iso_c_binding
use qg_fields
use qg_configs

implicit none
integer(c_int), intent(in) :: c_key_conf  !< Config structure
integer(c_int), intent(in) :: c_key_state !< Model fields

type(qg_config), pointer :: conf
type(qg_field),  pointer :: flds

call qg_config_registry%get(c_key_conf, conf)
call qg_field_registry%get(c_key_state,flds)

call propagate(flds, conf)

return
end subroutine c_qg_propagate

! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model - Adjoint
!> This routine is called from C++ to propagate the adjoint variables

subroutine c_qg_propagate_ad(c_key_conf, c_key_incr, c_key_traj) &
           bind(c,name='qg_propagate_ad_f90')

use iso_c_binding
use qg_fields
use qg_trajectories
use qg_configs
use qg_constants, only: bet
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf !< Config structure
integer(c_int), intent(in) :: c_key_incr !< Model fields
integer(c_int), intent(in) :: c_key_traj !< Trajectory structure

type(qg_config),     pointer :: conf
type(qg_field),      pointer :: flds
type(qg_trajectory), pointer :: traj

real(kind=kind_real), allocatable :: qnew(:,:,:), x_traj(:,:,:)
real(kind=kind_real), allocatable :: q_traj(:,:,:), u_traj(:,:,:), v_traj(:,:,:)
real(kind=kind_real), allocatable :: qn_traj(:,:), qs_traj(:,:)
real(kind=kind_real) :: xn_traj(2), xs_traj(2)

call qg_config_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_incr,flds)
call qg_traj_registry%get(c_key_traj,traj)

!--- workspace and trajectory
allocate(qnew(flds%nx,flds%ny,2))
allocate(x_traj(flds%nx,flds%ny,2))
allocate(qn_traj(flds%nx,2))
allocate(qs_traj(flds%nx,2))
allocate(q_traj(flds%nx,flds%ny,2))
allocate(u_traj(flds%nx,flds%ny,2))
allocate(v_traj(flds%nx,flds%ny,2))
call get_traj(traj,flds%nx,flds%ny,x_traj,xn_traj,xs_traj,qn_traj,qs_traj)

!--- generate trajectory values for potential vorticity and wind

call calc_pv(flds%nx,flds%ny,q_traj,x_traj,xn_traj,xs_traj, &
    &        conf%f1,conf%f2,conf%deltax,conf%deltay,bet,conf%rs)
call zonal_wind (u_traj,x_traj,xn_traj,xs_traj,flds%nx,flds%ny, conf%deltay)
call meridional_wind (v_traj,x_traj,flds%nx,flds%ny,conf%deltax)

! -- calculate potential vorticity and wind components

call meridional_wind_ad(flds%v,flds%x,flds%nx,flds%ny,conf%deltax)
call zonal_wind_ad(flds%u,flds%x,flds%nx,flds%ny,conf%deltay)
qnew(:,:,:) = flds%q(:,:,:)

!--- invert the potential vorticity to determine streamfunction

call invert_pv_ad(flds%x,qnew,flds%nx,flds%ny, &
                & conf%deltax,conf%deltay,conf%f1,conf%f2)

!--- advect the potential vorticity

call advect_pv_ad(qnew,flds%q,q_traj,qn_traj,qs_traj, &
    &             flds%u,u_traj,flds%v,v_traj,flds%nx,flds%ny, &
    &             conf%deltax,conf%deltay,conf%dt)

!--- clean-up
deallocate(qnew)
deallocate(x_traj)
deallocate(qn_traj)
deallocate(qs_traj)
deallocate(q_traj)
deallocate(u_traj)
deallocate(v_traj)

return
end subroutine c_qg_propagate_ad

! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model - Tangent Linear
!> This routine is called from C++ to propagate the increment

subroutine c_qg_propagate_tl(c_key_conf, c_key_incr, c_key_traj) &
           bind(c,name='qg_propagate_tl_f90')

use iso_c_binding
use qg_fields
use qg_trajectories
use qg_configs
use qg_constants, only: bet
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf !< Config structure
integer(c_int), intent(in) :: c_key_incr !< Model fields
integer(c_int), intent(in) :: c_key_traj !< Trajectory structure

type(qg_config),     pointer :: conf
type(qg_field),      pointer :: flds
type(qg_trajectory), pointer :: traj

real(kind_real), allocatable :: qnew(:,:,:), x_traj(:,:,:)
real(kind_real), allocatable :: q_traj(:,:,:), u_traj(:,:,:), v_traj(:,:,:)
real(kind_real), allocatable :: qn_traj(:,:), qs_traj(:,:)
real(kind_real) :: xn_traj(2), xs_traj(2)

call qg_config_registry%get(c_key_conf, conf)
call qg_field_registry%get(c_key_incr,flds)
call qg_traj_registry%get(c_key_traj,traj)

!--- workspace and trajectory
allocate(qnew(flds%nx,flds%ny,2))
allocate(x_traj(flds%nx,flds%ny,2))
allocate(qn_traj(flds%nx,2))
allocate(qs_traj(flds%nx,2))
allocate(q_traj(flds%nx,flds%ny,2))
allocate(u_traj(flds%nx,flds%ny,2))
allocate(v_traj(flds%nx,flds%ny,2))
call get_traj(traj,flds%nx,flds%ny,x_traj,xn_traj,xs_traj,qn_traj,qs_traj)

!--- generate trajectory values for potential vorticity and wind

call calc_pv(flds%nx,flds%ny,q_traj,x_traj,xn_traj,xs_traj, &
    &        conf%f1,conf%f2,conf%deltax,conf%deltay,bet,conf%rs)
call zonal_wind (u_traj,x_traj,xn_traj,xs_traj,flds%nx,flds%ny, conf%deltay)
call meridional_wind (v_traj,x_traj,flds%nx,flds%ny,conf%deltax)

!--- advect the potential vorticity

qnew(:,:,:)=0.0_kind_real
call advect_pv_tl(qnew,flds%q,q_traj,qn_traj,qs_traj, &
    &             flds%u,u_traj,flds%v,v_traj,flds%nx,flds%ny,&
    &             conf%deltax,conf%deltay,conf%dt)

!--- invert the potential vorticity to determine streamfunction

call invert_pv_tl(flds%x,qnew,flds%nx,flds%ny, &
                & conf%deltax,conf%deltay,conf%f1,conf%f2)

! -- calculate potential vorticity and wind components

flds%q(:,:,:) = qnew(:,:,:)
call zonal_wind_tl(flds%u,flds%x,flds%nx,flds%ny,conf%deltay)
call meridional_wind_tl(flds%v,flds%x,flds%nx,flds%ny,conf%deltax)

!--- clean-up
deallocate(qnew)
deallocate(x_traj)
deallocate(qn_traj)
deallocate(qs_traj)
deallocate(q_traj)
deallocate(u_traj)
deallocate(v_traj)

return
end subroutine c_qg_propagate_tl

! ------------------------------------------------------------------------------
!> Save trajectory and perform a timestep of the QG model

subroutine c_qg_prop_traj(c_key_conf, c_key_state, c_key_traj) bind(c,name='qg_prop_traj_f90')

use iso_c_binding
use qg_fields
use qg_configs
use qg_trajectories

implicit none
integer(c_int), intent(in)    :: c_key_conf  !< Config structure
integer(c_int), intent(in)    :: c_key_state !< Model fields
integer(c_int), intent(inout) :: c_key_traj  !< Trajectory structure

type(qg_config),     pointer :: conf
type(qg_field),      pointer :: flds
type(qg_trajectory), pointer :: traj

call qg_config_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_state,flds)

call qg_traj_registry%init()            
call qg_traj_registry%add(c_key_traj)
call qg_traj_registry%get(c_key_traj,traj)
call set_traj(traj,flds%nx,flds%ny, &
            & flds%x,flds%x_north,flds%x_south,flds%q_north,flds%q_south)

return
end subroutine c_qg_prop_traj

! ------------------------------------------------------------------------------
!> Wipe out trajectory of the QG model

subroutine c_qg_wipe_traj(c_key_traj) bind(c,name='qg_wipe_traj_f90')

use iso_c_binding
use qg_trajectories

implicit none
integer(c_int), intent(inout) :: c_key_traj  !< Trajectory structure
type(qg_trajectory), pointer :: traj

call qg_traj_registry%get(c_key_traj,traj)
call delete_traj(traj)
call qg_traj_registry%remove(c_key_traj)

return
end subroutine c_qg_wipe_traj

! ------------------------------------------------------------------------------
