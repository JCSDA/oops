! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_interp_mod

use fckit_log_module,only: fckit_log
use kinds
use qg_constants_mod
use qg_geom_mod
use qg_projection_mod
use qg_tools_mod

implicit none

private
public :: qg_interp_trilinear,qg_interp_trilinear_ad, &
        & qg_interp_bicubic,qg_interp_bicubic_tl,qg_interp_bicubic_ad, &
        & find_x_indices,find_y_indices
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Trilinear interpolation
subroutine qg_interp_trilinear(geom,lon,lat,z,gfld3d,val)

! Passed variables
type(qg_geom),intent(in) :: geom                              !< Geometry
real(kind_real),intent(in) :: lon                             !< Longitude
real(kind_real),intent(in) :: lat                             !< Latitude
real(kind_real),intent(in) :: z                               !< Altitude
real(kind_real),intent(in) :: gfld3d(geom%nx,geom%ny,geom%nz) !< 3D Field
real(kind_real),intent(out) :: val                            !< Value

! Local variables
integer :: jxm1,jxo,jxp1,jxp2
integer :: jym1,jyo,jyp1,jyp2
integer :: jzo,jzp1
real(kind_real) :: x,y
real(kind_real) :: ax,ay,az
real(kind_real) :: oo,op1,p1o,p1p1,o,p1

! Convert lon/lat to x/y
call lonlat_to_xy(lon,lat,x,y)

! Find indices
call find_x_indices(geom,x,jxm1,jxo,jxp1,jxp2,ax)
call find_y_indices(geom,y,jym1,jyo,jyp1,jyp2,ay)
call find_z_indices(geom,z,jzo,jzp1,az)

! Extrapolate along y if needed
if (jyo<1) then
  ay = ay+real(jyo-1,kind_real)
  jyo = 1
  jyp1 = 2
endif
if (jyp1>geom%ny) then
  ay = ay+real(jyp1-geom%ny,kind_real)
  jyo = geom%ny-1
  jyp1 = geom%ny
endif

! Interpolate along x
oo = (1.0-ax)*gfld3d(jxo,jyo,jzo)+ax*gfld3d(jxp1,jyo,jzo)
op1 = (1.0-ax)*gfld3d(jxo,jyo,jzp1)+ax*gfld3d(jxp1,jyo,jzp1)
p1o = (1.0-ax)*gfld3d(jxo,jyp1,jzo)+ax*gfld3d(jxp1,jyp1,jzo)
p1p1 = (1.0-ax)*gfld3d(jxo,jyp1,jzp1)+ax*gfld3d(jxp1,jyp1,jzp1)

! Interpolate along y
o = (1.0-ay)*oo+ay*p1o
p1 = (1.0-ay)*op1+ay*p1p1

! Interpolate along z
val = (1.0-az)*o+az*p1

end subroutine qg_interp_trilinear
! ------------------------------------------------------------------------------
!> Interpolation - adjoint
subroutine qg_interp_trilinear_ad(geom,lon,lat,z,val,gfld3d)

! Passed variables
type(qg_geom),intent(in) :: geom                                 !< Geometry
real(kind_real),intent(in) :: lon                                !< Longitude
real(kind_real),intent(in) :: lat                                !< Latitude
real(kind_real),intent(in) :: z                                  !< Altitude
real(kind_real),intent(in) :: val                                !< Value
real(kind_real),intent(inout) :: gfld3d(geom%nx,geom%ny,geom%nz) !< 3D field

! Local variables
integer :: jxm1,jxo,jxp1,jxp2
integer :: jym1,jyo,jyp1,jyp2
integer :: jzo,jzp1
real(kind_real) :: x,y
real(kind_real) :: ax,ay,az
real(kind_real) :: oo,op1,p1o,p1p1,o,p1

! Convert lon/lat to x/y
call lonlat_to_xy(lon,lat,x,y)

! Find indices
call find_x_indices(geom,x,jxm1,jxo,jxp1,jxp2,ax)
call find_y_indices(geom,y,jym1,jyo,jyp1,jyp2,ay)
call find_z_indices(geom,z,jzo,jzp1,az)

! Extrapolate along y if needed
if (jyo<1) then
  ay = ay+real(jyo-1,kind_real)
  jyo = 1
  jyp1 = 2
endif
if (jyp1>geom%ny) then
  ay = ay+real(jyp1-geom%ny,kind_real)
  jyo = geom%ny-1
  jyp1 = geom%ny
endif

! Interpolate along z
o = (1.0-az)*val
p1 = az*val

! Interpolate along y
oo = (1.0-ay)*o
p1o = ay*o
op1 = (1.0-ay)*p1
p1p1 = ay*p1

! Interpolate along x
gfld3d(jxo,jyo,jzo) = gfld3d(jxo,jyo,jzo)+(1.0-ax)*oo
gfld3d(jxp1,jyo,jzo) = gfld3d(jxp1,jyo,jzo)+ax*oo
gfld3d(jxo,jyo,jzp1) = gfld3d(jxo,jyo,jzp1)+(1.0-ax)*op1
gfld3d(jxp1,jyo,jzp1) = gfld3d(jxp1,jyo,jzp1)+ax*op1
gfld3d(jxo,jyp1,jzo) = gfld3d(jxo,jyp1,jzo)+(1.0-ax)*p1o
gfld3d(jxp1,jyp1,jzo) = gfld3d(jxp1,jyp1,jzo)+ax*p1o
gfld3d(jxo,jyp1,jzp1) = gfld3d(jxo,jyp1,jzp1)+(1.0-ax)*p1p1
gfld3d(jxp1,jyp1,jzp1) = gfld3d(jxp1,jyp1,jzp1)+ax*p1p1

end subroutine qg_interp_trilinear_ad
! ------------------------------------------------------------------------------
!> Bicubic horizontal interpolation
subroutine qg_interp_bicubic(geom,x,y,gfld2dext,val)

! Passed variables
type(qg_geom),intent(in) :: geom                              !< Geometry
real(kind_real),intent(in) :: x                               !< X value
real(kind_real),intent(in) :: y                               !< Y value
real(kind_real),intent(in) :: gfld2dext(geom%nx,-1:geom%ny+2) !< Extended 2D field
real(kind_real),intent(out) :: val                            !< Value

! Local variables
integer :: jxm1,jxo,jxp1,jxp2
integer :: jym1,jyo,jyp1,jyp2
real(kind_real) :: ax,ay,m1,o,p1,p2

! Find indices
call find_x_indices(geom,x,jxm1,jxo,jxp1,jxp2,ax)
call find_y_indices(geom,y,jym1,jyo,jyp1,jyp2,ay)

! Interpolation along x
call cubic(ax,gfld2dext(jxm1,jym1),gfld2dext(jxo,jym1),gfld2dext(jxp1,jym1),gfld2dext(jxp2,jym1),m1)
call cubic(ax,gfld2dext(jxm1,jyo),gfld2dext(jxo,jyo),gfld2dext(jxp1,jyo),gfld2dext(jxp2,jyo),o)
call cubic(ax,gfld2dext(jxm1,jyp1),gfld2dext(jxo,jyp1),gfld2dext(jxp1,jyp1),gfld2dext(jxp2,jyp1),p1)
call cubic(ax,gfld2dext(jxm1,jyp2),gfld2dext(jxo,jyp2),gfld2dext(jxp1,jyp2),gfld2dext(jxp2,jyp2),p2)

! Interpolation along y
call cubic(ay,m1,o,p1,p2,val)

end subroutine qg_interp_bicubic
! ------------------------------------------------------------------------------
!> Bicubic horizontal interpolation - tangent linear
subroutine qg_interp_bicubic_tl(geom,x,y,gfld2dext_traj,dx,dy,gfld2dext,val)

! Passed variables
type(qg_geom),intent(in) :: geom                                   !< Geometry
real(kind_real),intent(in) :: x                                    !< X value
real(kind_real),intent(in) :: y                                    !< Y value
real(kind_real),intent(in) :: gfld2dext_traj(geom%nx,-1:geom%ny+2) !< Extended 2D trajectory
real(kind_real),intent(in) :: dx                                   !< X perturbation
real(kind_real),intent(in) :: dy                                   !< Y perturbation
real(kind_real),intent(in) :: gfld2dext(geom%nx,-1:geom%ny+2)      !< Extended 2D perturbation
real(kind_real),intent(out) :: val                                 !< Value

! Local variables
integer :: jxm1,jxo,jxp1,jxp2
integer :: jym1,jyo,jyp1,jyp2
real(kind_real) :: ax_traj,ay_traj,m1_traj,o_traj,p1_traj,p2_traj,ax,ay,m1,o,p1,p2

! Find indices
call find_x_indices(geom,x,jxm1,jxo,jxp1,jxp2,ax_traj)
call find_y_indices(geom,y,jym1,jyo,jyp1,jyp2,ay_traj)

! Interpolation along x (trajectory)
call cubic(ax_traj,gfld2dext_traj(jxm1,jym1),gfld2dext_traj(jxo,jym1),gfld2dext_traj(jxp1,jym1),gfld2dext_traj(jxp2,jym1),m1_traj)
call cubic(ax_traj,gfld2dext_traj(jxm1,jyo),gfld2dext_traj(jxo,jyo),gfld2dext_traj(jxp1,jyo),gfld2dext_traj(jxp2,jyo),o_traj)
call cubic(ax_traj,gfld2dext_traj(jxm1,jyp1),gfld2dext_traj(jxo,jyp1),gfld2dext_traj(jxp1,jyp1),gfld2dext_traj(jxp2,jyp1),p1_traj)
call cubic(ax_traj,gfld2dext_traj(jxm1,jyp2),gfld2dext_traj(jxo,jyp2),gfld2dext_traj(jxp1,jyp2),gfld2dext_traj(jxp2,jyp2),p2_traj)

! Compute linearization coefficients
ax = dx/geom%deltax
ay = dy/geom%deltay

! Interpolation along x (perturbation)
call cubic_tl(ax_traj,gfld2dext_traj(jxm1,jym1),gfld2dext_traj(jxo,jym1),gfld2dext_traj(jxp1,jym1),gfld2dext_traj(jxp2,jym1), &
            & ax,gfld2dext(jxm1,jym1),gfld2dext(jxo,jym1),gfld2dext(jxp1,jym1),gfld2dext(jxp2,jym1),m1)
call cubic_tl(ax_traj,gfld2dext_traj(jxm1,jyo),gfld2dext_traj(jxo,jyo),gfld2dext_traj(jxp1,jyo),gfld2dext_traj(jxp2,jyo), &
            & ax,gfld2dext(jxm1,jyo),gfld2dext(jxo,jyo),gfld2dext(jxp1,jyo),gfld2dext(jxp2,jyo),o)
call cubic_tl(ax_traj,gfld2dext_traj(jxm1,jyp1),gfld2dext_traj(jxo,jyp1),gfld2dext_traj(jxp1,jyp1),gfld2dext_traj(jxp2,jyp1), &
            & ax,gfld2dext(jxm1,jyp1),gfld2dext(jxo,jyp1),gfld2dext(jxp1,jyp1),gfld2dext(jxp2,jyp1),p1)
call cubic_tl(ax_traj,gfld2dext_traj(jxm1,jyp2),gfld2dext_traj(jxo,jyp2),gfld2dext_traj(jxp1,jyp2),gfld2dext_traj(jxp2,jyp2), &
            & ax,gfld2dext(jxm1,jyp2),gfld2dext(jxo,jyp2),gfld2dext(jxp1,jyp2),gfld2dext(jxp2,jyp2),p2)

! Interpolation along y (perturbation)
call cubic_tl(ay_traj,m1_traj,o_traj,p1_traj,p2_traj,ay,m1,o,p1,p2,val)

end subroutine qg_interp_bicubic_tl
! ------------------------------------------------------------------------------
!> Bicubic horizontal interpolation - adjoint
subroutine qg_interp_bicubic_ad(geom,x,y,gfld2dext_traj,val,dx,dy,gfld2dext)

! Passed variables
type(qg_geom),intent(in) :: geom                                   !< Geometry
real(kind_real),intent(in) :: x                                    !< X value
real(kind_real),intent(in) :: y                                    !< Y value
real(kind_real),intent(in) :: gfld2dext_traj(geom%nx,-1:geom%ny+2) !< Extended 2D trajectory
real(kind_real),intent(in) :: val                                  !< Value
real(kind_real),intent(inout) :: dx                                !< X perturbation
real(kind_real),intent(inout) :: dy                                !< Y perturbation
real(kind_real),intent(inout) :: gfld2dext(geom%nx,-1:geom%ny+2)   !< Extended 2D perturbation

! Local variables
integer :: jxm1,jxo,jxp1,jxp2
integer :: jym1,jyo,jyp1,jyp2
real(kind_real) :: ax_traj,ay_traj,m1_traj,o_traj,p1_traj,p2_traj,ax,ay,m1,o,p1,p2

! Find indices
call find_x_indices(geom,x,jxm1,jxo,jxp1,jxp2,ax_traj)
call find_y_indices(geom,y,jym1,jyo,jyp1,jyp2,ay_traj)

! Interpolation along x (trajectory)
call cubic(ax_traj,gfld2dext_traj(jxm1,jym1),gfld2dext_traj(jxo,jym1),gfld2dext_traj(jxp1,jym1),gfld2dext_traj(jxp2,jym1),m1_traj)
call cubic(ax_traj,gfld2dext_traj(jxm1,jyo),gfld2dext_traj(jxo,jyo),gfld2dext_traj(jxp1,jyo),gfld2dext_traj(jxp2,jyo),o_traj)
call cubic(ax_traj,gfld2dext_traj(jxm1,jyp1),gfld2dext_traj(jxo,jyp1),gfld2dext_traj(jxp1,jyp1),gfld2dext_traj(jxp2,jyp1),p1_traj)
call cubic(ax_traj,gfld2dext_traj(jxm1,jyp2),gfld2dext_traj(jxo,jyp2),gfld2dext_traj(jxp1,jyp2),gfld2dext_traj(jxp2,jyp2),p2_traj)

! Initialization
m1 = 0.0
o = 0.0
p1 = 0.0
p2 = 0.0
ax = 0.0
ay = 0.0

! Interpolation along y (perturbation)
call cubic_ad(ay_traj,m1_traj,o_traj,p1_traj,p2_traj,val,ay,m1,o,p1,p2)

! Interpolation along x (perturbation)
call cubic_ad(ax_traj,gfld2dext_traj(jxm1,jym1),gfld2dext_traj(jxo,jym1),gfld2dext_traj(jxp1,jym1),gfld2dext_traj(jxp2,jym1), &
            & m1,ax,gfld2dext(jxm1,jym1),gfld2dext(jxo,jym1),gfld2dext(jxp1,jym1),gfld2dext(jxp2,jym1))
call cubic_ad(ax_traj,gfld2dext_traj(jxm1,jyo),gfld2dext_traj(jxo,jyo),gfld2dext_traj(jxp1,jyo),gfld2dext_traj(jxp2,jyo), &
            & o,ax,gfld2dext(jxm1,jyo),gfld2dext(jxo,jyo),gfld2dext(jxp1,jyo),gfld2dext(jxp2,jyo))
call cubic_ad(ax_traj,gfld2dext_traj(jxm1,jyp1),gfld2dext_traj(jxo,jyp1),gfld2dext_traj(jxp1,jyp1),gfld2dext_traj(jxp2,jyp1), &
            & p1,ax,gfld2dext(jxm1,jyp1),gfld2dext(jxo,jyp1),gfld2dext(jxp1,jyp1),gfld2dext(jxp2,jyp1))
call cubic_ad(ax_traj,gfld2dext_traj(jxm1,jyp2),gfld2dext_traj(jxo,jyp2),gfld2dext_traj(jxp1,jyp2),gfld2dext_traj(jxp2,jyp2), &
            & p2,ax,gfld2dext(jxm1,jyp2),gfld2dext(jxo,jyp2),gfld2dext(jxp1,jyp2),gfld2dext(jxp2,jyp2))

! Compute linearization coefficients
dx = dx+ax/geom%deltax
dy = dy+ay/geom%deltay

end subroutine qg_interp_bicubic_ad
! ------------------------------------------------------------------------------
!> Find x indices 
subroutine find_x_indices(geom,x,jxm1,jxo,jxp1,jxp2,ax)

! Passed variables
type(qg_geom),intent(in) :: geom  !< Geometry
real(kind_real),intent(in) :: x   !< X value
integer,intent(out) :: jxm1       !< Minus 1 index
integer,intent(out) :: jxo        !< Origin index
integer,intent(out) :: jxp1       !< Plus 1 index
integer,intent(out) :: jxp2       !< Plus 2 index
real(kind_real),intent(out) :: ax !< Coefficient

! Local variables
integer :: ix
real(kind_real) :: xx

! Initialization
if (x<0.0) then
  xx = x+domain_zonal
elseif (x>domain_zonal) then
  xx = x-domain_zonal
else
  xx = x
endif

! Find base index and weight
if ((0.0<=xx).and.(xx<geom%x(1))) then
  jxo = geom%nx
  ax = (domain_zonal-geom%x(jxo)+xx)/geom%deltax
elseif ((geom%x(geom%nx)<=xx).and.(xx<=domain_zonal)) then
  jxo = geom%nx
  ax = (xx-geom%x(jxo))/geom%deltax
else
  jxo = -1
  do ix=1,geom%nx-1
    if ((geom%x(ix)<=xx).and.(xx<geom%x(ix+1))) then
      jxo = ix
      exit
    end if
  end do
  if (jxo==-1) call abor1_ftn('find_x_indices: cannot find jxo')
  ax = (xx-geom%x(jxo))/geom%deltax
end if

! Define other indices
jxm1 = jxo-1
if (jxm1<1) jxm1 = jxm1+geom%nx
jxp1 = jxo+1
if (jxp1>geom%nx) jxp1 = jxp1-geom%nx
jxp2 = jxo+2
if (jxp2>geom%nx) jxp2 = jxp2-geom%nx

end subroutine find_x_indices
! ------------------------------------------------------------------------------
!> Find y indices 
subroutine find_y_indices(geom,y,jym1,jyo,jyp1,jyp2,ay)

! Passed variables
type(qg_geom),intent(in) :: geom  !< Geometry
real(kind_real),intent(in) :: y   !< Y value
integer,intent(out) :: jym1       !< Minus 1 index
integer,intent(out) :: jyo        !< Origin index
integer,intent(out) :: jyp1       !< Plus 1 index
integer,intent(out) :: jyp2       !< Plus 2 index
real(kind_real),intent(out) :: ay !< Coefficient

! Local variables
integer :: iy
real(kind_real) :: yy

! Initialization
yy = max(0.0_kind_real,min(y,domain_meridional))

! Find base index and weight
if ((0.0<=yy).and.(yy<geom%y(1))) then
  jyo = 0
  ay = yy/geom%deltay
elseif ((geom%y(geom%ny)<=yy).and.(yy<=domain_meridional)) then
  jyo = geom%ny
  ay = (yy-geom%y(jyo))/geom%deltay
else
  jyo = -1
  do iy=1,geom%ny-1
    if ((geom%y(iy)<=yy).and.(yy<geom%y(iy+1))) then
      jyo = iy
      exit
    end if
  end do
  if (jyo==-1) call abor1_ftn('find_y_indices: cannot find jyo')
  ay = (yy-geom%y(jyo))/geom%deltay
end if

! Define other indices
jym1 = jyo-1
jyp1 = jyo+1
jyp2 = jyo+2

end subroutine find_y_indices
! ------------------------------------------------------------------------------
!> Find z indices 
subroutine find_z_indices(geom,z,jzo,jzp1,az)

! Passed variables
type(qg_geom),intent(in) :: geom  !< Geometry
real(kind_real),intent(in) :: z   !< Z value
integer,intent(out) :: jzo        !< Origin index
integer,intent(out) :: jzp1       !< Plus 1 index
real(kind_real),intent(out) :: az !< Coefficient

! Local variables
integer :: iz
real(kind_real) :: zz

! Initialization
zz = max(0.0_kind_real,min(z,domain_depth))

! Find base index and weight
if ((0.0<=zz).and.(zz<geom%z(1))) then
  jzo = 1
elseif ((geom%z(geom%nz)<=zz).and.(zz<=domain_depth)) then
  jzo = geom%nz-1
else
  jzo = -1
  do iz=1,geom%nz-1
    if ((geom%z(iz)<=zz).and.(zz<geom%z(iz+1))) then
      jzo = iz
      exit
    end if
  end do
  if (jzo==-1) call abor1_ftn('find_z_indices: cannot find jzo')  
end if
az = (zz-geom%z(jzo))/(geom%z(jzo+1)-geom%z(jzo))

! Define other indices
jzp1 = jzo+1

end subroutine find_z_indices
! ------------------------------------------------------------------------------
!> Cubic interpolation
subroutine cubic(a,m1,o,p1,p2,res)

! Passed variables
real(kind_real),intent(in) :: a    !< Coefficient
real(kind_real),intent(in) :: m1   !< Minus 1 value
real(kind_real),intent(in) :: o    !< Origin value
real(kind_real),intent(in) :: p1   !< Plus 1 value
real(kind_real),intent(in) :: p2   !< Plus 2 value
real(kind_real),intent(out) :: res !< Result

res = -a*(a-1.0)*(a-2.0)/6.0*m1 &
    & +(a+1.0)*(a-1.0)*(a-2.0)/2.0*o &
    & -(a+1.0)*a*(a-2.0)/2.0*p1 &
    & +(a+1.0)*a*(a-1.0)/6.0*p2

end subroutine cubic
! ------------------------------------------------------------------------------
!> Cubic interpolation - tangent linear
subroutine cubic_tl(a_traj,m1_traj,o_traj,p1_traj,p2_traj,a,m1,o,p1,p2,res)

! Passed variables
real(kind_real),intent(in) :: a_traj  !< Coefficient trajectory
real(kind_real),intent(in) :: m1_traj !< Minus 1 trajectory
real(kind_real),intent(in) :: o_traj  !< Origin value
real(kind_real),intent(in) :: p1_traj !< Plus 1 value
real(kind_real),intent(in) :: p2_traj !< Plus 2 value
real(kind_real),intent(in) :: a       !< Coefficient perturbation
real(kind_real),intent(in) :: m1      !< Minus perturbation
real(kind_real),intent(in) :: o       !< Origin perturbation
real(kind_real),intent(in) :: p1      !< Plus 1 perturbation
real(kind_real),intent(in) :: p2      !< Plus 2 perturbation
real(kind_real),intent(out) :: res    !< Result

res = -a*(a_traj-1.0)*(a_traj-2.0)/6.0*m1_traj &
    & -a_traj*a*(a_traj-2.0)/6.0*m1_traj &
    & -a_traj*(a_traj-1.0)*a/6.0*m1_traj &
    & -a_traj*(a_traj-1.0)*(a_traj-2.0)/6.0*m1 &
    & +a*(a_traj-1.0)*(a_traj-2.0)/2.0*o_traj &
    & +(a_traj+1.0)*a*(a_traj-2.0)/2.0*o_traj &
    & +(a_traj+1.0)*(a_traj-1.0)*a/2.0*o_traj &
    & +(a_traj+1.0)*(a_traj-1.0)*(a_traj-2.0)/2.0*o &
    & -a*a_traj*(a_traj-2.0)/2.0*p1_traj &
    & -(a_traj+1.0)*a*(a_traj-2.0)/2.0*p1_traj &
    & -(a_traj+1.0)*a_traj*a/2.0*p1_traj &
    & -(a_traj+1.0)*a_traj*(a_traj-2.0)/2.0*p1 &
    & +a*a_traj*(a_traj-1.0)/6.0*p2_traj &
    & +(a_traj+1.0)*a*(a_traj-1.0)/6.0*p2_traj &
    & +(a_traj+1.0)*a_traj*a/6.0*p2_traj &
    & +(a_traj+1.0)*a_traj*(a_traj-1.0)/6.0*p2

end subroutine cubic_tl
! ------------------------------------------------------------------------------
!> Cubic interpolation - adjoint
subroutine cubic_ad(a_traj,m1_traj,o_traj,p1_traj,p2_traj,res,a,m1,o,p1,p2)

! Passed variables
real(kind_real),intent(in) :: a_traj  !< Coefficient trajectory
real(kind_real),intent(in) :: m1_traj !< Minus 1 trajectory
real(kind_real),intent(in) :: o_traj  !< Origin value
real(kind_real),intent(in) :: p1_traj !< Plus 1 value
real(kind_real),intent(in) :: p2_traj !< Plus 2 value
real(kind_real),intent(in) :: res     !< Result
real(kind_real),intent(inout) :: a    !< Coefficient perturbation
real(kind_real),intent(inout) :: m1   !< Minus perturbation
real(kind_real),intent(inout) :: o    !< Origin perturbation
real(kind_real),intent(inout) :: p1   !< Plus 1 perturbation
real(kind_real),intent(inout) :: p2   !< Plus 2 perturbation

a = a+(-(a_traj-1.0)*(a_traj-2.0)/6.0*m1_traj &
  & -a_traj*(a_traj-2.0)/6.0*m1_traj &
  & -a_traj*(a_traj-1.0)/6.0*m1_traj &
  & +(a_traj-1.0)*(a_traj-2.0)/2.0*o_traj &
  & +(a_traj+1.0)*(a_traj-2.0)/2.0*o_traj &
  & +(a_traj+1.0)*(a_traj-1.0)/2.0*o_traj &
  & -a_traj*(a_traj-2.0)/2.0*p1_traj &
  & -(a_traj+1.0)*(a_traj-2.0)/2.0*p1_traj &
  & -(a_traj+1.0)*a_traj/2.0*p1_traj &
  & +a_traj*(a_traj-1.0)/6.0*p2_traj &
  & +(a_traj+1.0)*(a_traj-1.0)/6.0*p2_traj &
  & +(a_traj+1.0)*a_traj/6.0*p2_traj) &
  & *res
m1 = m1-a_traj*(a_traj-1.0)*(a_traj-2.0)/6.0*res
o = o+(a_traj+1.0)*(a_traj-1.0)*(a_traj-2.0)/2.0*res
p1 = p1-(a_traj+1.0)*a_traj*(a_traj-2.0)/2.0*res
p2 = p2+(a_traj+1.0)*a_traj*(a_traj-1.0)/6.0*res

end subroutine cubic_ad
! ------------------------------------------------------------------------------
end module qg_interp_mod
