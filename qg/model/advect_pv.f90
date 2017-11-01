! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Advect potential vorticity

!> Potential vorticity is advected using cubic Lagrange interpolation
!! in the zonal direction, and then the meridional direction.
!!
!! Note that this is a first-order scheme. The upstream point is
!! determined solely from the wind at the arrival point.

subroutine advect_pv (qnew,q,q_north,q_south,u,v,nx,ny,deltax,deltay,dt)

!--- semi-Lagrangian advection of pv

use kinds

implicit none
real(kind_real), intent(out) :: qnew(nx,ny,2) !< Output potential vorticity
real(kind_real), intent(in)  :: q(nx,ny,2)    !< Input potential vorticity
real(kind_real), intent(in)  :: q_north(nx,2) !< PV on northern wall
real(kind_real), intent(in)  :: q_south(nx,2) !< PV on southern wall
real(kind_real), intent(in)  :: u(nx,ny,2)    !< Advecting zonal wind
real(kind_real), intent(in)  :: v(nx,ny,2)    !< Advecting meridional wind
integer, intent(in) :: nx          !< Zonal grid dimensions
integer, intent(in) :: ny          !< Meridional grid dimensions
real(kind_real), intent(in)  :: deltax        !< Zonal grid spacing (non-dimensional)
real(kind_real), intent(in)  :: deltay        !< Meridional grid spacing (non-dimensional)
real(kind_real), intent(in)  :: dt            !< Timestep (non-dimensional)

integer :: ii,jj,kk,ixm1,ix,ixp1,ixp2,jym1,jy,jyp1,jyp2
real(kind_real) :: ax,ay,qjm1,qj,qjp1,qjp2
real(kind_real), parameter :: one=1.0_kind_real
real(kind_real), parameter :: two=2.0_kind_real
real(kind_real), parameter :: half=0.5_kind_real
real(kind_real), parameter :: sixth=1.0_kind_real/6.0_kind_real

!--- advect q, returning qnew

do kk=1,2
  do jj=1,ny
    do ii=1,nx

!--- find the interpolation point

      ax = real(ii,kind_real) - u(ii,jj,kk)*dt/deltax
      ix = floor(ax)
      ax = ax-real(ix,kind_real)
      ixm1 = 1 + modulo(ix-2,nx)
      ixp1 = 1 + modulo(ix  ,nx)
      ixp2 = 1 + modulo(ix+1,nx)
      ix   = 1 + modulo(ix-1,nx)

      ay = real(jj,kind_real) - v(ii,jj,kk)*dt/deltay
      jy = floor(ay)
      ay = ay-real(jy,kind_real)
      jym1 = jy-1
      jyp1 = jy+1
      jyp2 = jy+2

!--- Lagrange interpolation in the zonal direction
      
      if (jym1 < 1) then
        qjm1 =  ax     *(ax-one)*(ax-two)*q_south(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_south(ix,  kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q_south(ixp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q_south(ixp2,kk)*(sixth)
      else if (jym1 > ny) then
        qjm1 =  ax     *(ax-one)*(ax-two)*q_north(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_north(ix,  kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q_north(ixp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q_north(ixp2,kk)*(sixth)
      else
        qjm1 =  ax     *(ax-one)*(ax-two)*q(ixm1,jym1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q(ix,  jym1,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q(ixp1,jym1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q(ixp2,jym1,kk)*(sixth)
      endif

      if (jy < 1) then
        qj   =  ax     *(ax-one)*(ax-two)*q_south(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_south(ix  ,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q_south(ixp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q_south(ixp2,kk)*(sixth)
      else if (jy > ny) then
        qj   =  ax     *(ax-one)*(ax-two)*q_north(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_north(ix  ,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q_north(ixp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q_north(ixp2,kk)*(sixth)
      else
        qj   =  ax     *(ax-one)*(ax-two)*q(ixm1,jy,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q(ix  ,jy,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q(ixp1,jy,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q(ixp2,jy,kk)*(sixth)
      endif

      if (jyp1 < 1) then
        qjp1 =  ax     *(ax-one)*(ax-two)*q_south(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_south(ix  ,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q_south(ixp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q_south(ixp2,kk)*(sixth)
      else if (jyp1 > ny) then
        qjp1 =  ax     *(ax-one)*(ax-two)*q_north(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_north(ix  ,kk)*half    + &
             & (ax+one)* ax     *(ax-two)*q_north(ixp1,kk)*(-half) + &
             & (ax+one)* ax     *(ax-one)*q_north(ixp2,kk)*(sixth)
      else
        qjp1 =  ax     *(ax-one)*(ax-two)*q(ixm1,jyp1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q(ix  ,jyp1,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q(ixp1,jyp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q(ixp2,jyp1,kk)*(sixth)
      endif

      if (jyp2 < 1) then
        qjp2 =  ax     *(ax-one)*(ax-two)*q_south(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_south(ix  ,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q_south(ixp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q_south(ixp2,kk)*(sixth)
      else if (jyp2 > ny) then
        qjp2 =  ax     *(ax-one)*(ax-two)*q_north(ixm1,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q_north(ix  ,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q_north(ixp1,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q_north(ixp2,kk)*(sixth)
      else
        qjp2 =  ax     *(ax-one)*(ax-two)*q(ixm1,jyp2,kk)*(-sixth) + &
             & (ax+one)*(ax-one)*(ax-two)*q(ix  ,jyp2,kk)*half      + &
             & (ax+one)* ax     *(ax-two)*q(ixp1,jyp2,kk)*(-half)   + &
             & (ax+one)* ax     *(ax-one)*q(ixp2,jyp2,kk)*(sixth)
      endif

!--- Lagrange interpolation in the meridional direction

      qnew(ii,jj,kk) =  ay     *(ay-one)*(ay-two)*(-sixth)*qjm1 + &
                     & (ay+one)*(ay-one)*(ay-two)*half     *qj   + &
                     & (ay+one)* ay     *(ay-two)*(-half)  *qjp1 + &
                     & (ay+one)* ay     *(ay-one)*(sixth) *qjp2
    enddo
  enddo
enddo

end subroutine advect_pv
