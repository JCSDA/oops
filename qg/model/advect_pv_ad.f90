! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Advect potential vorticity - Adjoint

subroutine advect_pv_ad (qnew,q,q_traj,q_north,q_south, &
          &              u,u_traj,v,v_traj,             &
          &              nx,ny,deltax,deltay,dt)

use kinds
implicit none

real(kind=kind_real), intent(inout) :: qnew(nx,ny,2)   !< Input  adjoint variable for PV
real(kind=kind_real), intent(inout) :: q(nx,ny,2)      !< Output adjoint variable for PV
real(kind=kind_real), intent(in)    :: q_traj(nx,ny,2) !< Trajectory potential vorticity
real(kind=kind_real), intent(in)    :: q_north(nx,2)   !< PV on northern wall
real(kind=kind_real), intent(in)    :: q_south(nx,2)   !< PV on southern wall
real(kind=kind_real), intent(out)   :: u(nx,ny,2)      !< Zonal wind adjoint variable
real(kind=kind_real), intent(in)    :: u_traj(nx,ny,2) !< Trajectory zonal wind
real(kind=kind_real), intent(out)   :: v(nx,ny,2)      !< Meridional wind adjoint variable
real(kind=kind_real), intent(in)    :: v_traj(nx,ny,2) !< Trajectory meridional wind
integer, intent(in) :: nx              !< Zonal grid dimension
integer, intent(in) :: ny              !< Meridional grid dimension
real(kind=kind_real), intent(in)    :: deltax          !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in)    :: deltay          !< Meridional grid spacing (non-dim)
real(kind=kind_real), intent(in)    :: dt              !< Timestep (non-dimensional)

integer :: ii,jj,kk,ixm1,ix,ixp1,ixp2,jym1,jy,jyp1,jyp2
real(kind=kind_real) :: ax_traj,ay_traj,qjm1_traj,qj_traj,qjp1_traj,qjp2_traj
real(kind=kind_real) :: ax,ay,qjm1,qj,qjp1,qjp2
real(kind_real), parameter :: one=1.0_kind_real
real(kind_real), parameter :: two=2.0_kind_real
real(kind_real), parameter :: half=0.5_kind_real
real(kind_real), parameter :: sixth=1.0_kind_real/6.0_kind_real

!--- advect

do kk=1,2
  do jj=1,ny
    do ii=1,nx
      q(ii,jj,kk) = 0.0_kind_real
    enddo
  enddo
enddo

do kk=1,2
  do jj=1,ny
    do ii=1,nx

!==== Trajectory calculations

!--- find the interpolation point (nonlinear)

      ax_traj = real(ii,kind_real) - u_traj(ii,jj,kk)*dt/deltax
      ix = floor(ax_traj)
      ax_traj = ax_traj-real(ix,kind_real)
      ixm1 = 1 + modulo(ix-2,nx)
      ixp1 = 1 + modulo(ix  ,nx)
      ixp2 = 1 + modulo(ix+1,nx)
      ix   = 1 + modulo(ix-1,nx)

      ay_traj = real(jj,kind_real) - v_traj(ii,jj,kk)*dt/deltay
      jy = floor(ay_traj)
      ay_traj = ay_traj-real(jy,kind_real)
      jym1 = jy-1
      jyp1 = jy+1
      jyp2 = jy+2

      if (jym1 < 1) then
        qjm1_traj =  &
            &  ax_traj      *(ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix,  kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            & (ax_traj+one) * ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)
      else if (jym1 > ny) then
        qjm1_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix,  kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)
      else
        qjm1_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jym1,kk)*(-sixth) +                             &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix,  jym1,kk)*half      +                             &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jym1,kk)*(-half)   +                             &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_traj(ixp2,jym1,kk)*(sixth)
      endif

      if (jy < 1) then
        qj_traj   =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)
      else if (jy > ny) then
        qj_traj   =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)
      else
        qj_traj   =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jy,kk)*(-sixth) +                               &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix  ,jy,kk)*half      +                               &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jy,kk)*(-half)   +                               &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_traj(ixp2,jy,kk)*(sixth)
      endif

      if (jyp1 < 1) then
        qjp1_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)
      else if (jyp1 > ny) then
        qjp1_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)
      else
        qjp1_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jyp1,kk)*(-sixth) +                             &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix  ,jyp1,kk)*half      +                             &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jyp1,kk)*(-half)   +                             &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_traj(ixp2,jyp1,kk)*(sixth)
      endif

      if (jyp2 < 1) then
        qjp2_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)
      else if (jyp2 > ny) then
        qjp2_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)
      else
        qjp2_traj =  &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jyp2,kk)*(-sixth) +                             &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix  ,jyp2,kk)*half      +                             &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jyp2,kk)*(-half)   +                             &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)                   &
            &*q_traj(ixp2,jyp2,kk)*(sixth)
      endif

!--- Lagrange interpolation in the meridional direction

!       qnew_traj(ii,jj,kk) =  &
!           &  ay_traj     *(ay_traj-one)*(ay_traj-two)*(-sixth)*qjm1_traj + &
!           & (ay_traj+one)*(ay_traj-one)*(ay_traj-two)*half     *qj_traj   + &
!           & (ay_traj+one)* ay_traj     *(ay_traj-two)*(-half)  *qjp1_traj + &
!           & (ay_traj+one)* ay_traj     *(ay_traj-one)*(sixth) *qjp2_traj

!==== Adjoint calculations

!--- Lagrange interpolation in the meridional direction

      ay   =  ( &
          &               (ay_traj-one)*(ay_traj-two)*(-sixth)*qjm1_traj + &
          &  ay_traj     *              (ay_traj-two)*(-sixth)*qjm1_traj + &
          &  ay_traj     *(ay_traj-one)*              (-sixth)*qjm1_traj + &
          &               (ay_traj-one)*(ay_traj-two)*half     *qj_traj   + &
          & (ay_traj+one)*              (ay_traj-two)*half     *qj_traj   + &
          & (ay_traj+one)*(ay_traj-one)*              half     *qj_traj   + &
          &                ay_traj     *(ay_traj-two)*(-half)  *qjp1_traj + &
          & (ay_traj+one)*              (ay_traj-two)*(-half)  *qjp1_traj + &
          & (ay_traj+one)* ay_traj     *              (-half)  *qjp1_traj + &
          &                ay_traj     *(ay_traj-one)*(sixth) *qjp2_traj + &
          & (ay_traj+one)*              (ay_traj-one)*(sixth) *qjp2_traj + &
          & (ay_traj+one)* ay_traj     *              (sixth) *qjp2_traj ) &
          &*qnew(ii,jj,kk)

      qjm1 = (  ay_traj     *(ay_traj-one)*(ay_traj-two)*(-sixth) ) &
          & *qnew(ii,jj,kk)

      qj   = ( (ay_traj+one)*(ay_traj-one)*(ay_traj-two)*half      ) &
          & *qnew(ii,jj,kk)

      qjp1 = ( (ay_traj+one)* ay_traj     *(ay_traj-two)*(-half)   ) &
          & *qnew(ii,jj,kk)

      qjp2 = ( (ay_traj+one)* ay_traj     *(ay_traj-one)*(sixth)  ) &
          & *qnew(ii,jj,kk)

      qnew(ii,jj,kk) = 0.0_kind_real

!--- Lagrange interpolation in the zonal direction

      if (jyp2 < 1) then
        ax = ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_south(ix  ,kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_south(ixp2,kk)*(sixth)                                 &
            &)*qjp2

        qjp2 =  0.0_kind_real
      else if (jyp2 > ny) then
        ax = ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_north(ix  ,kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp2,kk)*(sixth)                                 &
            &)*qjp2

        qjp2 = 0.0_kind_real
      else
        ax = ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jyp2,kk)*(-sixth) +                             &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_traj(ixm1,jyp2,kk)*(-sixth) +                             &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_traj(ixm1,jyp2,kk)*(-sixth) +                             &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix  ,jyp2,kk)*half      +                             &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ix  ,jyp2,kk)*half      +                             &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_traj(ix  ,jyp2,kk)*half      +                             &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jyp2,kk)*(-half)   +                             &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ixp1,jyp2,kk)*(-half)   +                             &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_traj(ixp1,jyp2,kk)*(-half)   +                             &
            &                 ax_traj     *(ax_traj-one)                   &
            &*q_traj(ixp2,jyp2,kk)*(sixth)   +                             &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            &*q_traj(ixp2,jyp2,kk)*(sixth)   +                             &
            &  (ax_traj+one)* ax_traj                                      &
            &*q_traj(ixp2,jyp2,kk)*(sixth)                                 &
            &)*qjp2

        q(ixm1,jyp2,kk)  = q(ixm1,jyp2,kk) + ( &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)*(-sixth)           &
            &)*qjp2

        q(ix  ,jyp2,kk) = q(ix  ,jyp2,kk) + ( &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)*half               &
            &)*qjp2

        q(ixp1,jyp2,kk) = q(ixp1,jyp2,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)*(-half)            &
            &)*qjp2

        q(ixp2,jyp2,kk) = q(ixp2,jyp2,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)*(sixth)           &
            &)*qjp2

        qjp2 = 0.0_kind_real
      endif

      if (jyp1 < 1) then
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_south(ix  ,kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_south(ixp2,kk)*(sixth)                                 &
            &)*qjp1

        qjp1 =  0.0_kind_real
      else if (jyp1 > ny) then
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_north(ix  ,kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp2,kk)*(sixth)                                 &
            &)*qjp1

        qjp1 = 0.0_kind_real
      else
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jyp1,kk)*(-sixth) +                             &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_traj(ixm1,jyp1,kk)*(-sixth) +                             &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_traj(ixm1,jyp1,kk)*(-sixth) +                             &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix  ,jyp1,kk)*half      +                             &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ix  ,jyp1,kk)*half      +                             &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_traj(ix  ,jyp1,kk)*half      +                             &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jyp1,kk)*(-half)   +                             &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ixp1,jyp1,kk)*(-half)   +                             &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_traj(ixp1,jyp1,kk)*(-half)   +                             &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_traj(ixp2,jyp1,kk)*(sixth)  +                             &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_traj(ixp2,jyp1,kk)*(sixth)  +                             &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_traj(ixp2,jyp1,kk)*(sixth)                                &
            &)*qjp1

        q(ixm1,jyp1,kk)  = q(ixm1,jyp1,kk) + ( &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)*(-sixth)           &
            &)*qjp1


        q(ix  ,jyp1,kk) = q(ix  ,jyp1,kk) + ( &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)*half               &
            &)*qjp1
        

        q(ixp1,jyp1,kk) = q(ixp1,jyp1,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)*(-half)            &
            &)*qjp1

        q(ixp2,jyp1,kk) = q(ixp2,jyp1,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)*(sixth)           &
            &)*qjp1
      endif

      if (jy < 1) then
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_south(ix  ,kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_south(ixp2,kk)*(sixth)                                 &
            &)*qj

        qj   =  0.0_kind_real
      else if (jy > ny) then
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ix  ,kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_north(ix  ,kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp2,kk)*(sixth)                                 &
            &)*qj

        qj   =  0.0_kind_real
      else
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jy,kk)*(-sixth) +                               &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_traj(ixm1,jy,kk)*(-sixth) +                               &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_traj(ixm1,jy,kk)*(-sixth) +                               &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix  ,jy,kk)*half      +                               &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ix  ,jy,kk)*half      +                               &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_traj(ix  ,jy,kk)*half      +                               &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jy,kk)*(-half)   +                               &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ixp1,jy,kk)*(-half)   +                               &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_traj(ixp1,jy,kk)*(-half)   +                               &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_traj(ixp2,jy,kk)*(sixth)  +                               &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_traj(ixp2,jy,kk)*(sixth)  +                               &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_traj(ixp2,jy,kk)*(sixth)                                  &
            &)*qj

        q(ixm1,jy,kk) = q(ixm1,jy,kk) + ( &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)*(-sixth)           &
            &)*qj

        q(ix  ,jy,kk) = q(ix  ,jy,kk) + ( &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)*half               &
            &)*qj

        q(ixp1,jy,kk) = q(ixp1,jy,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)*(-half)            &
            &)*qj

        q(ixp2,jy,kk) = q(ixp2,jy,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)*(sixth)           &
            &)*qj

        qj = 0.0_kind_real
      endif

      if (jym1 < 1) then
        ax = ax + ( &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj      *              (ax_traj-two)                   &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &  ax_traj      *(ax_traj-one)                                 &
            & *q_south(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_south(ix,  kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ix,  kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_south(ix,  kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_south(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            & (ax_traj+one) *              (ax_traj-one)                   &
            & *q_south(ixp2,kk)*(sixth)  +                              &
            & (ax_traj+one) * ax_traj                                      &
            & *q_south(ixp2,kk)*(sixth)                                 &
            &)*qjm1

        qjm1 =  0.0_kind_real
      else if (jym1 > ny) then
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_north(ixm1,kk)*(-sixth) +                              &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_north(ix,  kk)*half      +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ix,  kk)*half      +                              &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_north(ix,  kk)*half      +                              &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp1,kk)*(-half)   +                              &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_north(ixp2,kk)*(sixth)  +                              &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_north(ixp2,kk)*(sixth)                                 &
            &)*qjm1

        qjm1 =  0.0_kind_real
      else
        ax = ax + ( &
            &               (ax_traj-one)*(ax_traj-two)                    &
            & *q_traj(ixm1,jym1,kk)*(-sixth) +                             &
            &  ax_traj     *              (ax_traj-two)                    &
            & *q_traj(ixm1,jym1,kk)*(-sixth) +                             &
            &  ax_traj     *(ax_traj-one)                                  &
            & *q_traj(ixm1,jym1,kk)*(-sixth) +                             &
            &                (ax_traj-one)*(ax_traj-two)                   &
            & *q_traj(ix,  jym1,kk)*half      +                             &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ix,  jym1,kk)*half      +                             &
            &  (ax_traj+one)*(ax_traj-one)                                 &
            & *q_traj(ix,  jym1,kk)*half      +                             &
            &                 ax_traj     *(ax_traj-two)                   &
            & *q_traj(ixp1,jym1,kk)*(-half)   +                             &
            &  (ax_traj+one)*              (ax_traj-two)                   &
            & *q_traj(ixp1,jym1,kk)*(-half)   +                             &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_traj(ixp1,jym1,kk)*(-half)   +                             &
            &                 ax_traj     *(ax_traj-one)                   &
            & *q_traj(ixp2,jym1,kk)*(sixth)  +                             &
            &  (ax_traj+one)*              (ax_traj-one)                   &
            & *q_traj(ixp2,jym1,kk)*(sixth)  +                             &
            &  (ax_traj+one)* ax_traj                                      &
            & *q_traj(ixp2,jym1,kk)*(sixth)                                &
            &)*qjm1

        q(ixm1,jym1,kk) = q(ixm1,jym1,kk) + ( &
            &  ax_traj     *(ax_traj-one)*(ax_traj-two)*(-sixth)           &
            &)*qjm1

        q(ix,  jym1,kk) = q(ix,  jym1,kk) + ( &
            &  (ax_traj+one)*(ax_traj-one)*(ax_traj-two)*half               &
            &)*qjm1

        q(ixp1,jym1,kk) = q(ixp1,jym1,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-two)*(-half)            &
            &)*qjm1

        q(ixp2,jym1,kk) = q(ixp2,jym1,kk) + ( &
            &  (ax_traj+one)* ax_traj     *(ax_traj-one)*(sixth)           &
            &)*qjm1

        qjm1 =  0.0_kind_real
      endif

!--- find the interpolation point (adjoint)

      v(ii,jj,kk) = ay * (-dt/deltay)
      ay = 0.0_kind_real

      u(ii,jj,kk) = ax * (-dt/deltax)
      ax = 0.0_kind_real

    enddo
  enddo
enddo

end subroutine advect_pv_ad
