function det ( x1, y1, z1, x2, y2, z2, x0, y0, z0 )
!*****************************************************************************80
!
!!  Statement function:
!
!  DET(X1,...,Z0) >= 0 if and only if (X0,Y0,Z0) is in the
!  (closed) left hemisphere defined by the plane containing (0,0,0),
!  (X1,Y1,Z1), and (X2,Y2,Z2), where left is defined relative to an 
!  observer at (X1,Y1,Z1) facing (X2,Y2,Z2).
!
  use tools_kinds, only: kind_real
  implicit none

  real ( kind = kind_real ) x1
  real ( kind = kind_real ) y1
  real ( kind = kind_real ) z1
  real ( kind = kind_real ) x2
  real ( kind = kind_real ) y2
  real ( kind = kind_real ) z2
  real ( kind = kind_real ) x0
  real ( kind = kind_real ) y0
  real ( kind = kind_real ) z0
  real ( kind = kind_real ) det

  det = x0*(y1*z2-y2*z1) - y0*(x1*z2-x2*z1) + z0*(x1*y2-x2*y1)

end

function jrand ( n, ix, iy, iz )

!*****************************************************************************80
!
!! JRAND returns a random integer between 1 and N.
!
!  Discussion:
!
!   This function returns a uniformly distributed pseudorandom integer 
!   in the range 1 to N.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:  
!
!    Brian Wichmann, David Hill, 
!    An Efficient and Portable Pseudo-random Number Generator,
!    Applied Statistics, 
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum value to be returned.
!
!    Input/output, integer ( kind = 4 ) IX, IY, IZ = seeds initialized to 
!    values in the range 1 to 30,000 before the first call to JRAND, and 
!    not altered between subsequent calls (unless a sequence of random 
!    numbers is to be repeated by reinitializing the seeds).
!
!    Output, integer ( kind = 4 ) JRAND, a random integer in the range 1 to N.
!
!  Local parameters:
!
!    U = Pseudo-random number uniformly distributed in the interval (0,1).
!    X = Pseudo-random number in the range 0 to 3 whose fractional part is U.
!
  use tools_kinds, only: kind_real
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jrand
  integer ( kind = 4 ) n
  real ( kind = kind_real ) u
  real ( kind = kind_real ) x

  ix = mod ( 171 * ix, 30269 )
  iy = mod ( 172 * iy, 30307 )
  iz = mod ( 170 * iz, 30323 )

  x = ( real ( ix, kind = kind_real ) / 30269.0_kind_real ) &
    + ( real ( iy, kind = kind_real ) / 30307.0_kind_real ) &
    + ( real ( iz, kind = kind_real ) / 30323.0_kind_real )

  u = x - int ( x )
  jrand = int ( real ( n, kind = kind_real ) * u ) + 1

  return
end
function left ( x1, y1, z1, x2, y2, z2, x0, y0, z0 )

!*****************************************************************************80
!
!! LEFT determines whether a node is to the left of a plane through the origin.
!
!  Discussion:
!
!    This function determines whether node N0 is in the
!    (closed) left hemisphere defined by the plane containing
!    N1, N2, and the origin, where left is defined relative to
!    an observer at N1 facing N2.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = kind_real ) X1, Y1, Z1 = Coordinates of N1.
!
!    Input, real ( kind = kind_real ) X2, Y2, Z2 = Coordinates of N2.
!
!    Input, real ( kind = kind_real ) X0, Y0, Z0 = Coordinates of N0.
!
!    Output, logical LEFT = TRUE if and only if N0 is in the closed
!    left hemisphere.
!
  use tools_kinds, only: kind_real
  implicit none

  logical              left
  real ( kind = kind_real ) x0
  real ( kind = kind_real ) x1
  real ( kind = kind_real ) x2
  real ( kind = kind_real ) y0
  real ( kind = kind_real ) y1
  real ( kind = kind_real ) y2
  real ( kind = kind_real ) z0
  real ( kind = kind_real ) z1
  real ( kind = kind_real ) z2
!
!  LEFT = TRUE iff <N0,N1 X N2> = det(N0,N1,N2) >= 0.
!
  left = x0 * ( y1 * z2 - y2 * z1 ) &
       - y0 * ( x1 * z2 - x2 * z1 ) &
       + z0 * ( x1 * y2 - x2 * y1 ) >= 0.0_kind_real

  return
end
function lstptr ( lpl, nb, list, lptr )

!*****************************************************************************80
!
!! LSTPTR returns the index of NB in the adjacency list.
!
!  Discussion:
!
!    This function returns the index (LIST pointer) of NB in
!    the adjacency list for N0, where LPL = LEND(N0).
!
!    This function is identical to the similarly named function in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LPL, is LEND(N0).
!
!    Input, integer ( kind = 4 ) NB, index of the node whose pointer is to 
!    be returned.  NB must be connected to N0.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), the data 
!    structure defining the triangulation, created by TRMESH.
!
!    Output, integer ( kind = 4 ) LSTPTR, pointer such that LIST(LSTPTR) = NB or
!    LIST(LSTPTR) = -NB, unless NB is not a neighbor of N0, in which 
!    case LSTPTR = LPL.
!
!  Local parameters:
!
!    LP = LIST pointer
!    ND = Nodal index
!
  implicit none

  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nd

  lp = lptr(lpl)

  do

    nd = list(lp)

    if ( nd == nb ) then
      exit
    end if

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do

  lstptr = lp

  return
end
function swptst ( n1, n2, n3, n4, x, y, z )

!*****************************************************************************80
!
!! SWPTST decides whether to replace a diagonal arc by the other.
!
!  Discussion:
!
!    This function decides whether or not to replace a
!    diagonal arc in a quadrilateral with the other diagonal.
!    The decision will be to swap (SWPTST = TRUE) if and only
!    if N4 lies above the plane (in the half-space not contain-
!    ing the origin) defined by (N1,N2,N3), or equivalently, if
!    the projection of N4 onto this plane is interior to the
!    circumcircle of (N1,N2,N3).  The decision will be for no
!    swap if the quadrilateral is not strictly convex.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, N4, indexes of the four nodes 
!    defining the quadrilateral with N1 adjacent to N2, and (N1,N2,N3) in 
!    counterclockwise order.  The arc connecting N1 to N2 should be replaced 
!    by an arc connecting N3 to N4 if SWPTST = TRUE.  Refer to subroutine SWAP.
!
!    Input, real ( kind = kind_real ) X(N), Y(N), Z(N), the coordinates of the nodes. 
!
!    Output, logical SWPTST, TRUE if and only if the arc connecting N1
!    and N2 should be swapped for an arc connecting N3 and N4.
!
!  Local parameters:
!
!    DX1,DY1,DZ1 = Coordinates of N4->N1
!    DX2,DY2,DZ2 = Coordinates of N4->N2
!    DX3,DY3,DZ3 = Coordinates of N4->N3
!    X4,Y4,Z4 =    Coordinates of N4
!
  use tools_kinds, only: kind_real
  implicit none

  real ( kind = kind_real ) dx1
  real ( kind = kind_real ) dx2
  real ( kind = kind_real ) dx3
  real ( kind = kind_real ) dy1
  real ( kind = kind_real ) dy2
  real ( kind = kind_real ) dy3
  real ( kind = kind_real ) dz1
  real ( kind = kind_real ) dz2
  real ( kind = kind_real ) dz3
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  logical              swptst
  real ( kind = kind_real ) x(*)
  real ( kind = kind_real ) x4
  real ( kind = kind_real ) y(*)
  real ( kind = kind_real ) y4
  real ( kind = kind_real ) z(*)
  real ( kind = kind_real ) z4

  x4 = x(n4)
  y4 = y(n4)
  z4 = z(n4)
  dx1 = x(n1) - x4
  dx2 = x(n2) - x4
  dx3 = x(n3) - x4
  dy1 = y(n1) - y4
  dy2 = y(n2) - y4
  dy3 = y(n3) - y4
  dz1 = z(n1) - z4
  dz2 = z(n2) - z4
  dz3 = z(n3) - z4
!
!  N4 lies above the plane of (N1,N2,N3) iff N3 lies above
!  the plane of (N2,N1,N4) iff Det(N3-N4,N2-N4,N1-N4) =
!  (N3-N4,N2-N4 X N1-N4) > 0.
!
  swptst =  dx3 * ( dy2 * dz1 - dy1 * dz2 ) &
          - dy3 * ( dx2 * dz1 - dx1 * dz2 ) &
          + dz3 * ( dx2 * dy1 - dx1 * dy2 ) > 0.0_kind_real

  return
end

!----------------------------------------------------------------------
! Module: tools_stripack.f90
!> Purpose: stripack routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_stripack

use tools_kinds, only: kind_real

implicit none

private
public :: areas,trans,trfind,trlist,trmesh

contains

subroutine addnod ( nst, k, x, y, z, list, lptr, lend, lnew, ier )

!*****************************************************************************80
!
!! ADDNOD adds a node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds node K to a triangulation of the
!    convex hull of nodes 1, ..., K-1, producing a triangulation
!    of the convex hull of nodes 1, ..., K.
!
!    The algorithm consists of the following steps:  node K
!    is located relative to the triangulation (TRFIND), its
!    index is added to the data structure (INTADD or BDYADD),
!    and a sequence of swaps (SWPTST and SWAP) are applied to
!    the arcs opposite K so that all arcs incident on node K
!    and opposite node K are locally optimal (satisfy the circumcircle test).  
!
!    Thus, if a Delaunay triangulation of nodes 1 through K-1 is input, 
!    a Delaunay triangulation of nodes 1 through K will be output.
!
!  Modified:
!
!    15 May 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the index of a node at which TRFIND 
!    begins its search.  Search time depends on the proximity of this node to 
!    K.  If NST < 1, the search is begun at node K-1.
!
!    Input, integer ( kind = 4 ) K, the nodal index (index for X, Y, Z, and 
!    LEND) of the new node to be added.  4 <= K.
!
!    Input, real ( kind = kind_real ) X(K), Y(K), Z(K), the coordinates of the nodes.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(K), 
!    LNEW.  On input, the data structure associated with the triangulation of 
!    nodes 1 to K-1.  On output, the data has been updated to include node 
!    K.  The array lengths are assumed to be large enough to add node K. 
!    Refer to TRMESH.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!     0 if no errors were encountered.
!    -1 if K is outside its valid range on input.
!    -2 if all nodes (including K) are collinear (lie on a common geodesic).
!     L if nodes L and K coincide for some L < K.
!
!  Local parameters:
!
!    B1,B2,B3 = Unnormalized barycentric coordinates returned by TRFIND.
!    I1,I2,I3 = Vertex indexes of a triangle containing K
!    IN1 =      Vertex opposite K:  first neighbor of IO2
!               that precedes IO1.  IN1,IO1,IO2 are in
!               counterclockwise order.
!    IO1,IO2 =  Adjacent neighbors of K defining an arc to
!               be tested for a swap
!    IST =      Index of node at which TRFIND begins its search
!    KK =       Local copy of K
!    KM1 =      K-1
!    L =        Vertex index (I1, I2, or I3) returned in IER
!               if node K coincides with a vertex
!    LP =       LIST pointer
!    LPF =      LIST pointer to the first neighbor of K
!    LPO1 =     LIST pointer to IO1
!    LPO1S =    Saved value of LPO1
!    P =        Cartesian coordinates of node K
!
  implicit none

  integer ( kind = 4 ) k

  real ( kind = kind_real ) b1
  real ( kind = kind_real ) b2
  real ( kind = kind_real ) b3
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lend(k)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpo1
  integer ( kind = 4 ) lpo1s
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nst
  real ( kind = kind_real ) p(3)
  logical swptst
  real ( kind = kind_real ) x(k)
  real ( kind = kind_real ) y(k)
  real ( kind = kind_real ) z(k)

  kk = k

  if ( kk < 4 ) then
    ier = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADDNOD - Fatal error!'
    write ( *, '(a)' ) '  K < 4.'
    stop
  end if
!
!  Initialization:
!
  km1 = kk - 1
  ist = nst
  if ( ist < 1 ) then
    ist = km1
  end if

  p(1) = x(kk)
  p(2) = y(kk)
  p(3) = z(kk)
!
!  Find a triangle (I1,I2,I3) containing K or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed
!  from node K.
!
  call trfind ( ist, p, km1, x, y, z, list, lptr, lend, b1, b2, b3, &
    i1, i2, i3 )
!
!  Test for collinear or duplicate nodes.
!
  if ( i1 == 0 ) then
    ier = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADDNOD - Fatal error!'
    write ( *, '(a)' ) '  The nodes are coplanar.'
    stop
  end if

  if ( i3 /= 0 ) then

    l = i1

    if ( .not.(abs(p(1)-x(l))>0.0) .and. .not.(abs(p(2)-y(l))>0.0)  .and. .not.(abs(p(3)-z(l))>0.0) ) then
      ier = l
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDNOD - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Node ', l, ' is equal to node ', k
      stop
    end if

    l = i2

    if ( .not.(abs(p(1)-x(l))>0.0) .and. .not.(abs(p(2)-y(l))>0.0)  .and. .not.(abs(p(3)-z(l))>0.0) ) then
      ier = l
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDNOD - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Node ', l, ' is equal to node ', k
      stop
    end if

    l = i3
    if ( .not.(abs(p(1)-x(l))>0.0) .and. .not.(abs(p(2)-y(l))>0.0)  .and. .not.(abs(p(3)-z(l))>0.0) ) then
      ier = l
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDNOD - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Node ', l, ' is equal to node ', k
      stop
    end if

    call intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )

  else

    if ( i1 /= i2 ) then
      call bdyadd ( kk, i1,i2, list, lptr, lend, lnew )
    else
      call covsph ( kk, i1, list, lptr, lend, lnew )
    end if

  end if

  ier = 0
!
!  Initialize variables for optimization of the triangulation.
!
  lp = lend(kk)
  lpf = lptr(lp)
  io2 = list(lpf)
  lpo1 = lptr(lpf)
  io1 = abs ( list(lpo1) )
!
!  Begin loop: find the node opposite K.
!
  do

    lp = lstptr ( lend(io1), io2, list, lptr )

    if ( 0 <= list(lp) ) then

      lp = lptr(lp)
      in1 = abs ( list(lp) )
!
!  Swap test:  if a swap occurs, two new arcs are
!  opposite K and must be tested.
!
      lpo1s = lpo1

      if ( .not. swptst ( in1, kk, io1, io2, x, y, z ) ) then

        if ( lpo1 == lpf .or. list(lpo1) < 0 ) then
          exit
        end if

        io2 = io1
        lpo1 = lptr(lpo1)
        io1 = abs ( list(lpo1) )
        cycle

      end if

      call swap ( in1, kk, io1, io2, list, lptr, lend, lpo1 )
!
!  A swap is not possible because KK and IN1 are already
!  adjacent.  This error in SWPTST only occurs in the
!  neutral case and when there are nearly duplicate nodes.
!
      if ( lpo1 /= 0 ) then
        io1 = in1
        cycle
      end if

      lpo1 = lpo1s

    end if
!
!  No swap occurred.  Test for termination and reset IO2 and IO1.
!
    if ( lpo1 == lpf .or. list(lpo1) < 0 ) then
      exit
    end if

    io2 = io1
    lpo1 = lptr(lpo1)
    io1 = abs ( list(lpo1) )

  end do

  return
end
function areas ( v1, v2, v3 )

!*****************************************************************************80
!
!! AREAS computes the area of a spherical triangle on the unit sphere.
!
!  Discussion:
!
!    This function returns the area of a spherical triangle
!    on the unit sphere.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = kind_real ) V1(3), V2(3), V3(3), the Cartesian coordinates
!    of unit vectors (the three triangle vertices in any order).  These 
!    vectors, if nonzero, are implicitly scaled to have length 1.
!
!    Output, real ( kind = kind_real ) AREAS, the area of the spherical triangle 
!    defined by V1, V2, and V3, in the range 0 to 2*PI (the area of a
!    hemisphere).  AREAS = 0 (or 2*PI) if and only if V1, V2, and V3 lie in (or
!    close to) a plane containing the origin.
!
!  Local parameters:
!
!    A1,A2,A3 =    Interior angles of the spherical triangle.
!
!    CA1,CA2,CA3 = cos(A1), cos(A2), and cos(A3), respectively.
!
!    DV1,DV2,DV3 = copies of V1, V2, and V3.
!
!    I =           DO-loop index and index for Uij.
!
!    S12,S23,S31 = Sum of squared components of U12, U23, U31.
!
!    U12,U23,U31 = Unit normal vectors to the planes defined by
!                 pairs of triangle vertices.
!
  implicit none

  real ( kind = kind_real ) a1
  real ( kind = kind_real ) a2
  real ( kind = kind_real ) a3
  real ( kind = kind_real ) areas
  real ( kind = kind_real ) ca1
  real ( kind = kind_real ) ca2
  real ( kind = kind_real ) ca3
  real ( kind = kind_real ) dv1(3)
  real ( kind = kind_real ) dv2(3)
  real ( kind = kind_real ) dv3(3)
  real ( kind = kind_real ) s12
  real ( kind = kind_real ) s23
  real ( kind = kind_real ) s31
  real ( kind = kind_real ) u12(3)
  real ( kind = kind_real ) u23(3)
  real ( kind = kind_real ) u31(3)
  real ( kind = kind_real ) v1(3)
  real ( kind = kind_real ) v2(3)
  real ( kind = kind_real ) v3(3)

  dv1(1:3) = v1(1:3)
  dv2(1:3) = v2(1:3)
  dv3(1:3) = v3(1:3)
!
!  Compute cross products Uij = Vi X Vj.
!
  u12(1) = dv1(2) * dv2(3) - dv1(3) * dv2(2)
  u12(2) = dv1(3) * dv2(1) - dv1(1) * dv2(3)
  u12(3) = dv1(1) * dv2(2) - dv1(2) * dv2(1)

  u23(1) = dv2(2) * dv3(3) - dv2(3) * dv3(2)
  u23(2) = dv2(3) * dv3(1) - dv2(1) * dv3(3)
  u23(3) = dv2(1) * dv3(2) - dv2(2) * dv3(1)

  u31(1) = dv3(2) * dv1(3) - dv3(3) * dv1(2)
  u31(2) = dv3(3) * dv1(1) - dv3(1) * dv1(3)
  u31(3) = dv3(1) * dv1(2) - dv3(2) * dv1(1)
!
!  Normalize Uij to unit vectors.
!
  s12 = dot_product ( u12(1:3), u12(1:3) )
  s23 = dot_product ( u23(1:3), u23(1:3) )
  s31 = dot_product ( u31(1:3), u31(1:3) )
!
!  Test for a degenerate triangle associated with collinear vertices.
!
  if ( .not.(abs(s12)>0.0) .or. .not.(abs(s23)>0.0) .or. .not.(abs(s31)>0.0) ) then
    areas = 0.0_kind_real
    return
  end if

  s12 = sqrt ( s12 )
  s23 = sqrt ( s23 )
  s31 = sqrt ( s31 )

  u12(1:3) = u12(1:3) / s12
  u23(1:3) = u23(1:3) / s23
  u31(1:3) = u31(1:3) / s31
!
!  Compute interior angles Ai as the dihedral angles between planes:
!  CA1 = cos(A1) = -<U12,U31>
!  CA2 = cos(A2) = -<U23,U12>
!  CA3 = cos(A3) = -<U31,U23>
!
  ca1 = - dot_product ( u12(1:3), u31(1:3) )
  ca2 = - dot_product ( u23(1:3), u12(1:3) )
  ca3 = - dot_product ( u31(1:3), u23(1:3) )

  ca1 = max ( ca1, -1.0_kind_real )
  ca1 = min ( ca1, +1.0_kind_real )
  ca2 = max ( ca2, -1.0_kind_real )
  ca2 = min ( ca2, +1.0_kind_real )
  ca3 = max ( ca3, -1.0_kind_real )
  ca3 = min ( ca3, +1.0_kind_real )

  a1 = acos ( ca1 )
  a2 = acos ( ca2 )
  a3 = acos ( ca3 )
!
!  Compute AREAS = A1 + A2 + A3 - PI.
!
  areas = a1 + a2 + a3 - acos ( -1.0_kind_real )

  if ( areas < 0.0_kind_real ) then
    areas = 0.0_kind_real
  end if

  return
end
subroutine bdyadd ( kk, i1, i2, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! BDYADD adds a boundary node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds a boundary node to a triangulation
!    of a set of KK-1 points on the unit sphere.  The data
!    structure is updated with the insertion of node KK, but no
!    optimization is performed.
!
!    This routine is identical to the similarly named routine
!    in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK, the index of a node to be connected to 
!    the sequence of all visible boundary nodes.  1 <= KK and
!    KK must not be equal to I1 or I2.
!
!    Input, integer ( kind = 4 ) I1, the first (rightmost as viewed from KK) 
!    boundary node in the triangulation that is visible from
!    node KK (the line segment KK-I1 intersects no arcs.
!
!    Input, integer ( kind = 4 ) I2, the last (leftmost) boundary node that 
!    is visible from node KK.  I1 and I2 may be determined by TRFIND.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    LNEW, the triangulation data structure created by TRMESH.  
!    Nodes I1 and I2 must be included 
!    in the triangulation.  On output, the data structure is updated with
!    the addition of node KK.  Node KK is connected to I1, I2, and
!    all boundary nodes in between.
!
!  Local parameters:
!
!    K =     Local copy of KK
!    LP =    LIST pointer
!    LSAV =  LIST pointer
!    N1,N2 = Local copies of I1 and I2, respectively
!    NEXT =  Boundary node visible from K
!    NSAV =  Boundary node visible from K
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nsav

  k = kk
  n1 = i1
  n2 = i2
!
!  Add K as the last neighbor of N1.
!
  lp = lend(n1)
  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = -k
  lptr(lnew) = lsav
  lend(n1) = lnew
  lnew = lnew + 1
  next = -list(lp)
  list(lp) = next
  nsav = next
!
!  Loop on the remaining boundary nodes between N1 and N2,
!  adding K as the first neighbor.
!
  do

    lp = lend(next)
    call insert ( k, lp, list, lptr, lnew )

    if ( next == n2 ) then
      exit
    end if

    next = -list(lp)
    list(lp) = next

  end do
!
!  Add the boundary nodes between N1 and N2 as neighbors of node K.
!
  lsav = lnew
  list(lnew) = n1
  lptr(lnew) = lnew + 1
  lnew = lnew + 1
  next = nsav

  do

    if ( next == n2 ) then
      exit
    end if

    list(lnew) = next
    lptr(lnew) = lnew + 1
    lnew = lnew + 1
    lp = lend(next)
    next = list(lp)

  end do

  list(lnew) = -n2
  lptr(lnew) = lsav
  lend(k) = lnew
  lnew = lnew + 1

  return
end
subroutine covsph ( kk, n0, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! COVSPH connects an exterior node to boundary nodes, covering the sphere.
!
!  Discussion:
!
!    This subroutine connects an exterior node KK to all
!    boundary nodes of a triangulation of KK-1 points on the
!    unit sphere, producing a triangulation that covers the
!    sphere.  The data structure is updated with the addition
!    of node KK, but no optimization is performed.  All 
!    boundary nodes must be visible from node KK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK = Index of the node to be connected to the
!    set of all boundary nodes.  4 <= KK.
!
!    Input, integer ( kind = 4 ) N0 = Index of a boundary node (in the range
!    1 to KK-1).  N0 may be determined by TRFIND.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    LNEW, the triangulation data structure created by TRMESH.  Node N0 must
!    be included in the triangulation.  On output, updated with the addition 
!    of node KK as the last entry.  The updated triangulation contains no
!    boundary nodes.
!
!  Local parameters:
!
!    K =     Local copy of KK
!    LP =    LIST pointer
!    LSAV =  LIST pointer
!    NEXT =  Boundary node visible from K
!    NST =   Local copy of N0
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nst

  k = kk
  nst = n0
!
!  Traverse the boundary in clockwise order, inserting K as
!  the first neighbor of each boundary node, and converting
!  the boundary node to an interior node.
!
  next = nst

  do

    lp = lend(next)
    call insert ( k, lp, list, lptr, lnew )
    next = -list(lp)
    list(lp) = next

    if ( next == nst ) then
      exit
    end if

  end do
!
!  Traverse the boundary again, adding each node to K's adjacency list.
!
  lsav = lnew

  do

    lp = lend(next)
    list(lnew) = next
    lptr(lnew) = lnew + 1
    lnew = lnew + 1
    next = list(lp)

    if ( next == nst ) then
      exit
    end if

  end do

  lptr(lnew-1) = lsav
  lend(k) = lnew - 1

  return
end
subroutine insert ( k, lp, list, lptr, lnew )

!*****************************************************************************80
!
!! INSERT inserts K as a neighbor of N1.
!
!  Discussion:
!
!    This subroutine inserts K as a neighbor of N1 following
!    N2, where LP is the LIST pointer of N2 as a neighbor of
!    N1.  Note that, if N2 is the last neighbor of N1, K will
!    become the first neighbor (even if N1 is a boundary node).
!
!    This routine is identical to the similarly named routine in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the index of the node to be inserted.
!
!    Input, integer ( kind = 4 ) LP, the LIST pointer of N2 as a neighbor of N1.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LNEW, 
!    the data structure defining the triangulation, created by TRMESH.
!    On output, updated with the addition of node K.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav

  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = k
  lptr(lnew) = lsav
  lnew = lnew + 1

  return
end
function inside ( p, lv, xv, yv, zv, nv, listv, ier )

!*****************************************************************************80
!
!! INSIDE determines if a point is inside a polygonal region.
!
!  Discussion:
!
!    This function locates a point P relative to a polygonal
!    region R on the surface of the unit sphere, returning
!    INSIDE = TRUE if and only if P is contained in R.  R is
!    defined by a cyclically ordered sequence of vertices which
!    form a positively-oriented simple closed curve.  Adjacent
!    vertices need not be distinct but the curve must not be
!    self-intersecting.  Also, while polygon edges are by definition
!    restricted to a single hemisphere, R is not so
!    restricted.  Its interior is the region to the left as the
!    vertices are traversed in order.
!
!    The algorithm consists of selecting a point Q in R and
!    then finding all points at which the great circle defined
!    by P and Q intersects the boundary of R.  P lies inside R
!    if and only if there is an even number of intersection
!    points between Q and P.  Q is taken to be a point immediately
!    to the left of a directed boundary edge -- the first
!    one that results in no consistency-check failures.
!
!    If P is close to the polygon boundary, the problem is
!    ill-conditioned and the decision may be incorrect.  Also,
!    an incorrect decision may result from a poor choice of Q
!    (if, for example, a boundary edge lies on the great circle
!    defined by P and Q).  A more reliable result could be
!    obtained by a sequence of calls to INSIDE with the vertices
!    cyclically permuted before each call (to alter the
!    choice of Q).
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = kind_real ) P(3), the coordinates of the point (unit vector)
!    to be located.
!
!    Input, integer ( kind = 4 ) LV, the length of arrays XV, YV, and ZV.
!
!    Input, real ( kind = kind_real ) XV(LV), YV(LV), ZV(LV), the coordinates of unit
!    vectors (points on the unit sphere).  
!
!    Input, integer ( kind = 4 ) NV, the number of vertices in the polygon. 
!    3 <= NV <= LV.
!
!    Input, integer ( kind = 4 ) LISTV(NV), the indexes (for XV, YV, and ZV) 
!    of a cyclically-ordered (and CCW-ordered) sequence of vertices that
!    define R.  The last vertex (indexed by LISTV(NV)) is followed by the 
!    first (indexed by LISTV(1)).  LISTV entries must be in the range 1 to LV.
!
!    Output, logical INSIDE, TRUE if and only if P lies inside R unless
!    IER /= 0, in which case the value is not altered.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!    0, if no errors were encountered.
!    1, if LV or NV is outside its valid range.
!    2, if a LISTV entry is outside its valid range.
!    3, if the polygon boundary was found to be self-intersecting.  This 
!      error will not necessarily be detected.
!    4, if every choice of Q (one for each boundary edge) led to failure of 
!      some internal consistency check.  The most likely cause of this error 
!      is invalid input:  P = (0,0,0), a null or self-intersecting polygon, etc.
!
!  Local parameters:
!
!    B =         Intersection point between the boundary and
!                the great circle defined by P and Q.
!
!    BP,BQ =     <B,P> and <B,Q>, respectively, maximized over
!                intersection points B that lie between P and
!                Q (on the shorter arc) -- used to find the
!                closest intersection points to P and Q
!    CN =        Q X P = normal to the plane of P and Q
!    D =         Dot product <B,P> or <B,Q>
!    EPS =       Parameter used to define Q as the point whose
!                orthogonal distance to (the midpoint of)
!                boundary edge V1->V2 is approximately EPS/
!                (2*Cos(A/2)), where <V1,V2> = Cos(A).
!    EVEN =      TRUE iff an even number of intersection points
!                lie between P and Q (on the shorter arc)
!    I1,I2 =     Indexes (LISTV elements) of a pair of adjacent
!                boundary vertices (endpoints of a boundary
!                edge)
!    IERR =      Error flag for calls to INTRSC (not tested)
!    IMX =       Local copy of LV and maximum value of I1 and I2
!    K =         DO-loop index and LISTV index
!    K0 =        LISTV index of the first endpoint of the
!                boundary edge used to compute Q
!    LFT1,LFT2 = Logical variables associated with I1 and I2 in
!                the boundary traversal:  TRUE iff the vertex
!                is strictly to the left of Q->P (<V,CN> > 0)
!    N =         Local copy of NV
!    NI =        Number of intersections (between the boundary
!                curve and the great circle P-Q) encountered
!    PINR =      TRUE iff P is to the left of the directed
!                boundary edge associated with the closest
!                intersection point to P that lies between P
!                and Q (a left-to-right intersection as
!                viewed from Q), or there is no intersection
!                between P and Q (on the shorter arc)
!    PN,QN =     P X CN and CN X Q, respectively:  used to
!                locate intersections B relative to arc Q->P
!    Q =         (V1 + V2 + EPS*VN/VNRM)/QNRM, where V1->V2 is
!                the boundary edge indexed by LISTV(K0) ->
!                LISTV(K0+1)
!    QINR =      TRUE iff Q is to the left of the directed
!                boundary edge associated with the closest
!                intersection point to Q that lies between P
!                and Q (a right-to-left intersection as
!                viewed from Q), or there is no intersection
!                between P and Q (on the shorter arc)
!    QNRM =      Euclidean norm of V1+V2+EPS*VN/VNRM used to
!                compute (normalize) Q
!    V1,V2 =     Vertices indexed by I1 and I2 in the boundary
!                traversal
!    VN =        V1 X V2, where V1->V2 is the boundary edge
!                indexed by LISTV(K0) -> LISTV(K0+1)
!    VNRM =      Euclidean norm of VN
!
  implicit none

  integer ( kind = 4 ) lv
  integer ( kind = 4 ) nv

  real ( kind = kind_real ) b(3)
  real ( kind = kind_real ) bp
  real ( kind = kind_real ) bq
  real ( kind = kind_real ) cn(3)
  real ( kind = kind_real ) d
  real ( kind = kind_real ), parameter :: eps = 0.001_kind_real
  logical even
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) imx
  logical inside
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  logical lft1
  logical lft2
  integer ( kind = 4 ) listv(nv)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ni
  real ( kind = kind_real ) p(3)
  logical pinr
  real ( kind = kind_real ) pn(3)
  real ( kind = kind_real ) q(3)
  logical qinr
  real ( kind = kind_real ) qn(3)
  real ( kind = kind_real ) qnrm
  real ( kind = kind_real ) v1(3)
  real ( kind = kind_real ) v2(3)
  real ( kind = kind_real ) vn(3)
  real ( kind = kind_real ) vnrm
  real ( kind = kind_real ) xv(lv)
  real ( kind = kind_real ) yv(lv)
  real ( kind = kind_real ) zv(lv)
!
!  Store local parameters.
!
  imx = lv
  n = nv
!
!  Test for error 1.
!
  if ( n < 3 .or. imx < n ) then
    ier = 1
    return
  end if
!
!  Initialize K0.
!
  k0 = 0
  i1 = listv(1)

  if ( i1 < 1 .or. imx < i1 ) then
    ier = 2
    return
  end if
!
!  Increment K0 and set Q to a point immediately to the left
!  of the midpoint of edge V1->V2 = LISTV(K0)->LISTV(K0+1):
!  Q = (V1 + V2 + EPS*VN/VNRM)/QNRM, where VN = V1 X V2.
!
1 continue

  k0 = k0 + 1

  if ( n < k0 ) then
    ier = 4
    return
  end if

  i1 = listv(k0)

  if ( k0 < n ) then
    i2 = listv(k0+1)
  else
    i2 = listv(1)
  end if

  if ( i2 < 1 .or. imx < i2 ) then
    ier = 2
    return
  end if

  vn(1) = yv(i1) * zv(i2) - zv(i1) * yv(i2)
  vn(2) = zv(i1) * xv(i2) - xv(i1) * zv(i2)
  vn(3) = xv(i1) * yv(i2) - yv(i1) * xv(i2)
  vnrm = sqrt ( sum ( vn(1:3)**2 ) )

  if ( .not.(abs(vnrm)>0.0) ) then
    go to 1
  end if

  q(1) = xv(i1) + xv(i2) + eps * vn(1) / vnrm
  q(2) = yv(i1) + yv(i2) + eps * vn(2) / vnrm
  q(3) = zv(i1) + zv(i2) + eps * vn(3) / vnrm

  qnrm = sqrt ( sum ( q(1:3)**2 ) )

  q(1) = q(1) / qnrm
  q(2) = q(2) / qnrm
  q(3) = q(3) / qnrm
!
!  Compute CN = Q X P, PN = P X CN, and QN = CN X Q.
!
  cn(1) = q(2) * p(3) - q(3) * p(2)
  cn(2) = q(3) * p(1) - q(1) * p(3)
  cn(3) = q(1) * p(2) - q(2) * p(1)

  if ( .not.(abs(cn(1))>0.0) .and. .not.(abs(cn(2))>0.0)  .and. .not.(abs(cn(2))>0.0) ) then
    go to 1
  end if

  pn(1) = p(2) * cn(3) - p(3) * cn(2)
  pn(2) = p(3) * cn(1) - p(1) * cn(3)
  pn(3) = p(1) * cn(2) - p(2) * cn(1)
  qn(1) = cn(2) * q(3) - cn(3) * q(2)
  qn(2) = cn(3) * q(1) - cn(1) * q(3)
  qn(3) = cn(1) * q(2) - cn(2) * q(1)
!
!  Initialize parameters for the boundary traversal.
!
  ni = 0
  even = .true.
  bp = -2.0_kind_real
  bq = -2.0_kind_real
  pinr = .true.
  qinr = .true.
  i2 = listv(n)

  if ( i2 < 1 .or. imx < i2 ) then
    ier = 2
    return
  end if

  lft2 = 0.0_kind_real < cn(1) * xv(i2) + cn(2) * yv(i2) + cn(3) * zv(i2) 
!
!  Loop on boundary arcs I1->I2.
!
  do k = 1, n

    i1 = i2
    lft1 = lft2
    i2 = listv(k)

    if ( i2 < 1 .or. imx < i2 ) then
      ier = 2
      return
    end if

    lft2 = ( 0.0_kind_real < cn(1) * xv(i2) + cn(2) * yv(i2) + cn(3) * zv(i2) )

    if ( lft1 .eqv. lft2 ) then
      cycle
    end if
!
!  I1 and I2 are on opposite sides of Q->P.  Compute the
!  point of intersection B.
!
    ni = ni + 1
    v1(1) = xv(i1)
    v1(2) = yv(i1)
    v1(3) = zv(i1)
    v2(1) = xv(i2)
    v2(2) = yv(i2)
    v2(3) = zv(i2)

    call intrsc ( v1, v2, cn, b, ierr )
!
!  B is between Q and P (on the shorter arc) iff
!  B Forward Q->P and B Forward P->Q       iff
!  <B,QN> > 0 and 0 < <B,PN>.
!
    if ( 0.0_kind_real < dot_product ( b(1:3), qn(1:3) ) .and. &
         0.0_kind_real < dot_product ( b(1:3), pn(1:3) ) ) then
!
!  Update EVEN, BQ, QINR, BP, and PINR.
!
      even = .not. even
      d = dot_product ( b(1:3), q(1:3) )

      if ( bq < d ) then
        bq = d
        qinr = lft2
      end if

      d = dot_product ( b(1:3), p(1:3) )

      if ( bp < d ) then
        bp = d
        pinr = lft1
      end if

    end if

  end do
!
!  Test for consistency:  NI must be even and QINR must be TRUE.
!
  if ( ni /= 2 * ( ni / 2 ) .or. .not. qinr ) then
    go to 1
  end if
!
!  Test for error 3:  different values of PINR and EVEN.
!
  if ( pinr .neqv. even ) then
    ier = 3
    return
  end if

  ier = 0
  inside = even

  return
end
subroutine intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! INTADD adds an interior node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds an interior node to a triangulation
!    of a set of points on the unit sphere.  The data structure
!    is updated with the insertion of node KK into the triangle
!    whose vertices are I1, I2, and I3.  No optimization of the
!    triangulation is performed.
!
!    This routine is identical to the similarly named routine in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK, the index of the node to be inserted. 
!    1 <= KK and KK must not be equal to I1, I2, or I3.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, indexes of the 
!    counterclockwise-ordered sequence of vertices of a triangle which contains 
!    node KK.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), LNEW, 
!    the data structure defining the triangulation, created by TRMESH.  Triangle
!    (I1,I2,I3) must be included in the triangulation.
!    On output, updated with the addition of node KK.  KK
!    will be connected to nodes I1, I2, and I3.
!
!  Local parameters:
!
!    K =        Local copy of KK
!    LP =       LIST pointer
!    N1,N2,N3 = Local copies of I1, I2, and I3
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  k = kk
!
!  Initialization.
!
  n1 = i1
  n2 = i2
  n3 = i3
!
!  Add K as a neighbor of I1, I2, and I3.
!
  lp = lstptr ( lend(n1), n2, list, lptr )
  call insert ( k, lp, list, lptr, lnew )

  lp = lstptr ( lend(n2), n3, list, lptr )
  call insert ( k, lp, list, lptr, lnew )

  lp = lstptr ( lend(n3), n1, list, lptr )
  call insert ( k, lp, list, lptr, lnew )
!
!  Add I1, I2, and I3 as neighbors of K.
!
  list(lnew) = n1
  list(lnew+1) = n2
  list(lnew+2) = n3
  lptr(lnew) = lnew + 1
  lptr(lnew+1) = lnew + 2
  lptr(lnew+2) = lnew
  lend(k) = lnew + 2
  lnew = lnew + 3

  return
end
subroutine intrsc ( p1, p2, cn, p, ier )

!*****************************************************************************80
!
!! INTSRC finds the intersection of two great circles.
!
!  Discussion:
!
!    Given a great circle C and points P1 and P2 defining an
!    arc A on the surface of the unit sphere, where A is the
!    shorter of the two portions of the great circle C12
!    associated with P1 and P2, this subroutine returns the point
!    of intersection P between C and C12 that is closer to A.
!    Thus, if P1 and P2 lie in opposite hemispheres defined by
!    C, P is the point of intersection of C with A.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = kind_real ) P1(3), P2(3), the coordinates of unit vectors.
!
!    Input, real ( kind = kind_real ) CN(3), the coordinates of a nonzero vector 
!    which defines C as the intersection of the plane whose normal is CN 
!    with the unit sphere.  Thus, if C is to be the great circle defined
!    by P and Q, CN should be P X Q.
!
!    Output, real ( kind = kind_real ) P(3), point of intersection defined above 
!    unless IER is not 0, in which case P is not altered.
!
!    Output, integer ( kind = 4 ) IER, error indicator.
!    0, if no errors were encountered.
!    1, if <CN,P1> = <CN,P2>.  This occurs iff P1 = P2 or CN = 0 or there are
!      two intersection points at the same distance from A.
!    2, if P2 = -P1 and the definition of A is therefore ambiguous.
!
!  Local parameters:
!
!    D1 =  <CN,P1>
!    D2 =  <CN,P2>
!    I =   DO-loop index
!    PP =  P1 + T*(P2-P1) = Parametric representation of the
!          line defined by P1 and P2
!    PPN = Norm of PP
!    T =   D1/(D1-D2) = Parameter value chosen so that PP lies
!          in the plane of C
!
  implicit none

  real ( kind = kind_real ) cn(3)
  real ( kind = kind_real ) d1
  real ( kind = kind_real ) d2
  integer ( kind = 4 ) ier
  real ( kind = kind_real ) p(3)
  real ( kind = kind_real ) p1(3)
  real ( kind = kind_real ) p2(3)
  real ( kind = kind_real ) pp(3)
  real ( kind = kind_real ) ppn
  real ( kind = kind_real ) t

  d1 = dot_product ( cn(1:3), p1(1:3) )
  d2 = dot_product ( cn(1:3), p2(1:3) )

  if ( .not.(abs(d1-d2)>0.0) ) then
    ier = 1
    return
  end if
!
!  Solve for T such that <PP,CN> = 0 and compute PP and PPN.
!
  t = d1 / ( d1 - d2 )

  pp(1:3) = p1(1:3) + t * ( p2(1:3) - p1(1:3) )

  ppn = dot_product ( pp(1:3), pp(1:3) )
!
!  PPN = 0 iff PP = 0 iff P2 = -P1 (and T = .5).
!
  if ( .not.(abs(ppn)>0.0) ) then
    ier = 2
    return
  end if

  ppn = sqrt ( ppn )
!
!  Compute P = PP/PPN.
!
  p(1:3) = pp(1:3) / ppn

  ier = 0

  return
end
function nbcnt ( lpl, lptr )

!*****************************************************************************80
!
!! NBCNT returns the number of neighbors of a node.
!
!  Discussion:
!
!    This function returns the number of neighbors of a node
!    N0 in a triangulation created by TRMESH.
!
!    The number of neighbors also gives the order of the Voronoi
!    polygon containing the point.  Thus, a neighbor count of 6
!    means the node is contained in a 6-sided Voronoi region.
!
!    This function is identical to the similarly named function in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LPL = LIST pointer to the last neighbor of N0;
!    LPL = LEND(N0).
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), pointers associated with LIST.
!
!    Output, integer ( kind = 4 ) NBCNT, the number of neighbors of N0.
!
!  Local parameters:
!
!    K =  Counter for computing the number of neighbors.
!
!    LP = LIST pointer
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) nbcnt

  lp = lpl
  k = 1

  do

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

    k = k + 1

  end do

  nbcnt = k

  return
end
function nearnd ( p, ist, n, x, y, z, list, lptr, lend, al )

!*****************************************************************************80
!
!! NEARND returns the nearest node to a given point.
!
!  Discussion:
!
!    Given a point P on the surface of the unit sphere and a
!    Delaunay triangulation created by TRMESH, this
!    function returns the index of the nearest triangulation
!    node to P.
!
!    The algorithm consists of implicitly adding P to the
!    triangulation, finding the nearest neighbor to P, and
!    implicitly deleting P from the triangulation.  Thus, it
!    is based on the fact that, if P is a node in a Delaunay
!    triangulation, the nearest node to P is a neighbor of P.
!
!    For large values of N, this procedure will be faster than
!    the naive approach of computing the distance from P to every node.
!
!    Note that the number of candidates for NEARND (neighbors of P) 
!    is limited to LMAX defined in the PARAMETER statement below.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = kind_real ) P(3), the Cartesian coordinates of the point P to 
!    be located relative to the triangulation.  It is assumed 
!    that P(1)**2 + P(2)**2 + P(3)**2 = 1, that is, that the
!    point lies on the unit sphere.
!
!    Input, integer ( kind = 4 ) IST, the index of the node at which the search
!    is to begin.  The search time depends on the proximity of this 
!    node to P.  If no good candidate is known, any value between
!    1 and N will do.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    N must be at least 3.
!
!    Input, real ( kind = kind_real ) X(N), Y(N), Z(N), the Cartesian coordinates of
!    the nodes.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), 
!    the data structure defining the triangulation, created by TRMESH.
!
!    Output, real ( kind = kind_real ) AL, the arc length between P and node NEARND.
!    Because both points are on the unit sphere, this is also
!    the angular separation in radians.
!
!    Output, integer ( kind = 4 ) NEARND, the index of the nearest node to P.
!    NEARND will be 0 if N < 3 or the triangulation data structure 
!    is invalid.
!
!  Local parameters:
!
!    B1,B2,B3 =  Unnormalized barycentric coordinates returned by TRFIND
!    DS1 =       (Negative cosine of the) distance from P to N1
!    DSR =       (Negative cosine of the) distance from P to NR
!    DX1,..DZ3 = Components of vectors used by the swap test
!    I1,I2,I3 =  Nodal indexes of a triangle containing P, or
!                the rightmost (I1) and leftmost (I2) visible
!                boundary nodes as viewed from P
!    L =         Length of LISTP/LPTRP and number of neighbors of P
!    LMAX =      Maximum value of L
!    LISTP =     Indexes of the neighbors of P
!    LPTRP =     Array of pointers in 1-1 correspondence with LISTP elements
!    LP =        LIST pointer to a neighbor of N1 and LISTP pointer
!    LP1,LP2 =   LISTP indexes (pointers)
!    LPL =       Pointer to the last neighbor of N1
!    N1 =        Index of a node visible from P
!    N2 =        Index of an endpoint of an arc opposite P
!    N3 =        Index of the node opposite N1->N2
!    NN =        Local copy of N
!    NR =        Index of a candidate for the nearest node to P
!    NST =       Index of the node at which TRFIND begins the search
!
  implicit none

  integer ( kind = 4 ), parameter :: lmax = 25
  integer ( kind = 4 ) n

  real ( kind = kind_real ) al
  real ( kind = kind_real ) b1
  real ( kind = kind_real ) b2
  real ( kind = kind_real ) b3
  real ( kind = kind_real ) ds1
  real ( kind = kind_real ) dsr
  real ( kind = kind_real ) dx1
  real ( kind = kind_real ) dx2
  real ( kind = kind_real ) dx3
  real ( kind = kind_real ) dy1
  real ( kind = kind_real ) dy2
  real ( kind = kind_real ) dy3
  real ( kind = kind_real ) dz1
  real ( kind = kind_real ) dz2
  real ( kind = kind_real ) dz3
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) listp(lmax)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp1
  integer ( kind = 4 ) lp2
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lptrp(lmax)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nearnd
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nst
  real ( kind = kind_real ) p(3)
  real ( kind = kind_real ) x(n)
  real ( kind = kind_real ) y(n)
  real ( kind = kind_real ) z(n)

  nearnd = 0
  al = 0.0_kind_real
!
!  Store local parameters and test for N invalid.
!
  nn = n

  if ( nn < 3 ) then
    return
  end if

  nst = ist

  if ( nst < 1 .or. nn < nst ) then
    nst = 1
  end if
!
!  Find a triangle (I1,I2,I3) containing P, or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from P.
!
  call trfind ( nst, p, n, x, y, z, list, lptr, lend, b1, b2, b3, i1, i2, i3 )
!
!  Test for collinear nodes.
!
  if ( i1 == 0 ) then
    return
  end if
!
!  Store the linked list of 'neighbors' of P in LISTP and
!  LPTRP.  I1 is the first neighbor, and 0 is stored as
!  the last neighbor if P is not contained in a triangle.
!  L is the length of LISTP and LPTRP, and is limited to
!  LMAX.
!
  if ( i3 /= 0 ) then

    listp(1) = i1
    lptrp(1) = 2
    listp(2) = i2
    lptrp(2) = 3
    listp(3) = i3
    lptrp(3) = 1
    l = 3

  else

    n1 = i1
    l = 1
    lp1 = 2
    listp(l) = n1
    lptrp(l) = lp1
!
!  Loop on the ordered sequence of visible boundary nodes
!  N1 from I1 to I2.
!
    do

      lpl = lend(n1)
      n1 = -list(lpl)
      l = lp1
      lp1 = l+1
      listp(l) = n1
      lptrp(l) = lp1

      if ( n1 == i2 .or. lmax <= lp1 ) then
        exit
      end if

    end do

    l = lp1
    listp(l) = 0
    lptrp(l) = 1

  end if
!
!  Initialize variables for a loop on arcs N1-N2 opposite P
!  in which new 'neighbors' are 'swapped' in.  N1 follows
!  N2 as a neighbor of P, and LP1 and LP2 are the LISTP
!  indexes of N1 and N2.
!
  lp2 = 1
  n2 = i1
  lp1 = lptrp(1)
  n1 = listp(lp1)
!
!  Begin loop:  find the node N3 opposite N1->N2.
!
  do

    lp = lstptr ( lend(n1), n2, list, lptr )

    if ( 0 <= list(lp) ) then

      lp = lptr(lp)
      n3 = abs ( list(lp) )
!
!  Swap test:  Exit the loop if L = LMAX.
!
      if ( l == lmax ) then
        exit
      end if

      dx1 = x(n1) - p(1)
      dy1 = y(n1) - p(2)
      dz1 = z(n1) - p(3)

      dx2 = x(n2) - p(1)
      dy2 = y(n2) - p(2)
      dz2 = z(n2) - p(3)

      dx3 = x(n3) - p(1)
      dy3 = y(n3) - p(2)
      dz3 = z(n3) - p(3)
!
!  Swap:  Insert N3 following N2 in the adjacency list for P.
!  The two new arcs opposite P must be tested.
!
      if ( dx3 * ( dy2 * dz1 - dy1 * dz2 ) - &
           dy3 * ( dx2 * dz1 - dx1 * dz2 ) + &
           dz3 * ( dx2 * dy1 - dx1 * dy2 ) > 0.0_kind_real ) then

        l = l+1
        lptrp(lp2) = l
        listp(l) = n3
        lptrp(l) = lp1
        lp1 = l
        n1 = n3
        cycle

      end if

    end if
!
!  No swap:  Advance to the next arc and test for termination
!  on N1 = I1 (LP1 = 1) or N1 followed by 0.
!
    if ( lp1 == 1 ) then
      exit
    end if

    lp2 = lp1
    n2 = n1
    lp1 = lptrp(lp1)
    n1 = listp(lp1)

    if ( n1 == 0 ) then
      exit
    end if

  end do
!
!  Set NR and DSR to the index of the nearest node to P and
!  an increasing function (negative cosine) of its distance
!  from P, respectively.
!
  nr = i1
  dsr = -( x(nr) * p(1) + y(nr) * p(2) + z(nr) * p(3) )

  do lp = 2, l

    n1 = listp(lp)

    if ( n1 == 0 ) then
      cycle
    end if

    ds1 = -( x(n1) * p(1) + y(n1) * p(2) + z(n1) * p(3) )

    if ( ds1 < dsr ) then
      nr = n1
      dsr = ds1
    end if

  end do

  dsr = -dsr
  dsr = min ( dsr, 1.0_kind_real )

  al = acos ( dsr )
  nearnd = nr

  return
end
subroutine swap ( in1, in2, io1, io2, list, lptr, lend, lp21 )

!*****************************************************************************80
!
!! SWAP replaces the diagonal arc of a quadrilateral with the other diagonal.
!
!  Discussion:
!
!    Given a triangulation of a set of points on the unit
!    sphere, this subroutine replaces a diagonal arc in a
!    strictly convex quadrilateral (defined by a pair of adja-
!    cent triangles) with the other diagonal.  Equivalently, a
!    pair of adjacent triangles is replaced by another pair
!    having the same union.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IN1, IN2, IO1, IO2, nodal indexes of the 
!    vertices of the quadrilateral.  IO1-IO2 is replaced by IN1-IN2.  
!    (IO1,IO2,IN1) and (IO2,IO1,IN2) must be triangles on input.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    the data structure defining the triangulation, created by TRMESH.  
!    On output, updated with the swap; triangles (IO1,IO2,IN1) an (IO2,IO1,IN2) 
!    are replaced by (IN1,IN2,IO2) and (IN2,IN1,IO1) unless LP21 = 0.
!
!    Output, integer ( kind = 4 ) LP21, index of IN1 as a neighbor of IN2 after
!    the swap is performed unless IN1 and IN2 are adjacent on input, in which 
!    case LP21 = 0.
!
!  Local parameters:
!
!    LP, LPH, LPSAV = LIST pointers
!
  implicit none

  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpsav
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
!
!  Test for IN1 and IN2 adjacent.
!
  lp = lstptr ( lend(in1), in2, list, lptr )

  if ( abs ( list(lp) ) == in2 ) then
    lp21 = 0
    return
  end if
!
!  Delete IO2 as a neighbor of IO1.
!
  lp = lstptr ( lend(io1), in2, list, lptr )
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO2 is the last neighbor of IO1, make IN2 the last neighbor.
!
  if ( lend(io1) == lph ) then
    lend(io1) = lp
  end if
!
!  Insert IN2 as a neighbor of IN1 following IO1 using the hole created above.
!
  lp = lstptr ( lend(in1), io1, list, lptr )
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in2
  lptr(lph) = lpsav
!
!  Delete IO1 as a neighbor of IO2.
!
  lp = lstptr ( lend(io2), in1, list, lptr )
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO1 is the last neighbor of IO2, make IN1 the last neighbor.
!
  if ( lend(io2) == lph ) then
    lend(io2) = lp
  end if
!
!  Insert IN1 as a neighbor of IN2 following IO2.
!
  lp = lstptr ( lend(in2), io2, list, lptr )
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in1
  lptr(lph) = lpsav
  lp21 = lph

  return
end
subroutine trans ( n, rlat, rlon, x, y, z )

!*****************************************************************************80
!
!! TRANS transforms spherical coordinates to Cartesian coordinates.
!
!  Discussion:
!
!    This subroutine transforms spherical coordinates into
!    Cartesian coordinates on the unit sphere for input to
!    TRMESH.  Storage for X and Y may coincide with
!    storage for RLAT and RLON if the latter need not be saved.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes (points on the unit 
!    sphere) whose coordinates are to be transformed.
!
!    Input, real ( kind = kind_real ) RLAT(N), latitudes of the nodes in radians.
!
!    Input, real ( kind = kind_real ) RLON(N), longitudes of the nodes in radians.
!
!    Output, real ( kind = kind_real ) X(N), Y(N), Z(N), the coordinates in the 
!    range -1 to 1.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1 for I = 1 to N.
!
!  Local parameters:
!
!    COSPHI = cos(PHI)
!    I =      DO-loop index
!    NN =     Local copy of N
!    PHI =    Latitude
!    THETA =  Longitude
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = kind_real ) cosphi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nn
  real ( kind = kind_real ) phi
  real ( kind = kind_real ) rlat(n)
  real ( kind = kind_real ) rlon(n)
  real ( kind = kind_real ) theta
  real ( kind = kind_real ) x(n)
  real ( kind = kind_real ) y(n)
  real ( kind = kind_real ) z(n)

  nn = n

  do i = 1, nn
    phi = rlat(i)
    theta = rlon(i)
    cosphi = cos ( phi )
    x(i) = cosphi * cos ( theta )
    y(i) = cosphi * sin ( theta )
    z(i) = sin ( phi )
  end do

  return
end
subroutine trfind ( nst, p, n, x, y, z, list, lptr, lend, b1, b2, b3, i1, &
  i2, i3 )

!*****************************************************************************80
!
!! TRFIND locates a point relative to a triangulation.
!
!  Discussion:
!
!    This subroutine locates a point P relative to a triangulation
!    created by TRMESH.  If P is contained in
!    a triangle, the three vertex indexes and barycentric 
!    coordinates are returned.  Otherwise, the indexes of the
!    visible boundary nodes are returned.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, index of a node at which TRFIND begins
!    its search.  Search time depends on the proximity of this node to P.
!
!    Input, real ( kind = kind_real ) P(3), the x, y, and z coordinates (in that order)
!    of the point P to be located.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, real ( kind = kind_real ) X(N), Y(N), Z(N), the coordinates of the
!    triangulation nodes (unit vectors).
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), the 
!    data structure defining the triangulation, created by TRMESH.
!
!    Output, real ( kind = kind_real ) B1, B2, B3, the unnormalized barycentric
!    coordinates of the central projection of P onto the underlying planar
!    triangle if P is in the convex hull of the nodes.  These parameters 
!    are not altered if I1 = 0.
!
!    Output, integer ( kind = 4 ) I1, I2, I3, the counterclockwise-ordered 
!    vertex indexes of a triangle containing P if P is contained in a triangle.
!    If P is not in the convex hull of the nodes, I1 and I2 are the rightmost 
!    and leftmost (boundary) nodes that are visible from P, and I3 = 0.  (If 
!    all boundary nodes are visible from P, then I1 and I2 coincide.)
!    I1 = I2 = I3 = 0 if P and all of the nodes are coplanar (lie on a 
!    common great circle.
!
!  Local parameters:
!
!    EPS =      Machine precision
!    IX,IY,IZ = Integer seeds for JRAND
!    LP =       LIST pointer
!    N0,N1,N2 = Nodes in counterclockwise order defining a
!               cone (with vertex N0) containing P, or end-
!               points of a boundary edge such that P Right
!               N1->N2
!    N1S,N2S =  Initially-determined values of N1 and N2
!    N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
!    NEXT =     Candidate for I1 or I2 when P is exterior
!    NF,NL =    First and last neighbors of N0, or first
!               (rightmost) and last (leftmost) nodes
!               visible from P when P is exterior to the
!               triangulation
!    PTN1 =     Scalar product <P,N1>
!    PTN2 =     Scalar product <P,N2>
!    Q =        (N2 X N1) X N2  or  N1 X (N2 X N1) -- used in
!               the boundary traversal when P is exterior
!    S12 =      Scalar product <N1,N2>
!    TOL =      Tolerance (multiple of EPS) defining an upper
!               bound on the magnitude of a negative bary-
!               centric coordinate (B1 or B2) for P in a
!               triangle -- used to avoid an infinite number
!               of restarts with 0 <= B3 < EPS and B1 < 0 or
!               B2 < 0 but small in magnitude
!    XP,YP,ZP = Local variables containing P(1), P(2), and P(3)
!    X0,Y0,Z0 = Dummy arguments for DET
!    X1,Y1,Z1 = Dummy arguments for DET
!    X2,Y2,Z2 = Dummy arguments for DET
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = kind_real ) b1
  real ( kind = kind_real ) b2
  real ( kind = kind_real ) b3
  real ( kind = kind_real ) det
  real ( kind = kind_real ) eps
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ), save :: ix = 1
  integer ( kind = 4 ), save :: iy = 2
  integer ( kind = 4 ), save :: iz = 3
  integer ( kind = 4 ) jrand
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1s
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2s
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nst
  real ( kind = kind_real ) p(3)
  real ( kind = kind_real ) ptn1
  real ( kind = kind_real ) ptn2
  real ( kind = kind_real ) q(3)
  real ( kind = kind_real ) s12
  real ( kind = kind_real ) tol
  real ( kind = kind_real ) x(n)
  real ( kind = kind_real ) xp
  real ( kind = kind_real ) y(n)
  real ( kind = kind_real ) yp
  real ( kind = kind_real ) z(n)
  real ( kind = kind_real ) zp
!
!  Initialize variables.
!
  xp = p(1)
  yp = p(2)
  zp = p(3)
  n0 = nst

  if ( n0 < 1 .or. n < n0 ) then
    n0 = jrand ( n, ix, iy, iz )
  end if
!
!  Compute the relative machine precision EPS and TOL.
!
  eps = epsilon ( eps )
  tol = 100.0_kind_real * eps
!
!  Set NF and NL to the first and last neighbors of N0, and initialize N1 = NF.
!
2 continue

  lp = lend(n0)
  nl = list(lp)
  lp = lptr(lp)
  nf = list(lp)
  n1 = nf
!
!  Find a pair of adjacent neighbors N1,N2 of N0 that define
!  a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
!
  if ( 0 < nl ) then
!
!  N0 is an interior node.  Find N1.
!
3   continue

    if ( det ( x(n0),y(n0),z(n0),x(n1),y(n1),z(n1),xp,yp,zp ) < 0.0_kind_real ) then
      lp = lptr(lp)
      n1 = list(lp)
      if ( n1 == nl ) then
        go to 6
      end if
      go to 3
    end if

  else
!
!  N0 is a boundary node.  Test for P exterior.
!
    nl = -nl
!
!  Is P to the right of the boundary edge N0->NF?
!
    if ( det(x(n0),y(n0),z(n0),x(nf),y(nf),z(nf), xp,yp,zp) < 0.0_kind_real ) then
      n1 = n0
      n2 = nf
      go to 9
    end if
!
!  Is P to the right of the boundary edge NL->N0?
!
    if ( det(x(nl),y(nl),z(nl),x(n0),y(n0),z(n0),xp,yp,zp) < 0.0_kind_real ) then
      n1 = nl
      n2 = n0
      go to 9
    end if

  end if
!
!  P is to the left of arcs N0->N1 and NL->N0.  Set N2 to the
!  next neighbor of N0 (following N1).
!
4 continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    if ( det(x(n0),y(n0),z(n0),x(n2),y(n2),z(n2),xp,yp,zp) < 0.0_kind_real ) then
      go to 7
    end if

    n1 = n2

    if ( n1 /= nl ) then
      go to 4
    end if

  if ( det ( x(n0), y(n0), z(n0), x(nf), y(nf), z(nf), xp, yp, zp ) &
    < 0.0_kind_real ) then
    go to 6
  end if
!
!  P is left of or on arcs N0->NB for all neighbors NB
!  of N0.  Test for P = +/-N0.
!
  if ( abs ( x(n0 ) * xp + y(n0) * yp + z(n0) * zp) < 1.0_kind_real - 4.0_kind_real * eps ) then
!
!  All points are collinear iff P Left NB->N0 for all
!  neighbors NB of N0.  Search the neighbors of N0.
!  Note:  N1 = NL and LP points to NL.
!
    do

      if ( det(x(n1),y(n1),z(n1),x(n0),y(n0),z(n0),xp,yp,zp) < 0.0_kind_real ) then
        exit
      end if

      lp = lptr(lp)
      n1 = abs ( list(lp) )

      if ( n1 == nl ) then
        i1 = 0
        i2 = 0
        i3 = 0
        return
      end if

    end do

  end if
!
!  P is to the right of N1->N0, or P = +/-N0.  Set N0 to N1 and start over.
!
  n0 = n1
  go to 2
!
!  P is between arcs N0->N1 and N0->NF.
!
6 continue

  n2 = nf
!
!  P is contained in a wedge defined by geodesics N0-N1 and
!  N0-N2, where N1 is adjacent to N2.  Save N1 and N2 to
!  test for cycling.
!
7 continue

  n3 = n0
  n1s = n1
  n2s = n2
!
!  Top of edge-hopping loop:
!
8 continue

  b3 = det ( x(n1),y(n1),z(n1),x(n2),y(n2),z(n2),xp,yp,zp )

  if ( b3 < 0.0_kind_real ) then
!
!  Set N4 to the first neighbor of N2 following N1 (the
!  node opposite N2->N1) unless N1->N2 is a boundary arc.
!
    lp = lstptr ( lend(n2), n1, list, lptr )

    if ( list(lp) < 0 ) then
      go to 9
    end if

    lp = lptr(lp)
    n4 = abs ( list(lp) )
!
!  Define a new arc N1->N2 which intersects the geodesic N0-P.
!
    if ( det ( x(n0),y(n0),z(n0),x(n4),y(n4),z(n4),xp,yp,zp ) < 0.0_kind_real ) then
      n3 = n2
      n2 = n4
      n1s = n1
      if ( n2 /= n2s .and. n2 /= n0 ) then
        go to 8
      end if
    else
      n3 = n1
      n1 = n4
      n2s = n2
      if ( n1 /= n1s .and. n1 /= n0 ) then
        go to 8
      end if
    end if
!
!  The starting node N0 or edge N1-N2 was encountered
!  again, implying a cycle (infinite loop).  Restart
!  with N0 randomly selected.
!
    n0 = jrand ( n, ix, iy, iz )
    go to 2

  end if
!
!  P is in (N1,N2,N3) unless N0, N1, N2, and P are collinear
!  or P is close to -N0.
!
  if ( b3 >= eps ) then
!
!  B3 /= 0.
!
    b1 = det(x(n2),y(n2),z(n2),x(n3),y(n3),z(n3),xp,yp,zp)
    b2 = det(x(n3),y(n3),z(n3),x(n1),y(n1),z(n1),xp,yp,zp)
!
!  Restart with N0 randomly selected.
!
    if ( b1 < -tol .or. b2 < -tol ) then
      n0 = jrand ( n, ix, iy, iz )
      go to 2
    end if

  else
!
!  B3 = 0 and thus P lies on N1->N2. Compute
!  B1 = Det(P,N2 X N1,N2) and B2 = Det(P,N1,N2 X N1).
!
    b3 = 0.0_kind_real
    s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)
    ptn1 = xp * x(n1) + yp * y(n1) + zp * z(n1)
    ptn2 = xp * x(n2) + yp * y(n2) + zp * z(n2)
    b1 = ptn1 - s12 * ptn2
    b2 = ptn2 - s12 * ptn1
!
!  Restart with N0 randomly selected.
!
    if ( b1 < -tol .or. b2 < -tol ) then
      n0 = jrand ( n, ix, iy, iz )
      go to 2
    end if

  end if
!
!  P is in (N1,N2,N3).
!
  i1 = n1
  i2 = n2
  i3 = n3
  b1 = max ( b1, 0.0_kind_real )
  b2 = max ( b2, 0.0_kind_real )
  return
!
!  P Right N1->N2, where N1->N2 is a boundary edge.
!  Save N1 and N2, and set NL = 0 to indicate that
!  NL has not yet been found.
!
9 continue

  n1s = n1
  n2s = n2
  nl = 0
!
!  Counterclockwise Boundary Traversal:
!
10 continue

  lp = lend(n2)
  lp = lptr(lp)
  next = list(lp)

  if ( det(x(n2),y(n2),z(n2),x(next),y(next),z(next),xp,yp,zp) >= 0.0_kind_real ) then
!
!  N2 is the rightmost visible node if P Forward N2->N1
!  or NEXT Forward N2->N1.  Set Q to (N2 X N1) X N2.
!
    s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)

    q(1) = x(n1) - s12 * x(n2)
    q(2) = y(n1) - s12 * y(n2)
    q(3) = z(n1) - s12 * z(n2)

    if ( xp * q(1) + yp * q(2) + zp * q(3) >= 0.0_kind_real ) then
      go to 11
    end if

    if ( x(next) * q(1) + y(next) * q(2) + z(next) * q(3) >= 0.0_kind_real ) then
      go to 11
    end if
!
!  N1, N2, NEXT, and P are nearly collinear, and N2 is
!  the leftmost visible node.
!
    nl = n2
  end if
!
!  Bottom of counterclockwise loop:
!
  n1 = n2
  n2 = next

  if ( n2 /= n1s ) then
    go to 10
  end if
!
!  All boundary nodes are visible from P.
!
  i1 = n1s
  i2 = n1s
  i3 = 0
  return
!
!  N2 is the rightmost visible node.
!
11 continue

  nf = n2

  if ( nl == 0 ) then
!
!  Restore initial values of N1 and N2, and begin the search
!  for the leftmost visible node.
!
    n2 = n2s
    n1 = n1s
!
!  Clockwise Boundary Traversal:
!
12  continue

    lp = lend(n1)
    next = -list(lp)

    if ( 0.0_kind_real <= &
      det ( x(next), y(next), z(next), x(n1), y(n1), z(n1), xp, yp, zp )  ) then
!
!  N1 is the leftmost visible node if P or NEXT is
!  forward of N1->N2.  Compute Q = N1 X (N2 X N1).
!
      s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)
      q(1) = x(n2) - s12 * x(n1)
      q(2) = y(n2) - s12 * y(n1)
      q(3) = z(n2) - s12 * z(n1)

      if ( xp * q(1) + yp * q(2) + zp * q(3) >= 0.0_kind_real ) then
        go to 13
      end if

      if ( x(next) * q(1) + y(next) * q(2) + z(next) * q(3) >= 0.0_kind_real ) then
        go to 13
      end if
!
!  P, NEXT, N1, and N2 are nearly collinear and N1 is the rightmost 
!  visible node.
!
      nf = n1
    end if
!
!  Bottom of clockwise loop:
!
    n2 = n1
    n1 = next

    if ( n1 /= n1s ) then
      go to 12
    end if
!
!  All boundary nodes are visible from P.
!
    i1 = n1
    i2 = n1
    i3 = 0
    return
!
!  N1 is the leftmost visible node.
!
13   continue

    nl = n1

  end if
!
!  NF and NL have been found.
!
  i1 = nf
  i2 = nl
  i3 = 0

  return
end
subroutine trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )

!*****************************************************************************80
!
!! TRLIST converts a triangulation data structure to a triangle list.
!
!  Discussion:
!
!    This subroutine converts a triangulation data structure
!    from the linked list created by TRMESH to a triangle list.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), linked
!    list data structure defining the triangulation.  Refer to TRMESH.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows (entries per triangle)
!    reserved for the triangle list LTRI.  The value must be 6 if only the 
!    vertex indexes and neighboring triangle indexes are to be stored, or 9
!    if arc indexes are also to be assigned and stored.  Refer to LTRI.
!
!    Output, integer ( kind = 4 ) NT, the number of triangles in the 
!    triangulation unless IER /=0, in which case NT = 0.  NT = 2N-NB-2 if 
!    NB >= 3 or 2N-4 if NB = 0, where NB is the number of boundary nodes.
!
!    Output, integer ( kind = 4 ) LTRI(NROW,*).  The second dimension of LTRI
!    must be at least NT, where NT will be at most 2*N-4.  The J-th column 
!    contains the vertex nodal indexes (first three rows), neighboring triangle 
!    indexes (second three rows), and, if NROW = 9, arc indexes (last three
!    rows) associated with triangle J for J = 1,...,NT.  The vertices are
!    ordered counterclockwise with the first vertex taken to be the one 
!    with smallest index.  Thus, LTRI(2,J) and LTRI(3,J) are larger than
!    LTRI(1,J) and index adjacent neighbors of node LTRI(1,J).  For 
!    I = 1,2,3, LTRI(I+3,J) and LTRI(I+6,J) index the triangle and arc,
!    respectively, which are opposite (not shared by) node LTRI(I,J), with
!    LTRI(I+3,J) = 0 if LTRI(I+6,J) indexes a boundary arc.  Vertex indexes
!    range from 1 to N, triangle indexes from 0 to NT, and, if included, 
!    arc indexes from 1 to NA, where NA = 3N-NB-3 if NB >= 3 or 3N-6 if 
!    NB = 0.  The triangles are ordered on first (smallest) vertex indexes.
!
!    Output, integer ( kind = 4 ) IER, error indicator.
!    0, if no errors were encountered.
!    1, if N or NROW is outside its valid range on input.
!    2, if the triangulation data structure (LIST,LPTR,LEND) is invalid.  
!      Note, however, that these arrays are not completely tested for validity.
!
!  Local parameters:
!
!    ARCS =     Logical variable with value TRUE iff are
!               indexes are to be stored
!    I,J =      LTRI row indexes (1 to 3) associated with
!               triangles KT and KN, respectively
!    I1,I2,I3 = Nodal indexes of triangle KN
!    ISV =      Variable used to permute indexes I1,I2,I3
!    KA =       Arc index and number of currently stored arcs
!    KN =       Index of the triangle that shares arc I1-I2 with KT
!    KT =       Triangle index and number of currently stored triangles
!    LP =       LIST pointer
!    LP2 =      Pointer to N2 as a neighbor of N1
!    LPL =      Pointer to the last neighbor of I1
!    LPLN1 =    Pointer to the last neighbor of N1
!    N1,N2,N3 = Nodal indexes of triangle KT
!    NM2 =      N-2
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrow

  logical              arcs
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) isv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kn
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp2
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpln1
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) ltri(nrow,*)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) nm2
  integer ( kind = 4 ) nt
!
!  Test for invalid input parameters.
!
  if ( n < 3 .or. ( nrow /= 6 .and. nrow /= 9 ) ) then
    nt = 0
    ier = 1
    return
  end if
!
!  Initialize parameters for loop on triangles KT = (N1,N2,
!  N3), where N1 < N2 and N1 < N3.
!
!  ARCS = TRUE iff arc indexes are to be stored.
!  KA,KT = Numbers of currently stored arcs and triangles.
!  NM2 = Upper bound on candidates for N1.
!
  arcs = nrow == 9
  ka = 0
  kt = 0
  nm2 = n-2
!
!  Loop on nodes N1.
!
  do n1 = 1, nm2
!
!  Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!  to the last neighbor of N1, and LP2 points to N2.
!
    lpln1 = lend(n1)
    lp2 = lpln1

1   continue

      lp2 = lptr(lp2)
      n2 = list(lp2)
      lp = lptr(lp2)
      n3 = abs ( list(lp) )

      if ( n2 < n1 .or. n3 < n1 ) then
        go to 8
      end if
!
!  Add a new triangle KT = (N1,N2,N3).
!
      kt = kt + 1
      ltri(1,kt) = n1
      ltri(2,kt) = n2
      ltri(3,kt) = n3
!
!  Loop on triangle sides (I2,I1) with neighboring triangles
!  KN = (I1,I2,I3).
!
      do i = 1, 3

        if ( i == 1 ) then
          i1 = n3
          i2 = n2
        else if ( i == 2 ) then
          i1 = n1
          i2 = n3
        else
          i1 = n2
          i2 = n1
        end if
        j = 0
!
!  Set I3 to the neighbor of I1 that follows I2 unless
!  I2->I1 is a boundary arc.
!
        lpl = lend(i1)
        lp = lptr(lpl)

        do

          if ( list(lp) == i2 ) then
            go to 3
          end if

          lp = lptr(lp)

          if ( lp == lpl ) then
            exit
          end if

        end do
!
!  Invalid triangulation data structure:  I1 is a neighbor of
!  I2, but I2 is not a neighbor of I1.
!
        if ( abs ( list(lp) ) /= i2 ) then
          nt = 0
          ier = 2
          return
        end if
!
!  I2 is the last neighbor of I1.  Bypass the search for a neighboring
!  triangle if I2->I1 is a boundary arc.
!
        kn = 0

        if ( list(lp) < 0 ) then
          go to 6
        end if
!
!  I2->I1 is not a boundary arc, and LP points to I2 as
!  a neighbor of I1.
!
3       continue

        lp = lptr(lp)
        i3 = abs ( list(lp) )
!
!  Find J such that LTRI(J,KN) = I3 (not used if KT < KN),
!  and permute the vertex indexes of KN so that I1 is smallest.
!
        if ( i1 < i2 .and. i1 < i3 ) then
          j = 3
        else if ( i2 < i3 ) then
          j = 2
          isv = i1
          i1 = i2
          i2 = i3
          i3 = isv
        else
          j = 1
          isv = i1
          i1 = i3
          i3 = i2
          i2 = isv
        end if
!
!  Test for KT < KN (triangle index not yet assigned).
!
        if ( n1 < i1 ) then
          cycle
        end if
!
!  Find KN, if it exists, by searching the triangle list in
!  reverse order.
!
        do kn = kt-1, 1, -1
          if ( ltri(1,kn) == i1 .and. &
               ltri(2,kn) == i2 .and. &
               ltri(3,kn) == i3 ) then
            go to 5
          end if
        end do

        cycle
!
!  Store KT as a neighbor of KN.
!
5       continue

        ltri(j+3,kn) = kt
!
!  Store KN as a neighbor of KT, and add a new arc KA.
!
6       continue

        ltri(i+3,kt) = kn

        if ( arcs ) then
          ka = ka + 1
          ltri(i+6,kt) = ka
          if ( kn /= 0 ) then
            ltri(j+6,kn) = ka
          end if
        end if

        end do
!
!  Bottom of loop on triangles.
!
8     continue

      if ( lp2 /= lpln1 ) then
        go to 1
      end if

  end do

  nt = kt
  ier = 0

  return
end
subroutine trmesh ( n, x, y, z, list, lptr, lend, lnew, near, next, dist, ier )

!*****************************************************************************80
!
!! TRMESH creates a Delaunay triangulation on the unit sphere.
!
!  Discussion:
!
!    This subroutine creates a Delaunay triangulation of a
!    set of N arbitrarily distributed points, referred to as
!    nodes, on the surface of the unit sphere.  The Delaunay
!    triangulation is defined as a set of (spherical) triangles
!    with the following five properties:
!
!     1)  The triangle vertices are nodes.
!     2)  No triangle contains a node other than its vertices.
!     3)  The interiors of the triangles are pairwise disjoint.
!     4)  The union of triangles is the convex hull of the set
!           of nodes (the smallest convex set that contains
!           the nodes).  If the nodes are not contained in a
!           single hemisphere, their convex hull is the 
!           entire sphere and there are no boundary nodes.
!           Otherwise, there are at least three boundary nodes.
!     5)  The interior of the circumcircle of each triangle
!           contains no node.
!
!    The first four properties define a triangulation, and the
!    last property results in a triangulation which is as close
!    as possible to equiangular in a certain sense and which is
!    uniquely defined unless four or more nodes lie in a common
!    plane.  This property makes the triangulation well-suited
!    for solving closest-point problems and for triangle-based
!    interpolation.
!
!    Provided the nodes are randomly ordered, the algorithm
!    has expected time complexity O(N*log(N)) for most nodal
!    distributions.  Note, however, that the complexity may be
!    as high as O(N**2) if, for example, the nodes are ordered
!    on increasing latitude.
!
!    Spherical coordinates (latitude and longitude) may be
!    converted to Cartesian coordinates by Subroutine TRANS.
!
!    The following is a list of the software package modules
!    which a user may wish to call directly:
!
!    ADDNOD - Updates the triangulation by appending a new node.
!
!    AREAS  - Returns the area of a spherical triangle.
!
!    BNODES - Returns an array containing the indexes of the
!             boundary nodes (if any) in counterclockwise
!             order.  Counts of boundary nodes, triangles,
!             and arcs are also returned.
!
!    CIRCUM - Returns the circumcenter of a spherical triangle.
!
!    CRLIST - Returns the set of triangle circumcenters
!             (Voronoi vertices) and circumradii associated
!             with a triangulation.
!
!    DELARC - Deletes a boundary arc from a triangulation.
!
!    DELNOD - Updates the triangulation with a nodal deletion.
!
!    EDGE   - Forces an arbitrary pair of nodes to be connected
!             by an arc in the triangulation.
!
!    GETNP  - Determines the ordered sequence of L closest nodes
!             to a given node, along with the associated distances.
!
!    INSIDE - Locates a point relative to a polygon on the
!             surface of the sphere.
!
!    INTRSC - Returns the point of intersection between a
!             pair of great circle arcs.
!
!    JRAND  - Generates a uniformly distributed pseudo-random integer.
!
!    LEFT   - Locates a point relative to a great circle.
!
!    NEARND - Returns the index of the nearest node to an
!             arbitrary point, along with its squared
!             distance.
!
!    SCOORD - Converts a point from Cartesian coordinates to
!             spherical coordinates.
!
!    TRANS  - Transforms spherical coordinates into Cartesian
!             coordinates on the unit sphere for input to
!             Subroutine TRMESH.
!
!    TRLIST - Converts the triangulation data structure to a
!             triangle list more suitable for use in a finite
!             element code.
!
!    TRMESH - Creates a Delaunay triangulation of a set of
!             nodes.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, real ( kind = kind_real ) X(N), Y(N), Z(N), the coordinates of distinct
!    nodes.  (X(K),Y(K), Z(K)) is referred to as node K, and K is referred 
!    to as a nodal index.  It is required that X(K)**2 + Y(K)**2 + Z(K)**2 = 1
!    for all K.  The first three nodes must not be collinear (lie on a 
!    common great circle).
!
!    Output, integer ( kind = 4 ) LIST(6*(N-2)), nodal indexes which, along 
!    with LPTR, LEND, and LNEW, define the triangulation as a set of N 
!    adjacency lists; counterclockwise-ordered sequences of neighboring nodes 
!    such that the first and last neighbors of a boundary node are boundary 
!    nodes (the first neighbor of an interior node is arbitrary).  In order to 
!    distinguish between interior and boundary nodes, the last neighbor of 
!    each boundary node is represented by the negative of its index.
!
!    Output, integer ( kind = 4 ) LPTR(6*(N-2)), = Set of pointers (LIST 
!    indexes) in one-to-one correspondence with the elements of LIST.
!    LIST(LPTR(I)) indexes the node which follows LIST(I) in cyclical
!    counterclockwise order (the first neighbor follows the last neighbor).
!
!    Output, integer ( kind = 4 ) LEND(N), pointers to adjacency lists.  
!    LEND(K) points to the last neighbor of node K.  LIST(LEND(K)) < 0 if and 
!    only if K is a boundary node.
!
!    Output, integer ( kind = 4 ) LNEW, pointer to the first empty location 
!    in LIST and LPTR (list length plus one).  LIST, LPTR, LEND, and LNEW are 
!    not altered if IER < 0, and are incomplete if 0 < IER.
!
!    Workspace, integer ( kind = 4 ) NEAR(N), 
!    used to efficiently determine the nearest triangulation node to each
!    unprocessed node for use by ADDNOD.
!
!    Workspace, integer ( kind = 4 ) NEXT(N),
!    used to efficiently determine the nearest triangulation node to each
!    unprocessed node for use by ADDNOD.
!
!    Workspace, real ( kind = kind_real ) DIST(N), 
!    used to efficiently determine the nearest triangulation node to each
!    unprocessed node for use by ADDNOD.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!     0, if no errors were encountered.
!    -1, if N < 3 on input.
!    -2, if the first three nodes are collinear.
!     L, if nodes L and M coincide for some L < M.  The data structure 
!      represents a triangulation of nodes 1 to M-1 in this case.
!
!  Local parameters:
!
!    D =        (Negative cosine of) distance from node K to node I
!    D1,D2,D3 = Distances from node K to nodes 1, 2, and 3, respectively
!    I,J =      Nodal indexes
!    I0 =       Index of the node preceding I in a sequence of
!               unprocessed nodes:  I = NEXT(I0)
!    K =        Index of node to be added and DO-loop index: 3 < K
!    LP =       LIST index (pointer) of a neighbor of K
!    LPL =      Pointer to the last neighbor of K
!    NEXTI =    NEXT(I)
!    NN =       Local copy of N
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = kind_real ) d
  real ( kind = kind_real ) d1
  real ( kind = kind_real ) d2
  real ( kind = kind_real ) d3
  real ( kind = kind_real ) dist(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) near(n)
  integer ( kind = 4 ) next(n)
  integer ( kind = 4 ) nexti
  integer ( kind = 4 ) nn
  real ( kind = kind_real ) x(n)
  real ( kind = kind_real ) y(n)
  real ( kind = kind_real ) z(n)

  nn = n

  if ( nn < 3 ) then
    ier = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRMESH - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop
  end if
!
!  Store the first triangle in the linked list.
!
  if ( .not. left (x(1),y(1),z(1),x(2),y(2),z(2), &
                   x(3),y(3),z(3) ) ) then
!
!  The first triangle is (3,2,1) = (2,1,3) = (1,3,2).
!
    list(1) = 3
    lptr(1) = 2
    list(2) = -2
    lptr(2) = 1
    lend(1) = 2

    list(3) = 1
    lptr(3) = 4
    list(4) = -3
    lptr(4) = 3
    lend(2) = 4

    list(5) = 2
    lptr(5) = 6
    list(6) = -1
    lptr(6) = 5
    lend(3) = 6

  else if ( .not. left ( x(2),y(2),z(2),x(1),y(1),z(1),x(3),y(3),z(3) ) ) then
!
!  The first triangle is (1,2,3):  3 Strictly Left 1->2,
!  i.e., node 3 lies in the left hemisphere defined by arc 1->2.
!
    list(1) = 2
    lptr(1) = 2
    list(2) = -3
    lptr(2) = 1
    lend(1) = 2

    list(3) = 3
    lptr(3) = 4
    list(4) = -1
    lptr(4) = 3
    lend(2) = 4

    list(5) = 1
    lptr(5) = 6
    list(6) = -2
    lptr(6) = 5
    lend(3) = 6
!
!  The first three nodes are collinear.
!
  else

    ier = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRMESH - Fatal error!'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    write ( *, '(a)' ) '  Try reordering the data.'
    stop

  end if
!
!  Initialize LNEW and test for N = 3.
!
  lnew = 7

  if ( nn == 3 ) then
    ier = 0
    return
  end if
!
!  A nearest-node data structure (NEAR, NEXT, and DIST) is
!  used to obtain an expected-time (N*log(N)) incremental
!  algorithm by enabling constant search time for locating
!  each new node in the triangulation.
!
!  For each unprocessed node K, NEAR(K) is the index of the
!  triangulation node closest to K (used as the starting
!  point for the search in Subroutine TRFIND) and DIST(K)
!  is an increasing function of the arc length (angular
!  distance) between nodes K and NEAR(K):  -Cos(a) for arc
!  length a.
!
!  Since it is necessary to efficiently find the subset of
!  unprocessed nodes associated with each triangulation
!  node J (those that have J as their NEAR entries), the
!  subsets are stored in NEAR and NEXT as follows:  for
!  each node J in the triangulation, I = NEAR(J) is the
!  first unprocessed node in J's set (with I = 0 if the
!  set is empty), L = NEXT(I) (if 0 < I) is the second,
!  NEXT(L) (if 0 < L) is the third, etc.  The nodes in each
!  set are initially ordered by increasing indexes (which
!  maximizes efficiency) but that ordering is not main-
!  tained as the data structure is updated.
!
!  Initialize the data structure for the single triangle.
!
  near(1) = 0
  near(2) = 0
  near(3) = 0

  do k = nn, 4, -1

    d1 = -( x(k) * x(1) + y(k) * y(1) + z(k) * z(1) )
    d2 = -( x(k) * x(2) + y(k) * y(2) + z(k) * z(2) )
    d3 = -( x(k) * x(3) + y(k) * y(3) + z(k) * z(3) )

    if ( d1 <= d2 .and. d1 <= d3 ) then
      near(k) = 1
      dist(k) = d1
      next(k) = near(1)
      near(1) = k
    else if ( d2 <= d1 .and. d2 <= d3 ) then
      near(k) = 2
      dist(k) = d2
      next(k) = near(2)
      near(2) = k
    else
      near(k) = 3
      dist(k) = d3
      next(k) = near(3)
      near(3) = k
    end if

  end do
!
!  Add the remaining nodes.
!
  do k = 4, nn

    call addnod ( near(k), k, x, y, z, list, lptr, lend, lnew, ier )

    if ( ier /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRMESH - Fatal error!'
      write ( *, '(a,i8)' ) '  ADDNOD returned error code IER = ', ier
      stop
    end if
!
!  Remove K from the set of unprocessed nodes associated with NEAR(K).
!
    i = near(k)
    i0 = i

    if ( near(i) == k ) then

      near(i) = next(k)

    else

      i = near(i)

      do

        i0 = i
        i = next(i0)

        if ( i == k ) then
          exit
        end if

      end do

      next(i0) = next(k)

    end if

    near(k) = 0
!
!  Loop on neighbors J of node K.
!
    lpl = lend(k)
    lp = lpl

3   continue

    lp = lptr(lp)
    j = abs ( list(lp) )
!
!  Loop on elements I in the sequence of unprocessed nodes
!  associated with J:  K is a candidate for replacing J
!  as the nearest triangulation node to I.  The next value
!  of I in the sequence, NEXT(I), must be saved before I
!  is moved because it is altered by adding I to K's set.
!
    i = near(j)

    do

      if ( i == 0 ) then
        exit
      end if

      nexti = next(i)
!
!  Test for the distance from I to K less than the distance
!  from I to J.
!
      d = - ( x(i) * x(k) + y(i) * y(k) + z(i) * z(k) )
      if ( d < dist(i) ) then
!
!  Replace J by K as the nearest triangulation node to I:
!  update NEAR(I) and DIST(I), and remove I from J's set
!  of unprocessed nodes and add it to K's set.
!
        near(i) = k
        dist(i) = d

        if ( i == near(j) ) then
          near(j) = nexti
        else
          next(i0) = nexti
        end if

        next(i) = near(k)
        near(k) = i
      else
        i0 = i
      end if

      i = nexti

    end do
!
!  Bottom of loop on neighbors J.
!
    if ( lp /= lpl ) then
      go to 3
    end if

  end do

  return
end

end module tools_stripack
