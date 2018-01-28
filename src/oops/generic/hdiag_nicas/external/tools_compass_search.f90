!----------------------------------------------------------------------
! Module: tools_compass_search.f90
!> Purpose: compass search minimization routine
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_compass_search

use tools_kinds, only: kind_real
use type_mdata, only: mdatatype
implicit none

private
public :: compass_search

contains

subroutine compass_search ( mdata, func, m, x0, delta_tol, delta_init, &
  k_max, x, fx, k )

!*****************************************************************************80
!
!! COMPASS_SEARCH carries out a direct search minimization algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
!    Optimization by Direct Search: New Perspectives on Some Classical
!    and Modern Methods,
!    SIAM Review,
!    Volume 45, Number 3, 2003, pages 385-482.
!
!  Parameters:
!
!    Input, external real ( kind_real ) FUNCTION_HANDLE, the name of
!    a FORTRAN90 function which evaluates the function to be minimized, of the
!    form FUNCTION FUNCTION_HANDLE ( M, X ).
!
!    Input, integer M, the number of variables.
!
!    Input, real ( kind_real ) X0(M), a starting estimate for the minimizer.
!
!    Input, real ( kind_real ) DELTA_TOL, the smallest step size that is allowed.
!
!    Input, real ( kind_real ) DELTA_INIT, the starting stepsize.
!
!    Input, integer K_MAX, the maximum number of steps allowed.
!
!    Output, real ( kind_real ) X(M), the estimated minimizer.
!
!    Output, real ( kind_real ) FX, the function value at X.
!
!    Output, integer K, the number of steps taken.
!
  implicit none

  type(mdatatype),intent(inout) :: mdata !< Minimization data
  interface
    subroutine func(mdata,x,f)
    use tools_kinds, only: kind_real
    use type_mdata, only: mdatatype
    type(mdatatype),intent(in) :: mdata
    real(kind_real),intent(in) :: x(mdata%nx)
    real(kind_real),intent(out) :: f
    end subroutine
  end interface

  integer m

  logical decrease
  real ( kind_real ) delta
  real ( kind_real ) delta_init
  real ( kind_real ) delta_tol
  real ( kind_real ) fx
  real ( kind_real ) fxd
  integer i
  integer ii
  integer k
  integer k_max
  real ( kind_real ) s
  real ( kind_real ) x(m)
  real ( kind_real ) x0(m)
  real ( kind_real ) xd(m)

  k = 0
  x(1:m) = x0(1:m)
  call func(mdata,x,fx)

  if ( delta_tol <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
    write ( *, '(a)' ) '  DELTA_TOL <= 0.0.'
    write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
    stop
  end if

  if ( delta_init <= delta_tol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
    write ( *, '(a)' ) '  DELTA_INIT < DELTA_TOL.'
    write ( *, '(a,g14.6)' ) '  DELTA_INIT = ', delta_init
    write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
    stop
  end if

  delta = delta_init

  do while ( k < k_max )

    k = k + 1
!
!  For each coordinate direction I, seek a lower function value
!  by increasing or decreasing X(I) by DELTA.
!
    decrease = .false.
    s = + 1.0
    i = 1

    do ii = 1, 2 * m

      xd = x
      xd(i) = xd(i) + s * delta
      call func(mdata,xd,fxd)
!
!  As soon as a decrease is noticed, accept the new point.
!
      if ( fxd < fx ) then
        x = xd
        fx = fxd
        decrease = .true.
        exit
      end if

      s = - s
      if ( .not.(abs(s - 1.0)>0.0) ) then
        i = i + 1
      end if

    end do
!
!  If no decrease occurred, reduce DELTA.
!
    if ( .not. decrease ) then
      delta = delta / 2.0
      if ( delta < delta_tol ) then
        exit
      end if
    end if

  end do

  return
end

end module tools_compass_search

