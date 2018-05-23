!----------------------------------------------------------------------
! Module: tools_test
!> Purpose: test tools
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module tools_test

use tools_const, only: deg2rad
use tools_kinds,only: kind_real
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_kdtree, only: kdtree_type
use type_mpl, only: mpl
use type_nam, only: nam_type
use type_rng, only: rng

implicit none

private
public :: define_dirac,define_test_vectors

contains

!----------------------------------------------------------------------
! Subroutine: define_dirac
!> Purpose: define dirac indices
!----------------------------------------------------------------------
subroutine define_dirac(nam,geom,iprocdir,ic0adir,il0dir)

implicit none

! Passed variables
type(nam_type),intent(in) :: nam          !< Namelist
type(geom_type),intent(in) :: geom        !< Geometry
integer,intent(out) :: iprocdir(nam%ndir) !< Dirac processor index
integer,intent(out) :: ic0adir(nam%ndir)  !< Dirac gridpoint index
integer,intent(out) :: il0dir(nam%ndir)   !< Dirac level index

! Local variables
integer :: idir,il0,nn_index(1),ic0dir
real(kind_real) :: nn_dist(1)

do idir=1,nam%ndir
   ! Find nearest neighbor
   call geom%kdtree%find_nearest_neighbors(nam%londir(idir),nam%latdir(idir),1,nn_index,nn_dist)
   ic0dir = nn_index(1)
   iprocdir(idir) = geom%c0_to_proc(ic0dir)
   ic0adir(idir) = geom%c0_to_c0a(ic0dir)

   ! Find level index
   do il0=1,geom%nl0
      if (nam%levs(il0)==nam%levdir(idir)) il0dir(idir) = il0
   end do
end do

end subroutine define_dirac

!----------------------------------------------------------------------
! Subroutine: define_test_vectors
!> Purpose: define test vectors
!----------------------------------------------------------------------
subroutine define_test_vectors(nam,geom,ntest,fld)

! Passed variables
type(nam_type),intent(in) :: nam                                            !< Namelist
type(geom_type),intent(in) :: geom                                          !< Geometry
integer,intent(in) :: ntest                                                 !< Number of vectors
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest) !< Field

! Local variables
integer :: ic0dir(ntest),il0dir(ntest),ivdir(ntest),itsdir(ntest)
integer :: itest,ic0,iproc,ic0a

! Define dirac locations
if (mpl%main) then
   call rng%rand_integer(1,geom%nc0,ic0dir)
   call rng%rand_integer(1,geom%nl0,il0dir)
end if
call mpl%bcast(ic0dir)
call mpl%bcast(il0dir)
ivdir = 1
itsdir = 1

! Define test vectors
do itest=1,ntest
   fld(:,:,:,:,itest) = 0.0
   ic0 = ic0dir(itest)
   iproc = geom%c0_to_proc(ic0)
   if (iproc==mpl%myproc) then
      ic0a = geom%c0_to_c0a(ic0)
      fld(ic0a,il0dir(itest),ivdir(itest),itsdir(itest),itest) = 1.0
   end if
end do

end subroutine define_test_vectors

end module tools_test
