!----------------------------------------------------------------------
! Module: hdiag_localization.f90
!> Purpose: localization routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_localization

use hdiag_fit, only: compute_fit
use tools_display, only: msgwarning,msgerror,prog_init,prog_print
use tools_fit, only: ver_smooth
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr
use type_avg, only: avgtype
use type_curve, only: curvetype,curve_normalization
use type_hdata, only: hdatatype
use type_mpl, only: mpl

implicit none

interface compute_localization
  module procedure compute_localization
  module procedure compute_localization_local
end interface

private
public :: compute_localization

contains

!----------------------------------------------------------------------
! Subroutine: compute_localization
!> Purpose: compute localization
!----------------------------------------------------------------------
subroutine compute_localization(hdata,ib,avg,loc)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata  !< HDIAG data
integer,intent(in) :: ib             !< Block index
type(avgtype),intent(in) :: avg      !< Averaged statistics
type(curvetype),intent(inout) :: loc !< Localizations

! Local variables
integer :: il0,jl0r,jc3

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Compute raw localization
!$omp parallel do schedule(static) private(il0,jl0r,jc3)
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (isnotmsr(avg%m11asysq(jc3,jl0r,il0)).and.isnotmsr(avg%m11sq(jc3,jl0r,il0))) &
       & loc%raw(jc3,jl0r,il0) = avg%m11asysq(jc3,jl0r,il0)/avg%m11sq(jc3,jl0r,il0)
      end do
   end do
end do
!$omp end parallel do

! Normalize localization
call curve_normalization(hdata,ib,loc)

! Compute localization fits
if (bpar%fit_block(ib)) then
   ! Compute fit weight
   if (nam%fit_wgt) loc%fit_wgt = abs(avg%cor)

   ! Compute initial fit
   call compute_fit(hdata,ib,loc)
end if

! End associate
end associate

end subroutine compute_localization

!----------------------------------------------------------------------
! Subroutine: compute_localization_local
!> Purpose: compute localization, local
!----------------------------------------------------------------------
subroutine compute_localization_local(hdata,ib,avg,loc)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata              !< HDIAG data
integer,intent(in) :: ib                         !< Block index
type(avgtype),intent(in) :: avg(hdata%nc2a)      !< Averaged statistics
type(curvetype),intent(inout) :: loc(hdata%nc2a) !< Localizations

! Local variables
integer :: ic2a,progint
logical :: done(hdata%nc2a)

! Initialization
call prog_init(progint,done)

! Loop over points
do ic2a=1,hdata%nc2a
   ! Compute localization
   call compute_localization(hdata,ib,avg(ic2a),loc(ic2a))

   ! Print progression
   done(ic2a) = .true.
   call prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'

end subroutine compute_localization_local

end module hdiag_localization
