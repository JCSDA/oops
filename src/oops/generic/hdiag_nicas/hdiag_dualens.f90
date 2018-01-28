!----------------------------------------------------------------------
! Module: hdiag_dualens.f90
!> Purpose: dual-ensemble routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_dualens

use hdiag_fit, only: compute_fit
use tools_display, only: msgwarning
use tools_fit, only: ver_smooth
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr
use type_avg, only: avgtype
use type_curve, only: curvetype,curve_normalization
use type_hdata, only: hdatatype
implicit none

private
public :: compute_dualens

contains

!----------------------------------------------------------------------
! Subroutine: compute_dualens
!> Purpose: compute dual-ensemble hybridization formulae
!----------------------------------------------------------------------
subroutine compute_dualens(hdata,ib,avg,avg_lr,loc_deh,loc_deh_lr)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata         !< HDIAG data
integer,intent(in) :: ib                    !< Block index
type(avgtype),intent(in) :: avg             !< Averaged statistics
type(avgtype),intent(in) :: avg_lr          !< Low-resolution averaged statistics
type(curvetype),intent(inout) :: loc_deh    !< Adapted high-resolution localizations
type(curvetype),intent(inout) :: loc_deh_lr !< Adapted low-resolution localizations

! Local variables
integer :: il0,jl0r,jc3
real(kind_real) :: num(hdata%nam%nc3),num_lr(hdata%nam%nc3),den(hdata%nam%nc3)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Compute raw dual-ensemble hybridization
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      do jc3=1,bpar%nc3(ib)
         if (isnotmsr(avg%m11asysq(jc3,jl0r,il0)).and.isnotmsr(avg%m11sq(jc3,jl0r,il0)) &
       & .and.isnotmsr(avg_lr%m11lrm11asy(jc3,jl0r,il0)).and.isnotmsr(avg_lr%m11lrm11(jc3,jl0r,il0)) &
       & .and.isnotmsr(avg_lr%m11sq(jc3,jl0r,il0)).and.isnotmsr(avg_lr%m11lrm11asy(jc3,jl0r,il0))) then
            num(jc3) = avg%m11asysq(jc3,jl0r,il0)*avg_lr%m11sq(jc3,jl0r,il0) &
                    & -avg_lr%m11lrm11asy(jc3,jl0r,il0)*avg_lr%m11lrm11(jc3,jl0r,il0)
            num_lr(jc3) = avg_lr%m11lrm11asy(jc3,jl0r,il0)*avg%m11sq(jc3,jl0r,il0) &
                       & -avg%m11asysq(jc3,jl0r,il0)*avg_lr%m11lrm11(jc3,jl0r,il0)
            den(jc3) = avg%m11sq(jc3,jl0r,il0)*avg_lr%m11sq(jc3,jl0r,il0)-avg_lr%m11lrm11(jc3,jl0r,il0)**2
            if ((num(jc3)>0.0).and.(den(jc3)>0.0)) then
               loc_deh%raw(jc3,jl0r,il0) = num(jc3)/den(jc3)
               loc_deh_lr%raw(jc3,jl0r,il0) = num_lr(jc3)/den(jc3)
            end if
         end if
      end do
   end do
end do

! Normalize dual-ensemble hybridization
call curve_normalization(hdata,ib,loc_deh)
call curve_normalization(hdata,ib,loc_deh_lr)

! Compute dual-ensemble hybridization fits
if (bpar%fit_block(ib)) then
   ! Compute fit weight
   if (nam%fit_wgt) then
      loc_deh%fit_wgt = abs(avg%cor)
      loc_deh_lr%fit_wgt = abs(avg_lr%cor)
   end if

   ! Compute initial fit
   call compute_fit(hdata,ib,loc_deh)
   call compute_fit(hdata,ib,loc_deh_lr)
end if

! End associate
end associate

end subroutine compute_dualens

end module hdiag_dualens
