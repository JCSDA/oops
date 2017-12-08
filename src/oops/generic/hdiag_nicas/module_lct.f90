!----------------------------------------------------------------------
! Module: module_lct.f90
!> Purpose: LCT routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_lct

use model_interface, only: model_write
use module_fit_lct, only: compute_fit_lct
use omp_lib
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr
use type_lct, only: lcttype,lct_alloc,lct_dealloc
use type_mom, only: momtype
use type_mpl, only: mpl
use type_hdata, only: hdatatype
implicit none

private
public :: compute_lct

contains

!----------------------------------------------------------------------
! Subroutine: compute_lct
!> Purpose: compute LCT
!----------------------------------------------------------------------
subroutine compute_lct(hdata,mom,lct)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                                          !< HDIAG data
type(momtype),intent(in) :: mom(hdata%bpar%nb)                               !< Moments
type(lcttype),intent(out) :: lct(hdata%nam%nc1,hdata%geom%nl0,hdata%bpar%nb) !< LCT


! Local variables
integer :: ib,isub,il0,jl0,ic,ic1
real(kind_real) :: den

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

   ! Allocation
   do il0=1,geom%nl0
      do ic1=1,nam%nc1
         call lct_alloc(hdata,ib,lct(ic1,il0,ib))
      end do
   end do

   ! Compute correlations
   write(mpl%unit,'(a10,a)') '','Compute correlations'
   do isub=1,mom(ib)%nsub
      do jl0=1,geom%nl0
         do il0=1,bpar%nl0(ib)
            do ic=1,nam%nc
               do ic1=1,nam%nc1
                  den = mom(ib)%m2_1(ic1,ic,il0,jl0,isub)*mom(ib)%m2_2(ic1,ic,il0,jl0,isub)
                  if (den>0.0) then
                     lct(ic1,jl0,ib)%raw(ic,il0) = lct(ic1,jl0,ib)%raw(ic,il0)+mom(ib)%m11(ic1,ic,il0,jl0,isub)/sqrt(den)
                     lct(ic1,jl0,ib)%norm(ic,il0) = lct(ic1,jl0,ib)%norm(ic,il0)+1.0
                  end if
               end do
            end do
         end do
      end do
   end do
   do jl0=1,geom%nl0
      do il0=1,bpar%nl0(ib)
         do ic=1,nam%nc
            do ic1=1,nam%nc
               if (lct(ic1,jl0,ib)%norm(ic,il0)>0.0) lct(ic1,jl0,ib)%raw(ic,il0) = lct(ic1,jl0,ib)%raw(ic,il0) &
                                                                                 & /lct(ic1,jl0,ib)%norm(ic,il0)
            end do
         end do
      end do
   end do

   ! Compute LCT fit
   write(mpl%unit,'(a10,a)') '','Compute LCT fit'
   call compute_fit_lct(hdata,ib,lct(:,:,ib))
end do

! End associate
end associate

end subroutine compute_lct

end module module_lct
