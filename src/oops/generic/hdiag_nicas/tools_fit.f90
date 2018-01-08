!----------------------------------------------------------------------
! Module: tools_fit
!> Purpose: fit routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_fit

use tools_const, only:  gc99
use tools_display, only: msgerror,msgwarning
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isanynotmsr
use type_mpl, only: mpl
implicit none

integer,parameter :: niter = 50 !< Maximum number of iterations for the threshold definition

private
public :: fast_fit,ver_smooth

contains

!----------------------------------------------------------------------
! Subroutine: fast_fit
!> Purpose: compute a fast fit of a raw function
!----------------------------------------------------------------------
subroutine fast_fit(n,iz,dist,raw,fit_r)

implicit none

! Passed variables
integer,intent(in) :: n               !< Vector size
integer,intent(in) :: iz              !< Zero separation index
real(kind_real),intent(in) :: dist(n) !< Distance
real(kind_real),intent(in) :: raw(n)  !< Raw data
real(kind_real),intent(out) :: fit_r  !< Fast fit result

! Local variables
integer :: di,i,im,ip,iter
real(kind_real) :: th,thinv,dthinv,thtest
real(kind_real) :: fit_rm,fit_rp
logical :: valid

if (raw(iz)>0.0) then
   valid = .false.
   do i=1,n
     if (i/=iz) then
        if (raw(i)>0.0) valid = .true.
     end if
   end do

   fit_r = 0.0
   if (valid) then
      ! Define threshold and its inverse
      if (minval(raw/raw(iz),mask=(raw>0.0))<0.5) then
         ! Default threshold
         th = 0.5
         thinv = 0.33827292796125663
      else
         ! Curve-dependent threshold
         th = minval(raw/raw(iz),mask=(raw>0.0))+1.0e-6

         ! Find inverse threshold by dichotomy
         thinv = 0.5
         dthinv = 0.25
         do iter=1,niter
            thtest = gc99(thinv)
            if (th>thtest) then
               thinv = thinv-dthinv
            else
               thinv = thinv+dthinv
            end if
            dthinv = 0.5*dthinv
         end do
      end if

      ! Find support radius
      call msr(fit_rm)
      call msr(fit_rp)
      do di=1,n
         if (.not.isnotmsr(fit_rm)) then
            i = iz-di
            if (i>=1) then
               if (isnotmsr(raw(i)).and.(raw(i)<th*raw(iz))) then
                  im = i
                  ip = i+1
                  fit_rm = dist(im)+(dist(ip)-dist(im))*(th*raw(iz)-raw(im))/(raw(ip)-raw(im))
               end if
            end if
         end if
         if (.not.isnotmsr(fit_rp)) then
            i = iz+di
            if (i<=n) then
               if (isnotmsr(raw(i)).and.(raw(i)<th*raw(iz))) then
                  im = i-1
                  ip = i
                  fit_rp = dist(im)+(dist(ip)-dist(im))*(th*raw(iz)-raw(im))/(raw(ip)-raw(im))
               end if
            end if
         end if
      end do
      if (isnotmsr(fit_rm).and.isnotmsr(fit_rp)) then
         fit_r = 0.5*(fit_rm+fit_rp)
      elseif (isnotmsr(fit_rm)) then
         fit_r = fit_rm
      elseif (isnotmsr(fit_rp)) then
         fit_r = fit_rp
      end if

      ! Normalize
      if (isnotmsr(fit_r)) fit_r = fit_r/thinv
   end if
else
   call msr(fit_r)
end if

end subroutine fast_fit

!----------------------------------------------------------------------
! Subroutine: ver_smooth
!> Purpose: homogeneous vertical smoothing
!----------------------------------------------------------------------
subroutine ver_smooth(n,x,rv,profile)

implicit none

! Passed variables
integer,intent(in) :: n                     !< Vector size
real(kind_real),intent(in) :: x(n)          !< Coordinate
real(kind_real),intent(in) :: rv            !< Filtering support radius
real(kind_real),intent(inout) :: profile(n) !< Vertical profile

! Local variables
integer :: i,j
real(kind_real) :: kernel(n,n),distnorm,profile_init(n),norm

if ((rv>0.0).and.isanynotmsr(profile)) then
   ! Vertical smoothing kernel
   kernel = 0.0
   do i=1,n
      do j=1,n
         if (isnotmsr(profile(j))) then
            ! Gaspari-Cohn (1999) function
            distnorm = abs(x(j)-x(i))/rv
            kernel(i,j) = gc99(distnorm)
         end if
      end do
   end do

   ! Apply kernel
   profile_init = profile
   profile = 0.0
   do i=1,n
      norm = 0.0
      do j=1,n
         profile(i) = profile(i)+kernel(i,j)*profile_init(j)
         norm = norm+kernel(i,j)
      end do
      if (norm>0.0) then
         profile(i) = profile(i)/norm
      else
         call msr(profile(i))
      end if
   end do
end if

end subroutine ver_smooth

end module tools_fit
