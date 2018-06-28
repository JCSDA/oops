!----------------------------------------------------------------------
! Module: tools_fit
!> Purpose: fit routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module tools_fit

use tools_func, only:  gc99
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr,isanynotmsr,ismsi,ismsr
use type_mpl, only: mpl_type

implicit none

integer,parameter :: itermax = 10 !< Maximum number of iteration for the threshold definition

private
public :: fast_fit,ver_smooth,ver_fill

contains

!----------------------------------------------------------------------
! Subroutine: fast_fit
!> Purpose: compute a fast fit of a raw function
!----------------------------------------------------------------------
subroutine fast_fit(mpl,n,iz,dist,raw,fit_r)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl      !< MPI data
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

if (any(dist<0.0)) call mpl%abort('negative distance in fast_fit')

if (raw(iz)>0.0) then
   ! Initialization
   valid = .false.

   ! At least one positive point is needed
   do i=1,n
     if (i/=iz) then
        if (raw(i)>0.0) valid = .true.
     end if
   end do

   ! The zero-separation point must be the maximum
   do i=1,n
     if (i/=iz) then
        if (.not.(raw(i)<raw(iz))) valid = .false.
     end if
   end do

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
         do iter=1,itermax
            thtest = gc99(mpl,thinv)
            if (th>thtest) then
               thinv = thinv-dthinv
            else
               thinv = thinv+dthinv
            end if
            dthinv = 0.5*dthinv
         end do
      end if

      ! Find support radius, lower value
      call msr(fit_rm)
      ip = iz
      do di=1,n
         ! Check whether fit value has been found
         if (ismsr(fit_rm)) then
            ! Index
            im = iz-di

            ! Check index validity
            if (im>=1) then
               ! Check raw value validity
               if (raw(im)>0.0) then
                  ! Check whether threshold has been crossed
                  if (raw(im)<th*raw(iz)) then
                     ! Set fit value
                     fit_rm = dist(im)+(dist(ip)-dist(im))*(th*raw(iz)-raw(im))/(raw(ip)-raw(im))
                  else
                     ! Update index
                     ip = im
                  end if
               end if
            end if
         end if
      end do

      ! Find support radius, upper value
      call msr(fit_rp)
      im = iz
      do di=1,n
         ! Check whether fit value has been found
         if (ismsr(fit_rp)) then
            ! Index
            ip = iz+di

            ! Check index validity
            if (ip<=n) then
               ! Check raw value validity
               if (raw(ip)>0.0) then
                  ! Check whether threshold has been crossed
                  if (raw(ip)<th*raw(iz)) then
                     ! Set fit value
                     fit_rp = dist(im)+(dist(ip)-dist(im))*(th*raw(iz)-raw(im))/(raw(ip)-raw(im))
                  else
                     ! Update index
                     im = ip
                  end if
               end if
            end if
         end if
      end do

      ! Gather values
      if (isnotmsr(fit_rm).and.isnotmsr(fit_rp)) then
         fit_r = 0.5*(fit_rm+fit_rp)
      elseif (isnotmsr(fit_rm)) then
         fit_r = fit_rm
      elseif (isnotmsr(fit_rp)) then
         fit_r = fit_rp
      end if

      ! Check
      if (fit_r<0.0) then
         write(mpl%unit,*) iz,valid
         write(mpl%unit,*) raw
         write(mpl%unit,*) th,thinv
         write(mpl%unit,*) fit_rm,fit_rp,fit_r
         call mpl%abort('negative fit_r in fast_fit')
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
subroutine ver_smooth(mpl,n,x,rv,profile)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl            !< MPI data
integer,intent(in) :: n                     !< Vector size
real(kind_real),intent(in) :: x(n)          !< Coordinate
real(kind_real),intent(in) :: rv            !< Filtering support radius
real(kind_real),intent(inout) :: profile(n) !< Vertical profile

! Local variables
integer :: i,j
real(kind_real) :: kernel(n,n),distnorm,profile_init(n),norm

if (rv<0.0) call mpl%abort('negative filtering support radius in ver_smooth')

if ((rv>0.0).and.isanynotmsr(profile)) then
   ! Vertical smoothing kernel
   kernel = 0.0
   do i=1,n
      do j=1,n
         if (isnotmsr(profile(j))) then
            ! Gaspari-Cohn (1999) function
            distnorm = abs(x(j)-x(i))/rv
            kernel(i,j) = gc99(mpl,distnorm)
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

!----------------------------------------------------------------------
! Subroutine: ver_fill
!> Purpose: fill missing values
!----------------------------------------------------------------------
subroutine ver_fill(mpl,n,x,profile)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl            !< MPI data
integer,intent(in) :: n                     !< Vector size
real(kind_real),intent(in) :: x(n)          !< Coordinate
real(kind_real),intent(inout) :: profile(n) !< Vertical profile

! Local variables
integer :: i,j,iinf,isup
real(kind_real) :: profile_init(n)

if (isanynotmsr(profile)) then
   ! Initialization
   profile_init = profile
   call msi(iinf)

   do i=1,n
      if (isnotmsr(profile_init(i))) then
         ! Valid inferior point
         iinf = i
      else
         ! Look for a superior point
         call msi(isup)
         j = i+1
         do while ((j<=n).and.(ismsi(isup)))
            if (isnotmsr(profile_init(j))) isup = j
            j = j+1
         end do

         if (isnotmsi(iinf).and.isnotmsi(isup)) then
            ! Interpolation
            profile(i) = profile_init(iinf)+(x(i)-x(iinf))*(profile_init(isup)-profile_init(iinf))/(x(isup)-x(iinf))
         elseif (isnotmsi(isup)) then
            ! Extrapolation with nearest superior point
            profile(i) = profile(isup)
         elseif (isnotmsi(iinf)) then
            ! Extrapolation with nearest inferior point
            profile(i) = profile(iinf)
         else
            call mpl%abort('ver_fill failed')
         end if
      end if
   end do
end if

end subroutine ver_fill

end module tools_fit
