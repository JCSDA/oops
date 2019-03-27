!----------------------------------------------------------------------
! Module: tools_fit
! Purpose: fit-related tools
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_fit

use tools_func, only: gc99
use tools_kinds, only: kind_real
use tools_repro, only: inf,sup
use type_mpl, only: mpl_type

implicit none

integer,parameter :: itermax = 10 ! Maximum number of iteration for the threshold definition

private
public :: fast_fit,ver_smooth,ver_fill

contains

!----------------------------------------------------------------------
! Subroutine: fast_fit
! Purpose: fast fit length-scale estimation based on the value at mid-height
!----------------------------------------------------------------------
subroutine fast_fit(mpl,n,iz,dist,raw,fit_r)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: n               ! Vector size
integer,intent(in) :: iz              ! Zero separation index
real(kind_real),intent(in) :: dist(n) ! Distance
real(kind_real),intent(in) :: raw(n)  ! Raw data
real(kind_real),intent(out) :: fit_r  ! Fast fit result

! Local variables
integer :: di,i,im,ip,iter
real(kind_real) :: th,thinv,dthinv,thtest
real(kind_real) :: fit_rm,fit_rp,distmin
real(kind_real) :: raw_tmp(n)
logical :: valid
character(len=1024),parameter :: subr = 'fast_fit'

if (any(dist<0.0)) call mpl%abort(subr,'negative distance in fast_fit')

if (raw(iz)>0.0) then
   if (n>1) then
      ! Copy points that are lower than the zero-separation
      raw_tmp = mpl%msv%valr
      raw_tmp(iz) = raw(iz)
      do i=1,n
        if (i/=iz) then
           if (inf(raw(i),raw(iz))) raw_tmp(i) = raw(i)
        end if
      end do

      ! At least one positive point is needed
      valid = .false.
      do i=1,n
        if (i/=iz) then
           if (raw_tmp(i)>0.0) valid = .true.
        end if
      end do

      if (valid) then
         ! Define threshold and its inverse
         if (inf(minval(raw_tmp/raw_tmp(iz),mask=(raw_tmp>0.0)),0.5_kind_real)) then
            ! Default threshold
            th = 0.5
            thinv = 0.33827292796125663
         else
            ! Curve-dependent threshold
            th = minval(raw_tmp/raw_tmp(iz),mask=(raw_tmp>0.0))+1.0e-6

            ! Find inverse threshold by dichotomy
            thinv = 0.5
            dthinv = 0.25
            do iter=1,itermax
               thtest = gc99(mpl,thinv)
               if (sup(th,thtest)) then
                  thinv = thinv-dthinv
               else
                  thinv = thinv+dthinv
               end if
               dthinv = 0.5*dthinv
            end do
         end if

         ! Find support radius, lower value
         fit_rm = mpl%msv%valr
         ip = iz
         do di=1,n
            ! Check whether fit value has been found
            if (mpl%msv%isr(fit_rm)) then
               ! Index
               im = iz-di

               ! Check index validity
               if (im>=1) then
                  ! Check raw value validity
                  if (raw_tmp(im)>0.0) then
                     ! Check whether threshold has been crossed
                     if (inf(raw_tmp(im),th*raw_tmp(iz))) then
                        ! Set fit value
                        fit_rm = dist(im)+(dist(ip)-dist(im))*(th*raw_tmp(iz)-raw_tmp(im))/(raw_tmp(ip)-raw_tmp(im))
                     else
                        ! Update index
                        ip = im
                     end if
                  end if
               end if
            end if
         end do

         ! Find support radius, upper value
         fit_rp = mpl%msv%valr
         im = iz
         do di=1,n
            ! Check whether fit value has been found
            if (mpl%msv%isr(fit_rp)) then
               ! Index
               ip = iz+di

               ! Check index validity
               if (ip<=n) then
                  ! Check raw value validity
                  if (raw_tmp(ip)>0.0) then
                     ! Check whether threshold has been crossed
                     if (inf(raw_tmp(ip),th*raw_tmp(iz))) then
                        ! Set fit value
                        fit_rp = dist(im)+(dist(ip)-dist(im))*(th*raw_tmp(iz)-raw_tmp(im))/(raw_tmp(ip)-raw_tmp(im))
                     else
                        ! Update index
                        im = ip
                     end if
                  end if
               end if
            end if
         end do

         ! Gather values
         if (mpl%msv%isnotr(fit_rm).and.mpl%msv%isnotr(fit_rp)) then
            fit_r = 0.5*(fit_rm+fit_rp)
         elseif (mpl%msv%isnotr(fit_rm)) then
            fit_r = fit_rm
         elseif (mpl%msv%isnotr(fit_rp)) then
            fit_r = fit_rp
         end if

         ! Normalize
         if (mpl%msv%isnotr(fit_r)) fit_r = fit_r/thinv

         ! Check positivity
         if (inf(fit_r,0.0_kind_real)) then
            call mpl%warning(subr,'negative fit_r in fast_fit')
            fit_r = 0.0
         end if
      else
         ! Only the zero-separation point is valid, zero radius
         fit_r = 0.0
      end if

      ! Set minimum distance
      distmin = huge(1.0)
      if (iz>1) distmin = min(distmin,abs(dist(iz-1)-dist(iz)))
      if (iz<n) distmin = min(distmin,abs(dist(iz+1)-dist(iz)))
      fit_r = max(fit_r,distmin)
   else
      ! Only one point, zero radius
      fit_r = 0.0
   end if
else
   ! Zero-separation point is negative
   fit_r = mpl%msv%valr
end if

end subroutine fast_fit

!----------------------------------------------------------------------
! Subroutine: ver_smooth
! Purpose: homogeneous smoothing of a vertical profile
!----------------------------------------------------------------------
subroutine ver_smooth(mpl,n,x,rv,profile)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl         ! MPI data
integer,intent(in) :: n                     ! Vector size
real(kind_real),intent(in) :: x(n)          ! Coordinate
real(kind_real),intent(in) :: rv            ! Filtering support radius
real(kind_real),intent(inout) :: profile(n) ! Vertical profile

! Local variables
integer :: i,j
real(kind_real) :: kernel(n,n),distnorm,profile_init(n),norm
character(len=1024),parameter :: subr = 'ver_smooth'

if (rv<0.0) call mpl%abort(subr,'negative filtering support radius in ver_smooth')

if ((rv>0.0).and.mpl%msv%isanynotr(profile)) then
   ! Vertical smoothing kernel
   kernel = 0.0
   do i=1,n
      do j=1,n
         if (mpl%msv%isnotr(profile(j))) then
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
         profile(i) = mpl%msv%valr
      end if
   end do
end if

end subroutine ver_smooth

!----------------------------------------------------------------------
! Subroutine: ver_fill
! Purpose: missing values filling of a vertical profile
!----------------------------------------------------------------------
subroutine ver_fill(mpl,n,x,profile)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl         ! MPI data
integer,intent(in) :: n                     ! Vector size
real(kind_real),intent(in) :: x(n)          ! Coordinate
real(kind_real),intent(inout) :: profile(n) ! Vertical profile

! Local variables
integer :: i,j,iinf,isup
real(kind_real) :: profile_init(n)
character(len=1024),parameter :: subr = 'ver_fill'

if (mpl%msv%isanynotr(profile)) then
   ! Initialization
   profile_init = profile
   iinf = mpl%msv%vali

   do i=1,n
      if (mpl%msv%isnotr(profile_init(i))) then
         ! Valid inferior point
         iinf = i
      else
         ! Look for a superior point
         isup = mpl%msv%vali
         j = i+1
         do while ((j<=n).and.(mpl%msv%isi(isup)))
            if (mpl%msv%isnotr(profile_init(j))) isup = j
            j = j+1
         end do

         if (mpl%msv%isnoti(iinf).and.mpl%msv%isnoti(isup)) then
            ! Interpolation
            profile(i) = profile_init(iinf)+(x(i)-x(iinf))*(profile_init(isup)-profile_init(iinf))/(x(isup)-x(iinf))
         elseif (mpl%msv%isnoti(isup)) then
            ! Extrapolation with nearest superior point
            profile(i) = profile(isup)
         elseif (mpl%msv%isnoti(iinf)) then
            ! Extrapolation with nearest inferior point
            profile(i) = profile(iinf)
         else
            call mpl%abort(subr,'ver_fill failed')
         end if
      end if
   end do
end if

end subroutine ver_fill

end module tools_fit
