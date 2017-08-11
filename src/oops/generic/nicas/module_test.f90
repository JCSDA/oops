!----------------------------------------------------------------------
! Module: module_test.f90
!> Purpose: test NICAS
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_test

use model_interface, only: model_write
use module_apply_com, only: alpha_com_ab
use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_apply_nicas, only: apply_nicas_sqrt,apply_nicas_sqrt_ad,apply_nicas,apply_nicas_from_sqrt
use module_namelist, only: nam
use omp_lib
use tools_const, only: deg2rad,rad2deg,sphere_dist
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: isnotmsr,msr
use type_fields, only: fldtype,alphatype
use type_mpl, only: mpl
use type_sdata, only: sdatatype
use type_timer, only: timertype,timer_start,timer_end

implicit none

private
public :: test_adjoints,test_pos_def,test_mpi,test_dirac,test_perf

contains

!----------------------------------------------------------------------
! Subroutine: test_adjoints
!> Purpose: test adjoints accuracy
!----------------------------------------------------------------------
subroutine test_adjoints(sdata)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata !< Sampling data

! Local variables
integer :: iproc
real(kind_real) :: sum1,sum2
type(alphatype) :: alpha(sdata%nproc),alpha_save(sdata%nproc)
type(alphatype) :: alpha1(sdata%nproc),alpha1_save(sdata%nproc)
type(alphatype) :: alpha2(sdata%nproc),alpha2_save(sdata%nproc)
type(fldtype) :: fld(sdata%nproc),fld_save(sdata%nproc)
type(fldtype) :: fld1,fld1_save
type(fldtype) :: fld2,fld2_save

if (nam%nproc==0) then
   ! Allocation
   allocate(alpha(1)%val(sdata%ns))
   allocate(alpha_save(1)%val(sdata%ns))
   allocate(fld(1)%val(sdata%nc0,sdata%nl0))
   allocate(fld_save(1)%val(sdata%nc0,sdata%nl0))

   ! Initialization
   call random_number(alpha_save(1)%val)
   call random_number(fld_save(1)%val)
elseif (nam%nproc>0) then
   do iproc=1,sdata%nproc
      ! Allocation
      allocate(alpha(iproc)%valb(sdata%mpi(iproc)%nsb))
      allocate(alpha_save(iproc)%vala(sdata%mpi(iproc)%nsa))
      allocate(fld(iproc)%vala(sdata%mpi(iproc)%nc0a,sdata%nl0))
      allocate(fld_save(iproc)%vala(sdata%mpi(iproc)%nc0a,sdata%nl0))

      ! Initialization
      call random_number(alpha_save(iproc)%vala)
      call random_number(fld_save(iproc)%vala)
   end do

   ! Halo extension from zone A to zone B
   call alpha_com_AB(sdata,alpha_save)
end if

! Adjoint test
call interp(sdata,alpha_save,fld)
call interp_ad(sdata,fld_save,alpha)

! Print result
if (nam%nproc==0) then
   sum1 = sum(alpha(1)%val*alpha_save(1)%val)
   sum2 = sum(fld(1)%val*fld_save(1)%val)
elseif (nam%nproc>0) then
   sum1 = 0.0
   sum2 = 0.0
   do iproc=1,sdata%nproc
      sum1 = sum1+sum(alpha(iproc)%valb*alpha_save(iproc)%valb)
      sum2 = sum2+sum(fld(iproc)%vala*fld_save(iproc)%vala)
   end do
end if
write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','Interpolation adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

if (nam%nproc==0) then
   ! Allocation
   allocate(alpha1(1)%val(sdata%ns))
   allocate(alpha2(1)%val(sdata%ns))
   allocate(alpha1_save(1)%val(sdata%ns))
   allocate(alpha2_save(1)%val(sdata%ns))

   ! Initialization
   call random_number(alpha1_save(1)%val)
   call random_number(alpha2_save(1)%val)
   alpha1(1)%val = alpha1_save(1)%val
   alpha2(1)%val = alpha2_save(1)%val
elseif (nam%nproc>0) then
   do iproc=1,sdata%nproc
      ! Allocation
      allocate(alpha1(iproc)%vala(sdata%mpi(iproc)%nsa))
      allocate(alpha2(iproc)%vala(sdata%mpi(iproc)%nsa))
      allocate(alpha1_save(iproc)%vala(sdata%mpi(iproc)%nsa))
      allocate(alpha2_save(iproc)%vala(sdata%mpi(iproc)%nsa))

      ! Initialization
      call random_number(alpha1_save(iproc)%vala)
      call random_number(alpha2_save(iproc)%vala)
      alpha1(iproc)%vala = alpha1_save(iproc)%vala
      alpha2(iproc)%vala = alpha2_save(iproc)%vala
   end do

   ! Halo extension from zone A to zone B
   call alpha_com_AB(sdata,alpha1)
   call alpha_com_AB(sdata,alpha2)
   call alpha_com_AB(sdata,alpha1_save)
   call alpha_com_AB(sdata,alpha2_save)
end if

! Adjoint test
call convol(sdata,alpha1)
call convol(sdata,alpha2)

! Print result
if (nam%nproc==0) then
   sum1 = sum(alpha1(1)%val*alpha2_save(1)%val)
   sum2 = sum(alpha2(1)%val*alpha1_save(1)%val)
elseif (nam%nproc>0) then
   sum1 = 0.0
   sum2 = 0.0
   do iproc=1,sdata%nproc
      sum1 = sum1+sum(alpha1(iproc)%valb*alpha2_save(iproc)%valb)
      sum2 = sum2+sum(alpha2(iproc)%valb*alpha1_save(iproc)%valb)
   end do
end if
write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','Convolution adjoint test:   ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

if (nam%lsqrt) then
   ! Allocation
   allocate(fld1_save%val(sdata%nc0,sdata%nl0))

   ! Initialization
   call random_number(fld1_save%val)
   if (nam%nproc==0) then
      call random_number(alpha1_save(1)%val)
   elseif (nam%nproc>0) then
      do iproc=1,sdata%nproc
         call random_number(alpha1_save(iproc)%vala)
      end do

      ! Halo extension from zone A to zone B
      call alpha_com_AB(sdata,alpha1_save)
   end if

   ! Adjoint test
   call apply_nicas_sqrt(sdata,alpha1_save,fld1)
   call apply_nicas_sqrt_ad(sdata,fld1_save,alpha1)

   ! Print result
   if (nam%nproc==0) then
      sum1 = sum(alpha1(1)%val*alpha1_save(1)%val)
   elseif (nam%nproc>0) then
      sum1 = 0.0
      do iproc=1,sdata%nproc
         sum1 = sum1+sum(alpha1(iproc)%vala*alpha1_save(iproc)%vala)
      end do
   end if
   sum2 = sum(fld1%val*fld1_save%val)
   write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','NICAS square-root adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
else
   ! Allocation
   allocate(fld1%val(sdata%nc0,sdata%nl0))
   allocate(fld1_save%val(sdata%nc0,sdata%nl0))
   allocate(fld2%val(sdata%nc0,sdata%nl0))
   allocate(fld2_save%val(sdata%nc0,sdata%nl0))

   ! Initialization
   call random_number(fld1_save%val)
   call random_number(fld2_save%val)
   fld1%val = fld1_save%val
   fld2%val = fld2_save%val

   ! Adjoint test
   call apply_nicas(sdata,fld1)
   call apply_nicas(sdata,fld2)

   ! Print result
   sum1 = sum(fld1%val*fld2_save%val)
   sum2 = sum(fld2%val*fld1_save%val)
   write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','NICAS adjoint test:         ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
end if

end subroutine test_adjoints

!----------------------------------------------------------------------
! Subroutine: test_pos_def
!> Purpose: test positive_definiteness
!----------------------------------------------------------------------
subroutine test_pos_def(sdata)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata !< Sampling data

! Local variables
real(kind_real),parameter :: tol = 1.0e-3
integer,parameter :: nitermax = 1
integer :: iter
real(kind_real) :: norm,egvmax,egvmax_prev,egvmin,egvmin_prev
type(fldtype) :: fld,fld_prev

! Allocation
allocate(fld_prev%val(sdata%nc0,sdata%nl0))
allocate(fld%val(sdata%nc0,sdata%nl0))

! Power method to find the largest eigenvalue
call random_number(fld_prev%val)
norm = sum(fld_prev%val**2)
fld_prev%val = fld_prev%val/norm
egvmax_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld%val = fld_prev%val

   ! Apply C
   if (nam%lsqrt) then
      call apply_nicas_from_sqrt(sdata,fld)
   else
      call apply_nicas(sdata,fld)
   end if

   ! Compute Rayleigh quotient
   egvmax = sum(fld%val*fld_prev%val)/sum(fld_prev%val*fld_prev%val)

   ! Renormalize the vector
   norm = sum(fld%val**2)
   fld%val = fld%val/norm

   ! Exit test
   if (abs(egvmax-egvmax_prev)<tol) exit

   ! Update
   iter = iter+1
   fld_prev%val = fld%val
   egvmax_prev = egvmax
end do

! Power method to find the smallest eigenvalue
call random_number(fld_prev%val)
norm = sum(fld_prev%val**2)
egvmin_prev = huge(1.0)
fld_prev%val = fld_prev%val/norm
egvmin_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld%val = fld_prev%val

   ! Apply C
   if (nam%lsqrt) then
      call apply_nicas_from_sqrt(sdata,fld)
   else
      call apply_nicas(sdata,fld)
   end if
   fld%val = fld%val-egvmax*fld_prev%val

   ! Compute Rayleigh quotient
   egvmin = sum(fld%val*fld_prev%val)/sum(fld_prev%val*fld_prev%val)

   ! Renormalize the vector
   norm = sum(fld%val**2)
   fld%val = fld%val/norm

   ! Exit test
   if (egvmax+egvmin<-tol*egvmax) then
      write(mpl%unit,'(a7,a)') '','NICAS is not positive definite'
      exit
   end if

   ! Update
   iter = iter+1
   fld_prev%val = fld%val
   egvmin_prev = egvmin
end do

! Non conclusive test
if (iter==nitermax+1) write(mpl%unit,'(a7,a,e14.8,a,i4,a,e14.8)') '','NICAS seems to be positive definite: difference ', &
 & egvmax+egvmin,' after ',nitermax,' iterations for a tolerance ',tol

end subroutine test_pos_def

!----------------------------------------------------------------------
! Subroutine: test_mpi
!> Purpose: test single/mutli-procs equivalence
!----------------------------------------------------------------------
subroutine test_mpi(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: nproc_save,sdata_nproc_save
type(fldtype) :: flds,fldm

! Allocation
allocate(flds%val(sdata%nc0,sdata%nl0))
allocate(fldm%val(sdata%nc0,sdata%nl0))

! Initialization
call random_number(flds%val)
fldm%val = flds%val

! Multi-procs execution
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(sdata,fldm)
else
   call apply_nicas(sdata,fldm)
end if
! Single proc execution
nproc_save = nam%nproc
sdata_nproc_save = sdata%nproc
nam%nproc = 0
sdata%nproc = 1
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(sdata,flds)
else
   call apply_nicas(sdata,flds)
end if
nam%nproc = nproc_save
sdata%nproc = sdata_nproc_save

! Print difference
write(mpl%unit,'(a7,a,e14.8)') '','RMSE between single-proc and multi-procs execution: ', &
 & sqrt(sum((flds%val-fldm%val)**2)/float(sdata%nc0*sdata%nl0))

end subroutine test_mpi

!----------------------------------------------------------------------
! Subroutine: test_dirac
!> Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine test_dirac(sdata)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata !< Sampling data

! Local variables
integer :: ic0,idir,il0
integer,allocatable :: ic0_dir(:)
real(kind_real) :: dist
real(kind_real),allocatable :: distmin(:)
character(len=3) :: levchar
type(fldtype) :: dirac

! Find nearest point (brute force...)
allocate(distmin(nam%ndir))
allocate(ic0_dir(nam%ndir))
distmin = huge(1.0)
do ic0=1,sdata%nc0
   if (sdata%mask(ic0,1)) then
      do idir=1,nam%ndir
         call sphere_dist(nam%dirlon(idir)*deg2rad,nam%dirlat(idir)*deg2rad,sdata%lon(ic0),sdata%lat(ic0),dist)
         if (dist<distmin(idir)) then
            distmin(idir) = dist
            ic0_dir(idir) = ic0
         end if
      end do
   end if
end do

! Allocation
allocate(dirac%val(sdata%nc0,sdata%nl0))

do il0=1,sdata%nl0
   ! Remove files
   write(levchar,'(i3.3)') nam%levs(il0)
   call system('rm -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_dirac_'//levchar//'.nc')

   ! Generate diracs field
   dirac%val = 0.0
   do idir=1,nam%ndir
      dirac%val(ic0_dir(idir),il0) = 1.0
   end do
   do ic0=1,sdata%nc0
      if (.not.sdata%mask(ic0,il0)) call msr(dirac%val(ic0,il0))
   end do

   ! Apply NICAS method
   if (nam%lsqrt) then
      call apply_nicas_from_sqrt(sdata,dirac)
   else
      call apply_nicas(sdata,dirac)
   end if

   ! Write field
   call model_write(trim(nam%prefix)//'_dirac_'//levchar//'.nc','hr',sdata,dirac%val)

   ! Print results
   write(mpl%unit,'(a7,a)') '','Normalization at level '//levchar
   do idir=1,nam%ndir
      write(mpl%unit,'(a10,f6.1,a,f6.1,a,f10.7)') '',nam%dirlon(idir),' / ',nam%dirlat(idir),': ',dirac%val(ic0_dir(idir),il0)
   end do
   write(mpl%unit,'(a7,a,f10.7,a,f10.7)') '','Min - max at level '//levchar//': ', &
 & minval(dirac%val(:,il0),mask=sdata%mask(:,il0)),' - ',maxval(dirac%val(:,il0),mask=sdata%mask(:,il0))
end do

end subroutine test_dirac

!----------------------------------------------------------------------
! Subroutine: test_perf
!> Purpose: test NICAS performance
!----------------------------------------------------------------------
subroutine test_perf(sdata)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata !< Sampling data

! Local variables
type(fldtype) :: fld(sdata%nproc)
type(alphatype) :: alpha(sdata%nproc)
type(timertype) :: timer_interp_ad,timer_convol,timer_interp

if (nam%nproc==0) then
   ! Allocation
   allocate(fld(1)%val(sdata%nc0,sdata%nl0))
   allocate(alpha(1)%val(sdata%ns))

   ! Random initialization
   call random_number(fld(1)%val)

! Adjoint interpolation
call timer_start(timer_interp_ad)
call interp_ad(sdata,fld,alpha)
call timer_end(timer_interp_ad)

! Convolution
call timer_start(timer_convol)
call convol(sdata,alpha)
call timer_start(timer_convol)

! Interpolation
call timer_start(timer_interp)
call interp(sdata,alpha,fld)
call timer_end(timer_interp)

elseif (nam%nproc>0) then
   
end if

! Print results
write(mpl%unit,'(a7,a)') '','Performance results (elapsed time):'
write(mpl%unit,'(a10,a,f6.1,a)') '','Adjoint interpolation: ',timer_interp_ad%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Convolution          : ',timer_convol%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Interpolation        : ',timer_interp%elapsed,' s'

end subroutine test_perf
end module module_test
