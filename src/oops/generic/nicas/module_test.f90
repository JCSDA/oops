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
use module_apply_com, only: alpha_com_ab,fld_com_gl,fld_com_lg
use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_apply_nicas, only: apply_nicas,apply_nicas_from_sqrt
use module_apply_nicas_sqrt, only: apply_nicas_sqrt,apply_nicas_sqrt_ad
use module_namelist, only: nam
use omp_lib
use tools_const, only: deg2rad,rad2deg,sphere_dist
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr
use type_ctree, only: ctreetype,create_ctree,delete_ctree,find_nearest_neighbors
use type_fields, only: fldtype,alphatype
use type_mpl, only: mpl,mpl_bcast
use type_ndata, only: ndatatype,ndataloctype
use type_timer, only: timertype,timer_start,timer_end

implicit none

private
public :: test_adjoints,test_pos_def,test_mpi,test_dirac,test_perf

contains

!----------------------------------------------------------------------
! Subroutine: test_adjoints
!> Purpose: test adjoints accuracy
!----------------------------------------------------------------------
subroutine test_adjoints(ndata)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< Sampling data

! Local variables
integer :: iproc
real(kind_real) :: sum1,sum2
type(alphatype) :: alpha,alpha_save
type(alphatype) :: alpha1,alpha1_save
type(alphatype) :: alpha2,alpha2_save
type(fldtype) :: fld,fld_save
type(fldtype) :: fld1,fld1_save
type(fldtype) :: fld2,fld2_save

! Allocation
allocate(alpha%val(ndata%ns))
allocate(alpha_save%val(ndata%ns))
allocate(fld%val(ndata%nc0,ndata%nl0))
allocate(fld_save%val(ndata%nc0,ndata%nl0))

! Initialization
call random_number(alpha_save%val)
call random_number(fld_save%val)

! Adjoint test
call interp(ndata,alpha_save,fld)
call interp_ad(ndata,fld_save,alpha)

! Print result
sum1 = sum(alpha%val*alpha_save%val)
sum2 = sum(fld%val*fld_save%val)
write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','Interpolation adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Allocation
allocate(alpha1%val(ndata%ns))
allocate(alpha2%val(ndata%ns))
allocate(alpha1_save%val(ndata%ns))
allocate(alpha2_save%val(ndata%ns))

! Initialization
call random_number(alpha1_save%val)
call random_number(alpha2_save%val)
alpha1%val = alpha1_save%val
alpha2%val = alpha2_save%val

! Adjoint test
call convol(ndata,alpha1)
call convol(ndata,alpha2)

! Print result
sum1 = sum(alpha1%val*alpha2_save%val)
sum2 = sum(alpha2%val*alpha1_save%val)
write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','Convolution adjoint test:   ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

if (nam%lsqrt) then
   ! Allocation
   allocate(fld1_save%val(ndata%nc0,ndata%nl0))

   ! Initialization
   call random_number(fld1_save%val)
   call random_number(alpha1_save%val)

   ! Adjoint test
   call apply_nicas_sqrt(ndata,alpha1_save,fld1)
   call apply_nicas_sqrt_ad(ndata,fld1_save,alpha1)

   ! Print result
   sum1 = sum(alpha1%val*alpha1_save%val)
   sum2 = sum(fld1%val*fld1_save%val)
   write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','NICAS square-root adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
else
   ! Allocation
   allocate(fld1%val(ndata%nc0,ndata%nl0))
   allocate(fld1_save%val(ndata%nc0,ndata%nl0))
   allocate(fld2%val(ndata%nc0,ndata%nl0))
   allocate(fld2_save%val(ndata%nc0,ndata%nl0))

   ! Initialization
   call random_number(fld1_save%val)
   call random_number(fld2_save%val)
   fld1%val = fld1_save%val
   fld2%val = fld2_save%val

   ! Adjoint test
   call apply_nicas(ndata,fld1)
   call apply_nicas(ndata,fld2)

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
subroutine test_pos_def(ndata)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< Sampling data

! Local variables
real(kind_real),parameter :: tol = 1.0e-3
integer,parameter :: nitermax = 1
integer :: iter
real(kind_real) :: norm,egvmax,egvmax_prev,egvmin,egvmin_prev
type(fldtype) :: fld,fld_prev

! Allocation
allocate(fld_prev%val(ndata%nc0,ndata%nl0))
allocate(fld%val(ndata%nc0,ndata%nl0))

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
      call apply_nicas_from_sqrt(ndata,fld)
   else
      call apply_nicas(ndata,fld)
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
      call apply_nicas_from_sqrt(ndata,fld)
   else
      call apply_nicas(ndata,fld)
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
!> Purpose: test global/local equivalence
!----------------------------------------------------------------------
subroutine test_mpi(ndata,ndataloc)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata       !< Sampling data
type(ndataloctype),intent(in) :: ndataloc !< Sampling data, local

! Local variables
type(fldtype) :: fldg,fldl

! Allocation
if (mpl%main) then
   allocate(fldg%val(ndata%nc0,ndata%nl0))
   allocate(fldl%val(ndata%nc0,ndata%nl0))

   ! Initialization
   call random_number(fldg%val)
   fldl%val = fldg%val
end if

! Global to local
call fld_com_gl(ndata,ndataloc,fldl)

if (nam%lsqrt) then
   ! Global
   if (mpl%main) call apply_nicas_from_sqrt(ndata,fldg)

   ! Local
   call apply_nicas_from_sqrt(ndataloc,fldl)
else
   ! Global
   if (mpl%main) call apply_nicas(ndata,fldg)

   ! Local
   call apply_nicas(ndataloc,fldl)
end if

! Local to global
call fld_com_lg(ndata,ndataloc,fldl)

! Print difference
if (mpl%main) write(mpl%unit,'(a7,a,e14.8)') '','RMSE between single-proc and multi-procs execution: ', &
 & sqrt(sum((fldg%val-fldl%val)**2)/float(ndata%nc0*ndata%nl0))

end subroutine test_mpi

!----------------------------------------------------------------------
! Subroutine: test_dirac
!> Purpose: apply NICAS to Diracs
!----------------------------------------------------------------------
subroutine test_dirac(ndata,ndataloc)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata          !< Sampling data
type(ndataloctype),intent(inout) :: ndataloc !< Sampling data, local

! Local variables
integer :: idir,ic0_dir(nam%ndir),ic0a_dir(nam%ndir),il0_dir
type(fldtype) :: fld

! Define dirac
call define_dirac(ndata,ndataloc,ic0_dir,ic0a_dir,il0_dir)

! Allocation
allocate(fld%vala(ndataloc%nc0a,ndataloc%nl0))

! Generate diracs field
fld%vala = 0.0
do idir=1,nam%ndir
   if (isnotmsi(ic0a_dir(idir))) fld%vala(ic0a_dir(idir),il0_dir) = 1.0
end do

! Apply NICAS method
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(ndataloc,fld)
else
   call apply_nicas(ndataloc,fld)
end if

! Local to global
call fld_com_lg(ndata,ndataloc,fld)

if (mpl%main) then
   ! Write field
   call system('rm -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_dirac.nc')
   call model_write(trim(nam%prefix)//'_dirac.nc','hr',ndata,fld%val)

   ! Print results
   write(mpl%unit,'(a7,a)') '','Normalization:'
   do idir=1,nam%ndir
      write(mpl%unit,'(a10,f6.1,a,f6.1,a,f10.7)') '',nam%dirlon(idir),' / ',nam%dirlat(idir),': ',fld%val(ic0_dir(idir),il0_dir)
   end do
   write(mpl%unit,'(a7,a,f10.7,a,f10.7)') '','Min - max: ', &
 & minval(fld%val(:,il0_dir),mask=ndata%mask(:,il0_dir)),' - ',maxval(fld%val(:,il0_dir),mask=ndata%mask(:,il0_dir))
end if

end subroutine test_dirac

!----------------------------------------------------------------------
! Subroutine: define_dirac
!> Purpose: define field with Dirac functions
!----------------------------------------------------------------------
subroutine define_dirac(ndata,ndataloc,ic0_dir,ic0a_dir,il0_dir)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata       !< Sampling data
type(ndataloctype),intent(in) :: ndataloc !< Sampling data, local
integer,intent(out) :: ic0_dir(nam%ndir)
integer,intent(out) :: ic0a_dir(nam%ndir)
integer,intent(out) :: il0_dir

! Local variables
integer :: ic0,idir,il0,ic0a,iproc
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
type(ctreetype) :: ctree

if (mpl%main) then
   ! Compute cover tree
   allocate(mask_ctree(ndata%nc0))
   do ic0=1,ndata%nc0
      if (any(ndata%mask(ic0,:))) then
         mask_ctree(ic0) = 1
      else
         mask_ctree(ic0) = 0
      end if
   end do
   ctree = create_ctree(ndata%nc0,dble(ndata%lon),dble(ndata%lat),mask_ctree)

   ! Compute nearest neighbors
   do idir=1,nam%ndir
      call find_nearest_neighbors(ctree,dble(nam%dirlon(idir)*deg2rad), &
    & dble(nam%dirlat(idir)*deg2rad),1,ic0_dir(idir:idir),dum)
   end do
   
   ! Find level index
   do il0=1,ndata%nl0
      if (nam%levs(il0)==nam%dirlev) il0_dir = il0
   end do

   ! Release memory
   deallocate(mask_ctree)
   call delete_ctree(ctree)
end if

! Broadcast
call mpl_bcast(ic0_dir,mpl%ioproc)
call mpl_bcast(il0_dir,mpl%ioproc)

! Transfer to local
call msi(ic0a_dir)
do idir=1,nam%ndir
   ic0 = ic0_dir(idir)
   ic0a = ndata%ic0_to_ic0a(ic0)
   iproc = ndata%ic0_to_iproc(ic0)
   if (iproc==mpl%myproc) ic0a_dir(idir) = ic0a
end do

end subroutine define_dirac

!----------------------------------------------------------------------
! Subroutine: test_perf
!> Purpose: test NICAS performance
!----------------------------------------------------------------------
subroutine test_perf(ndata)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< Sampling data

! Local variables
type(fldtype) :: fld
type(alphatype) :: alpha
type(timertype) :: timer_interp_ad,timer_convol,timer_interp

! Allocation
allocate(fld%val(ndata%nc0,ndata%nl0))
allocate(alpha%val(ndata%ns))

! Random initialization
call random_number(fld%val)

! Adjoint interpolation
call timer_start(timer_interp_ad)
call interp_ad(ndata,fld,alpha)
call timer_end(timer_interp_ad)

! Convolution
call timer_start(timer_convol)
call convol(ndata,alpha)
call timer_start(timer_convol)

! Interpolation
call timer_start(timer_interp)
call interp(ndata,alpha,fld)
call timer_end(timer_interp)

! Print results
write(mpl%unit,'(a7,a)') '','Performance results (elapsed time):'
write(mpl%unit,'(a10,a,f6.1,a)') '','Adjoint interpolation: ',timer_interp_ad%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Convolution          : ',timer_convol%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Interpolation        : ',timer_interp%elapsed,' s'

end subroutine test_perf

end module module_test
