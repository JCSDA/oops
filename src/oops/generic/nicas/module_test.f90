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
use module_apply_convol, only: convol
use module_apply_interp, only: interp,interp_ad
use module_apply_nicas, only: apply_nicas,apply_nicas_from_sqrt
use module_apply_nicas_sqrt, only: apply_nicas_sqrt,apply_nicas_sqrt_ad
use module_namelist, only: namtype
use omp_lib
use tools_const, only: deg2rad,rad2deg,sphere_dist
use tools_dirac, only: setup_dirac_points
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr
use type_com, only: com_ext,com_red
use type_geom, only: fld_com_gl,fld_com_lg
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype
use type_randgen, only: rng,rand_real
use type_timer, only: timertype,timer_start,timer_end

implicit none

private
public :: test_adjoints,test_pos_def,test_mpi,test_dirac,test_perf

contains

!----------------------------------------------------------------------
! Subroutine: test_adjoints
!> Purpose: test adjoints accuracy
!----------------------------------------------------------------------
subroutine test_adjoints(nam,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(in) :: ndata !< Sampling data

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real) :: alpha(ndata%ns),alpha_save(ndata%ns)
real(kind_real) :: alpha1(ndata%ns),alpha1_save(ndata%ns)
real(kind_real) :: alpha2(ndata%ns),alpha2_save(ndata%ns)
real(kind_real) :: fld(ndata%nc0,ndata%nl0),fld_save(ndata%nc0,ndata%nl0)
real(kind_real) :: fld1(ndata%nc0,ndata%nl0),fld1_save(ndata%nc0,ndata%nl0)
real(kind_real) :: fld2(ndata%nc0,ndata%nl0),fld2_save(ndata%nc0,ndata%nl0)

! Initialization
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,alpha_save)
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,fld_save)

! Adjoint test
call interp(nam,ndata,alpha_save,fld)
call interp_ad(nam,ndata,fld_save,alpha)

! Print result
sum1 = sum(alpha*alpha_save)
sum2 = sum(fld*fld_save)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Initialization
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,alpha1_save)
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
call convol(ndata,alpha1)
call convol(ndata,alpha2)

! Print result
sum1 = sum(alpha1*alpha2_save)
sum2 = sum(alpha2*alpha1_save)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Convolution adjoint test:   ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Initialization
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,fld1_save)
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,fld2_save)
fld1 = fld1_save
fld2 = fld2_save

! Adjoint test
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(nam,ndata,fld1)
   call apply_nicas_from_sqrt(nam,ndata,fld2)
else
   call apply_nicas(nam,ndata,fld1)
   call apply_nicas(nam,ndata,fld2)
end if

! Print result
sum1 = sum(fld1*fld2_save)
sum2 = sum(fld2*fld1_save)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:         ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

end subroutine test_adjoints

!----------------------------------------------------------------------
! Subroutine: test_pos_def
!> Purpose: test positive_definiteness
!----------------------------------------------------------------------
subroutine test_pos_def(nam,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(in) :: ndata !< Sampling data

! Local variables
real(kind_real),parameter :: tol = 1.0e-3
integer,parameter :: nitermax = 1
integer :: iter
real(kind_real) :: norm,egvmax,egvmax_prev,egvmin,egvmin_prev
real(kind_real) :: fld(ndata%nc0,ndata%nl0),fld_prev(ndata%nc0,ndata%nl0)

! Power method to find the largest eigenvalue
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,fld_prev)
norm = sum(fld_prev**2)
fld_prev = fld_prev/norm
egvmax_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld = fld_prev

   ! Apply C
   if (nam%lsqrt) then
      call apply_nicas_from_sqrt(nam,ndata,fld)
   else
      call apply_nicas(nam,ndata,fld)
   end if

   ! Compute Rayleigh quotient
   egvmax = sum(fld*fld_prev)/sum(fld_prev*fld_prev)

   ! Renormalize the vector
   norm = sum(fld**2)
   fld = fld/norm

   ! Exit test
   if (abs(egvmax-egvmax_prev)<tol) exit

   ! Update
   iter = iter+1
   fld_prev = fld
   egvmax_prev = egvmax
end do

! Power method to find the smallest eigenvalue
call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,fld_prev)
norm = sum(fld_prev**2)
egvmin_prev = huge(1.0)
fld_prev = fld_prev/norm
egvmin_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld = fld_prev

   ! Apply C
   if (nam%lsqrt) then
      call apply_nicas_from_sqrt(nam,ndata,fld)
   else
      call apply_nicas(nam,ndata,fld)
   end if
   fld = fld-egvmax*fld_prev

   ! Compute Rayleigh quotient
   egvmin = sum(fld*fld_prev)/sum(fld_prev*fld_prev)

   ! Renormalize the vector
   norm = sum(fld**2)
   fld = fld/norm

   ! Exit test
   if (egvmax+egvmin<-tol*egvmax) then
      write(mpl%unit,'(a7,a)') '','NICAS is not positive definite'
      exit
   end if

   ! Update
   iter = iter+1
   fld_prev = fld
   egvmin_prev = egvmin
end do

! Non conclusive test
if (iter==nitermax+1) write(mpl%unit,'(a7,a,e15.8,a,i4,a,e15.8)') '','NICAS seems to be positive definite: difference ', &
 & egvmax+egvmin,' after ',nitermax,' iterations for a tolerance ',tol

end subroutine test_pos_def

!----------------------------------------------------------------------
! Subroutine: test_mpi
!> Purpose: test global/local equivalence
!----------------------------------------------------------------------
subroutine test_mpi(nam,ndata,ndataloc)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(in) :: ndata       !< Sampling data
type(ndataloctype),intent(in) :: ndataloc !< Sampling data, local

! Local variables
real(kind_real),allocatable :: fld(:,:),fldloc(:,:)

! Allocation
if (mpl%main) then
   ! Allocation
   allocate(fld(ndata%nc0,ndata%nl0))
   allocate(fldloc(ndata%nc0,ndata%nl0))

   ! Initialization
   call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,fld)
   fldloc = fld
end if

! Global to local
call fld_com_gl(ndata%geom,fldloc)

if (nam%lsqrt) then
   ! Global
   if (mpl%main) call apply_nicas_from_sqrt(nam,ndata,fld)

   ! Local
   call apply_nicas_from_sqrt(nam,ndataloc,fldloc)
else
   ! Global
   if (mpl%main) call apply_nicas(nam,ndata,fld)

   ! Local
   call apply_nicas(nam,ndataloc,fldloc)
end if

! Local to global
call fld_com_lg(ndata%geom,fldloc)

! Print difference
if (mpl%main) write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','RMSE for single-proc and multi-procs executions: ', &
 & sqrt(sum(fld**2)/float(ndata%nc0*ndata%nl0)),' / ',sqrt(sum(fldloc**2)/float(ndata%nc0*ndata%nl0)), &
 & ' / ',sqrt(sum((fld-fldloc)**2)/float(ndata%nc0*ndata%nl0))

end subroutine test_mpi

!----------------------------------------------------------------------
! Subroutine: test_dirac
!> Purpose: apply NICAS to Diracs
!----------------------------------------------------------------------
subroutine test_dirac(nam,ndata,ndataloc)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(in) :: ndata          !< Sampling data
type(ndataloctype),intent(inout) :: ndataloc !< Sampling data, local

! Local variables
integer :: il0,il0dir,ic0dir(nam%ndir),idir
real(kind_real),allocatable :: fld(:,:)

! Find level index
call msi(il0dir)
do il0=1,ndata%nl0
   if (nam%levs(il0)==nam%levdir) il0dir = il0
end do

if (mpl%main) then
   ! Setup dirac field
   call setup_dirac_points(ndata%nc0,ndata%geom%lon,ndata%geom%lat,ndata%geom%mask(:,il0dir),nam%ndir,nam%londir,nam%latdir,ic0dir)

   ! Allocation
   allocate(fld(ndata%nc0,ndata%nl0))

   ! Generate diracs field
   fld = 0.0
   do idir=1,nam%ndir
      fld(ic0dir(idir),il0dir) = 1.0
   end do
end if

! Global to local
call fld_com_gl(ndata%geom,fld)

! Apply NICAS method
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(nam,ndataloc,fld)
else
   call apply_nicas(nam,ndataloc,fld)
end if

! Local to global
call fld_com_lg(ndata%geom,fld)

if (mpl%main) then
   ! Write field
   call system('rm -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_dirac.nc')
   call model_write(nam,trim(nam%prefix)//'_dirac.nc','hr',ndata%geom,fld)

   ! Print results
   write(mpl%unit,'(a7,a)') '','Normalization:'
   do idir=1,nam%ndir
      write(mpl%unit,'(a10,f6.1,a,f6.1,a,f10.7)') '',nam%londir(idir),' / ',nam%latdir(idir),': ',fld(ic0dir(idir),il0dir)
   end do
   write(mpl%unit,'(a7,a,f10.7,a,f10.7)') '','Min - max: ', &
 & minval(fld(:,il0dir),mask=ndata%geom%mask(:,il0dir)),' - ',maxval(fld(:,il0dir),mask=ndata%geom%mask(:,il0dir))
end if

end subroutine test_dirac

!----------------------------------------------------------------------
! Subroutine: test_perf
!> Purpose: test NICAS performance
!----------------------------------------------------------------------
subroutine test_perf(nam,ndataloc)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data

! Local variables
real(kind_real) :: fld(ndataloc%nc0a,ndataloc%nl0)
real(kind_real),allocatable :: alpha(:),alpha_tmp(:)
type(timertype) :: timer_interp_ad,timer_com_1,timer_convol,timer_com_2,timer_interp

! Allocation
allocate(alpha(ndataloc%nsb))

! Random initialization
call rand_real(rng,0.0_kind_real,1.0_kind_real,.true.,fld)

! Adjoint interpolation
call timer_start(timer_interp_ad)
call interp_ad(nam,ndataloc,fld,alpha)
call timer_end(timer_interp_ad)

! Communication
call timer_start(timer_com_1)
if (nam%mpicom==1) then
   ! Allocation 
   allocate(alpha_tmp(ndataloc%nsb))

   ! Copy zone B
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndataloc%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone B into zone C
   alpha(ndataloc%isb_to_isc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
elseif (nam%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call com_red(ndataloc%AB,alpha)

   ! Allocation 
   allocate(alpha_tmp(ndataloc%nsb))

   ! Copy zone A
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndataloc%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone A into zone C
   alpha(ndataloc%isa_to_isc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
end if
call timer_end(timer_com_1)

! Convolution
call timer_start(timer_convol)
call convol(ndataloc,alpha)
call timer_start(timer_convol)

call timer_start(timer_com_2)
! Halo reduction from zone C to zone A
call com_red(ndataloc%AC,alpha)

! Halo extension from zone A to zone B
call com_ext(ndataloc%AB,alpha)
call timer_end(timer_com_2)

! Interpolation
call timer_start(timer_interp)
call interp(nam,ndataloc,alpha,fld)
call timer_end(timer_interp)

! Release memory
deallocate(alpha)

! Print results
write(mpl%unit,'(a7,a)') '','Performance results (elapsed time):'
write(mpl%unit,'(a10,a,f6.1,a)') '','Adjoint interpolation: ',timer_interp_ad%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Communication - 1    : ',timer_com_1%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Convolution          : ',timer_convol%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Communication - 2    : ',timer_com_2%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Interpolation        : ',timer_interp%elapsed,' s'

end subroutine test_perf

end module module_test
