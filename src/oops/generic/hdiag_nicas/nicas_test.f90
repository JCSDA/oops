!----------------------------------------------------------------------
! Module: nicas_test.f90
!> Purpose: test NICAS
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_test

use driver_hdiag, only: run_hdiag
use model_interface, only: model_write
use nicas_apply_bens, only: apply_bens
use nicas_apply_convol, only: apply_convol
use nicas_apply_interp, only: apply_interp,apply_interp_ad
use nicas_apply_localization, only: apply_localization,apply_localization_from_sqrt,randomize_localization
use nicas_apply_nicas, only: apply_nicas,apply_nicas_sqrt,apply_nicas_sqrt_ad,apply_nicas_from_sqrt
use omp_lib
use tools_const, only: deg2rad,rad2deg,sphere_dist
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr
use type_bdata, only: bdatatype
use type_bpar, only: bpartype
use type_com, only: com_ext,com_red
use type_ctree, only: find_nearest_neighbors
use type_geom, only: geomtype,fld_com_gl,fld_com_lg
use type_mpl, only: mpl,mpl_dot_prod
use type_nam, only: namtype
use type_ndata, only: ndatatype
use type_randgen, only: rand_real
use type_timer, only: timertype,timer_start,timer_end

implicit none

real(kind_real),parameter :: tol = 1.0e-3 !< Positive-definiteness test tolerance
integer,parameter :: nitermax = 50        !< Nunmber of iterations for the positive-definiteness test

private
public :: test_nicas_adjoints,test_nicas_pos_def,test_nicas_sqrt,test_nicas_dirac,test_nicas_perf
public :: test_loc_adjoint,test_loc_dirac,test_loc_ens_dirac,test_hdiag

contains

!----------------------------------------------------------------------
! Subroutine: test_nicas_adjoints
!> Purpose: test adjoints accuracy
!----------------------------------------------------------------------
subroutine test_nicas_adjoints(ndata)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< NICAS data

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real) :: alpha(ndata%nsb),alpha_save(ndata%nsb)
real(kind_real) :: alpha1(ndata%nsc),alpha1_save(ndata%nsc)
real(kind_real) :: alpha2(ndata%nsc),alpha2_save(ndata%nsc)
real(kind_real) :: fld(ndata%geom%nc0a,ndata%geom%nl0),fld_save(ndata%geom%nc0a,ndata%geom%nl0)
real(kind_real),allocatable :: fld1(:,:),fld1_save(:,:)
real(kind_real),allocatable :: fld2(:,:),fld2_save(:,:)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
call rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call apply_interp(geom,ndata,alpha_save,fld)
call apply_interp_ad(geom,ndata,fld_save,alpha)

! Print result
sum1 = sum(alpha*alpha_save)
sum2 = sum(fld*fld_save)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
call apply_convol(ndata,alpha1)
call apply_convol(ndata,alpha2)

! Print result
sum1 = sum(alpha1*alpha2_save)
sum2 = sum(alpha2*alpha1_save)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Convolution adjoint test:   ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

if (mpl%main) then
   ! Allocation
   allocate(fld1_save(geom%nc0,geom%nl0))
   allocate(fld2_save(geom%nc0,geom%nl0))

   ! Generate random field
   call rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
   call rand_real(0.0_kind_real,1.0_kind_real,fld2_save)
end if

! Global to local
call fld_com_gl(geom,fld1_save)
call fld_com_gl(geom,fld2_save)

! Allocation
allocate(fld1(ndata%geom%nc0a,ndata%geom%nl0))
allocate(fld2(ndata%geom%nc0a,ndata%geom%nl0))

! Adjoint test
fld1 = fld1_save
fld2 = fld2_save
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(geom,ndata,fld1)
   call apply_nicas_from_sqrt(geom,ndata,fld2)
else
   call apply_nicas(geom,ndata,fld1)
   call apply_nicas(geom,ndata,fld2)
end if

! Print result
call mpl_dot_prod(fld1,fld2_save,sum1)
call mpl_dot_prod(fld2,fld1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:         ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! End associate
end associate

end subroutine test_nicas_adjoints

!----------------------------------------------------------------------
! Subroutine: test_nicas_pos_def
!> Purpose: test positive_definiteness
!----------------------------------------------------------------------
subroutine test_nicas_pos_def(ndata)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< NICAS data

! Local variables
integer :: iter
real(kind_real) :: norm,num,den,egvmax,egvmax_prev,egvmin,egvmin_prev
real(kind_real) :: fld(ndata%geom%nc0a,ndata%geom%nl0),fld_prev(ndata%geom%nc0a,ndata%geom%nl0)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Power method to find the largest eigenvalue
call rand_real(0.0_kind_real,1.0_kind_real,fld_prev)
call mpl_dot_prod(fld_prev,fld_prev,norm)
fld_prev = fld_prev/norm
egvmax_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld = fld_prev

   ! Apply C
   if (nam%lsqrt) then
      call apply_nicas_from_sqrt(geom,ndata,fld)
   else
      call apply_nicas(geom,ndata,fld)
   end if

   ! Compute Rayleigh quotient
   call mpl_dot_prod(fld,fld_prev,num)
   call mpl_dot_prod(fld_prev,fld_prev,den)
   egvmax = num/den

   ! Renormalize the vector
   call mpl_dot_prod(fld,fld,norm)
   fld = fld/norm

   ! Exit test
   if (abs(egvmax-egvmax_prev)<tol) exit

   ! Update
   iter = iter+1
   fld_prev = fld
   egvmax_prev = egvmax
end do

! Power method to find the smallest eigenvalue
call rand_real(0.0_kind_real,1.0_kind_real,fld_prev)
call mpl_dot_prod(fld_prev,fld_prev,norm)
egvmin_prev = huge(1.0)
fld_prev = fld_prev/norm
egvmin_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld = fld_prev

   ! Apply C
   if (nam%lsqrt) then
      call apply_nicas_from_sqrt(geom,ndata,fld)
   else
      call apply_nicas(geom,ndata,fld)
   end if
   fld = fld-egvmax*fld_prev

   ! Compute Rayleigh quotient
   call mpl_dot_prod(fld,fld_prev,num)
   call mpl_dot_prod(fld_prev,fld_prev,den)
   egvmin = num/den

   ! Renormalize the vector
   call mpl_dot_prod(fld,fld,norm)
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

! End associate
end associate

end subroutine test_nicas_pos_def

!----------------------------------------------------------------------
! Subroutine: test_nicas_sqrt
!> Purpose: test full/square-root equivalence
!----------------------------------------------------------------------
subroutine test_nicas_sqrt(bdata,ndata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata !< B data
type(ndatatype),intent(in) :: ndata !< NICAS data

! Local variables
real(kind_real) :: fld(ndata%geom%nc0a,ndata%geom%nl0),fld_sqrt(ndata%geom%nc0a,ndata%geom%nl0)
type(ndatatype) :: ndata_other

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! TODO
call msgerror('not implemented yet')

! Print difference
write(mpl%unit,'(a7,a,f6.1,a)') '','Full / square-root error:',sqrt(sum((fld_sqrt-fld)**2)/sum(fld**2))*100.0,'%'

! End associate
end associate

end subroutine test_nicas_sqrt

!----------------------------------------------------------------------
! Subroutine: test_nicas_dirac
!> Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine test_nicas_dirac(nam,geom,blockname,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam          !< Namelist
type(geomtype),intent(in) :: geom        !< Geometry
character(len=*),intent(in) :: blockname !< Block name
type(ndatatype),intent(in) :: ndata      !< NICAS data

! Local variables
integer :: ic0,ic1,il0,il1,is
integer :: il0dir(nam%ndir),ic0dir(nam%ndir),idir
real(kind_real) :: dum(1)
real(kind_real),allocatable :: fld(:,:),fld_c1(:,:),fld_s(:,:)

if (mpl%main) then
   ! Allocation
   allocate(fld(geom%nc0,geom%nl0))

   ! Generate diracs field
   fld = 0.0
   do idir=1,nam%ndir
      ! Find level index
      do il0=1,geom%nl0
          if (nam%levs(il0)==nam%levdir(idir)) il0dir(idir) = il0
      end do

      ! Find nearest neighbor
      call find_nearest_neighbors(geom%ctree,dble(nam%londir(idir)*deg2rad),dble(nam%latdir(idir)*deg2rad),1,ic0dir(idir:idir),dum)

      ! Dirac value
      fld(ic0dir(idir),il0dir(idir)) = 1.0
   end do
end if

! Global to local
call fld_com_gl(geom,fld)

! Apply NICAS method
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(geom,ndata,fld)
else
   call apply_nicas(geom,ndata,fld)
end if

! Local to global
call fld_com_lg(geom,fld)

if (mpl%main) then
   ! Write field
   call model_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(blockname)//'_dirac',fld)

   ! Print results
   do idir=1,nam%ndir
      write(mpl%unit,'(a10,f6.1,a,f6.1,a,f10.7)') '',nam%londir(idir),' / ',nam%latdir(idir),': ',fld(ic0dir(idir),il0dir(idir))
   end do
   write(mpl%unit,'(a7,a,f10.7,a,f10.7)') '','Min - max: ', &
 & minval(fld(:,il0dir),mask=geom%mask(:,il0dir)),' - ',maxval(fld(:,il0dir),mask=geom%mask(:,il0dir))

   if (nam%new_param) then
      ! Write field at interpolation points
      allocate(fld_c1(geom%nc0,geom%nl0))
      call msr(fld_c1)
      do ic1=1,ndata%nc1
         ic0 = ndata%c1_to_c0(ic1)
         fld_c1(ic0,:) = fld(ic0,:)
      end do
      call model_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(blockname)//'_dirac_c1',fld_c1)

      ! Write field at subgrid points
      allocate(fld_s(geom%nc0,geom%nl0))
      call msr(fld_s)
      do is=1,ndata%ns
         ic1 = ndata%s_to_c1(is)
         ic0 = ndata%c1_to_c0(ic1)
         il1 = ndata%s_to_l1(is)
         il0 = ndata%l1_to_l0(il1)
         fld_s(ic0,il0) = fld(ic0,il0)
      end do
      call model_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(blockname)//'_dirac_s',fld_s)
   end if
end if

end subroutine test_nicas_dirac

!----------------------------------------------------------------------
! Subroutine: test_nicas_perf
!> Purpose: test NICAS performance
!----------------------------------------------------------------------
subroutine test_nicas_perf(nam,geom,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam     !< Namelist
type(geomtype),intent(in) :: geom   !< Geometry
type(ndatatype),intent(in) :: ndata !< NICAS data

! Local variables
real(kind_real),allocatable :: fld(:,:),alpha(:),alpha_tmp(:)
type(timertype) :: timer_interp_ad,timer_com_1,timer_convol,timer_com_2,timer_interp

! Allocation
allocate(alpha(ndata%nsb))

if (mpl%main) then
   ! Allocation
   allocate(fld(geom%nc0,geom%nl0))

   ! Generate random field
   call rand_real(0.0_kind_real,1.0_kind_real,fld)
end if

! Global to local
call fld_com_gl(geom,fld)

! Adjoint interpolation
call timer_start(timer_interp_ad)
call apply_interp_ad(geom,ndata,fld,alpha)
call timer_end(timer_interp_ad)

! Communication
call timer_start(timer_com_1)
if (ndata%mpicom==1) then
   ! Allocation
   allocate(alpha_tmp(ndata%nsb))

   ! Copy zone B
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndata%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone B into zone C
   alpha(ndata%sb_to_sc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
elseif (ndata%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call com_red(ndata%AB,alpha)

   ! Allocation
   allocate(alpha_tmp(ndata%nsb))

   ! Copy zone A
   alpha_tmp = alpha

   ! Reallocation
   deallocate(alpha)
   allocate(alpha(ndata%nsc))

   ! Initialize
   alpha = 0.0

   ! Copy zone A into zone C
   alpha(ndata%sa_to_sc) = alpha_tmp

   ! Release memory
   deallocate(alpha_tmp)
end if
call timer_end(timer_com_1)

! Convolution
call timer_start(timer_convol)
call apply_convol(ndata,alpha)
call timer_start(timer_convol)

call timer_start(timer_com_2)
! Halo reduction from zone C to zone A
call com_red(ndata%AC,alpha)

! Halo extension from zone A to zone B
call com_ext(ndata%AB,alpha)
call timer_end(timer_com_2)

! Interpolation
call timer_start(timer_interp)
call apply_interp(geom,ndata,alpha,fld)
call timer_end(timer_interp)

! Release memory
deallocate(alpha)

! Print results
write(mpl%unit,'(a10,a,f6.1,a)') '','Adjoint interpolation: ',timer_interp_ad%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Communication - 1    : ',timer_com_1%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Convolution          : ',timer_convol%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Communication - 2    : ',timer_com_2%elapsed,' s'
write(mpl%unit,'(a10,a,f6.1,a)') '','Interpolation        : ',timer_interp%elapsed,' s'

end subroutine test_nicas_perf

!----------------------------------------------------------------------
! Subroutine: test_loc_adjoint
!> Purpose: test localization adjoint
!----------------------------------------------------------------------
subroutine test_loc_adjoint(nam,geom,bpar,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                !< Namelist
type(geomtype),intent(in) :: geom              !< Geometry
type(bpartype),intent(in) :: bpar              !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1) !< NICAS data

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real),allocatable :: fld1(:,:,:,:),fld1_save(:,:,:,:)
real(kind_real),allocatable :: fld2(:,:,:,:),fld2_save(:,:,:,:)

if (mpl%main) then
   ! Allocation
   allocate(fld1_save(geom%nc0,geom%nl0,nam%nv,nam%nts))
   allocate(fld2_save(geom%nc0,geom%nl0,nam%nv,nam%nts))

   ! Generate random field
   call rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
   call rand_real(0.0_kind_real,1.0_kind_real,fld2_save)
end if

! Global to local
call fld_com_gl(nam,geom,fld1_save)
call fld_com_gl(nam,geom,fld2_save)

! Allocation
allocate(fld1(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld2(geom%nc0a,geom%nl0,nam%nv,nam%nts))

! Adjoint test
fld1 = fld1_save
fld2 = fld2_save
if (nam%lsqrt) then
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld1)
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld2)
else
   call apply_localization(nam,geom,bpar,ndata,fld1)
   call apply_localization(nam,geom,bpar,ndata,fld2)
end if

! Print result
call mpl_dot_prod(fld1,fld2_save,sum1)
call mpl_dot_prod(fld2,fld1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Localization adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

end subroutine test_loc_adjoint

!----------------------------------------------------------------------
! Subroutine: test_loc_dirac
!> Purpose: apply localization to diracs
!----------------------------------------------------------------------
subroutine test_loc_dirac(nam,geom,bpar,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                !< Namelist
type(geomtype),intent(in) :: geom              !< Geometry
type(bpartype),intent(in) :: bpar              !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1) !< NICAS data

! Local variables
integer :: il0,il0dir(nam%ndir),ic0dir(nam%ndir),idir,iv,its
real(kind_real) :: dum(1)
real(kind_real),allocatable :: dirac(:,:,:,:)
character(len=2) :: itschar

if (mpl%main) then
   ! Allocation
   allocate(dirac(geom%nc0,geom%nl0,nam%nv,nam%nts))

   ! Generate diracs field
   dirac = 0.0
   do idir=1,nam%ndir
      ! Find level index
      do il0=1,geom%nl0
          if (nam%levs(il0)==nam%levdir(idir)) il0dir(idir) = il0
      end do

      ! Find nearest neighbor
      call find_nearest_neighbors(geom%ctree,dble(nam%londir(idir)*deg2rad),dble(nam%latdir(idir)*deg2rad),1,ic0dir(idir:idir),dum)

      ! Dirac value
      dirac(ic0dir(idir),il0dir(idir),nam%ivdir(idir),nam%itsdir(idir)) = 1.0
   end do
end if

! Global to local
call fld_com_gl(nam,geom,dirac)

! Apply localization to dirac
if (nam%lsqrt) then
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,dirac)
else
   call apply_localization(nam,geom,bpar,ndata,dirac)
end if

! Local to global
call fld_com_lg(nam,geom,dirac)

if (mpl%main) then
   ! Write field
   do its=1,nam%nts
      write(itschar,'(i2.2)') its
      do iv=1,nam%nv
         call model_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(nam%varname(iv))//'_'//itschar,dirac(:,:,iv,its))
      end do
   end do
end if

end subroutine test_loc_dirac

!----------------------------------------------------------------------
! Subroutine: test_loc_ens_dirac
!> Purpose: apply localized ensemble to diracs
!----------------------------------------------------------------------
subroutine test_loc_ens_dirac(nam,geom,bpar,ndata,ens1)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                                                   !< Namelist
type(geomtype),intent(in) :: geom                                                 !< Geometry
type(bpartype),intent(in) :: bpar                                                 !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                                    !< NICAS data
real(kind_real),intent(in) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: il0,il0dir(nam%ndir),ic0dir(nam%ndir),idir,iv,its
real(kind_real) :: dum(1)
real(kind_real),allocatable :: dirac(:,:,:,:)
character(len=2) :: itschar

if (mpl%main) then
   ! Allocation
   allocate(dirac(geom%nc0,geom%nl0,nam%nv,nam%nts))

   ! Generate diracs field
   dirac = 0.0
   do idir=1,nam%ndir
      ! Find level index
      do il0=1,geom%nl0
          if (nam%levs(il0)==nam%levdir(idir)) il0dir(idir) = il0
      end do

      ! Find nearest neighbor
      call find_nearest_neighbors(geom%ctree,dble(nam%londir(idir)*deg2rad),dble(nam%latdir(idir)*deg2rad),1,ic0dir(idir:idir),dum)

      ! Dirac value
      dirac(ic0dir(idir),il0dir(idir),nam%ivdir(idir),nam%itsdir(idir)) = 1.0
   end do
end if

! Global to local
call fld_com_gl(nam,geom,dirac)

! Apply localized ensemble covariance
call apply_bens(nam,geom,bpar,ndata,ens1,dirac)

! Local to global
call fld_com_lg(nam,geom,dirac)

if (mpl%main) then
   ! Write field
   do its=1,nam%nts
      write(itschar,'(i2.2)') its
      do iv=1,nam%nv
         call model_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(nam%varname(iv))//'_'//itschar//'_Bens',dirac(:,:,iv,its))
      end do
   end do
end if

end subroutine test_loc_ens_dirac

!----------------------------------------------------------------------
! Subroutine: test_hdiag
!> Purpose: test hdiag with a randomization method
!----------------------------------------------------------------------
subroutine test_hdiag(nam,geom,bpar,bdata,ndata)

implicit none

! Passed variables
type(namtype),intent(inout) :: nam             !< Namelist variables
type(geomtype),intent(in) :: geom              !< Geometry
type(bpartype),intent(in) :: bpar              !< Block parameters
type(bdatatype),intent(in) :: bdata(bpar%nb+1) !< B data
type(ndatatype),intent(in) :: ndata(bpar%nb+1) !< NICAS data

! Local variables
integer,parameter :: ne = 150
integer :: ens1_ne,ens1_ne_offset,ens1_nsub,ib,il0
real(kind_real) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne)
logical :: new_hdiag,displ_diag,gau_approx,full_var
character(len=1024) :: prefix,method
type(bdatatype),allocatable :: bdata_test(:)

! Randomize ensemble
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Randomize ensemble'
call randomize_localization(nam,geom,bpar,ndata,ne,ens1)

! Copy sampling
call system('cp -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc ' &
 & //trim(nam%datadir)//'/'//trim(nam%prefix)//'_hdiag-test_sampling.nc')

! Save namelist variables
prefix = nam%prefix
method = nam%method
new_hdiag = nam%new_hdiag
displ_diag = nam%displ_diag
gau_approx = nam%gau_approx
full_var = nam%full_var
ens1_ne = nam%ens1_ne
ens1_ne_offset = nam%ens1_ne_offset
ens1_nsub = nam%ens1_nsub

! Set namelist variables
nam%prefix = trim(nam%prefix)//'_hdiag-test'
nam%method = 'cor'
nam%new_hdiag = .true.
nam%displ_diag = .false.
nam%gau_approx = .true.
nam%full_var = .false.
nam%ens1_ne = ne
nam%ens1_ne_offset = 0
nam%ens1_nsub = 1

! Call hdiag driver
call run_hdiag(nam,geom,bpar,bdata_test,ens1)

! Print scores
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- hdiag consistency results'
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      do il0=1,geom%nl0
         write(mpl%unit,'(a10,a7,i3,a4,a25,f6.1,a)') '','Level: ',nam%levs(il0),' ~> ','ensemble coefficient: ', &
       & sqrt(sum((bdata_test(ib)%coef_ens(:,il0)-bdata(ib)%coef_ens(:,il0))**2)/sum(bdata(ib)%coef_ens(:,il0)**2))*100.0,'%'
         write(mpl%unit,'(a49,f6.1,a)') 'horizontal length-scale: ', &
       & sqrt(sum((bdata_test(ib)%rh0(:,il0)-bdata(ib)%rh0(:,il0))**2)/sum(bdata(ib)%rh0(:,il0)**2))*100.0,'%'
         if (any(abs(bdata(ib)%rv0(:,il0))>0.0)) then
            write(mpl%unit,'(a49,f6.1,a)') 'vertical length-scale: ', &
          & sqrt(sum((bdata_test(ib)%rv0(:,il0)-bdata(ib)%rv0(:,il0))**2)/sum(bdata(ib)%rv0(:,il0)**2))*100.0,'%'
         end if
      end do
   end if
end do

! Reset namelist variables
nam%prefix = prefix
nam%method = method
nam%new_hdiag = new_hdiag
nam%displ_diag = displ_diag
nam%gau_approx = gau_approx
nam%full_var = full_var
nam%ens1_ne = ens1_ne
nam%ens1_ne_offset = ens1_ne_offset
nam%ens1_nsub = ens1_nsub

end subroutine test_hdiag

end module nicas_test
