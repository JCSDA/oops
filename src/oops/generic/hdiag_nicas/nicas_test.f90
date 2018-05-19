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
use nicas_apply_bens, only: apply_bens,apply_bens_noloc
use nicas_apply_convol, only: apply_convol
use nicas_apply_interp, only: apply_interp,apply_interp_s,apply_interp_v,apply_interp_h, &
                            & apply_interp_ad,apply_interp_h_ad,apply_interp_v_ad,apply_interp_s_ad
use nicas_apply_localization, only: apply_localization,apply_localization_from_sqrt,randomize_localization
use nicas_apply_nicas, only: apply_nicas,apply_nicas_sqrt,apply_nicas_sqrt_ad,apply_nicas_from_sqrt
use nicas_parameters, only: compute_parameters
use omp_lib
use tools_const, only: deg2rad,rad2deg,sphere_dist,reqkm
use tools_display, only: msgerror,vunitchar
use tools_kinds,only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr
use type_bdata, only: bdatatype,bdata_alloc,bdata_copy
use type_bpar, only: bpartype
use type_com, only: com_ext,com_red
use type_ctree, only: find_nearest_neighbors
use type_geom, only: geomtype,fld_com_gl,fld_com_lg,fld_write
use type_linop, only: apply_linop,apply_linop_ad
use type_mpl, only: mpl,mpl_dot_prod,mpl_bcast
use type_nam, only: namtype
use type_ndata, only: ndatatype,ndata_dealloc
use type_randgen, only: rand_real,rand_integer
use type_timer, only: timertype,timer_start,timer_end

implicit none

real(kind_real),parameter :: tol = 1.0e-3 !< Positive-definiteness test tolerance
integer,parameter :: nitermax = 50        !< Number of iterations for the positive-definiteness test
integer,parameter :: ne = 150             !< Ensemble size for randomization
integer,parameter :: nfac = 10            !< Number of length-scale factors
integer,parameter :: ntest = 100          !< Number of tests

private
public :: test_nicas_adjoint,test_nicas_pos_def,test_nicas_sqrt,test_nicas_dirac
public :: test_loc_adjoint,test_loc_sqrt,test_loc_dirac
public :: test_randomization,test_consistency,test_optimality

contains

!----------------------------------------------------------------------
! Subroutine: test_nicas_adjoint
!> Purpose: test NICAS adjoint accuracy
!----------------------------------------------------------------------
subroutine test_nicas_adjoint(ndata)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< NICAS data

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real),allocatable :: alpha(:),alpha_save(:),alpha1(:),alpha1_save(:),alpha2(:),alpha2_save(:)
real(kind_real),allocatable :: gamma(:,:),gamma_save(:,:),delta(:,:),delta_save(:,:)
real(kind_real),allocatable :: fld(:,:),fld_save(:,:),fld1(:,:),fld1_save(:,:),fld2(:,:),fld2_save(:,:)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Allocation
allocate(alpha(ndata%nsb))
allocate(alpha_save(ndata%nsb))
allocate(gamma(ndata%nc1b,ndata%nl1))
allocate(gamma_save(ndata%nc1b,ndata%nl1))
allocate(delta(ndata%nc1b,geom%nl0))
allocate(delta_save(ndata%nc1b,geom%nl0))
allocate(fld(geom%nc0a,ndata%geom%nl0))
allocate(fld_save(geom%nc0a,ndata%geom%nl0))

! Interpolation (subsampling)

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
call rand_real(0.0_kind_real,1.0_kind_real,gamma_save)

! Adjoint test
call apply_interp_s(ndata,alpha_save,gamma)
call apply_interp_s_ad(ndata,gamma_save,alpha)

! Print result
call mpl_dot_prod(alpha,alpha_save,sum1)
call mpl_dot_prod(gamma,gamma_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (subsampling): ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Interpolation (vertical)

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,gamma_save)
call rand_real(0.0_kind_real,1.0_kind_real,delta_save)

! Adjoint test
call apply_interp_v(geom,ndata,gamma_save,delta)
call apply_interp_v_ad(geom,ndata,delta_save,gamma)

! Print result
call mpl_dot_prod(gamma,gamma_save,sum1)
call mpl_dot_prod(delta,delta_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (vertical):    ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Interpolation (horizontal)

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,delta_save)
call rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call apply_interp_h(geom,ndata,delta_save,fld)
call apply_interp_h_ad(geom,ndata,fld_save,delta)

! Print result
call mpl_dot_prod(delta,delta_save,sum1)
call mpl_dot_prod(fld,fld_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (horizontal):  ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Interpolation (total)

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
call rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call apply_interp(geom,ndata,alpha_save,fld)
call apply_interp_ad(geom,ndata,fld_save,alpha)

! Print result
call mpl_dot_prod(alpha,alpha_save,sum1)
call mpl_dot_prod(fld,fld_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (total):       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Allocation
allocate(alpha1(ndata%nsc))
allocate(alpha1_save(ndata%nsc))
allocate(alpha2(ndata%nsc))
allocate(alpha2_save(ndata%nsc))

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
call apply_convol(ndata,alpha1)
call apply_convol(ndata,alpha2)

! Print result
call mpl_dot_prod(alpha1,alpha2_save,sum1)
call mpl_dot_prod(alpha2,alpha1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Convolution adjoint test:                 ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Allocation
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)
allocate(alpha1(ndata%nsb))
allocate(alpha1_save(ndata%nsb))
allocate(alpha2(ndata%nsa))
allocate(alpha2_save(ndata%nsa))

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
call com_red(ndata%AB,alpha1)
call com_ext(ndata%AB,alpha2)

! Print result
call mpl_dot_prod(alpha1,alpha2_save,sum1)
call mpl_dot_prod(alpha2,alpha1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AB adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Allocation
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)
allocate(alpha1(ndata%nsc))
allocate(alpha1_save(ndata%nsc))
allocate(alpha2(ndata%nsa))
allocate(alpha2_save(ndata%nsa))

! Initialization
call rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
call com_red(ndata%AC,alpha1)
call com_ext(ndata%AC,alpha2)

! Print result
call mpl_dot_prod(alpha1,alpha2_save,sum1)
call mpl_dot_prod(alpha2,alpha1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AC adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Allocation
allocate(fld1(geom%nc0a,ndata%geom%nl0))
allocate(fld2(geom%nc0a,ndata%geom%nl0))
allocate(fld1_save(geom%nc0,geom%nl0))
allocate(fld2_save(geom%nc0,geom%nl0))

! Generate random field
call rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rand_real(0.0_kind_real,1.0_kind_real,fld2_save)
fld1 = fld1_save
fld2 = fld2_save

! Adjoint test
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
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:                       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! End associate
end associate

end subroutine test_nicas_adjoint

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
subroutine test_nicas_sqrt(blockname,bdata,ndata)

implicit none

! Passed variables
character(len=*),intent(in) :: blockname !< Block name
type(bdatatype),intent(in) :: bdata      !< B data
type(ndatatype),intent(in) :: ndata      !< NICAS data

! Local variables
real(kind_real) :: fld(ndata%geom%nc0a,ndata%geom%nl0),fld_sqrt(ndata%geom%nc0a,ndata%geom%nl0)
type(ndatatype) :: ndata_other

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Set namelist and geometry
ndata_other%nam => ndata%nam
ndata_other%geom => ndata%geom

! Generate random field
call rand_real(-1.0_kind_real,1.0_kind_real,fld)
fld_sqrt = fld

! Apply NICAS, initial version
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(geom,ndata,fld_sqrt)
else
   call apply_nicas(geom,ndata,fld)
end if

! Switch lsqrt
nam%lsqrt = .not.nam%lsqrt

! Compute NICAS parameters
call compute_parameters(bdata,ndata_other)

! Apply NICAS, other version
if (nam%lsqrt) then
   call apply_nicas_from_sqrt(geom,ndata_other,fld_sqrt)
else
   ! Apply NICAS
   call apply_nicas(geom,ndata_other,fld)
end if

! Compute dirac
if (nam%check_dirac) call test_nicas_dirac(trim(blockname)//'_sqrt',ndata_other)

! Reset lsqrt value
nam%lsqrt = .not.nam%lsqrt

! Print difference
write(mpl%unit,'(a7,a,f6.1,a)') '','NICAS full / square-root error:',sqrt(sum((fld_sqrt-fld)**2)/sum(fld**2))*100.0,'%'

! End associate
end associate

end subroutine test_nicas_sqrt

!----------------------------------------------------------------------
! Subroutine: test_nicas_dirac
!> Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine test_nicas_dirac(blockname,ndata)

implicit none

! Passed variables
character(len=*),intent(in) :: blockname !< Block name
type(ndatatype),intent(in) :: ndata      !< NICAS data

! Local variables
integer :: ic0,ic1,il0,il1,is
integer :: il0dir(ndata%nam%ndir),ic0dir(ndata%nam%ndir),idir
real(kind_real),allocatable :: fld(:,:),fld_c1(:,:),fld_s(:,:)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

if (mpl%main) then
   ! Allocation
   allocate(fld(geom%nc0,geom%nl0))

   ! Find gridpoint and level indices
   call define_dirac(nam,geom,ic0dir,il0dir)

   ! Generate dirac field
   fld = 0.0
   do idir=1,nam%ndir
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
   call model_write(nam,geom,trim(nam%prefix)//'_dirac_gridded.nc',trim(blockname)//'_dirac',fld)
   call fld_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(blockname)//'_dirac',fld)

   ! Print results
   write(mpl%unit,'(a7,a)') '','Values at dirac points:'
   do idir=1,nam%ndir
      write(mpl%unit,'(a10,f6.1,a,f6.1,a,f10.7)') '',nam%londir(idir),' / ',nam%latdir(idir),': ',fld(ic0dir(idir),il0dir(idir))
   end do
   write(mpl%unit,'(a7,a)') '','Min - max: '
   do il0=1,geom%nl0
      write(mpl%unit,'(a10,a,i3,a,f10.7,a,f10.7)') '','Level ',nam%levs(il0),': ', &
    & minval(fld(:,il0),mask=geom%mask(:,il0)),' - ',maxval(fld(:,il0),mask=geom%mask(:,il0))
   end do

   if (nam%new_param) then
      ! Write field at interpolation points
      allocate(fld_c1(geom%nc0,geom%nl0))
      call msr(fld_c1)
      do ic1=1,ndata%nc1
         ic0 = ndata%c1_to_c0(ic1)
         fld_c1(ic0,:) = fld(ic0,:)
      end do
      call model_write(nam,geom,trim(nam%prefix)//'_dirac_gridded.nc',trim(blockname)//'_dirac_c1',fld_c1)
      call fld_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(blockname)//'_dirac_c1',fld_c1)

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
      call model_write(nam,geom,trim(nam%prefix)//'_dirac_gridded.nc',trim(blockname)//'_dirac_s',fld_s)
      call fld_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(blockname)//'_dirac_s',fld_s)
   end if
end if

! End associate
end associate

end subroutine test_nicas_dirac

!----------------------------------------------------------------------
! Subroutine: test_loc_adjoint
!> Purpose: test localization adjoint
!----------------------------------------------------------------------
subroutine test_loc_adjoint(nam,geom,bpar,ndata,ens1)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                                                            !< Namelist
type(geomtype),intent(in) :: geom                                                          !< Geometry
type(bpartype),intent(in) :: bpar                                                          !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                                             !< NICAS data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real),allocatable :: fld1_loc(:,:,:,:),fld1_bens(:,:,:,:),fld1_save(:,:,:,:)
real(kind_real),allocatable :: fld2_loc(:,:,:,:),fld2_bens(:,:,:,:),fld2_save(:,:,:,:)

! Allocation
allocate(fld1_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld2_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld1_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld2_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (present(ens1)) then
   allocate(fld1_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   allocate(fld2_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
end if

! Generate random field
call rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rand_real(0.0_kind_real,1.0_kind_real,fld2_save)

! Adjoint test
fld1_loc = fld1_save
fld2_loc = fld2_save
if (nam%lsqrt) then
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld1_loc)
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld2_loc)
else
   call apply_localization(nam,geom,bpar,ndata,fld1_loc)
   call apply_localization(nam,geom,bpar,ndata,fld2_loc)
end if
if (present(ens1)) then
   fld1_bens = fld1_save
   fld2_bens = fld2_save
   call apply_bens(nam,geom,bpar,ndata,ens1,fld1_bens)
   call apply_bens(nam,geom,bpar,ndata,ens1,fld2_bens)
end if

! Print result
call mpl_dot_prod(fld1_loc,fld2_save,sum1)
call mpl_dot_prod(fld2_loc,fld1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Localization adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
if (present(ens1)) then
   call mpl_dot_prod(fld1_bens,fld2_save,sum1)
   call mpl_dot_prod(fld2_bens,fld1_save,sum2)
   write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Ensemble B adjoint test:   ', &
    & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
end if

end subroutine test_loc_adjoint

!----------------------------------------------------------------------
! Subroutine: test_loc_sqrt
!> Purpose: test full/square-root equivalence
!----------------------------------------------------------------------
subroutine test_loc_sqrt(nam,geom,bpar,bdata,ndata,ens1)

implicit none

! Passed variables
type(namtype),intent(inout),target :: nam                                                  !< Namelist
type(geomtype),intent(in),target :: geom                                                   !< Geometry
type(bpartype),intent(in) :: bpar                                                          !< Block parameters
type(bdatatype),intent(in) :: bdata(bpar%nb+1)                                             !< B data
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                                             !< NICAS data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ib,ic0,ic0a,iv
real(kind_real),allocatable :: fld_loc(:,:,:,:),fld_loc_sqrt(:,:,:,:)
real(kind_real),allocatable :: fld_bens(:,:,:,:),fld_bens_sqrt(:,:,:,:)
character(len=1024) :: varname(nam%nv)
type(ndatatype) :: ndata_other(bpar%nb+1)

! Allocation
allocate(fld_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld_loc_sqrt(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (present(ens1)) then
   allocate(fld_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   allocate(fld_bens_sqrt(geom%nc0a,geom%nl0,nam%nv,nam%nts))
end if

! Generate random field
call rand_real(-1.0_kind_real,1.0_kind_real,fld_loc)
fld_loc_sqrt = fld_loc
if (present(ens1)) then
   call rand_real(-1.0_kind_real,1.0_kind_real,fld_bens)
   fld_bens_sqrt = fld_bens
end if

! Apply localization, initial version
if (nam%lsqrt) then
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld_loc_sqrt)
else
   call apply_localization(nam,geom,bpar,ndata,fld_loc)
end if
if (present(ens1)) then
   if (nam%lsqrt) then
      call apply_bens(nam,geom,bpar,ndata,ens1,fld_bens_sqrt)
   else
      call apply_bens(nam,geom,bpar,ndata,ens1,fld_bens)
   end if
end if

! Switch lsqrt
nam%lsqrt = .not.nam%lsqrt

! Prepare ndata, other version
do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      ! Set name, namelist and geometry
      ndata_other(ib)%nam => nam
      ndata_other(ib)%geom => geom
   end if

   if (bpar%nicas_block(ib)) then
      ! Compute NICAS parameters
      call compute_parameters(bdata(ib),ndata_other(ib))
   end if

   if (bpar%B_block(ib)) then
      ! Copy weights
      ndata_other(ib)%wgt = bdata(ib)%wgt
      if (bpar%nicas_block(ib)) then
         allocate(ndata_other(ib)%coef_ens(geom%nc0a,geom%nl0))
         do ic0=1,geom%nc0
            if (geom%c0_to_proc(ic0)==mpl%myproc) then
               ic0a = geom%c0_to_c0a(ic0)
               ndata_other(ib)%coef_ens(ic0a,:) = bdata(ib)%coef_ens(ic0,:)
            end if
         end do
      end if
   end if
end do

! Apply localization, other version
if (nam%lsqrt) then
   call apply_localization_from_sqrt(nam,geom,bpar,ndata_other,fld_loc_sqrt)
else
   call apply_localization(nam,geom,bpar,ndata_other,fld_loc)
end if
if (present(ens1)) then
   if (nam%lsqrt) then
      call apply_bens(nam,geom,bpar,ndata_other,ens1,fld_bens_sqrt)
   else
      call apply_bens(nam,geom,bpar,ndata_other,ens1,fld_bens)
   end if
end if

! Compute dirac
do iv=1,nam%nv
   varname(iv) = nam%varname(iv)
   nam%varname(iv) = trim(varname(iv))//'_sqrt'
end do
if (nam%check_dirac) call test_loc_dirac(nam,geom,bpar,ndata_other,ens1)
do iv=1,nam%nv
   nam%varname(iv) = varname(iv)
end do

! Reset lsqrt value
nam%lsqrt = .not.nam%lsqrt

! Print difference
write(mpl%unit,'(a7,a,f6.1,a)') '','Localization full / square-root error : ', &
 & sqrt(sum((fld_loc_sqrt-fld_loc)**2)/sum(fld_loc**2))*100.0,'%'
if (present(ens1)) write(mpl%unit,'(a7,a,f6.1,a)') '','Ensemble B full / square-root error:  ', &
 & sqrt(sum((fld_bens_sqrt-fld_bens)**2)/sum(fld_bens**2))*100.0,'%'

end subroutine test_loc_sqrt

!----------------------------------------------------------------------
! Subroutine: test_loc_dirac
!> Purpose: apply localization to diracs
!----------------------------------------------------------------------
subroutine test_loc_dirac(nam,geom,bpar,ndata,ens1)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                                                            !< Namelist
type(geomtype),intent(in) :: geom                                                          !< Geometry
type(bpartype),intent(in) :: bpar                                                          !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1)                                             !< NICAS data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: il0dir(nam%ndir),ic0dir(nam%ndir),idir,iv,its
real(kind_real),allocatable :: fld(:,:,:,:),fld_loc(:,:,:,:),fld_bens(:,:,:,:)
character(len=2) :: itschar

if (mpl%main) then
   ! Allocation
   allocate(fld(geom%nc0,geom%nl0,nam%nv,nam%nts))

   ! Find gridpoint and level indices
   call define_dirac(nam,geom,ic0dir,il0dir)

   ! Generate dirac field
   fld = 0.0
   do idir=1,nam%ndir
      fld(ic0dir(idir),il0dir(idir),nam%ivdir(idir),nam%itsdir(idir)) = 1.0
   end do
end if

! Global to local
call fld_com_gl(nam,geom,fld)

! Allocation
allocate(fld_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (present(ens1)) allocate(fld_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))

! Apply localization to dirac
fld_loc = fld
if (nam%lsqrt) then
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld_loc)
else
   call apply_localization(nam,geom,bpar,ndata,fld_loc)
end if

if (present(ens1)) then
   ! Apply localized ensemble covariance
   fld_bens = fld
   call apply_bens(nam,geom,bpar,ndata,ens1,fld_bens)
end if

! Local to global
call fld_com_lg(nam,geom,fld_loc)
if (present(ens1)) call fld_com_lg(nam,geom,fld_bens)

if (mpl%main) then
   ! Write field
   do its=1,nam%nts
      write(itschar,'(i2.2)') its
      do iv=1,nam%nv
         call model_write(nam,geom,trim(nam%prefix)//'_dirac_gridded.nc',trim(nam%varname(iv))//'_'//itschar,fld_loc(:,:,iv,its))
         call fld_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(nam%varname(iv))//'_'//itschar,fld_loc(:,:,iv,its))
         if (present(ens1)) then
            call model_write(nam,geom,trim(nam%prefix)//'_dirac_gridded.nc',trim(nam%varname(iv))//'_'//itschar//'_Bens', &
          & fld_bens(:,:,iv,its))
            call fld_write(nam,geom,trim(nam%prefix)//'_dirac.nc',trim(nam%varname(iv))//'_'//itschar//'_Bens', &
          & fld_bens(:,:,iv,its))
         end if
      end do
   end do
end if

end subroutine test_loc_dirac

!----------------------------------------------------------------------
! Subroutine: test_randomization
!> Purpose: test NICAS randomization method with respect to theoretical error statistics
!----------------------------------------------------------------------
subroutine test_randomization(nam,geom,bpar,ndata)

implicit none

! Passed variables
type(namtype),intent(inout),target :: nam      !< Namelist variables
type(geomtype),intent(in),target :: geom       !< Geometry
type(bpartype),intent(in) :: bpar              !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1) !< NICAS data

! Local variables
integer :: ifac,itest,nefac(nfac),ens1_ne,iv,its
real(kind_real) :: fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest),fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest)
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts),mse(ntest,nfac),mse_th(ntest,nfac)
real(kind_real),allocatable :: ens1(:,:,:,:,:),fld_tmp(:,:,:,:)
character(len=2) :: itschar
character(len=4) :: nechar,itestchar

! Define test vectors
write(mpl%unit,'(a4,a)') '','Define test vectors'
call define_test_vectors(nam,geom,fld_save)

! Apply localization to test vectors
write(mpl%unit,'(a4,a)') '','Apply localization to test vectors'
fld_ref = fld_save
do itest=1,ntest
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld_ref(:,:,:,:,itest))
end do

! Write first 10 test vectors
do itest=1,min(ntest,10)
   ! Allocation
   allocate(fld_tmp(geom%nc0a,geom%nl0,nam%nv,nam%nts))

   ! Copy
   fld_tmp = fld_ref(:,:,:,:,itest)

   ! Local to global
   call fld_com_lg(nam,geom,fld_tmp)

   if (mpl%main) then
      ! Write field
      write(itestchar,'(i4.4)') itest
      do its=1,nam%nts
         write(itschar,'(i2.2)') its
         do iv=1,nam%nv
            call model_write(nam,geom,trim(nam%prefix)//'_randomize_'//itestchar//'_gridded.nc', &
          & trim(nam%varname(iv))//'_ref_'//itschar,fld_tmp(:,:,iv,its))
         end do
      end do
   end if

   ! Release memory
   if (mpl%main) deallocate(fld_tmp)
end do

! Save namelist variables
ens1_ne = nam%ens1_ne

write(mpl%unit,'(a4,a)') '','Test randomization for various ensemble sizes:'
do ifac=1,nfac
   ! Ensemble size
   nefac(ifac) = max(int(5.0*float(ifac)/float(nfac)*float(ne)),3)
   nam%ens1_ne = nefac(ifac)
   write(nechar,'(i4.4)') nefac(ifac)

   ! Allocation
   allocate(ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nefac(ifac)))

   ! Randomize ensemble
   call randomize_localization(nam,geom,bpar,ndata,nefac(ifac),ens1)

   do itest=1,ntest
      ! Test localization
      fld = fld_save(:,:,:,:,itest)
      call apply_bens_noloc(nam,geom,bpar,ens1,fld)

      ! RMSE
      mse(itest,ifac) = sum((fld-fld_ref(:,:,:,:,itest))**2)
      mse_th(itest,ifac) = 1.0/float(nam%ens1_ne-1)*sum(1+fld_ref(:,:,:,:,itest)**2)

      ! Write first 10 test vectors
      if (itest<=min(ntest,10)) then
         ! Allocation
         allocate(fld_tmp(geom%nc0a,geom%nl0,nam%nv,nam%nts))

         ! Copy
         fld_tmp = fld

         ! Local to global
         call fld_com_lg(nam,geom,fld_tmp)

         if (mpl%main) then
            ! Write field
            write(itestchar,'(i4.4)') itest
            do its=1,nam%nts
               write(itschar,'(i2.2)') its
               do iv=1,nam%nv
                  call model_write(nam,geom,trim(nam%prefix)//'_randomize_'//itestchar//'_gridded.nc', &
                & trim(nam%varname(iv))//'_rand_'//nechar//'_'//itschar,fld_tmp(:,:,iv,its))
               end do
            end do
         end if

         ! Release memory
         if (mpl%main) deallocate(fld_tmp)
      end if
   end do

   ! Print scores
   write(mpl%unit,'(a7,a,i4,a,e15.8,a,e15.8)') '','Ensemble size ',nefac(ifac),', MSE (exp. / th.): ', &
 & sum(mse(:,ifac))/float(ntest),' / ',sum(mse_th(:,ifac))/float(ntest)

   ! Release memory
   deallocate(ens1)
end do

! Reset namelist variables
nam%ens1_ne = ens1_ne

end subroutine test_randomization

!----------------------------------------------------------------------
! Subroutine: test_consistency
!> Purpose: test HDIAG_NICAS consistency with a randomization method
!----------------------------------------------------------------------
subroutine test_consistency(nam,geom,bpar,bdata,ndata)

implicit none

! Passed variables
type(namtype),intent(inout) :: nam             !< Namelist variables
type(geomtype),intent(in) :: geom              !< Geometry
type(bpartype),intent(in) :: bpar              !< Block parameters
type(bdatatype),intent(in) :: bdata(bpar%nb+1) !< B data
type(ndatatype),intent(in) :: ndata(bpar%nb+1) !< NICAS data

! Local variables
integer :: ens1_ne,ens1_ne_offset,ens1_nsub,ib,il0
real(kind_real) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne)
character(len=1024) :: prefix,method
type(bdatatype),allocatable :: bdata_test(:)

! Randomize ensemble
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Randomize ensemble'
call randomize_localization(nam,geom,bpar,ndata,ne,ens1)

! Copy sampling
call system('cp -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc ' &
 & //trim(nam%datadir)//'/'//trim(nam%prefix)//'_consistency-test_sampling.nc')

! Save namelist variables
prefix = nam%prefix
method = nam%method
ens1_ne = nam%ens1_ne
ens1_ne_offset = nam%ens1_ne_offset
ens1_nsub = nam%ens1_nsub

! Set namelist variables
nam%prefix = trim(nam%prefix)//'_consistency-test'
nam%method = 'cor'
nam%ens1_ne = ne
nam%ens1_ne_offset = 0
nam%ens1_nsub = 1

! Call hdiag driver
call run_hdiag(nam,geom,bpar,bdata_test,ens1)

! Print scores
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- hdiag_nicas consistency results'
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      do il0=1,geom%nl0
         write(mpl%unit,'(a10,a7,i3,a4,a25,f6.1,a)') '','Level: ',nam%levs(il0),' ~> ','horizontal length-scale: ', &
       & sum(bdata_test(ib)%rh0(:,il0)-bdata(ib)%rh0(:,il0))/float(geom%nc0)*reqkm,' km'
         if (any(abs(bdata(ib)%rv0(:,il0))>0.0)) then
            write(mpl%unit,'(a49,f6.1,a)') 'vertical length-scale: ', &
          & sum(bdata_test(ib)%rv0(:,il0)-bdata(ib)%rv0(:,il0))/float(geom%nc0),' '//trim(vunitchar)
         end if
      end do
   end if
end do

! Reset namelist variables
nam%prefix = prefix
nam%method = method
nam%ens1_ne = ens1_ne
nam%ens1_ne_offset = ens1_ne_offset
nam%ens1_nsub = ens1_nsub

end subroutine test_consistency

!----------------------------------------------------------------------
! Subroutine: test_optimality
!> Purpose: test HDIAG localization optimality with a randomization method
!----------------------------------------------------------------------
subroutine test_optimality(nam,geom,bpar,ndata)

implicit none

! Passed variables
type(namtype),intent(inout),target :: nam      !< Namelist variables
type(geomtype),intent(in),target :: geom       !< Geometry
type(bpartype),intent(in) :: bpar              !< Block parameters
type(ndatatype),intent(in) :: ndata(bpar%nb+1) !< NICAS data

! Local variables
integer :: ib,ic0,ic0a,ifac,itest
real(kind_real) :: fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest),fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest)
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts),fac(nfac),mse(ntest,nfac)
real(kind_real) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne)
character(len=1024) :: prefix,method
type(bdatatype),allocatable :: bdata_save(:),bdata_test(:)
type(ndatatype),allocatable :: ndata_test(:)

! Define test vectors
write(mpl%unit,'(a4,a)') '','Define test vectors'
call define_test_vectors(nam,geom,fld_save)

! Apply localization to test vectors
write(mpl%unit,'(a4,a)') '','Apply localization to test vectors'
fld_ref = fld_save
do itest=1,ntest
   call apply_localization_from_sqrt(nam,geom,bpar,ndata,fld_ref(:,:,:,:,itest))
end do

! Randomize ensemble
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Randomize ensemble'
call randomize_localization(nam,geom,bpar,ndata,nam%ens1_ne,ens1)

! Copy sampling
call system('cp -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc ' &
 & //trim(nam%datadir)//'/'//trim(nam%prefix)//'_optimality-test_sampling.nc')

! Save namelist variables
prefix = nam%prefix
method = nam%method

! Set namelist variables
nam%prefix = trim(nam%prefix)//'_optimality-test'
nam%method = 'loc'

! Call hdiag driver
call run_hdiag(nam,geom,bpar,bdata_save,ens1)

! Allocation
call bdata_alloc(nam,geom,bpar,bdata_test)
allocate(ndata_test(bpar%nb+1))

! Copy bdata
call bdata_copy(bpar,bdata_save,bdata_test)

do ifac=1,nfac
   ! Multiplication factor
   fac(ifac) = 2.0*float(ifac)/float(nfac)

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,f4.2,a)') '--- Apply a multiplicative factor ',fac(ifac),' to length-scales'

   do ib=1,bpar%nb+1
      ! Set namelist and geometry
      ndata_test(ib)%nam => nam
      ndata_test(ib)%geom => geom

      if (bpar%nicas_block(ib)) then
         ! Length-scales multiplication
         bdata_test(ib)%rh0 = fac(ifac)*bdata_save(ib)%rh0
         bdata_test(ib)%rv0 = fac(ifac)*bdata_save(ib)%rv0
         if (trim(nam%strategy)=='specific_multivariate') then
            bdata_test(ib)%rh0s = fac(ifac)*bdata_save(ib)%rh0s
            bdata_test(ib)%rv0s = fac(ifac)*bdata_save(ib)%rv0s
         end if

         ! Compute NICAS parameters
         call compute_parameters(bdata_test(ib),ndata_test(ib))
      end if

      if (bpar%B_block(ib)) then
         ! Copy weights
         ndata_test(ib)%wgt = bdata_test(ib)%wgt
         if (bpar%nicas_block(ib)) then
            allocate(ndata_test(ib)%coef_ens(geom%nc0a,geom%nl0))
            do ic0=1,geom%nc0
               if (geom%c0_to_proc(ic0)==mpl%myproc) then
                  ic0a = geom%c0_to_c0a(ic0)
                  ndata_test(ib)%coef_ens(ic0a,:) = bdata_test(ib)%coef_ens(ic0,:)
               end if
            end do
         end if
      end if
   end do

   do itest=1,ntest
      ! Test localization
      fld = fld_save(:,:,:,:,itest)
      call apply_bens(nam,geom,bpar,ndata_test,ens1,fld)

      ! RMSE
      mse(itest,ifac) = sum((fld-fld_ref(:,:,:,:,itest))**2)
   end do


   ! Release memory
   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) call ndata_dealloc(ndata_test(ib))
   end do

   ! Print scores
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,f4.2,a,e15.8)') '--- Optimality results for a factor ',fac(ifac),', MSE: ',sum(mse(:,ifac))/float(ntest)
end do

! Print scores summary
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Optimality results summary'
do ifac=1,nfac
   write(mpl%unit,'(a7,a,f4.2,a,e15.8)') '','Factor ',fac(ifac),', MSE: ',sum(mse(:,ifac))/float(ntest)
end do

! Reset namelist variables
nam%prefix = prefix
nam%method = method

end subroutine test_optimality

!----------------------------------------------------------------------
! Subroutine: define_dirac
!> Purpose: define dirac indices
!----------------------------------------------------------------------
subroutine define_dirac(nam,geom,ic0dir,il0dir)

implicit none

! Passed variables
type(namtype),intent(in) :: nam         !< Namelist
type(geomtype),intent(in) :: geom       !< Geometry
integer,intent(out) :: ic0dir(nam%ndir) !< Dirac gridpoint indice
integer,intent(out) :: il0dir(nam%ndir) !< Dirac level indice

! Local variables
integer :: idir,il0,nn_index(1)
real(kind_real) :: nn_dist(1)

do idir=1,nam%ndir
   ! Find nearest neighbor
   call find_nearest_neighbors(geom%ctree,dble(nam%londir(idir)*deg2rad),dble(nam%latdir(idir)*deg2rad),1,nn_index,nn_dist)
   ic0dir(idir) = nn_index(1)

   ! Find level index
   do il0=1,geom%nl0
      if (nam%levs(il0)==nam%levdir(idir)) il0dir(idir) = il0
   end do
end do

end subroutine define_dirac

!----------------------------------------------------------------------
! Subroutine: define_test_vectors
!> Purpose: define test vectors
!----------------------------------------------------------------------
subroutine define_test_vectors(nam,geom,fld)

! Passed variables
type(namtype),intent(in) :: nam                                             !< Namelist
type(geomtype),intent(in) :: geom                                           !< Geometry
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest) !< Field

! Local variables
integer :: ic0dir(ntest),il0dir(ntest),ivdir(ntest),itsdir(ntest)
integer :: itest,ic0,iproc,ic0a

! Define dirac locations
if (mpl%main) then
   call rand_integer(1,geom%nc0,ic0dir)
   call rand_integer(1,geom%nl0,il0dir)
end if
call mpl_bcast(ic0dir,mpl%ioproc)
call mpl_bcast(il0dir,mpl%ioproc)
ivdir = 1
itsdir = 1

! Define test vectors
do itest=1,ntest
   fld(:,:,:,:,itest) = 0.0
   ic0 = ic0dir(itest)
   iproc = geom%c0_to_proc(ic0)
   if (iproc==mpl%myproc) then
      ic0a = geom%c0_to_c0a(ic0)
      fld(ic0a,il0dir(itest),ivdir(itest),itsdir(itest),itest) = 1.0
   end if
end do

end subroutine define_test_vectors

end module nicas_test
