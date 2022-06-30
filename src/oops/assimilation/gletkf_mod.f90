module letkf
!$$$  module documentation block
!
! module: letkf                        Update model state variables and
!                                      bias coefficients with the LETKF.
!
! prgmmr: ota              org: np23                   date: 2011-06-01
!         updates, optimizations by whitaker
!         adaptation to oops: frolov
!
! abstract: Updates the model state using the LETKF (Hunt et al 2007,
!  Physica D, 112-126). Uses 'gain form' of LETKF algorithm described
!  Bishop et al 2017 (https://doi.org/10.1175/MWR-D-17-0102.1).
!
!   2020-04-01 (frolov) adopted letkf_core routine from GSI to be used 
!               inside of the JEDI GETKFSolver.h
!
! attributes:
!   language: f95
!
!$$$

use iso_c_binding

implicit none

integer, parameter :: r_double=c_double, r_single=c_float, num_bytes_for_r_single=4
integer, parameter :: i_kind=c_int, r_kind=r_single

contains

subroutine letkf_core(nobsl,hxens,hxens_orig,dep,wts_ensmean,wts_ensperts,&
                      rdiaginv_loc,nanals,neigv,getkf_inflation,denkf,getkf,mult_infl)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    letkf_core
!
!   prgmmr: whitaker
!
! abstract:  LETKF core subroutine. Returns analysis weights
! for ensemble mean update and ensemble perturbation update (increments
! represented as a linear combination of prior ensemble perturbations).
! Uses 'gain form' LETKF, which works when ensemble used to estimate
! covariances not the same as ensemble being updated.  If getkf=false
! (no modulated ensemble model-space localization), then the 
! traditional form of the LETKF analysis weights for ensemble
! perturbations are returned (which represents posterior
! ensemble perturbations, not analysis increments, as a linear
! combination of prior ensemble perturbations).
!
! program history log:
!   2011-06-03  ota: created from miyoshi's LETKF core subroutine
!   2014-06-20  whitaker: Use LAPACK routine dsyev for eigenanalysis.
!   2016-02-01  whitaker: Use LAPACK dsyevr for eigenanalysis (faster
!               than dsyev in most cases). 
!   2018-07-01  whitaker: implement gain form of LETKF from Bishop et al 2017
!   (https://doi.org/10.1175/MWR-D-17-0102.1), allow for use of modulated 
!   ensemble vert localization (ensemble used to estimate posterior covariance
!   in ensemble space different than ensemble being updated). Add denkf,
!   getkf,getkf_inflation options.
!
!   input argument list:
!     nobsl    - number of observations in the local patch
!     hxens    - on input: first-guess modulated ensemble in observation space (Yb)
!                on output: overwritten with Yb * R**-1.
!     hxens_orig  - first-guess original ensembles in observation space.
!                   not used if getkf=false. 
!     dep      - nobsl observation departures (y-Hxmean)
!     rdiaginv_loc - inverse of diagonal element of observation error covariance
!                localized (inflated) based on distance to analysis point
!     nanals   - number of ensemble members (1st dimension of hxens)
!     neigv    - for modulated ensemble model-space localization, number
!                of eigenvectors of vertical localization (1 if not using
!                model space localization).  1st dimension of hxens_orig is 
!                nanals/neigv.
!     getkf_inflation - if true (and getkf=T,denkf=F), 
!                return posterior covariance matrix in
!                needed to compute getkf inflation (eqn 30 in Bishop et al
!                2017).
!     denkf - if true, use DEnKF approximation (implies getkf=T)
!             See Sakov and Oke 2008 https://doi.org/10.1111/j.1600-0870.2007.00299.x
!     getkf - if true, use gain formulation
!     mult_infl - multiplicative prior inflation.
!                 mult_infl == 1 is no inflation. mult_infl>1 prior is inflated. 
!
!   output argument list:
!
!     wts_ensmean - Factor used to compute ens mean analysis increment
!       by pre-multiplying with
!       model space ensemble perts. In notation from Bishop et al 2017,
!       wts_ensmean = C (Gamma + I)**-1 C^T (HZ)^ T R**-1/2 (y - Hxmean)
!       where HZ^T = Yb*R**-1/2 (YbRinvsqrt),
!       C are eigenvectors of (HZ)^T HZ and Gamma are eigenvalues
!       Has dimension (nanals) - increment is weighted average of ens
!       perts, wts_ensmean are weights. 
!
!     wts_ensperts  - if getkf=T same as above, but for computing increments to 
!       ensemble perturbations. From Bishop et al 2017 eqn 29
!       wts_ensperts = -C [ (I - (Gamma+I)**-1/2)*Gamma**-1 ] C^T (HZ)^T R**-1/2 Hxprime
!       Has dimension (nanals,nanals/neigv), analysis weights for each
!       member. Hxprime (hxens_orig) is the original, unmodulated
!       ensemble in observation space, HZ is the modulated ensemble in
!       ob space times R**-1/2. If denkf=T, wts_ensperts is approximated
!       as wts_ensperts = -0.5*C (Gamma + I)**-1 C^T (HZ)^ T R**-1/2 Hxprime
!       If getkf=F and denkf=F, then the original LETKF formulation is used with
!       wts_ensperts =
!       C (Gamma + I)**-1/2 C^T (square root of analysis error cov in ensemble space)
!       and these weights are applied to transform the background ensemble into an
!       analysis ensemble.  Note that modulated ensemble vertical localization
!       requires the gain form (getkf=T and/or denkf=T) since this form of the weights
!       requires that the background ensemble used to compute covariances is
!       the same ensemble being updated.
!    
!     paens - only allocated and returned
!       if getkf_inflation=T (and denkf=F).  In this case
!       paens is allocated dimension (nanals,nanals) and contains posterior 
!       covariance matrix in (modulated) ensemble space.
!
! attributes:
!   language:  f95
!   machine:
!
!$$$ end documentation block

implicit none
integer(i_kind), intent(in) :: nobsl,nanals,neigv
real(r_kind),dimension(nobsl),intent(in ) :: rdiaginv_loc
real(r_kind),dimension(nanals,nobsl),intent(inout)  :: hxens
real(r_single),dimension(nanals/neigv,nobsl),intent(in)  :: hxens_orig
real(r_single),dimension(nobsl),intent(in)  :: dep
real(r_single),dimension(nanals),intent(out)  :: wts_ensmean
real(r_single),dimension(nanals,nanals/neigv),intent(out)  :: wts_ensperts
real(r_single), intent(in) :: mult_infl
!real(r_single),dimension(:,:),allocatable, intent(inout) :: paens
! local variables.
real(r_kind),allocatable,dimension(:,:) :: work3,evecs
real(r_single),allocatable,dimension(:,:) :: swork2,pa,swork3,shxens
real(r_single),allocatable,dimension(:) :: swork1
real(r_kind),allocatable,dimension(:) :: rrloc,evals,gammapI,gamma_inv
real(r_kind) eps
integer(i_kind) :: nanal,ierr,lwork,liwork
!for LAPACK dsyevr
integer(i_kind) isuppz(2*nanals)
real(r_kind) vl,vu,normfact
integer(i_kind), allocatable, dimension(:) :: iwork
real(r_kind), dimension(:), allocatable :: work1
logical, intent(in) :: getkf_inflation,denkf,getkf

if (neigv < 1) then
  !print *,'neigv must be >=1 in letkf_core'
  call abor1_ftn("neigv must be >=1 in letkf_core")
endif

allocate(work3(nanals,nanals),evecs(nanals,nanals))
allocate(rrloc(nobsl),gammapI(nanals),evals(nanals),gamma_inv(nanals))
! for dsyevr
!allocate(iwork(10*nanals),work1(70*nanals))
! for dsyevd
allocate(iwork(3+5*nanals),work1(1+6*nanals+2*nanals*nanals))

! HZ^T = hxens sqrt(Rinv)
rrloc = rdiaginv_loc
eps = epsilon(0.0_r_single)
where (rrloc < eps) rrloc = eps
rrloc = sqrt(rrloc)
normfact = sqrt(real((nanals/neigv)-1,r_kind))
! normalize so dot product is covariance

!$OMP PARALLEL DO PRIVATE(nanal) shared(hxens, nobsl, rrloc, normfact)
do nanal=1,nanals
   hxens(nanal,1:nobsl) = hxens(nanal,1:nobsl) * &
   rrloc(1:nobsl)/normfact
end do
!$OMP END PARALLEL DO

! compute eigenvectors/eigenvalues of HZ^T HZ (left SV)
! (in Bishop paper HZ is nobsl, nanals, here is it nanals, nobsl)
lwork = size(work1); liwork = size(iwork)
if(r_kind == kind(1.d0)) then ! double precision 
   !work3 = matmul(hxens,transpose(hxens))
   call dgemm('n','t',nanals,nanals,nobsl,1.d0,hxens,nanals, &
               hxens,nanals,0.d0,work3,nanals)
   ! evecs contains eigenvectors of HZ^T HZ, or left singular vectors of HZ
   ! evals contains eigenvalues (singular values squared)
!   call dsyevr('V','A','L',nanals,work3,nanals,vl,vu,1,nanals,-1.d0,nanals,evals,evecs, &
!               nanals,isuppz,work1,lwork,iwork,liwork,ierr)
! use LAPACK dsyevd instead of dsyevr
   evecs = work3
   call dsyevd('V','L',nanals,evecs,nanals,evals,work1,lwork,iwork,liwork,ierr)
else ! single precision
   call sgemm('n','t',nanals,nanals,nobsl,1.e0,hxens,nanals, &
               hxens,nanals,0.e0,work3,nanals)
!   call ssyevr('V','A','L',nanals,work3,nanals,vl,vu,1,nanals,-1.e0,nanals,evals,evecs, &
!               nanals,isuppz,work1,lwork,iwork,liwork,ierr)

! use LAPACK dsyevd instead of dsyevr
   evecs = work3
   call ssyevd('V','L',nanals,evecs,nanals,evals,work1,lwork,iwork,liwork,ierr)
end if
if (ierr .ne. 0) print *,'warning: dsyev* failed, ierr=',ierr
deallocate(work1,iwork,work3) ! no longer needed
gamma_inv = 0.0_r_kind

!$OMP PARALLEL DO PRIVATE(nanal)
do nanal=1,nanals
   if (evals(nanal) > eps) then
       gamma_inv(nanal) = 1./evals(nanal)
   else
       evals(nanal) = 0.0_r_kind
   endif
enddo
!$OMP END PARALLEL DO

! gammapI used in calculation of posterior cov in ensemble space
gammapI = evals+1.0/mult_infl
deallocate(evals)

! create HZ^T R**-1/2 
allocate(shxens(nanals,nobsl))
!$OMP PARALLEL DO PRIVATE(nanal)
do nanal=1,nanals
   shxens(nanal,1:nobsl) = hxens(nanal,1:nobsl) * rrloc(1:nobsl)
end do
!$OMP END PARALLEL DO
deallocate(rrloc)

! compute factor to multiply with model space ensemble perturbations
! to compute analysis increment (for mean update), save in single precision.
! This is the factor C (Gamma + I)**-1 C^T (HZ)^ T R**-1/2 (y - HXmean)
! in Bishop paper (eqs 10-12).

allocate(swork3(nanals,nanals),swork2(nanals,nanals),pa(nanals,nanals))
!$OMP PARALLEL DO PRIVATE(nanal)
do nanal=1,nanals
   swork3(nanal,:) = evecs(nanal,:)/gammapI
   swork2(nanal,:) = evecs(nanal,:)
enddo
!$OMP END PARALLEL DO

! pa = C (Gamma + I)**-1 C^T (analysis error cov in ensemble space)
!pa = matmul(swork3,transpose(swork2))
call sgemm('n','t',nanals,nanals,nanals,1.e0,swork3,nanals,swork2,&
            nanals,0.e0,pa,nanals)
! work1 = (HZ)^ T R**-1/2 (y - HXmean)
! (nanals, nobsl) x (nobsl,) = (nanals,)
! in Bishop paper HZ is nobsl, nanals, here is it nanals, nobsl
allocate(swork1(nanals))
!$OMP PARALLEL DO PRIVATE(nanal)
do nanal=1,nanals
   swork1(nanal) = sum(shxens(nanal,:)*dep(:))
end do
!$OMP END PARALLEL DO
! wts_ensmean = C (Gamma + I)**-1 C^T (HZ)^ T R**-1/2 (y - HXmean)
! (nanals, nanals) x (nanals,) = (nanals,)
!$OMP PARALLEL DO PRIVATE(nanal)
do nanal=1,nanals
   wts_ensmean(nanal) = sum(pa(nanal,:)*swork1(:))/normfact
end do
!$OMP END PARALLEL DO

!if (.not. denkf .and. getkf_inflation) then
!   allocate(paens(nanals,nanals))
!   paens = pa/normfact**2
!endif
deallocate(swork1)

! compute factor to multiply with model space ensemble perturbations
! to compute analysis increment (for perturbation update), save in single precision.
! This is -C [ (I - (Gamma+I)**-1/2)*Gamma**-1 ] C^T (HZ)^T R**-1/2 HXprime
! in Bishop paper (eqn 29).
! For DEnKF factor is -0.5*C (Gamma + I)**-1 C^T (HZ)^ T R**-1/2 HXprime
! = -0.5 Pa (HZ)^ T R**-1/2 HXprime (Pa already computed)

if (getkf) then ! use Gain formulation for LETKF weights

if (denkf) then
   ! use Pa = C (Gamma + I)**-1 C^T (already computed)
   ! wts_ensperts = -0.5 Pa (HZ)^ T R**-1/2 HXprime
   pa = 0.5*pa
else
   gammapI = sqrt(1.0/gammapI)
!$OMP PARALLEL DO PRIVATE(nanal)
   do nanal=1,nanals
      swork3(nanal,:) = &
      evecs(nanal,:)*(1.-gammapI(:))*gamma_inv(:)
   enddo
!$OMP END PARALLEL DO
   ! swork2 still contains eigenvectors, over-write pa
   ! pa = C [ (I - (Gamma+I)**-1/2)*Gamma**-1 ] C^T
   !pa = matmul(swork3,transpose(swork2))
   call sgemm('n','t',nanals,nanals,nanals,1.e0,swork3,nanals,swork2,&
               nanals,0.e0,pa,nanals)
endif
deallocate(swork2,swork3)

! work2 = (HZ)^ T R**-1/2 HXprime
! (nanals, nobsl) x (nobsl, nanals/neigv) = (nanals, nanals/neigv)
! in Bishop paper HZ is nobsl, nanals, here is it nanals, nobsl
! HXprime in paper is nobsl, nanals/neigv here it is nanals/neigv, nobsl
allocate(swork2(nanals,nanals/neigv))
!swork2 = matmul(shxens,transpose(hxens_orig))
call sgemm('n','t',nanals,nanals/neigv,nobsl,1.e0,&
            shxens,nanals,hxens_orig,nanals/neigv,0.e0,swork2,nanals)
! wts_ensperts = -C [ (I - (Gamma+I)**-1/2)*Gamma**-1 ] C^T (HZ)^T R**-1/2 HXprime
! (nanals, nanals) x (nanals, nanals/eigv) = (nanals, nanals/neigv)
! if denkf, wts_ensperts = -0.5 C (Gamma + I)**-1 C^T (HZ)^T R**-1/2 HXprime
!wts_ensperts = -matmul(pa, swork2)/normfact
call sgemm('n','n',nanals,nanals/neigv,nanals,-1.e0,&
            pa,nanals,swork2,nanals,0.e0,wts_ensperts,nanals)
wts_ensperts = wts_ensperts/normfact

! clean up
deallocate(shxens,swork2,pa)

else  ! use original LETKF formulation (won't work if neigv != 1)

if (neigv > 1) then
  !print *,'neigv must be 1 in letkf_core if getkf=F'
  !stop 993
  call abor1_ftn("neigv must be 1 in letkf_core if getkf=F")
endif
! compute sqrt(Pa) - analysis weights
! (apply to prior ensemble to determine posterior ensemble,
!  not analysis increments as in Gain formulation)
! hxens_orig not used
! saves two matrix multiplications (nanals, nobsl) x (nobsl, nanals) and
! (nanals, nanals) x (nanals, nanals)
deallocate(shxens,pa)
gammapI = sqrt(1.0/gammapI)
!$OMP PARALLEL DO PRIVATE(nanal)
do nanal=1,nanals
   swork3(nanal,:) = evecs(nanal,:)*gammapI
enddo
!$OMP END PARALLEL DO
! swork2 already contains evecs
! wts_ensperts = 
! C (Gamma + I)**-1/2 C^T (square root of analysis error cov in ensemble space)
!wts_ensperts = matmul(swork3,transpose(swork2))
call sgemm('n','t',nanals,nanals,nanals,1.0,swork3,nanals,swork2,&
            nanals,0.e0,wts_ensperts,nanals)
deallocate(swork3,swork2)

endif

deallocate(evecs,gammapI,gamma_inv)


return
end subroutine letkf_core

end module letkf
