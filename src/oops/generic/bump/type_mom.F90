!----------------------------------------------------------------------
! Module: type_mom
! Purpose: moments derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mom

!$ use omp_lib
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr
use type_bpar, only: bpar_type
use type_com, only: com_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

! Moments derived type
type mom_type
   integer :: ne
   integer :: nsub
   type(mom_blk_type),allocatable :: blk(:) ! Moments blocks
contains
   procedure :: alloc => mom_alloc
   procedure :: compute => mom_compute
end type mom_type

private
public :: mom_type

contains

!----------------------------------------------------------------------
! Subroutine: mom_alloc
! Purpose: allocate moments
!----------------------------------------------------------------------
subroutine mom_alloc(mom,nam,geom,bpar,samp,ne,nsub)

implicit none

! Passed variables
class(mom_type),intent(inout) :: mom ! Moments
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(samp_type),intent(in) :: samp   ! Sampling
integer,intent(in) :: ne             ! Ensemble size
integer,intent(in) :: nsub           ! Number of sub-ensembles

! Local variables
integer :: ib

! Set attributes
mom%ne = ne
mom%nsub = nsub

! Allocation
allocate(mom%blk(bpar%nb))

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      ! Allocation
      mom%blk(ib)%ne = ne
      mom%blk(ib)%nsub = nsub
      allocate(mom%blk(ib)%m2_1(samp%nc1a,bpar%nc3(ib),geom%nl0,mom%blk(ib)%nsub))
      allocate(mom%blk(ib)%m2_2(samp%nc1a,bpar%nc3(ib),geom%nl0,mom%blk(ib)%nsub))
      allocate(mom%blk(ib)%m11(samp%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      if (.not.nam%gau_approx) allocate(mom%blk(ib)%m22(samp%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      if (nam%var_full) allocate(mom%blk(ib)%m2full(geom%nc0a,geom%nl0,mom%blk(ib)%nsub))

      ! Initialization
      mom%blk(ib)%m2_1 = 0.0
      mom%blk(ib)%m2_2 = 0.0
      mom%blk(ib)%m11 = 0.0
      if (.not.nam%gau_approx) mom%blk(ib)%m22 = 0.0
      if (nam%var_full) mom%blk(ib)%m2full = 0.0
   end if
end do

end subroutine mom_alloc

!----------------------------------------------------------------------
! Subroutine: mom_compute
! Purpose: compute centered moments (iterative formulae)
!----------------------------------------------------------------------
subroutine mom_compute(mom,mpl,nam,geom,bpar,samp,ens)

implicit none

! Passed variables
class(mom_type),intent(inout) :: mom ! Moments
type(mpl_type),intent(in) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(samp_type),intent(in) :: samp   ! Sampling
type(ens_type), intent(in) :: ens    ! Ensemble

! Local variables
integer :: ie,ie_sub,jc0,ic0c,jc0c,ic0,jl0r,jl0,il0,isub,jc3,ic1,ic1a,ib,jv,iv,jts,its
real(kind_real),allocatable :: fld_ext(:,:,:,:),fld_1(:,:,:),fld_2(:,:,:)
logical,allocatable :: mask_unpack(:,:)

! Allocation
call mom%alloc(nam,geom,bpar,samp,ens%ne,ens%nsub)

! Loop on sub-ensembles
do isub=1,ens%nsub
   if (ens%nsub==1) then
      write(mpl%info,'(a10,a)',advance='no') '','Full ensemble, member:'
   else
      write(mpl%info,'(a10,a,i4,a)',advance='no') '','Sub-ensemble ',isub,', member:'
   end if
   call flush(mpl%info)

   ! Compute centered moments
   do ie_sub=1,ens%ne/ens%nsub
      write(mpl%info,'(i4)',advance='no') ie_sub
      call flush(mpl%info)

      ! Full ensemble index
      ie = ie_sub+(isub-1)*ens%ne/ens%nsub

      ! Allocation
      allocate(fld_ext(samp%nc0c,geom%nl0,nam%nv,nam%nts))
      allocate(mask_unpack(samp%nc0c,geom%nl0))
      mask_unpack = .true.

      do ib=1,bpar%nb
         ! Indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         its = bpar%b_to_ts1(ib)
         jts = bpar%b_to_ts2(ib)

         ! Halo extension
         if ((iv==jv).and.(its==jts)) call samp%com_AC%ext(mpl,geom%nl0,ens%fld(:,:,iv,its,ie),fld_ext(:,:,iv,its))
      end do

      do ib=1,bpar%nb
         if (bpar%diag_block(ib)) then
            ! Allocation
            allocate(fld_1(samp%nc1a,bpar%nc3(ib),geom%nl0))
            allocate(fld_2(samp%nc1a,bpar%nc3(ib),geom%nl0))

            ! Initialization
            iv = bpar%b_to_v1(ib)
            jv = bpar%b_to_v2(ib)
            its = bpar%b_to_ts1(ib)
            jts = bpar%b_to_ts2(ib)

            ! Copy valid field points
            call msr(fld_1)
            call msr(fld_2)
            if ((iv/=jv).and.(its/=jts).and.nam%displ_diag) then
               ! Interpolate zero separation points
               !$omp parallel do schedule(static) private(il0)
               do il0=1,geom%nl0
                  call samp%d(il0,its)%apply(mpl,fld_ext(:,il0,iv,its),fld_1(:,1,il0))
                  call samp%d(il0,jts)%apply(mpl,fld_ext(:,il0,jv,jts),fld_2(:,1,il0))
               end do
               !$omp end parallel do
            else
               ! Copy all separations points
               !$omp parallel do schedule(static) private(il0,jc3,ic1a,ic1,ic0,jc0,ic0c,jc0c)
               do il0=1,geom%nl0
                  do jc3=1,bpar%nc3(ib)
                     do ic1a=1,samp%nc1a
                        ! Indices
                        ic1 = samp%c1a_to_c1(ic1a)

                        if (samp%c1l0_log(ic1,il0).and.samp%c1c3l0_log(ic1,jc3,il0)) then
                           ! Indices
                           ic0 = samp%c1_to_c0(ic1)
                           jc0 = samp%c1c3_to_c0(ic1,jc3)
                           ic0c = samp%c0_to_c0c(ic0)
                           jc0c = samp%c0_to_c0c(jc0)

                           ! Copy points
                           fld_1(ic1a,jc3,il0) = fld_ext(ic0c,il0,iv,its)
                           fld_2(ic1a,jc3,il0) = fld_ext(jc0c,il0,jv,jts)
                        end if
                     end do
                  end do
               end do
               !$omp end parallel do
            end if

            !$omp parallel do schedule(static) private(il0,jl0r,jl0)
            do il0=1,geom%nl0
               do jl0r=1,bpar%nl0r(ib)
                  jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)

                  ! Fourth-order moment
                  if (.not.nam%gau_approx) mom%blk(ib)%m22(:,:,jl0r,il0,isub) = mom%blk(ib)%m22(:,:,jl0r,il0,isub) &
                                                                              & +fld_1(:,:,il0)**2*fld_2(:,:,jl0)**2

                  ! Covariance
                  mom%blk(ib)%m11(:,:,jl0r,il0,isub) = mom%blk(ib)%m11(:,:,jl0r,il0,isub)+fld_1(:,:,il0)*fld_2(:,:,jl0)
               end do
            end do
            !$omp end parallel do

            ! Variances
            mom%blk(ib)%m2_1(:,:,:,isub) = mom%blk(ib)%m2_1(:,:,:,isub)+fld_1**2
            mom%blk(ib)%m2_2(:,:,:,isub) = mom%blk(ib)%m2_2(:,:,:,isub)+fld_2**2

            ! Full variance
            if (nam%var_full) mom%blk(ib)%m2full(:,:,isub) = mom%blk(ib)%m2full(:,:,isub)+ens%fld(:,:,iv,its,ie)**2

            ! Release memory
            deallocate(fld_1)
            deallocate(fld_2)
         end if
      end do

      ! Release memory
      deallocate(fld_ext)
      deallocate(mask_unpack)
   end do
   write(mpl%info,'(a)') ''
   call flush(mpl%info)
end do

! Normalize moments
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      mom%blk(ib)%m2_1 = mom%blk(ib)%m2_1/real(mom%ne/mom%nsub-1,kind_real)
      mom%blk(ib)%m2_2 = mom%blk(ib)%m2_2/real(mom%ne/mom%nsub-1,kind_real)
      mom%blk(ib)%m11 = mom%blk(ib)%m11/real(mom%ne/mom%nsub-1,kind_real)
      if (.not.nam%gau_approx) mom%blk(ib)%m22 = mom%blk(ib)%m22/real(mom%ne/mom%nsub,kind_real)
      if (nam%var_full) mom%blk(ib)%m2full = mom%blk(ib)%m2full/real(mom%ne/mom%nsub-1,kind_real)
   end if
end do

end subroutine mom_compute

end module type_mom
