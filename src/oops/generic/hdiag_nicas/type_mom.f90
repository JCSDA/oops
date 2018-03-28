!----------------------------------------------------------------------
! Module: type_mom
!> Purpose: moments derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_mom

use model_interface, only: model_read
use omp_lib
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr
use type_bpar, only: bpar_type
use type_com, only: com_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl
use type_hdata, only: hdata_type
use type_nam, only: nam_type

implicit none

! Moments derived type
type mom_type
   integer :: ne
   integer :: nsub
   type(mom_blk_type),allocatable :: blk(:) !< Moments blocks
contains
   procedure :: alloc => mom_alloc
   procedure :: compute => mom_compute
end type mom_type

private
public :: mom_type

contains

!----------------------------------------------------------------------
! Subroutine: mom_alloc
!> Purpose: allocate moments
!----------------------------------------------------------------------
subroutine mom_alloc(mom,nam,geom,bpar,hdata,ne,nsub)

implicit none

! Passed variables
class(mom_type),intent(inout) :: mom !< Moments
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data
integer,intent(in) :: ne             !< Ensemble size
integer,intent(in) :: nsub           !< Number of sub-ensembles

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
      allocate(mom%blk(ib)%m1_1(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      allocate(mom%blk(ib)%m2_1(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      allocate(mom%blk(ib)%m1_2(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      allocate(mom%blk(ib)%m2_2(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      allocate(mom%blk(ib)%m11(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      if (.not.nam%gau_approx) then
         allocate(mom%blk(ib)%m12(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
         allocate(mom%blk(ib)%m21(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
         allocate(mom%blk(ib)%m22(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,mom%blk(ib)%nsub))
      end if
      if (nam%full_var) then
         allocate(mom%blk(ib)%m1full(geom%nc0a,geom%nl0,mom%blk(ib)%nsub))
         allocate(mom%blk(ib)%m2full(geom%nc0a,geom%nl0,mom%blk(ib)%nsub))
      end if

      ! Initialization
      mom%blk(ib)%m1_1 = 0.0
      mom%blk(ib)%m2_1 = 0.0
      mom%blk(ib)%m1_2 = 0.0
      mom%blk(ib)%m2_2 = 0.0
      mom%blk(ib)%m11 = 0.0
      if (.not.nam%gau_approx) then
         mom%blk(ib)%m12 = 0.0
         mom%blk(ib)%m21 = 0.0
         mom%blk(ib)%m22 = 0.0
      end if
      if (nam%full_var) then
         mom%blk(ib)%m1full = 0.0
         mom%blk(ib)%m2full = 0.0
      end if
   end if
end do

end subroutine mom_alloc

!----------------------------------------------------------------------
! Subroutine: mom_compute
!> Purpose: compute centered moments (iterative formulae)
!----------------------------------------------------------------------
subroutine mom_compute(mom,nam,geom,bpar,hdata,filename,ens1)

implicit none

! Passed variables
class(mom_type),intent(inout) :: mom                                                       !< Moments
type(nam_type),intent(in) :: nam                                                           !< Namelist
type(geom_type),intent(in) :: geom                                                         !< Geometry
type(bpar_type),intent(in) :: bpar                                                         !< Block parameters
type(hdata_type),intent(in) :: hdata                                                       !< HDIAG data
character(len=*),intent(in) :: filename                                                    !< File name
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ne,ne_offset,nsub,ie,jc0,ic0c,jc0c,ic0,jl0r,jl0,il0,isub,jsub,jc3,ic1,ic1a,ib,jv,iv,jts,its
real(kind_real) :: fac1,fac2,fac3,fac4,fac5,fac6
real(kind_real),allocatable :: fld(:,:,:,:),fld_ext(:,:,:,:),fld_1(:,:,:,:),fld_2(:,:,:,:)
logical,allocatable :: mask_unpack(:,:)

! Setup
select case (trim(filename))
case ('ens1')
   ne = nam%ens1_ne
   ne_offset = nam%ens1_ne_offset
   nsub = nam%ens1_nsub
case ('ens2')
   ne = nam%ens2_ne
   ne_offset = nam%ens2_ne_offset
   nsub = nam%ens2_nsub
case default
   call msi(ne)
   call msi(ne_offset)
   call msi(nsub)
   call msgerror('wrong filename in ens_read')
end select

! Allocation
call mom%alloc(nam,geom,bpar,hdata,ne,nsub)
allocate(fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))

! Loop on sub-ensembles
do isub=1,nsub
   if (nsub==1) then
      write(mpl%unit,'(a10,a)',advance='no') '','Full ensemble, member:'
   else
      write(mpl%unit,'(a10,a,i4,a)',advance='no') '','Sub-ensemble ',isub,', member:'
   end if

   ! Compute centered moments iteratively
   do ie=1,ne/nsub
      write(mpl%unit,'(i4)',advance='no') ne_offset+ie

      ! Computation factors
      fac1 = 2.0/float(ie)
      fac2 = 1.0/float(ie**2)
      fac3 = float((ie-1)*(ie**2-3*ie+3))/float(ie**3)
      fac4 = 1.0/float(ie)
      fac5 = float((ie-1)*(ie-2))/float(ie**2)
      fac6 = float(ie-1)/float(ie)

      if (present(ens1)) then
         ! Copy field
         fld = ens1(:,:,:,:,ie+(isub-1)*nsub)
      else
         ! Load field
         if (nsub==1) then
            jsub = 0
         else
            jsub = isub
         end if
         call model_read(nam,geom,filename,ne_offset+ie,jsub,fld)
      end if

      ! Allocation
      allocate(fld_ext(hdata%nc0c,geom%nl0,nam%nv,nam%nts))
      allocate(mask_unpack(hdata%nc0c,geom%nl0))
      mask_unpack = .true.

      do ib=1,bpar%nb
         ! Indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         its = bpar%b_to_ts1(ib)
         jts = bpar%b_to_ts2(ib)

         ! Halo extension
         if ((iv==jv).and.(its==jts)) call hdata%com_AC%ext(geom%nl0,fld(:,:,iv,its),fld_ext(:,:,iv,its))
      end do

      do ib=1,bpar%nb
         if (bpar%diag_block(ib)) then
            ! Allocation
            allocate(fld_1(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))
            allocate(fld_2(hdata%nc1a,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))

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
                  call hdata%d(il0,its)%apply(fld_ext(:,il0,iv,its),fld_1(:,1,1,il0))
                  call hdata%d(il0,jts)%apply(fld_ext(:,il0,jv,jts),fld_2(:,1,1,il0))
               end do
               !$omp end parallel do
            else
               ! Copy all separations points
               !$omp parallel do schedule(static) private(il0,jl0r,jl0,jc3,ic1a,ic1,ic0,jc0,ic0c,jc0c)
               do il0=1,geom%nl0
                  do jl0r=1,bpar%nl0r(ib)
                     jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
                     do jc3=1,bpar%nc3(ib)
                        do ic1a=1,hdata%nc1a
                           ! Indices
                           ic1 = hdata%c1a_to_c1(ic1a)

                           if (hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,il0).and. &
                         & hdata%c1l0_log(ic1,jl0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) then
                              ! Indices
                              ic0 = hdata%c1_to_c0(ic1)
                              jc0 = hdata%c1c3_to_c0(ic1,jc3)
                              ic0c = hdata%c0_to_c0c(ic0)
                              jc0c = hdata%c0_to_c0c(jc0)

                              ! Copy points
                              fld_1(ic1a,jc3,jl0r,il0) = fld_ext(ic0c,il0,iv,its)
                              fld_2(ic1a,jc3,jl0r,il0) = fld_ext(jc0c,jl0,jv,jts)
                           end if
                        end do
                     end do
                  end do
               end do
               !$omp end parallel do
            end if

            !$omp parallel do schedule(static) private(il0)
            do il0=1,geom%nl0
               ! Remove means
               fld_1(:,:,:,il0) = fld_1(:,:,:,il0) - mom%blk(ib)%m1_1(:,:,:,il0,isub)
               fld_2(:,:,:,il0) = fld_2(:,:,:,il0) - mom%blk(ib)%m1_2(:,:,:,il0,isub)

               ! Update high-order moments
               if (ie>1) then
                  if (.not.nam%gau_approx) then
                       ! Fourth-order moment
                        mom%blk(ib)%m22(:,:,:,il0,isub) = mom%blk(ib)%m22(:,:,:,il0,isub) &
                                                    & -fac1*(mom%blk(ib)%m12(:,:,:,il0,isub)*fld_1(:,:,:,il0) &
                                                    & +mom%blk(ib)%m21(:,:,:,il0,isub)*fld_2(:,:,:,il0)) &
                                                    & +fac2*(4.0*mom%blk(ib)%m11(:,:,:,il0,isub)*fld_1(:,:,:,il0)*fld_2(:,:,:,il0) &
                                                    & +mom%blk(ib)%m2_1(:,:,:,il0,isub)*fld_2(:,:,:,il0)**2 &
                                                    & +mom%blk(ib)%m2_2(:,:,:,il0,isub)*fld_1(:,:,:,il0)**2) &
                                                    & +fac3*fld_1(:,:,:,il0)**2*fld_2(:,:,:,il0)**2

                        ! Third-order moments
                        mom%blk(ib)%m12(:,:,:,il0,isub) = mom%blk(ib)%m12(:,:,:,il0,isub) &
                                                    & -fac4*(2.0*mom%blk(ib)%m11(:,:,:,il0,isub)*fld_2(:,:,:,il0) &
                                                    & +mom%blk(ib)%m2_2(:,:,:,il0,isub)*fld_1(:,:,:,il0)) &
                                                    & +fac5*fld_2(:,:,:,il0)**2*fld_1(:,:,:,il0)

                        mom%blk(ib)%m21(:,:,:,il0,isub) = mom%blk(ib)%m21(:,:,:,il0,isub) &
                                                    & -fac4*(2.0*mom%blk(ib)%m11(:,:,:,il0,isub)*fld_1(:,:,:,il0) &
                                                    & +mom%blk(ib)%m2_1(:,:,:,il0,isub)*fld_2(:,:,:,il0)) &
                                                    & +fac5*fld_1(:,:,:,il0)**2*fld_2(:,:,:,il0)
                  end if

                  ! Covariance
                  mom%blk(ib)%m11(:,:,:,il0,isub) = mom%blk(ib)%m11(:,:,:,il0,isub)+fac6*fld_1(:,:,:,il0)*fld_2(:,:,:,il0)

                  ! Variances
                  mom%blk(ib)%m2_1(:,:,:,il0,isub) = mom%blk(ib)%m2_1(:,:,:,il0,isub)+fac6*fld_1(:,:,:,il0)**2
                  mom%blk(ib)%m2_2(:,:,:,il0,isub) = mom%blk(ib)%m2_2(:,:,:,il0,isub)+fac6*fld_2(:,:,:,il0)**2

                  ! Full variance
                  if (nam%full_var) mom%blk(ib)%m2full(:,il0,isub) = mom%blk(ib)%m2full(:,il0,isub) &
                                                               & +fac6*(fld(:,il0,iv,its)-mom%blk(ib)%m1full(:,il0,isub))**2
               end if

               ! Update means
               mom%blk(ib)%m1_1(:,:,:,il0,isub) = mom%blk(ib)%m1_1(:,:,:,il0,isub)+fac4*fld_1(:,:,:,il0)
               mom%blk(ib)%m1_2(:,:,:,il0,isub) = mom%blk(ib)%m1_2(:,:,:,il0,isub)+fac4*fld_2(:,:,:,il0)

               ! Full mean
               if (nam%full_var) mom%blk(ib)%m1full(:,il0,isub) = mom%blk(ib)%m1full(:,il0,isub) &
                                                            & +fac4*(fld(:,il0,iv,its)-mom%blk(ib)%m1full(:,il0,isub))
            end do
            !$omp end parallel do

            ! Release memory
            deallocate(fld_1)
            deallocate(fld_2)
         end if
      end do

      ! Release memory
      deallocate(fld_ext)
      deallocate(mask_unpack)
   end do
   write(mpl%unit,'(a)') ''

   ! Normalize statistics
   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         !$omp parallel do schedule(static) private(il0)
         do il0=1,geom%nl0
            mom%blk(ib)%m2_1(:,:,:,il0,isub) = mom%blk(ib)%m2_1(:,:,:,il0,isub)/float(ne/nsub-1)
            mom%blk(ib)%m2_2(:,:,:,il0,isub) = mom%blk(ib)%m2_2(:,:,:,il0,isub)/float(ne/nsub-1)
            mom%blk(ib)%m11(:,:,:,il0,isub) = mom%blk(ib)%m11(:,:,:,il0,isub)/float(ne/nsub-1)
            if (.not.nam%gau_approx) then
               mom%blk(ib)%m12(:,:,:,il0,isub) = mom%blk(ib)%m12(:,:,:,il0,isub)/float(ne/nsub)
               mom%blk(ib)%m21(:,:,:,il0,isub) = mom%blk(ib)%m21(:,:,:,il0,isub)/float(ne/nsub)
               mom%blk(ib)%m22(:,:,:,il0,isub) = mom%blk(ib)%m22(:,:,:,il0,isub)/float(ne/nsub)
            end if
            if (nam%full_var) mom%blk(ib)%m2full(:,il0,isub) = mom%blk(ib)%m2full(:,il0,isub)/float(ne/nsub-1)
         end do
         !$omp end parallel do
      end if
   end do
end do

end subroutine mom_compute

end module type_mom
