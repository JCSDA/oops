!----------------------------------------------------------------------
! Module: module_moments.f90
!> Purpose: moments routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_moments

use model_interface, only: model_read
use omp_lib
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr
use type_displ, only: displtype
use type_geom, only: fld_com_lg
use type_linop, only: apply_linop
use type_mom, only: momtype
use type_mpl, only: mpl,mpl_bcast
use type_hdata, only: hdatatype
implicit none

private
public :: compute_moments

contains

!----------------------------------------------------------------------
! Subroutine: compute_moments
!> Purpose: compute centered moments (iterative formulae)
!----------------------------------------------------------------------
subroutine compute_moments(hdata,filename,displ,mom,ens1)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                                                                                      !< HDIAG data
character(len=*),intent(in) :: filename                                                                                  !< File name
type(displtype),intent(in) :: displ                                                                                      !< Displacement
type(momtype),intent(inout) :: mom(hdata%bpar%nb)                                                                        !< Moments
real(kind_real),intent(in),optional :: ens1(hdata%geom%nc0a,hdata%geom%nl0,hdata%nam%nv,hdata%nam%nts,hdata%nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ne,ne_offset,nsub,ie,ic0,jc0,il0r,il0,jl0,isub,jsub,ic,ic1,ib,iv,jv,its,jts
real(kind_real) :: fac1,fac2,fac3,fac4,fac5,fac6
real(kind_real),allocatable :: fld(:,:,:,:),fld_1(:,:,:,:),fld_2(:,:,:,:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

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

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      ! Allocation
      mom(ib)%ne = ne
      mom(ib)%nsub = nsub
      allocate(mom(ib)%m1_1(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
      allocate(mom(ib)%m2_1(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
      allocate(mom(ib)%m1_2(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
      allocate(mom(ib)%m2_2(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
      allocate(mom(ib)%m11(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
      if (.not.nam%gau_approx) then
         allocate(mom(ib)%m12(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
         allocate(mom(ib)%m21(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
         allocate(mom(ib)%m22(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0,mom(ib)%nsub))
      end if
      if (nam%full_var) then
         allocate(mom(ib)%m1full(geom%nc0,geom%nl0,mom(ib)%nsub))
         allocate(mom(ib)%m2full(geom%nc0,geom%nl0,mom(ib)%nsub))
      end if

      ! Initialization
      mom(ib)%m1_1 = 0.0
      mom(ib)%m2_1 = 0.0
      mom(ib)%m1_2 = 0.0
      mom(ib)%m2_2 = 0.0
      mom(ib)%m11 = 0.0
      if (.not.nam%gau_approx) then
         mom(ib)%m12 = 0.0
         mom(ib)%m21 = 0.0
         mom(ib)%m22 = 0.0
      end if
      if (nam%full_var) then
         mom(ib)%m1full = 0.0
         mom(ib)%m2full = 0.0
      end if
   end if
end do

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
         ! Copy and broadcast field
         if (allocated(fld)) deallocate(fld)
         allocate(fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))
         fld = ens1(:,:,:,:,ie+(isub-1)*nsub)
         call fld_com_lg(nam,geom,fld)
         if (.not.mpl%main) allocate(fld(geom%nc0,geom%nl0,nam%nv,nam%nts))
         call mpl_bcast(fld,mpl%ioproc)
      else
         ! Load field
         if (allocated(fld)) deallocate(fld)
         allocate(fld(geom%nc0,geom%nl0,nam%nv,nam%nts))
         if (nsub==1) then
            jsub = 0
         else
            jsub = isub
         end if
         call model_read(nam,geom,filename,ie,jsub,fld)
      end if

      do ib=1,bpar%nb
         if (bpar%diag_block(ib)) then
            ! Allocation
            allocate(fld_1(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0))
            allocate(fld_2(nam%nc1,bpar%icmax(ib),bpar%nl0(ib),geom%nl0))

            ! Initialization
            iv = bpar%ib_to_iv(ib)
            jv = bpar%ib_to_jv(ib)
            its = bpar%ib_to_its(ib)
            jts = bpar%ib_to_jts(ib)

            ! Copy valid field points
            call msr(fld_1)
            call msr(fld_2)
            if ((iv==jv).and.(its==jts)) then
               ! Copy all separations points
               !$omp parallel do schedule(static) private(jl0,il0r,il0,ic,ic1,ic0,jc0)
               do jl0=1,geom%nl0
                  do il0r=1,bpar%nl0(ib)
                     il0 = bpar%il0rjl0ib_to_il0(il0r,jl0,ib)
                     do ic=1,bpar%icmax(ib)
                        do ic1=1,nam%nc1
                           if (hdata%ic1il0_log(ic1,jl0).and.hdata%ic1icil0_log(ic1,ic,il0)) then
                              ! Indices
                              ic0 = hdata%ic1icil0_to_ic0(ic1,ic,il0)
                              jc0 = hdata%ic1_to_ic0(ic1)

                              ! Copy points
                              fld_1(ic1,ic,il0r,jl0) = fld(jc0,jl0,jv,jts)
                              fld_2(ic1,ic,il0r,jl0) = fld(ic0,il0,iv,its)
                           end if
                        end do
                     end do
                  end do
               end do
               !$omp end parallel do
            else
               if (nam%displ_diag) then
                  ! Interpolate zero separation points
                  do jl0=1,geom%nl0
                     call apply_linop(displ%d(jl0,jts),fld(:,jl0,jv,jts),fld_1(:,1,1,jl0))
                     call apply_linop(displ%d(jl0,its),fld(:,jl0,iv,its),fld_2(:,1,1,jl0))
                  end do
               else
                  ! Copy zero separation points
                  !$omp parallel do schedule(static) private(jl0,ic1,jc0)
                  do jl0=1,geom%nl0
                     do ic1=1,nam%nc1
                        if (hdata%ic1il0_log(ic1,jl0)) then
                           ! Indices
                           jc0 = hdata%ic1_to_ic0(ic1)

                           ! Copy points
                           fld_1(ic1,1,1,jl0) = fld(jc0,jl0,jv,jts)
                           fld_2(ic1,1,1,jl0) = fld(jc0,jl0,iv,its)
                        end if
                     end do
                  end do
                  !$omp end parallel do
               end if
            end if

            !$omp parallel do schedule(static) private(jl0)
            do jl0=1,geom%nl0
               ! Remove means
               fld_1(:,:,:,jl0) = fld_1(:,:,:,jl0) - mom(ib)%m1_1(:,:,:,jl0,isub)
               fld_2(:,:,:,jl0) = fld_2(:,:,:,jl0) - mom(ib)%m1_2(:,:,:,jl0,isub)

               ! Update high-order moments
               if (ie>1) then
                  if (.not.nam%gau_approx) then
                       ! Fourth-order moment
                        mom(ib)%m22(:,:,:,jl0,isub) = mom(ib)%m22(:,:,:,jl0,isub) &
                                                  & -fac1*(mom(ib)%m12(:,:,:,jl0,isub)*fld_1(:,:,:,jl0) &
                                                  & +mom(ib)%m21(:,:,:,jl0,isub)*fld_2(:,:,:,jl0)) &
                                                  & +fac2*(4.0*mom(ib)%m11(:,:,:,jl0,isub)*fld_1(:,:,:,jl0)*fld_2(:,:,:,jl0) &
                                                  & +mom(ib)%m2_1(:,:,:,jl0,isub)*fld_2(:,:,:,jl0)**2 &
                                                  & +mom(ib)%m2_2(:,:,:,jl0,isub)*fld_1(:,:,:,jl0)**2) &
                                                  & +fac3*fld_1(:,:,:,jl0)**2*fld_2(:,:,:,jl0)**2

                        ! Third-order moments
                        mom(ib)%m12(:,:,:,jl0,isub) = mom(ib)%m12(:,:,:,jl0,isub) &
                                                  & -fac4*(2.0*mom(ib)%m11(:,:,:,jl0,isub)*fld_2(:,:,:,jl0) &
                                                  & +mom(ib)%m2_2(:,:,:,jl0,isub)*fld_1(:,:,:,jl0)) &
                                                  & +fac5*fld_2(:,:,:,jl0)**2*fld_1(:,:,:,jl0)

                        mom(ib)%m21(:,:,:,jl0,isub) = mom(ib)%m21(:,:,:,jl0,isub) &
                                                  & -fac4*(2.0*mom(ib)%m11(:,:,:,jl0,isub)*fld_1(:,:,:,jl0) &
                                                  & +mom(ib)%m2_1(:,:,:,jl0,isub)*fld_2(:,:,:,jl0)) &
                                                  & +fac5*fld_1(:,:,:,jl0)**2*fld_2(:,:,:,jl0)
                  end if

                  ! Covariance
                  mom(ib)%m11(:,:,:,jl0,isub) = mom(ib)%m11(:,:,:,jl0,isub)+fac6*fld_1(:,:,:,jl0)*fld_2(:,:,:,jl0)

                  ! Variances
                  mom(ib)%m2_1(:,:,:,jl0,isub) = mom(ib)%m2_1(:,:,:,jl0,isub)+fac6*fld_1(:,:,:,jl0)**2
                  mom(ib)%m2_2(:,:,:,jl0,isub) = mom(ib)%m2_2(:,:,:,jl0,isub)+fac6*fld_2(:,:,:,jl0)**2

                  ! Full variance
                  if (nam%full_var) mom(ib)%m2full(:,jl0,isub) = mom(ib)%m2full(:,jl0,isub) &
                                                             & +fac6*(fld(:,jl0,jv,jts)-mom(ib)%m1full(:,jl0,isub))**2
               end if

               ! Update means
               mom(ib)%m1_1(:,:,:,jl0,isub) = mom(ib)%m1_1(:,:,:,jl0,isub)+fac4*fld_1(:,:,:,jl0)
               mom(ib)%m1_2(:,:,:,jl0,isub) = mom(ib)%m1_2(:,:,:,jl0,isub)+fac4*fld_2(:,:,:,jl0)

               ! Full mean
               if (nam%full_var) mom(ib)%m1full(:,jl0,isub) = mom(ib)%m1full(:,jl0,isub) &
                                                            & +fac4*(fld(:,jl0,jv,jts)-mom(ib)%m1full(:,jl0,isub))
            end do
            !$omp end parallel do

            ! Release memory
            deallocate(fld_1)
            deallocate(fld_2)
         end if
      end do
   end do
   write(mpl%unit,'(a)') ''

   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         !$omp parallel do schedule(static) private(jl0)
         do jl0=1,geom%nl0
            ! Normalize
            mom(ib)%m2_1(:,:,:,jl0,isub) = mom(ib)%m2_1(:,:,:,jl0,isub)/float(ne/nsub-1)
            mom(ib)%m2_2(:,:,:,jl0,isub) = mom(ib)%m2_2(:,:,:,jl0,isub)/float(ne/nsub-1)
            mom(ib)%m11(:,:,:,jl0,isub) = mom(ib)%m11(:,:,:,jl0,isub)/float(ne/nsub-1)
            if (.not.nam%gau_approx) then
               mom(ib)%m12(:,:,:,jl0,isub) = mom(ib)%m12(:,:,:,jl0,isub)/float(ne/nsub)
               mom(ib)%m21(:,:,:,jl0,isub) = mom(ib)%m21(:,:,:,jl0,isub)/float(ne/nsub)
               mom(ib)%m22(:,:,:,jl0,isub) = mom(ib)%m22(:,:,:,jl0,isub)/float(ne/nsub)
            end if
            if (nam%full_var) mom(ib)%m2full(:,jl0,isub) = mom(ib)%m2full(:,jl0,isub)/float(ne/nsub-1)
        end do
        !$omp end parallel do
      end if
   end do
end do

! End associate
end associate

end subroutine compute_moments

end module module_moments
