!----------------------------------------------------------------------
! Module: nicas_parameters_sampling.f90
!> Purpose: compute NICAS parameters
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_parameters_sampling

use netcdf
use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: aqua,black,vunitchar,msgerror,msgwarning
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_bdata, only: bdatatype
use type_mpl, only: mpl,mpl_bcast
use type_ndata, only: ndatatype
use type_randgen, only: initialize_sampling,reseed_randgen

implicit none

integer,parameter :: nc1max = 15000 !< Maximum size of the Sc1 subset

private
public :: compute_sampling

contains

!----------------------------------------------------------------------
! Subroutine: compute_sampling
!> Purpose: compute NICAS sampling
!----------------------------------------------------------------------
subroutine compute_sampling(bdata,ndata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata    !< B data
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: il0,il0_prev,il1,ic0,ic1,ic2,is
integer :: mask_c0(ndata%geom%nc0)
integer,allocatable :: mask_c1(:),c2_to_c1(:)
real(kind_real) :: rh0s(ndata%geom%nc0,ndata%geom%nl0),rv0s(ndata%geom%nc0,ndata%geom%nl0)
real(kind_real) :: rh0minavg,rv1min,distnorm
real(kind_real) :: rh0min(ndata%geom%nc0)
logical :: inside

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Copy
if (allocated(bdata%rh0s).and.allocated(bdata%rv0s)) then
   ! Specific support radii for sampling
   rh0s = bdata%rh0s
   rv0s = bdata%rv0s
else
   ! Same support radii for sampling and convolution
   rh0s = bdata%rh0
   rv0s = bdata%rv0
end if

! Reset random numbers seed
if (trim(nam%strategy)=='specific_multivariate') call reseed_randgen()

! Compute support radii
write(mpl%unit,'(a10,a,a,f8.2,a,f8.2,a)') '','Average support radii (H/V): ', &
 & trim(aqua),sum(rh0s)/float(geom%nc0*geom%nl0)*reqkm,trim(black)//' km  / ' &
 & //trim(aqua),sum(rv0s)/float(geom%nc0*geom%nl0),trim(black)//' '//trim(vunitchar)

! Basic horizontal mesh defined with the minimum support radius
rh0min = minval(rh0s,dim=2)
rh0minavg = sum(rh0min,mask=isnotmsr(rh0min))/float(count(isnotmsr(rh0min)))
if (rh0minavg>0.0) then
   ndata%nc1 = floor(2.0*maxval(geom%area)*nam%resol**2/(sqrt(3.0)*rh0minavg**2))
else
   ndata%nc1 = geom%nc0
end if
ndata%nc1 = min(ndata%nc1,geom%nc0)
write(mpl%unit,'(a10,a,i8)') '','Estimated nc1 from horizontal support radius: ',ndata%nc1
if (ndata%nc1>nc1max) then
   call msgwarning('required nc1 larger than nc1max, resetting to nc1max')
   ndata%nc1 = nc1max
   write(mpl%unit,'(a10,a,f5.2)') '','Effective resolution: ',sqrt(float(ndata%nc1)*sqrt(3.0)*rh0minavg**2/(2.0*maxval(geom%area)))
end if
mask_c0 = 0
do ic0=1,geom%nc0
   if (any(geom%mask(ic0,:))) mask_c0(ic0) = 1
end do

! Compute subset
write(mpl%unit,'(a7,a)') '','Compute horizontal subset C1'
allocate(ndata%c1_to_c0(ndata%nc1))
if (mpl%main) call initialize_sampling(geom%nc0,dble(geom%lon),dble(geom%lat),mask_c0,rh0min,nam%ntry,nam%nrep, &
 & ndata%nc1,ndata%c1_to_c0)
call mpl_bcast(ndata%c1_to_c0,mpl%ioproc)

! Inverse conversion
allocate(ndata%c0_to_c1(geom%nc0))
call msi(ndata%c0_to_c1)
do ic1=1,ndata%nc1
   ic0 = ndata%c1_to_c0(ic1)
   ndata%c0_to_c1(ic0) = ic1
end do

! Vertical sampling
write(mpl%unit,'(a7,a)',advance='no') '','Compute vertical subset L1: '
allocate(ndata%llev(ndata%geom%nl0))
il0_prev = 1
do il0=1,geom%nl0
   ! Look for convolution levels
   if ((il0==1).or.(il0==geom%nl0)) then
      ! Keep first and last levels
      ndata%llev(il0) = .true.
   else
      ! Compute normalized distance with level il0_prev
      rv1min = sqrt(0.5*(minval(rv0s(ndata%c1_to_c0,il0))**2+minval(rv0s(ndata%c1_to_c0,il0_prev))**2))
      if (rv1min>0.0) then
         distnorm = abs(geom%vunit(il0)-geom%vunit(il0_prev))/rv1min
         ndata%llev(il0) = distnorm>1.0/nam%resol
      else
         ndata%llev(il0) = .true.
      end if
   end if

   ! Update
   if (ndata%llev(il0)) il0_prev = il0
end do
ndata%nl1 = count(ndata%llev)
allocate(ndata%l1_to_l0(ndata%nl1))
il1 = 0
do il0=1,geom%nl0
   if (ndata%llev(il0)) then
      write(mpl%unit,'(i3,a)',advance='no') nam%levs(il0),' '
      il1 = il1+1
      ndata%l1_to_l0(il1) = il0
   end if
end do
write(mpl%unit,'(a)') ''

! Find bottom and top for each point of S1
allocate(ndata%vbot(ndata%nc1))
allocate(ndata%vtop(ndata%nc1))
!$omp parallel do schedule(static) private(ic1,ic0,inside,il1,il0)
do ic1=1,ndata%nc1
   ic0 = ndata%c1_to_c0(ic1)
   inside = .false.
   ndata%vtop(ic1) = geom%nl0
   do il1=1,ndata%nl1
      il0 = ndata%l1_to_l0(il1)
      if (.not.inside.and.geom%mask(ic0,il0)) then
         ! Bottom level
         ndata%vbot(ic1) = il0
         inside = .true.
      end if
      if (inside.and.(.not.geom%mask(ic0,il0))) then
         ! Top level
         ndata%vtop(ic1) = il0
         inside = .false.
      end if
   end do
   if (ndata%vbot(ic1)>ndata%vtop(ic1)) call msgerror('non contiguous mask')
end do
!$omp end parallel do

! Inverse conversion
allocate(ndata%l0_to_l1(geom%nl0))
call msi(ndata%l0_to_l1)
do il1=1,ndata%nl1
   il0 = ndata%l1_to_l0(il1)
   ndata%l0_to_l1(il0) = il1
end do

! Horizontal subsampling
allocate(ndata%nc2(ndata%nl1))
allocate(ndata%c2mask(ndata%nc1,ndata%nl1))
allocate(mask_c1(ndata%nc1))
do il1=1,ndata%nl1
   write(mpl%unit,'(a7,a,i3,a)') '','Compute horizontal subset C2 (level ',il1,')'

   ! Compute nc2
   il0 = ndata%l1_to_l0(il1)
   ndata%nc2(il1) = floor(2.0*geom%area(il0)*nam%resol**2/(sqrt(3.0)*(sum(rh0s(ndata%c1_to_c0,ndata%l1_to_l0(il1)), &
                  & mask=isnotmsr(rh0s(ndata%c1_to_c0,ndata%l1_to_l0(il1)))) &
                  & /float(count(isnotmsr(rh0s(ndata%c1_to_c0,ndata%l1_to_l0(il1))))))**2))
   ndata%nc2(il1) = max(ndata%nc1/4,min(ndata%nc2(il1),ndata%nc1))

   if (ndata%nc2(il1)<ndata%nc1) then
      ! Allocation
      allocate(c2_to_c1(ndata%nc2(il1)))

      ! Mask
      mask_c1 = 0
      do ic1=1,ndata%nc1
         ic0 = ndata%c1_to_c0(ic1)
         if (geom%mask(ic0,il0)) mask_c1(ic1) = 1
      end do

      ! Compute subset
      if (mpl%main) call initialize_sampling(ndata%nc1,dble(geom%lon(ndata%c1_to_c0)), &
    & dble(geom%lat(ndata%c1_to_c0)),mask_c1,rh0s(ndata%c1_to_c0,il0),nam%ntry,nam%nrep, &
    & ndata%nc2(il1),c2_to_c1)
      call mpl_bcast(c2_to_c1,mpl%ioproc)

      ! Fill C2 mask
      ndata%c2mask(:,il1) = .false.
      do ic2=1,ndata%nc2(il1)
         ic1 = c2_to_c1(ic2)
         ndata%c2mask(ic1,il1) = .true.
      end do

      ! Release memory
      deallocate(c2_to_c1)
   else
      ! Fill C2 mask
      ndata%c2mask(:,il1) = .true.
   end if
end do

! Final conversions
ndata%ns = sum(ndata%nc2)
allocate(ndata%s_to_c1(ndata%ns))
allocate(ndata%s_to_l1(ndata%ns))
allocate(ndata%c1l1_to_s(ndata%nc1,ndata%nl1))
call msi(ndata%c1l1_to_s)
is = 0
do il1=1,ndata%nl1
   do ic1=1,ndata%nc1
      if (ndata%c2mask(ic1,il1)) then
         is = is+1
         ndata%s_to_c1(is) = ic1
         ndata%s_to_l1(is) = il1
         ndata%c1l1_to_s(ic1,il1) = is
      end if
   end do
end do

! End associate
end associate

end subroutine compute_sampling

end module nicas_parameters_sampling
