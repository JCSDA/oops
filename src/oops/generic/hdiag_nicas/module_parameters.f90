!----------------------------------------------------------------------
! Module: module_parameters.f90
!> Purpose: compute NICAS parameters
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_parameters

use module_parameters_convol, only: compute_convol_network,compute_convol_distance
use module_parameters_interp, only: compute_interp_h,compute_interp_v,compute_interp_s
use netcdf
use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: aqua,black,vunitchar,msgerror,msgwarning
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_bdata, only: bdatatype
use type_mpl, only: mpl
use type_ndata, only: ndatatype
use type_randgen, only: initialize_sampling

implicit none

integer,parameter :: nc1max = 15000

private
public :: compute_parameters

contains

!----------------------------------------------------------------------
! Subroutine: compute_parameters
!> Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine compute_parameters(bdata,ndata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata    !< B data
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: il0,il0_prev,il1,ic0,ic1,ic2,is,il0i
integer :: mask_ind_col(ndata%geom%nc0)
integer,allocatable :: mask_ind(:,:)
real(kind_real) :: rh0minavg,rv1min,distnorm
real(kind_real) :: rh0(ndata%geom%nc0,ndata%geom%nl0),rv0(ndata%geom%nc0,ndata%geom%nl0),rh0min(ndata%geom%nc0)
real(kind_real),allocatable :: rh1(:,:),rv1(:,:),rh2(:,:),rv2(:,:),rhs(:),rvs(:)
logical :: llev(ndata%geom%nl0),inside

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Compute adaptive sampling
write(mpl%unit,'(a7,a)') '','Compute adaptive sampling'

! Compute scaled support radii
rh0 = bdata%rh0
rv0 = bdata%rv0
write(mpl%unit,'(a10,a,a,f8.2,a,f8.2,a)') '','Average support radii (H/V): ', &
 & trim(aqua),sum(rh0)/float(geom%nc0*geom%nl0)*reqkm,trim(black)//' km  / ' &
 & //trim(aqua),sum(rv0)/float(geom%nc0*geom%nl0),trim(black)//' '//trim(vunitchar)

! Basic horizontal mesh defined with the minimum support radius
rh0min = minval(rh0,dim=2)
rh0minavg = sum(rh0min,mask=isnotmsr(rh0min))/float(count(isnotmsr(rh0min)))
if (rh0minavg>0.0) then
   ndata%nc1 = floor(2.0*maxval(geom%area)*nam%resol**2/(sqrt(3.0)*rh0minavg**2))
else
   ndata%nc1 = ndata%geom%nc0
end if
write(mpl%unit,'(a10,a,i8)') '','Estimated nc1 from horizontal support radius: ',ndata%nc1
if (ndata%nc1>nc1max) then
   call msgwarning('required nc1 larger than nc1max, resetting to nc1max')
   ndata%nc1 = nc1max
   write(mpl%unit,'(a10,a,f5.2)') '','Effective resolution: ',sqrt(float(ndata%nc1)*sqrt(3.0)*rh0minavg**2/(2.0*maxval(geom%area)))
end if
mask_ind_col = 0
do ic0=1,geom%nc0
   if (any(geom%mask(ic0,:))) mask_ind_col(ic0) = 1
end do

! Compute subset
write(mpl%unit,'(a7,a)') '','Compute horizontal subset C1'
allocate(ndata%ic1_to_ic0(ndata%nc1))
call initialize_sampling(geom%nc0,dble(geom%lon),dble(geom%lat),mask_ind_col,rh0min,nam%ntry,nam%nrep, &
 & ndata%nc1,ndata%ic1_to_ic0)

! Inverse conversion
allocate(ndata%ic0_to_ic1(geom%nc0))
call msi(ndata%ic0_to_ic1)
do ic1=1,ndata%nc1
   ic0 = ndata%ic1_to_ic0(ic1)
   ndata%ic0_to_ic1(ic0) = ic1
end do

! Adapt support radii
allocate(rh1(ndata%nc1,geom%nl0))
allocate(rv1(ndata%nc1,geom%nl0))
do ic1=1,ndata%nc1
   rh1(ic1,:) = rh0(ndata%ic1_to_ic0(ic1),:)
   rv1(ic1,:) = rv0(ndata%ic1_to_ic0(ic1),:)
end do

! Vertical sampling
write(mpl%unit,'(a7,a)',advance='no') '','Compute vertical subset L1: '
il0_prev = 1
do il0=1,geom%nl0
   ! Look for convolution levels
   if ((il0==1).or.(il0==geom%nl0)) then
      ! Keep first and last levels
      llev(il0) = .true.
   else
      ! Compute normalized distance with level il0_prev
      rv1min = sqrt(0.5*(minval(rv1(:,il0))**2+minval(rv1(:,il0_prev))**2))
      if (rv1min>0.0) then
         distnorm = abs(geom%vunit(il0)-geom%vunit(il0_prev))/rv1min
         llev(il0) = distnorm>1.0/nam%resol
      else
         llev(il0) = .true.
      end if
   end if

   ! Update
   if (llev(il0)) il0_prev = il0
end do
ndata%nl1 = count(llev)
allocate(ndata%il1_to_il0(ndata%nl1))
il1 = 0
do il0=1,geom%nl0
   if (llev(il0)) then
      write(mpl%unit,'(i3,a)',advance='no') nam%levs(il0),' '
      il1 = il1+1
      ndata%il1_to_il0(il1) = il0
   end if
end do
write(mpl%unit,'(a)') ''

! Find bottom and top for each point of S1
allocate(ndata%vbot(ndata%nc1))
allocate(ndata%vtop(ndata%nc1))
!$omp parallel do schedule(static) private(ic1,ic0,inside,il1,il0)
do ic1=1,ndata%nc1
   ic0 = ndata%ic1_to_ic0(ic1)
   inside = .false.
   ndata%vtop(ic1) = geom%nl0
   do il1=1,ndata%nl1
      il0 = ndata%il1_to_il0(il1)
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
allocate(ndata%il0_to_il1(geom%nl0))
call msi(ndata%il0_to_il1)
do il1=1,ndata%nl1
   il0 = ndata%il1_to_il0(il1)
   ndata%il0_to_il1(il0) = il1
end do

! Adapt support radii
allocate(rh2(ndata%nc1,ndata%nl1))
allocate(rv2(ndata%nc1,ndata%nl1))
do il1=1,ndata%nl1
   rh2(:,il1) = rh1(:,ndata%il1_to_il0(il1))
   rv2(:,il1) = rv1(:,ndata%il1_to_il0(il1))
end do

! Horizontal subsampling
allocate(ndata%nc2(ndata%nl1))
allocate(ndata%ic2il1_to_ic1(ndata%nc1,ndata%nl1))
allocate(mask_ind(ndata%nc1,ndata%nl1))
call msi(ndata%ic2il1_to_ic1)
mask_ind = 0
do il1=1,ndata%nl1
   il0 = ndata%il1_to_il0(il1)
   do ic1=1,ndata%nc1
      ic0 = ndata%ic1_to_ic0(ic1)
      if (geom%mask(ic0,il0)) mask_ind(ic1,il1) = 1
   end do
end do

do il1=1,ndata%nl1
   write(mpl%unit,'(a7,a,i3,a)') '','Compute horizontal subset C2 (level ',il1,')'
   il0 = ndata%il1_to_il0(il1)
   ndata%nc2(il1) = floor(2.0*geom%area(il0)*nam%resol**2/(sqrt(3.0)*(sum(rh2(:,il1),mask=isnotmsr(rh2(:,il1))) &
                  & /float(count(isnotmsr(rh2(:,il1)))))**2))
   ndata%nc2(il1) = min(ndata%nc2(il1),ndata%nc1)
   if (ndata%nc2(il1)<ndata%nc1) then
      ! Compute subset
      call initialize_sampling(ndata%nc1,dble(geom%lon(ndata%ic1_to_ic0)), &
    & dble(geom%lat(ndata%ic1_to_ic0)),mask_ind(:,il1),rh0(ndata%ic1_to_ic0,il0),nam%ntry,nam%nrep, &
    & ndata%nc2(il1),ndata%ic2il1_to_ic1(1:ndata%nc2(il1),il1))
   else
      do ic2=1,ndata%nc2(il1)
         ndata%ic2il1_to_ic1(ic2,il1) = ic2
      end do
   end if
end do

! Final conversions
ndata%ns = sum(ndata%nc2)
allocate(ndata%is_to_ic1(ndata%ns))
allocate(ndata%is_to_il1(ndata%ns))
allocate(ndata%is_to_ic2(ndata%ns))
allocate(ndata%ic0il0_to_is(geom%nc0,geom%nl0))
allocate(ndata%ic2il1_to_ic0(ndata%nc1,ndata%nl1))
allocate(ndata%ic2il1_to_is(ndata%nc1,ndata%nl1))
allocate(ndata%ic1il1_to_is(ndata%nc1,ndata%nl1))
call msi(ndata%ic0il0_to_is)
call msi(ndata%ic2il1_to_ic0)
call msi(ndata%ic2il1_to_is)
call msi(ndata%ic1il1_to_is)
is = 0
do il1=1,ndata%nl1
   do ic2=1,ndata%nc2(il1)
      is = is+1
      ic1 = ndata%ic2il1_to_ic1(ic2,il1)
      ic0 = ndata%ic1_to_ic0(ic1)
      il0 = ndata%il1_to_il0(il1)

      ndata%is_to_ic1(is) = ic1
      ndata%is_to_il1(is) = il1
      ndata%is_to_ic2(is) = ic2
      ndata%ic0il0_to_is(ic0,il0) = is
      ndata%ic2il1_to_ic0(ic2,il1) = ic0
      ndata%ic2il1_to_is(ic2,il1) = is
      ndata%ic1il1_to_is(ic1,il1) = is
   end do
end do

! Adapt support radii
allocate(rhs(ndata%ns))
allocate(rvs(ndata%ns))
do is=1,ndata%ns
   rhs(is) = rh2(ndata%is_to_ic1(is),ndata%is_to_il1(is))
   rvs(is) = rv2(ndata%is_to_ic1(is),ndata%is_to_il1(is))
end do

! Compute horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute horizontal interpolation data'
call compute_interp_h(ndata)

! Compute vertical interpolation data
write(mpl%unit,'(a7,a)') '','Compute vertical interpolation data'
call compute_interp_v(ndata,llev)

! Compute subsampling horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute subsampling horizontal interpolation data'
call compute_interp_s(ndata)

! Compute convolution data
write(mpl%unit,'(a7,a)') '','Compute convolution data'
if (nam%network) then
   call compute_convol_network(ndata,rh0,rv0)
else
   call compute_convol_distance(ndata,rhs,rvs)
end if

! Print results
write(mpl%unit,'(a7,a)') '','Global parameters'
write(mpl%unit,'(a10,a,i8)') '','nc0 =       ',geom%nc0
write(mpl%unit,'(a10,a,i8)') '','nl0 =       ',geom%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =       ',ndata%nc1
write(mpl%unit,'(a10,a,i8)') '','nl1 =       ',ndata%nl1
do il1=1,ndata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',ndata%nc2(il1)
end do
write(mpl%unit,'(a10,a,i8)') '','ns =        ',ndata%ns
do il0i=1,geom%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',ndata%h(il0i)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','v%n_s =     ',ndata%v%n_s
do il1=1,ndata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',ndata%s(il1)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','c%n_s =     ',ndata%c%n_s

! End associate
end associate

end subroutine compute_parameters

end module module_parameters
