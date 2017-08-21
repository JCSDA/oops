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

use model_interface, only: model_read
use module_namelist, only: nam
use module_parameters_convol, only: compute_convol_network,compute_convol_distance
use module_parameters_interp, only: compute_interp_h,compute_interp_v,compute_interp_s,interp_horiz
use netcdf
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: msgerror,msgwarning
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_mesh, only: meshtype,create_mesh
use type_mpl, only: mpl
use type_ndata, only: ndatatype
use type_randgen, only: initialize_sampling,rand_integer

implicit none

private
public :: compute_parameters

contains

!----------------------------------------------------------------------
! Subroutine: compute_parameters
!> Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine compute_parameters(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: il0,il0_prev,il1,ic0,ic1,ic2,i,is,il0i
integer :: mask_cell(ndata%nc0)
integer,allocatable :: ifors(:)
real(kind_real) :: distnorm
real(kind_real) :: rh0(ndata%nc0,ndata%nl0),rv0(ndata%nc0,ndata%nl0),rh0min(ndata%nc0)
real(kind_real),allocatable :: rh1(:,:),rv1(:,:),rh2(:,:),rv2(:,:),rhs(:),rvs(:)

! Compute grid mesh
write(mpl%unit,'(a7,a)') '','Compute grid mesh'
call compute_grid_mesh(ndata)

! Forced points in the subrid (TODO: rethink that)
ndata%nfor = 1
allocate(ndata%ifor(ndata%nfor))
do ic0=1,ndata%nc0
   if (all(ndata%mask(ic0,:))) then
      ndata%ifor(1) = ic0
      exit
   end if
end do

! Compute adaptive sampling
write(mpl%unit,'(a7,a)') '','Compute adaptive sampling'

! Read length-scales and transform to Gaspari-Cohn (1999) function parameters
if (trim(nam%Lbh_file)=='none') then
   do il0=1,ndata%nl0
      rh0(:,il0) = nam%Lbh(il0)
   end do
else
   call model_read(trim(nam%Lbh_file),'Lh',ndata,rh0)
end if
if (trim(nam%Lbv_file)=='none') then
   do il0=1,ndata%nl0
      rv0(:,il0) = nam%Lbv(il0)
   end do
else
   call model_read(trim(nam%Lbv_file),'Lv',ndata,rv0)
end if
write(mpl%unit,'(a10,a,e10.3,a,e10.3)') '','Average length-scales (H/V): ', &
 & sum(rh0)/float(ndata%nc0*ndata%nl0),' / ',sum(rv0)/float(ndata%nc0*ndata%nl0)

! Compute scaled support radii
rh0 = 2.0*rh0/(sqrt(0.3)*req)
rv0 = 2.0*rv0/sqrt(0.3)
write(mpl%unit,'(a10,a,e10.3,a,e10.3)') '','Average scaled support radii (H/V): ', &
 & sum(rh0)/float(ndata%nc0*ndata%nl0),' / ',sum(rv0)/float(ndata%nc0*ndata%nl0)

! Basic horizontal mesh defined with the minimum support radius
rh0min = minval(rh0,dim=2)
ndata%nc1 = floor(2.0*maxval(ndata%area)*nam%resol**2/(sqrt(3.0)*(sum(rh0min,mask=isnotmsr(rh0min)) &
          & /float(count(isnotmsr(rh0min))))**2))+ndata%nfor
if (ndata%nc1>ndata%nc0/2) then
   call msgwarning('required nc1 larger than nc0/2, resetting to nc0/2')
   ndata%nc1 = ndata%nc0/2
end if
mask_cell = 0
do ic0=1,ndata%nc0
   if (any(ndata%mask(ic0,:))) mask_cell(ic0) = 1
end do
allocate(ndata%ic1_to_ic0(ndata%nc1))

! Compute subset
write(mpl%unit,'(a7,a)') '','Compute horizontal subset C1'
call initialize_sampling(ndata%rng,ndata%nc0,dble(ndata%lon),dble(ndata%lat),mask_cell,rh0min,nam%ntry,nam%nrep, &
 & ndata%nc1,ndata%nfor,ndata%ifor,ndata%ic1_to_ic0)

! Inverse conversion
allocate(ndata%ic0_to_ic1(ndata%nc0))
call msi(ndata%ic0_to_ic1)
do ic1=1,ndata%nc1
   ic0 = ndata%ic1_to_ic0(ic1)
   ndata%ic0_to_ic1(ic0) = ic1
end do

! Adapt support radii
allocate(rh1(ndata%nc1,ndata%nl0))
allocate(rv1(ndata%nc1,ndata%nl0))
do ic1=1,ndata%nc1
   rh1(ic1,:) = rh0(ndata%ic1_to_ic0(ic1),:)
   rv1(ic1,:) = rv0(ndata%ic1_to_ic0(ic1),:)
end do

! Vertical sampling
write(mpl%unit,'(a7,a)',advance='no') '','Compute vertical subset L1: '
allocate(ndata%llev(ndata%nl0))
il0_prev = 1
do il0=1,ndata%nl0
   ! Look for convolution levels
   if ((il0==1).or.(il0==ndata%nl0)) then
      ! Keep first and last levels
      ndata%llev(il0) = .true.
   else
      ! Compute normalized distance with level il0_prev
      distnorm = abs(ndata%vunit(il0)-ndata%vunit(il0_prev))/sqrt(0.5*(minval(rv1(:,il0))**2+minval(rv1(:,il0_prev))**2))
      ndata%llev(il0) = distnorm>1.0/nam%resol
   end if

   ! Update
   if (ndata%llev(il0)) il0_prev = il0
end do
ndata%nl1 = count(ndata%llev)
allocate(ndata%il1_to_il0(ndata%nl1))
il1 = 0
do il0=1,ndata%nl0
   if (ndata%llev(il0)) then
      write(mpl%unit,'(i3,a)',advance='no') nam%levs(il0),' '
      il1 = il1+1
      ndata%il1_to_il0(il1) = il0
   end if
end do
write(mpl%unit,'(a)') ''

! Inverse conversion
allocate(ndata%il0_to_il1(ndata%nl0))
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
allocate(ifors(ndata%nfor))
do i=1,ndata%nfor
   ifors(i) = i
end do
do il1=1,ndata%nl1
   write(mpl%unit,'(a7,a,i3,a)') '','Compute horizontal subset C2 (level ',il1,')'
   il0 = ndata%il1_to_il0(il1)
   ndata%nc2(il1) = floor(2.0*ndata%area(il0)*nam%resol**2/(sqrt(3.0)*(sum(rh2(:,il1),mask=isnotmsr(rh2(:,il1))) &
                  & /float(count(isnotmsr(rh2(:,il1)))))**2))+ndata%nfor
   ndata%nc2(il1) = min(ndata%nc2(il1),ndata%nc1)
   if (ndata%nc2(il1)<ndata%nc1) then
      ! Compute subset

      call initialize_sampling(ndata%rng,ndata%nc1,dble(ndata%lon(ndata%ic1_to_ic0)), &
    & dble(ndata%lat(ndata%ic1_to_ic0)), mask_cell(ndata%ic1_to_ic0),rh0(ndata%ic1_to_ic0,il0),nam%ntry,nam%nrep, &
    & ndata%nc2(il1),ndata%nfor,ifors,ndata%ic2il1_to_ic1(1:ndata%nc2(il1),il1))
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
allocate(ndata%ic0il0_to_is(ndata%nc0,ndata%nl0))
allocate(ndata%ic2il1_to_ic0(ndata%nc1,ndata%nl1))
allocate(ndata%ic2il1_to_is(ndata%nc1,ndata%nl1))
allocate(ndata%ic1il1_to_is(ndata%nc1,ndata%nl1))
call msi(ndata%ic0il0_to_is)
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
call compute_interp_v(ndata)

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
write(mpl%unit,'(a10,a,i8)') '','nc0 =       ',ndata%nc0
write(mpl%unit,'(a10,a,i8)') '','nl0 =       ',ndata%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =       ',ndata%nc1
write(mpl%unit,'(a10,a,i8)') '','nl1 =       ',ndata%nl1
do il1=1,ndata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',ndata%nc2(il1)
end do
write(mpl%unit,'(a10,a,i8)') '','ns =        ',ndata%ns
do il0i=1,ndata%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',ndata%h(il0i)%n_s
end do
do il0i=1,ndata%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','v(',il0i,')%n_s = ',ndata%v(il0i)%n_s
end do
do il1=1,ndata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',ndata%s(il1)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','c%n_s =     ',ndata%c%n_s

end subroutine compute_parameters

!----------------------------------------------------------------------
! Subroutine: compute_grid_mesh
!> Purpose: compute grid mesh
!----------------------------------------------------------------------
subroutine compute_grid_mesh(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: nc0,info,ic0,jc0,kc0,i,ibnd,il0,nt,it
integer,allocatable :: ltri(:,:),ic0_bnd(:,:,:)
real(kind_real) :: area,areas,frac,latbnd(2),lonbnd(2),v1(3),v2(3)
logical :: init
type(meshtype) :: mesh

! Create mesh
if ((.not.all(ndata%area>0.0)).or.nam%mask_check.or.nam%network) &
 & call create_mesh(ndata%rng,ndata%nc0,ndata%lon,ndata%lat,.true.,mesh)

if ((.not.all(ndata%area>0.0))) then
   ! Allocation
   allocate(ltri(6,2*(mesh%nnr-2)))

   ! Create triangles list
   call trlist(mesh%nnr,mesh%list,mesh%lptr,mesh%lend,6,nt,ltri,info)

   ! Compute area
   ndata%area = 0.0
   do it=1,nt
      area = areas((/mesh%x(ltri(1,it)),mesh%y(ltri(1,it)),mesh%z(ltri(1,it))/), &
                 & (/mesh%x(ltri(2,it)),mesh%y(ltri(2,it)),mesh%z(ltri(2,it))/), &
                 & (/mesh%x(ltri(3,it)),mesh%y(ltri(3,it)),mesh%z(ltri(3,it))/))
      do il0=1,ndata%nl0
         frac = float(count(ndata%mask(mesh%order(ltri(1:3,it)),il0)))/3.0
         ndata%area(il0) = ndata%area(il0)+frac*area
      end do  
   end do

   ! Release memory
   deallocate(ltri)
end if

if (nam%mask_check) then
   ! Allocation
   allocate(ndata%nbnd(ndata%nl0))
   allocate(ic0_bnd(2,mesh%nnr,ndata%nl0))
   
   ! Find border points
   do il0=1,ndata%nl0
      ndata%nbnd(il0) = 0
      do ic0=1,mesh%nnr
         ! Check mask points only
         if (.not.ndata%mask(mesh%order(ic0),il0)) then
            i = mesh%lend(ic0)
            init = .true.
            do while ((i/=mesh%lend(ic0)).or.init)
               jc0 = abs(mesh%list(i))
               kc0 = abs(mesh%list(mesh%lptr(i)))
               if (.not.ndata%mask(mesh%order(jc0),il0).and.ndata%mask(mesh%order(kc0),il0)) then
                  ! Create a new boundary arc
                  ndata%nbnd(il0) = ndata%nbnd(il0)+1
                  if (ndata%nbnd(il0)>mesh%nnr) call msgerror('too many boundary arcs')
                  ic0_bnd(1,ndata%nbnd(il0),il0) = mesh%order(ic0)
                  ic0_bnd(2,ndata%nbnd(il0),il0) = mesh%order(jc0)
               end if
               i = mesh%lptr(i)
               init = .false.
            end do
         end if
      end do
   end do
   
   ! Allocation
   allocate(ndata%xbnd(2,maxval(ndata%nbnd),ndata%nl0))
   allocate(ndata%ybnd(2,maxval(ndata%nbnd),ndata%nl0))
   allocate(ndata%zbnd(2,maxval(ndata%nbnd),ndata%nl0))
   allocate(ndata%vbnd(3,maxval(ndata%nbnd),ndata%nl0))
   
   do il0=1,ndata%nl0
      ! Compute boundary arcs
      do ibnd=1,ndata%nbnd(il0)
         latbnd = ndata%lat(ic0_bnd(:,ibnd,il0))
         lonbnd = ndata%lon(ic0_bnd(:,ibnd,il0))
         call trans(2,latbnd,lonbnd,ndata%xbnd(:,ibnd,il0),ndata%ybnd(:,ibnd,il0),ndata%zbnd(:,ibnd,il0))
      end do
      do ibnd=1,ndata%nbnd(il0)
         v1 = (/ndata%xbnd(1,ibnd,il0),ndata%ybnd(1,ibnd,il0),ndata%zbnd(1,ibnd,il0)/)
         v2 = (/ndata%xbnd(2,ibnd,il0),ndata%ybnd(2,ibnd,il0),ndata%zbnd(2,ibnd,il0)/)
         call vector_product(v1,v2,ndata%vbnd(:,ibnd,il0))
      end do
   end do
end if


if (nam%network) then
   if (.not.allocated(ndata%net_nnb)) then
      ! Allocation
      allocate(ndata%net_nnb(ndata%nc0))
   
      ! Count neighbors
      ndata%net_nnb = 0
      do ic0=1,mesh%nnr
         i = mesh%lend(ic0)
         init = .true.
         do while ((i/=mesh%lend(ic0)).or.init)
            ndata%net_nnb(mesh%order(ic0)) = ndata%net_nnb(mesh%order(ic0))+1
            i = mesh%lptr(i)
            init = .false.
         end do
      end do
   
      ! Allocation
      allocate(ndata%net_inb(maxval(ndata%net_nnb),ndata%nc0))
   
      ! Find neighbors
      ndata%net_nnb = 0
      do ic0=1,mesh%nnr
         i = mesh%lend(ic0)
         init = .true.
         do while ((i/=mesh%lend(ic0)).or.init)
            ndata%net_nnb(mesh%order(ic0)) = ndata%net_nnb(mesh%order(ic0))+1
            ndata%net_inb(ndata%net_nnb(mesh%order(ic0)),mesh%order(ic0)) = mesh%order(abs(mesh%list(i)))
            i = mesh%lptr(i)
            init = .false.
         end do
      end do
   
      ! Copy neighbors for redudant points
      do ic0=1,ndata%nc0
         if (isnotmsi(mesh%redundant(ic0))) then
            ndata%net_nnb(ic0) = ndata%net_nnb(mesh%redundant(ic0))
            ndata%net_inb(:,ic0) = ndata%net_inb(:,mesh%redundant(ic0))
         end if
      end do
   end if

   ! Compute distances
   allocate(ndata%net_dnb(maxval(ndata%net_nnb),ndata%nc0))
   do ic0=1,ndata%nc0
      do i=1,ndata%net_nnb(ic0)
         call sphere_dist(ndata%lon(ic0),ndata%lat(ic0),ndata%lon(ndata%net_inb(i,ic0)), &
       & ndata%lat(ndata%net_inb(i,ic0)),ndata%net_dnb(i,ic0))
         ndata%net_dnb(i,ic0) = (ndata%net_dnb(i,ic0)/req)**2
      end do
   end do
end if

end subroutine compute_grid_mesh

end module module_parameters
