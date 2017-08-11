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

use model_interface, only: model_read,model_write
use module_namelist, only: nam
use module_parameters_convol, only: compute_convol_network,compute_convol_distance
use module_parameters_interp, only: compute_interp_h,compute_interp_v,compute_interp_s
use netcdf
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: msgerror,msgwarning
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_mpl, only: mpl,mpl_bcast
use type_randgen, only: initialize_sampling,rand_integer
use type_sdata, only: sdatatype
implicit none

private
public :: compute_parameters

contains

!----------------------------------------------------------------------
! Subroutine: compute_parameters
!> Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine compute_parameters(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0,il0_prev,il1,ic0,ic1,ic2,i,is,il0i
integer :: mask_cell(sdata%nc0)
integer :: ifors(sdata%nfor)
real(kind_real) :: distnorm
real(kind_real) :: rh0(sdata%nc0,sdata%nl0),rv0(sdata%nc0,sdata%nl0),rh0min(sdata%nc0)
real(kind_real),allocatable :: rh1(:,:),rv1(:,:),rh2(:,:),rv2(:,:),rhs(:),rvs(:)

! Forced points in the subrid (TODO: rethink that)
sdata%nfor = 0
allocate(sdata%ifor(sdata%nfor))

! Compute adaptive sampling
write(mpl%unit,'(a7,a)') '','Compute adaptive sampling'

! Read length-scales and transform to Gaspari-Cohn (1999) function parameters
if (trim(nam%Lbh_file)=='none') then
   do il0=1,sdata%nl0
      rh0(:,il0) = nam%Lbh(il0)
   end do
else
   call model_read(trim(nam%Lbh_file),'Lh',sdata,rh0)
end if
if (trim(nam%Lbv_file)=='none') then
   do il0=1,sdata%nl0
      rv0(:,il0) = nam%Lbv(il0)
   end do
else
   call model_read(trim(nam%Lbv_file),'Lv',sdata,rv0)
end if
write(mpl%unit,'(a10,a,e10.3,a,e10.3)') '','Average length-scales (H/V): ', &
 & sum(rh0)/float(sdata%nc0*sdata%nl0),' / ',sum(rv0)/float(sdata%nc0*sdata%nl0)

! Compute scaled support radii
rh0 = 2.0*rh0/(sqrt(0.3)*req)
rv0 = 2.0*rv0/sqrt(0.3)
write(mpl%unit,'(a10,a,e10.3,a,e10.3)') '','Average scaled support radii (H/V): ', &
 & sum(rh0)/float(sdata%nc0*sdata%nl0),' / ',sum(rv0)/float(sdata%nc0*sdata%nl0)

! Basic horizontal mesh defined with the minimum support radius
rh0min = minval(rh0,dim=2)
sdata%nc1 = floor(2.0*maxval(sdata%area)*nam%resol**2/(sqrt(3.0)*(sum(rh0min,mask=isnotmsr(rh0min)) &
          & /float(count(isnotmsr(rh0min))))**2))+sdata%nfor
if (sdata%nc1>sdata%nc0/2) then
   call msgwarning('required nc1 larger than nc0/2, resetting to nc0/2')
   sdata%nc1 = sdata%nc0/2
end if
mask_cell = 0
do ic0=1,sdata%nc0
   if (any(sdata%mask(ic0,:))) mask_cell(ic0) = 1
end do
allocate(sdata%ic1_to_ic0(sdata%nc1))

! Compute subset
write(mpl%unit,'(a7,a)') '','Compute horizontal subset C1'
call initialize_sampling(sdata%rng,sdata%nc0,dble(sdata%lon),dble(sdata%lat),mask_cell,rh0min,nam%ntry,nam%nrep, &
 & sdata%nc1,sdata%nfor,sdata%ifor,sdata%ic1_to_ic0)

! Inverse conversion
allocate(sdata%ic0_to_ic1(sdata%nc0))
call msi(sdata%ic0_to_ic1)
do ic1=1,sdata%nc1
   ic0 = sdata%ic1_to_ic0(ic1)
   sdata%ic0_to_ic1(ic0) = ic1
end do

! Adapt support radii
allocate(rh1(sdata%nc1,sdata%nl0))
allocate(rv1(sdata%nc1,sdata%nl0))
do ic1=1,sdata%nc1
   rh1(ic1,:) = rh0(sdata%ic1_to_ic0(ic1),:)
   rv1(ic1,:) = rv0(sdata%ic1_to_ic0(ic1),:)
end do

! Vertical sampling
write(mpl%unit,'(a7,a)',advance='no') '','Compute vertical subset L1: '
allocate(sdata%llev(sdata%nl0))
il0_prev = 1
do il0=1,sdata%nl0
   ! Look for convolution levels
   if ((il0==1).or.(il0==sdata%nl0)) then
      ! Keep first and last levels
      sdata%llev(il0) = .true.
   else
      ! Compute normalized distance with level il0_prev
      distnorm = abs(sdata%vunit(il0)-sdata%vunit(il0_prev))/sqrt(0.5*(minval(rv1(:,il0))**2+minval(rv1(:,il0_prev))**2))
      sdata%llev(il0) = distnorm>1.0/nam%resol
   end if

   ! Update
   if (sdata%llev(il0)) il0_prev = il0
end do
sdata%nl1 = count(sdata%llev)
allocate(sdata%il1_to_il0(sdata%nl1))
il1 = 0
do il0=1,sdata%nl0
   if (sdata%llev(il0)) then
      write(mpl%unit,'(i3,a)',advance='no') nam%levs(il0),' '
      il1 = il1+1
      sdata%il1_to_il0(il1) = il0
   end if
end do
write(mpl%unit,'(a)') ''

! Inverse conversion
allocate(sdata%il0_to_il1(sdata%nl0))
call msi(sdata%il0_to_il1)
do il1=1,sdata%nl1
   il0 = sdata%il1_to_il0(il1)
   sdata%il0_to_il1(il0) = il1
end do

! Adapt support radii
allocate(rh2(sdata%nc1,sdata%nl1))
allocate(rv2(sdata%nc1,sdata%nl1))
do il1=1,sdata%nl1
   rh2(:,il1) = rh1(:,sdata%il1_to_il0(il1))
   rv2(:,il1) = rv1(:,sdata%il1_to_il0(il1))
end do

! Horizontal subsampling
allocate(sdata%nc2(sdata%nl1))
allocate(sdata%ic2il1_to_ic1(sdata%nc1,sdata%nl1))
do il1=1,sdata%nl1
   write(mpl%unit,'(a7,a,i3,a)') '','Compute horizontal subset C2 (level ',il1,')'
   il0 = sdata%il1_to_il0(il1)
   sdata%nc2(il1) = floor(2.0*sdata%area(il0)*nam%resol**2/(sqrt(3.0)*(sum(rh2(:,il1),mask=isnotmsr(rh2(:,il1))) &
                  & /float(count(isnotmsr(rh2(:,il1)))))**2))+sdata%nfor
   sdata%nc2(il1) = min(sdata%nc2(il1),sdata%nc1)
   if (sdata%nc2(il1)<sdata%nc1) then
      ! Compute subset
      do i=1,sdata%nfor
         ifors(i) = i
      end do
      call initialize_sampling(sdata%rng,sdata%nc1,dble(sdata%lon(sdata%ic1_to_ic0)), &
    & dble(sdata%lat(sdata%ic1_to_ic0)), mask_cell(sdata%ic1_to_ic0),rh0(sdata%ic1_to_ic0,il0),nam%ntry,nam%nrep, &
    & sdata%nc2(il1),sdata%nfor,ifors,sdata%ic2il1_to_ic1(1:sdata%nc2(il1),il1))
   else
      do ic2=1,sdata%nc2(il1)
         sdata%ic2il1_to_ic1(ic2,il1) = ic2
      end do
   end if
end do

! Final conversions
sdata%ns = sum(sdata%nc2)
allocate(sdata%is_to_ic1(sdata%ns))
allocate(sdata%is_to_il1(sdata%ns))
allocate(sdata%is_to_ic2(sdata%ns))
allocate(sdata%ic0il0_to_is(sdata%nc0,sdata%nl0))
allocate(sdata%ic2il1_to_ic0(sdata%nc1,sdata%nl1))
allocate(sdata%ic2il1_to_is(sdata%nc1,sdata%nl1))
allocate(sdata%ic1il1_to_is(sdata%nc1,sdata%nl1))
call msi(sdata%ic0il0_to_is)
call msi(sdata%ic2il1_to_is)
call msi(sdata%ic1il1_to_is)
is = 0
do il1=1,sdata%nl1
   do ic2=1,sdata%nc2(il1)
      is = is+1
      ic1 = sdata%ic2il1_to_ic1(ic2,il1)
      ic0 = sdata%ic1_to_ic0(ic1)
      il0 = sdata%il1_to_il0(il1)

      sdata%is_to_ic1(is) = ic1
      sdata%is_to_il1(is) = il1
      sdata%is_to_ic2(is) = ic2
      sdata%ic0il0_to_is(ic0,il0) = is
      sdata%ic2il1_to_ic0(ic2,il1) = ic0
      sdata%ic2il1_to_is(ic2,il1) = is
      sdata%ic1il1_to_is(ic1,il1) = is
   end do
end do

! Adapt support radii
allocate(rhs(sdata%ns))
allocate(rvs(sdata%ns))
do is=1,sdata%ns
   rhs(is) = rh2(sdata%is_to_ic1(is),sdata%is_to_il1(is))
   rvs(is) = rv2(sdata%is_to_ic1(is),sdata%is_to_il1(is))
end do

! Compute grid mesh
write(mpl%unit,'(a7,a)') '','Compute grid mesh'
call compute_grid_mesh(sdata)

! Compute subgrid mesh
write(mpl%unit,'(a7,a)') '','Compute subgrid mesh'
call compute_subgrid_mesh(sdata)

! Compute horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute horizontal interpolation data'
call compute_interp_h(sdata)

! Compute vertical interpolation data
write(mpl%unit,'(a7,a)') '','Compute vertical interpolation data'
call compute_interp_v(sdata)

! Compute subsampling horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute subsampling horizontal interpolation data'
call compute_interp_s(sdata)

! Compute convolution data
write(mpl%unit,'(a7,a)') '','Compute convolution data'
if (nam%network) then
   call compute_convol_network(sdata,rh0,rv0)
else
   call compute_convol_distance(sdata,rhs,rvs)
end if

! Print results
write(mpl%unit,'(a7,a)') '','Global parameters'
write(mpl%unit,'(a10,a,i8)') '','nc0 =       ',sdata%nc0
write(mpl%unit,'(a10,a,i8)') '','nl0 =       ',sdata%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =       ',sdata%nc1
write(mpl%unit,'(a10,a,i8)') '','nl1 =       ',sdata%nl1
do il1=1,sdata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',sdata%nc2(il1)
end do
write(mpl%unit,'(a10,a,i8)') '','ns =        ',sdata%ns
do il0i=1,sdata%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',sdata%h(il0i)%n_s
end do
do il0i=1,sdata%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','v(',il0i,')%n_s = ',sdata%v(il0i)%n_s
end do
do il1=1,sdata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',sdata%s(il1)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','c%n_s =     ',sdata%c%n_s

end subroutine compute_parameters

!----------------------------------------------------------------------
! Subroutine: compute_grid_mesh
!> Purpose: compute grid mesh
!----------------------------------------------------------------------
subroutine compute_grid_mesh(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: nc0,lnew,info,ic0,jc0,kc0,i,ibnd,il0
integer :: redundant(sdata%nc0)
integer,allocatable :: order(:),list(:),lptr(:),lend(:),near(:),next(:),ic0_bnd(:,:,:)
real(kind_real) :: latbnd(2),lonbnd(2),v1(3),v2(3)
real(kind_real),allocatable :: x(:),y(:),z(:),dist(:)
logical :: init

! Look for redundant or masked points TODO : change that
call msi(redundant)
do ic0=1,sdata%nc0
   if (.not.isnotmsi(redundant(ic0))) then
      do jc0=ic0+1,sdata%nc0
         if ((abs(sdata%lon(ic0)-sdata%lon(jc0))<tiny(1.0)).and.(abs(sdata%lat(ic0)-sdata%lat(jc0))<tiny(1.0))) redundant(jc0) = ic0
      end do
   end if
end do
nc0 = count(.not.isnotmsi(redundant))

! Allocation
allocate(order(nc0))
allocate(list(6*(nc0-2)))
allocate(lptr(6*(nc0-2)))
allocate(lend(nc0))
allocate(near(nc0))
allocate(next(nc0))
allocate(x(nc0))
allocate(y(nc0))
allocate(z(nc0))
allocate(dist(nc0))

! Shuffle arrays (more efficient to compute the Delaunay triangulation)
ic0 = 0
do jc0=1,sdata%nc0
   if (.not.isnotmsi(redundant(jc0))) then
      ic0 = ic0+1
      order(ic0) = jc0
   end if
end do
do ic0=nc0,2,-1
   call rand_integer(sdata%rng,1,nc0,jc0)
   kc0 = order(jc0)
   order(jc0) = order(ic0)
   order(ic0) = kc0
end do

! Transform to cartesian coordinates
call trans(nc0,sdata%lat(order),sdata%lon(order),x,y,z)

! Create mesh
list = 0
call trmesh(nc0,x,y,z,list,lptr,lend,lnew,near,next,dist,info)

if (.not.allocated(sdata%grid_nnb)) then
   ! Allocation
   allocate(sdata%grid_nnb(sdata%nc0))

   ! Count neighbors
   sdata%grid_nnb = 0
   do ic0=1,nc0
      i = lend(ic0)
      init = .true.
      do while ((i/=lend(ic0)).or.init)
         sdata%grid_nnb(order(ic0)) = sdata%grid_nnb(order(ic0))+1
         i = lptr(i)
         init = .false.
      end do
   end do

   ! Allocation
   allocate(sdata%grid_inb(maxval(sdata%grid_nnb),sdata%nc0))

   ! Find neighbors
   sdata%grid_nnb = 0
   do ic0=1,nc0
      i = lend(ic0)
      init = .true.
      do while ((i/=lend(ic0)).or.init)
         sdata%grid_nnb(order(ic0)) = sdata%grid_nnb(order(ic0))+1
         sdata%grid_inb(sdata%grid_nnb(order(ic0)),order(ic0)) = order(abs(list(i)))
         i = lptr(i)
         init = .false.
      end do
   end do

   ! Copy neighbors for redudant points
   do ic0=1,sdata%nc0
      if (isnotmsi(redundant(ic0))) then
         sdata%grid_nnb(ic0) = sdata%grid_nnb(redundant(ic0))
         sdata%grid_inb(:,ic0) = sdata%grid_inb(:,redundant(ic0))
      end if
   end do
end if

! Compute distances
allocate(sdata%grid_dnb(maxval(sdata%grid_nnb),sdata%nc0))
do ic0=1,sdata%nc0
   do i=1,sdata%grid_nnb(ic0)
      call sphere_dist(sdata%lon(ic0),sdata%lat(ic0),sdata%lon(sdata%grid_inb(i,ic0)), &
    & sdata%lat(sdata%grid_inb(i,ic0)),sdata%grid_dnb(i,ic0))
      sdata%grid_dnb(i,ic0) = (sdata%grid_dnb(i,ic0)/req)**2
   end do
end do

! Allocation
allocate(sdata%nbnd(sdata%nl0))
allocate(ic0_bnd(2,nc0,sdata%nl0))

! Find border points
do il0=1,sdata%nl0
   sdata%nbnd(il0) = 0
   do ic0=1,nc0
      ! Check mask points only
      if (.not.sdata%mask(order(ic0),il0)) then
         i = lend(ic0)
         init = .true.
         do while ((i/=lend(ic0)).or.init)
            jc0 = abs(list(i))
            kc0 = abs(list(lptr(i)))
            if (.not.sdata%mask(order(jc0),il0).and.sdata%mask(order(kc0),il0)) then
               ! Create a new boundary arc
               sdata%nbnd(il0) = sdata%nbnd(il0)+1
               if (sdata%nbnd(il0)>nc0) call msgerror('too many boundary arcs')
               ic0_bnd(1,sdata%nbnd(il0),il0) = order(ic0)
               ic0_bnd(2,sdata%nbnd(il0),il0) = order(jc0)
            end if
            i = lptr(i)
            init = .false.
         end do
      end if
   end do
end do

! Allocation
allocate(sdata%xbnd(2,maxval(sdata%nbnd),sdata%nl0))
allocate(sdata%ybnd(2,maxval(sdata%nbnd),sdata%nl0))
allocate(sdata%zbnd(2,maxval(sdata%nbnd),sdata%nl0))
allocate(sdata%vbnd(3,maxval(sdata%nbnd),sdata%nl0))

do il0=1,sdata%nl0
   ! Compute boundary arcs
   do ibnd=1,sdata%nbnd(il0)
      latbnd = sdata%lat(ic0_bnd(:,ibnd,il0))
      lonbnd = sdata%lon(ic0_bnd(:,ibnd,il0))
      call trans(2,latbnd,lonbnd,sdata%xbnd(:,ibnd,il0),sdata%ybnd(:,ibnd,il0),sdata%zbnd(:,ibnd,il0))
   end do
   do ibnd=1,sdata%nbnd(il0)
      v1 = (/sdata%xbnd(1,ibnd,il0),sdata%ybnd(1,ibnd,il0),sdata%zbnd(1,ibnd,il0)/)
      v2 = (/sdata%xbnd(2,ibnd,il0),sdata%ybnd(2,ibnd,il0),sdata%zbnd(2,ibnd,il0)/)
      call vector_product(v1,v2,sdata%vbnd(:,ibnd,il0))
   end do
end do

end subroutine compute_grid_mesh

!----------------------------------------------------------------------
! Subroutine: compute_subgrid_mesh
!> Purpose: compute subgrid mesh
!----------------------------------------------------------------------
subroutine compute_subgrid_mesh(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: lnew,info,nt,ia,il1,it,i,i1,i2,ic1,namax
integer,allocatable :: list(:),lptr(:),lend(:),near(:),next(:),ltri(:,:)
real(kind_real),allocatable :: x(:),y(:),z(:),dist(:)
logical :: init

! Allocation
allocate(list(6*(sdata%nc1-2)))
allocate(lptr(6*(sdata%nc1-2)))
allocate(lend(sdata%nc1))
allocate(near(sdata%nc1))
allocate(next(sdata%nc1))
allocate(x(sdata%nc1))
allocate(y(sdata%nc1))
allocate(z(sdata%nc1))
allocate(dist(sdata%nc1))
allocate(ltri(9,2*(sdata%nc1-2)))

! Transform to cartesian coordinates
call trans(sdata%nc1,sdata%lat(sdata%ic1_to_ic0),sdata%lon(sdata%ic1_to_ic0),x,y,z)

! Create mesh
list = 0
call trmesh(sdata%nc1,x,y,z,list,lptr,lend,lnew,near,next,dist,info)

! Allocation
allocate(sdata%subgrid_nnb(sdata%nc1))

! Count neighbors
sdata%subgrid_nnb = 0
do ic1=1,sdata%nc1
   i = lend(ic1)
   init = .true.
   do while ((i/=lend(ic1)).or.init)
      sdata%subgrid_nnb(ic1) = sdata%subgrid_nnb(ic1)+1
      i = lptr(i)
      init = .false.
   end do
end do

! Allocation
allocate(sdata%subgrid_inb(maxval(sdata%subgrid_nnb),sdata%nc1))
allocate(sdata%subgrid_dnb(maxval(sdata%subgrid_nnb),sdata%nc1))

! Find neighbors
sdata%subgrid_nnb = 0
do ic1=1,sdata%nc1
   i = lend(ic1)
   init = .true.
   do while ((i/=lend(ic1)).or.init)
      sdata%subgrid_nnb(ic1) = sdata%subgrid_nnb(ic1)+1
      sdata%subgrid_inb(sdata%subgrid_nnb(ic1),ic1) = abs(list(i))
      i = lptr(i)
      init = .false.
   end do
end do

! Compute distances
do ic1=1,sdata%nc1
   do i=1,sdata%subgrid_nnb(ic1)
      call sphere_dist(sdata%lon(ic1),sdata%lat(ic1),sdata%lon(sdata%subgrid_inb(i,ic1)), &
    & sdata%lat(sdata%subgrid_inb(i,ic1)),sdata%subgrid_dnb(i,ic1))
      sdata%subgrid_dnb(i,ic1) = (sdata%subgrid_dnb(i,ic1)/req)**2
   end do
end do

! Create triangles list
call trlist(sdata%nc1,list,lptr,lend,9,nt,ltri,info)
namax = maxval(ltri(7:9,1:nt))

! Release memory
deallocate(list)
deallocate(lptr)
deallocate(lend)
deallocate(near)
deallocate(next)
deallocate(x)
deallocate(y)
deallocate(z)
deallocate(dist)
deallocate(ltri)

! Initialization
allocate(sdata%na(sdata%nl1))
allocate(sdata%larc(2,namax,sdata%nl1))
call msi(sdata%larc)

do il1=1,sdata%nl1
   ! Allocation
   allocate(list(6*(sdata%nc2(il1)-2)))
   allocate(lptr(6*(sdata%nc2(il1)-2)))
   allocate(lend(sdata%nc2(il1)))
   allocate(near(sdata%nc2(il1)))
   allocate(next(sdata%nc2(il1)))
   allocate(x(sdata%nc2(il1)))
   allocate(y(sdata%nc2(il1)))
   allocate(z(sdata%nc2(il1)))
   allocate(dist(sdata%nc2(il1)))
   allocate(ltri(9,2*(sdata%nc2(il1)-2)))

   ! Transform to cartesian coordinates
   call trans(sdata%nc2(il1),sdata%lat(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1)), &
 & sdata%lon(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1)),x,y,z)

   ! Create mesh
   list = 0
   call trmesh(sdata%nc2(il1),x,y,z,list,lptr,lend,lnew,near,next,dist,info)

   ! Create triangles list
   call trlist(sdata%nc2(il1),list,lptr,lend,9,nt,ltri,info)

   ! Copy arcs list
   sdata%na(il1) = maxval(ltri(7:9,1:nt))
   do ia=1,sdata%na(il1)
      it = 1
      do while (it<=nt)
         if (any(ltri(7:9,it)==ia)) exit
         it = it+1
      end do
      i = 1
      do while (i<=3)
         if (ltri(6+i,it)==ia) exit
         i = i+1
      end do
      i1 = mod(i+1,3)
      if (i1==0) i1 = 3
      i2 = mod(i+2,3)
      if (i2==0) i2 = 3
      sdata%larc(1,ia,il1) = ltri(i1,it)
      sdata%larc(2,ia,il1) = ltri(i2,it)
   end do

   ! Release memory
   deallocate(list)
   deallocate(lptr)
   deallocate(lend)
   deallocate(near)
   deallocate(next)
   deallocate(x)
   deallocate(y)
   deallocate(z)
   deallocate(dist)
   deallocate(ltri)
end do

end subroutine compute_subgrid_mesh

end module module_parameters
