!----------------------------------------------------------------------
! Module: nicas_parameters_convol.f90
!> Purpose: compute NICAS parameters (convolution)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_parameters_convol

use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product,gc99
use tools_display, only: prog_init,prog_print,msgerror
use tools_interp, only: check_arc
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use type_bdata, only: bdatatype
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy
use type_mesh, only: meshtype,create_mesh
use type_mpl, only: mpl
use type_nam, only: namtype
use type_ndata, only: ndatatype

implicit none

real(kind_real),parameter :: sqrt_fac = 0.721 !< Square-root factor (empirical)
real(kind_real),parameter :: S_inf = 1.0e-2   !< Minimum value for the convolution coefficients
real(kind_real),parameter :: deform = 0.0     !< Deformation coefficient (maximum absolute value: -0.318)

private
public :: compute_convol_network,compute_convol_distance

contains

!----------------------------------------------------------------------
! Subroutine: compute_convol_network
!> Purpose: compute convolution with a network approach
!----------------------------------------------------------------------
subroutine compute_convol_network(bdata,ndata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata    !< B data
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: n_s_max,progint,ithread,is,ic1,ic0,jl0,jl1,np,np_new,i,j,k,ip,jc0,jc1,il0,djl0,il1,kc0,dkl0,kl0,jp,js,isb,offset
integer :: c_n_s(mpl%nthread),c_nor_n_s(mpl%nthread)
integer,allocatable :: net_nnb(:),net_inb(:,:),plist(:,:),plist_new(:,:)
real(kind_real) :: dnb,disthsq,distvsq,rhsq,rvsq,distnorm,disttest,S_test
real(kind_real),allocatable :: net_dnb(:,:,:,:),dist(:,:)
logical :: init,add_to_front,valid_arc
logical,allocatable :: done(:),valid(:,:)
type(linoptype) :: c(mpl%nthread),c_nor(mpl%nthread)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Allocation
allocate(net_nnb(geom%nc0))

! Count neighbors
net_nnb = 0
do ic0=1,geom%nc0
   i = geom%mesh%lend(ic0)
   init = .true.
   do while ((i/=geom%mesh%lend(ic0)).or.init)
      jc0 = geom%mesh%order(ic0)
      net_nnb(jc0) = net_nnb(jc0)+1
      i = geom%mesh%lptr(i)
      init = .false.
   end do
end do

! Allocation
allocate(net_inb(maxval(net_nnb),geom%nc0))
allocate(net_dnb(maxval(net_nnb),-1:1,geom%nc0,geom%nl0))

! Find neighbors
net_nnb = 0
do ic0=1,geom%nc0
   i = geom%mesh%lend(ic0)
   init = .true.
   do while ((i/=geom%mesh%lend(ic0)).or.init)
      jc0 = geom%mesh%order(ic0)
      net_nnb(jc0) = net_nnb(jc0)+1
      net_inb(net_nnb(jc0),jc0) = geom%mesh%order(abs(geom%mesh%list(i)))
      i = geom%mesh%lptr(i)
      init = .false.
   end do
end do

! Compute distances
net_dnb = huge(1.0)
!$omp parallel do schedule(static) private(ic0,j,jc0,dnb,il0,valid_arc,djl0,jl0,disthsq,distvsq,rhsq,rvsq,distnorm)
do ic0=1,geom%nc0
   do j=1,net_nnb(ic0)
      ! Index
      jc0 = net_inb(j,ic0)

      if (any(geom%mask(jc0,:))) then
         ! True distance
         call sphere_dist(geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),dnb)

         do il0=1,geom%nl0
            ! Check mask bounds
            valid_arc = .true.
            if (nam%mask_check) call check_arc(geom,il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),valid_arc)

            if (valid_arc) then
               do djl0=-1,1
                  ! Index
                  jl0 = max(1,min(il0+djl0,geom%nl0))

                  if (geom%mask(ic0,il0).and.geom%mask(jc0,jl0)) then
                     ! Squared support radii
                     disthsq = dnb**2
                     distvsq = geom%distv(il0,jl0)**2
                     rhsq = 0.5*(bdata%rh0(ic0,il0)**2+bdata%rh0(jc0,jl0)**2)
                     rvsq = 0.5*(bdata%rv0(ic0,il0)**2+bdata%rv0(jc0,jl0)**2)
                     distnorm = 0.0
                     if (rhsq>0.0) then
                        distnorm = distnorm+disthsq/rhsq
                     elseif (disthsq>0.0) then
                        distnorm = distnorm+0.5*huge(0.0)
                     end if
                     if (rvsq>0.0) then
                        distnorm = distnorm+distvsq/rvsq
                     elseif (distvsq>0.0) then
                        distnorm = distnorm+0.5*huge(0.0)
                     end if
                     net_dnb(j,djl0,ic0,il0) = sqrt(distnorm)
                  end if
               end do
            end if
         end do
      end if
   end do
end do
!$omp end parallel do

! Allocation
n_s_max = 100*nint(float(geom%nc0*geom%nl0)/float(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   c_nor(ithread)%n_s = n_s_max
   call linop_alloc(c(ithread))
   call linop_alloc(c_nor(ithread))
end do
allocate(done(ndata%nsb))

! Compute weights
write(mpl%unit,'(a10,a)',advance='no') '','Compute weights: '
call prog_init(progint,done)
c_n_s = 0
c_nor_n_s = 0
!$omp parallel do schedule(static) private(isb,is,ithread,ic1,il1,ic0,il0,plist,plist_new,dist,valid,np,np_new,ip,jc0,jl0), &
!$omp&                             private(k,kc0,kl0,disttest,add_to_front,jp,js,jc1,jl1,distnorm,S_test)
do isb=1,ndata%nsb
   ! Indices
   is = ndata%sb_to_s(isb)
   ithread = omp_get_thread_num()+1
   ic1 = ndata%s_to_c1(is)
   il1 = ndata%s_to_l1(is)
   ic0 = ndata%c1_to_c0(ic1)
   il0 = ndata%l1_to_l0(il1)

   ! Allocation
   allocate(plist(geom%nc0*geom%nl0,2))
   allocate(plist_new(geom%nc0*geom%nl0,2))
   allocate(dist(geom%nc0,geom%nl0))
   allocate(valid(geom%nc0,geom%nl0))

   ! Initialize the front
   np = 1
   plist(1,1) = ic0
   plist(1,2) = il0
   dist = 1.0
   dist(ic0,il0) = 0.0
   valid = .false.
   valid(ic0,il0) = .true.

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         jc0 = plist(ip,1)
         jl0 = plist(ip,2)

         ! Loop over neighbors
         do k=1,net_nnb(jc0)
            kc0 = net_inb(k,jc0)
            do dkl0=-1,1
               kl0 = max(1,min(jl0+dkl0,geom%nl0))
               disttest = dist(jc0,jl0)+net_dnb(k,dkl0,jc0,jl0)
               if (disttest<1.0) then
                  ! Point is inside the support
                  if (disttest<dist(kc0,kl0)) then
                     ! Update distance
                     dist(kc0,kl0) = disttest
                     valid(kc0,kl0) = .true.

                     ! Check if the point should be added to the front (avoid duplicates)
                     add_to_front = .true.
                     do jp=1,np_new
                        if ((plist_new(jp,1)==kc0).and.(plist_new(jp,2)==kl0)) then
                           add_to_front = .false.
                           exit
                        end if
                     end do

                     if (add_to_front) then
                        ! Add point to the front
                        np_new = np_new+1
                        plist_new(np_new,1) = kc0
                        plist_new(np_new,2) = kl0
                     end if
                  end if
               end if
            end do
         end do
      end do

      ! Copy new front
      np = np_new
      plist(1:np,:) = plist_new(1:np,:)
   end do

   ! Count convolution operations
   do js=1,ndata%ns
      ! Indices
      jc1 = ndata%s_to_c1(js)
      jl1 = ndata%s_to_l1(js)
      jc0 = ndata%c1_to_c0(jc1)
      jl0 = ndata%l1_to_l0(jl1)

      if (valid(jc0,jl0)) then
         ! Normalized distance
         distnorm = dist(jc0,jl0)

         if (distnorm<1.0) then
            ! Distance deformation
            distnorm = distnorm+deform*sin(pi*distnorm)

            ! Square-root
            if (nam%lsqrt) distnorm = distnorm/sqrt_fac

            ! Gaspari-Cohn (1999) function
            S_test = gc99(distnorm)

            ! Store coefficient for convolution
            if (nam%mpicom==1) then
               if ((ndata%lcheck_sb(js).and.(is<js)).or.(.not.ndata%lcheck_sb(js))) &
             & call check_convol(is,js,S_test,c_n_s(ithread),c(ithread))
            elseif (nam%mpicom==2) then
               if (ndata%lcheck_sa(is).and.((ndata%lcheck_sa(js).and.(is<js)).or.(.not.ndata%lcheck_sa(js)))) &
             & call check_convol(is,js,S_test,c_n_s(ithread),c(ithread))
            end if

            ! Store coefficient for normalization
            if (nam%lsqrt) then
               if ((ndata%lcheck_sb(js).and.(is<js)).or.(.not.ndata%lcheck_sb(js))) &
            &  call check_convol(is,js,S_test,c_nor_n_s(ithread),c_nor(ithread))
            else
               if (ndata%lcheck_sb(js).and.(is<js)) call check_convol(is,js,S_test,c_nor_n_s(ithread),c_nor(ithread))
            end if
         end if
      end if
   end do

   ! Print progression
   done(isb) = .true.
   call prog_print(progint,done)

   ! Release memory
   deallocate(plist)
   deallocate(plist_new)
   deallocate(dist)
   deallocate(valid)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Initialize object
ndata%c%prefix = 'c'
ndata%c_nor%prefix = 'c_nor'

! Gather data
ndata%c%n_s = sum(c_n_s)
ndata%c_nor%n_s = sum(c_nor_n_s)
call linop_alloc(ndata%c)
call linop_alloc(ndata%c_nor)

! Gather convolution data from OpenMP threads
offset = 0
do ithread=1,mpl%nthread
   ndata%c%row(offset+1:offset+c_n_s(ithread)) = c(ithread)%row(1:c_n_s(ithread))
   ndata%c%col(offset+1:offset+c_n_s(ithread)) = c(ithread)%col(1:c_n_s(ithread))
   ndata%c%S(offset+1:offset+c_n_s(ithread)) = c(ithread)%S(1:c_n_s(ithread))
   offset = offset+c_n_s(ithread)
end do
offset = 0
do ithread=1,mpl%nthread
   ndata%c_nor%row(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%row(1:c_nor_n_s(ithread))
   ndata%c_nor%col(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%col(1:c_nor_n_s(ithread))
   ndata%c_nor%S(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%S(1:c_nor_n_s(ithread))
   offset = offset+c_nor_n_s(ithread)
end do

! Release memory
do ithread=1,mpl%nthread
   call linop_dealloc(c(ithread))
   call linop_dealloc(c_nor(ithread))
end do

! End associate
end associate

end subroutine compute_convol_network

!----------------------------------------------------------------------
! Subroutine: compute_convol_distance
!> Purpose: compute convolution with a distance approach
!----------------------------------------------------------------------
subroutine compute_convol_distance(bdata,ndata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata    !< B data
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: ms,n_s_max,progint,ithread,is,ic1,jl0,jc1,il1,jl1,il0,j,js,isb,ic1b,offset
integer :: c_n_s(mpl%nthread),c_nor_n_s(mpl%nthread)
integer,allocatable :: nn_index(:,:)
real(kind_real) :: disthsq,distvsq,rhsq,rvsq,distnorm,S_test
real(kind_real),allocatable :: nn_dist(:,:)
logical :: force,test,valid,submask(ndata%nc1,ndata%nl1)
logical,allocatable :: mask_ctree(:),done(:)
type(ctreetype) :: ctree
type(linoptype) :: c(mpl%nthread),c_nor(mpl%nthread)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Define submask
submask = .false.
do is=1,ndata%ns
   ic1 = ndata%s_to_c1(is)
   il1 = ndata%s_to_l1(is)
   submask(ic1,il1) = .true.
end do

! Compute cover tree
write(mpl%unit,'(a10,a)') '','Compute cover tree'
allocate(mask_ctree(ndata%nc1))
mask_ctree = .true.
ctree = create_ctree(ndata%nc1,dble(geom%lon(ndata%c1_to_c0)),dble(geom%lat(ndata%c1_to_c0)),mask_ctree)
deallocate(mask_ctree)

! Theoretical number of neighbors
ms = 0
do il1=1,ndata%nl1
   ms = max(ms,floor(0.25*float(ndata%nc1)*(sum(bdata%rh0(ndata%c1_to_c0,ndata%l1_to_l0(il1)),mask=ndata%c2mask(:,il1)) &
      & /float(ndata%nc2(il1)))**2))
end do

! Find all necessary neighbors
valid = .false.
do while (.not.valid)
   ! Check number of neighbors
   force = ms>ndata%nc1
   ms = max(1,min(ms,ndata%nc1))
   write(mpl%unit,'(a10,a,i8,a,f5.1,a)') '','Try with ',ms,' neighbors (',float(ms)/float(ndata%nc1)*100.0,'%)'

   ! Allocation
   if (allocated(nn_index)) deallocate(nn_index)
   if (allocated(nn_dist)) deallocate(nn_dist)
   allocate(nn_index(ms,ndata%nc1b))
   allocate(nn_dist(ms,ndata%nc1b))

   ! Loop over points
   test = .true.
   ic1b = 1
   do while ((test.or.force).and.(ic1b<=ndata%nc1b))
      ic1 = ndata%c1b_to_c1(ic1b)

      ! Compute nearest neighbors
      call find_nearest_neighbors(ctree,dble(geom%lon(ndata%c1_to_c0(ic1))), &
    & dble(geom%lat(ndata%c1_to_c0(ic1))),ms,nn_index(:,ic1b),nn_dist(:,ic1b))

      ! Loop over levels
      jc1 = nn_index(ms,ic1b)
      il1 = 1
      do while ((test.or.force).and.(il1<=ndata%nl1))
         rhsq = 0.5*(bdata%rh0(ndata%c1_to_c0(ic1),ndata%l1_to_l0(il1))**2+ &
               & bdata%rh0(ndata%c1_to_c0(jc1),ndata%l1_to_l0(il1))**2)
         distnorm = nn_dist(ms,ic1b)**2/rhsq
         if ((distnorm<1.0).and.(.not.force)) then
            ms = 2*ms
            test = .false.
         end if
         il1 = il1+1
      end do
      ic1b = ic1b+1
   end do
   if ((ic1b>ndata%nc1b).or.force) valid = .true.
end do

! Allocation
n_s_max = 100*nint(float(geom%nc0*geom%nl0)/float(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   c_nor(ithread)%n_s = n_s_max
   call linop_alloc(c(ithread))
   call linop_alloc(c_nor(ithread))
end do
allocate(done(ndata%nsb))

! Compute weights
write(mpl%unit,'(a10,a)',advance='no') '','Compute weights: '
call prog_init(progint,done)
c_n_s = 0
c_nor_n_s = 0
!$omp parallel do schedule(static) private(isb,is,ithread,ic1,ic1b,il1,il0,j,jc1,jl1,jl0,js,disthsq,distvsq,rhsq,rvsq), &
!$omp&                             private(distnorm,S_test)
do isb=1,ndata%nsb
   ! Indices
   is = ndata%sb_to_s(isb)
   ithread = omp_get_thread_num()+1
   ic1 = ndata%s_to_c1(is)
   ic1b = ndata%c1_to_c1b(ic1)
   il1 = ndata%s_to_l1(is)
   il0 = ndata%l1_to_l0(il1)

   ! Loop on nearest neighbors
   do j=1,ms
      jc1 = nn_index(j,ic1b)
      do jl1=1,ndata%nl1
         if (submask(jc1,jl1)) then
            jl0 = ndata%l1_to_l0(jl1)
            js = ndata%c1l1_to_s(jc1,jl1)

            ! Normalized distance
            disthsq = nn_dist(j,ic1b)**2
            distvsq = geom%distv(il0,jl0)**2
            rhsq = 0.5*(bdata%rh0(ndata%c1_to_c0(ndata%s_to_c1(is)),ndata%l1_to_l0(ndata%s_to_l1(is)))**2+ &
                 & bdata%rh0(ndata%c1_to_c0(ndata%s_to_c1(js)),ndata%l1_to_l0(ndata%s_to_l1(js)))**2)
            rvsq = 0.5*(bdata%rv0(ndata%c1_to_c0(ndata%s_to_c1(is)),ndata%l1_to_l0(ndata%s_to_l1(is)))**2+ &
                 & bdata%rv0(ndata%c1_to_c0(ndata%s_to_c1(js)),ndata%l1_to_l0(ndata%s_to_l1(js)))**2)
            distnorm = 0.0
            if (rhsq>0.0) then
               distnorm = distnorm+disthsq/rhsq
            elseif (disthsq>0.0) then
               distnorm = distnorm+0.5*huge(0.0)
            end if
            if (rvsq>0.0) then
               distnorm = distnorm+distvsq/rvsq
            elseif (distvsq>0.0) then
               distnorm = distnorm+0.5*huge(0.0)
            end if
            distnorm = sqrt(distnorm)

            if (distnorm<1.0) then
               ! Distance deformation
               distnorm = distnorm+deform*sin(pi*distnorm)

               ! Square-root
               if (nam%lsqrt) distnorm = distnorm/sqrt_fac

               ! Gaspari-Cohn (1999) function
               S_test = gc99(distnorm)

               ! Store coefficient for convolution
               if (nam%mpicom==1) then
                  if ((ndata%lcheck_sb(js).and.(is<js)).or.(.not.ndata%lcheck_sb(js))) &
                & call check_convol(is,js,S_test,c_n_s(ithread),c(ithread))
               elseif (nam%mpicom==2) then
                  if (ndata%lcheck_sa(is).and.((ndata%lcheck_sa(js).and.(is<js)).or.(.not.ndata%lcheck_sa(js)))) &
                & call check_convol(is,js,S_test,c_n_s(ithread),c(ithread))
               end if

               ! Store coefficient for normalization
               if (nam%lsqrt) then
                  if ((ndata%lcheck_sb(js).and.(is<js)).or.(.not.ndata%lcheck_sb(js))) &
               &  call check_convol(is,js,S_test,c_nor_n_s(ithread),c_nor(ithread))
               else
                  if (ndata%lcheck_sb(js).and.(is<js)) call check_convol(is,js,S_test,c_nor_n_s(ithread),c_nor(ithread))
               end if
            end if
         end if
      end do
   end do

   ! Print progression
   done(isb) = .true.
   call prog_print(progint,done)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Initialize object
ndata%c%prefix = 'c'
ndata%c_nor%prefix = 'c_nor'

! Gather data
ndata%c%n_s = sum(c_n_s)
ndata%c_nor%n_s = sum(c_nor_n_s)
call linop_alloc(ndata%c)
call linop_alloc(ndata%c_nor)

! Gather convolution data from OpenMP threads
offset = 0
do ithread=1,mpl%nthread
   ndata%c%row(offset+1:offset+c_n_s(ithread)) = c(ithread)%row(1:c_n_s(ithread))
   ndata%c%col(offset+1:offset+c_n_s(ithread)) = c(ithread)%col(1:c_n_s(ithread))
   ndata%c%S(offset+1:offset+c_n_s(ithread)) = c(ithread)%S(1:c_n_s(ithread))
   offset = offset+c_n_s(ithread)
end do
offset = 0
do ithread=1,mpl%nthread
   ndata%c_nor%row(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%row(1:c_nor_n_s(ithread))
   ndata%c_nor%col(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%col(1:c_nor_n_s(ithread))
   ndata%c_nor%S(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%S(1:c_nor_n_s(ithread))
   offset = offset+c_nor_n_s(ithread)
end do

! Release memory
do ithread=1,mpl%nthread
   call linop_dealloc(c(ithread))
   call linop_dealloc(c_nor(ithread))
end do

! End associate
end associate

end subroutine compute_convol_distance

!----------------------------------------------------------------------
! Subroutine: check_convol
!> Purpose: check convolution value and add it if necessary
!----------------------------------------------------------------------
subroutine check_convol(is,js,S_test,c_n_s,c)

implicit none

! Passed variables
integer,intent(in) :: is             !< First point index
integer,intent(in) :: js             !< Second point index
real(kind_real),intent(in) :: S_test !< Test interpolation convolution value
integer,intent(inout) :: c_n_s       !< Number of convolution operations
type(linoptype),intent(inout) :: c   !< Convolution data

! Local variables
type(linoptype) :: ctmp

if (S_test>S_inf) then
   c_n_s = c_n_s+1
   if (c_n_s>c%n_s) then
      ! Copy
      call linop_copy(c,ctmp)

      ! Reallocate larger linear operation
      call linop_dealloc(c)
      c%n_s = 2*ctmp%n_s
      call linop_alloc(c)

      ! Copy data
      c%row(1:ctmp%n_s) = ctmp%row
      c%col(1:ctmp%n_s) = ctmp%col
      c%S(1:ctmp%n_s) = ctmp%S

      ! Release memory
      call linop_dealloc(ctmp)
   end if

   ! New operation
   c%row(c_n_s) = is
   c%col(c_n_s) = js
   c%S(c_n_s) = S_test
end if

end subroutine check_convol

end module nicas_parameters_convol
