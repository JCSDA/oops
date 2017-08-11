!----------------------------------------------------------------------
! Module: module_parameters_convol.f90
!> Purpose: compute NICAS parameters (convolution)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_parameters_convol

use module_namelist, only: nam
use netcdf
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: prog_init,prog_print,msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_reorder
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_randgen, only: initialize_sampling,rand_integer
use type_sdata, only: sdatatype
implicit none

real(kind_real),parameter :: S_inf = 1.0e-2  !< Minimum value for the convolution coefficients
real(kind_real) :: deform = 0.0              !< Deformation coefficient (maximum absolute value: -0.318)

private
public :: compute_convol_network,compute_convol_distance

contains

!----------------------------------------------------------------------
! Subroutine: compute_convol_network
!> Purpose: compute convolution with a network approach
!----------------------------------------------------------------------
subroutine compute_convol_network(sdata,rh0,rv0)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata                 !< Sampling data
real(kind_real),intent(in) :: rh0(sdata%nc0,sdata%nl0) !< Scaled horizontal support radius
real(kind_real),intent(in) :: rv0(sdata%nc0,sdata%nl0) !< Scaled vertical support radius

! Local variables
integer :: n_s_max,progint,ithread,is,ic1,il1,ic0,il0,np,np_new,ip,jc0,jl0,kc0,kl0,jp,i,js
integer :: iproc,is_s(mpl%nproc),is_e(mpl%nproc),ns_loc(mpl%nproc),is_loc
integer :: c_n_s(mpl%nthread)
integer,allocatable :: plist(:,:),plist_new(:,:)
real(kind_real) :: distnorm,disttest,S_test
real(kind_real),allocatable :: dist(:,:)
logical :: add_to_front
logical,allocatable :: done(:),valid(:,:)
type(linoptype) :: c(mpl%nthread),ctmp(mpl%nthread)

! MPI splitting
do iproc=1,mpl%nproc
   is_s(iproc) = (iproc-1)*(sdata%ns/mpl%nproc+1)+1
   is_e(iproc) = min(iproc*(sdata%ns/mpl%nproc+1),sdata%ns)
   ns_loc(iproc) = is_e(iproc)-is_s(iproc)+1
end do

! Allocation
n_s_max = 100*nint(float(sdata%nc0*sdata%nl0)/float(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   call linop_alloc(c(ithread))
end do
allocate(done(ns_loc(mpl%myproc)))

! Compute weights
write(mpl%unit,'(a10,a)',advance='no') '','Compute weights: '
call prog_init(progint,done)
c_n_s = 0
!$omp parallel do private(is_loc,is,ithread,ic1,il1,ic0,il0,plist,plist_new,dist,valid,np,np_new), &
!$omp&            private(ip,jc0,jl0,i,kc0,kl0,distnorm,disttest,add_to_front,jp,js)
do is_loc=1,ns_loc(mpl%myproc)
   ! Indices
   is = is_s(mpl%myproc)+is_loc-1
   ithread = omp_get_thread_num()+1
   ic1 = sdata%is_to_ic1(is)
   il1 = sdata%is_to_il1(is)
   ic0 = sdata%ic1_to_ic0(ic1)
   il0 = sdata%il1_to_il0(il1)

   ! Allocation
   allocate(plist(sdata%nc0*sdata%nl0,2))
   allocate(plist_new(sdata%nc0*sdata%nl0,2))
   allocate(dist(sdata%nc0,sdata%nl0))
   allocate(valid(sdata%nc0,sdata%nl0))

   ! Initialize the front
   np = 1
   plist(1,1) = ic0
   plist(1,2) = il0
   dist = 2.0
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
         do i=1,sdata%grid_nnb(jc0)
            kc0 = sdata%grid_inb(i,jc0)
            do kl0=max(jl0-1,1),min(jl0+1,sdata%nl0)
               if (sdata%mask(kc0,kl0)) then
                  distnorm = sqrt(sdata%grid_dnb(i,jc0)/(0.5*(rh0(jc0,jl0)**2+rh0(kc0,kl0)**2)) &
                           & +abs(sdata%vunit(jl0)-sdata%vunit(kl0))/(0.5*(rv0(jc0,jl0)**2+rv0(kc0,kl0)**2)))
                  disttest = dist(jc0,jl0)+distnorm
                  if (disttest<1.0) then
                     ! Point is inside the support
                     if (disttest<dist(kc0,kl0)) then
                        ! Update distance
                        dist(kc0,kl0) = disttest
                        valid(kc0,kl0) = isnotmsi(sdata%ic0il0_to_is(kc0,kl0))

                        ! Check if the point should be added to the front (avoid duplicates)
                        add_to_front = .true.
                        do jp=1,np_new
                           if ((plist_new(jp,1)==kc0).and.(plist_new(jp,1)==kl0)) then
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
               end if
            end do
         end do
      end do

      ! Copy new front
      np = np_new
      plist(1:np,:) = plist_new(1:np,:)
   end do

   ! Count convolution operations
   do il0=1,sdata%nl0
      do ic0=1,sdata%nc0
         if (valid(ic0,il0)) then
            js = sdata%ic0il0_to_is(ic0,il0)

            if (is>js) then
               ! Distance deformation
               distnorm = dist(ic0,il0)+deform*sin(pi*dist(ic0,il0))

               ! Distance check bound
               if (.not.(distnorm>0.0)) call msgerror('negative normalized distance')

               ! Square-root
               if (nam%lsqrt) distnorm = distnorm*sqrt(2.0)

               ! Gaspari-Cohn (1999) function
               if (distnorm<tiny(1.0)) then
                  S_test = 1.0
               elseif (distnorm<0.5) then
                  S_test = 1.0-8.0*distnorm**5+8.0*distnorm**4+5.0*distnorm**3-20.0/3.0*distnorm**2
               else if (distnorm<1.0) then
                  S_test = 8.0/3.0*distnorm**5-8.0*distnorm**4+5.0*distnorm**3 &
                         & +20.0/3.0*distnorm**2-10.0*distnorm+4.0-1.0/(3.0*distnorm)
               else
                  S_test = 0.0
               end if

               ! Check convolution value
               if (S_test>S_inf) then
                  c_n_s(ithread) = c_n_s(ithread)+1
                  if (c_n_s(ithread)>c(ithread)%n_s) then
                     ! Allocation
                     ctmp(ithread)%n_s = c(ithread)%n_s
                     call linop_alloc(ctmp(ithread))

                     ! Copy data
                     ctmp(ithread)%row = c(ithread)%row
                     ctmp(ithread)%col = c(ithread)%col
                     ctmp(ithread)%S = c(ithread)%S

                     ! Reallocate larger linear operation
                     call linop_dealloc(c(ithread))
                     c(ithread)%n_s = 2*ctmp(ithread)%n_s
                     call linop_alloc(c(ithread))

                     ! Copy data
                     c(ithread)%row(1:ctmp(ithread)%n_s) = ctmp(ithread)%row
                     c(ithread)%col(1:ctmp(ithread)%n_s) = ctmp(ithread)%col
                     c(ithread)%S(1:ctmp(ithread)%n_s) = ctmp(ithread)%S

                     ! Release memory
                     call linop_dealloc(ctmp(ithread))
                  end if

                  ! New operation
                  c(ithread)%row(c_n_s(ithread)) = is
                  c(ithread)%col(c_n_s(ithread)) = js
                  c(ithread)%S(c_n_s(ithread)) = S_test
               end if
            end if
         end if
      end do
   end do

   ! Print progression
   done(is_loc) = .true.
   call prog_print(progint,done)

   ! Release memory
   deallocate(plist)
   deallocate(plist_new)
   deallocate(dist)
   deallocate(valid)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Gather data
sdata%c%prefix = 'c'
sdata%c%n_src = sdata%ns
sdata%c%n_dst = sdata%ns
call convol_gather_data(c_n_s,c,sdata%c)

! Release memory
do ithread=1,mpl%nthread
   call linop_dealloc(c(ithread))
end do

end subroutine compute_convol_network

!----------------------------------------------------------------------
! Subroutine: compute_convol_distance
!> Purpose: compute convolution with a distance approach
!----------------------------------------------------------------------
subroutine compute_convol_distance(sdata,rhs,rvs)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata             !< Sampling data
real(kind_real),intent(in) :: rhs(sdata%ns)        !< Scaled horizontal support radius
real(kind_real),intent(in) :: rvs(sdata%ns)        !< Scaled vertical support radius

! Local variables
integer :: ms,n_s_max,progint,ithread,is,ic1,il1,il0,jc1,jl1,jl0,js,i
integer :: iproc,is_s(mpl%nproc),is_e(mpl%nproc),ns_loc(mpl%nproc),is_loc
integer :: c_n_s(mpl%nthread)
integer,allocatable :: mask_ctree(:),nn_index(:,:)
real(kind_real) :: distnorm,S_test
real(kind_real),allocatable :: nn_dist(:,:)
logical :: submask(sdata%nc1,sdata%nl1)
logical,allocatable :: done(:)
type(ctreetype) :: ctree
type(linoptype) :: c(mpl%nthread),ctmp(mpl%nthread)

! Define submask
submask = .false.
do is=1,sdata%ns
   ic1 = sdata%is_to_ic1(is)
   il1 = sdata%is_to_il1(is)
   submask(ic1,il1) = .true.
end do

! Compute cover tree
write(mpl%unit,'(a10,a)') '','Compute cover tree'
allocate(mask_ctree(sdata%nc1))
mask_ctree = 1
ctree = create_ctree(sdata%nc1,dble(sdata%lon(sdata%ic1_to_ic0)),dble(sdata%lat(sdata%ic1_to_ic0)),mask_ctree)
deallocate(mask_ctree)

! Compute nearest neighbors
write(mpl%unit,'(a10,a)') '','Compute nearest neighbors'
ms = 10*min(floor(pi*nam%resol**2*(1.0-cos(minval(rhs)))/(sqrt(3.0)*minval(rhs)**2)),sdata%nc1)
ms = min(ms,sdata%nc1)
allocate(nn_index(ms,sdata%nc1))
allocate(nn_dist(ms,sdata%nc1))
do ic1=1,sdata%nc1
   call find_nearest_neighbors(ctree,dble(sdata%lon(sdata%ic1_to_ic0(ic1))), &
 & dble(sdata%lat(sdata%ic1_to_ic0(ic1))),ms,nn_index(:,ic1),nn_dist(:,ic1))
end do

! MPI splitting
do iproc=1,mpl%nproc
   is_s(iproc) = (iproc-1)*(sdata%ns/mpl%nproc+1)+1
   is_e(iproc) = min(iproc*(sdata%ns/mpl%nproc+1),sdata%ns)
   ns_loc(iproc) = is_e(iproc)-is_s(iproc)+1
end do

! Allocation
n_s_max = 100*nint(float(sdata%nc0*sdata%nl0)/float(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   call linop_alloc(c(ithread))
end do
allocate(done(ns_loc(mpl%myproc)))

! Compute weights
write(mpl%unit,'(a10,a)',advance='no') '','Compute weights: '
call prog_init(progint,done)
c_n_s = 0
!$omp parallel do private(is_loc,is,ithread,ic1,il1,il0,i,jc1,jl1,jl0,js,distnorm,S_test)
do is_loc=1,ns_loc(mpl%myproc)
   ! Indices
   is = is_s(mpl%myproc)+is_loc-1
   ithread = omp_get_thread_num()+1
   ic1 = sdata%is_to_ic1(is)
   il1 = sdata%is_to_il1(is)
   il0 = sdata%il1_to_il0(il1)

   ! Loop on nearest neighbors
   do i=1,ms
      jc1 = nn_index(i,ic1)
      do jl1=1,sdata%nl1
         if (submask(jc1,jl1)) then
            jl0 = sdata%il1_to_il0(jl1)
            js = sdata%ic1il1_to_is(jc1,jl1)
            ! Only half of the (symmetric) matrix coefficients should be stored
            if (is>js) then
               ! Normalized distance
               distnorm = sqrt(nn_dist(i,ic1)**2/(0.5*(rhs(is)**2+rhs(js)**2)) &
                        & +(sdata%vunit(il0)-sdata%vunit(jl0))**2/(0.5*(rvs(is)**2+rvs(js)**2)))

               if (distnorm<1.0) then
                  ! Distance deformation
                  distnorm = distnorm+deform*sin(pi*distnorm)

                  ! Square-root
                  if (nam%lsqrt) distnorm = distnorm*sqrt(2.0)

                  ! Gaspari-Cohn (1999) function
                  if (distnorm<0.5) then
                     S_test = 1.0-8.0*distnorm**5+8.0*distnorm**4+5.0*distnorm**3-20.0/3.0*distnorm**2
                  else if (distnorm<1.0) then
                     S_test = 8.0/3.0*distnorm**5-8.0*distnorm**4+5.0*distnorm**3 &
                            & +20.0/3.0*distnorm**2-10.0*distnorm+4.0-1.0/(3.0*distnorm)
                  else
                     S_test = 0.0
                  end if

                  ! Check convolution value
                  if (S_test>S_inf) then
                     c_n_s(ithread) = c_n_s(ithread)+1
                     if (c_n_s(ithread)>c(ithread)%n_s) then
                        ! Allocation
                        ctmp(ithread)%n_s = c(ithread)%n_s
                        call linop_alloc(ctmp(ithread))

                        ! Copy data
                        ctmp(ithread)%row = c(ithread)%row
                        ctmp(ithread)%col = c(ithread)%col
                        ctmp(ithread)%S = c(ithread)%S

                        ! Reallocate larger linear operation
                        call linop_dealloc(c(ithread))
                        c(ithread)%n_s = 2*ctmp(ithread)%n_s
                        call linop_alloc(c(ithread))

                        ! Copy data
                        c(ithread)%row(1:ctmp(ithread)%n_s) = ctmp(ithread)%row
                        c(ithread)%col(1:ctmp(ithread)%n_s) = ctmp(ithread)%col
                        c(ithread)%S(1:ctmp(ithread)%n_s) = ctmp(ithread)%S

                        ! Release memory
                        call linop_dealloc(ctmp(ithread))
                     end if

                     ! New operation
                     c(ithread)%row(c_n_s(ithread)) = is
                     c(ithread)%col(c_n_s(ithread)) = js
                     c(ithread)%S(c_n_s(ithread)) = S_test
                  end if
               end if
            end if
         end if
      end do
   end do

   ! Print progression
   done(is_loc) = .true.
   call prog_print(progint,done)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Gather data
sdata%c%prefix = 'c'
sdata%c%n_src = sdata%ns
sdata%c%n_dst = sdata%ns
call convol_gather_data(c_n_s,c,sdata%c)

! Release memory
call delete_ctree(ctree)
deallocate(nn_index)
deallocate(nn_dist)
do ithread=1,mpl%nthread
   call linop_dealloc(c(ithread))
end do

end subroutine compute_convol_distance

!----------------------------------------------------------------------
! Subroutine: convol_gather_data
!> Purpose: gather convolution data
!----------------------------------------------------------------------
subroutine convol_gather_data(c_n_s,cin,cout)

implicit none

! Passed variables
integer,intent(in) :: c_n_s(mpl%nthread)       !< Number of operations handled by each thread
type(linoptype),intent(in) :: cin(mpl%nthread) !< Linear operator for each thread
type(linoptype),intent(inout) :: cout          !< Gathered linear operator

! Local variables
integer :: ithread,offset,iproc,i_s
integer :: csum_n_sg(mpl%nproc)
type(linoptype) :: csum,cg(mpl%nproc)

! Allocation
csum%n_s = sum(c_n_s)
call linop_alloc(csum)

! Gather convolution data from OpenMP threads
offset = 0
do ithread=1,mpl%nthread
   csum%row(offset+1:offset+c_n_s(ithread)) = cin(ithread)%row(1:c_n_s(ithread))
   csum%col(offset+1:offset+c_n_s(ithread)) = cin(ithread)%col(1:c_n_s(ithread))
   csum%S(offset+1:offset+c_n_s(ithread)) = cin(ithread)%S(1:c_n_s(ithread))
   offset = offset+c_n_s(ithread)
end do

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         csum_n_sg(iproc) = csum%n_s
      else
         ! Receive data on ioproc
         call mpl_recv(csum_n_sg(iproc),iproc,mpl%tag)
      end if
   end do
else
   ! Send data to ioproc
   call mpl_send(csum%n_s,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl_bcast(csum_n_sg,mpl%ioproc)

! Allocation
cout%n_s = sum(csum_n_sg)
call linop_alloc(cout)

! Communication
if (mpl%main) then
   ! Allocation
   do iproc=1,mpl%nproc
      cg(iproc)%n_s = maxval(csum_n_sg)
      call linop_alloc(cg(iproc))
   end do

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         cg(iproc)%row(1:csum_n_sg(iproc)) = csum%row
         cg(iproc)%col(1:csum_n_sg(iproc)) = csum%col
         cg(iproc)%S(1:csum_n_sg(iproc)) = csum%S
      else
         ! Receive data on ioproc
         call mpl_recv(csum_n_sg(iproc),cg(iproc)%row(1:csum_n_sg(iproc)),iproc,mpl%tag)
         call mpl_recv(csum_n_sg(iproc),cg(iproc)%col(1:csum_n_sg(iproc)),iproc,mpl%tag+1)
         call mpl_recv(csum_n_sg(iproc),cg(iproc)%S(1:csum_n_sg(iproc)),iproc,mpl%tag+2)
      end if
   end do

   ! Format data
   offset = 0
   do iproc=1,mpl%nproc
      do i_s=1,csum_n_sg(iproc)
         cout%row(offset+i_s) = cg(iproc)%row(i_s)
         cout%col(offset+i_s) = cg(iproc)%col(i_s)
         cout%S(offset+i_s) = cg(iproc)%S(i_s)
      end do
      offset = offset+csum_n_sg(iproc)
   end do

   ! Release memory
   do iproc=1,mpl%nproc
      call linop_dealloc(cg(iproc))
   end do
else
   ! Send data to ioproc
   call mpl_send(csum%n_s,csum%row,mpl%ioproc,mpl%tag)
   call mpl_send(csum%n_s,csum%col,mpl%ioproc,mpl%tag+1)
   call mpl_send(csum%n_s,csum%S,mpl%ioproc,mpl%tag+2)
end if
mpl%tag = mpl%tag+3

! Broadcast data
call mpl_bcast(cout%row,mpl%ioproc)
call mpl_bcast(cout%col,mpl%ioproc)
call mpl_bcast(cout%S,mpl%ioproc)

! Reorder linear operator
call linop_reorder(cout)

end subroutine convol_gather_data

end module module_parameters_convol
