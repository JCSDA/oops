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

use module_namelist, only: namtype
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: prog_init,prog_print,msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_reorder
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_ndata, only: ndatatype
use type_randgen, only: initialize_sampling,rand_integer

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
subroutine compute_convol_network(nam,ndata,rh0,rv0)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(inout) :: ndata                 !< Sampling data
real(kind_real),intent(in) :: rh0(ndata%nc0,ndata%nl0) !< Scaled horizontal support radius
real(kind_real),intent(in) :: rv0(ndata%nc0,ndata%nl0) !< Scaled vertical support radius

! Local variables
integer :: n_s_max,progint,ithread,is,ic1,il1,ic0,il0,np,np_new,ip,jc0,jl0,kc0,kl0,jp,i,js
integer :: iproc,is_s(mpl%nproc),is_e(mpl%nproc),ns_loc(mpl%nproc),is_loc
integer :: c_n_s(mpl%nthread)
integer,allocatable :: plist(:,:),plist_new(:,:)
real(kind_real) :: distnorm,disttest,S_test
real(kind_real),allocatable :: dist(:,:)
logical :: add_to_front
logical,allocatable :: done(:),valid(:,:)
type(linoptype) :: c(mpl%nthread)

! MPI splitting
do iproc=1,mpl%nproc
   is_s(iproc) = (iproc-1)*(ndata%ns/mpl%nproc+1)+1
   is_e(iproc) = min(iproc*(ndata%ns/mpl%nproc+1),ndata%ns)
   ns_loc(iproc) = is_e(iproc)-is_s(iproc)+1
end do

! Allocation
n_s_max = 100*nint(float(ndata%nc0*ndata%nl0)/float(mpl%nthread*mpl%nproc))
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
   ic1 = ndata%is_to_ic1(is)
   il1 = ndata%is_to_il1(is)
   ic0 = ndata%ic1_to_ic0(ic1)
   il0 = ndata%il1_to_il0(il1)

   ! Allocation
   allocate(plist(ndata%nc0*ndata%nl0,2))
   allocate(plist_new(ndata%nc0*ndata%nl0,2))
   allocate(dist(ndata%nc0,ndata%nl0))
   allocate(valid(ndata%nc0,ndata%nl0))

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
         do i=1,ndata%geom%net_nnb(jc0)
            kc0 = ndata%geom%net_inb(i,jc0)
            do kl0=max(jl0-1,1),min(jl0+1,ndata%nl0)
               if (ndata%geom%mask(kc0,kl0)) then
                  distnorm = sqrt(ndata%geom%net_dnb(i,jc0)/(0.5*(rh0(jc0,jl0)**2+rh0(kc0,kl0)**2)) &
                           & +abs(ndata%geom%vunit(jl0)-ndata%geom%vunit(kl0))/(0.5*(rv0(jc0,jl0)**2+rv0(kc0,kl0)**2)))
                  disttest = dist(jc0,jl0)+distnorm
                  if (disttest<1.0) then
                     ! Point is inside the support
                     if (disttest<dist(kc0,kl0)) then
                        ! Update distance
                        dist(kc0,kl0) = disttest
                        valid(kc0,kl0) = isnotmsi(ndata%ic0il0_to_is(kc0,kl0))

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
   do il0=1,ndata%nl0
      do ic0=1,ndata%nc0
         if (valid(ic0,il0)) then
            js = ndata%ic0il0_to_is(ic0,il0)

            ! Only half of the (symmetric) matrix coefficients should be stored
            if (is>js) then
               ! Normalized distance
               distnorm = dist(ic0,il0)

               if (distnorm<1.0) then
                  ! Gaspari-Cohn (1999) function
                  S_test = gc99(nam,distnorm)

                  ! Check convolution value
                  call check_convol(is,js,S_test,c_n_s(ithread),c(ithread))
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

! Initialize object
ndata%c%prefix = 'c'
ndata%c%n_src = ndata%ns
ndata%c%n_dst = ndata%ns

! Gather data
call convol_gather_data(c_n_s,c,ndata%c)

! Release memory
do ithread=1,mpl%nthread
   call linop_dealloc(c(ithread))
end do

end subroutine compute_convol_network

!----------------------------------------------------------------------
! Subroutine: compute_convol_distance
!> Purpose: compute convolution with a distance approach
!----------------------------------------------------------------------
subroutine compute_convol_distance(nam,ndata,rhs,rvs)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(inout) :: ndata             !< Sampling data
real(kind_real),intent(in) :: rhs(ndata%ns)        !< Scaled horizontal support radius
real(kind_real),intent(in) :: rvs(ndata%ns)        !< Scaled vertical support radius

! Local variables
integer :: ms,n_s_max,progint,ithread,is,ic1,il1,il0,jc1,jl1,jl0,js,i,iproc
integer :: is_s(mpl%nproc),is_e(mpl%nproc),ns_loc(mpl%nproc),is_loc
integer :: ic1_s(mpl%nproc),ic1_e(mpl%nproc),nc1_loc(mpl%nproc),ic1_loc
integer :: c_n_s(mpl%nthread)
integer,allocatable :: mask_ctree(:),nn_index(:,:)
integer,allocatable :: rbuf_index(:),sbuf_index(:)
real(kind_real) :: distnorm,S_test
real(kind_real),allocatable :: nn_dist(:,:)
real(kind_real),allocatable :: rbuf_dist(:),sbuf_dist(:)
logical :: submask(ndata%nc1,ndata%nl1)
logical,allocatable :: done(:)
type(ctreetype) :: ctree
type(linoptype) :: c(mpl%nthread)

! Define submask
submask = .false.
do is=1,ndata%ns
   ic1 = ndata%is_to_ic1(is)
   il1 = ndata%is_to_il1(is)
   submask(ic1,il1) = .true.
end do

! Compute cover tree
write(mpl%unit,'(a10,a)') '','Compute cover tree'
allocate(mask_ctree(ndata%nc1))
mask_ctree = 1
ctree = create_ctree(ndata%nc1,dble(ndata%geom%lon(ndata%ic1_to_ic0)),dble(ndata%geom%lat(ndata%ic1_to_ic0)),mask_ctree)
deallocate(mask_ctree)

! Number of neighbors
ms = 10*min(floor(pi*nam%resol**2*(1.0-cos(minval(rhs)))/(sqrt(3.0)*minval(rhs)**2)),ndata%nc1)
ms = min(ms,ndata%nc1)

! MPI splitting
do iproc=1,mpl%nproc
   ic1_s(iproc) = (iproc-1)*(ndata%nc1/mpl%nproc+1)+1
   ic1_e(iproc) = min(iproc*(ndata%nc1/mpl%nproc+1),ndata%nc1)
   nc1_loc(iproc) = ic1_e(iproc)-ic1_s(iproc)+1
end do

! Allocation
allocate(nn_index(ms,ndata%nc1))
allocate(nn_dist(ms,ndata%nc1))

! Compute nearest neighbors
write(mpl%unit,'(a10,a)') '','Compute nearest neighbors'
do ic1_loc=1,nc1_loc(mpl%myproc)
   ic1 = ic1_s(mpl%myproc)+ic1_loc-1
   call find_nearest_neighbors(ctree,dble(ndata%geom%lon(ndata%ic1_to_ic0(ic1))), &
 & dble(ndata%geom%lat(ndata%ic1_to_ic0(ic1))),ms,nn_index(:,ic1),nn_dist(:,ic1))
end do

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc/=mpl%ioproc) then
         ! Allocation
         allocate(rbuf_index(ms*nc1_loc(iproc)))
         allocate(rbuf_dist(ms*nc1_loc(iproc)))

         ! Receive data on ioproc
         call mpl_recv(ms*nc1_loc(iproc),rbuf_index,iproc,mpl%tag)
         call mpl_recv(ms*nc1_loc(iproc),rbuf_dist,iproc,mpl%tag+1)

         ! Format data
         do ic1_loc=1,nc1_loc(iproc)
            ic1 = ic1_s(iproc)+ic1_loc-1
            nn_index(:,ic1) = rbuf_index((ic1_loc-1)*ms+1:ic1_loc*ms)
            nn_dist(:,ic1) = rbuf_dist((ic1_loc-1)*ms+1:ic1_loc*ms)
         end do
 
         ! Release memory
         deallocate(rbuf_index)
         deallocate(rbuf_dist)
      end if
   end do
else
   ! Allocation
   allocate(sbuf_index(ms*nc1_loc(mpl%myproc)))
   allocate(sbuf_dist(ms*nc1_loc(mpl%myproc)))

   ! Format data
   do ic1_loc=1,nc1_loc(mpl%myproc)
      ic1 = ic1_s(mpl%myproc)+ic1_loc-1
      sbuf_index((ic1_loc-1)*ms+1:ic1_loc*ms) = nn_index(:,ic1)
      sbuf_dist((ic1_loc-1)*ms+1:ic1_loc*ms) = nn_dist(:,ic1)
   end do

   ! Send data to ioproc
   call mpl_send(ms*nc1_loc(mpl%myproc),sbuf_index,mpl%ioproc,mpl%tag)
   call mpl_send(ms*nc1_loc(mpl%myproc),sbuf_dist,mpl%ioproc,mpl%tag+1)

   ! Release memory
   deallocate(sbuf_index)
   deallocate(sbuf_dist)
end if
mpl%tag = mpl%tag+2

! Broadcast
call mpl_bcast(nn_index,mpl%ioproc)
call mpl_bcast(nn_dist,mpl%ioproc)

! MPI splitting
do iproc=1,mpl%nproc
   is_s(iproc) = (iproc-1)*(ndata%ns/mpl%nproc+1)+1
   is_e(iproc) = min(iproc*(ndata%ns/mpl%nproc+1),ndata%ns)
   ns_loc(iproc) = is_e(iproc)-is_s(iproc)+1
end do

! Allocation
n_s_max = 100*nint(float(ndata%nc0*ndata%nl0)/float(mpl%nthread*mpl%nproc))
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
   ic1 = ndata%is_to_ic1(is)
   il1 = ndata%is_to_il1(is)
   il0 = ndata%il1_to_il0(il1)

   ! Loop on nearest neighbors
   do i=1,ms
      jc1 = nn_index(i,ic1)
      do jl1=1,ndata%nl1
         if (submask(jc1,jl1)) then
            jl0 = ndata%il1_to_il0(jl1)
            js = ndata%ic1il1_to_is(jc1,jl1)

            ! Only half of the (symmetric) matrix coefficients should be stored
            if (is>js) then
               ! Normalized distance
               distnorm = sqrt(nn_dist(i,ic1)**2/(0.5*(rhs(is)**2+rhs(js)**2)) &
                        & +(ndata%geom%vunit(il0)-ndata%geom%vunit(jl0))**2/(0.5*(rvs(is)**2+rvs(js)**2)))

               if (distnorm<1.0) then
                  ! Gaspari-Cohn (1999) function
                  S_test = gc99(nam,distnorm)

                  ! Check convolution value
                  call check_convol(is,js,S_test,c_n_s(ithread),c(ithread))
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

! Initialize object
ndata%c%prefix = 'c'
ndata%c%n_src = ndata%ns
ndata%c%n_dst = ndata%ns

! Gather data
call convol_gather_data(c_n_s,c,ndata%c)

! Release memory
call delete_ctree(ctree)
deallocate(nn_index)
deallocate(nn_dist)
do ithread=1,mpl%nthread
   call linop_dealloc(c(ithread))
end do

end subroutine compute_convol_distance

!----------------------------------------------------------------------
! Function: gc99
!> Purpose: Gaspari and Cohn (1999) function, with the support radius as a parameter
!----------------------------------------------------------------------
function gc99(nam,distnorm)

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
real(kind_real),intent(in) :: distnorm

! Returned variable
real(kind_real) :: gc99

! Local variable
real(kind_real) :: d

! Distance deformation
d = distnorm+deform*sin(pi*distnorm)

! Distance check bound
if (.not.(d>0.0)) call msgerror('negative normalized distance')

! Square-root
if (nam%lsqrt) d = d*sqrt(2.0)

if (d<0.5) then
   gc99 = 1.0-8.0*d**5+8.0*d**4+5.0*d**3-20.0/3.0*d**2
else if (d<1.0) then
   gc99 = 8.0/3.0*d**5-8.0*d**4+5.0*d**3+20.0/3.0*d**2-10.0*d+4.0-1.0/(3.0*d)
else
   gc99 = 0.0
end if

return

end function gc99

!----------------------------------------------------------------------
! Subroutine: check_convol
!> Purpose: check convolution value and add it if necessary
!----------------------------------------------------------------------
subroutine check_convol(is,js,S_test,c_n_s,c)

implicit none

! Passed variables
integer,intent(in) :: is
integer,intent(in) :: js
real(kind_real),intent(in) :: S_test
integer,intent(inout) :: c_n_s
type(linoptype),intent(inout) :: c         

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
integer :: ithread,offset,iproc
integer :: csum_n_sg(mpl%nproc)
type(linoptype) :: csum

write(mpl%unit,'(a10,a)') '','Gather convolution weights'

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
   offset = 0
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         cout%row(offset+1:offset+csum_n_sg(iproc)) = csum%row
         cout%col(offset+1:offset+csum_n_sg(iproc)) = csum%col
         cout%S(offset+1:offset+csum_n_sg(iproc)) = csum%S
      else
         ! Receive data on ioproc
         call mpl_recv(csum_n_sg(iproc),cout%row(offset+1:offset+csum_n_sg(iproc)),iproc,mpl%tag)
         call mpl_recv(csum_n_sg(iproc),cout%col(offset+1:offset+csum_n_sg(iproc)),iproc,mpl%tag+1)
         call mpl_recv(csum_n_sg(iproc),cout%S(offset+1:offset+csum_n_sg(iproc)),iproc,mpl%tag+2)
      end if

      !  Update offset
      offset = offset+csum_n_sg(iproc)
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
