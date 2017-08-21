!----------------------------------------------------------------------
! Module: module_parameters_interp.f90
!> Purpose: compute NICAS parameters (interpolation)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_parameters_interp

use module_namelist, only: nam
use netcdf
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_reorder
use type_mesh, only: meshtype,create_mesh,mesh_dealloc
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_ndata, only: ndatatype
use type_randgen, only: initialize_sampling,rand_integer

implicit none

real(kind_real),parameter :: S_inf = 1.0e-3 !< Minimum value for the interpolation coefficients

private
public :: compute_interp_h,compute_interp_v,compute_interp_s,interp_horiz

contains

!----------------------------------------------------------------------
! Subroutine: compute_interp_h
!> Purpose: compute basic horizontal interpolation
!----------------------------------------------------------------------
subroutine compute_interp_h(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ic0,ic1,i_s,ibnd,jc0,jc1,progint,il0i
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
real(kind_real) :: renorm(ndata%nc0)
real(kind_real),allocatable :: v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:),valid(:),validg(:,:),missing(:)
type(linoptype) :: hbase,htmp
type(ctreetype) :: ctree

! Compute interpolation
call interp_horiz(ndata,ndata%nc1,ndata%lon(ndata%ic1_to_ic0), &
 & ndata%lat(ndata%ic1_to_ic0),ndata%nc0,ndata%lon,ndata%lat,any(ndata%mask,dim=2),hbase)

! Allocation
allocate(valid(hbase%n_s))
allocate(ndata%h(ndata%nl0i))

do il0i=1,ndata%nl0i
   ! Copy
   call linop_copy(hbase,htmp)

   ! Check interpolation coefficient
   valid = .true.

   ! Check mask boundaries
   if (nam%mask_check) then
      call check_mask_bnd(ndata,htmp,valid,il0i,col_to_ic0=ndata%ic1_to_ic0)
   else
      write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(il0i)
   end if

   ! Renormalization
   renorm = 0.0
   do i_s=1,htmp%n_s
      if (valid(i_s)) renorm(htmp%row(i_s)) = renorm(htmp%row(i_s))+htmp%S(i_s)
   end do

   ! Initialize object
   ndata%h(il0i)%prefix = 'h'
   ndata%h(il0i)%n_src = ndata%nc1
   ndata%h(il0i)%n_dst = ndata%nc0
   ndata%h(il0i)%n_s = count(valid)
   call linop_alloc(ndata%h(il0i))
   ndata%h(il0i)%n_s = 0
   do i_s=1,htmp%n_s
      if (valid(i_s)) then
         ndata%h(il0i)%n_s = ndata%h(il0i)%n_s+1
         ndata%h(il0i)%row(ndata%h(il0i)%n_s) = htmp%row(i_s)
         ndata%h(il0i)%col(ndata%h(il0i)%n_s) = htmp%col(i_s)
         ndata%h(il0i)%S(ndata%h(il0i)%n_s) = htmp%S(i_s)/renorm(htmp%row(i_s))
      end if
   end do

   ! Release memory
   call linop_dealloc(htmp)
   deallocate(valid)

   ! Allocation
   allocate(missing(ndata%nc0))

   ! Count points that are not interpolated
   missing = .false.
   do ic0=1,ndata%nc0
      if (ndata%mask(ic0,il0i)) missing(ic0) = .true.
   end do
   do i_s=1,ndata%h(il0i)%n_s
      missing(ndata%h(il0i)%row(i_s)) = .false.
   end do
   if (count(missing)>0) then
      ! Copy
      call linop_copy(ndata%h(il0i),htmp)

      ! Reallocate interpolation
      call linop_dealloc(ndata%h(il0i))
      ndata%h(il0i)%n_s = ndata%h(il0i)%n_s+count(missing)
      call linop_alloc(ndata%h(il0i))

      ! Fill permanent arrays
      ndata%h(il0i)%row(1:htmp%n_s) = htmp%row
      ndata%h(il0i)%col(1:htmp%n_s) = htmp%col
      ndata%h(il0i)%S(1:htmp%n_s) = htmp%S

      ! Compute cover tree
      allocate(mask_ctree(ndata%nc1))
      mask_ctree = 0
      do ic1=1,ndata%nc1
         if (ndata%mask(ndata%ic1_to_ic0(ic1),il0i)) mask_ctree(ic1) = 1
      end do
      ctree = create_ctree(ndata%nc1,dble(ndata%lon(ndata%ic1_to_ic0)),dble(ndata%lat(ndata%ic1_to_ic0)),mask_ctree)
      deallocate(mask_ctree)

      ! Compute nearest neighbors
      do ic0=1,ndata%nc0
         if (missing(ic0)) then
            htmp%n_s = htmp%n_s+1
            ndata%h(il0i)%row(htmp%n_s) = ic0
            call find_nearest_neighbors(ctree,dble(ndata%lon(ic0)),dble(ndata%lat(ic0)),1,ndata%h(il0i)%col(htmp%n_s:htmp%n_s),dum)
            ndata%h(il0i)%S(htmp%n_s) = 1.0
         end if
      end do

      ! Release memory
      call linop_dealloc(htmp)
      call delete_ctree(ctree)
   end if

   ! Reorder linear operator
   call linop_reorder(ndata%h(il0i))

   ! Release memory
   deallocate(missing)
end do

! Release memory
call linop_dealloc(hbase)

end subroutine compute_interp_h

!----------------------------------------------------------------------
! Subroutine: compute_interp_v
!> Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine compute_interp_v(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: il0,jl0,il1,il0inf,il0sup,il0i,i_s,ic0,ic1
logical,allocatable :: valid(:)
type(linoptype) :: vbase,vtmp

! Linear interpolation
vbase%n_s = ndata%nl1
il0inf = 1
do il0=1,ndata%nl0
   if (ndata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         vbase%n_s = vbase%n_s+2
      end do
      il0inf = il0
   end if
end do
call linop_alloc(vbase)
do il1=1,ndata%nl1
   il0 = ndata%il1_to_il0(il1)
   vbase%row(il1) = il0
   vbase%col(il1) = il0
   vbase%S(il1) = 1.0
end do
vbase%n_s = ndata%nl1
il0inf = 1
do il0=1,ndata%nl0
   if (ndata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         vbase%n_s = vbase%n_s+1
         vbase%row(vbase%n_s) = jl0
         vbase%col(vbase%n_s) = il0inf
         vbase%S(vbase%n_s) = abs(ndata%vunit(il0sup)-ndata%vunit(jl0))/abs(ndata%vunit(il0sup)-ndata%vunit(il0inf))

         vbase%n_s = vbase%n_s+1
         vbase%row(vbase%n_s) = jl0
         vbase%col(vbase%n_s) = il0sup
         vbase%S(vbase%n_s) = abs(ndata%vunit(jl0)-ndata%vunit(il0inf))/abs(ndata%vunit(il0sup)-ndata%vunit(il0inf))
      end do
      il0inf = il0
   end if
end do

! Allocation
allocate(ndata%v(ndata%nl0i))
allocate(valid(vbase%n_s))

do il0i=1,ndata%nl0i
   ! Initialize vertical interpolation
   ndata%v(il0i)%prefix = 'v'
   ndata%v(il0i)%n_src = ndata%nl1
   ndata%v(il0i)%n_dst = ndata%nl0

   if (ndata%nl0i==1) then
      ! Copy basic vertical interpolation
      ndata%v(il0i)%n_s = vbase%n_s
      call linop_alloc(ndata%v(il0i))
      ndata%v(il0i)%row = vbase%row
      ndata%v(il0i)%col = vbase%col
      ndata%v(il0i)%S = vbase%S
   else
      ! Check valid operations
      do i_s=1,vbase%n_s
         valid(i_s) = (vbase%row(i_s)<=il0i).and.(vbase%col(i_s)<=il0i)
      end do

      ! Copy valid operations
      vtmp%n_s = count(valid)
      call linop_alloc(vtmp)
      vtmp%n_s = 0
      do i_s=1,vbase%n_s
         if (valid(i_s)) then
            vtmp%n_s = vtmp%n_s+1
            vtmp%row(vtmp%n_s) = vbase%row(i_s)
            vtmp%col(vtmp%n_s) = vbase%col(i_s)
            vtmp%S(vtmp%n_s) = vbase%S(i_s)
         end if
      end do

      ! Add missing levels to the interpolation
      ndata%v(il0i)%n_s = vtmp%n_s+(il0i-maxval(vtmp%col))
      call linop_alloc(ndata%v(il0i))
      ndata%v(il0i)%row(1:vtmp%n_s) = vtmp%row
      ndata%v(il0i)%col(1:vtmp%n_s) = vtmp%col
      ndata%v(il0i)%S(1:vtmp%n_s) = vtmp%S
      do jl0=maxval(vtmp%col)+1,il0i
         vtmp%n_s = vtmp%n_s+1
         ndata%v(il0i)%row(vtmp%n_s) = jl0
         ndata%v(il0i)%col(vtmp%n_s) = jl0
         ndata%v(il0i)%S(vtmp%n_s) = 1.0
      end do

      ! Release memory
      call linop_dealloc(vtmp)
   end if

   ! Conversion
   ndata%v(il0i)%col = ndata%il0_to_il1(ndata%v(il0i)%col)

   ! Reorder linear operator
   call linop_reorder(ndata%v(il0i))
end do

! Release memory
call linop_dealloc(vbase)
deallocate(valid)

! Find the bottom for each point of S1
allocate(ndata%vbot(ndata%nc1))
!$omp parallel do private(ic1,ic0,il0)
do ic1=1,ndata%nc1
   ic0 = ndata%ic1_to_ic0(ic1)
   il0 = 1
   do while (ndata%mask(ic0,il0).and.(il0<ndata%nl0))
      il0 = il0+1
   end do
   ndata%vbot(ic1) = min(il0,ndata%nl0i)
end do
!$omp end parallel do

end subroutine compute_interp_v

!----------------------------------------------------------------------
! Subroutine: compute_interp_s
!> Purpose: compute horizontal subsampling interpolation
!----------------------------------------------------------------------
subroutine compute_interp_s(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ic1,ic2,il1,i_s,ibnd,il0,ic0,jc2,jc1,jc0,progint
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
real(kind_real) :: renorm(ndata%nc1)
real(kind_real),allocatable :: v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:),valid(:),validg(:,:),missing(:)
type(linoptype) :: stmp
type(ctreetype) :: ctree

! Allocation
allocate(ndata%s(ndata%nl1))

do il1=1,ndata%nl1
   ! Initialize object
   ndata%s(il1)%prefix = 's'
   ndata%s(il1)%n_src = ndata%nc2(il1)
   ndata%s(il1)%n_dst = ndata%nc1

   if (ndata%nc2(il1)==ndata%nc1) then
      ! No interpolation
      ndata%s(il1)%n_s = ndata%nc2(il1)
      call linop_alloc(ndata%s(il1))
      do i_s=1,ndata%s(il1)%n_s
         ndata%s(il1)%row(i_s) = i_s
         ndata%s(il1)%col(i_s) = i_s
         ndata%s(il1)%S(i_s) = 1.0
      end do
   else
      ! Compute interpolation
      call interp_horiz(ndata,ndata%nc2(il1),ndata%lon(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1)), & 
    & ndata%lat(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1)),ndata%nc1,ndata%lon(ndata%ic1_to_ic0), &
    & ndata%lat(ndata%ic1_to_ic0),ndata%mask(ndata%ic1_to_ic0,ndata%il1_to_il0(il1)),stmp)
   
      ! Allocation
      allocate(valid(stmp%n_s))
   
      ! Check interpolation coefficient
      valid = .not.(stmp%S<S_inf)
   
      ! Check mask boundaries
      if (nam%mask_check) then
         call check_mask_bnd(ndata,stmp,valid,ndata%il1_to_il0(il1), &
       & row_to_ic0=ndata%ic1_to_ic0,col_to_ic0=ndata%ic1_to_ic0(ndata%ic2il1_to_ic1(1:ndata%nc2(il1),il1)))
      else
         write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(ndata%il1_to_il0(il1))
      end if
   
      ! Renormalization
      renorm = 0.0
      do i_s=1,stmp%n_s
         if (valid(i_s)) renorm(stmp%row(i_s)) = renorm(stmp%row(i_s))+stmp%S(i_s)
      end do
   
      ! Copy valid operations
      ndata%s(il1)%n_s = count(valid)
      call linop_alloc(ndata%s(il1))
      ndata%s(il1)%n_s = 0
      do i_s=1,stmp%n_s
         if (valid(i_s)) then
            ndata%s(il1)%n_s = ndata%s(il1)%n_s+1
            ndata%s(il1)%row(ndata%s(il1)%n_s) = stmp%row(i_s)
            ndata%s(il1)%col(ndata%s(il1)%n_s) = stmp%col(i_s)
            ndata%s(il1)%S(ndata%s(il1)%n_s) = stmp%S(i_s)/renorm(stmp%row(i_s))
         end if
      end do
   
      ! Release memory
      call linop_dealloc(stmp)
      deallocate(valid)
   
      ! Allocation
      allocate(missing(ndata%nc1))
   
      ! Count points that are not interpolated
      missing = .false.
      do ic1=1,ndata%nc1
         if (ndata%mask(ndata%ic1_to_ic0(ic1),ndata%il1_to_il0(il1))) missing(ic1) = .true.
      end do
      do i_s=1,ndata%s(il1)%n_s
         missing(ndata%s(il1)%row(i_s)) = .false.
      end do
      if (count(missing)>0) then
         ! Copy
         call linop_copy(ndata%s(il1),stmp)
   
         ! Reallocate permanent arrays
         call linop_dealloc(ndata%s(il1))
         ndata%s(il1)%n_s = ndata%s(il1)%n_s+count(missing)
         call linop_alloc(ndata%s(il1))
   
         ! Fill permanent arrays
         ndata%s(il1)%row(1:stmp%n_s) = stmp%row
         ndata%s(il1)%col(1:stmp%n_s) = stmp%col
         ndata%s(il1)%S(1:stmp%n_s) = stmp%S
   
         ! Compute cover tree
         allocate(mask_ctree(ndata%nc2(il1)))
         mask_ctree = 0
         do ic2=1,ndata%nc2(il1)
            if (ndata%mask(ndata%ic1_to_ic0(ndata%ic2il1_to_ic1(ic2,il1)),ndata%il1_to_il0(il1))) mask_ctree(ic2) = 1
         end do
         ctree = create_ctree(ndata%nc2(il1),dble(ndata%lon(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1))), &
       & dble(ndata%lat(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1))),mask_ctree)
         deallocate(mask_ctree)
   
         ! Compute nearest neighbors
         do ic1=1,ndata%nc1
            if (missing(ic1)) then
               stmp%n_s = stmp%n_s+1
               ndata%s(il1)%row(stmp%n_s) = ic1
               call find_nearest_neighbors(ctree,dble(ndata%lon(ndata%ic1_to_ic0(ic1))), &
             & dble(ndata%lat(ndata%ic1_to_ic0(ic1))),1,ndata%s(il1)%col(stmp%n_s:stmp%n_s),dum)
               ndata%s(il1)%S(stmp%n_s) = 1.0
            end if
         end do
   
         ! Release memory
         call linop_dealloc(stmp)
         call delete_ctree(ctree)
      end if
   
      ! Reorder linear operator
      call linop_reorder(ndata%s(il1))
   
      ! Release memory
      deallocate(missing)
   end if
end do

end subroutine compute_interp_s

!----------------------------------------------------------------------
! Subroutine: interp_horiz
!> Purpose: compute horizontal interpolation
!----------------------------------------------------------------------
subroutine interp_horiz(ndata,n_src,lon_src,lat_src,n_dst,lon_dst,lat_dst,mask_dst,interp)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata
integer,intent(in) :: n_src
real(kind_real),intent(in) :: lon_src(n_src)
real(kind_real),intent(in) :: lat_src(n_src)
integer,intent(in) :: n_dst
real(kind_real),intent(in) :: lon_dst(n_dst)
real(kind_real),intent(in) :: lat_dst(n_dst)
logical,intent(in) :: mask_dst(n_dst)
type(linoptype),intent(inout) :: interp

! Local variables
integer :: info,i_src,j_src,k_src,i,ibnd,il0,nt,i_dst,inn(1),n_s,i_s,offset,progint
integer :: ib(3)
integer :: iproc,i_dst_s(mpl%nproc),i_dst_e(mpl%nproc),n_dst_loc(mpl%nproc),i_dst_loc,n_sg(mpl%nproc)
integer,allocatable :: mask_ctree(:),row(:),col(:)
real(kind_real) :: dist(1),p(3),b(3)
real(kind_real),allocatable :: S(:)
logical :: init
logical,allocatable :: done(:)
type(ctreetype) :: ctree
type(meshtype) :: mesh

! Create mesh
call create_mesh(ndata%rng,n_src,lon_src,lat_src,.false.,mesh)

! Compute cover tree
allocate(mask_ctree(mesh%nnr))
mask_ctree = 1
ctree = create_ctree(mesh%nnr,dble(lon_src(mesh%order)),dble(lat_src(mesh%order)),mask_ctree)

! MPI splitting
do iproc=1,mpl%nproc
   i_dst_s(iproc) = (iproc-1)*(n_dst/mpl%nproc+1)+1
   i_dst_e(iproc) = min(iproc*(n_dst/mpl%nproc+1),n_dst)
   n_dst_loc(iproc) = i_dst_e(iproc)-i_dst_s(iproc)+1
end do

! Allocation
allocate(row(3*n_dst_loc(mpl%myproc)))
allocate(col(3*n_dst_loc(mpl%myproc)))
allocate(S(3*n_dst_loc(mpl%myproc)))
allocate(done(n_dst_loc(mpl%myproc)))

! Compute interpolation
write(mpl%unit,'(a10,a)',advance='no') '','Compute interpolation: '
call prog_init(progint,done)
n_s = 0
do i_dst_loc=1,n_dst_loc(mpl%myproc)
   ! Indices
   i_dst = i_dst_s(mpl%myproc)+i_dst_loc-1

   if (mask_dst(i_dst)) then
      ! Find nearest neighbor
      call find_nearest_neighbors(ctree,dble(lon_dst(i_dst)), &
    & dble(lat_dst(i_dst)),1,inn,dist)

      if (dist(1)<tiny(1.0)) then
         ! Subset point
         n_s = n_s+1
         row(n_s) = i_dst
         col(n_s) = mesh%order(inn(1))
         S(n_s) = 1.0
      else
         ! Transform to cartesian coordinates
         call trans(1,lat_dst(i_dst),lon_dst(i_dst),p(1),p(2),p(3))

         ! Compute barycentric coordinates
         call trfind(inn(1),dble(p),mesh%nnr,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend, &
       & b(1),b(2),b(3),ib(1),ib(2),ib(3))

         if (all(ib>0)) then
            ! Valid interpolation
            if (sum(b)>0.0) then
               ! Normalize barycentric coordinates
               b = b/sum(b)

               ! Add interpolation elements
               do i=1,3
                  if (b(i)>S_inf) then
                     n_s = n_s+1
                     row(n_s) = i_dst
                     col(n_s) = mesh%order(ib(i))
                     S(n_s) = b(i)
                  end if
               end do
            end if
         end if
      end if
   end if

   done(i_dst_loc) = .true.
   call prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         n_sg(iproc) = n_s
      else
         ! Receive data on ioproc
         call mpl_recv(n_sg(iproc),iproc,mpl%tag)
      end if
   end do
else
   ! Send data to ioproc
   call mpl_send(n_s,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl_bcast(n_sg,mpl%ioproc)

! Allocation
interp%n_s = sum(n_sg)
interp%n_src = n_src
interp%n_dst = n_dst
call linop_alloc(interp)

! Communication
if (mpl%main) then
   offset = 0
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         interp%row(offset+1:offset+n_sg(iproc)) = row(1:n_sg(iproc))
         interp%col(offset+1:offset+n_sg(iproc)) = col(1:n_sg(iproc))
         interp%S(offset+1:offset+n_sg(iproc)) = S(1:n_sg(iproc))
      else
         ! Receive data on ioproc
         call mpl_recv(n_sg(iproc),interp%row(offset+1:offset+n_sg(iproc)),iproc,mpl%tag)
         call mpl_recv(n_sg(iproc),interp%col(offset+1:offset+n_sg(iproc)),iproc,mpl%tag+1)
         call mpl_recv(n_sg(iproc),interp%S(offset+1:offset+n_sg(iproc)),iproc,mpl%tag+2)
      end if

      ! Update offset
      offset = offset+n_sg(iproc)
   end do
else
   ! Send data to ioproc
   call mpl_send(n_s,row(1:n_s),mpl%ioproc,mpl%tag)
   call mpl_send(n_s,col(1:n_s),mpl%ioproc,mpl%tag+1)
   call mpl_send(n_s,S(1:n_s),mpl%ioproc,mpl%tag+2)
end if
mpl%tag = mpl%tag+3

! Broadcast data
call mpl_bcast(interp%row,mpl%ioproc)
call mpl_bcast(interp%col,mpl%ioproc)
call mpl_bcast(interp%S,mpl%ioproc)

! Release memory
call delete_ctree(ctree)
call mesh_dealloc(mesh)

end subroutine interp_horiz

!----------------------------------------------------------------------
! Subroutine: check_mask_bnd
!> Purpose: check mask boundaries for interpolations
!----------------------------------------------------------------------
subroutine check_mask_bnd(ndata,interp,valid,il0,row_to_ic0,col_to_ic0)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata          !< Sampling data
type(linoptype),intent(inout) :: interp          !< Interpolation
logical,intent(inout) :: valid(interp%n_s)   !< Valid vector
integer,intent(in),optional :: row_to_ic0(interp%n_dst) !< Conversion from row to ic0 (identity if missing)
integer,intent(in),optional :: col_to_ic0(interp%n_src) !< Conversion from col to ic0 (identity if missing)

! Local variables
integer :: ic0,ic1,i_s,ibnd,jc0,jc1,progint,il0
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
real(kind_real) :: renorm(ndata%nc0)
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:),missing(:)

! MPI splitting
do iproc=1,mpl%nproc
   i_s_s(iproc) = (iproc-1)*(interp%n_s/mpl%nproc+1)+1
   i_s_e(iproc) = min(iproc*(interp%n_s/mpl%nproc+1),interp%n_s)
   n_s_loc(iproc) = i_s_e(iproc)-i_s_s(iproc)+1
end do

! Allocation
allocate(done(n_s_loc(mpl%myproc)))

! Check that interpolations are not crossing mask boundaries
write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),': '
call prog_init(progint,done)
!$omp parallel do private(i_s_loc,i_s,x,y,z,v1,v2,va,vp,t,ic0,jc1,jc0)
do i_s_loc=1,n_s_loc(mpl%myproc)
   ! Indices
   i_s = i_s_s(mpl%myproc)+i_s_loc-1

   if (valid(i_s)) then
      ! Allocation
      allocate(x(2))
      allocate(y(2))
      allocate(z(2))
      allocate(v1(3))
      allocate(v2(3))
      allocate(va(3))
      allocate(vp(3))
      allocate(t(4))

      ! Indices
      if (present(row_to_ic0)) then
         ic0 = row_to_ic0(interp%row(i_s))
      else
         ic0 = interp%row(i_s)
      end if
      if (present(col_to_ic0)) then
         jc0 = col_to_ic0(interp%col(i_s))
      else
         jc0 = interp%col(i_s)
      end if

      ! Transform to cartesian coordinates
      call trans(2,ndata%lat((/ic0,jc0/)),ndata%lon((/ic0,jc0/)),x,y,z)

      ! Compute arc orthogonal vector
      v1 = (/x(1),y(1),z(1)/)
      v2 = (/x(2),y(2),z(2)/)
      call vector_product(v1,v2,va)

      ! Check if arc is crossing boundary arcs
      do ibnd=1,ndata%nbnd(il0)
         call vector_product(va,ndata%vbnd(:,ibnd,il0),vp)
         v1 = (/x(1),y(1),z(1)/)
         call vector_triple_product(v1,va,vp,t(1))
         v1 = (/x(2),y(2),z(2)/)
         call vector_triple_product(v1,va,vp,t(2))
         v1 = (/ndata%xbnd(1,ibnd,il0),ndata%ybnd(1,ibnd,il0),ndata%zbnd(1,ibnd,il0)/)
         call vector_triple_product(v1,ndata%vbnd(:,ibnd,il0),vp,t(3))
         v1 = (/ndata%xbnd(2,ibnd,il0),ndata%ybnd(2,ibnd,il0),ndata%zbnd(2,ibnd,il0)/)
         call vector_triple_product(v1,ndata%vbnd(:,ibnd,il0),vp,t(4))
         t(1) = -t(1)
         t(3) = -t(3)
         if (all(t>0).or.(all(t<0))) then
            valid(i_s) = .false.
            exit
         end if
      end do

      ! Memory release
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(v1)
      deallocate(v2)
      deallocate(va)
      deallocate(vp)
      deallocate(t)
   end if

   ! Print progression
   done(i_s_loc) = .true.
   call prog_print(progint,done)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Broadcast
do iproc=1,mpl%nproc
   call mpl_bcast(valid(i_s_s(iproc):i_s_e(iproc)),iproc)
end do

! Release memory
deallocate(done)

end subroutine check_mask_bnd

end module module_parameters_interp
