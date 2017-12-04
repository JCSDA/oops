!----------------------------------------------------------------------
! Module: tools_interp.f90
!> Purpose: horizontal interpolation tools
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_interp

use omp_lib
use tools_const, only: vector_product,vector_triple_product
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_stripack, only: trans,trmesh,trfind
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_reorder
use type_mesh, only: meshtype,create_mesh,mesh_dealloc
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send

implicit none

real(kind_real),parameter :: S_inf = 1.0e-3 !< Minimum value for the interpolation coefficients

interface compute_interp_bilin
  module procedure compute_interp_bilin_from_lat_lon
  module procedure compute_interp_bilin_from_mesh_ctree
end interface

private
public :: compute_interp_bilin,compute_grid_interp_bilin,check_mask_bnd

contains

!----------------------------------------------------------------------
! Subroutine: compute_interp_bilin_from_lat_lon
!> Purpose: compute horizontal interpolation from source latitude/longitude
!----------------------------------------------------------------------
subroutine compute_interp_bilin_from_lat_lon(n_src,lon_src,lat_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp)

implicit none

! Passed variables
integer,intent(in) :: n_src                  !< Source size
real(kind_real),intent(in) :: lon_src(n_src) !< Source longitudes
real(kind_real),intent(in) :: lat_src(n_src) !< Source latitudes
logical,intent(in) :: mask_src(n_src)        !< Source mask
integer,intent(in) :: n_dst                  !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) !< Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst) !< Destination latitudes
logical,intent(in) :: mask_dst(n_dst)        !< Destination mask
type(linoptype),intent(inout) :: interp      !< Interpolation data

! Local variables
logical,allocatable :: mask_ctree(:)
type(ctreetype) :: ctree
type(meshtype) :: mesh

! Create mesh
call create_mesh(n_src,lon_src,lat_src,.false.,mesh)

! Compute cover tree
allocate(mask_ctree(mesh%nnr))
mask_ctree = .true.
ctree = create_ctree(mesh%nnr,dble(lon_src(mesh%order)),dble(lat_src(mesh%order)),mask_ctree)
deallocate(mask_ctree)

! Compute interpolation
call compute_interp_bilin_from_mesh_ctree(mesh,ctree,n_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp)

! Release memory
call delete_ctree(ctree)
call mesh_dealloc(mesh)

end subroutine compute_interp_bilin_from_lat_lon

!----------------------------------------------------------------------
! Subroutine: compute_interp_bilin_from_mesh_ctree
!> Purpose: compute horizontal interpolation from source mesh and ctree
!----------------------------------------------------------------------
subroutine compute_interp_bilin_from_mesh_ctree(mesh,ctree,n_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp)

implicit none

! Passed variables
type(meshtype),intent(in) :: mesh            !< Mesh
type(ctreetype),intent(in) :: ctree          !< Cover tree
integer,intent(in) :: n_src                  !< Source size
logical,intent(in) :: mask_src(n_src)        !< Source mask
integer,intent(in) :: n_dst                  !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) !< Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst) !< Destination latitudes
logical,intent(in) :: mask_dst(n_dst)        !< Destination mask
type(linoptype),intent(inout) :: interp      !< Interpolation data

! Local variables
integer :: i,i_dst,inn(1),n_s,offset,progint
integer :: ib(3)
integer :: iproc,i_dst_s(mpl%nproc),i_dst_e(mpl%nproc),n_dst_loc(mpl%nproc),i_dst_loc,n_sg(mpl%nproc)
integer,allocatable :: row(:),col(:)
real(kind_real) :: dist(1),p(3),b(3)
real(kind_real),allocatable :: S(:)
logical,allocatable :: done(:)

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

      if (abs(dist(1))>0.0) then
         ! Transform to cartesian coordinates
         call trans(1,lat_dst(i_dst),lon_dst(i_dst),p(1),p(2),p(3))

         ! Compute barycentric coordinates
         call trfind(mesh%order_inv(inn(1)),dble(p),mesh%nnr,mesh%x,mesh%y,mesh%z,mesh%list,mesh%lptr,mesh%lend, &
       & b(1),b(2),b(3),ib(1),ib(2),ib(3))

         if (all(ib>0)) then
            if (all(mask_src(mesh%order(ib)))) then
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
      else
         ! Subsampled point
         n_s = n_s+1
         row(n_s) = i_dst
         col(n_s) = mesh%order(inn(1))
         S(n_s) = 1.0
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

end subroutine compute_interp_bilin_from_mesh_ctree

!----------------------------------------------------------------------
! Subroutine: compute_grid_interp_bilin
!> Purpose: compute grid horizontal interpolation
!----------------------------------------------------------------------
subroutine compute_grid_interp_bilin(geom,nc1,ic1_to_ic0,mask_check,vbot,vtop,h)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom             !< Geometry
integer,intent(in) :: nc1                     !< Subset Sc1 size
integer,intent(in) :: ic1_to_ic0(nc1)         !< Subset Sc1 to subset Sc0
logical,intent(in) :: mask_check              !< Mask check key
integer,intent(in) :: vbot(nc1)               !< Bottom level
integer,intent(in) :: vtop(nc1)               !< Top level
type(linoptype),intent(inout) :: h(geom%nl0i) !< Horizontal interpolation data

! Local variables
integer :: ic0,ic1,i_s,il0i
real(kind_real) :: dum(1)
real(kind_real) :: renorm(geom%nc0)
logical :: test_nc0(geom%nc0),test_nc1(nc1)
logical,allocatable :: mask_extra(:),valid(:),missing(:)
type(linoptype) :: hbase,htmp
type(ctreetype) :: ctree

! Compute interpolation
call compute_interp_bilin(nc1,geom%lon(ic1_to_ic0),geom%lat(ic1_to_ic0), any(geom%mask(ic1_to_ic0,:),dim=2), &
 & geom%nc0,geom%lon,geom%lat,any(geom%mask,dim=2),hbase)

! Allocation
allocate(valid(hbase%n_s))
allocate(mask_extra(nc1))

do il0i=1,geom%nl0i
   ! Copy
   call linop_copy(hbase,htmp)
   valid = .true.

   ! Check mask boundaries
   if (mask_check) then
      write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Sublevel ',il0i,': '
      call check_mask_bnd(geom,htmp,valid,il0i,col_to_ic0=ic1_to_ic0)
   else
      write(mpl%unit,'(a10,a,i3)') '','Sublevel ',il0i
   end if

   if (geom%nl0i>1) then
      ! Check extrapolated points because of masks vertical variations (unsampled levels)
      mask_extra = .false.
      do ic1=1,nc1
         ic0 = ic1_to_ic0(ic1)
         if (geom%mask(ic0,il0i).and.((il0i<vbot(ic1)).or.(il0i>vtop(ic1)))) mask_extra(ic1) = .true.
      end do

      ! Remove operations for extrapolated points
      do i_s=1,htmp%n_s
         if (valid(i_s)) then
            ic1 = hbase%col(i_s)
            if (mask_extra(ic1)) valid(i_s) = .false.
         end if
      end do
      if (count(mask_extra)>0) write(mpl%unit,'(a10,a,i5)') '','Extrapolated points: ',count(mask_extra)
   end if

   ! Renormalization
   renorm = 0.0
   do i_s=1,htmp%n_s
      if (valid(i_s)) renorm(htmp%row(i_s)) = renorm(htmp%row(i_s))+htmp%S(i_s)
   end do

   ! Initialize object
   h(il0i)%prefix = 'h'
   h(il0i)%n_src = nc1
   h(il0i)%n_dst = geom%nc0
   h(il0i)%n_s = count(valid)
   call linop_alloc(h(il0i))
   h(il0i)%n_s = 0
   do i_s=1,htmp%n_s
      if (valid(i_s)) then
         h(il0i)%n_s = h(il0i)%n_s+1
         h(il0i)%row(h(il0i)%n_s) = htmp%row(i_s)
         h(il0i)%col(h(il0i)%n_s) = htmp%col(i_s)
         h(il0i)%S(h(il0i)%n_s) = htmp%S(i_s)/renorm(htmp%row(i_s))
      end if
   end do

   ! Release memory
   call linop_dealloc(htmp)

   ! Allocation
   allocate(missing(geom%nc0))

   ! Count points that are not interpolated
   missing = .false.
   do ic0=1,geom%nc0
      if (geom%mask(ic0,il0i)) missing(ic0) = .true.
   end do
   do i_s=1,h(il0i)%n_s
      missing(h(il0i)%row(i_s)) = .false.
   end do
   if (count(missing)>0) then
      ! Copy
      call linop_copy(h(il0i),htmp)

      ! Reallocate interpolation
      call linop_dealloc(h(il0i))
      h(il0i)%n_s = h(il0i)%n_s+count(missing)
      call linop_alloc(h(il0i))

      ! Fill permanent arrays
      h(il0i)%row(1:htmp%n_s) = htmp%row
      h(il0i)%col(1:htmp%n_s) = htmp%col
      h(il0i)%S(1:htmp%n_s) = htmp%S

      ! Compute cover tree
      ctree = create_ctree(nc1,dble(geom%lon(ic1_to_ic0)),dble(geom%lat(ic1_to_ic0)), &
            & geom%mask(ic1_to_ic0,il0i).and.(.not.mask_extra))

      ! Compute nearest neighbors
      do ic0=1,geom%nc0
         if (missing(ic0)) then
            htmp%n_s = htmp%n_s+1
            h(il0i)%row(htmp%n_s) = ic0
            call find_nearest_neighbors(ctree,dble(geom%lon(ic0)),dble(geom%lat(ic0)),1,h(il0i)%col(htmp%n_s:htmp%n_s),dum)
            h(il0i)%S(htmp%n_s) = 1.0
         end if
      end do

      ! Release memory
      call linop_dealloc(htmp)
      call delete_ctree(ctree)
   end if

   ! Reorder linear operator
   call linop_reorder(h(il0i))

   ! Release memory
   deallocate(missing)

   ! Test interpolation
   test_nc0 = geom%mask(:,min(il0i,geom%nl0i))
   test_nc1 = .true.
   do i_s=1,h(il0i)%n_s
      test_nc0(h(il0i)%row(i_s)) = .false.
      test_nc1(h(il0i)%col(i_s)) = .false.
   end do
   if (any(test_nc0)) call msgerror('error with the grid interpolation row')
   if (any(test_nc1)) call msgerror('error with the grid interpolation col')
end do

! Release memory
call linop_dealloc(hbase)
deallocate(valid)
deallocate(mask_extra)

end subroutine compute_grid_interp_bilin

!----------------------------------------------------------------------
! Subroutine: check_mask_bnd
!> Purpose: check mask boundaries for interpolations
!----------------------------------------------------------------------
subroutine check_mask_bnd(geom,interp,valid,il0,row_to_ic0,col_to_ic0)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                       !< Geometry
type(linoptype),intent(inout) :: interp                 !< Interpolation data
logical,intent(inout) :: valid(interp%n_s)              !< Valid points
integer,intent(in),optional :: row_to_ic0(interp%n_dst) !< Conversion from row to ic0 (identity if missing)
integer,intent(in),optional :: col_to_ic0(interp%n_src) !< Conversion from col to ic0 (identity if missing)

! Local variables
integer :: ic0,i_s,ibnd,jc0,jc1,progint,il0
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:)

! MPI splitting
do iproc=1,mpl%nproc
   i_s_s(iproc) = (iproc-1)*(interp%n_s/mpl%nproc+1)+1
   i_s_e(iproc) = min(iproc*(interp%n_s/mpl%nproc+1),interp%n_s)
   n_s_loc(iproc) = i_s_e(iproc)-i_s_s(iproc)+1
end do

! Allocation
allocate(done(n_s_loc(mpl%myproc)))

! Check that interpolations are not crossing mask boundaries
call prog_init(progint,done)
!$omp parallel do schedule(static) private(i_s_loc,i_s,x,y,z,v1,v2,va,vp,t,ic0,jc1,jc0)
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
      call trans(2,geom%lat((/ic0,jc0/)),geom%lon((/ic0,jc0/)),x,y,z)

      ! Compute arc orthogonal vector
      v1 = (/x(1),y(1),z(1)/)
      v2 = (/x(2),y(2),z(2)/)
      call vector_product(v1,v2,va)

      ! Check if arc is crossing boundary arcs
      do ibnd=1,geom%nbnd(il0)
         call vector_product(va,geom%vbnd(:,ibnd,il0),vp)
         v1 = (/x(1),y(1),z(1)/)
         call vector_triple_product(v1,va,vp,t(1))
         v1 = (/x(2),y(2),z(2)/)
         call vector_triple_product(v1,va,vp,t(2))
         v1 = (/geom%xbnd(1,ibnd,il0),geom%ybnd(1,ibnd,il0),geom%zbnd(1,ibnd,il0)/)
         call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(3))
         v1 = (/geom%xbnd(2,ibnd,il0),geom%ybnd(2,ibnd,il0),geom%zbnd(2,ibnd,il0)/)
         call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(4))
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

end module tools_interp
