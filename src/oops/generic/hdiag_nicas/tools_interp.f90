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
use tools_stripack, only: trans
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_reorder
use type_mesh, only: meshtype,create_mesh,mesh_dealloc,copy_mesh,barycentric,addnode,polygon
use type_mpl, only: mpl,mpl_bcast,mpl_allgather,mpl_split,mpl_send,mpl_recv

implicit none

integer,parameter :: nnatmax = 40           !< Maximum number of natural neighbors
real(kind_real),parameter :: S_inf = 1.0e-3 !< Minimum value for the interpolation coefficients

interface compute_interp
  module procedure compute_interp_from_lat_lon
  module procedure compute_interp_from_mesh_ctree
end interface

private
public :: compute_interp,compute_grid_interp,check_mask_bnd,check_arc,interp_missing

contains

!----------------------------------------------------------------------
! Subroutine: compute_interp_from_lat_lon
!> Purpose: compute horizontal interpolation from source latitude/longitude
!----------------------------------------------------------------------
subroutine compute_interp_from_lat_lon(n_src,lon_src,lat_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp_type,interp)

implicit none

! Passed variables
integer,intent(in) :: n_src                   !< Source size
real(kind_real),intent(in) :: lon_src(n_src)  !< Source longitudes
real(kind_real),intent(in) :: lat_src(n_src)  !< Source latitudes
logical,intent(in) :: mask_src(n_src)         !< Source mask
integer,intent(in) :: n_dst                   !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst)  !< Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst)  !< Destination latitudes
logical,intent(in) :: mask_dst(n_dst)         !< Destination mask
character(len=1024),intent(in) :: interp_type !< Interpolation type
type(linoptype),intent(inout) :: interp       !< Interpolation data

! Local variables
integer :: n_src_eff,i_src,i_src_eff
integer,allocatable :: src_eff_to_src(:)
logical,allocatable :: mask_ctree(:),mask_src_eff(:)
type(ctreetype) :: ctree
type(meshtype) :: mesh

! Count non-missing source points
n_src_eff = count(mask_src)

! Allocation
allocate(src_eff_to_src(n_src_eff))
allocate(mask_src_eff(n_src_eff))

! Conversion
i_src_eff = 0
do i_src=1,n_src
   if (mask_src(i_src)) then
      i_src_eff = i_src_eff+1
      src_eff_to_src(i_src_eff) = i_src
   end if
end do

! Create mesh
call create_mesh(n_src_eff,lon_src(src_eff_to_src),lat_src(src_eff_to_src),.false.,mesh)

! Compute cover tree
allocate(mask_ctree(n_src_eff))
mask_ctree = .true.
ctree = create_ctree(n_src_eff,dble(lon_src(src_eff_to_src)),dble(lat_src(src_eff_to_src)),mask_ctree)
deallocate(mask_ctree)

! Compute interpolation
mask_src_eff = .true.
call compute_interp_from_mesh_ctree(mesh,ctree,n_src_eff,mask_src_eff,n_dst,lon_dst,lat_dst,mask_dst,interp_type,interp)

! Effective points conversion
interp%n_src = n_src
interp%col = src_eff_to_src(interp%col)

! Release memory
call delete_ctree(ctree)
call mesh_dealloc(mesh)

end subroutine compute_interp_from_lat_lon

!----------------------------------------------------------------------
! Subroutine: compute_interp_from_mesh_ctree
!> Purpose: compute horizontal interpolation from source mesh and ctree
!----------------------------------------------------------------------
subroutine compute_interp_from_mesh_ctree(mesh,ctree,n_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp_type,interp)

implicit none

! Passed variables
type(meshtype),intent(in) :: mesh             !< Mesh
type(ctreetype),intent(in) :: ctree           !< Cover tree
integer,intent(in) :: n_src                   !< Source size
logical,intent(in) :: mask_src(n_src)         !< Source mask
integer,intent(in) :: n_dst                   !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst)  !< Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst)  !< Destination latitudes
logical,intent(in) :: mask_dst(n_dst)         !< Destination mask
character(len=1024),intent(in) :: interp_type !< Interpolation type
type(linoptype),intent(inout) :: interp       !< Interpolation data

! Local variables
integer :: i,i_src,i_dst,nn_index(1),n_s,progint,ib(3),nnat,inat,np,iproc,offset
integer :: i_dst_s(mpl%nproc),i_dst_e(mpl%nproc),n_dst_loc(mpl%nproc),i_dst_loc,n_sg(mpl%nproc)
integer,allocatable :: natis(:),row(:),col(:)
real(kind_real) :: nn_dist(1),b(3)
real(kind_real),allocatable :: area_polygon(:),area_polygon_new(:),natwgt(:),S(:)
logical :: loop
logical,allocatable :: done(:)
type(meshtype) :: meshnew

! MPI splitting
call mpl_split(n_dst,i_dst_s,i_dst_e,n_dst_loc)

! Allocation
call msi(np)
if (trim(interp_type)=='bilin') then
   ! Bilinear interpolation
   np = 3
elseif (trim(interp_type)=='natural') then
   ! Natural neighbors
   np = nnatmax
   allocate(area_polygon(mesh%nnr))
   allocate(area_polygon_new(nnatmax))
   allocate(natis(mesh%nnr))
   allocate(natwgt(nnatmax))
else
   call msgerror('wrong interpolation type')
end if
allocate(row(np*n_dst_loc(mpl%myproc)))
allocate(col(np*n_dst_loc(mpl%myproc)))
allocate(S(np*n_dst_loc(mpl%myproc)))
allocate(done(n_dst_loc(mpl%myproc)))

if (trim(interp_type)=='natural') then
   ! Compute polygons areas
   do i_src=1,mesh%nnr
      natis(i_src) = i_src
   end do
   call polygon(mesh,mesh%nnr,natis,area_polygon)
end if

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
    & dble(lat_dst(i_dst)),1,nn_index,nn_dist)

      if (abs(nn_dist(1))>0.0) then
         ! Compute barycentric coordinates
         call barycentric(lon_dst(i_dst),lat_dst(i_dst),mesh,nn_index(1),b,ib)

         if (all(ib>0)) then
            if (all(mask_src(mesh%order(ib)))) then
               ! Valid interpolation
               if (trim(interp_type)=='bilin') then
                  ! Bilinear interpolation
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
               elseif (trim(interp_type)=='natural') then
                  ! Natural neighbors interpolation

                  ! Copy mesh
                  call copy_mesh(mesh,meshnew)

                  ! Add a node
                  call addnode(lon_dst(i_dst),lat_dst(i_dst),meshnew)

                  ! Find natural neighbors
                  i_src = meshnew%lend(meshnew%nnr)
                  loop = .true.
                  nnat = 0
                  do while (loop)
                     nnat = nnat+1
                     natis(nnat) = abs(meshnew%list(i_src))
                     i_src = meshnew%lptr(i_src)
                     loop = (i_src/=meshnew%lend(meshnew%nnr))
                  end do

                  if (all(mask_src(natis(1:nnat)))) then
                     ! Compute natural neighbors polygons areas
                     call polygon(meshnew,nnat,natis(1:nnat),area_polygon_new(1:nnat))

                     ! Compute weight
                     natwgt(1:nnat) = area_polygon(natis(1:nnat))-area_polygon_new(1:nnat)

                     ! Add interpolation elements
                     do inat=1,nnat
                        n_s = n_s+1
                        row(n_s) = i_dst
                        col(n_s) = mesh%order(natis(inat))
                        S(n_s) = natwgt(inat)/sum(natwgt(1:nnat))
                     end do
                  end if
               end if
            end if
         end if
      else
         ! Subsampled point
         n_s = n_s+1
         row(n_s) = i_dst
         col(n_s) = nn_index(1)
         S(n_s) = 1.0
      end if
   end if

   done(i_dst_loc) = .true.
   call prog_print(progint,done)
end do
write(mpl%unit,'(a)') '100%'

! Communication
call mpl_allgather(1,(/n_s/),n_sg)

! Allocation
interp%n_s = sum(n_sg)
interp%n_src = n_src
interp%n_dst = n_dst
call linop_alloc(interp)

! Communication
if (mpl%main) then
   offset = 0
   do iproc=1,mpl%nproc
      if (n_sg(iproc)>0) then
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
      end if

      ! Update offset
      offset = offset+n_sg(iproc)
   end do
else
   if (n_s>0) then
      ! Send data to ioproc
      call mpl_send(n_s,row(1:n_s),mpl%ioproc,mpl%tag)
      call mpl_send(n_s,col(1:n_s),mpl%ioproc,mpl%tag+1)
      call mpl_send(n_s,S(1:n_s),mpl%ioproc,mpl%tag+2)
   end if
end if
mpl%tag = mpl%tag+3

! Broadcast data
call mpl_bcast(interp%row,mpl%ioproc)
call mpl_bcast(interp%col,mpl%ioproc)
call mpl_bcast(interp%S,mpl%ioproc)

end subroutine compute_interp_from_mesh_ctree

!----------------------------------------------------------------------
! Subroutine: compute_grid_interp
!> Purpose: compute grid horizontal interpolation
!----------------------------------------------------------------------
subroutine compute_grid_interp(geom,nc1,c1_to_c0,mask_check,vbot,vtop,interp_type,h)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom             !< Geometry
integer,intent(in) :: nc1                     !< Subset Sc1 size
integer,intent(in) :: c1_to_c0(nc1)           !< Subset Sc1 to subset Sc0
logical,intent(in) :: mask_check              !< Mask check key
integer,intent(in) :: vbot(nc1)               !< Bottom level
integer,intent(in) :: vtop(nc1)               !< Top level
character(len=1024),intent(in) :: interp_type !< Interpolation type
type(linoptype),intent(inout) :: h(geom%nl0i) !< Horizontal interpolation data

! Local variables
integer :: ic0,ic1,i_s,il0i
real(kind_real) :: renorm(geom%nc0)
logical :: test_c0(geom%nc0)
logical,allocatable :: mask_extra(:),valid(:),missing(:),test_c1(:)
type(linoptype) :: hbase,htmp

! Compute interpolation
call compute_interp(nc1,geom%lon(c1_to_c0),geom%lat(c1_to_c0), any(geom%mask(c1_to_c0,:),dim=2), &
 & geom%nc0,geom%lon,geom%lat,any(geom%mask,dim=2),interp_type,hbase)

! Allocation
allocate(valid(hbase%n_s))
allocate(mask_extra(nc1))
allocate(test_c1(nc1))

do il0i=1,geom%nl0i
   ! Copy
   call linop_copy(hbase,htmp)
   valid = .true.

   ! Check mask boundaries
   if (mask_check) then
      write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Sublevel ',il0i,': '
      call check_mask_bnd(geom,htmp,valid,il0i,col_to_ic0=c1_to_c0)
   else
      write(mpl%unit,'(a10,a,i3)') '','Sublevel ',il0i
   end if

   if (geom%nl0i>1) then
      ! Check extrapolated points because of masks vertical variations (unsampled levels)
      mask_extra = .false.
      do ic1=1,nc1
         ic0 = c1_to_c0(ic1)
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
   else
      mask_extra = .false.
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
   call interp_missing(geom%nc0,geom%lon,geom%lat,geom%mask(:,il0i),interp_type,h(il0i))

   ! Test interpolation
   test_c0 = geom%mask(:,min(il0i,geom%nl0i))
   test_c1 = .true.
   do i_s=1,h(il0i)%n_s
      test_c0(h(il0i)%row(i_s)) = .false.
      test_c1(h(il0i)%col(i_s)) = .false.
   end do
   if (any(test_c0)) call msgerror('error with the grid interpolation row')
   if (any(test_c1)) call msgerror('error with the grid interpolation col')
end do

! Release memory
call linop_dealloc(hbase)
deallocate(valid)
deallocate(mask_extra)
deallocate(test_c1)

end subroutine compute_grid_interp

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
integer :: ic0,i_s,jc0,jc1,progint,il0,iproc
integer :: i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:)

! MPI splitting
call mpl_split(interp%n_s,i_s_s,i_s_e,n_s_loc)

! Allocation
allocate(done(n_s_loc(mpl%myproc)))

! Check that interpolations are not crossing mask boundaries
call prog_init(progint,done)
!$omp parallel do schedule(static) private(i_s_loc,i_s,x,y,z,v1,v2,va,vp,t,ic0,jc1,jc0)
do i_s_loc=1,n_s_loc(mpl%myproc)
   ! Indices
   i_s = i_s_s(mpl%myproc)+i_s_loc-1

   if (valid(i_s)) then
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

      ! Check if arc is crossing boundary arcs
      call check_arc(geom,il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),valid(i_s))
   end if

   ! Print progression
   done(i_s_loc) = .true.
   call prog_print(progint,done)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (n_s_loc(iproc)>0) then
         if (iproc/=mpl%ioproc) then
            ! Receive data on ioproc
            call mpl_recv(n_s_loc(iproc),valid(i_s_s(iproc):i_s_e(iproc)),iproc,mpl%tag)
         end if
      end if
   end do
else
   if (n_s_loc(mpl%myproc)>0) then
      ! Send data to ioproc
      call mpl_send(n_s_loc(mpl%myproc),valid(i_s_s(mpl%myproc):i_s_e(mpl%myproc)),mpl%ioproc,mpl%tag)
   end if
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl_bcast(valid,mpl%ioproc)

! Release memory
deallocate(done)

end subroutine check_mask_bnd


!----------------------------------------------------------------------
! Subroutine: check_arc
!> Purpose: check if an arc is crossing boundaries
!----------------------------------------------------------------------
subroutine check_arc(geom,il0,lon_s,lat_s,lon_e,lat_e,valid)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom   !< Geometry
integer,intent(in) :: il0           !< Level
real(kind_real),intent(in) :: lon_s !< First point longitude
real(kind_real),intent(in) :: lat_s !< First point latitude
real(kind_real),intent(in) :: lon_e !< Second point longitude
real(kind_real),intent(in) :: lat_e !< Second point latitude
logical,intent(out) :: valid        !< True for valid arcs

! Local variables
integer :: ibnd
real(kind_real) :: x(2),y(2),z(2),v1(3),v2(3),va(3),vp(3),t(4)

! Transform to cartesian coordinates
call trans(2,(/lat_s,lat_e/),(/lon_s,lon_e/),x,y,z)

! Compute arc orthogonal vector
v1 = (/x(1),y(1),z(1)/)
v2 = (/x(2),y(2),z(2)/)
call vector_product(v1,v2,va)

! Check if arc is crossing boundary arcs
valid = .true.
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
      valid = .false.
      exit
   end if
end do

end subroutine check_arc

!----------------------------------------------------------------------
! Subroutine: interp_missing
!> Purpose: deal with missing interpolation points
!----------------------------------------------------------------------
subroutine interp_missing(n_dst,lon_dst,lat_dst,mask_dst,interp_type,interp)

implicit none

! Passed variables
integer,intent(in) :: n_dst                   !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst)  !< Destination longitude
real(kind_real),intent(in) :: lat_dst(n_dst)  !< Destination latitude
logical,intent(in) :: mask_dst(n_dst)         !< Destination mask
character(len=1024),intent(in) :: interp_type !< Interpolation type
type(linoptype),intent(inout) :: interp       !< Interpolation data

! Local variables
integer :: i_dst,i_s
integer :: nn(1)
real(kind_real) :: dum(1)
logical :: missing(n_dst),lmask(n_dst),found
type(linoptype) :: interp_tmp
type(ctreetype) :: ctree

! Find missing points
missing = .false.
do i_dst=1,n_dst
   if (mask_dst(i_dst)) missing(i_dst) = .true.
end do
do i_s=1,interp%n_s
   missing(interp%row(i_s)) = .false.
end do

if (count(missing)>0) then
   ! Allocate temporary interpolation
   if (trim(interp_type)=='bilin') then
      interp_tmp%n_s = interp%n_s+3*count(missing)
   elseif (trim(interp_type)=='natural') then
      interp_tmp%n_s = interp%n_s+40*count(missing)
   end if
   call linop_alloc(interp_tmp)

   ! Fill arrays
   interp_tmp%row(1:interp%n_s) = interp%row
   interp_tmp%col(1:interp%n_s) = interp%col
   interp_tmp%S(1:interp%n_s) = interp%S

   ! Reset size
   interp_tmp%n_s = interp%n_s

   ! Mask
   lmask = mask_dst.and.(.not.missing)

   ! Compute cover tree
   ctree = create_ctree(n_dst,dble(lon_dst),dble(lat_dst),lmask)

   do i_dst=1,n_dst
      if (missing(i_dst)) then
         ! Compute nearest neighbor
         call find_nearest_neighbors(ctree,dble(lon_dst(i_dst)),dble(lat_dst(i_dst)),1,nn,dum)

         ! Copy data
         found = .false.
         do i_s=1,interp%n_s
            if (interp%row(i_s)==nn(1)) then
               found = .true.
               interp_tmp%n_s = interp_tmp%n_s+1
               interp_tmp%row(interp_tmp%n_s) = i_dst
               interp_tmp%col(interp_tmp%n_s) = interp%col(i_s)
               interp_tmp%S(interp_tmp%n_s) = interp%S(i_s)
            end if
         end do
         if (.not.found) call msgerror('missing point not found')
      end if
   end do

   ! Reallocate interpolation
   call linop_dealloc(interp)
   interp%n_s = interp_tmp%n_s

   ! Fill arrays
   interp%row = interp_tmp%row(1:interp%n_s)
   interp%col = interp_tmp%col(1:interp%n_s)
   interp%S = interp_tmp%S(1:interp%n_s)

   ! Release memory
   call linop_dealloc(interp_tmp)
   call delete_ctree(ctree)
end if

end subroutine interp_missing

end module tools_interp
