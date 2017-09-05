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
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_stripack, only: trans,trmesh,trfind
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_reorder
use type_mesh, only: meshtype,create_mesh,mesh_dealloc
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_randgen, only: randgentype

implicit none

real(kind_real),parameter :: S_inf = 1.0e-3 !< Minimum value for the interpolation coefficients

interface interp_horiz
  module procedure interp_horiz_from_lat_lon
  module procedure interp_horiz_from_mesh_ctree
end interface

private
public :: interp_horiz

contains

!----------------------------------------------------------------------
! Subroutine: interp_horiz_from_lat_lon
!> Purpose: compute horizontal interpolation from source latitude/longitude
!----------------------------------------------------------------------
subroutine interp_horiz_from_lat_lon(rng,n_src,lon_src,lat_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp)

implicit none

! Passed variables
type(randgentype),intent(in) :: rng
integer,intent(in) :: n_src
real(kind_real),intent(in) :: lon_src(n_src)
real(kind_real),intent(in) :: lat_src(n_src)
logical,intent(in) :: mask_src(n_src)
integer,intent(in) :: n_dst
real(kind_real),intent(in) :: lon_dst(n_dst)
real(kind_real),intent(in) :: lat_dst(n_dst)
logical,intent(in) :: mask_dst(n_dst)
type(linoptype),intent(inout) :: interp

! Local variables
integer,allocatable :: mask_ctree(:)
type(ctreetype) :: ctree
type(meshtype) :: mesh

! Create mesh
call create_mesh(rng,n_src,lon_src,lat_src,.false.,mesh)

! Compute cover tree
allocate(mask_ctree(mesh%nnr))
mask_ctree = 1
ctree = create_ctree(mesh%nnr,dble(lon_src(mesh%order)),dble(lat_src(mesh%order)),mask_ctree)
deallocate(mask_ctree)

! Compute interpolation
call interp_horiz_from_mesh_ctree(mesh,ctree,n_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp)

! Release memory
call delete_ctree(ctree)
call mesh_dealloc(mesh)

end subroutine interp_horiz_from_lat_lon

!----------------------------------------------------------------------
! Subroutine: interp_horiz_from_mesh_ctree
!> Purpose: compute horizontal interpolation from source mesh and ctree
!----------------------------------------------------------------------
subroutine interp_horiz_from_mesh_ctree(mesh,ctree,n_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp)

implicit none

! Passed variables
type(meshtype),intent(in) :: mesh
type(ctreetype),intent(in) :: ctree
integer,intent(in) :: n_src
logical,intent(in) :: mask_src(n_src)
integer,intent(in) :: n_dst
real(kind_real),intent(in) :: lon_dst(n_dst)
real(kind_real),intent(in) :: lat_dst(n_dst)
logical,intent(in) :: mask_dst(n_dst)
type(linoptype),intent(inout) :: interp

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

      if (.not.(abs(dist(1))>0.0)) then
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

end subroutine interp_horiz_from_mesh_ctree

end module tools_interp
