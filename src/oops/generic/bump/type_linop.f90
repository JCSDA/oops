!----------------------------------------------------------------------
! Module: type_linop
!> Purpose: linear operator derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_linop

use netcdf
!$ use omp_lib
use tools_display, only: msgerror,msgwarning,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use tools_qsort, only: qsort
use type_ctree, only: ctree_type
use type_geom, only: geom_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl

implicit none

! Linear operator derived type
type linop_type
   character(len=1024) :: prefix            !< Operator prefix (for I/O)
   integer :: n_src                         !< Source vector size
   integer :: n_dst                         !< Destination vector size
   integer :: n_s                           !< Operator size
   integer,allocatable :: row(:)            !< Output indices
   integer,allocatable :: col(:)            !< Input indices
   real(kind_real),allocatable :: S(:)      !< Coefficients
   integer :: nvec                          !< Size of the vector of linear operators with similar row and col
   real(kind_real),allocatable :: Svec(:,:) !< Coefficients of the vector of linear operators with similar row and col
contains
   procedure :: alloc => linop_alloc
   procedure :: dealloc => linop_dealloc
   procedure :: copy => linop_copy
   procedure :: reorder => linop_reorder
   procedure :: read => linop_read
   procedure :: write => linop_write
   procedure :: apply => linop_apply
   procedure :: apply_ad => linop_apply_ad
   procedure :: apply_sym => linop_apply_sym
   procedure :: add_op => linop_add_op
   procedure :: linop_interp_from_lat_lon
   procedure :: linop_interp_from_mesh_ctree
   procedure :: linop_interp_grid
   generic :: interp => linop_interp_from_lat_lon,linop_interp_from_mesh_ctree,linop_interp_grid
   procedure :: interp_check_mask => linop_interp_check_mask
   procedure :: interp_missing => linop_interp_missing
end type linop_type

logical :: check_data = .false.   !< Activate data check for all linear operations
integer,parameter :: nnatmax = 40 !< Maximum number of natural neighbors
real,parameter :: S_inf = 1.0e-2  !< Minimum interpolation coefficient

private
public :: linop_type

contains

!----------------------------------------------------------------------
! Subroutine: linop_alloc
!> Purpose: linear operator object allocation
!----------------------------------------------------------------------
subroutine linop_alloc(linop,nvec)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop !< Linear operator
integer,intent(in),optional :: nvec      !< Size of the vector of linear operators with similar row and col

! Allocation
allocate(linop%row(linop%n_s))
allocate(linop%col(linop%n_s))
if (present(nvec)) then
   linop%nvec = nvec
   allocate(linop%Svec(linop%n_s,linop%nvec))
else
   allocate(linop%S(linop%n_s))
end if

! Initialization
call msi(linop%row)
call msi(linop%col)
if (present(nvec)) then
   call msr(linop%Svec)
else
   call msi(linop%nvec)
   call msr(linop%S)
end if

end subroutine linop_alloc

!----------------------------------------------------------------------
! Subroutine: linop_dealloc
!> Purpose: linear operator object deallocation
!----------------------------------------------------------------------
subroutine linop_dealloc(linop)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop !< Linear operator

! Release memory
if (allocated(linop%row)) deallocate(linop%row)
if (allocated(linop%col)) deallocate(linop%col)
if (allocated(linop%S)) deallocate(linop%S)
if (allocated(linop%Svec)) deallocate(linop%Svec)

end subroutine linop_dealloc

!----------------------------------------------------------------------
! Function: linop_copy
!> Purpose: linear operator object copy
!----------------------------------------------------------------------
type(linop_type) function linop_copy(linop)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop !< Linear operator

! Copy attributes
linop_copy%prefix = trim(linop%prefix)
linop_copy%n_src = linop%n_src
linop_copy%n_dst = linop%n_dst
linop_copy%n_s = linop%n_s

! Deallocation
call linop_copy%dealloc

! Allocation
if (isnotmsi(linop%nvec)) then
   call linop_copy%alloc(linop%nvec)
else
   call linop_copy%alloc
end if

! Copy data
linop_copy%row = linop%row
linop_copy%col = linop%col
if (isnotmsi(linop_copy%nvec)) then
   linop_copy%Svec = linop%Svec
else
   linop_copy%S = linop%S
end if

end function linop_copy

!----------------------------------------------------------------------
! Subroutine: linop_reorder
!> Purpose: reorder linear operator
!----------------------------------------------------------------------
subroutine linop_reorder(linop)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop !< Linear operator

! Local variables
integer :: row,i_s_s,i_s_e,n_s,i_s
integer,allocatable :: order(:)

if (linop%n_s<1000000) then
   ! Sort with respect to row
   allocate(order(linop%n_s))
   call qsort(linop%n_s,linop%row,order)

   ! Sort col and S
   linop%col = linop%col(order)
   if (isnotmsi(linop%nvec)) then
      linop%Svec = linop%Svec(order,:)
   else
      linop%S = linop%S(order)
   end if
   deallocate(order)

   ! Sort with respect to col for each row
   row = minval(linop%row)
   i_s_s = 1
   call msi(i_s_e)
   do i_s=1,linop%n_s
      if (linop%row(i_s)==row) then
         i_s_e = i_s
      else
         n_s = i_s_e-i_s_s+1
         allocate(order(n_s))
         call qsort(n_s,linop%col(i_s_s:i_s_e),order)
         order = order+i_s_s-1
         if (isnotmsi(linop%nvec)) then
            linop%Svec(i_s_s:i_s_e,:) = linop%Svec(order,:)
         else
            linop%S(i_s_s:i_s_e) = linop%S(order)
         end if
         deallocate(order)
         i_s_s = i_s+1
         row = linop%row(i_s)
      end if
   end do
else
   call msgwarning('linear operator is too big, impossible to reorder')
end if

end subroutine linop_reorder


!----------------------------------------------------------------------
! Subroutine: linop_read
!> Purpose: read linear operator from a NetCDF file
!----------------------------------------------------------------------
subroutine linop_read(linop,ncid)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop !< Linear operator
integer,intent(in) :: ncid               !< NetCDF file ID

! Local variables
integer :: info,nvec
integer :: n_s_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_read'

! Get operator size
info = nf90_inq_dimid(ncid,trim(linop%prefix)//'_n_s',n_s_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n_s_id,len=linop%n_s))
else
   linop%n_s = 0
end if

if (linop%n_s>0) then
   ! Get source/destination dimensions
   call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(linop%prefix)//'_n_src',linop%n_src))
   call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(linop%prefix)//'_n_dst',linop%n_dst))
   call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(linop%prefix)//'_nvec',nvec))

   ! Allocation
   if (isnotmsi(nvec)) then
      call linop%alloc(nvec)
   else
      call linop%alloc
   end if

   ! Get variables id
   call ncerr(subr,nf90_inq_varid(ncid,trim(linop%prefix)//'_row',row_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(linop%prefix)//'_col',col_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(linop%prefix)//'_S',S_id))

   ! Get variables
   call ncerr(subr,nf90_get_var(ncid,row_id,linop%row))
   call ncerr(subr,nf90_get_var(ncid,col_id,linop%col))
   if (isnotmsi(linop%nvec)) then
      call ncerr(subr,nf90_get_var(ncid,S_id,linop%Svec))
   else
      call ncerr(subr,nf90_get_var(ncid,S_id,linop%S))
   end if
end if

end subroutine linop_read

!----------------------------------------------------------------------
! Subroutine: linop_write
!> Purpose: write linear operator to a NetCDF file
!----------------------------------------------------------------------
subroutine linop_write(linop,ncid)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop !< Linear operator
integer,intent(in) :: ncid            !< NetCDF file ID

! Local variables
integer :: n_s_id,nvec_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_write'

if (linop%n_s>0) then
   ! Start definition mode
   call ncerr(subr,nf90_redef(ncid))

   ! Write source/destination dimensions
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_src',linop%n_src))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_dst',linop%n_dst))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_nvec',linop%nvec))

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,trim(linop%prefix)//'_n_s',linop%n_s,n_s_id))
   if (isnotmsi(linop%nvec)) call ncerr(subr,nf90_def_dim(ncid,trim(linop%prefix)//'_nvec',linop%nvec,nvec_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_row',nf90_int,(/n_s_id/),row_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_col',nf90_int,(/n_s_id/),col_id))
   if (isnotmsi(linop%nvec)) then
      call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_S',ncfloat,(/n_s_id,nvec_id/),S_id))
   else
      call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_S',ncfloat,(/n_s_id/),S_id))
   end if

   ! End definition mode
   call ncerr(subr,nf90_enddef(ncid))

   ! Put variables
   call ncerr(subr,nf90_put_var(ncid,row_id,linop%row))
   call ncerr(subr,nf90_put_var(ncid,col_id,linop%col))
   if (isnotmsi(linop%nvec)) then
      call ncerr(subr,nf90_put_var(ncid,S_id,linop%Svec))
   else
      call ncerr(subr,nf90_put_var(ncid,S_id,linop%S))
   end if
end if

end subroutine linop_write

!----------------------------------------------------------------------
! Subroutine: linop_apply
!> Purpose: apply linear operator
!----------------------------------------------------------------------
subroutine linop_apply(linop,fld_src,fld_dst,ivec)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop               !< Linear operator
real(kind_real),intent(in) :: fld_src(linop%n_src)  !< Source vector
real(kind_real),intent(out) :: fld_dst(linop%n_dst) !< Destination vector
integer,intent(in),optional :: ivec                 !< Index of the vector of linear operators with similar row and col

! Local variables
integer :: i_s,i_dst
logical :: missing(linop%n_dst)

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call msgerror('row>n_dst for linear operation '//trim(linop%prefix))
   if (present(ivec)) then
      if (any(isnan(linop%Svec))) call msgerror('NaN in Svec for linear operation '//trim(linop%prefix))
   else
      if (any(isnan(linop%S))) call msgerror('NaN in S for linear operation '//trim(linop%prefix))
   end if

   ! Check input
   if (any(fld_src>huge(1.0))) call msgerror('Overflowing number in fld_src for linear operation '//trim(linop%prefix))
   if (any(isnan(fld_src))) call msgerror('NaN in fld_src for linear operation '//trim(linop%prefix))
   if (any(.not.isnotmsr(fld_src))) call msgerror('Missing value in fld_src for linear operation '//trim(linop%prefix))
end if

! Initialization
fld_dst = 0.0
missing = .true.

! Apply weights
do i_s=1,linop%n_s
   if (present(ivec)) then
      fld_dst(linop%row(i_s)) = fld_dst(linop%row(i_s))+linop%Svec(i_s,ivec)*fld_src(linop%col(i_s))
   else
      fld_dst(linop%row(i_s)) = fld_dst(linop%row(i_s))+linop%S(i_s)*fld_src(linop%col(i_s))
   end if
   missing(linop%row(i_s)) = .false.
end do

! Missing destination values
do i_dst=1,linop%n_dst
   if (missing(i_dst)) call msr(fld_dst(i_dst))
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld_dst))) call msgerror('NaN in fld_dst for linear operation '//trim(linop%prefix))
end if

end subroutine linop_apply

!----------------------------------------------------------------------
! Subroutine: linop_apply_ad
!> Purpose: apply linear operator, adjoint
!----------------------------------------------------------------------
subroutine linop_apply_ad(linop,fld_dst,fld_src,ivec)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop               !< Linear operator
real(kind_real),intent(in) :: fld_dst(linop%n_dst)  !< Destination vector
real(kind_real),intent(out) :: fld_src(linop%n_src) !< Source vector
integer,intent(in),optional :: ivec                 !< Index of the vector of linear operators with similar row and col

! Local variables
integer :: i_s

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for adjoint linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call msgerror('row>n_dst for adjoint linear operation '//trim(linop%prefix))
   if (present(ivec)) then
      if (any(isnan(linop%Svec))) call msgerror('NaN in Svec for adjoint linear operation '//trim(linop%prefix))
   else
      if (any(isnan(linop%S))) call msgerror('NaN in S for adjoint linear operation '//trim(linop%prefix))
   end if

   ! Check input
   if (any(fld_dst>huge(1.0))) call msgerror('Overflowing number in fld_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(isnan(fld_dst))) call msgerror('NaN in fld_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(.not.isnotmsr(fld_dst))) call msgerror('Missing value in fld_dst for adjoint linear operation '//trim(linop%prefix))
end if

! Initialization
fld_src = 0.0

! Apply weights
do i_s=1,linop%n_s
   if (present(ivec)) then
      fld_src(linop%col(i_s)) = fld_src(linop%col(i_s))+linop%Svec(i_s,ivec)*fld_dst(linop%row(i_s))
   else
      fld_src(linop%col(i_s)) = fld_src(linop%col(i_s))+linop%S(i_s)*fld_dst(linop%row(i_s))
   end if
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld_src))) call msgerror('NaN in fld_src for adjoint linear operation '//trim(linop%prefix))
end if

end subroutine linop_apply_ad

!----------------------------------------------------------------------
! Subroutine: linop_apply_sym
!> Purpose: apply linear operator, symmetric
!----------------------------------------------------------------------
subroutine linop_apply_sym(linop,fld,ivec)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop             !< Linear operator
real(kind_real),intent(inout) :: fld(linop%n_src) !< Source/destination vector
integer,intent(in),optional :: ivec                 !< Index of the vector of linear operators with similar row and col

! Local variables
integer :: i_s,ithread
real(kind_real) :: fld_arr(linop%n_dst,mpl%nthread)

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for symmetric linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_src) call msgerror('row>n_dst for symmetric linear operation '//trim(linop%prefix))
   if (present(ivec)) then
      if (any(isnan(linop%Svec))) call msgerror('NaN in Svec for symmetric linear operation '//trim(linop%prefix))
   else
      if (any(isnan(linop%S))) call msgerror('NaN in S for symmetric linear operation '//trim(linop%prefix))
   end if

   ! Check input
   if (any(fld>huge(1.0))) call msgerror('Overflowing number in fld for symmetric linear operation '//trim(linop%prefix))
   if (any(isnan(fld))) call msgerror('NaN in fld for symmetric linear operation '//trim(linop%prefix))
   if (any(.not.isnotmsr(fld))) call msgerror('Missing value in fld for symmetric linear operation '//trim(linop%prefix))
end if

! Apply weights
fld_arr = 0.0
!$omp parallel do schedule(static) private(i_s,ithread)
do i_s=1,linop%n_s
   ithread = omp_get_thread_num()+1
   if (present(ivec)) then
      fld_arr(linop%row(i_s),ithread) = fld_arr(linop%row(i_s),ithread)+linop%Svec(i_s,ivec)*fld(linop%col(i_s))
      fld_arr(linop%col(i_s),ithread) = fld_arr(linop%col(i_s),ithread)+linop%Svec(i_s,ivec)*fld(linop%row(i_s))
   else
      fld_arr(linop%row(i_s),ithread) = fld_arr(linop%row(i_s),ithread)+linop%S(i_s)*fld(linop%col(i_s))
      fld_arr(linop%col(i_s),ithread) = fld_arr(linop%col(i_s),ithread)+linop%S(i_s)*fld(linop%row(i_s))
   end if
end do
!$omp end parallel do

! Sum over threads
do ithread=1,mpl%nthread
   fld = fld+fld_arr(:,ithread)
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld))) call msgerror('NaN in fld for symmetric linear operation '//trim(linop%prefix))
end if

end subroutine linop_apply_sym

!----------------------------------------------------------------------
! Subroutine: linop_add_op
!> Purpose: add operation
!----------------------------------------------------------------------
subroutine linop_add_op(linop,n_s,row,col,S)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop !< Linear operator
integer,intent(inout) :: n_s             !< Number of operations
integer,intent(in) :: row                !< Row index
integer,intent(in) :: col                !< Column index
real(kind_real),intent(in) :: S          !< Value

! Local variables
type(linop_type) :: linop_tmp

! Update
n_s = n_s+1
if (n_s>linop%n_s) then
   ! Copy
   linop_tmp = linop%copy()

   ! Reallocate larger linear operation
   call linop%dealloc
   linop%n_s = 2*linop_tmp%n_s
   call linop%alloc

   ! Copy data
   linop%row(1:linop_tmp%n_s) = linop_tmp%row
   linop%col(1:linop_tmp%n_s) = linop_tmp%col
   linop%S(1:linop_tmp%n_s) = linop_tmp%S

   ! Release memory
   call linop_tmp%dealloc
end if

! New operation
linop%row(n_s) = row
linop%col(n_s) = col
linop%S(n_s) = S

end subroutine linop_add_op

!----------------------------------------------------------------------
! Subroutine: linop_interp_from_lat_lon
!> Purpose: compute horizontal interpolation from source latitude/longitude
!----------------------------------------------------------------------
subroutine linop_interp_from_lat_lon(linop,n_src,lon_src,lat_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop     !< Linear operator
integer,intent(in) :: n_src                  !< Source size
real(kind_real),intent(in) :: lon_src(n_src) !< Source longitudes
real(kind_real),intent(in) :: lat_src(n_src) !< Source latitudes
logical,intent(in) :: mask_src(n_src)        !< Source mask
integer,intent(in) :: n_dst                  !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) !< Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst) !< Destination latitudes
logical,intent(in) :: mask_dst(n_dst)        !< Destination mask
character(len=*),intent(in) :: interp_type   !< Interpolation type

! Local variables
integer :: n_src_eff,i_src,i_src_eff
integer,allocatable :: src_eff_to_src(:)
logical,allocatable :: mask_ctree(:),mask_src_eff(:)
type(ctree_type) :: ctree
type(mesh_type) :: mesh

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
call mesh%create(n_src_eff,lon_src(src_eff_to_src),lat_src(src_eff_to_src))

! Compute cover tree
allocate(mask_ctree(n_src_eff))
mask_ctree = .true.
call ctree%create(n_src_eff,lon_src(src_eff_to_src),lat_src(src_eff_to_src),mask_ctree)
deallocate(mask_ctree)

! Compute interpolation
mask_src_eff = .true.
call linop%interp(mesh,ctree,n_src_eff,mask_src_eff,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

! Effective points conversion
linop%n_src = n_src
linop%col = src_eff_to_src(linop%col)

! Release memory
call ctree%delete
call mesh%dealloc

end subroutine linop_interp_from_lat_lon

!----------------------------------------------------------------------
! Subroutine: linop_interp_from_mesh_ctree
!> Purpose: compute horizontal interpolation from source mesh and ctree
!----------------------------------------------------------------------
subroutine linop_interp_from_mesh_ctree(linop,mesh,ctree,n_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop     !< Linear operator
type(mesh_type),intent(in) :: mesh           !< Mesh
type(ctree_type),intent(in) :: ctree         !< Cover tree
integer,intent(in) :: n_src                  !< Source size
logical,intent(in) :: mask_src(n_src)        !< Source mask
integer,intent(in) :: n_dst                  !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) !< Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst) !< Destination latitudes
logical,intent(in) :: mask_dst(n_dst)        !< Destination mask
character(len=*),intent(in) :: interp_type   !< Interpolation type

! Local variables
integer :: i,i_src,i_dst,nn_index(1),n_s,progint,ib(3),nnat,inat,np,iproc,offset,i_s
integer :: i_dst_s(mpl%nproc),i_dst_e(mpl%nproc),n_dst_loc(mpl%nproc),i_dst_loc,proc_to_n_s(mpl%nproc)
integer,allocatable :: natis(:),row(:),col(:)
real(kind_real) :: nn_dist(1),b(3)
real(kind_real),allocatable :: area_polygon(:),area_polygon_new(:),natwgt(:),S(:)
logical :: loop
logical,allocatable :: done(:),missing(:)
type(mesh_type) :: meshnew

! MPI splitting
call mpl%split(n_dst,i_dst_s,i_dst_e,n_dst_loc)

! Allocation
call msi(np)
if (trim(interp_type)=='bilin') then
   ! Bilinear interpolation
   np = 3
elseif (trim(interp_type)=='natural') then
   ! Natural neighbors
   np = nnatmax
   allocate(area_polygon(mesh%n))
   allocate(area_polygon_new(nnatmax))
   allocate(natis(mesh%n))
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
   do i_src=1,mesh%n
      natis(i_src) = i_src
   end do
   call mesh%polygon(mesh%n,natis,area_polygon)
end if

! Compute interpolation
write(mpl%unit,'(a10,a)',advance='no') '','Compute interpolation: '
call flush(mpl%unit)
call prog_init(progint,done)
n_s = 0
do i_dst_loc=1,n_dst_loc(mpl%myproc)
   ! Indices
   i_dst = i_dst_s(mpl%myproc)+i_dst_loc-1

   if (mask_dst(i_dst)) then
      ! Find nearest neighbor
      call ctree%find_nearest_neighbors(lon_dst(i_dst),lat_dst(i_dst),1,nn_index,nn_dist)

      if (abs(nn_dist(1))>0.0) then
         ! Compute barycentric coordinates
         call mesh%barycentric(lon_dst(i_dst),lat_dst(i_dst),nn_index(1),b,ib)
         if (sum(b)>0.0) b = b/sum(b)
         do i=1,3
            b(i) = nint(b(i)/S_inf)*S_inf
         end do
         if (sum(b)>0.0) b = b/sum(b)

         if (all(ib>0)) then
            if (all(mask_src(mesh%order(ib)))) then
               ! Valid interpolation
               if (trim(interp_type)=='bilin') then
                  ! Bilinear interpolation
                  do i=1,3
                     if (b(i)>0.0) then
                        n_s = n_s+1
                        row(n_s) = i_dst
                        col(n_s) = mesh%order(ib(i))
                        S(n_s) = b(i)
                     end if
                  end do
               elseif (trim(interp_type)=='natural') then
                  ! Natural neighbors interpolation

                  ! Copy mesh
                  meshnew = mesh%copy()

                  ! Add a node
                  call meshnew%addnode(lon_dst(i_dst),lat_dst(i_dst))

                  ! Find natural neighbors
                  i_src = meshnew%lend(meshnew%n)
                  loop = .true.
                  nnat = 0
                  do while (loop)
                     nnat = nnat+1
                     natis(nnat) = abs(meshnew%list(i_src))
                     i_src = meshnew%lptr(i_src)
                     loop = (i_src/=meshnew%lend(meshnew%n))
                  end do

                  if (all(mask_src(natis(1:nnat)))) then
                     ! Compute natural neighbors polygons areas
                     call meshnew%polygon(nnat,natis(1:nnat),area_polygon_new(1:nnat))

                     ! Compute weight
                     natwgt = 0.0
                     natwgt(1:nnat) = area_polygon(natis(1:nnat))-area_polygon_new(1:nnat)
                     if (sum(natwgt(1:nnat))>0.0) natwgt(1:nnat) = natwgt(1:nnat)/sum(natwgt(1:nnat))
                     do inat=1,nnat
                        natwgt(inat) = nint(natwgt(inat)/S_inf)*S_inf
                     end do
                     if (sum(natwgt(1:nnat))>0.0) natwgt(1:nnat) = natwgt(1:nnat)/sum(natwgt(1:nnat))

                     ! Add interpolation element
                     do inat=1,nnat
                        if (natwgt(inat)>0.0) then
                           n_s = n_s+1
                           row(n_s) = i_dst
                           col(n_s) = mesh%order(natis(inat))
                           S(n_s) = natwgt(inat)
                        end if
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
call flush(mpl%unit)

! Communication
call mpl%allgather(1,(/n_s/),proc_to_n_s)

! Allocation
linop%n_s = sum(proc_to_n_s)
linop%n_src = n_src
linop%n_dst = n_dst
call linop%alloc

! Communication
if (mpl%main) then
   offset = 0
   do iproc=1,mpl%nproc
      if (proc_to_n_s(iproc)>0) then
         if (iproc==mpl%ioproc) then
            ! Copy data
            linop%row(offset+1:offset+proc_to_n_s(iproc)) = row(1:proc_to_n_s(iproc))
            linop%col(offset+1:offset+proc_to_n_s(iproc)) = col(1:proc_to_n_s(iproc))
            linop%S(offset+1:offset+proc_to_n_s(iproc)) = S(1:proc_to_n_s(iproc))
         else
            ! Receive data on ioproc
            call mpl%recv(proc_to_n_s(iproc),linop%row(offset+1:offset+proc_to_n_s(iproc)),iproc,mpl%tag)
            call mpl%recv(proc_to_n_s(iproc),linop%col(offset+1:offset+proc_to_n_s(iproc)),iproc,mpl%tag+1)
            call mpl%recv(proc_to_n_s(iproc),linop%S(offset+1:offset+proc_to_n_s(iproc)),iproc,mpl%tag+2)
         end if
      end if

      ! Update offset
      offset = offset+proc_to_n_s(iproc)
   end do
else
   if (n_s>0) then
      ! Send data to ioproc
      call mpl%send(n_s,row(1:n_s),mpl%ioproc,mpl%tag)
      call mpl%send(n_s,col(1:n_s),mpl%ioproc,mpl%tag+1)
      call mpl%send(n_s,S(1:n_s),mpl%ioproc,mpl%tag+2)
   end if
end if
mpl%tag = mpl%tag+3

! Broadcast data
call mpl%bcast(linop%row)
call mpl%bcast(linop%col)
call mpl%bcast(linop%S)

! Deal with missing points
call linop%interp_missing(n_dst,lon_dst,lat_dst,mask_dst,interp_type)

! Check interpolation
allocate(missing(n_dst))
missing = .false.
do i_dst=1,n_dst
   if (mask_dst(i_dst)) missing(i_dst) = .true.
end do
do i_s=1,linop%n_s
   missing(linop%row(i_s)) = .false.
end do
if (any(missing)) call msgerror('missing destination points in interp_from_mesh_ctree')

end subroutine linop_interp_from_mesh_ctree

!----------------------------------------------------------------------
! Subroutine: linop_interp_grid
!> Purpose: compute horizontal grid interpolation
!----------------------------------------------------------------------
subroutine linop_interp_grid(linop,geom,il0i,nc1,c1_to_c0,mask_check,vbot,vtop,interp_type,interp_base)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop      !< Linear operator
type(geom_type),intent(in) :: geom            !< Geometry
integer,intent(in) :: il0i                    !< Level
integer,intent(in) :: nc1                     !< Subset Sc1 size
integer,intent(in) :: c1_to_c0(nc1)           !< Subset Sc1 to subset Sc0
logical,intent(in) :: mask_check              !< Mask check key
integer,intent(in) :: vbot(nc1)               !< Bottom level
integer,intent(in) :: vtop(nc1)               !< Top level
character(len=*),intent(in) :: interp_type    !< Interpolation type
type(linop_type),intent(inout) :: interp_base !< Linear operator (base interpolation)

! Local variables
integer :: ic0,ic1,i_s
real(kind_real) :: renorm(geom%nc0)
logical :: test_c0(geom%nc0)
logical,allocatable :: mask_extra(:),valid(:)

if (.not.allocated(interp_base%row)) then
   ! Compute base interpolation
   call interp_base%interp(nc1,geom%lon(c1_to_c0),geom%lat(c1_to_c0), any(geom%mask(c1_to_c0,:),dim=2), &
 & geom%nc0,geom%lon,geom%lat,any(geom%mask,dim=2),interp_type)
end if

! Allocation
allocate(valid(interp_base%n_s))
allocate(mask_extra(nc1))

! Initialization
valid = .true.

! Check mask boundaries
if (mask_check) then
   write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Sublevel ',il0i,': '
   call flush(mpl%unit)
   call interp_base%interp_check_mask(geom,valid,il0i,col_to_ic0=c1_to_c0)
else
   write(mpl%unit,'(a10,a,i3)') '','Sublevel ',il0i
   call flush(mpl%unit)
end if

if (geom%nl0i>1) then
   ! Check extrapolated points because of masks vertical variations (unsampled levels)
   mask_extra = .false.
   do ic1=1,nc1
      ic0 = c1_to_c0(ic1)
      if (geom%mask(ic0,il0i).and.((il0i<vbot(ic1)).or.(il0i>vtop(ic1)))) mask_extra(ic1) = .true.
   end do

   ! Remove operations for extrapolated points
   do i_s=1,interp_base%n_s
      if (valid(i_s)) then
         ic1 = interp_base%col(i_s)
         if (mask_extra(ic1)) valid(i_s) = .false.
      end if
   end do
   if (count(mask_extra)>0) then
      write(mpl%unit,'(a10,a,i5)') '','Extrapolated points: ',count(mask_extra)
      call flush(mpl%unit)
   end if
else
   mask_extra = .false.
end if

! Renormalization
renorm = 0.0
do i_s=1,interp_base%n_s
   if (valid(i_s)) renorm(interp_base%row(i_s)) = renorm(interp_base%row(i_s))+interp_base%S(i_s)
end do

! Initialize object
linop%n_src = nc1
linop%n_dst = geom%nc0
linop%n_s = count(valid)
call linop%alloc
linop%n_s = 0
do i_s=1,interp_base%n_s
   if (valid(i_s)) then
      linop%n_s = linop%n_s+1
      linop%row(linop%n_s) = interp_base%row(i_s)
      linop%col(linop%n_s) = interp_base%col(i_s)
      linop%S(linop%n_s) = interp_base%S(i_s)/renorm(interp_base%row(i_s))
   end if
end do

! Release memory
call interp_base%dealloc

! Deal with missing points
call linop%interp_missing(geom%nc0,geom%lon,geom%lat,geom%mask(:,il0i),interp_type)

! Check interpolation
test_c0 = geom%mask(:,min(il0i,geom%nl0i))
do i_s=1,linop%n_s
   test_c0(linop%row(i_s)) = .false.
end do
if (any(test_c0)) call msgerror('error with the grid interpolation row')

! Release memory
deallocate(valid)
deallocate(mask_extra)

end subroutine linop_interp_grid

!----------------------------------------------------------------------
! Subroutine: interp_check_mask
!> Purpose: check mask boundaries for interpolations
!----------------------------------------------------------------------
subroutine linop_interp_check_mask(linop,geom,valid,il0,row_to_ic0,col_to_ic0)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop               !< Linear operator
type(geom_type),intent(in) :: geom                     !< Geometry
logical,intent(inout) :: valid(linop%n_s)              !< Valid points
integer,intent(in),optional :: row_to_ic0(linop%n_dst) !< Conversion from row to ic0 (identity if missing)
integer,intent(in),optional :: col_to_ic0(linop%n_src) !< Conversion from col to ic0 (identity if missing)

! Local variables
integer :: ic0,i_s,jc0,jc1,progint,il0,iproc
integer :: i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:)

! MPI splitting
call mpl%split(linop%n_s,i_s_s,i_s_e,n_s_loc)

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
         ic0 = row_to_ic0(linop%row(i_s))
      else
         ic0 = linop%row(i_s)
      end if
      if (present(col_to_ic0)) then
         jc0 = col_to_ic0(linop%col(i_s))
      else
         jc0 = linop%col(i_s)
      end if

      ! Check if arc is crossing boundary arcs
      call geom%check_arc(il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),valid(i_s))
   end if

   ! Print progression
   done(i_s_loc) = .true.
   call prog_print(progint,done)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'
call flush(mpl%unit)

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (n_s_loc(iproc)>0) then
         if (iproc/=mpl%ioproc) then
            ! Receive data on ioproc
            call mpl%recv(n_s_loc(iproc),valid(i_s_s(iproc):i_s_e(iproc)),iproc,mpl%tag)
         end if
      end if
   end do
else
   if (n_s_loc(mpl%myproc)>0) then
      ! Send data to ioproc
      call mpl%send(n_s_loc(mpl%myproc),valid(i_s_s(mpl%myproc):i_s_e(mpl%myproc)),mpl%ioproc,mpl%tag)
   end if
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl%bcast(valid)

! Release memory
deallocate(done)

end subroutine linop_interp_check_mask

!----------------------------------------------------------------------
! Subroutine: linop_interp_missing
!> Purpose: deal with missing interpolation points
!----------------------------------------------------------------------
subroutine linop_interp_missing(linop,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop     !< Linear operator (interpolation)
integer,intent(in) :: n_dst                  !< Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) !< Destination longitude
real(kind_real),intent(in) :: lat_dst(n_dst) !< Destination latitude
logical,intent(in) :: mask_dst(n_dst)        !< Destination mask
character(len=*),intent(in) :: interp_type   !< Interpolation type

! Local variables
integer :: i_dst,i_s
integer :: nn(1)
real(kind_real) :: dum(1)
logical :: missing(n_dst),lmask(n_dst),found
type(linop_type) :: interp_tmp
type(ctree_type) :: ctree

! Find missing points
missing = .false.
do i_dst=1,n_dst
   if (mask_dst(i_dst)) missing(i_dst) = .true.
end do
do i_s=1,linop%n_s
   missing(linop%row(i_s)) = .false.
end do

if (count(missing)>0) then
   ! Allocate temporary interpolation
   if (trim(interp_type)=='bilin') then
      interp_tmp%n_s = linop%n_s+3*count(missing)
   elseif (trim(interp_type)=='natural') then
      interp_tmp%n_s = linop%n_s+40*count(missing)
   else
      call msgerror('wrong interpolation')
   end if
   call interp_tmp%alloc

   ! Fill arrays
   interp_tmp%row(1:linop%n_s) = linop%row
   interp_tmp%col(1:linop%n_s) = linop%col
   interp_tmp%S(1:linop%n_s) = linop%S

   ! Reset size
   interp_tmp%n_s = linop%n_s

   ! Mask
   lmask = mask_dst.and.(.not.missing)

   ! Compute cover tree
   call ctree%create(n_dst,lon_dst,lat_dst,lmask)

   do i_dst=1,n_dst
      if (missing(i_dst)) then
         ! Compute nearest neighbor
         call ctree%find_nearest_neighbors(lon_dst(i_dst),lat_dst(i_dst),1,nn,dum)

         ! Copy data
         found = .false.
         do i_s=1,linop%n_s
            if (linop%row(i_s)==nn(1)) then
               found = .true.
               interp_tmp%n_s = interp_tmp%n_s+1
               interp_tmp%row(interp_tmp%n_s) = i_dst
               interp_tmp%col(interp_tmp%n_s) = linop%col(i_s)
               interp_tmp%S(interp_tmp%n_s) = linop%S(i_s)
            end if
         end do
         if (.not.found) call msgerror('missing point not found')
      end if
   end do

   ! Reallocate interpolation
   call linop%dealloc
   linop%n_s = interp_tmp%n_s
   call linop%alloc

   ! Fill arrays
   linop%row = interp_tmp%row(1:linop%n_s)
   linop%col = interp_tmp%col(1:linop%n_s)
   linop%S = interp_tmp%S(1:linop%n_s)

   ! Release memory
   call interp_tmp%dealloc
   call ctree%delete
end if

end subroutine linop_interp_missing

end module type_linop
