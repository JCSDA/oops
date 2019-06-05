!----------------------------------------------------------------------
! Module: type_linop
! Purpose: linear operator derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_linop

use netcdf
!$ use omp_lib
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: inf
use type_geom, only: geom_type
use type_tree, only: tree_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_rng, only: rng_type

implicit none

logical,parameter :: check_data = .false.             ! Activate data check for all linear operations
integer,parameter :: reorder_max = 1000000            ! Maximum size of linear operation to allow reordering
integer,parameter :: nnatmax = 40                     ! Maximum number of natural neighbors
real(kind_real),parameter :: S_inf = 1.0e-2_kind_real ! Minimum interpolation coefficient

! Linear operator derived type
type linop_type
   character(len=1024) :: prefix            ! Operator prefix (for I/O)
   integer :: n_src                         ! Source vector size
   integer :: n_dst                         ! Destination vector size
   integer :: n_s                           ! Operator size
   integer,allocatable :: row(:)            ! Output indices
   integer,allocatable :: col(:)            ! Input indices
   real(kind_real),allocatable :: S(:)      ! Coefficients
   integer :: nvec                          ! Size of the vector of linear operators with similar row and col
   real(kind_real),allocatable :: Svec(:,:) ! Coefficients of the vector of linear operators with similar row and col
contains
   procedure :: alloc => linop_alloc
   procedure :: dealloc => linop_dealloc
   procedure :: copy => linop_copy
   procedure :: read => linop_read
   procedure :: write => linop_write
   procedure :: reorder => linop_reorder
   procedure :: apply => linop_apply
   procedure :: apply_ad => linop_apply_ad
   procedure :: apply_sym => linop_apply_sym
   procedure :: add_op => linop_add_op
   procedure :: gather => linop_gather
   procedure :: linop_interp_from_lat_lon
   procedure :: linop_interp_from_mesh_tree
   procedure :: linop_interp_grid
   generic :: interp => linop_interp_from_lat_lon,linop_interp_from_mesh_tree,linop_interp_grid
   procedure :: check_mask => linop_check_mask
   procedure :: interp_missing => linop_interp_missing
end type linop_type

private
public :: linop_type

contains

!----------------------------------------------------------------------
! Subroutine: linop_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine linop_alloc(linop,nvec)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop ! Linear operator
integer,intent(in),optional :: nvec      ! Size of the vector of linear operators with similar row and col

! Vector size
if (present(nvec)) then
   linop%nvec = nvec
else
   linop%nvec = 0
end if

! Allocation
allocate(linop%row(linop%n_s))
allocate(linop%col(linop%n_s))
if (linop%nvec>0) then
   allocate(linop%Svec(linop%n_s,linop%nvec))
else
   allocate(linop%S(linop%n_s))
end if

end subroutine linop_alloc

!----------------------------------------------------------------------
! Subroutine: linop_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine linop_dealloc(linop)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop ! Linear operator

! Release memory
if (allocated(linop%row)) deallocate(linop%row)
if (allocated(linop%col)) deallocate(linop%col)
if (allocated(linop%S)) deallocate(linop%S)
if (allocated(linop%Svec)) deallocate(linop%Svec)

end subroutine linop_dealloc

!----------------------------------------------------------------------
! Subroutine: linop_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine linop_copy(linop_out,linop_in)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop_out ! Output linear operator
type(linop_type),intent(in) :: linop_in      ! Input linear operator

! Release memory
call linop_out%dealloc

! Copy attributes
linop_out%prefix = linop_in%prefix
linop_out%n_src = linop_in%n_src
linop_out%n_dst = linop_in%n_dst
linop_out%n_s = linop_in%n_s

! Allocation
call linop_out%alloc(linop_in%nvec)

! Copy data
if (linop_in%n_s>0) then
   linop_out%row = linop_in%row
   linop_out%col = linop_in%col
   if (linop_out%nvec>0) then
      linop_out%Svec = linop_in%Svec
   else
      linop_out%S = linop_in%S
   end if
end if

end subroutine linop_copy

!----------------------------------------------------------------------
! Subroutine: linop_read
! Purpose: read
!----------------------------------------------------------------------
subroutine linop_read(linop,mpl,ncid)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop ! Linear operator
type(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: ncid               ! NetCDF file ID

! Local variables
integer :: info,nvec
integer :: n_s_id,row_id,col_id,S_id,Svec_id
character(len=1024),parameter :: subr = 'linop_read'

! Get operator size
info = nf90_inq_dimid(ncid,trim(linop%prefix)//'_n_s',n_s_id)
if (info==nf90_noerr) then
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,n_s_id,len=linop%n_s))
else
   linop%n_s = 0
end if

! Get source/destination dimensions
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,trim(linop%prefix)//'_n_src',linop%n_src))
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,trim(linop%prefix)//'_n_dst',linop%n_dst))
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,trim(linop%prefix)//'_nvec',nvec))

! Allocation
call linop%alloc(nvec)

if (linop%n_s>0) then
   ! Get variables id
   call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(linop%prefix)//'_row',row_id))
   call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(linop%prefix)//'_col',col_id))
   if (linop%nvec>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(linop%prefix)//'_Svec',Svec_id))
   else
      call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(linop%prefix)//'_S',S_id))
   end if

   ! Get variables
   call mpl%ncerr(subr,nf90_get_var(ncid,row_id,linop%row))
   call mpl%ncerr(subr,nf90_get_var(ncid,col_id,linop%col))
   if (linop%nvec>0) then
      call mpl%ncerr(subr,nf90_get_var(ncid,Svec_id,linop%Svec))
   else
      call mpl%ncerr(subr,nf90_get_var(ncid,S_id,linop%S))
   end if
end if

end subroutine linop_read

!----------------------------------------------------------------------
! Subroutine: linop_write
! Purpose: write
!----------------------------------------------------------------------
subroutine linop_write(linop,mpl,ncid)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop ! Linear operator
type(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid            ! NetCDF file ID

! Local variables
integer :: n_s_id,nvec_id,row_id,col_id,S_id,Svec_id
character(len=1024),parameter :: subr = 'linop_write'

! Start definition mode
call mpl%ncerr(subr,nf90_redef(ncid))

! Write source/destination dimensions
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_src',linop%n_src))
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_dst',linop%n_dst))
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_nvec',linop%nvec))

if (linop%n_s>0) then
   ! Define dimensions
   call mpl%ncerr(subr,nf90_def_dim(ncid,trim(linop%prefix)//'_n_s',linop%n_s,n_s_id))
   if (linop%nvec>0) call mpl%ncerr(subr,nf90_def_dim(ncid,trim(linop%prefix)//'_nvec', &
 & linop%nvec,nvec_id))

   ! Define variables
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_row',nf90_int,(/n_s_id/),row_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_col',nf90_int,(/n_s_id/),col_id))
   if (linop%nvec>0) then
      call mpl%ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_Svec',nc_kind_real,(/n_s_id,nvec_id/),Svec_id))
   else
      call mpl%ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_S',nc_kind_real,(/n_s_id/),S_id))
   end if

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Put variables
   call mpl%ncerr(subr,nf90_put_var(ncid,row_id,linop%row))
   call mpl%ncerr(subr,nf90_put_var(ncid,col_id,linop%col))
   if (linop%nvec>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,Svec_id,linop%Svec))
   else
      call mpl%ncerr(subr,nf90_put_var(ncid,S_id,linop%S))
   end if
else
   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))
end if

end subroutine linop_write

!----------------------------------------------------------------------
! Subroutine: linop_reorder
! Purpose: reorder linear operator
!----------------------------------------------------------------------
subroutine linop_reorder(linop,mpl)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop ! Linear operator
type(mpl_type),intent(inout) :: mpl      ! MPI data

! Local variables
integer :: row,i_s_s,i_s_e,n_s,i_s
integer,allocatable :: order(:)
character(len=1024),parameter :: subr = 'linop_reorder'

if ((linop%n_s>0).and.(linop%n_s<reorder_max)) then
   ! Sort with respect to row
   allocate(order(linop%n_s))
   call qsort(linop%n_s,linop%row,order)

   ! Sort col and S
   linop%col = linop%col(order)
   if (linop%nvec>0) then
      linop%Svec = linop%Svec(order,:)
   else
      linop%S = linop%S(order)
   end if
   deallocate(order)

   ! Sort with respect to col for each row
   row = minval(linop%row)
   i_s_s = 1
   i_s_e = mpl%msv%vali
   do i_s=1,linop%n_s
      if (linop%row(i_s)==row) then
         i_s_e = i_s
      else
         n_s = i_s_e-i_s_s+1
         allocate(order(n_s))
         call qsort(n_s,linop%col(i_s_s:i_s_e),order)
         order = order+i_s_s-1
         if (linop%nvec>0) then
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
   call mpl%warning(subr,'linear operator is too big, impossible to reorder')
end if

end subroutine linop_reorder

!----------------------------------------------------------------------
! Subroutine: linop_apply
! Purpose: apply linear operator
!----------------------------------------------------------------------
subroutine linop_apply(linop,mpl,fld_src,fld_dst,ivec,mssrc,msdst)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop               ! Linear operator
type(mpl_type),intent(inout) :: mpl                 ! MPI data
real(kind_real),intent(in) :: fld_src(linop%n_src)  ! Source vector
real(kind_real),intent(out) :: fld_dst(linop%n_dst) ! Destination vector
integer,intent(in),optional :: ivec                 ! Index of the vector of linear operators with similar row and col
logical,intent(in),optional :: mssrc                ! Check for missing source
logical,intent(in),optional :: msdst                ! Check for missing destination

! Local variables
integer :: i_s,i_dst
logical :: lmssrc,lmsdst,valid
logical,allocatable :: missing(:)
character(len=1024),parameter :: subr = 'linop_apply'

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call mpl%abort(subr,'col<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call mpl%abort(subr,'col>n_src for linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call mpl%abort(subr,'row<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call mpl%abort(subr,'row>n_dst for linear operation '//trim(linop%prefix))
   if (present(ivec)) then
      if (any(isnan(linop%Svec))) call mpl%abort(subr,'NaN in Svec for linear operation '//trim(linop%prefix))
   else
      if (any(isnan(linop%S))) call mpl%abort(subr,'NaN in S for linear operation '//trim(linop%prefix))
   end if

   ! Check input
   if (any(fld_src>huge(1.0))) call mpl%abort(subr,'Overflowing number in fld_src for linear operation '//trim(linop%prefix))
   if (any(isnan(fld_src))) call mpl%abort(subr,'NaN in fld_src for linear operation '//trim(linop%prefix))
end if

! Initialization
fld_dst = 0.0
lmssrc = .false.
if (present(mssrc)) lmssrc = mssrc
lmsdst = .true.
if (present(msdst)) lmsdst = msdst
if (lmsdst) then
   allocate(missing(linop%n_dst))
   missing = .true.
end if

! Apply weights
do i_s=1,linop%n_s
   if (lmssrc) then
      ! Check for missing source (WARNING: source-dependent => no adjoint)
      valid = mpl%msv%isnotr(fld_src(linop%col(i_s)))
   else
      ! Source independent
      valid = .true.
   end if

   if (valid) then
      if (present(ivec)) then
         fld_dst(linop%row(i_s)) = fld_dst(linop%row(i_s))+linop%Svec(i_s,ivec)*fld_src(linop%col(i_s))
      else
         fld_dst(linop%row(i_s)) = fld_dst(linop%row(i_s))+linop%S(i_s)*fld_src(linop%col(i_s))
      end if

      ! Check for missing destination
      if (lmsdst) missing(linop%row(i_s)) = .false.
   end if
end do

if (lmsdst) then
   ! Missing destination values
   do i_dst=1,linop%n_dst
      if (missing(i_dst)) fld_dst(i_dst) = mpl%msv%valr
   end do

   ! Release memory
   deallocate(missing)
end if

if (check_data) then
   ! Check output
   if (any(isnan(fld_dst))) call mpl%abort(subr,'NaN in fld_dst for linear operation '//trim(linop%prefix))
end if

end subroutine linop_apply

!----------------------------------------------------------------------
! Subroutine: linop_apply_ad
! Purpose: apply linear operator, adjoint
!----------------------------------------------------------------------
subroutine linop_apply_ad(linop,mpl,fld_dst,fld_src,ivec)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop               ! Linear operator
type(mpl_type),intent(inout) :: mpl                 ! MPI data
real(kind_real),intent(in) :: fld_dst(linop%n_dst)  ! Destination vector
real(kind_real),intent(out) :: fld_src(linop%n_src) ! Source vector
integer,intent(in),optional :: ivec                 ! Index of the vector of linear operators with similar row and col

! Local variables
integer :: i_s
character(len=1024),parameter :: subr = 'linop_apply_ad'

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call mpl%abort(subr,'col<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call mpl%abort(subr,'col>n_src for adjoint linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call mpl%abort(subr,'row<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call mpl%abort(subr,'row>n_dst for adjoint linear operation '//trim(linop%prefix))
   if (present(ivec)) then
      if (any(isnan(linop%Svec))) call mpl%abort(subr,'NaN in Svec for adjoint linear operation '//trim(linop%prefix))
   else
      if (any(isnan(linop%S))) call mpl%abort(subr,'NaN in S for adjoint linear operation '//trim(linop%prefix))
   end if

   ! Check input
   if (any(fld_dst>huge(1.0))) &
 & call mpl%abort(subr,'Overflowing number in fld_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(isnan(fld_dst))) call mpl%abort(subr,'NaN in fld_dst for adjoint linear operation '//trim(linop%prefix))
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
   if (any(isnan(fld_src))) call mpl%abort(subr,'NaN in fld_src for adjoint linear operation '//trim(linop%prefix))
end if

end subroutine linop_apply_ad

!----------------------------------------------------------------------
! Subroutine: linop_apply_sym
! Purpose: apply linear operator, symmetric
!----------------------------------------------------------------------
subroutine linop_apply_sym(linop,mpl,fld,ivec)

implicit none

! Passed variables
class(linop_type),intent(in) :: linop             ! Linear operator
type(mpl_type),intent(inout) :: mpl               ! MPI data
real(kind_real),intent(inout) :: fld(linop%n_src) ! Source/destination vector
integer,intent(in),optional :: ivec               ! Index of the vector of linear operators with similar row and col

! Local variables
integer :: i_s,ithread
real(kind_real) :: fld_arr(linop%n_dst,mpl%nthread)
character(len=1024),parameter :: subr = 'linop_apply_sym'

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call mpl%abort(subr,'col<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call mpl%abort(subr,'col>n_src for symmetric linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call mpl%abort(subr,'row<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_src) call mpl%abort(subr,'row>n_dst for symmetric linear operation '//trim(linop%prefix))
   if (present(ivec)) then
      if (any(isnan(linop%Svec))) call mpl%abort(subr,'NaN in Svec for symmetric linear operation '//trim(linop%prefix))
   else
      if (any(isnan(linop%S))) call mpl%abort(subr,'NaN in S for symmetric linear operation '//trim(linop%prefix))
   end if

   ! Check input
   if (any(fld>huge(1.0))) call mpl%abort(subr,'Overflowing number in fld for symmetric linear operation '//trim(linop%prefix))
   if (any(isnan(fld))) call mpl%abort(subr,'NaN in fld for symmetric linear operation '//trim(linop%prefix))
end if

! Apply weights
fld_arr = 0.0
!$omp parallel do schedule(static) private(i_s,ithread)
do i_s=1,linop%n_s
   ithread = 1
!$ ithread = omp_get_thread_num()+1
   if (present(ivec)) then
      fld_arr(linop%row(i_s),ithread) = fld_arr(linop%row(i_s),ithread)+linop%Svec(i_s,ivec)*fld(linop%col(i_s))
      if (linop%col(i_s)/=linop%row(i_s)) fld_arr(linop%col(i_s),ithread) = fld_arr(linop%col(i_s),ithread) &
                                                                          & +linop%Svec(i_s,ivec)*fld(linop%row(i_s))
   else
      fld_arr(linop%row(i_s),ithread) = fld_arr(linop%row(i_s),ithread)+linop%S(i_s)*fld(linop%col(i_s))
      if (linop%col(i_s)/=linop%row(i_s)) fld_arr(linop%col(i_s),ithread) = fld_arr(linop%col(i_s),ithread) &
                                                                          & +linop%S(i_s)*fld(linop%row(i_s))
   end if
end do
!$omp end parallel do

! Sum over threads
fld = 0.0
do ithread=1,mpl%nthread
   fld = fld+fld_arr(:,ithread)
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld))) call mpl%abort(subr,'NaN in fld for symmetric linear operation '//trim(linop%prefix))
end if

end subroutine linop_apply_sym

!----------------------------------------------------------------------
! Subroutine: linop_add_op
! Purpose: add operation
!----------------------------------------------------------------------
subroutine linop_add_op(linop,n_s,row,col,S)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop ! Linear operators
integer,intent(inout) :: n_s             ! Number of operations
integer,intent(in) :: row                ! Row index
integer,intent(in) :: col                ! Column index
real(kind_real),intent(in) :: S          ! Value

! Local variables
type(linop_type) :: linop_tmp

! Update
n_s = n_s+1
if (n_s>linop%n_s) then
   ! Copy
   call linop_tmp%copy(linop)

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
! Subroutine: linop_gather
! Purpose: gather data from OpenMP threads
!----------------------------------------------------------------------
subroutine linop_gather(linop,mpl,n_s_arr,linop_arr)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop              ! Linear operator
type(mpl_type),intent(inout) :: mpl                   ! MPI data
integer,intent(in) :: n_s_arr(mpl%nthread)            ! Number of operations
type(linop_type),intent(in) :: linop_arr(mpl%nthread) ! Linear operator array

! Local variables
integer :: ithread,offset

! Total number of operations
linop%n_s = sum(n_s_arr)

! Allocation
call linop%alloc

! Gather data
offset = 0
do ithread=1,mpl%nthread
   linop%row(offset+1:offset+n_s_arr(ithread)) = linop_arr(ithread)%row(1:n_s_arr(ithread))
   linop%col(offset+1:offset+n_s_arr(ithread)) = linop_arr(ithread)%col(1:n_s_arr(ithread))
   linop%S(offset+1:offset+n_s_arr(ithread)) = linop_arr(ithread)%S(1:n_s_arr(ithread))
   offset = offset+n_s_arr(ithread)
end do

end subroutine linop_gather

!----------------------------------------------------------------------
! Subroutine: linop_interp_from_lat_lon
! Purpose: compute horizontal interpolation from source latitude/longitude
!----------------------------------------------------------------------
subroutine linop_interp_from_lat_lon(linop,mpl,rng,n_src,lon_src,lat_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop     ! Linear operator
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(rng_type),intent(inout) :: rng          ! Random number generator
integer,intent(in) :: n_src                  ! Source size
real(kind_real),intent(in) :: lon_src(n_src) ! Source longitudes
real(kind_real),intent(in) :: lat_src(n_src) ! Source latitudes
logical,intent(in) :: mask_src(n_src)        ! Source mask
integer,intent(in) :: n_dst                  ! Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) ! Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst) ! Destination latitudes
logical,intent(in) :: mask_dst(n_dst)        ! Destination mask
character(len=*),intent(in) :: interp_type   ! Interpolation type

! Local variables
integer :: n_src_eff,i_src,i_src_eff
integer,allocatable :: src_eff_to_src(:)
real(kind_real),allocatable :: lon_src_eff(:),lat_src_eff(:)
logical,allocatable :: mask_src_eff(:)
type(tree_type) :: tree
type(mesh_type) :: mesh

! Count non-missing source points
n_src_eff = count(mask_src)

! Allocation
allocate(src_eff_to_src(n_src_eff))
allocate(lon_src_eff(n_src_eff))
allocate(lat_src_eff(n_src_eff))
allocate(mask_src_eff(n_src_eff))

! Conversion
i_src_eff = 0
do i_src=1,n_src
   if (mask_src(i_src)) then
      i_src_eff = i_src_eff+1
      src_eff_to_src(i_src_eff) = i_src
   end if
end do
lon_src_eff = lon_src(src_eff_to_src)
lat_src_eff = lat_src(src_eff_to_src)
mask_src_eff = .true.

! Allocation
call mesh%alloc(n_src_eff)
call tree%alloc(mpl,n_src_eff)

! Initialization
call mesh%init(mpl,rng,lon_src_eff,lat_src_eff)
call tree%init(lon_src_eff,lat_src_eff)

! Compute interpolation
call linop%interp(mpl,mesh,tree,n_src_eff,mask_src_eff,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

! Effective points conversion
linop%n_src = n_src
linop%col = src_eff_to_src(linop%col)

! Release memory
deallocate(src_eff_to_src)
deallocate(lon_src_eff)
deallocate(lat_src_eff)
deallocate(mask_src_eff)
call mesh%dealloc
call tree%dealloc

end subroutine linop_interp_from_lat_lon

!----------------------------------------------------------------------
! Subroutine: linop_interp_from_mesh_tree
! Purpose: compute horizontal interpolation from source mesh and tree
!----------------------------------------------------------------------
subroutine linop_interp_from_mesh_tree(linop,mpl,mesh,tree,n_src,mask_src,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop     ! Linear operator
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(mesh_type),intent(in) :: mesh           ! Mesh
type(tree_type),intent(in) :: tree       ! KD-tree
integer,intent(in) :: n_src                  ! Source size
logical,intent(in) :: mask_src(n_src)        ! Source mask
integer,intent(in) :: n_dst                  ! Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) ! Destination longitudes
real(kind_real),intent(in) :: lat_dst(n_dst) ! Destination latitudes
logical,intent(in) :: mask_dst(n_dst)        ! Destination mask
character(len=*),intent(in) :: interp_type   ! Interpolation type

! Local variables
integer :: i,i_src,i_dst,nn_index(1),n_s,ib(3),nnat,inat,np,iproc,offset,i_s
integer :: n_dst_loc(0:mpl%nproc),i_dst_loc,proc_to_n_s(mpl%nproc),n_s_loc(0:mpl%nproc)
integer,allocatable :: natis(:),row(:),col(:)
real(kind_real) :: nn_dist(1),b(3)
real(kind_real),allocatable :: area_polygon(:),area_polygon_new(:),natwgt(:),S(:)
logical :: loop
logical,allocatable :: missing(:)
character(len=1024),parameter :: subr = 'linop_interp_from_mesh_tree'
type(mesh_type) :: meshnew

! MPI splitting
call mpl%split(n_dst,n_dst_loc)

! Allocation
np = mpl%msv%vali
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
   call mpl%abort(subr,'wrong interpolation type')
end if
allocate(row(np*n_dst_loc(mpl%myproc)))
allocate(col(np*n_dst_loc(mpl%myproc)))
allocate(S(np*n_dst_loc(mpl%myproc)))

if (trim(interp_type)=='natural') then
   ! Compute polygons areas
   do i_src=1,mesh%n
      natis(i_src) = i_src
   end do
   call mesh%polygon(mesh%n,natis,area_polygon)
end if

! Compute interpolation
write(mpl%info,'(a10,a)') '','Compute interpolation: '
call mpl%flush(.false.)
call mpl%prog_init(n_dst_loc(mpl%myproc))
n_s = 0
do i_dst_loc=1,n_dst_loc(mpl%myproc)
   ! Indices
   i_dst = sum(n_dst_loc(0:mpl%myproc-1))+i_dst_loc

   if (mask_dst(i_dst)) then
      ! Find nearest neighbor
      call tree%find_nearest_neighbors(lon_dst(i_dst),lat_dst(i_dst),1,nn_index,nn_dist)

      if (abs(nn_dist(1))>0.0) then
         ! Compute barycentric coordinates
         call mesh%barycentric(mpl,lon_dst(i_dst),lat_dst(i_dst),nn_index(1),b,ib)
         if (sum(b)>0.0) b = b/sum(b)
         do i=1,3
            if (inf(b(i),S_inf)) b(i) = 0.0
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
                  call meshnew%copy(mesh)

                  ! Add a node
                  call meshnew%addnode(mpl,lon_dst(i_dst),lat_dst(i_dst))

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
                        if (inf(natwgt(inat),S_inf)) natwgt(inat) = 0.0
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

   ! Update
   call mpl%prog_print(i_dst_loc)
end do
call mpl%prog_final

! Communication
call mpl%f_comm%allgather(n_s,proc_to_n_s)

! Allocation
linop%n_s = sum(proc_to_n_s)
linop%n_src = n_src
linop%n_dst = n_dst
call linop%alloc

! Copy data
offset = 0
do iproc=1,mpl%nproc
   if (proc_to_n_s(iproc)>0) then
      if (iproc==mpl%myproc) then
         linop%row(offset+1:offset+proc_to_n_s(iproc)) = row(1:proc_to_n_s(iproc))
         linop%col(offset+1:offset+proc_to_n_s(iproc)) = col(1:proc_to_n_s(iproc))
         linop%S(offset+1:offset+proc_to_n_s(iproc)) = S(1:proc_to_n_s(iproc))
      end if
   end if
   offset = offset+proc_to_n_s(iproc)
end do

! MPI sharing
n_s_loc(0) = 0
n_s_loc(1:mpl%nproc) = proc_to_n_s
call mpl%share(linop%n_s,n_s_loc,linop%row)
call mpl%share(linop%n_s,n_s_loc,linop%col)
call mpl%share(linop%n_s,n_s_loc,linop%S)

! Deal with missing points
call linop%interp_missing(mpl,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

! Check interpolation
allocate(missing(n_dst))
missing = .false.
do i_dst=1,n_dst
   if (mask_dst(i_dst)) missing(i_dst) = .true.
end do
do i_s=1,linop%n_s
   missing(linop%row(i_s)) = .false.
end do
if (any(missing)) call mpl%abort(subr,'missing destination points in interp_from_mesh_tree')

! Release memory
if (trim(interp_type)=='natural') then
   deallocate(area_polygon)
   deallocate(area_polygon_new)
   deallocate(natis)
   deallocate(natwgt)
end if
deallocate(row)
deallocate(col)
deallocate(S)
deallocate(missing)

end subroutine linop_interp_from_mesh_tree

!----------------------------------------------------------------------
! Subroutine: linop_interp_grid
! Purpose: compute horizontal grid interpolation
!----------------------------------------------------------------------
subroutine linop_interp_grid(linop,mpl,rng,geom,il0i,nc1,c1_to_c0,mask_check,vbot,vtop,interp_type,interp_base)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop      ! Linear operator
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(rng_type),intent(inout) :: rng           ! Random number generator
type(geom_type),intent(in) :: geom            ! Geometry
integer,intent(in) :: il0i                    ! Level
integer,intent(in) :: nc1                     ! Subset Sc1 size
integer,intent(in) :: c1_to_c0(nc1)           ! Subset Sc1 to subset Sc0
logical,intent(in) :: mask_check              ! Mask check key
integer,intent(in) :: vbot(nc1)               ! Bottom level
integer,intent(in) :: vtop(nc1)               ! Top level
character(len=*),intent(in) :: interp_type    ! Interpolation type
type(linop_type),intent(inout) :: interp_base ! Linear operator (base interpolation)

! Local variables
integer :: ic0,ic1,jc0,jc1,i_s
real(kind_real) :: renorm(geom%nc0)
real(kind_real),allocatable :: lon_c1(:),lat_c1(:),lon_col(:),lat_col(:)
logical :: test_c0(geom%nc0),mask_extra(nc1)
logical,allocatable :: mask_c1(:),valid(:)
character(len=1024),parameter :: subr = 'linop_interp_grid'

if (.not.allocated(interp_base%row)) then
   ! Allocation
   allocate(lon_c1(nc1))
   allocate(lat_c1(nc1))
   allocate(mask_c1(nc1))

   ! Initialization
   lon_c1 = geom%lon(c1_to_c0)
   lat_c1 = geom%lat(c1_to_c0)
   mask_c1 = geom%mask_hor_c0(c1_to_c0)

   ! Compute base interpolation
   call interp_base%interp(mpl,rng,nc1,lon_c1,lat_c1,mask_c1,geom%nc0,geom%lon,geom%lat,geom%mask_hor_c0,interp_type)

   ! Release memory
   deallocate(lon_c1)
   deallocate(lat_c1)
   deallocate(mask_c1)
end if

! Allocation
allocate(valid(interp_base%n_s))

! Check mask
do i_s=1,interp_base%n_s
   ic0 = interp_base%row(i_s)
   jc1 = interp_base%col(i_s)
   jc0 = c1_to_c0(jc1)
   valid(i_s) = geom%mask_c0(ic0,il0i).and.geom%mask_c0(jc0,il0i)
end do

if (mask_check) then
   write(mpl%info,'(a10,a,i3,a)') '','Sublevel ',il0i,': '
   call mpl%flush(.false.)

   ! Allocation
   allocate(lon_col(nc1))
   allocate(lat_col(nc1))

   ! Initialization
   lon_col = geom%lon(c1_to_c0)
   lat_col = geom%lat(c1_to_c0)

   ! Check mask boundaries
   call interp_base%check_mask(mpl,geom,valid,il0i,lon_col=lon_col,lat_col=lat_col)

   ! Release memory
   deallocate(lon_col)
   deallocate(lat_col)
end if

if (geom%nl0i>1) then
   ! Check extrapolated points because of masks vertical variations (unsampled levels)
   mask_extra = .false.
   do ic1=1,nc1
      ic0 = c1_to_c0(ic1)
      if (geom%mask_c0(ic0,il0i).and.((il0i<vbot(ic1)).or.(il0i>vtop(ic1)))) mask_extra(ic1) = .true.
   end do

   ! Remove operations for extrapolated points
   do i_s=1,interp_base%n_s
      if (valid(i_s)) then
         ic1 = interp_base%col(i_s)
         if (mask_extra(ic1)) valid(i_s) = .false.
      end if
   end do
   if (count(mask_extra)>0) then
      write(mpl%info,'(a10,a,i5)') '','Extrapolated points: ',count(mask_extra)
      call mpl%flush
   end if
else
   mask_extra = .false.
end if

! Renormalization
renorm = 0.0
do i_s=1,interp_base%n_s
   if (valid(i_s)) renorm(interp_base%row(i_s)) = renorm(interp_base%row(i_s))+interp_base%S(i_s)
end do

! Allocation
linop%n_src = nc1
linop%n_dst = geom%nc0
linop%n_s = count(valid)
call linop%alloc

! Initialization
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
call linop%interp_missing(mpl,geom%nc0,geom%lon,geom%lat,geom%mask_c0(:,il0i),interp_type)

! Check interpolation
test_c0 = geom%mask_c0(:,il0i)
do i_s=1,linop%n_s
   test_c0(linop%row(i_s)) = .false.
end do
if (any(test_c0)) call mpl%abort(subr,'error with the grid interpolation row')

! Release memory
deallocate(valid)

end subroutine linop_interp_grid

!----------------------------------------------------------------------
! Subroutine: linop_check_mask
! Purpose: check mask boundaries for linear operators
!----------------------------------------------------------------------
subroutine linop_check_mask(linop,mpl,geom,valid,il0,lon_row,lat_row,lon_col,lat_col)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop                    ! Linear operator
type(mpl_type),intent(inout) :: mpl                         ! MPI data
type(geom_type),intent(in) :: geom                          ! Geometry
logical,intent(inout) :: valid(linop%n_s)                   ! Valid points
integer,intent(in) :: il0                                   ! Level
real(kind_real),intent(in),optional :: lon_row(linop%n_dst) ! Row longitude (Sc0 longitude if missing)
real(kind_real),intent(in),optional :: lat_row(linop%n_dst) ! Row latitude (Sc0 latitude if missing)
real(kind_real),intent(in),optional :: lon_col(linop%n_src) ! Column longitude (Sc0 longitude if missing)
real(kind_real),intent(in),optional :: lat_col(linop%n_src) ! Column latitude (Sc0 latitude if missing)

! Local variables
integer :: i_s
integer :: n_s_loc(0:mpl%nproc),i_s_loc
real(kind_real) :: llon_row,llat_row,llon_col,llat_col

! MPI splitting
call mpl%split(linop%n_s,n_s_loc)

! Check that interpolations are not crossing mask boundaries
call mpl%prog_init(n_s_loc(mpl%myproc))
!$omp parallel do schedule(static) private(i_s_loc,i_s,llon_row,llat_row,llon_col,llat_col)
do i_s_loc=1,n_s_loc(mpl%myproc)
   ! Indices
   i_s = sum(n_s_loc(0:mpl%myproc-1))+i_s_loc

   if (valid(i_s)) then
      ! Row lon/lat
      if (present(lon_row).and.present(lat_row)) then
         llon_row = lon_row(linop%row(i_s))
         llat_row = lat_row(linop%row(i_s))
      else
         llon_row = geom%lon(linop%row(i_s))
         llat_row = geom%lat(linop%row(i_s))
      end if

      ! Column lon/lat
      if (present(lon_col).and.present(lat_col)) then
         llon_col = lon_col(linop%col(i_s))
         llat_col = lat_col(linop%col(i_s))
      else
         llon_col = geom%lon(linop%col(i_s))
         llat_col = geom%lat(linop%col(i_s))
      end if

      ! Check if arc is crossing boundary arcs
      call geom%check_arc(mpl,il0,llon_row,llat_row,llon_col,llat_col,valid(i_s))
   end if

   ! Update
   call mpl%prog_print(i_s_loc)
end do
!$omp end parallel do
call mpl%prog_final

! MPI sharing
call mpl%share(linop%n_s,n_s_loc,valid)

end subroutine linop_check_mask

!----------------------------------------------------------------------
! Subroutine: linop_interp_missing
! Purpose: deal with missing interpolation points
!----------------------------------------------------------------------
subroutine linop_interp_missing(linop,mpl,n_dst,lon_dst,lat_dst,mask_dst,interp_type)

implicit none

! Passed variables
class(linop_type),intent(inout) :: linop     ! Linear operator (interpolation)
type(mpl_type),intent(inout) :: mpl          ! MPI data
integer,intent(in) :: n_dst                  ! Destination size
real(kind_real),intent(in) :: lon_dst(n_dst) ! Destination longitude
real(kind_real),intent(in) :: lat_dst(n_dst) ! Destination latitude
logical,intent(in) :: mask_dst(n_dst)        ! Destination mask
character(len=*),intent(in) :: interp_type   ! Interpolation type

! Local variables
integer :: i_dst,i_s
integer :: nn(1)
logical :: missing(n_dst),lmask(n_dst),found
character(len=1024),parameter :: subr = 'linop_interp_missing'
type(linop_type) :: interp_tmp
type(tree_type) :: tree

! Find missing points
missing = .false.
do i_dst=1,n_dst
   if (mask_dst(i_dst)) missing(i_dst) = .true.
end do
do i_s=1,linop%n_s
   missing(linop%row(i_s)) = .false.
end do

if (count(missing)>0) then
   write(mpl%info,'(a10,a,i6,a)') '','Deal with ',count(missing),' missing interpolation points'
   call mpl%flush

   ! Allocate temporary interpolation
   if (trim(interp_type)=='bilin') then
      interp_tmp%n_s = linop%n_s+3*count(missing)
   elseif (trim(interp_type)=='natural') then
      interp_tmp%n_s = linop%n_s+40*count(missing)
   else
      call mpl%abort(subr,'wrong interpolation')
   end if
   call interp_tmp%alloc

   ! Fill arrays
   interp_tmp%row(1:linop%n_s) = linop%row
   interp_tmp%col(1:linop%n_s) = linop%col
   interp_tmp%S(1:linop%n_s) = linop%S

   ! Reset size
   interp_tmp%n_s = linop%n_s

   ! Allocation
   lmask = mask_dst.and.(.not.missing)
   call tree%alloc(mpl,n_dst,mask=lmask)

   ! Initialization
   call tree%init(lon_dst,lat_dst)

   do i_dst=1,n_dst
      if (missing(i_dst)) then
         ! Compute nearest neighbor
         call tree%find_nearest_neighbors(lon_dst(i_dst),lat_dst(i_dst),1,nn)

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
         if (.not.found) call mpl%abort(subr,'missing point not found')
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
   call tree%dealloc
   call interp_tmp%dealloc
end if

end subroutine linop_interp_missing

end module type_linop
