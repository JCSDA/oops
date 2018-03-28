!----------------------------------------------------------------------
! Module: type_linop
!> Purpose: linear operator derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_linop

use netcdf
use omp_lib
use tools_display, only: msgerror,msgwarning
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr
use tools_nc, only: ncfloat,ncerr
use tools_qsort, only: qsort
use type_mpl, only: mpl

implicit none

! Linear operator derived type
type linoptype
   character(len=1024) :: prefix       !< Operator prefix (for I/O)
   integer :: n_src                    !< Source vector size
   integer :: n_dst                    !< Destination vector size
   integer :: n_s                      !< Operator size
   integer,allocatable :: row(:)       !< Output indices
   integer,allocatable :: col(:)       !< Input indices
   real(kind_real),allocatable :: S(:) !< Coefficients
end type linoptype

interface linop_read
  module procedure linop_read_0d
  module procedure linop_read_1d
  module procedure linop_read_2d
end interface
interface linop_write
  module procedure linop_write_0d
  module procedure linop_write_1d
  module procedure linop_write_2d
end interface

logical :: check_data = .false. !< Activate data check for all linear operations

private
public :: linoptype
public :: linop_alloc,linop_dealloc,linop_copy,linop_reorder, &
 & apply_linop,apply_linop_ad,apply_linop_sym,linop_read,linop_write

contains

!----------------------------------------------------------------------
! Subroutine: linop_alloc
!> Purpose: linear operator object allocation
!----------------------------------------------------------------------
subroutine linop_alloc(linop)

implicit none

! Passed variables
type(linoptype),intent(inout) :: linop !< Linear operator

! Allocation
allocate(linop%row(linop%n_s))
allocate(linop%col(linop%n_s))
allocate(linop%S(linop%n_s))

! Initialization
call msi(linop%row)
call msi(linop%col)
call msr(linop%S)

end subroutine linop_alloc

!----------------------------------------------------------------------
! Subroutine: linop_dealloc
!> Purpose: linear operator object deallocation
!----------------------------------------------------------------------
subroutine linop_dealloc(linop)

implicit none

! Passed variables
type(linoptype),intent(inout) :: linop !< Linear operator

! Release memory
if (allocated(linop%row)) deallocate(linop%row)
if (allocated(linop%col)) deallocate(linop%col)
if (allocated(linop%S)) deallocate(linop%S)

end subroutine linop_dealloc

!----------------------------------------------------------------------
! Subroutine: linop_copy
!> Purpose: linear operator object copy
!----------------------------------------------------------------------
subroutine linop_copy(linop_in,linop_out)

implicit none

! Passed variables
type(linoptype),intent(in) :: linop_in     !< Input linear operator
type(linoptype),intent(inout) :: linop_out !< Output linear operator

! Copy attributes
linop_out%prefix = trim(linop_in%prefix)
linop_out%n_src = linop_in%n_src
linop_out%n_dst = linop_in%n_dst
linop_out%n_s = linop_in%n_s

! Deallocation
call linop_dealloc(linop_out)

! Allocation
call linop_alloc(linop_out)

! Copy data
linop_out%row = linop_in%row
linop_out%col = linop_in%col
linop_out%S = linop_in%S

end subroutine linop_copy

!----------------------------------------------------------------------
! Subroutine: linop_reorder
!> Purpose: reorder linear operator
!----------------------------------------------------------------------
subroutine linop_reorder(linop)

implicit none

! Passed variables
type(linoptype),intent(inout) :: linop !< Linear operator

! Local variables
integer :: row,i_s_s,i_s_e,n_s,i_s
integer,allocatable :: order(:)

if (linop%n_s<1000000) then
   ! Sort with respect to row
   allocate(order(linop%n_s))
   call qsort(linop%n_s,linop%row,order)

   ! Sort col and S
   linop%col = linop%col(order)
   linop%S = linop%S(order)
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
         linop%S(i_s_s:i_s_e) = linop%S(order)
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
! Subroutine: apply_linop
!> Purpose: apply linear operator
!----------------------------------------------------------------------
subroutine apply_linop(linop,fld_src,fld_dst)

implicit none

! Passed variables
type(linoptype),intent(in) :: linop                 !< Linear operator
real(kind_real),intent(in) :: fld_src(linop%n_src)  !< Source vector
real(kind_real),intent(out) :: fld_dst(linop%n_dst) !< Destination vector

! Local variables
integer :: i_s

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call msgerror('row>n_dst for linear operation '//trim(linop%prefix))
   if (any(isnan(linop%S))) call msgerror('NaN in S for linear operation '//trim(linop%prefix))

   ! Check input
   if (any(fld_src>huge(1.0))) call msgerror('Overflowing number in fld_src for linear operation '//trim(linop%prefix))
   if (any(isnan(fld_src))) call msgerror('NaN in fld_src for linear operation '//trim(linop%prefix))
   if (any(.not.isnotmsr(fld_src))) call msgerror('Missing value in fld_src for linear operation '//trim(linop%prefix))
end if

! Initialization
fld_dst = 0.0

! Apply weights
do i_s=1,linop%n_s
   fld_dst(linop%row(i_s)) = fld_dst(linop%row(i_s))+linop%S(i_s)*fld_src(linop%col(i_s))
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld_dst))) call msgerror('NaN in fld_dst for linear operation '//trim(linop%prefix))
end if

end subroutine apply_linop

!----------------------------------------------------------------------
! Subroutine: apply_linop_ad
!> Purpose: apply linear operator, adjoint
!----------------------------------------------------------------------
subroutine apply_linop_ad(linop,fld_dst,fld_src)

implicit none

! Passed variables
type(linoptype),intent(in) :: linop                 !< Linear operator
real(kind_real),intent(in) :: fld_dst(linop%n_dst)  !< Destination vector
real(kind_real),intent(out) :: fld_src(linop%n_src) !< Source vector

! Local variables
integer :: i_s

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for adjoint linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call msgerror('row>n_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(isnan(linop%S))) call msgerror('NaN in S for adjoint linear operation '//trim(linop%prefix))

   ! Check input
   if (any(fld_dst>huge(1.0))) call msgerror('Overflowing number in fld_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(isnan(fld_dst))) call msgerror('NaN in fld_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(.not.isnotmsr(fld_dst))) call msgerror('Missing value in fld_dst for adjoint linear operation '//trim(linop%prefix))
end if

! Initialization
fld_src = 0.0

! Apply weights
do i_s=1,linop%n_s
   fld_src(linop%col(i_s)) = fld_src(linop%col(i_s))+linop%S(i_s)*fld_dst(linop%row(i_s))
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld_src))) call msgerror('NaN in fld_src for adjoint linear operation '//trim(linop%prefix))
end if

end subroutine apply_linop_ad

!----------------------------------------------------------------------
! Subroutine: apply_linop_sym
!> Purpose: apply linear operator, symmetric
!----------------------------------------------------------------------
subroutine apply_linop_sym(linop,fld)

implicit none

! Passed variables
type(linoptype),intent(in) :: linop               !< Linear operator
real(kind_real),intent(inout) :: fld(linop%n_src) !< Source/destination vector

! Local variables
integer :: i_s,ithread
real(kind_real) :: fld_arr(linop%n_dst,mpl%nthread)

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for symmetric linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_src) call msgerror('row>n_dst for symmetric linear operation '//trim(linop%prefix))
   if (any(isnan(linop%S))) call msgerror('NaN in S for symmetric linear operation '//trim(linop%prefix))

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
   fld_arr(linop%row(i_s),ithread) = fld_arr(linop%row(i_s),ithread)+linop%S(i_s)*fld(linop%col(i_s))
   fld_arr(linop%col(i_s),ithread) = fld_arr(linop%col(i_s),ithread)+linop%S(i_s)*fld(linop%row(i_s))
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

end subroutine apply_linop_sym

!----------------------------------------------------------------------
! Subroutine: linop_read_0d
!> Purpose: read single linear operator from a NetCDF file
!----------------------------------------------------------------------
subroutine linop_read_0d(ncid,prefix,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: prefix  !< Linear operator prefix
type(linoptype),intent(inout) :: linop !< Linear operator

! Local variables
integer :: info
integer :: n_s_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_read_0d'

! Copy prefix
linop%prefix = trim(prefix)

! Get operator size
info = nf90_inq_dimid(ncid,trim(prefix)//'_n_s',n_s_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n_s_id,len=linop%n_s))
else
   linop%n_s = 0
end if

if (linop%n_s>0) then
   ! Allocation
   call linop_alloc(linop)

   ! Get variables id
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_row',row_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_col',col_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_S',S_id))

   ! Get source/destination dimensions
   call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(prefix)//'_n_src',linop%n_src))
   call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(prefix)//'_n_dst',linop%n_dst))

   ! Get variables
   call ncerr(subr,nf90_get_var(ncid,row_id,linop%row))
   call ncerr(subr,nf90_get_var(ncid,col_id,linop%col))
   call ncerr(subr,nf90_get_var(ncid,S_id,linop%S))
end if

end subroutine linop_read_0d

!----------------------------------------------------------------------
! Subroutine: linop_read_1d
!> Purpose: read array of linear operators from a NetCDF file, 1D
!----------------------------------------------------------------------
subroutine linop_read_1d(ncid,prefix,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid                            !< NetCDF file id
character(len=*),intent(in) :: prefix                 !< Linear operator prefix
type(linoptype),allocatable,intent(inout) :: linop(:) !< Linear operators

! Local variables
integer :: info,n1,n_s_max,i1
integer :: n_s_max_id,n1_id,n_s_id,n_src_id,n_dst_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_read_1d'

! Get maximum operator size
info = nf90_inq_dimid(ncid,trim(prefix)//'_n_s_max',n_s_max_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n_s_max_id,len=n_s_max))
else
   n_s_max = 0
end if

! Get array size
info = nf90_inq_dimid(ncid,trim(prefix)//'_n1',n1_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n1_id,len=n1))
else
   n1 = 0
end if

if ((n1>0).and.(n_s_max>0)) then
   ! Allocation
   allocate(linop(n1))

   ! Get variables id
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_n_s',n_s_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_n_src',n_src_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_n_dst',n_dst_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_row',row_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_col',col_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_S',S_id))

   do i1=1,n1
      ! Copy prefix
      linop(i1)%prefix = trim(prefix)

      ! Get operator size
      call ncerr(subr,nf90_get_var(ncid,n_s_id,linop(i1)%n_s,(/i1/)))
      call ncerr(subr,nf90_get_var(ncid,n_src_id,linop(i1)%n_src,(/i1/)))
      call ncerr(subr,nf90_get_var(ncid,n_dst_id,linop(i1)%n_dst,(/i1/)))

      if (linop(i1)%n_s>0) then
         ! Allocation
         call linop_alloc(linop(i1))

         ! Get variables
         call ncerr(subr,nf90_get_var(ncid,row_id,linop(i1)%row,(/1,i1/),(/linop(i1)%n_s,1/)))
         call ncerr(subr,nf90_get_var(ncid,col_id,linop(i1)%col,(/1,i1/),(/linop(i1)%n_s,1/)))
         call ncerr(subr,nf90_get_var(ncid,S_id,linop(i1)%S,(/1,i1/),(/linop(i1)%n_s,1/)))
      end if
   end do
end if

end subroutine linop_read_1d

!----------------------------------------------------------------------
! Subroutine: linop_read_2d
!> Purpose: read array of linear operators from a NetCDF file, 2D
!----------------------------------------------------------------------
subroutine linop_read_2d(ncid,prefix,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid                              !< NetCDF file id
character(len=*),intent(in) :: prefix                   !< Linear operator prefix
type(linoptype),allocatable,intent(inout) :: linop(:,:) !< Linear operators

! Local variables
integer :: info,n1,n2,n_s_max,i1,i2
integer :: n_s_max_id,n1_id,n2_id,n_s_id,n_src_id,n_dst_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_read_2d'

! Get maximum operator size
info = nf90_inq_dimid(ncid,trim(prefix)//'_n_s_max',n_s_max_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n_s_max_id,len=n_s_max))
else
   n_s_max = 0
end if

! Get array size
info = nf90_inq_dimid(ncid,trim(prefix)//'_n1',n1_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n1_id,len=n1))
else
   n1 = 0
end if
info = nf90_inq_dimid(ncid,trim(prefix)//'_n2',n2_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n2_id,len=n2))
else
   n2 = 0
end if

if ((n1>0).and.(n2>0).and.(n_s_max>0)) then
   ! Allocation
   allocate(linop(n1,n2))

   ! Get variables id
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_n_s',n_s_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_n_src',n_src_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_n_dst',n_dst_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_row',row_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_col',col_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_S',S_id))

   do i2=1,n2
      do i1=1,n1
         ! Copy prefix
         linop(i1,i2)%prefix = trim(prefix)

         ! Get operator size
         call ncerr(subr,nf90_get_var(ncid,n_s_id,linop(i1,i2)%n_s,(/i1,i2/)))
         call ncerr(subr,nf90_get_var(ncid,n_src_id,linop(i1,i2)%n_src,(/i1,i2/)))
         call ncerr(subr,nf90_get_var(ncid,n_dst_id,linop(i1,i2)%n_dst,(/i1,i2/)))

         if (linop(i1,i2)%n_s>0) then
            ! Allocation
            call linop_alloc(linop(i1,i2))

            ! Get variables
            call ncerr(subr,nf90_get_var(ncid,row_id,linop(i1,i2)%row,(/1,i1,i2/),(/linop(i1,i2)%n_s,1,1/)))
            call ncerr(subr,nf90_get_var(ncid,col_id,linop(i1,i2)%col,(/1,i1,i2/),(/linop(i1,i2)%n_s,1,1/)))
            call ncerr(subr,nf90_get_var(ncid,S_id,linop(i1,i2)%S,(/1,i1,i2/),(/linop(i1,i2)%n_s,1,1/)))
         end if
      end do
   end do
end if

end subroutine linop_read_2d

!----------------------------------------------------------------------
! Subroutine: linop_write_0d
!> Purpose: write single linear operator to a NetCDF file
!----------------------------------------------------------------------
subroutine linop_write_0d(ncid,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid          !< NetCDF file id
type(linoptype),intent(in) :: linop !< Linear operator

! Local variables
integer :: n_s_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_write_0d'

if (linop%n_s>0) then
   ! Start definition mode
   call ncerr(subr,nf90_redef(ncid))

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,trim(linop%prefix)//'_n_s',linop%n_s,n_s_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_row',nf90_int,(/n_s_id/),row_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_col',nf90_int,(/n_s_id/),col_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_S',ncfloat,(/n_s_id/),S_id))

   ! Write source/destination dimensions
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_src',linop%n_src))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_dst',linop%n_dst))

   ! End definition mode
   call ncerr(subr,nf90_enddef(ncid))

   ! Put variables
   call ncerr(subr,nf90_put_var(ncid,row_id,linop%row))
   call ncerr(subr,nf90_put_var(ncid,col_id,linop%col))
   call ncerr(subr,nf90_put_var(ncid,S_id,linop%S))
end if

end subroutine linop_write_0d

!----------------------------------------------------------------------
! Subroutine: linop_write_1d
!> Purpose: write array of linear operators to a NetCDF file, 1D
!----------------------------------------------------------------------
subroutine linop_write_1d(ncid,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
type(linoptype),intent(in) :: linop(:) !< Linear operator

! Local variables
integer :: n1,i1,n_s_max
integer :: n_s_max_id,n1_id,n_s_id,n_src_id,n_dst_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_write_1d'

! Array size
n1 = size(linop)

! Maximum operator size
n_s_max = 0
do i1=1,n1
   n_s_max = max(n_s_max,linop(i1)%n_s)
end do

if ((n1>0).and.(n_s_max>0)) then
   ! Start definition mode
   call ncerr(subr,nf90_redef(ncid))

   ! Define dimension
   call ncerr(subr,nf90_def_dim(ncid,trim(linop(1)%prefix)//'_n_s_max',n_s_max,n_s_max_id))
   call ncerr(subr,nf90_def_dim(ncid,trim(linop(1)%prefix)//'_n1',n1,n1_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_n_s',nf90_int,(/n1_id/),n_s_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_n_src',nf90_int,(/n1_id/),n_src_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_n_dst',nf90_int,(/n1_id/),n_dst_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_row',nf90_int,(/n_s_max_id,n1_id/),row_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_col',nf90_int,(/n_s_max_id,n1_id/),col_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_S',ncfloat,(/n_s_max_id,n1_id/),S_id))

   ! End definition mode
   call ncerr(subr,nf90_enddef(ncid))

   do i1=1,n1
      ! Put variables
      call ncerr(subr,nf90_put_var(ncid,n_s_id,linop(i1)%n_s,(/i1/)))
      call ncerr(subr,nf90_put_var(ncid,n_src_id,linop(i1)%n_src,(/i1/)))
      call ncerr(subr,nf90_put_var(ncid,n_dst_id,linop(i1)%n_dst,(/i1/)))
      call ncerr(subr,nf90_put_var(ncid,row_id,linop(i1)%row,(/1,i1/),(/linop(i1)%n_s,1/)))
      call ncerr(subr,nf90_put_var(ncid,col_id,linop(i1)%col,(/1,i1/),(/linop(i1)%n_s,1/)))
      call ncerr(subr,nf90_put_var(ncid,S_id,linop(i1)%S,(/1,i1/),(/linop(i1)%n_s,1/)))
   end do
end if

end subroutine linop_write_1d

!----------------------------------------------------------------------
! Subroutine: linop_write_2d
!> Purpose: write array of linear operators to a NetCDF file, 2D
!----------------------------------------------------------------------
subroutine linop_write_2d(ncid,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid               !< NetCDF file id
type(linoptype),intent(in) :: linop(:,:) !< Linear operator

! Local variables
integer :: n1,n2,i1,i2,n_s_max
integer :: n_s_max_id,n1_id,n2_id,n_s_id,n_src_id,n_dst_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_write_2d'

! Array size
n1 = size(linop,1)
n2 = size(linop,2)

! Maximum operator size
n_s_max = 0
do i2=1,n2
   do i1=1,n1
      n_s_max = max(n_s_max,linop(i1,i2)%n_s)
   end do
end do

if ((n1>0).and.(n2>0).and.(n_s_max>0)) then
   ! Start definition mode
   call ncerr(subr,nf90_redef(ncid))

   ! Define dimension
   call ncerr(subr,nf90_def_dim(ncid,trim(linop(1,1)%prefix)//'_n_s_max',n_s_max,n_s_max_id))
   call ncerr(subr,nf90_def_dim(ncid,trim(linop(1,1)%prefix)//'_n1',n1,n1_id))
   call ncerr(subr,nf90_def_dim(ncid,trim(linop(1,1)%prefix)//'_n2',n2,n2_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1,1)%prefix)//'_n_s',nf90_int,(/n1_id,n2_id/),n_s_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1,1)%prefix)//'_n_src',nf90_int,(/n1_id,n2_id/),n_src_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1,1)%prefix)//'_n_dst',nf90_int,(/n1_id,n2_id/),n_dst_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1,1)%prefix)//'_row',nf90_int,(/n_s_max_id,n1_id,n2_id/),row_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1,1)%prefix)//'_col',nf90_int,(/n_s_max_id,n1_id,n2_id/),col_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1,1)%prefix)//'_S',ncfloat,(/n_s_max_id,n1_id,n2_id/),S_id))

   ! End definition mode
   call ncerr(subr,nf90_enddef(ncid))

   do i2=1,n2
      do i1=1,n1
         ! Put variables
         call ncerr(subr,nf90_put_var(ncid,n_s_id,linop(i1,i2)%n_s,(/i1,i2/)))
         call ncerr(subr,nf90_put_var(ncid,n_src_id,linop(i1,i2)%n_src,(/i1,i2/)))
         call ncerr(subr,nf90_put_var(ncid,n_dst_id,linop(i1,i2)%n_dst,(/i1,i2/)))
         call ncerr(subr,nf90_put_var(ncid,row_id,linop(i1,i2)%row,(/1,i1,i2/),(/linop(i1,i2)%n_s,1,1/)))
         call ncerr(subr,nf90_put_var(ncid,col_id,linop(i1,i2)%col,(/1,i1,i2/),(/linop(i1,i2)%n_s,1,1/)))
         call ncerr(subr,nf90_put_var(ncid,S_id,linop(i1,i2)%S,(/1,i1,i2/),(/linop(i1,i2)%n_s,1,1/)))
      end do
   end do
end if

end subroutine linop_write_2d

end module type_linop
