!----------------------------------------------------------------------
! Module: type_mpl
! Purpose: MPI parameters derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mpl

use iso_c_binding
use iso_fortran_env, only : output_unit
use fckit_log_module, only: fckit_log
use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_kinds, only: kind_real
use type_msv, only: msv_type

implicit none

integer,parameter :: lunit_min=10   ! Minimum unit number
integer,parameter :: lunit_max=1000 ! Maximum unit number
integer,parameter :: ddis = 5       ! Progression display step

type mpl_type
   ! MPI parameters
   integer :: nproc                 ! Number of MPI tasks
   integer :: myproc                ! MPI task index
   integer :: ioproc                ! Main task index
   logical :: main                  ! Main task logical
   integer :: tag                   ! MPI tag
   integer :: nthread               ! Number of OpenMP threads

   ! fckit communicator
   type(fckit_mpi_comm) :: f_comm   ! fckit communicator

   ! Missing values
   type(msv_type) :: msv            ! Missing values

   ! Display parameters
   character(len=1024) :: verbosity ! Verbosity level
   character(len=1024) :: info      ! Info buffer
   character(len=1024) :: trace     ! Trace buffer
   character(len=1024) :: stats     ! Stats buffer
   character(len=1024) :: test      ! Test buffer
   integer :: lunit                 ! Listing unit

   ! Display colors
   character(len=1024) :: black     ! Black color code
   character(len=1024) :: green     ! Green color code
   character(len=1024) :: peach     ! Peach color code
   character(len=1024) :: aqua      ! Aqua color code
   character(len=1024) :: purple    ! Purple color code
   character(len=1024) :: err       ! Error color code
   character(len=1024) :: wng       ! Warning color code

   ! Progression print
   integer :: nprog                 ! Progression array size
   integer :: progint               ! Progression integer
   logical,allocatable :: done(:)   ! Progression array
contains
   procedure :: newunit => mpl_newunit
   procedure :: init => mpl_init
   procedure :: final => mpl_final
   procedure :: flush => mpl_flush
   procedure :: abort => mpl_abort
   procedure :: warning => mpl_warning
   procedure :: prog_init => mpl_prog_init
   procedure :: prog_print => mpl_prog_print
   procedure :: prog_final => mpl_prog_final
   procedure :: ncdimcheck => mpl_ncdimcheck
   procedure :: ncerr => mpl_ncerr
   procedure :: update_tag => mpl_update_tag
   procedure :: bcast => mpl_bcast_string_1d
   procedure :: mpl_dot_prod_1d
   procedure :: mpl_dot_prod_2d
   procedure :: mpl_dot_prod_3d
   procedure :: mpl_dot_prod_4d
   generic :: dot_prod => mpl_dot_prod_1d,mpl_dot_prod_2d,mpl_dot_prod_3d,mpl_dot_prod_4d
   procedure :: split => mpl_split
   procedure :: mpl_share_integer_1d
   procedure :: mpl_share_integer_2d
   procedure :: mpl_share_real_1d
   procedure :: mpl_share_real_4d
   procedure :: mpl_share_logical_1d
   procedure :: mpl_share_logical_3d
   procedure :: mpl_share_logical_4d
   generic :: share => mpl_share_integer_1d,mpl_share_integer_2d,mpl_share_real_1d,mpl_share_real_4d, &
            & mpl_share_logical_1d,mpl_share_logical_3d,mpl_share_logical_4d
   procedure :: glb_to_loc_index => mpl_glb_to_loc_index
   procedure :: mpl_glb_to_loc_real_1d
   procedure :: mpl_glb_to_loc_real_2d
   generic :: glb_to_loc => mpl_glb_to_loc_real_1d,mpl_glb_to_loc_real_2d
   procedure :: mpl_loc_to_glb_real_1d
   procedure :: mpl_loc_to_glb_real_2d
   procedure :: mpl_loc_to_glb_logical_2d
   generic :: loc_to_glb => mpl_loc_to_glb_real_1d,mpl_loc_to_glb_real_2d,mpl_loc_to_glb_logical_2d
   procedure :: mpl_write_integer
   procedure :: mpl_write_integer_array
   procedure :: mpl_write_real
   procedure :: mpl_write_real_array
   procedure :: mpl_write_logical
   procedure :: mpl_write_logical_array
   procedure :: mpl_write_string
   procedure :: mpl_write_string_array
   generic :: write => mpl_write_integer,mpl_write_integer_array,mpl_write_real,mpl_write_real_array, &
            & mpl_write_logical,mpl_write_logical_array,mpl_write_string,mpl_write_string_array
end type mpl_type

private
public :: mpl_type

contains

!----------------------------------------------------------------------
! Subroutine: mpl_newunit
! Purpose: find a free unit
!----------------------------------------------------------------------
subroutine mpl_newunit(mpl,lunit)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(out) :: lunit         ! New unit

! Local variables
integer :: lun
logical :: lopened
character(len=1024),parameter :: subr = 'mpl_newunit'

! Loop over possible units
do lun=lunit_min,lunit_max
   inquire(unit=lun,opened=lopened)
   if (.not.lopened) then
      lunit=lun
      exit
   end if
end do

! Check
if (lopened) call mpl%abort(subr,'cannot find a free unit')

end subroutine mpl_newunit

!----------------------------------------------------------------------
! Subroutine: mpl_init
! Purpose: initialize MPL object
!----------------------------------------------------------------------
subroutine mpl_init(mpl)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data

! Get MPI communicator
mpl%f_comm = fckit_mpi_comm()

! Get MPI size
mpl%nproc = mpl%f_comm%size()

! Get MPI rank
mpl%myproc = mpl%f_comm%rank()+1

! Define main task
mpl%ioproc = 1
mpl%main = (mpl%myproc==mpl%ioproc)

! Time-based tag
if (mpl%main) then
   call system_clock(count=mpl%tag)
   call mpl%update_tag(0)
end if
call mpl%f_comm%broadcast(mpl%tag,mpl%ioproc-1)

! Set max number of OpenMP threads
mpl%nthread = 1
!$ mpl%nthread = omp_get_max_threads()
!$ call omp_set_num_threads(mpl%nthread)

! Set log at 'no message'
mpl%info = 'no_message'
mpl%trace = 'no_message'
mpl%stats = 'no_message'
mpl%test = 'no_message'

! Set default listing
mpl%black = ' '
mpl%green = ' '
mpl%peach = ' '
mpl%aqua = ' '
mpl%purple = ' '
mpl%err = ' '
mpl%wng = ' '
mpl%verbosity = 'all'
mpl%lunit = mpl%msv%vali

end subroutine mpl_init

!----------------------------------------------------------------------
! Subroutine: mpl_final
! Purpose: finalize MPI
!----------------------------------------------------------------------
subroutine mpl_final(mpl)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data

! Release memory
if (allocated(mpl%done)) deallocate(mpl%done)

! Finalize MPI communicator
call mpl%f_comm%final

end subroutine mpl_final

!----------------------------------------------------------------------
! Subroutine: mpl_flush
! Purpose: flush listings
!----------------------------------------------------------------------
subroutine mpl_flush(mpl,advance_flag)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl        ! MPI data
logical,intent(in),optional :: advance_flag ! Advance flag

! Local variables
logical :: ladvance_flag

if ((trim(mpl%verbosity)=='all').or.((trim(mpl%verbosity)=='main').and.mpl%main)) then
   ! Set advance flag
   ladvance_flag = .true.
   if (present(advance_flag)) ladvance_flag = advance_flag

   ! Check info message
   if (trim(mpl%info)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info(trim(mpl%info))
         else
            write(mpl%lunit,'(a)') trim(mpl%info)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info(trim(mpl%info),newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') trim(mpl%info)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%info = 'no_message'
   end if

   ! Check trace message
   if (trim(mpl%trace)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info('OOPS_TRACE: '//trim(mpl%trace))
         else
            write(mpl%lunit,'(a)') 'OOPS_TRACE: '//trim(mpl%trace)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info('OOPS_TRACE: '//trim(mpl%trace),newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') 'OOPS_TRACE: '//trim(mpl%trace)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%trace = 'no_message'
   end if

   ! Check stats message
   if (trim(mpl%stats)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info('OOPS_STATS: '//trim(mpl%stats))
         else
            write(mpl%lunit,'(a)') 'OOPS_STATS: '//trim(mpl%stats)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info('OOPS_STATS: '//trim(mpl%stats),newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') 'OOPS_STATS: '//trim(mpl%stats)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%stats = 'no_message'
   end if

   ! Check test message
   if (trim(mpl%test)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info('Test     : '//trim(mpl%test))
         else
            write(mpl%lunit,'(a)') 'Test     : '//trim(mpl%test)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%isi(mpl%lunit)) then
            call fckit_log%info('Test     : '//trim(mpl%test),newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') 'Test     : '//trim(mpl%test)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%test = 'no_message'
   end if
end if

end subroutine mpl_flush

!----------------------------------------------------------------------
! Subroutine: mpl_abort
! Purpose: clean MPI abort
!----------------------------------------------------------------------
subroutine mpl_abort(mpl,subr,message)

implicit none

! Passed variable
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
character(len=*),intent(in) :: message ! Message

! Write message
write(mpl%info,'(a)') trim(mpl%err)//'!!! Error in '//trim(subr)//': '//trim(message)//trim(mpl%black)

! Flush listing
call mpl%flush

! Write standard output message
write(output_unit,'(a,i4.4,a)') '!!! ABORT in '//trim(subr)//' on task #',mpl%myproc,': '//trim(message)
call flush(output_unit)

! Abort MPI
call mpl%f_comm%abort(1)


end subroutine mpl_abort

!----------------------------------------------------------------------
! Subroutine: mpl_warning
! Purpose: print warning message
!----------------------------------------------------------------------
subroutine mpl_warning(mpl,subr,message)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
character(len=*),intent(in) :: message ! Message

! Print warning message
write(mpl%info,'(a)') trim(mpl%wng)//'!!! Warning in '//trim(subr)//': '//trim(message)//trim(mpl%black)
call mpl%flush

end subroutine mpl_warning

!----------------------------------------------------------------------
! Subroutine: mpl_prog_init
! Purpose: initialize progression display
!----------------------------------------------------------------------
subroutine mpl_prog_init(mpl,nprog)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in) :: nprog          ! Array size

! Print message
write(mpl%info,'(a)') ' 0%'
call mpl%flush(.false.)

! Allocation
allocate(mpl%done(nprog))

! Initialization
mpl%nprog = nprog
mpl%progint = ddis
mpl%done = .false.

end subroutine mpl_prog_init

!----------------------------------------------------------------------
! Subroutine: mpl_prog_print
! Purpose: print progression display
!----------------------------------------------------------------------
subroutine mpl_prog_print(mpl,i)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in),optional :: i     ! Index

! Local variables
integer :: ithread
real(kind_real) :: prog

! Update progression array
if (present(i)) mpl%done(i) = .true.

! Print message
prog = 100.0*real(count(mpl%done),kind_real)/real(mpl%nprog,kind_real)
ithread = 0
!$ ithread = omp_get_thread_num()
if ((int(prog)>mpl%progint).and.(ithread==0)) then
   if (mpl%progint<100) then
      if (mpl%progint<10) then
         write(mpl%info,'(i2,a)') mpl%progint,'% '
      else
         write(mpl%info,'(i3,a)') mpl%progint,'% '
      end if
      call mpl%flush(.false.)
   end if
   mpl%progint = mpl%progint+ddis
end if

end subroutine mpl_prog_print

!----------------------------------------------------------------------
! Subroutine: mpl_prog_final
! Purpose: finalize progression display
!----------------------------------------------------------------------
subroutine mpl_prog_final(mpl,advance_flag)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl        ! MPI data
logical,intent(in),optional :: advance_flag ! Advance flag

! Local variables
logical :: ladvance_flag

! Set advance flag
ladvance_flag = .true.
if (present(advance_flag)) ladvance_flag = advance_flag

! Print message
write(mpl%info,'(a)') ' 100%'
call mpl%flush(ladvance_flag)

! Release memory
deallocate(mpl%done)

end subroutine mpl_prog_final

!----------------------------------------------------------------------
! Subroutine: mpl_ncdimcheck
! Purpose: check if NetCDF file dimension exists and has the right size
!----------------------------------------------------------------------
function mpl_ncdimcheck(mpl,subr,ncid,dimname,dimsize,mode) result (dimid)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
integer,intent(in) :: ncid             ! NetCDF file ID
character(len=*),intent(in) :: dimname ! Dimension name
integer,intent(in) :: dimsize          ! Dimension size
logical,intent(in) :: mode             ! Mode (.true.: create and return id if missing, abort if present but wrong size / .false.: return id if present and has right size, else return missing value)

! Result
integer :: dimid                       ! NetCDF dimension ID

! Local variables
integer :: info,dimsize_test

! End definition mode
if (mode) call mpl%ncerr(subr,nf90_enddef(ncid))

! Get dimension ID
info = nf90_inq_dimid(ncid,trim(dimname),dimid)

if (info==nf90_noerr) then
   ! Get dimension size
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,dimid,len=dimsize_test))

   if (dimsize_test==dimsize) then
      ! Definition mode
      if (mode) call mpl%ncerr(subr,nf90_redef(ncid))
   else
      ! Wrong dimension
      if (mode) then
         ! Abort
         call mpl%abort(subr,'dimension '//trim(dimname)//' has a different size in file')
      else
         ! Return missing value
         dimid = mpl%msv%vali
      end if
   end if
else
   if (mode) then
      ! Definition mode
      call mpl%ncerr(subr,nf90_redef(ncid))

      ! Create dimension
      call mpl%ncerr(subr,nf90_def_dim(ncid,trim(dimname),dimsize,dimid))
   else
      ! Return missing value
      dimid = mpl%msv%vali
   end if
end if

end function mpl_ncdimcheck

!----------------------------------------------------------------------
! Subroutine: mpl_ncerr
! Purpose: handle NetCDF error
!----------------------------------------------------------------------
subroutine mpl_ncerr(mpl,subr,info)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
character(len=*),intent(in) :: subr  ! Calling subroutine
integer,intent(in) :: info           ! Info index

! Check status
if (info/=nf90_noerr) call mpl%abort(subr,trim(nf90_strerror(info)))

end subroutine mpl_ncerr

!----------------------------------------------------------------------
! Subroutine: mpl_update_tag
! Purpose: update MPL tag
!----------------------------------------------------------------------
subroutine mpl_update_tag(mpl,add)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in) :: add            ! Tag update incrememnt

! Update tag
mpl%tag = mpl%tag+add

! Apply bounds (between 1 and 10000)
mpl%tag = mod(mpl%tag,10000)
mpl%tag = max(mpl%tag,1)

end subroutine mpl_update_tag

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_string_1d
! Purpose: broadcast 1d string array
!----------------------------------------------------------------------
subroutine mpl_bcast_string_1d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl                  ! MPI data
character(len=*),dimension(:),intent(inout) :: var ! Logical array, 1d
integer,intent(in) :: root                         ! Root task

! Local variable
integer :: i

! Broadcast one string at a time
do i=1,size(var)
   call mpl%f_comm%broadcast(var(i),root)
end do

end subroutine mpl_bcast_string_1d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_1d
! Purpose: global dot product over local fields, 1d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_1d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl     ! MPI data
real(kind_real),intent(in) :: fld1(:) ! Field 1
real(kind_real),intent(in) :: fld2(:) ! Field 2
real(kind_real),intent(out) :: dp     ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnotr(fld1).and.mpl%msv%isnotr(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%ioproc-1)

end subroutine mpl_dot_prod_1d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_2d
! Purpose: global dot product over local fields, 2d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_2d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl       ! MPI data
real(kind_real),intent(in) :: fld1(:,:) ! Field 1
real(kind_real),intent(in) :: fld2(:,:) ! Field 2
real(kind_real),intent(out) :: dp       ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnotr(fld1).and.mpl%msv%isnotr(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%ioproc-1)

end subroutine mpl_dot_prod_2d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_3d
! Purpose: global dot product over local fields, 3d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_3d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl         ! MPI data
real(kind_real),intent(in) :: fld1(:,:,:) ! Field 1
real(kind_real),intent(in) :: fld2(:,:,:) ! Field 2
real(kind_real),intent(out) :: dp         ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnotr(fld1).and.mpl%msv%isnotr(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%ioproc-1)

end subroutine mpl_dot_prod_3d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_4d
! Purpose: global dot product over local fields, 4d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_4d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl           ! MPI data
real(kind_real),intent(in) :: fld1(:,:,:,:) ! Field 1
real(kind_real),intent(in) :: fld2(:,:,:,:) ! Field 2
real(kind_real),intent(out) :: dp           ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnotr(fld1).and.mpl%msv%isnotr(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%ioproc-1)

end subroutine mpl_dot_prod_4d

!----------------------------------------------------------------------
! Subroutine: mpl_split
! Purpose: split array over different MPI tasks
!----------------------------------------------------------------------
subroutine mpl_split(mpl,n,n_loc)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl         ! MPI data
integer,intent(in) :: n                   ! Total array size
integer,intent(out) :: n_loc(0:mpl%nproc) ! Local array size

! Local variable
integer :: iproc,nres,delta
integer :: i_s,i_e

! MPI splitting
n_loc(0) = 0
nres = n
i_e = 0
do iproc=1,mpl%nproc
   i_s = i_e+1
   delta = n/mpl%nproc
   if (nres>(mpl%nproc-iproc+1)*delta) delta = delta+1
   i_e = i_s+delta-1
   n_loc(iproc) = delta
   nres = nres-delta
end do

end subroutine mpl_split

!----------------------------------------------------------------------
! Subroutine: mpl_share_integer_1d
! Purpose: share integer array over different MPI tasks, 1d
!----------------------------------------------------------------------
subroutine mpl_share_integer_1d(mpl,n1,n1_loc,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n1                  ! Global array size
integer,intent(in) :: n1_loc(0:mpl%nproc) ! Local array size
integer,intent(inout) :: var(n1)          ! Integer array, 1d

! Local variable
integer :: iproc
integer :: sendcount,recvcounts(mpl%nproc),displs(mpl%nproc)
integer,allocatable :: sbuf(:)

! Dimensions
sendcount = n1_loc(mpl%myproc)
recvcounts = n1_loc(1:mpl%nproc)
do iproc=1,mpl%nproc
   displs(iproc) = sum(n1_loc(0:iproc-1))
end do

! Allocation
allocate(sbuf(sendcount))

! Prepare buffer
sbuf = var(sum(n1_loc(0:mpl%myproc-1))+1:sum(n1_loc(0:mpl%myproc)))

! Gather data
call mpl%f_comm%allgather(sbuf,var,sendcount,recvcounts,displs)

end subroutine mpl_share_integer_1d

!----------------------------------------------------------------------
! Subroutine: mpl_share_integer_2d
! Purpose: share integer array over different MPI tasks, 2d
!----------------------------------------------------------------------
subroutine mpl_share_integer_2d(mpl,n1,n2,n1_loc,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n1                  ! Global array size 1
integer,intent(in) :: n2                  ! Global array size 2
integer,intent(in) :: n1_loc(0:mpl%nproc) ! Local array size 1
integer,intent(inout) :: var(n1,n2)       ! Integer array, 2d

! Local variable
integer :: iproc,i1,i2,i1_loc,i
integer :: sendcount,recvcounts(mpl%nproc),displs(mpl%nproc)
integer,allocatable :: rbuf(:),sbuf(:)

! Dimensions
sendcount = n1_loc(mpl%myproc)*n2
recvcounts = n1_loc(1:mpl%nproc)*n2
do iproc=1,mpl%nproc
   displs(iproc) = sum(n1_loc(0:iproc-1))*n2
end do

! Allocation
allocate(sbuf(sendcount))
allocate(rbuf(n1*n2))

! Prepare buffer
i = 0
do i2=1,n2
   do i1_loc=1,n1_loc(mpl%myproc)
      i1 = sum(n1_loc(0:mpl%myproc-1))+i1_loc
      i = i+1
      sbuf(i) = var(i1,i2)
   end do
end do

! Gather data
call mpl%f_comm%allgather(sbuf,rbuf,sendcount,recvcounts,displs)

! Format data
i = 0
do iproc=1,mpl%nproc
   do i2=1,n2
      do i1_loc=1,n1_loc(iproc)
         i1 = sum(n1_loc(0:iproc-1))+i1_loc
         i = i+1
         var(i1,i2) = rbuf(i)
      end do
   end do
end do

end subroutine mpl_share_integer_2d

!----------------------------------------------------------------------
! Subroutine: mpl_share_real_1d
! Purpose: share real array over different MPI tasks, 1d
!----------------------------------------------------------------------
subroutine mpl_share_real_1d(mpl,n1,n1_loc,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n1                  ! Global array size
integer,intent(in) :: n1_loc(0:mpl%nproc) ! Local array size
real(kind_real),intent(inout) :: var(n1)  ! Real array, 1d

! Local variable
integer :: iproc
integer :: sendcount,recvcounts(mpl%nproc),displs(mpl%nproc)
real(kind_real),allocatable :: sbuf(:)

! Dimensions
sendcount = n1_loc(mpl%myproc)
recvcounts = n1_loc(1:mpl%nproc)
do iproc=1,mpl%nproc
   displs(iproc) = sum(n1_loc(0:iproc-1))
end do

! Allocation
allocate(sbuf(sendcount))

! Prepare buffer
sbuf = var(sum(n1_loc(0:mpl%myproc-1))+1:sum(n1_loc(0:mpl%myproc)))

! Gather data
call mpl%f_comm%allgather(sbuf,var,sendcount,recvcounts,displs)

end subroutine mpl_share_real_1d

!----------------------------------------------------------------------
! Subroutine: mpl_share_real_4d
! Purpose: share real array over different MPI tasks, 4d
!----------------------------------------------------------------------
subroutine mpl_share_real_4d(mpl,n1,n2,n3,n4,n1_loc,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl              ! MPI data
integer,intent(in) :: n1                          ! Array size 1
integer,intent(in) :: n2                          ! Array size 2
integer,intent(in) :: n3                          ! Array size 3
integer,intent(in) :: n4                          ! Array size 4
integer,intent(in) :: n1_loc(0:mpl%nproc)         ! Local array size 1
real(kind_real),intent(inout) :: var(n1,n2,n3,n4) ! Real array, 4d

! Local variable
integer :: iproc,i1,i2,i3,i4,i1_loc,i
integer :: sendcount,recvcounts(mpl%nproc),displs(mpl%nproc)
real(kind_real),allocatable :: rbuf(:),sbuf(:)

! Dimensions
sendcount = n1_loc(mpl%myproc)*n2*n3*n4
recvcounts = n1_loc(1:mpl%nproc)*n2*n3*n4
do iproc=1,mpl%nproc
   displs(iproc) = sum(n1_loc(0:iproc-1))*n2*n3*n4
end do

! Allocation
allocate(sbuf(sendcount))
allocate(rbuf(n1*n2*n3*n4))

! Prepare buffer
i = 0
do i4=1,n4
   do i3=1,n3
      do i2=1,n2
         do i1_loc=1,n1_loc(mpl%myproc)
            i1 = sum(n1_loc(0:mpl%myproc-1))+i1_loc
            i = i+1
            sbuf(i) = var(i1,i2,i3,i4)
         end do
      end do
   end do
end do

! Gather data
call mpl%f_comm%allgather(sbuf,rbuf,sendcount,recvcounts,displs)

! Format data
i = 0
do iproc=1,mpl%nproc
   do i4=1,n4
      do i3=1,n3
         do i2=1,n2
            do i1_loc=1,n1_loc(iproc)
               i1 = sum(n1_loc(0:iproc-1))+i1_loc
               i = i+1
               var(i1,i2,i3,i4) = rbuf(i)
            end do
         end do
      end do
   end do
end do

end subroutine mpl_share_real_4d

!----------------------------------------------------------------------
! Subroutine: mpl_share_logical_1d
! Purpose: share logical array over different MPI tasks, 1d
!----------------------------------------------------------------------
subroutine mpl_share_logical_1d(mpl,n1,n1_loc,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n1                  ! Global array size
integer,intent(in) :: n1_loc(0:mpl%nproc) ! Local array size
logical,intent(inout) :: var(n1)          ! Logical array, 1d

! Local variable
integer :: iproc,i1_loc,i1
integer :: sendcount,recvcounts(mpl%nproc),displs(mpl%nproc)
integer,allocatable :: sbuf(:),rbuf(:)
character(len=1024),parameter :: subr = 'mpl_share_logical_1d'

! Dimensions
sendcount = n1_loc(mpl%myproc)
recvcounts = n1_loc(1:mpl%nproc)
do iproc=1,mpl%nproc
   displs(iproc) = sum(n1_loc(0:iproc-1))
end do

! Allocation
allocate(sbuf(sendcount))
allocate(rbuf(n1))

! Prepare buffer
do i1_loc=1,n1_loc(mpl%myproc)
   i1 = sum(n1_loc(0:mpl%myproc-1))+i1_loc
   if (var(i1)) then
      sbuf(i1_loc) = 1
   else
      sbuf(i1_loc) = 0
   end if
end do

! Gather data
call mpl%f_comm%allgather(sbuf,rbuf,sendcount,recvcounts,displs)

! Format data
do i1=1,n1
   if (rbuf(i1)==0) then
      var(i1) = .false.
   elseif (rbuf(i1)==1) then
      var(i1) = .true.
   else
      call mpl%abort(subr,'wrong value for received buffer in mpl_share_logical_1d')
   end if
end do

end subroutine mpl_share_logical_1d

!----------------------------------------------------------------------
! Subroutine: mpl_share_logical_3d
! Purpose: share logical array over different MPI tasks, 3d
!----------------------------------------------------------------------
subroutine mpl_share_logical_3d(mpl,n1,n2,n3,n1_loc,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n1                  ! Array size 1
integer,intent(in) :: n2                  ! Array size 2
integer,intent(in) :: n3                  ! Array size 3
integer,intent(in) :: n1_loc(0:mpl%nproc) ! Local array size 1
logical,intent(inout) :: var(n1,n2,n3)    ! Logical array, 3d

! Local variable
integer :: iproc,i1,i2,i3,i1_loc,i
integer :: sendcount,recvcounts(mpl%nproc),displs(mpl%nproc)
integer,allocatable :: rbuf(:),sbuf(:)
character(len=1024),parameter :: subr = 'mpl_share_logical_3d'

! Dimensions
sendcount = n1_loc(mpl%myproc)*n2*n3
recvcounts = n1_loc(1:mpl%nproc)*n2*n3
do iproc=1,mpl%nproc
   displs(iproc) = sum(n1_loc(0:iproc-1))*n2*n3
end do

! Allocation
allocate(sbuf(sendcount))
allocate(rbuf(n1*n2*n3))

! Prepare buffer
i = 0
do i3=1,n3
   do i2=1,n2
      do i1_loc=1,n1_loc(mpl%myproc)
         i1 = sum(n1_loc(0:mpl%myproc-1))+i1_loc
         i = i+1
         if (var(i1,i2,i3)) then
            sbuf(i) = 1
         else
            sbuf(i) = 0
         end if
      end do
   end do
end do

! Gather data
call mpl%f_comm%allgather(sbuf,rbuf,sendcount,recvcounts,displs)

! Format data
i = 0
do iproc=1,mpl%nproc
   do i3=1,n3
      do i2=1,n2
         do i1_loc=1,n1_loc(iproc)
            i1 = sum(n1_loc(0:iproc-1))+i1_loc
            i = i+1
            if (rbuf(i)==0) then
               var(i1,i2,i3) = .false.
            elseif (rbuf(i)==1) then
               var(i1,i2,i3) = .true.
            else
               call mpl%abort(subr,'wrong value for received buffer in mpl_share_logical_3d')
            end if
         end do
      end do
   end do
end do

end subroutine mpl_share_logical_3d

!----------------------------------------------------------------------
! Subroutine: mpl_share_logical_4d
! Purpose: share logical array over different MPI tasks, 4d
!----------------------------------------------------------------------
subroutine mpl_share_logical_4d(mpl,n1,n2,n3,n4,n1_loc,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n1                  ! Array size 1
integer,intent(in) :: n2                  ! Array size 2
integer,intent(in) :: n3                  ! Array size 3
integer,intent(in) :: n4                  ! Array size 4
integer,intent(in) :: n1_loc(0:mpl%nproc) ! Local array size 1
logical,intent(inout) :: var(n1,n2,n3,n4) ! Logical array, 4d

! Local variable
integer :: iproc,i1,i2,i3,i4,i1_loc,i
integer :: sendcount,recvcounts(mpl%nproc),displs(mpl%nproc)
integer,allocatable :: rbuf(:),sbuf(:)
character(len=1024),parameter :: subr = 'mpl_share_logical_4d'

! Dimensions
sendcount = n1_loc(mpl%myproc)*n2*n3*n4
recvcounts = n1_loc(1:mpl%nproc)*n2*n3*n4
do iproc=1,mpl%nproc
   displs(iproc) = sum(n1_loc(0:iproc-1))*n2*n3*n4
end do

! Allocation
allocate(sbuf(sendcount))
allocate(rbuf(n1*n2*n3*n4))

! Prepare buffer
i = 0
do i4=1,n4
   do i3=1,n3
      do i2=1,n2
         do i1_loc=1,n1_loc(mpl%myproc)
            i1 = sum(n1_loc(0:mpl%myproc-1))+i1_loc
            i = i+1
            if (var(i1,i2,i3,i4)) then
               sbuf(i) = 1
            else
               sbuf(i) = 0
            end if
         end do
      end do
   end do
end do

! Gather data
call mpl%f_comm%allgather(sbuf,rbuf,sendcount,recvcounts,displs)

! Format data
i = 0
do iproc=1,mpl%nproc
   do i4=1,n4
      do i3=1,n3
         do i2=1,n2
            do i1_loc=1,n1_loc(iproc)
               i1 = sum(n1_loc(0:iproc-1))+i1_loc
               i = i+1
               if (rbuf(i)==0) then
                  var(i1,i2,i3,i4) = .false.
               elseif (rbuf(i)==1) then
                  var(i1,i2,i3,i4) = .true.
               else
                  call mpl%abort(subr,'wrong value for received buffer in mpl_share_logical_4d')
               end if
            end do
         end do
      end do
   end do
end do

end subroutine mpl_share_logical_4d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_index
! Purpose: communicate global index to local index
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_index(mpl,n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl               ! MPI data
integer,intent(in) :: n_loc                        ! Local dimension
integer,intent(in) :: loc_to_glb(n_loc)            ! Local to global index
integer,intent(in) :: n_glb                        ! Global dimension
integer,intent(out) :: glb_to_loc(n_glb)           ! Global to local index
integer,intent(out),optional :: glb_to_proc(n_glb) ! Global to processor

! Local variables
integer :: iproc,i_loc,n_loc_tmp
integer,allocatable :: loc_to_glb_tmp(:)
type(fckit_mpi_status) :: status

if (mpl%main) then
   ! Initialization
   glb_to_loc = mpl%msv%vali
   if (present(glb_to_proc)) glb_to_proc = mpl%msv%vali

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         n_loc_tmp = n_loc
      else
         ! Receive dimension on ioproc
         call mpl%f_comm%receive(n_loc_tmp,iproc-1,mpl%tag,status)
      end if

      ! Allocation
      allocate(loc_to_glb_tmp(n_loc_tmp))

      if (iproc==mpl%ioproc) then
         ! Copy data
         loc_to_glb_tmp = loc_to_glb
      else
         ! Receive data on ioproc
         call mpl%f_comm%receive(loc_to_glb_tmp,iproc-1,mpl%tag+1,status)
      end if

      ! Fill glb_to_loc and glb_to_proc if required
      do i_loc=1,n_loc_tmp
         glb_to_loc(loc_to_glb_tmp(i_loc)) = i_loc
         if (present(glb_to_proc)) glb_to_proc(loc_to_glb_tmp(i_loc)) = iproc
      end do

      ! Release memory
      deallocate(loc_to_glb_tmp)
   end do
else
   ! Send dimensions to ioproc
   call mpl%f_comm%send(n_loc,mpl%ioproc-1,mpl%tag)

   ! Send data to ioproc
   call mpl%f_comm%send(loc_to_glb,mpl%ioproc-1,mpl%tag+1)
end if
call mpl%update_tag(2)

! Broadcast
call mpl%f_comm%broadcast(glb_to_loc,mpl%ioproc-1)
if (present(glb_to_proc)) call mpl%f_comm%broadcast(glb_to_proc,mpl%ioproc-1)

end subroutine mpl_glb_to_loc_index

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_real_1d
! Purpose: global to local, 1d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_real_1d(mpl,n_glb,glb_to_proc,glb_to_loc,glb,n_loc,loc)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n_glb               ! Global array size
integer,intent(in) :: glb_to_proc(n_glb)  ! Global index to task index
integer,intent(in) :: glb_to_loc(n_glb)   ! Global index to local index
real(kind_real),intent(in) :: glb(:)      ! Global array
integer,intent(in) :: n_loc               ! Local array size
real(kind_real),intent(out) :: loc(n_loc) ! Local array
type(fckit_mpi_status) :: status

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
real(kind_real),allocatable :: sbuf(:)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_real_1d'

! Check global array size
if (mpl%main) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_glb_to_loc_real_1d')
end if

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(sbuf(n_loc_tmp))

      ! Prepare buffers
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            sbuf(i_loc) = glb(i_glb)
         end if
      end do

      if (iproc==mpl%ioproc) then
         ! Copy data
         loc = sbuf
      else
         ! Send data to iproc
         call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Receive data from ioproc
   call mpl%f_comm%receive(loc,mpl%ioproc-1,mpl%tag,status)
end if
call mpl%update_tag(1)

end subroutine mpl_glb_to_loc_real_1d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_real_2d
! Purpose: global to local, 2d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_real_2d(mpl,nl,n_glb,glb_to_proc,glb_to_loc,glb,n_loc,loc)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl         ! MPI data
integer,intent(in) :: nl                     ! Number of levels
integer,intent(in) :: n_glb                  ! Global array size
integer,intent(in) :: glb_to_proc(n_glb)     ! Global index to task index
integer,intent(in) :: glb_to_loc(n_glb)      ! Global index to local index
real(kind_real),intent(in) :: glb(:,:)       ! Global array
integer,intent(in) :: n_loc                  ! Local array size
real(kind_real),intent(out) :: loc(n_loc,nl) ! Local array

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il
real(kind_real),allocatable :: sbuf(:),rbuf(:)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_real_2d'
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main) then
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_glb_to_loc_real_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_glb_to_loc_real_2d')
end if

! Allocation
allocate(rbuf(n_loc*nl))

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(sbuf(n_loc_tmp*nl))

      ! Prepare buffers
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            do il=1,nl
               sbuf((il-1)*n_loc_tmp+i_loc) = glb(i_glb,il)
            end do
         end if
      end do

      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf
      else
         ! Send data to iproc
         call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Receive data from ioproc
   call mpl%f_comm%receive(rbuf,mpl%ioproc-1,mpl%tag,status)
end if
call mpl%update_tag(1)

! Unpack buffer
do il=1,nl
   do i_loc=1,n_loc
      loc(i_loc,il) = rbuf((il-1)*n_loc+i_loc)
   end do
end do

! Release memory
deallocate(rbuf)

end subroutine mpl_glb_to_loc_real_2d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_real_1d
! Purpose: local to global, 1d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_real_1d(mpl,n_loc,loc,n_glb,glb_to_proc,glb_to_loc,bcast,glb)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
integer,intent(in) :: n_loc              ! Local array size
real(kind_real),intent(in) :: loc(n_loc) ! Local array
integer,intent(in) :: n_glb              ! Global array size
integer,intent(in) :: glb_to_proc(n_glb) ! Global index to task index
integer,intent(in) :: glb_to_loc(n_glb)  ! Global index to local index
logical,intent(in) :: bcast              ! Broadcast option
real(kind_real),intent(out) :: glb(:)    ! Global array

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
real(kind_real),allocatable :: rbuf(:)
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_real_1d'
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main.or.bcast) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_loc_to_glb_real_1d')
end if

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp))

      if (iproc==mpl%ioproc) then
          ! Copy data
          rbuf = loc
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            glb(i_glb) = rbuf(i_loc)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to ioproc
   call mpl%f_comm%send(loc,mpl%ioproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (bcast) call mpl%f_comm%broadcast(glb,mpl%ioproc-1)

end subroutine mpl_loc_to_glb_real_1d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_real_2d
! Purpose: local to global, 2d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_real_2d(mpl,nl,n_loc,loc,n_glb,glb_to_proc,glb_to_loc,bcast,glb)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl        ! MPI data
integer,intent(in) :: nl                    ! Number of levels
integer,intent(in) :: n_loc                 ! Local array size
real(kind_real),intent(in) :: loc(n_loc,nl) ! Local array
integer,intent(in) :: n_glb                 ! Global array size
integer,intent(in) :: glb_to_proc(n_glb)    ! Global index to task index
integer,intent(in) :: glb_to_loc(n_glb)     ! Global index to local index
logical,intent(in) :: bcast                 ! Broadcast option
real(kind_real),intent(out) :: glb(:,:)     ! Global array

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il
real(kind_real),allocatable :: rbuf(:),sbuf(:)
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_real_2d'
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main.or.bcast) then
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_loc_to_glb_real_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_loc_to_glb_real_1d')
end if

! Allocation
allocate(sbuf(n_loc*nl))

! Prepare buffer
do il=1,nl
   do i_loc=1,n_loc
      sbuf((il-1)*n_loc+i_loc) = loc(i_loc,il)
   end do
end do

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp*nl))

      if (iproc==mpl%ioproc) then
          ! Copy data
          rbuf = sbuf
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            do il=1,nl
               glb(i_glb,il) = rbuf((il-1)*n_loc_tmp+i_loc)
            end do
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to ioproc
   call mpl%f_comm%send(sbuf,mpl%ioproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (bcast) call mpl%f_comm%broadcast(glb,mpl%ioproc-1)

! Release memory
deallocate(sbuf)

end subroutine mpl_loc_to_glb_real_2d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_logical_2d
! Purpose: local to global for a logical, 2d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_logical_2d(mpl,nl,n_loc,loc,n_glb,glb_to_proc,glb_to_loc,bcast,glb)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
integer,intent(in) :: nl                 ! Number of levels
integer,intent(in) :: n_loc              ! Local array size
logical,intent(in) :: loc(n_loc,nl)      ! Local array
integer,intent(in) :: n_glb              ! Global array size
integer,intent(in) :: glb_to_proc(n_glb) ! Global index to task index
integer,intent(in) :: glb_to_loc(n_glb)  ! Global index to local index
logical,intent(in) :: bcast              ! Broadcast option
logical,intent(out) :: glb(:,:)          ! Global array

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il
logical,allocatable :: rbuf(:),sbuf(:)
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_logical_2d'
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main.or.bcast) then
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_loc_to_glb_real_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_loc_to_glb_real_1d')
end if

! Allocation
allocate(sbuf(n_loc*nl))

! Prepare buffer
do il=1,nl
   do i_loc=1,n_loc
      sbuf((il-1)*n_loc+i_loc) = loc(i_loc,il)
   end do
end do

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp*nl))

      if (iproc==mpl%ioproc) then
          ! Copy data
          rbuf = sbuf
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            do il=1,nl
               glb(i_glb,il) = rbuf((il-1)*n_loc_tmp+i_loc)
            end do
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to ioproc
   call mpl%f_comm%send(sbuf,mpl%ioproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (bcast) call mpl%f_comm%broadcast(glb,mpl%ioproc-1)

! Release memory
deallocate(sbuf)

end subroutine mpl_loc_to_glb_logical_2d

!----------------------------------------------------------------------
! Subroutine: mpl_write_integer
! Purpose: write integer into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_integer(mpl,ncid,varname,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
integer,intent(in) :: var              ! Integer

! Local variables
integer :: delta
character(len=1024) :: str
character(len=1024),parameter :: subr = 'mpl_write_integer'

if (mpl%msv%isi(ncid)) then
   ! Write integer into a log file
   delta = 2
   if (var<0) delta = delta+1
   if (abs(var)>0) then
      write(str,'(a,i4.4,a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a,i', &
       & floor(log(abs(real(var,kind_real)))/log(10.0))+delta,',a)'
   else
      write(str,'(a,i4.4,a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a,i',delta,',a)'
   end if
   write(mpl%info,str) '',trim(varname),'',':',var
   call mpl%flush
else
   ! Write integer into a NetCDF file
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))
end if

end subroutine mpl_write_integer

!----------------------------------------------------------------------
! Subroutine: mpl_write_integer_array
! Purpose: write integer array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_integer_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
integer,intent(in) :: n                ! Integer array size
integer,intent(in) :: var(n)           ! Integer array

! Local variables
integer :: i,delta
character(len=1024) :: str,fullstr
character(len=1024),parameter :: subr = 'mpl_write_integer_array'

if (mpl%msv%isi(ncid)) then
   ! Write integer array into a log file
   write(str,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a)'
   write(mpl%info,str) '',trim(varname),'',':'
   call mpl%flush(.false.)
   do i=1,n
      delta = 2
      if (var(i)<0) delta = delta+1
      if (abs(var(i))>0) then
         write(str,'(a,i4.4,a)') '(i',floor(log(abs(real(var(i),kind_real)))/log(10.0))+delta,',a)'
      else
         write(str,'(a,i4.4,a)') '(i',delta,',a)'
      end if
      write(mpl%info,str) var(i),','
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write integer array into a NetCDF file
   if (n>0) then
      write(fullstr,'(i3.3)') var(1)
      do i=2,n
         write(str,'(i3.3)') var(i)
         fullstr = trim(fullstr)//':'//trim(str)
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
   end if
end if

end subroutine mpl_write_integer_array

!----------------------------------------------------------------------
! Subroutine: mpl_write_real
! Purpose: write real into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_real(mpl,ncid,varname,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
real(kind_real),intent(in) :: var      ! Real

! Local variables
integer :: delta
character(len=1024) :: str
character(len=1024),parameter :: subr = 'mpl_write_real'

if (mpl%msv%isi(ncid)) then
   ! Write real into a log file
   delta = 10
   if (var<0.0) delta = delta+1
   write(str,'(a,i4.4,a,i4.4,a,i2,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a,e',delta,'.3,a)'
   write(mpl%info,str) '',trim(varname),'',':',var
   call mpl%flush
else
   ! Write real into a NetCDF file
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))
end if

end subroutine mpl_write_real

!----------------------------------------------------------------------
! Subroutine: mpl_write_real_array
! Purpose: write real array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_real_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
integer,intent(in) :: n                ! Real array size
real(kind_real),intent(in) :: var(n)   ! Real array

! Local variables
integer :: i,delta
character(len=1024) :: str,fullstr
character(len=1024),parameter :: subr = 'mpl_write_real_array'

if (mpl%msv%isi(ncid)) then
   ! Write real array into a log file
   write(str,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a)'
   write(mpl%info,str) '',trim(varname),'',':'
   call mpl%flush(.false.)
   do i=1,n
      delta = 10
      if (var(i)<0.0) delta = delta+1
      write(str,'(a,i2,a)') '(e',delta,'.3,a)'
      write(mpl%info,str) var(i),','
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write real array into a NetCDF file
   if (n>0) then
      write(fullstr,'(e10.3)') var(1)
      do i=2,n
         write(str,'(e10.3)') var(i)
         fullstr = trim(fullstr)//':'//trim(str)
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
   end if
end if

end subroutine mpl_write_real_array

!----------------------------------------------------------------------
! Subroutine: mpl_write_logical
! Purpose: write logical into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_logical(mpl,ncid,varname,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
logical,intent(in) :: var              ! Logical

! Local variables
character(len=1024) :: str
character(len=1024),parameter :: subr = 'mpl_write_logical'

if (mpl%msv%isi(ncid)) then
   ! Write logical into a log file
   write(str,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a,l2,a)'
   write(mpl%info,str) '',trim(varname),'',':',var
   call mpl%flush
else
   ! Write logical into a NetCDF file
   if (var) then
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.true.'))
   else
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.false.'))
   end if
end if

end subroutine mpl_write_logical

!----------------------------------------------------------------------
! Subroutine: mpl_write_logical_array
! Purpose: write logical array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_logical_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
integer,intent(in) :: n                ! Real array size
logical,intent(in) :: var(n)           ! Logical array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024),parameter :: subr = 'mpl_write_logical_array'

if (mpl%msv%isi(ncid)) then
   ! Write logical array into a log file
   write(str,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a)'
   write(mpl%info,str) '',trim(varname),'',':'
   call mpl%flush(.false.)
   do i=1,n
      write(str,'(a)') '(l2,a)'
      write(mpl%info,str) var(i),','
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write real array into a NetCDF file
   if (n>0) then
      if (var(1)) then
         write(fullstr,'(a6)') '.true.'
      else
         write(fullstr,'(a7)') '.false.'
      end if
      do i=2,n
         if (var(i)) then
            write(str,'(a6)') '.true.'
         else
            write(str,'(a7)') '.false.'
         end if
         fullstr = trim(fullstr)//':'//trim(str)
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
   end if
end if

end subroutine mpl_write_logical_array

!----------------------------------------------------------------------
! Subroutine: mpl_write_string
! Purpose: write string into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_string(mpl,ncid,varname,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
character(len=*),intent(in) :: var     ! String

! Local variables
character(len=1024) :: str
character(len=1024),parameter :: subr = 'mpl_write_string'

if (mpl%msv%isi(ncid)) then
   ! Write string into a log file
   if (len_trim(var)>0) then
      write(str,'(a,i4.4,a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a,a1,a',len_trim(var),')'
      write(mpl%info,str) '',trim(varname),'',':','',trim(var)
   else
      write(str,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a)'
      write(mpl%info,str) '',trim(varname),'',':'
   end if
   call mpl%flush
else
   ! Write string into a NetCDF file
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(var)))
end if

end subroutine mpl_write_string

!----------------------------------------------------------------------
! Subroutine: mpl_write_string_array
! Purpose: write string array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_string_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: varname ! Variable name
integer,intent(in) :: n                ! String array size
character(len=*),intent(in) :: var(n)  ! String array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024),parameter :: subr = 'mpl_write_string_array'

if (mpl%msv%isi(ncid)) then
   ! Write string array into a log file
   write(str,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(varname),',a',25-len_trim(varname),',a)'
   write(mpl%info,str) '',trim(varname),'',':'
   call mpl%flush(.false.)
   do i=1,n
      if (len_trim(var(i))>0) then
         write(str,'(a,i4.4,a)') '(a1,a',len_trim(var(i)),',a)'
         write(mpl%info,str) '',trim(var(i)),','
      else
         write(mpl%info,'(a1,a)') '',','
      end if
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write string array into a NetCDF file
   if (n>0) then
      fullstr = trim(var(1))
      do i=2,n
         fullstr = trim(fullstr)//':'//trim(var(i))
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
   end if
end if

end subroutine mpl_write_string_array

end module type_mpl
