!----------------------------------------------------------------------
! Module: type_mpl
! Purpose: MPI parameters derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mpl

use iso_fortran_env, only : output_unit
use iso_c_binding
use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_kinds, only: kind_real
use type_msv, only: msv_type

implicit none

integer,parameter :: lunit_min=10                  ! Minimum unit number
integer,parameter :: lunit_max=1000                ! Maximum unit number
integer,parameter :: ddis = 5                      ! Progression display step

type mpl_type
   ! MPI parameters
   integer :: nproc                 ! Number of MPI tasks
   integer :: myproc                ! MPI task index
   integer :: ioproc                ! Main task index
   logical :: main                  ! Main task logical
   character(len=1024) :: verbosity ! Verbosity level
   character(len=1024) :: info      ! Info buffer
   character(len=1024) :: test      ! Test buffer
   integer :: info_unit             ! Info listing unit
   integer :: test_unit             ! Test listing unit
   integer :: tag                   ! MPI tag
   integer :: nthread               ! Number of OpenMP threads

   type(fckit_mpi_comm) :: f_comm   ! Interface to fckit

   ! Progression print
   integer :: nprog                 ! Progression array size
   integer :: progint               ! Progression integer
   logical,allocatable :: done(:)   ! Progression array

   ! Display colors
   character(len=1024) :: black     ! Black color code
   character(len=1024) :: green     ! Green color code
   character(len=1024) :: peach     ! Peach color code
   character(len=1024) :: aqua      ! Aqua color code
   character(len=1024) :: purple    ! Purple color code
   character(len=1024) :: err       ! Error color code
   character(len=1024) :: wng       ! Warning color code

   character(len=1024) :: vunitchar ! Vertical unit

   type(msv_type) :: msv            ! Missing values
contains
   procedure :: newunit => mpl_newunit
   procedure :: init => mpl_init
   procedure :: final => mpl_final
   procedure :: init_listing => mpl_init_listing
   procedure :: flush => mpl_flush
   procedure :: close_listing => mpl_close_listing
   procedure :: abort => mpl_abort
   procedure :: warning => mpl_warning
   procedure :: prog_init => mpl_prog_init
   procedure :: prog_print => mpl_prog_print
   procedure :: prog_final => mpl_prog_final
   procedure :: ncerr => mpl_ncerr
   procedure :: update_tag => mpl_update_tag
   procedure :: bcast => mpl_bcast_string_1d
   procedure :: mpl_dot_prod_1d
   procedure :: mpl_dot_prod_2d
   procedure :: mpl_dot_prod_3d
   procedure :: mpl_dot_prod_4d
   generic :: dot_prod => mpl_dot_prod_1d,mpl_dot_prod_2d,mpl_dot_prod_3d,mpl_dot_prod_4d
   procedure :: split => mpl_split
   procedure :: glb_to_loc_index => mpl_glb_to_loc_index
   procedure :: mpl_glb_to_loc_real_1d
   procedure :: mpl_glb_to_loc_real_2d
   generic :: glb_to_loc => mpl_glb_to_loc_real_1d,mpl_glb_to_loc_real_2d
   procedure :: mpl_loc_to_glb_real_1d
   procedure :: mpl_loc_to_glb_real_2d
   procedure :: mpl_loc_to_glb_logical_2d
   generic :: loc_to_glb => mpl_loc_to_glb_real_1d,mpl_loc_to_glb_real_2d,mpl_loc_to_glb_logical_2d
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

! Loop over possible units
do lun=lunit_min,lunit_max
   inquire(unit=lun,opened=lopened)
   if (.not.lopened) then
      lunit=lun
      exit
   end if
end do

! Check
if (lopened) call mpl%abort('cannot find a free unit')

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

end subroutine mpl_init

!----------------------------------------------------------------------
! Subroutine: mpl_final
! Purpose: finalize MPI
!----------------------------------------------------------------------
subroutine mpl_final(mpl)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data

! Finalize MPI communicator
call mpl%f_comm%final

end subroutine mpl_final

!----------------------------------------------------------------------
! Subroutine: mpl_init_listing
! Purpose: initialize listings
!----------------------------------------------------------------------
subroutine mpl_init_listing(mpl,prefix,model,verbosity,colorlog,logpres,lunit)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
character(len=*),intent(in) :: prefix    ! Output prefix
character(len=*),intent(in) :: model     ! Model
character(len=*),intent(in) :: verbosity ! Verbosity level
logical,intent(in) :: colorlog           ! Color listing flag
logical,intent(in) :: logpres            ! Vertical unit flag
integer,intent(in),optional :: lunit     ! Main listing unit

! Local variables
integer :: iproc
character(len=1024) :: filename

! Set verbosity level
mpl%verbosity = trim(verbosity)

! Setup display colors
if (colorlog) then
   mpl%black = char(27)//'[0;0m'
   mpl%green = char(27)//'[0;32m'
   mpl%peach = char(27)//'[1;91m'
   mpl%aqua = char(27)//'[1;36m'
   mpl%purple = char(27)//'[1;35m'
   mpl%err = char(27)//'[0;37;41;1m'
   mpl%wng = char(27)//'[0;37;42;1m'
else
   mpl%black = ' '
   mpl%green = ' '
   mpl%peach = ' '
   mpl%aqua = ' '
   mpl%purple = ' '
   mpl%err = ' '
   mpl%wng = ' '
end if

! Vertical unit
if (trim(model)=='online') then
   mpl%vunitchar = 'vert. unit'
else
   if (logpres) then
      mpl%vunitchar = 'log(Pa)'
   else
      mpl%vunitchar = 'lev.'
   end if
end if

! Define info and test units, and open files
do iproc=1,mpl%nproc
   if ((trim(mpl%verbosity)=='all').or.((trim(mpl%verbosity)=='main').and.(iproc==mpl%ioproc))) then
      ! Info listing
      if (mpl%main.and.present(lunit)) then
         ! Specific listing unit
         mpl%info_unit = lunit
      else
         ! Deal with each proc sequentially
         if (iproc==mpl%myproc) then
            ! Find a free unit
            call mpl%newunit(mpl%info_unit)

            ! Open listing file
            write(filename,'(a,i4.4)') trim(prefix)//'.info.',mpl%myproc-1
            open(unit=mpl%info_unit,file=trim(filename),action='write',status='replace')
         end if
      end if

      ! Test listing
      if (iproc==mpl%myproc) then
         ! Find a free unit
         call mpl%newunit(mpl%test_unit)

         ! Open listing file
         write(filename,'(a,i4.4)') trim(prefix)//'.test.',mpl%myproc-1
         open(unit=mpl%test_unit,file=trim(filename),action='write',status='replace')
      end if

      ! Wait
      call mpl%f_comm%barrier
   end if
end do

! Initialize message
mpl%info = 'no_message'
mpl%test = 'no_message'

end subroutine mpl_init_listing

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
         write(mpl%info_unit,'(a)') trim(mpl%info)
      else
         write(mpl%info_unit,'(a)',advance='no') trim(mpl%info)
      end if
      call flush(mpl%info_unit)

      ! Set at 'no message'
      mpl%info = 'no_message'
   end if

   ! Check test message
   if (trim(mpl%test)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         write(mpl%test_unit,'(a)') trim(mpl%test)
      else
         write(mpl%test_unit,'(a)',advance='no') trim(mpl%test)
      end if
      call flush(mpl%test_unit)

      ! Set at 'no message'
      mpl%test = 'no_message'
   end if
end if

end subroutine mpl_flush

!----------------------------------------------------------------------
! Subroutine: mpl_close_listing
! Purpose: close listings
!----------------------------------------------------------------------
subroutine mpl_close_listing(mpl)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl   ! MPI data

if ((trim(mpl%verbosity)=='all').or.((trim(mpl%verbosity)=='main').and.mpl%main)) then
   ! Close info listing
   close(unit=mpl%info_unit)

   ! Close test listing
   close(unit=mpl%test_unit)
end if

end subroutine mpl_close_listing

!----------------------------------------------------------------------
! Subroutine: mpl_abort
! Purpose: clean MPI abort
!----------------------------------------------------------------------
subroutine mpl_abort(mpl,message)

implicit none

! Passed variable
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: message ! Message

! Write message
write(mpl%info,'(a)') trim(mpl%err)//'!!! Error: '//trim(message)//trim(mpl%black)

! Flush listing
call mpl%flush

! Write standard output message
write(output_unit,'(a,i4.4,a)') '!!! ABORT on task #',mpl%myproc,': '//trim(message)
call flush(output_unit)

! Abort MPI
call mpl%f_comm%abort(1)

end subroutine mpl_abort

!----------------------------------------------------------------------
! Subroutine: mpl_warning
! Purpose: print warning message
!----------------------------------------------------------------------
subroutine mpl_warning(mpl,message)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: message ! Message

! Print warning message
write(mpl%info,'(a)') trim(mpl%wng)//'!!! Warning: '//trim(message)//trim(mpl%black)
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
real(kind_real) :: prog

! Update progression array
if (present(i)) mpl%done(i) = .true.

! Print message
prog = 100.0*real(count(mpl%done),kind_real)/real(mpl%nprog,kind_real)
if (int(prog)>mpl%progint) then
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
if (info/=nf90_noerr) call mpl%abort('in '//trim(subr)//': '//trim(nf90_strerror(info)))

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
dp_loc(1) = sum(fld1*fld2)

! Allreduce
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
dp_loc(1) = sum(fld1*fld2)

! Allreduce
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
dp_loc(1) = sum(fld1*fld2)

! Allreduce
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
dp_loc(1) = sum(fld1*fld2)

! Allreduce
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%ioproc-1)

end subroutine mpl_dot_prod_4d

!----------------------------------------------------------------------
! Subroutine: mpl_split
! Purpose: split array over different MPI tasks
!----------------------------------------------------------------------
subroutine mpl_split(mpl,n,i_s,i_e,n_loc)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl       ! MPI data
integer,intent(in) :: n                 ! Total array size
integer,intent(out) :: i_s(mpl%nproc)   ! Index start
integer,intent(out) :: i_e(mpl%nproc)   ! Index end
integer,intent(out) :: n_loc(mpl%nproc) ! Local array size

! Local variable
integer :: iproc,nres,delta

! MPI splitting
nres = n
do iproc=1,mpl%nproc
   if (iproc==1) then
      i_s(iproc) = 1
   else
      i_s(iproc) = i_e(iproc-1)+1
   end if
   delta = n/mpl%nproc
   if (nres>(mpl%nproc-iproc+1)*delta) delta = delta+1
   i_e(iproc) = i_s(iproc)+delta-1
   n_loc(iproc) = delta
   nres = nres-delta
end do

end subroutine mpl_split

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_index
! Purpose: communicate global index to local index
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_index(mpl,n_loc,loc_to_glb,n_glb,glb_to_loc)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
integer,intent(in) :: n_loc              ! Local dimension
integer,intent(in) :: loc_to_glb(n_loc)  ! Local to global index
integer,intent(in) :: n_glb              ! Global dimension
integer,intent(out) :: glb_to_loc(n_glb) ! Global to local index

! Local variables
integer :: iproc,i_loc,n_loc_tmp
integer,allocatable :: loc_to_glb_tmp(:)
type(fckit_mpi_status) :: status

if (mpl%main) then
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

      ! Fill glb_to_loc
      do i_loc=1,n_loc_tmp
         glb_to_loc(loc_to_glb_tmp(i_loc)) = i_loc
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

! Check global array size
if (mpl%main) then
   if (size(glb)/=n_glb) call mpl%abort('wrong dimension for the global array in mpl_glb_to_loc_real_1d')
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
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main) then
   if (size(glb,1)/=n_glb) call mpl%abort('wrong first dimension for the global array in mpl_glb_to_loc_real_2d')
   if (size(glb,2)/=nl) call mpl%abort('wrong second dimension for the global array in mpl_glb_to_loc_real_2d')
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
type(fckit_mpi_status) :: status

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
real(kind_real),allocatable :: rbuf(:)

! Check global array size
if (mpl%main.or.bcast) then
   if (size(glb)/=n_glb) call mpl%abort('wrong dimension for the global array in mpl_loc_to_glb_real_1d')
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
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main.or.bcast) then
   if (size(glb,1)/=n_glb) call mpl%abort('wrong first dimension for the global array in mpl_loc_to_glb_real_2d')
   if (size(glb,2)/=nl) call mpl%abort('wrong second dimension for the global array in mpl_loc_to_glb_real_1d')
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
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main.or.bcast) then
   if (size(glb,1)/=n_glb) call mpl%abort('wrong first dimension for the global array in mpl_loc_to_glb_real_2d')
   if (size(glb,2)/=nl) call mpl%abort('wrong second dimension for the global array in mpl_loc_to_glb_real_1d')
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

end module type_mpl
