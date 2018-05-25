!----------------------------------------------------------------------
! Module: type_mpl
!> Purpose: MPI parameters derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_mpl

use iso_c_binding
use mpi
!$ use omp_lib
use tools_const, only: msvali,msvalr
use tools_kinds, only: kind_real

implicit none

type mpl_type
   ! MPL parameters
   integer :: mpi_comm   !< MPI communicator
   integer :: nproc      !< Number of MPI tasks
   integer :: myproc     !< MPI task index
   integer :: ioproc = 1 !< Main task index
   logical :: main       !< Main task logical
   integer :: unit       !< Listing unit
   integer :: rtype      !< MPI real type
   integer :: tag        !< MPI tag
   integer :: nthread    !< Number of OpenMP threads
contains
   procedure :: newunit => mpl_newunit
   procedure :: check => mpl_check
   procedure :: init => mpl_init
   procedure :: init_listing => mpl_init_listing
   procedure :: abort => mpl_abort
   procedure :: mpl_bcast_integer
   procedure :: mpl_bcast_integer_array_1d
   procedure :: mpl_bcast_integer_array_2d
   procedure :: mpl_bcast_real
   procedure :: mpl_bcast_real_array_1d
   procedure :: mpl_bcast_real_array_2d
   procedure :: mpl_bcast_real_array_3d
   procedure :: mpl_bcast_real_array_4d
   procedure :: mpl_bcast_real_array_5d
   procedure :: mpl_bcast_real_array_6d
   procedure :: mpl_bcast_logical
   procedure :: mpl_bcast_logical_array_1d
   procedure :: mpl_bcast_logical_array_2d
   procedure :: mpl_bcast_logical_array_3d
   procedure :: mpl_bcast_string
   procedure :: mpl_bcast_string_array_1d
   generic :: bcast => mpl_bcast_integer,mpl_bcast_integer_array_1d,mpl_bcast_integer_array_2d,mpl_bcast_real, &
                     & mpl_bcast_real_array_1d,mpl_bcast_real_array_2d,mpl_bcast_real_array_3d,mpl_bcast_real_array_4d, &
                     & mpl_bcast_real_array_5d,mpl_bcast_real_array_6d,mpl_bcast_logical,mpl_bcast_logical_array_1d, &
                     & mpl_bcast_logical_array_2d,mpl_bcast_logical_array_3d,mpl_bcast_string,mpl_bcast_string_array_1d
   procedure :: mpl_recv_integer
   procedure :: mpl_recv_integer_array_1d
   procedure :: mpl_recv_real_array_1d
   procedure :: mpl_recv_logical_array_1d
   generic :: recv => mpl_recv_integer,mpl_recv_integer_array_1d,mpl_recv_real_array_1d,mpl_recv_logical_array_1d
   procedure :: mpl_send_integer
   procedure :: mpl_send_integer_array_1d
   procedure :: mpl_send_real_array_1d
   procedure :: mpl_send_logical_array_1d
   generic :: send => mpl_send_integer,mpl_send_integer_array_1d,mpl_send_real_array_1d,mpl_send_logical_array_1d
   procedure :: mpl_gatherv_real
   generic :: gatherv => mpl_gatherv_real
   procedure :: mpl_scatterv_real
   generic :: scatterv => mpl_scatterv_real
   procedure :: mpl_allgather_integer
   procedure :: mpl_allgather_real
   procedure :: mpl_allgather_logical
   generic :: allgather => mpl_allgather_integer,mpl_allgather_real,mpl_allgather_logical
   procedure :: mpl_allgatherv_real
   generic :: allgatherv => mpl_allgatherv_real
   procedure :: mpl_alltoallv_real
   generic :: alltoallv => mpl_alltoallv_real
   procedure :: mpl_allreduce_sum_integer
   procedure :: mpl_allreduce_sum_real
   procedure :: mpl_allreduce_sum_real_array_1d
   generic :: allreduce_sum => mpl_allreduce_sum_integer,mpl_allreduce_sum_real,mpl_allreduce_sum_real_array_1d
   procedure :: mpl_allreduce_min_real
   procedure :: mpl_allreduce_min_real_array_1d
   generic :: allreduce_min => mpl_allreduce_min_real,mpl_allreduce_min_real_array_1d
   procedure :: mpl_allreduce_max_real
   procedure :: mpl_allreduce_max_real_array_1d
   generic :: allreduce_max => mpl_allreduce_max_real,mpl_allreduce_max_real_array_1d
   procedure :: mpl_dot_prod_1d
   procedure :: mpl_dot_prod_2d
   procedure :: mpl_dot_prod_3d
   procedure :: mpl_dot_prod_4d
   generic :: dot_prod => mpl_dot_prod_1d,mpl_dot_prod_2d,mpl_dot_prod_3d,mpl_dot_prod_4d
   procedure :: split => mpl_split
end type mpl_type

type(mpl_type) :: mpl

integer,parameter :: lunit_min=10   !< Minimum unit number
integer,parameter :: lunit_max=1000 !< Maximum unit number

private
public :: mpl

contains

!----------------------------------------------------------------------
! Subroutine: mpl_newunit
!> Purpose: find a free unit
!----------------------------------------------------------------------
subroutine mpl_newunit(mpl,lunit)

implicit none

! Passed variables
class(mpl_type) :: mpl       !< MPL object
integer,intent(out) :: lunit !< New unit

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
! Subroutine: mpl_check
!> Purpose: check MPI error
!----------------------------------------------------------------------
subroutine mpl_check(mpl,info)

implicit none

! Passed variables
class(mpl_type) :: mpl     !< MPL object
integer,intent(in) :: info !< Error index

! Local variables
integer :: len,info_loc
character(len=mpi_max_error_string) :: message

if (info/=mpi_success) then
   ! Get string
   call mpi_error_string(info,message,len,info_loc)

   ! Abort MPI
   call mpl%abort(message(1:len))
end if

end subroutine mpl_check

!----------------------------------------------------------------------
! Subroutine: mpl_init
!> Purpose: start MPI
!----------------------------------------------------------------------
subroutine mpl_init(mpl,mpi_comm)

implicit none

! Passed variables
class(mpl_type) :: mpl                 !< MPL object
integer,intent(in) :: mpi_comm         !< MPI communicator

! Local variables
integer :: info

! Copy MPI communicator
mpl%mpi_comm = mpi_comm

! Get MPI size
call mpi_comm_size(mpl%mpi_comm,mpl%nproc,info)
call mpl%check(info)

! Get MPI rank
call mpi_comm_rank(mpl%mpi_comm,mpl%myproc,info)
call mpl%check(info)
mpl%myproc = mpl%myproc+1

! Define main task
mpl%main = (mpl%myproc==mpl%ioproc)

! Define real type for MPI
if (kind_real==4) then
   mpl%rtype = mpi_real
elseif (kind_real==8) then
   mpl%rtype = mpi_double
else
   call mpl%abort('Unknown real kind for MPI')
end if

! Initialize tag
mpl%tag = 4321

! Set max number of OpenMP threads
mpl%nthread = 1
!$ mpl%nthread = omp_get_max_threads()
!$ call omp_set_num_threads(mpl%nthread)

end subroutine mpl_init

!----------------------------------------------------------------------
! Subroutine: mpl_init_listing
!> Purpose: start MPI
!----------------------------------------------------------------------
subroutine mpl_init_listing(mpl,prefix,listing)

implicit none

! Passed variables
class(mpl_type) :: mpl                 !< MPL object
character(len=*),intent(in) :: prefix  !< Output prefix
integer,intent(in),optional :: listing !< Main listing unit

! Local variables
integer :: iproc,info
character(len=4) :: myprocchar

! Define unit and open file
do iproc=1,mpl%nproc
   if (mpl%main.and.present(listing)) then
      ! Specific listing unit
      mpl%unit = listing
   else
      ! Deal with each proc sequentially
      if (iproc==mpl%myproc) then
         ! Find a free unit
         call mpl%newunit(mpl%unit)

         ! Open listing file
         write(myprocchar,'(i4.4)') mpl%myproc-1
         open(unit=mpl%unit,file=trim(prefix)//'.out.'//myprocchar,action='write',status='replace')
      end if
   end if

   ! Wait
   call mpi_barrier(mpl%mpi_comm,info)

   ! Check
   call mpl%check(info)
end do

end subroutine mpl_init_listing

!----------------------------------------------------------------------
! Subroutine: mpl_abort
!> Purpose: clean MPI abort
!----------------------------------------------------------------------
subroutine mpl_abort(mpl,message)

implicit none

! Passed variable
class(mpl_type) :: mpl                 !< MPL object
character(len=*),intent(in) :: message !< Message

! Local variables
integer :: info

! Write message
write(mpl%unit,'(a)') trim(message)
call flush(mpl%unit)

! Finalize MPI
call mpi_finalize(info)
call mpl%check(info)

! Stop
stop

end subroutine mpl_abort

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_integer
!> Purpose: broadcast integer
!----------------------------------------------------------------------
subroutine mpl_bcast_integer(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl              !< MPL object
integer,intent(in) :: var           !< Integer
integer,intent(in),optional :: root !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,1,mpi_integer,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_integer

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_integer_array_1d
!> Purpose: broadcast 1d integer array
!----------------------------------------------------------------------
subroutine mpl_bcast_integer_array_1d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                 !< MPL object
integer,dimension(:),intent(in) :: var !< Integer array, 1d
integer,intent(in),optional :: root    !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpi_integer,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_integer_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_integer_array_2d
!> Purpose: broadcast 2d integer array
!----------------------------------------------------------------------
subroutine mpl_bcast_integer_array_2d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                   !< MPL object
integer,dimension(:,:),intent(in) :: var !< Integer array, 2d
integer,intent(in),optional :: root      !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpi_integer,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_integer_array_2d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real
!> Purpose: broadcast real
!----------------------------------------------------------------------
subroutine mpl_bcast_real(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl              !< MPL object
real(kind_real),intent(in) :: var   !< Real
integer,intent(in),optional :: root !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,1,mpl%rtype,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_real

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_1d
!> Purpose: broadcast 1d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_1d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                         !< MPL object
real(kind_real),dimension(:),intent(in) :: var !< Real array, 1d
integer,intent(in),optional :: root            !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_2d
!> Purpose: broadcast 2d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_2d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                           !< MPL object
real(kind_real),dimension(:,:),intent(in) :: var !< Real array, 2d
integer,intent(in),optional :: root              !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_real_array_2d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_3d
!> Purpose: broadcast 3d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_3d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                             !< MPL object
real(kind_real),dimension(:,:,:),intent(in) :: var !< Real array, 3d
integer,intent(in),optional :: root                !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_real_array_3d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_4d
!> Purpose: broadcast 4d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_4d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                               !< MPL object
real(kind_real),dimension(:,:,:,:),intent(in) :: var !< Real array, 4d
integer,intent(in),optional :: root                  !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_real_array_4d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_5d
!> Purpose: broadcast 5d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_5d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                                 !< MPL object
real(kind_real),dimension(:,:,:,:,:),intent(in) :: var !< Real array, 5d
integer,intent(in),optional :: root                    !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_real_array_5d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_6d
!> Purpose: broadcast 6d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_6d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                                   !< MPL object
real(kind_real),dimension(:,:,:,:,:,:),intent(in) :: var !< Real array, 6d
integer,intent(in),optional :: root                      !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_real_array_6d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_logical
!> Purpose: broadcast logical
!----------------------------------------------------------------------
subroutine mpl_bcast_logical(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl              !< MPL object
logical,intent(in) :: var           !< Logical
integer,intent(in),optional :: root !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,1,mpi_logical,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_logical

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_logical_array_1d
!> Purpose: broadcast 1d logical array
!----------------------------------------------------------------------
subroutine mpl_bcast_logical_array_1d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                 !< MPL object
logical,dimension(:),intent(in) :: var !< Logical array, 1d
integer,intent(in),optional :: root    !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpi_logical,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_logical_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_logical_array_2d
!> Purpose: broadcast 2d logical array
!----------------------------------------------------------------------
subroutine mpl_bcast_logical_array_2d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                   !< MPL object
logical,dimension(:,:),intent(in) :: var !< Logical array, 1d
integer,intent(in),optional :: root      !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpi_logical,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_logical_array_2d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_logical_array_3d
!> Purpose: broadcast 3d logical array
!----------------------------------------------------------------------
subroutine mpl_bcast_logical_array_3d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                     !< MPL object
logical,dimension(:,:,:),intent(in) :: var !< Logical array, 1d
integer,intent(in),optional :: root        !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,size(var),mpi_logical,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_logical_array_3d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_string
!> Purpose: broadcast string
!----------------------------------------------------------------------
subroutine mpl_bcast_string(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl              !< MPL object
character(len=*),intent(in) :: var  !< String
integer,intent(in),optional :: root !< Root task

! Local variable
integer :: info,mpi_root

! Find root
if (present(root)) then
   mpi_root = root-1
else
   mpi_root = mpl%ioproc-1
end if

! Broadcast
call mpi_bcast(var,len(var),mpi_character,mpi_root,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_bcast_string

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_string_array_1d
!> Purpose: broadcast 1d string array
!----------------------------------------------------------------------
subroutine mpl_bcast_string_array_1d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type) :: mpl                          !< MPL object
character(len=*),dimension(:),intent(in) :: var !< Logical array, 1d
integer,intent(in),optional :: root             !< Root task

! Local variable
integer :: i

! Broadcast one string at a time
do i=1,size(var)
   if (present(root)) then
      call mpl%bcast(var(i),root)
   else
      call mpl%bcast(var(i))
   end if
end do

end subroutine mpl_bcast_string_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_recv_integer
!> Purpose: receive integer
!----------------------------------------------------------------------
subroutine mpl_recv_integer(mpl,var,src,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl     !< MPL object
integer,intent(out) :: var !< Integer
integer,intent(in) :: src  !< Source task
integer,intent(in) :: tag  !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,1,mpi_integer,src-1,tag,mpl%mpi_comm,status,info)

! Check
call mpl%check(info)

end subroutine mpl_recv_integer

!----------------------------------------------------------------------
! Subroutine: mpl_recv_integer_array_1d
!> Purpose: receive 1d integer array
!----------------------------------------------------------------------
subroutine mpl_recv_integer_array_1d(mpl,n,var,src,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl        !< MPL object
integer,intent(in) :: n       !< Array size
integer,intent(out) :: var(n) !< Integer array, 1d
integer,intent(in) :: src     !< Source task
integer,intent(in) :: tag     !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,n,mpi_integer,src-1,tag,mpl%mpi_comm,status,info)

! Check
call mpl%check(info)

end subroutine mpl_recv_integer_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_recv_real_array_1d
!> Purpose: receive 1d real array
!----------------------------------------------------------------------
subroutine mpl_recv_real_array_1d(mpl,n,var,src,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl                !< MPL object
integer,intent(in) :: n               !< Array size
real(kind_real),intent(out) :: var(n) !< Real array, 1d
integer,intent(in) :: src             !< Source task
integer,intent(in) :: tag             !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,n,mpl%rtype,src-1,tag,mpl%mpi_comm,status,info)

! Check
call mpl%check(info)

end subroutine mpl_recv_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_recv_logical_array_1d
!> Purpose: receive 1d logical array
!----------------------------------------------------------------------
subroutine mpl_recv_logical_array_1d(mpl,n,var,src,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl        !< MPL object
integer,intent(in) :: n       !< Array size
logical,intent(out) :: var(n) !< Logical array, 1d
integer,intent(in) :: src     !< Source task
integer,intent(in) :: tag     !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,n,mpi_logical,src-1,tag,mpl%mpi_comm,status,info)

! Check
call mpl%check(info)

end subroutine mpl_recv_logical_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_send_integer
!> Purpose: send integer
!----------------------------------------------------------------------
subroutine mpl_send_integer(mpl,var,dst,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl    !< MPL object
integer,intent(in) :: var !< Integer
integer,intent(in) :: dst !< Destination task
integer,intent(in) :: tag !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,1,mpi_integer,dst-1,tag,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_send_integer

!----------------------------------------------------------------------
! Subroutine: mpl_send_integer_array_1d
!> Purpose: send 1d integer array
!----------------------------------------------------------------------
subroutine mpl_send_integer_array_1d(mpl,n,var,dst,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl       !< MPL object
integer,intent(in) :: n      !< Array size
integer,intent(in) :: var(n) !< Integer array, 1d
integer,intent(in) :: dst    !< Destination task
integer,intent(in) :: tag    !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,n,mpi_integer,dst-1,tag,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_send_integer_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_send_integer_array_1d
!> Purpose: send 1d real array
!----------------------------------------------------------------------
subroutine mpl_send_real_array_1d(mpl,n,var,dst,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl               !< MPL object
integer,intent(in) :: n              !< Array size
real(kind_real),intent(in) :: var(n) !< Real array, 1d
integer,intent(in) :: dst            !< Destination task
integer,intent(in) :: tag            !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,n,mpl%rtype,dst-1,tag,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_send_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_send_logical_array_1d
!> Purpose: send 1d logical array
!----------------------------------------------------------------------
subroutine mpl_send_logical_array_1d(mpl,n,var,dst,tag)

implicit none

! Passed variables
class(mpl_type) :: mpl       !< MPL object
integer,intent(in) :: n      !< Array size
logical,intent(in) :: var(n) !< Logical array, 1d
integer,intent(in) :: dst    !< Destination task
integer,intent(in) :: tag    !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,n,mpi_logical,dst-1,tag,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_send_logical_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_gatherv_real
!> Purpose: gatherv for a real array
!----------------------------------------------------------------------
subroutine mpl_gatherv_real(mpl,ns,sbuf,rcounts,nr,rbuf)

implicit none

! Passed variables
class(mpl_type) :: mpl                   !< MPL object
integer,intent(in) :: ns                 !< Sent buffer size
real(kind_real),intent(in) :: sbuf(ns)   !< Sent buffer
integer,intent(in) :: rcounts(mpl%nproc) !< Received counts
integer,intent(in) :: nr                 !< Received buffer size
real(kind_real),intent(out) :: rbuf(nr)  !< Received buffer

! Local variable
integer :: info,displ(mpl%nproc),iproc

! Check
if (sum(rcounts)/=nr) call mpl%abort('inconsistency in mpl_gatherv_real')

! Define displacement
displ(1) = 0
do iproc=2,mpl%nproc
   displ(iproc) = displ(iproc-1)+rcounts(iproc-1)
end do

! Gatherv
call mpi_gatherv(sbuf,ns,mpl%rtype,rbuf,rcounts,displ,mpl%rtype,mpl%ioproc-1,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_gatherv_real

!----------------------------------------------------------------------
! Subroutine: mpl_scatterv_real
!> Purpose: scatterv for a real array
!----------------------------------------------------------------------
subroutine mpl_scatterv_real(mpl,scounts,ns,sbuf,nr,rbuf)

implicit none

! Passed variables
class(mpl_type) :: mpl                   !< MPL object
integer,intent(in) :: scounts(mpl%nproc) !< Sent counts
integer,intent(in) :: ns                 !< Sent buffer size
real(kind_real),intent(in) :: sbuf(ns)   !< Sent buffer
integer,intent(in) :: nr                 !< Received buffer size
real(kind_real),intent(out) :: rbuf(nr)  !< Received buffer

! Local variable
integer :: info,displ(mpl%nproc),iproc

! Check
if (sum(scounts)/=ns) call mpl%abort('inconsistency in mpl_scatterv_real')

! Define displacement
displ(1) = 0
do iproc=2,mpl%nproc
   displ(iproc) = displ(iproc-1)+scounts(iproc-1)
end do

! Scatterv
call mpi_scatterv(sbuf,scounts,displ,mpl%rtype,rbuf,nr,mpl%rtype,mpl%ioproc-1,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_scatterv_real

!----------------------------------------------------------------------
! Subroutine: mpl_allgather_integer
!> Purpose: allgather for a integer array
!----------------------------------------------------------------------
subroutine mpl_allgather_integer(mpl,ns,sbuf,rbuf)

implicit none

! Passed variables
class(mpl_type) :: mpl                     !< MPL object
integer,intent(in) :: ns                   !< Sent buffer size
integer,intent(in) :: sbuf(ns)             !< Sent buffer
integer,intent(out) :: rbuf(mpl%myproc*ns) !< Received buffer

! Local variable
integer :: info

! Allgather
call mpi_allgather(sbuf,ns,mpi_integer,rbuf,ns,mpi_integer,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_allgather_integer

!----------------------------------------------------------------------
! Subroutine: mpl_allgather_real
!> Purpose: allgather for a real array
!----------------------------------------------------------------------
subroutine mpl_allgather_real(mpl,ns,sbuf,rbuf)

implicit none

! Passed variables
class(mpl_type) :: mpl                             !< MPL object
integer,intent(in) :: ns                           !< Sent buffer size
real(kind_real),intent(in) :: sbuf(ns)             !< Sent buffer
real(kind_real),intent(out) :: rbuf(mpl%myproc*ns) !< Received buffer

! Local variable
integer :: info

! Allgather
call mpi_allgather(sbuf,ns,mpl%rtype,rbuf,ns,mpl%rtype,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_allgather_real

!----------------------------------------------------------------------
! Subroutine: mpl_allgather_logical
!> Purpose: allgather for a logical array
!----------------------------------------------------------------------
subroutine mpl_allgather_logical(mpl,ns,sbuf,rbuf)

implicit none

! Passed variables
class(mpl_type) :: mpl                     !< MPL object
integer,intent(in) :: ns                   !< Sent buffer size
logical,intent(in) :: sbuf(ns)             !< Sent buffer
logical,intent(out) :: rbuf(mpl%myproc*ns) !< Received buffer

! Local variable
integer :: info

! Allgather
call mpi_allgather(sbuf,ns,mpi_logical,rbuf,ns,mpi_logical,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_allgather_logical

!----------------------------------------------------------------------
! Subroutine: mpl_allgatherv_real
!> Purpose: allgatherv for a real array
!----------------------------------------------------------------------
subroutine mpl_allgatherv_real(mpl,ns,sbuf,rcounts,nr,rbuf)

implicit none

! Passed variables
class(mpl_type) :: mpl                   !< MPL object
integer,intent(in) :: ns                 !< Sent buffer size
real(kind_real),intent(in) :: sbuf(ns)   !< Sent buffer
integer,intent(in) :: rcounts(mpl%nproc) !< Received counts
integer,intent(in) :: nr                 !< Received buffer size
real(kind_real),intent(out) :: rbuf(nr)  !< Received buffer

! Local variable
integer :: info,displ(mpl%nproc),iproc

! Check
if (sum(rcounts)/=nr) call mpl%abort('inconsistency in mpl_allgatherv_real')

! Define displacement
displ(1) = 0
do iproc=2,mpl%nproc
   displ(iproc) = displ(iproc-1)+rcounts(iproc-1)
end do

! Allgatherv
call mpi_allgatherv(sbuf,ns,mpl%rtype,rbuf,rcounts,displ,mpl%rtype,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_allgatherv_real


!----------------------------------------------------------------------
! Subroutine: mpl_alltoallv_real
!> Purpose: alltoallv for a real array
!----------------------------------------------------------------------
subroutine mpl_alltoallv_real(mpl,ns,sbuf,scounts,sdispl,nr,rbuf,rcounts,rdispl)

implicit none

! Passed variables
class(mpl_type) :: mpl                   !< MPL object
integer,intent(in) :: ns                 !< Sent buffer size
real(kind_real),intent(in) :: sbuf(ns)   !< Sent buffer
integer,intent(in) :: scounts(mpl%nproc) !< Sending counts
integer,intent(in) :: sdispl(mpl%nproc)  !< Sending displacement
integer,intent(in) :: nr                 !< Received buffer size
real(kind_real),intent(out) :: rbuf(nr)  !< Received buffer
integer,intent(in) :: rcounts(mpl%nproc) !< Receiving counts
integer,intent(in) :: rdispl(mpl%nproc)  !< Receiving displacement

! Local variable
integer :: info

! Alltoallv
call mpi_alltoallv(sbuf,scounts,sdispl,mpl%rtype,rbuf,rcounts,rdispl,mpl%rtype,mpl%mpi_comm,info)

! Check
call mpl%check(info)

end subroutine mpl_alltoallv_real

!----------------------------------------------------------------------
! Subroutine: mpl_allreduce_sum_integer
!> Purpose: allreduce sum for an integer
!----------------------------------------------------------------------
subroutine mpl_allreduce_sum_integer(mpl,var_in,var_out)

implicit none

! Passed variables
class(mpl_type) :: mpl         !< MPL object
integer,intent(in) :: var_in   !< Input integer
integer,intent(out) :: var_out !< Output integer

! Local variable
integer :: info
integer :: sbuf(1),rbuf(1)

! Check for missing values
if (abs(var_in-msvali)>0) then
   sbuf(1) = var_in
else
   sbuf(1) = 0
end if

! Allreduce
call mpi_allreduce(sbuf,rbuf,1,mpi_integer,mpi_sum,mpl%mpi_comm,info)
var_out = rbuf(1)

! Check
call mpl%check(info)

end subroutine mpl_allreduce_sum_integer

!----------------------------------------------------------------------
! Subroutine: mpl_allreduce_sum_real
!> Purpose: allreduce sum for a real number
!----------------------------------------------------------------------
subroutine mpl_allreduce_sum_real(mpl,var_in,var_out)

implicit none

! Passed variables
class(mpl_type) :: mpl                 !< MPL object
real(kind_real),intent(in) :: var_in   !< Input real
real(kind_real),intent(out) :: var_out !< Output real

! Local variable
integer :: info
real(kind_real) :: sbuf(1),rbuf(1)

! Check for missing values
if (abs(var_in-msvalr)>0.0) then
   sbuf(1) = var_in
else
   sbuf(1) = 0.0
end if

! Allreduce
call mpi_allreduce(sbuf,rbuf,1,mpl%rtype,mpi_sum,mpl%mpi_comm,info)
var_out = rbuf(1)

! Check
call mpl%check(info)

end subroutine mpl_allreduce_sum_real

!----------------------------------------------------------------------
! Subroutine: mpl_allreduce_sum_real_array_1d
!> Purpose: allreduce sum for a real array, 1d
!----------------------------------------------------------------------
subroutine mpl_allreduce_sum_real_array_1d(mpl,var_in,var_out)

implicit none

! Passed variables
class(mpl_type) :: mpl                    !< MPL object
real(kind_real),intent(in) :: var_in(:)   !< Input real
real(kind_real),intent(out) :: var_out(:) !< Output real

! Local variable
integer :: i,info
real(kind_real) :: sbuf(size(var_in)),rbuf(size(var_in))

! Check for missing values
do i=1,size(var_in)
   if (abs(var_in(i)-msvalr)>0.0) then
      sbuf(i) = var_in(i)
   else
      sbuf(i) = 0.0
   end if
end do

! Allreduce
call mpi_allreduce(sbuf,rbuf,size(var_in),mpl%rtype,mpi_sum,mpl%mpi_comm,info)
var_out = rbuf

! Check
call mpl%check(info)

end subroutine mpl_allreduce_sum_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_allreduce_min_real
!> Purpose: allreduce min for a real number
!----------------------------------------------------------------------
subroutine mpl_allreduce_min_real(mpl,var_in,var_out)

implicit none

! Passed variables
class(mpl_type) :: mpl                 !< MPL object
real(kind_real),intent(in) :: var_in   !< Input real
real(kind_real),intent(out) :: var_out !< Output real

! Local variable
integer :: info
real(kind_real) :: sbuf(1),rbuf(1)

! Check for missing values
if (abs(var_in-msvalr)>0.0) then
   sbuf(1) = var_in
else
   sbuf(1) = huge(1.0)
end if

! Allreduce
call mpi_allreduce(sbuf,rbuf,1,mpl%rtype,mpi_min,mpl%mpi_comm,info)
var_out = rbuf(1)

! Check
call mpl%check(info)

end subroutine mpl_allreduce_min_real

!----------------------------------------------------------------------
! Subroutine: mpl_allreduce_min_real_array_1d
!> Purpose: allreduce min for a real array, 1d
!----------------------------------------------------------------------
subroutine mpl_allreduce_min_real_array_1d(mpl,var_in,var_out)

implicit none

! Passed variables
class(mpl_type) :: mpl                    !< MPL object
real(kind_real),intent(in) :: var_in(:)   !< Input real
real(kind_real),intent(out) :: var_out(:) !< Output real

! Local variable
integer :: i,info
real(kind_real) :: sbuf(size(var_in)),rbuf(size(var_in))

! Check for missing values
do i=1,size(var_in)
   if (abs(var_in(i)-msvalr)>0.0) then
      sbuf(i) = var_in(i)
   else
      sbuf(i) = 0.0
   end if
end do

! Allreduce
call mpi_allreduce(sbuf,rbuf,size(var_in),mpl%rtype,mpi_min,mpl%mpi_comm,info)
var_out = rbuf

! Check
call mpl%check(info)

end subroutine mpl_allreduce_min_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_allreduce_max_real
!> Purpose: allreduce max for a real number
!----------------------------------------------------------------------
subroutine mpl_allreduce_max_real(mpl,var_in,var_out)

implicit none

! Passed variables
class(mpl_type) :: mpl                 !< MPL object
real(kind_real),intent(in) :: var_in   !< Input real
real(kind_real),intent(out) :: var_out !< Output real

! Local variable
integer :: info
real(kind_real) :: sbuf(1),rbuf(1)

! Check for missing values
if (abs(var_in-msvalr)>0.0) then
   sbuf(1) = var_in
else
   sbuf(1) = -huge(1.0)
end if

! Allreduce
call mpi_allreduce(sbuf,rbuf,1,mpl%rtype,mpi_max,mpl%mpi_comm,info)
var_out = rbuf(1)

! Check
call mpl%check(info)

end subroutine mpl_allreduce_max_real

!----------------------------------------------------------------------
! Subroutine: mpl_allreduce_max_real_array_1d
!> Purpose: allreduce max for a real array, 1d
!----------------------------------------------------------------------
subroutine mpl_allreduce_max_real_array_1d(mpl,var_in,var_out)

implicit none

! Passed variables
class(mpl_type) :: mpl                    !< MPL object
real(kind_real),intent(in) :: var_in(:)   !< Input real
real(kind_real),intent(out) :: var_out(:) !< Output real

! Local variable
integer :: i,info
real(kind_real) :: sbuf(size(var_in)),rbuf(size(var_in))

! Check for missing values
do i=1,size(var_in)
   if (abs(var_in(i)-msvalr)>0.0) then
      sbuf(i) = var_in(i)
   else
      sbuf(i) = 0.0
   end if
end do

! Allreduce
call mpi_allreduce(sbuf,rbuf,size(var_in),mpl%rtype,mpi_max,mpl%mpi_comm,info)
var_out = rbuf

! Check
call mpl%check(info)

end subroutine mpl_allreduce_max_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_1d
!> Purpose: global dot product over local fields, 1d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_1d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type) :: mpl                !< MPL object
real(kind_real),intent(in) :: fld1(:) !< Field 1
real(kind_real),intent(in) :: fld2(:) !< Field 2
real(kind_real),intent(out) :: dp     !< Global dot product

! Local variable
integer :: info
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = sum(fld1*fld2)

! Allreduce
call mpi_allreduce(dp_loc,dp_out,1,mpl%rtype,mpi_sum,mpl%mpi_comm,info)
dp = dp_out(1)

! Check
call mpl%check(info)

! Broadcast
call mpl%bcast(dp)

! Check
call mpl%check(info)

end subroutine mpl_dot_prod_1d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_2d
!> Purpose: global dot product over local fields, 2d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_2d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type) :: mpl                  !< MPL object
real(kind_real),intent(in) :: fld1(:,:) !< Field 1
real(kind_real),intent(in) :: fld2(:,:) !< Field 2
real(kind_real),intent(out) :: dp       !< Global dot product

! Local variable
integer :: info
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = sum(fld1*fld2)

! Allreduce
call mpi_allreduce(dp_loc,dp_out,1,mpl%rtype,mpi_sum,mpl%mpi_comm,info)
dp = dp_out(1)

! Check
call mpl%check(info)

! Broadcast
call mpl%bcast(dp)

! Check
call mpl%check(info)

end subroutine mpl_dot_prod_2d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_3d
!> Purpose: global dot product over local fields, 3d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_3d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type) :: mpl                    !< MPL object
real(kind_real),intent(in) :: fld1(:,:,:) !< Field 1
real(kind_real),intent(in) :: fld2(:,:,:) !< Field 2
real(kind_real),intent(out) :: dp         !< Global dot product

! Local variable
integer :: info
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = sum(fld1*fld2)

! Allreduce
call mpi_allreduce(dp_loc,dp_out,1,mpl%rtype,mpi_sum,mpl%mpi_comm,info)
dp = dp_out(1)

! Check
call mpl%check(info)

! Broadcast
call mpl%bcast(dp)

! Check
call mpl%check(info)

end subroutine mpl_dot_prod_3d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_4d
!> Purpose: global dot product over local fields, 4d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_4d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type) :: mpl                      !< MPL object
real(kind_real),intent(in) :: fld1(:,:,:,:) !< Field 1
real(kind_real),intent(in) :: fld2(:,:,:,:) !< Field 2
real(kind_real),intent(out) :: dp           !< Global dot product

! Local variable
integer :: info
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = sum(fld1*fld2)

! Allreduce
call mpi_allreduce(dp_loc,dp_out,1,mpl%rtype,mpi_sum,mpl%mpi_comm,info)
dp = dp_out(1)

! Check
call mpl%check(info)

! Broadcast
call mpl%bcast(dp)

! Check
call mpl%check(info)

end subroutine mpl_dot_prod_4d

!----------------------------------------------------------------------
! Subroutine: mpl_split
!> Purpose: split array over different MPI tasks
!----------------------------------------------------------------------
subroutine mpl_split(mpl,n,i_s,i_e,n_loc)

implicit none

! Passed variables
class(mpl_type) :: mpl                  !< MPL object
integer,intent(in) :: n                 !< Total array size
integer,intent(out) :: i_s(mpl%nproc)   !< Index start
integer,intent(out) :: i_e(mpl%nproc)   !< Index end
integer,intent(out) :: n_loc(mpl%nproc) !< Local array size

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

end module type_mpl
