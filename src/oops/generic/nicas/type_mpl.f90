!----------------------------------------------------------------------
! Module: type_mpl
!> Purpose: MPI parameters derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_mpl

use iso_c_binding
use mpi
use omp_lib
use tools_kinds, only: kind_real
implicit none

type mpltype
   ! MPL parameters
   integer :: nproc      !< Number of MPI tasks
   integer :: myproc     !< MPI task index
   integer :: ioproc = 1 !< Main task index
   logical :: main       !< Main task logical
   integer :: unit       !< Listing unit
   integer :: rtype      !< MPI real type
   integer :: tag        !< MPI tag
   integer :: nthread    !< Number of OpenMP threads
end type mpltype

type(mpltype) :: mpl

interface mpl_bcast
  module procedure mpl_bcast_integer
  module procedure mpl_bcast_integer_array_1d
  module procedure mpl_bcast_integer_array_2d
  module procedure mpl_bcast_real
  module procedure mpl_bcast_real_array_1d
  module procedure mpl_bcast_real_array_2d
  module procedure mpl_bcast_logical
  module procedure mpl_bcast_logical_array_1d
  module procedure mpl_bcast_string
end interface

interface mpl_recv
  module procedure mpl_recv_integer
  module procedure mpl_recv_integer_array_1d
  module procedure mpl_recv_real_array_1d
  module procedure mpl_recv_logical_array_1d
end interface

interface mpl_send
  module procedure mpl_send_integer
  module procedure mpl_send_integer_array_1d
  module procedure mpl_send_real_array_1d
  module procedure mpl_send_logical_array_1d
end interface

private
public :: mpl
public :: mpl_start,mpl_end,mpl_abort,mpl_barrier,mpl_bcast,mpl_recv,mpl_send

contains

!----------------------------------------------------------------------
! Subroutine: mpl_check
!> Purpose: check MPI error
!----------------------------------------------------------------------
subroutine mpl_check(info)

implicit none

! Passed variables
integer,intent(in) :: info !< Error index

! Local variables
integer :: len,ierr
character(len=mpi_max_error_string) :: message

if (info/=mpi_success) then
   ! Get string
   call mpi_error_string(info,message,len,ierr)

   ! Abort MPI
   call mpl_abort(message(1:len))
endif

end subroutine mpl_check

!----------------------------------------------------------------------
! Subroutine: mpl_start
!> Purpose: start MPI
!----------------------------------------------------------------------
subroutine mpl_start() bind(c, name='mpl_start_f90')

implicit none

! Local variables
integer :: info
logical :: init

! Check whether MPI has been initialized
call mpi_initialized(init,info)
call mpl_check(info)

if (.not.init) then
   ! Initialize MPI
   call mpi_init(info)
   call mpl_check(info)
end if

! Get MPI size
call mpi_comm_size(mpi_comm_world,mpl%nproc,info)
call mpl_check(info)

! Get MPI rank
call mpi_comm_rank(mpi_comm_world,mpl%myproc,info)
call mpl_check(info)
mpl%myproc = mpl%myproc+1

! Define main task
mpl%main = (mpl%myproc==mpl%ioproc)

! Define unit and open file
mpl%unit = 5+mpl%myproc

! Define real type for bcast
if (kind_real==4) then
   mpl%rtype = mpi_real
elseif (kind_real==8) then
   mpl%rtype = mpi_double
else
   call mpl_abort('unknown real kind for mpl_bcast interface')
end if

! Initialize tag
mpl%tag = 4321

! Set max number of OpenMP threads
mpl%nthread = omp_get_max_threads()
call omp_set_num_threads(mpl%nthread)

end subroutine mpl_start

!----------------------------------------------------------------------
! Subroutine: mpl_end
!> Purpose: end MPI
!----------------------------------------------------------------------
subroutine mpl_end() bind(c, name='mpl_end_f90')

implicit none

! Local variables
integer :: info

! Finalize MPI
call mpi_finalize(info)

! Check
call mpl_check(info)

end subroutine mpl_end

!----------------------------------------------------------------------
! Subroutine: mpl_abort
!> Purpose: clean MPI abort
!----------------------------------------------------------------------
subroutine mpl_abort(message)

implicit none

! Passed variables
character(len=*),intent(in) :: message !< Message

! Flush
call flush(mpl%unit)

! Write message
write(mpl%unit,'(a)') trim(message)

! Finalize MPI
call mpl_end

! Stop
stop

end subroutine mpl_abort

!----------------------------------------------------------------------
! Subroutine: mpl_barrier
!> Purpose: MPI barrier
!----------------------------------------------------------------------
subroutine mpl_barrier

implicit none

! Local variable
integer :: info

! Wait
call mpi_barrier(mpi_comm_world,info)

! Check
call mpl_check(info)

end subroutine mpl_barrier

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_integer
!> Purpose: broadcast integer
!----------------------------------------------------------------------
subroutine mpl_bcast_integer(var,root)

implicit none

! Passed variables
integer,intent(in) :: var  !< Integer
integer,intent(in) :: root !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,1,mpi_integer,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_integer

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_integer_array_1d
!> Purpose: broadcast 1d integer array
!----------------------------------------------------------------------
subroutine mpl_bcast_integer_array_1d(var,root)

implicit none

! Passed variables
integer,dimension(:),intent(in) :: var !< Integer array, 1d
integer,intent(in) :: root             !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,size(var),mpi_integer,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_integer_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_integer_array_2d
!> Purpose: broadcast 2d integer array
!----------------------------------------------------------------------
subroutine mpl_bcast_integer_array_2d(var,root)

implicit none

! Passed variables
integer,dimension(:,:),intent(in) :: var !< Integer array, 2d
integer,intent(in) :: root               !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,size(var),mpi_integer,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_integer_array_2d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real
!> Purpose: broadcast real
!----------------------------------------------------------------------
subroutine mpl_bcast_real(var,root)

implicit none

! Passed variables
real(kind_real),intent(in) :: var !< Real
integer,intent(in) :: root        !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,1,mpl%rtype,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_real

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_1d
!> Purpose: broadcast 1d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_1d(var,root)

implicit none

! Passed variables
real(kind_real),dimension(:),intent(in) :: var !< Real array, 1d
integer,intent(in) :: root                     !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_real_array_2d
!> Purpose: broadcast 2d real array
!----------------------------------------------------------------------
subroutine mpl_bcast_real_array_2d(var,root)

implicit none

! Passed variables
real(kind_real),dimension(:,:),intent(in) :: var !< Real array, 2d
integer,intent(in) :: root                       !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,size(var),mpl%rtype,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_real_array_2d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_logical
!> Purpose: broadcast logical
!----------------------------------------------------------------------
subroutine mpl_bcast_logical(var,root)

implicit none

! Passed variables
logical,intent(in) :: var  !< Logical
integer,intent(in) :: root !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,1,mpi_logical,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_logical

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_logical_array_1d
!> Purpose: broadcast 1d logical array
!----------------------------------------------------------------------
subroutine mpl_bcast_logical_array_1d(var,root)

implicit none

! Passed variables
logical,dimension(:),intent(in) :: var !< Logical array, 1d
integer,intent(in) :: root             !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,size(var),mpi_logical,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_logical_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_bcast_string
!> Purpose: broadcast string
!----------------------------------------------------------------------
subroutine mpl_bcast_string(var,root)

implicit none

! Passed variables
character(len=*),intent(in) :: var  !< String
integer,intent(in) :: root          !< Root task

! Local variable
integer :: info

! Broadcast
call mpi_bcast(var,len(var),mpi_character,root-1,mpi_comm_world,info)

! Check
call mpl_check(info)

! Wait
call mpl_barrier

end subroutine mpl_bcast_string

!----------------------------------------------------------------------
! Subroutine: mpl_recv_integer
!> Purpose: receive integer
!----------------------------------------------------------------------
subroutine mpl_recv_integer(var,src,tag)

implicit none

! Passed variables
integer,intent(out) :: var !< Integer
integer,intent(in) :: src  !< Source task
integer,intent(in) :: tag  !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,1,mpi_integer,src-1,tag,mpi_comm_world,status,info)

! Check
call mpl_check(info)

end subroutine mpl_recv_integer

!----------------------------------------------------------------------
! Subroutine: mpl_recv_integer_array_1d
!> Purpose: receive 1d integer array
!----------------------------------------------------------------------
subroutine mpl_recv_integer_array_1d(n,var,src,tag)

implicit none

! Passed variables
integer,intent(in) :: n       !< Array size
integer,intent(out) :: var(n) !< Integer array, 1d
integer,intent(in) :: src     !< Source task
integer,intent(in) :: tag     !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,n,mpi_integer,src-1,tag,mpi_comm_world,status,info)

! Check
call mpl_check(info)

end subroutine mpl_recv_integer_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_recv_real_array_1d
!> Purpose: receive 1d real array
!----------------------------------------------------------------------
subroutine mpl_recv_real_array_1d(n,var,src,tag)

implicit none

! Passed variables
integer,intent(in) :: n               !< Array size
real(kind_real),intent(out) :: var(n) !< Real array, 1d
integer,intent(in) :: src             !< Source task
integer,intent(in) :: tag             !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,n,mpl%rtype,src-1,tag,mpi_comm_world,status,info)

! Check
call mpl_check(info)

end subroutine mpl_recv_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_recv_logical_array_1d
!> Purpose: receive 1d logical array
!----------------------------------------------------------------------
subroutine mpl_recv_logical_array_1d(n,var,src,tag)

implicit none

! Passed variables
integer,intent(in) :: n       !< Array size
logical,intent(out) :: var(n) !< Logical array, 1d
integer,intent(in) :: src     !< Source task
integer,intent(in) :: tag     !< Tag

! Local variable
integer :: info
integer,dimension(mpi_status_size) :: status

! Receive
call mpi_recv(var,n,mpi_logical,src-1,tag,mpi_comm_world,status,info)

! Check
call mpl_check(info)

end subroutine mpl_recv_logical_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_send_integer
!> Purpose: send integer
!----------------------------------------------------------------------
subroutine mpl_send_integer(var,dst,tag)

implicit none

! Passed variables
integer,intent(in) :: var !< Integer
integer,intent(in) :: dst !< Destination task
integer,intent(in) :: tag !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,1,mpi_integer,dst-1,tag,mpi_comm_world,info)

! Check
call mpl_check(info)

end subroutine mpl_send_integer

!----------------------------------------------------------------------
! Subroutine: mpl_send_integer_array_1d
!> Purpose: send 1d integer array
!----------------------------------------------------------------------
subroutine mpl_send_integer_array_1d(n,var,dst,tag)

implicit none

! Passed variables
integer,intent(in) :: n      !< Array size
integer,intent(in) :: var(n) !< Integer array, 1d
integer,intent(in) :: dst    !< Destination task
integer,intent(in) :: tag    !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,n,mpi_integer,dst-1,tag,mpi_comm_world,info)

! Check
call mpl_check(info)

end subroutine mpl_send_integer_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_send_integer_array_1d
!> Purpose: send 1d real array
!----------------------------------------------------------------------
subroutine mpl_send_real_array_1d(n,var,dst,tag)

implicit none

! Passed variables
integer,intent(in) :: n              !< Array size
real(kind_real),intent(in) :: var(n) !< Real array, 1d
integer,intent(in) :: dst            !< Destination task
integer,intent(in) :: tag            !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,n,mpl%rtype,dst-1,tag,mpi_comm_world,info)

! Check
call mpl_check(info)

end subroutine mpl_send_real_array_1d

!----------------------------------------------------------------------
! Subroutine: mpl_send_logical_array_1d
!> Purpose: send 1d logical array
!----------------------------------------------------------------------
subroutine mpl_send_logical_array_1d(n,var,dst,tag)

implicit none

! Passed variables
integer,intent(in) :: n      !< Array size
logical,intent(in) :: var(n) !< Logical array, 1d
integer,intent(in) :: dst    !< Destination task
integer,intent(in) :: tag    !< Tag

! Local variable
integer :: info

! Send
call mpi_send(var,n,mpi_logical,dst-1,tag,mpi_comm_world,info)

! Check
call mpl_check(info)

end subroutine mpl_send_logical_array_1d

end module type_mpl
