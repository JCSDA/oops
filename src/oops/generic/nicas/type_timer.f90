!----------------------------------------------------------------------
! Module: type_timer
!> Purpose: timer data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_timer

use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_mpl, only: mpl

implicit none

! Timer data derived type
type timertype
   real(kind_real) :: cpu_time_start !< CPU time start
   real(kind_real) :: cpu_time_end   !< CPU time end
   integer :: count_rate             !< Count rate
   integer :: count_max              !< Count maximum
   integer :: system_clock_start     !< System clock start
   integer :: system_clock_end       !< System clock end
   real(kind_real) :: cpu            !< CPU time
   real(kind_real) :: elapsed        !< Elapsed time
   character(len=9) :: VmPeak        !< Peak virtual memory size
   character(len=9) :: VmHWM         !< Peak resident set size
   character(len=9) :: rchar         !< Read data volume
   character(len=9) :: wchar         !< Written data volume
end type timertype

! I/O unit
integer,parameter :: iunit = 100 !< I/O unit

private
public :: timertype,timer_start,timer_end,timer_display

contains

!----------------------------------------------------------------------
! Subroutine: timer_start
!> Purpose: initialize timer
!----------------------------------------------------------------------
subroutine timer_start(timer)

implicit none

! Passed variables
type(timertype),intent(inout) :: timer !< Timer data

! Execution times  initialization
call system_clock(count_rate=timer%count_rate,count_max=timer%count_max)
call system_clock(count=timer%system_clock_start)
call cpu_time(timer%cpu_time_start)

end subroutine timer_start

!----------------------------------------------------------------------
! Subroutine: timer_end
!> Purpose: finalize timer
!----------------------------------------------------------------------
subroutine timer_end(timer)

implicit none

! Passed variables
type(timertype),intent(inout) :: timer !< Timer data

! Execution times calculation
call system_clock(count=timer%system_clock_end)
call cpu_time(timer%cpu_time_end)
timer%cpu = timer%cpu_time_end-timer%cpu_time_start
if (timer%system_clock_end<timer%system_clock_start) then
  timer%elapsed = float(timer%system_clock_end-timer%system_clock_start+timer%count_max)/float(timer%count_rate)
else
  timer%elapsed = float(timer%system_clock_end-timer%system_clock_start)/float(timer%count_rate)
end if

end subroutine timer_end

!----------------------------------------------------------------------
! Subroutine: timer_display
!> Purpose: display timer
!----------------------------------------------------------------------
subroutine timer_display(timer)

implicit none

! Passed variables
type(timertype),intent(inout) :: timer !< Timer data

! Local variables
integer :: ierr,get_pid
real(kind_real) :: VmPeak,VmHWM,rchar,wchar
logical :: isfile
character(len=8) :: pidchar
character(len=1024) :: filename,line

! Execution times calculation
call timer_end(timer)

! Maximum memory usage
timer%VmPeak = 'unknown'
timer%VmHWM = 'unknown'
call fgetpid(get_pid)
write(pidchar,'(i8)') get_pid
filename = '/proc/'//trim(adjustl(pidchar))//'/status'
inquire(file=filename,exist=isfile)
if (isfile) then
   open(unit=iunit,file=filename,action='read')
   do
      read(iunit,'(a)',iostat=ierr) line
      if (ierr>0) then
         call msgerror('cannot read /proc/pid/status')
      elseif (ierr<0) then
         exit
      else
         if (line(1:7)=='VmPeak:') then
            read(line(8:16),*) VmPeak
            timer%VmPeak = ''
            if (VmPeak<1e3) then
               write(timer%VmPeak(2:6),'(f5.1)') VmPeak
               timer%VmPeak(7:9) = ' kb'
            elseif (VmPeak<1e6) then
               write(timer%VmPeak(2:6),'(f5.1)') 1.0e-3*VmPeak
               timer%VmPeak(7:9) = ' Mb'
            elseif (VmPeak<1e9) then
               write(timer%VmPeak(2:6),'(f5.1)') 1.0e-6*VmPeak
               timer%VmPeak(7:9) = ' Gb'
            elseif (VmPeak<1e12) then
               write(timer%VmPeak(2:6),'(f5.1)') 1.0e-9*VmPeak
               timer%VmPeak(7:9) = ' Tb'
            end if
         end if
         if (line(1:6)=='VmHWM:') then
            read(line(8:16),*) VmHWM
            timer%VmHWM = ''
            if (VmHWM<1e3) then
               write(timer%VmHWM(2:6),'(f5.1)') VmHWM
               timer%VmHWM(7:9) = ' kb'
            elseif (VmHWM<1e6) then
               write(timer%VmHWM(2:6),'(f5.1)') 1.0e-3*VmHWM
               timer%VmHWM(7:9) = ' Mb'
            elseif (VmHWM<1e9) then
               write(timer%VmHWM(2:6),'(f5.1)') 1.0e-6*VmHWM
               timer%VmHWM(7:9) = ' Gb'
            elseif (VmHWM<1e12) then
               write(timer%VmHWM(2:6),'(f5.1)') 1.0e-9*VmHWM
               timer%VmHWM(7:9) = ' Tb'
            end if
         end if
      end if
   end do
   close(iunit)
end if

! I/O info
timer%rchar = 'unknown'
timer%wchar = 'unknown'
call fgetpid(get_pid)
write(pidchar,'(i8)') get_pid
filename = '/proc/'//trim(adjustl(pidchar))//'/io'
if (isfile) then
   open(unit=iunit,file=filename,action='read')
   do
      read(iunit,'(a)',iostat=ierr) line
      if (ierr>0) then
         call msgerror('cannot read /proc/pid/io')
      elseif (ierr<0) then
         exit
      else
         if (line(1:6)=='rchar:') then
            read(line(7:),*) rchar
            timer%rchar = ''
            if (rchar<1e3) then
               write(timer%rchar(2:6),'(f5.1)') rchar
               timer%rchar(7:9) = '  b'
            elseif (rchar<1e6) then
               write(timer%rchar(2:6),'(f5.1)') 1.0e-3*rchar
               timer%rchar(7:9) = ' kb'
            elseif (rchar<1e9) then
               write(timer%rchar(2:6),'(f5.1)') 1.0e-6*rchar
               timer%rchar(7:9) = ' Mb'
            elseif (rchar<1e12) then
               write(timer%rchar(2:6),'(f5.1)') 1.0e-9*rchar
               timer%rchar(7:9) = ' Gb'
            end if
         end if
         if (line(1:6)=='wchar:') then
            read(line(7:),*) wchar
            timer%wchar = ''
            if (wchar<1e3) then
               write(timer%wchar(2:6),'(f5.1)') wchar
               timer%wchar(7:9) = '  b'
            elseif (wchar<1e6) then
               write(timer%wchar(2:6),'(f5.1)') 1.0e-3*wchar
               timer%wchar(7:9) = ' kb'
            elseif (wchar<1e9) then
               write(timer%wchar(2:6),'(f5.1)') 1.0e-6*wchar
               timer%wchar(7:9) = ' Mb'
            elseif (wchar<1e12) then
               write(timer%wchar(2:6),'(f5.1)') 1.0e-9*wchar
               timer%wchar(7:9) = ' Gb'
            end if
         end if
      end if
   end do
   close(iunit)
endif

! Display
write(mpl%unit,'(a,f8.3,a)') '--- CPU time:                   ',timer%cpu,' s'
write(mpl%unit,'(a,f8.3,a)') '--- Elapsed time:               ',timer%elapsed,' s'
if (timer%elapsed>0.0) write(mpl%unit,'(a,f8.3)') '--- CPU/elapsed ratio:          ',timer%cpu/timer%elapsed
write(mpl%unit,'(a)') '--- Peak virtual mem. size:       '//timer%VmPeak
write(mpl%unit,'(a)') '--- Peak resident set size:       '//timer%VmHWM
write(mpl%unit,'(a)') '--- Read data volume:             '//timer%rchar
write(mpl%unit,'(a)') '--- Written data volume:          '//timer%wchar

end subroutine timer_display

end module type_timer
