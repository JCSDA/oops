!----------------------------------------------------------------------
! Module: type_obsop
! Purpose: observation operator data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_obsop

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max,fckit_mpi_status
use netcdf
use tools_const, only: pi,deg2rad,rad2deg,req,reqkm
use tools_func, only: sphere_dist
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: rth
use type_com, only: com_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Observation operator data derived type
type obsop_type
   ! Observations
   integer :: nobs                          ! Number of observations
   real(kind_real),allocatable :: lonobs(:) ! Observations longitudes
   real(kind_real),allocatable :: latobs(:) ! Observations latitudes
   integer,allocatable :: obsa_to_obs(:)    ! Local to global observation

   ! Required data to apply an observation operator

   ! Number of points
   integer :: nc0b                          ! Halo B size

   ! Number of observations
   integer :: nobsa                         ! Local number of observations

   ! Interpolation data
   type(linop_type) :: h                    ! Interpolation data

   ! Communication data
   type(com_type) :: com                    ! Communication data
contains
   procedure :: partial_dealloc => obsop_partial_dealloc
   procedure :: dealloc => obsop_dealloc
   procedure :: read => obsop_read
   procedure :: write => obsop_write
   procedure :: from => obsop_from
   procedure :: run_obsop => obsop_run_obsop
   procedure :: run_obsop_tests => obsop_run_obsop_tests
   procedure :: apply => obsop_apply
   procedure :: apply_ad => obsop_apply_ad
   procedure :: test_adjoint => obsop_test_adjoint
   procedure :: test_accuracy => obsop_test_accuracy
end type obsop_type

private
public :: obsop_type

contains

!----------------------------------------------------------------------
! Subroutine: obsop_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine obsop_partial_dealloc(obsop)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data

! Release memory
if (allocated(obsop%lonobs)) deallocate(obsop%lonobs)
if (allocated(obsop%latobs)) deallocate(obsop%latobs)

end subroutine obsop_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: obsop_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine obsop_dealloc(obsop)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data

! Release memory
call obsop%partial_dealloc
if (allocated(obsop%obsa_to_obs)) deallocate(obsop%obsa_to_obs)
call obsop%h%dealloc
call obsop%com%dealloc

end subroutine obsop_dealloc

!----------------------------------------------------------------------
! Subroutine: obsop_read
! Purpose: read observations locations
!----------------------------------------------------------------------
subroutine obsop_read(obsop,mpl,nam)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: ncid
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'obsop_read'

! Create file
write(filename,'(a,a,i4.4,a,i4.4)') trim(nam%prefix),'_obs_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

! Get attributes
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'nc0b',obsop%nc0b))
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'nobsa',obsop%nobsa))

! Read interpolation
obsop%h%prefix = 'o'
call obsop%h%read(mpl,ncid)

! Read communication
call obsop%com%read(mpl,ncid,'com')

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

end subroutine obsop_read

!----------------------------------------------------------------------
! Subroutine: obsop_write
! Purpose: write observations locations
!----------------------------------------------------------------------
subroutine obsop_write(obsop,mpl,nam)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: ncid
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'obsop_write'

! Create file
write(filename,'(a,a,i4.4,a,i4.4)') trim(nam%prefix),'_obs_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%write(mpl,ncid)

! Write attributes
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nc0b',obsop%nc0b))
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nobsa',obsop%nobsa))

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write interpolation
call obsop%h%write(mpl,ncid)

! Write communication
call obsop%com%write(mpl,ncid)

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

end subroutine obsop_write

!----------------------------------------------------------------------
! Subroutine: obsop_from
! Purpose: copy observation operator data
!----------------------------------------------------------------------
subroutine obsop_from(obsop,nobsa,lonobs,latobs)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop    ! Observation operator data
integer,intent(in) :: nobsa                 ! Number of observations
real(kind_real),intent(in) :: lonobs(nobsa) ! Observations longitudes (in degrees)
real(kind_real),intent(in) :: latobs(nobsa) ! Observations latitudes (in degrees)

! Get size
obsop%nobsa = nobsa

! Release memory
call obsop%dealloc

! Allocation
allocate(obsop%lonobs(obsop%nobsa))
allocate(obsop%latobs(obsop%nobsa))

if (obsop%nobsa>0) then
   ! Copy
   obsop%lonobs = lonobs*deg2rad
   obsop%latobs = latobs*deg2rad
end if

end subroutine obsop_from

!----------------------------------------------------------------------
! Subroutine: obsop_run_obsop
! Purpose: observation operator driver
!----------------------------------------------------------------------
subroutine obsop_run_obsop(obsop,mpl,rng,nam,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
integer :: offset,iobs,jobs,iobsa,iproc,i_s,ic0,ic0b,i,ic0a,delta,nres,ind(1),lunit,nobs_eff,nobsa_eff
integer :: nn_index(1),imin(1),imax(1),nmoves,imoves
integer :: proc_to_nobsa(mpl%nproc),proc_to_nobsa_eff(mpl%nproc),nobs_to_move(mpl%nproc),nobs_to_move_tmp(mpl%nproc)
integer :: c0_to_c0b(geom%nc0),c0a_to_c0b(geom%nc0a)
integer,allocatable :: nop(:),iop(:),srcproc(:,:),srcic0(:,:),order(:),obs_moved(:,:)
integer,allocatable :: obs_to_proc(:),obs_to_obsa(:),c0b_to_c0(:)
real(kind_real) :: nn_dist(1),N_max,C_max
real(kind_real),allocatable :: lonobs(:),latobs(:),list(:)
logical :: maskobsa(obsop%nobsa),lcheck_nc0b(geom%nc0)
logical,allocatable :: maskobs(:)
character(len=1024),parameter :: subr = 'obsop_run_obsop'
type(fckit_mpi_status) :: status
type(linop_type) :: hfull

! Check whether observations are inside the mesh
do iobsa=1,obsop%nobsa
   call geom%mesh%inside(mpl,obsop%lonobs(iobsa),obsop%latobs(iobsa),maskobsa(iobsa))
   if (.not.maskobsa(iobsa)) then
      ! Check for very close points
      call geom%tree%find_nearest_neighbors(obsop%lonobs(iobsa),obsop%latobs(iobsa),1,nn_index,nn_dist)
      if (nn_dist(1)<rth*req) maskobsa(iobsa) = .true.
   end if
end do
nobsa_eff = count(maskobsa)

! Get global number of observations
call mpl%f_comm%allgather(obsop%nobsa,proc_to_nobsa)
call mpl%f_comm%allgather(nobsa_eff,proc_to_nobsa_eff)
obsop%nobs = sum(proc_to_nobsa)
nobs_eff = sum(proc_to_nobsa_eff)

! Print input
write(mpl%info,'(a7,a)') '','Number of observations / valid observations per MPI task on input:'
call mpl%flush
do iproc=1,mpl%nproc
   write(mpl%info,'(a10,a,i3,a,i8,a,i8)') '','Task ',iproc,': ',proc_to_nobsa(iproc),' / ',proc_to_nobsa_eff(iproc)
   call mpl%flush
end do
write(mpl%info,'(a10,a,i8,a,i8)') '','Total   : ',obsop%nobs,' / ',nobs_eff
call mpl%flush

! Allocation
allocate(lonobs(obsop%nobs))
allocate(latobs(obsop%nobs))
allocate(maskobs(obsop%nobs))

! Get observations coordinates on main task
if (mpl%main) then
   offset = 0
   do iproc=1,mpl%nproc
      if (proc_to_nobsa(iproc)>0) then
         if (iproc==mpl%ioproc) then
            ! Copy data
            lonobs(offset+1:offset+proc_to_nobsa(iproc)) = obsop%lonobs
            latobs(offset+1:offset+proc_to_nobsa(iproc)) = obsop%latobs
            maskobs(offset+1:offset+proc_to_nobsa(iproc)) = maskobsa
         else
            ! Receive data on ioproc
            call mpl%f_comm%receive(lonobs(offset+1:offset+proc_to_nobsa(iproc)),iproc-1,mpl%tag,status)
            call mpl%f_comm%receive(latobs(offset+1:offset+proc_to_nobsa(iproc)),iproc-1,mpl%tag+1,status)
            call mpl%f_comm%receive(maskobs(offset+1:offset+proc_to_nobsa(iproc)),iproc-1,mpl%tag+2,status)
         end if
      end if

      ! Update offset
      offset = offset+proc_to_nobsa(iproc)
   end do
else
   if (obsop%nobsa>0) then
      ! Send data to ioproc
      call mpl%f_comm%send(obsop%lonobs,mpl%ioproc-1,mpl%tag)
      call mpl%f_comm%send(obsop%latobs,mpl%ioproc-1,mpl%tag+1)
      call mpl%f_comm%send(maskobsa,mpl%ioproc-1,mpl%tag+2)
   end if
end if
call mpl%update_tag(3)

! Broadcast data
call mpl%f_comm%broadcast(lonobs,mpl%ioproc-1)
call mpl%f_comm%broadcast(latobs,mpl%ioproc-1)
call mpl%f_comm%broadcast(maskobs,mpl%ioproc-1)

! Allocation
allocate(nop(obsop%nobs))
allocate(iop(obsop%nobs))
allocate(obs_to_proc(obsop%nobs))
allocate(obs_to_obsa(obsop%nobs))

! Compute interpolation
hfull%prefix = 'o'
write(mpl%info,'(a7,a)') '','Single level:'
call mpl%flush
call hfull%interp(mpl,geom%mesh,geom%tree,geom%nc0,geom%mask_hor_c0,obsop%nobs,lonobs,latobs,maskobs,nam%obsop_interp)

! Count interpolation points
nop = 0
do i_s=1,hfull%n_s
   iobs = hfull%row(i_s)
   nop(iobs) = nop(iobs)+1
end do

! Allocation
allocate(srcproc(maxval(nop),obsop%nobs))
allocate(srcic0(maxval(nop),obsop%nobs))

! Find grid points origin
iop = 0
srcproc = mpl%msv%vali
srcic0 = mpl%msv%vali
do i_s=1,hfull%n_s
   ic0 = hfull%col(i_s)
   iproc = geom%c0_to_proc(ic0)
   iobs = hfull%row(i_s)
   iop(iobs) = iop(iobs)+1
   srcproc(iop(iobs),iobs) = iproc
   srcic0(iop(iobs),iobs) = ic0
end do

! Generate observation distribution on processors
select case (trim(nam%obsdis))
case('')
   ! Observations distributed as provided by user
   write(mpl%info,'(a7,a)') '','Observations distributed as provided by user'
   call mpl%flush

   ! Fill obs_to_proc
   iobs = 0
   do iproc=1,mpl%nproc
      do iobsa=1,proc_to_nobsa(iproc)
         iobs = iobs+1
         obs_to_proc(iobs) = iproc
      end do
   end do
case('random')
   ! Observations randomly distributed
   write(mpl%info,'(a7,a)') '','Observations randomly distributed'
   call mpl%flush

   if (mpl%main) then
      ! Allocation
      allocate(list(obsop%nobs))
      allocate(order(obsop%nobs))

      ! Generate random order
      call rng%rand_real(0.0_kind_real,1.0_kind_real,list)
      call qsort(obsop%nobs,list,order)

      ! Split observations
      iproc = 1
      do iobs=1,obsop%nobs
         jobs = order(iobs)
         obs_to_proc(jobs) = iproc
         iproc = iproc+1
         if (iproc>mpl%nproc) iproc = 1
      end do

      ! Release memory
      deallocate(list)
      deallocate(order)
   end if

   ! Broadcast
   call mpl%f_comm%broadcast(obs_to_proc,mpl%ioproc-1)
case ('local','adjusted')
   ! Grid-based repartition
   do iobs=1,obsop%nobs
      if (maskobs(iobs)) then
         ! Set observation proc
         if (mpl%msv%isnoti(srcproc(2,iobs)).and.(srcproc(2,iobs)==srcproc(3,iobs))) then
            ! Set to second point proc
            obs_to_proc(iobs) = srcproc(2,iobs)
         else
            ! Set to first point proc
            obs_to_proc(iobs) = srcproc(1,iobs)
         end if
      else
         ! Missing obs on main task
         obs_to_proc(iobs) = mpl%ioproc
      end if
   end do

   if (trim(nam%obsdis)=='adjusted') then
      ! Observations distribution adjusted
      write(mpl%info,'(a7,a)') '','Observations distribution adjusted'
      call mpl%flush

      ! Proc to nobsa
      proc_to_nobsa = 0
      do iobs=1,obsop%nobs
         if (maskobs(iobs)) then
            ! Concerned proc
            iproc = obs_to_proc(iobs)

            ! Number of observations per proc
            proc_to_nobsa(iproc) = proc_to_nobsa(iproc)+1
         end if
      end do

      ! Target nobsa, nobsa to move
      nres = nobs_eff
      do iproc=1,mpl%nproc
         delta = nobs_eff/mpl%nproc
         if (nres>(mpl%nproc-iproc+1)*delta) delta = delta+1
         nobs_to_move(iproc) = delta
         nres = nres-delta
      end do
      nobs_to_move = int(real(proc_to_nobsa-nobs_to_move,kind_real))
      if (sum(nobs_to_move)>0) then
         ind = maxloc(nobs_to_move)
      elseif (sum(nobs_to_move)<0) then
         ind = minloc(nobs_to_move)
      else
         ind = 1
      end if
      nobs_to_move(ind(1)) = nobs_to_move(ind(1))-sum(nobs_to_move)

      ! Select observations to move
      allocate(obs_moved(maxval(abs(nobs_to_move)),mpl%nproc))
      obs_moved = mpl%msv%vali
      nobs_to_move_tmp = nobs_to_move
      do iobs=1,obsop%nobs
         if (maskobs(iobs)) then
            iproc = obs_to_proc(iobs)
            if (nobs_to_move_tmp(iproc)>0) then
               ! Move this observation from iproc
               obs_moved(nobs_to_move_tmp(iproc),iproc) = iobs
               nobs_to_move_tmp(iproc) = nobs_to_move_tmp(iproc)-1
            end if
         end if
      end do

      ! Select destination proc
      do while (any(nobs_to_move<0))
         imin = minloc(nobs_to_move)
         imax = maxloc(nobs_to_move)
         nmoves = min(nobs_to_move(imax(1)),-nobs_to_move(imin(1)))
         do imoves=1,nmoves
            iobs = obs_moved(nobs_to_move(imax(1)),imax(1))
            obs_moved(nobs_to_move(imax(1)),imax(1)) = mpl%msv%vali
            obs_to_proc(iobs) = imin(1)
            nobs_to_move(imax(1)) = nobs_to_move(imax(1))-1
            nobs_to_move(imin(1)) = nobs_to_move(imin(1))+1
         end do
      end do

      ! Release memory
      deallocate(obs_moved)
   else
      ! Observations distribution based on grid distribution
      write(mpl%info,'(a7,a)') '','Observations distribution based on grid distribution'
      call mpl%flush
   end if
case default
   call mpl%abort(subr,'wrong obsdis')
end select

! Release memory
deallocate(obsop%lonobs)
deallocate(obsop%latobs)

! Allocation
obsop%nobsa = count(obs_to_proc==mpl%myproc)
allocate(obsop%obsa_to_obs(obsop%nobsa))
allocate(obsop%lonobs(obsop%nobsa))
allocate(obsop%latobs(obsop%nobsa))

! Fill proc_to_nobsa, obs_to_obsa and obsa_to_obs
proc_to_nobsa = 0
proc_to_nobsa_eff = 0
do iobs=1,obsop%nobs
   ! Concerned proc
   iproc = obs_to_proc(iobs)

   ! Number of observations per proc
   proc_to_nobsa(iproc) = proc_to_nobsa(iproc)+1
   if (maskobs(iobs)) proc_to_nobsa_eff(iproc) = proc_to_nobsa_eff(iproc)+1

   ! Observations local index
   iobsa = proc_to_nobsa(iproc)
   obs_to_obsa(iobs) = iobsa
   if (iproc==mpl%myproc) obsop%obsa_to_obs(iobsa) = iobs

   ! Fill lonobs/latobs
   if (iproc==mpl%myproc) then
      obsop%lonobs(iobsa) = lonobs(iobs)
      obsop%latobs(iobsa) = latobs(iobs)
   end if
end do

! Count number of local interpolation operations
obsop%h%n_s = 0
do i_s=1,hfull%n_s
   iobs = hfull%row(i_s)
   iproc = obs_to_proc(iobs)
   if (iproc==mpl%myproc) obsop%h%n_s = obsop%h%n_s+1
end do

! Count halo points
lcheck_nc0b = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   if (geom%c0_to_proc(ic0)==mpl%myproc) lcheck_nc0b(ic0) = .true.
end do
do iobs=1,obsop%nobs
   iproc = obs_to_proc(iobs)
   if (iproc==mpl%myproc) then
      do i=1,iop(iobs)
         ic0 = srcic0(i,iobs)
         lcheck_nc0b(ic0) = .true.
      end do
   end if
end do
obsop%nc0b = count(lcheck_nc0b)

! Define halo points
allocate(c0b_to_c0(obsop%nc0b))
c0_to_c0b = mpl%msv%vali
ic0b = 0
do ic0=1,geom%nc0
   if (lcheck_nc0b(ic0)) then
      ic0b = ic0b+1
      c0b_to_c0(ic0b) = ic0
      c0_to_c0b(ic0) = ic0b
   end if
end do

! Split interpolation data
obsop%h%prefix = 'o'
obsop%h%n_src = obsop%nc0b
obsop%h%n_dst = obsop%nobsa
call obsop%h%alloc
obsop%h%n_s = 0
do i_s=1,hfull%n_s
   iobs = hfull%row(i_s)
   ic0 = hfull%col(i_s)
   iproc = obs_to_proc(iobs)
   if (iproc==mpl%myproc) then
      obsop%h%n_s = obsop%h%n_s+1
      obsop%h%row(obsop%h%n_s) = obs_to_obsa(iobs)
      obsop%h%col(obsop%h%n_s) = c0_to_c0b(ic0)
      obsop%h%S(obsop%h%n_s) = hfull%S(i_s)
   end if
end do

! Inter-halo conversions
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0b = c0_to_c0b(ic0)
   c0a_to_c0b(ic0a) = ic0b
end do

! Setup communications
call obsop%com%setup(mpl,'com',geom%nc0,geom%nc0a,obsop%nc0b,c0b_to_c0,c0a_to_c0b,geom%c0_to_proc,geom%c0_to_c0a)

! Compute scores
call mpl%f_comm%allreduce(real(obsop%com%nhalo,kind_real),C_max,fckit_mpi_max())
C_max = C_max/(3.0*real(nobs_eff,kind_real)/real(mpl%nproc,kind_real))
N_max = real(maxval(proc_to_nobsa_eff),kind_real)/(real(nobs_eff,kind_real)/real(mpl%nproc,kind_real))

! Print output
write(mpl%info,'(a7,a)') '','Number of observations / valid observations per MPI task on output:'
call mpl%flush
do iproc=1,mpl%nproc
   write(mpl%info,'(a10,a,i3,a,i8,a,i8)') '','Task ',iproc,': ',proc_to_nobsa(iproc),' / ',proc_to_nobsa_eff(iproc)
   call mpl%flush
end do
write(mpl%info,'(a7,a,f5.1,a)') '','Observation repartition imbalance: ',100.0*real(maxval(proc_to_nobsa_eff) &
 & -minval(proc_to_nobsa_eff),kind_real)/(real(sum(proc_to_nobsa_eff),kind_real)/real(mpl%nproc,kind_real)),' %'
call mpl%flush
write(mpl%info,'(a7,a,i8,a,i8,a,i8)') '','Number of grid points / halo size / number of received values: ', &
 & obsop%com%nred,' / ',obsop%com%next,' / ',obsop%com%nhalo
call mpl%flush
write(mpl%info,'(a7,a,f10.2,a,f10.2)') '','Scores (N_max / C_max):',N_max,' / ',C_max
call mpl%flush

! Write observation operator
if (nam%write_obsop) call obsop%write(mpl,nam)

if (mpl%main.and.nam%write_obsop) then
   ! Write scores
   if (mpl%nproc>1) then
      call mpl%newunit(lunit)
      open(unit=lunit,file=trim(nam%datadir)//'/'//trim(nam%prefix)//'_obs_scores_'//trim(nam%obsdis)//'.dat',status='replace')
      write(lunit,*) N_max,C_max
      close(unit=lunit)
   end if
end if

! Release memory
deallocate(lonobs)
deallocate(latobs)
deallocate(maskobs)
deallocate(nop)
deallocate(iop)
deallocate(obs_to_proc)
deallocate(obs_to_obsa)
deallocate(srcproc)
deallocate(srcic0)
deallocate(c0b_to_c0)

end subroutine obsop_run_obsop

!----------------------------------------------------------------------
! Subroutine: obsop_run_obsop_tests
! Purpose: observation operator tests driver
!----------------------------------------------------------------------
subroutine obsop_run_obsop_tests(obsop,mpl,rng,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(geom_type),intent(in) :: geom       ! Geometry

! Test adjoints
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Test observation operator adjoint'
call mpl%flush
call obsop%test_adjoint(mpl,rng,geom)

if (allocated(obsop%obsa_to_obs)) then
   ! Test precision
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test observation operator precision'
   call mpl%flush
   call obsop%test_accuracy(mpl,geom)
end if

end subroutine obsop_run_obsop_tests

!----------------------------------------------------------------------
! Subroutine: obsop_apply
! Purpose: observation operator interpolation
!----------------------------------------------------------------------
subroutine obsop_apply(obsop,mpl,geom,fld,obs)

implicit none

! Passed variables
class(obsop_type),intent(in) :: obsop                    ! Observation operator data
type(mpl_type),intent(inout) :: mpl                      ! MPI data
type(geom_type),intent(in) :: geom                       ! Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0)    ! Field
real(kind_real),intent(out) :: obs(obsop%nobsa,geom%nl0) ! Observations columns

! Local variables
integer :: il0
real(kind_real) :: fld_ext(obsop%nc0b,geom%nl0)

! Halo extension
call obsop%com%ext(mpl,geom%nl0,fld,fld_ext)

if (obsop%nobsa>0) then
   ! Horizontal interpolation
   !$omp parallel do schedule(static) private(il0)
   do il0=1,geom%nl0
      call obsop%h%apply(mpl,fld_ext(:,il0),obs(:,il0),msdst=.true.)
   end do
   !$omp end parallel do
end if

end subroutine obsop_apply

!----------------------------------------------------------------------
! Subroutine: obsop_apply_ad
! Purpose: observation operator interpolation adjoint
!----------------------------------------------------------------------
subroutine obsop_apply_ad(obsop,mpl,geom,obs,fld)

implicit none

! Passed variables
class(obsop_type),intent(in) :: obsop                   ! Observation operator data
type(mpl_type),intent(inout) :: mpl                     ! MPI data
type(geom_type),intent(in) :: geom                      ! Geometry
real(kind_real),intent(in) :: obs(obsop%nobsa,geom%nl0) ! Observations columns
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0)  ! Field

! Local variables
integer :: il0
real(kind_real) :: fld_ext(obsop%nc0b,geom%nl0)

if (obsop%nobsa>0) then
   ! Horizontal interpolation
   !$omp parallel do schedule(static) private(il0)
   do il0=1,geom%nl0
      call obsop%h%apply_ad(mpl,obs(:,il0),fld_ext(:,il0))
   end do
   !$omp end parallel do
else
   ! No observation on this task
   fld_ext = 0.0
end if

! Halo reduction
call obsop%com%red(mpl,geom%nl0,fld_ext,fld)

end subroutine obsop_apply_ad

!----------------------------------------------------------------------
! Subroutine: obsop_test_adjoint
! Purpose: test observation operator adjoints accuracy
!----------------------------------------------------------------------
subroutine obsop_test_adjoint(obsop,mpl,rng,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
real(kind_real) :: sum1,sum2_loc,sum2
real(kind_real) :: fld(geom%nc0a,geom%nl0),fld_save(geom%nc0a,geom%nl0)
real(kind_real) :: yobs(obsop%nobsa,geom%nl0),yobs_save(obsop%nobsa,geom%nl0)

! Generate random fields
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)
if (obsop%nobsa>0) call rng%rand_real(0.0_kind_real,1.0_kind_real,yobs_save)

! Apply direct and adjoint obsservation operators
call obsop%apply(mpl,geom,fld_save,yobs)
call obsop%apply_ad(mpl,geom,yobs_save,fld)

! Compute adjoint test
call mpl%dot_prod(fld,fld_save,sum1)
if (obsop%nobsa>0) then
   sum2_loc = sum(yobs*yobs_save,mask=mpl%msv%isnotr(yobs))
else
   sum2_loc = 0.0
end if
call mpl%f_comm%allreduce(sum2_loc,sum2,fckit_mpi_sum())

! Print results
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Observation operator adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

end subroutine obsop_test_adjoint

!----------------------------------------------------------------------
! Subroutine: obsop_test_accuracy
! Purpose: test observation operator accuracy
!----------------------------------------------------------------------
subroutine obsop_test_accuracy(obsop,mpl,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
integer :: ic0a,ic0,iobsa
integer :: iprocmax(1),iobsamax(1)
real(kind_real) :: lonmax,latmax,ylonmax,ylatmax
real(kind_real) :: norm,distmin,distmax,distsum
real(kind_real) :: norm_tot,distmin_tot,proc_to_distmax(mpl%nproc),distsum_tot
real(kind_real) :: lon(geom%nc0a,geom%nl0),lat(geom%nc0a,geom%nl0)
real(kind_real) :: ylon(obsop%nobsa,geom%nl0),ylat(obsop%nobsa,geom%nl0)
real(kind_real) :: dist(obsop%nobsa)
character(len=1024),parameter :: subr = 'obsop_test_accuracy'

! Initialization
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   lon(ic0a,:) = geom%lon(ic0)
   lat(ic0a,:) = geom%lat(ic0)
end do

! Apply obsop
call obsop%apply(mpl,geom,lon,ylon)
call obsop%apply(mpl,geom,lat,ylat)

if (obsop%nobsa>0) then
   ! Remove points close to the longitude discontinuity and to the poles
   dist = mpl%msv%valr
   do iobsa=1,obsop%nobsa
      if (mpl%msv%isnotr(ylon(iobsa,1)).and.mpl%msv%isnotr(ylat(iobsa,1))) then
         if ((abs(obsop%lonobs(iobsa))<0.8*pi).and.(abs(obsop%latobs(iobsa))<0.4*pi)) then
            call sphere_dist(ylon(iobsa,1),ylat(iobsa,1),obsop%lonobs(iobsa),obsop%latobs(iobsa),dist(iobsa))
            dist(iobsa) = dist(iobsa)*reqkm
         end if
      end if
   end do
   norm = real(count(mpl%msv%isnotr(dist)),kind_real)
   if (norm>0) then
      distmin = minval(dist,mask=mpl%msv%isnotr(dist))
      distmax = maxval(dist,mask=mpl%msv%isnotr(dist))
      distsum = sum(dist,mask=mpl%msv%isnotr(dist))
   else
      distmin = huge(1.0)
      distmax = 0.0
      distsum = 0.0
   end if
else
   ! No observation on this task
   norm = 0
   distmin = huge(1.0)
   distmax = 0.0
   distsum = 0.0
end if

! Gather results
call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
call mpl%f_comm%allreduce(distmin,distmin_tot,fckit_mpi_min())
call mpl%f_comm%allgather(distmax,proc_to_distmax)
call mpl%f_comm%allreduce(distsum,distsum_tot,fckit_mpi_sum())

! Maximum error detail
iprocmax = maxloc(proc_to_distmax)
if (iprocmax(1)==mpl%myproc) then
   iobsamax = maxloc(dist,mask=mpl%msv%isnotr(dist))
   lonmax = obsop%lonobs(iobsamax(1))
   latmax = obsop%latobs(iobsamax(1))
   ylonmax = ylon(iobsamax(1),1)
   ylatmax = ylat(iobsamax(1),1)
end if

! Broadcast results
call mpl%f_comm%broadcast(lonmax,iprocmax(1)-1)
call mpl%f_comm%broadcast(latmax,iprocmax(1)-1)
call mpl%f_comm%broadcast(ylonmax,iprocmax(1)-1)
call mpl%f_comm%broadcast(ylatmax,iprocmax(1)-1)

! Print results
if (norm_tot>0.0) then
   write(mpl%info,'(a7,a,f10.2,a,f10.2,a,f10.2,a)') '','Interpolation error (min/mean/max): ',distmin_tot, &
 & ' km / ',distsum_tot/norm_tot,' km / ',maxval(proc_to_distmax),' km'
   call mpl%flush
   write(mpl%info,'(a7,a)') '','Max. interpolation error location (lon/lat): '
   call mpl%flush
   write(mpl%info,'(a10,a14,f10.2,a,f10.2,a)') '','Observation:  ',lonmax*rad2deg,' deg. / ' ,latmax*rad2deg,' deg.'
   call mpl%flush
   write(mpl%info,'(a10,a14,f10.2,a,f10.2,a)') '','Interpolation:',ylonmax*rad2deg,' deg. / ' ,ylatmax*rad2deg,' deg.'
   call mpl%flush
else
   call mpl%abort(subr,'all observations are out of the test windows')
end if

end subroutine obsop_test_accuracy

end module type_obsop
