subroutine split_interp_linop(hinterp_glb, hinterp_loc, kobs, iobs_to_iproc, kpts, ipts_to_iproc, nhalo, obscom)

use tools_kinds, only: kind_real
use tools_missing, only: msi
use type_com, only: comtype,com_setup,com_bcast
use type_linop, only: linoptype,linop_alloc
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast

implicit none

! Passed variables
type(linoptype), intent(in) :: hinterp_glb
type(linoptype), intent(inout) :: hinterp_loc
integer, intent(in) :: kobs  ! Total number of obs
integer, intent(in) :: iobs_to_iproc(kobs)  ! Index of task that owns each obs
integer, intent(in) :: kpts  ! Total number of grid points
integer, intent(in) :: ipts_to_iproc(kpts)  ! Index of task that owns each grid point
integer, intent(out) :: nhalo  ! Number of grid points local with halo
type(comtype), intent(inout) :: obscom

! Local variables
integer :: iobs,jobs,iobsa,iproc,nobsa,i_s,ic0,ic0b,i,jproc,ic0a,nc0a,nc0b,nptsloc
integer :: ic0_to_ic0b(kpts), ic0_to_ic0a(kpts)
integer :: iop(kobs),srcic0(3,kobs)
integer,allocatable :: ic0b_to_ic0(:),ic0b_to_ic0_copy(:)
logical :: lcheck_nc0b(kpts)
type(comtype) :: comobs(mpl%nproc)
integer :: iobs_to_iobsa(kobs)
integer :: iproc_to_nobsa(mpl%nproc), indx(mpl%nproc)

! Find grid points origin
iop = 0
call msi(srcic0)
do i_s=1,hinterp_glb%n_s
   ic0 = hinterp_glb%col(i_s)
   iproc = ipts_to_iproc(ic0)
   iobs = hinterp_glb%row(i_s)
   iop(iobs) = iop(iobs)+1
   srcic0(iop(iobs),iobs) = ic0
end do

! Local observations
iproc_to_nobsa(:) = 0
do iobs=1,kobs
   ! Concerned proc
   iproc = iobs_to_iproc(iobs)

   ! Number of observations per proc
   iproc_to_nobsa(iproc) = iproc_to_nobsa(iproc)+1

   ! Observations local index
   iobs_to_iobsa(iobs) = iproc_to_nobsa(iproc)
end do

! Count number of local interpolation operations
do i_s=1,hinterp_glb%n_s
   iobs = hinterp_glb%row(i_s)
   iproc = iobs_to_iproc(iobs)
   if (iproc==mpl%myproc) hinterp_loc%n_s = hinterp_loc%n_s+1
end do

! Count halo points
indx(:)=0
lcheck_nc0b = .false.
do ic0=1,kpts
   jproc = ipts_to_iproc(ic0)
   if (jproc==mpl%myproc) lcheck_nc0b(ic0) = .true.
   indx(jproc)=indx(jproc)+1
   ic0_to_ic0a(ic0)=indx(jproc)
end do
do iobs=1,kobs
   jproc = iobs_to_iproc(iobs)
   if (jproc==mpl%myproc) then
      do i=1,iop(iobs)
         ic0 = srcic0(i,iobs)
         lcheck_nc0b(ic0) = .true.
      end do
   end if
end do
nptsloc = count(ipts_to_iproc==mpl%myproc)
nhalo = count(lcheck_nc0b)

! Define halo points
allocate(ic0b_to_ic0(nhalo))
call msi(ic0_to_ic0b)
ic0b = 0
do ic0=1,kpts
   if (lcheck_nc0b(ic0)) then
      ic0b = ic0b+1
      ic0_to_ic0b(ic0) = ic0b
      ic0b_to_ic0(ic0b) = ic0
   end if
end do

! Split interpolation data
hinterp_loc%prefix = 'o'
hinterp_loc%n_src = nhalo
hinterp_loc%n_dst = iproc_to_nobsa(mpl%myproc)
call linop_alloc(hinterp_loc)
hinterp_loc%n_s = 0
do i_s=1,hinterp_glb%n_s
   iobs = hinterp_glb%row(i_s)
   jproc = iobs_to_iproc(iobs)
   if (iproc==jproc) then
      hinterp_loc%n_s = hinterp_loc%n_s+1
      hinterp_loc%row(hinterp_loc%n_s) = iobs_to_iobsa(iobs)
      hinterp_loc%col(hinterp_loc%n_s) = ic0_to_ic0b(hinterp_glb%col(i_s))
      hinterp_loc%S(hinterp_loc%n_s) = hinterp_glb%S(i_s)
   end if
end do

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = nptsloc
         nc0b = nhalo
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
         call mpl_recv(nc0b,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(ic0b_to_ic0_copy(nc0b))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         ic0b_to_ic0_copy = ic0b_to_ic0
      else
         ! Receive data on ioproc
         call mpl_recv(nc0b,ic0b_to_ic0_copy,iproc,mpl%tag+2)
      end if

      ! Allocation
      comobs(iproc)%nred = nc0a
      comobs(iproc)%next = nc0b
      allocate(comobs(iproc)%iext_to_iproc(comobs(iproc)%next))
      allocate(comobs(iproc)%iext_to_ired(comobs(iproc)%next))
      allocate(comobs(iproc)%ired_to_iext(comobs(iproc)%nred))

      ! Define halo origin
      do ic0b=1,nc0b
         ic0 = ic0b_to_ic0_copy(ic0b)
         ic0a = ic0_to_ic0a(ic0)
         comobs(iproc)%iext_to_iproc(ic0b) = ipts_to_iproc(ic0)
         comobs(iproc)%iext_to_ired(ic0b) = ic0a
         comobs(iproc)%ired_to_iext(ic0a) = ic0b
      end do

      ! Release memory
      deallocate(ic0b_to_ic0_copy)
   end do

   ! Communications setup
   call com_setup(mpl%nproc,comobs)
else
   ! Send dimensions to ioproc
   call mpl_send(nptsloc,mpl%ioproc,mpl%tag)
   call mpl_send(nhalo,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl_send(nhalo,ic0b_to_ic0,mpl%ioproc,mpl%tag+2)
end if
mpl%tag = mpl%tag+3
deallocate(ic0b_to_ic0)

! Communications broadcast
obscom%prefix = 'o'
call com_bcast(mpl%nproc,comobs,obscom)

! Print results
write(mpl%unit,'(a7,a)') '','Number of observations per MPI task:'
do iproc=1,mpl%nproc
   write(mpl%unit,'(a10,a,i3,a,i8)') '','Task ',iproc,': ',count(iobs_to_iproc==iproc)
end do
write(mpl%unit,'(a7,a,f5.1,a)') '','Observation repartition imbalance: ', &
 & 100.0*float(maxval(iproc_to_nobsa)-minval(iproc_to_nobsa))/(float(sum(iproc_to_nobsa))/float(mpl%nproc)),' %'
write(mpl%unit,'(a7,a)') '','Number of grid points, halo size and number of received values per MPI task:'
do iproc=1,mpl%nproc
   write(mpl%unit,'(a10,a,i3,a,i8,a,i8,a,i8)') '','Task ',iproc,': ', &
 & comobs(iproc)%nred,' / ',comobs(iproc)%next,' / ',comobs(iproc)%nhalo
end do

end subroutine split_interp_linop
