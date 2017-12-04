!----------------------------------------------------------------------
! Module: module_mpi_obsop.f90
!> Purpose: compute observation operator MPI distribution
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_mpi_obsop

use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi
use tools_qsort, only: qsort
use type_com, only: comtype,com_setup,com_bcast
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_alloc
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_bcast
use type_nam, only: namtype
use type_odata, only: odatatype,odataloctype
use type_randgen, only: rand_real

implicit none

private
public :: compute_mpi_obsop

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi_obsop
!> Purpose: compute observation operator MPI distribution
!----------------------------------------------------------------------
subroutine compute_mpi_obsop(odata,odataloc)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata       !< Observation operator data
type(odataloctype),intent(inout) :: odataloc !< Observation operator data, local

! Local variables
integer :: iobs,jobs,iobsa,iproc,nobsa,i_s,ic0,ic0b,i,jproc,ic0a,nc0a,nc0b
integer :: ic0_to_ic0b(odata%geom%nc0)
integer,allocatable :: iop(:),srcproc(:,:),srcic0(:,:),order(:),ic0b_to_ic0(:),ic0b_to_ic0_copy(:)
real(kind_real),allocatable :: list(:)
logical :: lcheck_nc0b(odata%geom%nc0)
type(comtype) :: comobs(mpl%nproc)

! Associate
associate(nam=>odata%nam,geom=>odata%geom)

! Allocation
allocate(iop(odata%nobs))
allocate(srcproc(3,odata%nobs))
allocate(srcic0(3,odata%nobs))
allocate(order(odata%nobs))
allocate(odata%iobs_to_iproc(odata%nobs))
allocate(odata%iobs_to_iobsa(odata%nobs))
allocate(odata%iproc_to_nobsa(mpl%nproc))

! Find grid points origin
iop = 0
call msi(srcproc)
call msi(srcic0)
do i_s=1,odata%interp%n_s
   ic0 = odata%interp%col(i_s)
   iproc = geom%ic0_to_iproc(ic0)
   iobs = odata%interp%row(i_s)
   iop(iobs) = iop(iobs)+1
   srcproc(iop(iobs),iobs) = iproc
   srcic0(iop(iobs),iobs) = ic0
end do

! Generate observation distribution on processors
if (.false.) then
   ! Random repartition
   if (mpl%main) call rand_real(0.0_kind_real,1.0_kind_real,list)
   call mpl_bcast(list,mpl%ioproc)
   call qsort(odata%nobs,list,order)
   nobsa = odata%nobs/mpl%nproc
   if (nobsa*mpl%nproc<odata%nobs) nobsa = nobsa+1
   iproc = 1
   iobsa = 1
   do iobs=1,odata%nobs
      jobs = order(iobs)
      odata%iobs_to_iproc(jobs) = iproc
      iobsa = iobsa+1
      if (iobsa>nobsa) then
         iproc = iproc+1
         iobsa = 1
      end if
   end do
else
   ! Source grid-based repartition
   do iobs=1,odata%nobs
      ! Set observation proc
      if (srcproc(2,iobs)==srcproc(3,iobs)) then
         ! Set to second point proc
         odata%iobs_to_iproc(iobs) = srcproc(2,iobs)
      else
         ! Set to first point proc
         odata%iobs_to_iproc(iobs) = srcproc(1,iobs)
      end if
   end do
end if

! Local observations
odata%iproc_to_nobsa = 0
do iobs=1,odata%nobs
   ! Concerned proc
   iproc = odata%iobs_to_iproc(iobs)

   ! Number of observations per proc
   odata%iproc_to_nobsa(iproc) = odata%iproc_to_nobsa(iproc)+1

   ! Observations local index
   odata%iobs_to_iobsa(iobs) = odata%iproc_to_nobsa(iproc)
end do

! Count number of local interpolation operations
do i_s=1,odata%interp%n_s
   iobs = odata%interp%row(i_s)
   iproc = odata%iobs_to_iproc(iobs)
   if (iproc==mpl%myproc) odataloc%interp%n_s = odataloc%interp%n_s+1
end do

! Count halo points
lcheck_nc0b = .false.
do ic0=1,geom%nc0
   jproc = geom%ic0_to_iproc(ic0)
   if (jproc==mpl%myproc) lcheck_nc0b(ic0) = .true.
end do
do iobs=1,odata%nobs
   jproc = odata%iobs_to_iproc(iobs)
   if (jproc==mpl%myproc) then
      do i=1,iop(iobs)
         ic0 = srcic0(i,iobs)
         lcheck_nc0b(ic0) = .true.
      end do
   end if
end do
odataloc%nc0a = count(geom%ic0_to_iproc==mpl%myproc)
odataloc%nc0b = count(lcheck_nc0b)
odataloc%nobsa = odata%iproc_to_nobsa(mpl%myproc)

! Define halo points
allocate(ic0b_to_ic0(odataloc%nc0b))
call msi(ic0_to_ic0b)
ic0b = 0
do ic0=1,geom%nc0
   if (lcheck_nc0b(ic0)) then
      ic0b = ic0b+1
      ic0_to_ic0b(ic0) = ic0b
      ic0b_to_ic0(ic0b) = ic0
   end if
end do

! Split interpolation data
odataloc%interp%prefix = 'o'
odataloc%interp%n_src = odataloc%nc0b
odataloc%interp%n_dst = odataloc%nobsa
call linop_alloc(odataloc%interp)
odataloc%interp%n_s = 0
do i_s=1,odata%interp%n_s
   iobs = odata%interp%row(i_s)
   jproc = odata%iobs_to_iproc(iobs)
   if (iproc==jproc) then
      odataloc%interp%n_s = odataloc%interp%n_s+1
      odataloc%interp%row(odataloc%interp%n_s) = odata%iobs_to_iobsa(iobs)
      odataloc%interp%col(odataloc%interp%n_s) = ic0_to_ic0b(odata%interp%col(i_s))
      odataloc%interp%S(odataloc%interp%n_s) = odata%interp%S(i_s)
   end if
end do

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = odataloc%nc0a
         nc0b = odataloc%nc0b
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
         ic0a = geom%ic0_to_ic0a(ic0)
         comobs(iproc)%iext_to_iproc(ic0b) = geom%ic0_to_iproc(ic0)
         comobs(iproc)%iext_to_ired(ic0b) = ic0a
         comobs(iproc)%ired_to_iext(ic0a) = ic0b
      end do

      ! Release memory
      deallocate(ic0b_to_ic0)
   end do

   ! Communications setup
   call com_setup(mpl%nproc,comobs)
else
   ! Send dimensions to ioproc
   call mpl_send(odataloc%nc0a,mpl%ioproc,mpl%tag)
   call mpl_send(odataloc%nc0b,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl_send(odataloc%nc0b,ic0b_to_ic0,mpl%ioproc,mpl%tag+2)
end if
mpl%tag = mpl%tag+3

! Communications broadcast
odataloc%com%prefix = 'o'
call com_bcast(mpl%nproc,comobs,odataloc%com)

! Print results
write(mpl%unit,'(a7,a)') '','Number of observations per MPI task:'
do iproc=1,mpl%nproc
   write(mpl%unit,'(a10,a,i3,a,i8)') '','Task ',iproc,': ',count(odata%iobs_to_iproc==iproc)
end do
write(mpl%unit,'(a7,a,f5.1,a)') '','Observation repartition imbalance: ', &
 & 100.0*float(maxval(odata%iproc_to_nobsa)-minval(odata%iproc_to_nobsa))/(float(sum(odata%iproc_to_nobsa))/float(mpl%nproc)),' %'
write(mpl%unit,'(a7,a)') '','Number of grid points, halo size and number of received values per MPI task:'
do iproc=1,mpl%nproc
   write(mpl%unit,'(a10,a,i3,a,i8,a,i8,a,i8)') '','Task ',iproc,': ', &
 & comobs(iproc)%nred,' / ',comobs(iproc)%next,' / ',comobs(iproc)%nhalo
end do

! End associate
end associate

end subroutine compute_mpi_obsop

end module module_mpi_obsop
