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

use module_namelist, only: namtype
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi
use tools_qsort, only: qsort
use type_com, only: comtype,com_copy,com_setup
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_alloc
use type_mpl, only: mpl
use type_odata, only: odatatype,odataloctype
use type_randgen, only: rng,rand_real

implicit none

private
public :: compute_mpi_obsop

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi_obsop
!> Purpose: compute observation operator MPI distribution
!----------------------------------------------------------------------
subroutine compute_mpi_obsop(nam,odata,odataloc)

implicit none

! Passed variables
type(namtype),intent(in) :: nam
type(odatatype),intent(inout) :: odata
type(odataloctype),intent(inout) :: odataloc

! Local variables
integer :: iobs,jobs,iobsa,iproc,nobsa,i_s,ic0,ic0b,i,jproc,ic0a,nc0a,nc0b
integer :: ic0_to_ic0b(odata%nc0)
integer,allocatable :: iop(:),srcproc(:,:),srcic0(:,:),order(:)
real(kind_real),allocatable :: list(:)
logical :: lcheck_nc0b(odata%nc0)
type(comtype) :: comobs(nam%nproc)

! Allocation
allocate(iop(odata%nobs))
allocate(srcproc(3,odata%nobs))
allocate(srcic0(3,odata%nobs))
allocate(order(odata%nobs))
allocate(odata%iobs_to_iproc(odata%nobs))
allocate(odata%iobs_to_iobsa(odata%nobs))
allocate(odata%iproc_to_nobsa(nam%nproc))

! Find grid points origin
iop = 0
call msi(srcproc)
call msi(srcic0) 
do i_s=1,odata%interp%n_s
   ic0 = odata%interp%col(i_s)
   iproc = odata%geom%ic0_to_iproc(ic0)
   iobs = odata%interp%row(i_s)
   iop(iobs) = iop(iobs)+1
   srcproc(iop(iobs),iobs) = iproc
   srcic0(iop(iobs),iobs) = ic0
end do

! Generate observation distribution on processors
if (.false.) then
   ! Random repartition
   call rand_real(rng,0.0_kind_real,1.0_kind_real,.true.,list)
   call qsort(odata%nobs,list,order)
   nobsa = odata%nobs/nam%nproc
   if (nobsa*nam%nproc<odata%nobs) nobsa = nobsa+1 
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

do iproc=1,nam%nproc
   ! Count halo points
   lcheck_nc0b = .false.
   do ic0=1,odata%nc0
      jproc = odata%geom%ic0_to_iproc(ic0)
      if (iproc==jproc) lcheck_nc0b(ic0) = .true.
   end do
   do iobs=1,odata%nobs
      jproc = odata%iobs_to_iproc(iobs)
      if (iproc==jproc) then
         do i=1,iop(iobs)
            ic0 = srcic0(i,iobs)
            lcheck_nc0b(ic0) = .true.
         end do
      end if
   end do
   nc0a = count(odata%geom%ic0_to_iproc==iproc)
   nc0b = count(lcheck_nc0b)

   ! Communication
   comobs(iproc)%prefix = 'comobs'
   comobs(iproc)%nred = nc0a
   comobs(iproc)%next = nc0b

   ! Allocation
   allocate(comobs(iproc)%iext_to_iproc(comobs(iproc)%next))
   allocate(comobs(iproc)%iext_to_ired(comobs(iproc)%next))
   allocate(comobs(iproc)%ired_to_iext(comobs(iproc)%nred))

   ! Define halo origin
   call msi(ic0_to_ic0b)
   ic0b = 0
   do ic0=1,odata%nc0
      if (lcheck_nc0b(ic0)) then
         ic0b = ic0b+1
         comobs(iproc)%iext_to_iproc(ic0b) = odata%geom%ic0_to_iproc(ic0)
         ic0a = odata%geom%ic0_to_ic0a(ic0)
         comobs(iproc)%iext_to_ired(ic0b) = ic0a
         jproc = odata%geom%ic0_to_iproc(ic0)
         if (iproc==jproc) comobs(iproc)%ired_to_iext(ic0a) = ic0b
         ic0_to_ic0b(ic0) = ic0b
      end if
   end do

   if (iproc==mpl%myproc) then
      ! Parameters
      odataloc%nc0a = nc0a
      odataloc%nc0b = nc0b
      odataloc%nl0 = odata%nl0
      odataloc%nobsa = odata%iproc_to_nobsa(iproc)

      ! Split interpolation data
      odataloc%interp%prefix = 'o'
      odataloc%interp%n_src = comobs(iproc)%next
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
   end if
end do

! Communications setup
call com_setup(nam,comobs)

! Communications copy
call com_copy(nam,comobs(mpl%myproc),odataloc%com)
odataloc%com%prefix = 'o'

! Print results
write(mpl%unit,'(a7,a)') '','Number of observations per MPI task:'
do iproc=1,nam%nproc
   write(mpl%unit,'(a10,a,i3,a,i8)') '','Task ',iproc,': ',count(odata%iobs_to_iproc==iproc)
end do
write(mpl%unit,'(a7,a,f5.1,a)') '','Observation repartition imbalance: ', &
 & 100.0*float(maxval(odata%iproc_to_nobsa)-minval(odata%iproc_to_nobsa))/(float(sum(odata%iproc_to_nobsa))/float(nam%nproc)),' %'
write(mpl%unit,'(a7,a)') '','Number of grid points, halo size and number of received values per MPI task:'
do iproc=1,nam%nproc
   write(mpl%unit,'(a10,a,i3,a,i8,a,i8,a,i8)') '','Task ',iproc,': ', &
 & comobs(iproc)%nred,' / ',comobs(iproc)%next,' / ',comobs(iproc)%nhalo
end do

end subroutine compute_mpi_obsop

end module module_mpi_obsop
