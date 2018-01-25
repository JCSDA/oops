!----------------------------------------------------------------------
! Module: nicas_parameters_obsop.f90
!> Purpose: compute observation operator interpolation
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module obsop_parameters

use tools_const, only: rad2deg
use tools_display, only: msgerror,newunit
use tools_interp, only: compute_interp
use tools_kinds, only: kind_real
use tools_missing, only: msi,isnotmsi
use tools_qsort, only: qsort
use type_com, only: comtype,com_setup,com_bcast,com_dealloc
use type_ctree, only: ctreetype,create_ctree
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_alloc,linop_reorder,linop_dealloc
use type_mesh, only: create_mesh
use type_mpl, only: mpl,mpl_allgather,mpl_send,mpl_recv,mpl_bcast
use type_nam, only: namtype
use type_odata, only: odatatype
use type_randgen, only: rand_real

implicit none

private
public :: compute_parameters

contains

!----------------------------------------------------------------------
! Subroutine: compute_parameters
!> Purpose: compute observation operator interpolation
!----------------------------------------------------------------------
subroutine compute_parameters(odata,local)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata !< Observation operator data
logical,intent(in),optional :: local   !< Local input data (default = .false.)

! Local variables
integer :: offset,iobs,jobs,iobsa,iproc,jproc,nobsa,i_s,ic0,ic0b,i,ic0a,nc0a,nc0b,nhalo,delta,nres,ind(1),lunit
integer,allocatable :: nop(:),iop(:),srcproc(:,:),srcic0(:,:),order(:),nobsa_to_move(:)
integer,allocatable :: c0_to_c0a(:),c0a_to_c0(:),c0b_to_c0(:),c0a_to_c0b(:)
real(kind_real),allocatable :: lonobs(:),latobs(:),list(:)
logical :: llocal
logical,allocatable :: maskobs(:),lcheck_nc0b(:)
type(comtype) :: comobs(mpl%nproc)

! Associate
associate(geom=>odata%geom,nam=>odata%nam)

! Default: global input data
llocal = .false.
if (present(local)) llocal = local

if (llocal) then
   ! Allocation
   allocate(odata%proc_to_nobsa(mpl%nproc))

   ! Get global number of observations
   call mpl_allgather(1,(/odata%nobsa/),odata%proc_to_nobsa)
   odata%nobs = sum(odata%proc_to_nobsa)

   ! Allocation
   allocate(lonobs(odata%nobs))
   allocate(latobs(odata%nobs))

   ! Get observations coordinates
   if (mpl%main) then
      offset = 0
      do iproc=1,mpl%nproc
         if (odata%proc_to_nobsa(iproc)>0) then
            if (iproc==mpl%ioproc) then
               ! Copy data
               lonobs(offset+1:offset+odata%proc_to_nobsa(iproc)) = odata%lonobs
               latobs(offset+1:offset+odata%proc_to_nobsa(iproc)) = odata%latobs
            else
               ! Receive data on ioproc
               call mpl_recv(odata%proc_to_nobsa(iproc),lonobs(offset+1:offset+odata%proc_to_nobsa(iproc)),iproc,mpl%tag)
               call mpl_recv(odata%proc_to_nobsa(iproc),latobs(offset+1:offset+odata%proc_to_nobsa(iproc)),iproc,mpl%tag+1)
            end if
         end if

         ! Update offset
         offset = offset+odata%proc_to_nobsa(iproc)
      end do
   else
      if (odata%nobsa>0) then
         ! Send data to ioproc
         call mpl_send(odata%nobsa,odata%lonobs,mpl%ioproc,mpl%tag)
         call mpl_send(odata%nobsa,odata%latobs,mpl%ioproc,mpl%tag+1)
      end if
   end if
   mpl%tag = mpl%tag+2

   ! Broadcast data
   call mpl_bcast(lonobs,mpl%ioproc)
   call mpl_bcast(latobs,mpl%ioproc)
else
   ! Allocation
   allocate(lonobs(odata%nobs))
   allocate(latobs(odata%nobs))

   ! Copy coordinates
   lonobs = odata%lonobs
   latobs = odata%latobs
end if

! Allocation
allocate(maskobs(odata%nobs))
allocate(odata%obs_to_proc(odata%nobs))
allocate(odata%obs_to_obsa(odata%nobs))

! Compute interpolation
odata%hfull%prefix = 'o'
write(mpl%unit,'(a7,a)') '','Single level:'
maskobs = .true.
call compute_interp(geom%mesh,geom%ctree,geom%nc0,any(geom%mask,dim=2),odata%nobs,lonobs,latobs, &
 & maskobs,nam%obsop_interp,odata%hfull)

if (llocal) then
   ! Fill obs_to_proc and obs_to_obsa
   iobs = 0
   do iproc=1,mpl%nproc
      do iobsa=1,odata%proc_to_nobsa(iproc)
         iobs = iobs+1
         odata%obs_to_proc(iobs) = iproc
         odata%obs_to_obsa(iobs) = iobsa
      end do
   end do
else
   ! Split observations between processors

   ! Allocation
   allocate(nop(odata%nobs))
   allocate(iop(odata%nobs))
   allocate(odata%proc_to_nobsa(mpl%nproc))
   nop = 0
   do i_s=1,odata%hfull%n_s
      iobs = odata%hfull%row(i_s)
      nop(iobs) = nop(iobs)+1
   end do
   allocate(srcproc(maxval(nop),odata%nobs))
   allocate(srcic0(maxval(nop),odata%nobs))
   
   ! Find grid points origin
   iop = 0
   call msi(srcproc)
   call msi(srcic0)
   do i_s=1,odata%hfull%n_s
      ic0 = odata%hfull%col(i_s)
      iproc = geom%c0_to_proc(ic0)
      iobs = odata%hfull%row(i_s)
      iop(iobs) = iop(iobs)+1
      srcproc(iop(iobs),iobs) = iproc
      srcic0(iop(iobs),iobs) = ic0
   end do
   
   ! Generate observation distribution on processors
   if (nam%obsdis<0.0) then
      ! Allocation
      allocate(list(odata%nobs))
      allocate(order(odata%nobs))
   
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
         odata%obs_to_proc(jobs) = iproc
         iobsa = iobsa+1
         if (iobsa>nobsa) then
            iproc = iproc+1
            iobsa = 1
         end if
      end do
   
      ! Release memory
      deallocate(list)
      deallocate(order)
   else
      ! Source grid-based repartition
      do iobs=1,odata%nobs
         ! Set observation proc
         if (srcproc(2,iobs)==srcproc(3,iobs)) then
            ! Set to second point proc
            odata%obs_to_proc(iobs) = srcproc(2,iobs)
         else
            ! Set to first point proc
            odata%obs_to_proc(iobs) = srcproc(1,iobs)
         end if
      end do
   
      odata%proc_to_nobsa = 0
      do iobs=1,odata%nobs
         ! Concerned proc
         iproc = odata%obs_to_proc(iobs)
   
         ! Number of observations per proc
         odata%proc_to_nobsa(iproc) = odata%proc_to_nobsa(iproc)+1
      end do
   
      ! Allocation
      allocate(nobsa_to_move(mpl%nproc))
   
      ! Target nobsa, nobsa to move
      nres = odata%nobs
      do iproc=1,mpl%nproc
         delta = odata%nobs/mpl%nproc
         if (nres>(mpl%nproc-iproc+1)*delta) delta = delta+1
         nobsa_to_move(iproc) = delta
         nres = nres-delta
      end do
      nobsa_to_move = int(nam%obsdis*float(odata%proc_to_nobsa-nobsa_to_move))
      if (sum(nobsa_to_move)>0) then
         ind = maxloc(nobsa_to_move)
      elseif (sum(nobsa_to_move)<0) then
         ind = minloc(nobsa_to_move)
      else
         ind = 1
      end if
      nobsa_to_move(ind(1)) = nobsa_to_move(ind(1))-sum(nobsa_to_move)
   
      ! Move observations between processors
      do iobs=1,odata%nobs
         iproc = odata%obs_to_proc(iobs)
         if (nobsa_to_move(iproc)>0) then
            ! Move this observation from iproc
            do jproc=1,mpl%nproc
               if (nobsa_to_move(jproc)<0) then
                  ! Move this observation to jproc
                  odata%obs_to_proc(iobs) = jproc
                  nobsa_to_move(iproc) = nobsa_to_move(iproc)-1
                  nobsa_to_move(jproc) = nobsa_to_move(jproc)+1
                  exit
               end if
            end do
         end if
      end do
   
      ! Release memory
      deallocate(nobsa_to_move)
   end if
   
   ! Local observations
   odata%proc_to_nobsa = 0
   do iobs=1,odata%nobs
      ! Concerned proc
      iproc = odata%obs_to_proc(iobs)
   
      ! Number of observations per proc
      odata%proc_to_nobsa(iproc) = odata%proc_to_nobsa(iproc)+1
   
      ! Observations local index
      odata%obs_to_obsa(iobs) = odata%proc_to_nobsa(iproc)
   end do
end if

! Allocation
allocate(lcheck_nc0b(geom%nc0))

! Count number of local interpolation operations
odata%h%n_s = 0
do i_s=1,odata%hfull%n_s
   iobs = odata%hfull%row(i_s)
   iproc = odata%obs_to_proc(iobs)
   if (iproc==mpl%myproc) odata%h%n_s = odata%h%n_s+1
end do

! Count halo points
lcheck_nc0b = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   if (geom%c0_to_proc(ic0)==mpl%myproc) lcheck_nc0b(ic0) = .true.
end do
do iobs=1,odata%nobs
   iproc = odata%obs_to_proc(iobs)
   if (iproc==mpl%myproc) then
      do i=1,iop(iobs)
         ic0 = srcic0(i,iobs)
         lcheck_nc0b(ic0) = .true.
      end do
   end if
end do
odata%nc0a = count(geom%c0_to_proc==mpl%myproc)
odata%nc0b = count(lcheck_nc0b)
odata%nobsa = odata%proc_to_nobsa(mpl%myproc)

! Define halo points
allocate(odata%c0b_to_c0(odata%nc0b))
allocate(odata%c0_to_c0b(geom%nc0))
call msi(odata%c0_to_c0b)
ic0b = 0
do ic0=1,geom%nc0
   if (lcheck_nc0b(ic0)) then
      ic0b = ic0b+1
      odata%c0b_to_c0(ic0b) = ic0
      odata%c0_to_c0b(ic0) = ic0b
   end if
end do

! Split interpolation data
odata%h%prefix = 'o'
odata%h%n_src = odata%nc0b
odata%h%n_dst = odata%nobsa
call linop_alloc(odata%h)
odata%h%n_s = 0
do i_s=1,odata%hfull%n_s
   iobs = odata%hfull%row(i_s)
   ic0 = odata%hfull%col(i_s)
   iproc = odata%obs_to_proc(iobs)
   if (iproc==mpl%myproc) then
      odata%h%n_s = odata%h%n_s+1
      odata%h%row(odata%h%n_s) = odata%obs_to_obsa(iobs)
      odata%h%col(odata%h%n_s) = odata%c0_to_c0b(ic0)
      odata%h%S(odata%h%n_s) = odata%hfull%S(i_s)
   end if
end do

! Inter-halo conversions
allocate(odata%c0a_to_c0b(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0b = odata%c0_to_c0b(ic0)
   odata%c0a_to_c0b(ic0a) = ic0b
end do

! Get global distribution of the subgrid on ioproc
if (mpl%main) then
   ! Allocation
   allocate(c0_to_c0a(geom%nc0))

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         nc0a = geom%nc0a
      else
         ! Receive dimension on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
      end if

      ! Allocation
      allocate(c0a_to_c0(nc0a))

      if (iproc==mpl%ioproc) then
         ! Copy data
         c0a_to_c0 = geom%c0a_to_c0
      else
         ! Receive data on ioproc
         call mpl_recv(nc0a,c0a_to_c0,iproc,mpl%tag+1)
      end if

      ! Fill c0_to_c0a
      do ic0a=1,nc0a
         c0_to_c0a(c0a_to_c0(ic0a)) = ic0a
      end do

      ! Release memory
      deallocate(c0a_to_c0)
   end do
else
   ! Send dimensions to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)

   ! Send data to ioproc
   call mpl_send(geom%nc0a,geom%c0a_to_c0,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = odata%nc0a
         nc0b = odata%nc0b
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
         call mpl_recv(nc0b,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(c0b_to_c0(nc0b))
      allocate(c0a_to_c0b(nc0a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0b_to_c0 = odata%c0b_to_c0
         c0a_to_c0b = odata%c0a_to_c0b
      else
         ! Receive data on ioproc
         call mpl_recv(nc0b,c0b_to_c0,iproc,mpl%tag+2)
         call mpl_recv(nc0a,c0a_to_c0b,iproc,mpl%tag+3)
      end if

      ! Allocation
      comobs(iproc)%nred = nc0a
      comobs(iproc)%next = nc0b
      allocate(comobs(iproc)%ext_to_proc(comobs(iproc)%next))
      allocate(comobs(iproc)%ext_to_red(comobs(iproc)%next))
      allocate(comobs(iproc)%red_to_ext(comobs(iproc)%nred))

      ! Define halo origin
      do ic0b=1,nc0b
         ic0 = c0b_to_c0(ic0b)
         comobs(iproc)%ext_to_proc(ic0b) = geom%c0_to_proc(ic0)
         ic0a = c0_to_c0a(ic0)
         comobs(iproc)%ext_to_red(ic0b) = ic0a
      end do
      comobs(iproc)%red_to_ext = c0a_to_c0b

      ! Release memory
      deallocate(c0b_to_c0)
      deallocate(c0a_to_c0b)
   end do

   ! Communications setup
   call com_setup(comobs)

   ! Release memory
   deallocate(c0_to_c0a)
else
   ! Send dimensions to ioproc
   call mpl_send(odata%nc0a,mpl%ioproc,mpl%tag)
   call mpl_send(odata%nc0b,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl_send(odata%nc0b,odata%c0b_to_c0,mpl%ioproc,mpl%tag+2)
   call mpl_send(geom%nc0a,odata%c0a_to_c0b,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4

! Communications broadcast
odata%com%prefix = 'o'
call com_bcast(comobs,odata%com)

! Write data
if (mpl%main) then
   lunit = newunit()
   open(unit=lunit,file=trim(nam%datadir)//'/'//trim(nam%prefix)//'_obs_out.dat',status='replace')
   do iobs=1,odata%nobs
      write(lunit,*) odata%lonobs(iobs)*rad2deg,odata%latobs(iobs)*rad2deg
   end do
   close(unit=lunit)

   if (mpl%nproc>1) then
      lunit = newunit()
      open(unit=lunit,file=trim(nam%datadir)//'/'//trim(nam%prefix)//'_obs_comm.dat',status='replace')
      do iproc=1,mpl%nproc
         write(lunit,*) odata%proc_to_nobsa(iproc),comobs(iproc)%nhalo,comobs(iproc)%nexcl
      end do
      close(unit=lunit)
   end if
end if
if (mpl%main) then
   do iproc=1,mpl%nproc
      call com_dealloc(comobs(iproc))
   end do
end if

! Print results
write(mpl%unit,'(a7,a)') '','Number of observations per MPI task:'
do iproc=1,mpl%nproc
   write(mpl%unit,'(a10,a,i3,a,i8)') '','Task ',iproc,': ',odata%proc_to_nobsa(iproc)
end do
write(mpl%unit,'(a7,a,f5.1,a)') '','Observation repartition imbalance: ', &
 & 100.0*float(maxval(odata%proc_to_nobsa)-minval(odata%proc_to_nobsa))/(float(sum(odata%proc_to_nobsa))/float(mpl%nproc)),' %'
write(mpl%unit,'(a7,a)') '','Number of grid points, halo size and number of received values per MPI task:'
if (mpl%main) then
   do iproc=1,mpl%nproc
      write(mpl%unit,'(a10,a,i3,a,i8,a,i8,a,i8)') '','Task ',iproc,': ', &
    & comobs(iproc)%nred,' / ',comobs(iproc)%next,' / ',comobs(iproc)%nhalo
   end do
   nhalo = 0
   do iproc=1,mpl%nproc
      nhalo = nhalo+comobs(iproc)%nhalo
   end do
   write(mpl%unit,'(a7,a,f5.1,a)') '','Ratio between communications and observations: ',100.0*float(nhalo)/float(odata%nobs),' %'
else
   write(mpl%unit,'(a10,a,i3,a,i8,a,i8,a,i8)') '','Task ',mpl%myproc,': ', &
 & odata%com%nred,' / ',odata%com%next,' / ',odata%com%nhalo
end if

! End associate
end associate

end subroutine compute_parameters

end module obsop_parameters
