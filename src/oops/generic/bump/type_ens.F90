!----------------------------------------------------------------------
! Module: type_ens
! Purpose: ensemble derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_ens

use fckit_mpi_module, only: fckit_mpi_sum
use tools_const, only: deg2rad,rad2deg
use tools_func, only: sphere_dist
use tools_kinds, only: kind_real
use type_geom, only: geom_type
use type_io, only: io_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Ensemble derived type
type ens_type
   ! Attributes
   integer :: ne                                  ! Ensemble size
   integer :: nsub                                ! Number of sub-ensembles
   logical :: allocated                           ! Allocation flag

   ! Data
   real(kind_real),allocatable :: fld(:,:,:,:,:)  ! Ensemble perturbation
   real(kind_real),allocatable :: mean(:,:,:,:,:) ! Ensemble mean
contains
   procedure :: set_att => ens_set_att
   procedure :: alloc => ens_alloc
   procedure :: dealloc => ens_dealloc
   procedure :: copy => ens_copy
   procedure :: remove_mean => ens_remove_mean
   procedure :: from => ens_from
   procedure :: apply_bens => ens_apply_bens
   procedure :: cortrack => ens_cortrack
end type ens_type

private
public :: ens_type

contains

!----------------------------------------------------------------------
! Subroutine: ens_set_att
! Purpose: set attributes
!----------------------------------------------------------------------
subroutine ens_set_att(ens,ne,nsub)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble
integer,intent(in) :: ne             ! Ensemble size
integer,intent(in) :: nsub           ! Number of sub-ensembles

! Copy attributes
ens%ne = ne
ens%nsub = nsub
ens%allocated = .false.

end subroutine ens_set_att

!----------------------------------------------------------------------
! Subroutine: ens_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine ens_alloc(ens,nam,geom,ne,nsub)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
integer,intent(in) :: ne             ! Ensemble size
integer,intent(in) :: nsub           ! Number of sub-ensembles

! Copy attributes
call ens%set_att(ne,nsub)

! Allocation
if (ne>0) then
   allocate(ens%fld(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne))
   allocate(ens%mean(geom%nc0a,geom%nl0,nam%nv,nam%nts,nsub))
   ens%allocated = .true.
end if

end subroutine ens_alloc

!----------------------------------------------------------------------
! Subroutine: ens_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine ens_dealloc(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble

! Release memory
if (allocated(ens%fld)) deallocate(ens%fld)
if (allocated(ens%mean)) deallocate(ens%mean)
ens%allocated = .false.

end subroutine ens_dealloc

!----------------------------------------------------------------------
! Subroutine: ens_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine ens_copy(ens_out,ens_in)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens_out ! Output ensemble
type(ens_type),intent(in) :: ens_in      ! Input ensemble

! Copy data
if (allocated(ens_in%fld)) ens_out%fld = ens_in%fld
if (allocated(ens_in%mean)) ens_out%mean = ens_in%mean

end subroutine ens_copy

!----------------------------------------------------------------------
! Subroutine: ens_remove_mean
! Purpose: remove ensemble mean
!----------------------------------------------------------------------
subroutine ens_remove_mean(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble

! Local variables
integer :: isub,ie_sub,ie

if (ens%allocated) then
   ! Loop over sub-ensembles
   do isub=1,ens%nsub
      ! Compute mean
      ens%mean(:,:,:,:,isub) = 0.0
      do ie_sub=1,ens%ne/ens%nsub
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub
         ens%mean(:,:,:,:,isub) = ens%mean(:,:,:,:,isub)+ens%fld(:,:,:,:,ie)
      end do
      ens%mean(:,:,:,:,isub) = ens%mean(:,:,:,:,isub)/(ens%ne/ens%nsub)

      ! Remove mean
      do ie_sub=1,ens%ne/ens%nsub
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub
         ens%fld(:,:,:,:,ie) = ens%fld(:,:,:,:,ie)-ens%mean(:,:,:,:,isub)
      end do
   end do
end if

end subroutine ens_remove_mean

!----------------------------------------------------------------------
! Subroutine: ens_from
! Purpose: copy ensemble array into ensemble data
!----------------------------------------------------------------------
subroutine ens_from(ens,nam,geom,ne,ens_mga)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                        ! Ensemble
type(nam_type),intent(in) :: nam                                            ! Namelist
type(geom_type),intent(in) :: geom                                          ! Geometry
integer,intent(in) :: ne                                                    ! Ensemble size
real(kind_real),intent(in) :: ens_mga(geom%nmga,geom%nl0,nam%nv,nam%nts,ne) ! Ensemble on model grid, halo A

! Local variables
integer :: ie,its,iv,il0

! Allocation
call ens%alloc(nam,geom,ne,1)

if (ens%allocated) then
   ! Copy
   do ie=1,ens%ne
      do its=1,nam%nts
         do iv=1,nam%nv
            do il0=1,geom%nl0
               ens%fld(:,il0,iv,its,ie) = ens_mga(geom%c0a_to_mga,il0,iv,its,ie)
            end do
         end do
      end do
   end do

   ! Remove mean
   call ens%remove_mean
end if

end subroutine ens_from

!----------------------------------------------------------------------
! Subroutine: ens_apply_bens
! Purpose: apply raw ensemble covariance
!----------------------------------------------------------------------
subroutine ens_apply_bens(ens,mpl,nam,geom,fld)

implicit none

! Passed variables
class(ens_type),intent(in) :: ens                                       ! Ensemble
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variable
integer :: ie,ic0a,il0,iv,its
real(kind_real) :: alpha,norm
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: pert(geom%nc0a,geom%nl0,nam%nv,nam%nts)

! Initialization
fld_copy = fld

! Apply localized ensemble covariance formula
fld = 0.0
norm = sqrt(real(nam%ens1_ne-1,kind_real))
do ie=1,nam%ens1_ne
   ! Compute perturbation
   !$omp parallel do schedule(static) private(its,iv,il0,ic0a)
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask_c0a(ic0a,il0)) pert(ic0a,il0,iv,its) = ens%fld(ic0a,il0,iv,its,ie)/norm
            end do
         end do
      end do
   end do
   !$omp end parallel do

   ! Dot product
   call mpl%dot_prod(pert,fld_copy,alpha)

   ! Schur product
   !$omp parallel do schedule(static) private(its,iv,il0,ic0a)
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask_c0a(ic0a,il0)) fld(ic0a,il0,iv,its) = fld(ic0a,il0,iv,its)+alpha*pert(ic0a,il0,iv,its)
            end do
         end do
      end do
   end do
   !$omp end parallel do
end do

end subroutine ens_apply_bens

!----------------------------------------------------------------------
! Subroutine: ens_cortrack
! Purpose: correlation tracker
!----------------------------------------------------------------------
subroutine ens_cortrack(ens,mpl,rng,nam,geom,io)

implicit none

! Passed variables
class(ens_type),intent(in) :: ens   ! Ensemble
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(io_type),intent(in) :: io      ! I/O

! Local variable
integer :: ic0a,ic0,il0,jc0a,jc0,ie,its,ind(2)
integer :: proc_to_ic0a(mpl%nproc),iproc(1),nn_index(1)
real(kind_real) :: dist,proc_to_val(mpl%nproc),lon,lat,val,var_dirac
real(kind_real) :: u(geom%nc0a,geom%nl0,nam%nts),v(geom%nc0a,geom%nl0,nam%nts),ffsq(geom%nc0a,geom%nl0,nam%nts)
real(kind_real) :: var(geom%nc0a,geom%nl0,nam%nts),dirac(geom%nc0a,geom%nl0,nam%nts),cor(geom%nc0a,geom%nl0,nam%nv,nam%nts)
character(len=2) :: timeslotchar
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'ens_cortrack'

! File name
filename = trim(nam%prefix)//'_cortrack'
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)

! Compute variance
write(mpl%info,'(a7,a)') '','Compute variance'
call mpl%flush
var = 0.0
do ie=1,ens%ne
   var = var+ens%fld(:,:,1,:,ie)**2
end do
var = var/real(ens%ne-ens%nsub,kind_real)

if (nam%nv==3) then
   ! Compute wind speed squared (variables 2 and 3 should be u and v)
   write(mpl%info,'(a7,a)') '','Compute wind speed'
   call mpl%flush
   u = sum(ens%mean(:,:,2,:,:),dim=4)/real(ens%nsub,kind_real)
   v = sum(ens%mean(:,:,3,:,:),dim=4)/real(ens%nsub,kind_real)
   ffsq = u**2+v**2

   ! Write wind
   write(mpl%info,'(a7,a)') '','Write wind'
   call mpl%flush
   do its=1,nam%nts
      write(timeslotchar,'(i2.2)') nam%timeslot(its)
      call io%fld_write(mpl,nam,geom,filename,'u_'//timeslotchar,u(:,:,its))
      call io%fld_write(mpl,nam,geom,filename,'v_'//timeslotchar,v(:,:,its))
   end do

   if (.true.) then
      ! Find local maximum value and index
      write(mpl%info,'(a7,a)') '','Dirac point based on maximum wind speed'
      call mpl%flush
      val = maxval(ffsq(:,:,1))
      ind = maxloc(ffsq(:,:,1))
      ic0a = ind(1)
      il0 = ind(2)
      call mpl%f_comm%allgather(val,proc_to_val)
      call mpl%f_comm%allgather(ic0a,proc_to_ic0a)
      iproc = maxloc(proc_to_val)
   else
      ! Find local minimum value and index
      write(mpl%info,'(a7,a)') '','Dirac point based on minimum wind speed'
      call mpl%flush
      val = minval(ffsq(:,:,1))
      ind = minloc(ffsq(:,:,1))
      ic0a = ind(1)
      il0 = ind(2)
      call mpl%f_comm%allgather(val,proc_to_val)
      call mpl%f_comm%allgather(ic0a,proc_to_ic0a)
      iproc = minloc(proc_to_val)
   end if
   call mpl%f_comm%broadcast(ic0a,iproc(1)-1)
   call mpl%f_comm%broadcast(il0,iproc(1)-1)
else
   if (.false.) then
      ! Define lon/lat at the domain center
      lon = 0.5*(minval(geom%lon)+maxval(geom%lon))
      lat = 0.5*(minval(geom%lat)+maxval(geom%lat))

      ! Find nearest neighbor
      call geom%tree%find_nearest_neighbors(lon,lat,1,nn_index)
      ic0 = nn_index(1)
   else
      ! Random point
      if (mpl%main) call rng%rand_integer(1,geom%nc0,ic0)
      call mpl%f_comm%broadcast(ic0,mpl%ioproc-1)
   end if
   if (mpl%main) call rng%rand_integer(1,geom%nl0,il0)
   call mpl%f_comm%broadcast(il0,mpl%ioproc-1)

   ! Broadcast dirac point
   ic0a = mpl%msv%vali
   iproc(1) = geom%c0_to_proc(ic0)
   if (mpl%myproc==iproc(1)) ic0a = geom%c0_to_c0a(ic0)
   call mpl%f_comm%broadcast(ic0a,iproc(1)-1)
   lon = geom%lon(ic0)
   lat = geom%lat(ic0)
   write(mpl%info,'(a7,a,i2,a,e10.3,a,e10.3)') '','Timeslot ',nam%timeslot(1),' ~> lon / lat:       ', &
 & lon*rad2deg,' / ',lat*rad2deg
   call mpl%flush
end if

! Generate dirac field
dirac = 0.0
if (mpl%myproc==iproc(1)) dirac(ic0a,il0,1) = 1.0

! Apply raw ensemble covariance
cor = 0.0
cor(:,:,1,:) = dirac
call ens%apply_bens(mpl,nam,geom,cor)

! Normalize
call mpl%f_comm%allreduce(sum(dirac*var),var_dirac,fckit_mpi_sum())
cor(:,:,1,:) = cor(:,:,1,:)/sqrt(var*var_dirac)

! Write correlation
do its=1,nam%nts
   write(timeslotchar,'(i2.2)') nam%timeslot(its)
   call io%fld_write(mpl,nam,geom,filename,'cor_'//timeslotchar,cor(:,:,1,its))
end do

! Correlation tracker
write(timeslotchar,'(i2.2)') nam%timeslot(1)
call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_0',cor(:,:,1,1))
call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_1',cor(:,:,1,2))
do its=2,nam%nts-1
   ! Locate correlation maximum and index
   val = 0.0
   ic0a = mpl%msv%vali
   do jc0a=1,geom%nc0a
      jc0 = geom%c0a_to_c0(jc0a)
      call sphere_dist(lon,lat,geom%lon(jc0),geom%lat(jc0),dist)
      if (dist<nam%adv_rad) then
         if (cor(jc0a,il0,1,its)>val) then
            val = cor(jc0a,il0,1,its)
            ic0a = jc0a
         end if
      end if
   end do
   if (mpl%msv%isi(ic0a)) call mpl%abort(subr,'cannot locate correlation maximum and index')
   call mpl%f_comm%allgather(val,proc_to_val)
   call mpl%f_comm%allgather(ic0a,proc_to_ic0a)
   iproc = maxloc(proc_to_val)
   if (mpl%myproc==iproc(1)) then
      ic0a = proc_to_ic0a(iproc(1))
      ic0 = geom%c0a_to_c0(ic0a)
      val = cor(ic0a,il0,1,its)
   end if
   call mpl%f_comm%broadcast(ic0,iproc(1)-1)
   lon = geom%lon(ic0)
   lat = geom%lat(ic0)
   call mpl%f_comm%broadcast(val,iproc(1)-1)
   write(mpl%info,'(a7,a,i2,a,e10.3,a,e10.3,a,f6.1)') '','Timeslot ',nam%timeslot(its),' ~> lon / lat / val: ', &
 & lon*rad2deg,' / ',lat*rad2deg,' / ',val
   call mpl%flush

   ! Generate dirac field
   dirac = 0.0
   if (mpl%myproc==iproc(1)) dirac(ic0a,il0,its) = 1.0

   ! Apply raw ensemble covariance
   cor = 0.0
   cor(:,:,1,:) = dirac
   call ens%apply_bens(mpl,nam,geom,cor)

   ! Normalize
   call mpl%f_comm%allreduce(sum(dirac*var),var_dirac,fckit_mpi_sum())
   cor(:,:,1,:) = cor(:,:,1,:)/sqrt(var*var_dirac)

   ! Write correlation tracker
   write(timeslotchar,'(i2.2)') nam%timeslot(its)
   call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_0',cor(:,:,1,its))
   call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_1',cor(:,:,1,its+1))
end do

end subroutine ens_cortrack

end module type_ens
