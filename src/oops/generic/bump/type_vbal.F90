!----------------------------------------------------------------------
! Module: type_vbal
! Purpose: vertical balance derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_vbal

use fckit_mpi_module, only: fckit_mpi_sum
!$ use omp_lib
use netcdf
use tools_func, only: syminv
use tools_kinds, only: kind_real,nc_kind_real
use tools_repro, only: infeq
use type_bpar, only: bpar_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type
use type_samp, only: samp_type
use type_vbal_blk, only: vbal_blk_type

implicit none

! Vertical balance derived type
type vbal_type
   type(samp_type) :: samp                     ! Sampling
   integer :: np                               ! Maximum number of neighbors
   integer :: nc2b                             ! Subset Sc2 size, halo B
   logical :: allocated                        ! Allocation flag
   integer,allocatable :: h_n_s(:,:)           ! Number of neighbors for the horizontal interpolation
   integer,allocatable :: h_c2b(:,:,:)         ! Index of neighbors for the horizontal interpolation
   real(kind_real),allocatable :: h_S(:,:,:)   ! Weight of neighbors for the horizontal interpolation
   type(vbal_blk_type),allocatable :: blk(:,:) ! Vertical balance blocks
contains
   procedure :: alloc => vbal_alloc
   procedure :: partial_dealloc => vbal_partial_dealloc
   procedure :: dealloc => vbal_dealloc
   procedure :: read => vbal_read
   procedure :: write => vbal_write
   procedure :: run_vbal => vbal_run_vbal
   procedure :: run_vbal_tests => vbal_run_vbal_tests
   procedure :: apply => vbal_apply
   procedure :: apply_inv => vbal_apply_inv
   procedure :: apply_ad => vbal_apply_ad
   procedure :: apply_inv_ad => vbal_apply_inv_ad
   procedure :: test_inverse => vbal_test_inverse
   procedure :: test_adjoint => vbal_test_adjoint
end type vbal_type

logical,parameter :: diag_auto = .true.   ! Diagonal auto-covariance for the inversion
real(kind_real),parameter :: var_th = 0.8 ! Variance threshold to truncate the vertical auto-covariance spectrum

private
public :: vbal_type

contains

!----------------------------------------------------------------------
! Subroutine: vbal_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine vbal_alloc(vbal,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(vbal_type),intent(inout) :: vbal ! Vertical balance
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters

! Local variables
integer :: iv,jv
character(len=1024),parameter :: subr = 'vbal_alloc'

! Find number of neighbors
vbal%np = 0
if (trim(nam%diag_interp)=='bilin') then
   ! Bilinear interpolation
   vbal%np = 3
elseif (trim(nam%diag_interp)=='natural') then
   ! Natural neighbors
   vbal%np = 40
else
   call mpl%abort(subr,'wrong interpolation type')
end if

! Allocation
allocate(vbal%h_n_s(geom%nc0a,geom%nl0i))
allocate(vbal%h_c2b(vbal%np,geom%nc0a,geom%nl0i))
allocate(vbal%h_S(vbal%np,geom%nc0a,geom%nl0i))
allocate(vbal%blk(nam%nv,nam%nv))
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         call vbal%blk(iv,jv)%alloc(nam,geom,vbal%samp%nc2b,iv,jv)
      end if
   end do
end do

! Update allocation flag
vbal%allocated = .true.

end subroutine vbal_alloc

!----------------------------------------------------------------------
! Subroutine: vbal_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine vbal_partial_dealloc(vbal)

implicit none

! Passed variables
class(vbal_type),intent(inout) :: vbal ! Vertical balance

! Local variables
integer :: iv,jv

! Release memory
if (allocated(vbal%blk)) then
   do jv=1,size(vbal%blk,2)
      do iv=1,size(vbal%blk,1)
         call vbal%blk(iv,jv)%partial_dealloc
      end do
   end do
end if

end subroutine vbal_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: vbal_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine vbal_dealloc(vbal)

implicit none

! Passed variables
class(vbal_type),intent(inout) :: vbal ! Vertical balance

! Local variables
integer :: iv,jv

! Release memory
if (allocated(vbal%h_n_s)) deallocate(vbal%h_n_s)
if (allocated(vbal%h_c2b)) deallocate(vbal%h_c2b)
if (allocated(vbal%h_S)) deallocate(vbal%h_S)
if (allocated(vbal%blk)) then
   do jv=1,size(vbal%blk,2)
      do iv=1,size(vbal%blk,1)
         call vbal%blk(iv,jv)%dealloc
      end do
   end do
   deallocate(vbal%blk)
end if

! Update allocation flag
vbal%allocated = .false.

end subroutine vbal_dealloc

!----------------------------------------------------------------------
! Subroutine: vbal_read
! Purpose: read
!----------------------------------------------------------------------
subroutine vbal_read(vbal,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(vbal_type),intent(inout) :: vbal ! Vertical balance
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters

! Local variables
integer :: iv,jv
integer :: ncid,np_id,nc0a_id,nc2b_id,nl0i_id,nl0_1_id,nl0_2_id,h_n_s_id,h_c2b_id,h_S_id,reg_id(nam%nv,nam%nv)
integer :: nc0a_test,nl0i_test,nl0_1_test,nl0_2_test
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'vbal_read'

! Open file
write(filename,'(a,a,i4.4,a,i4.4)') trim(nam%prefix),'_vbal_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

! Get dimensions
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'np',np_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,np_id,len=vbal%np))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc0a',nc0a_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc0a_id,len=nc0a_test))
if (nc0a_test/=geom%nc0a) call mpl%abort(subr,'wrong dimension when reading vbal')
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc2b',nc2b_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc2b_id,len=vbal%samp%nc2b))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nl0i',nl0i_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nl0i_id,len=nl0i_test))
if (nl0i_test/=geom%nl0i) call mpl%abort(subr,'wrong dimension when reading vbal')
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nl0_1',nl0_1_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nl0_1_id,len=nl0_1_test))
if (nl0_1_test/=geom%nl0) call mpl%abort(subr,'wrong dimension when reading vbal')
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nl0_2',nl0_2_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nl0_2_id,len=nl0_2_test))
if (nl0_2_test/=geom%nl0) call mpl%abort(subr,'wrong dimension when reading vbal')

! Allocation
allocate(vbal%h_n_s(geom%nc0a,geom%nl0i))
allocate(vbal%h_c2b(vbal%np,geom%nc0a,geom%nl0i))
allocate(vbal%h_S(vbal%np,geom%nc0a,geom%nl0i))
allocate(vbal%blk(nam%nv,nam%nv))
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         call vbal%blk(iv,jv)%alloc(nam,geom,vbal%samp%nc2b,iv,jv)
      end if
   end do
end do

! Get variable id
call mpl%ncerr(subr,nf90_inq_varid(ncid,'h_n_s',h_n_s_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'h_c2b',h_c2b_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'h_S',h_S_id))
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(vbal%blk(iv,jv)%name)//'_reg',reg_id(iv,jv)))
   end do
end do

! Read data
call mpl%ncerr(subr,nf90_get_var(ncid,h_n_s_id,vbal%h_n_s))
call mpl%ncerr(subr,nf90_get_var(ncid,h_c2b_id,vbal%h_c2b))
call mpl%ncerr(subr,nf90_get_var(ncid,h_S_id,vbal%h_S))
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) call mpl%ncerr(subr,nf90_get_var(ncid,reg_id(iv,jv),vbal%blk(iv,jv)%reg))
   end do
end do

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

end subroutine vbal_read

!----------------------------------------------------------------------
! Subroutine: vbal_write
! Purpose: write
!----------------------------------------------------------------------
subroutine vbal_write(vbal,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(vbal_type),intent(inout) :: vbal ! Vertical balance
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters

! Local variables
integer :: iv,jv
integer :: ncid,np_id,nc0a_id,nc2b_id,nl0i_id,nl0_1_id,nl0_2_id,h_n_s_id,h_c2b_id,h_S_id
integer :: reg_id(nam%nv,nam%nv),auto_id(nam%nv,nam%nv),cross_id(nam%nv,nam%nv),auto_inv_id(nam%nv,nam%nv)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'vbal_write'

! Create file
write(filename,'(a,a,i4.4,a,i4.4)') trim(nam%prefix),'_vbal_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%write(mpl,ncid)

! Define dimensions
call mpl%ncerr(subr,nf90_def_dim(ncid,'np',vbal%np,np_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0a',geom%nc0a,nc0a_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2b',vbal%samp%nc2b,nc2b_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0i',geom%nl0i,nl0i_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0_1',geom%nl0,nl0_1_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0_2',geom%nl0,nl0_2_id))

! Define variables
call mpl%ncerr(subr,nf90_def_var(ncid,'h_n_s',nf90_int,(/nc0a_id,nl0i_id/),h_n_s_id))
call mpl%ncerr(subr,nf90_put_att(ncid,h_n_s_id,'_FillValue',mpl%msv%vali))
call mpl%ncerr(subr,nf90_def_var(ncid,'h_c2b',nf90_int,(/np_id,nc0a_id,nl0i_id/),h_c2b_id))
call mpl%ncerr(subr,nf90_put_att(ncid,h_c2b_id,'_FillValue',mpl%msv%vali))
call mpl%ncerr(subr,nf90_def_var(ncid,'h_S',nc_kind_real,(/np_id,nc0a_id,nl0i_id/),h_S_id))
call mpl%ncerr(subr,nf90_put_att(ncid,h_S_id,'_FillValue',mpl%msv%valr))
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(vbal%blk(iv,jv)%name)//'_auto',nc_kind_real,(/nc2b_id,nl0_1_id,nl0_2_id/), &
       & auto_id(iv,jv)))
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(vbal%blk(iv,jv)%name)//'_cross',nc_kind_real,(/nc2b_id,nl0_1_id,nl0_2_id/), &
       & cross_id(iv,jv)))
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(vbal%blk(iv,jv)%name)//'_auto_inv',nc_kind_real,(/nc2b_id,nl0_1_id,nl0_2_id/), &
       & auto_inv_id(iv,jv)))
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(vbal%blk(iv,jv)%name)//'_reg',nc_kind_real,(/nc2b_id,nl0_1_id,nl0_2_id/), &
       & reg_id(iv,jv)))
      end if
   end do
end do

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write variables
call mpl%ncerr(subr,nf90_put_var(ncid,h_n_s_id,vbal%h_n_s))
call mpl%ncerr(subr,nf90_put_var(ncid,h_c2b_id,vbal%h_c2b))
call mpl%ncerr(subr,nf90_put_var(ncid,h_S_id,vbal%h_S))
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         call mpl%ncerr(subr,nf90_put_var(ncid,reg_id(iv,jv),vbal%blk(iv,jv)%reg))
         call mpl%ncerr(subr,nf90_put_var(ncid,auto_id(iv,jv),vbal%blk(iv,jv)%auto))
         call mpl%ncerr(subr,nf90_put_var(ncid,cross_id(iv,jv),vbal%blk(iv,jv)%cross))
         call mpl%ncerr(subr,nf90_put_var(ncid,auto_inv_id(iv,jv),vbal%blk(iv,jv)%auto_inv))
      end if
   end do
end do

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

end subroutine vbal_write

!----------------------------------------------------------------------
! Subroutine: vbal_run_vbal
! Purpose: compute vertical balance
!----------------------------------------------------------------------
subroutine vbal_run_vbal(vbal,mpl,rng,nam,geom,bpar,io,ens,ensu)

implicit none

! Passed variables
class(vbal_type),intent(inout) :: vbal ! Vertical balance
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O
type(ens_type), intent(in) :: ens      ! Ensemble
type(ens_type),intent(inout) :: ensu   ! Unbalanced ensemble

! Local variables
integer :: il0i,i_s,ic0a,ic2b,ic2,ie,ie_sub,ic0,jl0,il0,isub,ic1,ic1a,iv,jv,offset,nc1a
real(kind_real) :: fld(geom%nc0a,geom%nl0)
real(kind_real) :: auto_avg(nam%nc2,geom%nl0,geom%nl0),cross_avg(nam%nc2,geom%nl0,geom%nl0),auto_inv(geom%nl0,geom%nl0)
real(kind_real) :: sbuf(2*nam%nc2*geom%nl0**2),rbuf(2*nam%nc2*geom%nl0**2)
real(kind_real),allocatable :: list_auto(:),list_cross(:),auto_avg_tmp(:,:)
real(kind_real),allocatable :: fld_1(:,:),fld_2(:,:),auto(:,:,:,:),cross(:,:,:,:)
logical :: valid,mask_unpack(geom%nl0,geom%nl0)

! Setup sampling
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'
call mpl%flush
call vbal%samp%setup_sampling(mpl,rng,nam,geom,bpar,io,ens)

! Compute MPI distribution, halo A
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute MPI distribution, halos A'
call mpl%flush
call vbal%samp%compute_mpi_a(mpl,nam,geom)

! Compute MPI distribution, halos A-B
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute MPI distribution, halos A-B'
call mpl%flush
call vbal%samp%compute_mpi_ab(mpl,nam,geom)

! Allocation
call ensu%alloc(nam,geom,ens%ne,ens%nsub)

! Copy ensemble
call ensu%copy(ens)

! Allocation
allocate(fld_1(vbal%samp%nc1a,geom%nl0))
allocate(fld_2(vbal%samp%nc1a,geom%nl0))
allocate(auto(vbal%samp%nc1a,geom%nl0,geom%nl0,ens%nsub))
allocate(cross(vbal%samp%nc1a,geom%nl0,geom%nl0,ens%nsub))
if (.not.diag_auto) allocate(auto_avg_tmp(geom%nl0,geom%nl0))
call vbal%alloc(mpl,nam,geom,bpar)

! Initialization
mask_unpack = .true.
vbal%h_n_s = 0
vbal%h_c2b = mpl%msv%vali
vbal%h_S = mpl%msv%valr

! Get interpolation coefficients
do il0i=1,geom%nl0i
   do i_s=1,vbal%samp%h(il0i)%n_s
      ic0a = vbal%samp%h(il0i)%row(i_s)
      vbal%h_n_s(ic0a,il0i) = vbal%h_n_s(ic0a,il0i)+1
      vbal%h_c2b(vbal%h_n_s(ic0a,il0i),ic0a,il0i) = vbal%samp%h(il0i)%col(i_s)
      vbal%h_S(vbal%h_n_s(ic0a,il0i),ic0a,il0i) = vbal%samp%h(il0i)%S(i_s)
   end do
end do

do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         ! Initialization
         write(mpl%info,'(a7,a)') '','Unbalancing: '//trim(nam%varname(iv))//' with respect to unbalanced '//trim(nam%varname(jv))
         call mpl%flush
         auto = 0.0
         cross = 0.0

         ! Loop on sub-ensembles
         do isub=1,ensu%nsub
            if (ensu%nsub==1) then
               write(mpl%info,'(a10,a)') '','Full ensemble, member:'
               call mpl%flush(.false.)
            else
               write(mpl%info,'(a10,a,i4,a)') '','Sub-ensemble ',isub,', member:'
               call mpl%flush(.false.)
            end if

            ! Compute centered moments iteratively
            do ie_sub=1,ensu%ne/ensu%nsub
               write(mpl%info,'(i4)') ie_sub
               call mpl%flush(.false.)

               ! Full ensemble index
               ie = ie_sub+(isub-1)*ensu%ne/ensu%nsub

               ! Copy all separations points
               !$omp parallel do schedule(static) private(il0,ic1a,ic1,ic0,ic0a)
               do il0=1,geom%nl0
                  do ic1a=1,vbal%samp%nc1a
                     ! Indice
                     ic1 = vbal%samp%c1a_to_c1(ic1a)

                     if (vbal%samp%c1l0_log(ic1,il0)) then
                        ! Indice
                        ic0 = vbal%samp%c1_to_c0(ic1)
                        ic0a = geom%c0_to_c0a(ic0)

                        ! Copy points
                        fld_1(ic1a,il0) = ensu%fld(ic0a,il0,iv,1,ie)
                        fld_2(ic1a,il0) = ensu%fld(ic0a,il0,jv,1,ie)
                     end if
                  end do
               end do
               !$omp end parallel do

               !$omp parallel do schedule(static) private(il0,jl0,ic1a,ic1)
               do il0=1,geom%nl0
                  do jl0=1,geom%nl0
                     do ic1a=1,vbal%samp%nc1a
                        ! Indice
                        ic1 = vbal%samp%c1a_to_c1(ic1a)

                        ! Covariance
                        if (vbal%samp%c1l0_log(ic1,il0).and.vbal%samp%c1l0_log(ic1,jl0)) then
                           auto(ic1a,jl0,il0,isub) = auto(ic1a,jl0,il0,isub)+fld_2(ic1a,il0)*fld_2(ic1a,jl0)
                           cross(ic1a,jl0,il0,isub) = cross(ic1a,jl0,il0,isub)+fld_1(ic1a,il0)*fld_2(ic1a,jl0)
                        end if
                     end do
                  end do
               end do
               !$omp end parallel do
            end do
            write(mpl%info,'(a)') ''
            call mpl%flush
         end do

         ! Average covariances
         write(mpl%info,'(a10,a)') '','Average covariances: '
         call mpl%flush(.false.)
         call mpl%prog_init(nam%nc2)
         do ic2=1,nam%nc2
            !$omp parallel do schedule(static) private(il0,jl0,nc1a,ic1a,ic1,valid,isub), &
            !$omp&                             firstprivate(list_auto,list_cross)
            do il0=1,geom%nl0
               do jl0=1,geom%nl0
                  ! Allocation
                  allocate(list_auto(vbal%samp%nc1a))
                  allocate(list_cross(vbal%samp%nc1a))

                  ! Fill lists
                  nc1a = 0
                  do ic1a=1,vbal%samp%nc1a
                     ! Index
                     ic1 = vbal%samp%c1a_to_c1(ic1a)

                     ! Check validity
                     valid = vbal%samp%c1l0_log(ic1,il0).and.vbal%samp%c1l0_log(ic1,jl0).and.vbal%samp%vbal_mask(ic1a,ic2)

                     if (valid) then
                        ! Update
                         nc1a = nc1a+1

                        ! Averages for diagnostics
                        list_auto(nc1a) = sum(auto(ic1a,jl0,il0,:))/real(ensu%nsub,kind_real)
                        list_cross(nc1a) = sum(cross(ic1a,jl0,il0,:))/real(ensu%nsub,kind_real)
                     end if
                  end do

                  ! Average
                  if (nc1a>0) then
                     auto_avg(ic2,jl0,il0) = sum(list_auto(1:nc1a))
                     cross_avg(ic2,jl0,il0) = sum(list_cross(1:nc1a))
                  else
                     auto_avg(ic2,jl0,il0) = 0.0
                     cross_avg(ic2,jl0,il0) = 0.0
                  end if

                  ! Release memory
                  deallocate(list_auto)
                  deallocate(list_cross)
               end do
             end do
            !$omp end parallel do

            ! Update
            call mpl%prog_print(ic2)
         end do
         call mpl%prog_final

         ! Gather data
         write(mpl%info,'(a10,a)') '','Gather data'
         call mpl%flush
         if (mpl%nproc>1) then
            ! Pack data
            offset = 0
            do ic2=1,nam%nc2
               sbuf(offset+1:offset+geom%nl0**2) = pack(auto_avg(ic2,:,:),.true.)
               offset = offset+geom%nl0**2
               sbuf(offset+1:offset+geom%nl0**2) = pack(cross_avg(ic2,:,:),.true.)
               offset = offset+geom%nl0**2
            end do

            ! Reduce data
            call mpl%f_comm%allreduce(sbuf,rbuf,fckit_mpi_sum())

            ! Unpack data
            offset = 0
            do ic2=1,nam%nc2
               auto_avg(ic2,:,:) = unpack(rbuf(offset+1:offset+geom%nl0**2),mask_unpack,auto_avg(ic2,:,:))
               offset = offset+geom%nl0**2
               cross_avg(ic2,:,:) = unpack(rbuf(offset+1:offset+geom%nl0**2),mask_unpack,cross_avg(ic2,:,:))
               offset = offset+geom%nl0**2
            end do
         end if

         ! Compute regressions
         write(mpl%info,'(a10,a)') '','Compute regressions: '
         call mpl%flush(.false.)
         call mpl%prog_init(vbal%samp%nc2b)
         do ic2b=1,vbal%samp%nc2b
            ! Global index
            ic2 = vbal%samp%c2b_to_c2(ic2b)

            if (diag_auto) then
               ! Diagonal inversion
               auto_inv = 0.0
               do il0=1,geom%nl0
                  if (auto_avg(ic2,il0,il0)>0.0) auto_inv(il0,il0) = 1.0/auto_avg(ic2,il0,il0)
               end do
            else
               ! Inverse the vertical auto-covariance
               auto_avg_tmp = auto_avg(ic2,:,:)
               call syminv(mpl,geom%nl0,auto_avg_tmp,auto_inv)
            end if

            ! Compute the regression
            vbal%blk(iv,jv)%auto(ic2b,:,:) = auto_avg(ic2,:,:)
            vbal%blk(iv,jv)%cross(ic2b,:,:) = cross_avg(ic2,:,:)
            vbal%blk(iv,jv)%auto_inv(ic2b,:,:) = auto_inv
            vbal%blk(iv,jv)%reg(ic2b,:,:) = matmul(cross_avg(ic2,:,:),auto_inv)

            ! Update
            call mpl%prog_print(ic2b)
         end do
         call mpl%prog_final
      end if
   end do

   ! Unbalance ensemble
   if (any(bpar%vbal_block(iv,1:iv-1))) then
      write(mpl%info,'(a10,a)') '','Unbalance ensemble members: '
      call mpl%flush(.false.)
      do ie=1,ensu%ne
         write(mpl%info,'(i4)') ie
         call mpl%flush(.false.)
         do jv=1,iv-1
            if (bpar%vbal_block(iv,jv)) then
               fld = ensu%fld(:,:,jv,1,ie)
               call vbal%blk(iv,jv)%apply(geom,vbal%np,vbal%h_n_s,vbal%h_c2b,vbal%h_S,fld)
               ensu%fld(:,:,iv,1,ie) = ensu%fld(:,:,iv,1,ie)-fld
            end if
         end do
      end do
      write(mpl%info,'(a)') ''
      call mpl%flush
   end if
end do

! Write balance operator
if (nam%write_vbal) call vbal%write(mpl,nam,geom,bpar)

! Release memory
deallocate(fld_1)
deallocate(fld_2)
deallocate(auto)
deallocate(cross)
if (.not.diag_auto) deallocate(auto_avg_tmp)

end subroutine vbal_run_vbal

!----------------------------------------------------------------------
! Subroutine: vbal_run_vbal_tests
! Purpose: compute vertical balance tests
!----------------------------------------------------------------------
subroutine vbal_run_vbal_tests(vbal,mpl,rng,nam,geom,bpar)

implicit none

! Passed variables
class(vbal_type),intent(inout) :: vbal ! Vertical balance
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters

! Test inverse
call vbal%test_inverse(mpl,rng,nam,geom,bpar)

! Test adjoint
call vbal%test_adjoint(mpl,rng,nam,geom,bpar)

end subroutine vbal_run_vbal_tests

!----------------------------------------------------------------------
! Subroutine: vbal_apply
! Purpose: apply vertical balance
!----------------------------------------------------------------------
subroutine vbal_apply(vbal,nam,geom,bpar,fld)

implicit none

! Passed variables
class(vbal_type),intent(in) :: vbal                             ! Vertical balance
type(nam_type),intent(in) :: nam                                ! Namelist
type(geom_type),intent(in) :: geom                              ! Geometry
type(bpar_type),intent(in) :: bpar                              ! Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv) ! Source/destination vector

! Local variables
integer :: iv,jv
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0),fld_out(geom%nc0a,geom%nl0,nam%nv)

! Initialization
fld_out = fld

! Add balance component
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         fld_tmp = fld(:,:,jv)
         call vbal%blk(iv,jv)%apply(geom,vbal%np,vbal%h_n_s,vbal%h_c2b,vbal%h_S,fld_tmp)
         fld_out(:,:,iv) = fld_out(:,:,iv)+fld_tmp
      end if
   end do
end do

! Final copy
fld = fld_out

end subroutine vbal_apply

!----------------------------------------------------------------------
! Subroutine: vbal_apply_inv
! Purpose: apply inverse vertical balance
!----------------------------------------------------------------------
subroutine vbal_apply_inv(vbal,nam,geom,bpar,fld)

implicit none

! Passed variables
class(vbal_type),intent(in) :: vbal                             ! Vertical balance
type(nam_type),intent(in) :: nam                                ! Namelist
type(geom_type),intent(in) :: geom                              ! Geometry
type(bpar_type),intent(in) :: bpar                              ! Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv) ! Source/destination vector

! Local variables
integer :: iv,jv
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0),fld_out(geom%nc0a,geom%nl0,nam%nv)

! Initialization
fld_out = fld

! Remove balance component
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         fld_tmp = fld_out(:,:,jv)
         call vbal%blk(iv,jv)%apply(geom,vbal%np,vbal%h_n_s,vbal%h_c2b,vbal%h_S,fld_tmp)
         fld_out(:,:,iv) = fld_out(:,:,iv)-fld_tmp
      end if
   end do
end do

! Final copy
fld = fld_out

end subroutine vbal_apply_inv

!----------------------------------------------------------------------
! Subroutine: vbal_apply_ad
! Purpose: apply adjoint vertical balance
!----------------------------------------------------------------------
subroutine vbal_apply_ad(vbal,nam,geom,bpar,fld)

implicit none

! Passed variables
class(vbal_type),intent(in) :: vbal                             ! Vertical balance
type(nam_type),intent(in) :: nam                                ! Namelist
type(geom_type),intent(in) :: geom                              ! Geometry
type(bpar_type),intent(in) :: bpar                              ! Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv) ! Source/destination vector

! Local variables
integer :: iv,jv
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0),fld_out(geom%nc0a,geom%nl0,nam%nv)

! Initialization
fld_out = fld

! Add balance component
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         fld_tmp = fld(:,:,iv)
         call vbal%blk(iv,jv)%apply_ad(geom,vbal%np,vbal%h_n_s,vbal%h_c2b,vbal%h_S,fld_tmp)
         fld_out(:,:,jv) = fld_out(:,:,jv)+fld_tmp
      end if
   end do
end do

! Final copy
fld = fld_out

end subroutine vbal_apply_ad

!----------------------------------------------------------------------
! Subroutine: vbal_apply_inv_ad
! Purpose: apply inverse adjoint vertical balance
!----------------------------------------------------------------------
subroutine vbal_apply_inv_ad(vbal,nam,geom,bpar,fld)

implicit none

! Passed variables
class(vbal_type),intent(in) :: vbal                             ! Vertical balance
type(nam_type),intent(in) :: nam                                ! Namelist
type(geom_type),intent(in) :: geom                              ! Geometry
type(bpar_type),intent(in) :: bpar                              ! Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv) ! Source/destination vector

! Local variables
integer :: iv,jv
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0),fld_out(geom%nc0a,geom%nl0,nam%nv)

! Initialization
fld_out = fld

! Remove balance component
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         fld_tmp = fld_out(:,:,iv)
         call vbal%blk(iv,jv)%apply_ad(geom,vbal%np,vbal%h_n_s,vbal%h_c2b,vbal%h_S,fld_tmp)
         fld_out(:,:,jv) = fld_out(:,:,jv)-fld_tmp
      end if
   end do
end do

! Final copy
fld = fld_out

end subroutine vbal_apply_inv_ad

!----------------------------------------------------------------------
! Subroutine: vbal_test_inverse
! Purpose: test vertical balance inverse
!----------------------------------------------------------------------
subroutine vbal_test_inverse(vbal,mpl,rng,nam,geom,bpar)

implicit none

! Passed variables
class(vbal_type),intent(in) :: vbal ! Vertical balance
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(bpar_type),intent(in) :: bpar  ! Block parameters

! Local variables
real(kind_real) :: mse,mse_tot
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv),fld_save(geom%nc0a,geom%nl0,nam%nv)

! Generate random field
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Direct / inverse
fld = fld_save
call vbal%apply(nam,geom,bpar,fld)
call vbal%apply_inv(nam,geom,bpar,fld)
mse = sum((fld-fld_save)**2)
call mpl%f_comm%allreduce(mse,mse_tot,fckit_mpi_sum())
write(mpl%info,'(a7,a,e15.8)') '','Vertical balance direct/inverse test:  ',mse_tot
call mpl%flush

! Inverse / direct
fld = fld_save
call vbal%apply_inv(nam,geom,bpar,fld)
call vbal%apply(nam,geom,bpar,fld)
mse = sum((fld-fld_save)**2)
call mpl%f_comm%allreduce(mse,mse_tot,fckit_mpi_sum())
write(mpl%info,'(a7,a,e15.8)') '','Vertical balance inverse/direct test:  ',mse_tot
call mpl%flush

! Direct / inverse, adjoint
fld = fld_save
call vbal%apply_ad(nam,geom,bpar,fld)
call vbal%apply_inv_ad(nam,geom,bpar,fld)
mse = sum((fld-fld_save)**2)
call mpl%f_comm%allreduce(mse,mse_tot,fckit_mpi_sum())
write(mpl%info,'(a7,a,e15.8)') '','Vertical balance direct/inverse (adjoint) test:  ',mse_tot
call mpl%flush

! Inverse / direct
fld = fld_save
call vbal%apply_inv_ad(nam,geom,bpar,fld)
call vbal%apply_ad(nam,geom,bpar,fld)
mse = sum((fld-fld_save)**2)
call mpl%f_comm%allreduce(mse,mse_tot,fckit_mpi_sum())
write(mpl%info,'(a7,a,e15.8)') '','Vertical balance inverse/direct (adjoint) test:  ',mse_tot
call mpl%flush

end subroutine vbal_test_inverse

!----------------------------------------------------------------------
! Subroutine: vbal_test_adjoint
! Purpose: test vertical balance adjoint
!----------------------------------------------------------------------
subroutine vbal_test_adjoint(vbal,mpl,rng,nam,geom,bpar)

implicit none

! Passed variables
class(vbal_type),intent(in) :: vbal ! Vertical balance
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(bpar_type),intent(in) :: bpar  ! Block parameters

! Local variables
integer :: iv,jv
real(kind_real) :: sum1,sum2
real(kind_real) :: fld1_blk(geom%nc0a,geom%nl0,nam%nv),fld1_dir(geom%nc0a,geom%nl0,nam%nv)
real(kind_real) :: fld1_inv(geom%nc0a,geom%nl0,nam%nv),fld1_save(geom%nc0a,geom%nl0,nam%nv)
real(kind_real) :: fld2_blk(geom%nc0a,geom%nl0,nam%nv),fld2_dir(geom%nc0a,geom%nl0,nam%nv)
real(kind_real) :: fld2_inv(geom%nc0a,geom%nl0,nam%nv),fld2_save(geom%nc0a,geom%nl0,nam%nv)

! Generate random field
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld2_save)

! Block adjoint test
fld1_blk = fld1_save
fld2_blk = fld2_save
do iv=1,nam%nv
   do jv=1,nam%nv
      if (bpar%vbal_block(iv,jv)) then
         call vbal%blk(iv,jv)%apply(geom,vbal%np,vbal%h_n_s,vbal%h_c2b,vbal%h_S,fld1_blk(:,:,iv))
         call vbal%blk(iv,jv)%apply_ad(geom,vbal%np,vbal%h_n_s,vbal%h_c2b,vbal%h_S,fld2_blk(:,:,iv))
         call mpl%dot_prod(fld1_blk(:,:,iv),fld2_save(:,:,iv),sum1)
         call mpl%dot_prod(fld2_blk(:,:,iv),fld1_save(:,:,iv),sum2)
         write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Vertical balance block adjoint test:  ', &
         & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
         call mpl%flush
      end if
   end do
end do

! Direct adjoint test
fld1_dir = fld1_save
fld2_dir = fld2_save
call vbal%apply(nam,geom,bpar,fld1_dir)
call vbal%apply_ad(nam,geom,bpar,fld2_dir)

! Inverse adjoint test
fld1_inv = fld1_save
fld2_inv = fld2_save
call vbal%apply_inv(nam,geom,bpar,fld1_inv)
call vbal%apply_inv_ad(nam,geom,bpar,fld2_inv)

! Print result
call mpl%dot_prod(fld1_dir,fld2_save,sum1)
call mpl%dot_prod(fld2_dir,fld1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Vertical balance direct adjoint test:  ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush
call mpl%dot_prod(fld1_inv,fld2_save,sum1)
call mpl%dot_prod(fld2_inv,fld1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Vertical balance inverse adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

end subroutine vbal_test_adjoint

end module type_vbal
