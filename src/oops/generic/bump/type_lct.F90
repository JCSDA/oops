!----------------------------------------------------------------------
! Module: type_lct
! Purpose: LCT data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_lct

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use tools_const, only: reqkm,pi
use tools_func, only: gau2gc,fit_lct,lct_d2h,check_cond,Dmin
use tools_kinds, only: kind_real
use type_bpar, only: bpar_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_lct_blk, only: lct_blk_type
use type_mom, only: mom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type
use type_rng, only: rng_type

implicit none

logical,parameter :: write_cor = .true. ! Write raw and fitted correlations

! LCT data derived type
type lct_type
   type(samp_type) :: samp                  ! Sampling
   type(mom_type) :: mom                    ! Moments
   type(lct_blk_type),allocatable :: blk(:) ! LCT blocks
   logical :: allocated                     ! Allocation flag
contains
   procedure :: alloc => lct_alloc
   procedure :: dealloc => lct_dealloc
   procedure :: run_lct => lct_run_lct
   procedure :: compute => lct_compute
   procedure :: filter => lct_filter
   procedure :: rmse => lct_rmse
   procedure :: interp => lct_interp
   procedure :: write => lct_write
   procedure :: write_cor => lct_write_cor
end type lct_type

private
public :: lct_type

contains

!----------------------------------------------------------------------
! Subroutine: lct_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine lct_alloc(lct,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib

! Allocation
allocate(lct%blk(bpar%nb))
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) call lct%blk(ib)%alloc(nam,geom,bpar,lct%samp,ib)
end do

! Update allocation flag
lct%allocated = .true.

end subroutine lct_alloc

!----------------------------------------------------------------------
! Subroutine: lct_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine lct_dealloc(lct)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT

! Local variables
integer :: ib

! Release memory
if (allocated(lct%blk)) then
   do ib=1,size(lct%blk)
      call lct%blk(ib)%dealloc
   end do
   deallocate(lct%blk)
end if

! Update allocation flag
lct%allocated = .false.

end subroutine lct_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_run_lct
! Purpose: LCT driver
!----------------------------------------------------------------------
subroutine lct_run_lct(lct,mpl,rng,nam,geom,bpar,io,ens)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(rng_type),intent(inout) :: rng  ! Random number generator
type(nam_type),intent(inout) :: nam  ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(io_type),intent(in) :: io       ! I/O
type(ens_type),intent(in) :: ens     ! Ensemble

! Setup sampling
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'
call mpl%flush

! Set artificially small local radius
nam%local_rad = 1.0e-12

! Setup sampling
call lct%samp%setup_sampling(mpl,rng,nam,geom,bpar,io,ens)

! Compute MPI distribution, halo A
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute MPI distribution, halos A'
call mpl%flush
call lct%samp%compute_mpi_a(mpl,nam,geom)

! Compute MPI distribution, halos A-B
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute MPI distribution, halos A-B'
call mpl%flush
call lct%samp%compute_mpi_ab(mpl,nam,geom)

! Compute MPI distribution, halo C
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute MPI distribution, halo C'
call mpl%flush
call lct%samp%compute_mpi_c(mpl,nam,geom)

if (nam%diag_rhflt>0.0) then
   ! Compute MPI distribution, halo F
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute MPI distribution, halo F'
   call mpl%flush
   call lct%samp%compute_mpi_f(mpl,nam)
end if

if (nam%new_mom) then
   ! Compute sample moments
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute sample moments'
   call mpl%flush
   call lct%mom%compute(mpl,nam,geom,bpar,lct%samp,ens,'mom')
else
   ! Load sample moments
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Load sample moments'
   call mpl%flush
   call lct%mom%read(mpl,nam,geom,bpar,lct%samp,ens,'mom')
end if

! Compute LCT
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute LCT'
call mpl%flush
call lct%compute(mpl,nam,geom,bpar)

! Filter LCT
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Filter LCT'
call mpl%flush
call lct%filter(mpl,nam,geom,bpar)

! LCT RMSE
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- LCT RMSE'
call mpl%flush
call lct%rmse(mpl,nam,geom,bpar)

! Interpolate LCT
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Interpolate LCT'
call mpl%flush
call lct%interp(mpl,nam,geom,bpar)

if (nam%write_lct) then
   ! Write LCT
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Write LCT'
   call mpl%flush
   call lct%write(mpl,nam,geom,bpar,io)

   if (write_cor) then
      ! Write correlation and LCT fit
      write(mpl%info,'(a)') '-------------------------------------------------------------------'
      call mpl%flush
      write(mpl%info,'(a)') '--- Write correlation and LCT fit'
      call mpl%flush
      call lct%write_cor(mpl,nam,geom,bpar,io)
   end if
end if

end subroutine lct_run_lct

!----------------------------------------------------------------------
! Subroutine: lct_compute
! Purpose: compute LCT
!----------------------------------------------------------------------
subroutine lct_compute(lct,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib

! Allocation
call lct%alloc(nam,geom,bpar)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Compute correlation
      write(mpl%info,'(a10,a)') '','Compute correlation'
      call mpl%flush
      call lct%blk(ib)%correlation(nam,geom,bpar,lct%samp,lct%mom%blk(ib))

      ! Compute LCT fit
      write(mpl%info,'(a10,a)') '','Compute LCT fit'
      call mpl%flush
      call lct%blk(ib)%fitting(mpl,nam,geom,bpar,lct%samp)
   end if
end do

end subroutine lct_compute

!----------------------------------------------------------------------
! Subroutine: lct_filter
! Purpose: filter LCT
!----------------------------------------------------------------------
subroutine lct_filter(lct,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: ib,il0,jl0,jl0r,ic1a,ic1,ic0,jc3,jc0,icomp,iscales,nmsr,nmsr_tot
real(kind_real) :: fld_c1a(lct%samp%nc1a)
real(kind_real),allocatable :: fld_filt_c1a(:),dx(:,:),dy(:,:),dz(:,:)
logical :: valid,mask_c1a(lct%samp%nc1a,geom%nl0)
logical,allocatable :: dmask(:,:)

! Allocation
if (nam%diag_rhflt>0.0) allocate(fld_filt_c1a(lct%samp%nc1a))

! Define mask
do il0=1,geom%nl0
   do ic1a=1,lct%samp%nc1a
      ic1 = lct%samp%c1a_to_c1(ic1a)
      mask_c1a(ic1a,il0) = lct%samp%c1l0_log(ic1,il0)
   end do
end do

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Allocation
      if (nam%diag_rhflt>0.0) then
         allocate(dx(nam%nc3,bpar%nl0r(ib)))
         allocate(dy(nam%nc3,bpar%nl0r(ib)))
         allocate(dz(nam%nc3,bpar%nl0r(ib)))
         allocate(dmask(nam%nc3,bpar%nl0r(ib)))
      end if

      do il0=1,geom%nl0
         ! Count missing LCT
         nmsr = 0
         do ic1a=1,lct%samp%nc1a
            if (mask_c1a(ic1a,il0).and.(.not.(all(mpl%msv%isnotr(lct%blk(ib)%coef(:,ic1a,il0))) &
         & .and.all(mpl%msv%isnotr(lct%blk(ib)%coef(:,ic1a,il0)))))) nmsr = nmsr+1
         end do
         call mpl%f_comm%allreduce(nmsr,nmsr_tot,fckit_mpi_sum())
         write(mpl%info,'(a10,a,i3,a,i8,a)') '','Level',nam%levs(il0),': ',nmsr_tot,' missing points'
         call mpl%flush(.false.)

         do iscales=1,lct%blk(ib)%nscales
            do icomp=1,4+1
               ! Copy
               if (icomp<=4) then
                  fld_c1a = lct%blk(ib)%D(icomp,iscales,:,il0)
               else
                  fld_c1a = lct%blk(ib)%coef(iscales,:,il0)
               end if

               if (nam%diag_rhflt>0.0) then
                  ! Copy
                  fld_filt_c1a = fld_c1a

                  ! Filter
                  call lct%samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,fld_filt_c1a)
                  call lct%samp%diag_filter(mpl,nam,'gc99',nam%diag_rhflt,fld_filt_c1a)
               end if

               ! Fill missing values
               if (nmsr_tot>0) then
                  call lct%samp%diag_fill(mpl,nam,fld_c1a)
                  if (nam%diag_rhflt>0.0) call lct%samp%diag_fill(mpl,nam,fld_filt_c1a)
               end if

               ! Copy
               if (icomp<=4) then
                  lct%blk(ib)%D(icomp,iscales,:,il0) = fld_c1a
                  if (nam%diag_rhflt>0.0) lct%blk(ib)%D_filt(icomp,iscales,:,il0) = fld_filt_c1a
               else
                  lct%blk(ib)%coef(iscales,:,il0) = fld_c1a
                  if (nam%diag_rhflt>0.0) lct%blk(ib)%coef_filt(iscales,:,il0) = fld_filt_c1a
               end if
            end do
         end do

         ! Count missing LCT
         nmsr = 0
         do ic1a=1,lct%samp%nc1a
            if (mask_c1a(ic1a,il0).and.(.not.(all(mpl%msv%isnotr(lct%blk(ib)%coef(:,ic1a,il0))) &
         & .and.all(mpl%msv%isnotr(lct%blk(ib)%coef(:,ic1a,il0)))))) nmsr = nmsr+1
         end do
         call mpl%f_comm%allreduce(nmsr,nmsr_tot,fckit_mpi_sum())
         write(mpl%info,'(a,i8,a)') ' ~> ',nmsr_tot,' missing points'
         call mpl%flush(.false.)
         if (nam%diag_rhflt>0.0) then
            write(mpl%info,'(a,f10.2,a)') ', filtering at ',nam%diag_rhflt*reqkm,' km'
            call mpl%flush
         else
            write(mpl%info,'(a,f10.2,a)') ', no filtering'
            call mpl%flush
         end if

         ! Copy variances
         lct%blk(ib)%m2flt = lct%blk(ib)%m2

         if (nam%diag_rhflt>0.0) then
            do ic1a=1,lct%samp%nc1a
               ! Global index
               ic1 = lct%samp%c1a_to_c1(ic1a)

               if (lct%samp%c1l0_log(ic1,il0)) then
                  ! Compute deltas
                  ic0 = lct%samp%c1_to_c0(ic1)
                  do jl0r=1,bpar%nl0r(ib)
                     jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
                     do jc3=1,nam%nc3
                        dmask(jc3,jl0r) = lct%samp%c1l0_log(ic1,il0).and.lct%samp%c1c3l0_log(ic1,jc3,jl0)
                        if (dmask(jc3,jl0r)) then
                           jc0 = lct%samp%c1c3_to_c0(ic1,jc3)
                           call geom%compute_deltas(ic0,il0,jc0,jl0,dx(jc3,jl0r),dy(jc3,jl0r),dz(jc3,jl0r))
                        end if
                     end do
                  end do

                  ! Check tensor validity
                  valid = .true.
                  do iscales=1,lct%blk(ib)%nscales
                     if (valid) then
                        call check_cond(lct%blk(ib)%D_filt(1,iscales,ic1a,il0),lct%blk(ib)%D_filt(2,iscales,ic1a,il0), &
                      & lct%blk(ib)%D_filt(4,iscales,ic1a,il0),valid)
                        if (bpar%nl0r(ib)>1) valid = valid.and.(lct%blk(ib)%D_filt(3,iscales,ic1a,il0)>0.0)
                        valid = valid.and.(lct%blk(ib)%coef_filt(iscales,ic1a,il0)>0.0)
                        if (lct%blk(ib)%nscales>1) valid = valid.and.(lct%blk(ib)%coef_filt(iscales,ic1a,il0)<1.0)
                     end if
                  end do
                  if (valid) then
                     ! Rebuild fit
                     call fit_lct(mpl,nam%nc3,bpar%nl0r(ib),dx,dy,dz,dmask,lct%blk(ib)%nscales, &
                   & lct%blk(ib)%D_filt(:,:,ic1a,il0),lct%blk(ib)%coef_filt(:,ic1a,il0),lct%blk(ib)%fit_filt(:,:,ic1a,il0))
                  else
                     ! Missing values
                     lct%blk(ib)%fit_filt(:,:,ic1a,il0) = mpl%msv%valr
                  end if
               end if
            end do

            ! Filter variances
            call lct%samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,lct%blk(ib)%m2flt)
            call lct%samp%diag_filter(mpl,nam,'gc99',nam%diag_rhflt,lct%blk(ib)%m2flt)
         end if
      end do

      ! Release memory
      if (nam%diag_rhflt>0.0) then
         deallocate(dx)
         deallocate(dy)
         deallocate(dz)
         deallocate(dmask)
      end if
   end if
end do

! Release memory
if (nam%diag_rhflt>0.0) deallocate(fld_filt_c1a)

end subroutine lct_filter

!----------------------------------------------------------------------
! Subroutine: lct_rmse
! Purpose: compute LCT fit RMSE
!----------------------------------------------------------------------
subroutine lct_rmse(lct,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(in) :: lct       ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: ib,il0,jl0r,jl0,ic1a,ic1,jc3
real(kind_real) :: rmse,norm,rmse_tot,norm_tot
real(kind_real) :: rmse_filt,norm_filt,rmse_filt_tot,norm_filt_tot

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Compute RMSE
      rmse = 0.0
      norm = 0.0
      if (nam%diag_rhflt>0.0) then
         rmse_filt = 0.0
         norm_filt = 0.0
      end if
      do il0=1,geom%nl0
         do ic1a=1,lct%samp%nc1a
            ic1 = lct%samp%c1a_to_c1(ic1a)
            do jl0r=1,bpar%nl0r(ib)
               jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
               do jc3=1,nam%nc3
                  if (lct%samp%c1l0_log(ic1,il0).and.lct%samp%c1c3l0_log(ic1,jc3,jl0)) then
                     if (mpl%msv%isnotr(lct%blk(ib)%fit(jc3,jl0r,ic1a,il0))) then
                        rmse = rmse+(lct%blk(ib)%fit(jc3,jl0r,ic1a,il0)-lct%blk(ib)%raw(jc3,jl0r,ic1a,il0))**2
                        norm = norm+1.0
                     end if
                     if (nam%diag_rhflt>0.0) then
                        if (mpl%msv%isnotr(lct%blk(ib)%fit_filt(jc3,jl0r,ic1a,il0))) then
                           rmse_filt = rmse_filt+(lct%blk(ib)%fit_filt(jc3,jl0r,ic1a,il0)-lct%blk(ib)%raw(jc3,jl0r,ic1a,il0))**2
                           norm_filt = norm_filt+1.0
                        end if
                     end if
                  end if
               end do
            end do
         end do
      end do
      call mpl%f_comm%allreduce(rmse,rmse_tot,fckit_mpi_sum())
      call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
      if (norm_tot>0.0) rmse_tot = sqrt(rmse_tot/norm_tot)
      write(mpl%info,'(a10,a,e15.8,a,i8,a)') '','LCT fit RMSE:          ',rmse_tot,' for ',int(norm_tot),' diagnostic points'
      call mpl%flush
      if (nam%diag_rhflt>0.0) then
         call mpl%f_comm%allreduce(rmse_filt,rmse_filt_tot,fckit_mpi_sum())
         call mpl%f_comm%allreduce(norm_filt,norm_filt_tot,fckit_mpi_sum())
         if (norm_filt_tot>0.0) rmse_filt_tot = sqrt(rmse_filt_tot/norm_filt_tot)
         write(mpl%info,'(a10,a,e15.8,a,i8,a)') '','LCT filtered fit RMSE: ',rmse_filt_tot,' for ',int(norm_tot), &
       & ' diagnostic points'
         call mpl%flush
      end if
   end if
end do

end subroutine lct_rmse

!----------------------------------------------------------------------
! Subroutine: lct_interp
! Purpose: interpolate LCT
!----------------------------------------------------------------------
subroutine lct_interp(lct,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: ib,il0,il0i,ic1a,ic1,icomp,ic0a,iscales
real(kind_real) :: det,Lavg_tot,norm_tot
real(kind_real) :: fld_c1a(lct%samp%nc1a,geom%nl0,2*4+2),fld_c1b(lct%samp%nc2b,geom%nl0),fld(geom%nc0a,geom%nl0,2*4+3)
real(kind_real),allocatable :: D(:,:,:,:),coef(:,:,:)
logical :: mask_c1a(lct%samp%nc1a,geom%nl0)
character(len=1024),parameter :: subr = 'lct_interp'

! Define mask
do il0=1,geom%nl0
   do ic1a=1,lct%samp%nc1a
      ic1 = lct%samp%c1a_to_c1(ic1a)
      mask_c1a(ic1a,il0) = lct%samp%c1l0_log(ic1,il0)
   end do
end do

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Allocation
      allocate(D(4,lct%blk(ib)%nscales,lct%samp%nc1a,geom%nl0))
      allocate(coef(lct%blk(ib)%nscales,lct%samp%nc1a,geom%nl0))

      ! Initialization
      if (nam%diag_rhflt>0.0) then
         D = lct%blk(ib)%D_filt
         coef = lct%blk(ib)%coef_filt
      else
         D = lct%blk(ib)%D
         coef = lct%blk(ib)%coef
      end if

      do iscales=1,lct%blk(ib)%nscales
         write(mpl%info,'(a10,a,i2)') '','Scale: ',iscales
         call mpl%flush

         ! Initialization
         fld_c1a = mpl%msv%valr
         fld = mpl%msv%valr

         ! Copy and inverse diffusion tensor
         write(mpl%info,'(a13,a)') '','Copy and inverse diffusion tensor'
         call mpl%flush
         do il0=1,geom%nl0
            do ic1a=1,lct%samp%nc1a
               ic1 = lct%samp%c1a_to_c1(ic1a)
               if (mask_c1a(ic1a,il0)) then
                  ! Ensure positive-definiteness of D
                  D(1,iscales,ic1a,il0) = max(Dmin,D(1,iscales,ic1a,il0))
                  D(2,iscales,ic1a,il0) = max(Dmin,D(2,iscales,ic1a,il0))
                  if (bpar%nl0r(ib)>1) D(3,iscales,ic1a,il0) = max(Dmin,D(3,iscales,ic1a,il0))
                  D(4,iscales,ic1a,il0) = max(-1.0_kind_real+Dmin,min(D(4,iscales,ic1a,il0),1.0_kind_real-Dmin))

                  ! Copy diffusion tensor
                  fld_c1a(ic1a,il0,1) = D(1,iscales,ic1a,il0)
                  fld_c1a(ic1a,il0,2) = D(2,iscales,ic1a,il0)
                  fld_c1a(ic1a,il0,3) = D(3,iscales,ic1a,il0)
                  fld_c1a(ic1a,il0,4) = sqrt(D(1,iscales,ic1a,il0)*D(2,iscales,ic1a,il0))*D(4,iscales,ic1a,il0)

                  ! Inverse diffusion tensor
                  call lct_d2h(mpl,fld_c1a(ic1a,il0,1),fld_c1a(ic1a,il0,2),fld_c1a(ic1a,il0,3),fld_c1a(ic1a,il0,4), &
                & fld_c1a(ic1a,il0,4+1),fld_c1a(ic1a,il0,4+2),fld_c1a(ic1a,il0,4+3),fld_c1a(ic1a,il0,4+4))

                  ! Copy coefficient
                  fld_c1a(ic1a,il0,2*4+1) = coef(iscales,ic1a,il0)

                  ! Copy filtered variances
                  fld_c1a(ic1a,il0,2*4+2) = lct%blk(ib)%m2flt(ic1a,il0)
               end if
            end do
         end do

         ! Interpolate components
         write(mpl%info,'(a13,a)') '','Interpolate components'
         call mpl%flush
         do icomp=1,2*4+2
            call lct%samp%com_AB%ext(mpl,geom%nl0,fld_c1a(:,:,icomp),fld_c1b)
            do il0=1,geom%nl0
               il0i = min(il0,geom%nl0i)
               call lct%samp%h(il0i)%apply(mpl,fld_c1b(:,il0),fld(:,il0,icomp))
            end do
         end do

         ! Compute horizontal length-scale and equivalent support radius
         write(mpl%info,'(a13,a)') '','Compute horizontal length-scale and equivalent support radius:'
         call mpl%flush
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask_c0a(ic0a,il0)) then
                  ! Length-scale = D determinant^{1/4}
                  det = fld(ic0a,il0,1)*fld(ic0a,il0,2)-fld(ic0a,il0,4)**2
                  if (det>0.0) then
                     fld(ic0a,il0,2*4+3) = sqrt(sqrt(det))
                  else
                     call mpl%abort(subr,'non-valid horizontal diffusion tensor determinant, grid c0')
                  end if
               end if
            end do
            call mpl%f_comm%allreduce(sum(fld(:,il0,2*4+3),mpl%msv%isnotr(fld(:,il0,2*3+3))),Lavg_tot,fckit_mpi_sum())
            call mpl%f_comm%allreduce(real(count(mpl%msv%isnotr(fld(:,il0,2*4+3))),kind_real),norm_tot,fckit_mpi_sum())
            if (norm_tot>0.0) then
               write(mpl%info,'(a16,a,i3,a,f10.2,a,f10.2,a)') '','Level',nam%levs(il0),' ~> ', &
             & Lavg_tot/norm_tot*reqkm,' km / ',Lavg_tot/norm_tot*gau2gc*reqkm,' km'
               call mpl%flush
            end if
         end do

         ! Copy output values
         lct%blk(ib)%D11(:,:,iscales) = fld(:,:,1)
         lct%blk(ib)%D22(:,:,iscales) = fld(:,:,2)
         lct%blk(ib)%D33(:,:,iscales) = fld(:,:,3)
         lct%blk(ib)%D12(:,:,iscales) = fld(:,:,4)
         lct%blk(ib)%H11(:,:,iscales) = fld(:,:,4+1)
         lct%blk(ib)%H22(:,:,iscales) = fld(:,:,4+2)
         lct%blk(ib)%H33(:,:,iscales) = fld(:,:,4+3)
         lct%blk(ib)%H12(:,:,iscales) = fld(:,:,4+4)
         lct%blk(ib)%Dcoef(:,:,iscales) = fld(:,:,2*4+1)
         lct%blk(ib)%coef_ens = fld(:,:,2*4+2)
         lct%blk(ib)%DLh(:,:,iscales) = fld(:,:,2*4+3)
      end do

      ! Allocation
      deallocate(D)
      deallocate(coef)
   end if
end do

end subroutine lct_interp

!----------------------------------------------------------------------
! Subroutine: lct_write
! Purpose: write LCT
!----------------------------------------------------------------------
subroutine lct_write(lct,mpl,nam,geom,bpar,io)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters
type(io_type),intent(in) :: io          ! I/O

! Local variables
integer :: ib,iv,iscales
character(len=1) :: iscaleschar
character(len=1024) :: filename

! Set file name
filename = trim(nam%prefix)//'_lct'
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)

! Write vertical unit
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      iv = bpar%b_to_v2(ib)
      do iscales=1,lct%blk(ib)%nscales
         ! Write fields
         write(iscaleschar,'(i1)') iscales
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D11_'//iscaleschar,lct%blk(ib)%D11(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D22_'//iscaleschar,lct%blk(ib)%D22(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D33_'//iscaleschar,lct%blk(ib)%D33(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D12_'//iscaleschar,lct%blk(ib)%D12(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H11_'//iscaleschar,lct%blk(ib)%H11(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H22_'//iscaleschar,lct%blk(ib)%H22(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H33_'//iscaleschar,lct%blk(ib)%H33(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H12_'//iscaleschar,lct%blk(ib)%H12(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_coef_'//iscaleschar,lct%blk(ib)%Dcoef(:,:,iscales))
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_Lh_'//iscaleschar,lct%blk(ib)%DLh(:,:,iscales))
      end do
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_coef_ens',lct%blk(ib)%coef_ens)
   end if
end do

end subroutine lct_write

!----------------------------------------------------------------------
! Subroutine: lct_write_cor
! Purpose: write correlation and LCT fit
!----------------------------------------------------------------------
subroutine lct_write_cor(lct,mpl,nam,geom,bpar,io)

implicit none

! Passed variables
class(lct_type),intent(in) :: lct       ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters
type(io_type),intent(in) :: io          ! I/O

! Local variables
integer :: ib,iv,il0,jl0r,jl0,ic1a,ic1,jc3,i,iproc,ic0,nf
real(kind_real),allocatable :: fld_c0a(:,:,:),fld_c0(:,:,:),sbuf(:),rbuf(:)
logical :: valid
logical :: free(geom%nc0,geom%nl0)
character(len=1024) :: filename
type(fckit_mpi_status) :: status

! Number of fields
if (nam%diag_rhflt>0.0) then
   nf = 3
else
   nf = 2
end if

! Allocation
allocate(fld_c0a(geom%nc0a,geom%nl0,nf))

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Allocation
      if (mpl%main) allocate(rbuf(nam%nc3*bpar%nl0r(ib)*nf))
      allocate(fld_c0(geom%nc0,geom%nl0,nf))

      ! Select level
      il0 = 1

      ! Prepare field
      if (mpl%main) fld_c0 = mpl%msv%valr
      free = .true.
      do ic1=1,nam%nc1
         ! Select tensor to plot
         valid  = .true.
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               if (valid.and.lct%samp%c1l0_log(ic1,il0).and.lct%samp%c1c3l0_log(ic1,jc3,jl0)) &
            &  valid = valid.and.free(lct%samp%c1c3_to_c0(ic1,jc3),jl0)
            end do
         end do

         if (valid) then
            ! Remember that the footprint is not free anymore
            do jl0r=1,bpar%nl0r(ib)
               jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
               do jc3=1,nam%nc3
                  free(lct%samp%c1c3_to_c0(ic1,jc3),jl0) = .false.
               end do
            end do

            ! Find processor
            iproc = lct%samp%c2_to_proc(ic1)
            if (iproc==mpl%myproc) then
               ! Allocate buffer
               allocate(sbuf(nam%nc3*bpar%nl0r(ib)*nf))

               ! Prepare buffer
               sbuf = mpl%msv%valr
               ic1a = lct%samp%c1_to_c1a(ic1)
               i = 1
               do jl0r=1,bpar%nl0r(ib)
                  jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
                  do jc3=1,nam%nc3
                     if (lct%samp%c1l0_log(ic1,il0).and.lct%samp%c1c3l0_log(ic1,jc3,jl0)) then
                        sbuf(i) = lct%blk(ib)%raw(jc3,jl0r,ic1a,il0)
                        sbuf(i+1) = lct%blk(ib)%fit(jc3,jl0r,ic1a,il0)
                        if (nf==3) sbuf(i+2) = lct%blk(ib)%fit_filt(jc3,jl0r,ic1a,il0)
                     end if
                     i = i+nf
                  end do
               end do
            end if

            if (mpl%main) then
               if (iproc==mpl%ioproc) then
                  ! Copy
                  rbuf = sbuf
               else
                  ! Receive data
                  call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
               end if

               ! Fill field
               i = 1
               do jl0r=1,bpar%nl0r(ib)
                  jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
                  do jc3=1,nam%nc3
                     if (lct%samp%c1l0_log(ic1,il0).and.lct%samp%c1c3l0_log(ic1,jc3,jl0)) then
                        ic0 = lct%samp%c1c3_to_c0(ic1,jc3)
                        fld_c0(ic0,jl0,1) = rbuf(i)
                        fld_c0(ic0,jl0,2) = rbuf(i+1)
                        if (nf==3) fld_c0(ic0,jl0,3) = rbuf(i+2)
                     end if
                     i = i+nf
                  end do
               end do
            else
               ! Send data
               if (iproc==mpl%myproc) call mpl%f_comm%send(sbuf,mpl%ioproc-1,mpl%tag)
            end if
            call mpl%update_tag(1)

            ! Release memory
            if (iproc==mpl%myproc) deallocate(sbuf)
         end if
      end do

      ! Global to local
      call mpl%glb_to_loc(geom%nl0,geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,fld_c0(:,:,1),geom%nc0a,fld_c0a(:,:,1))
      call mpl%glb_to_loc(geom%nl0,geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,fld_c0(:,:,2),geom%nc0a,fld_c0a(:,:,2))
      if (nf==3) call mpl%glb_to_loc(geom%nl0,geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,fld_c0(:,:,3),geom%nc0a,fld_c0a(:,:,3))

      ! Write LCT diagnostics
      write(mpl%info,'(a10,a)') '','Write LCT diagnostics'
      call mpl%flush
      filename = trim(nam%prefix)//'_lct'
      iv = bpar%b_to_v2(ib)
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_raw',fld_c0a(:,:,1))
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_fit',fld_c0a(:,:,2))
      if (nf==3) call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_fit_filt',fld_c0a(:,:,3))

      ! Release memory
      if (mpl%main) deallocate(rbuf)
      deallocate(fld_c0)
   end if
end do

! Release memory
deallocate(fld_c0a)

end subroutine lct_write_cor

end module type_lct
