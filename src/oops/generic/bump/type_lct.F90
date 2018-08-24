!----------------------------------------------------------------------
! Module: type_lct
!> Purpose: LCT data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_lct

use tools_const, only: req,reqkm,pi
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr
use type_bpar, only: bpar_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_io, only: io_type
use type_lct_blk, only: lct_blk_type
use type_mom, only: mom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

logical,parameter :: write_cor = .true.              !< Write raw and fitted correlations

! LCT data derived type
type lct_type
   type(lct_blk_type),allocatable :: blk(:) !< LCT blocks
   logical :: allocated                     !< Allocation flag
contains
   procedure :: alloc => lct_alloc
   procedure :: dealloc => lct_dealloc
   procedure :: run_lct => lct_run_lct
   procedure :: compute => lct_compute
   procedure :: filter => lct_filter
   procedure :: rmse => lct_rmse
   procedure :: write => lct_write
   procedure :: write_cor => lct_write_cor
end type lct_type

private
public :: lct_type

contains

!----------------------------------------------------------------------
! Subroutine: lct_alloc
!> Purpose: LCT data allocation
!----------------------------------------------------------------------
subroutine lct_alloc(lct,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct !< LCT
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data

! Local variables
integer :: ib

! Allocation
allocate(lct%blk(bpar%nb))
do ib=1,bpar%nb
   call lct%blk(ib)%alloc(nam,geom,bpar,hdata,ib)
end do

! Update allocation flag
lct%allocated = .true.

end subroutine lct_alloc

!----------------------------------------------------------------------
! Subroutine: lct_dealloc
!> Purpose: LCT data deallocation
!----------------------------------------------------------------------
subroutine lct_dealloc(lct,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct !< LCT
type(bpar_type),intent(in) :: bpar   !< Block parameters

! Local variables
integer :: ib

! Release memory
if (allocated(lct%blk)) then
   do ib=1,bpar%nb
      call lct%blk(ib)%dealloc
   end do
   deallocate(lct%blk)
end if

! Update allocation flag
lct%allocated = .false.

end subroutine lct_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_run_lct
!> Purpose: LCT driver
!----------------------------------------------------------------------
subroutine lct_run_lct(lct,mpl,rng,nam,geom,bpar,io,ens)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct !< LCT
type(mpl_type),intent(inout) :: mpl  !< MPI data
type(rng_type),intent(inout) :: rng  !< Random number generator
type(nam_type),intent(inout) :: nam  !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(io_type),intent(in) :: io       !< I/O
type(ens_type),intent(in) :: ens     !< Ensemble

! Local variables
type(hdata_type) :: hdata
type(mom_type) :: mom

! Setup sampling
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'
call flush(mpl%unit)

! Set artificially small local radius
nam%local_rad = 1.0e-12

! Setup sampling
call hdata%setup_sampling(mpl,rng,nam,geom,io)

! Compute MPI distribution, halo A
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A'
call flush(mpl%unit)
call hdata%compute_mpi_a(mpl,nam,geom)

! Compute MPI distribution, halos A-B
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute MPI distribution, halos A-B'
call flush(mpl%unit)
call hdata%compute_mpi_ab(mpl,nam,geom)

! Compute MPI distribution, halo C
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute MPI distribution, halo C'
call flush(mpl%unit)
call hdata%compute_mpi_c(mpl,nam,geom)

! Compute sample moments
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute sample moments'
call flush(mpl%unit)
call mom%compute(mpl,nam,geom,bpar,hdata,ens)

! Compute LCT
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute LCT'
call flush(mpl%unit)
call lct%compute(mpl,nam,geom,bpar,hdata,mom)

! Filter LCT
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Filter LCT'
call flush(mpl%unit)
call lct%filter(mpl,nam,geom,bpar,hdata)

! LCT RMSE
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- LCT RMSE'
call flush(mpl%unit)
call lct%rmse(mpl,nam,geom,bpar,hdata)

! Write LCT
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Write LCT'
call flush(mpl%unit)
call lct%write(mpl,nam,geom,bpar,io,hdata)

if (write_cor) then
   ! Write correlation and LCT fit
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write correlation and LCT fit'
   call flush(mpl%unit)
   call lct%write_cor(mpl,nam,geom,bpar,io,hdata)
end if

end subroutine lct_run_lct

!----------------------------------------------------------------------
! Subroutine: lct_compute
!> Purpose: compute LCT
!----------------------------------------------------------------------
subroutine lct_compute(lct,mpl,nam,geom,bpar,hdata,mom)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct !< LCT
type(mpl_type),intent(in) :: mpl     !< MPI data
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data
type(mom_type),intent(in) :: mom     !< Moments

! Local variables
integer :: ib

! Allocation
call lct%alloc(nam,geom,bpar,hdata)

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
   call flush(mpl%unit)

   ! Compute correlation
   write(mpl%unit,'(a10,a)') '','Compute correlation'
   call flush(mpl%unit)
   call lct%blk(ib)%correlation(nam,geom,bpar,hdata,mom%blk(ib))

   ! Compute LCT fit
   write(mpl%unit,'(a10,a)') '','Compute LCT fit'
   call flush(mpl%unit)
   call lct%blk(ib)%fitting(mpl,nam,geom,bpar,hdata)
end do

end subroutine lct_compute

!----------------------------------------------------------------------
! Subroutine: lct_filter
!> Purpose: filter LCT
!----------------------------------------------------------------------
subroutine lct_filter(lct,mpl,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    !< LCT
type(mpl_type),intent(in) :: mpl        !< MPI data
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,il0,ic1a,ic1,icomp,iscales,offset,nmsr,nmsr_tot
real(kind_real) :: fld_c1a(hdata%nc1a)
logical :: mask_c1a(hdata%nc1a,geom%nl0)

! Define mask
do il0=1,geom%nl0
   do ic1a=1,hdata%nc1a
      ic1 = hdata%c1a_to_c1(ic1a)
      mask_c1a(ic1a,il0) = hdata%c1l0_log(ic1,il0)
   end do
end do

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
   call flush(mpl%unit)

   do il0=1,geom%nl0
      ! Count missing LCT
      nmsr = 0
      do ic1a=1,hdata%nc1a
         if (mask_c1a(ic1a,il0).and.(.not.(all(isnotmsr(lct%blk(ib)%coef(:,ic1a,il0))) &
      & .and.all(isnotmsr(lct%blk(ib)%coef(:,ic1a,il0)))))) nmsr = nmsr+1
      end do
      call mpl%allreduce_sum(nmsr,nmsr_tot)
      write(mpl%unit,'(a10,a,i3,a,i8,a)',advance='no') '','Level',nam%levs(il0),': ',nmsr_tot,' missing points'

      if (nmsr_tot>0) then
         offset = 0
         do iscales=1,lct%blk(ib)%nscales
            do icomp=1,lct%blk(ib)%ncomp(iscales)+1
               ! Copy
               if (icomp<=lct%blk(ib)%ncomp(iscales)) then
                  fld_c1a = lct%blk(ib)%D(offset+icomp,:,il0)
               else
                  fld_c1a = lct%blk(ib)%coef(iscales,:,il0)
               end if

               ! Fill missing values
               call hdata%diag_fill(mpl,nam,geom,il0,fld_c1a)

               ! Copy
               if (icomp<=lct%blk(ib)%ncomp(iscales)) then
                  lct%blk(ib)%D(offset+icomp,:,il0) = fld_c1a
               else
                  lct%blk(ib)%coef(iscales,:,il0) = fld_c1a
               end if
            end do

            ! Update offset
            offset = offset+lct%blk(ib)%ncomp(iscales)
         end do
      end if

      ! Count missing LCT
      nmsr = 0
      do ic1a=1,hdata%nc1a
         if (mask_c1a(ic1a,il0).and.(.not.(all(isnotmsr(lct%blk(ib)%coef(:,ic1a,il0))) &
      & .and.all(isnotmsr(lct%blk(ib)%coef(:,ic1a,il0)))))) nmsr = nmsr+1
      end do
      call mpl%allreduce_sum(nmsr,nmsr_tot)
      write(mpl%unit,'(a,i8,a)') ' ~> ',nmsr_tot,' missing points'
   end do
end do

end subroutine lct_filter

!----------------------------------------------------------------------
! Subroutine: lct_rmse
!> Purpose: compute LCT fit RMSE
!----------------------------------------------------------------------
subroutine lct_rmse(lct,mpl,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_type),intent(in) :: lct       !< LCT
type(mpl_type),intent(in) :: mpl        !< MPI data
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,il0,jl0r,jl0,ic1a,ic1,jc3
real(kind_real) :: rmse,norm,rmse_tot,norm_tot

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
   call flush(mpl%unit)

   ! Compute RMSE
   rmse = 0.0
   norm = 0.0
   do il0=1,geom%nl0
      do ic1a=1,hdata%nc1a
         ic1 = hdata%c1a_to_c1(ic1a)
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               if (hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) then
                  if (isnotmsr(lct%blk(ib)%fit(jc3,jl0r,ic1a,il0))) then
                     rmse = rmse+(lct%blk(ib)%fit(jc3,jl0r,ic1a,il0)-lct%blk(ib)%raw(jc3,jl0r,ic1a,il0))**2
                     norm = norm+1.0
                  end if
               end if
            end do
         end do
      end do
   end do
   call mpl%allreduce_sum(rmse,rmse_tot)
   call mpl%allreduce_sum(norm,norm_tot)
   if (norm_tot>0.0) rmse_tot = sqrt(rmse_tot/norm_tot)
   write(mpl%unit,'(a10,a,e15.8,a,i8,a)') '','LCT diag RMSE: ',rmse_tot,' for ',int(norm_tot),' diagnostic points'
   call flush(mpl%unit)
end do

end subroutine lct_rmse

!----------------------------------------------------------------------
! Subroutine: lct_write
!> Purpose: interpolate and write LCT
!----------------------------------------------------------------------
subroutine lct_write(lct,mpl,nam,geom,bpar,io,hdata)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    !< LCT
type(mpl_type),intent(inout) :: mpl     !< MPI data
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(io_type),intent(in) :: io          !< I/O
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,iv,il0,il0i,ic1a,ic1,icomp,ic0a,ic0,iscales,offset
real(kind_real) :: det,diag_prod,Lavg,Lavg_tot,norm,norm_tot
real(kind_real),allocatable :: fld_c1a(:,:,:),fld_c1b(:,:),fld(:,:,:)
logical :: valid_coef,mask_c1a(hdata%nc1a,geom%nl0)
character(len=1) :: iscaleschar
character(len=1024) :: filename


! Define mask
do il0=1,geom%nl0
   do ic1a=1,hdata%nc1a
      ic1 = hdata%c1a_to_c1(ic1a)
      mask_c1a(ic1a,il0) = hdata%c1l0_log(ic1,il0)
   end do
end do

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
   call flush(mpl%unit)

   ! Initialization
   offset = 0

   do iscales=1,lct%blk(ib)%nscales
      write(mpl%unit,'(a10,a,i2)') '','Scale: ',iscales

      ! Allocation
      allocate(fld_c1a(hdata%nc1a,geom%nl0,lct%blk(ib)%ncomp(iscales)+1))
      allocate(fld_c1b(hdata%nc2b,geom%nl0))
      allocate(fld(geom%nc0a,geom%nl0,lct%blk(ib)%ncomp(iscales)+2))

      ! Initialization
      call msr(fld_c1a)
      call msr(fld)

      ! Invert LCT to get DT
      write(mpl%unit,'(a13,a)') '','Invert LCT to get DT '
      call flush(mpl%unit)
      do il0=1,geom%nl0
         do ic1a=1,hdata%nc1a
            ic1 = hdata%c1a_to_c1(ic1a)
            if (mask_c1a(ic1a,il0)) then
               ! Check D determinant
               diag_prod = lct%blk(ib)%D(offset+1,ic1a,il0)*lct%blk(ib)%D(offset+2,ic1a,il0)
               if (bpar%nl0r(ib)>1) det = det*lct%blk(ib)%D(offset+3,ic1a,il0)
               if (lct%blk(ib)%ncomp(iscales)==3) then
                  det = diag_prod
               else
                  det = diag_prod*(1.0-lct%blk(ib)%D(offset+4,ic1a,il0)**2)
               end if
               if (bpar%nl0r(ib)>1) det = det*lct%blk(ib)%D(offset+3,ic1a,il0)

               ! Check coefficient
               valid_coef = (lct%blk(ib)%coef(iscales,ic1a,il0)>0.0)
               if (lct%blk(ib)%nscales>1) valid_coef = valid_coef.and.(lct%blk(ib)%coef(iscales,ic1a,il0)<1.0)

               if ((det>0.0).and.valid_coef) then
                  ! Copy diffusion tensor
                  fld_c1a(ic1a,il0,1) = lct%blk(ib)%D(offset+1,ic1a,il0)
                  fld_c1a(ic1a,il0,2) = lct%blk(ib)%D(offset+2,ic1a,il0)
                  if (bpar%nl0r(ib)==1) then
                     fld_c1a(ic1a,il0,3) = 0.0
                  else
                     fld_c1a(ic1a,il0,3) = lct%blk(ib)%D(offset+3,ic1a,il0)
                  end if
                  if (lct%blk(ib)%ncomp(iscales)==4) fld_c1a(ic1a,il0,4) = lct%blk(ib)%D(offset+4,ic1a,il0)

                  ! Copy coefficient
                  fld_c1a(ic1a,il0,lct%blk(ib)%ncomp(iscales)+1) = lct%blk(ib)%coef(iscales,ic1a,il0)
               else
                  call mpl%abort('non-valid LCT, grid c1')
               end if
            end if
         end do
      end do

      ! Interpolate DT
      write(mpl%unit,'(a13,a)') '','Interpolate DT'
      call flush(mpl%unit)
      do icomp=1,lct%blk(ib)%ncomp(iscales)+1
         call hdata%com_AB%ext(mpl,geom%nl0,fld_c1a(:,:,icomp),fld_c1b)
         do il0=1,geom%nl0
            il0i = min(il0,geom%nl0i)
            call hdata%h(il0i)%apply(mpl,fld_c1b(:,il0),fld(:,il0,icomp))
         end do
      end do

      ! Compute horizontal length-scale
      write(mpl%unit,'(a13,a)') '','Compute horizontal length-scale:'
      call flush(mpl%unit)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            ic0 = geom%c0a_to_c0(ic0a)
            if (geom%mask(ic0,il0)) then
               ! Check D determinant
               if (lct%blk(ib)%ncomp(iscales)==3) then
                  det = fld(ic0a,il0,1)*fld(ic0a,il0,2)
               else
                  det = fld(ic0a,il0,1)*fld(ic0a,il0,2)*(1.0-fld(ic0a,il0,4)**2)
               end if
               if (det>0.0) then
                  ! Length-scale = D determinant^{1/4}
                  fld(ic0a,il0,lct%blk(ib)%ncomp(iscales)+2) = sqrt(sqrt(det))
               else
                  call mpl%abort('non-valid horizontal determinant in LCT, grid c0')
               end if

               ! Check coefficient
               valid_coef = (fld(ic0a,il0,lct%blk(ib)%ncomp(iscales)+1)>0.0)
               if (lct%blk(ib)%nscales>1) valid_coef = valid_coef.and.(fld(ic0a,il0,lct%blk(ib)%ncomp(iscales)+1)<1.0)
               if (.not.valid_coef) call mpl%abort('non-valid coefficient in LCT, grid c0')
            end if
         end do
         Lavg = sum(fld(:,il0,lct%blk(ib)%ncomp(iscales)+2),isnotmsr(fld(:,il0,lct%blk(ib)%ncomp(iscales)+2)))
         norm = real(count(isnotmsr(fld(:,il0,lct%blk(ib)%ncomp(iscales)+2))),kind_real)
         call mpl%allreduce_sum(Lavg,Lavg_tot)
         call mpl%allreduce_sum(norm,norm_tot)
         if (norm_tot>0.0) write(mpl%unit,'(a13,a,i3,a,f10.2,a)') '','Level',nam%levs(il0),' ~> ',Lavg_tot/norm_tot*reqkm,' km'
      end do

      ! Copy to LCT
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            ic0 = geom%c0a_to_c0(ic0a)
            if (geom%mask(ic0,il0)) then
               lct%blk(ib)%D11(ic0a,il0,iscales) = fld(ic0a,il0,1)
               lct%blk(ib)%D22(ic0a,il0,iscales) = fld(ic0a,il0,2)
               lct%blk(ib)%D33(ic0a,il0,iscales) = fld(ic0a,il0,3)
               if (lct%blk(ib)%ncomp(iscales)==4) then
                  lct%blk(ib)%D12(ic0a,il0,iscales) = sqrt(fld(ic0a,il0,1)*fld(ic0a,il0,2))*fld(ic0a,il0,4)
               else
                  lct%blk(ib)%D12(ic0a,il0,iscales) = 0.0
               end if
               lct%blk(ib)%Dcoef(ic0a,il0,iscales) = fld(ic0a,il0,lct%blk(ib)%ncomp(iscales)+1)
               lct%blk(ib)%DLh(ic0a,il0,iscales) = fld(ic0a,il0,lct%blk(ib)%ncomp(iscales)+2)
            end if
         end do
      end do

      ! Write LCT
      write(mpl%unit,'(a13,a)') '','Write LCT'
      call flush(mpl%unit)
      filename = trim(nam%prefix)//'_lct'
      iv = bpar%b_to_v2(ib)
      write(iscaleschar,'(i1)') iscales
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D11_'//iscaleschar,lct%blk(ib)%D11(:,:,iscales))
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D22_'//iscaleschar,lct%blk(ib)%D22(:,:,iscales))
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D33_'//iscaleschar,lct%blk(ib)%D33(:,:,iscales))
      if (lct%blk(ib)%ncomp(iscales)==4) call io%fld_write(mpl,nam,geom,filename, &
    & trim(nam%varname(iv))//'_D12_'//iscaleschar,fld(:,:,4)*req**2)
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_coef_'//iscaleschar,lct%blk(ib)%Dcoef(:,:,iscales))
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_Lh_'//iscaleschar,lct%blk(ib)%DLh(:,:,iscales))

      ! Update offset
      offset = offset+lct%blk(ib)%ncomp(iscales)

      ! Release memory
      deallocate(fld_c1a)
      deallocate(fld_c1b)
      deallocate(fld)
   end do
end do

end subroutine lct_write

!----------------------------------------------------------------------
! Subroutine: lct_write_cor
!> Purpose: write correlation and LCT fit
!----------------------------------------------------------------------
subroutine lct_write_cor(lct,mpl,nam,geom,bpar,io,hdata)

implicit none

! Passed variables
class(lct_type),intent(in) :: lct       !< LCT
type(mpl_type),intent(inout) :: mpl     !< MPI data
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(io_type),intent(in) :: io          !< I/O
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,iv,il0,jl0r,jl0,ic1a,ic1,jc3,i,iproc,ic0
real(kind_real) :: fld_glb(geom%nc0,geom%nl0,2),fld(geom%nc0a,geom%nl0,2)
real(kind_real),allocatable :: sbuf(:),rbuf(:)
logical :: valid
logical :: free(geom%nc0,geom%nl0)
character(len=1024) :: filename

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
   call flush(mpl%unit)

   ! Allocation
   if (mpl%main) allocate(rbuf(nam%nc3*bpar%nl0r(ib)*2))

   ! Select level
   il0 = 1

   ! Prepare field
   call msr(fld_glb)
   free = .true.
   do ic1=1,nam%nc1
      ! Select tensor to plot
      valid  = .true.
      do jl0r=1,bpar%nl0r(ib)
         jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
         do jc3=1,nam%nc3
            if (valid.and.hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) &
         &  valid = valid.and.free(hdata%c1c3_to_c0(ic1,jc3),jl0)
         end do
      end do

      if (valid) then
         ! Remember that the footprint is not free anymore
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               free(hdata%c1c3_to_c0(ic1,jc3),jl0) = .false.
            end do
         end do

         ! Find processor
         iproc = hdata%c2_to_proc(ic1)
         if (iproc==mpl%myproc) then
            ! Allocate buffer
            allocate(sbuf(nam%nc3*bpar%nl0r(ib)*2))

            ! Prepare buffer
            call msr(sbuf)
            ic1a = hdata%c1_to_c1a(ic1)
            i = 1
            do jl0r=1,bpar%nl0r(ib)
               jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
               do jc3=1,nam%nc3
                  if (hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) then
                     sbuf(i) = lct%blk(ib)%raw(jc3,jl0r,ic1a,il0)
                     sbuf(i+1) = lct%blk(ib)%fit(jc3,jl0r,ic1a,il0)
                  end if
                  i = i+2
               end do
            end do
         end if

         if (mpl%main) then
            if (iproc==mpl%ioproc) then
               ! Copy
               rbuf = sbuf
            else
               ! Receive data
               call mpl%recv(nam%nc3*bpar%nl0r(ib)*2,rbuf,iproc,mpl%tag)
            end if

            ! Fill field
            i = 1
            do jl0r=1,bpar%nl0r(ib)
               jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
               do jc3=1,nam%nc3
                  if (hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) then
                     ic0 = hdata%c1c3_to_c0(ic1,jc3)
                     fld_glb(ic0,jl0,1) = rbuf(i)
                     fld_glb(ic0,jl0,2) = rbuf(i+1)
                  end if
                  i = i+2
               end do
            end do
         else
            ! Send data
            if (iproc==mpl%myproc) call mpl%send(nam%nc3*bpar%nl0r(ib)*2,sbuf,mpl%ioproc,mpl%tag)
         end if
         call mpl%update_tag(1)

         ! Release memory
         if (iproc==mpl%myproc) deallocate(sbuf)
      end if
   end do

   ! Global to local
   do il0=1,geom%nl0
      call mpl%glb_to_loc(geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,fld_glb(:,il0,1),geom%nc0a,fld(:,il0,1))
      call mpl%glb_to_loc(geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,fld_glb(:,il0,2),geom%nc0a,fld(:,il0,2))
   end do

   ! Write LCT diagnostics
   write(mpl%unit,'(a10,a)') '','Write LCT diagnostics'
   call flush(mpl%unit)
   filename = trim(nam%prefix)//'_lct'
   iv = bpar%b_to_v2(ib)
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_raw',fld(:,:,1))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_fit',fld(:,:,2))

   ! Release memory
   if (mpl%main) deallocate(rbuf)
end do

end subroutine lct_write_cor

end module type_lct
