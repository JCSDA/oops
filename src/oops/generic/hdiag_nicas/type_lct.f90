!----------------------------------------------------------------------
! Module: type_lct
!> Purpose: LCT data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_lct

use model_interface, only: model_write
use tools_const, only: req,reqkm
use tools_display, only: prog_init,prog_print,msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_lct_blk, only: lct_blk_type
use type_mom, only: mom_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! LCT data derived type
type lct_type
   ! Attributes
   integer :: nscales                           !< Number of LCT scales
   integer,allocatable :: ncomp(:)              !< Number of LCT components

   ! Data
   type(lct_blk_type),allocatable :: blk(:,:,:) !< LCT blocks
contains
   procedure :: alloc => lct_alloc
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
!> Purpose: lct object allocation
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
integer :: iscales,ib,il0,ic1a

! Number of scales and components
lct%nscales = nam%lct_nscales
allocate(lct%ncomp(lct%nscales))
do iscales=1,lct%nscales
   if (nam%lct_diag(iscales)) then
      lct%ncomp(iscales) = 3
   else
      lct%ncomp(iscales) = 4
   end if
end do

! Allocation
allocate(lct%blk(hdata%nc1a,geom%nl0,bpar%nb))
do ib=1,bpar%nb
   do il0=1,geom%nl0
      do ic1a=1,hdata%nc1a
         call lct%blk(ic1a,il0,ib)%alloc(nam,bpar,ic1a,il0,ib)
      end do
   end do
end do

end subroutine lct_alloc

!----------------------------------------------------------------------
! Subroutine: lct_compute
!> Purpose: compute LCT
!----------------------------------------------------------------------
subroutine lct_compute(lct,nam,geom,bpar,hdata,mom)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct !< LCT
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
type(bpar_type),intent(in) :: bpar   !< Block parameters
type(hdata_type),intent(in) :: hdata !< HDIAG data
type(mom_type),intent(in) :: mom     !< Moments

! Local variables
integer :: ib,il0,ic1a,progint,missing_tot
logical,allocatable :: missing(:),done(:)

! Allocation
call lct%alloc(nam,geom,bpar,hdata)
allocate(missing(hdata%nc1a))
allocate(done(hdata%nc1a))

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

   do il0=1,geom%nl0
      write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),':'

      ! Initialization
      call prog_init(progint,done)

      do ic1a=1,hdata%nc1a
         ! Compute correlation
         call lct%blk(ic1a,il0,ib)%correlation(nam,bpar,mom%blk(ib))

         ! Compute LCT fit
         call lct%blk(ic1a,il0,ib)%fitting(nam,geom,bpar,hdata)
         missing(ic1a) = .not.isnotmsr(lct%blk(ic1a,il0,ib)%H(1))

         ! Update
         done(ic1a) = .true.
         call prog_print(progint,done)
      end do

      ! Gather missing points
      call mpl%allreduce_sum(count(missing),missing_tot)
      write(mpl%unit,'(a,i6,a)') '100% (',missing_tot,' missing points)'
   end do
end do

end subroutine lct_compute

!----------------------------------------------------------------------
! Subroutine: lct_filter
!> Purpose: filter LCT
!----------------------------------------------------------------------
subroutine lct_filter(lct,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    !< LCT
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,il0,ic1a,ic1,icomp,ic0,iscales,offset
real(kind_real) :: fac
real(kind_real),allocatable :: fld_c1a(:,:,:),fld_c1a_tmp(:,:)
logical :: valid,proc_to_valid(mpl%nproc)

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

   ! Initialization
   offset = 0

   do iscales=1,lct%nscales
      ! Allocation
      allocate(fld_c1a(hdata%nc1a,geom%nl0,lct%ncomp(iscales)+1))
      allocate(fld_c1a_tmp(hdata%nc1a,lct%ncomp(iscales)+1))

      ! Copy
      do il0=1,geom%nl0
         do ic1a=1,hdata%nc1a
            fld_c1a(ic1a,il0,1:lct%ncomp(iscales)) = lct%blk(ic1a,il0,ib)%H(offset+1:offset+lct%ncomp(iscales))
            fld_c1a(ic1a,il0,lct%ncomp(iscales)+1) = lct%blk(ic1a,il0,ib)%coef(iscales)
         end do
      end do

      do il0=1,geom%nl0
         write(mpl%unit,'(a10,a,i3,a,i3,a)') '','Scale ',iscales,', level',nam%levs(il0),': '

         ! Initialization
         fld_c1a_tmp = fld_c1a(:,il0,:)
         fac = 1.0

         ! Check invalid points
         valid = .true.
         do ic1a=1,hdata%nc1a
            ic1 = hdata%c1a_to_c1(ic1a)
            ic0 = hdata%c1_to_c0(ic1)
            if (geom%mask(ic0,il0).and.(.not.isallnotmsr(fld_c1a_tmp(ic1a,:)))) valid = .false.
         end do

         ! Gather validity results
         call mpl%allgather(1,(/valid/),proc_to_valid)

         do while (.not.all(proc_to_valid))
            ! Copy
            fld_c1a_tmp = fld_c1a(:,il0,:)

            ! Filter LCT
            write(mpl%unit,'(a13,a,f9.2,a)') '','Filter LCT with radius ',fac*nam%diag_rhflt*reqkm,' km'
            do icomp=1,lct%ncomp(iscales)+1
               call hdata%diag_filter(geom,il0,'median',fac*nam%diag_rhflt,fld_c1a_tmp(:,icomp))
               call hdata%diag_filter(geom,il0,'average',fac*nam%diag_rhflt,fld_c1a_tmp(:,icomp))
            end do

            ! Check invalid points
            valid = .true.
            do ic1a=1,hdata%nc1a
               ic1 = hdata%c1a_to_c1(ic1a)
               ic0 = hdata%c1_to_c0(ic1)
               if (geom%mask(ic0,il0).and.(.not.isallnotmsr(fld_c1a_tmp(ic1a,:)))) valid = .false.
            end do

            ! Gather validity results
            call mpl%allgather(1,(/valid/),proc_to_valid)

            ! Update fac (increase smoothing)
            fac = 2.0*fac
         end do

         ! Copy
         fld_c1a(:,il0,:) = fld_c1a_tmp
      end do

      ! Copy
      do il0=1,geom%nl0
         do ic1a=1,hdata%nc1a
            lct%blk(ic1a,il0,ib)%H(offset+1:offset+lct%ncomp(iscales)) = fld_c1a(ic1a,il0,1:lct%ncomp(iscales))
            lct%blk(ic1a,il0,ib)%coef(iscales) = fld_c1a(ic1a,il0,lct%ncomp(iscales)+1)
         end do
      end do

      ! Release memory
      deallocate(fld_c1a)
      deallocate(fld_c1a_tmp)
   end do
end do

end subroutine lct_filter

!----------------------------------------------------------------------
! Subroutine: lct_rmse
!> Purpose: compute LCT fit RMSE
!----------------------------------------------------------------------
subroutine lct_rmse(lct,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_type),intent(in) :: lct       !< LCT
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,il0,jl0r,jl0,ic1a,ic1,jc3
real(kind_real) :: rmse,norm,rmse_tot,norm_tot

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
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
                  if (isnotmsr(lct%blk(ic1a,il0,ib)%fit(jc3,jl0))) then
                     rmse = rmse+(lct%blk(ic1a,il0,ib)%fit(jc3,jl0)-lct%blk(ic1a,il0,ib)%raw(jc3,jl0))**2
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
end do

end subroutine lct_rmse

!----------------------------------------------------------------------
! Subroutine: lct_write
!> Purpose: interpolate and write LCT
!----------------------------------------------------------------------
subroutine lct_write(lct,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_type),intent(in) :: lct       !< LCT
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,iv,il0,il0i,ic1a,icomp,ic0a,ic0,iscales,offset
real(kind_real) :: det
real(kind_real),allocatable :: fld_c1a(:,:,:),fld_c1b(:,:),fld(:,:,:)
character(len=1) :: iscaleschar
character(len=1024) :: filename

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

   ! Initialization
   offset = 0

   do iscales=1,lct%nscales
      ! Allocation
      allocate(fld_c1a(hdata%nc1a,geom%nl0,lct%ncomp(iscales)+1))
      allocate(fld_c1b(hdata%nc2b,geom%nl0))
      allocate(fld(geom%nc0a,geom%nl0,lct%ncomp(iscales)+2))

      ! Copy
      do il0=1,geom%nl0
         do ic1a=1,hdata%nc1a
            fld_c1a(ic1a,il0,1:lct%ncomp(iscales)) = lct%blk(ic1a,il0,ib)%H(offset+1:offset+lct%ncomp(iscales))
            fld_c1a(ic1a,il0,lct%ncomp(iscales)+1) = lct%blk(ic1a,il0,ib)%coef(iscales)
         end do
      end do

      ! Interpolate LCT
      write(mpl%unit,'(a10,a)') '','Interpolate LCT'
      do icomp=1,lct%ncomp(iscales)+1
         call hdata%com_AB%ext(geom%nl0,fld_c1a(:,:,icomp),fld_c1b)
         do il0=1,geom%nl0
            il0i = min(il0,geom%nl0i)
            call hdata%h(il0i)%apply(fld_c1b,fld(:,:,icomp))
         end do
      end do

      ! Compute horizontal length-scale
      write(mpl%unit,'(a10,a)') '','Compute horizontal length-scale'
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            ic0 = geom%c0a_to_c0(ic0a)
            if (geom%mask(ic0,il0)) then
               ! Compute determinant
               if (lct%ncomp(iscales)==3) then
                  det = fld(ic0a,il0,1)*fld(ic0a,il0,2)
               else
                  det = fld(ic0a,il0,1)*fld(ic0a,il0,2)-fld(ic0a,il0,4)**2
               end if

               ! Length-scale = determinant^{1/4}
               if (det>0.0) then
                  fld(ic0a,il0,lct%ncomp(iscales)+2) = 1.0/sqrt(sqrt(det))
               else
                  call msgerror('non-positive determinant in LCT')
               end if
            end if
         end do
      end do

      ! Write gridded LCT
      write(mpl%unit,'(a10,a)') '','Write gridded LCT'
      filename = trim(nam%prefix)//'_lct_gridded.nc'
      iv = bpar%b_to_v2(ib)
      write(iscaleschar,'(i1)') iscales
      call model_write(nam,geom,filename,trim(nam%varname(iv))//'_H11_'//iscaleschar,fld(:,:,1)/req**2)
      call model_write(nam,geom,filename,trim(nam%varname(iv))//'_H22_'//iscaleschar,fld(:,:,2)/req**2)
      call model_write(nam,geom,filename,trim(nam%varname(iv))//'_H33_'//iscaleschar,fld(:,:,3))
      if (lct%ncomp(iscales)==4) call model_write(nam,geom,filename,trim(nam%varname(iv))//'_Hc12_'//iscaleschar,fld(:,:,4))
      call model_write(nam,geom,filename,trim(nam%varname(iv))//'_coef_'//iscaleschar,fld(:,:,lct%ncomp(iscales)+1))
      call model_write(nam,geom,filename,trim(nam%varname(iv))//'_Lh_'//iscaleschar,fld(:,:,lct%ncomp(iscales)+2)*reqkm)

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
subroutine lct_write_cor(lct,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(lct_type),intent(in) :: lct       !< LCT
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(bpar_type),intent(in) :: bpar      !< Block parameters
type(hdata_type),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ib,iv,il0,jl0r,jl0,ic1a,ic1,jc3,i,iproc,jproc,ic0
real(kind_real) :: fld_glb(geom%nc0,geom%nl0,2),fld(geom%nc0a,geom%nl0,2)
real(kind_real),allocatable :: sbuf(:),rbuf(:)
logical :: valid
logical :: free(geom%nc0,geom%nl0)
character(len=1024) :: filename

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

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
                     sbuf(i) = lct%blk(ic1a,il0,ib)%raw(jc3,jl0)
                     sbuf(i+1) = lct%blk(ic1a,il0,ib)%fit(jc3,jl0)
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
         mpl%tag = mpl%tag+1

         ! Release memory
         if (iproc==mpl%myproc) deallocate(sbuf)
      end if
   end do

   ! Global to local
   call geom%fld_com_gl(fld_glb(:,:,1),fld(:,:,1))
   call geom%fld_com_gl(fld_glb(:,:,2),fld(:,:,2))

   ! Write LCT diagnostics
   write(mpl%unit,'(a10,a)') '','Write LCT diagnostics'
   filename = trim(nam%prefix)//'_lct_gridded.nc'
   iv = bpar%b_to_v2(ib)
   call model_write(nam,geom,filename,trim(nam%varname(iv))//'_raw',fld(:,:,1))
   call model_write(nam,geom,filename,trim(nam%varname(iv))//'_fit',fld(:,:,2))

   ! Release memory
   if (mpl%main) deallocate(rbuf)
end do

end subroutine lct_write_cor

end module type_lct
