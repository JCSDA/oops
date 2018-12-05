!----------------------------------------------------------------------
! Module: type_cmat
! Purpose: C matrix derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_cmat

use fckit_mpi_module, only: fckit_mpi_sum
use netcdf
use tools_const, only: rad2deg,reqkm,req
use tools_func, only: gau2gc
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsi,isnotmsr,isallnotmsr,isanynotmsr
use type_avg, only: avg_type
use type_bpar, only: bpar_type
use type_cmat_blk, only: cmat_blk_type
use type_diag, only: diag_type
use type_displ, only: displ_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_hdiag, only: hdiag_type
use type_io, only: io_type
use type_lct, only: lct_type
use type_mom, only: mom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type
use type_samp, only: samp_type

implicit none

! C matrix derived type
type cmat_type
   character(len=1024) :: prefix             ! Prefix
   type(cmat_blk_type),allocatable :: blk(:) ! C matrix blocks
   logical :: allocated                      ! Allocation flag
contains
   procedure :: cmat_alloc
   procedure :: cmat_alloc_blk
   generic :: alloc => cmat_alloc,cmat_alloc_blk
   procedure :: dealloc => cmat_dealloc
   procedure :: copy => cmat_copy
   procedure :: read => cmat_read
   procedure :: write => cmat_write
   procedure :: from_hdiag => cmat_from_hdiag
   procedure :: from_lct => cmat_from_lct
   procedure :: from_nam => cmat_from_nam
   procedure :: from_oops => cmat_from_oops
   procedure :: setup_sampling => cmat_setup_sampling
end type cmat_type

private
public :: cmat_type

contains

!----------------------------------------------------------------------
! Subroutine: cmat_alloc
! Purpose: C matrix allocation
!----------------------------------------------------------------------
subroutine cmat_alloc(cmat,bpar,prefix)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat    ! C matrix
type(bpar_type),intent(in) :: bpar        ! Block parameters
character(len=*),intent(in) :: prefix     ! Prefix

! Local variables
integer :: ib

! Copy prefix
cmat%prefix = prefix

! Allocation
if (.not.allocated(cmat%blk)) allocate(cmat%blk(bpar%nbe))

! Set block name
do ib=1,bpar%nbe
   cmat%blk(ib)%name = trim(prefix)//'_'//trim(bpar%blockname(ib))
end do

end subroutine cmat_alloc

!----------------------------------------------------------------------
! Subroutine: cmat_alloc_blk
! Purpose: C matrix block allocation
!----------------------------------------------------------------------
subroutine cmat_alloc_blk(cmat,nam,geom,bpar)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat    ! C matrix
type(nam_type),intent(in) :: nam          ! Namelist
type(geom_type),intent(in) :: geom        ! Geometry
type(bpar_type),intent(in) :: bpar        ! Block parameters

! Local variables
integer :: ib

! Allocation
do ib=1,bpar%nbe
   cmat%blk(ib)%ib = ib
   call cmat%blk(ib)%alloc(nam,geom,bpar)
end do

! Update allocation flag
cmat%allocated = .true.

end subroutine cmat_alloc_blk

!----------------------------------------------------------------------
! Subroutine: cmat_dealloc
! Purpose: C matrix allocation
!----------------------------------------------------------------------
subroutine cmat_dealloc(cmat,bpar)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat ! C matrix
type(bpar_type),intent(in) :: bpar     ! Block parameters

! Local variables
integer :: ib

! Release memory
if (allocated(cmat%blk)) then
   do ib=1,bpar%nbe
      call cmat%blk(ib)%dealloc
   end do
   deallocate(cmat%blk)
end if

! Update allocation flag
cmat%allocated = .false.

end subroutine cmat_dealloc

!----------------------------------------------------------------------
! Function: cmat_copy
! Purpose: C matrix copy
!----------------------------------------------------------------------
type(cmat_type) function cmat_copy(cmat,nam,geom,bpar)

implicit none

! Passed variables
class(cmat_type),intent(in) :: cmat       ! C matrix
type(nam_type),intent(in) :: nam          ! Namelist
type(geom_type),intent(in) :: geom        ! Geometry
type(bpar_type),intent(in) :: bpar        ! Block parameters

! Local variables
integer :: ib

! Allocation
call cmat_copy%alloc(bpar,trim(cmat%prefix))

! Copy attributes
do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      cmat_copy%blk(ib)%double_fit = cmat%blk(ib)%double_fit
      cmat_copy%blk(ib)%anisotropic = cmat%blk(ib)%anisotropic
   end if
end do

! Allocation
call cmat_copy%alloc(nam,geom,bpar)

! Copy
do ib=1,bpar%nbe
   if (allocated(cmat%blk(ib)%coef_ens)) cmat_copy%blk(ib)%coef_ens = cmat%blk(ib)%coef_ens
   if (allocated(cmat%blk(ib)%coef_sta)) cmat_copy%blk(ib)%coef_sta = cmat%blk(ib)%coef_sta
   if (allocated(cmat%blk(ib)%rh)) cmat_copy%blk(ib)%rh = cmat%blk(ib)%rh
   if (allocated(cmat%blk(ib)%rv)) cmat_copy%blk(ib)%rv = cmat%blk(ib)%rv
   if (allocated(cmat%blk(ib)%rv_rfac)) cmat_copy%blk(ib)%rv_rfac = cmat%blk(ib)%rv_rfac
   if (allocated(cmat%blk(ib)%rv_coef)) cmat_copy%blk(ib)%rv_coef = cmat%blk(ib)%rv_coef
   if (allocated(cmat%blk(ib)%rhs)) cmat_copy%blk(ib)%rhs = cmat%blk(ib)%rhs
   if (allocated(cmat%blk(ib)%rvs)) cmat_copy%blk(ib)%rvs = cmat%blk(ib)%rvs
   if (allocated(cmat%blk(ib)%H11)) cmat_copy%blk(ib)%H11 = cmat%blk(ib)%H11
   if (allocated(cmat%blk(ib)%H22)) cmat_copy%blk(ib)%H22 = cmat%blk(ib)%H22
   if (allocated(cmat%blk(ib)%H33)) cmat_copy%blk(ib)%H33 = cmat%blk(ib)%H33
   if (allocated(cmat%blk(ib)%H12)) cmat_copy%blk(ib)%H12 = cmat%blk(ib)%H12
   if (allocated(cmat%blk(ib)%Hcoef)) cmat_copy%blk(ib)%Hcoef = cmat%blk(ib)%Hcoef
   if (allocated(cmat%blk(ib)%displ_lon)) cmat_copy%blk(ib)%displ_lon = cmat%blk(ib)%displ_lon
   if (allocated(cmat%blk(ib)%displ_lat)) cmat_copy%blk(ib)%displ_lat = cmat%blk(ib)%displ_lat
end do

end function cmat_copy

!----------------------------------------------------------------------
! Subroutine: cmat_read
! Purpose: read C matrix
!----------------------------------------------------------------------
subroutine cmat_read(cmat,mpl,nam,geom,bpar,io)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat ! C matrix
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O

! Local variables
integer :: ib,ncid,double_fit,anisotropic
character(len=1024) :: filename
character(len=1024) :: subr = 'cmat_read'

! Allocation
call cmat%alloc(bpar,'cmat')

do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      if (mpl%main) then
         ! Set filename
         filename = trim(nam%prefix)//'_'//trim(cmat%blk(ib)%name)

         ! Read attributes
         call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))
         call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'double_fit',double_fit))
         if (double_fit==1) then
            cmat%blk(ib)%double_fit = .true.
         else
            cmat%blk(ib)%double_fit = .false.
         end if
         call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'anisotropic',anisotropic))
         if (anisotropic==1) then
            cmat%blk(ib)%anisotropic = .true.
         else
            cmat%blk(ib)%anisotropic = .false.
         end if
         call mpl%ncerr(subr,nf90_close(ncid))
      end if

      ! Broadcast
      call mpl%f_comm%broadcast(cmat%blk(ib)%double_fit,mpl%ioproc-1)
      call mpl%f_comm%broadcast(cmat%blk(ib)%anisotropic,mpl%ioproc-1)
   end if
end do

! Allocation
call cmat%alloc(nam,geom,bpar)

do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      ! Set filename
      filename = trim(nam%prefix)//'_'//trim(cmat%blk(ib)%name)

      ! Read fields
      if (bpar%nicas_block(ib)) then
         call io%fld_read(mpl,nam,geom,filename,'coef_ens',cmat%blk(ib)%coef_ens)
         call io%fld_read(mpl,nam,geom,filename,'coef_sta',cmat%blk(ib)%coef_sta)
         call io%fld_read(mpl,nam,geom,filename,'rh',cmat%blk(ib)%rh)
         call io%fld_read(mpl,nam,geom,filename,'rv',cmat%blk(ib)%rv)
         if (cmat%blk(ib)%double_fit) then
            call io%fld_read(mpl,nam,geom,filename,'rv_rfac',cmat%blk(ib)%rv_rfac)
            call io%fld_read(mpl,nam,geom,filename,'rv_coef',cmat%blk(ib)%rv_coef)
         end if
         call io%fld_read(mpl,nam,geom,filename,'rhs',cmat%blk(ib)%rhs)
         call io%fld_read(mpl,nam,geom,filename,'rvs',cmat%blk(ib)%rvs)
         if (cmat%blk(ib)%anisotropic) then
            call io%fld_read(mpl,nam,geom,filename,'H11',cmat%blk(ib)%H11)
            call io%fld_read(mpl,nam,geom,filename,'H22',cmat%blk(ib)%H22)
            call io%fld_read(mpl,nam,geom,filename,'H33',cmat%blk(ib)%H33)
            call io%fld_read(mpl,nam,geom,filename,'H12',cmat%blk(ib)%H12)
         end if
      end if
      if ((ib==bpar%nbe).and.nam%displ_diag) then
         call io%fld_read(mpl,nam,geom,filename,'displ_lon',cmat%blk(ib)%displ_lon)
         call io%fld_read(mpl,nam,geom,filename,'displ_lat',cmat%blk(ib)%displ_lat)
      end if
   end if
end do

end subroutine cmat_read

!----------------------------------------------------------------------
! Subroutine: cmat_write
! Purpose: write C matrix
!----------------------------------------------------------------------
subroutine cmat_write(cmat,mpl,nam,geom,bpar,io)

implicit none

! Passed variables
class(cmat_type),intent(in) :: cmat ! C matrix
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(bpar_type),intent(in) :: bpar  ! Block parameters
type(io_type),intent(in) :: io      ! I/O

! Local variables
integer :: ib,ncid
character(len=1024) :: filename
character(len=1024) :: subr = 'cmat_write'

do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      ! Set filename
      filename = trim(nam%prefix)//'_'//trim(cmat%blk(ib)%name)

      ! Write fields
      if (bpar%nicas_block(ib)) then
         call io%fld_write(mpl,nam,geom,filename,'coef_ens',cmat%blk(ib)%coef_ens)
         call io%fld_write(mpl,nam,geom,filename,'coef_sta',cmat%blk(ib)%coef_sta)
         call io%fld_write(mpl,nam,geom,filename,'rh',cmat%blk(ib)%rh)
         call io%fld_write(mpl,nam,geom,filename,'rv',cmat%blk(ib)%rv)
         if (cmat%blk(ib)%double_fit) then
            call io%fld_write(mpl,nam,geom,filename,'rv_rfac',cmat%blk(ib)%rv_rfac)
            call io%fld_write(mpl,nam,geom,filename,'rv_coef',cmat%blk(ib)%rv_coef)
         end if
         call io%fld_write(mpl,nam,geom,filename,'rhs',cmat%blk(ib)%rhs)
         call io%fld_write(mpl,nam,geom,filename,'rvs',cmat%blk(ib)%rvs)
         if (cmat%blk(ib)%anisotropic) then
            call io%fld_write(mpl,nam,geom,filename,'H11',cmat%blk(ib)%H11)
            call io%fld_write(mpl,nam,geom,filename,'H22',cmat%blk(ib)%H22)
            call io%fld_write(mpl,nam,geom,filename,'H33',cmat%blk(ib)%H33)
            call io%fld_write(mpl,nam,geom,filename,'H12',cmat%blk(ib)%H12)
         end if
      end if
      if ((ib==bpar%nbe).and.nam%displ_diag) then
         call io%fld_write(mpl,nam,geom,filename,'displ_lon',cmat%blk(ib)%displ_lon)
         call io%fld_write(mpl,nam,geom,filename,'displ_lat',cmat%blk(ib)%displ_lat)
      end if

      if (mpl%main) then
         ! Write attributes
         call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_write,ncid))
         call mpl%ncerr(subr,nf90_redef(ncid))
         if (cmat%blk(ib)%double_fit) then
            call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'double_fit',1))
         else
            call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'double_fit',0))
         end if
         if (cmat%blk(ib)%anisotropic) then
            call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'anisotropic',1))
         else
            call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'anisotropic',0))
         end if
         call mpl%ncerr(subr,nf90_close(ncid))
      end if
   end if
end do

end subroutine cmat_write

!----------------------------------------------------------------------
! Subroutine: cmat_from_hdiag
! Purpose: transform HDIAG into C matrix
!----------------------------------------------------------------------
subroutine cmat_from_hdiag(cmat,mpl,nam,geom,bpar,hdiag)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat ! C matrix
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(hdiag_type),intent(in) :: hdiag   ! Hybrid diagnostics

! Local variables
integer :: ib,n,i,il0,il0i,ic2a,its
real(kind_real) :: fld_c2a(hdiag%samp%nc2a,geom%nl0,6),fld_c2b(hdiag%samp%nc2b,geom%nl0),fld_c0a(geom%nc0a,geom%nl0,6)

! Allocation
call cmat%alloc(bpar,'cmat')

! Copy attributes
do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      select case (trim(nam%method))
      case ('cor')
         cmat%blk(ib)%double_fit = hdiag%cor_1%blk(0,ib)%double_fit
      case ('loc_norm','loc')
         cmat%blk(ib)%double_fit = hdiag%loc_1%blk(0,ib)%double_fit
      case ('hyb-avg','hyb-rnd')
         cmat%blk(ib)%double_fit = hdiag%loc_2%blk(0,ib)%double_fit
      case ('dual-ens')
         call mpl%abort('dual-ens not ready yet for C matrix')
      case default
         call mpl%abort('cmat not implemented yet for this method')
      end select
      cmat%blk(ib)%anisotropic = .false.
   end if
end do

! Allocation
call cmat%alloc(nam,geom,bpar)

! Convolution parameters
do ib=1,bpar%nbe
   if (bpar%B_block(ib)) then
      if (bpar%nicas_block(ib)) then
         if (nam%local_diag) then
            ! Copy data
            n = 4
            if (cmat%blk(ib)%double_fit) n = n+2
            do ic2a=1,hdiag%samp%nc2a
               select case (trim(nam%method))
               case ('cor')
                  fld_c2a(ic2a,:,1) = hdiag%cor_1%blk(ic2a,ib)%raw_coef_ens
                  fld_c2a(ic2a,:,2) = 0.0
                  fld_c2a(ic2a,:,3) = hdiag%cor_1%blk(ic2a,ib)%fit_rh
                  fld_c2a(ic2a,:,4) = hdiag%cor_1%blk(ic2a,ib)%fit_rv
                  if (cmat%blk(ib)%double_fit) then
                     fld_c2a(ic2a,:,5) = hdiag%cor_1%blk(ic2a,ib)%fit_rv_rfac
                     fld_c2a(ic2a,:,6) = hdiag%cor_1%blk(ic2a,ib)%fit_rv_coef
                  end if
               case ('loc_norm','loc')
                  fld_c2a(ic2a,:,1) = hdiag%loc_1%blk(ic2a,ib)%raw_coef_ens
                  fld_c2a(ic2a,:,2) = 0.0
                  fld_c2a(ic2a,:,3) = hdiag%loc_1%blk(ic2a,ib)%fit_rh
                  fld_c2a(ic2a,:,4) = hdiag%loc_1%blk(ic2a,ib)%fit_rv
                  if (cmat%blk(ib)%double_fit) then
                     fld_c2a(ic2a,:,5) = hdiag%loc_1%blk(ic2a,ib)%fit_rv_rfac
                     fld_c2a(ic2a,:,6) = hdiag%loc_1%blk(ic2a,ib)%fit_rv_coef
                  end if
               case ('hyb-avg','hyb-rnd')
                  fld_c2a(ic2a,:,1) = hdiag%loc_2%blk(ic2a,ib)%raw_coef_ens
                  fld_c2a(ic2a,:,2) = hdiag%loc_2%blk(ic2a,ib)%raw_coef_sta
                  fld_c2a(ic2a,:,3) = hdiag%loc_2%blk(ic2a,ib)%fit_rh
                  fld_c2a(ic2a,:,4) = hdiag%loc_2%blk(ic2a,ib)%fit_rv
                  if (cmat%blk(ib)%double_fit) then
                     fld_c2a(ic2a,:,5) = hdiag%loc_2%blk(ic2a,ib)%fit_rv_rfac
                     fld_c2a(ic2a,:,6) = hdiag%loc_2%blk(ic2a,ib)%fit_rv_coef
                  end if
               case ('dual-ens')
                  call mpl%abort('dual-ens not ready yet for C matrix')
               case default
                  call mpl%abort('cmat not implemented yet for this method')
               end select
            end do

            do i=1,n
               ! Fill missing values
               do il0=1,geom%nl0
                  call hdiag%samp%diag_fill(mpl,nam,geom,il0,fld_c2a(:,il0,i))
               end do

               ! Interpolate
               call hdiag%samp%com_AB%ext(mpl,geom%nl0,fld_c2a(:,:,i),fld_c2b)
               do il0=1,geom%nl0
                  il0i = min(il0,geom%nl0i)
                  call hdiag%samp%h(il0i)%apply(mpl,fld_c2b(:,il0),fld_c0a(:,il0,i))
               end do
            end do

            ! Copy to C matrix
            cmat%blk(ib)%coef_ens = fld_c0a(:,:,1)
            call mpl%f_comm%allreduce(sum(cmat%blk(ib)%coef_ens,mask=geom%mask_c0a),cmat%blk(ib)%wgt,fckit_mpi_sum())
            cmat%blk(ib)%wgt = cmat%blk(ib)%wgt/real(count(geom%mask_c0),kind_real)
            cmat%blk(ib)%coef_sta = fld_c0a(:,:,2)
            cmat%blk(ib)%rh = fld_c0a(:,:,3)
            cmat%blk(ib)%rv = fld_c0a(:,:,4)
            if (cmat%blk(ib)%double_fit) then
               cmat%blk(ib)%rv_rfac = fld_c0a(:,:,5)
               cmat%blk(ib)%rv_coef = fld_c0a(:,:,6)
            end if
         else
            ! Copy to C matrix
            do il0=1,geom%nl0
               ! Copy data
               select case (trim(nam%method))
               case ('cor')
                  cmat%blk(ib)%coef_ens(:,il0) = hdiag%cor_1%blk(0,ib)%raw_coef_ens(il0)
                  cmat%blk(ib)%wgt = sum(hdiag%cor_1%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
                  cmat%blk(ib)%coef_sta(:,il0) = 0.0
                  cmat%blk(ib)%rh(:,il0) = hdiag%cor_1%blk(0,ib)%fit_rh(il0)
                  cmat%blk(ib)%rv(:,il0) = hdiag%cor_1%blk(0,ib)%fit_rv(il0)
                  if (cmat%blk(ib)%double_fit) then
                     cmat%blk(ib)%rv_rfac(:,il0) = hdiag%cor_1%blk(0,ib)%fit_rv_rfac(il0)
                     cmat%blk(ib)%rv_coef(:,il0) = hdiag%cor_1%blk(0,ib)%fit_rv_coef(il0)
                  end if
               case ('loc_norm','loc')
                  cmat%blk(ib)%coef_ens(:,il0) = hdiag%loc_1%blk(0,ib)%raw_coef_ens(il0)
                  cmat%blk(ib)%wgt = sum(hdiag%loc_1%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
                  cmat%blk(ib)%coef_sta(:,il0) = 0.0
                  cmat%blk(ib)%rh(:,il0) = hdiag%loc_1%blk(0,ib)%fit_rh(il0)
                  cmat%blk(ib)%rv(:,il0) = hdiag%loc_1%blk(0,ib)%fit_rv(il0)
                  if (cmat%blk(ib)%double_fit) then
                     cmat%blk(ib)%rv_rfac(:,il0) = hdiag%loc_1%blk(0,ib)%fit_rv_rfac(il0)
                     cmat%blk(ib)%rv_coef(:,il0) = hdiag%loc_1%blk(0,ib)%fit_rv_coef(il0)
                  end if
               case ('hyb-avg','hyb-rnd')
                  cmat%blk(ib)%coef_ens(:,il0) = hdiag%loc_2%blk(0,ib)%raw_coef_ens(il0)
                  cmat%blk(ib)%wgt = sum(hdiag%loc_2%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
                  cmat%blk(ib)%coef_sta(:,il0) = hdiag%loc_2%blk(0,ib)%raw_coef_sta
                  cmat%blk(ib)%rh(:,il0) = hdiag%loc_2%blk(0,ib)%fit_rh(il0)
                  cmat%blk(ib)%rv(:,il0) = hdiag%loc_2%blk(0,ib)%fit_rv(il0)
                  if (cmat%blk(ib)%double_fit) then
                     cmat%blk(ib)%rv_rfac(:,il0) = hdiag%loc_2%blk(0,ib)%fit_rv_rfac(il0)
                     cmat%blk(ib)%rv_coef(:,il0) = hdiag%loc_2%blk(0,ib)%fit_rv_coef(il0)
                  end if
               case ('dual-ens')
                  call mpl%abort('dual-ens not ready yet for C matrix')
               case default
                  call mpl%abort('cmat not implemented yet for this method')
               end select
            end do
         end if
      else
         ! Define weight only
         select case (trim(nam%method))
         case ('cor')
            cmat%blk(ib)%wgt = sum(hdiag%cor_1%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
         case ('loc_norm','loc')
            cmat%blk(ib)%wgt = sum(hdiag%loc_1%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
         case ('hyb-avg','hyb-rnd')
            cmat%blk(ib)%wgt = sum(hdiag%loc_2%blk(0,ib)%raw_coef_ens)/real(geom%nl0,kind_real)
         case ('dual-ens')
            call mpl%abort('dual-ens not ready yet for C matrix')
         case default
            call mpl%abort('cmat not implemented yet for this method')
         end select
      end if
   end if
end do

! Displacement
if (nam%displ_diag) then
   do its=2,nam%nts
      cmat%blk(bpar%nbe)%displ_lon(:,:,its) = hdiag%samp%displ_lon(:,:,its)
      cmat%blk(bpar%nbe)%displ_lat(:,:,its) = hdiag%samp%displ_lat(:,:,its)
   end do
end if

end subroutine cmat_from_hdiag

!----------------------------------------------------------------------
! Subroutine: cmat_from_lct
! Purpose: copy LCT into C matrix
!----------------------------------------------------------------------
subroutine cmat_from_lct(cmat,mpl,nam,geom,bpar,lct)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat ! C matrix
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(lct_type),intent(in) :: lct       ! LCT

! Local variables
integer :: ib,iv,jv,its,jts,iscales,il0,ic0a
real(kind_real) :: tr,det,diff

! Allocation
call cmat%alloc(bpar,'cmat')

! Copy attributes
do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      cmat%blk(ib)%double_fit = .false.
      cmat%blk(ib)%anisotropic = .true.
   end if
end do

! Allocation
call cmat%alloc(nam,geom,bpar)

! Convolution parameters
do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      ! Indices
      iv = bpar%b_to_v1(ib)
      jv = bpar%b_to_v2(ib)
      its = bpar%b_to_ts1(ib)
      jts = bpar%b_to_ts2(ib)
      if ((iv/=jv).or.(its/=jts)) call mpl%abort('only diagonal blocks for cmat_from_lct')

      if (lct%blk(ib)%nscales>1) call mpl%warning('only the first scale is used to define cmat from LCT')
      iscales = 1

      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%mask_c0a(ic0a,il0)) then
               ! Copy LCT
               cmat%blk(ib)%H11(ic0a,il0) = lct%blk(ib)%H11(ic0a,il0,iscales)
               cmat%blk(ib)%H22(ic0a,il0) = lct%blk(ib)%H22(ic0a,il0,iscales)
               cmat%blk(ib)%H33(ic0a,il0) = lct%blk(ib)%H33(ic0a,il0,iscales)
               cmat%blk(ib)%H12(ic0a,il0) = lct%blk(ib)%H12(ic0a,il0,iscales)

               ! Copy scale coefficient
               cmat%blk(ib)%Hcoef(ic0a,il0) = lct%blk(ib)%Dcoef(ic0a,il0,iscales)

               ! Compute support radii from the largest scale
               tr = cmat%blk(ib)%H11(ic0a,il0)+cmat%blk(ib)%H22(ic0a,il0)
               det = cmat%blk(ib)%H11(ic0a,il0)*cmat%blk(ib)%H22(ic0a,il0)-cmat%blk(ib)%H12(ic0a,il0)**2
               diff = 0.25*(cmat%blk(ib)%H11(ic0a,il0)-cmat%blk(ib)%H22(ic0a,il0))**2+cmat%blk(ib)%H12(ic0a,il0)**2
               if ((det>0.0).and..not.(diff<0.0)) then
                  if (0.5*tr>sqrt(diff)) then
                     cmat%blk(ib)%rh(ic0a,il0) = gau2gc/(sqrt(0.5*tr-sqrt(diff)))
                  else
                     write(mpl%info,*) 0.5*tr,sqrt(diff)
                     call mpl%abort('non positive-definite LCT in cmat_from_lct (eigenvalue)')
                  end if
               else
                  write(mpl%info,*) cmat%blk(ib)%H11(ic0a,il0),cmat%blk(ib)%H22(ic0a,il0), &
 & cmat%blk(ib)%H12(ic0a,il0)
                  call mpl%abort('non positive-definite LCT in cmat_from_lct (determinant)')
               end if
               if (cmat%blk(ib)%H33(ic0a,il0)>0.0) then
                  cmat%blk(ib)%rv(ic0a,il0) = gau2gc/(sqrt(cmat%blk(ib)%H33(ic0a,il0)))
               else
                  cmat%blk(ib)%rv(ic0a,il0) = 0.0
               end if
            end if
         end do
      end do

      ! Set coefficients
      cmat%blk(ib)%coef_ens = 1.0
      cmat%blk(ib)%coef_sta = 0.0
      cmat%blk(ib)%wgt = 1.0
   end if
end do

end subroutine cmat_from_lct

!----------------------------------------------------------------------
! Subroutine: cmat_from_nam
! Purpose: copy radii into C matrix
!----------------------------------------------------------------------
subroutine cmat_from_nam(cmat,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat ! C matrix
type(mpl_type),intent(in) :: mpl       ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters

! Local variables
integer :: ib,iv,jv,its,jts

write(mpl%info,'(a)') '-------------------------------------------------------------------'
write(mpl%info,'(a)') '--- Copy namelist radii into C matrix'
call flush(mpl%info)

! Allocation
call cmat%alloc(bpar,'cmat')

! Copy attributes
do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      cmat%blk(ib)%double_fit = .false.
      cmat%blk(ib)%anisotropic = .false.
   end if
end do

! Allocation
call cmat%alloc(nam,geom,bpar)

! Convolution parameters
do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      ! Indices
      iv = bpar%b_to_v1(ib)
      jv = bpar%b_to_v2(ib)
      its = bpar%b_to_ts1(ib)
      jts = bpar%b_to_ts2(ib)
      if ((iv/=jv).or.(its/=jts)) call mpl%abort('only diagonal blocks for cmat_from_nam')

      ! Copy support radii
      cmat%blk(ib)%rh = nam%rh
      cmat%blk(ib)%rhs = nam%rh
      cmat%blk(ib)%rv = nam%rv
      cmat%blk(ib)%rvs = nam%rv

      ! Set coefficients
      cmat%blk(ib)%coef_ens = 1.0
      cmat%blk(ib)%coef_sta = 0.0
      cmat%blk(ib)%wgt = 1.0
   end if
end do

end subroutine cmat_from_nam

!----------------------------------------------------------------------
! Subroutine: cmat_from_oops
! Purpose: copy C matrix from OOPS
!----------------------------------------------------------------------
subroutine cmat_from_oops(cmat,mpl,geom,bpar)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat ! C matrix
type(mpl_type),intent(in) :: mpl       ! MPI data
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters

! Local variables
integer :: ib

do ib=1,bpar%nbe
   if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
      if (allocated(cmat%blk(ib)%oops_coef_ens)) then
         write(mpl%info,'(a7,a,a)') '','Ensemble coefficient copied from OOPS for block ',trim(bpar%blockname(ib))
         cmat%blk(ib)%coef_ens = cmat%blk(ib)%oops_coef_ens
         call mpl%f_comm%allreduce(sum(cmat%blk(ib)%coef_ens,mask=geom%mask_c0a),cmat%blk(ib)%wgt,fckit_mpi_sum())
         cmat%blk(ib)%wgt = cmat%blk(ib)%wgt/real(count(geom%mask_c0),kind_real)
      end if
      if (allocated(cmat%blk(ib)%oops_coef_sta)) then
         write(mpl%info,'(a7,a,a)') '','Static coefficient copied from OOPS for block ',trim(bpar%blockname(ib))
         cmat%blk(ib)%coef_sta = cmat%blk(ib)%oops_coef_sta
      end if
      if (allocated(cmat%blk(ib)%oops_rh)) then
         write(mpl%info,'(a7,a,a)') '','Horizontal fit support radius copied from OOPS for block ',trim(bpar%blockname(ib))
         cmat%blk(ib)%rh = cmat%blk(ib)%oops_rh
      end if
      if (allocated(cmat%blk(ib)%oops_rv)) then
         write(mpl%info,'(a7,a,a)') '','Vertical fit support radius copied from OOPS for block ',trim(bpar%blockname(ib))
         cmat%blk(ib)%rv = cmat%blk(ib)%oops_rv
      end if
      if (allocated(cmat%blk(ib)%oops_rv_rfac)) then
         write(mpl%info,'(a7,a,a)') '','Vertical fit factor copied from OOPS for block ',trim(bpar%blockname(ib))
         cmat%blk(ib)%rv_rfac = cmat%blk(ib)%oops_rv_rfac
      end if
      if (allocated(cmat%blk(ib)%oops_rv_coef)) then
         write(mpl%info,'(a7,a,a)') '','Vertical fit coefficient copied from OOPS for block ',trim(bpar%blockname(ib))
         cmat%blk(ib)%rv_coef = cmat%blk(ib)%oops_rv_coef
      end if
   end if
end do

end subroutine cmat_from_oops

!----------------------------------------------------------------------
! Subroutine: cmat_setup_sampling
! Purpose: setup C matrix sampling
!----------------------------------------------------------------------
subroutine cmat_setup_sampling(cmat,nam,geom,bpar)

implicit none

! Passed variables
class(cmat_type),intent(inout) :: cmat    ! C matrix
type(nam_type),intent(in) :: nam          ! Namelist
type(geom_type),intent(in) :: geom        ! Geometry
type(bpar_type),intent(in) :: bpar        ! Block parameters

! Local variables
integer :: ib,il0,ic0a

! Sampling parameters
if (trim(nam%strategy)=='specific_multivariate') then
   ! Initialization
   cmat%blk(bpar%nbe)%rhs = huge(1.0)
   cmat%blk(bpar%nbe)%rvs = huge(1.0)
   do ib=1,bpar%nb
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         ! Get minimum
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               cmat%blk(bpar%nbe)%rhs(ic0a,il0) = min(cmat%blk(bpar%nbe)%rhs(ic0a,il0),cmat%blk(ib)%rh(ic0a,il0))
               cmat%blk(bpar%nbe)%rvs(ic0a,il0) = min(cmat%blk(bpar%nbe)%rvs(ic0a,il0),cmat%blk(ib)%rv(ic0a,il0))
            end do
         end do
      end if
   end do
else
   ! Copy
   do ib=1,bpar%nbe
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         cmat%blk(ib)%rhs = cmat%blk(ib)%rh
         cmat%blk(ib)%rvs = cmat%blk(ib)%rv
      end if
   end do
end if

end subroutine cmat_setup_sampling

end module type_cmat
