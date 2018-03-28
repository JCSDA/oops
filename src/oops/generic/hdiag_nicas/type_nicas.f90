!----------------------------------------------------------------------
! Module: type_nicas
!> Purpose: NICAS data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_nicas

use driver_hdiag, only: run_hdiag
use model_interface, only: model_write
use netcdf
use tools_const, only: rad2deg,reqkm
use tools_display, only: msgerror,vunitchar
use tools_func, only: cholesky
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr
use tools_nc, only: ncerr,ncfloat
use tools_test, only: define_dirac,define_test_vectors
use type_bpar, only: bpar_type
use type_cmat, only: cmat_type
use type_com, only: com_type
use type_cv, only: cv_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_nicas_blk, only: nicas_blk_type
use type_mpl, only: mpl
use type_nam, only: nam_type
use type_rng, only: rng
use yomhook, only: lhook,dr_hook

implicit none

! NICAS derived type
type nicas_type
   character(len=1024) :: prefix              !< Prefix
   type(nicas_blk_type),allocatable :: blk(:) !< NICAS data blocks
contains
   procedure :: alloc => nicas_alloc
   procedure :: read => nicas_read
   procedure :: write => nicas_write
   procedure :: write_mpi_summary => nicas_write_mpi_summary
   procedure :: alloc_cv => nicas_alloc_cv
   procedure :: random_cv => nicas_random_cv
   procedure :: apply => nicas_apply
   procedure :: apply_from_sqrt => nicas_apply_from_sqrt
   procedure :: apply_sqrt => nicas_apply_sqrt
   procedure :: apply_sqrt_ad => nicas_apply_sqrt_ad
   procedure :: randomize => nicas_randomize
   procedure :: apply_bens => nicas_apply_bens
   procedure :: apply_bens_noloc => nicas_apply_bens_noloc
   procedure :: test_adjoint => nicas_test_adjoint
   procedure :: test_sqrt => nicas_test_sqrt
   procedure :: test_dirac => nicas_test_dirac
   procedure :: test_randomization => nicas_test_randomization
   procedure :: test_consistency => nicas_test_consistency
   procedure :: test_optimality => nicas_test_optimality
end type nicas_type

logical :: lcoef_ens = .false.     !< Apply ensemble coefficient (will reduce variance)
integer,parameter :: ne_rand = 150 !< Ensemble size for randomization
integer,parameter :: nfac = 10     !< Number of length-scale factors
integer,parameter :: ntest = 100   !< Number of tests

private
public :: nicas_blk_type,nicas_type

contains

!----------------------------------------------------------------------
! Subroutine: nicas_alloc
!> Purpose: nicas object allocation
!----------------------------------------------------------------------
subroutine nicas_alloc(nicas,nam,bpar,prefix)

! Passed variables
class(nicas_type),intent(out) :: nicas !< NICAS data
type(nam_type),intent(in) :: nam       !< Namelist
type(bpar_type),intent(in) :: bpar     !< Block parameters
character(len=*),intent(in) :: prefix  !< Prefix

! Local variable
integer :: ib

! Copy prefix
nicas%prefix = prefix

! Allocation
allocate(nicas%blk(bpar%nb+1))

! Set name
do ib=1,bpar%nb+1
   nicas%blk(ib)%ib = ib
   write(nicas%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(prefix)//'_',nam%mpicom,'_',mpl%nproc,'-',mpl%myproc, &
 & '_',trim(bpar%blockname(ib))
end do

end subroutine nicas_alloc

!----------------------------------------------------------------------
! Subroutine: nicas_read
!> Purpose: read nicas object
!----------------------------------------------------------------------
subroutine nicas_read(nicas,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas !< NICAS data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
type(bpar_type),intent(in) :: bpar       !< Block parameters

! Local variables
integer :: ib,il0i,il1,its,il0
integer :: info,nl0_test,nc0a_test
integer :: ncid,nl0_id,nc0a_id,nc1b_id,nl1_id,nsa_id,nsb_id,nc0d_id,nc0dinv_id
integer :: sb_to_c1b_id,sb_to_l1_id
integer :: sa_to_sc_id,sb_to_sc_id
integer :: norm_id,coef_ens_id
character(len=1024) :: filename
character(len=1024) :: subr = 'nicas_read'

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      ! Open file and get dimensions
      filename = trim(nam%prefix)//'_'//trim(nicas%blk(ib)%name)//'.nc'
      call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))

      ! Get dimensions
      call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
      call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=nl0_test))
      if (nl0_test/=geom%nl0) call msgerror('wrong dimension when reading nicas')
      if (bpar%nicas_block(ib)) then
         call ncerr(subr,nf90_inq_dimid(ncid,'nc0a',nc0a_id))
         call ncerr(subr,nf90_inquire_dimension(ncid,nc0a_id,len=nc0a_test))
         if (nc0a_test/=geom%nc0a) call msgerror('wrong dimension when reading nicas')
         info = nf90_inq_dimid(ncid,'nc1b',nc1b_id)
         if (info==nf90_noerr) then
            call ncerr(subr,nf90_inquire_dimension(ncid,nc1b_id,len=nicas%blk(ib)%nc1b))
         else
            nicas%blk(ib)%nc1b = 0
         end if
         call ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
         call ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=nicas%blk(ib)%nl1))
         info = nf90_inq_dimid(ncid,'nsa',nsa_id)
         if (info==nf90_noerr) then
            call ncerr(subr,nf90_inquire_dimension(ncid,nsa_id,len=nicas%blk(ib)%nsa))
         else
            nicas%blk(ib)%nsa = 0
         end if
         info = nf90_inq_dimid(ncid,'nsb',nsb_id)
         if (info==nf90_noerr) then
            call ncerr(subr,nf90_inquire_dimension(ncid,nsb_id,len=nicas%blk(ib)%nsb))
         else
            nicas%blk(ib)%nsb = 0
         end if
         call ncerr(subr,nf90_get_att(ncid,nf90_global,'nsc',nicas%blk(ib)%nsc))
         call ncerr(subr,nf90_get_att(ncid,nf90_global,'mpicom',nicas%blk(ib)%mpicom))
      end if
      if ((ib==bpar%nb+1).and.(abs(nam%advmode)==1)) then
         call ncerr(subr,nf90_inq_dimid(ncid,'nc0d',nc0d_id))
         call ncerr(subr,nf90_inquire_dimension(ncid,nc0d_id,len=nicas%blk(ib)%nc0d))
         call ncerr(subr,nf90_inq_dimid(ncid,'nc0dinv',nc0dinv_id))
         call ncerr(subr,nf90_inquire_dimension(ncid,nc0dinv_id,len=nicas%blk(ib)%nc0dinv))
      end if

      ! Allocation
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) allocate(nicas%blk(ib)%sb_to_c1b(nicas%blk(ib)%nsb))
         if (nicas%blk(ib)%nsb>0) allocate(nicas%blk(ib)%sb_to_l1(nicas%blk(ib)%nsb))
         if (nicas%blk(ib)%nsa>0) allocate(nicas%blk(ib)%sa_to_sc(nicas%blk(ib)%nsa))
         if (nicas%blk(ib)%nsb>0) allocate(nicas%blk(ib)%sb_to_sc(nicas%blk(ib)%nsb))
         allocate(nicas%blk(ib)%norm(geom%nc0a,geom%nl0))
         allocate(nicas%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
         allocate(nicas%blk(ib)%h(geom%nl0i))
         allocate(nicas%blk(ib)%s(nicas%blk(ib)%nl1))
      end if
      if ((ib==bpar%nb+1).and.(abs(nam%advmode)==1)) then
         allocate(nicas%blk(ib)%d(geom%nl0,2:nam%nts))
         allocate(nicas%blk(ib)%dinv(geom%nl0,2:nam%nts))
      end if

      ! Get variable id
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_c1b',sb_to_c1b_id))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_l1',sb_to_l1_id))
         if (nicas%blk(ib)%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'sa_to_sc',sa_to_sc_id))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_sc',sb_to_sc_id))
         call ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
         call ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))
      end if

      ! Read data
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_get_var(ncid,sb_to_c1b_id,nicas%blk(ib)%sb_to_c1b))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_get_var(ncid,sb_to_l1_id,nicas%blk(ib)%sb_to_l1))
         if (nicas%blk(ib)%nsa>0) call ncerr(subr,nf90_get_var(ncid,sa_to_sc_id,nicas%blk(ib)%sa_to_sc))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_get_var(ncid,sb_to_sc_id,nicas%blk(ib)%sb_to_sc))
         call ncerr(subr,nf90_get_var(ncid,norm_id,nicas%blk(ib)%norm))
         call ncerr(subr,nf90_get_var(ncid,coef_ens_id,nicas%blk(ib)%coef_ens))
         call nicas%blk(ib)%com_AB%read(ncid,'com_AB')
         call nicas%blk(ib)%com_AC%read(ncid,'com_AC')
         nicas%blk(ib)%c%prefix = 'c'
         call nicas%blk(ib)%c%read(ncid)
         do il0i=1,geom%nl0i
            write(nicas%blk(ib)%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
            call nicas%blk(ib)%h(il0i)%read(ncid)
         end do
         nicas%blk(ib)%v%prefix = 'v'
         call nicas%blk(ib)%v%read(ncid)
         do il1=1,nicas%blk(ib)%nl1
            write(nicas%blk(ib)%s(il1)%prefix,'(a,i3.3)') 's_',il1
            call nicas%blk(ib)%s(il1)%read(ncid)
         end do
      end if
      if ((ib==bpar%nb+1).and.(abs(nam%advmode)==1)) then
         call nicas%blk(ib)%com_AD%read(ncid,'com_AD')
         call nicas%blk(ib)%com_ADinv%read(ncid,'com_ADinv')
         do its=2,nam%nts
            do il0=1,geom%nl0
               write(nicas%blk(ib)%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
               call nicas%blk(ib)%d(il0,its)%read(ncid)
               write(nicas%blk(ib)%dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'dinv_',il0,'_',its
               call nicas%blk(ib)%dinv(il0,its)%read(ncid)
            end do
         end do
      end if

      ! Read main weight
      call ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',nicas%blk(ib)%wgt))

      ! Close file
      call ncerr(subr,nf90_close(ncid))
   end if
end do

end subroutine nicas_read

!----------------------------------------------------------------------
! Subroutine: nicas_write
!> Purpose: write nicas object
!----------------------------------------------------------------------
subroutine nicas_write(nicas,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters

! Local variables
integer :: ib,il0i,il1,its,il0
integer :: ncid,nl0_id,nc0a_id,nc1b_id,nl1_id,nsa_id,nsb_id,nc0d_id,nc0dinv_id
integer :: sb_to_c1b_id,sb_to_l1_id
integer :: sa_to_sc_id,sb_to_sc_id
integer :: norm_id,coef_ens_id
character(len=1024) :: filename
character(len=1024) :: subr = 'nicas_write'

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      ! Create file
      filename = trim(nam%prefix)//'_'//trim(nicas%blk(ib)%name)//'.nc'
      call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

      ! Write namelist parameters
      call nam%ncwrite(ncid)

      ! Define dimensions
      call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
      if (bpar%nicas_block(ib)) then
         call ncerr(subr,nf90_def_dim(ncid,'nc0a',geom%nc0a,nc0a_id))
         call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',geom%nl0i))
         if (nicas%blk(ib)%nc1b>0) call ncerr(subr,nf90_def_dim(ncid,'nc1b',nicas%blk(ib)%nc1b,nc1b_id))
         call ncerr(subr,nf90_def_dim(ncid,'nl1',nicas%blk(ib)%nl1,nl1_id))
         if (nicas%blk(ib)%nsa>0) call ncerr(subr,nf90_def_dim(ncid,'nsa',nicas%blk(ib)%nsa,nsa_id))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_def_dim(ncid,'nsb',nicas%blk(ib)%nsb,nsb_id))
         call ncerr(subr,nf90_put_att(ncid,nf90_global,'nsc',nicas%blk(ib)%nsc))
         call ncerr(subr,nf90_put_att(ncid,nf90_global,'mpicom',nicas%blk(ib)%mpicom))
      end if
      if ((ib==bpar%nb+1).and.nam%displ_diag) then
         call ncerr(subr,nf90_def_dim(ncid,'nc0d',nicas%blk(ib)%nc0d,nc0d_id))
         call ncerr(subr,nf90_def_dim(ncid,'nc0dinv',nicas%blk(ib)%nc0dinv,nc0dinv_id))
      end if

      ! Write main weight
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',nicas%blk(ib)%wgt))

      ! Define variables
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'sb_to_c1b',nf90_int,(/nsb_id/),sb_to_c1b_id))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'sb_to_l1',nf90_int,(/nsb_id/),sb_to_l1_id))
         if (nicas%blk(ib)%nsa>0) call ncerr(subr,nf90_def_var(ncid,'sa_to_sc',nf90_int,(/nsa_id/),sa_to_sc_id))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_def_var(ncid,'sb_to_sc',nf90_int,(/nsb_id/),sb_to_sc_id))
         call ncerr(subr,nf90_def_var(ncid,'norm',ncfloat,(/nc0a_id,nl0_id/),norm_id))
         call ncerr(subr,nf90_def_var(ncid,'coef_ens',ncfloat,(/nc0a_id,nl0_id/),coef_ens_id))

         if (nicas%blk(ib)%nsa>0) call ncerr(subr,nf90_put_att(ncid,sa_to_sc_id,'_FillValue',msvali))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_put_att(ncid,sb_to_sc_id,'_FillValue',msvali))
         call ncerr(subr,nf90_put_att(ncid,norm_id,'_FillValue',msvalr))
         call ncerr(subr,nf90_put_att(ncid,coef_ens_id,'_FillValue',msvalr))
      end if

      ! End definition mode
      call ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_put_var(ncid,sb_to_c1b_id,nicas%blk(ib)%sb_to_c1b))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_put_var(ncid,sb_to_l1_id,nicas%blk(ib)%sb_to_l1))
         if (nicas%blk(ib)%nsa>0) call ncerr(subr,nf90_put_var(ncid,sa_to_sc_id,nicas%blk(ib)%sa_to_sc))
         if (nicas%blk(ib)%nsb>0) call ncerr(subr,nf90_put_var(ncid,sb_to_sc_id,nicas%blk(ib)%sb_to_sc))
         call ncerr(subr,nf90_put_var(ncid,norm_id,nicas%blk(ib)%norm))
         call ncerr(subr,nf90_put_var(ncid,coef_ens_id,nicas%blk(ib)%coef_ens))
         call nicas%blk(ib)%com_AB%write(ncid)
         call nicas%blk(ib)%com_AC%write(ncid)
         call nicas%blk(ib)%c%write(ncid)
         do il0i=1,geom%nl0i
            call nicas%blk(ib)%h(il0i)%write(ncid)
         end do
         call nicas%blk(ib)%v%write(ncid)
         do il1=1,nicas%blk(ib)%nl1
            call nicas%blk(ib)%s(il1)%write(ncid)
         end do
      end if
      if ((ib==bpar%nb+1).and.nam%displ_diag) then
         call nicas%blk(ib)%com_AD%write(ncid)
         call nicas%blk(ib)%com_ADinv%write(ncid)
         do its=2,nam%nts
            do il0=1,geom%nl0
               call nicas%blk(ib)%d(il0,its)%write(ncid)
               call nicas%blk(ib)%dinv(il0,its)%write(ncid)
            end do
         end do
      end if

      ! Close file
      call ncerr(subr,nf90_close(ncid))
   end if
end do

end subroutine nicas_write

!----------------------------------------------------------------------
! Subroutine: nicas_write_mpi_summary
!> Purpose: write nicas object
!----------------------------------------------------------------------
subroutine nicas_write_mpi_summary(nicas,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(nam_type),intent(in) :: nam      !< Namelist
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters

! Local variables
integer :: ib,is,ic1,il1
integer :: ncid,nc0_id,nc1_id,nl1_id
integer :: lon_id,lat_id,c0_to_proc_id,c1_to_c0_id,l1_to_l0_id,lcheck_id
real(kind_real),allocatable :: lcheck(:,:)
character(len=1024) :: filename
character(len=1024) :: subr = 'nicas_write_mpi_summary'

do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      ! Allocation
      allocate(lcheck(nicas%blk(ib)%nc1,nicas%blk(ib)%nl1))

      ! Create summary file
      filename = trim(nam%prefix)//'_'//trim(nicas%blk(ib)%name)//'_summary.nc'
      call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

      ! Write namelist parameters
      call nam%ncwrite(ncid)

      ! Define dimensions
      call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
      call ncerr(subr,nf90_def_dim(ncid,'nc1',nicas%blk(ib)%nc1,nc1_id))
      call ncerr(subr,nf90_def_dim(ncid,'nl1',nicas%blk(ib)%nl1,nl1_id))

      ! Define variables
      call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
      call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
      call ncerr(subr,nf90_def_var(ncid,'c0_to_proc',nf90_int,(/nc0_id/),c0_to_proc_id))
      call ncerr(subr,nf90_def_var(ncid,'c1_to_c0',nf90_int,(/nc1_id/),c1_to_c0_id))
      call ncerr(subr,nf90_def_var(ncid,'l1_to_l0',nf90_int,(/nl1_id/),l1_to_l0_id))
      call ncerr(subr,nf90_def_var(ncid,'lcheck',ncfloat,(/nc1_id,nl1_id/),lcheck_id))

      ! End definition mode
      call ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
      call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
      call ncerr(subr,nf90_put_var(ncid,c0_to_proc_id,geom%c0_to_proc))
      call ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,nicas%blk(ib)%c1_to_c0))
      call ncerr(subr,nf90_put_var(ncid,l1_to_l0_id,nicas%blk(ib)%l1_to_l0))
      call msr(lcheck)
      do is=1,nicas%blk(ib)%ns
         ic1 = nicas%blk(ib)%s_to_c1(is)
         il1 = nicas%blk(ib)%s_to_l1(is)
         if (nicas%blk(ib)%lcheck_sa(is)) then
            lcheck(ic1,il1) = 1.0
         elseif (nicas%blk(ib)%lcheck_sb(is)) then
            lcheck(ic1,il1) = 2.0
         elseif (nicas%blk(ib)%lcheck_sc(is)) then
            lcheck(ic1,il1) = 3.0
         else
            lcheck(ic1,il1) = 4.0
         end if
      end do
      call ncerr(subr,nf90_put_var(ncid,lcheck_id,lcheck))

      ! Close summary file
      call ncerr(subr,nf90_close(ncid))

      ! Release memory
      deallocate(lcheck)
   end if
end do

end subroutine nicas_write_mpi_summary

!----------------------------------------------------------------------
! Subroutine: nicas_alloc_cv
!> Purpose: control vector object allocation
!----------------------------------------------------------------------
subroutine nicas_alloc_cv(nicas,bpar,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(cv_type),intent(inout) :: cv     !< Control vector

! Local variables
integer :: ib

do ib=1,bpar%nb+1
   if (bpar%cv_block(ib)) then
      ! Allocation
      if (bpar%nicas_block(ib)) then
         allocate(cv%blk(ib)%alpha(nicas%blk(ib)%nsa))
      else
         allocate(cv%blk(ib)%alpha(nicas%blk(1)%nsa))
      end if

      ! Initialization
      call msr(cv%blk(ib)%alpha)
   end if
end do

end subroutine nicas_alloc_cv

!----------------------------------------------------------------------
! Subroutine: nicas_random_cv
!> Purpose: generate a random control vector
!----------------------------------------------------------------------
subroutine nicas_random_cv(nicas,bpar,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(cv_type),intent(out) :: cv       !< Control vector

! Local variables
integer :: ib

! Allocation
call nicas%alloc_cv(bpar,cv)

! Random initialization
do ib=1,bpar%nb+1
   if (bpar%cv_block(ib)) call rng%rand_gau(cv%blk(ib)%alpha)
end do

end subroutine nicas_random_cv

!----------------------------------------------------------------------
! Subroutine: nicas_apply
!> Purpose: apply localization
!----------------------------------------------------------------------
subroutine nicas_apply(nicas,nam,geom,bpar,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   !< NICAS data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
type(bpar_type),target,intent(in) :: bpar                               !< Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ib,its,iv,jv,il0,ic0a
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:)

if (lhook) call dr_hook('nicas_apply',0,zhook_handle)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Adjoint advection
   if (nam%advmode==1) call nicas%blk(bpar%nb+1)%apply_adv_ad(nam,geom,fld)

   ! Sum product over variables and timeslots
   fld_3d = 0.0
   !$omp parallel do schedule(static) private(il0,ic0a,its,iv)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         do its=1,nam%nts
            do iv=1,nam%nv
               fld_3d(ic0a,il0) = fld_3d(ic0a,il0)+fld(ic0a,il0,iv,its)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   if (lcoef_ens) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nb+1)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
   end if

   ! Apply common localization
   call nicas%blk(bpar%nb+1)%apply(geom,fld_3d)

   if (lcoef_ens) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nb+1)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
   end if

   ! Build final vector
   !$omp parallel do schedule(static) private(il0,ic0a,its,iv)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         do its=1,nam%nts
            do iv=1,nam%nv
               fld(ic0a,il0,iv,its) = fld_3d(ic0a,il0)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   ! Advection
   if (nam%advmode==1) call nicas%blk(bpar%nb+1)%apply_adv(nam,geom,fld)
case ('specific_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call nicas%blk(ib)%apply(geom,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('specific_multivariate')
   call msgerror('specific multivariate strategy should not be called from apply_localization (lsqrt required)')
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))

   ! Copy weights
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         if (isnotmsr(nicas%blk(ib)%wgt)) then
            wgt(iv,jv) = nicas%blk(ib)%wgt
         else
            wgt(iv,jv) = 0.0
         end if
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do iv=1,nam%nv
      ! Apply common ensemble coefficient square-root
      if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(bpar%nb+1)%coef_ens)

      ! Apply common localization
      call nicas%blk(bpar%nb+1)%apply(geom,fld_4d(:,:,iv))

      ! Apply common ensemble coefficient square-root
      if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(bpar%nb+1)%coef_ens)
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=1,nam%nv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d_tmp
   end do
end select

if (lhook) call dr_hook('nicas_apply',1,zhook_handle)

end subroutine nicas_apply

!----------------------------------------------------------------------
! Subroutine: nicas_apply_from_sqrt
!> Purpose: apply localization from square-root
!----------------------------------------------------------------------
subroutine nicas_apply_from_sqrt(nicas,nam,geom,bpar,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   !< NICAS data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
type(bpar_type),target,intent(in) :: bpar                               !< Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
real(kind_real) :: zhook_handle
type(cv_type) :: cv

if (lhook) call dr_hook('nicas_apply_from_sqrt',0,zhook_handle)

! Apply square-root adjoint
call nicas%apply_sqrt_ad(nam,geom,bpar,fld,cv)

! Apply square-root
call nicas%apply_sqrt(nam,geom,bpar,cv,fld)

if (lhook) call dr_hook('nicas_apply_from_sqrt',1,zhook_handle)

end subroutine nicas_apply_from_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_apply_sqrt
!> Purpose: apply localization square-root
!----------------------------------------------------------------------
subroutine nicas_apply_sqrt(nicas,nam,geom,bpar,cv,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                 !< NICAS data
type(nam_type),target,intent(in) :: nam                               !< Namelist
type(geom_type),target,intent(in) :: geom                             !< Geometry
type(bpar_type),target,intent(in) :: bpar                             !< Block parameters
type(cv_type),intent(in) :: cv                                        !< Control variable
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ib,its,iv,jv,i
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),a(:),u(:)

if (lhook) call dr_hook('nicas_apply_sqrt',0,zhook_handle)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Apply common localization
   call nicas%blk(bpar%nb+1)%apply_sqrt(geom,cv%blk(bpar%nb+1)%alpha,fld_3d)

   ! Apply common ensemble coefficient square-root
   if (lcoef_ens) fld_3d = fld_3d*sqrt(nicas%blk(bpar%nb+1)%coef_ens)

   ! Build final vector
   do its=1,nam%nts
      do iv=1,nam%nv
         fld(:,:,iv,its) = fld_3d
      end do
   end do
case ('specific_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific localization (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt(geom,cv%blk(ib)%alpha,fld_4d(:,:,iv))

         ! Apply specific ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('specific_multivariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific localization (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt(geom,cv%blk(bpar%nb+1)%alpha,fld_4d(:,:,iv))

         ! Apply specific ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(ib)%coef_ens)
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))
   allocate(a((nam%nv*(nam%nv+1))/2))
   allocate(u((nam%nv*(nam%nv+1))/2))

   ! Copy weights
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         if (isnotmsr(nicas%blk(ib)%wgt)) then
            wgt(iv,jv) = nicas%blk(ib)%wgt
         else
            wgt(iv,jv) = 0.0
         end if
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Cholesky decomposition
   i = 0
   do iv=1,nam%nv
      do jv=1,iv
         i = i+1
         a(i) = wgt(iv,jv)
      end do
   end do
   call cholesky(nam%nv,(nam%nv*(nam%nv+1))/2,a,u)
   i = 0
   wgt = 0.0
   do iv=1,nam%nv
      do jv=1,iv
         i = i+1
         wgt(iv,jv) = u(i)
      end do
   end do

   do ib=1,bpar%nb
      if (bpar%cv_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific localization (same for all timeslots)
         call nicas%blk(bpar%nb+1)%apply_sqrt(geom,cv%blk(ib)%alpha,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(bpar%nb+1)%coef_ens)
      end if
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=1,iv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d_tmp
   end do
end select

if (lhook) call dr_hook('nicas_apply_sqrt',1,zhook_handle)

end subroutine nicas_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_apply_sqrt_ad
!> Purpose: apply localization square-root, adjoint
!----------------------------------------------------------------------
subroutine nicas_apply_sqrt_ad(nicas,nam,geom,bpar,fld,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                !< NICAS data
type(nam_type),target,intent(in) :: nam                              !< Namelist
type(geom_type),target,intent(in) :: geom                            !< Geometry
type(bpar_type),target,intent(in) :: bpar                            !< Block parameters
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field
type(cv_type),intent(out) :: cv                                      !< Control variable

! Local variable
integer :: ib,its,iv,jv,i
real(kind_real) :: zhook_handle
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),a(:),u(:)
type(cv_type) :: cv_tmp

if (lhook) call dr_hook('nicas_apply_sqrt_ad',0,zhook_handle)

! Allocation
call nicas%alloc_cv(bpar,cv)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Sum product over variables and timeslots
   fld_3d = 0.0
   do its=1,nam%nts
      do iv=1,nam%nv
         fld_3d = fld_3d+fld(:,:,iv,its)
      end do
   end do

   ! Apply common ensemble coefficient square-root
   if (lcoef_ens) fld_3d = fld_3d*sqrt(nicas%blk(bpar%nb+1)%coef_ens)

   ! Apply common localization
   call nicas%blk(bpar%nb+1)%apply_sqrt_ad(geom,fld_3d,cv%blk(bpar%nb+1)%alpha)
case ('specific_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt_ad(geom,fld_4d(:,:,iv),cv%blk(ib)%alpha)
      end if
   end do
case ('specific_multivariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   call nicas%alloc_cv(bpar,cv_tmp)

   ! Initialization
   cv%blk(bpar%nb+1)%alpha = 0.0

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d(:,:,iv) = fld_4d(:,:,iv)*sqrt(nicas%blk(ib)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt_ad(geom,fld_4d(:,:,iv),cv_tmp%blk(bpar%nb+1)%alpha)

         ! Sum control variable
         cv%blk(bpar%nb+1)%alpha = cv%blk(bpar%nb+1)%alpha+cv_tmp%blk(bpar%nb+1)%alpha
      end if
   end do
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))
   allocate(a((nam%nv*(nam%nv+1))/2))
   allocate(u((nam%nv*(nam%nv+1))/2))

   ! Copy weights
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         if (isnotmsr(nicas%blk(ib)%wgt)) then
            wgt(iv,jv) = nicas%blk(ib)%wgt
         else
            wgt(iv,jv) = 0.0
         end if
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Cholesky decomposition
   i = 0
   do iv=1,nam%nv
      do jv=1,iv
         i = i+1
         a(i) = wgt(iv,jv)
      end do
   end do
   call cholesky(nam%nv,(nam%nv*(nam%nv+1))/2,a,u)
   i = 0
   wgt = 0.0
   do jv=1,nam%nv
      do iv=1,jv
         i = i+1
         wgt(iv,jv) = u(i)
      end do
   end do

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=iv,nam%nv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   do ib=1,bpar%nb
      if (bpar%cv_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         if (lcoef_ens) fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)*sqrt(nicas%blk(bpar%nb+1)%coef_ens)

         ! Apply specific localization (same for all timeslots)
         call nicas%blk(bpar%nb+1)%apply_sqrt_ad(geom,fld_4d_tmp(:,:,iv),cv%blk(ib)%alpha)
      end if
   end do
end select

if (lhook) call dr_hook('nicas_apply_sqrt_ad',1,zhook_handle)

end subroutine nicas_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: nicas_randomize
!> Purpose: randomize localization from square-root
!----------------------------------------------------------------------
subroutine nicas_randomize(nicas,nam,geom,bpar,ne,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                    !< NICAS data
type(nam_type),target,intent(in) :: nam                                  !< Namelist
type(geom_type),target,intent(in) :: geom                                !< Geometry
type(bpar_type),target,intent(in) :: bpar                                !< Blocal parameters
integer,intent(in) :: ne                                                 !< Number of members
real(kind_real),intent(out) :: ens(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne) !< Ensemble

! Local variable
integer :: ie
real(kind_real) :: mean(geom%nc0a,geom%nl0,nam%nv,nam%nts),std(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: zhook_handle
type(cv_type) :: cv_ens(ne)

if (lhook) call dr_hook('nicas_randomize',0,zhook_handle)

do ie=1,ne
   ! Generate random control vector
   call nicas%random_cv(bpar,cv_ens(ie))

   ! Apply square-root
   call nicas%apply_sqrt(nam,geom,bpar,cv_ens(ie),ens(:,:,:,:,ie))
end do

! Normalize ensemble
mean = sum(ens,dim=5)/float(ne)
do ie=1,ne
   ens(:,:,:,:,ie) = ens(:,:,:,:,ie)-mean
end do
std = sqrt(sum(ens**2,dim=5)/float(ne-1))
do ie=1,ne
   ens(:,:,:,:,ie) = ens(:,:,:,:,ie)/std
end do

if (lhook) call dr_hook('nicas_randomize',1,zhook_handle)

end subroutine nicas_randomize

!----------------------------------------------------------------------
! Subroutine: nicas_apply_bens
!> Purpose: apply localized ensemble covariance
!----------------------------------------------------------------------
subroutine nicas_apply_bens(nicas,nam,geom,bpar,ens1,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                             !< NICAS data
type(nam_type),target,intent(in) :: nam                                           !< Namelist
type(geom_type),target,intent(in) :: geom                                         !< Geometry
type(bpar_type),target,intent(in) :: bpar                                         !< Blocal parameters
real(kind_real),intent(in) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts)           !< Field

! Local variable
integer :: ie
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv,nam%nts),fld_tmp(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: mean(geom%nc0a,geom%nl0,nam%nv,nam%nts),pert(geom%nc0a,geom%nl0,nam%nv,nam%nts)

! Compute mean
mean = sum(ens1,dim=5)/float(nam%ens1_ne)

! Copy field
fld_copy = fld

! Adjoint advection
if (nam%advmode==-1) call nicas%blk(bpar%nb+1)%apply_adv_ad(nam,geom,fld_copy)

! Apply localized ensemble covariance formula
fld = 0.0
do ie=1,nam%ens1_ne
   ! Compute perturbation
   pert = (ens1(:,:,:,:,ie)-mean)/sqrt(float(nam%ens1_ne-1))

   ! Inverse advection
   if (nam%advmode==-1) call nicas%blk(bpar%nb+1)%apply_adv_inv(nam,geom,pert)

   ! Schur product
   fld_tmp = pert*fld_copy

   ! Apply localization
   if (nam%lsqrt) then
      call nicas%apply_from_sqrt(nam,geom,bpar,fld_tmp)
   else
      call nicas%apply(nam,geom,bpar,fld_tmp)
   end if

   ! Schur product
   fld = fld+fld_tmp*pert
end do

! Advection
if (nam%advmode==-1) call nicas%blk(bpar%nb+1)%apply_adv(nam,geom,fld)

end subroutine nicas_apply_bens

!----------------------------------------------------------------------
! Subroutine: nicas_apply_bens_noloc
!> Purpose: apply ensemble covariance, without localization
!----------------------------------------------------------------------
subroutine nicas_apply_bens_noloc(nicas,nam,geom,ens1,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                             !< NICAS data
type(nam_type),target,intent(in) :: nam                                           !< Namelist
type(geom_type),target,intent(in) :: geom                                         !< Geometry
real(kind_real),intent(in) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts)           !< Field

! Local variable
integer :: ie
real(kind_real) :: alpha
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: mean(geom%nc0a,geom%nl0,nam%nv,nam%nts),pert(geom%nc0a,geom%nl0,nam%nv,nam%nts)

! Compute mean
mean = sum(ens1,dim=5)/float(nam%ens1_ne)

! Initialization
fld_copy = fld

! Apply localized ensemble covariance formula
fld = 0.0
do ie=1,nam%ens1_ne
   ! Compute perturbation
   pert = (ens1(:,:,:,:,ie)-mean)/sqrt(float(nam%ens1_ne-1))

   ! Dot product
   call mpl%dot_prod(pert,fld_copy,alpha)

   ! Schur product
   fld = fld+alpha*pert
end do

end subroutine nicas_apply_bens_noloc

!----------------------------------------------------------------------
! Subroutine: nicas_test_adjoint
!> Purpose: test localization adjoint
!----------------------------------------------------------------------
subroutine nicas_test_adjoint(nicas,nam,geom,bpar,ens1)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                                      !< NICAS data
type(nam_type),intent(in) :: nam                                                           !< Namelist
type(geom_type),intent(in) :: geom                                                         !< Geometry
type(bpar_type),intent(in) :: bpar                                                         !< Block parameters
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real),allocatable :: fld1_loc(:,:,:,:),fld1_adv(:,:,:,:),fld1_bens(:,:,:,:),fld1_save(:,:,:,:)
real(kind_real),allocatable :: fld2_loc(:,:,:,:),fld2_adv(:,:,:,:),fld2_bens(:,:,:,:),fld2_save(:,:,:,:)

! Allocation
allocate(fld1_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld2_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld1_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld2_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (present(ens1)) then
   allocate(fld1_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   allocate(fld2_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
end if

! Generate random field
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld2_save)

! Adjoint test
fld1_loc = fld1_save
fld2_loc = fld2_save
if (nam%lsqrt) then
   call nicas%apply_from_sqrt(nam,geom,bpar,fld1_loc)
   call nicas%apply_from_sqrt(nam,geom,bpar,fld2_loc)
else
   call nicas%apply(nam,geom,bpar,fld1_loc)
   call nicas%apply(nam,geom,bpar,fld2_loc)
end if
if (abs(nam%advmode)==1) then
   fld1_adv = fld1_save
   fld2_adv = fld2_save
   call nicas%blk(bpar%nb+1)%apply_adv(nam,geom,fld1_adv)
   call nicas%blk(bpar%nb+1)%apply_adv_ad(nam,geom,fld2_adv)
end if
if (present(ens1)) then
   fld1_bens = fld1_save
   fld2_bens = fld2_save
   call nicas%apply_bens(nam,geom,bpar,ens1,fld1_bens)
   call nicas%apply_bens(nam,geom,bpar,ens1,fld2_bens)
end if

! Print result
call mpl%dot_prod(fld1_loc,fld2_save,sum1)
call mpl%dot_prod(fld2_loc,fld1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Localization adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
if (abs(nam%advmode)==1) then
   call mpl%dot_prod(fld1_adv,fld2_save,sum1)
   call mpl%dot_prod(fld2_adv,fld1_save,sum2)
   write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Advection adjoint test:    ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
end if
if (present(ens1)) then
   call mpl%dot_prod(fld1_bens,fld2_save,sum1)
   call mpl%dot_prod(fld2_bens,fld1_save,sum2)
   write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Ensemble B adjoint test:   ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
end if

end subroutine nicas_test_adjoint

!----------------------------------------------------------------------
! Subroutine: nicas_test_sqrt
!> Purpose: test full/square-root equivalence
!----------------------------------------------------------------------
subroutine nicas_test_sqrt(nicas,nam,geom,bpar,cmat,ens1)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                                      !< NICAS data
type(nam_type),intent(inout),target :: nam                                                 !< Namelist
type(geom_type),intent(in),target :: geom                                                  !< Geometry
type(bpar_type),intent(in) :: bpar                                                         !< Block parameters
type(cmat_type),intent(in) :: cmat                                                         !< C matrix data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ib,ic0,ic0a,iv
real(kind_real),allocatable :: fld_loc(:,:,:,:),fld_loc_sqrt(:,:,:,:)
real(kind_real),allocatable :: fld_bens(:,:,:,:),fld_bens_sqrt(:,:,:,:)
character(len=1024) :: varname(nam%nv)
type(nicas_type) :: nicas_other

! Allocation
allocate(fld_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld_loc_sqrt(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (present(ens1)) then
   allocate(fld_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   allocate(fld_bens_sqrt(geom%nc0a,geom%nl0,nam%nv,nam%nts))
end if

! Generate random field
call rng%rand_real(-1.0_kind_real,1.0_kind_real,fld_loc)
fld_loc_sqrt = fld_loc
if (present(ens1)) then
   call rng%rand_real(-1.0_kind_real,1.0_kind_real,fld_bens)
   fld_bens_sqrt = fld_bens
end if

! Apply localization, initial version
if (nam%lsqrt) then
   call nicas%apply_from_sqrt(nam,geom,bpar,fld_loc_sqrt)
else
   call nicas%apply(nam,geom,bpar,fld_loc)
end if
if (present(ens1)) then
   if (nam%lsqrt) then
      call nicas%apply_bens(nam,geom,bpar,ens1,fld_bens_sqrt)
   else
      call nicas%apply_bens(nam,geom,bpar,ens1,fld_bens)
   end if
end if

! Switch lsqrt
nam%lsqrt = .not.nam%lsqrt

! Allocation
call nicas_other%alloc(nam,bpar,'nicas_other')

! Prepare nicas, other version
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      ! Compute NICAS parameters
      call nicas_other%blk(ib)%compute_parameters(nam,geom,cmat%blk(ib))
   end if

   if (bpar%B_block(ib)) then
      ! Copy weights
      nicas_other%blk(ib)%wgt = cmat%blk(ib)%wgt
      if (bpar%nicas_block(ib)) then
         allocate(nicas_other%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
         do ic0=1,geom%nc0
            if (geom%c0_to_proc(ic0)==mpl%myproc) then
               ic0a = geom%c0_to_c0a(ic0)
               nicas_other%blk(ib)%coef_ens(ic0a,:) = cmat%blk(ib)%coef_ens(ic0,:)
            end if
         end do
      end if
   end if
end do

! Apply localization, other version
if (nam%lsqrt) then
   call nicas_other%apply_from_sqrt(nam,geom,bpar,fld_loc_sqrt)
else
   call nicas_other%apply(nam,geom,bpar,fld_loc)
end if
if (present(ens1)) then
   if (nam%lsqrt) then
      call nicas_other%apply_bens(nam,geom,bpar,ens1,fld_bens_sqrt)
   else
      call nicas_other%apply_bens(nam,geom,bpar,ens1,fld_bens)
   end if
end if

! Compute dirac
do iv=1,nam%nv
   varname(iv) = nam%varname(iv)
   nam%varname(iv) = trim(varname(iv))//'_sqrt'
end do
if (nam%check_dirac) call nicas_other%test_dirac(nam,geom,bpar,ens1)
do iv=1,nam%nv
   nam%varname(iv) = varname(iv)
end do

! Reset lsqrt value
nam%lsqrt = .not.nam%lsqrt

! Print difference
write(mpl%unit,'(a7,a,f6.1,a)') '','Localization full / square-root error : ', &
 & sqrt(sum((fld_loc_sqrt-fld_loc)**2)/sum(fld_loc**2))*100.0,'%'
if (present(ens1)) write(mpl%unit,'(a7,a,f6.1,a)') '','Ensemble B full / square-root error:  ', &
 & sqrt(sum((fld_bens_sqrt-fld_bens)**2)/sum(fld_bens**2))*100.0,'%'

end subroutine nicas_test_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_test_dirac
!> Purpose: apply localization to diracs
!----------------------------------------------------------------------
subroutine nicas_test_dirac(nicas,nam,geom,bpar,ens1)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                                      !< NICAS data
type(nam_type),intent(in) :: nam                                                           !< Namelist
type(geom_type),intent(in) :: geom                                                         !< Geometry
type(bpar_type),intent(in) :: bpar                                                         !< Block parameters
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: iprocdir(nam%ndir),ic0adir(nam%ndir),il0dir(nam%ndir),idir,iv,its
real(kind_real),allocatable :: fld(:,:,:,:),fld_loc(:,:,:,:),fld_bens(:,:,:,:)
character(len=2) :: itschar
character(len=1024) :: filename

! Allocation
allocate(fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))

! Find gridpoint and level indices
call define_dirac(nam,geom,iprocdir,ic0adir,il0dir)

! Generate dirac field
fld = 0.0
do idir=1,nam%ndir
   if (iprocdir(idir)==mpl%myproc) fld(ic0adir(idir),il0dir(idir),nam%ivdir(idir),nam%itsdir(idir)) = 1.0
end do

! Allocation
allocate(fld_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (present(ens1)) allocate(fld_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))

! Apply localization to dirac
fld_loc = fld
if (nam%lsqrt) then
   call nicas%apply_from_sqrt(nam,geom,bpar,fld_loc)
else
   call nicas%apply(nam,geom,bpar,fld_loc)
end if

if (present(ens1)) then
   ! Apply localized ensemble covariance
   fld_bens = fld
   call nicas%apply_bens(nam,geom,bpar,ens1,fld_bens)
end if

! Write field
filename = trim(nam%prefix)//'_dirac.nc'
do its=1,nam%nts
   write(itschar,'(i2.2)') its
   do iv=1,nam%nv
      call geom%fld_write(nam,trim(nam%prefix)//'_dirac.nc',trim(nam%varname(iv))//'_'//itschar,fld_loc(:,:,iv,its))
      if (present(ens1)) call geom%fld_write(nam,trim(nam%prefix)//'_dirac.nc',trim(nam%varname(iv))//'_'//itschar//'_Bens', &
       & fld_bens(:,:,iv,its))
   end do
end do

! Write gridded field
filename = trim(nam%prefix)//'_dirac_gridded.nc'
do its=1,nam%nts
   write(itschar,'(i2.2)') its
   do iv=1,nam%nv
      call model_write(nam,geom,trim(nam%prefix)//'_dirac_gridded.nc',trim(nam%varname(iv))//'_'//itschar,fld_loc(:,:,iv,its))
      if (present(ens1)) call model_write(nam,geom,trim(nam%prefix)//'_dirac_gridded.nc', &
       & trim(nam%varname(iv))//'_'//itschar//'_Bens',fld_bens(:,:,iv,its))
   end do
end do

end subroutine nicas_test_dirac

!----------------------------------------------------------------------
! Subroutine: nicas_test_randomization
!> Purpose: test NICAS randomization method with respect to theoretical error statistics
!----------------------------------------------------------------------
subroutine nicas_test_randomization(nicas,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(nam_type),intent(inout) :: nam   !< Namelist variables
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters

! Local variables
integer :: ifac,itest,nefac(nfac),ens1_ne,iv,its
real(kind_real) :: fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest),fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest)
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts),mse(ntest,nfac),mse_th(ntest,nfac)
real(kind_real),allocatable :: ens1(:,:,:,:,:)
character(len=2) :: itschar
character(len=4) :: nechar,itestchar
character(len=1024) :: filename

! Define test vectors
write(mpl%unit,'(a4,a)') '','Define test vectors'
call define_test_vectors(nam,geom,ntest,fld_save)

! Apply localization to test vectors
write(mpl%unit,'(a4,a)') '','Apply localization to test vectors'
fld_ref = fld_save
do itest=1,ntest
   call nicas%apply_from_sqrt(nam,geom,bpar,fld_ref(:,:,:,:,itest))
end do

! Write first 10 test vectors
do itest=1,min(ntest,10)
   ! Write field
   write(itestchar,'(i4.4)') itest
   filename = trim(nam%prefix)//'_randomize_'//itestchar//'_gridded.nc'
   do its=1,nam%nts
      write(itschar,'(i2.2)') its
      do iv=1,nam%nv
         call model_write(nam,geom,filename,trim(nam%varname(iv))//'_ref_'//itschar,fld_ref(:,:,iv,its,itest))
      end do
   end do
end do

! Save namelist variables
ens1_ne = nam%ens1_ne

write(mpl%unit,'(a4,a)') '','Test randomization for various ensemble sizes:'
do ifac=1,nfac
   ! Ensemble size
   nefac(ifac) = max(int(5.0*float(ifac)/float(nfac)*float(ne_rand)),3)
   nam%ens1_ne = nefac(ifac)
   write(nechar,'(i4.4)') nefac(ifac)

   ! Allocation
   allocate(ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nefac(ifac)))

   ! Randomize ensemble
   call nicas%randomize(nam,geom,bpar,nefac(ifac),ens1)

   do itest=1,ntest
      ! Test localization
      fld = fld_save(:,:,:,:,itest)
      call nicas%apply_bens_noloc(nam,geom,ens1,fld)

      ! RMSE
      mse(itest,ifac) = sum((fld-fld_ref(:,:,:,:,itest))**2)
      mse_th(itest,ifac) = 1.0/float(nam%ens1_ne-1)*sum(1+fld_ref(:,:,:,:,itest)**2)

      ! Write first 10 test vectors
      if (itest<=min(ntest,10)) then
         ! Write field
         write(itestchar,'(i4.4)') itest
         filename = trim(nam%prefix)//'_randomize_'//itestchar//'_gridded.nc'
         do its=1,nam%nts
            write(itschar,'(i2.2)') its
            do iv=1,nam%nv
               call model_write(nam,geom,filename,trim(nam%varname(iv))//'_rand_'//nechar//'_'//itschar,fld(:,:,iv,its))
            end do
         end do
      end if
   end do

   ! Print scores
   write(mpl%unit,'(a7,a,i4,a,e15.8,a,e15.8)') '','Ensemble size ',nefac(ifac),', MSE (exp. / th.): ', &
 & sum(mse(:,ifac))/float(ntest),' / ',sum(mse_th(:,ifac))/float(ntest)

   ! Release memory
   deallocate(ens1)
end do

! Reset namelist variables
nam%ens1_ne = ens1_ne

end subroutine nicas_test_randomization

!----------------------------------------------------------------------
! Subroutine: nicas_test_consistency
!> Purpose: test HDIAG_NICAS consistency with a randomization method
!----------------------------------------------------------------------
subroutine nicas_test_consistency(nicas,nam,geom,bpar,cmat)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(nam_type),intent(inout) :: nam   !< Namelist variables
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(cmat_type),intent(in) :: cmat    !< C matrix data

! Local variables
integer :: ens1_ne,ens1_ne_offset,ens1_nsub,ib,il0
real(kind_real) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne_rand)
character(len=1024) :: prefix,method
type(cmat_type) :: cmat_test

! Randomize ensemble
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Randomize ensemble'
call nicas%randomize(nam,geom,bpar,ne_rand,ens1)

! Copy sampling
call system('cp -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc ' &
 & //trim(nam%datadir)//'/'//trim(nam%prefix)//'_consistency-test_sampling.nc')

! Save namelist variables
prefix = nam%prefix
method = nam%method
ens1_ne = nam%ens1_ne
ens1_ne_offset = nam%ens1_ne_offset
ens1_nsub = nam%ens1_nsub

! Set namelist variables
nam%prefix = trim(nam%prefix)//'_consistency-test'
nam%method = 'cor'
nam%ens1_ne = ne_rand
nam%ens1_ne_offset = 0
nam%ens1_nsub = 1

! Call hdiag driver
call run_hdiag(nam,geom,bpar,cmat_test,ens1)

! Print scores
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- hdiag_nicas consistency results'
do ib=1,bpar%nb+1
   if (bpar%nicas_block(ib)) then
      write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      do il0=1,geom%nl0
         write(mpl%unit,'(a10,a7,i3,a4,a25,f6.1,a)') '','Level: ',nam%levs(il0),' ~> ','horizontal length-scale: ', &
       & sum(cmat_test%blk(ib)%rh0(:,il0)-cmat%blk(ib)%rh0(:,il0))/float(geom%nc0)*reqkm,' km'
         if (any(abs(cmat%blk(ib)%rv0(:,il0))>0.0)) then
            write(mpl%unit,'(a49,f6.1,a)') 'vertical length-scale: ', &
          & sum(cmat_test%blk(ib)%rv0(:,il0)-cmat%blk(ib)%rv0(:,il0))/float(geom%nc0),' '//trim(vunitchar)
         end if
      end do
   end if
end do

! Reset namelist variables
nam%prefix = prefix
nam%method = method
nam%ens1_ne = ens1_ne
nam%ens1_ne_offset = ens1_ne_offset
nam%ens1_nsub = ens1_nsub

end subroutine nicas_test_consistency

!----------------------------------------------------------------------
! Subroutine: nicas_test_optimality
!> Purpose: test HDIAG localization optimality with a randomization method
!----------------------------------------------------------------------
subroutine nicas_test_optimality(nicas,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(nam_type),intent(inout) :: nam   !< Namelist variables
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters

! Local variables
integer :: ib,ic0,ic0a,ifac,itest
real(kind_real) :: fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest),fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest)
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts),fac(nfac),mse(ntest,nfac)
real(kind_real) :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne)
character(len=1024) :: prefix,method
type(cmat_type) :: cmat_save,cmat_test
type(nicas_type) :: nicas_test

! Define test vectors
write(mpl%unit,'(a4,a)') '','Define test vectors'
call define_test_vectors(nam,geom,ntest,fld_save)

! Apply localization to test vectors
write(mpl%unit,'(a4,a)') '','Apply localization to test vectors'
fld_ref = fld_save
do itest=1,ntest
   call nicas%apply_from_sqrt(nam,geom,bpar,fld_ref(:,:,:,:,itest))
end do

! Randomize ensemble
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Randomize ensemble'
call nicas%randomize(nam,geom,bpar,nam%ens1_ne,ens1)

! Copy sampling
call system('cp -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc ' &
 & //trim(nam%datadir)//'/'//trim(nam%prefix)//'_optimality-test_sampling.nc')

! Save namelist variables
prefix = nam%prefix
method = nam%method

! Set namelist variables
nam%prefix = trim(nam%prefix)//'_optimality-test'
nam%method = 'loc'

! Allocation
call nicas_test%alloc(nam,bpar,'nicas_test')
call cmat_save%alloc(nam,geom,bpar,'cmat_save')

! Call hdiag driver
call run_hdiag(nam,geom,bpar,cmat_save,ens1)

! Copy cmat
cmat_test = cmat_save%copy(nam,geom,bpar)

do ifac=1,nfac
   ! Multiplication factor
   fac(ifac) = 2.0*float(ifac)/float(nfac)

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,f4.2,a)') '--- Apply a multiplicative factor ',fac(ifac),' to length-scales'

   do ib=1,bpar%nb+1
      if (bpar%nicas_block(ib)) then
         ! Length-scales multiplication
         cmat_test%blk(ib)%rh0 = fac(ifac)*cmat_save%blk(ib)%rh0
         cmat_test%blk(ib)%rv0 = fac(ifac)*cmat_save%blk(ib)%rv0
         if (trim(nam%strategy)=='specific_multivariate') then
            cmat_test%blk(ib)%rh0s = fac(ifac)*cmat_save%blk(ib)%rh0s
            cmat_test%blk(ib)%rv0s = fac(ifac)*cmat_save%blk(ib)%rv0s
         end if

         ! Compute NICAS parameters
         call nicas_test%blk(ib)%compute_parameters(nam,geom,cmat_test%blk(ib))
      end if

      if (bpar%B_block(ib)) then
         ! Copy weights
         nicas_test%blk(ib)%wgt = cmat_test%blk(ib)%wgt
         if (bpar%nicas_block(ib)) then
            allocate(nicas_test%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
            do ic0=1,geom%nc0
               if (geom%c0_to_proc(ic0)==mpl%myproc) then
                  ic0a = geom%c0_to_c0a(ic0)
                  nicas_test%blk(ib)%coef_ens(ic0a,:) = cmat_test%blk(ib)%coef_ens(ic0,:)
               end if
            end do
         end if
      end if
   end do

   do itest=1,ntest
      ! Test localization
      fld = fld_save(:,:,:,:,itest)
      call nicas_test%apply_bens(nam,geom,bpar,ens1,fld)

      ! RMSE
      mse(itest,ifac) = sum((fld-fld_ref(:,:,:,:,itest))**2)
   end do

   ! Print scores
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,f4.2,a,e15.8)') '--- Optimality results for a factor ',fac(ifac),', MSE: ',sum(mse(:,ifac))/float(ntest)
end do

! Print scores summary
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Optimality results summary'
do ifac=1,nfac
   write(mpl%unit,'(a7,a,f4.2,a,e15.8)') '','Factor ',fac(ifac),', MSE: ',sum(mse(:,ifac))/float(ntest)
end do

! Reset namelist variables
nam%prefix = prefix
nam%method = method

end subroutine nicas_test_optimality

end module type_nicas
