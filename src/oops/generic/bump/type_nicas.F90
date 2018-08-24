!----------------------------------------------------------------------
! Module: type_nicas
!> Purpose: NICAS data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_nicas

use netcdf
use tools_const, only: rad2deg,reqkm,pi,msvali,msvalr
use tools_func, only: cholesky
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr
use tools_nc, only: ncfloat
use tools_test, only: define_dirac,define_test_vectors
use type_bpar, only: bpar_type
use type_cmat, only: cmat_type
use type_com, only: com_type
use type_cv, only: cv_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_linop, only: linop_type
use type_nicas_blk, only: nicas_blk_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

integer,parameter :: ne_rand = 150          !< Ensemble size for randomization
integer,parameter :: nfac = 10              !< Number of length-scale factors
integer,parameter :: ntest = 100            !< Number of tests
logical,parameter :: pos_def_test = .false. !< Positive-definiteness test

! NICAS derived type
type nicas_type
   character(len=1024) :: prefix              !< Prefix
   type(nicas_blk_type),allocatable :: blk(:) !< NICAS data blocks
   logical :: allocated                       !< Allocation flag
contains
   procedure :: alloc => nicas_alloc
   procedure :: dealloc => nicas_dealloc
   procedure :: read => nicas_read
   procedure :: write => nicas_write
   procedure :: write_mpi_summary => nicas_write_mpi_summary
   procedure :: run_nicas => nicas_run_nicas
   procedure :: run_nicas_tests => nicas_run_nicas_tests
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

private
public :: nicas_type

contains

!----------------------------------------------------------------------
! Subroutine: nicas_alloc
!> Purpose: NICAS data allocation
!----------------------------------------------------------------------
subroutine nicas_alloc(nicas,mpl,nam,bpar,prefix)

! Passed variables
class(nicas_type),intent(inout) :: nicas !< NICAS data
type(mpl_type),intent(in) :: mpl         !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
type(bpar_type),intent(in) :: bpar       !< Block parameters
character(len=*),intent(in) :: prefix    !< Prefix

! Local variable
integer :: ib

! Copy prefix
nicas%prefix = prefix

! Allocation
allocate(nicas%blk(bpar%nbe))

! Set name
do ib=1,bpar%nbe
   nicas%blk(ib)%ib = ib
   if (nam%lsqrt) then
      write(nicas%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(prefix)//'-',nam%mpicom,'-sqrt_',mpl%nproc,'-',mpl%myproc, &
    & '_',trim(bpar%blockname(ib))
   else
      write(nicas%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(prefix)//'-',nam%mpicom,'_',mpl%nproc,'-',mpl%myproc, &
    & '_',trim(bpar%blockname(ib))
   end if
end do

! Update allocation flag
nicas%allocated = .true.

end subroutine nicas_alloc

!----------------------------------------------------------------------
! Subroutine: nicas_dealloc
!> Purpose: NICAS data deallocation
!----------------------------------------------------------------------
subroutine nicas_dealloc(nicas,nam,geom,bpar)

! Passed variables
class(nicas_type),intent(inout) :: nicas !< NICAS data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
type(bpar_type),intent(in) :: bpar       !< Block parameters

! Local variable
integer :: ib

! Release memory
if (allocated(nicas%blk)) then
   do ib=1,bpar%nbe
      call nicas%blk(ib)%dealloc(nam,geom)
   end do
   deallocate(nicas%blk)
end if

! Update allocation flag
nicas%allocated = .false.

end subroutine nicas_dealloc

!----------------------------------------------------------------------
! Subroutine: nicas_read
!> Purpose: read NICAS data
!----------------------------------------------------------------------
subroutine nicas_read(nicas,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas !< NICAS data
type(mpl_type),intent(in) :: mpl         !< MPI data
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

! Allocation
call nicas%alloc(mpl,nam,bpar,'nicas')

do ib=1,bpar%nbe
   if (bpar%B_block(ib)) then
      ! Open file and get dimensions
      filename = trim(nam%prefix)//'_'//trim(nicas%blk(ib)%name)//'.nc'
      call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))

      ! Get dimensions
      call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=nl0_test))
      if (nl0_test/=geom%nl0) call mpl%abort('wrong dimension when reading nicas')
      if (bpar%nicas_block(ib)) then
         call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc0a',nc0a_id))
         call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc0a_id,len=nc0a_test))
         if (nc0a_test/=geom%nc0a) call mpl%abort('wrong dimension when reading nicas')
         info = nf90_inq_dimid(ncid,'nc1b',nc1b_id)
         if (info==nf90_noerr) then
            call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc1b_id,len=nicas%blk(ib)%nc1b))
         else
            nicas%blk(ib)%nc1b = 0
         end if
         call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
         call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=nicas%blk(ib)%nl1))
         info = nf90_inq_dimid(ncid,'nsa',nsa_id)
         if (info==nf90_noerr) then
            call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nsa_id,len=nicas%blk(ib)%nsa))
         else
            nicas%blk(ib)%nsa = 0
         end if
         info = nf90_inq_dimid(ncid,'nsb',nsb_id)
         if (info==nf90_noerr) then
            call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nsb_id,len=nicas%blk(ib)%nsb))
         else
            nicas%blk(ib)%nsb = 0
         end if
         call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'nsc',nicas%blk(ib)%nsc))
      end if
      if ((ib==bpar%nbe).and.(abs(nam%advmode)==1)) then
         call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc0d',nc0d_id))
         call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc0d_id,len=nicas%blk(ib)%nc0d))
         call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc0dinv',nc0dinv_id))
         call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc0dinv_id,len=nicas%blk(ib)%nc0dinv))
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
      if ((ib==bpar%nbe).and.(abs(nam%advmode)==1)) then
         allocate(nicas%blk(ib)%d(geom%nl0,2:nam%nts))
         allocate(nicas%blk(ib)%dinv(geom%nl0,2:nam%nts))
      end if

      ! Get variable id
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_c1b',sb_to_c1b_id))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_l1',sb_to_l1_id))
         if (nicas%blk(ib)%nsa>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,'sa_to_sc',sa_to_sc_id))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_sc',sb_to_sc_id))
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))
      end if

      ! Read data
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_c1b_id,nicas%blk(ib)%sb_to_c1b))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_l1_id,nicas%blk(ib)%sb_to_l1))
         if (nicas%blk(ib)%nsa>0) call mpl%ncerr(subr,nf90_get_var(ncid,sa_to_sc_id,nicas%blk(ib)%sa_to_sc))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_sc_id,nicas%blk(ib)%sb_to_sc))
         call mpl%ncerr(subr,nf90_get_var(ncid,norm_id,nicas%blk(ib)%norm))
         call mpl%ncerr(subr,nf90_get_var(ncid,coef_ens_id,nicas%blk(ib)%coef_ens))
         call nicas%blk(ib)%com_AB%read(mpl,ncid,'com_AB')
         call nicas%blk(ib)%com_AC%read(mpl,ncid,'com_AC')
         nicas%blk(ib)%c%prefix = 'c'
         call nicas%blk(ib)%c%read(mpl,ncid)
         do il0i=1,geom%nl0i
            write(nicas%blk(ib)%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
            call nicas%blk(ib)%h(il0i)%read(mpl,ncid)
         end do
         nicas%blk(ib)%v%prefix = 'v'
         call nicas%blk(ib)%v%read(mpl,ncid)
         do il1=1,nicas%blk(ib)%nl1
            write(nicas%blk(ib)%s(il1)%prefix,'(a,i3.3)') 's_',il1
            call nicas%blk(ib)%s(il1)%read(mpl,ncid)
         end do
      end if
      if ((ib==bpar%nbe).and.(abs(nam%advmode)==1)) then
         call nicas%blk(ib)%com_AD%read(mpl,ncid,'com_AD')
         call nicas%blk(ib)%com_ADinv%read(mpl,ncid,'com_ADinv')
         do its=2,nam%nts
            do il0=1,geom%nl0
               write(nicas%blk(ib)%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
               call nicas%blk(ib)%d(il0,its)%read(mpl,ncid)
               write(nicas%blk(ib)%dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'dinv_',il0,'_',its
               call nicas%blk(ib)%dinv(il0,its)%read(mpl,ncid)
            end do
         end do
      end if

      ! Read main weight
      call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',nicas%blk(ib)%wgt))

      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   end if
end do

end subroutine nicas_read

!----------------------------------------------------------------------
! Subroutine: nicas_write
!> Purpose: write NICAS data
!----------------------------------------------------------------------
subroutine nicas_write(nicas,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(mpl_type),intent(in) :: mpl      !< MPI data
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

do ib=1,bpar%nbe
   if (bpar%B_block(ib)) then
      ! Create file
      filename = trim(nam%prefix)//'_'//trim(nicas%blk(ib)%name)//'.nc'
      call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

      ! Write namelist parameters
      call nam%ncwrite(mpl,ncid)

      ! Define dimensions
      call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
      if (bpar%nicas_block(ib)) then
         call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0a',geom%nc0a,nc0a_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',geom%nl0i))
         if (nicas%blk(ib)%nc1b>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1b',nicas%blk(ib)%nc1b,nc1b_id))
         call mpl%ncerr(subr,nf90_def_dim(ncid,'nl1',nicas%blk(ib)%nl1,nl1_id))
         if (nicas%blk(ib)%nsa>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nsa',nicas%blk(ib)%nsa,nsa_id))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nsb',nicas%blk(ib)%nsb,nsb_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nsc',nicas%blk(ib)%nsc))
      end if
      if ((ib==bpar%nbe).and.nam%displ_diag) then
         call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0d',nicas%blk(ib)%nc0d,nc0d_id))
         call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0dinv',nicas%blk(ib)%nc0dinv,nc0dinv_id))
      end if

      ! Write main weight
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',nicas%blk(ib)%wgt))

      ! Define variables
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_def_var(ncid,'sb_to_c1b',nf90_int,(/nsb_id/),sb_to_c1b_id))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_def_var(ncid,'sb_to_l1',nf90_int,(/nsb_id/),sb_to_l1_id))
         if (nicas%blk(ib)%nsa>0) call mpl%ncerr(subr,nf90_def_var(ncid,'sa_to_sc',nf90_int,(/nsa_id/),sa_to_sc_id))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_def_var(ncid,'sb_to_sc',nf90_int,(/nsb_id/),sb_to_sc_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'norm',ncfloat,(/nc0a_id,nl0_id/),norm_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'coef_ens',ncfloat,(/nc0a_id,nl0_id/),coef_ens_id))

         if (nicas%blk(ib)%nsa>0) call mpl%ncerr(subr,nf90_put_att(ncid,sa_to_sc_id,'_FillValue',msvali))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_put_att(ncid,sb_to_sc_id,'_FillValue',msvali))
         call mpl%ncerr(subr,nf90_put_att(ncid,norm_id,'_FillValue',msvalr))
         call mpl%ncerr(subr,nf90_put_att(ncid,coef_ens_id,'_FillValue',msvalr))
      end if

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      if (bpar%nicas_block(ib)) then
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_c1b_id,nicas%blk(ib)%sb_to_c1b))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_l1_id,nicas%blk(ib)%sb_to_l1))
         if (nicas%blk(ib)%nsa>0) call mpl%ncerr(subr,nf90_put_var(ncid,sa_to_sc_id,nicas%blk(ib)%sa_to_sc))
         if (nicas%blk(ib)%nsb>0) call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_sc_id,nicas%blk(ib)%sb_to_sc))
         call mpl%ncerr(subr,nf90_put_var(ncid,norm_id,nicas%blk(ib)%norm))
         call mpl%ncerr(subr,nf90_put_var(ncid,coef_ens_id,nicas%blk(ib)%coef_ens))
         call nicas%blk(ib)%com_AB%write(mpl,ncid)
         call nicas%blk(ib)%com_AC%write(mpl,ncid)
         call nicas%blk(ib)%c%write(mpl,ncid)
         do il0i=1,geom%nl0i
            call nicas%blk(ib)%h(il0i)%write(mpl,ncid)
         end do
         call nicas%blk(ib)%v%write(mpl,ncid)
         do il1=1,nicas%blk(ib)%nl1
            call nicas%blk(ib)%s(il1)%write(mpl,ncid)
         end do
      end if
      if ((ib==bpar%nbe).and.nam%displ_diag) then
         call nicas%blk(ib)%com_AD%write(mpl,ncid)
         call nicas%blk(ib)%com_ADinv%write(mpl,ncid)
         do its=2,nam%nts
            do il0=1,geom%nl0
               call nicas%blk(ib)%d(il0,its)%write(mpl,ncid)
               call nicas%blk(ib)%dinv(il0,its)%write(mpl,ncid)
            end do
         end do
      end if

      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   end if
end do

end subroutine nicas_write

!----------------------------------------------------------------------
! Subroutine: nicas_write_mpi_summary
!> Purpose: write NICAS MPI related data summary
!----------------------------------------------------------------------
subroutine nicas_write_mpi_summary(nicas,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(mpl_type),intent(in) :: mpl      !< MPI data
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

do ib=1,bpar%nbe
   if (bpar%nicas_block(ib)) then
      ! Allocation
      allocate(lcheck(nicas%blk(ib)%nc1,nicas%blk(ib)%nl1))

      ! Create summary file
      filename = trim(nam%prefix)//'_'//trim(nicas%blk(ib)%name)//'_summary.nc'
      call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

      ! Write namelist parameters
      call nam%ncwrite(mpl,ncid)

      ! Define dimensions
      call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
      call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1',nicas%blk(ib)%nc1,nc1_id))
      call mpl%ncerr(subr,nf90_def_dim(ncid,'nl1',nicas%blk(ib)%nl1,nl1_id))

      ! Define variables
      call mpl%ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'c0_to_proc',nf90_int,(/nc0_id/),c0_to_proc_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'c1_to_c0',nf90_int,(/nc1_id/),c1_to_c0_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'l1_to_l0',nf90_int,(/nl1_id/),l1_to_l0_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lcheck',ncfloat,(/nc1_id,nl1_id/),lcheck_id))

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))

      ! Write variables
      call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
      call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
      call mpl%ncerr(subr,nf90_put_var(ncid,c0_to_proc_id,geom%c0_to_proc))
      call mpl%ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,nicas%blk(ib)%c1_to_c0))
      call mpl%ncerr(subr,nf90_put_var(ncid,l1_to_l0_id,nicas%blk(ib)%l1_to_l0))
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
      call mpl%ncerr(subr,nf90_put_var(ncid,lcheck_id,lcheck))

      ! Close summary file
      call mpl%ncerr(subr,nf90_close(ncid))

      ! Release memory
      deallocate(lcheck)
   end if
end do

end subroutine nicas_write_mpi_summary

!----------------------------------------------------------------------
! Subroutine: nicas_run_nicas
!> Purpose: NICAS driver
!----------------------------------------------------------------------
subroutine nicas_run_nicas(nicas,mpl,rng,nam,geom,bpar,cmat)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas  !< NICAS data
type(mpl_type),intent(inout) :: mpl       !< MPI data
type(rng_type),intent(inout) :: rng       !< Random number generator
type(nam_type),intent(inout) :: nam       !< Namelist
type(geom_type),intent(inout) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar        !< Block parameters
type(cmat_type),intent(in) :: cmat        !< C matrix data

! Local variables
integer :: ib

! Allocation
call nicas%alloc(mpl,nam,bpar,'nicas')

! Compute NICAS parameters
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute NICAS parameters'
call flush(mpl%unit)

do ib=1,bpar%nbe
   if (bpar%nicas_block(ib).or.((ib==bpar%nbe).and.nam%displ_diag)) then
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))
      call flush(mpl%unit)
   end if

   ! NICAS parameters
   if (bpar%nicas_block(ib)) call nicas%blk(ib)%compute_parameters(mpl,rng,nam,geom,cmat%blk(ib))

   ! Advection
   if ((ib==bpar%nbe).and.nam%displ_diag) call nicas%blk(ib)%compute_adv(mpl,rng,nam,geom,cmat%blk(ib))

   ! Coefficient
   if (bpar%B_block(ib)) then
      ! Copy weights
      nicas%blk(ib)%wgt = cmat%blk(ib)%wgt
      if (bpar%nicas_block(ib)) then
         allocate(nicas%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
         nicas%blk(ib)%coef_ens = cmat%blk(ib)%coef_ens
      end if
   end if
end do

! Write NICAS parameters
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Write NICAS parameters'
call flush(mpl%unit)
call nicas%write(mpl,nam,geom,bpar)

! Write NICAS MPI summary
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Write NICAS MPI summary'
call flush(mpl%unit)
call nicas%write_mpi_summary(mpl,nam,geom,bpar)

end subroutine nicas_run_nicas

!----------------------------------------------------------------------
! Subroutine: nicas_run_nicas_tests
!> Purpose: NICAS tests driver
!----------------------------------------------------------------------
subroutine nicas_run_nicas_tests(nicas,mpl,rng,nam,geom,bpar,io,cmat,ens)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas  !< NICAS data
type(mpl_type),intent(inout) :: mpl       !< MPI data
type(rng_type),intent(inout) :: rng       !< Random number generator
type(nam_type),intent(inout) :: nam       !< Namelist
type(geom_type),intent(inout) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar        !< Block parameters
type(io_type),intent(in) :: io            !< I/O
type(cmat_type),intent(in) :: cmat        !< C matrix data
type(ens_type),intent(in) :: ens          !< Ensemble

! Local variables
integer :: ib

if (nam%check_adjoints) then
   ! Test adjoint
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS adjoint'
   call flush(mpl%unit)

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))
         call flush(mpl%unit)
         call nicas%blk(ib)%test_adjoint(mpl,rng,nam,geom)
      end if
   end do

   ! Test NICAS adjoint
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS adjoint'
   call flush(mpl%unit)
   call nicas%test_adjoint(mpl,rng,nam,geom,bpar,ens)
end if

if (nam%check_pos_def) then
   ! Test NICAS positive definiteness
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS positive definiteness'
   call flush(mpl%unit)

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))
         call flush(mpl%unit)
         call nicas%blk(ib)%test_pos_def(mpl,rng,nam,geom)
      end if
   end do
end if

if (nam%check_sqrt) then
   ! Test NICAS full/square-root equivalence
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS full/square-root equivalence'
   call flush(mpl%unit)

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))
         call flush(mpl%unit)
         call nicas%blk(ib)%test_sqrt(mpl,rng,nam,geom,bpar,io,cmat%blk(ib))
      end if
   end do

   ! Test NICAS full/square-root equivalence
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS full/square-root equivalence'
   call flush(mpl%unit)
   call nicas%test_sqrt(mpl,rng,nam,geom,bpar,io,cmat,ens)
end if

if (nam%check_dirac) then
   ! Apply NICAS to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'
   call flush(mpl%unit)

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Block: '//trim(bpar%blockname(ib))
         call flush(mpl%unit)
         call nicas%blk(ib)%test_dirac(mpl,nam,geom,bpar,io)
      end if
   end do

   ! Apply NICAS to diracs
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'
   call flush(mpl%unit)
   call nicas%test_dirac(mpl,nam,geom,bpar,io,ens)
end if

if (nam%check_randomization) then
   ! Test NICAS randomization
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS randomization'
   call flush(mpl%unit)
   call nicas%test_randomization(mpl,rng,nam,geom,bpar,io)
end if

if (nam%check_consistency) then
   ! Test HDIAG-NICAS consistency
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test HDIAG-NICAS consistency'
   call flush(mpl%unit)
   call nicas%test_consistency(mpl,rng,nam,geom,bpar,io,cmat)
end if

if (nam%check_optimality) then
   ! Test HDIAG optimality
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test HDIAG optimality'
   call flush(mpl%unit)
   call nicas%test_optimality(mpl,rng,nam,geom,bpar,io)
end if

end subroutine nicas_run_nicas_tests

!----------------------------------------------------------------------
! Subroutine: nicas_alloc_cv
!> Purpose: control vector allocation
!----------------------------------------------------------------------
subroutine nicas_alloc_cv(nicas,bpar,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(cv_type),intent(inout) :: cv     !< Control vector

! Local variables
integer :: ib

! Allocation
allocate(cv%blk(bpar%nbe))

do ib=1,bpar%nbe
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
subroutine nicas_random_cv(nicas,rng,bpar,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(rng_type),intent(inout) :: rng   !< Random number generator
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(cv_type),intent(out) :: cv       !< Control vector

! Local variables
integer :: ib

! Allocation
call nicas%alloc_cv(bpar,cv)

! Random initialization
do ib=1,bpar%nbe
   if (bpar%cv_block(ib)) call rng%rand_gau(cv%blk(ib)%alpha)
end do

end subroutine nicas_random_cv

!----------------------------------------------------------------------
! Subroutine: nicas_apply
!> Purpose: apply NICAS
!----------------------------------------------------------------------
subroutine nicas_apply(nicas,mpl,nam,geom,bpar,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   !< NICAS data
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
type(bpar_type),target,intent(in) :: bpar                               !< Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ib,its,iv,jv,il0,ic0a
real(kind_real) :: prod,prod_tot
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:)
real(kind_real),allocatable :: fld_save(:,:,:,:)

if (pos_def_test) then
   ! Save field for positive-definiteness test
   allocate(fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   fld_save = fld
end if

! Adjoint advection
if (nam%advmode==1) call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

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

   ! Apply common ensemble coefficient square-root
   !$omp parallel do schedule(static) private(il0,ic0a)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
      end do
   end do
   !$omp end parallel do

   ! Apply common NICAS
   call nicas%blk(bpar%nbe)%apply(mpl,nam,geom,fld_3d)

   ! Apply common ensemble coefficient square-root
   !$omp parallel do schedule(static) private(il0,ic0a)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
      end do
   end do
   !$omp end parallel do

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
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do

         ! Apply specific NICAS (same for all timeslots)
         call nicas%blk(ib)%apply(mpl,nam,geom,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do
case ('specific_multivariate')
   call mpl%abort('specific multivariate strategy should not be called from apply_NICAS (lsqrt required)')
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
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                         & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do

      ! Apply common NICAS
      call nicas%blk(bpar%nbe)%apply(mpl,nam,geom,fld_4d(:,:,iv))

      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                         & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
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

if (pos_def_test) then
   ! Positive-definiteness test
   prod = sum(fld_save*fld)
   call mpl%allreduce_sum(prod,prod_tot)
   if (prod_tot<0.0) call mpl%abort('negative result in nicas_apply')
end if

! Advection
if (nam%advmode==1) call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld)

end subroutine nicas_apply

!----------------------------------------------------------------------
! Subroutine: nicas_apply_from_sqrt
!> Purpose: apply NICAS from square-root
!----------------------------------------------------------------------
subroutine nicas_apply_from_sqrt(nicas,mpl,nam,geom,bpar,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   !< NICAS data
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
type(bpar_type),target,intent(in) :: bpar                               !< Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
real(kind_real) :: prod,prod_tot
real(kind_real),allocatable :: fld_save(:,:,:,:)
type(cv_type) :: cv

if (pos_def_test) then
   ! Save field for positivity test
   allocate(fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   fld_save = fld
end if

! Apply square-root adjoint
call nicas%apply_sqrt_ad(mpl,nam,geom,bpar,fld,cv)

! Apply square-root
call nicas%apply_sqrt(mpl,nam,geom,bpar,cv,fld)

if (pos_def_test) then
   ! Positivity test
   prod = sum(fld_save*fld)
   call mpl%allreduce_sum(prod,prod_tot)
   if (prod_tot<0.0) call mpl%abort('negative result in nicas_apply')
end if

end subroutine nicas_apply_from_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_apply_sqrt
!> Purpose: apply NICAS square-root
!----------------------------------------------------------------------
subroutine nicas_apply_sqrt(nicas,mpl,nam,geom,bpar,cv,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                 !< NICAS data
type(mpl_type),intent(in) :: mpl                                      !< MPI data
type(nam_type),target,intent(in) :: nam                               !< Namelist
type(geom_type),target,intent(in) :: geom                             !< Geometry
type(bpar_type),target,intent(in) :: bpar                             !< Block parameters
type(cv_type),intent(in) :: cv                                        !< Control variable
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ib,its,iv,jv,i,ic0a,il0
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),a(:),u(:)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Apply common NICAS
   call nicas%blk(bpar%nbe)%apply_sqrt(mpl,geom,cv%blk(bpar%nbe)%alpha,fld_3d)

   ! Apply common ensemble coefficient square-root
   !$omp parallel do schedule(static) private(il0,ic0a)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
      end do
   end do
   !$omp end parallel do

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

         ! Apply specific NICAS (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt(mpl,geom,cv%blk(ib)%alpha,fld_4d(:,:,iv))

         ! Apply specific ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
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

         ! Apply specific NICAS (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt(mpl,geom,cv%blk(bpar%nbe)%alpha,fld_4d(:,:,iv))

         ! Apply specific ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
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
   call cholesky(mpl,nam%nv,nam%nv*(nam%nv+1)/2,a,u)
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

         ! Apply specific NICAS (same for all timeslots)
         call nicas%blk(bpar%nbe)%apply_sqrt(mpl,geom,cv%blk(ib)%alpha,fld_4d(:,:,iv))

         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
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

! Advection
if (nam%advmode==1) call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld)

end subroutine nicas_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_apply_sqrt_ad
!> Purpose: apply NICAS square-root, adjoint
!----------------------------------------------------------------------
subroutine nicas_apply_sqrt_ad(nicas,mpl,nam,geom,bpar,fld,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                !< NICAS data
type(mpl_type),intent(in) :: mpl                                     !< MPI data
type(nam_type),target,intent(in) :: nam                              !< Namelist
type(geom_type),target,intent(in) :: geom                            !< Geometry
type(bpar_type),target,intent(in) :: bpar                            !< Block parameters
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field
type(cv_type),intent(out) :: cv                                      !< Control variable

! Local variable
integer :: ib,its,iv,jv,i,ic0a,il0
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:),fld_5d(:,:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),a(:),u(:)
type(cv_type) :: cv_tmp

! Allocation
allocate(fld_5d(geom%nc0a,geom%nl0,nam%nv,nam%nts))
call nicas%alloc_cv(bpar,cv)

! Copy
fld_5d = fld

! Adjoint advection
if (nam%advmode==1) call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld_5d)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Sum product over variables and timeslots
   fld_3d = 0.0
   do its=1,nam%nts
      do iv=1,nam%nv
         fld_3d = fld_3d+fld_5d(:,:,iv,its)
      end do
   end do

   ! Apply common ensemble coefficient square-root
   !$omp parallel do schedule(static) private(il0,ic0a)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
      end do
   end do
   !$omp end parallel do

   ! Apply common NICAS
   call nicas%blk(bpar%nbe)%apply_sqrt_ad(mpl,geom,fld_3d,cv%blk(bpar%nbe)%alpha)
case ('specific_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld_5d(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do

         ! Apply specific NICAS (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt_ad(mpl,geom,fld_4d(:,:,iv),cv%blk(ib)%alpha)
      end if
   end do
case ('specific_multivariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   call nicas%alloc_cv(bpar,cv_tmp)

   ! Initialization
   cv%blk(bpar%nbe)%alpha = 0.0

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld_5d(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do

         ! Apply specific NICAS (same for all timeslots)
         call nicas%blk(ib)%apply_sqrt_ad(mpl,geom,fld_4d(:,:,iv),cv_tmp%blk(bpar%nbe)%alpha)

         ! Sum control variable
         cv%blk(bpar%nbe)%alpha = cv%blk(bpar%nbe)%alpha+cv_tmp%blk(bpar%nbe)%alpha
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
   call cholesky(mpl,nam%nv,nam%nv*(nam%nv+1)/2,a,u)
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
      fld_4d = fld_4d+fld_5d(:,:,:,its)
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
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) fld_4d_tmp(ic0a,il0,iv) = fld_4d_tmp(ic0a,il0,iv) &
                                                                            & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do

         ! Apply specific NICAS (same for all timeslots)
         call nicas%blk(bpar%nbe)%apply_sqrt_ad(mpl,geom,fld_4d_tmp(:,:,iv),cv%blk(ib)%alpha)
      end if
   end do
end select

end subroutine nicas_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: nicas_randomize
!> Purpose: randomize NICAS from square-root
!----------------------------------------------------------------------
subroutine nicas_randomize(nicas,mpl,rng,nam,geom,bpar,ne,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas      !< NICAS data
type(mpl_type),intent(in) :: mpl           !< MPI data
type(rng_type),intent(inout) :: rng        !< Random number generator
type(nam_type),target,intent(in) :: nam    !< Namelist
type(geom_type),target,intent(in) :: geom  !< Geometry
type(bpar_type),target,intent(in) :: bpar  !< Blocal parameters
integer,intent(in) :: ne                   !< Number of members
type(ens_type),intent(out) :: ens          !< Ensemble

! Local variable
integer :: ie,ic0a,il0,its,iv
real(kind_real) :: norm,mean(geom%nc0a,geom%nl0,nam%nv,nam%nts),std(geom%nc0a,geom%nl0,nam%nv,nam%nts)
type(cv_type) :: cv_ens(ne)

! Allocation
call ens%alloc(nam,geom,ne,1)

do ie=1,ne
   ! Generate random control vector
   call nicas%random_cv(rng,bpar,cv_ens(ie))

   ! Apply square-root
   call nicas%apply_sqrt(mpl,nam,geom,bpar,cv_ens(ie),ens%fld(:,:,:,:,ie))
end do

! Remove mean
mean = sum(ens%fld,dim=5)/real(ne,kind_real)
do ie=1,ne
   ens%fld(:,:,:,:,ie) = ens%fld(:,:,:,:,ie)-mean
end do

! Compute standard deviation
norm = real(ne-1,kind_real)
!$omp parallel do schedule(static) private(its,iv,il0,ic0a)
do its=1,nam%nts
   do iv=1,nam%nv
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%mask(geom%c0a_to_c0(ic0a),il0)) std(ic0a,il0,iv,its) = sqrt(sum(ens%fld(ic0a,il0,iv,its,:)**2)) &
                                                                          & /norm
         end do
      end do
   end do
end do
!$omp end parallel do

! Normalize perturbations
!$omp parallel do schedule(static) private(its,iv,il0,ic0a)
do its=1,nam%nts
   do iv=1,nam%nv
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%mask(geom%c0a_to_c0(ic0a),il0)) ens%fld(ic0a,il0,iv,its,:) = ens%fld(ic0a,il0,iv,its,:) &
                                                                                & /std(ic0a,il0,iv,its)
         end do
      end do
   end do
end do
!$omp end parallel do

end subroutine nicas_randomize

!----------------------------------------------------------------------
! Subroutine: nicas_apply_bens
!> Purpose: apply localized ensemble covariance
!----------------------------------------------------------------------
subroutine nicas_apply_bens(nicas,mpl,nam,geom,bpar,ens,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   !< NICAS data
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
type(bpar_type),target,intent(in) :: bpar                               !< Blocal parameters
type(ens_type),intent(in) :: ens                                        !< Ensemble
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ie
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv,nam%nts),fld_tmp(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: pert(geom%nc0a,geom%nl0,nam%nv,nam%nts)

! Copy field
fld_copy = fld

! Adjoint advection
if (nam%advmode==-1) call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld_copy)

! Apply localized ensemble covariance formula
fld = 0.0
do ie=1,nam%ens1_ne
   ! Compute perturbation
   pert = ens%fld(:,:,:,:,ie)

   ! Inverse advection
   if (nam%advmode==-1) call nicas%blk(bpar%nbe)%apply_adv_inv(mpl,nam,geom,pert)

   ! Schur product
   fld_tmp = pert*fld_copy

   ! Apply NICAS
   if (nam%lsqrt) then
      call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_tmp)
   else
      call nicas%apply(mpl,nam,geom,bpar,fld_tmp)
   end if

   ! Schur product
   fld = fld+fld_tmp*pert

   ! Normalization
   fld = fld/real(nam%ens1_ne-1,kind_real)
end do

! Advection
if (nam%advmode==-1) call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld)

end subroutine nicas_apply_bens

!----------------------------------------------------------------------
! Subroutine: nicas_apply_bens_noloc
!> Purpose: apply ensemble covariance, without localization
!----------------------------------------------------------------------
subroutine nicas_apply_bens_noloc(nicas,mpl,nam,geom,ens,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   !< NICAS data
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
type(ens_type),intent(in) :: ens                                        !< Ensemble
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variable
integer :: ie,ic0a,il0,iv,its
real(kind_real) :: alpha,norm
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: pert(geom%nc0a,geom%nl0,nam%nv,nam%nts)
character(len=1024) :: dum

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
               if (geom%mask(geom%c0a_to_c0(ic0a),il0)) pert(ic0a,il0,iv,its) = ens%fld(ic0a,il0,iv,its,ie)/norm
            end do
         end do
      end do
   end do
   !$omp end parallel do

   ! Dot product
   call mpl%dot_prod(pert,fld_copy,alpha)

   ! Schur product
   fld = fld+alpha*pert
end do

! To avoid compilation warnings
dum = nicas%prefix

end subroutine nicas_apply_bens_noloc

!----------------------------------------------------------------------
! Subroutine: nicas_test_adjoint
!> Purpose: test NICAS adjoint
!----------------------------------------------------------------------
subroutine nicas_test_adjoint(nicas,mpl,rng,nam,geom,bpar,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas     !< NICAS data
type(mpl_type),intent(in) :: mpl          !< MPI data
type(rng_type),intent(inout) :: rng       !< Random number generator
type(nam_type),intent(in) :: nam          !< Namelist
type(geom_type),intent(in) :: geom        !< Geometry
type(bpar_type),intent(in) :: bpar        !< Block parameters
type(ens_type),intent(in) :: ens          !< Ensemble

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real),allocatable :: fld1_loc(:,:,:,:),fld1_adv(:,:,:,:),fld1_bens(:,:,:,:),fld1_save(:,:,:,:)
real(kind_real),allocatable :: fld2_loc(:,:,:,:),fld2_adv(:,:,:,:),fld2_bens(:,:,:,:),fld2_save(:,:,:,:)

! Allocation
allocate(fld1_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld2_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld1_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld2_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (ens%ne>0) then
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
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld1_loc)
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld2_loc)
else
   call nicas%apply(mpl,nam,geom,bpar,fld1_loc)
   call nicas%apply(mpl,nam,geom,bpar,fld2_loc)
end if
if (abs(nam%advmode)==1) then
   fld1_adv = fld1_save
   fld2_adv = fld2_save
   call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld1_adv)
   call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld2_adv)
end if
if (ens%ne>0) then
   fld1_bens = fld1_save
   fld2_bens = fld2_save
   call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld1_bens)
   call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld2_bens)
end if

! Print result
call mpl%dot_prod(fld1_loc,fld2_save,sum1)
call mpl%dot_prod(fld2_loc,fld1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%unit)
if (abs(nam%advmode)==1) then
   call mpl%dot_prod(fld1_adv,fld2_save,sum1)
   call mpl%dot_prod(fld2_adv,fld1_save,sum2)
   write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Advection adjoint test:    ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
   call flush(mpl%unit)
end if
if (ens%ne>0) then
   call mpl%dot_prod(fld1_bens,fld2_save,sum1)
   call mpl%dot_prod(fld2_bens,fld1_save,sum2)
   write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Ensemble B adjoint test:   ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
   call flush(mpl%unit)
end if

end subroutine nicas_test_adjoint

!----------------------------------------------------------------------
! Subroutine: nicas_test_sqrt
!> Purpose: test full/square-root equivalence
!----------------------------------------------------------------------
subroutine nicas_test_sqrt(nicas,mpl,rng,nam,geom,bpar,io,cmat,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas      !< NICAS data
type(mpl_type),intent(inout) :: mpl        !< MPI data
type(rng_type),intent(inout) :: rng        !< Random number generator
type(nam_type),intent(inout),target :: nam !< Namelist
type(geom_type),intent(in),target :: geom  !< Geometry
type(bpar_type),intent(in) :: bpar         !< Block parameters
type(io_type),intent(in) :: io             !< I/O
type(cmat_type),intent(in) :: cmat         !< C matrix data
type(ens_type),intent(in) :: ens           !< Ensemble

! Local variables
integer :: ib,iv
real(kind_real),allocatable :: fld_loc(:,:,:,:),fld_loc_sqrt(:,:,:,:)
real(kind_real),allocatable :: fld_bens(:,:,:,:),fld_bens_sqrt(:,:,:,:)
character(len=1024) :: varname(nam%nv)
type(nicas_type) :: nicas_other

! Allocation
allocate(fld_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(fld_loc_sqrt(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (ens%ne>0) then
   allocate(fld_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   allocate(fld_bens_sqrt(geom%nc0a,geom%nl0,nam%nv,nam%nts))
end if

! Generate random field
call rng%rand_real(-1.0_kind_real,1.0_kind_real,fld_loc)
fld_loc_sqrt = fld_loc
if (ens%ne>0) then
   call rng%rand_real(-1.0_kind_real,1.0_kind_real,fld_bens)
   fld_bens_sqrt = fld_bens
end if

! Apply NICAS, initial version
if (nam%lsqrt) then
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_loc_sqrt)
else
   call nicas%apply(mpl,nam,geom,bpar,fld_loc)
end if
if (ens%ne>0) then
   if (nam%lsqrt) then
      call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld_bens_sqrt)
   else
      call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld_bens)
   end if
end if

! Switch lsqrt
nam%lsqrt = .not.nam%lsqrt

! Allocation
call nicas_other%alloc(mpl,nam,bpar,'nicas_other')

! Prepare nicas, other version
do ib=1,bpar%nbe
   if (bpar%nicas_block(ib)) then
      ! Compute NICAS parameters
      call nicas_other%blk(ib)%compute_parameters(mpl,rng,nam,geom,cmat%blk(ib))
   end if

   if (bpar%B_block(ib)) then
      ! Copy weights
      nicas_other%blk(ib)%wgt = cmat%blk(ib)%wgt
      if (bpar%nicas_block(ib)) then
         allocate(nicas_other%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
         nicas_other%blk(ib)%coef_ens = cmat%blk(ib)%coef_ens
      end if
   end if
end do

! Apply NICAS, other version
if (nam%lsqrt) then
   call nicas_other%apply_from_sqrt(mpl,nam,geom,bpar,fld_loc_sqrt)
else
   call nicas_other%apply(mpl,nam,geom,bpar,fld_loc)
end if
if (ens%ne>0) then
   if (nam%lsqrt) then
      call nicas_other%apply_bens(mpl,nam,geom,bpar,ens,fld_bens_sqrt)
   else
      call nicas_other%apply_bens(mpl,nam,geom,bpar,ens,fld_bens)
   end if
end if

! Compute dirac
do iv=1,nam%nv
   varname(iv) = nam%varname(iv)
   nam%varname(iv) = trim(varname(iv))//'_sqrt'
end do
if (nam%check_dirac) call nicas_other%test_dirac(mpl,nam,geom,bpar,io,ens)
do iv=1,nam%nv
   nam%varname(iv) = varname(iv)
end do

! Reset lsqrt value
nam%lsqrt = .not.nam%lsqrt

! Print difference
write(mpl%unit,'(a7,a,f6.1,a)') '','NICAS full / square-root error : ', &
 & sqrt(sum((fld_loc_sqrt-fld_loc)**2)/sum(fld_loc**2))*100.0,'%'
if (ens%ne>0) write(mpl%unit,'(a7,a,f6.1,a)') '','Ensemble B full / square-root error:  ', &
 & sqrt(sum((fld_bens_sqrt-fld_bens)**2)/sum(fld_bens**2))*100.0,'%'
call flush(mpl%unit)

end subroutine nicas_test_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_test_dirac
!> Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine nicas_test_dirac(nicas,mpl,nam,geom,bpar,io,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas     !< NICAS data
type(mpl_type),intent(inout) :: mpl       !< MPI data
type(nam_type),intent(in) :: nam          !< Namelist
type(geom_type),intent(in) :: geom        !< Geometry
type(bpar_type),intent(in) :: bpar        !< Block parameters
type(io_type),intent(in) :: io            !< I/O
type(ens_type),intent(in) :: ens          !< Ensemble

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
if (ens%ne>0) allocate(fld_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))

! Apply NICAS to dirac
fld_loc = fld
if (nam%lsqrt) then
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_loc)
else
   call nicas%apply(mpl,nam,geom,bpar,fld_loc)
end if

if (ens%ne>0) then
   ! Apply localized ensemble covariance
   fld_bens = fld
   call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld_bens)
end if

! Write field
filename = trim(nam%prefix)//'_dirac'
do its=1,nam%nts
   write(itschar,'(i2.2)') its
   do iv=1,nam%nv
      call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_'//itschar,fld_loc(:,:,iv,its))
      if (ens%ne>0) call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_'//itschar//'_Bens', &
       & fld_bens(:,:,iv,its))
   end do
end do

end subroutine nicas_test_dirac

!----------------------------------------------------------------------
! Subroutine: nicas_test_randomization
!> Purpose: test NICAS randomization method with respect to theoretical error statistics
!----------------------------------------------------------------------
subroutine nicas_test_randomization(nicas,mpl,rng,nam,geom,bpar,io)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(mpl_type),intent(inout) :: mpl   !< MPI data
type(rng_type),intent(inout) :: rng   !< Random number generator
type(nam_type),intent(inout) :: nam   !< Namelist variables
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(io_type),intent(in) :: io        !< I/O

! Local variables
integer :: ifac,itest,nefac(nfac),ens1_ne,iv,its
real(kind_real) :: fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest),fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest)
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts),mse(ntest,nfac),mse_th(ntest,nfac)
character(len=2) :: itschar
character(len=4) :: nechar,itestchar
character(len=1024) :: filename
type(ens_type) :: ens

! Define test vectors
write(mpl%unit,'(a4,a)') '','Define test vectors'
call flush(mpl%unit)
call define_test_vectors(mpl,rng,nam,geom,ntest,fld_save)

! Apply NICAS to test vectors
write(mpl%unit,'(a4,a)') '','Apply NICAS to test vectors'
call flush(mpl%unit)
fld_ref = fld_save
do itest=1,ntest
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_ref(:,:,:,:,itest))
end do

! Write first 10 test vectors
do itest=1,min(ntest,10)
   ! Write field
   write(itestchar,'(i4.4)') itest
   filename = trim(nam%prefix)//'_randomize_'//itestchar
   do its=1,nam%nts
      write(itschar,'(i2.2)') its
      do iv=1,nam%nv
         call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_ref_'//itschar,fld_ref(:,:,iv,its,itest))
      end do
   end do
end do

! Save namelist variables
ens1_ne = nam%ens1_ne

write(mpl%unit,'(a4,a)') '','Test randomization for various ensemble sizes:'
call flush(mpl%unit)
do ifac=1,nfac
   ! Ensemble size
   nefac(ifac) = max(int(5.0*real(ifac,kind_real)/real(nfac,kind_real)*real(ne_rand,kind_real)),3)
   nam%ens1_ne = nefac(ifac)
   write(nechar,'(i4.4)') nefac(ifac)

   ! Randomize ensemble
   call nicas%randomize(mpl,rng,nam,geom,bpar,nefac(ifac),ens)

   do itest=1,ntest
      ! Test NICAS
      fld = fld_save(:,:,:,:,itest)
      call nicas%apply_bens_noloc(mpl,nam,geom,ens,fld)

      ! RMSE
      mse(itest,ifac) = sum((fld-fld_ref(:,:,:,:,itest))**2)
      mse_th(itest,ifac) = 1.0/real(nam%ens1_ne-1,kind_real)*sum(1+fld_ref(:,:,:,:,itest)**2)

      ! Write first 10 test vectors
      if (itest<=min(ntest,10)) then
         ! Write field
         write(itestchar,'(i4.4)') itest
         filename = trim(nam%prefix)//'_randomize_'//itestchar
         do its=1,nam%nts
            write(itschar,'(i2.2)') its
            do iv=1,nam%nv
               call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_rand_'//nechar//'_'//itschar,fld(:,:,iv,its))
            end do
         end do
      end if
   end do

   ! Print scores
   write(mpl%unit,'(a7,a,i4,a,e15.8,a,e15.8)') '','Ensemble size ',nefac(ifac),', MSE (exp. / th.): ', &
 & sum(mse(:,ifac))/real(ntest,kind_real),' / ',sum(mse_th(:,ifac))/real(ntest,kind_real)
   call flush(mpl%unit)

   ! Release memory
   call ens%dealloc
end do

! Reset namelist variables
nam%ens1_ne = ens1_ne

end subroutine nicas_test_randomization

!----------------------------------------------------------------------
! Subroutine: nicas_test_consistency
!> Purpose: test HDIAG-NICAS consistency with a randomization method
!----------------------------------------------------------------------
subroutine nicas_test_consistency(nicas,mpl,rng,nam,geom,bpar,io,cmat)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(mpl_type),intent(inout) :: mpl   !< MPI data
type(rng_type),intent(inout) :: rng   !< Random number generator
type(nam_type),intent(inout) :: nam   !< Namelist variables
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(io_type),intent(in) :: io        !< I/O
type(cmat_type),intent(in) :: cmat    !< C matrix data

! Local variables
integer :: ens1_ne,ens1_ne_offset,ens1_nsub,ib,il0
real(kind_real) :: rh_c0_sum,rv_c0_sum,norm
character(len=1024) :: prefix,method
type(cmat_type) :: cmat_test
type(ens_type) :: ens

! Randomize ensemble
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Randomize ensemble'
call flush(mpl%unit)
call nicas%randomize(mpl,rng,nam,geom,bpar,ne_rand,ens)

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
call cmat_test%run_hdiag(mpl,rng,nam,geom,bpar,io,ens)

! Print scores
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- HDIAG/NICAS consistency results'
do ib=1,bpar%nbe
   if (bpar%nicas_block(ib)) then
      write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      do il0=1,geom%nl0
         call mpl%allreduce_sum(sum(cmat_test%blk(ib)%rh(:,il0)-cmat%blk(ib)%rh(:,il0),geom%mask(geom%c0a_to_c0,il0)), &
       & rh_c0_sum)
         call mpl%allreduce_sum(sum(cmat_test%blk(ib)%rv(:,il0)-cmat%blk(ib)%rv(:,il0),geom%mask(geom%c0a_to_c0,il0)), &
       & rv_c0_sum)
         call mpl%allreduce_sum(real(count(geom%mask(geom%c0a_to_c0,il0)),kind_real),norm)
         write(mpl%unit,'(a10,a7,i3,a4,a25,f6.1,a)') '','Level: ',nam%levs(il0),' ~> ','horizontal length-scale: ', &
       & rh_c0_sum/norm*reqkm,' km'
         if (any(abs(cmat%blk(ib)%rv(:,il0))>0.0)) then
            write(mpl%unit,'(a49,f6.1,a)') 'vertical length-scale: ',rh_c0_sum/norm,' '//trim(mpl%vunitchar)
         end if
      end do
   end if
end do
call flush(mpl%unit)

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
subroutine nicas_test_optimality(nicas,mpl,rng,nam,geom,bpar,io)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas !< NICAS data
type(mpl_type),intent(inout) :: mpl   !< MPI data
type(rng_type),intent(inout) :: rng   !< Random number generator
type(nam_type),intent(inout) :: nam   !< Namelist variables
type(geom_type),intent(in) :: geom    !< Geometry
type(bpar_type),intent(in) :: bpar    !< Block parameters
type(io_type),intent(in) :: io        !< I/O

! Local variables
integer :: ib,ifac,itest
real(kind_real) :: fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest),fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest)
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts),fac(nfac),mse(ntest,nfac)
character(len=1024) :: prefix,method
type(cmat_type) :: cmat_save,cmat_test
type(ens_type) :: ens
type(nicas_type) :: nicas_test

! Define test vectors
write(mpl%unit,'(a4,a)') '','Define test vectors'
call flush(mpl%unit)
call define_test_vectors(mpl,rng,nam,geom,ntest,fld_save)

! Apply NICAS to test vectors
write(mpl%unit,'(a4,a)') '','Apply NICAS to test vectors'
call flush(mpl%unit)
fld_ref = fld_save
do itest=1,ntest
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_ref(:,:,:,:,itest))
end do

! Randomize ensemble
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Randomize ensemble'
call flush(mpl%unit)
call nicas%randomize(mpl,rng,nam,geom,bpar,nam%ens1_ne,ens)

! Copy sampling
call system('cp -f '//trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc ' &
 & //trim(nam%datadir)//'/'//trim(nam%prefix)//'_optimality-test_sampling.nc')

! Save namelist variables
prefix = nam%prefix
method = nam%method

! Set namelist variables
nam%prefix = trim(nam%prefix)//'_optimality-test'
nam%method = 'loc_norm'

! Allocation
call nicas_test%alloc(mpl,nam,bpar,'nicas_test')
call cmat_save%alloc(nam,geom,bpar,'cmat_save')

! Call hdiag driver
call cmat_save%run_hdiag(mpl,rng,nam,geom,bpar,io,ens)

! Copy cmat
cmat_test = cmat_save%copy(nam,geom,bpar)

do ifac=1,nfac
   ! Multiplication factor
   fac(ifac) = 2.0*real(ifac,kind_real)/real(nfac,kind_real)

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,f4.2,a)') '--- Apply a multiplicative factor ',fac(ifac),' to length-scales'
   call flush(mpl%unit)

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         ! Length-scales multiplication
         cmat_test%blk(ib)%rh = fac(ifac)*cmat_save%blk(ib)%rh
         cmat_test%blk(ib)%rv = fac(ifac)*cmat_save%blk(ib)%rv
         if (trim(nam%strategy)=='specific_multivariate') then
            cmat_test%blk(ib)%rhs = fac(ifac)*cmat_save%blk(ib)%rhs
            cmat_test%blk(ib)%rvs = fac(ifac)*cmat_save%blk(ib)%rvs
         end if

         ! Compute NICAS parameters
         call nicas_test%blk(ib)%compute_parameters(mpl,rng,nam,geom,cmat_test%blk(ib))
      end if

      if (bpar%B_block(ib)) then
         ! Copy weights
         nicas_test%blk(ib)%wgt = cmat_test%blk(ib)%wgt
         if (bpar%nicas_block(ib)) then
            allocate(nicas_test%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
            nicas_test%blk(ib)%coef_ens = cmat_test%blk(ib)%coef_ens
         end if
      end if
   end do

   do itest=1,ntest
      ! Test NICAS
      fld = fld_save(:,:,:,:,itest)
      call nicas_test%apply_bens(mpl,nam,geom,bpar,ens,fld)

      ! RMSE
      mse(itest,ifac) = sum((fld-fld_ref(:,:,:,:,itest))**2)
   end do

   ! Print scores
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,f4.2,a,e15.8)') '--- Optimality results for a factor ',fac(ifac),', MSE: ', &
 & sum(mse(:,ifac))/real(ntest,kind_real)
   call flush(mpl%unit)
end do

! Print scores summary
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Optimality results summary'
do ifac=1,nfac
   write(mpl%unit,'(a7,a,f4.2,a,e15.8)') '','Factor ',fac(ifac),', MSE: ',sum(mse(:,ifac))/real(ntest,kind_real)
end do
call flush(mpl%unit)

! Reset namelist variables
nam%prefix = prefix
nam%method = method

end subroutine nicas_test_optimality

end module type_nicas
