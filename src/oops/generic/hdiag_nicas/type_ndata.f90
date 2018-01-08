!----------------------------------------------------------------------
! Module: type_ndata
!> Purpose: sampling data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_ndata

use netcdf
use tools_const, only: rad2deg
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr
use tools_nc, only: ncerr,ncfloat
use type_com, only: comtype,com_dealloc,com_read,com_write
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_read,linop_write
use type_mpl, only: mpl
use type_nam, only: namtype,namncwrite

implicit none

! Sampling data derived type
type ndatatype
   ! Block name
   character(len=1024) :: cname                 !< Block name

   ! Namelist
   type(namtype),pointer :: nam                 !< Namelist

   ! Geometry
   type(geomtype),pointer :: geom               !< Geometry

   ! Specific geometry
   integer :: nc1                               !< Number of points in subset Sc1
   integer,allocatable :: vbot(:)               !< Bottom level in grid Gh
   integer,allocatable :: vtop(:)               !< Top level in grid Gh
   integer,allocatable :: nc2(:)                !< Number of points in subset Sc2
   integer :: ns                                !< Number of subgrid nodes

   ! Full interpolations
   type(linoptype),allocatable :: hfull(:)      !< Horizontal interpolation
   type(linoptype),allocatable :: sfull(:)      !< Subsample interpolation

   ! Other data

   ! Level selection
   logical,allocatable :: llev(:)               !< Level selection

   ! Parameters/normalization conversion
   integer,allocatable :: s_to_c1(:)          !< Subgrid to subset Sc1
   integer,allocatable :: s_to_l1(:)          !< Subgrid to subset Sl1
   integer,allocatable :: s_to_c2(:)          !< Subgrid to subset sc2
   integer,allocatable :: c1_to_c0(:)         !< Subset Sc1 to subset Sc0
   integer,allocatable :: l1_to_l0(:)         !< Subset Sl1 to subset Sl0
   integer,allocatable :: c2l1_to_c0(:,:)    !< Grid Gs to subset Sc0
   integer,allocatable :: c0_to_c1(:)         !< Subset Sc0 to subset Sc1
   integer,allocatable :: l0_to_l1(:)         !< Subset Sl0 to subset Sl1
   integer,allocatable :: c0l0_to_s(:,:)     !< Grid Gf to subgrid
   integer,allocatable :: c2l1_to_s(:,:)     !< Grid Gs to subgrid
   integer,allocatable :: c2l1_to_c1(:,:)    !< Grid Gs to subset Sc1
   integer,allocatable :: c1l1_to_s(:,:)     !< Grid Gv to subgrid

   ! MPI distribution
   integer :: nc1a                              !< Number of points in subset Sc1 on halo A
   logical,allocatable :: lcheck_c1a(:)        !<
   logical,allocatable :: lcheck_c1b(:)        !<
   logical,allocatable :: lcheck_c1bb(:)        !<
   logical,allocatable :: lcheck_c2b(:,:)      !<
   logical,allocatable :: lcheck_sa(:)         !<
   logical,allocatable :: lcheck_sb(:)         !<
   logical,allocatable :: lcheck_sc(:)         !<
   logical,allocatable :: lcheck_sc_nor(:)         !<
   logical,allocatable :: lcheck_h(:,:)         !<
   logical,allocatable :: lcheck_s(:,:)         !<
   integer,allocatable :: c1a_to_c1(:)        !<
   integer,allocatable :: c1_to_c1a(:)        !<
   integer,allocatable :: c1b_to_c1(:)        !<
   integer,allocatable :: c1_to_c1b(:)        !<
   integer,allocatable :: c1bb_to_c1(:)       !<
   integer,allocatable :: c1_to_c1bb(:)       !<
   integer,allocatable :: c2l1_to_c2b(:,:)   !<
   integer,allocatable :: c2bl1_to_sb(:,:)   !<
   integer,allocatable :: interph_lg(:,:)   !<
   integer,allocatable :: interps_lg(:,:)   !<
   integer,allocatable :: sc_to_sb(:)         !< Subgrid, halo C to halo B
   integer,allocatable :: sa_to_s(:)          !<
   integer,allocatable :: s_to_sa(:)          !<
   integer,allocatable :: sb_to_s(:)          !<
   integer,allocatable :: s_to_sb(:)          !<
   integer,allocatable :: sc_to_s(:)          !<
   integer,allocatable :: s_to_sc(:)          !<

   ! Specific for normalization computation
   integer :: nsc_nor                           !< Number of subgrid nodes on halo C
   integer,allocatable :: sc_nor_to_s(:)          !<
   integer,allocatable :: s_to_sc_nor(:)          !<
   integer,allocatable :: sb_to_sc_nor(:)         !<
   type(linoptype) :: c_nor                     !< Convolution (with more coefficients to compute the normalization)

   ! Illustration
   integer,allocatable :: halo(:)               !< Halo points for illustration

   ! Required data to apply a localization

   ! Number of points
   integer :: nc1b                              !< Number of points in subset Sc1 on halo B
   integer :: nc1bb                             !< Number of points in subset Sc1 on halo B, extended for subgrid nodes on halo B
   integer :: nl1                               !< Number of levels in subset Sl1
   integer,allocatable :: nc2b(:)               !< Number of points in subset Sc2 on halo B
   integer :: nsa                               !< Number of subgrid nodes on halo A
   integer :: nsb                               !< Number of subgrid nodes on halo B
   integer :: nsc                               !< Number of subgrid nodes on halo C

   ! Inter-halo conversions
   integer,allocatable :: sa_to_sb(:)         !< Subgrid, halo A to halo B
   integer,allocatable :: sa_to_sc(:)         !< Subgrid, halo A to halo C
   integer,allocatable :: sb_to_sc(:)         !< Subgrid, halo B to halo C

   ! Linear operators
   type(linoptype) :: c                         !< Convolution
   type(linoptype),allocatable :: h(:)          !< Horizontal interpolation
   type(linoptype) :: v                         !< Vertical interpolation
   type(linoptype),allocatable :: s(:)          !< Subsample interpolation

   ! Copy conversions
   integer,allocatable :: sb_to_c2b(:)        !< Subgrid to subset Sc2 on halo B
   integer,allocatable :: sb_to_l1(:)         !< Subgrid to subset Sl1 on halo B

   ! Normalization
   real(kind_real),allocatable :: norm(:,:)     !< Normalization factor

   ! Localization weights
   real(kind_real),allocatable :: coef_ens(:,:) !< Ensemble coefficient square-root
   real(kind_real) :: wgt                       !< Main weight

   ! Communications
   type(comtype) :: AB                          !< Communication between halos A and B
   type(comtype) :: AC                          !< Communication between halos A and C
   integer :: mpicom                            !< Number of communication steps

   ! Transforms
   real(kind_real),allocatable :: trans(:,:)    !< Direct transform
   real(kind_real),allocatable :: transinv(:,:) !< Inverse transform
end type ndatatype

private
public :: ndatatype
public :: ndata_dealloc,ndata_read,ndata_write,ndata_write_mpi_summary

contains

!----------------------------------------------------------------------
! Subroutine: ndata_dealloc
!> Purpose: deallocate ndata object
!----------------------------------------------------------------------
subroutine ndata_dealloc(ndata,nicas_block)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data
logical,intent(in) :: nicas_block      !< NICAS block key

! Local variables
integer :: il0i,il1

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! TODO : partial deallocation

if (.false.) then
! Release memory
deallocate(ndata%coef_ens)
if (nicas_block) then
   deallocate(ndata%vbot)
   deallocate(ndata%vtop)
   do il0i=1,geom%nl0i
      call linop_dealloc(ndata%h(il0i))
   end do
   deallocate(ndata%h)
   call linop_dealloc(ndata%v)
   do il1=1,ndata%nl1
      call linop_dealloc(ndata%s(il1))
   end do
   deallocate(ndata%s)
   deallocate(ndata%norm)
   deallocate(ndata%s_to_c1)
   deallocate(ndata%s_to_l1)
   deallocate(ndata%s_to_c2)
   deallocate(ndata%c1_to_c0)
   deallocate(ndata%l1_to_l0)
   deallocate(ndata%c2l1_to_c0)
   deallocate(ndata%c0_to_c1)
   deallocate(ndata%l0_to_l1)
   deallocate(ndata%c0l0_to_s)
   deallocate(ndata%c2l1_to_s)
   deallocate(ndata%c2l1_to_c1)
   deallocate(ndata%c1l1_to_s)
   if (allocated(ndata%halo)) deallocate(ndata%halo)
end if
end if

! End associate
end associate

end subroutine ndata_dealloc

!----------------------------------------------------------------------
! Subroutine: ndata_read
!> Purpose: read ndata object
!----------------------------------------------------------------------
subroutine ndata_read(nam,geom,ndata,nicas_block,auto_block)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam       !< Namelist
type(geomtype),target,intent(inout) :: geom  !< Geometry
type(ndatatype),intent(inout) :: ndata       !< NICAS data
logical,intent(in) :: nicas_block            !< NICAS block key
logical,intent(in) :: auto_block             !< Autocovariance block key

! Local variables
integer :: ncid,info
integer :: nc0a_id,nc1b_id,nl1_id,nsa_id,nsb_id
integer :: nc2b_id,sb_to_c2b_id,sb_to_l1_id
integer :: sa_to_sb_id,sa_to_sc_id,sb_to_sc_id
integer :: norm_id,coef_ens_id,trans_id,transinv_id
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_read'

! Open file and get dimensions
filename = trim(nam%prefix)//'_'//trim(ndata%cname)//'.nc'
call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))

! Read main weight
call ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',ndata%wgt))

if (nicas_block) then
   ! Get dimensions
   call ncerr(subr,nf90_inq_dimid(ncid,'nc0a',nc0a_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,nc0a_id,len=geom%nc0a))
   info = nf90_inq_dimid(ncid,'nc1b',nc1b_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nc1b_id,len=ndata%nc1b))
   else
      ndata%nc1b = 0
   end if
   call ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=ndata%nl1))
   info = nf90_inq_dimid(ncid,'nsa',nsa_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nsa_id,len=ndata%nsa))
   else
      ndata%nsa = 0
   end if
   info = nf90_inq_dimid(ncid,'nsb',nsb_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nsb_id,len=ndata%nsb))
   else
      ndata%nsb = 0
   end if
   call ncerr(subr,nf90_get_att(ncid,nf90_global,'nsc',ndata%nsc))
   call ncerr(subr,nf90_get_att(ncid,nf90_global,'mpicom',ndata%mpicom))

   ! Allocation
   allocate(ndata%nc2b(ndata%nl1))
   if (ndata%nsb>0) allocate(ndata%sb_to_c2b(ndata%nsb))
   if (ndata%nsb>0) allocate(ndata%sb_to_l1(ndata%nsb))
   if (ndata%nsa>0) allocate(ndata%sa_to_sb(ndata%nsa))
   if (ndata%nsa>0) allocate(ndata%sa_to_sc(ndata%nsa))
   if (ndata%nsb>0) allocate(ndata%sb_to_sc(ndata%nsb))
   allocate(ndata%norm(geom%nc0a,geom%nl0))
   allocate(ndata%coef_ens(geom%nc0a,geom%nl0))

   ! Get variable id
   call ncerr(subr,nf90_inq_varid(ncid,'nc2b',nc2b_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_c2b',sb_to_c2b_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_l1',sb_to_l1_id))
   if (ndata%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'sa_to_sb',sa_to_sb_id))
   if (ndata%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'sa_to_sc',sa_to_sc_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_sc',sb_to_sc_id))
   call ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
   call ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))

   ! Read data
   call ncerr(subr,nf90_get_var(ncid,nc2b_id,ndata%nc2b))
   if (ndata%nsb>0) call ncerr(subr,nf90_get_var(ncid,sb_to_c2b_id,ndata%sb_to_c2b))
   if (ndata%nsb>0) call ncerr(subr,nf90_get_var(ncid,sb_to_l1_id,ndata%sb_to_l1))
   if (ndata%nsa>0) call ncerr(subr,nf90_get_var(ncid,sa_to_sb_id,ndata%sa_to_sb))
   if (ndata%nsa>0) call ncerr(subr,nf90_get_var(ncid,sa_to_sc_id,ndata%sa_to_sc))
   if (ndata%nsb>0) call ncerr(subr,nf90_get_var(ncid,sb_to_sc_id,ndata%sb_to_sc))
   call ncerr(subr,nf90_get_var(ncid,norm_id,ndata%norm))
   call ncerr(subr,nf90_get_var(ncid,coef_ens_id,ndata%coef_ens))

   ! Read communications
   call com_read(ncid,'AB',ndata%AB)
   call com_read(ncid,'AC',ndata%AC)

   ! Read linear operators
   call linop_read(ncid,'c',ndata%c)
   call linop_read(ncid,'h',ndata%h)
   call linop_read(ncid,'v',ndata%v)
   call linop_read(ncid,'s',ndata%s)
end if

if (nam%transform.and.auto_block) then
   ! Allocation
   allocate(ndata%trans(geom%nc0,geom%nl0))
   allocate(ndata%transinv(geom%nc0,geom%nl0))

   ! Get variable id
   call ncerr(subr,nf90_inq_varid(ncid,'trans',trans_id))
   call ncerr(subr,nf90_inq_varid(ncid,'transinv',transinv_id))

   ! Read data
   call ncerr(subr,nf90_get_var(ncid,trans_id,ndata%trans))
   call ncerr(subr,nf90_get_var(ncid,transinv_id,ndata%transinv))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine ndata_read

!----------------------------------------------------------------------
! Subroutine: ndata_write
!> Purpose: write ndata object
!----------------------------------------------------------------------
subroutine ndata_write(nam,geom,ndata,nicas_block,auto_block)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam    !< Namelist
type(geomtype),target,intent(in) :: geom  !< Geometry
type(ndatatype),intent(in) :: ndata       !< NICAS data
logical,intent(in) :: nicas_block         !< NICAS block key
logical,intent(in) :: auto_block          !< Autocovariance block key

! Local variables
integer :: ncid
integer :: nc0a_id,nl0_id,nc1b_id,nl1_id,nsa_id,nsb_id,nl0_2_id
integer :: nc2b_id,sb_to_c2b_id,sb_to_l1_id
integer :: sa_to_sb_id,sa_to_sc_id,sb_to_sc_id
integer :: norm_id,coef_ens_id,trans_id,transinv_id
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_write'

! Create file
filename = trim(nam%prefix)//'_'//trim(ndata%cname)//'.nc'
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call namncwrite(nam,ncid)

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
if (nicas_block) then
   call ncerr(subr,nf90_def_dim(ncid,'nc0a',geom%nc0a,nc0a_id))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',geom%nl0i))
   if (ndata%nc1b>0) call ncerr(subr,nf90_def_dim(ncid,'nc1b',ndata%nc1b,nc1b_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl1',ndata%nl1,nl1_id))
   if (ndata%nsa>0) call ncerr(subr,nf90_def_dim(ncid,'nsa',ndata%nsa,nsa_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_def_dim(ncid,'nsb',ndata%nsb,nsb_id))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'nsc',ndata%nsc))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'mpicom',ndata%mpicom))
end if
if (nam%transform.and.auto_block) call ncerr(subr,nf90_def_dim(ncid,'nl0_2',geom%nl0,nl0_2_id))

! Write main weight
call ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',ndata%wgt))

! Define variables
if (nicas_block) then
   call ncerr(subr,nf90_def_var(ncid,'nc2b',nf90_int,(/nl1_id/),nc2b_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_def_var(ncid,'sb_to_c2b',nf90_int,(/nsb_id/),sb_to_c2b_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_def_var(ncid,'sb_to_l1',nf90_int,(/nsb_id/),sb_to_l1_id))
   if (ndata%nsa>0) call ncerr(subr,nf90_def_var(ncid,'sa_to_sb',nf90_int,(/nsa_id/),sa_to_sb_id))
   if (ndata%nsa>0) call ncerr(subr,nf90_def_var(ncid,'sa_to_sc',nf90_int,(/nsa_id/),sa_to_sc_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_def_var(ncid,'sb_to_sc',nf90_int,(/nsb_id/),sb_to_sc_id))
   call ncerr(subr,nf90_def_var(ncid,'norm',ncfloat,(/nc0a_id,nl0_id/),norm_id))
   call ncerr(subr,nf90_def_var(ncid,'coef_ens',ncfloat,(/nc0a_id,nl0_id/),coef_ens_id))

   if (ndata%nsa>0) call ncerr(subr,nf90_put_att(ncid,sa_to_sb_id,'_FillValue',msvali))
   if (ndata%nsa>0) call ncerr(subr,nf90_put_att(ncid,sa_to_sc_id,'_FillValue',msvali))
   if (ndata%nsb>0) call ncerr(subr,nf90_put_att(ncid,sb_to_sc_id,'_FillValue',msvali))
   call ncerr(subr,nf90_put_att(ncid,norm_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_put_att(ncid,coef_ens_id,'_FillValue',msvalr))
end if
if (nam%transform.and.auto_block) then
   call ncerr(subr,nf90_def_var(ncid,'trans',ncfloat,(/nl0_id,nl0_2_id/),trans_id))
   call ncerr(subr,nf90_put_att(ncid,trans_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'transinv',ncfloat,(/nl0_id,nl0_2_id/),transinv_id))
   call ncerr(subr,nf90_put_att(ncid,transinv_id,'_FillValue',msvalr))
end if

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write variables
if (nicas_block) then
   call ncerr(subr,nf90_put_var(ncid,nc2b_id,ndata%nc2b))
   if (ndata%nsb>0) call ncerr(subr,nf90_put_var(ncid,sb_to_c2b_id,ndata%sb_to_c2b))
   if (ndata%nsb>0) call ncerr(subr,nf90_put_var(ncid,sb_to_l1_id,ndata%sb_to_l1))
   if (ndata%nsa>0) call ncerr(subr,nf90_put_var(ncid,sa_to_sb_id,ndata%sa_to_sb))
   if (ndata%nsa>0) call ncerr(subr,nf90_put_var(ncid,sa_to_sc_id,ndata%sa_to_sc))
   if (ndata%nsb>0) call ncerr(subr,nf90_put_var(ncid,sb_to_sc_id,ndata%sb_to_sc))
   call ncerr(subr,nf90_put_var(ncid,norm_id,ndata%norm))
   call ncerr(subr,nf90_put_var(ncid,coef_ens_id,ndata%coef_ens))

   ! Write communications
   call com_write(ncid,ndata%AB)
   call com_write(ncid,ndata%AC)

   ! Write linear operators
   call linop_write(ncid,ndata%c)
   call linop_write(ncid,ndata%h)
   call linop_write(ncid,ndata%v)
   call linop_write(ncid,ndata%s)
end if
if (nam%transform.and.auto_block) then
   call ncerr(subr,nf90_put_var(ncid,trans_id,ndata%trans))
   call ncerr(subr,nf90_put_var(ncid,transinv_id,ndata%transinv))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine ndata_write

!----------------------------------------------------------------------
! Subroutine: ndata_write_mpi_summary
!> Purpose: write ndata object
!----------------------------------------------------------------------
subroutine ndata_write_mpi_summary(ndata)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata !< NICAS data

! Local variables
integer :: ncid
integer :: nc0_id
integer :: lon_id,lat_id,c0_to_proc_id,halo_id
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_write_mpi_summary'

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Create summary file
filename = trim(nam%prefix)//'_'//trim(ndata%cname)//'_summary.nc'
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call namncwrite(nam,ncid)

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))

! Define variables
call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
call ncerr(subr,nf90_def_var(ncid,'c0_to_proc',nf90_int,(/nc0_id/),c0_to_proc_id))
call ncerr(subr,nf90_def_var(ncid,'halo',nf90_int,(/nc0_id/),halo_id))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write variables
call ncerr(subr,nf90_put_var(ncid,lon_id,ndata%geom%lon*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_id,ndata%geom%lat*rad2deg))
call ncerr(subr,nf90_put_var(ncid,c0_to_proc_id,ndata%geom%c0_to_proc))
call ncerr(subr,nf90_put_var(ncid,halo_id,ndata%halo))

! Close summary file
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine ndata_write_mpi_summary

end module type_ndata
