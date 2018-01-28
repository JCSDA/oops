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
   integer,allocatable :: s_to_c1(:)            !< Subgrid to subset Sc1
   integer,allocatable :: s_to_l1(:)            !< Subgrid to subset Sl1
   integer,allocatable :: c1_to_c0(:)           !< Subset Sc1 to subset Sc0
   integer,allocatable :: l1_to_l0(:)           !< Subset Sl1 to subset Sl0
   integer,allocatable :: c0_to_c1(:)           !< Subset Sc0 to subset Sc1
   integer,allocatable :: l0_to_l1(:)           !< Subset Sl0 to subset Sl1
   logical,allocatable :: c2mask(:,:)           !< Mask from G1 to subgrid
   integer,allocatable :: c1l1_to_s(:,:)        !< Grid Gv to subgrid

   ! MPI distribution
   integer :: nc1a                              !< Number of points in subset Sc1 on halo A
   logical,allocatable :: lcheck_c1a(:)         !< Detection of halo A on subset Sc1
   logical,allocatable :: lcheck_c1b(:)         !< Detection of halo B on subset Sc1
   logical,allocatable :: lcheck_sa(:)          !< Detection of halo A on subgrid
   logical,allocatable :: lcheck_sb(:)          !< Detection of halo B on subgrid
   logical,allocatable :: lcheck_sc(:)          !< Detection of halo C on subgrid
   logical,allocatable :: lcheck_sc_nor(:)      !< Detection of halo A on subgrid
   logical,allocatable :: lcheck_h(:,:)         !< Detection of horizontal interpolation coefficients
   logical,allocatable :: lcheck_s(:,:)         !< Detection of subsampling interpolation coefficients
   integer,allocatable :: c1a_to_c1(:)          !< Subset Sc1, halo A to global
   integer,allocatable :: c1_to_c1a(:)          !< Subset Sc1, global to halo A
   integer,allocatable :: c1b_to_c1(:)          !< Subset Sc1, halo B to global
   integer,allocatable :: c1_to_c1b(:)          !< Subset Sc1, global to halo B
   integer,allocatable :: c1bl1_to_sb(:,:)      !< Halo B, subset Sc1 to subgrid
   integer,allocatable :: interph_lg(:,:)       !< Local to global for horizontal interpolation
   integer,allocatable :: interps_lg(:,:)       !< Local to global for subsampling interpolation
   integer,allocatable :: sc_to_sb(:)           !< Subgrid, halo C to halo B
   integer,allocatable :: sa_to_s(:)            !< Subgrid, halo A to global
   integer,allocatable :: s_to_sa(:)            !< Subgrid, global to halo A
   integer,allocatable :: sb_to_s(:)            !< Subgrid, halo B to global
   integer,allocatable :: s_to_sb(:)            !< Subgrid, global to halo B
   integer,allocatable :: sc_to_s(:)            !< Subgrid, halo C to global
   integer,allocatable :: s_to_sc(:)            !< Subgrid, global to halo C

   ! Extended data for normalization computation
   integer :: nsc_nor                           !< Number of subgrid nodes on halo C (extended for normalization)
   integer,allocatable :: sc_nor_to_s(:)        !< Subgrid, halo C to global (extended for normalization)
   integer,allocatable :: s_to_sc_nor(:)        !< Subgrid, global to halo C (extended for normalization)
   integer,allocatable :: sb_to_sc_nor(:)       !< Subgrid, halo B to halo C (extended for normalization)
   type(linoptype) :: c_nor                     !< Convolution (extended for normalization)

   ! Required data to apply a localization

   ! Number of points
   integer :: nc1b                              !< Number of points in subset Sc1 on halo B
   integer :: nl1                               !< Number of levels in subset Sl1
   integer :: nsa                               !< Number of subgrid nodes on halo A
   integer :: nsb                               !< Number of subgrid nodes on halo B
   integer :: nsc                               !< Number of subgrid nodes on halo C

   ! Inter-halo conversions
   integer,allocatable :: sa_to_sb(:)           !< Subgrid, halo A to halo B
   integer,allocatable :: sa_to_sc(:)           !< Subgrid, halo A to halo C
   integer,allocatable :: sb_to_sc(:)           !< Subgrid, halo B to halo C

   ! Linear operators
   type(linoptype) :: c                         !< Convolution
   type(linoptype),allocatable :: h(:)          !< Horizontal interpolation
   type(linoptype) :: v                         !< Vertical interpolation
   type(linoptype),allocatable :: s(:)          !< Subsample interpolation

   ! Copy conversions
   integer,allocatable :: sb_to_c1b(:)          !< Subgrid to subset Sc1 on halo B
   integer,allocatable :: sb_to_l1(:)           !< Subgrid to subset Sl1 on halo B

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
!> Purpose: ndata object deallocation
!----------------------------------------------------------------------
subroutine ndata_dealloc(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: il0,il1

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Release memory
if (allocated(ndata%vbot)) deallocate(ndata%vbot)
if (allocated(ndata%vtop)) deallocate(ndata%vtop)
if (allocated(ndata%nc2)) deallocate(ndata%nc2)
if (allocated(ndata%hfull)) then
   do il0=1,geom%nl0
      call linop_dealloc(ndata%hfull(il0))
   end do
   deallocate(ndata%hfull)
end if
if (allocated(ndata%sfull)) then
   do il1=1,ndata%nl1
      call linop_dealloc(ndata%sfull(il1))
   end do
   deallocate(ndata%sfull)
end if
if (allocated(ndata%llev)) deallocate(ndata%llev)
if (allocated(ndata%s_to_c1)) deallocate(ndata%s_to_c1)
if (allocated(ndata%s_to_l1)) deallocate(ndata%s_to_l1)
if (allocated(ndata%c1_to_c0)) deallocate(ndata%c1_to_c0)
if (allocated(ndata%l1_to_l0)) deallocate(ndata%l1_to_l0)
if (allocated(ndata%c0_to_c1)) deallocate(ndata%c0_to_c1)
if (allocated(ndata%l0_to_l1)) deallocate(ndata%l0_to_l1)
if (allocated(ndata%c2mask)) deallocate(ndata%c2mask)
if (allocated(ndata%c1l1_to_s)) deallocate(ndata%c1l1_to_s)
if (allocated(ndata%lcheck_c1a)) deallocate(ndata%lcheck_c1a)
if (allocated(ndata%lcheck_c1b)) deallocate(ndata%lcheck_c1b)
if (allocated(ndata%lcheck_sa)) deallocate(ndata%lcheck_sa)
if (allocated(ndata%lcheck_sb)) deallocate(ndata%lcheck_sb)
if (allocated(ndata%lcheck_sc)) deallocate(ndata%lcheck_sc)
if (allocated(ndata%lcheck_sc_nor)) deallocate(ndata%lcheck_sc_nor)
if (allocated(ndata%lcheck_h)) deallocate(ndata%lcheck_h)
if (allocated(ndata%lcheck_s)) deallocate(ndata%lcheck_s)
if (allocated(ndata%c1a_to_c1)) deallocate(ndata%c1a_to_c1)
if (allocated(ndata%c1_to_c1a)) deallocate(ndata%c1_to_c1a)
if (allocated(ndata%c1b_to_c1)) deallocate(ndata%c1b_to_c1)
if (allocated(ndata%c1_to_c1b)) deallocate(ndata%c1_to_c1b)
if (allocated(ndata%c1bl1_to_sb)) deallocate(ndata%c1bl1_to_sb)
if (allocated(ndata%interph_lg)) deallocate(ndata%interph_lg)
if (allocated(ndata%interps_lg)) deallocate(ndata%interps_lg)
if (allocated(ndata%sc_to_sb)) deallocate(ndata%sc_to_sb)
if (allocated(ndata%sa_to_s)) deallocate(ndata%sa_to_s)
if (allocated(ndata%s_to_sa)) deallocate(ndata%s_to_sa)
if (allocated(ndata%sb_to_s)) deallocate(ndata%sb_to_s)
if (allocated(ndata%s_to_sb)) deallocate(ndata%s_to_sb)
if (allocated(ndata%sc_to_s)) deallocate(ndata%sc_to_s)
if (allocated(ndata%s_to_sc)) deallocate(ndata%s_to_sc)
if (allocated(ndata%sc_nor_to_s)) deallocate(ndata%sc_nor_to_s)
if (allocated(ndata%s_to_sc_nor)) deallocate(ndata%s_to_sc_nor)
if (allocated(ndata%sb_to_sc_nor)) deallocate(ndata%sb_to_sc_nor)
call linop_dealloc(ndata%c_nor)
if (allocated(ndata%sa_to_sb)) deallocate(ndata%sa_to_sb)
if (allocated(ndata%sa_to_sc)) deallocate(ndata%sa_to_sc)
if (allocated(ndata%sb_to_sc)) deallocate(ndata%sb_to_sc)
call linop_dealloc(ndata%c)
if (allocated(ndata%h)) then
   do il0=1,geom%nl0
      call linop_dealloc(ndata%h(il0))
   end do
   deallocate(ndata%h)
end if
call linop_dealloc(ndata%v)
if (allocated(ndata%s)) then
   do il1=1,ndata%nl1
      call linop_dealloc(ndata%s(il1))
   end do
   deallocate(ndata%s)
end if
if (allocated(ndata%sb_to_c1b)) deallocate(ndata%sb_to_c1b)
if (allocated(ndata%sb_to_l1)) deallocate(ndata%sb_to_l1)
if (allocated(ndata%norm)) deallocate(ndata%norm)
if (allocated(ndata%coef_ens)) deallocate(ndata%coef_ens)
if (allocated(ndata%trans)) deallocate(ndata%trans)
call com_dealloc(ndata%AB)
call com_dealloc(ndata%AC)
if (allocated(ndata%transinv)) deallocate(ndata%transinv)
if (allocated(ndata%sb_to_c1b)) deallocate(ndata%sb_to_c1b)

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
integer :: sb_to_c1b_id,sb_to_l1_id
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
   if (ndata%nsb>0) allocate(ndata%sb_to_c1b(ndata%nsb))
   if (ndata%nsb>0) allocate(ndata%sb_to_l1(ndata%nsb))
   if (ndata%nsa>0) allocate(ndata%sa_to_sb(ndata%nsa))
   if (ndata%nsa>0) allocate(ndata%sa_to_sc(ndata%nsa))
   if (ndata%nsb>0) allocate(ndata%sb_to_sc(ndata%nsb))
   allocate(ndata%norm(geom%nc0a,geom%nl0))
   allocate(ndata%coef_ens(geom%nc0a,geom%nl0))

   ! Get variable id
   if (ndata%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_c1b',sb_to_c1b_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_l1',sb_to_l1_id))
   if (ndata%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'sa_to_sb',sa_to_sb_id))
   if (ndata%nsa>0) call ncerr(subr,nf90_inq_varid(ncid,'sa_to_sc',sa_to_sc_id))
   if (ndata%nsb>0) call ncerr(subr,nf90_inq_varid(ncid,'sb_to_sc',sb_to_sc_id))
   call ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
   call ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))

   ! Read data
   if (ndata%nsb>0) call ncerr(subr,nf90_get_var(ncid,sb_to_c1b_id,ndata%sb_to_c1b))
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
integer :: sb_to_c1b_id,sb_to_l1_id
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
   if (ndata%nsb>0) call ncerr(subr,nf90_def_var(ncid,'sb_to_c1b',nf90_int,(/nsb_id/),sb_to_c1b_id))
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
   if (ndata%nsb>0) call ncerr(subr,nf90_put_var(ncid,sb_to_c1b_id,ndata%sb_to_c1b))
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
integer :: ncid,nc0_id,nc1_id,nl1_id
integer :: lon_id,lat_id,c0_to_proc_id,c1_to_c0_id,l1_to_l0_id,lcheck_id
integer :: is,ic1,il1
real(kind_real) :: lcheck(ndata%nc1,ndata%nl1)
character(len=1024) :: filename
character(len=1024) :: subr = 'ndata_write_mpi_summary'

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Create summary file
filename = trim(nam%prefix)//'_'//trim(ndata%cname)//'_summary.nc'
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call namncwrite(nam,ncid)

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
call ncerr(subr,nf90_def_dim(ncid,'nc1',ndata%nc1,nc1_id))
call ncerr(subr,nf90_def_dim(ncid,'nl1',ndata%nl1,nl1_id))

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
call ncerr(subr,nf90_put_var(ncid,lon_id,ndata%geom%lon*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_id,ndata%geom%lat*rad2deg))
call ncerr(subr,nf90_put_var(ncid,c0_to_proc_id,ndata%geom%c0_to_proc))
call ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,ndata%c1_to_c0))
call ncerr(subr,nf90_put_var(ncid,l1_to_l0_id,ndata%l1_to_l0))
call msr(lcheck)
do is=1,ndata%ns
   ic1 = ndata%s_to_c1(is)
   il1 = ndata%s_to_l1(is)
   if (ndata%lcheck_sa(is)) then
      lcheck(ic1,il1) = 1.0
   elseif (ndata%lcheck_sb(is)) then
      lcheck(ic1,il1) = 2.0
   elseif (ndata%lcheck_sc(is)) then
      lcheck(ic1,il1) = 3.0
   else
      lcheck(ic1,il1) = 4.0
   end if
end do
call ncerr(subr,nf90_put_var(ncid,lcheck_id,lcheck))

! Close summary file
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine ndata_write_mpi_summary

end module type_ndata
