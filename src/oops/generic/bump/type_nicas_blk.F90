!----------------------------------------------------------------------
! Module: type_nicas_blk
!> Purpose: NICAS data block derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_nicas_blk

use netcdf
!$ use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg,msvalr
use tools_func, only: sphere_dist,vector_product,vector_triple_product,gc99
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat
use tools_qsort, only: qsort
use tools_repro, only: supeq
use tools_test, only: define_dirac
use type_bpar, only: bpar_type
use type_cmat_blk, only: cmat_blk_type
use type_com, only: com_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_kdtree, only: kdtree_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

integer,parameter :: nc1max = 15000                       !< Maximum size of the Sc1 subset
logical,parameter :: write_grids = .false.                !< Write Sc1 and Sc2 subsets lon/lat
real(kind_real),parameter :: sqrt_r = 0.721_kind_real     !< Square-root factor (empirical)
real(kind_real),parameter :: sqrt_r_dble = 0.86_kind_real !< Square-root factor (empirical)
real(kind_real),parameter :: sqrt_rfac = 0.9_kind_real    !< Square-root factor (empirical)
real(kind_real),parameter :: sqrt_coef = 0.54_kind_real   !< Square-root factor (empirical)
real(kind_real),parameter :: S_inf = 1.0e-2_kind_real     !< Minimum value for the convolution coefficients
real(kind_real),parameter :: tol = 1.0e-3_kind_real       !< Positive-definiteness test tolerance
integer,parameter :: nitermax = 50                        !< Number of iterations for the positive-definiteness test
logical,parameter :: test_no_point = .false.              !< Test NICAS with no subgrid point on the last MPI task

! Ball data derived type
type balldata_type
   integer :: nbd                        !< Number of values
   integer,allocatable :: bd_to_c1(:)    !< Ball data index to subset Sc1
   integer,allocatable :: bd_to_l1(:)    !< Ball data index to subset Sl1
   real(kind_real),allocatable :: val(:) !< Values
contains
   procedure :: alloc => balldata_alloc
   procedure :: dealloc => balldata_dealloc
   procedure :: pack => balldata_pack
end type balldata_type

! NICAS block derived type
type nicas_blk_type
   ! Block index and name
   integer :: ib                                !< Block index
   character(len=1024) :: name                  !< Name
   logical :: double_fit                        !< Double fit

   ! Specific geometry
   integer :: nc1                               !< Number of points in subset Sc1
   integer,allocatable :: vbot(:)               !< Bottom level in grid Gh
   integer,allocatable :: vtop(:)               !< Top level in grid Gh
   integer,allocatable :: nc2(:)                !< Number of points in subset Sc2
   integer :: ns                                !< Number of subgrid nodes

   ! Full interpolations
   type(linop_type),allocatable :: hfull(:)     !< Horizontal interpolation
   type(linop_type) :: vfull                    !< Vertical interpolation
   type(linop_type),allocatable :: sfull(:)     !< Subsample interpolation

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
   logical,allocatable :: mask_c1(:,:)          !< Mask from subset C1 to subgrid
   logical,allocatable :: mask_c2(:,:)          !< Mask from subset C2 to subgrid
   integer,allocatable :: c1l1_to_s(:,:)        !< Grid Gv to subgrid
   integer,allocatable :: c1_to_proc(:)         !< Subset Sc1 to local task
   integer,allocatable :: s_to_proc(:)          !< Subgrid to local task


   ! MPI distribution
   integer :: nc1a                              !< Number of points in subset Sc1 on halo A
   integer :: nc1bb                             !< Number of points in subset Sc1 on halo B (extended)
   integer :: nsbb                              !< Number of points in subgrid on halo B (extended)
   logical,allocatable :: lcheck_sa(:)          !< Detection of halo A on subgrid
   logical,allocatable :: lcheck_sb(:)          !< Detection of halo B on subgrid
   logical,allocatable :: lcheck_sc(:)          !< Detection of halo C on subgrid
   integer,allocatable :: c1a_to_c1(:)          !< Subset Sc1, halo A to global
   integer,allocatable :: c1_to_c1a(:)          !< Subset Sc1, global to halo A
   integer,allocatable :: c1b_to_c1(:)          !< Subset Sc1, halo B to global
   integer,allocatable :: c1_to_c1b(:)          !< Subset Sc1, global to halo B
   integer,allocatable :: c1bl1_to_sb(:,:)      !< Halo B, subset Sc1 to subgrid
   integer,allocatable :: interph_lg(:,:)       !< Local to global for horizontal interpolation
   integer,allocatable :: interps_lg(:,:)       !< Local to global for subsampling interpolation
   integer,allocatable :: c1a_to_c0a(:)         !< Halo A, subset Sc1 to subset Sc0
   integer,allocatable :: sa_to_c0a(:)          !< Halo A, subgrid to subset Sc0
   integer,allocatable :: sa_to_l0(:)           !< Halo A, subgrid to subset Sl0
   integer,allocatable :: sa_to_sb(:)           !< Subgrid, halo A to halo B
   integer,allocatable :: sc_to_sb(:)           !< Subgrid, halo C to halo B
   integer,allocatable :: sa_to_s(:)            !< Subgrid, halo A to global
   integer,allocatable :: s_to_sa(:)            !< Subgrid, global to halo A
   integer,allocatable :: sb_to_s(:)            !< Subgrid, halo B to global
   integer,allocatable :: s_to_sb(:)            !< Subgrid, global to halo B
   integer,allocatable :: sc_to_s(:)            !< Subgrid, halo C to global
   integer,allocatable :: s_to_sc(:)            !< Subgrid, global to halo C
   integer,allocatable :: sbb_to_s(:)           !< Subgrid, halo B (extended) to global
   integer,allocatable :: c1_to_c1bb(:)         !< Subset Sc1, global to halo B (extended)
   integer,allocatable :: c1bb_to_c1(:)         !< Subset Sc1, halo B (extended) to global

   ! Convolution parameters
   real(kind_real) :: rhmax                        !<
   real(kind_real),allocatable :: rh_c1(:,:)       !<
   real(kind_real),allocatable :: rv_c1(:,:)       !<
   type(balldata_type),allocatable :: distnorm(:)  !<
   type(balldata_type),allocatable :: distnormv(:) !<
   type(balldata_type),allocatable :: rfac(:)      !<
   type(balldata_type),allocatable :: coef(:)      !<

   ! Extended data for normalization computation
   integer :: nsc_nor                           !< Number of subgrid nodes on halo C (extended for normalization)
   integer,allocatable :: sc_nor_to_s(:)        !< Subgrid, halo C to global (extended for normalization)
   integer,allocatable :: s_to_sc_nor(:)        !< Subgrid, global to halo C (extended for normalization)
   integer,allocatable :: sb_to_sc_nor(:)       !< Subgrid, halo B to halo C (extended for normalization)
   type(linop_type) :: c_nor                    !< Convolution (extended for normalization)

   ! KD-tree
   type(kdtree_type) :: kdtree                  !< KD-tree

   ! Required data to apply a localization

   ! Number of points
   integer :: nc1b                              !< Number of points in subset Sc1 on halo B
   integer :: nl1                               !< Number of levels in subset Sl1
   integer :: nsa                               !< Number of subgrid nodes on halo A
   integer :: nsb                               !< Number of subgrid nodes on halo B
   integer :: nsc                               !< Number of subgrid nodes on halo C
   integer :: nc0d                              !< Number of points in subset Sc1 on halo D
   integer :: nc0dinv                           !< Number of points in subset Sc1 on halo Dinv

   ! Inter-halo conversions
   integer,allocatable :: sa_to_sc(:)           !< Subgrid, halo A to halo C
   integer,allocatable :: sb_to_sc(:)           !< Subgrid, halo B to halo C

   ! Linear operators
   type(linop_type) :: c                        !< Convolution
   type(linop_type),allocatable :: h(:)         !< Horizontal interpolation
   type(linop_type) :: v                        !< Vertical interpolation
   type(linop_type),allocatable :: s(:)         !< Subsample interpolation
   type(linop_type),allocatable :: d(:,:)       !< Advection
   type(linop_type),allocatable :: dinv(:,:)    !< Inverse advection

   ! Copy conversions
   integer,allocatable :: sb_to_c1b(:)          !< Subgrid to subset Sc1 on halo B
   integer,allocatable :: sb_to_l1(:)           !< Subgrid to subset Sl1 on halo B

   ! Normalization
   real(kind_real),allocatable :: norm(:,:)     !< Normalization factor

   ! Localization weights
   real(kind_real),allocatable :: coef_ens(:,:) !< Ensemble coefficient square-root
   real(kind_real) :: wgt                       !< Main weight

   ! Communications
   type(com_type) :: com_AB                     !< Communication between halos A and B
   type(com_type) :: com_AC                     !< Communication between halos A and C
   type(com_type) :: com_AD                     !< Communication between halos A and D
   type(com_type) :: com_ADinv                  !< Communication between halos A and Dinv
contains
   procedure :: dealloc => nicas_blk_dealloc
   procedure :: compute_parameters => nicas_blk_compute_parameters
   procedure :: compute_sampling => nicas_blk_compute_sampling
   procedure :: compute_interp_h => nicas_blk_compute_interp_h
   procedure :: compute_interp_v => nicas_blk_compute_interp_v
   procedure :: compute_interp_s => nicas_blk_compute_interp_s
   procedure :: compute_mpi_ab => nicas_blk_compute_mpi_ab
   procedure :: compute_convol => nicas_blk_compute_convol
   procedure :: compute_convol_network => nicas_blk_compute_convol_network
   procedure :: compute_convol_distance => nicas_blk_compute_convol_distance
   procedure :: compute_convol_weights => nicas_blk_compute_convol_weights
   procedure :: compute_mpi_c => nicas_blk_compute_mpi_c
   procedure :: compute_normalization => nicas_blk_compute_normalization
   procedure :: compute_adv => nicas_blk_compute_adv
   procedure :: apply => nicas_blk_apply
   procedure :: apply_from_sqrt => nicas_blk_apply_from_sqrt
   procedure :: apply_sqrt => nicas_blk_apply_sqrt
   procedure :: apply_sqrt_ad => nicas_blk_apply_sqrt_ad
   procedure :: apply_interp => nicas_blk_apply_interp
   procedure :: apply_interp_ad => nicas_blk_apply_interp_ad
   procedure :: apply_interp_h => nicas_blk_apply_interp_h
   procedure :: apply_interp_h_ad => nicas_blk_apply_interp_h_ad
   procedure :: apply_interp_v => nicas_blk_apply_interp_v
   procedure :: apply_interp_v_ad => nicas_blk_apply_interp_v_ad
   procedure :: apply_interp_s => nicas_blk_apply_interp_s
   procedure :: apply_interp_s_ad => nicas_blk_apply_interp_s_ad
   procedure :: apply_convol => nicas_blk_apply_convol
   procedure :: apply_adv => nicas_blk_apply_adv
   procedure :: apply_adv_ad => nicas_blk_apply_adv_ad
   procedure :: apply_adv_inv => nicas_blk_apply_adv_inv
   procedure :: test_adjoint => nicas_blk_test_adjoint
   procedure :: test_pos_def => nicas_blk_test_pos_def
   procedure :: test_sqrt => nicas_blk_test_sqrt
   procedure :: test_dirac => nicas_blk_test_dirac
end type nicas_blk_type

private
public :: nicas_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: balldata_alloc
!> Purpose: ball data allocation
!----------------------------------------------------------------------
subroutine balldata_alloc(balldata)

implicit none

! Passed variables
class(balldata_type),intent(inout) :: balldata !< Ball data

! Allocation
allocate(balldata%bd_to_c1(balldata%nbd))
allocate(balldata%bd_to_l1(balldata%nbd))
allocate(balldata%val(balldata%nbd))

end subroutine balldata_alloc

!----------------------------------------------------------------------
! Subroutine: balldata_dealloc
!> Purpose: ball data deallocation
!----------------------------------------------------------------------
subroutine balldata_dealloc(balldata)

implicit none

! Passed variables
class(balldata_type),intent(inout) :: balldata !< Ball data

! Release memory
if (allocated(balldata%bd_to_c1)) deallocate(balldata%bd_to_c1)
if (allocated(balldata%bd_to_l1)) deallocate(balldata%bd_to_l1)
if (allocated(balldata%val)) deallocate(balldata%val)

end subroutine balldata_dealloc

!----------------------------------------------------------------------
! Subroutine: balldata_pack
!> Purpose: pack data into balldata object
!----------------------------------------------------------------------
subroutine balldata_pack(balldata,nc1,nl1,val)

implicit none

! Passed variables
class(balldata_type),intent(inout) :: balldata !< Ball data
integer,intent(in) :: nc1                      !< TODO
integer,intent(in) :: nl1                      !< TODO
real(kind_real),intent(in) :: val(nc1,nl1)     !< TODO

! Local variables
integer :: ibd,ic1,il1

! Count non-missing values
balldata%nbd = count(isnotmsr(val))

! Allocation
call balldata%alloc

! Pack data
ibd = 0
do il1=1,nl1
   do ic1=1,nc1
      if (isnotmsr(val(ic1,il1))) then
         ibd = ibd+1
         balldata%bd_to_c1(ibd) = ic1
         balldata%bd_to_l1(ibd) = il1
         balldata%val(ibd) = val(ic1,il1)
      end if
   end do
end do

end subroutine balldata_pack

!----------------------------------------------------------------------
! Subroutine: nicas_blk_dealloc
!> Purpose: NICAS block data deallocation
!----------------------------------------------------------------------
subroutine nicas_blk_dealloc(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),target,intent(in) :: nam          !< Namelist
type(geom_type),target,intent(in) :: geom        !< Geometry

! Local variables
integer :: il0,il1,its,isbb

! Release memory
if (allocated(nicas_blk%vbot)) deallocate(nicas_blk%vbot)
if (allocated(nicas_blk%vtop)) deallocate(nicas_blk%vtop)
if (allocated(nicas_blk%nc2)) deallocate(nicas_blk%nc2)
if (allocated(nicas_blk%hfull)) then
   do il0=1,geom%nl0i
      call nicas_blk%hfull(il0)%dealloc
   end do
   deallocate(nicas_blk%hfull)
end if
call nicas_blk%vfull%dealloc
if (allocated(nicas_blk%sfull)) then
   do il1=1,nicas_blk%nl1
      call nicas_blk%sfull(il1)%dealloc
   end do
   deallocate(nicas_blk%sfull)
end if
if (allocated(nicas_blk%llev)) deallocate(nicas_blk%llev)
if (allocated(nicas_blk%s_to_c1)) deallocate(nicas_blk%s_to_c1)
if (allocated(nicas_blk%s_to_l1)) deallocate(nicas_blk%s_to_l1)
if (allocated(nicas_blk%c1_to_c0)) deallocate(nicas_blk%c1_to_c0)
if (allocated(nicas_blk%l1_to_l0)) deallocate(nicas_blk%l1_to_l0)
if (allocated(nicas_blk%c0_to_c1)) deallocate(nicas_blk%c0_to_c1)
if (allocated(nicas_blk%l0_to_l1)) deallocate(nicas_blk%l0_to_l1)
if (allocated(nicas_blk%mask_c1)) deallocate(nicas_blk%mask_c1)
if (allocated(nicas_blk%mask_c2)) deallocate(nicas_blk%mask_c2)
if (allocated(nicas_blk%c1l1_to_s)) deallocate(nicas_blk%c1l1_to_s)
if (allocated(nicas_blk%c1_to_proc)) deallocate(nicas_blk%c1_to_proc)
if (allocated(nicas_blk%s_to_proc)) deallocate(nicas_blk%s_to_proc)
if (allocated(nicas_blk%lcheck_sa)) deallocate(nicas_blk%lcheck_sa)
if (allocated(nicas_blk%lcheck_sb)) deallocate(nicas_blk%lcheck_sb)
if (allocated(nicas_blk%lcheck_sc)) deallocate(nicas_blk%lcheck_sc)
if (allocated(nicas_blk%c1a_to_c1)) deallocate(nicas_blk%c1a_to_c1)
if (allocated(nicas_blk%c1_to_c1a)) deallocate(nicas_blk%c1_to_c1a)
if (allocated(nicas_blk%c1b_to_c1)) deallocate(nicas_blk%c1b_to_c1)
if (allocated(nicas_blk%c1_to_c1b)) deallocate(nicas_blk%c1_to_c1b)
if (allocated(nicas_blk%c1bl1_to_sb)) deallocate(nicas_blk%c1bl1_to_sb)
if (allocated(nicas_blk%interph_lg)) deallocate(nicas_blk%interph_lg)
if (allocated(nicas_blk%interps_lg)) deallocate(nicas_blk%interps_lg)
if (allocated(nicas_blk%c1a_to_c0a)) deallocate(nicas_blk%c1a_to_c0a)
if (allocated(nicas_blk%sa_to_c0a)) deallocate(nicas_blk%sa_to_c0a)
if (allocated(nicas_blk%sa_to_l0)) deallocate(nicas_blk%sa_to_l0)
if (allocated(nicas_blk%sa_to_sb)) deallocate(nicas_blk%sa_to_sb)
if (allocated(nicas_blk%sc_to_sb)) deallocate(nicas_blk%sc_to_sb)
if (allocated(nicas_blk%sa_to_s)) deallocate(nicas_blk%sa_to_s)
if (allocated(nicas_blk%s_to_sa)) deallocate(nicas_blk%s_to_sa)
if (allocated(nicas_blk%sb_to_s)) deallocate(nicas_blk%sb_to_s)
if (allocated(nicas_blk%s_to_sb)) deallocate(nicas_blk%s_to_sb)
if (allocated(nicas_blk%sc_to_s)) deallocate(nicas_blk%sc_to_s)
if (allocated(nicas_blk%s_to_sc)) deallocate(nicas_blk%s_to_sc)
if (allocated(nicas_blk%sbb_to_s)) deallocate(nicas_blk%sbb_to_s)
if (allocated(nicas_blk%c1_to_c1bb)) deallocate(nicas_blk%c1_to_c1bb)
if (allocated(nicas_blk%c1bb_to_c1)) deallocate(nicas_blk%c1bb_to_c1)
if (allocated(nicas_blk%rh_c1)) deallocate(nicas_blk%rh_c1)
if (allocated(nicas_blk%rv_c1)) deallocate(nicas_blk%rv_c1)
if (allocated(nicas_blk%distnorm)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%distnorm(isbb)%dealloc
   end do
   deallocate(nicas_blk%distnorm)
end if
if (allocated(nicas_blk%distnormv)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%distnormv(isbb)%dealloc
   end do
   deallocate(nicas_blk%distnormv)
end if
if (allocated(nicas_blk%rfac)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%rfac(isbb)%dealloc
   end do
   deallocate(nicas_blk%rfac)
end if
if (allocated(nicas_blk%coef)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%coef(isbb)%dealloc
   end do
   deallocate(nicas_blk%coef)
end if
if (allocated(nicas_blk%sc_nor_to_s)) deallocate(nicas_blk%sc_nor_to_s)
if (allocated(nicas_blk%s_to_sc_nor)) deallocate(nicas_blk%s_to_sc_nor)
if (allocated(nicas_blk%sb_to_sc_nor)) deallocate(nicas_blk%sb_to_sc_nor)
call nicas_blk%c_nor%dealloc
call nicas_blk%kdtree%dealloc
if (allocated(nicas_blk%sa_to_sc)) deallocate(nicas_blk%sa_to_sc)
if (allocated(nicas_blk%sb_to_sc)) deallocate(nicas_blk%sb_to_sc)
call nicas_blk%c%dealloc
if (allocated(nicas_blk%h)) then
   do il0=1,geom%nl0i
      call nicas_blk%h(il0)%dealloc
   end do
   deallocate(nicas_blk%h)
end if
call nicas_blk%v%dealloc
if (allocated(nicas_blk%s)) then
   do il1=1,nicas_blk%nl1
     call nicas_blk%s(il1)%dealloc
   end do
   deallocate(nicas_blk%s)
end if
if (allocated(nicas_blk%d)) then
   do its=2,nam%nts
      do il0=1,geom%nl0
        call nicas_blk%d(il0,its)%dealloc
        call nicas_blk%dinv(il0,its)%dealloc
      end do
   end do
   deallocate(nicas_blk%d)
   deallocate(nicas_blk%dinv)
end if
if (allocated(nicas_blk%sb_to_c1b)) deallocate(nicas_blk%sb_to_c1b)
if (allocated(nicas_blk%sb_to_l1)) deallocate(nicas_blk%sb_to_l1)
if (allocated(nicas_blk%norm)) deallocate(nicas_blk%norm)
if (allocated(nicas_blk%coef_ens)) deallocate(nicas_blk%coef_ens)
if (allocated(nicas_blk%sb_to_c1b)) deallocate(nicas_blk%sb_to_c1b)
call nicas_blk%com_AB%dealloc
call nicas_blk%com_AC%dealloc
call nicas_blk%com_AD%dealloc
call nicas_blk%com_ADinv%dealloc

end subroutine nicas_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_parameters
!> Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine nicas_blk_compute_parameters(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng              !< Random number generator
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: il0i,il1

! Compute adaptive sampling
write(mpl%info,'(a7,a)') '','Compute adaptive sampling'
call flush(mpl%info)
call nicas_blk%compute_sampling(mpl,rng,nam,geom,cmat_blk)

! Compute horizontal interpolation data
write(mpl%info,'(a7,a)') '','Compute horizontal interpolation data'
call flush(mpl%info)
call nicas_blk%compute_interp_h(mpl,rng,nam,geom)

! Compute vertical interpolation data
write(mpl%info,'(a7,a)') '','Compute vertical interpolation data'
call flush(mpl%info)
call nicas_blk%compute_interp_v(geom)

! Compute subsampling horizontal interpolation data
write(mpl%info,'(a7,a)') '','Compute subsampling horizontal interpolation data'
call flush(mpl%info)
call nicas_blk%compute_interp_s(mpl,rng,nam,geom)

! Compute MPI distribution, halos A-B
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halos A-B'
call flush(mpl%info)
call nicas_blk%compute_mpi_ab(mpl,geom)

! Compute convolution data
write(mpl%info,'(a7,a)') '','Compute convolution data'
call flush(mpl%info)
call nicas_blk%compute_convol(mpl,rng,nam,geom,cmat_blk)

! Compute MPI distribution, halo C
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo C'
call flush(mpl%info)
call nicas_blk%compute_mpi_c(mpl,geom)

! Compute normalization
write(mpl%info,'(a7,a)') '','Compute normalization'
call flush(mpl%info)
call nicas_blk%compute_normalization(mpl,nam,geom)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%info,'(a10,a,i8)') '','nc0 =        ',geom%nc0
write(mpl%info,'(a10,a,i8)') '','nc0a =       ',geom%nc0a
write(mpl%info,'(a10,a,i8)') '','nl0 =        ',geom%nl0
write(mpl%info,'(a10,a,i8)') '','nc1 =        ',nicas_blk%nc1
write(mpl%info,'(a10,a,i8)') '','nc1a =       ',nicas_blk%nc1a
write(mpl%info,'(a10,a,i8)') '','nc1b =       ',nicas_blk%nc1b
write(mpl%info,'(a10,a,i8)') '','nl1 =        ',nicas_blk%nl1
do il1=1,nicas_blk%nl1
   write(mpl%info,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',nicas_blk%nc2(il1)
end do
write(mpl%info,'(a10,a,i8)') '','ns =         ',nicas_blk%ns
write(mpl%info,'(a10,a,i8)') '','nsa =        ',nicas_blk%nsa
write(mpl%info,'(a10,a,i8)') '','nsb =        ',nicas_blk%nsb
write(mpl%info,'(a10,a,i8)') '','nsc =        ',nicas_blk%nsc
do il0i=1,geom%nl0i
   write(mpl%info,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',nicas_blk%h(il0i)%n_s
end do
write(mpl%info,'(a10,a,i8)') '','v%n_s =      ',nicas_blk%v%n_s
do il1=1,nicas_blk%nl1
   write(mpl%info,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',nicas_blk%s(il1)%n_s
end do
write(mpl%info,'(a10,a,i8)') '','c%n_s =      ',nicas_blk%c%n_s
write(mpl%info,'(a10,a,i8)') '','c_nor%n_s =  ',nicas_blk%c_nor%n_s
call flush(mpl%info)

end subroutine nicas_blk_compute_parameters

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_sampling
!> Purpose: compute NICAS sampling
!----------------------------------------------------------------------
subroutine nicas_blk_compute_sampling(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng              !< Random number generator
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: il0,il0_prev,il1,ic0,ic1,ic2,is,ic0a,iproc
integer :: ncid,nc1_id,nl1_id,lon_c1_id,lat_c1_id,mask_c2_id
integer,allocatable :: c2_to_c1(:)
real(kind_real) :: rhs_sum(geom%nl0),rhs_avg(geom%nl0),rvs_sum(geom%nl0),rvs_avg(geom%nl0),norm(geom%nl0),distnorm(geom%nc0a)
real(kind_real) :: distnormmin,rv,rhs_minavg
real(kind_real),allocatable :: rhs_min(:),rhs_min_glb(:),rhs_c0(:),rhs_c1(:)
real(kind_real),allocatable :: lon_c1(:),lat_c1(:),mask_c2_real(:,:)
logical :: inside
logical,allocatable :: mask_hor_c0(:),mask_c1(:)
character(len=1024) :: filename
character(len=1024) :: subr = 'nicas_blk_compute_sampling'

! Allocation
allocate(rhs_min(geom%nc0a))
if (mpl%main) allocate(rhs_min_glb(geom%nc0))
allocate(nicas_blk%llev(geom%nl0))

! Reset random numbers seed
if (trim(nam%strategy)=='specific_multivariate') call rng%reseed(mpl)

! Compute support radii
norm = 1.0/real(geom%nc0_mask,kind_real)
rhs_sum = sum(cmat_blk%rhs,dim=1,mask=geom%mask_c0a)
call mpl%allreduce_sum(rhs_sum,rhs_avg)
rhs_avg = rhs_avg*norm
rvs_sum = sum(cmat_blk%rvs,dim=1,mask=geom%mask_c0a)
call mpl%allreduce_sum(rvs_sum,rvs_avg)
rvs_avg = rvs_avg*norm
write(mpl%info,'(a10,a)') '','Average support radii (H/V): '
do il0=1,geom%nl0
   write(mpl%info,'(a13,a,i3,a,f10.2,a,f10.2,a)') '','Level ',nam%levs(il0),': '//trim(mpl%aqua),rhs_avg(il0)*reqkm, &
 & trim(mpl%black)//' km  / '//trim(mpl%aqua),rvs_avg(il0),trim(mpl%black)//' '//trim(mpl%vunitchar)
end do
call flush(mpl%info)

! Basic horizontal mesh defined with the minimum support radius
norm(1) = 1.0/real(count(geom%mask_hor_c0),kind_real)
rhs_min = huge(1.0)
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   do il0=1,geom%nl0
      if (geom%mask_c0(ic0,il0)) then
         rhs_min(ic0a) = min(cmat_blk%rhs(ic0a,il0),rhs_min(ic0a))
      end if
   end do
end do
call mpl%allreduce_sum(sum(rhs_min,mask=geom%mask_hor_c0a),rhs_minavg)
rhs_minavg = rhs_minavg*norm(1)
if (rhs_minavg>0.0) then
   nicas_blk%nc1 = floor(2.0*maxval(geom%area)*nam%resol**2/(sqrt(3.0)*rhs_minavg**2))
else
   nicas_blk%nc1 = geom%nc0
end if
nicas_blk%nc1 = min(nicas_blk%nc1,geom%nc0)
write(mpl%info,'(a10,a,i8)') '','Estimated nc1 from horizontal support radius: ',nicas_blk%nc1
call flush(mpl%info)
if (nicas_blk%nc1>nc1max) then
   call mpl%warning('required nc1 larger than nc1max, resetting to nc1max')
   nicas_blk%nc1 = nc1max
   write(mpl%info,'(a10,a,f5.2)') '','Effective horizontal resolution: ',sqrt(real(nicas_blk%nc1,kind_real)*sqrt(3.0) &
 & *rhs_minavg**2/(2.0*maxval(geom%area)))
   call flush(mpl%info)
end if
if (nicas_blk%nc1<3) call mpl%abort('nicas_blk%nc1 lower than 3')

! Allocation
allocate(nicas_blk%c1_to_c0(nicas_blk%nc1))

! Communication
call mpl%loc_to_glb(geom%nc0a,rhs_min,geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.false.,rhs_min_glb)

! Compute subset
write(mpl%info,'(a7,a)',advance='no') '','Compute horizontal subset C1: '
call flush(mpl%info)

! Allocation
allocate(mask_hor_c0(geom%nc0))
mask_hor_c0 = geom%mask_hor_c0

if (test_no_point) then
   ! Mask points on the last MPI task
   if (mpl%nproc==1) call mpl%abort('at least 2 MPI tasks required for test_no_point')
   do ic0=1,geom%nc0
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%nproc) mask_hor_c0(ic0) = .false.
   end do
end if

! Compute subsampling
call rng%initialize_sampling(mpl,geom%nc0,geom%lon,geom%lat,mask_hor_c0,rhs_min_glb,nam%ntry,nam%nrep, &
 & nicas_blk%nc1,nicas_blk%c1_to_c0)
nicas_blk%c1_to_proc = geom%c0_to_proc(nicas_blk%c1_to_c0)

! Inverse conversion
allocate(nicas_blk%c0_to_c1(geom%nc0))
call msi(nicas_blk%c0_to_c1)
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   nicas_blk%c0_to_c1(ic0) = ic1
end do

! Vertical sampling
write(mpl%info,'(a7,a)',advance='no') '','Compute vertical subset L1: '
call flush(mpl%info)
il0_prev = 1
do il0=1,geom%nl0
   ! Look for convolution levels
   if ((il0==1).or.(il0==geom%nl0)) then
      ! Keep first and last levels
      nicas_blk%llev(il0) = .true.
   else
      ! Compute minimum normalized distance with level il0_prev
      distnorm = huge(1.0)
      do ic0a=1,geom%nc0a
         ic0 = geom%c0a_to_c0(ic0a)
         if (geom%mask_c0(ic0,il0)) then
            rv = sqrt(0.5*(cmat_blk%rvs(ic0a,il0)**2+cmat_blk%rvs(ic0a,il0_prev)**2))
            if (rv>0.0) distnorm(ic0a) = abs(geom%vunit(ic0,il0)-geom%vunit(ic0,il0_prev))/rv
         end if
      end do
      call mpl%allreduce_min(minval(distnorm),distnormmin)
      nicas_blk%llev(il0) = distnormmin>1.0/nam%resol
   end if

   ! Update
   if (nicas_blk%llev(il0)) il0_prev = il0
end do
nicas_blk%nl1 = count(nicas_blk%llev)
allocate(nicas_blk%l1_to_l0(nicas_blk%nl1))
il1 = 0
do il0=1,geom%nl0
   if (nicas_blk%llev(il0)) then
      write(mpl%info,'(i3,a)',advance='no') nam%levs(il0),' '
      call flush(mpl%info)
      il1 = il1+1
      nicas_blk%l1_to_l0(il1) = il0
   end if
end do
write(mpl%info,'(a)') ''
call flush(mpl%info)

! Find bottom and top for each point of S1
allocate(nicas_blk%vbot(nicas_blk%nc1))
allocate(nicas_blk%vtop(nicas_blk%nc1))
!$omp parallel do schedule(static) private(ic1,ic0,inside,il1,il0)
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   inside = .false.
   nicas_blk%vtop(ic1) = geom%nl0
   do il1=1,nicas_blk%nl1
      il0 = nicas_blk%l1_to_l0(il1)
      if (.not.inside.and.geom%mask_c0(ic0,il0)) then
         ! Bottom level
         nicas_blk%vbot(ic1) = il0
         inside = .true.
      end if
      if (inside.and.(.not.geom%mask_c0(ic0,il0))) then
         ! Top level
         nicas_blk%vtop(ic1) = il0
         inside = .false.
      end if
   end do
   if (nicas_blk%vbot(ic1)>nicas_blk%vtop(ic1)) call mpl%abort('non contiguous mask')
end do
!$omp end parallel do

! Inverse conversion
allocate(nicas_blk%l0_to_l1(geom%nl0))
call msi(nicas_blk%l0_to_l1)
do il1=1,nicas_blk%nl1
   il0 = nicas_blk%l1_to_l0(il1)
   nicas_blk%l0_to_l1(il0) = il1
end do

! Allocation
allocate(nicas_blk%nc2(nicas_blk%nl1))
allocate(nicas_blk%mask_c2(nicas_blk%nc1,nicas_blk%nl1))
allocate(lon_c1(nicas_blk%nc1))
allocate(lat_c1(nicas_blk%nc1))
allocate(mask_c1(nicas_blk%nc1))
allocate(rhs_c0(geom%nc0))
allocate(rhs_c1(nicas_blk%nc1))

! Horizontal subsampling
do il1=1,nicas_blk%nl1
   write(mpl%info,'(a7,a,i3,a)',advance='no') '','Compute horizontal subset C2 (level ',il1,'): '
   call flush(mpl%info)

   ! Compute nc2
   il0 = nicas_blk%l1_to_l0(il1)
   nicas_blk%nc2(il1) = floor(2.0*geom%area(il0)*nam%resol**2/(sqrt(3.0)*rhs_avg(il0)**2))
   nicas_blk%nc2(il1) = min(nicas_blk%nc2(il1),nicas_blk%nc1)
   if (nicas_blk%nc2(il1)<3) call mpl%abort('nicas_blk%nc2 lower than 3')

   if (nicas_blk%nc2(il1)<nicas_blk%nc1) then
      ! Allocation
      allocate(c2_to_c1(nicas_blk%nc2(il1)))

      ! Compute subset
      call mpl%loc_to_glb(geom%nc0a,cmat_blk%rhs(:,il0),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.false.,rhs_c0)

      ! Initialization
      lon_c1 = geom%lon(nicas_blk%c1_to_c0)
      lat_c1 = geom%lat(nicas_blk%c1_to_c0)
      mask_c1 = geom%mask_c0(nicas_blk%c1_to_c0,il0)
      rhs_c1 = rhs_c0(nicas_blk%c1_to_c0)

      ! Initialize sampling
      call rng%initialize_sampling(mpl,nicas_blk%nc1,lon_c1,lat_c1,mask_c1,rhs_c1,nam%ntry,nam%nrep,nicas_blk%nc2(il1),c2_to_c1)

      ! Fill C2 mask
      nicas_blk%mask_c2(:,il1) = .false.
      do ic2=1,nicas_blk%nc2(il1)
         ic1 = c2_to_c1(ic2)
         nicas_blk%mask_c2(ic1,il1) = .true.
      end do

      ! Release memory
      deallocate(c2_to_c1)
   else
      if (mpl%main) then
         write(mpl%info,'(a)') 'use all C1 points'
         call flush(mpl%info)
      end if

      ! Fill C2 mask
      nicas_blk%mask_c2(:,il1) = .true.
   end if
end do

! Size
nicas_blk%ns = sum(nicas_blk%nc2)

! Conversions
allocate(nicas_blk%s_to_c1(nicas_blk%ns))
allocate(nicas_blk%s_to_l1(nicas_blk%ns))
allocate(nicas_blk%s_to_proc(nicas_blk%ns))
is = 0
do il1=1,nicas_blk%nl1
   do ic1=1,nicas_blk%nc1
      if (nicas_blk%mask_c2(ic1,il1)) then
         is = is+1
         nicas_blk%s_to_c1(is) = ic1
         nicas_blk%s_to_l1(is) = il1
         nicas_blk%s_to_proc(is) = nicas_blk%c1_to_proc(ic1)
      end if
   end do
end do

! Conversions
allocate(nicas_blk%c1l1_to_s(nicas_blk%nc1,nicas_blk%nl1))
call msi(nicas_blk%c1l1_to_s)
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)
   nicas_blk%c1l1_to_s(ic1,il1) = is
end do

! Write grids
if (mpl%main.and.write_grids) then
   ! Allocation
   allocate(lon_c1(nicas_blk%nc1))
   allocate(lat_c1(nicas_blk%nc1))
   allocate(mask_c2_real(nicas_blk%nc1,nicas_blk%nl1))

   ! Copy data
   lon_c1 = geom%lon(nicas_blk%c1_to_c0)*rad2deg
   lat_c1 = geom%lat(nicas_blk%c1_to_c0)*rad2deg
   call msr(mask_c2_real)
   do il1=1,nicas_blk%nl1
      do ic1=1,nicas_blk%nc1
         if (nicas_blk%mask_c2(ic1,il1)) mask_c2_real(ic1,il1) = 1.0
      end do
   end do

   ! Create file
   filename = trim(nam%prefix)//'_'//trim(nicas_blk%name)//'_grids.nc'
   call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Define dimensions
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1',nicas_blk%nc1,nc1_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nl1',nicas_blk%nl1,nl1_id))

   ! Define variables
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_c1',ncfloat,(/nc1_id/),lon_c1_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_c1',ncfloat,(/nc1_id/),lat_c1_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'mask_c2',ncfloat,(/nc1_id,nl1_id/),mask_c2_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,mask_c2_id,'_FillValue',msvalr))

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write variables
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c1_id,lon_c1))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c1_id,lat_c1))
   call mpl%ncerr(subr,nf90_put_var(ncid,mask_c2_id,mask_c2_real))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(lon_c1)
   deallocate(lat_c1)
   deallocate(mask_c2_real)
end if

end subroutine nicas_blk_compute_sampling

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_interp_h
!> Purpose: compute basic horizontal interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_compute_interp_h(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng              !< Random number generator
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il0i
type(linop_type) :: hbase

! Allocation
allocate(nicas_blk%hfull(geom%nl0i))

do il0i=1,geom%nl0i
   ! Compute grid interpolation
   call nicas_blk%hfull(il0i)%interp(mpl,rng,geom,il0i,nicas_blk%nc1,nicas_blk%c1_to_c0,nam%mask_check, &
 & nicas_blk%vbot,nicas_blk%vtop,nam%nicas_interp,hbase)
end do

end subroutine nicas_blk_compute_interp_h

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_interp_v
!> Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_compute_interp_v(nicas_blk,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: ic1,ic0,jl0,il0,jl1,il0inf,il0sup

! Initialize vertical interpolation
nicas_blk%vfull%prefix = 'vfull'
nicas_blk%vfull%n_src = nicas_blk%nl1
nicas_blk%vfull%n_dst = geom%nl0

! Linear interpolation
nicas_blk%vfull%n_s = nicas_blk%nl1
il0inf = 1
do jl0=1,geom%nl0
   if (nicas_blk%llev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         nicas_blk%vfull%n_s = nicas_blk%vfull%n_s+2
      end do
      il0inf = jl0
   end if
end do
call nicas_blk%vfull%alloc(nicas_blk%nc1)

do jl1=1,nicas_blk%nl1
   jl0 = nicas_blk%l1_to_l0(jl1)
   nicas_blk%vfull%row(jl1) = jl0
   nicas_blk%vfull%col(jl1) = jl0
   do ic1=1,nicas_blk%nc1
      nicas_blk%vfull%Svec(jl1,ic1) = 1.0
   end do
end do

nicas_blk%vfull%n_s = nicas_blk%nl1
il0inf = 1
do jl0=1,geom%nl0
   if (nicas_blk%llev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         nicas_blk%vfull%n_s = nicas_blk%vfull%n_s+1
         nicas_blk%vfull%row(nicas_blk%vfull%n_s) = il0
         nicas_blk%vfull%col(nicas_blk%vfull%n_s) = il0inf
         do ic1=1,nicas_blk%nc1
            ic0 = nicas_blk%c1_to_c0(ic1)
            nicas_blk%vfull%Svec(nicas_blk%vfull%n_s,ic1) = abs(geom%vunit(ic0,il0sup)-geom%vunit(ic0,il0)) &
                                                          & /abs(geom%vunit(ic0,il0sup)-geom%vunit(ic0,il0inf))
         end do

         nicas_blk%vfull%n_s = nicas_blk%vfull%n_s+1
         nicas_blk%vfull%row(nicas_blk%vfull%n_s) = il0
         nicas_blk%vfull%col(nicas_blk%vfull%n_s) = il0sup
         do ic1=1,nicas_blk%nc1
            ic0 = nicas_blk%c1_to_c0(ic1)
            nicas_blk%vfull%Svec(nicas_blk%vfull%n_s,ic1) = abs(geom%vunit(ic0,il0)-geom%vunit(ic0,il0inf)) &
                                                          & /abs(geom%vunit(ic0,il0sup)-geom%vunit(ic0,il0inf))
         end do
      end do
      il0inf = jl0
   end if
end do

! Conversion
nicas_blk%vfull%col = nicas_blk%l0_to_l1(nicas_blk%vfull%col)

! Deallocate selected levels
deallocate(nicas_blk%llev)

end subroutine nicas_blk_compute_interp_v

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_interp_s
!> Purpose: compute horizontal subsampling interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_compute_interp_s(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng              !< Random number generator
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il1,i_s
real(kind_real) :: renorm(nicas_blk%nc1)
real(kind_real) :: lon_c1(nicas_blk%nc1),lat_c1(nicas_blk%nc1)
logical :: mask_src(nicas_blk%nc1),mask_dst(nicas_blk%nc1)
logical,allocatable :: valid(:)
type(linop_type) :: stmp

! Allocation
allocate(nicas_blk%sfull(nicas_blk%nl1))

do il1=1,nicas_blk%nl1
   ! Initialize NICAS block data
   write(nicas_blk%sfull(il1)%prefix,'(a,i3.3)') 's_',il1
   nicas_blk%sfull(il1)%n_src = nicas_blk%nc1
   nicas_blk%sfull(il1)%n_dst = nicas_blk%nc1

   if (nicas_blk%nc2(il1)==nicas_blk%nc1) then
      ! No interpolation
      nicas_blk%sfull(il1)%n_s = nicas_blk%nc1
      call nicas_blk%sfull(il1)%alloc
      do i_s=1,nicas_blk%sfull(il1)%n_s
         nicas_blk%sfull(il1)%row(i_s) = i_s
         nicas_blk%sfull(il1)%col(i_s) = i_s
         nicas_blk%sfull(il1)%S(i_s) = 1.0
      end do
   else
      ! Initialization
      lon_c1 = geom%lon(nicas_blk%c1_to_c0)
      lat_c1 = geom%lat(nicas_blk%c1_to_c0)
      mask_src = nicas_blk%mask_c2(:,il1)
      mask_dst = .true.

      ! Compute interpolation
      call stmp%interp(mpl,rng,nicas_blk%nc1,lon_c1,lat_c1,mask_src,nicas_blk%nc1,lon_c1,lat_c1,mask_dst,nam%nicas_interp)

      ! Allocation
      allocate(valid(stmp%n_s))
      valid = .true.

      ! Check mask boundaries
      if (nam%mask_check) then
         write(mpl%info,'(a10,a,i3,a)',advance='no') '','Sublevel ',il1,': '
         call flush(mpl%info)
         call stmp%interp_check_mask(mpl,geom,valid,nicas_blk%l1_to_l0(il1),row_to_ic0=nicas_blk%c1_to_c0, &
       & col_to_ic0=nicas_blk%c1_to_c0)
      else
         write(mpl%info,'(a10,a,i3)') '','Sublevel ',il1
         call flush(mpl%info)
      end if

      ! Renormalization
      renorm = 0.0
      do i_s=1,stmp%n_s
         if (valid(i_s)) renorm(stmp%row(i_s)) = renorm(stmp%row(i_s))+stmp%S(i_s)
      end do

      ! Copy valid operations
      nicas_blk%sfull(il1)%n_s = count(valid)
      call nicas_blk%sfull(il1)%alloc
      nicas_blk%sfull(il1)%n_s = 0
      do i_s=1,stmp%n_s
         if (valid(i_s)) then
            nicas_blk%sfull(il1)%n_s = nicas_blk%sfull(il1)%n_s+1
            nicas_blk%sfull(il1)%row(nicas_blk%sfull(il1)%n_s) = stmp%row(i_s)
            nicas_blk%sfull(il1)%col(nicas_blk%sfull(il1)%n_s) = stmp%col(i_s)
            nicas_blk%sfull(il1)%S(nicas_blk%sfull(il1)%n_s) = stmp%S(i_s)/renorm(stmp%row(i_s))
         end if
      end do

      ! Release memory
      call stmp%dealloc
      deallocate(valid)

      ! Count points that are not interpolated
      call nicas_blk%sfull(il1)%interp_missing(mpl,nicas_blk%nc1,lon_c1,lat_c1,mask_dst,nam%nicas_interp)
   end if
end do

end subroutine nicas_blk_compute_interp_s

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_ab
!> Purpose: compute NICAS MPI distribution, halos A-B
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_ab(nicas_blk,mpl,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il0i,ic0,ic0a,iproc,ic1,jc1,ic1a,ic1b,il0,il1,isa,isb,i_s,i_s_loc,is,js,h_n_s_max,s_n_s_max,h_n_s_max_loc,s_n_s_max_loc
integer,allocatable :: s_to_proc(:),interph_lg(:,:),interps_lg(:,:)
integer,allocatable :: proc_to_nc1a(:),proc_to_nsa(:)
logical,allocatable :: lcheck_c1a(:),lcheck_c1b_h(:),lcheck_c1b(:),lcheck_h(:,:),lcheck_s(:,:)

! Allocation
h_n_s_max = 0
do il0i=1,geom%nl0i
   h_n_s_max = max(h_n_s_max,nicas_blk%hfull(il0i)%n_s)
end do
s_n_s_max = 0
do il1=1,nicas_blk%nl1
   s_n_s_max = max(s_n_s_max,nicas_blk%sfull(il1)%n_s)
end do
allocate(nicas_blk%h(geom%nl0i))
allocate(nicas_blk%s(nicas_blk%nl1))
allocate(lcheck_c1a(nicas_blk%nc1))
allocate(lcheck_c1b_h(nicas_blk%nc1))
allocate(lcheck_c1b(nicas_blk%nc1))
allocate(nicas_blk%lcheck_sa(nicas_blk%ns))
allocate(nicas_blk%lcheck_sb(nicas_blk%ns))
allocate(lcheck_h(h_n_s_max,geom%nl0i))
allocate(lcheck_s(s_n_s_max,nicas_blk%nl1))
allocate(proc_to_nc1a(mpl%nproc))
allocate(proc_to_nsa(mpl%nproc))
allocate(s_to_proc(nicas_blk%ns))

! Halo definitions

! Halo A
lcheck_c1a = .false.
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) lcheck_c1a(ic1) = .true.
end do
nicas_blk%lcheck_sa = .false.
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   ic0 = nicas_blk%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) nicas_blk%lcheck_sa(is) = .true.
end do

! Halo B

! Horizontal interpolation
lcheck_h = .false.
lcheck_c1b_h = .false.
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%hfull(il0i)%n_s
      ic0 = nicas_blk%hfull(il0i)%row(i_s)
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%myproc) then
         jc1 = nicas_blk%hfull(il0i)%col(i_s)
         lcheck_h(i_s,il0i) = .true.
         lcheck_c1b_h(jc1) = .true.
      end if
   end do
end do

! Subsampling horizontal interpolation
nicas_blk%lcheck_sb = .false.
lcheck_s = .false.
lcheck_c1b = lcheck_c1b_h
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%sfull(il1)%n_s
      ic1 = nicas_blk%sfull(il1)%row(i_s)
      if (lcheck_c1b_h(ic1)) then
         jc1 = nicas_blk%sfull(il1)%col(i_s)
         js = nicas_blk%c1l1_to_s(jc1,il1)
         lcheck_c1b(jc1) = .true.
         nicas_blk%lcheck_sb(js) = .true.
         lcheck_s(i_s,il1) = .true.
      end if
   end do
end do

! Check halos consistency
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is).and.(.not.nicas_blk%lcheck_sb(is))) call mpl%abort('point in halo A but not in halo B')
end do

! Sizes
nicas_blk%nc1a = count(lcheck_c1a)
call mpl%allgather(nicas_blk%nc1a,proc_to_nc1a)
nicas_blk%nsa = count(nicas_blk%lcheck_sa)
call mpl%allgather(nicas_blk%nsa,proc_to_nsa)
do il0i=1,geom%nl0i
   nicas_blk%h(il0i)%n_s = count(lcheck_h(:,il0i))
end do
nicas_blk%nc1b = count(lcheck_c1b)
nicas_blk%nsb = count(nicas_blk%lcheck_sb)
do il1=1,nicas_blk%nl1
   nicas_blk%s(il1)%n_s = count(lcheck_s(:,il1))
end do

! Global <-> local conversions for fields

! Halo A
allocate(nicas_blk%c1a_to_c1(nicas_blk%nc1a))
allocate(nicas_blk%c1_to_c1a(nicas_blk%nc1))
ic1a = 0
do ic1=1,nicas_blk%nc1
   if (lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      nicas_blk%c1a_to_c1(ic1a) = ic1
   end if
end do
call mpl%glb_to_loc_index(nicas_blk%nc1a,nicas_blk%c1a_to_c1,nicas_blk%nc1,nicas_blk%c1_to_c1a)

allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
allocate(nicas_blk%s_to_sa(nicas_blk%ns))
isa = 0
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is)) then
      isa = isa+1
      nicas_blk%sa_to_s(isa) = is
   end if
end do
call mpl%glb_to_loc_index(nicas_blk%nsa,nicas_blk%sa_to_s,nicas_blk%ns,nicas_blk%s_to_sa)

! Halo B
allocate(nicas_blk%c1b_to_c1(nicas_blk%nc1b))
allocate(nicas_blk%c1_to_c1b(nicas_blk%nc1))
call msi(nicas_blk%c1_to_c1b)
ic1b = 0
do ic1=1,nicas_blk%nc1
   if (lcheck_c1b(ic1)) then
      ic1b = ic1b+1
      nicas_blk%c1b_to_c1(ic1b) = ic1
      nicas_blk%c1_to_c1b(ic1) = ic1b
   end if
end do

allocate(nicas_blk%sb_to_s(nicas_blk%nsb))
allocate(nicas_blk%s_to_sb(nicas_blk%ns))
call msi(nicas_blk%s_to_sb)
isb = 0
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sb(is)) then
      isb = isb+1
      nicas_blk%sb_to_s(isb) = is
      nicas_blk%s_to_sb(is) = isb
   end if
end do

! Inter-halo conversions
allocate(nicas_blk%sa_to_sb(nicas_blk%nsa))
do isa=1,nicas_blk%nsa
   is = nicas_blk%sa_to_s(isa)
   isb = nicas_blk%s_to_sb(is)
   nicas_blk%sa_to_sb(isa) = isb
end do

! Global <-> local conversions for data
h_n_s_max_loc = 0
do il0i=1,geom%nl0i
   h_n_s_max_loc = max(h_n_s_max_loc,nicas_blk%h(il0i)%n_s)
end do
allocate(interph_lg(h_n_s_max_loc,geom%nl0i))
do il0i=1,geom%nl0i
   i_s_loc = 0
   do i_s=1,nicas_blk%hfull(il0i)%n_s
      if (lcheck_h(i_s,il0i)) then
         i_s_loc = i_s_loc+1
         interph_lg(i_s_loc,il0i) = i_s
      end if
   end do
end do
s_n_s_max_loc = 0
do il1=1,nicas_blk%nl1
   s_n_s_max_loc = max(s_n_s_max_loc,nicas_blk%s(il1)%n_s)
end do
allocate(interps_lg(s_n_s_max_loc,nicas_blk%nl1))
do il1=1,nicas_blk%nl1
   i_s_loc = 0
   do i_s=1,nicas_blk%sfull(il1)%n_s
      if (lcheck_s(i_s,il1)) then
         i_s_loc = i_s_loc+1
         interps_lg(i_s_loc,il1) = i_s
      end if
   end do
end do

! Local data

! Horizontal interpolation
do il0i=1,geom%nl0i
   write(nicas_blk%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
   nicas_blk%h(il0i)%n_src = nicas_blk%nc1b
   nicas_blk%h(il0i)%n_dst = geom%nc0a
   call nicas_blk%h(il0i)%alloc
   do i_s_loc=1,nicas_blk%h(il0i)%n_s
      i_s = interph_lg(i_s_loc,il0i)
      nicas_blk%h(il0i)%row(i_s_loc) = geom%c0_to_c0a(nicas_blk%hfull(il0i)%row(i_s))
      nicas_blk%h(il0i)%col(i_s_loc) = nicas_blk%c1_to_c1b(nicas_blk%hfull(il0i)%col(i_s))
      nicas_blk%h(il0i)%S(i_s_loc) = nicas_blk%hfull(il0i)%S(i_s)
   end do
   call nicas_blk%h(il0i)%reorder(mpl)
end do

! Vertical interpolation
nicas_blk%v%prefix = 'v'
nicas_blk%v%n_src = nicas_blk%vfull%n_src
nicas_blk%v%n_dst = nicas_blk%vfull%n_dst
nicas_blk%v%n_s = nicas_blk%vfull%n_s
call nicas_blk%v%alloc(nicas_blk%nc1b)
do i_s=1,nicas_blk%v%n_s
   nicas_blk%v%row(i_s) = nicas_blk%vfull%row(i_s)
   nicas_blk%v%col(i_s) = nicas_blk%vfull%col(i_s)
   do ic1b=1,nicas_blk%nc1b
      ic1 = nicas_blk%c1b_to_c1(ic1b)
      nicas_blk%v%Svec(i_s,ic1b) = nicas_blk%vfull%Svec(i_s,ic1)
   end do
end do
call nicas_blk%v%reorder(mpl)

! Subsampling horizontal interpolation
do il1=1,nicas_blk%nl1
   write(nicas_blk%s(il1)%prefix,'(a,i3.3)') 's_',il1
   nicas_blk%s(il1)%n_src = nicas_blk%nc1b
   nicas_blk%s(il1)%n_dst = nicas_blk%nc1b
   call nicas_blk%s(il1)%alloc
   do i_s_loc=1,nicas_blk%s(il1)%n_s
      i_s = interps_lg(i_s_loc,il1)
      nicas_blk%s(il1)%row(i_s_loc) = nicas_blk%c1_to_c1b(nicas_blk%sfull(il1)%row(i_s))
      nicas_blk%s(il1)%col(i_s_loc) = nicas_blk%c1_to_c1b(nicas_blk%sfull(il1)%col(i_s))
      nicas_blk%s(il1)%S(i_s_loc) = nicas_blk%sfull(il1)%S(i_s)
   end do
   call nicas_blk%s(il1)%reorder(mpl)
end do

! Conversions
allocate(nicas_blk%c1a_to_c0a(nicas_blk%nc1a))
allocate(nicas_blk%sa_to_c0a(nicas_blk%nsa))
allocate(nicas_blk%sa_to_l0(nicas_blk%nsa))
allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
allocate(nicas_blk%c1bl1_to_sb(nicas_blk%nc1b,nicas_blk%nl1))
call msi(nicas_blk%c1bl1_to_sb)
do ic1a=1,nicas_blk%nc1a
   ic1 = nicas_blk%c1a_to_c1(ic1a)
   ic0 = nicas_blk%c1_to_c0(ic1)
   ic0a = geom%c0_to_c0a(ic0)
   nicas_blk%c1a_to_c0a(ic1a) = ic0a
end do
do isa=1,nicas_blk%nsa
   is = nicas_blk%sa_to_s(isa)
   ic1 = nicas_blk%s_to_c1(is)
   ic0 = nicas_blk%c1_to_c0(ic1)
   ic0a = geom%c0_to_c0a(ic0)
   il1 = nicas_blk%s_to_l1(is)
   il0 = nicas_blk%l1_to_l0(il1)
   nicas_blk%sa_to_c0a(isa) = ic0a
   nicas_blk%sa_to_l0(isa) = il0
end do
do isb=1,nicas_blk%nsb
   is = nicas_blk%sb_to_s(isb)
   il1 = nicas_blk%s_to_l1(is)
   ic1 = nicas_blk%s_to_c1(is)
   ic1b = nicas_blk%c1_to_c1b(ic1)
   nicas_blk%sb_to_c1b(isb) = ic1b
   nicas_blk%sb_to_l1(isb) = il1
   nicas_blk%c1bl1_to_sb(ic1b,il1) = isb
end do

! Setup communications
s_to_proc = geom%c0_to_proc(nicas_blk%c1_to_c0(nicas_blk%s_to_c1))
call nicas_blk%com_AB%setup(mpl,'com_AB',nicas_blk%ns,nicas_blk%nsa,nicas_blk%nsb,nicas_blk%sb_to_s,nicas_blk%sa_to_sb, &
 & s_to_proc,nicas_blk%s_to_sa)

end subroutine nicas_blk_compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol
!> Purpose: compute convolution
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng              !< Random number generator
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: n_s_max,ithread,is,ic1,jc1,il1,il0,j,js,isb,ic1b,ic0,ic0a,ic1a,i_s,jc,kc,ks,jbd,jc0,jl0,jl1,ic1bb,isbb
integer :: c_n_s(mpl%nthread)
integer,allocatable :: nn(:),nn_index(:),inec(:),c_ind(:,:)
real(kind_real) :: distvsq,rvsq
real(kind_real),allocatable :: lon_c1(:),lat_c1(:),nn_dist(:)
real(kind_real),allocatable :: rh_c1a(:,:),rv_c1a(:,:),rv_rfac_c1a(:,:),rv_coef_c1a(:,:)
real(kind_real),allocatable :: rv_rfac_c1(:,:),rv_coef_c1(:,:)
real(kind_real),allocatable :: distnormv(:,:),rfac(:,:),coef(:,:)
real(kind_real),allocatable :: c_S(:,:),c_S_conv(:)
logical :: add_op
logical,allocatable :: lcheck_c1bb(:)
type(linop_type) :: ctmp,c(mpl%nthread)

! Associate
associate(ib=>nicas_blk%ib)

! Set double-fit parameter
nicas_blk%double_fit = cmat_blk%double_fit

! Allocation
allocate(lon_c1(nicas_blk%nc1))
allocate(lat_c1(nicas_blk%nc1))
allocate(nicas_blk%mask_c1(nicas_blk%nc1,nicas_blk%nl1))
allocate(lcheck_c1bb(nicas_blk%nc1))

! Initialization
lon_c1 = geom%lon(nicas_blk%c1_to_c0)
lat_c1 = geom%lat(nicas_blk%c1_to_c0)
nicas_blk%mask_c1 = .false.
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)
   nicas_blk%mask_c1(ic1,il1) = .true.
end do

! Compute KD-tree
write(mpl%info,'(a10,a)') '','Compute KD-tree'
call flush(mpl%info)
call nicas_blk%kdtree%create(mpl,nicas_blk%nc1,lon_c1,lat_c1)

! Find largest possible radius
call mpl%allreduce_max(maxval(cmat_blk%rh),nicas_blk%rhmax)
if (nicas_blk%double_fit) then
   nicas_blk%rhmax = nicas_blk%rhmax*sqrt_r_dble
else
   nicas_blk%rhmax = nicas_blk%rhmax*sqrt_r
end if

if (nam%lsqrt) then
   ! Copy
   nicas_blk%nc1bb = nicas_blk%nc1b
else
   write(mpl%info,'(a10,a)',advance='no') '','Define extended halo: '
   call flush(mpl%info)

   ! Allocation
   allocate(nn(nicas_blk%nc1b))

   do ic1b=1,nicas_blk%nc1b
      ! Indices
      ic1 = nicas_blk%c1b_to_c1(ic1b)
      ic0 = nicas_blk%c1_to_c0(ic1)

      ! Count nearest neighbors
      call nicas_blk%kdtree%count_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nicas_blk%rhmax,nn(ic1b))
   end do

   ! Initialization
   call mpl%prog_init(nicas_blk%nc1b)
   lcheck_c1bb = .false.

   do ic1b=1,nicas_blk%nc1b
      ! Indices
      ic1 = nicas_blk%c1b_to_c1(ic1b)
      ic0 = nicas_blk%c1_to_c0(ic1)

      ! Allocation
      allocate(nn_index(nn(ic1b)))
      allocate(nn_dist(nn(ic1b)))

      ! Find nearest neighbors
      call nicas_blk%kdtree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nn(ic1b),nn_index,nn_dist)

      ! Fill mask
      do j=1,nn(ic1b)
         jc1 = nn_index(j)
         lcheck_c1bb(jc1) = .true.
      end do

      ! Release memory
      deallocate(nn_index)
      deallocate(nn_dist)

      ! Update
      call mpl%prog_print(ic1b)
   end do
   write(mpl%info,'(a)') '100%'
   call flush(mpl%info)

   ! Halo size
   nicas_blk%nc1bb = count(lcheck_c1bb)
   write(mpl%info,'(a10,a,i6,a,i6)') '','Halo sizes nc1b / nc1bb: ',nicas_blk%nc1b,' / ',nicas_blk%nc1bb
end if

! Allocation
allocate(nicas_blk%c1bb_to_c1(nicas_blk%nc1bb))
allocate(nicas_blk%c1_to_c1bb(nicas_blk%nc1))

if (nam%lsqrt) then
   ! Copy
   nicas_blk%c1bb_to_c1 = nicas_blk%c1b_to_c1
   nicas_blk%c1_to_c1bb = nicas_blk%c1_to_c1b
   nicas_blk%nsbb = nicas_blk%nsb
else
   ! Global <-> local conversions for fields
   call msi(nicas_blk%c1_to_c1bb)
   ic1bb = 0
   do ic1=1,nicas_blk%nc1
      if (lcheck_c1bb(ic1)) then
         ic1bb = ic1bb+1
         nicas_blk%c1bb_to_c1(ic1bb) = ic1
         nicas_blk%c1_to_c1bb(ic1) = ic1bb
      end if
   end do

   ! Count points in extended halo
   nicas_blk%nsbb = 0
   do is=1,nicas_blk%ns
      ic1 = nicas_blk%s_to_c1(is)
      if (lcheck_c1bb(ic1)) nicas_blk%nsbb = nicas_blk%nsbb+1
   end do
   write(mpl%info,'(a10,a,i6,a,i6)') '','Halo sizes nsb / nsbb:   ',nicas_blk%nsb,' / ',nicas_blk%nsbb
end if

! Allocation
allocate(nicas_blk%sbb_to_s(nicas_blk%nsbb))

if (nam%lsqrt) then
   ! Copy
   nicas_blk%sbb_to_s = nicas_blk%sb_to_s
else
   ! Global <-> local conversions for fields
   isbb = 0
   do is=1,nicas_blk%ns
      ic1 = nicas_blk%s_to_c1(is)
      if (lcheck_c1bb(ic1)) then
         isbb = isbb+1
         nicas_blk%sbb_to_s(isbb) = is
      end if
   end do
end if

! Compute horizontal and vertical parameters
write(mpl%info,'(a10,a)') '','Compute horizontal and vertical parameters'
call flush(mpl%info)

! Allocation
allocate(rh_c1a(nicas_blk%nc1a,nicas_blk%nl1))
allocate(rv_c1a(nicas_blk%nc1a,nicas_blk%nl1))
allocate(nicas_blk%rh_c1(nicas_blk%nc1,nicas_blk%nl1))
allocate(nicas_blk%rv_c1(nicas_blk%nc1,nicas_blk%nl1))

! Copy and rescale
write(mpl%info,'(a13,a)') '','Copy and rescale'
call flush(mpl%info)
do il1=1,nicas_blk%nl1
   do ic1a=1,nicas_blk%nc1a
      ! Copy
      ic0a = nicas_blk%c1a_to_c0a(ic1a)
      il0 = nicas_blk%l1_to_l0(il1)
      rh_c1a(ic1a,il1) = cmat_blk%rh(ic0a,il0)
      rv_c1a(ic1a,il1) = cmat_blk%rv(ic0a,il0)

      ! Square-root rescaling
      rh_c1a(ic1a,il1) = rh_c1a(ic1a,il1)*sqrt_r
      if (nicas_blk%double_fit) then
         rv_c1a(ic1a,il1) = rv_c1a(ic1a,il1)*sqrt_r_dble
      else
         rv_c1a(ic1a,il1) = rv_c1a(ic1a,il1)*sqrt_r
      end if
   end do
end do

! Communication
write(mpl%info,'(a13,a)') '','Communication'
call flush(mpl%info)
call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rh_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%rh_c1)
call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rv_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%rv_c1)

! Allocation
allocate(nicas_blk%distnorm(nicas_blk%nsbb))

! Compute distances
if (nam%network) then
   call nicas_blk%compute_convol_network(mpl,rng,nam,geom)
else
   call nicas_blk%compute_convol_distance(mpl,geom)
end if

! Release memory
deallocate(nicas_blk%rh_c1)
deallocate(nicas_blk%rv_c1)
deallocate(nicas_blk%mask_c1)

if (nicas_blk%double_fit) then
   ! Compute double-fit parameters
   write(mpl%info,'(a10,a)') '','Compute double-fit parameters'
   call flush(mpl%info)

   ! Allocation
   allocate(rv_rfac_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(rv_coef_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(rv_rfac_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(rv_coef_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(nicas_blk%distnormv(nicas_blk%nsbb))
   allocate(nicas_blk%rfac(nicas_blk%nsbb))
   allocate(nicas_blk%coef(nicas_blk%nsbb))

   ! Copy and rescale
   write(mpl%info,'(a13,a)') '','Copy and rescale'
   call flush(mpl%info)
   do il1=1,nicas_blk%nl1
      do ic1a=1,nicas_blk%nc1a
         ! Copy
         ic0a = nicas_blk%c1a_to_c0a(ic1a)
         il0 = nicas_blk%l1_to_l0(il1)
         rv_rfac_c1a(ic1a,il1) = cmat_blk%rv_rfac(ic0a,il0)
         rv_coef_c1a(ic1a,il1) = cmat_blk%rv_coef(ic0a,il0)

         ! Square-root rescaling
         rv_rfac_c1a(ic1a,il1) = rv_rfac_c1a(ic1a,il1)*sqrt_rfac
         rv_coef_c1a(ic1a,il1) = rv_coef_c1a(ic1a,il1)*sqrt_coef
      end do
   end do

   ! Communication
   write(mpl%info,'(a13,a)') '','Communication'
   call flush(mpl%info)
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rv_rfac_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & rv_rfac_c1)
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rv_coef_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & rv_coef_c1)

   ! Release memory
   deallocate(rv_rfac_c1a)
   deallocate(rv_coef_c1a)

   !$omp parallel do schedule(static) private(isbb,is,ic1,il1,ic0,il0,jbd,jc1,jl1,jc0,jl0,distvsq,rvsq), &
   !$omp&                             firstprivate(distnormv,rfac,coef)
   do isbb=1,nicas_blk%nsbb
      ! Indices
      is = nicas_blk%sbb_to_s(isbb)
      ic1 = nicas_blk%s_to_c1(is)
      il1 = nicas_blk%s_to_l1(is)
      ic0 = nicas_blk%c1_to_c0(ic1)
      il0 = nicas_blk%l1_to_l0(il1)

      ! Allocation
      allocate(distnormv(nicas_blk%nc1,nicas_blk%nl1))
      allocate(rfac(nicas_blk%nc1,nicas_blk%nl1))
      allocate(coef(nicas_blk%nc1,nicas_blk%nl1))

      ! Initialization
      call msr(distnormv)
      call msr(rfac)
      call msr(coef)

      do jbd=1,nicas_blk%distnorm(isbb)%nbd
         ! Indices
         jc1 = nicas_blk%distnorm(isbb)%bd_to_c1(jbd)
         jl1 = nicas_blk%distnorm(isbb)%bd_to_l1(jbd)
         jc0 = nicas_blk%c1_to_c0(jc1)
         jl0 = nicas_blk%l1_to_l0(jl1)

         ! Vertical distance
         distvsq = (geom%vunit(ic0,il0)-geom%vunit(jc0,jl0))**2
         rvsq = 0.5*(nicas_blk%rv_c1(ic1,il1)**2+nicas_blk%rv_c1(jc1,jl1)**2)
         if (rvsq>0.0) then
            distnormv(jc1,jl1) = sqrt(distvsq/rvsq)
         elseif (distvsq>0.0) then
            distnormv(jc1,jl1) = 0.5*huge(0.0)
         end if
         rfac(ic1,il1) = sqrt(rv_rfac_c1(ic1,il1)*rv_rfac_c1(jc1,jl1))
         coef(ic1,il1) = sqrt(rv_coef_c1(ic1,il1)*rv_coef_c1(jc1,jl1))
      end do
      
      ! Pack data
      call nicas_blk%distnormv(isbb)%pack(nicas_blk%nc1,nicas_blk%nl1,distnormv)
      call nicas_blk%rfac(isbb)%pack(nicas_blk%nc1,nicas_blk%nl1,rfac)
      call nicas_blk%coef(isbb)%pack(nicas_blk%nc1,nicas_blk%nl1,coef)

      ! Release memory
      deallocate(distnormv)
      deallocate(rfac)
      deallocate(coef)
   end do
   !$omp end parallel do

   ! Release memory
   deallocate(rv_rfac_c1)
   deallocate(rv_coef_c1)
end if

! Compute weights
call nicas_blk%compute_convol_weights(mpl,nam,geom,ctmp)

! Release memory
do isbb=1,nicas_blk%nsbb
   call nicas_blk%distnorm(isbb)%dealloc
   if (nicas_blk%double_fit) then
      call nicas_blk%distnormv(isbb)%dealloc
      call nicas_blk%rfac(isbb)%dealloc
      call nicas_blk%coef(isbb)%dealloc
   end if
end do
deallocate(nicas_blk%distnorm)
if (nicas_blk%double_fit) then
   deallocate(nicas_blk%distnormv)
   deallocate(nicas_blk%rfac)
   deallocate(nicas_blk%coef)
end if

if (nam%lsqrt) then
   ! Copy
   nicas_blk%c = ctmp%copy()
else
   ! Compute convolution inverse mapping
   allocate(inec(nicas_blk%ns))
   inec = 0
   do i_s=1,ctmp%n_s
      is = ctmp%col(i_s)
      inec(is) = inec(is)+1
   end do
   allocate(c_ind(maxval(inec),nicas_blk%ns))
   allocate(c_S(maxval(inec),nicas_blk%ns))
   call msi(c_ind)
   call msr(c_S)
   inec = 0
   do i_s=1,ctmp%n_s
      is = ctmp%col(i_s)
      js = ctmp%row(i_s)
      inec(is) = inec(is)+1
      c_ind(inec(is),is) = js
      c_S(inec(is),is) = ctmp%S(i_s)
   end do

   ! Initialization
   write(mpl%info,'(a10,a)',advance='no') '','Second pass:     '
   call flush(mpl%info)
   call mpl%prog_init(nicas_blk%nsb)
   n_s_max = 100*nint(real(geom%nc0*geom%nl0)/real(mpl%nthread*mpl%nproc))
   c_n_s = 0
   do ithread=1,mpl%nthread
      c(ithread)%n_s = n_s_max
      call c(ithread)%alloc
      call msi(c(ithread)%row)
      call msi(c(ithread)%col)
      call msr(c(ithread)%S)
   end do

   ! Apply convolution
   !$omp parallel do schedule(static) private(isb,is,ithread,jc,js,kc,ks,add_op) firstprivate(c_S_conv)
   do isb=1,nicas_blk%nsb
      ! Indices
      is = nicas_blk%sb_to_s(isb)
      ithread = 1
!$    ithread = omp_get_thread_num()+1

      ! Allocation
      allocate(c_S_conv(nicas_blk%ns))

      ! Initialization
      c_S_conv = 0.0

      ! Loop twice over points
      do jc=1,inec(is)
         js = c_ind(jc,is)
         do kc=1,inec(js)
            ks = c_ind(kc,js)
            c_S_conv(ks) = c_S_conv(ks)+c_S(jc,is)*c_S(kc,js)
         end do
      end do

      ! Store coefficient for convolution
      do js=1,nicas_blk%ns
         add_op = .false.
         if (nam%mpicom==1) then
            add_op = (nicas_blk%lcheck_sb(js).and.(is<=js)).or.(.not.nicas_blk%lcheck_sb(js))
         elseif (nam%mpicom==2) then
            add_op = nicas_blk%lcheck_sa(is).and.((nicas_blk%lcheck_sa(js).and.(is<=js)) &
                   & .or.(.not.nicas_blk%lcheck_sa(js)))
         end if
         if (add_op) call c(ithread)%add_op(c_n_s(ithread),is,js,c_S_conv(js))
      end do

      ! Release memory
      deallocate(c_S_conv)

      ! Update
      call mpl%prog_print(isb)
   end do
   !$omp end parallel do
   write(mpl%info,'(a)') '100%'
   call flush(mpl%info)

   ! Gather data from threads
   call nicas_blk%c%gather(mpl,c_n_s,c)

   ! Release memory
   do ithread=1,mpl%nthread
      call c(ithread)%dealloc
   end do
end if

! Set prefix
nicas_blk%c%prefix = 'c'
nicas_blk%c_nor%prefix = 'c_nor'

! End associate
end associate

end subroutine nicas_blk_compute_convol

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_network
!> Purpose: compute convolution with a network approach
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_network(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng              !< Random number generator
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: net_nnbmax,is,ic0,ic1,il0,jl1,np,np_new,i,j,k,ip,kc1,jc1,il1,dkl1,kl1,jp,isbb,djl1,inr,jc0,jl0
integer,allocatable :: net_nnb(:),net_inb(:,:),plist(:,:),plist_new(:,:)
real(kind_real) :: distnorm_network,disttest
real(kind_real) :: dnb,disthsq,distvsq,rhsq,rvsq
real(kind_real),allocatable :: distnorm(:,:),net_dnb(:,:,:,:)
logical :: init,valid_arc,add_to_front
type(mesh_type) :: mesh

! Allocation
allocate(net_nnb(nicas_blk%nc1))

! Create mesh
write(mpl%info,'(a10,a)') '','Create mesh'
call flush(mpl%info)
call mesh%create(mpl,rng,nicas_blk%nc1,geom%lon(nicas_blk%c1_to_c0),geom%lat(nicas_blk%c1_to_c0))

! Count neighbors
write(mpl%info,'(a10,a)') '','Count neighbors'
call flush(mpl%info)
net_nnb = 0
do ic1=1,nicas_blk%nc1
   inr = mesh%order_inv(ic1)
   i = mesh%lend(inr)
   init = .true.
   do while ((i/=mesh%lend(inr)).or.init)
      net_nnb(ic1) = net_nnb(ic1)+1
      i = mesh%lptr(i)
      init = .false.
   end do
end do

! Allocation
net_nnbmax = maxval(net_nnb)
allocate(net_inb(net_nnbmax,nicas_blk%nc1))
allocate(net_dnb(net_nnbmax,-1:1,nicas_blk%nc1,nicas_blk%nl1))

! Find mesh neighbors
write(mpl%info,'(a10,a)') '','Find mesh neighbors'
call flush(mpl%info)
net_nnb = 0
do ic1=1,nicas_blk%nc1
   inr = mesh%order_inv(ic1)
   i = mesh%lend(inr)
   init = .true.
   do while ((i/=mesh%lend(inr)).or.init)
      net_nnb(ic1) = net_nnb(ic1)+1
      net_inb(net_nnb(ic1),ic1) = mesh%order(abs(mesh%list(i)))
      i = mesh%lptr(i)
      init = .false.
   end do
end do

! Compute mesh edges distances
write(mpl%info,'(a10,a)',advance='no') '','Compute mesh edges distances: '
call flush(mpl%info)
call mpl%prog_init(nicas_blk%nc1)
net_dnb = 1.0
!$omp parallel do schedule(static) private(ic1,j,ic0,jc1,jc0,dnb,il1,il0,valid_arc,djl1,jl1,jl0,disthsq,distvsq,rhsq,rvsq), &
!$omp&                             private(distnorm_network)
do ic1=1,nicas_blk%nc1
   do j=1,net_nnb(ic1)
      ! Indices
      ic0 = nicas_blk%c1_to_c0(ic1)
      jc1 = net_inb(j,ic1)
      jc0 = nicas_blk%c1_to_c0(jc1)

      if (geom%mask_hor_c0(jc0)) then
         ! True distance
         call sphere_dist(geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),dnb)

         do il1=1,nicas_blk%nl1
            ! Check mask bounds
            il0 = nicas_blk%l1_to_l0(il1)
            valid_arc = .true.
            if (nam%mask_check) call geom%check_arc(il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),valid_arc)

            if (valid_arc) then
               do djl1=-1,1
                  ! Index
                  jl1 = max(1,min(il1+djl1,nicas_blk%nl1))
                  jl0 = nicas_blk%l1_to_l0(jl1)

                  if (geom%mask_c0(ic0,il0).and.geom%mask_c0(jc0,jl0)) then
                     ! Squared support radii
                     disthsq = dnb**2
                     distvsq = (geom%vunit(ic0,il0)-geom%vunit(jc0,jl0))**2
                     rhsq = 0.5*(nicas_blk%rh_c1(ic1,il1)**2+nicas_blk%rh_c1(jc1,jl1)**2)
                     rvsq = 0.5*(nicas_blk%rv_c1(ic1,il1)**2+nicas_blk%rv_c1(jc1,jl1)**2)
                     distnorm_network = 0.0
                     if (rhsq>0.0) then
                        distnorm_network = distnorm_network+disthsq/rhsq
                     elseif (disthsq>0.0) then
                        distnorm_network = distnorm_network+0.5*huge(0.0)
                     end if
                     if (rvsq>0.0) then
                        distnorm_network = distnorm_network+distvsq/rvsq
                     elseif (distvsq>0.0) then
                        distnorm_network = distnorm_network+0.5*huge(0.0)
                     end if
                     net_dnb(j,djl1,ic1,il1) = sqrt(distnorm_network)
                  end if
               end do
            end if
         end do
      end if
   end do

   ! Update
   call mpl%prog_print(ic1)
end do
!$omp end parallel do
write(mpl%info,'(a)') '100%'
call flush(mpl%info)

! Compute distances
write(mpl%info,'(a10,a)',advance='no') '','Compute distances: '
call flush(mpl%info)
call mpl%prog_init(nicas_blk%nsbb)
!$omp parallel do schedule(static) private(isbb,is,ic1,il1,np,np_new,jc1,jl1,k,kc1,dkl1,kl1,disttest,add_to_front,jp), &
!$omp&                             firstprivate(distnorm,plist,plist_new)
do isbb=1,nicas_blk%nsbb
   ! Indices
   is = nicas_blk%sbb_to_s(isbb)
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)

   ! Allocation
   allocate(distnorm(nicas_blk%nc1,nicas_blk%nl1))
   allocate(plist(nicas_blk%nc1*nicas_blk%nl1,2))
   allocate(plist_new(nicas_blk%nc1*nicas_blk%nl1,2))

   ! Initialize the front
   np = 1
   plist(1,1) = ic1
   plist(1,2) = il1
   distnorm = 1.0
   distnorm(ic1,il1) = 0.0

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         jc1 = plist(ip,1)
         jl1 = plist(ip,2)

         ! Loop over neighbors
         do k=1,net_nnb(jc1)
            kc1 = net_inb(k,jc1)
            do dkl1=-1,1
               kl1 = max(1,min(jl1+dkl1,nicas_blk%nl1))
               if (nicas_blk%mask_c1(kc1,kl1)) then
                  disttest = distnorm(jc1,jl1)+net_dnb(k,dkl1,jc1,jl1)
                  if (disttest<1.0) then
                     ! Point is inside the support
                     if (disttest<distnorm(kc1,kl1)) then
                        ! Update distance
                        distnorm(kc1,kl1) = disttest

                        ! Check if the point should be added to the front (avoid duplicates)
                        add_to_front = .true.
                        do jp=1,np_new
                           if ((plist_new(jp,1)==kc1).and.(plist_new(jp,2)==kl1)) then
                              add_to_front = .false.
                              exit
                           end if
                        end do

                        if (add_to_front) then
                           ! Add point to the front
                           np_new = np_new+1
                           plist_new(np_new,1) = kc1
                           plist_new(np_new,2) = kl1
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end do

      ! Copy new front
      np = np_new
      plist(1:np,:) = plist_new(1:np,:)
   end do

   ! Pack data
   do il1=1,nicas_blk%nl1
      do ic1=1,nicas_blk%nc1
         if (supeq(distnorm(ic1,il1),1.0_kind_real)) call msr(distnorm(ic1,il1))
      end do
   end do
   call nicas_blk%distnorm(isbb)%pack(nicas_blk%nc1,nicas_blk%nl1,distnorm)

   ! Release memory
   deallocate(distnorm)
   deallocate(plist)
   deallocate(plist_new)

   ! Update
   call mpl%prog_print(isbb)
end do
write(mpl%info,'(a)') '100%'
call flush(mpl%info)

end subroutine nicas_blk_compute_convol_network

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_distance
!> Purpose: compute convolution with a distance approach
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_distance(nicas_blk,mpl,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk                      !< NICAS data block
type(mpl_type),intent(inout) :: mpl                                   !< MPI data
type(geom_type),intent(in) :: geom                                    !< Geometry

! Local variables
integer :: nnmax,is,ic1,jc1,il1,il0,j,js,ic0,jc0,jl0,jl1
integer :: ic1bb,isbb
integer,allocatable :: nn(:),nn_index(:,:)
real(kind_real) :: disthsq,distvsq,rhsq,rvsq
real(kind_real),allocatable :: distnorm(:,:),nn_dist(:,:)

! Allocation
allocate(nn(nicas_blk%nc1bb))

! Count nearest neighbors
write(mpl%info,'(a10,a)') '','Count nearest neighbors'
call flush(mpl%info)
do ic1bb=1,nicas_blk%nc1bb
   ! Indices
   ic1 = nicas_blk%c1bb_to_c1(ic1bb)
   ic0 = nicas_blk%c1_to_c0(ic1)

   ! Count nearest neighbors
   call nicas_blk%kdtree%count_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nicas_blk%rhmax,nn(ic1bb))
end do
nnmax = maxval(nn)
write(mpl%info,'(a10,a,i6,a,f5.1,a)') '','Number of neighbors to find: ',nnmax,' (', &
& real(nnmax,kind_real)/real(nicas_blk%nc1,kind_real)*100.0,'%)'

! Allocation
allocate(nn_index(nnmax,nicas_blk%nc1bb))
allocate(nn_dist(nnmax,nicas_blk%nc1bb))

! Find nearest neighbors
write(mpl%info,'(a10,a)',advance='no') '','Find nearest neighbors: '
call flush(mpl%info)
call mpl%prog_init(nicas_blk%nc1bb)
do ic1bb=1,nicas_blk%nc1bb
   ! Indices
   ic1 = nicas_blk%c1bb_to_c1(ic1bb)
   ic0 = nicas_blk%c1_to_c0(ic1)

   ! Find nearest neighbors
   call nicas_blk%kdtree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0), &
 & nn(ic1bb),nn_index(1:nn(ic1bb),ic1bb),nn_dist(1:nn(ic1bb),ic1bb))

   ! Update
   call mpl%prog_print(ic1bb)
end do
write(mpl%info,'(a)') '100%'
call flush(mpl%info)

! Compute distances
write(mpl%info,'(a10,a)',advance='no') '','Compute distances: '
call flush(mpl%info)
call mpl%prog_init(nicas_blk%nsbb)
!$omp parallel do schedule(static) private(isbb,is,ic1,ic1bb,ic0,il1,il0,j,jl1,jl0,js,jc1,jc0,disthsq,distvsq,rhsq,rvsq), &
!$omp&                             firstprivate(distnorm)
do isbb=1,nicas_blk%nsbb
   ! Indices
   is = nicas_blk%sbb_to_s(isbb)
   ic1 = nicas_blk%s_to_c1(is)
   ic1bb = nicas_blk%c1_to_c1bb(ic1)
   ic0 = nicas_blk%c1_to_c0(ic1)
   il1 = nicas_blk%s_to_l1(is)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Allocation
   allocate(distnorm(nicas_blk%nc1,nicas_blk%nl1))

   ! Initialization
   call msr(distnorm)

   ! Loop on nearest neighbors
   do j=1,nn(ic1bb)
      jc1 = nn_index(j,ic1bb)
      do jl1=1,nicas_blk%nl1
         if (nicas_blk%mask_c1(jc1,jl1)) then
            jc0 = nicas_blk%c1_to_c0(jc1)
            jl0 = nicas_blk%l1_to_l0(jl1)
            js = nicas_blk%c1l1_to_s(jc1,jl1)

            ! Normalized distance
            disthsq = nn_dist(j,ic1bb)**2
            distvsq = (geom%vunit(ic0,il0)-geom%vunit(jc0,jl0))**2
            rhsq = 0.5*(nicas_blk%rh_c1(ic1,il1)**2+nicas_blk%rh_c1(jc1,jl1)**2)
            rvsq = 0.5*(nicas_blk%rv_c1(ic1,il1)**2+nicas_blk%rv_c1(jc1,jl1)**2)
            if (rhsq>0.0) then
               disthsq = disthsq/rhsq
            elseif (disthsq>0.0) then
               disthsq = 0.5*huge(0.0)
            end if
            if (rvsq>0.0) then
               distvsq = distvsq/rvsq
            elseif (distvsq>0.0) then
               distvsq = 0.5*huge(0.0)
            end if
            distnorm(jc1,jl1) = sqrt(disthsq+distvsq)
         end if
      end do
   end do

   ! Pack data
   do il1=1,nicas_blk%nl1
      do ic1=1,nicas_blk%nc1
         if (supeq(distnorm(ic1,il1),1.0_kind_real)) call msr(distnorm(ic1,il1))
      end do
   end do
   call nicas_blk%distnorm(isbb)%pack(nicas_blk%nc1,nicas_blk%nl1,distnorm)

   ! Release memory
   deallocate(distnorm)

   ! Update
   call mpl%prog_print(isbb)
end do
!$omp end parallel do
write(mpl%info,'(a)') '100%'
call flush(mpl%info)

end subroutine nicas_blk_compute_convol_distance

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_weights
!> Purpose: compute convolution weights
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_weights(nicas_blk,mpl,nam,geom,ctmp)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk      !< NICAS data block
type(mpl_type),intent(inout) :: mpl                   !< MPI data
type(nam_type),intent(in) :: nam                      !< Namelist
type(geom_type),intent(in) :: geom                    !< Geometry
type(linop_type),intent(inout) :: ctmp                !< Convolution operator

! Local variables
integer :: n_s_max,ithread,is,ic1,jc1,il1,jl1,il0,jbd,js,ic0,ic1bb,isbb
integer :: c_n_s(mpl%nthread),c_nor_n_s(mpl%nthread)
real(kind_real) :: disth,S_test
logical :: add_op
type(linop_type) :: c(mpl%nthread),c_nor(mpl%nthread)

! Allocation
n_s_max = 10*nint(real(geom%nc0*geom%nl0)/real(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   call c(ithread)%alloc
   c_nor(ithread)%n_s = n_s_max
   call c_nor(ithread)%alloc
end do

! Initialization
write(mpl%info,'(a10,a)',advance='no') '','Compute weights: '
call flush(mpl%info)
call mpl%prog_init(nicas_blk%nsbb)
c_n_s = 0
c_nor_n_s = 0

! Compute weights
!$omp parallel do schedule(static) private(isbb,is,ithread,ic1,ic1bb,ic0,il1,il0,jbd,jc1,jl1,js,disth,S_test,add_op)
do isbb=1,nicas_blk%nsbb
   ! Indices
   is = nicas_blk%sbb_to_s(isbb)
   ithread = 1
!$ ithread = omp_get_thread_num()+1
   ic1 = nicas_blk%s_to_c1(is)
   ic1bb = nicas_blk%c1_to_c1bb(ic1)
   ic0 = nicas_blk%c1_to_c0(ic1)
   il1 = nicas_blk%s_to_l1(is)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Count convolution operations
   do jbd=1,nicas_blk%distnorm(isbb)%nbd
      ! Indices
      jc1 = nicas_blk%distnorm(isbb)%bd_to_c1(jbd)
      jl1 = nicas_blk%distnorm(isbb)%bd_to_l1(jbd)
      js = nicas_blk%c1l1_to_s(jc1,jl1)

      if (nicas_blk%double_fit) then
         ! Double Gaspari-Cohn (1999) function
         disth = sqrt(nicas_blk%distnorm(isbb)%val(jbd)**2-nicas_blk%distnormv(isbb)%val(jbd)**2)
         S_test = gc99(mpl,disth)*((1.0+nicas_blk%coef(isbb)%val(jbd))*gc99(mpl,nicas_blk%distnormv(isbb)%val(jbd) &
                & /nicas_blk%rfac(isbb)%val(jbd))-nicas_blk%coef(isbb)%val(jbd)*gc99(mpl,nicas_blk%distnormv(isbb)%val(jbd)))
      else
         ! Gaspari-Cohn (1999) function
         S_test = gc99(mpl,nicas_blk%distnorm(isbb)%val(jbd))
      end if

      if (abs(S_test)>abs(S_inf)) then
         ! Store coefficient for convolution
         if (nam%lsqrt) then
            add_op = .false.
            if (nam%mpicom==1) then
               add_op = (nicas_blk%lcheck_sb(js).and.(is<=js)).or.(.not.nicas_blk%lcheck_sb(js))
            elseif (nam%mpicom==2) then
               add_op = nicas_blk%lcheck_sa(is).and.((nicas_blk%lcheck_sa(js).and.(is<=js)) &
                     & .or.(.not.nicas_blk%lcheck_sa(js)))
            end if
         else
            add_op = .true.
         end if
         if (add_op) call c(ithread)%add_op(c_n_s(ithread),is,js,S_test)

         ! Store coefficient for normalization
         add_op = nicas_blk%lcheck_sb(is).and.((nicas_blk%lcheck_sb(js).and.(is<=js)).or.(.not.nicas_blk%lcheck_sb(js)))
         if (add_op) call c_nor(ithread)%add_op(c_nor_n_s(ithread),is,js,S_test)
      end if
   end do

   ! Update
   call mpl%prog_print(isbb)
end do
!$omp end parallel do
write(mpl%info,'(a)') '100%'
call flush(mpl%info)

! Gather data from threads
call ctmp%gather(mpl,c_n_s,c)
call nicas_blk%c_nor%gather(mpl,c_nor_n_s,c_nor)

! Release memory
do ithread=1,mpl%nthread
   call c(ithread)%dealloc
   call c_nor(ithread)%dealloc
end do

end subroutine nicas_blk_compute_convol_weights

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_c
!> Purpose: compute NICAS MPI distribution, halo C
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_c(nicas_blk,mpl,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: isa,isb,isc,i_s,is,js
integer,allocatable :: s_to_proc(:)
logical,allocatable :: lcheck_sc_nor(:)

! Allocation
allocate(nicas_blk%lcheck_sc(nicas_blk%ns))
allocate(lcheck_sc_nor(nicas_blk%ns))
allocate(s_to_proc(nicas_blk%ns))

! Halo definitions

! Halo C
nicas_blk%lcheck_sc = nicas_blk%lcheck_sb
do i_s=1,nicas_blk%c%n_s
   is = nicas_blk%c%row(i_s)
   js = nicas_blk%c%col(i_s)
   nicas_blk%lcheck_sc(is) = .true.
   nicas_blk%lcheck_sc(js) = .true.
end do
nicas_blk%nsc = count(nicas_blk%lcheck_sc)
lcheck_sc_nor = nicas_blk%lcheck_sb
do i_s=1,nicas_blk%c_nor%n_s
   is = nicas_blk%c_nor%row(i_s)
   js = nicas_blk%c_nor%col(i_s)
   lcheck_sc_nor(is) = .true.
   lcheck_sc_nor(js) = .true.
end do
nicas_blk%nsc_nor = count(lcheck_sc_nor)

! Check halos consistency
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is).and.(.not.nicas_blk%lcheck_sc(is))) call mpl%abort('point in halo A but not in halo C')
   if (nicas_blk%lcheck_sb(is).and.(.not.nicas_blk%lcheck_sc(is))) call mpl%abort('point in halo B but not in halo C')
end do

! Global <-> local conversions for fields

! Halo C
allocate(nicas_blk%sc_to_s(nicas_blk%nsc))
allocate(nicas_blk%s_to_sc(nicas_blk%ns))
call msi(nicas_blk%s_to_sc)
isc = 0
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sc(is)) then
      isc = isc+1
      nicas_blk%sc_to_s(isc) = is
      nicas_blk%s_to_sc(is) = isc
   end if
end do

allocate(nicas_blk%sc_nor_to_s(nicas_blk%nsc_nor))
allocate(nicas_blk%s_to_sc_nor(nicas_blk%ns))
call msi(nicas_blk%s_to_sc_nor)
isc = 0
do is=1,nicas_blk%ns
   if (lcheck_sc_nor(is)) then
      isc = isc+1
      nicas_blk%sc_nor_to_s(isc) = is
      nicas_blk%s_to_sc_nor(is) = isc
   end if
end do

! Inter-halo conversions
allocate(nicas_blk%sa_to_sc(nicas_blk%nsa))
do isa=1,nicas_blk%nsa
   is = nicas_blk%sa_to_s(isa)
   isc = nicas_blk%s_to_sc(is)
   nicas_blk%sa_to_sc(isa) = isc
end do
allocate(nicas_blk%sb_to_sc(nicas_blk%nsb))
allocate(nicas_blk%sc_to_sb(nicas_blk%nsc))
allocate(nicas_blk%sb_to_sc_nor(nicas_blk%nsb))
call msi(nicas_blk%sc_to_sb)
do isb=1,nicas_blk%nsb
   is = nicas_blk%sb_to_s(isb)
   isc = nicas_blk%s_to_sc(is)
   nicas_blk%sb_to_sc(isb) = isc
   nicas_blk%sc_to_sb(isc) = isb
   isc = nicas_blk%s_to_sc_nor(is)
   nicas_blk%sb_to_sc_nor(isb) = isc
end do

! Local data

! Convolution
nicas_blk%c%n_src = nicas_blk%nsc
nicas_blk%c%n_dst = nicas_blk%nsc
do i_s=1,nicas_blk%c%n_s
   nicas_blk%c%row(i_s) = nicas_blk%s_to_sc(nicas_blk%c%row(i_s))
   nicas_blk%c%col(i_s) = nicas_blk%s_to_sc(nicas_blk%c%col(i_s))
end do
call nicas_blk%c%reorder(mpl)

! Convolution for normalization
nicas_blk%c_nor%n_src = nicas_blk%nsc
nicas_blk%c_nor%n_dst = nicas_blk%nsc
do i_s=1,nicas_blk%c_nor%n_s
   nicas_blk%c_nor%row(i_s) = nicas_blk%s_to_sc_nor(nicas_blk%c_nor%row(i_s))
   nicas_blk%c_nor%col(i_s) = nicas_blk%s_to_sc_nor(nicas_blk%c_nor%col(i_s))
end do
call nicas_blk%c_nor%reorder(mpl)

! Setup communications
s_to_proc = geom%c0_to_proc(nicas_blk%c1_to_c0(nicas_blk%s_to_c1))
call nicas_blk%com_AC%setup(mpl,'com_AC',nicas_blk%ns,nicas_blk%nsa,nicas_blk%nsc,nicas_blk%sc_to_s,nicas_blk%sa_to_sc, &
 & s_to_proc,nicas_blk%s_to_sa)

end subroutine nicas_blk_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_normalization
!> Purpose: compute normalization
!----------------------------------------------------------------------
subroutine nicas_blk_compute_normalization(nicas_blk,mpl,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il0i,i_s,ic1,ic1b,jc1b,is,js,isc,jsb,jsc,ic0,ic0a,il0,il1,ih,jv,nlr,ilr,ic,isc_add
integer,allocatable :: ineh(:,:),inev(:),ines(:,:),inec(:),order(:),isc_list(:)
integer,allocatable :: h_col(:,:,:),v_col(:,:),s_col(:,:,:),c_ind(:,:)
real(kind_real) :: S_add
real(kind_real),allocatable :: h_S(:,:,:),v_S(:,:,:),s_S(:,:,:),c_S(:,:)
real(kind_real),allocatable :: list(:),S_list(:),S_list_tmp(:)

! Compute horizontal interpolation inverse mapping
allocate(ineh(geom%nc0a,geom%nl0i))
ineh = 0
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%h(il0i)%n_s
      ic0a = nicas_blk%h(il0i)%row(i_s)
      ineh(ic0a,il0i) = ineh(ic0a,il0i)+1
   end do
end do
allocate(h_col(maxval(ineh),geom%nc0a,geom%nl0i))
allocate(h_S(maxval(ineh),geom%nc0a,geom%nl0i))
call msi(h_col)
call msr(h_S)
ineh = 0
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%h(il0i)%n_s
      ic0a = nicas_blk%h(il0i)%row(i_s)
      ineh(ic0a,il0i) = ineh(ic0a,il0i)+1
      h_col(ineh(ic0a,il0i),ic0a,il0i) = nicas_blk%h(il0i)%col(i_s)
      h_S(ineh(ic0a,il0i),ic0a,il0i) = nicas_blk%h(il0i)%S(i_s)
   end do
end do

! Compute vertical interpolation inverse mapping
allocate(inev(geom%nl0))
inev = 0
do i_s=1,nicas_blk%v%n_s
   il0 = nicas_blk%v%row(i_s)
   inev(il0) = inev(il0)+1
end do
allocate(v_col(maxval(inev),geom%nl0))
allocate(v_S(maxval(inev),nicas_blk%nc1b,geom%nl0))
call msi(v_col)
call msr(v_S)
inev = 0
do i_s=1,nicas_blk%v%n_s
   il0 = nicas_blk%v%row(i_s)
   inev(il0) = inev(il0)+1
   v_col(inev(il0),il0) = nicas_blk%v%col(i_s)
   do ic1b=1,nicas_blk%nc1b
      v_S(inev(il0),ic1b,il0) = nicas_blk%v%Svec(i_s,ic1b)
   end do
end do

! Compute subsampling interpolation inverse mapping
allocate(ines(nicas_blk%nc1b,nicas_blk%nl1))
ines = 0
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%s(il1)%n_s
      ic1b = nicas_blk%s(il1)%row(i_s)
      ines(ic1b,il1) = ines(ic1b,il1)+1
   end do
end do
allocate(s_col(maxval(ines),nicas_blk%nc1b,nicas_blk%nl1))
allocate(s_S(maxval(ines),nicas_blk%nc1b,nicas_blk%nl1))
call msi(s_col)
call msr(s_S)
ines = 0
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%s(il1)%n_s
      ic1b = nicas_blk%s(il1)%row(i_s)
      ines(ic1b,il1) = ines(ic1b,il1)+1
      jc1b = nicas_blk%s(il1)%col(i_s)
      jsb = nicas_blk%c1bl1_to_sb(jc1b,il1)
      jsc = nicas_blk%sb_to_sc_nor(jsb)
      s_col(ines(ic1b,il1),ic1b,il1) = jsc
      s_S(ines(ic1b,il1),ic1b,il1) = nicas_blk%s(il1)%S(i_s)
   end do
end do

! Compute convolution inverse mapping
allocate(inec(nicas_blk%nsc_nor))
inec = 0
do i_s=1,nicas_blk%c_nor%n_s
   isc = nicas_blk%c_nor%col(i_s)
   jsc = nicas_blk%c_nor%row(i_s)
   if (jsc/=isc) then
      is = nicas_blk%sc_nor_to_s(isc)
      inec(isc) = inec(isc)+1
      js = nicas_blk%sc_nor_to_s(jsc)
      inec(jsc) = inec(jsc)+1
   end if
end do
allocate(c_ind(maxval(inec),nicas_blk%nsc_nor))
allocate(c_S(maxval(inec),nicas_blk%nsc_nor))
call msi(c_ind)
call msr(c_S)
inec = 0
do i_s=1,nicas_blk%c_nor%n_s
   isc = nicas_blk%c_nor%col(i_s)
   jsc = nicas_blk%c_nor%row(i_s)
   if (jsc/=isc) then
      is = nicas_blk%sc_nor_to_s(isc)
      inec(isc) = inec(isc)+1
      c_ind(inec(isc),isc) = jsc
      c_S(inec(isc),isc) = nicas_blk%c_nor%S(i_s)
      js = nicas_blk%sc_nor_to_s(jsc)
      inec(jsc) = inec(jsc)+1
      c_ind(inec(jsc),jsc) = isc
      c_S(inec(jsc),jsc) = nicas_blk%c_nor%S(i_s)
   end if
end do

! Re-order indices
do isc=1,nicas_blk%nsc_nor
   ! Allocation
   allocate(order(inec(isc)))
   allocate(list(inec(isc)))

   ! Copy
   list = c_ind(1:inec(isc),isc)

   ! Order
   call qsort(inec(isc),list,order)

   ! Re-order
   c_ind(1:inec(isc),isc) = c_ind(order(1:inec(isc)),isc)
   c_S(1:inec(isc),isc) = c_S(order(1:inec(isc)),isc)

   ! Release memory
   deallocate(order)
   deallocate(list)
end do

! Allocation
allocate(nicas_blk%norm(geom%nc0a,geom%nl0))
call msr(nicas_blk%norm)

! Compute normalization weights
do il0=1,geom%nl0
   il0i = min(il0,geom%nl0i)
   write(mpl%info,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),': '
   call flush(mpl%info)
   call mpl%prog_init(geom%nc0a)

   !$omp parallel do schedule(static) private(ic0a,ic0,nlr,isc_add,S_add,ih,ic1b,ic1,jv,il1,is,ilr,ic,isc,jsc), &
   !$omp&                             firstprivate(isc_list,S_list,S_list_tmp)
   do ic0a=1,geom%nc0a
      ! Index
      ic0 = geom%c0a_to_c0(ic0a)

      if (geom%mask_c0(ic0,il0)) then
         ! Allocation
         allocate(isc_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         allocate(S_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         allocate(S_list_tmp(nicas_blk%nsc_nor))

         ! Initialization
         isc_list = 0
         S_list = 0.0

         ! Adjoint interpolation
         nlr = 0
         do ih=1,ineh(ic0a,il0i)
            ic1b = h_col(ih,ic0a,il0i)
            ic1 = nicas_blk%c1b_to_c1(ic1b)
            do jv=1,inev(il0)
               il1 = v_col(jv,il0)
               if ((nicas_blk%l1_to_l0(il1)>=nicas_blk%vbot(ic1)).and.(nicas_blk%l1_to_l0(il1)<=nicas_blk%vtop(ic1))) then
                  do is=1,ines(ic1b,il1)
                     isc_add = s_col(is,ic1b,il1)
                     S_add = h_S(ih,ic0a,il0i)*v_S(jv,ic1b,il0)*s_S(is,ic1b,il1)
                     if (nlr==0) then
                        ilr = 1
                        nlr = 1
                     else
                        do ilr=1,nlr
                           if (isc_add==isc_list(ilr)) exit
                        end do
                        if (ilr==nlr+1) nlr = nlr+1
                     end if
                     isc_list(ilr) = isc_add
                     S_list(ilr) = S_list(ilr)+S_add
                  end do
               end if
            end do
         end do

         ! Initialization
         S_list_tmp = 0.0
         do ilr=1,nlr
            isc = isc_list(ilr)
            S_list_tmp(isc) = S_list(ilr)
         end do

         ! Convolution
         do ilr=1,nlr
            isc = isc_list(ilr)
            do ic=1,inec(isc)
               jsc = c_ind(ic,isc)
               S_list_tmp(jsc) = S_list_tmp(jsc)+c_S(ic,isc)*S_list(ilr)
            end do
         end do

         ! Sum of squared values
         nicas_blk%norm(ic0a,il0) = sum(S_list_tmp**2)

         ! Normalization factor
         nicas_blk%norm(ic0a,il0) = 1.0/sqrt(nicas_blk%norm(ic0a,il0))

         ! Update
         call mpl%prog_print(ic0a)

         ! Release memory
         deallocate(isc_list)
         deallocate(S_list)
         deallocate(S_list_tmp)
      end if
   end do
   !$omp end parallel do
   write(mpl%info,'(a)') '100%'
   call flush(mpl%info)
end do

end subroutine nicas_blk_compute_normalization

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_adv
!> Purpose: compute advection
!----------------------------------------------------------------------
subroutine nicas_blk_compute_adv(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng              !< Random number generator
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: its,il0,ic0,ic0a,i_s,i_s_loc,iproc,jc0
integer :: ic0d,d_n_s_max,d_n_s_max_loc
integer :: ic0dinv,dinv_n_s_max,dinv_n_s_max_loc
integer,allocatable :: c0d_to_c0(:),c0_to_c0d(:),c0a_to_c0d(:),interpd_lg(:,:,:)
integer,allocatable :: c0dinv_to_c0(:),c0_to_c0dinv(:),c0a_to_c0dinv(:),interpdinv_lg(:,:,:)
real(kind_real) :: displ_lon(geom%nc0,geom%nl0),displ_lat(geom%nc0,geom%nl0)
logical :: mask_c0(geom%nc0)
logical,allocatable :: lcheck_c0d(:),lcheck_d(:,:,:)
logical,allocatable :: lcheck_c0dinv(:),lcheck_dinv(:,:,:)
type(linop_type),allocatable :: dfull(:,:)
type(linop_type),allocatable :: dinvfull(:,:)

write(mpl%info,'(a7,a)') '','Compute advection'
call flush(mpl%info)

! Allocation
allocate(dfull(geom%nl0,2:nam%nts))
allocate(dinvfull(geom%nl0,2:nam%nts))

! Initialization
mask_c0 = .true.

do its=2,nam%nts
   ! Local to global
   call mpl%loc_to_glb(geom%nl0,geom%nc0a,cmat_blk%displ_lon(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,displ_lon)
   call mpl%loc_to_glb(geom%nl0,geom%nc0a,cmat_blk%displ_lat(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,displ_lat)

   do il0=1,geom%nl0
      ! Compute direct interpolation
      call dfull(il0,its)%interp(mpl,rng,geom%nc0,displ_lon(:,il0),displ_lat(:,il0),mask_c0,geom%nc0,geom%lon,geom%lat,mask_c0, &
    & nam%diag_interp)

      ! Compute inverse interpolation
      call dinvfull(il0,its)%interp(mpl,rng,geom%nc0,geom%lon,geom%lat,mask_c0,geom%nc0,displ_lon(:,il0),displ_lat(:,il0),mask_c0, &
    & nam%diag_interp)
   end do
end do

! Allocation
d_n_s_max = 0
dinv_n_s_max = 0
do its=2,nam%nts
   do il0=1,geom%nl0
      d_n_s_max = max(d_n_s_max,dfull(il0,its)%n_s)
      dinv_n_s_max = max(dinv_n_s_max,dinvfull(il0,its)%n_s)
   end do
end do
allocate(lcheck_c0d(geom%nc0))
allocate(lcheck_c0dinv(geom%nc0))
allocate(lcheck_d(d_n_s_max,geom%nl0,2:nam%nts))
allocate(lcheck_dinv(dinv_n_s_max,geom%nl0,2:nam%nts))
allocate(nicas_blk%d(geom%nl0,2:nam%nts))
allocate(nicas_blk%dinv(geom%nl0,2:nam%nts))

! Halo definitions

! Halo A
lcheck_c0d = .false.
lcheck_c0dinv = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   lcheck_c0d(ic0) = .true.
   lcheck_c0dinv(ic0) = .true.
end do

! Halo D
lcheck_d = .false.
lcheck_dinv = .false.
do its=2,nam%nts
   do il0=1,geom%nl0
      do i_s=1,dfull(il0,its)%n_s
         ic0 = dfull(il0,its)%row(i_s)
         iproc = geom%c0_to_proc(ic0)
         if (iproc==mpl%myproc) then
            jc0 = dfull(il0,its)%col(i_s)
            lcheck_d(i_s,il0,its) = .true.
            lcheck_c0d(jc0) = .true.
         end if
      end do
      do i_s=1,dinvfull(il0,its)%n_s
         ic0 = dinvfull(il0,its)%row(i_s)
         iproc = geom%c0_to_proc(ic0)
         if (iproc==mpl%myproc) then
            jc0 = dinvfull(il0,its)%col(i_s)
            lcheck_dinv(i_s,il0,its) = .true.
            lcheck_c0dinv(jc0) = .true.
         end if
      end do
   end do
end do

! Halo sizes
do its=2,nam%nts
   do il0=1,geom%nl0
      nicas_blk%d(il0,its)%n_s = count(lcheck_d(:,il0,its))
      nicas_blk%dinv(il0,its)%n_s = count(lcheck_dinv(:,il0,its))
   end do
end do
nicas_blk%nc0d = count(lcheck_c0d)
nicas_blk%nc0dinv = count(lcheck_c0dinv)

! Global <-> local conversions for fields
allocate(c0d_to_c0(nicas_blk%nc0d))
allocate(c0dinv_to_c0(nicas_blk%nc0dinv))
allocate(c0_to_c0d(geom%nc0))
allocate(c0_to_c0dinv(geom%nc0))
call msi(c0_to_c0d)
call msi(c0_to_c0dinv)
ic0d = 0
ic0dinv = 0
do ic0=1,geom%nc0
   if (lcheck_c0d(ic0)) then
      ic0d = ic0d+1
      c0d_to_c0(ic0d) = ic0
      c0_to_c0d(ic0) = ic0d
   end if
   if (lcheck_c0dinv(ic0)) then
      ic0dinv = ic0dinv+1
      c0dinv_to_c0(ic0dinv) = ic0
      c0_to_c0dinv(ic0) = ic0dinv
   end if
end do

! Inter-halo conversions
allocate(c0a_to_c0d(geom%nc0a))
allocate(c0a_to_c0dinv(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0d = c0_to_c0d(ic0)
   ic0dinv = c0_to_c0dinv(ic0)
   c0a_to_c0d(ic0a) = ic0d
   c0a_to_c0dinv(ic0a) = ic0dinv
end do

! Global <-> local conversions for data
d_n_s_max_loc = 0
dinv_n_s_max_loc = 0
do its=2,nam%nts
   do il0=1,geom%nl0
      d_n_s_max_loc = max(d_n_s_max_loc,nicas_blk%d(il0,its)%n_s)
      dinv_n_s_max_loc = max(dinv_n_s_max_loc,nicas_blk%dinv(il0,its)%n_s)
   end do
end do
allocate(interpd_lg(d_n_s_max_loc,geom%nl0,2:nam%nts))
allocate(interpdinv_lg(dinv_n_s_max_loc,geom%nl0,2:nam%nts))
do its=2,nam%nts
   do il0=1,geom%nl0
      i_s_loc = 0
      do i_s=1,dfull(il0,its)%n_s
         if (lcheck_d(i_s,il0,its)) then
            i_s_loc = i_s_loc+1
            interpd_lg(i_s_loc,il0,its) = i_s
         end if
      end do
      i_s_loc = 0
      do i_s=1,dinvfull(il0,its)%n_s
         if (lcheck_dinv(i_s,il0,its)) then
            i_s_loc = i_s_loc+1
            interpdinv_lg(i_s_loc,il0,its) = i_s
         end if
      end do
   end do
end do

! Local data
do its=2,nam%nts
   do il0=1,geom%nl0
      write(nicas_blk%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
      write(nicas_blk%dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'dinv_',il0,'_',its
      nicas_blk%d(il0,its)%n_src = nicas_blk%nc0d
      nicas_blk%dinv(il0,its)%n_src = nicas_blk%nc0dinv
      nicas_blk%d(il0,its)%n_dst = geom%nc0a
      nicas_blk%dinv(il0,its)%n_dst = geom%nc0a
      call nicas_blk%d(il0,its)%alloc
      call nicas_blk%dinv(il0,its)%alloc
      do i_s_loc=1,nicas_blk%d(il0,its)%n_s
         i_s = interpd_lg(i_s_loc,il0,its)
         nicas_blk%d(il0,its)%row(i_s_loc) = geom%c0_to_c0a(dfull(il0,its)%row(i_s))
         nicas_blk%d(il0,its)%col(i_s_loc) = c0_to_c0d(dfull(il0,its)%col(i_s))
         nicas_blk%d(il0,its)%S(i_s_loc) = dfull(il0,its)%S(i_s)
      end do
      do i_s_loc=1,nicas_blk%dinv(il0,its)%n_s
         i_s = interpdinv_lg(i_s_loc,il0,its)
         nicas_blk%dinv(il0,its)%row(i_s_loc) = geom%c0_to_c0a(dinvfull(il0,its)%row(i_s))
         nicas_blk%dinv(il0,its)%col(i_s_loc) = c0_to_c0dinv(dinvfull(il0,its)%col(i_s))
         nicas_blk%dinv(il0,its)%S(i_s_loc) = dinvfull(il0,its)%S(i_s)
      end do
      call nicas_blk%d(il0,its)%reorder(mpl)
      call nicas_blk%dinv(il0,its)%reorder(mpl)
   end do
end do

! Setup communications
call nicas_blk%com_AD%setup(mpl,'com_AD',geom%nc0,geom%nc0a,nicas_blk%nc0d,c0d_to_c0,c0a_to_c0d,geom%c0_to_proc,geom%c0_to_c0a)
call nicas_blk%com_ADinv%setup(mpl,'com_ADinv',geom%nc0,geom%nc0a,nicas_blk%nc0dinv,c0dinv_to_c0,c0a_to_c0dinv, &
 & geom%c0_to_proc,geom%c0_to_c0a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%info,'(a10,a,i8)') '','nc0d =       ',nicas_blk%nc0d
write(mpl%info,'(a10,a,i8)') '','nc0dinv =    ',nicas_blk%nc0dinv
do its=2,nam%nts
   do il0=1,geom%nl0
      write(mpl%info,'(a10,a,i3,a,i2,a,i8)') '','d(',il0,',',its,')%n_s =    ',nicas_blk%d(il0,its)%n_s
      write(mpl%info,'(a10,a,i3,a,i2,a,i8)') '','dinv(',il0,',',its,')%n_s = ',nicas_blk%dinv(il0,its)%n_s
   end do
end do
call flush(mpl%info)

end subroutine nicas_blk_compute_adv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply
!> Purpose: apply NICAS method
!----------------------------------------------------------------------
subroutine nicas_blk_apply(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk            !< NICAS data block
type(mpl_type),intent(in) :: mpl                         !< MPI data
type(nam_type),intent(in) :: nam                         !< Namelist
type(geom_type),intent(in) :: geom                       !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: alpha_a(nicas_blk%nsa),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

! Adjoint interpolation
call nicas_blk%apply_interp_ad(mpl,geom,fld,alpha_b)

! Communication
if (nam%mpicom==1) then
   ! Initialization
   alpha_c = 0.0

   ! Copy zone B into zone C
   alpha_c(nicas_blk%sb_to_sc) = alpha_b
elseif (nam%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call nicas_blk%com_AB%red(mpl,alpha_b,alpha_a)

   ! Initialization
   alpha_c = 0.0

   ! Copy zone A into zone C
   alpha_c(nicas_blk%sa_to_sc) = alpha_a
end if

! Convolution
call nicas_blk%apply_convol(mpl,alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(mpl,alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(mpl,geom,alpha_b,fld)

end subroutine nicas_blk_apply

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_from_sqrt
!> Purpose: apply NICAS method from its square-root formulation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_from_sqrt(nicas_blk,mpl,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk            !< NICAS data block
type(mpl_type),intent(in) :: mpl                         !< MPI data
type(geom_type),intent(in) :: geom                       !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: alpha(nicas_blk%nsa)

! Apply square-root adjoint
call nicas_blk%apply_sqrt_ad(mpl,geom,fld,alpha)

! Apply square-root
call nicas_blk%apply_sqrt(mpl,geom,alpha,fld)

end subroutine nicas_blk_apply_from_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_sqrt
!> Purpose: apply NICAS method square-root
!----------------------------------------------------------------------
subroutine nicas_blk_apply_sqrt(nicas_blk,mpl,geom,alpha,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk          !< NICAS data block
type(mpl_type),intent(in) :: mpl                       !< MPI data
type(geom_type),intent(in) :: geom                     !< Geometry
real(kind_real),intent(in) :: alpha(nicas_blk%nsa)     !< Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variable
real(kind_real) :: alpha_a(nicas_blk%nsa),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

! Initialization
alpha_c = 0.0

! Copy zone A into zone C
alpha_c(nicas_blk%sa_to_sc) = alpha

! Convolution
call nicas_blk%apply_convol(mpl,alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(mpl,alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(mpl,geom,alpha_b,fld)

end subroutine nicas_blk_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_sqrt_ad
!> Purpose: apply NICAS method square-root adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_sqrt_ad(nicas_blk,mpl,geom,fld,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         !< NICAS data block
type(mpl_type),intent(in) :: mpl                      !< MPI data
type(geom_type),intent(in) :: geom                    !< Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field
real(kind_real),intent(out) :: alpha(nicas_blk%nsa)   !< Subgrid field

! Local variable
real(kind_real) :: alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

! Adjoint interpolation
call nicas_blk%apply_interp_ad(mpl,geom,fld,alpha_b)

! Halo reduction from zone B to zone A
call nicas_blk%com_AB%red(mpl,alpha_b,alpha)

! Initialization
alpha_c = 0.0

! Copy zone A into zone C
alpha_c(nicas_blk%sa_to_sc) = alpha

! Convolution
call nicas_blk%apply_convol(mpl,alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha)

end subroutine nicas_blk_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp
!> Purpose: apply interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp(nicas_blk,mpl,geom,alpha,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk          !< NICAS data block
type(mpl_type),intent(in) :: mpl                       !< MPI data
type(geom_type),intent(in) :: geom                     !< Geometry
real(kind_real),intent(in) :: alpha(nicas_blk%nsb)     !< Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: gamma(nicas_blk%nc1b,nicas_blk%nl1),delta(nicas_blk%nc1b,geom%nl0)

! Subsampling horizontal interpolation
call nicas_blk%apply_interp_s(mpl,alpha,gamma)

! Vertical interpolation
call nicas_blk%apply_interp_v(mpl,geom,gamma,delta)

! Horizontal interpolation
call nicas_blk%apply_interp_h(mpl,geom,delta,fld)

! Normalization
fld = fld*nicas_blk%norm

end subroutine nicas_blk_apply_interp

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_ad
!> Purpose: apply interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_ad(nicas_blk,mpl,geom,fld,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         !< NICAS data block
type(mpl_type),intent(in) :: mpl                      !< MPI data
type(geom_type),intent(in) :: geom                    !< Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field
real(kind_real),intent(out) :: alpha(nicas_blk%nsb)   !< Subgrid field

! Local variables
real(kind_real) :: gamma(nicas_blk%nc1b,nicas_blk%nl1),delta(nicas_blk%nc1b,geom%nl0)
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0)

! Normalization
fld_tmp = fld*nicas_blk%norm

! Horizontal interpolation
call nicas_blk%apply_interp_h_ad(mpl,geom,fld_tmp,delta)

! Vertical interpolation
call nicas_blk%apply_interp_v_ad(mpl,geom,delta,gamma)

! Subsampling horizontal interpolation
call nicas_blk%apply_interp_s_ad(mpl,gamma,alpha)

end subroutine nicas_blk_apply_interp_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_h
!> Purpose: apply horizontal interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_h(nicas_blk,mpl,geom,delta,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                !< NICAS data block
type(mpl_type),intent(in) :: mpl                             !< MPI data
type(geom_type),intent(in) :: geom                           !< Geometry
real(kind_real),intent(in) :: delta(nicas_blk%nc1b,geom%nl0) !< Subset Sc1 field, full levels
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0)       !< Field

! Local variables
integer :: il0

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call nicas_blk%h(min(il0,geom%nl0i))%apply(mpl,delta(:,il0),fld(:,il0),msdst=.false.)
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_h

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_h_ad
!> Purpose: apply horizontal interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_h_ad(nicas_blk,mpl,geom,fld,delta)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                 !< NICAS data block
type(mpl_type),intent(in) :: mpl                              !< MPI data
type(geom_type),intent(in) :: geom                            !< Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0)         !< Field
real(kind_real),intent(out) :: delta(nicas_blk%nc1b,geom%nl0) !< Subset Sc1 field, full levels

! Local variables
integer :: il0

!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call nicas_blk%h(min(il0,geom%nl0i))%apply_ad(mpl,fld(:,il0),delta(:,il0))
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_h_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_v
!> Purpose: apply vertical interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_v(nicas_blk,mpl,geom,gamma,delta)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                     !< NICAS data block
type(mpl_type),intent(in) :: mpl                                  !< MPI data
type(geom_type),intent(in) :: geom                                !< Geometry
real(kind_real),intent(in) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) !< Subset Sc1 field, limited levels
real(kind_real),intent(out) :: delta(nicas_blk%nc1b,geom%nl0)     !< Subset Sc1 field, full levels

! Local variables
integer :: ic1b
real(kind_real),allocatable :: gamma_tmp(:),delta_tmp(:)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b) firstprivate(gamma_tmp,delta_tmp)
do ic1b=1,nicas_blk%nc1b
   ! Allocation
   allocate(gamma_tmp(nicas_blk%nl1))
   allocate(delta_tmp(geom%nl0))

   ! Copy data
   gamma_tmp = gamma(ic1b,:)

   ! Apply interpolation
   call nicas_blk%v%apply(mpl,gamma_tmp,delta_tmp,ivec=ic1b,msdst=.false.)

   ! Copy data
   delta(ic1b,:) = delta_tmp

   ! Release memory
   deallocate(gamma_tmp)
   deallocate(delta_tmp)
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_v

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_v_ad
!> Purpose: apply vertical interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_v_ad(nicas_blk,mpl,geom,delta,gamma)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                      !< NICAS data block
type(mpl_type),intent(in) :: mpl                                   !< MPI data
type(geom_type),intent(in) :: geom                                 !< Geometry
real(kind_real),intent(in) :: delta(nicas_blk%nc1b,geom%nl0)       !< Subset Sc1 field, full levels
real(kind_real),intent(out) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) !< Subset Sc1 field, limited levels

! Local variables
integer :: ic1b
real(kind_real),allocatable :: gamma_tmp(:),delta_tmp(:)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b) firstprivate(gamma_tmp,delta_tmp)
do ic1b=1,nicas_blk%nc1b
   ! Allocation
   allocate(gamma_tmp(nicas_blk%nl1))
   allocate(delta_tmp(geom%nl0))

   ! Copy data
   delta_tmp = delta(ic1b,:)

   ! Apply interpolation
   call nicas_blk%v%apply_ad(mpl,delta_tmp,gamma_tmp,ivec=ic1b)

   ! Copy data
   gamma(ic1b,:) = gamma_tmp

   ! Release memory
   deallocate(gamma_tmp)
   deallocate(delta_tmp)
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_v_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_s
!> Purpose: apply subsampling interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_s(nicas_blk,mpl,alpha,gamma)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                      !< NICAS data block
type(mpl_type),intent(in) :: mpl                                   !< MPI data
real(kind_real),intent(in) :: alpha(nicas_blk%nsb)                 !< Subgrid field
real(kind_real),intent(out) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) !< Subset Sc1 field, limited levels

! Local variables
integer :: isb,il1
real(kind_real) :: beta(nicas_blk%nc1b,nicas_blk%nl1)

! Initialization
beta = 0.0

! Copy
!$omp parallel do schedule(static) private(isb)
do isb=1,nicas_blk%nsb
   beta(nicas_blk%sb_to_c1b(isb),nicas_blk%sb_to_l1(isb)) = alpha(isb)
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,nicas_blk%nl1
   call nicas_blk%s(il1)%apply(mpl,beta(:,il1),gamma(:,il1),msdst=.false.)
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_s

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_s_ad
!> Purpose: apply subsampling interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_s_ad(nicas_blk,mpl,gamma,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                     !< NICAS data block
type(mpl_type),intent(in) :: mpl                                  !< MPI data
real(kind_real),intent(in) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) !< Subset Sc1 field, limited levels
real(kind_real),intent(out) :: alpha(nicas_blk%nsb)               !< Subgrid field

! Local variables
integer :: il1,isb
real(kind_real) :: beta(nicas_blk%nc1b,nicas_blk%nl1)

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,nicas_blk%nl1
   call nicas_blk%s(il1)%apply_ad(mpl,gamma(:,il1),beta(:,il1))
end do
!$omp end parallel do

! Copy
!$omp parallel do schedule(static) private(isb)
do isb=1,nicas_blk%nsb
   alpha(isb) = beta(nicas_blk%sb_to_c1b(isb),nicas_blk%sb_to_l1(isb))
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_s_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_convol
!> Purpose: apply convolution
!----------------------------------------------------------------------
subroutine nicas_blk_apply_convol(nicas_blk,mpl,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         !< NICAS data block
type(mpl_type),intent(in) :: mpl                      !< MPI data
real(kind_real),intent(inout) :: alpha(nicas_blk%nsc) !< Subgrid field

! Apply linear operator, symmetric
call nicas_blk%c%apply_sym(mpl,alpha)

end subroutine nicas_blk_apply_convol

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv
!> Purpose: apply advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           !< NICAS data block
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_d(nicas_blk%nc0d,geom%nl0)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Halo extension from zone A to zone D
      call nicas_blk%com_AD%ext(mpl,geom%nl0,fld(:,:,iv,its),fld_d)

      ! Interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call nicas_blk%d(il0,its)%apply(mpl,fld_d(:,il0),fld(:,il0,iv,its),msdst=.false.)
      end do
      !$omp end parallel do
   end do
end do

end subroutine nicas_blk_apply_adv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv_ad
!> Purpose: apply advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv_ad(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           !< NICAS data block
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_d(nicas_blk%nc0d,geom%nl0)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Adjoint interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call nicas_blk%d(il0,its)%apply_ad(mpl,fld(:,il0,iv,its),fld_d(:,il0))
      end do
      !$omp end parallel do

      ! Halo reduction from zone D to zone A
      call nicas_blk%com_AD%red(mpl,geom%nl0,fld_d,fld(:,:,iv,its))
   end do
end do

end subroutine nicas_blk_apply_adv_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv_inv
!> Purpose: apply inverse advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv_inv(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           !< NICAS data block
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_dinv(nicas_blk%nc0dinv,geom%nl0)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Halo extension from zone A to zone Dinv
      call nicas_blk%com_ADinv%ext(mpl,geom%nl0,fld(:,:,iv,its),fld_dinv)

      ! Interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call nicas_blk%dinv(il0,its)%apply(mpl,fld_dinv(:,il0),fld(:,il0,iv,its),msdst=.false.)
      end do
      !$omp end parallel do
   end do
end do

end subroutine nicas_blk_apply_adv_inv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_adjoint
!> Purpose: test NICAS adjoint accuracy
!----------------------------------------------------------------------
subroutine nicas_blk_test_adjoint(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl           !< MPI data
type(rng_type),intent(inout) :: rng           !< Random number generator
type(nam_type),intent(in) :: nam              !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry

! Local variables
integer :: isb
real(kind_real) :: sum1,sum2
real(kind_real),allocatable :: alpha(:),alpha_save(:),alpha1(:),alpha1_save(:),alpha2(:),alpha2_save(:)
real(kind_real),allocatable :: gamma(:,:),gamma_save(:,:),delta(:,:),delta_save(:,:)
real(kind_real),allocatable :: fld(:,:),fld_save(:,:),fld1(:,:),fld1_save(:,:),fld2(:,:),fld2_save(:,:)

! Allocation
allocate(alpha(nicas_blk%nsb))
allocate(alpha_save(nicas_blk%nsb))
allocate(gamma(nicas_blk%nc1b,nicas_blk%nl1))
allocate(gamma_save(nicas_blk%nc1b,nicas_blk%nl1))
allocate(delta(nicas_blk%nc1b,geom%nl0))
allocate(delta_save(nicas_blk%nc1b,geom%nl0))
allocate(fld(geom%nc0a,geom%nl0))
allocate(fld_save(geom%nc0a,geom%nl0))

! Interpolation (subsampling)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
gamma_save = 0
do isb=1,nicas_blk%nsb
   call rng%rand_real(0.0_kind_real,1.0_kind_real,gamma_save(nicas_blk%sb_to_c1b(isb),nicas_blk%sb_to_l1(isb)))
end do

! Adjoint test
call nicas_blk%apply_interp_s(mpl,alpha_save,gamma)
call nicas_blk%apply_interp_s_ad(mpl,gamma_save,alpha)

! Print result
call mpl%dot_prod(alpha,alpha_save,sum1)
call mpl%dot_prod(gamma,gamma_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (subsampling): ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

! Interpolation (vertical)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,gamma_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,delta_save)

! Adjoint test
call nicas_blk%apply_interp_v(mpl,geom,gamma_save,delta)
call nicas_blk%apply_interp_v_ad(mpl,geom,delta_save,gamma)

! Print result
call mpl%dot_prod(gamma,gamma_save,sum1)
call mpl%dot_prod(delta,delta_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (vertical):    ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

! Interpolation (horizontal)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,delta_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call nicas_blk%apply_interp_h(mpl,geom,delta_save,fld)
call nicas_blk%apply_interp_h_ad(mpl,geom,fld_save,delta)

! Print result
call mpl%dot_prod(delta,delta_save,sum1)
call mpl%dot_prod(fld,fld_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (horizontal):  ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

! Interpolation (total)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call nicas_blk%apply_interp(mpl,geom,alpha_save,fld)
call nicas_blk%apply_interp_ad(mpl,geom,fld_save,alpha)

! Print result
call mpl%dot_prod(alpha,alpha_save,sum1)
call mpl%dot_prod(fld,fld_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (total):       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

! Allocation
allocate(alpha1(nicas_blk%nsc))
allocate(alpha1_save(nicas_blk%nsc))
allocate(alpha2(nicas_blk%nsc))
allocate(alpha2_save(nicas_blk%nsc))

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
call nicas_blk%apply_convol(mpl,alpha1)
call nicas_blk%apply_convol(mpl,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Convolution adjoint test:                 ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

! Allocation
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)
allocate(alpha1(nicas_blk%nsa))
allocate(alpha1_save(nicas_blk%nsb))
allocate(alpha2(nicas_blk%nsb))
allocate(alpha2_save(nicas_blk%nsa))

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)

! Adjoint test
call nicas_blk%com_AB%red(mpl,alpha1_save,alpha1)
call nicas_blk%com_AB%ext(mpl,alpha2_save,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AB adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

! Allocation
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)
allocate(alpha1(nicas_blk%nsa))
allocate(alpha1_save(nicas_blk%nsc))
allocate(alpha2(nicas_blk%nsc))
allocate(alpha2_save(nicas_blk%nsa))

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)

! Adjoint test
call nicas_blk%com_AC%red(mpl,alpha1_save,alpha1)
call nicas_blk%com_AC%ext(mpl,alpha2_save,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AC adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

! Allocation
allocate(fld1(geom%nc0a,geom%nl0))
allocate(fld2(geom%nc0a,geom%nl0))
allocate(fld1_save(geom%nc0,geom%nl0))
allocate(fld2_save(geom%nc0,geom%nl0))

! Generate random field
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld2_save)
fld1 = fld1_save
fld2 = fld2_save

! Adjoint test
if (nam%lsqrt) then
   call nicas_blk%apply_from_sqrt(mpl,geom,fld1)
   call nicas_blk%apply_from_sqrt(mpl,geom,fld2)
else
   call nicas_blk%apply(mpl,nam,geom,fld1)
   call nicas_blk%apply(mpl,nam,geom,fld2)
end if

! Print result
call mpl%dot_prod(fld1,fld2_save,sum1)
call mpl%dot_prod(fld2,fld1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:                       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call flush(mpl%info)

end subroutine nicas_blk_test_adjoint

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_pos_def
!> Purpose: test positive_definiteness
!----------------------------------------------------------------------
subroutine nicas_blk_test_pos_def(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
type(mpl_type),intent(in) :: mpl              !< MPI data
type(rng_type),intent(inout) :: rng           !< Random number generator
type(nam_type),intent(in) :: nam              !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry

! Local variables
integer :: iter
real(kind_real) :: norm,num,den,egvmax,egvmax_prev,egvmin,egvmin_prev
real(kind_real) :: fld(geom%nc0a,geom%nl0),fld_prev(geom%nc0a,geom%nl0)

! Power method to find the largest eigenvalue
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_prev)
call mpl%dot_prod(fld_prev,fld_prev,norm)
fld_prev = fld_prev/norm
egvmax_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld = fld_prev

   ! Apply NICAS
   if (nam%lsqrt) then
      call nicas_blk%apply_from_sqrt(mpl,geom,fld)
   else
      call nicas_blk%apply(mpl,nam,geom,fld)
   end if

   ! Compute Rayleigh quotient
   call mpl%dot_prod(fld,fld_prev,num)
   call mpl%dot_prod(fld_prev,fld_prev,den)
   egvmax = num/den

   ! Renormalize the vector
   call mpl%dot_prod(fld,fld,norm)
   fld = fld/norm

   ! Exit test
   if (abs(egvmax-egvmax_prev)<tol) exit

   ! Update
   iter = iter+1
   fld_prev = fld
   egvmax_prev = egvmax
end do

! Power method to find the smallest eigenvalue
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_prev)
call mpl%dot_prod(fld_prev,fld_prev,norm)
egvmin_prev = huge(1.0)
fld_prev = fld_prev/norm
egvmin_prev = huge(1.0)
iter = 1
do while (iter<=nitermax)
   ! Copy vector
   fld = fld_prev

   ! Apply C
   if (nam%lsqrt) then
      call nicas_blk%apply_from_sqrt(mpl,geom,fld)
   else
      call nicas_blk%apply(mpl,nam,geom,fld)
   end if
   fld = fld-egvmax*fld_prev

   ! Compute Rayleigh quotient
   call mpl%dot_prod(fld,fld_prev,num)
   call mpl%dot_prod(fld_prev,fld_prev,den)
   egvmin = num/den

   ! Renormalize the vector
   call mpl%dot_prod(fld,fld,norm)
   fld = fld/norm

   ! Exit test
   if (egvmax+egvmin<-tol*egvmax) then
      write(mpl%info,'(a7,a)') '','NICAS is not positive definite'
      call flush(mpl%info)
      exit
   end if

   ! Update
   iter = iter+1
   fld_prev = fld
   egvmin_prev = egvmin
end do

! Non conclusive test
if (iter==nitermax+1) write(mpl%info,'(a7,a,e15.8,a,i4,a,e15.8)') '','NICAS seems to be positive definite: difference ', &
 & egvmax+egvmin,' after ',nitermax,' iterations for a tolerance ',tol
call flush(mpl%info)

end subroutine nicas_blk_test_pos_def

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_sqrt
!> Purpose: test full/square-root equivalence
!----------------------------------------------------------------------
subroutine nicas_blk_test_sqrt(nicas_blk,mpl,rng,nam,geom,bpar,io,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl           !< MPI data
type(rng_type),intent(inout) :: rng           !< Random number generator
type(nam_type),intent(inout) :: nam           !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry
type(bpar_type),intent(in) :: bpar            !< Block parameters
type(io_type),intent(in) :: io                !< I/O
type(cmat_blk_type),intent(in) :: cmat_blk    !< C matrix data block

! Local variables
real(kind_real) :: fld(geom%nc0a,geom%nl0),fld_sqrt(geom%nc0a,geom%nl0)
type(nicas_blk_type) :: nicas_blk_other

! Associate
associate(ib=>nicas_blk%ib)

! Generate random field
call rng%rand_real(-1.0_kind_real,1.0_kind_real,fld)
fld_sqrt = fld

! Apply NICAS, initial version
if (nam%lsqrt) then
   call nicas_blk%apply_from_sqrt(mpl,geom,fld_sqrt)
else
   call nicas_blk%apply(mpl,nam,geom,fld)
end if

! Switch lsqrt
nam%lsqrt = .not.nam%lsqrt

! Compute NICAS parameters
call nicas_blk_other%compute_parameters(mpl,rng,nam,geom,cmat_blk)

! Apply NICAS, other version
if (nam%lsqrt) then
   call nicas_blk_other%apply_from_sqrt(mpl,geom,fld_sqrt)
else
   ! Apply NICAS
   call nicas_blk_other%apply(mpl,nam,geom,fld)
end if

! Compute dirac
if (nam%check_dirac) call nicas_blk_other%test_dirac(mpl,nam,geom,bpar,io)

! Reset lsqrt value
nam%lsqrt = .not.nam%lsqrt

! Print difference
write(mpl%info,'(a7,a,f6.1,a)') '','NICAS full / square-root error:',sqrt(sum((fld_sqrt-fld)**2)/sum(fld**2))*100.0,'%'
call flush(mpl%info)

! End associate
end associate

end subroutine nicas_blk_test_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_dirac
!> Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine nicas_blk_test_dirac(nicas_blk,mpl,nam,geom,bpar,io)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
type(mpl_type),intent(inout) :: mpl           !< MPI data
type(nam_type),intent(in) :: nam              !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry
type(bpar_type),intent(in) :: bpar            !< Block parameters
type(io_type),intent(in) :: io                !< I/O

! Local variables
integer :: il0
integer :: ic0adir(nam%ndir),iprocdir(nam%ndir),il0dir(nam%ndir),idir
real(kind_real) :: val,valmin,valmax,valmin_tot,valmax_tot
real(kind_real) :: fld(geom%nc0a,geom%nl0)
character(len=1024) :: suffix,filename

! Associate
associate(ib=>nicas_blk%ib)

! Find processor, gridpoint and level indices
call define_dirac(nam,geom,iprocdir,ic0adir,il0dir)

! Generate dirac field
fld = 0.0
do idir=1,nam%ndir
   if (iprocdir(idir)==mpl%myproc) fld(ic0adir(idir),il0dir(idir)) = 1.0
end do

! Apply NICAS method
if (nam%lsqrt) then
   call nicas_blk%apply_from_sqrt(mpl,geom,fld)
else
   call nicas_blk%apply(mpl,nam,geom,fld)
end if

if (nam%lsqrt) then
   suffix = '_sqrt'
else
   suffix = ''
end if

! Write field
filename = trim(nam%prefix)//'_dirac'
call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_dirac'//trim(suffix),fld)

! Print results
write(mpl%info,'(a7,a)') '','Values at dirac points:'
call flush(mpl%info)
do idir=1,nam%ndir
   if (iprocdir(idir)==mpl%myproc) val = fld(ic0adir(idir),il0dir(idir))
   call mpl%bcast(val,iprocdir(idir))
   write(mpl%info,'(a10,f6.1,a,f6.1,a,f10.7)') '',nam%londir(idir)*rad2deg,' / ',nam%latdir(idir)*rad2deg,': ',val
   call flush(mpl%info)
end do
write(mpl%info,'(a7,a)') '','Min - max: '
call flush(mpl%info)
do il0=1,geom%nl0
   valmin = minval(fld(:,il0),mask=geom%mask_c0a(:,il0))
   valmax = maxval(fld(:,il0),mask=geom%mask_c0a(:,il0))
   call mpl%allreduce_min(valmin,valmin_tot)
   call mpl%allreduce_max(valmax,valmax_tot)
   write(mpl%info,'(a10,a,i3,a,f10.7,a,f10.7)') '','Level ',nam%levs(il0),': ',valmin_tot,' - ',valmax_tot
   call flush(mpl%info)
end do

! End associate
end associate

end subroutine nicas_blk_test_dirac

end module type_nicas_blk
