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
module type_nicas_blk

use model_interface, only: model_write
use netcdf
use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg
use tools_display, only: prog_init,prog_print,msgerror,msgwarning,aqua,black,vunitchar
use tools_func, only: sphere_dist,vector_product,vector_triple_product,gc99
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_qsort, only: qsort
use tools_test, only: define_dirac
use type_bpar, only: bpar_type
use type_cmat_blk, only: cmat_blk_type
use type_com, only: com_type
use type_ctree, only: ctree_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mpl, only: mpl
use type_nam, only: nam_type
use type_rng, only: rng
use yomhook, only: lhook,dr_hook

implicit none

! NICAS block derived type
type nicas_blk_type
   ! Block index and name
   integer :: ib                                !< Block index
   character(len=1024) :: name                  !< Name

   ! Specific geometry
   integer :: nc1                               !< Number of points in subset Sc1
   integer,allocatable :: vbot(:)               !< Bottom level in grid Gh
   integer,allocatable :: vtop(:)               !< Top level in grid Gh
   integer,allocatable :: nc2(:)                !< Number of points in subset Sc2
   integer :: ns                                !< Number of subgrid nodes

   ! Full interpolations
   type(linop_type),allocatable :: hfull(:)     !< Horizontal interpolation
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
   integer,allocatable :: sa_to_sb(:)           !< Subgrid, halo A to halo B
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
   type(linop_type) :: c_nor                    !< Convolution (extended for normalization)

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
   integer :: mpicom                            !< Number of communication steps
contains
   procedure :: dealloc => nicas_blk_dealloc
   procedure :: compute_parameters => nicas_blk_compute_parameters
   procedure :: compute_sampling => nicas_blk_compute_sampling
   procedure :: compute_interp_h => nicas_blk_compute_interp_h
   procedure :: compute_interp_v => nicas_blk_compute_interp_v
   procedure :: compute_interp_s => nicas_blk_compute_interp_s
   procedure :: compute_mpi_ab => nicas_blk_compute_mpi_ab
   procedure :: compute_convol_network => nicas_blk_compute_convol_network
   procedure :: compute_convol_distance => nicas_blk_compute_convol_distance
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

integer,parameter :: nc1max = 15000           !< Maximum size of the Sc1 subset
real(kind_real),parameter :: sqrt_fac = 0.721 !< Square-root factor (empirical)
real(kind_real),parameter :: S_inf = 1.0e-2   !< Minimum value for the convolution coefficients
real(kind_real),parameter :: deform = 0.0     !< Deformation coefficient (maximum absolute value: -0.318)
real(kind_real),parameter :: tol = 1.0e-3     !< Positive-definiteness test tolerance
integer,parameter :: nitermax = 50            !< Number of iterations for the positive-definiteness test

private
public :: nicas_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: nicas_blk_dealloc
!> Purpose: nicas_blk object deallocation
!----------------------------------------------------------------------
subroutine nicas_blk_dealloc(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),target,intent(in) :: nam          !< Namelist
type(geom_type),target,intent(in) :: geom        !< Geometry

! Local variables
integer :: il0,il1,its

! Release memory
if (allocated(nicas_blk%vbot)) deallocate(nicas_blk%vbot)
if (allocated(nicas_blk%vtop)) deallocate(nicas_blk%vtop)
if (allocated(nicas_blk%nc2)) deallocate(nicas_blk%nc2)
if (allocated(nicas_blk%hfull)) then
   do il0=1,geom%nl0
      call nicas_blk%hfull(il0)%dealloc
   end do
   deallocate(nicas_blk%hfull)
end if
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
if (allocated(nicas_blk%c2mask)) deallocate(nicas_blk%c2mask)
if (allocated(nicas_blk%c1l1_to_s)) deallocate(nicas_blk%c1l1_to_s)
if (allocated(nicas_blk%lcheck_c1a)) deallocate(nicas_blk%lcheck_c1a)
if (allocated(nicas_blk%lcheck_c1b)) deallocate(nicas_blk%lcheck_c1b)
if (allocated(nicas_blk%lcheck_sa)) deallocate(nicas_blk%lcheck_sa)
if (allocated(nicas_blk%lcheck_sb)) deallocate(nicas_blk%lcheck_sb)
if (allocated(nicas_blk%lcheck_sc)) deallocate(nicas_blk%lcheck_sc)
if (allocated(nicas_blk%lcheck_sc_nor)) deallocate(nicas_blk%lcheck_sc_nor)
if (allocated(nicas_blk%lcheck_h)) deallocate(nicas_blk%lcheck_h)
if (allocated(nicas_blk%lcheck_s)) deallocate(nicas_blk%lcheck_s)
if (allocated(nicas_blk%c1a_to_c1)) deallocate(nicas_blk%c1a_to_c1)
if (allocated(nicas_blk%c1_to_c1a)) deallocate(nicas_blk%c1_to_c1a)
if (allocated(nicas_blk%c1b_to_c1)) deallocate(nicas_blk%c1b_to_c1)
if (allocated(nicas_blk%c1_to_c1b)) deallocate(nicas_blk%c1_to_c1b)
if (allocated(nicas_blk%c1bl1_to_sb)) deallocate(nicas_blk%c1bl1_to_sb)
if (allocated(nicas_blk%interph_lg)) deallocate(nicas_blk%interph_lg)
if (allocated(nicas_blk%interps_lg)) deallocate(nicas_blk%interps_lg)
if (allocated(nicas_blk%sa_to_sb)) deallocate(nicas_blk%sa_to_sb)
if (allocated(nicas_blk%sc_to_sb)) deallocate(nicas_blk%sc_to_sb)
if (allocated(nicas_blk%sa_to_s)) deallocate(nicas_blk%sa_to_s)
if (allocated(nicas_blk%s_to_sa)) deallocate(nicas_blk%s_to_sa)
if (allocated(nicas_blk%sb_to_s)) deallocate(nicas_blk%sb_to_s)
if (allocated(nicas_blk%s_to_sb)) deallocate(nicas_blk%s_to_sb)
if (allocated(nicas_blk%sc_to_s)) deallocate(nicas_blk%sc_to_s)
if (allocated(nicas_blk%s_to_sc)) deallocate(nicas_blk%s_to_sc)
if (allocated(nicas_blk%sc_nor_to_s)) deallocate(nicas_blk%sc_nor_to_s)
if (allocated(nicas_blk%s_to_sc_nor)) deallocate(nicas_blk%s_to_sc_nor)
if (allocated(nicas_blk%sb_to_sc_nor)) deallocate(nicas_blk%sb_to_sc_nor)
call nicas_blk%c_nor%dealloc
if (allocated(nicas_blk%sa_to_sc)) deallocate(nicas_blk%sa_to_sc)
if (allocated(nicas_blk%sb_to_sc)) deallocate(nicas_blk%sb_to_sc)
call nicas_blk%c%dealloc
if (allocated(nicas_blk%h)) then
   do il0=1,geom%nl0
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

end subroutine nicas_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_parameters
!> Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine nicas_blk_compute_parameters(nicas_blk,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: il0i,il1

! Compute adaptive sampling
write(mpl%unit,'(a7,a)') '','Compute adaptive sampling'
call nicas_blk%compute_sampling(nam,geom,cmat_blk)

! Compute horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute horizontal interpolation data'
call nicas_blk%compute_interp_h(nam,geom)

! Compute vertical interpolation data
write(mpl%unit,'(a7,a)') '','Compute vertical interpolation data'
call nicas_blk%compute_interp_v(geom)

! Compute subsampling horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute subsampling horizontal interpolation data'
call nicas_blk%compute_interp_s(nam,geom)

! Compute MPI distribution, halos A-B
write(mpl%unit,'(a7,a)') '','Compute MPI distribution, halos A-B'
call nicas_blk%compute_mpi_ab(geom)

! Compute convolution data
write(mpl%unit,'(a7,a)') '','Compute convolution data'
if (nam%network) then
   call nicas_blk%compute_convol_network(nam,geom,cmat_blk)
else
   call nicas_blk%compute_convol_distance(nam,geom,cmat_blk)
end if

! Compute MPI distribution, halo C
write(mpl%unit,'(a7,a)') '','Compute MPI distribution, halo C'
call nicas_blk%compute_mpi_c(nam,geom)

! Compute normalization
write(mpl%unit,'(a7,a)') '','Compute normalization'
call nicas_blk%compute_normalization(nam,geom)

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0 =        ',geom%nc0
write(mpl%unit,'(a10,a,i8)') '','nc0a =       ',geom%nc0a
write(mpl%unit,'(a10,a,i8)') '','nl0 =        ',geom%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =        ',nicas_blk%nc1
write(mpl%unit,'(a10,a,i8)') '','nc1a =       ',nicas_blk%nc1a
write(mpl%unit,'(a10,a,i8)') '','nc1b =       ',nicas_blk%nc1b
write(mpl%unit,'(a10,a,i8)') '','nl1 =        ',nicas_blk%nl1
do il1=1,nicas_blk%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',nicas_blk%nc2(il1)
end do
write(mpl%unit,'(a10,a,i8)') '','ns =         ',nicas_blk%ns
write(mpl%unit,'(a10,a,i8)') '','nsa =        ',nicas_blk%nsa
write(mpl%unit,'(a10,a,i8)') '','nsb =        ',nicas_blk%nsb
write(mpl%unit,'(a10,a,i8)') '','nsc =        ',nicas_blk%nsc
do il0i=1,geom%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',nicas_blk%h(il0i)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','v%n_s =      ',nicas_blk%v%n_s
do il1=1,nicas_blk%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',nicas_blk%s(il1)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','c%n_s =      ',nicas_blk%c%n_s

end subroutine nicas_blk_compute_parameters

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_sampling
!> Purpose: compute NICAS sampling
!----------------------------------------------------------------------
subroutine nicas_blk_compute_sampling(nicas_blk,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: il0,il0_prev,il1,ic0,ic1,ic2,is
integer,allocatable :: mask_c0(:),mask_c1(:),c2_to_c1(:)
real(kind_real) :: rh0sminavg,rv1min,rh0savg,distnorm
real(kind_real),allocatable :: rh0smin(:),rh0s_c1(:)
logical :: inside

! Allocation
allocate(rh0smin(geom%nc0))
allocate(mask_c0(geom%nc0))

! Reset random numbers seed
if (trim(nam%strategy)=='specific_multivariate') call rng%reseed

! Compute support radii
write(mpl%unit,'(a10,a,a,f8.2,a,f8.2,a)') '','Average support radii (H/V): ', &
 & trim(aqua),sum(cmat_blk%rh0s)/float(geom%nc0*geom%nl0)*reqkm,trim(black)//' km  / ' &
 & //trim(aqua),sum(cmat_blk%rv0s)/float(geom%nc0*geom%nl0),trim(black)//' '//trim(vunitchar)

! Basic horizontal mesh defined with the minimum support radius
rh0smin = huge(1.0)
mask_c0 = 0
do ic0=1,geom%nc0
   do il0=1,geom%nl0
      if (geom%mask(ic0,il0)) then
         rh0smin(ic0) = min(cmat_blk%rh0s(ic0,il0),rh0smin(ic0))
         mask_c0(ic0) = 1
      end if
   end do
end do
rh0sminavg = sum(rh0smin,mask=(mask_c0==1))/float(sum(mask_c0))
if (rh0sminavg>0.0) then
   nicas_blk%nc1 = floor(2.0*maxval(geom%area)*nam%resol**2/(sqrt(3.0)*rh0sminavg**2))
else
   nicas_blk%nc1 = geom%nc0
end if
nicas_blk%nc1 = min(nicas_blk%nc1,geom%nc0)
write(mpl%unit,'(a10,a,i8)') '','Estimated nc1 from horizontal support radius: ',nicas_blk%nc1
if (nicas_blk%nc1>nc1max) then
   call msgwarning('required nc1 larger than nc1max, resetting to nc1max')
   nicas_blk%nc1 = nc1max
   write(mpl%unit,'(a10,a,f5.2)') '','Effective resolution: ',sqrt(float(nicas_blk%nc1)*sqrt(3.0)*rh0sminavg**2 &
 & /(2.0*maxval(geom%area)))
end if
mask_c0 = 0
do ic0=1,geom%nc0
   if (any(geom%mask(ic0,:))) mask_c0(ic0) = 1
end do

! Compute subset
write(mpl%unit,'(a7,a)') '','Compute horizontal subset C1'
allocate(nicas_blk%c1_to_c0(nicas_blk%nc1))
if (mpl%main) call rng%initialize_sampling(geom%nc0,dble(geom%lon),dble(geom%lat),mask_c0,rh0smin,nam%ntry,nam%nrep, &
 & nicas_blk%nc1,nicas_blk%c1_to_c0)
call mpl%bcast(nicas_blk%c1_to_c0,mpl%ioproc)

! Inverse conversion
allocate(nicas_blk%c0_to_c1(geom%nc0))
call msi(nicas_blk%c0_to_c1)
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   nicas_blk%c0_to_c1(ic0) = ic1
end do

! Vertical sampling
write(mpl%unit,'(a7,a)',advance='no') '','Compute vertical subset L1: '
allocate(nicas_blk%llev(geom%nl0))
il0_prev = 1
do il0=1,geom%nl0
   ! Look for convolution levels
   if ((il0==1).or.(il0==geom%nl0)) then
      ! Keep first and last levels
      nicas_blk%llev(il0) = .true.
   else
      ! Compute normalized distance with level il0_prev
      rv1min = sqrt(0.5*(minval(cmat_blk%rv0s(nicas_blk%c1_to_c0,il0))**2+minval(cmat_blk%rv0s(nicas_blk%c1_to_c0,il0_prev))**2))
      if (rv1min>0.0) then
         distnorm = abs(geom%vunit(il0)-geom%vunit(il0_prev))/rv1min
         nicas_blk%llev(il0) = distnorm>1.0/nam%resol
      else
         nicas_blk%llev(il0) = .true.
      end if
   end if

   ! Update
   if (nicas_blk%llev(il0)) il0_prev = il0
end do
nicas_blk%nl1 = count(nicas_blk%llev)
allocate(nicas_blk%l1_to_l0(nicas_blk%nl1))
il1 = 0
do il0=1,geom%nl0
   if (nicas_blk%llev(il0)) then
      write(mpl%unit,'(i3,a)',advance='no') nam%levs(il0),' '
      il1 = il1+1
      nicas_blk%l1_to_l0(il1) = il0
   end if
end do
write(mpl%unit,'(a)') ''

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
      if (.not.inside.and.geom%mask(ic0,il0)) then
         ! Bottom level
         nicas_blk%vbot(ic1) = il0
         inside = .true.
      end if
      if (inside.and.(.not.geom%mask(ic0,il0))) then
         ! Top level
         nicas_blk%vtop(ic1) = il0
         inside = .false.
      end if
   end do
   if (nicas_blk%vbot(ic1)>nicas_blk%vtop(ic1)) call msgerror('non contiguous mask')
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
allocate(nicas_blk%c2mask(nicas_blk%nc1,nicas_blk%nl1))
allocate(rh0s_c1(nicas_blk%nc1))
allocate(mask_c1(nicas_blk%nc1))

! Horizontal subsampling
do il1=1,nicas_blk%nl1
   write(mpl%unit,'(a7,a,i3,a)') '','Compute horizontal subset C2 (level ',il1,')'

   ! Compute nc2
   il0 = nicas_blk%l1_to_l0(il1)
   mask_c1 = 0
   do ic1=1,nicas_blk%nc1
      ic0 = nicas_blk%c1_to_c0(ic1)
      if (geom%mask(ic0,il0)) then
         rh0s_c1(ic1) = cmat_blk%rh0s(ic0,il0)
         mask_c1(ic1) = 1
      end if
   end do
   rh0savg = sum(rh0s_c1,mask=(mask_c1==1))/float(sum(mask_c1))
   nicas_blk%nc2(il1) = floor(2.0*geom%area(il0)*nam%resol**2/(sqrt(3.0)*rh0savg**2))
   nicas_blk%nc2(il1) = max(nicas_blk%nc1/4,min(nicas_blk%nc2(il1),nicas_blk%nc1))

   if (nicas_blk%nc2(il1)<nicas_blk%nc1) then
      ! Allocation
      allocate(c2_to_c1(nicas_blk%nc2(il1)))

      ! Mask
      mask_c1 = 0
      do ic1=1,nicas_blk%nc1
         ic0 = nicas_blk%c1_to_c0(ic1)
         if (geom%mask(ic0,il0)) mask_c1(ic1) = 1
      end do

      ! Compute subset
      if (mpl%main) call rng%initialize_sampling(nicas_blk%nc1,dble(geom%lon(nicas_blk%c1_to_c0)), &
    & dble(geom%lat(nicas_blk%c1_to_c0)),mask_c1,rh0s_c1,nam%ntry,nam%nrep,nicas_blk%nc2(il1),c2_to_c1)
      call mpl%bcast(c2_to_c1,mpl%ioproc)

      ! Fill C2 mask
      nicas_blk%c2mask(:,il1) = .false.
      do ic2=1,nicas_blk%nc2(il1)
         ic1 = c2_to_c1(ic2)
         nicas_blk%c2mask(ic1,il1) = .true.
      end do

      ! Release memory
      deallocate(c2_to_c1)
   else
      ! Fill C2 mask
      nicas_blk%c2mask(:,il1) = .true.
   end if
end do

! Final conversions
nicas_blk%ns = sum(nicas_blk%nc2)
allocate(nicas_blk%s_to_c1(nicas_blk%ns))
allocate(nicas_blk%s_to_l1(nicas_blk%ns))
allocate(nicas_blk%c1l1_to_s(nicas_blk%nc1,nicas_blk%nl1))
call msi(nicas_blk%c1l1_to_s)
is = 0
do il1=1,nicas_blk%nl1
   do ic1=1,nicas_blk%nc1
      if (nicas_blk%c2mask(ic1,il1)) then
         is = is+1
         nicas_blk%s_to_c1(is) = ic1
         nicas_blk%s_to_l1(is) = il1
         nicas_blk%c1l1_to_s(ic1,il1) = is
      end if
   end do
end do

end subroutine nicas_blk_compute_sampling

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_interp_h
!> Purpose: compute basic horizontal interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_compute_interp_h(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il0i
type(linop_type) :: hbase

! Allocation
allocate(nicas_blk%hfull(geom%nl0i))

do il0i=1,geom%nl0i
   ! Compute grid interpolation
   call nicas_blk%hfull(il0i)%interp(geom,il0i,nicas_blk%nc1,nicas_blk%c1_to_c0,nam%mask_check, &
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
integer :: jl0,il0,jl1,il0inf,il0sup

! Initialize vertical interpolation
nicas_blk%v%prefix = 'v'
nicas_blk%v%n_src = nicas_blk%nl1
nicas_blk%v%n_dst = geom%nl0

! Linear interpolation
nicas_blk%v%n_s = nicas_blk%nl1
il0inf = 1
do jl0=1,geom%nl0
   if (nicas_blk%llev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         nicas_blk%v%n_s = nicas_blk%v%n_s+2
      end do
      il0inf = jl0
   end if
end do
call nicas_blk%v%alloc
do jl1=1,nicas_blk%nl1
   jl0 = nicas_blk%l1_to_l0(jl1)
   nicas_blk%v%row(jl1) = jl0
   nicas_blk%v%col(jl1) = jl0
   nicas_blk%v%S(jl1) = 1.0
end do
nicas_blk%v%n_s = nicas_blk%nl1
il0inf = 1
do jl0=1,geom%nl0
   if (nicas_blk%llev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         nicas_blk%v%n_s = nicas_blk%v%n_s+1
         nicas_blk%v%row(nicas_blk%v%n_s) = il0
         nicas_blk%v%col(nicas_blk%v%n_s) = il0inf
         nicas_blk%v%S(nicas_blk%v%n_s) = abs(geom%vunit(il0sup)-geom%vunit(il0)) &
                            & /abs(geom%vunit(il0sup)-geom%vunit(il0inf))

         nicas_blk%v%n_s = nicas_blk%v%n_s+1
         nicas_blk%v%row(nicas_blk%v%n_s) = il0
         nicas_blk%v%col(nicas_blk%v%n_s) = il0sup
         nicas_blk%v%S(nicas_blk%v%n_s) = abs(geom%vunit(il0)-geom%vunit(il0inf)) &
                            & /abs(geom%vunit(il0sup)-geom%vunit(il0inf))
      end do
      il0inf = jl0
   end if
end do

! Conversion
nicas_blk%v%col = nicas_blk%l0_to_l1(nicas_blk%v%col)

! Deallocate selected levels
deallocate(nicas_blk%llev)

end subroutine nicas_blk_compute_interp_v

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_interp_s
!> Purpose: compute horizontal subsampling interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_compute_interp_s(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il1,i_s
real(kind_real) :: renorm(nicas_blk%nc1)
logical :: mask_src(nicas_blk%nc1),mask_dst(nicas_blk%nc1)
logical,allocatable :: valid(:)
type(linop_type) :: stmp

! Allocation
allocate(nicas_blk%sfull(nicas_blk%nl1))

do il1=1,nicas_blk%nl1
   ! Initialize object
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
      ! Mask
      mask_src = nicas_blk%c2mask(:,il1)
      mask_dst = .true.

      ! Compute interpolation
      call stmp%interp(nicas_blk%nc1,geom%lon(nicas_blk%c1_to_c0),geom%lat(nicas_blk%c1_to_c0),mask_src, &
    & nicas_blk%nc1,geom%lon(nicas_blk%c1_to_c0),geom%lat(nicas_blk%c1_to_c0),mask_dst,nam%nicas_interp)

      ! Allocation
      allocate(valid(stmp%n_s))
      valid = .true.

      ! Check mask boundaries
      if (nam%mask_check) then
         write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Sublevel ',il1,': '
         call stmp%interp_check_mask(geom,valid,nicas_blk%l1_to_l0(il1),row_to_ic0=nicas_blk%c1_to_c0, &
       & col_to_ic0=nicas_blk%c1_to_c0)
      else
         write(mpl%unit,'(a10,a,i3)') '','Sublevel ',il1
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
      call nicas_blk%sfull(il1)%interp_missing(nicas_blk%nc1,geom%lon(nicas_blk%c1_to_c0),geom%lat(nicas_blk%c1_to_c0),mask_dst, &
    & nam%nicas_interp)
   end if
end do

end subroutine nicas_blk_compute_interp_s

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_ab
!> Purpose: compute NICAS MPI distribution, halos A-B
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_ab(nicas_blk,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il0i,ic0,iproc,ic1,jc1,ic1a,ic1b,il0,il1,isa,isb,i_s,i_s_loc,is,js,h_n_s_max,s_n_s_max,h_n_s_max_loc,s_n_s_max_loc
integer :: nsa,nsb
integer,allocatable :: interph_lg(:,:),interps_lg(:,:)
logical :: lcheck_c1b_h(nicas_blk%nc1)
integer,allocatable :: sa_to_s(:),sb_to_s(:),sa_to_sb(:)
type(com_type) :: com_AB(mpl%nproc)

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
allocate(nicas_blk%lcheck_c1a(nicas_blk%nc1))
allocate(nicas_blk%lcheck_c1b(nicas_blk%nc1))
allocate(nicas_blk%lcheck_sa(nicas_blk%ns))
allocate(nicas_blk%lcheck_sb(nicas_blk%ns))
allocate(nicas_blk%lcheck_h(h_n_s_max,geom%nl0i))
allocate(nicas_blk%lcheck_s(s_n_s_max,nicas_blk%nl1))

! Halo definitions

! Halo A
nicas_blk%lcheck_c1a = .false.
nicas_blk%lcheck_sa = .false.
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   ic0 = nicas_blk%c1_to_c0(ic1)
   il1 = nicas_blk%s_to_l1(is)
   il0 = nicas_blk%l1_to_l0(il1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) then
      nicas_blk%lcheck_c1a(ic1) = .true.
      nicas_blk%lcheck_sa(is) = .true.
   end if
end do

! Halo B

! Horizontal interpolation
nicas_blk%lcheck_h = .false.
lcheck_c1b_h = .false.
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%hfull(il0i)%n_s
      ic0 = nicas_blk%hfull(il0i)%row(i_s)
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%myproc) then
         jc1 = nicas_blk%hfull(il0i)%col(i_s)
         nicas_blk%lcheck_h(i_s,il0i) = .true.
         lcheck_c1b_h(jc1) = .true.
      end if
   end do
end do

! Subsampling horizontal interpolation
nicas_blk%lcheck_sb = .false.
nicas_blk%lcheck_s = .false.
nicas_blk%lcheck_c1b = lcheck_c1b_h
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%sfull(il1)%n_s
      ic1 = nicas_blk%sfull(il1)%row(i_s)
      if (lcheck_c1b_h(ic1)) then
         jc1 = nicas_blk%sfull(il1)%col(i_s)
         js = nicas_blk%c1l1_to_s(jc1,il1)
         nicas_blk%lcheck_c1b(jc1) = .true.
         nicas_blk%lcheck_sb(js) = .true.
         nicas_blk%lcheck_s(i_s,il1) = .true.
      end if
   end do
end do

! Check halos consistency
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is).and.(.not.nicas_blk%lcheck_sb(is))) call msgerror('point in halo A but not in halo B')
end do

! Sizes
nicas_blk%nc1a = count(nicas_blk%lcheck_c1a)
nicas_blk%nsa = count(nicas_blk%lcheck_sa)
do il0i=1,geom%nl0i
   nicas_blk%h(il0i)%n_s = count(nicas_blk%lcheck_h(:,il0i))
end do
nicas_blk%nc1b = count(nicas_blk%lcheck_c1b)
nicas_blk%nsb = count(nicas_blk%lcheck_sb)
do il1=1,nicas_blk%nl1
   nicas_blk%s(il1)%n_s = count(nicas_blk%lcheck_s(:,il1))
end do

! Global <-> local conversions for fields

! Halo A
allocate(nicas_blk%c1a_to_c1(nicas_blk%nc1a))
allocate(nicas_blk%c1_to_c1a(nicas_blk%nc1))
call msi(nicas_blk%c1_to_c1a)
ic1a = 0
do ic1=1,nicas_blk%nc1
   if (nicas_blk%lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      nicas_blk%c1a_to_c1(ic1a) = ic1
      nicas_blk%c1_to_c1a(ic1) = ic1a
   end if
end do

allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
isa = 0
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is)) then
      isa = isa+1
      nicas_blk%sa_to_s(isa) = is
   end if
end do

! Halo B
allocate(nicas_blk%c1b_to_c1(nicas_blk%nc1b))
allocate(nicas_blk%c1_to_c1b(nicas_blk%nc1))
call msi(nicas_blk%c1_to_c1b)
ic1b = 0
do ic1=1,nicas_blk%nc1
   if (nicas_blk%lcheck_c1b(ic1)) then
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
      if (nicas_blk%lcheck_h(i_s,il0i)) then
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
      if (nicas_blk%lcheck_s(i_s,il1)) then
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
   call nicas_blk%h(il0i)%reorder
end do

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
   call nicas_blk%s(il1)%reorder
end do

! Conversions
allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
allocate(nicas_blk%c1bl1_to_sb(nicas_blk%nc1b,nicas_blk%nl1))
call msi(nicas_blk%c1bl1_to_sb)
do isb=1,nicas_blk%nsb
   is = nicas_blk%sb_to_s(isb)
   il1 = nicas_blk%s_to_l1(is)
   ic1 = nicas_blk%s_to_c1(is)
   ic1b = nicas_blk%c1_to_c1b(ic1)
   nicas_blk%sb_to_c1b(isb) = ic1b
   nicas_blk%sb_to_l1(isb) = il1
   nicas_blk%c1bl1_to_sb(ic1b,il1) = isb
end do

! Get global distribution of the subgrid on ioproc
if (mpl%main) then
   ! Allocation
   allocate(nicas_blk%s_to_sa(nicas_blk%ns))

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         nsa = nicas_blk%nsa
      else
         ! Receive dimension on ioproc
         call mpl%recv(nsa,iproc,mpl%tag)
      end if

      ! Allocation
      allocate(sa_to_s(nsa))

      if (iproc==mpl%ioproc) then
         ! Copy data
         sa_to_s = nicas_blk%sa_to_s
      else
         ! Receive data on ioproc
         call mpl%recv(nsa,sa_to_s,iproc,mpl%tag+1)
      end if

      ! Fill s_to_sa
      do isa=1,nsa
         nicas_blk%s_to_sa(sa_to_s(isa)) = isa
      end do

      ! Release memory
      deallocate(sa_to_s)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(nicas_blk%nsa,mpl%ioproc,mpl%tag)

   ! Send data to ioproc
   call mpl%send(nicas_blk%nsa,nicas_blk%sa_to_s,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nsa = nicas_blk%nsa
         nsb = nicas_blk%nsb
      else
         ! Receive dimensions on ioproc
         call mpl%recv(nsa,iproc,mpl%tag)
         call mpl%recv(nsb,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(sb_to_s(nsb))
      allocate(sa_to_sb(nsa))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         sb_to_s = nicas_blk%sb_to_s
         sa_to_sb = nicas_blk%sa_to_sb
      else
         ! Receive data on ioproc
         call mpl%recv(nsb,sb_to_s,iproc,mpl%tag+2)
         call mpl%recv(nsa,sa_to_sb,iproc,mpl%tag+3)
      end if

      ! Allocation
      com_AB(iproc)%nred = nsa
      com_AB(iproc)%next = nsb
      allocate(com_AB(iproc)%ext_to_proc(com_AB(iproc)%next))
      allocate(com_AB(iproc)%ext_to_red(com_AB(iproc)%next))
      allocate(com_AB(iproc)%red_to_ext(com_AB(iproc)%nred))

      ! AB communication
      do isb=1,nsb
         is = sb_to_s(isb)
         ic1 = nicas_blk%s_to_c1(is)
         ic0 = nicas_blk%c1_to_c0(ic1)
         com_AB(iproc)%ext_to_proc(isb) = geom%c0_to_proc(ic0)
         isa = nicas_blk%s_to_sa(is)
         com_AB(iproc)%ext_to_red(isb) = isa
      end do
      com_AB(iproc)%red_to_ext = sa_to_sb

      ! Release memory
      deallocate(sb_to_s)
      deallocate(sa_to_sb)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(nicas_blk%nsa,mpl%ioproc,mpl%tag)
   call mpl%send(nicas_blk%nsb,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl%send(nicas_blk%nsb,nicas_blk%sb_to_s,mpl%ioproc,mpl%tag+2)
   call mpl%send(nicas_blk%nsa,nicas_blk%sa_to_sb,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4
call nicas_blk%com_AB%setup(com_AB,'com_AB')

end subroutine nicas_blk_compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_network
!> Purpose: compute convolution with a network approach
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_network(nicas_blk,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: n_s_max,progint,ithread,is,ic1,ic0,inr,jl0,jl1,np,np_new,i,j,k,ip,jc0,jc1,il0,djl0,il1,kc0,dkl0,kl0,jp,js,isb,offset
integer :: c_n_s(mpl%nthread),c_nor_n_s(mpl%nthread)
integer,allocatable :: net_nnb(:),net_inb(:,:),plist(:,:),plist_new(:,:)
real(kind_real) :: dnb,disthsq,distvsq,rhsq,rvsq,distnorm,disttest,S_test
real(kind_real),allocatable :: net_dnb(:,:,:,:),dist(:,:)
logical :: init,add_to_front,valid_arc
logical,allocatable :: done(:),valid(:,:)
type(linop_type) :: c(mpl%nthread),c_nor(mpl%nthread)

! Count neighbors
net_nnb = 0
do ic0=1,geom%nc0
   inr = geom%mesh%order_inv(ic0)
   i = geom%mesh%lend(inr)
   init = .true.
   do while ((i/=geom%mesh%lend(inr)).or.init)
      net_nnb(ic0) = net_nnb(ic0)+1
      i = geom%mesh%lptr(i)
      init = .false.
   end do
end do

! Allocation
allocate(net_inb(maxval(net_nnb),geom%nc0))
allocate(net_dnb(maxval(net_nnb),-1:1,geom%nc0,geom%nl0))

! Find neighbors
net_nnb = 0
do ic0=1,geom%nc0
   inr = geom%mesh%order_inv(ic0)
   i = geom%mesh%lend(inr)
   init = .true.
   do while ((i/=geom%mesh%lend(inr)).or.init)
      net_nnb(ic0) = net_nnb(ic0)+1
      net_inb(net_nnb(ic0),ic0) = geom%mesh%order(abs(geom%mesh%list(i)))
      i = geom%mesh%lptr(i)
      init = .false.
   end do
end do

! Compute distances
net_dnb = huge(1.0)
!$omp parallel do schedule(static) private(ic0,j,jc0,dnb,il0,valid_arc,djl0,jl0,disthsq,distvsq,rhsq,rvsq,distnorm)
do ic0=1,geom%nc0
   do j=1,net_nnb(ic0)
      ! Index
      jc0 = net_inb(j,ic0)

      if (any(geom%mask(jc0,:))) then
         ! True distance
         call sphere_dist(geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),dnb)

         do il0=1,geom%nl0
            ! Check mask bounds
            valid_arc = .true.
            if (nam%mask_check) call geom%check_arc(il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),valid_arc)

            if (valid_arc) then
               do djl0=-1,1
                  ! Index
                  jl0 = max(1,min(il0+djl0,geom%nl0))

                  if (geom%mask(ic0,il0).and.geom%mask(jc0,jl0)) then
                     ! Squared support radii
                     disthsq = dnb**2
                     distvsq = geom%distv(il0,jl0)**2
                     rhsq = 0.5*(cmat_blk%rh0(ic0,il0)**2+cmat_blk%rh0(jc0,jl0)**2)
                     rvsq = 0.5*(cmat_blk%rv0(ic0,il0)**2+cmat_blk%rv0(jc0,jl0)**2)
                     distnorm = 0.0
                     if (rhsq>0.0) then
                        distnorm = distnorm+disthsq/rhsq
                     elseif (disthsq>0.0) then
                        distnorm = distnorm+0.5*huge(0.0)
                     end if
                     if (rvsq>0.0) then
                        distnorm = distnorm+distvsq/rvsq
                     elseif (distvsq>0.0) then
                        distnorm = distnorm+0.5*huge(0.0)
                     end if
                     net_dnb(j,djl0,ic0,il0) = sqrt(distnorm)
                  end if
               end do
            end if
         end do
      end if
   end do
end do
!$omp end parallel do

! Allocation
n_s_max = 100*nint(float(geom%nc0*geom%nl0)/float(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   c_nor(ithread)%n_s = n_s_max
   call c(ithread)%alloc
   call c_nor(ithread)%alloc
end do
allocate(done(nicas_blk%nsb))

! Compute weights
write(mpl%unit,'(a10,a)',advance='no') '','Compute weights: '
call prog_init(progint,done)
c_n_s = 0
c_nor_n_s = 0
!$omp parallel do schedule(static) private(isb,is,ithread,ic1,il1,ic0,il0,np,np_new,ip,jc0,jl0,k,kc0,kl0,disttest,add_to_front), &
!$omp&                             private(jp,js,jc1,jl1,distnorm,S_test) firstprivate(plist,plist_new,dist,valid)
do isb=1,nicas_blk%nsb
   ! Indices
   is = nicas_blk%sb_to_s(isb)
   ithread = omp_get_thread_num()+1
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)
   ic0 = nicas_blk%c1_to_c0(ic1)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Allocation
   allocate(plist(geom%nc0*geom%nl0,2))
   allocate(plist_new(geom%nc0*geom%nl0,2))
   allocate(dist(geom%nc0,geom%nl0))
   allocate(valid(geom%nc0,geom%nl0))

   ! Initialize the front
   np = 1
   plist(1,1) = ic0
   plist(1,2) = il0
   dist = 1.0
   dist(ic0,il0) = 0.0
   valid = .false.
   valid(ic0,il0) = .true.

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         jc0 = plist(ip,1)
         jl0 = plist(ip,2)

         ! Loop over neighbors
         do k=1,net_nnb(jc0)
            kc0 = net_inb(k,jc0)
            do dkl0=-1,1
               kl0 = max(1,min(jl0+dkl0,geom%nl0))
               disttest = dist(jc0,jl0)+net_dnb(k,dkl0,jc0,jl0)
               if (disttest<1.0) then
                  ! Point is inside the support
                  if (disttest<dist(kc0,kl0)) then
                     ! Update distance
                     dist(kc0,kl0) = disttest
                     valid(kc0,kl0) = .true.

                     ! Check if the point should be added to the front (avoid duplicates)
                     add_to_front = .true.
                     do jp=1,np_new
                        if ((plist_new(jp,1)==kc0).and.(plist_new(jp,2)==kl0)) then
                           add_to_front = .false.
                           exit
                        end if
                     end do

                     if (add_to_front) then
                        ! Add point to the front
                        np_new = np_new+1
                        plist_new(np_new,1) = kc0
                        plist_new(np_new,2) = kl0
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

   ! Count convolution operations
   do js=1,nicas_blk%ns
      ! Indices
      jc1 = nicas_blk%s_to_c1(js)
      jl1 = nicas_blk%s_to_l1(js)
      jc0 = nicas_blk%c1_to_c0(jc1)
      jl0 = nicas_blk%l1_to_l0(jl1)

      if (valid(jc0,jl0)) then
         ! Normalized distance
         distnorm = dist(jc0,jl0)

         if (distnorm<1.0) then
            ! Distance deformation
            distnorm = distnorm+deform*sin(pi*distnorm)

            ! Square-root
            if (nam%lsqrt) distnorm = distnorm/sqrt_fac

            ! Gaspari-Cohn (1999) function
            S_test = gc99(distnorm)

            if (S_test>S_inf) then
               ! Store coefficient for convolution
               if (nam%mpicom==1) then
                  if ((nicas_blk%lcheck_sb(js).and.(is<js)).or.(.not.nicas_blk%lcheck_sb(js))) &
                & call c(ithread)%add_op(c_n_s(ithread),is,js,S_test)
               elseif (nam%mpicom==2) then
                  if (nicas_blk%lcheck_sa(is).and.((nicas_blk%lcheck_sa(js).and.(is<js)).or.(.not.nicas_blk%lcheck_sa(js)))) &
                & call c(ithread)%add_op(c_n_s(ithread),is,js,S_test)
               end if

               ! Store coefficient for normalization
               if (nam%lsqrt) then
                  if ((nicas_blk%lcheck_sb(js).and.(is<js)).or.(.not.nicas_blk%lcheck_sb(js))) &
               &  call c_nor(ithread)%add_op(c_nor_n_s(ithread),is,js,S_test)
               else
                  if (nicas_blk%lcheck_sb(js).and.(is<js)) call c_nor(ithread)%add_op(c_nor_n_s(ithread),is,js,S_test)
               end if
            end if
         end if
      end if
   end do

   ! Print progression
   done(isb) = .true.
   call prog_print(progint,done)

   ! Release memory
   deallocate(plist)
   deallocate(plist_new)
   deallocate(dist)
   deallocate(valid)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Initialize object
nicas_blk%c%prefix = 'c'
nicas_blk%c_nor%prefix = 'c_nor'

! Gather data
nicas_blk%c%n_s = sum(c_n_s)
nicas_blk%c_nor%n_s = sum(c_nor_n_s)
call nicas_blk%c%alloc
call nicas_blk%c_nor%alloc

! Gather convolution data from OpenMP threads
offset = 0
do ithread=1,mpl%nthread
   nicas_blk%c%row(offset+1:offset+c_n_s(ithread)) = c(ithread)%row(1:c_n_s(ithread))
   nicas_blk%c%col(offset+1:offset+c_n_s(ithread)) = c(ithread)%col(1:c_n_s(ithread))
   nicas_blk%c%S(offset+1:offset+c_n_s(ithread)) = c(ithread)%S(1:c_n_s(ithread))
   offset = offset+c_n_s(ithread)
end do
offset = 0
do ithread=1,mpl%nthread
   nicas_blk%c_nor%row(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%row(1:c_nor_n_s(ithread))
   nicas_blk%c_nor%col(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%col(1:c_nor_n_s(ithread))
   nicas_blk%c_nor%S(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%S(1:c_nor_n_s(ithread))
   offset = offset+c_nor_n_s(ithread)
end do

! Release memory
do ithread=1,mpl%nthread
   call c(ithread)%dealloc
   call c_nor(ithread)%dealloc
end do

end subroutine nicas_blk_compute_convol_network

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_distance
!> Purpose: compute convolution with a distance approach
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_distance(nicas_blk,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: ms,n_s_max,progint,ithread,is,ic1,jl0,jc1,il1,jl1,il0,j,js,isb,ic1b,offset
integer :: c_n_s(mpl%nthread),c_nor_n_s(mpl%nthread)
integer,allocatable :: nn_index(:,:)
real(kind_real) :: disthsq,distvsq,rhsq,rvsq,distnorm,S_test
real(kind_real),allocatable :: nn_dist(:,:)
logical :: force,test,valid,submask(nicas_blk%nc1,nicas_blk%nl1)
logical,allocatable :: mask_ctree(:),done(:)
type(ctree_type) :: ctree
type(linop_type) :: c(mpl%nthread),c_nor(mpl%nthread)

! Define submask
submask = .false.
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)
   submask(ic1,il1) = .true.
end do

! Compute cover tree
write(mpl%unit,'(a10,a)') '','Compute cover tree'
allocate(mask_ctree(nicas_blk%nc1))
mask_ctree = .true.
call ctree%create(nicas_blk%nc1,dble(geom%lon(nicas_blk%c1_to_c0)),dble(geom%lat(nicas_blk%c1_to_c0)),mask_ctree)
deallocate(mask_ctree)

! Theoretical number of neighbors
ms = 0
do il1=1,nicas_blk%nl1
   ms = max(ms,floor(0.25*float(nicas_blk%nc1)*(sum(cmat_blk%rh0(nicas_blk%c1_to_c0,nicas_blk%l1_to_l0(il1)), &
      & mask=nicas_blk%c2mask(:,il1))/float(nicas_blk%nc2(il1)))**2))
end do

! Find all necessary neighbors
valid = .false.
do while (.not.valid)
   ! Check number of neighbors
   force = ms>nicas_blk%nc1
   ms = max(1,min(ms,nicas_blk%nc1))
   write(mpl%unit,'(a10,a,i8,a,f5.1,a)') '','Try with ',ms,' neighbors (',float(ms)/float(nicas_blk%nc1)*100.0,'%)'

   ! Allocation
   if (allocated(nn_index)) deallocate(nn_index)
   if (allocated(nn_dist)) deallocate(nn_dist)
   allocate(nn_index(ms,nicas_blk%nc1b))
   allocate(nn_dist(ms,nicas_blk%nc1b))

   ! Loop over points
   test = .true.
   ic1b = 1
   do while ((test.or.force).and.(ic1b<=nicas_blk%nc1b))
      ic1 = nicas_blk%c1b_to_c1(ic1b)

      ! Compute nearest neighbors
      call ctree%find_nearest_neighbors(dble(geom%lon(nicas_blk%c1_to_c0(ic1))), &
    & dble(geom%lat(nicas_blk%c1_to_c0(ic1))),ms,nn_index(:,ic1b),nn_dist(:,ic1b))

      ! Loop over levels
      jc1 = nn_index(ms,ic1b)
      il1 = 1
      do while ((test.or.force).and.(il1<=nicas_blk%nl1))
         rhsq = 0.5*(cmat_blk%rh0(nicas_blk%c1_to_c0(ic1),nicas_blk%l1_to_l0(il1))**2+ &
               & cmat_blk%rh0(nicas_blk%c1_to_c0(jc1),nicas_blk%l1_to_l0(il1))**2)
         distnorm = nn_dist(ms,ic1b)**2/rhsq
         if ((distnorm<1.0).and.(.not.force)) then
            ms = 2*ms
            test = .false.
         end if
         il1 = il1+1
      end do
      ic1b = ic1b+1
   end do
   if ((ic1b>nicas_blk%nc1b).or.force) valid = .true.
end do

! Allocation
n_s_max = 100*nint(float(geom%nc0*geom%nl0)/float(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   c_nor(ithread)%n_s = n_s_max
   call c(ithread)%alloc
   call c_nor(ithread)%alloc
end do
allocate(done(nicas_blk%nsb))

! Compute weights
write(mpl%unit,'(a10,a)',advance='no') '','Compute weights: '
call prog_init(progint,done)
c_n_s = 0
c_nor_n_s = 0
!$omp parallel do schedule(static) private(isb,is,ithread,ic1,ic1b,il1,il0,j,jc1,jl1,jl0,js,disthsq,distvsq,rhsq,rvsq), &
!$omp&                             private(distnorm,S_test)
do isb=1,nicas_blk%nsb
   ! Indices
   is = nicas_blk%sb_to_s(isb)
   ithread = omp_get_thread_num()+1
   ic1 = nicas_blk%s_to_c1(is)
   ic1b = nicas_blk%c1_to_c1b(ic1)
   il1 = nicas_blk%s_to_l1(is)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Loop on nearest neighbors
   do j=1,ms
      jc1 = nn_index(j,ic1b)
      do jl1=1,nicas_blk%nl1
         if (submask(jc1,jl1)) then
            jl0 = nicas_blk%l1_to_l0(jl1)
            js = nicas_blk%c1l1_to_s(jc1,jl1)

            ! Normalized distance
            disthsq = nn_dist(j,ic1b)**2
            distvsq = geom%distv(il0,jl0)**2
            rhsq = 0.5*(cmat_blk%rh0(nicas_blk%c1_to_c0(nicas_blk%s_to_c1(is)),nicas_blk%l1_to_l0(nicas_blk%s_to_l1(is)))**2+ &
                 & cmat_blk%rh0(nicas_blk%c1_to_c0(nicas_blk%s_to_c1(js)),nicas_blk%l1_to_l0(nicas_blk%s_to_l1(js)))**2)
            rvsq = 0.5*(cmat_blk%rv0(nicas_blk%c1_to_c0(nicas_blk%s_to_c1(is)),nicas_blk%l1_to_l0(nicas_blk%s_to_l1(is)))**2+ &
                 & cmat_blk%rv0(nicas_blk%c1_to_c0(nicas_blk%s_to_c1(js)),nicas_blk%l1_to_l0(nicas_blk%s_to_l1(js)))**2)
            distnorm = 0.0
            if (rhsq>0.0) then
               distnorm = distnorm+disthsq/rhsq
            elseif (disthsq>0.0) then
               distnorm = distnorm+0.5*huge(0.0)
            end if
            if (rvsq>0.0) then
               distnorm = distnorm+distvsq/rvsq
            elseif (distvsq>0.0) then
               distnorm = distnorm+0.5*huge(0.0)
            end if
            distnorm = sqrt(distnorm)

            if (distnorm<1.0) then
               ! Distance deformation
               distnorm = distnorm+deform*sin(pi*distnorm)

               ! Square-root
               if (nam%lsqrt) distnorm = distnorm/sqrt_fac

               ! Gaspari-Cohn (1999) function
               S_test = gc99(distnorm)

               if (S_test>S_inf) then
                  ! Store coefficient for convolution
                  if (nam%mpicom==1) then
                     if ((nicas_blk%lcheck_sb(js).and.(is<js)).or.(.not.nicas_blk%lcheck_sb(js))) &
                   & call c(ithread)%add_op(c_n_s(ithread),is,js,S_test)
                  elseif (nam%mpicom==2) then
                     if (nicas_blk%lcheck_sa(is).and.((nicas_blk%lcheck_sa(js).and.(is<js)).or.(.not.nicas_blk%lcheck_sa(js)))) &
                   & call c(ithread)%add_op(c_n_s(ithread),is,js,S_test)
                  end if

                  ! Store coefficient for normalization
                  if (nam%lsqrt) then
                     if ((nicas_blk%lcheck_sb(js).and.(is<js)).or.(.not.nicas_blk%lcheck_sb(js))) &
                  &  call c_nor(ithread)%add_op(c_nor_n_s(ithread),is,js,S_test)
                  else
                     if (nicas_blk%lcheck_sb(js).and.(is<js)) call c_nor(ithread)%add_op(c_nor_n_s(ithread),is,js,S_test)
                  end if
               end if
            end if
         end if
      end do
   end do

   ! Print progression
   done(isb) = .true.
   call prog_print(progint,done)
end do
!$omp end parallel do
write(mpl%unit,'(a)') '100%'

! Initialize object
nicas_blk%c%prefix = 'c'
nicas_blk%c_nor%prefix = 'c_nor'

! Gather data
nicas_blk%c%n_s = sum(c_n_s)
nicas_blk%c_nor%n_s = sum(c_nor_n_s)
call nicas_blk%c%alloc
call nicas_blk%c_nor%alloc

! Gather convolution data from OpenMP threads
offset = 0
do ithread=1,mpl%nthread
   nicas_blk%c%row(offset+1:offset+c_n_s(ithread)) = c(ithread)%row(1:c_n_s(ithread))
   nicas_blk%c%col(offset+1:offset+c_n_s(ithread)) = c(ithread)%col(1:c_n_s(ithread))
   nicas_blk%c%S(offset+1:offset+c_n_s(ithread)) = c(ithread)%S(1:c_n_s(ithread))
   offset = offset+c_n_s(ithread)
end do
offset = 0
do ithread=1,mpl%nthread
   nicas_blk%c_nor%row(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%row(1:c_nor_n_s(ithread))
   nicas_blk%c_nor%col(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%col(1:c_nor_n_s(ithread))
   nicas_blk%c_nor%S(offset+1:offset+c_nor_n_s(ithread)) = c_nor(ithread)%S(1:c_nor_n_s(ithread))
   offset = offset+c_nor_n_s(ithread)
end do

! Release memory
do ithread=1,mpl%nthread
   call c(ithread)%dealloc
   call c_nor(ithread)%dealloc
end do

end subroutine nicas_blk_compute_convol_distance

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_c
!> Purpose: compute NICAS MPI distribution, halo C
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_c(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: iproc,ic0,ic1,isa,isb,isc,i_s,is,js,nsa,nsc
integer,allocatable :: sc_to_s(:),sa_to_sc(:)
type(com_type) :: com_AC(mpl%nproc)

! Allocation
allocate(nicas_blk%lcheck_sc(nicas_blk%ns))
allocate(nicas_blk%lcheck_sc_nor(nicas_blk%ns))

! Halo definitions

! Halo C
if (nam%mpicom==1) then
   ! 1 communication step
   nicas_blk%lcheck_sc = nicas_blk%lcheck_sb
   do i_s=1,nicas_blk%c%n_s
      is = nicas_blk%c%row(i_s)
      js = nicas_blk%c%col(i_s)
      nicas_blk%lcheck_sc(is) = .true.
      nicas_blk%lcheck_sc(js) = .true.
   end do
elseif (nam%mpicom==2) then
   ! 2 communication steps
   nicas_blk%lcheck_sc = nicas_blk%lcheck_sb
   do i_s=1,nicas_blk%c%n_s
      is = nicas_blk%c%row(i_s)
      js = nicas_blk%c%col(i_s)
      nicas_blk%lcheck_sc(is) = .true.
      nicas_blk%lcheck_sc(js) = .true.
   end do
end if
nicas_blk%nsc = count(nicas_blk%lcheck_sc)
nicas_blk%lcheck_sc_nor = nicas_blk%lcheck_sb
do i_s=1,nicas_blk%c_nor%n_s
   is = nicas_blk%c_nor%row(i_s)
   js = nicas_blk%c_nor%col(i_s)
   nicas_blk%lcheck_sc_nor(is) = .true.
   nicas_blk%lcheck_sc_nor(js) = .true.
end do
nicas_blk%nsc_nor = count(nicas_blk%lcheck_sc_nor)

! Check halos consistency
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is).and.(.not.nicas_blk%lcheck_sc(is))) call msgerror('point in halo A but not in halo C')
   if (nicas_blk%lcheck_sb(is).and.(.not.nicas_blk%lcheck_sc(is))) call msgerror('point in halo B but not in halo C')
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
   if (nicas_blk%lcheck_sc_nor(is)) then
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
call nicas_blk%c%reorder

! Convolution for normalization
nicas_blk%c_nor%n_src = nicas_blk%nsc
nicas_blk%c_nor%n_dst = nicas_blk%nsc
do i_s=1,nicas_blk%c_nor%n_s
   nicas_blk%c_nor%row(i_s) = nicas_blk%s_to_sc_nor(nicas_blk%c_nor%row(i_s))
   nicas_blk%c_nor%col(i_s) = nicas_blk%s_to_sc_nor(nicas_blk%c_nor%col(i_s))
end do
call nicas_blk%c_nor%reorder

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nsa = nicas_blk%nsa
         nsc = nicas_blk%nsc
      else
         ! Receive dimensions on ioproc
         call mpl%recv(nsa,iproc,mpl%tag)
         call mpl%recv(nsc,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(sc_to_s(nsc))
      allocate(sa_to_sc(nsa))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         sc_to_s = nicas_blk%sc_to_s
         sa_to_sc = nicas_blk%sa_to_sc
      else
         ! Receive data on ioproc
         call mpl%recv(nsc,sc_to_s,iproc,mpl%tag+2)
         call mpl%recv(nsa,sa_to_sc,iproc,mpl%tag+3)
      end if

      ! Allocation
      com_AC(iproc)%nred = nsa
      com_AC(iproc)%next = nsc
      allocate(com_AC(iproc)%ext_to_proc(com_AC(iproc)%next))
      allocate(com_AC(iproc)%ext_to_red(com_AC(iproc)%next))
      allocate(com_AC(iproc)%red_to_ext(com_AC(iproc)%nred))

      ! AC communication
      do isc=1,nsc
         is = sc_to_s(isc)
         ic1 = nicas_blk%s_to_c1(is)
         ic0 = nicas_blk%c1_to_c0(ic1)
         com_AC(iproc)%ext_to_proc(isc) = geom%c0_to_proc(ic0)
         isa = nicas_blk%s_to_sa(is)
         com_AC(iproc)%ext_to_red(isc) = isa
      end do
      com_AC(iproc)%red_to_ext = sa_to_sc

      ! Release memory
      deallocate(sc_to_s)
      deallocate(sa_to_sc)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(nicas_blk%nsa,mpl%ioproc,mpl%tag)
   call mpl%send(nicas_blk%nsc,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl%send(nicas_blk%nsc,nicas_blk%sc_to_s,mpl%ioproc,mpl%tag+2)
   call mpl%send(nicas_blk%nsa,nicas_blk%sa_to_sc,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4
call nicas_blk%com_AC%setup(com_AC,'com_AC')

! Copy mpicom
nicas_blk%mpicom = nam%mpicom

end subroutine nicas_blk_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_normalization
!> Purpose: compute normalization
!----------------------------------------------------------------------
subroutine nicas_blk_compute_normalization(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry

! Local variables
integer :: il0i,i_s,ic1,ic1b,jc1b,is,js,isc,jsb,jsc,ic0,ic0a,il0,il1,ih,jv,nlr,ilr,jlr,ic,isc_add,progint
integer,allocatable :: ineh(:,:),inev(:),ines(:,:),inec(:),order(:),isc_list(:),order2(:)
integer,allocatable :: h_col(:,:,:),v_col(:,:),s_col(:,:,:),c_ind(:,:)
real(kind_real) :: S_add
real(kind_real),allocatable :: h_S(:,:,:),v_S(:,:),s_S(:,:,:),c_S(:,:)
real(kind_real),allocatable :: list(:),S_list(:),S_list_tmp(:)
logical :: conv
logical,allocatable :: done(:)

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
allocate(v_S(maxval(inev),geom%nl0))
inev = 0
do i_s=1,nicas_blk%v%n_s
   il0 = nicas_blk%v%row(i_s)
   inev(il0) = inev(il0)+1
   v_col(inev(il0),il0) = nicas_blk%v%col(i_s)
   v_S(inev(il0),il0) = nicas_blk%v%S(i_s)
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
   is = nicas_blk%sc_nor_to_s(isc)
   js = nicas_blk%sc_nor_to_s(jsc)
   inec(isc) = inec(isc)+1
   inec(jsc) = inec(jsc)+1
end do
allocate(c_ind(maxval(inec),nicas_blk%nsc_nor))
allocate(c_S(maxval(inec),nicas_blk%nsc_nor))
call msi(c_ind)
call msr(c_S)
inec = 0
do i_s=1,nicas_blk%c_nor%n_s
   isc = nicas_blk%c_nor%col(i_s)
   jsc = nicas_blk%c_nor%row(i_s)
   is = nicas_blk%sc_nor_to_s(isc)
   js = nicas_blk%sc_nor_to_s(jsc)
   inec(isc) = inec(isc)+1
   c_ind(inec(isc),isc) = jsc
   c_S(inec(isc),isc) = nicas_blk%c_nor%S(i_s)
   inec(jsc) = inec(jsc)+1
   c_ind(inec(jsc),jsc) = isc
   c_S(inec(jsc),jsc) = nicas_blk%c_nor%S(i_s)
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
allocate(done(geom%nc0a))
call msr(nicas_blk%norm)

! Compute normalization weights
do il0=1,geom%nl0
   il0i = min(il0,geom%nl0i)
   write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),': '
   call prog_init(progint,done)

   !$omp parallel do schedule(static) private(ic0a,ic0,nlr,isc_add,S_add,ih,ic1b,ic1,jv,il1,is,ilr,conv,ic,jlr,isc,jsc), &
   !$omp&                             firstprivate(isc_list,order2,S_list,S_list_tmp)
   do ic0a=1,geom%nc0a
      ! Index
      ic0 = geom%c0a_to_c0(ic0a)

      if (geom%mask(ic0,il0)) then
         ! Allocation
         allocate(isc_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         allocate(S_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         if (nam%lsqrt) then
            allocate(S_list_tmp(nicas_blk%nsc_nor))
         else
            allocate(order2(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
            allocate(S_list_tmp(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         end if

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
                     S_add = h_S(ih,ic0a,il0i)*v_S(jv,il0)*s_S(is,ic1b,il1)
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

         if (nam%lsqrt) then
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
         else
            ! Sort arrays
            call qsort(nlr,isc_list(1:nlr),order2(1:nlr))
            S_list(1:nlr) = S_list(order2(1:nlr))

            ! Convolution
            S_list_tmp(1:nlr) = S_list(1:nlr)
            do ilr=1,nlr
               isc = isc_list(ilr)
               conv = .true.
               ic = 1
               do jlr=ilr+1,nlr
                  jsc = isc_list(jlr)
                  if (ic<=inec(isc)) then
                     do while (c_ind(ic,isc)/=jsc)
                        ic = ic+1
                        if (ic>inec(isc)) then
                           conv = .false.
                           exit
                        end if
                     end do
                  else
                     conv = .false.
                  end if
                  if (conv) then
                     S_list(ilr) = S_list(ilr)+c_S(ic,isc)*S_list_tmp(jlr)
                     S_list(jlr) = S_list(jlr)+c_S(ic,isc)*S_list_tmp(ilr)
                  end if
               end do
            end do

            ! Normalization factor
            if (nlr>0) nicas_blk%norm(ic0a,il0) = 1.0/sqrt(sum(S_list(1:nlr)*S_list_tmp(1:nlr)))
         end if

         ! Print progression
         done(ic0a) = .true.
         call prog_print(progint,done)

          ! Release memory
         deallocate(isc_list)
         if (.not.nam%lsqrt) deallocate(order2)
         deallocate(S_list)
         deallocate(S_list_tmp)
      end if
   end do
   !$omp end parallel do
   write(mpl%unit,'(a)') '100%'
end do

! Release memory
deallocate(done)

end subroutine nicas_blk_compute_normalization

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_adv
!> Purpose: compute advection
!----------------------------------------------------------------------
subroutine nicas_blk_compute_adv(nicas_blk,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam                 !< Namelist
type(geom_type),intent(in) :: geom               !< Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       !< C matrix data block

! Local variables
integer :: its,il0,ic0,ic0a,nc0a,i_s,i_s_loc,iproc,jc0
integer :: ic0d,nc0d,d_n_s_max,d_n_s_max_loc
integer :: ic0dinv,nc0dinv,dinv_n_s_max,dinv_n_s_max_loc
integer,allocatable :: c0d_to_c0(:),c0_to_c0d(:),c0a_to_c0d(:),interpd_lg(:,:,:)
integer,allocatable :: c0dinv_to_c0(:),c0_to_c0dinv(:),c0a_to_c0dinv(:),interpdinv_lg(:,:,:)
integer,allocatable :: c0d_to_c0_copy(:),c0a_to_c0d_copy(:)
integer,allocatable :: c0dinv_to_c0_copy(:),c0a_to_c0dinv_copy(:)
logical :: mask_c0(geom%nc0)
logical,allocatable :: lcheck_c0d(:),lcheck_d(:,:,:)
logical,allocatable :: lcheck_c0dinv(:),lcheck_dinv(:,:,:)
type(linop_type),allocatable :: dfull(:,:)
type(linop_type),allocatable :: dinvfull(:,:)
type(com_type) :: com_AD(mpl%nproc)
type(com_type) :: com_ADinv(mpl%nproc)

write(mpl%unit,'(a7,a)') '','Compute advection'

! Allocation
allocate(dfull(geom%nl0,2:nam%nts))
allocate(dinvfull(geom%nl0,2:nam%nts))

! Initialization
mask_c0 = .true.

do its=2,nam%nts
   do il0=1,geom%nl0
      ! Compute direct interpolation
      call dfull(il0,its)%interp(geom%nc0,cmat_blk%displ_lon(:,il0,its),cmat_blk%displ_lat(:,il0,its),mask_c0,geom%nc0, &
    & geom%lon,geom%lat,mask_c0,nam%diag_interp)

      ! Compute inverse interpolation
      call dinvfull(il0,its)%interp(geom%nc0,geom%lon,geom%lat,mask_c0,geom%nc0,cmat_blk%displ_lon(:,il0,its), &
    & cmat_blk%displ_lat(:,il0,its),mask_c0,nam%diag_interp)
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
      call nicas_blk%d(il0,its)%reorder
      call nicas_blk%dinv(il0,its)%reorder
   end do
end do

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = geom%nc0a
         nc0d = nicas_blk%nc0d
         nc0dinv = nicas_blk%nc0dinv
      else
         ! Receive dimensions on ioproc
         call mpl%recv(nc0a,iproc,mpl%tag)
         call mpl%recv(nc0d,iproc,mpl%tag+1)
         call mpl%recv(nc0dinv,iproc,mpl%tag+2)
      end if

      ! Allocation
      allocate(c0d_to_c0_copy(nc0d))
      allocate(c0dinv_to_c0_copy(nc0dinv))
      allocate(c0a_to_c0d_copy(nc0a))
      allocate(c0a_to_c0dinv_copy(nc0a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0d_to_c0_copy = c0d_to_c0
         c0dinv_to_c0_copy = c0dinv_to_c0
         c0a_to_c0d_copy = c0a_to_c0d
         c0a_to_c0dinv_copy = c0a_to_c0dinv
      else
         ! Receive data on ioproc
         call mpl%recv(nc0d,c0d_to_c0_copy,iproc,mpl%tag+3)
         call mpl%recv(nc0dinv,c0dinv_to_c0_copy,iproc,mpl%tag+4)
         call mpl%recv(nc0a,c0a_to_c0d_copy,iproc,mpl%tag+5)
         call mpl%recv(nc0a,c0a_to_c0dinv_copy,iproc,mpl%tag+6)
      end if

      ! Allocation
      com_AD(iproc)%nred = nc0a
      com_ADinv(iproc)%nred = nc0a
      com_AD(iproc)%next = nc0d
      com_ADinv(iproc)%next = nc0dinv
      allocate(com_AD(iproc)%ext_to_proc(com_AD(iproc)%next))
      allocate(com_ADinv(iproc)%ext_to_proc(com_ADinv(iproc)%next))
      allocate(com_AD(iproc)%ext_to_red(com_AD(iproc)%next))
      allocate(com_ADinv(iproc)%ext_to_red(com_ADinv(iproc)%next))
      allocate(com_AD(iproc)%red_to_ext(com_AD(iproc)%nred))
      allocate(com_ADinv(iproc)%red_to_ext(com_ADinv(iproc)%nred))

      ! AD communication
      do ic0d=1,nc0d
         ic0 = c0d_to_c0_copy(ic0d)
         com_AD(iproc)%ext_to_proc(ic0d) = geom%c0_to_proc(ic0)
         ic0a = geom%c0_to_c0a(ic0)
         com_AD(iproc)%ext_to_red(ic0d) = ic0a
      end do
      com_AD(iproc)%red_to_ext = c0a_to_c0d_copy

      ! ADinv communication
      do ic0dinv=1,nc0dinv
         ic0 = c0dinv_to_c0_copy(ic0dinv)
         com_ADinv(iproc)%ext_to_proc(ic0dinv) = geom%c0_to_proc(ic0)
         ic0a = geom%c0_to_c0a(ic0)
         com_ADinv(iproc)%ext_to_red(ic0dinv) = ic0a
      end do
      com_ADinv(iproc)%red_to_ext = c0a_to_c0dinv_copy

      ! Release memory
      deallocate(c0d_to_c0_copy)
      deallocate(c0dinv_to_c0_copy)
      deallocate(c0a_to_c0d_copy)
      deallocate(c0a_to_c0dinv_copy)
   end do
else
   ! Send dimensions to ioproc
   call mpl%send(geom%nc0a,mpl%ioproc,mpl%tag)
   call mpl%send(nicas_blk%nc0d,mpl%ioproc,mpl%tag+1)
   call mpl%send(nicas_blk%nc0dinv,mpl%ioproc,mpl%tag+2)

   ! Send data to ioproc
   call mpl%send(nicas_blk%nc0d,c0d_to_c0,mpl%ioproc,mpl%tag+3)
   call mpl%send(nicas_blk%nc0dinv,c0dinv_to_c0,mpl%ioproc,mpl%tag+4)
   call mpl%send(geom%nc0a,c0a_to_c0d,mpl%ioproc,mpl%tag+5)
   call mpl%send(geom%nc0a,c0a_to_c0dinv,mpl%ioproc,mpl%tag+6)
end if
mpl%tag = mpl%tag+7
call nicas_blk%com_AD%setup(com_AD,'com_AD')
call nicas_blk%com_ADinv%setup(com_ADinv,'com_ADinv')

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0d =       ',nicas_blk%nc0d
write(mpl%unit,'(a10,a,i8)') '','nc0dinv =    ',nicas_blk%nc0dinv
do its=2,nam%nts
   do il0=1,geom%nl0
      write(mpl%unit,'(a10,a,i3,a,i2,a,i8)') '','d(',il0,',',its,')%n_s =    ',nicas_blk%d(il0,its)%n_s
      write(mpl%unit,'(a10,a,i3,a,i2,a,i8)') '','dinv(',il0,',',its,')%n_s = ',nicas_blk%dinv(il0,its)%n_s
   end do
end do

end subroutine nicas_blk_compute_adv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply
!> Purpose: apply NICAS method
!----------------------------------------------------------------------
subroutine nicas_blk_apply(nicas_blk,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk            !< NICAS data block
type(geom_type),intent(in) :: geom                       !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: zhook_handle
real(kind_real) :: alpha_a(nicas_blk%nsa),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

if (lhook) call dr_hook('nicas_blk_apply',0,zhook_handle)

! Adjoint interpolation
call nicas_blk%apply_interp_ad(geom,fld,alpha_b)

! Communication
if (nicas_blk%mpicom==1) then
   ! Initialization
   alpha_c = 0.0

   ! Copy zone B into zone C
   alpha_c(nicas_blk%sb_to_sc) = alpha_b
elseif (nicas_blk%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call nicas_blk%com_AB%red(alpha_b,alpha_a)

   ! Initialization
   alpha_c = 0.0

   ! Copy zone A into zone C
   alpha_c(nicas_blk%sa_to_sc) = alpha_a
end if

! Convolution
call nicas_blk%apply_convol(alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(geom,alpha_b,fld)

if (lhook) call dr_hook('nicas_blk_apply',1,zhook_handle)

end subroutine nicas_blk_apply

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_from_sqrt
!> Purpose: apply NICAS method from its square-root formulation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_from_sqrt(nicas_blk,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk            !< NICAS data block
type(geom_type),intent(in) :: geom                       !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: zhook_handle
real(kind_real) :: alpha(nicas_blk%nsa)

if (lhook) call dr_hook('nicas_blk_apply_from_sqrt',0,zhook_handle)

! Apply square-root adjoint
call nicas_blk%apply_sqrt_ad(geom,fld,alpha)

! Apply square-root
call nicas_blk%apply_sqrt(geom,alpha,fld)

if (lhook) call dr_hook('nicas_blk_apply_from_sqrt',1,zhook_handle)

end subroutine nicas_blk_apply_from_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_sqrt
!> Purpose: apply NICAS method square-root
!----------------------------------------------------------------------
subroutine nicas_blk_apply_sqrt(nicas_blk,geom,alpha,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk          !< NICAS data block
type(geom_type),intent(in) :: geom                     !< Geometry
real(kind_real),intent(in) :: alpha(nicas_blk%nsa)     !< Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variable
real(kind_real) :: zhook_handle
real(kind_real) :: alpha_a(nicas_blk%nsa),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

if (lhook) call dr_hook('nicas_blk_apply_sqrt',0,zhook_handle)

! Initialization
alpha_c = 0.0

! Copy zone A into zone C
alpha_c(nicas_blk%sa_to_sc) = alpha

! Convolution
call nicas_blk%apply_convol(alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(geom,alpha_b,fld)

if (lhook) call dr_hook('nicas_blk_apply_sqrt',1,zhook_handle)

end subroutine nicas_blk_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_sqrt_ad
!> Purpose: apply NICAS method square-root adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_sqrt_ad(nicas_blk,geom,fld,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         !< NICAS data block
type(geom_type),intent(in) :: geom                    !< Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field
real(kind_real),intent(out) :: alpha(nicas_blk%nsa)   !< Subgrid field

! Local variable
real(kind_real) :: zhook_handle
real(kind_real) :: alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

if (lhook) call dr_hook('nicas_blk_apply_sqrt_ad',0,zhook_handle)

! Adjoint interpolation
call nicas_blk%apply_interp_ad(geom,fld,alpha_b)

! Halo reduction from zone B to zone A
call nicas_blk%com_AB%red(alpha_b,alpha)

! Initialization
alpha_c = 0.0

! Copy zone A into zone C
alpha_c(nicas_blk%sa_to_sc) = alpha

! Convolution
call nicas_blk%apply_convol(alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(alpha_c,alpha)

if (lhook) call dr_hook('nicas_blk_apply_sqrt_ad',1,zhook_handle)

end subroutine nicas_blk_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp
!> Purpose: apply interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp(nicas_blk,geom,alpha,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk          !< NICAS data block
type(geom_type),intent(in) :: geom                     !< Geometry
real(kind_real),intent(in) :: alpha(nicas_blk%nsb)     !< Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: gamma(nicas_blk%nc1b,nicas_blk%nl1),delta(nicas_blk%nc1b,geom%nl0)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('nicas_blk_apply_interp',0,zhook_handle)

! Subsampling horizontal interpolation
call nicas_blk%apply_interp_s(alpha,gamma)

! Vertical interpolation
call nicas_blk%apply_interp_v(geom,gamma,delta)

! Horizontal interpolation
call nicas_blk%apply_interp_h(geom,delta,fld)

! Normalization
fld = fld*nicas_blk%norm

if (lhook) call dr_hook('nicas_blk_apply_interp',1,zhook_handle)

end subroutine nicas_blk_apply_interp

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_ad
!> Purpose: apply interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_ad(nicas_blk,geom,fld,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         !< NICAS data block
type(geom_type),intent(in) :: geom                    !< Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field
real(kind_real),intent(out) :: alpha(nicas_blk%nsb)   !< Subgrid field

! Local variables
real(kind_real) :: gamma(nicas_blk%nc1b,nicas_blk%nl1),delta(nicas_blk%nc1b,geom%nl0)
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('nicas_blk_apply_interp_ad',0,zhook_handle)

! Normalization
fld_tmp = fld*nicas_blk%norm

! Horizontal interpolation
call nicas_blk%apply_interp_h_ad(geom,fld_tmp,delta)

! Vertical interpolation
call nicas_blk%apply_interp_v_ad(geom,delta,gamma)

! Subsampling horizontal interpolation
call nicas_blk%apply_interp_s_ad(gamma,alpha)

if (lhook) call dr_hook('nicas_blk_apply_interp_ad',1,zhook_handle)

end subroutine nicas_blk_apply_interp_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_h
!> Purpose: apply horizontal interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_h(nicas_blk,geom,delta,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                !< NICAS data block
type(geom_type),intent(in) :: geom                           !< Geometry
real(kind_real),intent(in) :: delta(nicas_blk%nc1b,geom%nl0) !< Subset Sc1 field, full levels
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0)       !< Field

! Local variables
integer :: il0

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call nicas_blk%h(min(il0,geom%nl0i))%apply(delta(:,il0),fld(:,il0))
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_h

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_h_ad
!> Purpose: apply horizontal interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_h_ad(nicas_blk,geom,fld,delta)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                 !< NICAS data block
type(geom_type),intent(in) :: geom                            !< Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0)         !< Field
real(kind_real),intent(out) :: delta(nicas_blk%nc1b,geom%nl0) !< Subset Sc1 field, full levels

! Local variables
integer :: il0

!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call nicas_blk%h(min(il0,geom%nl0i))%apply_ad(fld(:,il0),delta(:,il0))
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_h_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_v
!> Purpose: apply vertical interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_v(nicas_blk,geom,gamma,delta)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                     !< NICAS data block
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
   call nicas_blk%v%apply(gamma_tmp,delta_tmp)

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
subroutine nicas_blk_apply_interp_v_ad(nicas_blk,geom,delta,gamma)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                      !< NICAS data block
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
   call nicas_blk%v%apply_ad(delta_tmp,gamma_tmp)

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
subroutine nicas_blk_apply_interp_s(nicas_blk,alpha,gamma)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                      !< NICAS data block
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
   call nicas_blk%s(il1)%apply(beta(:,il1),gamma(:,il1))
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_s

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_s_ad
!> Purpose: apply subsampling interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_s_ad(nicas_blk,gamma,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                     !< NICAS data block
real(kind_real),intent(in) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) !< Subset Sc1 field, limited levels
real(kind_real),intent(out) :: alpha(nicas_blk%nsb)               !< Subgrid field

! Local variables
integer :: il1,isb
real(kind_real) :: beta(nicas_blk%nc1b,nicas_blk%nl1)

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,nicas_blk%nl1
   call nicas_blk%s(il1)%apply_ad(gamma(:,il1),beta(:,il1))
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
subroutine nicas_blk_apply_convol(nicas_blk,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         !< NICAS data block
real(kind_real),intent(inout) :: alpha(nicas_blk%nsc) !< Subgrid field

! Local variables
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('nicas_blk_apply_convol',0,zhook_handle)

! Apply linear operator, symmetric
call nicas_blk%c%apply_sym(alpha)

if (lhook) call dr_hook('nicas_blk_apply_convol',1,zhook_handle)

end subroutine nicas_blk_apply_convol

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv
!> Purpose: apply advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv(nicas_blk,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           !< NICAS data block
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_d(nicas_blk%nc0d,geom%nl0)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('nicas_blk_apply_adv',0,zhook_handle)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Halo extension from zone A to zone D
      call nicas_blk%com_AD%ext(geom%nl0,fld(:,:,iv,its),fld_d)

      ! Interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call nicas_blk%d(il0,its)%apply(fld_d(:,il0),fld(:,il0,iv,its))
      end do
      !$omp end parallel do
   end do
end do

if (lhook) call dr_hook('nicas_blk_apply_adv',1,zhook_handle)

end subroutine nicas_blk_apply_adv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv_ad
!> Purpose: apply advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv_ad(nicas_blk,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           !< NICAS data block
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_d(nicas_blk%nc0d,geom%nl0)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('nicas_blk_apply_adv_ad',0,zhook_handle)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Adjoint interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call nicas_blk%d(il0,its)%apply_ad(fld(:,il0,iv,its),fld_d(:,il0))
      end do
      !$omp end parallel do

      ! Halo reduction from zone D to zone A
      call nicas_blk%com_AD%red(geom%nl0,fld_d,fld(:,:,iv,its))
   end do
end do

if (lhook) call dr_hook('nicas_blk_apply_adv_ad',1,zhook_handle)

end subroutine nicas_blk_apply_adv_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv_inv
!> Purpose: apply inverse advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv_inv(nicas_blk,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           !< NICAS data block
type(nam_type),target,intent(in) :: nam                                 !< Namelist
type(geom_type),target,intent(in) :: geom                               !< Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_dinv(nicas_blk%nc0dinv,geom%nl0)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('nicas_blk_apply_adv_inv',0,zhook_handle)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Halo extension from zone A to zone Dinv
      call nicas_blk%com_ADinv%ext(geom%nl0,fld(:,:,iv,its),fld_dinv)

      ! Interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         call nicas_blk%dinv(il0,its)%apply(fld_dinv(:,il0),fld(:,il0,iv,its))
      end do
      !$omp end parallel do
   end do
end do

if (lhook) call dr_hook('nicas_blk_apply_adv_inv',1,zhook_handle)

end subroutine nicas_blk_apply_adv_inv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_adjoint
!> Purpose: test NICAS adjoint accuracy
!----------------------------------------------------------------------
subroutine nicas_blk_test_adjoint(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam              !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry

! Local variables
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
call rng%rand_real(0.0_kind_real,1.0_kind_real,gamma_save)

! Adjoint test
call nicas_blk%apply_interp_s(alpha_save,gamma)
call nicas_blk%apply_interp_s_ad(gamma_save,alpha)

! Print result
call mpl%dot_prod(alpha,alpha_save,sum1)
call mpl%dot_prod(gamma,gamma_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (subsampling): ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Interpolation (vertical)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,gamma_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,delta_save)

! Adjoint test
call nicas_blk%apply_interp_v(geom,gamma_save,delta)
call nicas_blk%apply_interp_v_ad(geom,delta_save,gamma)

! Print result
call mpl%dot_prod(gamma,gamma_save,sum1)
call mpl%dot_prod(delta,delta_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (vertical):    ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Interpolation (horizontal)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,delta_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call nicas_blk%apply_interp_h(geom,delta_save,fld)
call nicas_blk%apply_interp_h_ad(geom,fld_save,delta)

! Print result
call mpl%dot_prod(delta,delta_save,sum1)
call mpl%dot_prod(fld,fld_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (horizontal):  ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! Interpolation (total)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call nicas_blk%apply_interp(geom,alpha_save,fld)
call nicas_blk%apply_interp_ad(geom,fld_save,alpha)

! Print result
call mpl%dot_prod(alpha,alpha_save,sum1)
call mpl%dot_prod(fld,fld_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (total):       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

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
call nicas_blk%apply_convol(alpha1)
call nicas_blk%apply_convol(alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Convolution adjoint test:                 ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

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
call nicas_blk%com_AB%red(alpha1_save,alpha1)
call nicas_blk%com_AB%ext(alpha2_save,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AB adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

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
call nicas_blk%com_AC%red(alpha1_save,alpha1)
call nicas_blk%com_AC%ext(alpha2_save,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AC adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

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
   call nicas_blk%apply_from_sqrt(geom,fld1)
   call nicas_blk%apply_from_sqrt(geom,fld2)
else
   call nicas_blk%apply(geom,fld1)
   call nicas_blk%apply(geom,fld2)
end if

! Print result
call mpl%dot_prod(fld1,fld2_save,sum1)
call mpl%dot_prod(fld2,fld1_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:                       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

end subroutine nicas_blk_test_adjoint

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_pos_def
!> Purpose: test positive_definiteness
!----------------------------------------------------------------------
subroutine nicas_blk_test_pos_def(nicas_blk,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
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
      call nicas_blk%apply_from_sqrt(geom,fld)
   else
      call nicas_blk%apply(geom,fld)
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
      call nicas_blk%apply_from_sqrt(geom,fld)
   else
      call nicas_blk%apply(geom,fld)
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
      write(mpl%unit,'(a7,a)') '','NICAS is not positive definite'
      exit
   end if

   ! Update
   iter = iter+1
   fld_prev = fld
   egvmin_prev = egvmin
end do

! Non conclusive test
if (iter==nitermax+1) write(mpl%unit,'(a7,a,e15.8,a,i4,a,e15.8)') '','NICAS seems to be positive definite: difference ', &
 & egvmax+egvmin,' after ',nitermax,' iterations for a tolerance ',tol

end subroutine nicas_blk_test_pos_def

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_sqrt
!> Purpose: test full/square-root equivalence
!----------------------------------------------------------------------
subroutine nicas_blk_test_sqrt(nicas_blk,nam,geom,bpar,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
type(nam_type),intent(inout) :: nam           !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry
type(bpar_type),intent(in) :: bpar            !< Block parameters
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
   call nicas_blk%apply_from_sqrt(geom,fld_sqrt)
else
   call nicas_blk%apply(geom,fld)
end if

! Switch lsqrt
nam%lsqrt = .not.nam%lsqrt

! Compute NICAS parameters
call nicas_blk_other%compute_parameters(nam,geom,cmat_blk)

! Apply NICAS, other version
if (nam%lsqrt) then
   call nicas_blk_other%apply_from_sqrt(geom,fld_sqrt)
else
   ! Apply NICAS
   call nicas_blk_other%apply(geom,fld)
end if

! Compute dirac
if (nam%check_dirac) call nicas_blk_other%test_dirac(nam,geom,bpar)

! Reset lsqrt value
nam%lsqrt = .not.nam%lsqrt

! Print difference
write(mpl%unit,'(a7,a,f6.1,a)') '','NICAS full / square-root error:',sqrt(sum((fld_sqrt-fld)**2)/sum(fld**2))*100.0,'%'

! End associate
end associate

end subroutine nicas_blk_test_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_dirac
!> Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine nicas_blk_test_dirac(nicas_blk,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk !< NICAS data block
type(nam_type),intent(in) :: nam              !< Namelist
type(geom_type),intent(in) :: geom            !< Geometry
type(bpar_type),intent(in) :: bpar            !< Block parameters

! Local variables
integer :: ic0,ic0a,ic1,il0,il1,is,iproc
integer :: ic0adir(nam%ndir),iprocdir(nam%ndir),il0dir(nam%ndir),idir
real(kind_real) :: val,valmin,valmax,valmin_tot,valmax_tot
real(kind_real) :: fld(geom%nc0a,geom%nl0),fld_c1(geom%nc0a,geom%nl0),fld_s(geom%nc0a,geom%nl0)
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
   call nicas_blk%apply_from_sqrt(geom,fld)
else
   call nicas_blk%apply(geom,fld)
end if

if (nam%lsqrt) then
   suffix = '_sqrt'
else
   suffix = ''
end if

! Write field
filename = trim(nam%prefix)//'_dirac.nc'
call geom%fld_write(nam,trim(nam%prefix)//'_dirac.nc',trim(bpar%blockname(ib))//'_dirac'//trim(suffix),fld)

! Write gridded field
filename = trim(nam%prefix)//'_dirac_gridded.nc'
call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_dirac'//trim(suffix),fld)

! Print results
write(mpl%unit,'(a7,a)') '','Values at dirac points:'
do idir=1,nam%ndir
   if (iprocdir(idir)==mpl%myproc) val = fld(ic0adir(idir),il0dir(idir))
   call mpl%bcast(val,iprocdir(idir))
   write(mpl%unit,'(a10,f6.1,a,f6.1,a,f10.7)') '',nam%londir(idir),' / ',nam%latdir(idir),': ',val
end do
write(mpl%unit,'(a7,a)') '','Min - max: '
do il0=1,geom%nl0
   valmin = minval(fld(:,il0),mask=geom%mask(geom%c0a_to_c0,il0))
   valmax = maxval(fld(:,il0),mask=geom%mask(geom%c0a_to_c0,il0))
   call mpl%allreduce_min(valmin,valmin_tot)
   call mpl%allreduce_max(valmax,valmax_tot)
   write(mpl%unit,'(a10,a,i3,a,f10.7,a,f10.7)') '','Level ',nam%levs(il0),': ',valmin_tot,' - ',valmax_tot
end do

if (nam%new_param) then
   ! Write gridded field at interpolation points
   call msr(fld_c1)
   do ic1=1,nicas_blk%nc1
      ic0 = nicas_blk%c1_to_c0(ic1)
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%myproc) then
         ic0a = geom%c0_to_c0a(ic0)
         fld_c1(ic0a,:) = fld(ic0a,:)
      end if
   end do
   call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_dirac_c1'//trim(suffix),fld_c1)

   ! Write field at subgrid points
   call msr(fld_s)
   do is=1,nicas_blk%ns
      ic1 = nicas_blk%s_to_c1(is)
      ic0 = nicas_blk%c1_to_c0(ic1)
      il1 = nicas_blk%s_to_l1(is)
      il0 = nicas_blk%l1_to_l0(il1)
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%myproc) then
         ic0a = geom%c0_to_c0a(ic0)
         fld_s(ic0a,il0) = fld(ic0a,il0)
      end if
   end do
   call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_dirac_s'//trim(suffix),fld_s)
end if

! End associate
end associate

end subroutine nicas_blk_test_dirac

end module type_nicas_blk
