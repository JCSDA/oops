!----------------------------------------------------------------------
! Module: type_adv
! Purpose: advection derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_adv

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_const, only: req,reqkm,rad2deg,deg2rad
use tools_func, only: lonlatmod,sphere_dist,reduce_arc,lonlat2xyz,xyz2lonlat,vector_product
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: inf,sup
use type_bpar, only: bpar_type
use type_com, only: com_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type
use type_rng, only: rng_type

implicit none

logical,parameter :: cor_tracker = .false.                 ! Tracker method
real(kind_real),parameter :: max_std_ratio = 2.5_kind_real ! Minimum ratio between correlation maximum and local standard deviation

! Advection derived type
type adv_type
   integer :: niter                                   ! Number of stored iterations
   real(kind_real),allocatable :: search_rad(:,:)     ! Search radius
   real(kind_real),allocatable :: lon_c2a(:,:)        ! Longitude origin
   real(kind_real),allocatable :: lat_c2a(:,:)        ! Latitude origin
   real(kind_real),allocatable :: lon_c2a_raw(:,:,:)  ! Raw advected longitude
   real(kind_real),allocatable :: lat_c2a_raw(:,:,:)  ! Raw advected latitude
   real(kind_real),allocatable :: dist_c2a_raw(:,:,:) ! Raw advection distance
   real(kind_real),allocatable :: dist_raw(:,:)       ! Averaged raw advection distance
   real(kind_real),allocatable :: valid_raw(:,:)      ! Averaged raw advection validity
   real(kind_real),allocatable :: lon_c2a_flt(:,:,:)  ! Filtered advected longitude
   real(kind_real),allocatable :: lat_c2a_flt(:,:,:)  ! Filtered advected latitude
   real(kind_real),allocatable :: dist_c2a_flt(:,:,:) ! Filtered advection distance
   real(kind_real),allocatable :: dist_flt(:,:)       ! Averaged filtered advection distance
   real(kind_real),allocatable :: valid_flt(:,:)      ! Averaged filtered advection validity
   real(kind_real),allocatable :: rhflt(:,:)          ! Optimal filtering support radius
   real(kind_real),allocatable :: cor_loc(:,:,:,:)    ! Local correlation
   real(kind_real),allocatable :: cor_max(:,:,:,:)    ! Maximum correlation
   real(kind_real),allocatable :: cor_adv(:,:,:,:)    ! Advected correlation
contains
   procedure :: alloc => adv_alloc
   procedure :: dealloc => adv_dealloc
   procedure :: compute => adv_compute
   procedure :: compute_raw => adv_compute_raw
   procedure :: filter => adv_filter
   procedure :: interp => adv_interp
   procedure :: test => adv_test
   procedure :: write => adv_write
end type adv_type

private
public :: adv_type

contains

!----------------------------------------------------------------------
! Subroutine: adv_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine adv_alloc(adv,nam,geom,samp)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv ! Advection
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(samp_type),intent(in) :: samp   ! Sampling

! Allocation
allocate(adv%search_rad(geom%nc0,geom%nl0))
allocate(adv%lon_c2a(samp%nc2a,geom%nl0))
allocate(adv%lat_c2a(samp%nc2a,geom%nl0))
allocate(adv%lon_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%lat_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_raw(geom%nl0,2:nam%nts))
allocate(adv%valid_raw(geom%nl0,2:nam%nts))
allocate(adv%lon_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%lat_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_flt(geom%nl0,2:nam%nts))
allocate(adv%valid_flt(geom%nl0,2:nam%nts))
allocate(adv%rhflt(geom%nl0,2:nam%nts))
allocate(adv%cor_loc(geom%nc0a,geom%nl0,nam%nv,2:nam%nts))
allocate(adv%cor_max(geom%nc0a,geom%nl0,nam%nv,2:nam%nts))
allocate(adv%cor_adv(geom%nc0a,geom%nl0,nam%nv,2:nam%nts))

end subroutine adv_alloc

!----------------------------------------------------------------------
! Subroutine: adv_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine adv_dealloc(adv)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv ! Advection

! Release memory
if (allocated(adv%search_rad)) deallocate(adv%search_rad)
if (allocated(adv%lon_c2a)) deallocate(adv%lon_c2a)
if (allocated(adv%lat_c2a)) deallocate(adv%lat_c2a)
if (allocated(adv%lon_c2a_raw)) deallocate(adv%lon_c2a_raw)
if (allocated(adv%lat_c2a_raw)) deallocate(adv%lat_c2a_raw)
if (allocated(adv%dist_c2a_raw)) deallocate(adv%dist_c2a_raw)
if (allocated(adv%dist_raw)) deallocate(adv%dist_raw)
if (allocated(adv%valid_raw)) deallocate(adv%valid_raw)
if (allocated(adv%lon_c2a_flt)) deallocate(adv%lon_c2a_flt)
if (allocated(adv%lat_c2a_flt)) deallocate(adv%lat_c2a_flt)
if (allocated(adv%dist_c2a_flt)) deallocate(adv%dist_c2a_flt)
if (allocated(adv%dist_flt)) deallocate(adv%dist_flt)
if (allocated(adv%valid_flt)) deallocate(adv%valid_flt)
if (allocated(adv%rhflt)) deallocate(adv%rhflt)
if (allocated(adv%cor_loc)) deallocate(adv%cor_loc)
if (allocated(adv%cor_max)) deallocate(adv%cor_max)
if (allocated(adv%cor_adv)) deallocate(adv%cor_adv)

end subroutine adv_dealloc

!----------------------------------------------------------------------
! Subroutine: adv_compute
! Purpose: compute advection
!----------------------------------------------------------------------
subroutine adv_compute(adv,mpl,rng,nam,geom,samp,ens)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv  ! Advection
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(rng_type),intent(inout) :: rng   ! Random number generator
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(samp_type),intent(inout) :: samp ! Sampling
type(ens_type), intent(in) :: ens     ! Ensemble

! Allocation
call adv%alloc(nam,geom,samp)

! Compute raw advection
call adv%compute_raw(mpl,nam,geom,samp,ens)

if (nam%adv_niter>0) then
   ! Filter advection
   call adv%filter(mpl,nam,geom,samp)
else
   ! Copy
   adv%lon_c2a_flt = adv%lon_c2a_raw
   adv%lat_c2a_flt = adv%lat_c2a_raw
   adv%dist_c2a_flt = adv%dist_c2a_raw
   adv%dist_flt = adv%dist_raw
   adv%valid_flt = adv%valid_raw
   adv%rhflt = mpl%msv%valr
end if

! Interpolate advection
call adv%interp(mpl,nam,geom,samp)

! Test advection efficiency
call adv%test(mpl,rng,nam,geom,samp,ens)

end subroutine adv_compute

!----------------------------------------------------------------------
! Subroutine: adv_compute_raw
! Purpose: compute raw advection
!----------------------------------------------------------------------
subroutine adv_compute_raw(adv,mpl,nam,geom,samp,ens)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv  ! Advection
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(samp_type),intent(inout) :: samp ! Sampling
type(ens_type), intent(in) :: ens     ! Ensemble

! Local variables
integer :: ic0,ic2,ic2a,jn,jc0,il0,il0i,isub,iv,its,ie,ie_sub,ic0a,nc0d,ic0d,jc0d,jnmax,nnmax
integer :: nn(samp%nc2a,geom%nl0),ic0_rac(samp%nc2a,geom%nl0)
integer :: c0_to_c0d(geom%nc0),c0a_to_c0d(geom%nc0a)
integer,allocatable :: jc0_ra(:,:,:),c0d_to_c0(:)
real(kind_real) :: m2m2,fld_1,fld_2,cov,cor_avg_max
real(kind_real) :: mean,stddev,dist_sum,norm,norm_tot
real(kind_real) :: lon_c2(nam%nc2),lat_c2(nam%nc2)
real(kind_real),allocatable :: fld_ext_1(:,:),fld_ext_2(:,:)
real(kind_real),allocatable :: m2_1(:,:,:,:,:),m2_2(:,:,:,:,:)
real(kind_real),allocatable :: m11(:,:,:,:,:)
real(kind_real),allocatable :: cor(:,:),cor_avg(:)
logical :: valid_c2(nam%nc2)
logical,allocatable :: lcheck_c0d(:),mask_nn(:,:,:),mask_nn_tmp(:)
character(len=1024),parameter :: subr = 'adv_compute_raw'
type(com_type) :: com_AD
type(mesh_type) :: mesh

write(mpl%info,'(a7,a)') '','Compute raw advection'
call mpl%flush

! Initialization
do il0=1,geom%nl0
   do ic2a=1,samp%nc2a
      ic2 = samp%c2a_to_c2(ic2a)
      ic0 = samp%c2_to_c0(ic2)
      adv%lon_c2a(ic2a,il0) = geom%lon(ic0)
      adv%lat_c2a(ic2a,il0) = geom%lat(ic0)
   end do
end do
adv%cor_max = mpl%msv%valr

! Advection for timeslot 1
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      samp%adv_lon(ic0a,il0,1) = geom%lon(ic0)
      samp%adv_lat(ic0a,il0,1) = geom%lat(ic0)
   end do
end do

do its=2,nam%nts
   write(mpl%info,'(a10,a,i2)') '','Timeslot ',its
   call mpl%flush

   ! Find research area center
   write(mpl%info,'(a13,a)') '','Find research area center'
   call mpl%flush
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%mask_c2a(ic2a,il0)) then
            ! Research area center
            if (cor_tracker.and.(its>2)) then
               ! Find closer point to previous timeslot advected mesh point
               call geom%tree%find_nearest_neighbors(adv%lon_c2a_flt(ic2a,il0,its-1),adv%lat_c2a_flt(ic2a,il0,its-1),1, &
             & ic0_rac(ic2a:ic2a,il0))
            else
               ! Set at origin mesh point
               ic2 = samp%c2a_to_c2(ic2a)
               ic0_rac(ic2a,il0) = samp%c2_to_c0(ic2)
            end if
         end if
      end do
   end do

   ! Define search radius
   write(mpl%info,'(a13,a)') '','Define search radius'
   call mpl%flush
   do il0=1,geom%nl0
      il0i = min(il0,geom%nl0i)
      do ic0=1,geom%nc0
         adv%search_rad(ic0,il0) = min(nam%adv_rad,min(geom%mdist(ic0,il0i),geom%mesh%bdist(geom%mesh%order_inv(ic0))))
      end do
   end do

   ! Count nearest neighbors
   write(mpl%info,'(a13,a)') '','Count nearest neighbors:'
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc2a*geom%nl0)
   nn = mpl%msv%vali
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%mask_c2a(ic2a,il0)) then
            ! Index
            ic0 = ic0_rac(ic2a,il0)

            ! Count nearest neighbors
            call geom%tree%count_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),adv%search_rad(ic0,il0),nn(ic2a,il0))

            ! Update
            call mpl%prog_print((il0-1)*samp%nc2a+ic2a)
         end if
      end do
   end do
   call mpl%prog_final

   ! Get maximum
   nnmax = maxval(nn)

   ! Allocation
   allocate(jc0_ra(nnmax,samp%nc2a,geom%nl0))
   allocate(mask_nn(nnmax,samp%nc2a,geom%nl0))
   allocate(m2_1(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m2_2(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m11(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))

   ! Find nearest neighbors
   write(mpl%info,'(a13,a)') '','Find nearest neighbors:'
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc2a*geom%nl0)
   lcheck_c0d = samp%lcheck_c0a
   jc0_ra = mpl%msv%vali
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%mask_c2a(ic2a,il0)) then
            ! Index
            ic0 = ic0_rac(ic2a,il0)

            ! Find nearest neighbors
            call geom%tree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nn(ic2a,il0),jc0_ra(1:nn(ic2a,il0),ic2a,il0))

            ! Check points
            do jn=1,nn(ic2a,il0)
               jc0 = jc0_ra(jn,ic2a,il0)
               lcheck_c0d(jc0) = .true.
            end do

            ! Update
            call mpl%prog_print((il0-1)*samp%nc2a+ic2a)
         end if
      end do
   end do
   call mpl%prog_final

   ! Halo D size
   nc0d = count(lcheck_c0d)

   ! Allocation
   allocate(c0d_to_c0(nc0d))
   allocate(fld_ext_1(nc0d,geom%nl0))
   allocate(fld_ext_2(nc0d,geom%nl0))

   ! Global <-> local conversions for fields
   c0_to_c0d = mpl%msv%vali
   ic0d = 0
   do ic0=1,geom%nc0
      if (lcheck_c0d(ic0)) then
         ic0d = ic0d+1
         c0d_to_c0(ic0d) = ic0
         c0_to_c0d(ic0) = ic0d
      end if
   end do

   ! Inter-halo conversions
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      ic0d = c0_to_c0d(ic0)
      c0a_to_c0d(ic0a) = ic0d
   end do

   ! Setup communications
   call com_AD%setup(mpl,'com_AD',geom%nc0,geom%nc0a,nc0d,c0d_to_c0,c0a_to_c0d,geom%c0_to_proc,geom%c0_to_c0a)

   ! Mask
   mask_nn = .false.
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%mask_c2a(ic2a,il0)) then
            ic0 = ic0_rac(ic2a,il0)
            if (geom%mask_c0(ic0,il0)) then
               do jn=1,nn(ic2a,il0)
                  jc0 = jc0_ra(jn,ic2a,il0)
                  mask_nn(jn,ic2a,il0) = geom%mask_c0(jc0,il0)
               end do
            end if
         end if
      end do
   end do

   ! Initialization
   m11 = 0.0
   m2_1 = 0.0
   m2_2 = 0.0

   ! Compute moments
   write(mpl%info,'(a13,a)') '','Compute moments'
   call mpl%flush
   do isub=1,ens%nsub
      if (ens%nsub==1) then
         write(mpl%info,'(a16,a)') '','Full ensemble, member:'
         call mpl%flush(.false.)
      else
         write(mpl%info,'(a16,a,i4,a)') '','Sub-ensemble ',isub,', member:'
         call mpl%flush(.false.)
      end if

      ! Compute covariances
      do ie_sub=1,ens%ne/ens%nsub
         write(mpl%info,'(i4)') ie_sub
         call mpl%flush(.false.)

         ! Full ensemble index
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub

         do iv=1,nam%nv
            ! Halo extension
            call com_AD%ext(mpl,geom%nl0,ens%fld(:,:,iv,1,ie),fld_ext_1)
            call com_AD%ext(mpl,geom%nl0,ens%fld(:,:,iv,its,ie),fld_ext_2)

            !$omp parallel do schedule(static) private(il0,ic2a,jn,ic0,jc0,ic0d,jc0d,fld_1,fld_2)
            do il0=1,geom%nl0
               do ic2a=1,samp%nc2a
                  if (samp%mask_c2a(ic2a,il0)) then
                     do jn=1,nn(ic2a,il0)
                        if (mask_nn(jn,ic2a,il0)) then
                           ! Indices
                           ic0 = ic0_rac(ic2a,il0)
                           jc0 = jc0_ra(jn,ic2a,il0)

                           ! Halo D indices
                           ic0d = c0_to_c0d(ic0)
                           jc0d = c0_to_c0d(jc0)

                           ! Copy points
                           fld_1 = fld_ext_1(ic0d,il0)
                           fld_2 = fld_ext_2(jc0d,il0)

                           ! Covariance
                           m11(jn,ic2a,il0,iv,isub) = m11(jn,ic2a,il0,iv,isub)+fld_1*fld_2

                           ! Variances
                           m2_1(jn,ic2a,il0,iv,isub) = m2_1(jn,ic2a,il0,iv,isub)+fld_1**2
                           m2_2(jn,ic2a,il0,iv,isub) = m2_2(jn,ic2a,il0,iv,isub)+fld_2**2
                        end if
                     end do
                  end if
               end do
            end do
            !$omp end parallel do
         end do
      end do
      write(mpl%info,'(a)') ''
      call mpl%flush
   end do

   ! Find correlation propagation
   write(mpl%info,'(a13,a)') '','Find correlation propagation'
   call mpl%flush

   ! Initialization
   adv%lon_c2a_raw(:,:,its) = adv%lon_c2a
   adv%lat_c2a_raw(:,:,its) = adv%lat_c2a
   adv%dist_c2a_raw(:,:,its) = 0.0

   do il0=1,geom%nl0
      !$omp parallel do schedule(static) private(ic2a,jn,iv,cov,m2m2,jnmax,cor_avg_max,jc0,ic2,ic0,ic0a,mean,stddev), &
      !$omp&                             firstprivate(cor,cor_avg,mask_nn_tmp)
      do ic2a=1,samp%nc2a
         if (samp%mask_c2a(ic2a,il0)) then
            ! Allocation
            allocate(cor(nn(ic2a,il0),nam%nv))
            allocate(cor_avg(nn(ic2a,il0)))
            allocate(mask_nn_tmp(nn(ic2a,il0)))

            ! Index
            ic0 = ic0_rac(ic2a,il0)

            ! Initialization
            mask_nn_tmp = .false.
            jnmax = mpl%msv%vali
            cor_avg_max = 0.0

            do jn=1,nn(ic2a,il0)
               if (mask_nn(jn,ic2a,il0)) then
                  ! Compute covariance and correlation
                  do iv=1,nam%nv
                     ! Covariance
                     cov = sum(m11(jn,ic2a,il0,iv,:))/real(ens%nsub,kind_real)

                     ! Correlation
                     m2m2 = sum(m2_1(jn,ic2a,il0,iv,:))*sum(m2_2(jn,ic2a,il0,iv,:))/real(ens%nsub**2,kind_real)
                     if (m2m2>0.0) then
                        cor(jn,iv) = cov/sqrt(m2m2)
                     else
                        cor(jn,iv) = mpl%msv%valr
                     end if
                  end do

                  ! Average correlation
                  if (mpl%msv%isanynotr(cor(jn,:))) then
                     cor_avg(jn) = sum(cor(jn,:),mask=mpl%msv%isnotr(cor(jn,:)))/real(count(mpl%msv%isnotr(cor(jn,:))),kind_real)
                  else
                     cor_avg(jn) = mpl%msv%valr
                  end if

                  ! Set mask
                  mask_nn_tmp(jn) = (cor_avg(jn)>0.0).and.mpl%msv%isnotr(cor_avg(jn))
               end if

               ! Locate the maximum correlation
               if (mask_nn_tmp(jn)) then
                  if (sup(cor_avg(jn),cor_avg_max)) then
                     jnmax = jn
                     cor_avg_max = cor_avg(jn)
                  end if
               end if
            end do

            ! Save maximum correlation
            if (mpl%msv%isnoti(jnmax)) then
               ic0a = geom%c0_to_c0a(ic0)
               adv%cor_max(ic0a,il0,:,its) = cor(jnmax,:)
            end if

            ! Apply threshold on standard-deviation
            if (count(mask_nn_tmp)>2) then
               mean = sum(cor_avg,mask=mask_nn_tmp)/real(count(mask_nn_tmp),kind_real)
               stddev = sqrt(sum((cor_avg-mean)**2,mask=mask_nn_tmp)/real(count(mask_nn_tmp)-1,kind_real))
            else
               stddev = 0.0
            end if
            if (cor_avg_max<max_std_ratio*stddev) jnmax = mpl%msv%vali

            if (mpl%msv%isnoti(jnmax)) then
               ! Indices
               jc0 = jc0_ra(jnmax,ic2a,il0)

               ! Save advection
               adv%lon_c2a_raw(ic2a,il0,its) = geom%lon(jc0)
               adv%lat_c2a_raw(ic2a,il0,its) = geom%lat(jc0)

               ! Compute distance
               call sphere_dist(adv%lon_c2a(ic2a,il0),adv%lat_c2a(ic2a,il0),adv%lon_c2a_raw(ic2a,il0,its), &
             & adv%lat_c2a_raw(ic2a,il0,its),adv%dist_c2a_raw(ic2a,il0,its))
               if (adv%dist_c2a_raw(ic2a,il0,its)>adv%search_rad(ic0,il0)) call mpl%abort(subr,'search radius not respected')
            end if

            ! Release memory
            deallocate(cor)
            deallocate(cor_avg)
            deallocate(mask_nn_tmp)
         end if
      end do
      !$omp end parallel do

      ! Check raw mesh
      call mpl%loc_to_glb(samp%nc2a,adv%lon_c2a_raw(:,il0,its),nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lon_c2)
      call mpl%loc_to_glb(samp%nc2a,adv%lat_c2a_raw(:,il0,its),nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lat_c2)
      if (mpl%main) then
         call mesh%copy(samp%mesh)
         call mesh%store(mpl,lon_c2,lat_c2)
         call mesh%check(mpl,valid_c2)
      end if
      call mpl%f_comm%broadcast(valid_c2,mpl%ioproc-1)
      adv%valid_raw(il0,its) = real(count(valid_c2.and.samp%mesh%valid.and.samp%mask_c2(:,il0)),kind_real) &
                           & /real(count(samp%mesh%valid.and.samp%mask_c2(:,il0)),kind_real)

      ! Average distance
      dist_sum = sum(adv%dist_c2a_raw(:,il0,its),mask=mpl%msv%isnotr(adv%dist_c2a_raw(:,il0,its)))
      norm = real(count(mpl%msv%isnotr(adv%dist_c2a_raw(:,il0,its))),kind_real)
      call mpl%f_comm%allreduce(dist_sum,adv%dist_raw(il0,its),fckit_mpi_sum())
      call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
      adv%dist_raw(il0,its) = adv%dist_raw(il0,its)/norm_tot
   end do

   ! Release memory
   deallocate(jc0_ra)
   deallocate(mask_nn)
   deallocate(m2_1)
   deallocate(m2_2)
   deallocate(m11)
   deallocate(c0d_to_c0)
   call com_AD%dealloc
   deallocate(fld_ext_1)
   deallocate(fld_ext_2)
end do

end subroutine adv_compute_raw

!----------------------------------------------------------------------
! Subroutine: adv_filter
! Purpose: filter advection
!----------------------------------------------------------------------
subroutine adv_filter(adv,mpl,nam,geom,samp)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv  ! Advection
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(samp_type),intent(inout) :: samp ! Sampling

! Local variables
integer :: ic0,ic2,ic2a,il0,il0i,its,iter
real(kind_real) :: dist_sum,norm,norm_tot,valid_flt,dist_flt,rhflt,drhflt
real(kind_real) :: lon_c2a(samp%nc2a),lat_c2a(samp%nc2a),dist_c2a(samp%nc2a)
real(kind_real) :: x_ori(samp%nc2a),y_ori(samp%nc2a),z_ori(samp%nc2a)
real(kind_real) :: dx_ini(samp%nc2a),dy_ini(samp%nc2a),dz_ini(samp%nc2a)
real(kind_real) :: dx(samp%nc2a),dy(samp%nc2a),dz(samp%nc2a),dd(samp%nc2a)
real(kind_real) :: lon_c2(nam%nc2),lat_c2(nam%nc2)
logical :: dichotomy,convergence
logical :: valid_c2(nam%nc2)
character(len=1024),parameter :: subr = 'adv_filter'
type(mesh_type) :: mesh

write(mpl%info,'(a7,a)') '','Filter advection'
call mpl%flush

! Initialization
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic0 = samp%c2_to_c0(ic2)
   call lonlat2xyz(mpl,geom%lon(ic0),geom%lat(ic0),x_ori(ic2a),y_ori(ic2a),z_ori(ic2a))
end do

do its=2,nam%nts
   write(mpl%info,'(a10,a,i2)') '','Timeslot ',its
   call mpl%flush

   do il0=1,geom%nl0
      write(mpl%info,'(a13,a,i3)') '','Level ',nam%levs(il0)
      call mpl%flush

      ! Convert to cartesian coordinates
      do ic2a=1,samp%nc2a
         if (samp%mask_c2a(ic2a,il0)) then
            call lonlat2xyz(mpl,adv%lon_c2a_raw(ic2a,il0,its),adv%lat_c2a_raw(ic2a,il0,its),dx_ini(ic2a),dy_ini(ic2a),dz_ini(ic2a))
            if (mpl%msv%isnotr(dx_ini(ic2a))) dx_ini(ic2a) = dx_ini(ic2a)-x_ori(ic2a)
            if (mpl%msv%isnotr(dy_ini(ic2a))) dy_ini(ic2a) = dy_ini(ic2a)-y_ori(ic2a)
            if (mpl%msv%isnotr(dz_ini(ic2a))) dz_ini(ic2a) = dz_ini(ic2a)-z_ori(ic2a)
         else
            dx_ini(ic2a) = mpl%msv%valr
            dy_ini(ic2a) = mpl%msv%valr
            dz_ini(ic2a) = mpl%msv%valr
         end if
      end do

      ! Dichotomy initialization
      convergence = .true.
      dichotomy = .false.
      rhflt = nam%adv_rhflt
      drhflt = rhflt

      do iter=1,nam%adv_niter
         ! Copy advection
         dx = dx_ini
         dy = dy_ini
         dz = dz_ini

         ! Median filter to remove extreme values
         dd = dx**2+dy**2+dz**2
         call samp%diag_filter(mpl,nam,'median',rhflt,dx,dd)
         call samp%diag_filter(mpl,nam,'median',rhflt,dy,dd)
         call samp%diag_filter(mpl,nam,'median',rhflt,dz,dd)

         ! Average filter to smooth values
         call samp%diag_filter(mpl,nam,'gc99',rhflt,dx)
         call samp%diag_filter(mpl,nam,'gc99',rhflt,dy)
         call samp%diag_filter(mpl,nam,'gc99',rhflt,dz)

         ! Back to spherical coordinates
         dx = x_ori+dx
         dy = y_ori+dy
         dz = z_ori+dz
         do ic2a=1,samp%nc2a
            if (samp%mask_c2a(ic2a,il0)) then
               call xyz2lonlat(mpl,dx(ic2a),dy(ic2a),dz(ic2a),lon_c2a(ic2a),lat_c2a(ic2a))
            else
               lon_c2a(ic2a) = adv%lon_c2a(ic2a,il0)
               lat_c2a(ic2a) = adv%lat_c2a(ic2a,il0)
            end if
         end do

         ! Reduce distance with respect to the boundary
         il0i = min(il0,geom%nl0i)
         do ic2a=1,samp%nc2a
            ic2 = samp%c2a_to_c2(ic2a)
            ic0 = samp%c2_to_c0(ic2)
            call reduce_arc(adv%lon_c2a(ic2a,il0),adv%lat_c2a(ic2a,il0),lon_c2a(ic2a),lat_c2a(ic2a), &
          & adv%search_rad(ic0,il0i),dist_c2a(ic2a))
         end do

         ! Check mesh
         call mpl%loc_to_glb(samp%nc2a,lon_c2a,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lon_c2)
         call mpl%loc_to_glb(samp%nc2a,lat_c2a,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lat_c2)
         if (mpl%main) then
            call mesh%copy(samp%mesh)
            call mesh%store(mpl,lon_c2,lat_c2)
            call mesh%check(mpl,valid_c2)
         end if
         call mpl%f_comm%broadcast(valid_c2,mpl%ioproc-1)
         valid_flt = real(count(valid_c2.and.samp%mesh%valid.and.samp%mask_c2(:,il0)),kind_real) &
                   & /real(count(samp%mesh%valid.and.samp%mask_c2(:,il0)),kind_real)

         ! Average distance
         dist_sum = sum(dist_c2a,mask=mpl%msv%isnotr(dist_c2a))
         norm = real(count(mpl%msv%isnotr(dist_c2a)),kind_real)
         call mpl%f_comm%allreduce(dist_sum,dist_flt,fckit_mpi_sum())
         call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
         dist_flt = dist_flt/norm_tot

         ! Print result
         write(mpl%info,'(a16,a,i2,a,f10.2,a,f7.3,a)') '','Iteration ',iter,': rhflt = ', &
       & rhflt*reqkm,' km, validity = ',100.0*valid_flt,'%'
         call mpl%flush

         ! Update support radius
         if (inf(valid_flt,1.0_kind_real)) then
            ! Increase filtering support radius
            if (dichotomy) then
               drhflt = 0.5*drhflt
               if (iter<nam%adv_niter) rhflt = rhflt+drhflt
            else
               convergence = .false.
               if (iter<nam%adv_niter) rhflt = rhflt+drhflt
               drhflt = 2.0*drhflt
            end if
         else
            ! Convergence
            convergence = .true.

            ! Change dichotomy status
            if (.not.dichotomy) then
               dichotomy = .true.
               drhflt = 0.5*drhflt
            end if

            ! Decrease filtering support radius
            drhflt = 0.5*drhflt
            if (iter<nam%adv_niter) rhflt = rhflt-drhflt

            ! Save values
            adv%lon_c2a_flt(:,il0,its) = lon_c2a
            adv%lat_c2a_flt(:,il0,its) = lat_c2a
            adv%dist_c2a_flt(:,il0,its) = dist_c2a
            adv%valid_flt(il0,its) = valid_flt
            adv%dist_flt(il0,its) = dist_flt
            adv%rhflt(il0,its) = rhflt
         end if
      end do

      ! Check convergence
      if (.not.convergence) call mpl%abort(subr,'convergence failed in adv%compute')

      ! Print results
      write(mpl%info,'(a16,a,f10.2,a)') '','Optimal filtering support radius: ',adv%rhflt(il0,its)*reqkm,' km'
      call mpl%flush
      write(mpl%info,'(a16,a,f5.1,a,f5.1,a)') '','Valid points:             ',100.0*adv%valid_raw(il0,its),'%  ~> ', &
    & 100.0*adv%valid_flt(il0,its),'%'
      call mpl%flush
      write(mpl%info,'(a16,a,f10.1,a,f10.1,a)') '','Advection distance: ',adv%dist_raw(il0,its)*reqkm,' km ~> ', &
    & adv%dist_flt(il0,its)*reqkm,' km'
      call mpl%flush
   end do
end do

end subroutine adv_filter

!----------------------------------------------------------------------
! Subroutine: adv_interp
! Purpose: interpolate advection
!----------------------------------------------------------------------
subroutine adv_interp(adv,mpl,nam,geom,samp)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv  ! Advection
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(samp_type),intent(inout) :: samp ! Sampling

! Local variables
integer :: ic0,ic2,ic2a,il0,il0i,its,ic0a
real(kind_real) :: reduced_dist
real(kind_real) :: x_ori(samp%nc2a),y_ori(samp%nc2a),z_ori(samp%nc2a)
real(kind_real) :: dx(samp%nc2a),dy(samp%nc2a),dz(samp%nc2a)
real(kind_real) :: dx_c2b(samp%nc2b),dy_c2b(samp%nc2b),dz_c2b(samp%nc2b)
real(kind_real) :: dx_c0a(geom%nc0a),dy_c0a(geom%nc0a),dz_c0a(geom%nc0a)
real(kind_real) :: x_ori_c0a(geom%nc0a),y_ori_c0a(geom%nc0a),z_ori_c0a(geom%nc0a)

write(mpl%info,'(a7,a)') '','Interpolate advection'
call mpl%flush

! Initialization
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic0 = samp%c2_to_c0(ic2)
   call lonlat2xyz(mpl,geom%lon(ic0),geom%lat(ic0),x_ori(ic2a),y_ori(ic2a),z_ori(ic2a))
end do
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   call lonlat2xyz(mpl,geom%lon(ic0),geom%lat(ic0),x_ori_c0a(ic0a),y_ori_c0a(ic0a),z_ori_c0a(ic0a))
end do

do its=2,nam%nts
   do il0=1,geom%nl0
      ! Convert to cartesian coordinates
      do ic2a=1,samp%nc2a
         if (samp%mask_c2a(ic2a,il0)) then
            call lonlat2xyz(mpl,adv%lon_c2a_flt(ic2a,il0,its),adv%lat_c2a_flt(ic2a,il0,its),dx(ic2a),dy(ic2a),dz(ic2a))
            if (mpl%msv%isnotr(dx(ic2a))) dx(ic2a) = dx(ic2a)-x_ori(ic2a)
            if (mpl%msv%isnotr(dy(ic2a))) dy(ic2a) = dy(ic2a)-y_ori(ic2a)
            if (mpl%msv%isnotr(dz(ic2a))) dz(ic2a) = dz(ic2a)-z_ori(ic2a)
         else
            dx(ic2a) = 0.0
            dy(ic2a) = 0.0
            dz(ic2a) = 0.0
         end if
      end do

      ! Advection interpolation
      il0i = min(il0,geom%nl0i)
      call samp%com_AB%ext(mpl,dx,dx_c2b)
      call samp%com_AB%ext(mpl,dy,dy_c2b)
      call samp%com_AB%ext(mpl,dz,dz_c2b)
      call samp%h(il0i)%apply(mpl,dx_c2b,dx_c0a)
      call samp%h(il0i)%apply(mpl,dy_c2b,dy_c0a)
      call samp%h(il0i)%apply(mpl,dz_c2b,dz_c0a)

      ! Back to spherical coordinates
      dx_c0a = x_ori_c0a+dx_c0a
      dy_c0a = y_ori_c0a+dy_c0a
      dz_c0a = z_ori_c0a+dz_c0a
      do ic0a=1,geom%nc0a
         call xyz2lonlat(mpl,dx_c0a(ic0a),dy_c0a(ic0a),dz_c0a(ic0a),samp%adv_lon(ic0a,il0,its),samp%adv_lat(ic0a,il0,its))
      end do

      ! Reduce distance with respect to the boundary
      il0i = min(il0,geom%nl0i)
      do ic0a=1,geom%nc0a
         if (geom%mask_c0a(ic0a,il0)) then
            ic0 = geom%c0a_to_c0(ic0a)
            call reduce_arc(geom%lon(ic0),geom%lat(ic0),samp%adv_lon(ic0a,il0,its),samp%adv_lat(ic0a,il0,its), &
          & adv%search_rad(ic0,il0i),reduced_dist)
         end if
      end do
   end do
end do

end subroutine adv_interp

!----------------------------------------------------------------------
! Subroutine: adv_test
! Purpose: test advection efficiency
!----------------------------------------------------------------------
subroutine adv_test(adv,mpl,rng,nam,geom,samp,ens)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv  ! Advection
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(rng_type),intent(inout) :: rng   ! Random number generator
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(samp_type),intent(inout) :: samp ! Sampling
type(ens_type), intent(in) :: ens     ! Ensemble

! Local variables
integer :: its,il0,ic0,ic0a,i_s,i_s_loc,iproc,jc0,nc0d,ie,iv
integer :: ic0d,d_n_s_max,d_n_s_max_loc
integer :: c0_to_c0d(geom%nc0),c0a_to_c0d(geom%nc0a)
integer,allocatable :: c0d_to_c0(:),interpd_lg(:,:,:)
real(kind_real) :: score_loc,score_max,score_adv,norm
real(kind_real) :: score_loc_tot,score_max_tot,score_adv_tot,norm_tot
real(kind_real) :: adv_lon(geom%nc0,geom%nl0),adv_lat(geom%nc0,geom%nl0)
real(kind_real) :: m2_1(geom%nc0a,geom%nl0),m2_2_loc(geom%nc0a,geom%nl0),m11_loc(geom%nc0a,geom%nl0)
real(kind_real) :: m2_2_adv(geom%nc0a,geom%nl0),m11_adv(geom%nc0a,geom%nl0)
real(kind_real) :: fld(geom%nc0a,geom%nl0)
real(kind_real),allocatable :: fld_d(:,:)
logical :: mask_c0(geom%nc0),lcheck_c0d(geom%nc0)
logical,allocatable :: lcheck_d(:,:,:)
type(com_type) :: com_ADinv
type(linop_type) :: dinvfull(geom%nl0,2:nam%nts),dinv(geom%nl0,2:nam%nts)

write(mpl%info,'(a7,a)') '','Test advection'
call mpl%flush

! Initialization
mask_c0 = .true.

do its=2,nam%nts
   ! Local to global
   call mpl%loc_to_glb(geom%nl0,geom%nc0a,samp%adv_lon(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,adv_lon)
   call mpl%loc_to_glb(geom%nl0,geom%nc0a,samp%adv_lat(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,adv_lat)

   do il0=1,geom%nl0
      ! Compute direct interpolation
      call dinvfull(il0,its)%interp(mpl,rng,geom%nc0,geom%lon,geom%lat,mask_c0,geom%nc0,adv_lon(:,il0),adv_lat(:,il0),mask_c0, &
    & nam%diag_interp)
   end do
end do

! Allocation
d_n_s_max = 0
do its=2,nam%nts
   do il0=1,geom%nl0
      d_n_s_max = max(d_n_s_max,dinvfull(il0,its)%n_s)
   end do
end do
allocate(lcheck_d(d_n_s_max,geom%nl0,2:nam%nts))

! Halo definitions

! Halo A
lcheck_c0d = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   lcheck_c0d(ic0) = .true.
end do

! Halo D
lcheck_d = .false.
do its=2,nam%nts
   do il0=1,geom%nl0
      do i_s=1,dinvfull(il0,its)%n_s
         ic0 = dinvfull(il0,its)%row(i_s)
         iproc = geom%c0_to_proc(ic0)
         if (iproc==mpl%myproc) then
            jc0 = dinvfull(il0,its)%col(i_s)
            lcheck_d(i_s,il0,its) = .true.
            lcheck_c0d(jc0) = .true.
         end if
      end do
   end do
end do

! Halo sizes
do its=2,nam%nts
   do il0=1,geom%nl0
      dinv(il0,its)%n_s = count(lcheck_d(:,il0,its))
   end do
end do
nc0d = count(lcheck_c0d)

! Global <-> local conversions for fields
allocate(c0d_to_c0(nc0d))
c0_to_c0d = mpl%msv%vali
ic0d = 0
do ic0=1,geom%nc0
   if (lcheck_c0d(ic0)) then
      ic0d = ic0d+1
      c0d_to_c0(ic0d) = ic0
      c0_to_c0d(ic0) = ic0d
   end if
end do

! Inter-halo conversions
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0d = c0_to_c0d(ic0)
   c0a_to_c0d(ic0a) = ic0d
end do

! Global <-> local conversions for data
d_n_s_max_loc = 0
do its=2,nam%nts
   do il0=1,geom%nl0
      d_n_s_max_loc = max(d_n_s_max_loc,dinv(il0,its)%n_s)
   end do
end do
allocate(interpd_lg(d_n_s_max_loc,geom%nl0,2:nam%nts))
do its=2,nam%nts
   do il0=1,geom%nl0
      i_s_loc = 0
      do i_s=1,dinvfull(il0,its)%n_s
         if (lcheck_d(i_s,il0,its)) then
            i_s_loc = i_s_loc+1
            interpd_lg(i_s_loc,il0,its) = i_s
         end if
      end do
   end do
end do

! Local data
do its=2,nam%nts
   do il0=1,geom%nl0
      write(dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
      dinv(il0,its)%n_src = nc0d
      dinv(il0,its)%n_dst = geom%nc0a
      call dinv(il0,its)%alloc
      do i_s_loc=1,dinv(il0,its)%n_s
         i_s = interpd_lg(i_s_loc,il0,its)
         dinv(il0,its)%row(i_s_loc) = geom%c0_to_c0a(dinvfull(il0,its)%row(i_s))
         dinv(il0,its)%col(i_s_loc) = c0_to_c0d(dinvfull(il0,its)%col(i_s))
         dinv(il0,its)%S(i_s_loc) = dinvfull(il0,its)%S(i_s)
      end do
      call dinv(il0,its)%reorder(mpl)
   end do
end do

! Setup communications
call com_ADinv%setup(mpl,'com_ADinv',geom%nc0,geom%nc0a,nc0d,c0d_to_c0,c0a_to_c0d,geom%c0_to_proc,geom%c0_to_c0a)

! Allocation
allocate(fld_d(nc0d,geom%nl0))

! Compute correlation at zero separation
do its=2,nam%nts
   do iv=1,nam%nv
      ! Initialization
      m2_1 = 0.0
      m2_2_loc = 0.0
      m2_2_adv = 0.0
      m11_loc = 0.0
      m11_adv = 0.0

      do ie=1,ens%ne
         ! Halo extension from zone A to zone D
         call com_ADinv%ext(mpl,geom%nl0,ens%fld(:,:,iv,its,ie),fld_d)

         ! Interpolation
         !$omp parallel do schedule(static) private(il0)
         do il0=1,geom%nl0
            call dinv(il0,its)%apply(mpl,fld_d(:,il0),fld(:,il0),msdst=.false.)
         end do
         !$omp end parallel do

         ! Compute variance/covariance
         m2_1 = m2_1+ens%fld(:,:,iv,1,ie)**2
         m2_2_loc = m2_2_loc+ens%fld(:,:,iv,its,ie)**2
         m2_2_adv = m2_2_adv+fld**2
         m11_loc = m11_loc+ens%fld(:,:,iv,1,ie)*ens%fld(:,:,iv,its,ie)
         m11_adv = m11_adv+ens%fld(:,:,iv,1,ie)*fld
      end do

      ! Compute correlation
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if ((m2_1(ic0a,il0)>0.0).and.(m2_2_loc(ic0a,il0)>0.0)) then
               adv%cor_loc(ic0a,il0,iv,its) = m11_loc(ic0a,il0)/sqrt(m2_1(ic0a,il0)*m2_2_loc(ic0a,il0))
            else
               adv%cor_loc(ic0a,il0,iv,its) = mpl%msv%valr
            end if
            if ((m2_1(ic0a,il0)>0.0).and.(m2_2_adv(ic0a,il0)>0.0)) then
               adv%cor_adv(ic0a,il0,iv,its) = m11_adv(ic0a,il0)/sqrt(m2_1(ic0a,il0)*m2_2_adv(ic0a,il0))
            else
               adv%cor_adv(ic0a,il0,iv,its) = mpl%msv%valr
            end if
         end do
      end do
   end do
end do

! Print results
do its=2,nam%nts
   score_loc = sum(adv%cor_loc(:,:,:,its),mask=mpl%msv%isnotr(adv%cor_max(:,:,:,its)))
   score_max = sum(adv%cor_max(:,:,:,its),mask=mpl%msv%isnotr(adv%cor_max(:,:,:,its)))
   score_adv = sum(adv%cor_adv(:,:,:,its),mask=mpl%msv%isnotr(adv%cor_max(:,:,:,its)))
   norm = real(count(mpl%msv%isnotr(adv%cor_max(:,:,:,its))),kind_real)
   call mpl%f_comm%allreduce(score_loc,score_loc_tot,fckit_mpi_sum())
   call mpl%f_comm%allreduce(score_max,score_max_tot,fckit_mpi_sum())
   call mpl%f_comm%allreduce(score_adv,score_adv_tot,fckit_mpi_sum())
   call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
   score_loc = score_loc/norm
   score_max = score_max/norm
   score_adv = score_adv/norm
   write(mpl%info,'(a10,a,i2,a,f6.3,a,f6.3,a,f6.3)') '','Timeslot ',its,': ',score_loc,' ~> ',score_adv,' < ',score_max
   call mpl%flush
end do

end subroutine adv_test

!----------------------------------------------------------------------
! Subroutine: adv_write
! Purpose: write advection
!----------------------------------------------------------------------
subroutine adv_write(adv,mpl,nam,geom,bpar,io,samp)

implicit none

! Passed variables
class(adv_type),intent(in) :: adv   ! Advection
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(bpar_type),intent(in) :: bpar  ! Block parameters
type(io_type),intent(in) :: io      ! I/O
type(samp_type),intent(in) :: samp  ! Sampling

! Local variables
integer :: ncid,nc2_id,nl0_id,nts_id,na_id,vunit_id,larc_s_id,larc_e_id
integer :: lon_c2_id,lat_c2_id,lon_c2_raw_id,lat_c2_raw_id,dist_c2_raw_id,valid_raw_id,dist_raw_id
integer :: lon_c2_flt_id,lat_c2_flt_id,dist_c2_flt_id,valid_flt_id,dist_flt_id,rhflt_id
integer :: iproc,its,il0,ic2a,ic2,i,ib,iv,jv,jts
integer,allocatable :: c2a_to_c2(:)
real(kind_real),allocatable :: sbuf(:),rbuf(:),lon_c2(:,:),lat_c2(:,:)
real(kind_real),allocatable :: lon_c2_raw(:,:,:),lat_c2_raw(:,:,:),dist_c2_raw(:,:,:)
real(kind_real),allocatable :: lon_c2_flt(:,:,:),lat_c2_flt(:,:,:),dist_c2_flt(:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'adv_write'
type(fckit_mpi_status) :: status

! Allocation
allocate(sbuf(samp%nc2a*geom%nl0*(2+(nam%nts-1)*6)))

! Prepare buffer
i = 1
do il0=1,geom%nl0
   do ic2a=1,samp%nc2a
      sbuf(i) = adv%lon_c2a(ic2a,il0)*rad2deg
      i = i+1
      sbuf(i) = adv%lat_c2a(ic2a,il0)*rad2deg
      i = i+1
      do its=2,nam%nts
         sbuf(i) = adv%lon_c2a_raw(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = adv%lat_c2a_raw(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = adv%dist_c2a_raw(ic2a,il0,its)*reqkm
         i = i+1
         sbuf(i) = adv%lon_c2a_flt(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = adv%lat_c2a_flt(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = adv%dist_c2a_flt(ic2a,il0,its)*reqkm
         i = i+1
      end do
   end do
end do

if (mpl%main) then
   ! Allocation
   allocate(lon_c2(nam%nc2,geom%nl0))
   allocate(lat_c2(nam%nc2,geom%nl0))
   allocate(lon_c2_raw(nam%nc2,geom%nl0,nam%nts-1))
   allocate(lat_c2_raw(nam%nc2,geom%nl0,nam%nts-1))
   allocate(dist_c2_raw(nam%nc2,geom%nl0,nam%nts-1))
   allocate(lon_c2_flt(nam%nc2,geom%nl0,nam%nts-1))
   allocate(lat_c2_flt(nam%nc2,geom%nl0,nam%nts-1))
   allocate(dist_c2_flt(nam%nc2,geom%nl0,nam%nts-1))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(c2a_to_c2(samp%proc_to_nc2a(iproc)))
      allocate(rbuf(samp%proc_to_nc2a(iproc)*geom%nl0*(2+(nam%nts-1)*6)))

      if (iproc==mpl%ioproc) then
         ! Copy buffer
         c2a_to_c2 = samp%c2a_to_c2
         rbuf = sbuf
      else
         ! Receive data on ioproc
         call mpl%f_comm%receive(c2a_to_c2,iproc-1,mpl%tag,status)
         call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag+1,status)
      end if

      ! Write data
      i = 1
      do il0=1,geom%nl0
         do ic2a=1,samp%proc_to_nc2a(iproc)
            ic2 = c2a_to_c2(ic2a)
            lon_c2(ic2,il0) = rbuf(i)
            i = i+1
            lat_c2(ic2,il0) = rbuf(i)
            i = i+1
            do its=2,nam%nts
               lon_c2_raw(ic2,il0,its-1) = rbuf(i)
               i = i+1
               lat_c2_raw(ic2,il0,its-1) = rbuf(i)
               i = i+1
               dist_c2_raw(ic2,il0,its-1) = rbuf(i)
               i = i+1
               lon_c2_flt(ic2,il0,its-1) = rbuf(i)
               i = i+1
               lat_c2_flt(ic2,il0,its-1) = rbuf(i)
               i = i+1
               dist_c2_flt(ic2,il0,its-1) = rbuf(i)
               i = i+1
            end do
         end do
      end do

      ! Release memory
      deallocate(c2a_to_c2)
      deallocate(rbuf)
   end do
else
   ! Send data to ioproc
   call mpl%f_comm%send(samp%c2a_to_c2,mpl%ioproc-1,mpl%tag)
   call mpl%f_comm%send(sbuf,mpl%ioproc-1,mpl%tag+1)
end if
call mpl%update_tag(2)

! Release memory
deallocate(sbuf)

! Define filename
filename = trim(nam%prefix)//'_adv_diag'

if (mpl%main) then
   ! Create file
   call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call nam%write(mpl,ncid)

   ! Define dimensions
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2',nam%nc2,nc2_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nts',nam%nts-1,nts_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'na',samp%mesh%na,na_id))

   ! Define variables
   call mpl%ncerr(subr,nf90_def_var(ncid,'vunit',nc_kind_real,(/nc2_id,nl0_id/),vunit_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_c2',nc_kind_real,(/nc2_id,nl0_id/),lon_c2_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_c2_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_c2',nc_kind_real,(/nc2_id,nl0_id/),lat_c2_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_c2_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_c2_raw',nc_kind_real,(/nc2_id,nl0_id,nts_id/),lon_c2_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_c2_raw_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_c2_raw',nc_kind_real,(/nc2_id,nl0_id,nts_id/),lat_c2_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_c2_raw_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'dist_c2_raw',nc_kind_real,(/nc2_id,nl0_id,nts_id/),dist_c2_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,dist_c2_raw_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'valid_raw',nc_kind_real,(/nl0_id,nts_id/),valid_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,valid_raw_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'dist_raw',nc_kind_real,(/nl0_id,nts_id/),dist_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,dist_raw_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_c2_flt',nc_kind_real,(/nc2_id,nl0_id,nts_id/),lon_c2_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_c2_flt_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_c2_flt',nc_kind_real,(/nc2_id,nl0_id,nts_id/),lat_c2_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_c2_flt_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'dist_c2_flt',nc_kind_real,(/nc2_id,nl0_id,nts_id/),dist_c2_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,dist_c2_flt_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'valid_flt',nc_kind_real,(/nl0_id,nts_id/),valid_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,valid_flt_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'dist_flt',nc_kind_real,(/nl0_id,nts_id/),dist_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,dist_flt_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'rhflt',nc_kind_real,(/nl0_id,nts_id/),rhflt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,rhflt_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'larc_s',nf90_int,(/na_id/),larc_s_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'larc_e',nf90_int,(/na_id/),larc_e_id))

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write data
   call mpl%ncerr(subr,nf90_put_var(ncid,vunit_id,geom%vunit_c0(samp%c2_to_c0,:)))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_id,lon_c2))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_id,lat_c2))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_raw_id,lon_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_raw_id,lat_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_c2_raw_id,dist_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,valid_raw_id,adv%valid_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_raw_id,adv%dist_raw*reqkm))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_flt_id,lon_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_flt_id,lat_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_c2_flt_id,dist_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,valid_flt_id,adv%valid_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_flt_id,adv%dist_flt*reqkm))
   call mpl%ncerr(subr,nf90_put_var(ncid,rhflt_id,adv%rhflt*reqkm))
   call mpl%ncerr(subr,nf90_put_var(ncid,larc_s_id,samp%mesh%order(samp%mesh%larc(1,:))))
   call mpl%ncerr(subr,nf90_put_var(ncid,larc_e_id,samp%mesh%order(samp%mesh%larc(2,:))))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(lon_c2)
   deallocate(lat_c2)
   deallocate(lon_c2_raw)
   deallocate(lat_c2_raw)
   deallocate(dist_c2_raw)
   deallocate(lon_c2_flt)
   deallocate(lat_c2_flt)
   deallocate(dist_c2_flt)
end if

! Write test
do ib=1,bpar%nb
   iv = bpar%b_to_v1(ib)
   jv = bpar%b_to_v2(ib)
   its = bpar%b_to_ts1(ib)
   jts = bpar%b_to_ts2(ib)
   if ((iv>1).and.(iv==jv).and.(its==jts)) then
      call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_cor_loc',adv%cor_loc(:,:,iv,its))
      call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_cor_max',adv%cor_adv(:,:,iv,its))
      call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_cor_adv',adv%cor_max(:,:,iv,its))
   end if
end do

end subroutine adv_write

end module type_adv
