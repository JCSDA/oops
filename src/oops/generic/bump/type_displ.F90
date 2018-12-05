!----------------------------------------------------------------------
! Module: type_displ
! Purpose: displacement data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_displ

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_const, only: req,reqkm,rad2deg,deg2rad,msvalr
use tools_func, only: lonlatmod,sphere_dist,reduce_arc,vector_product
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi,isnotmsr,isallnotmsr,isanynotmsr
use tools_nc, only: ncfloat
use tools_qsort, only: qsort
use tools_stripack, only: trans,scoord
use type_com, only: com_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

logical,parameter :: cor_tracker = .true.                  ! Tracker method
real(kind_real),parameter :: max_std_ratio = 2.5_kind_real ! Minimum ratio between correlation maximum and local standard deviation

! Displacement data derived type
type displ_type
   integer :: niter                                   ! Number of stored iterations
   real(kind_real),allocatable :: dist(:,:,:)         ! Displacement distance
   real(kind_real),allocatable :: valid(:,:,:)        ! Displacement validity
   real(kind_real),allocatable :: rhflt(:,:,:)        ! Displacement filtering support radius

   real(kind_real),allocatable :: lon_c2a(:,:)        ! Longitude origin
   real(kind_real),allocatable :: lat_c2a(:,:)        ! Latitude origin
   real(kind_real),allocatable :: lon_c2a_raw(:,:,:)  ! Raw displaced longitude
   real(kind_real),allocatable :: lat_c2a_raw(:,:,:)  ! Raw displaced latitude
   real(kind_real),allocatable :: dist_c2a_raw(:,:,:) ! Raw displacement distance
   real(kind_real),allocatable :: lon_c2a_flt(:,:,:)  ! Filtered displaced longitude
   real(kind_real),allocatable :: lat_c2a_flt(:,:,:)  ! Filtered displaced latitude
   real(kind_real),allocatable :: dist_c2a_flt(:,:,:) ! Displacement distance, filtered
contains
   procedure :: alloc => displ_alloc
   procedure :: dealloc => displ_dealloc
   procedure :: compute => displ_compute
   procedure :: write => displ_write
end type displ_type

private
public :: displ_type

contains

!----------------------------------------------------------------------
! Subroutine: displ_alloc
! Purpose: displacement data allocation
!----------------------------------------------------------------------
subroutine displ_alloc(displ,nam,geom,samp)

implicit none

! Passed variables
class(displ_type),intent(inout) :: displ ! Displacement data
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry
type(samp_type),intent(in) :: samp       ! Sampling

! Allocation
allocate(displ%dist(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%valid(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%rhflt(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%lon_c2a(samp%nc2a,geom%nl0))
allocate(displ%lat_c2a(samp%nc2a,geom%nl0))
allocate(displ%lon_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(displ%lat_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(displ%dist_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(displ%lon_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(displ%lat_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(displ%dist_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))

! Initialization
call msr(displ%dist)
call msr(displ%valid)
call msr(displ%rhflt)
call msr(displ%lon_c2a)
call msr(displ%lat_c2a)
call msr(displ%lon_c2a_raw)
call msr(displ%lat_c2a_raw)
call msr(displ%dist_c2a_raw)
call msr(displ%lon_c2a_flt)
call msr(displ%lat_c2a_flt)
call msr(displ%dist_c2a_flt)

end subroutine displ_alloc

!----------------------------------------------------------------------
! Subroutine: displ_dealloc
! Purpose: displacement data deallocation
!----------------------------------------------------------------------
subroutine displ_dealloc(displ)

implicit none

! Passed variables
class(displ_type),intent(inout) :: displ ! Displacement data

! Deallocation
if (allocated(displ%dist)) deallocate(displ%dist)
if (allocated(displ%valid)) deallocate(displ%valid)
if (allocated(displ%rhflt)) deallocate(displ%rhflt)
if (allocated(displ%lon_c2a)) deallocate(displ%lon_c2a)
if (allocated(displ%lat_c2a)) deallocate(displ%lat_c2a)
if (allocated(displ%lon_c2a_raw)) deallocate(displ%lon_c2a_raw)
if (allocated(displ%lat_c2a_raw)) deallocate(displ%lat_c2a_raw)
if (allocated(displ%dist_c2a_raw)) deallocate(displ%dist_c2a_raw)
if (allocated(displ%lon_c2a_flt)) deallocate(displ%lon_c2a_flt)
if (allocated(displ%lat_c2a_flt)) deallocate(displ%lat_c2a_flt)
if (allocated(displ%dist_c2a_flt)) deallocate(displ%dist_c2a_flt)

end subroutine displ_dealloc

!----------------------------------------------------------------------
! Subroutine: displ_compute
! Purpose: compute correlation maximum displacement
!----------------------------------------------------------------------
subroutine displ_compute(displ,mpl,nam,geom,samp,ens)

implicit none

! Passed variables
class(displ_type),intent(inout) :: displ ! Displacement data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry
type(samp_type),intent(inout) :: samp    ! Sampling
type(ens_type), intent(in) :: ens        ! Ensemble

! Local variables
integer :: ic0,ic1,ic2,ic2a,jn,jc0,il0,il0i,isub,iv,its,ie,ie_sub,iter,ic0a,nc0d,ic0d,jc0d,jnmax,nnmax
integer :: ic0_ori(samp%nc2a),nn(samp%nc2a,geom%nl0)
integer,allocatable :: ic0_rac(:,:),jc0_ra(:,:,:),c0d_to_c0(:),c0_to_c0d(:),c0a_to_c0d(:)
real(kind_real) :: m2m2,fld_1,fld_2,cov,drhflt,dum,cmax
real(kind_real) :: mean,stddev,norm_tot,distsum_tot
real(kind_real),allocatable :: nn_dist(:)
real(kind_real),allocatable :: fld_ext_1(:,:,:),fld_ext_2(:,:,:)
real(kind_real),allocatable :: m1_1(:,:,:,:,:),m2_1(:,:,:,:,:)
real(kind_real),allocatable :: m1_2(:,:,:,:,:),m2_2(:,:,:,:,:)
real(kind_real),allocatable :: m11(:,:,:,:,:)
real(kind_real),allocatable :: cor(:),cor_avg(:)
real(kind_real) :: dlon_c0a(geom%nc0a),dlat_c0a(geom%nc0a)
real(kind_real) :: dlon_c2a(samp%nc2a),dlat_c2a(samp%nc2a),dist_c2a(samp%nc2a)
real(kind_real) :: dlon_c2b(samp%nc2b),dlat_c2b(samp%nc2b)
real(kind_real) :: lon_c2a_ori(samp%nc2a,geom%nl0),lat_c2a_ori(samp%nc2a,geom%nl0)
real(kind_real) :: lon_c2a(samp%nc2a),lat_c2a(samp%nc2a)
real(kind_real),allocatable :: lon_c2(:),lat_c2(:),valid_c2(:)
real(kind_real) :: x_ori(samp%nc2a),y_ori(samp%nc2a),z_ori(samp%nc2a)
real(kind_real) :: dx_ini(samp%nc2a),dy_ini(samp%nc2a),dz_ini(samp%nc2a)
real(kind_real) :: dx(samp%nc2a),dy(samp%nc2a),dz(samp%nc2a)
logical :: dichotomy,convergence
logical :: mask_c2a(samp%nc2a,geom%nl0),mask_c2(nam%nc2,geom%nl0)
logical,allocatable :: lcheck_c0d(:),mask_nn(:,:,:),mask_nn_tmp(:)
type(com_type) :: com_AD
type(mesh_type) :: mesh

! Allocation
call displ%alloc(nam,geom,samp)
allocate(ic0_rac(samp%nc2a,geom%nl0))
allocate(lon_c2(nam%nc2))
allocate(lat_c2(nam%nc2))
allocate(valid_c2(nam%nc2))
allocate(c0_to_c0d(geom%nc0))
allocate(c0a_to_c0d(geom%nc0a))

! Initialization
do il0=1,geom%nl0
   do ic2=1,nam%nc2
      ic1 = samp%c2_to_c1(ic2)
      mask_c2(ic2,il0) = samp%c1l0_log(ic1,il0)
   end do
   do ic2a=1,samp%nc2a
      ic2 = samp%c2a_to_c2(ic2a)
      mask_c2a(ic2a,il0) = mask_c2(ic2,il0)
      if (mask_c2a(ic2a,il0)) then
         ic0 = samp%c2_to_c0(ic2)
         lon_c2a_ori(ic2a,il0) = geom%lon(ic0)
         lat_c2a_ori(ic2a,il0) = geom%lat(ic0)
      end if
   end do
end do

! Copy
displ%lon_c2a = lon_c2a_ori
displ%lat_c2a = lat_c2a_ori

do its=2,nam%nts
   write(mpl%info,'(a7,a,i2)') '','Timeslot ',its

   ! Find origin point and research area center
   write(mpl%info,'(a10,a)') '','Find origin point and research area center'
   do ic2a=1,samp%nc2a
      ! Origin point
      ic2 = samp%c2a_to_c2(ic2a)
      ic0_ori(ic2a) = samp%c2_to_c0(ic2)

      if (its>2) then
         ! Research area center
         allocate(nn_dist(1))
         do il0=1,geom%nl0
            call geom%kdtree%find_nearest_neighbors(displ%lon_c2a_flt(ic2a,il0,its-1),displ%lat_c2a_flt(ic2a,il0,its-1),1, &
          & ic0_rac(ic2a:ic2a,il0),nn_dist)
         end do
         deallocate(nn_dist)
      else
         ic0_rac(ic2a,:) = ic0_ori(ic2a)
      end if
   end do

   ! Count nearest neighbors
   write(mpl%info,'(a10,a)',advance='no') '','Count nearest neighbors:'
   call mpl%prog_init(samp%nc2a*geom%nl0)
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         ! Index
         ic0 = ic0_rac(ic2a,il0)

         ! Count nearest neighbors
         call geom%kdtree%count_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nam%displ_rad,nn(ic2a,il0))
         nn(ic2a,il0) = max(1,nn(ic2a,il0))

         ! Update
         call mpl%prog_print((il0-1)*samp%nc2a+ic2a)
      end do
   end do
   write(mpl%info,'(a)') '100%'
   call flush(mpl%info)

   ! Get maximum
   nnmax = maxval(nn)

   ! Allocation
   allocate(jc0_ra(nnmax,samp%nc2a,geom%nl0))
   allocate(nn_dist(nnmax))
   allocate(mask_nn(nnmax,samp%nc2a,geom%nl0))
   allocate(m1_1(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m2_1(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m1_2(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m2_2(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m11(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))

   ! Find nearest neighbors
   write(mpl%info,'(a10,a)',advance='no') '','Find nearest neighbors:'
   call mpl%prog_init(samp%nc2a*geom%nl0)
   lcheck_c0d = samp%lcheck_c0a
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         ! Index
         ic0 = ic0_rac(ic2a,il0)

         ! Find nearest neighbors
         call geom%kdtree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nn(ic2a,il0), &
       & jc0_ra(1:nn(ic2a,il0),ic2a,il0),nn_dist(1:nn(ic2a,il0)))

         ! Check points
         do jn=1,nn(ic2a,il0)
            jc0 = jc0_ra(jn,ic2a,il0)
            lcheck_c0d(jc0) = .true.
         end do

         ! Update
         call mpl%prog_print((il0-1)*samp%nc2a+ic2a)
      end do
   end do
   write(mpl%info,'(a)') '100%'
   call flush(mpl%info)

   ! Halo D size
   nc0d = count(lcheck_c0d)

   ! Allocation
   allocate(c0d_to_c0(nc0d))
   allocate(fld_ext_1(nc0d,geom%nl0,nam%nv))
   allocate(fld_ext_2(nc0d,geom%nl0,nam%nv))

   ! Global <-> local conversions for fields
   call msi(c0_to_c0d)
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
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         ic0 = ic0_rac(ic2a,il0)
         if (geom%mask_c0(ic0,il0)) then
            do jn=1,nn(ic2a,il0)
               jc0 = jc0_ra(jn,ic2a,il0)
               mask_nn(jn,ic2a,il0) = geom%mask_c0(jc0,il0)
            end do
         else
            mask_nn(:,ic2a,il0) = .false.
         end if
      end do
   end do

   ! Initialize
   m11 = 0.0
   m2_1 = 0.0
   m2_2 = 0.0

   ! Compute moments
   write(mpl%info,'(a10,a)') '','Compute moments'
   call flush(mpl%info)
   do isub=1,ens%nsub
      if (ens%nsub==1) then
         write(mpl%info,'(a13,a)',advance='no') '','Full ensemble, member:'
      else
         write(mpl%info,'(a13,a,i4,a)',advance='no') '','Sub-ensemble ',isub,', member:'
      end if
      call flush(mpl%info)

      ! Compute covariances
      do ie_sub=1,ens%ne/ens%nsub
         write(mpl%info,'(i4)',advance='no') ie_sub
         call flush(mpl%info)

         ! Full ensemble index
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub

         do iv=1,nam%nv
            ! Halo extension
            call com_AD%ext(mpl,geom%nl0,ens%fld(:,:,iv,its-1,ie),fld_ext_1(:,:,iv))
            call com_AD%ext(mpl,geom%nl0,ens%fld(:,:,iv,its,ie),fld_ext_2(:,:,iv))

            !$omp parallel do schedule(static) private(il0,ic2a,jn,ic0,jc0,ic0d,jc0d,fld_1,fld_2)
            do il0=1,geom%nl0
               do ic2a=1,samp%nc2a
                  if (mask_c2a(ic2a,il0)) then
                     do jn=1,nn(ic2a,il0)
                        if (mask_nn(jn,ic2a,il0)) then
                           ! Indices
                           if (cor_tracker) then
                              ic0 = ic0_rac(ic2a,il0)
                           else
                              ic0 = ic0_ori(ic2a)
                           end if
                           jc0 = jc0_ra(jn,ic2a,il0)

                           ! Halo D indices
                           ic0d = c0_to_c0d(ic0)
                           jc0d = c0_to_c0d(jc0)

                           ! Copy points
                           fld_1 = fld_ext_1(ic0d,il0,iv)
                           fld_2 = fld_ext_2(jc0d,il0,iv)

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
      call flush(mpl%info)
   end do

   ! Find correlation propagation
   write(mpl%info,'(a10,a)') '','Find correlation propagation'
   call flush(mpl%info)

   do il0=1,geom%nl0
      write(mpl%info,'(a13,a,i3)') '','Level ',nam%levs(il0)
      call flush(mpl%info)

      ! Number of points
      if (its==2) call mpl%f_comm%allreduce(real(count(mask_c2a(:,il0)),kind_real),norm_tot,fckit_mpi_sum())

      !$omp parallel do schedule(static) private(ic2a,jn,iv,cov,m2m2,mean,stddev,jnmax,cmax,jc0d,jc0,ic2), &
      !$omp&                             firstprivate(cor,cor_avg,mask_nn_tmp)
      do ic2a=1,samp%nc2a
         if (mask_c2a(ic2a,il0)) then
            ! Allocation
            allocate(cor(nam%nv))
            allocate(cor_avg(nn(ic2a,il0)))
            allocate(mask_nn_tmp(nn(ic2a,il0)))

            do jn=1,nn(ic2a,il0)
               if (mask_nn(jn,ic2a,il0)) then
                  ! Compute covariance and correlation
                  do iv=1,nam%nv
                     ! Covariance
                     cov = sum(m11(jn,ic2a,il0,iv,:))/real(ens%nsub,kind_real)

                     ! Correlation
                     m2m2 = sum(m2_1(jn,ic2a,il0,iv,:))*sum(m2_2(jn,ic2a,il0,iv,:))/real(ens%nsub**2,kind_real)
                     if (m2m2>0.0) then
                        cor(iv) = cov/sqrt(m2m2)
                     else
                        call msr(cor(iv))
                     end if
                  end do

                  ! Average correlation
                  if (isanynotmsr(cor)) then
                     cor_avg(jn) = sum(cor,mask=isnotmsr(cor))/real(count(isnotmsr(cor)),kind_real)
                  else
                     call msr(cor_avg)
                  end if
               end if
            end do
            mask_nn_tmp = mask_nn(1:nn(ic2a,il0),ic2a,il0).and.isnotmsr(cor_avg)

            if (count(mask_nn_tmp)>2) then
               ! Compute the local standard-deviation
               mean = sum(cor_avg,mask=mask_nn_tmp)/real(count(mask_nn_tmp),kind_real)
               stddev = sqrt(sum((cor_avg-mean)**2,mask=mask_nn_tmp)/real(count(mask_nn_tmp)-1,kind_real))
            else
               stddev = 0.0
            end if

            ! Locate the maximum correlation
            call msi(jnmax)
            cmax = 0.0
            do jn=1,nn(ic2a,il0)
               if (mask_nn_tmp(jn)) then
                  if (cor_avg(jn)>max_std_ratio*stddev) then
                     jnmax = jn
                     cmax = cor_avg(jn)
                  end if
               end if
            end do
            if (isnotmsi(jnmax)) then
               ! Indices
               jn = jnmax
               jc0 = jc0_ra(jn,ic2a,il0)

               ! Compute displacement and distance
               dlon_c2a(ic2a) = geom%lon(jc0)-lon_c2a_ori(ic2a,il0)
               dlat_c2a(ic2a) = geom%lat(jc0)-lat_c2a_ori(ic2a,il0)
               call lonlatmod(dlon_c2a(ic2a),dlat_c2a(ic2a))
               call sphere_dist(lon_c2a_ori(ic2a,il0),lat_c2a_ori(ic2a,il0),geom%lon(jc0),geom%lat(jc0),dist_c2a(ic2a))
            else
               call msr(dlon_c2a(ic2a))
               call msr(dlat_c2a(ic2a))
               call msr(dist_c2a(ic2a))
            end if

            ! Release memory
            deallocate(cor)
            deallocate(cor_avg)
            deallocate(mask_nn_tmp)
         end if
      end do
      !$omp end parallel do

      ! Copy lon/lat
      do ic2a=1,samp%nc2a
         if (isnotmsr(dlon_c2a(ic2a)).and.isnotmsr(dlat_c2a(ic2a))) then
            lon_c2a(ic2a) = lon_c2a_ori(ic2a,il0)+dlon_c2a(ic2a)
            lat_c2a(ic2a) = lat_c2a_ori(ic2a,il0)+dlat_c2a(ic2a)
            call lonlatmod(lon_c2a(ic2a),lat_c2a(ic2a))
         else
            call msr(lon_c2a(ic2a))
            call msr(lat_c2a(ic2a))
         end if
      end do

      ! Check raw mesh
      call mpl%loc_to_glb(samp%nc2a,lon_c2a,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lon_c2)
      call mpl%loc_to_glb(samp%nc2a,lat_c2a,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lat_c2)
      if (mpl%main) then
         mesh = samp%mesh%copy()
         call mesh%trans(lon_c2,lat_c2)
         call mesh%check(valid_c2)
         displ%valid(0,il0,its) = sum(valid_c2,mask=mask_c2(:,il0))/real(count((mask_c2(:,il0))),kind_real)
      end if
      call mpl%f_comm%broadcast(displ%valid(0,il0,its),mpl%ioproc-1)
      displ%rhflt(0,il0,its) = 0.0

      ! Average distance
      call mpl%f_comm%allreduce(sum(dist_c2a,mask=mask_c2a(:,il0)),distsum_tot,fckit_mpi_sum())
      displ%dist(0,il0,its) = distsum_tot/norm_tot

      ! Copy
      displ%lon_c2a_raw(:,il0,its) = lon_c2a
      displ%lat_c2a_raw(:,il0,its) = lat_c2a
      displ%dist_c2a_raw(:,il0,its) = dist_c2a
   end do

   ! Filter displacement
   write(mpl%info,'(a10,a)') '','Filter displacement'
   call flush(mpl%info)

   do il0=1,geom%nl0
      write(mpl%info,'(a13,a,i3)') '','Level ',nam%levs(il0)
      call flush(mpl%info)

      if (nam%displ_niter>0) then
         ! Convert to cartesian coordinates
         call trans(samp%nc2a,lat_c2a_ori(:,il0),lon_c2a_ori(:,il0),x_ori,y_ori,z_ori)
         call trans(samp%nc2a,lat_c2a,lon_c2a,dx_ini,dy_ini,dz_ini)

         ! Dichotomy initialization
         do ic2a=1,samp%nc2a
            if (isnotmsr(dx_ini(ic2a))) dx_ini(ic2a) = dx_ini(ic2a)-x_ori(ic2a)
            if (isnotmsr(dy_ini(ic2a))) dy_ini(ic2a) = dy_ini(ic2a)-y_ori(ic2a)
            if (isnotmsr(dz_ini(ic2a))) dz_ini(ic2a) = dz_ini(ic2a)-z_ori(ic2a)
         end do
         convergence = .true.
         dichotomy = .false.
         displ%rhflt(1,il0,its) = nam%displ_rhflt
         drhflt = displ%rhflt(1,il0,its)

         do iter=1,nam%displ_niter
            ! Copy increment
            dx = dx_ini
            dy = dy_ini
            dz = dz_ini

            ! Median filter to remove extreme values
            call samp%diag_filter(mpl,nam,geom,il0,'median',displ%rhflt(iter,il0,its),dx)
            call samp%diag_filter(mpl,nam,geom,il0,'median',displ%rhflt(iter,il0,its),dy)
            call samp%diag_filter(mpl,nam,geom,il0,'median',displ%rhflt(iter,il0,its),dz)

            ! Average filter to smooth displacement
            call samp%diag_filter(mpl,nam,geom,il0,'gc99',displ%rhflt(iter,il0,its),dx)
            call samp%diag_filter(mpl,nam,geom,il0,'gc99',displ%rhflt(iter,il0,its),dy)
            call samp%diag_filter(mpl,nam,geom,il0,'gc99',displ%rhflt(iter,il0,its),dz)

            ! Back to spherical coordinates
            dx = dx+x_ori
            dy = dy+y_ori
            dz = dz+z_ori
            do ic2a=1,samp%nc2a
               call scoord(dx(ic2a),dy(ic2a),dz(ic2a),lat_c2a(ic2a),lon_c2a(ic2a),dum)
            end do

            ! Reduce distance with respect to boundary
            do ic2a=1,samp%nc2a
               if (mask_c2a(ic2a,il0)) then
                  ic2 = samp%c2a_to_c2(ic2a)
                  call reduce_arc(lon_c2a_ori(ic2a,il0),lat_c2a_ori(ic2a,il0),lon_c2a(ic2a),lat_c2a(ic2a),samp%mesh%bdist(ic2), &
                & dist_c2a(ic2a))
                  dlon_c2a(ic2a) = lon_c2a(ic2a)-lon_c2a_ori(ic2a,il0)
                  dlat_c2a(ic2a) = lat_c2a(ic2a)-lat_c2a_ori(ic2a,il0)
               end if
            end do

            ! Copy lon/lat
            do ic2a=1,samp%nc2a
               lon_c2a(ic2a) = lon_c2a_ori(ic2a,il0)+dlon_c2a(ic2a)
               lat_c2a(ic2a) = lat_c2a_ori(ic2a,il0)+dlat_c2a(ic2a)
               call lonlatmod(lon_c2a(ic2a),lat_c2a(ic2a))
            end do

            ! Check mesh
            call mpl%loc_to_glb(samp%nc2a,lon_c2a,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lon_c2)
            call mpl%loc_to_glb(samp%nc2a,lat_c2a,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.false.,lat_c2)
            if (mpl%main) then
               mesh = samp%mesh%copy()
               call mesh%trans(lon_c2,lat_c2)
               call mesh%check(valid_c2)
               displ%valid(iter,il0,its) = sum(valid_c2,mask=mask_c2(:,il0))/real(count((mask_c2(:,il0))),kind_real)
            end if
            call mpl%f_comm%broadcast(displ%valid(iter,il0,its),mpl%ioproc-1)

            ! Compute distances
            do ic2a=1,samp%nc2a
               if (mask_c2a(ic2a,il0)) call sphere_dist(lon_c2a_ori(ic2a,il0),lat_c2a_ori(ic2a,il0), &
             & lon_c2a(ic2a),lat_c2a(ic2a),dist_c2a(ic2a))
            end do

            ! Average distance
            call mpl%f_comm%allreduce(sum(dist_c2a,mask=mask_c2a(:,il0)),distsum_tot,fckit_mpi_sum())
            displ%dist(iter,il0,its) = distsum_tot/norm_tot

            ! Print results
            write(mpl%info,'(a13,a,i2,a,f10.2,a,f6.2,a,f6.2,a,f7.2,a)') '','Iteration ',iter,': rhflt = ', &
          & displ%rhflt(iter,il0,its)*reqkm,' km, valid points: ',100.0*displ%valid(0,il0,its),'% ~> ', &
          & 100.0*displ%valid(iter,il0,its),'%, average displacement = ',displ%dist(iter,il0,its)*reqkm,' km'
            call flush(mpl%info)

            ! Update support radius
            if (displ%valid(iter,il0,its)<1.0-nam%displ_tol) then
               ! Increase filtering support radius
               if (dichotomy) then
                  drhflt = 0.5*drhflt
                  if (iter<nam%displ_niter) displ%rhflt(iter+1,il0,its) = displ%rhflt(iter,il0,its)+drhflt
               else
                  convergence = .false.
                  if (iter<nam%displ_niter) displ%rhflt(iter+1,il0,its) = displ%rhflt(iter,il0,its)+drhflt
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
               if (iter<nam%displ_niter) displ%rhflt(iter+1,il0,its) = displ%rhflt(iter,il0,its)-drhflt

               ! Save values
               displ%lon_c2a_flt(:,il0,its) = lon_c2a
               displ%lat_c2a_flt(:,il0,its) = lat_c2a
               displ%dist_c2a_flt(:,il0,its) = dist_c2a
            end if
         end do

         ! Check convergence
         if (.not.convergence) call mpl%abort('iterative filtering failed')
      else
         ! Print results
         write(mpl%info,'(a10,a22,f10.2,a,f6.2,a,f7.2,a)') '','Raw displacement: rhflt = ', &
       & displ%rhflt(0,il0,its)*reqkm,' km, valid points: ',100.0*displ%valid(0,il0,its),'%, average displacement = ', &
       & displ%dist(0,il0,its)*reqkm,' km'
         call flush(mpl%info)

         ! Copy
         displ%lon_c2a_flt(:,il0,its) = displ%lon_c2a_raw(:,il0,its)
         displ%lat_c2a_flt(:,il0,its) = displ%lat_c2a_raw(:,il0,its)
         displ%dist_c2a_flt(:,il0,its) = displ%dist_c2a_raw(:,il0,its)
      end if

      ! Displacement interpolation
      do ic2a=1,samp%nc2a
         call lonlatmod(dlon_c2a(ic2a),dlat_c2a(ic2a))
      end do
      il0i = min(il0,geom%nl0i)
      call samp%com_AB%ext(mpl,dlon_c2a,dlon_c2b)
      call samp%com_AB%ext(mpl,dlat_c2a,dlat_c2b)
      call samp%h(il0i)%apply(mpl,dlon_c2b,dlon_c0a)
      call samp%h(il0i)%apply(mpl,dlat_c2b,dlat_c0a)

      ! Displaced grid
      do ic0a=1,geom%nc0a
         ic0 = geom%c0a_to_c0(ic0a)
         samp%displ_lon(ic0a,il0,its) = geom%lon(ic0)+dlon_c0a(ic0a)
         samp%displ_lat(ic0a,il0,its) = geom%lat(ic0)+dlat_c0a(ic0a)
         call lonlatmod(samp%displ_lon(ic0a,il0,its),samp%displ_lat(ic0a,il0,its))
      end do
   end do

   ! Release memory
   deallocate(jc0_ra)
   deallocate(nn_dist)
   deallocate(mask_nn)
   deallocate(m1_1)
   deallocate(m2_1)
   deallocate(m1_2)
   deallocate(m2_2)
   deallocate(m11)
   deallocate(c0d_to_c0)
   call com_AD%dealloc
   deallocate(fld_ext_1)
   deallocate(fld_ext_2)
end do

! Displaced grid for timeslot 1
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      samp%displ_lon(ic0a,il0,1) = geom%lon(ic0)
      samp%displ_lat(ic0a,il0,1) = geom%lat(ic0)
   end do
end do

end subroutine displ_compute

!----------------------------------------------------------------------
! Subroutine: displ_write
! Purpose: write displacement data
!----------------------------------------------------------------------
subroutine displ_write(displ,mpl,nam,geom,samp,filename)

implicit none

! Passed variables
class(displ_type),intent(in) :: displ   ! Displacement data
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(samp_type),intent(in) :: samp      ! Sampling
character(len=*),intent(in) :: filename ! File name

! Local variables
integer :: ncid,nc2_id,nl0_id,nts_id,na_id,displ_niter_id,vunit_id,valid_id,dist_id,rhflt_id,larc_s_id,larc_e_id
integer :: lon_c2_id,lat_c2_id,lon_c2_raw_id,lat_c2_raw_id,dist_c2_raw_id,lon_c2_flt_id,lat_c2_flt_id,dist_c2_flt_id
integer :: iproc,its,il0,ic2a,ic2,i
integer,allocatable :: c2a_to_c2(:)
real(kind_real),allocatable :: sbuf(:),rbuf(:),lon_c2(:,:),lat_c2(:,:)
real(kind_real),allocatable :: lon_c2_raw(:,:,:),lat_c2_raw(:,:,:),dist_c2_raw(:,:,:)
real(kind_real),allocatable :: lon_c2_flt(:,:,:),lat_c2_flt(:,:,:),dist_c2_flt(:,:,:)
character(len=1024) :: subr = 'displ_write'
type(fckit_mpi_status) :: status

! Allocation
allocate(sbuf(samp%nc2a*geom%nl0*(2+(nam%nts-1)*6)))

! Prepare buffer
i = 1
do il0=1,geom%nl0
   do ic2a=1,samp%nc2a
      sbuf(i) = displ%lon_c2a(ic2a,il0)*rad2deg
      i = i+1
      sbuf(i) = displ%lat_c2a(ic2a,il0)*rad2deg
      i = i+1
      do its=2,nam%nts
         sbuf(i) = displ%lon_c2a_raw(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = displ%lat_c2a_raw(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = displ%dist_c2a_raw(ic2a,il0,its)*reqkm
         i = i+1
         sbuf(i) = displ%lon_c2a_flt(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = displ%lat_c2a_flt(ic2a,il0,its)*rad2deg
         i = i+1
         sbuf(i) = displ%dist_c2a_flt(ic2a,il0,its)*reqkm
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

if (mpl%main) then
   ! Create file
   call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call nam%ncwrite(mpl,ncid)

   ! Define dimensions
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2',nam%nc2,nc2_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nts',nam%nts-1,nts_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'niter',nam%displ_niter+1,displ_niter_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'na',samp%mesh%na,na_id))

   ! Define variables
   call mpl%ncerr(subr,nf90_def_var(ncid,'vunit',ncfloat,(/nc2_id,nl0_id/),vunit_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'valid',ncfloat,(/displ_niter_id,nl0_id,nts_id/),valid_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,valid_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'dist',ncfloat,(/displ_niter_id,nl0_id,nts_id/),dist_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,dist_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'rhflt',ncfloat,(/displ_niter_id,nl0_id,nts_id/),rhflt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,rhflt_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_c2',ncfloat,(/nc2_id,nl0_id/),lon_c2_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_c2_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_c2',ncfloat,(/nc2_id,nl0_id/),lat_c2_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_c2_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),lon_c2_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_c2_raw_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),lat_c2_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_c2_raw_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'dist_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),dist_c2_raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,dist_c2_raw_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),lon_c2_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_c2_flt_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),lat_c2_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_c2_flt_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'dist_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),dist_c2_flt_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,dist_c2_flt_id,'_FillValue',msvalr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'larc_s',nf90_int,(/na_id/),larc_s_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'larc_e',nf90_int,(/na_id/),larc_e_id))

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write data
   call mpl%ncerr(subr,nf90_put_var(ncid,vunit_id,geom%vunit(samp%c2_to_c0,:)))
   call mpl%ncerr(subr,nf90_put_var(ncid,valid_id,displ%valid))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_id,displ%dist*reqkm))
   call mpl%ncerr(subr,nf90_put_var(ncid,rhflt_id,displ%rhflt*reqkm))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_id,lon_c2))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_id,lat_c2))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_raw_id,lon_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_raw_id,lat_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_c2_raw_id,dist_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_flt_id,lon_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_flt_id,lat_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_c2_flt_id,dist_c2_flt))
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

end subroutine displ_write

end module type_displ
