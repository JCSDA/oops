!----------------------------------------------------------------------
! Module: type_displ
!> Purpose: displacement data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_displ

use netcdf
!$ use omp_lib
use tools_const, only: req,reqkm,rad2deg,deg2rad,msvalr
use tools_display, only: msgerror,prog_init,prog_print
use tools_func, only: lonlatmod,sphere_dist,reduce_arc,vector_product
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_nc, only: ncfloat,ncerr
use tools_qsort, only: qsort
use tools_stripack, only: trans,scoord
use type_com, only: com_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_hdata, only: hdata_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! Displacement data derived type
type displ_type
   integer :: niter                                   !< Number of stored iterations
   real(kind_real),allocatable :: dist(:,:,:)         !< Displacement distance
   real(kind_real),allocatable :: valid(:,:,:)        !< Displacement validity
   real(kind_real),allocatable :: rhflt(:,:,:)        !< Displacement filtering support radius

   real(kind_real),allocatable :: lon_c2a(:,:)        !< Longitude origin
   real(kind_real),allocatable :: lat_c2a(:,:)        !< Latitude origin
   real(kind_real),allocatable :: lon_c2a_raw(:,:,:)  !< Raw displaced longitude
   real(kind_real),allocatable :: lat_c2a_raw(:,:,:)  !< Raw displaced latitude
   real(kind_real),allocatable :: dist_c2a_raw(:,:,:) !< Raw displacement distance
   real(kind_real),allocatable :: lon_c2a_flt(:,:,:)  !< Filtered displaced longitude
   real(kind_real),allocatable :: lat_c2a_flt(:,:,:)  !< Filtered displaced latitude
   real(kind_real),allocatable :: dist_c2a_flt(:,:,:) !< Displacement distance, filtered
contains
   procedure :: alloc => displ_alloc
   procedure :: dealloc => displ_dealloc
   procedure :: compute => displ_compute
   procedure :: write => displ_write
end type displ_type

character(len=1024) :: displ_method = 'cor_center_mass' !< Displacement computation method
real(kind_real),parameter :: cor_th = 0.2_kind_real     !< Correlation threshold

private
public :: displ_type

contains

!----------------------------------------------------------------------
! Subroutine: displ_alloc
!> Purpose: displacement data allocation
!----------------------------------------------------------------------
subroutine displ_alloc(displ,nam,geom,hdata)

implicit none

! Passed variables
class(displ_type),intent(inout) :: displ !< Displacement data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
type(hdata_type),intent(in) :: hdata     !< HDIAG data

! Allocation
allocate(displ%dist(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%valid(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%rhflt(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%lon_c2a(hdata%nc2a,geom%nl0))
allocate(displ%lat_c2a(hdata%nc2a,geom%nl0))
allocate(displ%lon_c2a_raw(hdata%nc2a,geom%nl0,2:nam%nts))
allocate(displ%lat_c2a_raw(hdata%nc2a,geom%nl0,2:nam%nts))
allocate(displ%dist_c2a_raw(hdata%nc2a,geom%nl0,2:nam%nts))
allocate(displ%lon_c2a_flt(hdata%nc2a,geom%nl0,2:nam%nts))
allocate(displ%lat_c2a_flt(hdata%nc2a,geom%nl0,2:nam%nts))
allocate(displ%dist_c2a_flt(hdata%nc2a,geom%nl0,2:nam%nts))

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
!> Purpose: displacement data deallocation
!----------------------------------------------------------------------
subroutine displ_dealloc(displ)

implicit none

! Passed variables
class(displ_type),intent(inout) :: displ !< Displacement data

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
!> Purpose: compute correlation maximum displacement
!----------------------------------------------------------------------
subroutine displ_compute(displ,nam,geom,hdata,ens)

implicit none

! Passed variables
class(displ_type),intent(inout) :: displ !< Displacement data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
type(hdata_type),intent(inout) :: hdata  !< HDIAG data
type(ens_type), intent(in) :: ens        !< Ensemble

! Local variables
integer :: ic0,ic1,ic2,ic2a,jc0,jc1,il0,il0i,isub,iv,its,ie,ie_sub,iter,ic0a,jc0d
integer,allocatable :: ind(:)
real(kind_real) :: fac4,fac6,m11_avg,m2m2_avg,fld_1,fld_2,drhflt,dum,distsum,norm
real(kind_real) :: lon_target,lat_target,rad_target,x_cm,y_cm,z_cm,n_cm
real(kind_real) :: norm_tot,distsum_tot
real(kind_real),allocatable :: fld_ext(:,:,:,:)
real(kind_real),allocatable :: m1_1(:,:,:,:,:,:),m2_1(:,:,:,:,:,:)
real(kind_real),allocatable :: m1_2(:,:,:,:,:,:),m2_2(:,:,:,:,:,:)
real(kind_real),allocatable :: m11(:,:,:,:,:,:)
real(kind_real),allocatable :: cor(:),cor_avg(:)
real(kind_real),allocatable :: x(:),y(:),z(:)
real(kind_real) :: dlon_c0a(geom%nc0a),dlat_c0a(geom%nc0a)
real(kind_real) :: dlon_c2a(hdata%nc2a),dlat_c2a(hdata%nc2a),dist_c2a(hdata%nc2a)
real(kind_real) :: dlon_c2b(hdata%nc2b),dlat_c2b(hdata%nc2b)
real(kind_real) :: lon_c2a_ori(hdata%nc2a,geom%nl0),lat_c2a_ori(hdata%nc2a,geom%nl0)
real(kind_real) :: lon_c2a(hdata%nc2a),lat_c2a(hdata%nc2a)
real(kind_real) :: lon_c2(hdata%nc2),lat_c2(hdata%nc2),valid_c2(hdata%nc2)
real(kind_real) :: x_ori(hdata%nc2a),y_ori(hdata%nc2a),z_ori(hdata%nc2a)
real(kind_real) :: dx_ini(hdata%nc2a),dy_ini(hdata%nc2a),dz_ini(hdata%nc2a)
real(kind_real) :: dx(hdata%nc2a),dy(hdata%nc2a),dz(hdata%nc2a)
logical :: dichotomy,convergence
logical :: mask_c2a(hdata%nc2a,geom%nl0),mask_c2(hdata%nc2,geom%nl0)
type(mesh_type) :: mesh

! Allocation
call displ%alloc(nam,geom,hdata)
allocate(fld_ext(hdata%nc0d,geom%nl0,nam%nv,2:nam%nts))
allocate(m1_1(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,ens%nsub))
allocate(m2_1(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,ens%nsub))
allocate(m1_2(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,ens%nsub))
allocate(m2_2(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,ens%nsub))
allocate(m11(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,ens%nsub))

! Initialization
m1_1 = 0.0
m2_1 = 0.0
m1_2 = 0.0
m2_2 = 0.0
m11 = 0.0

! Initialization
do il0=1,geom%nl0
   do ic2=1,hdata%nc2
      ic1 = hdata%c2_to_c1(ic2)
      mask_c2(ic2,il0) = hdata%c1l0_log(ic1,il0)
   end do
   do ic2a=1,hdata%nc2a
      ic2 = hdata%c2a_to_c2(ic2a)
      mask_c2a(ic2a,il0) = mask_c2(ic2,il0)
      if (mask_c2a(ic2a,il0)) then
         ic0 = hdata%c2_to_c0(ic2)
         lon_c2a_ori(ic2a,il0) = geom%lon(ic0)
         lat_c2a_ori(ic2a,il0) = geom%lat(ic0)
      end if
   end do
end do

! Copy
displ%lon_c2a = lon_c2a_ori
displ%lat_c2a = lat_c2a_ori

! Compute moments
write(mpl%unit,'(a7,a)') '','Compute moments'
call flush(mpl%unit)
do isub=1,ens%nsub
   if (ens%nsub==1) then
      write(mpl%unit,'(a10,a)',advance='no') '','Full ensemble, member:'
   else
      write(mpl%unit,'(a10,a,i4,a)',advance='no') '','Sub-ensemble ',isub,', member:'
   end if
   call flush(mpl%unit)

   ! Compute centered moments iteratively
   do ie_sub=1,ens%ne/ens%nsub
      write(mpl%unit,'(i4)',advance='no') ie_sub
      call flush(mpl%unit)

      ! Full ensemble index
      ie = ie_sub+(isub-1)*ens%ne/ens%nsub

      ! Computation factors
      fac4 = 1.0/real(ie_sub,kind_real)
      fac6 = real(ie_sub-1,kind_real)/real(ie_sub,kind_real)

      do its=2,nam%nts
         do iv=1,nam%nv
            ! Halo extension
            call hdata%com_AD%ext(geom%nl0,ens%fld(:,:,iv,its,ie),fld_ext(:,:,iv,its))
         end do
      end do

      do its=2,nam%nts
         do iv=1,nam%nv
            !$omp parallel do schedule(static) private(il0,ic2a,ic2,ic1,jc1,ic0,jc0,ic0a,jc0d,fld_1,fld_2)
            do il0=1,geom%nl0
               do ic2a=1,hdata%nc2a
                  ic2 = hdata%c2a_to_c2(ic2a)
                  ic1 = hdata%c2_to_c1(ic2)
                  if (hdata%c1l0_log(ic1,il0)) then
                     do jc1=1,nam%nc1
                        if (hdata%displ_mask(jc1,ic2,min(il0,geom%nl0i))) then
                           ! Indices
                           ic0 = hdata%c2_to_c0(ic2)
                           jc0 = hdata%c1_to_c0(jc1)
                           ic0a = geom%c0_to_c0a(ic0)
                           jc0d = hdata%c0_to_c0d(jc0)

                           ! Copy points
                           fld_1 = ens%fld(ic0a,il0,iv,1,ie)
                           fld_2 = fld_ext(jc0d,il0,iv,its)

                           ! Remove means
                           fld_1 = fld_1 - m1_1(jc1,ic2a,il0,iv,its,isub)
                           fld_2 = fld_2 - m1_2(jc1,ic2a,il0,iv,its,isub)

                           ! Update high-order moments
                           if (ie_sub>1) then
                              ! Covariance
                              m11(jc1,ic2a,il0,iv,its,isub) = m11(jc1,ic2a,il0,iv,its,isub)+fac6*fld_1*fld_2

                              ! Variances
                              m2_1(jc1,ic2a,il0,iv,its,isub) = m2_1(jc1,ic2a,il0,iv,its,isub)+fac6*fld_1**2
                              m2_2(jc1,ic2a,il0,iv,its,isub) = m2_2(jc1,ic2a,il0,iv,its,isub)+fac6*fld_2**2
                           end if

                           ! Update means
                           m1_1(jc1,ic2a,il0,iv,its,isub) = m1_1(jc1,ic2a,il0,iv,its,isub)+fac4*fld_1
                           m1_2(jc1,ic2a,il0,iv,its,isub) = m1_2(jc1,ic2a,il0,iv,its,isub)+fac4*fld_2
                        end if
                     end do
                  end if
               end do
            end do
            !$omp end parallel do
         end do
      end do
   end do
   write(mpl%unit,'(a)') ''
   call flush(mpl%unit)
end do

! Find correlation propagation
write(mpl%unit,'(a7,a)') '','Find correlation propagation'
call flush(mpl%unit)

do its=2,nam%nts
   do il0=1,geom%nl0
      write(mpl%unit,'(a10,a,i2,a,i3)') '','Timeslot ',its,' - level ',nam%levs(il0)
      call flush(mpl%unit)

      ! Number of points
      norm = real(count(mask_c2a(:,il0)),kind_real)
      call mpl%allreduce_sum(norm,norm_tot)

      !$omp parallel do schedule(static) private(ic2a,ic2,jc1,jc0,iv,m11_avg,m2m2_avg,lon_target,lat_target,rad_target), &
      !$omp&                             private(x_cm,y_cm,z_cm,n_cm) firstprivate(cor,cor_avg,ind,x,y,z)
      do ic2a=1,hdata%nc2a
         ic2 = hdata%c2a_to_c2(ic2a)
         if (mask_c2a(ic2a,il0)) then
            ! Allocation
            allocate(cor(nam%nv))
            allocate(cor_avg(nam%nc1))

            do jc1=1,nam%nc1
               ! Initialization
               call msr(cor_avg(jc1))

               if (hdata%displ_mask(jc1,ic2,min(il0,geom%nl0i))) then
                  ! Compute correlation for each variable
                  do iv=1,nam%nv
                     ! Correlation
                     m11_avg = sum(m11(jc1,ic2a,il0,iv,its,:))/real(ens%nsub,kind_real)
                     m2m2_avg = sum(m2_1(jc1,ic2a,il0,iv,its,:))*sum(m2_2(jc1,ic2a,il0,iv,its,:))/real(ens%nsub**2,kind_real)
                     if (m2m2_avg>0.0) then
                        cor(iv) = m11_avg/sqrt(m2m2_avg)
                     else
                        call msr(cor(iv))
                     end if
                  end do

                  ! Average correlations
                  if (isanynotmsr(cor)) then
                     cor_avg(jc1) = sum(cor,mask=isnotmsr(cor))/real(count(isnotmsr(cor)),kind_real)
                  else
                     call msgerror('average correlation contains missing values only')
                  end if
               end if
            end do

            select case (trim(displ_method))
            case ('cor_max')
               ! Locate the maximum correlation, with a correlation threshold

               ! Allocation
               allocate(ind(1))

               if (maxval(cor_avg)>cor_th) then
                  ! Find maximum
                  ind = maxloc(cor_avg)
                  jc1 = ind(1)
                  jc0 = hdata%c1_to_c0(jc1)
                  lon_target = geom%lon(jc0)
                  lat_target = geom%lon(jc0)
               else
                  lon_target = 0.0
                  lat_target = 0.0
               end if

               ! Release memory
               deallocate(ind)
            case ('cor_center_mass')
               ! Locate the correlation center of mass, with a correlation threshold

               ! Allocation
               allocate(x(1))
               allocate(y(1))
               allocate(z(1))

               ! Initialization
               x_cm = 0.0
               y_cm = 0.0
               z_cm = 0.0
               n_cm = 0.0

               do jc1=1,nam%nc1
                  if (cor_avg(jc1)>max(cor_th,0.5*maxval(cor_avg))) then
                     ! Compute cartesian coordinates
                     jc0 = hdata%c1_to_c0(jc1)
                     call trans(1,geom%lat(jc0),geom%lon(jc0),x,y,z)

                     ! Update center of mass
                     x_cm = x_cm+cor_avg(jc1)*x(1)
                     y_cm = y_cm+cor_avg(jc1)*y(1)
                     z_cm = z_cm+cor_avg(jc1)*z(1)
                     n_cm = n_cm+cor_avg(jc1)
                  end if
               end do

               if (n_cm>0.0) then
                  ! Final center of mass
                  x_cm = x_cm/n_cm
                  y_cm = y_cm/n_cm
                  z_cm = z_cm/n_cm

                  ! Back to spherical coordinates
                  call scoord(x_cm,y_cm,z_cm,lat_target,lon_target,rad_target)
               else
                  lon_target = 0.0
                  lat_target = 0.0
               end if

               ! Release memory
               deallocate(x)
               deallocate(y)
               deallocate(z)
            end select

            ! Compute displacement and distance
            dlon_c2a(ic2a) = lon_target-lon_c2a_ori(ic2a,il0)
            dlat_c2a(ic2a) = lat_target-lat_c2a_ori(ic2a,il0)
            call lonlatmod(dlon_c2a(ic2a),dlat_c2a(ic2a))
            call sphere_dist(lon_c2a_ori(ic2a,il0),lat_c2a_ori(ic2a,il0),lon_target,lat_target,dist_c2a(ic2a))

            ! Release memory
            deallocate(cor)
            deallocate(cor_avg)
         end if
      end do
      !$omp end parallel do

      ! Copy lon/lat
      do ic2a=1,hdata%nc2a
         lon_c2a(ic2a) = lon_c2a_ori(ic2a,il0)+dlon_c2a(ic2a)
         lat_c2a(ic2a) = lat_c2a_ori(ic2a,il0)+dlat_c2a(ic2a)
         call lonlatmod(lon_c2a(ic2a),lat_c2a(ic2a))
      end do

      ! Check raw mesh
      call mpl%gatherv(hdata%nc2a,lon_c2a,hdata%proc_to_nc2a,hdata%nc2,lon_c2)
      call mpl%gatherv(hdata%nc2a,lat_c2a,hdata%proc_to_nc2a,hdata%nc2,lat_c2)
      if (mpl%main) then
         mesh = hdata%mesh%copy()
         call mesh%trans(lon_c2,lat_c2)
         call mesh%check(valid_c2)
         displ%valid(0,il0,its) = sum(valid_c2,mask=mask_c2(:,il0))/real(count((mask_c2(:,il0))),kind_real)
      end if
      call mpl%bcast(displ%valid(0,il0,its))
      displ%rhflt(0,il0,its) = 0.0

      ! Average distance
      distsum = sum(dist_c2a,mask=mask_c2a(:,il0))
      call mpl%allreduce_sum(distsum,distsum_tot)
      displ%dist(0,il0,its) = distsum_tot/norm_tot

      ! Copy
      displ%lon_c2a_raw(:,il0,its) = lon_c2a
      displ%lat_c2a_raw(:,il0,its) = lat_c2a
      displ%dist_c2a_raw(:,il0,its) = dist_c2a

      if (nam%displ_niter>0) then
         ! Filter displacement

         ! Compute raw displacement in cartesian coordinates
         call trans(hdata%nc2a,lat_c2a_ori(:,il0),lon_c2a_ori(:,il0),x_ori,y_ori,z_ori)
         call trans(hdata%nc2a,lat_c2a,lon_c2a,dx_ini,dy_ini,dz_ini)
         dx_ini = dx_ini-x_ori
         dy_ini = dy_ini-y_ori
         dz_ini = dz_ini-z_ori

         ! Iterative filtering
         convergence = .true.
         dichotomy = .false.

         ! Dichotomy initialization
         displ%rhflt(1,il0,its) = nam%displ_rhflt
         drhflt = displ%rhflt(1,il0,its)

         do iter=1,nam%displ_niter
            ! Copy increment
            dx = dx_ini
            dy = dy_ini
            dz = dz_ini

            ! Median filter to remove extreme values
            call hdata%diag_filter(geom,il0,'median',displ%rhflt(iter,il0,its),dx)
            call hdata%diag_filter(geom,il0,'median',displ%rhflt(iter,il0,its),dy)
            call hdata%diag_filter(geom,il0,'median',displ%rhflt(iter,il0,its),dz)

            ! Average filter to smooth displacement
            call hdata%diag_filter(geom,il0,'gc99',displ%rhflt(iter,il0,its),dx)
            call hdata%diag_filter(geom,il0,'gc99',displ%rhflt(iter,il0,its),dy)
            call hdata%diag_filter(geom,il0,'gc99',displ%rhflt(iter,il0,its),dz)

            ! Back to spherical coordinates
            dx = dx+x_ori
            dy = dy+y_ori
            dz = dz+z_ori
            do ic2a=1,hdata%nc2a
               call scoord(dx(ic2a),dy(ic2a),dz(ic2a),lat_c2a(ic2a),lon_c2a(ic2a),dum)
            end do

            ! Reduce distance with respect to boundary
            do ic2a=1,hdata%nc2a
               if (mask_c2a(ic2a,il0)) then
                  ic2 = hdata%c2a_to_c2(ic2a)
                  call reduce_arc(lon_c2a_ori(ic2a,il0),lat_c2a_ori(ic2a,il0),lon_c2a(ic2a),lat_c2a(ic2a),hdata%mesh%bdist(ic2), &
                & dist_c2a(ic2a))
                  dlon_c2a(ic2a) = lon_c2a(ic2a)-lon_c2a_ori(ic2a,il0)
                  dlat_c2a(ic2a) = lat_c2a(ic2a)-lat_c2a_ori(ic2a,il0)
               end if
            end do

            ! Copy lon/lat
            do ic2a=1,hdata%nc2a
               lon_c2a(ic2a) = lon_c2a_ori(ic2a,il0)+dlon_c2a(ic2a)
               lat_c2a(ic2a) = lat_c2a_ori(ic2a,il0)+dlat_c2a(ic2a)
               call lonlatmod(lon_c2a(ic2a),lat_c2a(ic2a))
            end do

            ! Check mesh
            call mpl%gatherv(hdata%nc2a,lon_c2a,hdata%proc_to_nc2a,hdata%nc2,lon_c2)
            call mpl%gatherv(hdata%nc2a,lat_c2a,hdata%proc_to_nc2a,hdata%nc2,lat_c2)
            if (mpl%main) then
               mesh = hdata%mesh%copy()
               call mesh%trans(lon_c2,lat_c2)
               call mesh%check(valid_c2)
               displ%valid(iter,il0,its) = sum(valid_c2,mask=mask_c2(:,il0))/real(count((mask_c2(:,il0))),kind_real)
            end if
            call mpl%bcast(displ%valid(iter,il0,its))

            ! Compute distances
            do ic2a=1,hdata%nc2a
               if (mask_c2a(ic2a,il0)) call sphere_dist(lon_c2a_ori(ic2a,il0),lat_c2a_ori(ic2a,il0), &
             & lon_c2a(ic2a),lat_c2a(ic2a),dist_c2a(ic2a))
            end do

            ! Average distance
            distsum = sum(dist_c2a,mask=mask_c2a(:,il0))
            call mpl%allreduce_sum(distsum,distsum_tot)
            displ%dist(iter,il0,its) = distsum_tot/norm_tot

            ! Print results
            write(mpl%unit,'(a13,a,i2,a,f10.2,a,f6.2,a,f6.2,a,f7.2,a)') '','Iteration ',iter,': rhflt = ', &
          & displ%rhflt(iter,il0,its)*reqkm,' km, valid points: ',100.0*displ%valid(0,il0,its),'% ~> ', &
          & 100.0*displ%valid(iter,il0,its),'%, average displacement = ',displ%dist(iter,il0,its)*reqkm,' km'
            call flush(mpl%unit)

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
            end if
         end do

         ! Copy
         displ%lon_c2a_flt(:,il0,its) = lon_c2a
         displ%lat_c2a_flt(:,il0,its) = lat_c2a
         displ%dist_c2a_flt(:,il0,its) = dist_c2a

         ! Check convergence
         if (.not.convergence) call msgerror('iterative filtering failed')
      else
         ! Print results
         write(mpl%unit,'(a10,a22,f10.2,a,f6.2,a,f7.2,a)') '','Raw displacement: rhflt = ', &
       & displ%rhflt(0,il0,its)*reqkm,' km, valid points: ',100.0*displ%valid(0,il0,its),'%, average displacement = ', &
       & displ%dist(0,il0,its)*reqkm,' km'
         call flush(mpl%unit)
      end if

      ! Displacement interpolation
      do ic2a=1,hdata%nc2a
         call lonlatmod(dlon_c2a(ic2a),dlat_c2a(ic2a))
      end do
      il0i = min(il0,geom%nl0i)
      call hdata%com_AB%ext(dlon_c2a,dlon_c2b)
      call hdata%com_AB%ext(dlat_c2a,dlat_c2b)
      call hdata%h(il0i)%apply(dlon_c2b,dlon_c0a)
      call hdata%h(il0i)%apply(dlat_c2b,dlat_c0a)

      ! Displaced grid
      do ic0a=1,geom%nc0a
         ic0 = geom%c0a_to_c0(ic0a)
         hdata%displ_lon(ic0a,il0,its) = geom%lon(ic0)+dlon_c0a(ic0a)
         hdata%displ_lat(ic0a,il0,its) = geom%lat(ic0)+dlat_c0a(ic0a)
         call lonlatmod(hdata%displ_lon(ic0a,il0,its),hdata%displ_lat(ic0a,il0,its))
      end do
   end do
end do

! Displaced grid for timeslot 1
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      hdata%displ_lon(ic0a,il0,1) = geom%lon(ic0)
      hdata%displ_lat(ic0a,il0,1) = geom%lat(ic0)
   end do
end do

end subroutine displ_compute

!----------------------------------------------------------------------
! Subroutine: displ_write
!> Purpose: write displacement data
!----------------------------------------------------------------------
subroutine displ_write(displ,nam,geom,hdata,filename)

implicit none

! Passed variables
class(displ_type),intent(in) :: displ   !< Displacement data
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
type(hdata_type),intent(in) :: hdata    !< HDIAG data
character(len=*),intent(in) :: filename !< File name

! Local variables
integer :: ncid,nc2_id,nl0_id,nts_id,displ_niter_id,vunit_id,valid_id,dist_id,rhflt_id
integer :: lon_c2_id,lat_c2_id,lon_c2_raw_id,lat_c2_raw_id,dist_c2_raw_id,lon_c2_flt_id,lat_c2_flt_id,dist_c2_flt_id
integer :: iproc,its,il0,ic2a,ic2,i
integer,allocatable :: c2a_to_c2(:)
real(kind_real),allocatable :: sbuf(:),rbuf(:),lon_c2(:,:),lat_c2(:,:)
real(kind_real),allocatable :: lon_c2_raw(:,:,:),lat_c2_raw(:,:,:),dist_c2_raw(:,:,:)
real(kind_real),allocatable :: lon_c2_flt(:,:,:),lat_c2_flt(:,:,:),dist_c2_flt(:,:,:)
character(len=1024) :: subr = 'displ_write'

! Allocation
allocate(sbuf(hdata%nc2a*geom%nl0*(2+(nam%nts-1)*6)))

! Prepare buffer
i = 1
do il0=1,geom%nl0
   do ic2a=1,hdata%nc2a
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
   allocate(lon_c2(hdata%nc2,geom%nl0))
   allocate(lat_c2(hdata%nc2,geom%nl0))
   allocate(lon_c2_raw(hdata%nc2,geom%nl0,nam%nts-1))
   allocate(lat_c2_raw(hdata%nc2,geom%nl0,nam%nts-1))
   allocate(dist_c2_raw(hdata%nc2,geom%nl0,nam%nts-1))
   allocate(lon_c2_flt(hdata%nc2,geom%nl0,nam%nts-1))
   allocate(lat_c2_flt(hdata%nc2,geom%nl0,nam%nts-1))
   allocate(dist_c2_flt(hdata%nc2,geom%nl0,nam%nts-1))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(c2a_to_c2(hdata%proc_to_nc2a(iproc)))
      allocate(rbuf(hdata%proc_to_nc2a(iproc)*geom%nl0*(2+(nam%nts-1)*6)))

      if (iproc==mpl%ioproc) then
         ! Copy buffer
         c2a_to_c2 = hdata%c2a_to_c2
         rbuf = sbuf
      else
         ! Receive data on ioproc
         call mpl%recv(hdata%proc_to_nc2a(iproc),c2a_to_c2,iproc,mpl%tag)
         call mpl%recv(hdata%proc_to_nc2a(iproc)*geom%nl0*(2+(nam%nts-1)*6),rbuf,iproc,mpl%tag+1)
      end if

      ! Write data
      i = 1
      do il0=1,geom%nl0
         do ic2a=1,hdata%proc_to_nc2a(iproc)
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
   call mpl%send(hdata%nc2a,hdata%c2a_to_c2,mpl%ioproc,mpl%tag)
   call mpl%send(hdata%nc2a*geom%nl0*(2+(nam%nts-1)*6),sbuf,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Release memory
deallocate(sbuf)

if (mpl%main) then
   ! Create file
   call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call nam%ncwrite(ncid)

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,'nc2',hdata%nc2,nc2_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
   call ncerr(subr,nf90_def_dim(ncid,'nts',nam%nts-1,nts_id))
   call ncerr(subr,nf90_def_dim(ncid,'niter',nam%displ_niter+1,displ_niter_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,'vunit',ncfloat,(/nc2_id,nl0_id/),vunit_id))
   call ncerr(subr,nf90_def_var(ncid,'valid',ncfloat,(/displ_niter_id,nl0_id,nts_id/),valid_id))
   call ncerr(subr,nf90_put_att(ncid,valid_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'dist',ncfloat,(/displ_niter_id,nl0_id,nts_id/),dist_id))
   call ncerr(subr,nf90_put_att(ncid,dist_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'rhflt',ncfloat,(/displ_niter_id,nl0_id,nts_id/),rhflt_id))
   call ncerr(subr,nf90_put_att(ncid,rhflt_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'lon_c2',ncfloat,(/nc2_id,nl0_id/),lon_c2_id))
   call ncerr(subr,nf90_put_att(ncid,lon_c2_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'lat_c2',ncfloat,(/nc2_id,nl0_id/),lat_c2_id))
   call ncerr(subr,nf90_put_att(ncid,lat_c2_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'lon_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),lon_c2_raw_id))
   call ncerr(subr,nf90_put_att(ncid,lon_c2_raw_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'lat_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),lat_c2_raw_id))
   call ncerr(subr,nf90_put_att(ncid,lat_c2_raw_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'dist_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),dist_c2_raw_id))
   call ncerr(subr,nf90_put_att(ncid,dist_c2_raw_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'lon_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),lon_c2_flt_id))
   call ncerr(subr,nf90_put_att(ncid,lon_c2_flt_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'lat_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),lat_c2_flt_id))
   call ncerr(subr,nf90_put_att(ncid,lat_c2_flt_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'dist_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),dist_c2_flt_id))
   call ncerr(subr,nf90_put_att(ncid,dist_c2_flt_id,'_FillValue',msvalr))

   ! End definition mode
   call ncerr(subr,nf90_enddef(ncid))

   ! Write data
   call ncerr(subr,nf90_put_var(ncid,vunit_id,geom%vunit(hdata%c2_to_c0,:)))
   call ncerr(subr,nf90_put_var(ncid,valid_id,displ%valid))
   call ncerr(subr,nf90_put_var(ncid,dist_id,displ%dist*reqkm))
   call ncerr(subr,nf90_put_var(ncid,rhflt_id,displ%rhflt*reqkm))
   call ncerr(subr,nf90_put_var(ncid,lon_c2_id,lon_c2))
   call ncerr(subr,nf90_put_var(ncid,lat_c2_id,lat_c2))
   call ncerr(subr,nf90_put_var(ncid,lon_c2_raw_id,lon_c2_raw))
   call ncerr(subr,nf90_put_var(ncid,lat_c2_raw_id,lat_c2_raw))
   call ncerr(subr,nf90_put_var(ncid,dist_c2_raw_id,dist_c2_raw))
   call ncerr(subr,nf90_put_var(ncid,lon_c2_flt_id,lon_c2_flt))
   call ncerr(subr,nf90_put_var(ncid,lat_c2_flt_id,lat_c2_flt))
   call ncerr(subr,nf90_put_var(ncid,dist_c2_flt_id,dist_c2_flt))

   ! Close file
   call ncerr(subr,nf90_close(ncid))

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
