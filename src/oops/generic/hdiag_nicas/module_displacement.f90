!----------------------------------------------------------------------
! Module: module_displacement.f90
!> Purpose: displacement computation routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_displacement

use model_interface, only: model_read
use module_diag_tools, only: diag_filter
use omp_lib
use tools_const, only: req,reqkm,rad2deg,deg2rad,lonmod,sphere_dist,reduce_arc,vector_product
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_interp, only: compute_interp_bilin
use tools_missing, only: msr,isnotmsr,isallnotmsr
use tools_qsort, only: qsort
use type_ctree, only: find_nearest_neighbors
use type_displ, only: displtype,displ_alloc
use type_geom, only: fld_com_lg
use type_linop, only: apply_linop
use type_mesh, only: check_mesh
use type_mpl, only: mpl,mpl_bcast
use type_hdata, only: hdatatype
implicit none

real(kind_real),parameter :: cor_th = 0.0     !< Correlation threshold
logical,parameter :: common_sampling = .true. !< Common sampling for all variables

private
public :: compute_displacement

contains

!----------------------------------------------------------------------
! Subroutine: compute_displacement
!> Purpose: compute correlation maximum displacement
!----------------------------------------------------------------------
subroutine compute_displacement(hdata,displ,ens1)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                                                                                      !< HDIAG data
type(displtype),intent(inout) :: displ                                                                                   !< Displacement data
real(kind_real),intent(in),optional :: ens1(hdata%geom%nc0a,hdata%geom%nl0,hdata%nam%nv,hdata%nam%nts,hdata%nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ne,ne_offset,nsub,ic0,ic1,ic2,jc0,il0,isub,jsub,its,iv,ie,iter
integer,allocatable :: order(:)
real(kind_real) :: fac4,fac6,m11_avg,m2m2_avg,drhflt,lon_tmp,lat_tmp
real(kind_real),allocatable :: fld(:,:,:,:),fld_1(:,:,:),fld_2(:,:,:)
real(kind_real),allocatable :: m1_1(:,:,:,:,:,:),m2_1(:,:,:,:,:,:)
real(kind_real),allocatable :: m1_2(:,:,:,:,:,:),m2_2(:,:,:,:,:,:)
real(kind_real),allocatable :: m11(:,:,:,:,:,:)
real(kind_real),allocatable :: cor(:),cor_avg(:)
real(kind_real) :: lon(hdata%nc2,hdata%geom%nl0),lat(hdata%nc2,hdata%geom%nl0)
real(kind_real) :: dlon(hdata%nc2),dlat(hdata%nc2),dist(hdata%nc2)
real(kind_real) :: lon_nc2(hdata%nc2),lat_nc2(hdata%nc2),valid(hdata%nc2)
real(kind_real) :: dlon_ini(hdata%nc2),dlat_ini(hdata%nc2)
real(kind_real) :: lon_nc1(hdata%nam%nc1),lat_nc1(hdata%nam%nc1)
logical :: dichotomy,convergence

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Setup
ne = nam%ens1_ne
ne_offset = nam%ens1_ne_offset
nsub = nam%ens1_nsub

! Allocation
call displ_alloc(hdata,displ)
allocate(fld_1(nam%nc1,hdata%nc2,geom%nl0))
allocate(fld_2(nam%nc1,hdata%nc2,geom%nl0))
allocate(m1_1(nam%nc1,hdata%nc2,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m2_1(nam%nc1,hdata%nc2,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m1_2(nam%nc1,hdata%nc2,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m2_2(nam%nc1,hdata%nc2,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m11(nam%nc1,hdata%nc2,geom%nl0,nam%nv,2:nam%nts,nsub))

! Initialization
m1_1 = 0.0
m2_1 = 0.0
m1_2 = 0.0
m2_2 = 0.0
m11 = 0.0

! Initial point
do il0=1,geom%nl0
   do ic2=1,hdata%nc2
      if (hdata%ic1il0_log(hdata%ic2_to_ic1(ic2),il0)) then
         lon(ic2,il0) = geom%lon(hdata%ic2_to_ic0(ic2))
         lat(ic2,il0) = geom%lat(hdata%ic2_to_ic0(ic2))
      end if
   end do
end do

! Compute moments
write(mpl%unit,'(a7,a)') '','Compute moments'

! Initialization

! Loop on sub-ensembles
do isub=1,nsub
   if (nsub==1) then
      write(mpl%unit,'(a10,a)',advance='no') '','Full ensemble, member:'
   else
      write(mpl%unit,'(a10,a,i4,a)',advance='no') '','Sub-ensemble ',isub,', member:'
   end if

   ! Compute centered moments iteratively
   do ie=1,ne/nsub
      write(mpl%unit,'(i4)',advance='no') ne_offset+ie

      ! Computation factors
      fac4 = 1.0/float(ie)
      fac6 = float(ie-1)/float(ie)

      if (present(ens1)) then
         ! Copy and broadcast field
         if (allocated(fld)) deallocate(fld)
         allocate(fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))
         fld = ens1(:,:,:,:,ie+(isub-1)*nsub)
         call fld_com_lg(nam,geom,fld)
         call mpl_bcast(fld,mpl%ioproc)
      else
         ! Load field
         if (allocated(fld)) deallocate(fld)
         allocate(fld(geom%nc0,geom%nl0,nam%nv,nam%nts))
         if (nsub==1) then
            jsub = 0
         else
            jsub = isub
         end if
         call model_read(nam,geom,'ens1',ie,jsub,fld)
      end if

      do its=2,nam%nts
         do iv=1,nam%nv
            ! Initialization
            call msr(fld_1)
            call msr(fld_2)

            ! Copy valid field points
            do il0=1,geom%nl0
               do ic2=1,hdata%nc2
                  if (hdata%ic1il0_log(hdata%ic2_to_ic1(ic2),il0)) then
                     !$omp parallel do schedule(static) private(ic1,ic0,jc0)
                     do ic1=1,nam%nc1
                        if (hdata%displ_mask(ic1,ic2,min(il0,geom%nl0i))) then
                           ! Indices
                           ic0 = hdata%ic1_to_ic0(ic1)
                           jc0 = hdata%ic2_to_ic0(ic2)

                           ! Copy points
                           fld_1(ic1,ic2,il0) = fld(jc0,il0,iv,1)
                           fld_2(ic1,ic2,il0) = fld(ic0,il0,iv,its)
                        end if
                     end do
                     !$omp end parallel do
                  end if
               end do
            end do

            ! Remove means
            fld_1 = fld_1 - m1_1(:,:,:,iv,its,isub)
            fld_2 = fld_2 - m1_2(:,:,:,iv,its,isub)

            ! Update high-order moments
            if (ie>1) then
               ! Covariance
               m11(:,:,:,iv,its,isub) = m11(:,:,:,iv,its,isub)+fac6*fld_1*fld_2

               ! Variances
               m2_1(:,:,:,iv,its,isub) = m2_1(:,:,:,iv,its,isub)+fac6*fld_1**2
               m2_2(:,:,:,iv,its,isub) = m2_2(:,:,:,iv,its,isub)+fac6*fld_2**2
            end if

            ! Update means
            m1_1(:,:,:,iv,its,isub) = m1_1(:,:,:,iv,its,isub)+fac4*fld_1
            m1_2(:,:,:,iv,its,isub) = m1_2(:,:,:,iv,its,isub)+fac4*fld_2
         end do
      end do
   end do
   write(mpl%unit,'(a)') ''

   do its=2,nam%nts
      do iv=1,nam%nv
         ! Normalize
         m2_1(:,:,:,iv,its,isub) = m2_1(:,:,:,iv,its,isub)/float(ne/nsub-1)
         m2_2(:,:,:,iv,its,isub) = m2_2(:,:,:,iv,its,isub)/float(ne/nsub-1)
         m11(:,:,:,iv,its,isub) = m11(:,:,:,iv,its,isub)/float(ne/nsub-1)
      end do
   end do
end do

! Find correlation maximum propagation
write(mpl%unit,'(a7,a)') '','Find correlation maximum propagation'

do its=2,nam%nts
   do il0=1,geom%nl0
      write(mpl%unit,'(a10,a,i2,a,i3)') '','Timeslot ',its,' - level ',nam%levs(il0)
      !$omp parallel do schedule(static) private(ic2,cor,cor_avg,order,ic1,iv,m11_avg,m2m2_avg)
      do ic2=1,hdata%nc2
         if (hdata%ic1il0_log(hdata%ic2_to_ic1(ic2),il0)) then
            ! Allocation
            allocate(cor(nam%nv))
            allocate(cor_avg(nam%nc1))
            allocate(order(nam%nc1))

            do ic1=1,nam%nc1
               if (hdata%displ_mask(ic1,ic2,min(il0,geom%nl0i))) then
                  ! Compute correlation for each variable
                  do iv=1,nam%nv
                     ! Correlation
                     m11_avg = sum(m11(ic1,ic2,il0,iv,its,:))/float(nsub)
                     m2m2_avg = sum(m2_1(ic1,ic2,il0,iv,its,:))*sum(m2_2(ic1,ic2,il0,iv,its,:))/float(nsub**2)
                     if (m2m2_avg>0.0) then
                        cor(iv) = m11_avg/sqrt(m2m2_avg)
                     else
                        call msr(cor(iv))
                     end if
                  end do

                  ! Average correlations
                  cor_avg(ic1) = sum(cor,mask=isnotmsr(cor))/float(count(isnotmsr(cor)))
               else
                  call msr(cor_avg(ic1))
               end if
            end do

            ! Sort correlations
            call qsort(nam%nc1,cor_avg,order)

            ! Locate the maximum correlation, with a correlation threshold
            if (cor_avg(nam%nc1)>cor_th) then
               dlon(ic2) = lonmod(geom%lon(hdata%ic1_to_ic0(order(nam%nc1)))-lon(ic2,il0))
               dlat(ic2) = geom%lat(hdata%ic1_to_ic0(order(nam%nc1)))-lat(ic2,il0)
               call sphere_dist(lon(ic2,il0),lat(ic2,il0),geom%lon(hdata%ic1_to_ic0(order(nam%nc1))), &
             & geom%lat(hdata%ic1_to_ic0(order(nam%nc1))),dist(ic2))
            end if

            ! Release memory
            deallocate(cor)
            deallocate(cor_avg)
            deallocate(order)
         end if
      end do
      !$omp end parallel do

      ! Average distance
      displ%dist(0,il0,its) = sum(dist,mask=hdata%ic1il0_log(hdata%ic2_to_ic1,il0)) &
                            & /count(hdata%ic1il0_log(hdata%ic2_to_ic1,il0))

      ! Check raw mesh
      do ic2=1,hdata%nc2
         lon_nc2(ic2) = lonmod(lon(ic2,il0)+dlon(ic2))
         lat_nc2(ic2) = lat(ic2,il0)+dlat(ic2)
      end do
      call check_mesh(hdata%nc2,lon_nc2,lat_nc2,hdata%nt,hdata%ltri,valid)
      displ%valid(0,il0,its) = sum(valid,mask=hdata%ic1il0_log(hdata%ic2_to_ic1,il0)) &
                             & /count(hdata%ic1il0_log(hdata%ic2_to_ic1,il0))
      displ%rhflt(0,il0,its) = 0.0

      ! Raw displacement interpolation
      call apply_linop(hdata%h(min(il0,geom%nl0i)),dlon,displ%dlon_raw(:,il0,its))
      call apply_linop(hdata%h(min(il0,geom%nl0i)),dlat,displ%dlat_raw(:,il0,its))

      if (nam%displ_niter>0) then
         ! Filter displacement

         ! Allocation
         allocate(order(hdata%nc2))

         ! Compute displacement
         dlon_ini = dlon*cos(lat(:,il0))
         dlat_ini = dlat

         ! Iterative filtering
         convergence = .true.
         drhflt = 0.0
         dichotomy = .false.

         ! Dichotomy initialization
         displ%rhflt(1,il0,its) = nam%displ_rhflt

         do iter=1,nam%displ_niter
            ! Copy increment
            dlon = dlon_ini
            dlat = dlat_ini

            ! Median filter to remove extreme values
            call diag_filter(hdata,il0,'median',displ%rhflt(iter,il0,its),dlon)
            call diag_filter(hdata,il0,'median',displ%rhflt(iter,il0,its),dlat)

            ! Average filter to smooth displacement
            call diag_filter(hdata,il0,'average',displ%rhflt(iter,il0,its),dlon)
            call diag_filter(hdata,il0,'average',displ%rhflt(iter,il0,its),dlat)

            ! Compute displaced location
            dlon = dlon/cos(lat(:,il0))

            ! Reduce distance with respect to boundary
            do ic2=1,hdata%nc2
               if (hdata%ic1il0_log(hdata%ic2_to_ic1(ic2),il0)) then
                  lon_tmp = lonmod(lon(ic2,il0)+dlon(ic2))
                  lat_tmp = lat(ic2,il0)+dlat(ic2)
                  call reduce_arc(lon(ic2,il0),lat(ic2,il0),lon_tmp,lat_tmp,hdata%bdist(ic2),dist(ic2))
                  dlon(ic2) = lonmod(lon_tmp-lon(ic2,il0))
                  dlat(ic2) = lat_tmp-lat(ic2,il0)
               end if
            end do

            ! Check mesh
            do ic2=1,hdata%nc2
               lon_nc2(ic2) = lonmod(lon(ic2,il0)+dlon(ic2))
               lat_nc2(ic2) = lat(ic2,il0)+dlat(ic2)
            end do
            call check_mesh(hdata%nc2,lon_nc2,lat_nc2,hdata%nt,hdata%ltri,valid)
            displ%valid(iter,il0,its) = sum(valid,mask=hdata%ic1il0_log(hdata%ic2_to_ic1,il0)) &
                                      & /count(hdata%ic1il0_log(hdata%ic2_to_ic1,il0))

            ! Compute distances
            do ic2=1,hdata%nc2
               if (hdata%ic1il0_log(hdata%ic2_to_ic1(ic2),il0)) call sphere_dist(lon(ic2,il0),lat(ic2,il0), &
             & lonmod(lon(ic2,il0)+dlon(ic2)),lat(ic2,il0)+dlat(ic2),dist(ic2))
            end do
            displ%dist(iter,il0,its) = sum(dist,mask=hdata%ic1il0_log(hdata%ic2_to_ic1,il0)) &
                                     & /count(hdata%ic1il0_log(hdata%ic2_to_ic1,il0))

            ! Print results
            write(mpl%unit,'(a13,a,i2,a,f8.2,a,f6.2,a,f6.2,a,f7.2,a)') '','Iteration ',iter,': rhflt = ', &
          & displ%rhflt(iter,il0,its)*reqkm,' km, valid points: ',100.0*displ%valid(0,il0,its),'% ~> ', &
          & 100.0*displ%valid(iter,il0,its),'%, average displacement = ',displ%dist(iter,il0,its)*reqkm,' km'

            ! Update support radius
            if (displ%valid(iter,il0,its)<1.0-nam%displ_tol) then
               ! Increase filtering support radius
               if (dichotomy) then
                  if (iter<nam%displ_niter) displ%rhflt(iter+1,il0,its) = displ%rhflt(iter,il0,its)+drhflt
               else
                  ! No convergence
                  convergence = .false.
                  if (iter<nam%displ_niter) displ%rhflt(iter+1,il0,its) = 2.0*displ%rhflt(iter,il0,its)
               end if
            else
               ! Convergence
               convergence = .true.

               ! Check dichotomy status
               if (.not.dichotomy) then
                  dichotomy = .true.
                  drhflt = 0.5*displ%rhflt(iter,il0,its)
               end if

               ! Decrease filtering support radius
               if (iter<nam%displ_niter) displ%rhflt(iter+1,il0,its) = displ%rhflt(iter,il0,its)-drhflt
            end if
            if (dichotomy) drhflt = 0.5*drhflt
         end do

         ! Check convergence
         if (.not.convergence) call msgerror('iterative filtering failed')

         ! Deallocation
         deallocate(order)
      end if

      ! Filtered displacement interpolation
      call apply_linop(hdata%h(min(il0,geom%nl0i)),dlon,displ%dlon_flt(:,il0,its))
      call apply_linop(hdata%h(min(il0,geom%nl0i)),dlat,displ%dlat_flt(:,il0,its))

      ! Initialize sampling interpolation
      displ%d(il0,its)%prefix = 'd'
      displ%d(il0,its)%n_src = geom%nc0
      displ%d(il0,its)%n_dst = nam%nc1

      ! Compute sampling interpolation
      do ic1=1,nam%nc1
         lon_nc1(ic1) = lonmod(geom%lon(hdata%ic1_to_ic0(ic1))+displ%dlon_flt(hdata%ic1_to_ic0(ic1),il0,its))
         lat_nc1(ic1) = geom%lat(hdata%ic1_to_ic0(ic1))+displ%dlat_flt(hdata%ic1_to_ic0(ic1),il0,its)
      end do
      call compute_interp_bilin(geom%mesh,geom%ctree(min(il0,geom%nl0i)),geom%nc0,geom%mask(:,il0), &
    & nam%nc1,lon_nc1,lat_nc1,hdata%ic1il0_log(:,il0),displ%d(il0,its))
   end do
end do

! Compute dummy sampling interpolation (timeslot 1)
do il0=1,geom%nl0
   ! Initialize interpolation
   displ%d(il0,1)%prefix = 'd'
   displ%d(il0,1)%n_src = geom%nc0
   displ%d(il0,1)%n_dst = nam%nc1

   ! Compute interpolation
   call compute_interp_bilin(geom%mesh,geom%ctree(min(il0,geom%nl0i)),geom%nc0,geom%mask(:,il0), &
 & nam%nc1,geom%lon(hdata%ic1_to_ic0),geom%lat(hdata%ic1_to_ic0),hdata%ic1il0_log(:,il0),displ%d(il0,1))
end do

! Output units
displ%dlon_raw = displ%dlon_raw*rad2deg
displ%dlat_raw = displ%dlat_raw*rad2deg
displ%dlon_flt = displ%dlon_flt*rad2deg
displ%dlat_flt = displ%dlat_flt*rad2deg
displ%dist = displ%dist*reqkm
displ%rhflt = displ%rhflt*reqkm

! End associate
end associate

end subroutine compute_displacement

end module module_displacement
