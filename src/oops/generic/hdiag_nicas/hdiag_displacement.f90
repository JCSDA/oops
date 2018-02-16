!----------------------------------------------------------------------
! Module: hdiag_displacement.f90
!> Purpose: displacement computation routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_displacement

use model_interface, only: model_read
use omp_lib
use tools_const, only: req,reqkm,rad2deg,deg2rad,lonlatmod,sphere_dist,reduce_arc,vector_product
use hdiag_tools, only: diag_filter,diag_com_lg
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_qsort, only: qsort
use tools_stripack, only: trans,scoord
use type_com, only: com_ext
use type_displ, only: displtype,displ_alloc
use type_geom, only: fld_com_gl
use type_linop, only: apply_linop
use type_mesh, only: check_mesh
use type_mpl, only: mpl,mpl_bcast
use type_hdata, only: hdatatype
implicit none

real(kind_real),parameter :: cor_th = 0.2     !< Correlation threshold

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
integer :: ne,ne_offset,nsub,ic0,ic1,ic2,ic2a,jc0,jc1,il0,isub,jsub,iv,its,ie,iter,ic0a,jc0d
integer,allocatable :: order(:)
real(kind_real) :: fac4,fac6,m11_avg,m2m2_avg,fld_1,fld_2,drhflt,dum
real(kind_real),allocatable :: fld(:,:,:,:),fld_halo(:,:,:,:),fld_com(:,:)
real(kind_real),allocatable :: m1_1(:,:,:,:,:,:),m2_1(:,:,:,:,:,:)
real(kind_real),allocatable :: m1_2(:,:,:,:,:,:),m2_2(:,:,:,:,:,:)
real(kind_real),allocatable :: m11(:,:,:,:,:,:)
real(kind_real),allocatable :: cor(:),cor_avg(:)
real(kind_real),allocatable :: dlon_c2(:),dlat_c2(:),dist_c2(:)
real(kind_real) :: lon_c2_ori(hdata%nc2,hdata%geom%nl0),lat_c2_ori(hdata%nc2,hdata%geom%nl0)
real(kind_real) :: lon_c2(hdata%nc2),lat_c2(hdata%nc2),valid(hdata%nc2)
real(kind_real) :: x_ori(hdata%nc2),y_ori(hdata%nc2),z_ori(hdata%nc2)
real(kind_real) :: dx_ini(hdata%nc2),dy_ini(hdata%nc2),dz_ini(hdata%nc2)
real(kind_real) :: dx(hdata%nc2),dy(hdata%nc2),dz(hdata%nc2)
real(kind_real) :: dlon_c0(hdata%geom%nc0),dlat_c0(hdata%geom%nc0)
logical :: dichotomy,convergence

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Setup
ne = nam%ens1_ne
ne_offset = nam%ens1_ne_offset
nsub = nam%ens1_nsub

! Allocation
call displ_alloc(hdata,displ)
allocate(fld_halo(hdata%nc0d,geom%nl0,nam%nv,2:nam%nts))
allocate(m1_1(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m2_1(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m1_2(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m2_2(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,nsub))
allocate(m11(nam%nc1,hdata%nc2a,geom%nl0,nam%nv,2:nam%nts,nsub))

! Initialization
m1_1 = 0.0
m2_1 = 0.0
m1_2 = 0.0
m2_2 = 0.0
m11 = 0.0

! Initial point
do il0=1,geom%nl0
   do ic2=1,hdata%nc2
      ic1 = hdata%c2_to_c1(ic2)
      ic0 = hdata%c2_to_c0(ic2)
      if (hdata%c1l0_log(ic1,il0)) then
         lon_c2_ori(ic2,il0) = geom%lon(ic0)
         lat_c2_ori(ic2,il0) = geom%lat(ic0)
      end if
   end do
end do

! Copy
displ%lon_c2 = lon_c2_ori
displ%lat_c2 = lat_c2_ori

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
         ! Copy field
         if (allocated(fld)) deallocate(fld)
         allocate(fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))
         fld = ens1(:,:,:,:,ie+(isub-1)*nsub)
      else
         ! Load field
         if (mpl%main) then
            if (allocated(fld)) deallocate(fld)
            allocate(fld(geom%nc0,geom%nl0,nam%nv,nam%nts))
            if (nsub==1) then
               jsub = 0
            else
               jsub = isub
            end if
            call model_read(nam,geom,'ens1',ne_offset+ie,jsub,fld)
         end if
         call fld_com_gl(nam,geom,fld)
      end if

      do its=2,nam%nts
         do iv=1,nam%nv
            ! Allocation
            allocate(fld_com(geom%nc0a,geom%nl0))

            ! Copy points
            fld_com = fld(:,:,iv,its)

            ! Halo extension
            call com_ext(hdata%AD,fld_com)

            ! Copy points
            fld_halo(:,:,iv,its) = fld_com

            ! Release memory
            deallocate(fld_com)
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
                           fld_1 = fld(ic0a,il0,iv,1)
                           fld_2 = fld_halo(jc0d,il0,iv,its)

                           ! Remove means
                           fld_1 = fld_1 - m1_1(jc1,ic2a,il0,iv,its,isub)
                           fld_2 = fld_2 - m1_2(jc1,ic2a,il0,iv,its,isub)

                           ! Update high-order moments
                           if (ie>1) then
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
end do

! Find correlation maximum propagation
write(mpl%unit,'(a7,a)') '','Find correlation maximum propagation'

do its=2,nam%nts
   do il0=1,geom%nl0
      write(mpl%unit,'(a10,a,i2,a,i3)') '','Timeslot ',its,' - level ',nam%levs(il0)

      ! Allocation
      allocate(dlon_c2(hdata%nc2a))
      allocate(dlat_c2(hdata%nc2a))
      allocate(dist_c2(hdata%nc2a))

      !$omp parallel do schedule(static) private(ic2a,ic2,ic1,cor,cor_avg,order,jc1,jc0,iv,m11_avg,m2m2_avg)
      do ic2a=1,hdata%nc2a
         ic2 = hdata%c2a_to_c2(ic2a)
         ic1 = hdata%c2_to_c1(ic2)
         if (hdata%c1l0_log(ic1,il0)) then
            ! Allocation
            allocate(cor(nam%nv))
            allocate(cor_avg(nam%nc1))
            allocate(order(nam%nc1))

            do jc1=1,nam%nc1
               ! Initialization
               call msr(cor_avg(jc1))

               if (hdata%displ_mask(jc1,ic2,min(il0,geom%nl0i))) then
                  ! Compute correlation for each variable
                  do iv=1,nam%nv
                     ! Correlation
                     m11_avg = sum(m11(jc1,ic2a,il0,iv,its,:))/float(nsub)
                     m2m2_avg = sum(m2_1(jc1,ic2a,il0,iv,its,:))*sum(m2_2(jc1,ic2a,il0,iv,its,:))/float(nsub**2)
                     if (m2m2_avg>0.0) then
                        cor(iv) = m11_avg/sqrt(m2m2_avg)
                     else
                        call msr(cor(iv))
                     end if
                  end do

                  ! Average correlations
                  if (isanynotmsr(cor)) then
                     cor_avg(jc1) = sum(cor,mask=isnotmsr(cor))/float(count(isnotmsr(cor)))
                  else
                     call msgerror('average correlation contains missing values only')
                  end if
               end if
            end do

            ! Sort correlations
            call qsort(nam%nc1,cor_avg,order)

            ! Locate the maximum correlation, with a correlation threshold
            if (cor_avg(nam%nc1)>cor_th) then
               jc1 = order(nam%nc1)
               jc0 = hdata%c1_to_c0(jc1)
               dlon_c2(ic2a) = geom%lon(jc0)-lon_c2_ori(ic2,il0)
               dlat_c2(ic2a) = geom%lat(jc0)-lat_c2_ori(ic2,il0)
               call lonlatmod(dlon_c2(ic2a),dlat_c2(ic2a))
               call sphere_dist(lon_c2_ori(ic2,il0),lat_c2_ori(ic2,il0),geom%lon(jc0),geom%lat(jc0),dist_c2(ic2a))
            else
               dlon_c2(ic2a) = 0.0
               dlat_c2(ic2a) = 0.0
               dist_c2(ic2a) = 0.0
            end if

            ! Release memory
            deallocate(cor)
            deallocate(cor_avg)
            deallocate(order)
         end if
      end do
      !$omp end parallel do

      ! Gather data
      call diag_com_lg(hdata,dlon_c2)
      call diag_com_lg(hdata,dlat_c2)
      call diag_com_lg(hdata,dist_c2)

      if (.not.mpl%main) then
         ! Allocation
         allocate(dlon_c2(hdata%nc2))
         allocate(dlat_c2(hdata%nc2))
         allocate(dist_c2(hdata%nc2))
      end if

      ! Broadcast data
      call mpl_bcast(dlon_c2,mpl%ioproc)
      call mpl_bcast(dlat_c2,mpl%ioproc)
      call mpl_bcast(dist_c2,mpl%ioproc)

      ! Average distance
      displ%dist(0,il0,its) = sum(dist_c2,mask=hdata%c1l0_log(hdata%c2_to_c1,il0)) &
                            & /count(hdata%c1l0_log(hdata%c2_to_c1,il0))

      ! Check raw mesh
      do ic2=1,hdata%nc2
         lon_c2(ic2) = lon_c2_ori(ic2,il0)+dlon_c2(ic2)
         lat_c2(ic2) = lat_c2_ori(ic2,il0)+dlat_c2(ic2)
         call lonlatmod(lon_c2(ic2),lat_c2(ic2))
      end do
      call check_mesh(hdata%nc2,lon_c2,lat_c2,hdata%nt,hdata%ltri,valid)
      displ%valid(0,il0,its) = sum(valid,mask=hdata%c1l0_log(hdata%c2_to_c1,il0)) &
                             & /count(hdata%c1l0_log(hdata%c2_to_c1,il0))
      displ%rhflt(0,il0,its) = 0.0

      ! Copy
      displ%lon_c2_raw(:,il0,its) = lon_c2
      displ%lat_c2_raw(:,il0,its) = lat_c2
      displ%dist_c2_raw(:,il0,its) = dist_c2

      if (nam%displ_niter>0) then
         ! Filter displacement

         ! Compute raw displacement in cartesian coordinates
         call trans(hdata%nc2,lat_c2_ori(:,il0),lon_c2_ori(:,il0),x_ori,y_ori,z_ori)
         call trans(hdata%nc2,lat_c2,lon_c2,dx_ini,dy_ini,dz_ini)
         dx_ini = dx_ini-x_ori
         dy_ini = dy_ini-y_ori
         dz_ini = dz_ini-z_ori

         ! Iterative filtering
         convergence = .true.
         drhflt = 0.0
         dichotomy = .false.

         ! Dichotomy initialization
         displ%rhflt(1,il0,its) = nam%displ_rhflt

         do iter=1,nam%displ_niter
            ! Copy increment
            dx = dx_ini
            dy = dy_ini
            dz = dz_ini

            ! Median filter to remove extreme values
            call diag_filter(hdata,il0,'median',displ%rhflt(iter,il0,its),dx)
            call diag_filter(hdata,il0,'median',displ%rhflt(iter,il0,its),dy)
            call diag_filter(hdata,il0,'median',displ%rhflt(iter,il0,its),dz)

            ! Average filter to smooth displacement
            call diag_filter(hdata,il0,'gc99',displ%rhflt(iter,il0,its),dx)
            call diag_filter(hdata,il0,'gc99',displ%rhflt(iter,il0,its),dy)
            call diag_filter(hdata,il0,'gc99',displ%rhflt(iter,il0,its),dz)

            ! Back to spherical coordinates
            dx = dx+x_ori
            dy = dy+y_ori
            dz = dz+z_ori
            do ic2=1,hdata%nc2
               call scoord(dx(ic2),dy(ic2),dz(ic2),lat_c2(ic2),lon_c2(ic2),dum)
            end do

            ! Reduce distance with respect to boundary
            do ic2=1,hdata%nc2
               if (hdata%c1l0_log(hdata%c2_to_c1(ic2),il0)) then
                  call reduce_arc(lon_c2_ori(ic2,il0),lat_c2_ori(ic2,il0),lon_c2(ic2),lat_c2(ic2),hdata%bdist(ic2),dist_c2(ic2))
                  dlon_c2(ic2) = lon_c2(ic2)-lon_c2_ori(ic2,il0)
                  dlat_c2(ic2) = lat_c2(ic2)-lat_c2_ori(ic2,il0)
               end if
            end do

            ! Check mesh
            do ic2=1,hdata%nc2
               lon_c2(ic2) = lon_c2_ori(ic2,il0)+dlon_c2(ic2)
               lat_c2(ic2) = lat_c2_ori(ic2,il0)+dlat_c2(ic2)
               call lonlatmod(lon_c2(ic2),lat_c2(ic2))
            end do
            call check_mesh(hdata%nc2,lon_c2,lat_c2,hdata%nt,hdata%ltri,valid)
            displ%valid(iter,il0,its) = sum(valid,mask=hdata%c1l0_log(hdata%c2_to_c1,il0)) &
                                      & /count(hdata%c1l0_log(hdata%c2_to_c1,il0))

            ! Compute distances
            do ic2=1,hdata%nc2
               if (hdata%c1l0_log(hdata%c2_to_c1(ic2),il0)) call sphere_dist(lon_c2_ori(ic2,il0),lat_c2_ori(ic2,il0), &
             & lon_c2(ic2),lat_c2(ic2),dist_c2(ic2))
            end do
            displ%dist(iter,il0,its) = sum(dist_c2,mask=hdata%c1l0_log(hdata%c2_to_c1,il0)) &
                                     & /count(hdata%c1l0_log(hdata%c2_to_c1,il0))

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

         ! Copy
         displ%lon_c2_flt(:,il0,its) = lon_c2
         displ%lat_c2_flt(:,il0,its) = lat_c2
         displ%dist_c2_flt(:,il0,its) = dist_c2

         ! Check convergence
         if (.not.convergence) call msgerror('iterative filtering failed')
      else
         ! Print results
         write(mpl%unit,'(a10,a22,f8.2,a,f6.2,a,f7.2,a)') '','Raw displacement: rhflt = ', &
       & displ%rhflt(0,il0,its)*reqkm,' km, valid points: ',100.0*displ%valid(0,il0,its),'%, average displacement = ', &
       & displ%dist(0,il0,its)*reqkm,' km'
      end if

      ! Displacement interpolation
      do ic2=1,hdata%nc2
         call lonlatmod(dlon_c2(ic2),dlat_c2(ic2))
      end do
      call apply_linop(hdata%h(min(il0,geom%nl0i)),dlon_c2,dlon_c0)
      call apply_linop(hdata%h(min(il0,geom%nl0i)),dlat_c2,dlat_c0)

      ! Displaced grid
      displ%lon_c0_flt(:,il0,its) = geom%lon+dlon_c0
      displ%lat_c0_flt(:,il0,its) = geom%lat+dlat_c0
      do ic0=1,geom%nc0
         call lonlatmod(displ%lon_c0_flt(ic0,il0,its),displ%lat_c0_flt(ic0,il0,its))
      end do

      ! Release memory
      deallocate(dlon_c2)
      deallocate(dlat_c2)
      deallocate(dist_c2)
   end do
end do

! Displaced grid for timeslot 1
do il0=1,geom%nl0
   displ%lon_c0_flt(:,il0,1) = geom%lon
   displ%lat_c0_flt(:,il0,1) = geom%lat
end do

! End associate
end associate

end subroutine compute_displacement

end module hdiag_displacement
