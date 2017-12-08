!----------------------------------------------------------------------
! Module: module_sampling.f90
!> Purpose: sampling routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_sampling

use omp_lib
use tools_const, only: pi,req,deg2rad,sphere_dist,vector_product,vector_triple_product
use tools_display, only: prog_init,prog_print,msgerror,msgwarning,black,green,peach
use tools_interp, only: compute_grid_interp_bilin
use tools_kinds, only: kind_real
use tools_missing, only: msi,isnotmsi
use tools_stripack, only: trans,trmesh,trlist,bnodes,scoord
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_hdata, only: hdatatype,hdata_alloc,hdata_read,hdata_write
use type_linop, only: linop_alloc
use type_mpl, only: mpl,mpl_bcast
use type_randgen, only: rand_integer,initialize_sampling
implicit none

integer,parameter :: irmax = 10000 !< Maximum number of random number draws

private
public :: setup_sampling,compute_sampling_zs,compute_sampling_lct,compute_sampling_mesh

contains

!----------------------------------------------------------------------
! Subroutine: setup_sampling
!> Purpose: setup sampling
!----------------------------------------------------------------------
subroutine setup_sampling(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: info,il0,ic0,ic1,ic2,ildw,ic,i_s,il0i,jc1
integer :: mask_ind(hdata%nam%nc1)
integer,allocatable :: vbot(:),vtop(:),nn_nc1_index(:)
real(kind_real) :: rh0(hdata%geom%nc0,hdata%geom%nl0),dum(1)
real(kind_real),allocatable :: nn_nc1_dist(:)
type(ctreetype) :: ctree_diag

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Check subsampling size
if (nam%nc1>maxval(count(geom%mask,dim=1))) then
   call msgwarning('nc1 is too large for then mask, reset nc1 to the largest possible value')
   nam%nc1 = maxval(count(geom%mask,dim=1))
end if

! Define nc2
if (nam%new_lct) then
   hdata%nc2 = nam%nc1
else
   if (nam%local_diag.or.nam%displ_diag) then
      hdata%nc2 = int(2.0*maxval(geom%area)/(sqrt(3.0)*(nam%local_rad)**2))
      write(mpl%unit,'(a7,a,i8)') '','Estimated nc2 from local diagnostic radius: ',hdata%nc2
      hdata%nc2 = min(hdata%nc2,nam%nc1)
      write(mpl%unit,'(a7,a,i8)') '','Final nc2: ',hdata%nc2
   end if
end if

! Allocation
call hdata_alloc(hdata)

! Read or compute sampling data
info = 1
if (nam%sam_read) info = hdata_read(hdata)
if (info==1) then
   ! Compute zero-separation sampling
   call compute_sampling_zs(hdata)

   if (nam%new_lct) then
      ! Compute LCT sampling
      call compute_sampling_lct(hdata)
   else
      ! Compute positive separation sampling
      call compute_sampling_ps(hdata)
   end if
end if

if (nam%local_diag.or.nam%displ_diag) then
   if ((info==1).or.(info==2)) then
      ! Define subsampling
      mask_ind = 1
      rh0 = 1.0
      call initialize_sampling(nam%nc1,dble(geom%lon(hdata%ic1_to_ic0)),dble(geom%lat(hdata%ic1_to_ic0)),mask_ind, &
    & rh0,nam%ntry,nam%nrep,hdata%nc2,hdata%ic2_to_ic1)
      hdata%ic2_to_ic0 = hdata%ic1_to_ic0(hdata%ic2_to_ic1)
   end if

   if ((info==1).or.(info==2).or.(info==3).or.(info==4)) then
      if (trim(nam%flt_type)/='none') then
         do il0=1,geom%nl0
            if ((il0==1).or.(geom%nl0i>1)) then
               write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(il0)
               ctree_diag = create_ctree(hdata%nc2,geom%lon(hdata%ic2_to_ic0),geom%lat(hdata%ic2_to_ic0), &
                          & hdata%ic1il0_log(hdata%ic2_to_ic1,il0))
               do ic2=1,hdata%nc2
                  ic1 = hdata%ic2_to_ic1(ic2)
                  ic0 = hdata%ic2_to_ic0(ic2)
                  if (hdata%ic1il0_log(ic1,il0)) call find_nearest_neighbors(ctree_diag,geom%lon(ic0),geom%lat(ic0), &
                   & hdata%nc2,hdata%nn_nc2_index(:,ic2,il0),hdata%nn_nc2_dist(:,ic2,il0))
               end do
               call delete_ctree(ctree_diag)
            end if
         end do
      end if
   end if

   ! Compute sampling mesh
   write(mpl%unit,'(a7,a)') '','Compute sampling mesh'
   call compute_sampling_mesh(hdata)

   if ((info==1).or.(info==2).or.(info==3)) then
      ! Allocation
      allocate(nn_nc1_index(nam%nc1))
      allocate(nn_nc1_dist(nam%nc1))
      allocate(vbot(hdata%nc2))
      allocate(vtop(hdata%nc2))

      ! Compute nearest neighbors
      write(mpl%unit,'(a7,a)') '','Compute nearest neighbors'
      do il0i=1,geom%nl0i
         write(mpl%unit,'(a10,a,i3)') '','Independent level ',il0i
         ctree_diag = create_ctree(nam%nc1,geom%lon(hdata%ic1_to_ic0),geom%lat(hdata%ic1_to_ic0),hdata%ic1il0_log(:,il0i))
         do ic2=1,hdata%nc2
            ic0 = hdata%ic2_to_ic0(ic2)
            if (hdata%ic1il0_log(hdata%ic2_to_ic1(ic2),il0i)) then
               ! Find nearest neighbors
               call find_nearest_neighbors(ctree_diag,geom%lon(ic0),geom%lat(ic0), &
             & nam%nc1,nn_nc1_index,nn_nc1_dist)

               do ic1=1,nam%nc1
                  jc1 = nn_nc1_index(ic1)
                  hdata%local_mask(jc1,ic2,il0i) = (ic1==1).or.(nn_nc1_dist(ic1)<min(nam%local_rad,hdata%bdist(ic2)))
                  hdata%displ_mask(jc1,ic2,il0i) = (ic1==1).or.(nn_nc1_dist(ic1)<min(nam%displ_rad,hdata%bdist(ic2)))
               end do
            end if
         end do
         call delete_ctree(ctree_diag)
      end do

      ! Initialize vbot and vtop
      vbot = 1
      vtop = geom%nl0

      ! Compute grid interpolation to Sc0
      call compute_grid_interp_bilin(geom,hdata%nc2,hdata%ic2_to_ic0,nam%mask_check,vbot,vtop,hdata%h)

      ! Compute grid interpolation to Sc1
      do il0i=1,geom%nl0i
         hdata%s(il0i)%prefix = 's'
         hdata%s(il0i)%n_src = hdata%nc2
         hdata%s(il0i)%n_dst = nam%nc1
         hdata%s(il0i)%n_s = 0
         do i_s=1,hdata%h(il0i)%n_s
            do ic1=1,nam%nc1
               if (hdata%ic1_to_ic0(ic1)==hdata%h(il0i)%row(i_s)) hdata%s(il0i)%n_s = hdata%s(il0i)%n_s+1
            end do
         end do
         call linop_alloc(hdata%s(il0i))
         hdata%s(il0i)%n_s = 0
         do i_s=1,hdata%h(il0i)%n_s
            do ic1=1,nam%nc1
               if (hdata%ic1_to_ic0(ic1)==hdata%h(il0i)%row(i_s)) then
                  hdata%s(il0i)%n_s = hdata%s(il0i)%n_s+1
                  hdata%s(il0i)%row(hdata%s(il0i)%n_s) = ic1
                  hdata%s(il0i)%col(hdata%s(il0i)%n_s) = hdata%h(il0i)%col(i_s)
                  hdata%s(il0i)%S(hdata%s(il0i)%n_s) = hdata%h(il0i)%S(i_s)
               end if
            end do
         end do
      end do

      ! Release memory
      deallocate(nn_nc1_index)
      deallocate(nn_nc1_dist)
      deallocate(vbot)
      deallocate(vtop)
   end if
end if

! Write sampling data
if (nam%sam_write.and.mpl%main) call hdata_write(hdata)

! Compute nearest neighbors for local diagnostics output
if (nam%local_diag.and.(nam%nldwv>0)) then
   write(mpl%unit,'(a7,a)') '','Compute nearest neighbors for local diagnostics output'
   allocate(hdata%nn_ldwv_index(nam%nldwv))
   ctree_diag = create_ctree(hdata%nc2,geom%lon(hdata%ic2_to_ic0), &
                geom%lat(hdata%ic2_to_ic0),hdata%ic1il0_log(hdata%ic2_to_ic1,1))
   do ildw=1,nam%nldwv
      call find_nearest_neighbors(ctree_diag,nam%lon_ldwv(ildw)*deg2rad,nam%lat_ldwv(ildw)*deg2rad, &
    & 1,hdata%nn_ldwv_index(ildw:ildw),dum)
   end do
   call delete_ctree(ctree_diag)
end if

! Print results
write(mpl%unit,'(a7,a)') '','Sampling efficiency (%):'
do il0=1,geom%nl0
   write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),' ~> '
   do ic=1,nam%nc
      if (count(hdata%ic1icil0_log(:,ic,il0))>=nam%nc1/2) then
         ! Sucessful sampling
         write(mpl%unit,'(a,i3,a)',advance='no') trim(green), &
       & int(100.0*float(count(hdata%ic1icil0_log(:,ic,il0)))/float(nam%nc1)),trim(black)
      else
         ! Insufficient sampling
         write(mpl%unit,'(a,i3,a)',advance='no') trim(peach), &
       & int(100.0*float(count(hdata%ic1icil0_log(:,ic,il0)))/float(nam%nc1)),trim(black)
      end if
      if (ic<nam%nc) write(mpl%unit,'(a)',advance='no') '-'
   end do
   write(mpl%unit,'(a)') ' '
end do

! End associate
end associate

end subroutine setup_sampling

!----------------------------------------------------------------------
! Subroutine: compute_sampling_zs
!> Purpose: compute zero-separation sampling
!----------------------------------------------------------------------
subroutine compute_sampling_zs(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ic0,ic1,il0
integer :: mask_ind_col(hdata%geom%nc0)
real(kind_real) :: rh0(hdata%geom%nc0)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Initialize mask
mask_ind_col = 0
do ic0=1,geom%nc0
   if (any(geom%mask(ic0,:))) mask_ind_col(ic0) = 1
end do

! Initialize support radius to 1.0
rh0 = 1.0

! Compute subset
write(mpl%unit,'(a7,a)') '','Compute horizontal subset C1'
if (nam%nc1<maxval(count(geom%mask,dim=1))) then
   call initialize_sampling(geom%nc0,dble(geom%lon),dble(geom%lat),mask_ind_col,rh0,nam%ntry,nam%nrep, &
 & nam%nc1,hdata%ic1_to_ic0)
else
   ic1 = 0
   do ic0=1,geom%nc0
      if (any(geom%mask(ic0,:))) then
         ic1 = ic1+1
         hdata%ic1_to_ic0(ic1) = ic0
      end if
   end do
end if

! Copy for all levels
do il0=1,geom%nl0
   hdata%ic1il0_log(:,il0) = geom%mask(hdata%ic1_to_ic0,il0)
   hdata%ic1icil0_to_ic0(:,1,il0) = hdata%ic1_to_ic0
   hdata%ic1icil0_log(:,1,il0) = geom%mask(hdata%ic1_to_ic0,il0)
end do

! End associate
end associate

end subroutine compute_sampling_zs

!----------------------------------------------------------------------
! Subroutine: compute_sampling_ps
!> Purpose: compute positive separation sampling
!----------------------------------------------------------------------
subroutine compute_sampling_ps(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: il0,jl0,irmaxloc,progint,ic,ic1,ir,ipt,jpt,i,nvc0,ic0,ivc0,icinf,icsup,ictest,ibnd
integer,allocatable :: vipt(:)
real(kind_real) :: d
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical :: found,valid,done(hdata%nam%nc*hdata%nam%nc1)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

if (nam%nc>1) then
   write(mpl%unit,'(a7,a)') '','Compute positive separation sampling'

   do il0=1,geom%nl0
      if ((il0>1).and.(geom%nl0i==1).and.all(hdata%ic1icil0_to_ic0(:,1,il0)==hdata%ic1icil0_to_ic0(:,1,1))) then
         ! Copy sampling
         hdata%ic1icil0_log(:,:,il0) = hdata%ic1icil0_log(:,:,1)
         hdata%ic1icil0_to_ic0(:,:,il0) = hdata%ic1icil0_to_ic0(:,:,1)
      else
         ! New sampling
         write(mpl%unit,'(a10,a,i2,a,i3,a)',advance='no') '','Level ',nam%levs(il0),' ~> new sampling:'

         ! Initialize
         do ic=1,nam%nc
            if (ic/=1) then
               hdata%ic1icil0_log(:,ic,il0) = .false.
               call msi(hdata%ic1icil0_to_ic0(:,ic,il0))
            end if
         end do

         ! Define valid nodes vector
         nvc0 = count(geom%mask(:,il0))
         allocate(vipt(nvc0))
         ivc0 = 0
         do ic0=1,geom%nc0
            if (geom%mask(ic0,il0)) then
               ivc0 = ivc0+1
               vipt(ivc0) = ic0
            end if
         end do

         ! Sample classes of positive separation
         call prog_init(progint)
         ir = 0
         irmaxloc = irmax
         do while (.not.all(hdata%ic1icil0_log(:,:,il0)).and.(nvc0>1).and.(ir<=irmaxloc))
            ! Try a random point
            if (mpl%main) call rand_integer(1,nvc0,i)
            call mpl_bcast(i,mpl%ioproc)
            ir = ir+1
            jpt = vipt(i)

            !$omp parallel do schedule(static) private(ic1,ipt,d,icinf,icsup,found,ictest,ic,valid,x,y,z,v1,v2,va,ibnd,vp,t)
            do ic1=1,nam%nc1
               ! Allocation
               allocate(x(2))
               allocate(y(2))
               allocate(z(2))
               allocate(v1(3))
               allocate(v2(3))
               allocate(va(3))
               allocate(vp(3))
               allocate(t(4))

               ! Check if there is a valid first point
               if (hdata%ic1il0_log(ic1,il0)) then
                  ! Compute the distance
                  ipt = hdata%ic1icil0_to_ic0(ic1,1,il0)
                  call sphere_dist(geom%lon(ipt),geom%lat(ipt),geom%lon(jpt),geom%lat(jpt),d)

                  ! Find the class (dichotomy method)
                  if ((d>0.0).and.(d<(float(nam%nc)-0.5)*nam%dc)) then
                     ic = 1
                     icinf = 1
                     icsup = nam%nc
                     found = .false.
                     do while (.not.found)
                        ! New value
                        ictest = (icsup+icinf)/2

                        ! Update
                        if (d<(float(ictest)-0.5)*nam%dc) icsup = ictest
                        if (d>(float(ictest)-0.5)*nam%dc) icinf = ictest

                        ! Exit test
                        if (icsup==icinf+1) then
                           if (abs((float(icinf)-0.5)*nam%dc-d)<abs((float(icsup)-0.5)*nam%dc-d)) then
                              ic = icinf
                           else
                              ic = icsup
                           end if
                           found = .true.
                        end if
                     end do

                     ! Find if this class has not been aready filled
                     if ((ic/=1).and.(.not.hdata%ic1icil0_log(ic1,ic,il0))) then
                        valid = .true.
                        if (nam%mask_check) then
                           ! Transform to cartesian coordinates
                           call trans(2,geom%lat((/ipt,jpt/)),geom%lon((/ipt,jpt/)),x,y,z)

                           ! Compute arc orthogonal vector
                           v1 = (/x(1),y(1),z(1)/)
                           v2 = (/x(2),y(2),z(2)/)
                           call vector_product(v1,v2,va)

                           ! Check if arc is crossing boundary arcs
                           do ibnd=1,geom%nbnd(il0)
                              call vector_product(va,geom%vbnd(:,ibnd,il0),vp)
                              v1 = (/x(1),y(1),z(1)/)
                              call vector_triple_product(v1,va,vp,t(1))
                              v1 = (/x(2),y(2),z(2)/)
                              call vector_triple_product(v1,va,vp,t(2))
                              v1 = (/geom%xbnd(1,ibnd,il0),geom%ybnd(1,ibnd,il0),geom%zbnd(1,ibnd,il0)/)
                              call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(3))
                              v1 = (/geom%xbnd(2,ibnd,il0),geom%ybnd(2,ibnd,il0),geom%zbnd(2,ibnd,il0)/)
                              call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(4))
                              t(1) = -t(1)
                              t(3) = -t(3)
                              if (all(t>0).or.(all(t<0))) then
                                 valid = .false.
                                 exit
                              end if
                           end do
                        end if

                        if (valid) then
                           ! Fill hdata
                           hdata%ic1icil0_log(ic1,ic,il0) = .true.
                           hdata%ic1icil0_to_ic0(ic1,ic,il0) = jpt
                        end if
                     end if
                  end if
               end if

               ! Release memory
               deallocate(x)
               deallocate(y)
               deallocate(z)
               deallocate(v1)
               deallocate(v2)
               deallocate(va)
               deallocate(vp)
               deallocate(t)
            end do
            !$omp end parallel do

            ! Update valid nodes vector
            vipt(i) = vipt(nvc0)
            nvc0 = nvc0-1

            ! Print progression
            done = pack(hdata%ic1icil0_log(:,:,il0),mask=.true.)
            call prog_print(progint,done)
         end do
         write(mpl%unit,'(a)') '100%'

         ! Release memory
         deallocate(vipt)
      end if
   end do
elseif (nam%nc==1) then
   ! Special case with one class
   do ic1=1,nam%nc1
      ! Check location validity
      if (isnotmsi(hdata%ic1_to_ic0(ic1))) then
         hdata%ic1icil0_log(ic1,1,1) = .true.
         hdata%ic1icil0_to_ic0(ic1,1,1) = ic1
      end if
   end do
end if

! Define sampling weights
write(mpl%unit,'(a7,a)') '','Define sampling weights'
do jl0=1,geom%nl0
   do il0=1,geom%nl0
      do ic=1,nam%nc
         do ic1=1,nam%nc1
            ! Absolute weight of each couple
            if (hdata%ic1il0_log(ic1,jl0).and.hdata%ic1icil0_log(ic1,ic,il0)) then
               hdata%swgt(ic1,ic,il0,jl0) = 1.0
            else
               hdata%swgt(ic1,ic,il0,jl0) = 0.0
             end if
         end do

         ! Relative weight
         if (sum(hdata%swgt(:,ic,il0,jl0))>0) hdata%swgt(:,ic,il0,jl0) = hdata%swgt(:,ic,il0,jl0)/sum(hdata%swgt(:,ic,il0,jl0))
      end do
   end do
end do

! End associate
end associate

end subroutine compute_sampling_ps

!----------------------------------------------------------------------
! Subroutine: compute_sampling_lct
!> Purpose: compute LCT sampling
!----------------------------------------------------------------------
subroutine compute_sampling_lct(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: il0,ic1,ic0,jc0,ibnd,ic
integer :: nn(hdata%nam%nc)
real(kind_real) :: dum(hdata%nam%nc)
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

write(mpl%unit,'(a7,a)') '','Compute LCT sampling'

do il0=1,geom%nl0
   do ic1=1,nam%nc1
      ! Check location validity
      if (isnotmsi(hdata%ic1_to_ic0(ic1))) then
         ! Find neighbors
         call find_nearest_neighbors(geom%ctree(min(il0,geom%nl0i)),dble(geom%lon(hdata%ic1_to_ic0(ic1))), &
       & dble(geom%lat(hdata%ic1_to_ic0(ic1))),nam%nc,nn,dum)
         hdata%ic1icil0_to_ic0(ic1,:,il0) = nn
         hdata%ic1icil0_log(ic1,:,il0) = .true.

         if (nam%mask_check) then
            ! Check that great circle to neighbors is not crossing mask boundaries
            !$omp parallel do schedule(static) private(ic,x,y,z,v1,v2,va,vp,t,ic0,jc0)
            do ic=1,nam%nc
               ! Allocation
               allocate(x(2))
               allocate(y(2))
               allocate(z(2))
               allocate(v1(3))
               allocate(v2(3))
               allocate(va(3))
               allocate(vp(3))
               allocate(t(4))

               ! Indices
               ic0 = hdata%ic1_to_ic0(ic1)
               jc0 = hdata%ic1icil0_to_ic0(ic1,ic,il0)

               ! Transform to cartesian coordinates
               call trans(2,geom%lat((/ic0,jc0/)),geom%lon((/ic0,jc0/)),x,y,z)

               ! Compute arc orthogonal vector
               v1 = (/x(1),y(1),z(1)/)
               v2 = (/x(2),y(2),z(2)/)
               call vector_product(v1,v2,va)

               ! Check if arc is crossing boundary arcs
               do ibnd=1,geom%nbnd(il0)
                  call vector_product(va,geom%vbnd(:,ibnd,il0),vp)
                  v1 = (/x(1),y(1),z(1)/)
                  call vector_triple_product(v1,va,vp,t(1))
                  v1 = (/x(2),y(2),z(2)/)
                  call vector_triple_product(v1,va,vp,t(2))
                  v1 = (/geom%xbnd(1,ibnd,il0),geom%ybnd(1,ibnd,il0),geom%zbnd(1,ibnd,il0)/)
                  call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(3))
                  v1 = (/geom%xbnd(2,ibnd,il0),geom%ybnd(2,ibnd,il0),geom%zbnd(2,ibnd,il0)/)
                  call vector_triple_product(v1,geom%vbnd(:,ibnd,il0),vp,t(4))
                  t(1) = -t(1)
                  t(3) = -t(3)
                  if (all(t>0).or.(all(t<0))) then
                     call msi(hdata%ic1icil0_to_ic0(ic1,ic,il0))
                     hdata%ic1icil0_log(ic1,ic,il0) = .false.
                     exit
                  end if
               end do

               ! Memory release
               deallocate(x)
               deallocate(y)
               deallocate(z)
               deallocate(v1)
               deallocate(v2)
               deallocate(va)
               deallocate(vp)
               deallocate(t)
            end do
            !$omp end parallel do
         end if
      end if
   end do
end do

! End associate
end associate

end subroutine compute_sampling_lct

!----------------------------------------------------------------------
! Subroutine: compute_sampling_mesh
!> Purpose: compute sampling mesh
!----------------------------------------------------------------------
subroutine compute_sampling_mesh(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: info,ic2
integer :: list(6*(hdata%nc2-2)),lptr(6*(hdata%nc2-2)),lend(hdata%nc2),lnew,near(hdata%nc2)
integer :: next(hdata%nc2),ltri(9,2*(hdata%nc2-2))
integer :: na,ia,it,i,i1,i2,nodes(hdata%nc2),nb,natmp,nttmp,nab,larcb(2,3*(hdata%nc2-2)),iab
real(kind_real) :: x(hdata%nc2),y(hdata%nc2),z(hdata%nc2),dist(hdata%nc2),dist_12
real(kind_real) :: v1(3),v2(3),vp(3),v(3),vf(3),vt(3),tlat,tlon,trad,dist_t1,dist_t2

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Transform to cartesian coordinates
call trans(hdata%nc2,geom%lat(hdata%ic2_to_ic0),geom%lon(hdata%ic2_to_ic0),x,y,z)

! Create mesh
list = 0
call trmesh(hdata%nc2,x,y,z,list,lptr,lend,lnew,near,next,dist,info)

! Create triangles list
call trlist(hdata%nc2,list,lptr,lend,9,hdata%nt,ltri,info)

! Copy triangle list
allocate(hdata%ltri(3,hdata%nt))
hdata%ltri(:,1:hdata%nt) = ltri(1:3,1:hdata%nt)

! Copy arcs list
na = maxval(ltri(7:9,1:hdata%nt))
allocate(hdata%larc(2,na))
do ia=1,na
   it = 1
   do while (it<=hdata%nt)
      if (any(ltri(7:9,it)==ia)) exit
      it = it+1
   end do
   i = 1
   do while (i<=3)
      if (ltri(6+i,it)==ia) exit
      i = i+1
   end do
   i1 = mod(i+1,3)
   if (i1==0) i1 = 3
   i2 = mod(i+2,3)
   if (i2==0) i2 = 3
   hdata%larc(1,ia) = ltri(i1,it)
   hdata%larc(2,ia) = ltri(i2,it)
end do

! Find boundary nodes
allocate(hdata%bdist(hdata%nc2))
call bnodes(hdata%nc2,list,lptr,lend,nodes,nb,natmp,nttmp)
if (nb>0) then
   ! Find boundary arcs
   nab = 0
   do ia=1,na
      if (any(nodes==hdata%larc(1,ia)).and.any(nodes==hdata%larc(2,ia))) then
         nab = nab+1
         larcb(:,nab) = hdata%larc(:,ia)
      end if
   end do

   ! Find minimal distance to a boundary arc
   hdata%bdist = huge(1.0)
   do iab=1,nab
      ! Distance
      call sphere_dist(geom%lon(hdata%ic2_to_ic0(larcb(1,iab))),geom%lat(hdata%ic2_to_ic0(larcb(1,iab))), &
    & geom%lon(hdata%ic2_to_ic0(larcb(2,iab))),geom%lat(hdata%ic2_to_ic0(larcb(2,iab))),dist_12)

      ! Vectors
      v1 = (/x(larcb(1,iab)),y(larcb(1,iab)),z(larcb(1,iab))/)
      v2 = (/x(larcb(2,iab)),y(larcb(2,iab)),z(larcb(2,iab))/)

      ! Compute normal vector to the boundary arc plane
      call vector_product(v1,v2,vp)

      ! Compute the shortest distance from each point to the boundary arc great-circle
      do ic2=1,hdata%nc2
         ! Vector
         v = (/x(ic2),y(ic2),z(ic2)/)

         ! Vector products
         call vector_product(v,vp,vf)
         call vector_product(vp,vf,vt)

         ! Back to spherical coordinates
         call scoord(vt(1),vt(2),vt(3),tlat,tlon,trad)

         ! Check whether T is on the arc
         call sphere_dist(tlon,tlat,geom%lon(hdata%ic2_to_ic0(larcb(1,iab))),geom%lat(hdata%ic2_to_ic0(larcb(1,iab))),dist_t1)
         call sphere_dist(tlon,tlat,geom%lon(hdata%ic2_to_ic0(larcb(2,iab))),geom%lat(hdata%ic2_to_ic0(larcb(2,iab))),dist_t2)
         if ((dist_t1<dist_12).and.(dist_t2<dist_12)) then
            ! T is on the arc
            call sphere_dist(geom%lon(hdata%ic2_to_ic0(ic2)),geom%lat(hdata%ic2_to_ic0(ic2)),tlon,tlat,dist_t1)
            hdata%bdist(ic2) = min(hdata%bdist(ic2),dist_t1)
         else
            ! T is not on the arc
            call sphere_dist(geom%lon(hdata%ic2_to_ic0(ic2)),geom%lat(hdata%ic2_to_ic0(ic2)), &
          & geom%lon(hdata%ic2_to_ic0(larcb(1,iab))),geom%lat(hdata%ic2_to_ic0(larcb(1,iab))),dist_t1)
            call sphere_dist(geom%lon(hdata%ic2_to_ic0(ic2)),geom%lat(hdata%ic2_to_ic0(ic2)), &
          & geom%lon(hdata%ic2_to_ic0(larcb(2,iab))),geom%lat(hdata%ic2_to_ic0(larcb(2,iab))),dist_t2)
            hdata%bdist(ic2) = min(hdata%bdist(ic2),min(dist_t1,dist_t2))
         end if
      end do
   end do
   hdata%bdist = hdata%bdist/req
else
   hdata%bdist = huge(1.0)
end if

! End associate
end associate

end subroutine compute_sampling_mesh

end module module_sampling
