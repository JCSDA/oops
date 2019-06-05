!----------------------------------------------------------------------
! Module: type_rng
! Purpose: random numbers generator derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_rng

use iso_fortran_env, only: int64
use tools_const, only: pi
use tools_func, only: sphere_dist,gc99
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use tools_repro, only: inf,sup,infeq
use type_tree, only: tree_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

integer,parameter :: default_seed = 140587            ! Default seed
integer(kind=int64),parameter :: a = 1103515245_int64 ! Linear congruential multiplier
integer(kind=int64),parameter :: c = 12345_int64      ! Linear congruential offset
integer(kind=int64),parameter :: m = 2147483648_int64 ! Linear congruential modulo
logical,parameter :: nn_stats = .false.               ! Compute and print subsampling statistics

type rng_type
   integer(kind=int64) :: seed
contains
   procedure :: init => rng_init
   procedure :: reseed => rng_reseed
   procedure :: lcg => rng_lcg
   procedure :: rng_rand_integer_0d
   procedure :: rng_rand_integer_1d
   generic :: rand_integer => rng_rand_integer_0d,rng_rand_integer_1d
   procedure :: rng_rand_real_0d
   procedure :: rng_rand_real_1d
   procedure :: rng_rand_real_2d
   procedure :: rng_rand_real_3d
   procedure :: rng_rand_real_4d
   procedure :: rng_rand_real_5d
   generic :: rand_real => rng_rand_real_0d,rng_rand_real_1d,rng_rand_real_2d,rng_rand_real_3d,rng_rand_real_4d,rng_rand_real_5d
   procedure :: rng_rand_gau_1d
   procedure :: rng_rand_gau_5d
   generic :: rand_gau => rng_rand_gau_1d,rng_rand_gau_5d
   procedure :: initialize_sampling => rng_initialize_sampling
end type rng_type

private
public :: rng_type

contains

!----------------------------------------------------------------------
! Subroutine: rng_init
! Purpose: initialize the random number generator
!----------------------------------------------------------------------
subroutine rng_init(rng,mpl,nam)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist variables

! Local variable
integer :: seed

! Set seed
if (nam%default_seed) then
   ! Default seed
   seed = default_seed
else
   ! Time-based seed
   call system_clock(count=seed)
end if

! Different seed for each task
seed = seed+mpl%myproc

! Long integer
rng%seed = int(seed,kind=int64)

! Print result
if (nam%default_seed) then
   write(mpl%info,'(a7,a)') '','Linear congruential generator initialized with a default seed'
   call mpl%flush
else
   write(mpl%info,'(a7,a)') '','Linear congruential generator initialized'
   call mpl%flush
end if

end subroutine rng_init

!----------------------------------------------------------------------
! Subroutine: rng_reseed
! Purpose: re-seed the random number generator
!----------------------------------------------------------------------
subroutine rng_reseed(rng,mpl)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng ! Random number generator
type(mpl_type),intent(inout) :: mpl  ! MPI data

! Local variable
integer :: seed

! Default seed
seed = default_seed

! Different seed for each task
seed = seed+mpl%myproc

! Long integer
rng%seed = int(seed,kind=int64)

end subroutine rng_reseed

!----------------------------------------------------------------------
! Subroutine: rng_lcg
! Purpose: linear congruential generator
!----------------------------------------------------------------------
subroutine rng_lcg(rng,x)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(out) :: x             ! Random number between 0 and 1

! Update seed
rng%seed = mod(a*rng%seed+c,m)

! Random number
x = real(rng%seed,kind_real)/real(m-1,kind_real)

end subroutine rng_lcg

!----------------------------------------------------------------------
! Subroutine: rng_rand_integer_0d
! Purpose: generate a random integer, 0d
!----------------------------------------------------------------------
subroutine rng_rand_integer_0d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
integer,intent(in) :: binf           ! Lower bound
integer,intent(in) :: bsup           ! Upper bound
integer,intent(out) :: ir            ! Random integer

! Local variables
real(kind_real) :: x

! Generate random number between 0 and 1
call rng%lcg(x)

! Adapt range
x = x*real(bsup-binf+1,kind_real)

! Add offset
ir = binf+int(x)

end subroutine rng_rand_integer_0d

!----------------------------------------------------------------------
! Subroutine: rng_rand_integer_1d
! Purpose: generate a random integer, 1d
!----------------------------------------------------------------------
subroutine rng_rand_integer_1d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
integer,intent(in) :: binf           ! Lower bound
integer,intent(in) :: bsup           ! Upper bound
integer,intent(out) :: ir(:)         ! Random integer

! Local variables
integer :: i

do i=1,size(ir)
   call rng%rand_integer(binf,bsup,ir(i))
end do

end subroutine rng_rand_integer_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_0d
! Purpose: generate a random real, 0d
!----------------------------------------------------------------------
subroutine rng_rand_real_0d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(in) :: binf   ! Lower bound
real(kind_real),intent(in) :: bsup   ! Upper bound
real(kind_real),intent(out) :: rr    ! Random integer

! Local variables
real(kind_real) :: x

! Generate random number between 0 and 1
call rng%lcg(x)

! Adapt range
x = x*(bsup-binf)

! Add offset
rr = binf+x

end subroutine rng_rand_real_0d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_1d
! Purpose: generate a random real, 1d
!----------------------------------------------------------------------
subroutine rng_rand_real_1d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(in) :: binf   ! Lower bound
real(kind_real),intent(in) :: bsup   ! Upper bound
real(kind_real),intent(out) :: rr(:) ! Random integer

! Local variables
integer :: i

do i=1,size(rr)
   call rng%rand_real(binf,bsup,rr(i))
end do

end subroutine rng_rand_real_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_2d
! Purpose: generate a random real, 2d
!----------------------------------------------------------------------
subroutine rng_rand_real_2d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng    ! Random number generator
real(kind_real),intent(in) :: binf      ! Lower bound
real(kind_real),intent(in) :: bsup      ! Upper bound
real(kind_real),intent(out) :: rr(:,:)  ! Random integer

! Local variables
integer :: i,j

do i=1,size(rr,2)
   do j=1,size(rr,1)
   call rng%rand_real(binf,bsup,rr(j,i))
   end do
end do

end subroutine rng_rand_real_2d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_3d
! Purpose: generate a random real, 3d
!----------------------------------------------------------------------
subroutine rng_rand_real_3d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng     ! Random number generator
real(kind_real),intent(in) :: binf       ! Lower bound
real(kind_real),intent(in) :: bsup       ! Upper bound
real(kind_real),intent(out) :: rr(:,:,:) ! Random integer

! Local variables
integer :: i,j,k

do i=1,size(rr,3)
   do j=1,size(rr,2)
      do k=1,size(rr,1)
         call rng%rand_real(binf,bsup,rr(k,j,i))
      end do
   end do
end do

end subroutine rng_rand_real_3d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_4d
! Purpose: generate a random real, 4d
!----------------------------------------------------------------------
subroutine rng_rand_real_4d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng       ! Random number generator
real(kind_real),intent(in) :: binf         ! Lower bound
real(kind_real),intent(in) :: bsup         ! Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:) ! Random integer

! Local variables
integer :: i,j,k,l

do i=1,size(rr,4)
   do j=1,size(rr,3)
      do k=1,size(rr,2)
         do l=1,size(rr,1)
            call rng%rand_real(binf,bsup,rr(l,k,j,i))
         end do
      end do
   end do
end do

end subroutine rng_rand_real_4d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_5d
! Purpose: generate a random real, 5d
!----------------------------------------------------------------------
subroutine rng_rand_real_5d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         ! Random number generator
real(kind_real),intent(in) :: binf           ! Lower bound
real(kind_real),intent(in) :: bsup           ! Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:,:) ! Random integer

! Local variables
integer :: i,j,k,l,m

do i=1,size(rr,5)
   do j=1,size(rr,4)
      do k=1,size(rr,3)
         do l=1,size(rr,3)
            do m=1,size(rr,1)
               call rng%rand_real(binf,bsup,rr(m,l,k,j,i))
            end do
         end do
      end do
   end do
end do

end subroutine rng_rand_real_5d

!----------------------------------------------------------------------
! Subroutine: rng_rand_gau_1d
! Purpose: generate random Gaussian deviates, 1d
!----------------------------------------------------------------------
subroutine rng_rand_gau_1d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(out) :: rr(:) ! Random integer

! Local variables
integer :: i,iset
real(kind_real) :: gasdev,fac,gset,rsq,v1,v2

! Normal distribution
iset = 0
do i=1,size(rr,1)
   if( iset==0) then
      rsq = 0.0
      do while((rsq>=1.0).or.(rsq<=0.0))
         call rng%rand_real(0.0_kind_real,1.0_kind_real,v1)
         v1 = 2.0*v1-1.0
         call rng%rand_real(0.0_kind_real,1.0_kind_real,v2)
         v2 = 2.0*v2-1.0
         rsq = v1**2+v2**2
      end do
      fac = sqrt(-2.0*log(rsq)/rsq)
      gset = v1*fac
      gasdev = v2*fac
      iset = 1
   else
      gasdev = gset
      iset = 0
   end if
   rr(i) = gasdev
end do

end subroutine rng_rand_gau_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_gau_5d
! Purpose: generate random Gaussian deviates, 5d
!----------------------------------------------------------------------
subroutine rng_rand_gau_5d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         ! Random number generator
real(kind_real),intent(out) :: rr(:,:,:,:,:) ! Random integer

! Local variables
integer :: i,j,k,l

do i=1,size(rr,5)
   do j=1,size(rr,4)
      do k=1,size(rr,3)
         do l=1,size(rr,2)
            call rng%rand_gau(rr(:,l,k,j,i))
         end do
      end do
   end do
end do

end subroutine rng_rand_gau_5d

!----------------------------------------------------------------------
! Subroutine: rng_initialize_sampling
! Purpose: intialize sampling
!----------------------------------------------------------------------
subroutine rng_initialize_sampling(rng,mpl,n,lon,lat,mask,nfor,for,rh,ntry,nrep,ns,ihor,fast)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
type(mpl_type),intent(inout) :: mpl  ! MPI data
integer,intent(in) :: n              ! Number of points
real(kind_real),intent(in) :: lon(n) ! Longitudes
real(kind_real),intent(in) :: lat(n) ! Latitudes
logical,intent(in) :: mask(n)        ! Mask
integer,intent(in) :: nfor           ! Number of forced points (included into the subsampling)
integer,intent(in) :: for(nfor)      ! Forced points
real(kind_real),intent(in) :: rh(n)  ! Horizontal support radius
integer,intent(in) :: ntry           ! Number of tries
integer,intent(in) :: nrep           ! Number of replacements
integer,intent(in) :: ns             ! Number of samplings points
integer,intent(out) :: ihor(ns)      ! Horizontal sampling index
logical,intent(in),optional :: fast  ! Fast sampling flag

! Local variables
integer :: system_clock_start,system_clock_end,count_rate,count_max
integer :: is,js,ifor,i,irep,irmax,itry,irval,irvalmin,irvalmax,i_red,ir,ismin,nval,nrep_eff,nn_index(2)
integer,allocatable :: val_to_full(:),ihor_tmp(:)
real(kind_real) :: elapsed
real(kind_real) :: d,distmax,distmin,nn_dist(2)
real(kind_real) :: nn_sdist_min,nn_sdist_max,nn_sdist_avg,nn_sdist_std
real(kind_real) :: cdf_norm,rr
real(kind_real),allocatable :: cdf(:)
real(kind_real),allocatable :: lon_rep(:),lat_rep(:),dist(:)
real(kind_real),allocatable :: sdist(:,:),nn_sdist(:)
logical :: lfast
logical,allocatable :: lmask(:),smask(:),rmask(:)
character(len=1024),parameter :: subr = 'rng_initialize_sampling'
type(tree_type) :: tree

if (mpl%main) then
   ! Check forced points
   do ifor=1,nfor
      if (.not.mask(for(ifor))) call mpl%abort(subr,'a forced point is out of the mask')
   end do
   if (ns<nfor) call mpl%abort(subr,'ns lower than the number of forced points')

   ! Check mask size
   nval = count(mask)
   if (nval==0) then
       call mpl%abort(subr,'empty mask in initialize sampling')
   elseif (nval<ns) then
      call mpl%abort(subr,'ns greater than mask size in initialize_sampling')
   elseif (nval==ns) then
      write(mpl%info,'(a)') ' all points are used'
      call mpl%flush

      ! Allocation
      allocate(lmask(n))

      ! Initialization
      lmask = mask
      is = 0

      ! Forced points
      do ifor=1,nfor
         is = is+1
         ir = for(ifor)
         ihor(ifor) = ir
         lmask(ir) = .false.
      end do

      ! Other points
      do i=1,n
         if (lmask(i)) then
            is = is+1
            ihor(is) = i
         end if
      end do
   else
      if (nn_stats) then
         ! Save initial time
         call system_clock(count=system_clock_start)
      end if

      ! Allocation
      nrep_eff = min(nrep,nval-ns)
      allocate(ihor_tmp(ns+nrep_eff))
      allocate(lmask(n))
      allocate(smask(n))
      allocate(val_to_full(nval))

      ! Initialization
      ihor_tmp = mpl%msv%vali
      lmask = mask
      smask = .false.
      val_to_full = mpl%msv%vali
      i_red = 0
      do i=1,n
         if (lmask(i)) then
            i_red = i_red+1
            val_to_full(i_red) = i
         end if
      end do
      call mpl%prog_init(ns+nrep_eff)
      lfast = .false.
      if (present(fast)) lfast = fast

      ! Forced points
      do ifor=1,nfor
         ir = for(ifor)
         irval = mpl%msv%vali
         do i=1,nval
            if (val_to_full(i)==ir) then
               irval = i
               exit
            end if
         end do
         if (mpl%msv%isi(irval)) call mpl%abort(subr,'cannot find irval in initialize_sampling')
         ihor_tmp(ifor) = ir
         lmask(ir) = .false.
         smask(ir) = .true.

         ! Shift valid points array
         if (irval<nval) val_to_full(irval:nval-1) = val_to_full(irval+1:nval)
         nval = nval-1
      end do

      if (lfast) then
         ! Allocation
         allocate(cdf(nval))

         ! Define sampling with cumulative distribution function
         cdf(1) = 0.0
         do i_red=2,nval
             i = val_to_full(i_red)
             cdf(i_red) = cdf(i_red-1)+1.0/rh(i)**2
         end do
         cdf_norm = 1.0/cdf(nval)
         cdf(1:nval) = cdf(1:nval)*cdf_norm

         do is=1+nfor,ns+nrep_eff
            ! Generate random number
            call rng%rand_real(0.0_kind_real,1.0_kind_real,rr)

            ! Dichotomy to find the value
            irvalmin = 1
            irvalmax = nval
            do while (irvalmax-irvalmin>1)
               irval = (irvalmin+irvalmax)/2
               if ((cdf(irvalmin)-rr)*(cdf(irval)-rr)>0.0) then
                  irvalmin = irval
               else
                  irvalmax = irval
               end if
            end do

            ! New sampling point
            ir = val_to_full(irval)
            ihor_tmp(is) = ir
            lmask(ir) = .false.

            ! Shift valid points array
            if (irval<nval) then
               cdf(irval:nval-1) = cdf(irval+1:nval)
               val_to_full(irval:nval-1) = val_to_full(irval+1:nval)
            end if
            nval = nval-1

            ! Renormalize cdf
            cdf_norm = 1.0/cdf(nval)
            cdf(1:nval) = cdf(1:nval)*cdf_norm

            ! Update
            call mpl%prog_print(is)
         end do
         call mpl%prog_final(.false.)

         ! Release memory
         deallocate(cdf)
      else
         ! Define sampling with KD-tree
         do is=1+nfor,ns+nrep_eff
            if (is>2) then
               ! Allocation
               call tree%alloc(mpl,n,mask=smask)

               ! Initialization
               call tree%init(lon,lat)
            end if

            ! Initialization
            distmax = 0.0
            irmax = 0
            irvalmax = 0
            itry = 1

            ! Find a new point
            do itry=1,ntry
               ! Generate a random index among valid points
               call rng%rand_integer(1,nval,irval)
               ir = val_to_full(irval)

               ! Check point validity
               if (is==1) then
                  ! Accept point
                  irvalmax = irval
                  irmax = ir
               else
                  if (is==2) then
                     ! Compute distance
                     call sphere_dist(lon(ir),lat(ir),lon(ihor_tmp(1)),lat(ihor_tmp(1)),d)
                  else
                     ! Find nearest neighbor distance
                     call tree%find_nearest_neighbors(lon(ir),lat(ir),1,nn_index(1:1),nn_dist(1:1))
                     d = nn_dist(1)**2/(rh(ir)**2+rh(nn_index(1))**2)
                  end if

                  ! Check distance
                  if (sup(d,distmax)) then
                     distmax = d
                     irvalmax = irval
                     irmax = ir
                  end if
               end if
            end do

            ! Delete tree
            if (is>2) call tree%dealloc

            ! Add point to sampling
            if (irmax>0) then
               ! New sampling point
               ihor_tmp(is) = irmax
               lmask(irmax) = .false.
               smask(irmax) = .true.

               ! Shift valid points array
               if (irvalmax<nval) val_to_full(irvalmax:nval-1) = val_to_full(irvalmax+1:nval)
               nval = nval-1
            end if

            ! Update
            call mpl%prog_print(is)
         end do
         call mpl%prog_final(.false.)

         ! Release memory
         deallocate(smask)
      end if

      if (nrep_eff>0) then
         ! Continue printing
         write(mpl%info,'(a)') ' => '
         call mpl%flush(.false.)

         ! Allocation
         allocate(rmask(ns+nrep_eff))
         allocate(lon_rep(ns+nrep_eff))
         allocate(lat_rep(ns+nrep_eff))
         allocate(dist(ns+nrep_eff))

         ! Initialization
         rmask = .true.
         do is=1,ns+nrep_eff
            lon_rep(is) = lon(ihor_tmp(is))
            lat_rep(is) = lat(ihor_tmp(is))
         end do
         dist = mpl%msv%valr
         call mpl%prog_init(nrep_eff)

         ! Remove closest points
         do irep=1,nrep_eff
            ! Allocation
            call tree%alloc(mpl,ns+nrep_eff,mask=rmask)

            ! Initialization
            call tree%init(lon_rep,lat_rep)

            ! Get minimum distance
            do is=1+nfor,ns+nrep_eff
               if (rmask(is)) then
                  ! Find nearest neighbor distance
                  call tree%find_nearest_neighbors(lon(ihor_tmp(is)),lat(ihor_tmp(is)),2,nn_index,nn_dist)
                  if (nn_index(1)==is) then
                     dist(is) = nn_dist(2)
                  elseif (nn_index(2)==is) then
                     dist(is) = nn_dist(1)
                  else
                     call mpl%abort(subr,'wrong index in replacement')
                  end if
                  dist(is) = dist(is)**2/(rh(ihor_tmp(nn_index(1)))**2+rh(ihor_tmp(nn_index(2)))**2)
               end if
            end do

            ! Delete tree
            call tree%dealloc

            ! Remove worst point
            distmin = huge(1.0)
            ismin = mpl%msv%vali
            do is=1+nfor,ns+nrep_eff
               if (rmask(is)) then
                  if (inf(dist(is),distmin)) then
                     ismin = is
                     distmin = dist(is)
                  end if
               end if
            end do
            rmask(ismin) = .false.

             ! Update
            call mpl%prog_print(irep)
         end do
         call mpl%prog_final

         ! Copy ihor
         js = 0
         do is=1,ns+nrep_eff
            if (rmask(is)) then
               js = js+1
               ihor(js) = ihor_tmp(is)
            end if
         end do

         ! Release memory
         deallocate(rmask)
         deallocate(lon_rep)
         deallocate(lat_rep)
         deallocate(dist)
      else
         ! Stop printing
         write(mpl%info,'(a)') ''
         call mpl%flush

         ! Copy ihor
         ihor = ihor_tmp
      end if

      if (nn_stats) then
         ! Save final time
         call system_clock(count=system_clock_end)
         call system_clock(count_rate=count_rate,count_max=count_max)
         if (system_clock_end<system_clock_start) then
            elapsed = real(system_clock_end-system_clock_start+count_max,kind_real) &
                    & /real(count_rate,kind_real)
         else
            elapsed = real(system_clock_end-system_clock_start,kind_real) &
                    & /real(count_rate,kind_real)
         end if

         ! Allocation
         allocate(sdist(ns,ns))
         allocate(nn_sdist(ns))

         ! Compute normalized distances between sampling points
         do is=1,ns
            sdist(is,is) = huge(1.0)
            do js=1,is-1
               call sphere_dist(lon(ihor(is)),lat(ihor(is)),lon(ihor(js)),lat(ihor(js)),sdist(is,js))
               sdist(is,js) = sdist(is,js)/sqrt(rh(ihor(is))**2+rh(ihor(js))**2)
               sdist(js,is) = sdist(is,js)
            end do
         end do

         ! Find nearest neighbor normalized distance
         do is=1,ns
            nn_sdist(is) = minval(sdist(:,is))
         end do

         ! Compute statistics
         nn_sdist_min = minval(nn_sdist)
         nn_sdist_max = maxval(nn_sdist)
         nn_sdist_avg = sum(nn_sdist)/real(ns,kind_real)
         nn_sdist_std = sqrt(sum((nn_sdist-nn_sdist_avg)**2)/real(ns-1,kind_real))

         ! Print statistics
         write(mpl%info,'(a10,a)') '','Nearest neighbor normalized distance statistics:'
         call mpl%flush
         write(mpl%info,'(a13,a,e9.2)') '','Minimum: ',nn_sdist_min
         call mpl%flush
         write(mpl%info,'(a13,a,e9.2)') '','Maximum: ',nn_sdist_max
         call mpl%flush
         write(mpl%info,'(a13,a,e9.2)') '','Average: ',nn_sdist_avg
         call mpl%flush
         write(mpl%info,'(a13,a,e9.2)') '','Std-dev: ',nn_sdist_std
         call mpl%flush
         write(mpl%info,'(a10,a,f8.3,a)') '','Elapsed time to compute the subsampling: ',elapsed,' s'
         call mpl%flush

         ! Release memory
         deallocate(sdist)
         deallocate(nn_sdist)
      end if

      ! Release memory
      deallocate(ihor_tmp)
      deallocate(lmask)
      deallocate(val_to_full)
   end if
else
   write(mpl%info,'(a)') ''
   call mpl%flush
end if

! Broadcast
call mpl%f_comm%broadcast(ihor,mpl%ioproc-1)

end subroutine rng_initialize_sampling

end module type_rng
