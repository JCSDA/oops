!----------------------------------------------------------------------
! Module: module_normalization.f90
!> Purpose: normalization routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_normalization

use module_namelist, only: nam
use omp_lib
use tools_display, only: msgerror,ddis,prog_init,prog_print
use type_fields, only: fldtype,alphatype
use tools_kinds,only: kind_real
use tools_missing, only: msr,isnotmsi,msi
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send,mpl_barrier
use type_sdata, only: sdatatype
implicit none

private
public :: compute_normalization

contains

!----------------------------------------------------------------------
! Subroutine: compute_normalization
!> Purpose: compute normalization
!----------------------------------------------------------------------
subroutine compute_normalization(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0i,i_s,ic1,jc2,is,js,ic0,il0,il1,ih,iv,nlr,ilr,jlr,ic,progint,is_add
integer :: iproc,ic0_s(mpl%nproc),ic0_e(mpl%nproc),nc0_loc(mpl%nproc),ic0_loc,ibuf
integer,allocatable :: ineh(:,:),inev(:,:),ines(:,:),inec(:),order(:),is_list(:)
integer,allocatable :: interp_h(:,:,:),interp_v(:,:,:),interp_s(:,:,:),convol_c(:,:)
real(kind_real) :: S_add
real(kind_real),allocatable :: interp_h_S(:,:,:),interp_v_S(:,:,:),interp_s_S(:,:,:),convol_c_S(:,:)
real(kind_real),allocatable :: S_list(:),S_list_tmp(:),norm(:),normg(:,:)
logical :: conv
logical,allocatable :: done(:),valid_list_tmp(:)

! Compute horizontal interpolation inverse mapping
allocate(ineh(sdata%nc0,sdata%nl0i))
ineh = 0
do il0i=1,sdata%nl0i
   do i_s=1,sdata%h(il0i)%n_s
      ic0 = sdata%h(il0i)%row(i_s)
      ineh(ic0,il0i) = ineh(ic0,il0i)+1
   end do
end do
allocate(interp_h(maxval(ineh),sdata%nc0,sdata%nl0i))
allocate(interp_h_S(maxval(ineh),sdata%nc0,sdata%nl0i))
ineh = 0
do il0i=1,sdata%nl0i
   do i_s=1,sdata%h(il0i)%n_s
      ic0 = sdata%h(il0i)%row(i_s)
      ineh(ic0,il0i) = ineh(ic0,il0i)+1
      interp_h(ineh(ic0,il0i),ic0,il0i) = sdata%h(il0i)%col(i_s)
      interp_h_S(ineh(ic0,il0i),ic0,il0i) = sdata%h(il0i)%S(i_s)
   end do
end do

! Compute vertical interpolation inverse mapping
allocate(inev(sdata%nl0,sdata%nl0i))
inev = 0
do il0i=1,sdata%nl0i
   do i_s=1,sdata%v(il0i)%n_s
      il0 = sdata%v(il0i)%row(i_s)
      inev(il0,il0i) = inev(il0,il0i)+1
   end do
end do
allocate(interp_v(maxval(inev),sdata%nl0,sdata%nl0i))
allocate(interp_v_S(maxval(inev),sdata%nl0,sdata%nl0i))
inev = 0
do il0i=1,sdata%nl0i
   do i_s=1,sdata%v(il0i)%n_s
      il0 = sdata%v(il0i)%row(i_s)
      inev(il0,il0i) = inev(il0,il0i)+1
      interp_v(inev(il0,il0i),il0,il0i) = sdata%v(il0i)%col(i_s)
      interp_v_S(inev(il0,il0i),il0,il0i) = sdata%v(il0i)%S(i_s)
   end do
end do

! Compute subsampling interpolation inverse mapping
allocate(ines(sdata%nc1,sdata%nl1))
ines = 0
do il1=1,sdata%nl1
   do i_s=1,sdata%s(il1)%n_s
      ic1 = sdata%s(il1)%row(i_s)
      ines(ic1,il1) = ines(ic1,il1)+1
   end do
end do
allocate(interp_s(maxval(ines),sdata%nc1,sdata%nl1))
allocate(interp_s_S(maxval(ines),sdata%nc1,sdata%nl1))
ines = 0
do il1=1,sdata%nl1
   do i_s=1,sdata%s(il1)%n_s
      ic1 = sdata%s(il1)%row(i_s)
      ines(ic1,il1) = ines(ic1,il1)+1
      jc2 = sdata%s(il1)%col(i_s)
      interp_s(ines(ic1,il1),ic1,il1) = sdata%ic2il1_to_is(jc2,il1)
      interp_s_S(ines(ic1,il1),ic1,il1) = sdata%s(il1)%S(i_s)
   end do
end do

! Compute convolution inverse mapping
allocate(inec(sdata%ns))
inec = 0
do i_s=1,sdata%c%n_s
   is = sdata%c%col(i_s)
   js = sdata%c%row(i_s)

   inec(is) = inec(is)+1
   inec(js) = inec(js)+1
end do
allocate(convol_c(maxval(inec),sdata%ns))
allocate(convol_c_S(maxval(inec),sdata%ns))
call msi(convol_c)
call msr(convol_c_S)
inec = 0
do i_s=1,sdata%c%n_s
   is = sdata%c%col(i_s)
   js = sdata%c%row(i_s)

   inec(is) = inec(is)+1
   convol_c(inec(is),is) = sdata%c%row(i_s)
   convol_c_S(inec(is),is) = sdata%c%S(i_s)

   inec(js) = inec(js)+1
   convol_c(inec(js),js) = sdata%c%col(i_s)
   convol_c_S(inec(js),js) = sdata%c%S(i_s)
end do

! Re-order indices
do is=1,sdata%ns
   allocate(order(inec(is)))
   call qsort(inec(is),convol_c(1:inec(is),is),order)
   convol_c_S(1:inec(is),is) = convol_c_S(order(1:inec(is)),is)
   deallocate(order)
end do

if (nam%lsqrt) then
   ! Compute internal normalization weights
   write(mpl%unit,'(a7,a)') '','Compute internal normalization weights'

   ! Allocation
   allocate(sdata%norm_sqrt%val(sdata%ns))

   do is=1,sdata%ns
      ! Allocation
      allocate(S_list_tmp(sdata%ns))
      allocate(valid_list_tmp(sdata%ns))

      ! Convolution
      S_list_tmp = 0.0
      valid_list_tmp = .false.
      do ic=1,inec(is)
         js = convol_c(ic,is)
         S_list_tmp(js) = S_list_tmp(js)+convol_c_S(ic,is)
         valid_list_tmp(js) = .true.
      end do

      ! Sum of squared values
      sdata%norm_sqrt%val(is) = 1.0
      do js=1,sdata%ns
         if (valid_list_tmp(js)) sdata%norm_sqrt%val(is) = sdata%norm_sqrt%val(is)+S_list_tmp(js)**2
      end do

      ! Normalization factor
      sdata%norm_sqrt%val(is) = 1.0/sqrt(sdata%norm_sqrt%val(is))

      ! Release memory
      deallocate(S_list_tmp)
      deallocate(valid_list_tmp)
   end do
end if

! MPI splitting
do iproc=1,mpl%nproc
   ic0_s(iproc) = (iproc-1)*(sdata%nc0/mpl%nproc+1)+1
   ic0_e(iproc) = min(iproc*(sdata%nc0/mpl%nproc+1),sdata%nc0)
   nc0_loc(iproc) = ic0_e(iproc)-ic0_s(iproc)+1
end do

! Allocation
allocate(done(nc0_loc(mpl%myproc)*sdata%nl0))
allocate(norm(nc0_loc(mpl%myproc)*sdata%nl0))
call msr(norm)

! Compute normalization weights
write(mpl%unit,'(a7,a)',advance='no') '','Compute normalization weights: '
call prog_init(progint,done)
do il0=1,sdata%nl0
   il0i = min(il0,sdata%nl0i)

   !$omp parallel do private(ic0_loc,ic0,ibuf,is_list,order,S_list,S_list_tmp,valid_list_tmp,nlr) &
   !$omp&            private(is_add,S_add,ih,ic1,iv,il1,is,ilr,conv,ic,jlr,js)
   do ic0_loc=1,nc0_loc(mpl%myproc)
      ! MPI offset
      ic0 = ic0_s(mpl%myproc)+ic0_loc-1
      ibuf = (il0-1)*nc0_loc(mpl%myproc)+ic0_loc

      if (sdata%mask(ic0,il0)) then
         ! Allocation
         allocate(is_list(ineh(ic0,il0i)*maxval(inev(il0,:))*maxval(ines)))
         allocate(order(ineh(ic0,il0i)*maxval(inev(il0,:))*maxval(ines)))
         allocate(S_list(ineh(ic0,il0i)*maxval(inev(il0,:))*maxval(ines)))
         if (nam%lsqrt) then
            allocate(S_list_tmp(sdata%ns))
            allocate(valid_list_tmp(sdata%ns))
         else
            allocate(S_list_tmp(ineh(ic0,il0i)*inev(il0,il0i)*maxval(ines)))
         end if

         ! Initialization
         is_list = 0
         S_list = 0.0

         ! Adjoint interpolation
         nlr = 0
         do ih=1,ineh(ic0,il0i)
            ic1 = interp_h(ih,ic0,il0i)
            do iv=1,inev(il0,il0i)
               il1 = interp_v(iv,il0,il0i)
               do is=1,ines(ic1,il1)
                  is_add = interp_s(is,ic1,il1)
                  S_add = interp_h_S(ih,ic0,il0i)*interp_v_S(iv,il0,il0i)*interp_s_S(is,ic1,il1)
                  if (nlr==0) then
                     ilr = 1
                     nlr = 1
                  else
                     do ilr=1,nlr
                        if (is_add==is_list(ilr)) exit
                     end do
                     if (ilr==nlr+1) nlr = nlr+1
                  end if
                  is_list(ilr) = is_add
                  S_list(ilr) = S_list(ilr)+S_add
               end do
            end do
         end do

         if (nam%lsqrt) then
            ! Internal normalization
            do ilr=1,nlr
               is = is_list(ilr)
               S_list(ilr) = S_list(ilr)*sdata%norm_sqrt%val(is)
            end do

            ! Initialization
            S_list_tmp = 0.0
            valid_list_tmp = .false.
            do ilr=1,nlr
               is = is_list(ilr)
               S_list_tmp(is) = S_list(ilr)
               valid_list_tmp(is) = .true.
            end do

            ! Convolution
            do ilr=1,nlr
               is = is_list(ilr)
               do ic=1,inec(is)
                  js = convol_c(ic,is)
                  S_list_tmp(js) = S_list_tmp(js)+convol_c_S(ic,is)*S_list(ilr)
                  valid_list_tmp(js) = .true.
               end do
            end do

            ! Sum of squared values
            norm(ibuf) = 0.0
            do js=1,sdata%ns
               if (valid_list_tmp(js)) then
                  norm(ibuf) = norm(ibuf)+S_list_tmp(js)**2
               end if
            end do

            ! Normalization factor
            norm(ibuf) = 1.0/sqrt(norm(ibuf))
         else
            ! Sort arrays
            call qsort(nlr,is_list(1:nlr),order(1:nlr))
            S_list(1:nlr) = S_list(order(1:nlr))

            ! Convolution
            S_list_tmp(1:nlr) = S_list(1:nlr)
            do ilr=1,nlr
               is = is_list(ilr)
               conv = .true.
               ic = 1
               do jlr=ilr+1,nlr
                  js = is_list(jlr)
                  do while (convol_c(ic,is)/=js)
                     ic = ic+1
                     if (ic>inec(is)) then
                        conv = .false.
                        exit
                     end if
                  end do
                  if (conv) then
                     S_list(ilr) = S_list(ilr)+convol_c_S(ic,is)*S_list_tmp(jlr)
                     S_list(jlr) = S_list(jlr)+convol_c_S(ic,is)*S_list_tmp(ilr)
                  end if
               end do
            end do

            ! Normalization factor
            if (nlr>0) norm(ibuf) = 1.0/sqrt(sum(S_list(1:nlr)*S_list_tmp(1:nlr)))
         end if

         ! Print progression
         done(ibuf) = .true.
         call prog_print(progint,done)

         ! Release memory
         deallocate(is_list)
         deallocate(order)
         deallocate(S_list)
         deallocate(S_list_tmp)
         if (nam%lsqrt) deallocate(valid_list_tmp)
      end if
   end do
   !$omp end parallel do
end do
write(mpl%unit,'(a)') '100%'

! Allocation
allocate(sdata%norm%val(sdata%nc0,sdata%nl0))

! Communication
if (mpl%main) then
   ! Allocation
   allocate(normg(maxval(nc0_loc)*sdata%nl0,mpl%nproc))

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         normg(1:nc0_loc(iproc)*sdata%nl0,iproc) = norm
      else
         ! Receive data on ioproc
         call mpl_recv(nc0_loc(iproc)*sdata%nl0,normg(1:nc0_loc(iproc)*sdata%nl0,iproc),iproc,mpl%tag)
      end if
   end do

   ! Format data
   do iproc=1,mpl%nproc
      do il0=1,sdata%nl0
         do ic0_loc=1,nc0_loc(iproc)
            ic0 = ic0_s(iproc)+ic0_loc-1
            ibuf = (il0-1)*nc0_loc(iproc)+ic0_loc
            sdata%norm%val(ic0,il0) = normg(ibuf,iproc)
         end do
      end do
   end do

   ! Release memory
   deallocate(normg)
else
   ! Send data to ioproc
   call mpl_send(nc0_loc(mpl%myproc)*sdata%nl0,norm,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Broadcast data
call mpl_bcast(sdata%norm%val,mpl%ioproc)

! Release memory
deallocate(norm)
deallocate(done)

end subroutine compute_normalization

end module module_normalization
