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

use module_namelist, only: namtype
use omp_lib
use tools_display, only: msgerror,ddis,prog_init,prog_print
use tools_kinds,only: kind_real
use tools_missing, only: msr,isnotmsi,msi
use tools_qsort, only: qsort
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send,mpl_barrier
use type_ndata, only: ndatatype

implicit none

private
public :: compute_normalization

contains

!----------------------------------------------------------------------
! Subroutine: compute_normalization
!> Purpose: compute normalization
!----------------------------------------------------------------------
subroutine compute_normalization(nam,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: il0i,i_s,ic1,jc2,is,js,ic0,il0,il1,ih,iv,nlr,ilr,jlr,ic,progint,is_add
integer :: iproc,ic0_s(mpl%nproc),ic0_e(mpl%nproc),nc0_loc(mpl%nproc),ic0_loc
integer,allocatable :: ineh(:,:),inev(:),ines(:,:),inec(:),order(:),is_list(:)
integer,allocatable :: interp_h(:,:,:),interp_v(:,:),interp_s(:,:,:),convol_c(:,:)
real(kind_real) :: S_add
real(kind_real),allocatable :: interp_h_S(:,:,:),interp_v_S(:,:),interp_s_S(:,:,:),convol_c_S(:,:)
real(kind_real),allocatable :: S_list(:),S_list_tmp(:)
logical :: conv
logical,allocatable :: done(:),valid_list_tmp(:)

! Compute horizontal interpolation inverse mapping
allocate(ineh(ndata%nc0,ndata%nl0i))
ineh = 0
do il0i=1,ndata%nl0i
   do i_s=1,ndata%h(il0i)%n_s
      ic0 = ndata%h(il0i)%row(i_s)
      ineh(ic0,il0i) = ineh(ic0,il0i)+1
   end do
end do
allocate(interp_h(maxval(ineh),ndata%nc0,ndata%nl0i))
allocate(interp_h_S(maxval(ineh),ndata%nc0,ndata%nl0i))
ineh = 0
do il0i=1,ndata%nl0i
   do i_s=1,ndata%h(il0i)%n_s
      ic0 = ndata%h(il0i)%row(i_s)
      ineh(ic0,il0i) = ineh(ic0,il0i)+1
      interp_h(ineh(ic0,il0i),ic0,il0i) = ndata%h(il0i)%col(i_s)
      interp_h_S(ineh(ic0,il0i),ic0,il0i) = ndata%h(il0i)%S(i_s)
   end do
end do

! Compute vertical interpolation inverse mapping
allocate(inev(ndata%nl0))
inev = 0
do i_s=1,ndata%v%n_s
   il0 = ndata%v%row(i_s)
   inev(il0) = inev(il0)+1
end do
allocate(interp_v(maxval(inev),ndata%nl0))
allocate(interp_v_S(maxval(inev),ndata%nl0))
inev = 0
do i_s=1,ndata%v%n_s
   il0 = ndata%v%row(i_s)
   inev(il0) = inev(il0)+1
   interp_v(inev(il0),il0) = ndata%v%col(i_s)
   interp_v_S(inev(il0),il0) = ndata%v%S(i_s)
end do

! Compute subsampling interpolation inverse mapping
allocate(ines(ndata%nc1,ndata%nl1))
ines = 0
do il1=1,ndata%nl1
   do i_s=1,ndata%s(il1)%n_s
      ic1 = ndata%s(il1)%row(i_s)
      ines(ic1,il1) = ines(ic1,il1)+1
   end do
end do
allocate(interp_s(maxval(ines),ndata%nc1,ndata%nl1))
allocate(interp_s_S(maxval(ines),ndata%nc1,ndata%nl1))
ines = 0
do il1=1,ndata%nl1
   do i_s=1,ndata%s(il1)%n_s
      ic1 = ndata%s(il1)%row(i_s)
      ines(ic1,il1) = ines(ic1,il1)+1
      jc2 = ndata%s(il1)%col(i_s)
      interp_s(ines(ic1,il1),ic1,il1) = ndata%ic2il1_to_is(jc2,il1)
      interp_s_S(ines(ic1,il1),ic1,il1) = ndata%s(il1)%S(i_s)
   end do
end do

! Compute convolution inverse mapping
allocate(inec(ndata%ns))
inec = 0
do i_s=1,ndata%c%n_s
   is = ndata%c%col(i_s)
   js = ndata%c%row(i_s)

   inec(is) = inec(is)+1
   inec(js) = inec(js)+1
end do
allocate(convol_c(maxval(inec),ndata%ns))
allocate(convol_c_S(maxval(inec),ndata%ns))
call msi(convol_c)
call msr(convol_c_S)
inec = 0
do i_s=1,ndata%c%n_s
   is = ndata%c%col(i_s)
   js = ndata%c%row(i_s)

   inec(is) = inec(is)+1
   convol_c(inec(is),is) = ndata%c%row(i_s)
   convol_c_S(inec(is),is) = ndata%c%S(i_s)

   inec(js) = inec(js)+1
   convol_c(inec(js),js) = ndata%c%col(i_s)
   convol_c_S(inec(js),js) = ndata%c%S(i_s)
end do

! Re-order indices
do is=1,ndata%ns
   allocate(order(inec(is)))
   call qsort(inec(is),convol_c(1:inec(is),is),order)
   convol_c_S(1:inec(is),is) = convol_c_S(order(1:inec(is)),is)
   deallocate(order)
end do

if (nam%lsqrt) then
   ! Compute internal normalization weights
   write(mpl%unit,'(a7,a)') '','Compute internal normalization weights'

   ! Allocation
   allocate(ndata%norm_sqrt(ndata%ns))

   do is=1,ndata%ns
      ! Allocation
      allocate(S_list_tmp(ndata%ns))
      allocate(valid_list_tmp(ndata%ns))

      ! Convolution
      S_list_tmp = 0.0
      valid_list_tmp = .false.
      do ic=1,inec(is)
         js = convol_c(ic,is)
         S_list_tmp(js) = S_list_tmp(js)+convol_c_S(ic,is)
         valid_list_tmp(js) = .true.
      end do

      ! Sum of squared values
      ndata%norm_sqrt(is) = 1.0
      do js=1,ndata%ns
         if (valid_list_tmp(js)) ndata%norm_sqrt(is) = ndata%norm_sqrt(is)+S_list_tmp(js)**2
      end do

      ! Normalization factor
      ndata%norm_sqrt(is) = 1.0/sqrt(ndata%norm_sqrt(is))

      ! Release memory
      deallocate(S_list_tmp)
      deallocate(valid_list_tmp)
   end do
end if

! MPI splitting
do iproc=1,mpl%nproc
   ic0_s(iproc) = (iproc-1)*(ndata%nc0/mpl%nproc+1)+1
   ic0_e(iproc) = min(iproc*(ndata%nc0/mpl%nproc+1),ndata%nc0)
   nc0_loc(iproc) = ic0_e(iproc)-ic0_s(iproc)+1
end do

! Allocation
allocate(ndata%norm(ndata%nc0,ndata%nl0))
allocate(done(nc0_loc(mpl%myproc)*ndata%nl0))
call msr(ndata%norm)

! Compute normalization weights
write(mpl%unit,'(a7,a)',advance='no') '','Compute normalization weights: '
call prog_init(progint,done)
do il0=1,ndata%nl0
   il0i = min(il0,ndata%nl0i)

   !$omp parallel do private(ic0_loc,ic0,is_list,order,S_list,S_list_tmp,valid_list_tmp,nlr) &
   !$omp&            private(is_add,S_add,ih,ic1,iv,il1,is,ilr,conv,ic,jlr,js)
   do ic0_loc=1,nc0_loc(mpl%myproc)
      ! MPI offset
      ic0 = ic0_s(mpl%myproc)+ic0_loc-1

      if (ndata%geom%mask(ic0,il0)) then
         ! Allocation
         allocate(is_list(ineh(ic0,il0i)*inev(il0)*maxval(ines)))
         allocate(order(ineh(ic0,il0i)*inev(il0)*maxval(ines)))
         allocate(S_list(ineh(ic0,il0i)*inev(il0)*maxval(ines)))
         if (nam%lsqrt) then
            allocate(S_list_tmp(ndata%ns))
            allocate(valid_list_tmp(ndata%ns))
         else
            allocate(S_list_tmp(ineh(ic0,il0i)*inev(il0)*maxval(ines)))
         end if

         ! Initialization
         is_list = 0
         S_list = 0.0

         ! Adjoint interpolation
         nlr = 0
         do ih=1,ineh(ic0,il0i)
            ic1 = interp_h(ih,ic0,il0i)
            do iv=1,inev(il0)
               il1 = interp_v(iv,il0)
               if ((ndata%il1_to_il0(il1)>=ndata%vbot(ic1)).and.(ndata%il1_to_il0(il1)<=ndata%vtop(ic1))) then
                  do is=1,ines(ic1,il1)
                     is_add = interp_s(is,ic1,il1)
                     S_add = interp_h_S(ih,ic0,il0i)*interp_v_S(iv,il0)*interp_s_S(is,ic1,il1)
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
               end if
            end do
         end do

         if (nam%lsqrt) then
            ! Internal normalization
            do ilr=1,nlr
               is = is_list(ilr)
               S_list(ilr) = S_list(ilr)*ndata%norm_sqrt(is)
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
            ndata%norm(ic0,il0) = 0.0
            do js=1,ndata%ns
               if (valid_list_tmp(js)) then
                  ndata%norm(ic0,il0) = ndata%norm(ic0,il0)+S_list_tmp(js)**2
               end if
            end do

            ! Normalization factor
            ndata%norm(ic0,il0) = 1.0/sqrt(ndata%norm(ic0,il0))
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
            if (nlr>0) ndata%norm(ic0,il0) = 1.0/sqrt(sum(S_list(1:nlr)*S_list_tmp(1:nlr)))
         end if

         ! Print progression
         done((il0-1)*nc0_loc(mpl%myproc)+ic0_loc) = .true.
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

   ! Communication
   if (mpl%main) then 
      do iproc=1,mpl%nproc
         if (iproc/=mpl%ioproc) then
            ! Receive data on ioproc
            call mpl_recv(nc0_loc(iproc),ndata%norm(ic0_s(iproc):ic0_e(iproc),il0),iproc,mpl%tag)
         end if
      end do
   else
      ! Send data to ioproc
      call mpl_send(nc0_loc(mpl%myproc),ndata%norm(ic0_s(mpl%myproc):ic0_e(mpl%myproc),il0),mpl%ioproc,mpl%tag)
   end if
   mpl%tag = mpl%tag+1

   ! Broadcast data
   call mpl_bcast(ndata%norm(:,il0),mpl%ioproc)
end do
write(mpl%unit,'(a)') '100%'

! Release memory
deallocate(done)

end subroutine compute_normalization

end module module_normalization
