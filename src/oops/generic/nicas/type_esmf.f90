!----------------------------------------------------------------------
! Module: type_esmf
!> Purpose: ESMF derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-B license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_esmf

use iso_c_binding
use esmf
use netcdf
use tools_const, only: pi,rad2deg
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsi
use tools_nc, only: ncerr
use type_linop, only: linoptype,linop_alloc
use type_ndata, only: ndatatype
use type_randgen, only: rand_integer

implicit none

! ESMF parameters
type(esmf_typekind_flag) :: esmffloat !< ESMF type for floats

interface esmf_create_field 
   module procedure esmf_create_field_from_mesh
   module procedure esmf_create_field_from_grid
end interface

private
public :: esmffloat,esmf_start,esmf_end,esmf_check,esmf_create_field,esmf_create_interp

contains

!----------------------------------------------------------------------
! Subroutine: esmf_start
!> Purpose: initialize esmf
!----------------------------------------------------------------------
subroutine esmf_start() bind(c, name='esmf_start_f90')

implicit none

! Local variables
integer :: info

! Start ESMF
call esmf_initialize(defaultcalkind=esmf_calkind_gregorian,rc=info)
call esmf_check(info)

! Check float size
if (kind_real==4) then
   esmffloat = esmf_typekind_r4
elseif (kind_real==8) then
   esmffloat = esmf_typekind_r8
else
   call msgerror('unknown real kind for ESMF floats')
end if

end subroutine esmf_start

!----------------------------------------------------------------------
! Subroutine: esmf_end
!> Purpose: finalize esmf
!----------------------------------------------------------------------
subroutine esmf_end() bind(c, name='esmf_end_f90')

implicit none

! Finalize ESMF
call esmf_finalize(endflag=esmf_end_keepmpi)

end subroutine esmf_end

!----------------------------------------------------------------------
! Subroutine: esmf_check
!> Purpose: check esmf
!----------------------------------------------------------------------
subroutine esmf_check(info)

implicit none

! Passed variables
integer,intent(in) :: info !< Returned code

! Check ESMF code returned
if (info/=esmf_success) then
   call esmf_finalize(endflag=esmf_end_keepmpi)
   call msgerror('ESMF problem, see log file for details')
end if

end subroutine esmf_check

!----------------------------------------------------------------------
! Subroutine: esmf_create_field_from_mesh
!> Purpose: create an ESMF field with an appropriate mesh
!----------------------------------------------------------------------
subroutine esmf_create_field_from_mesh(ndata,n,lon,lat,meshmask,field)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata   !< Sampling data
integer,intent(in) :: n               !< Number of nodes
real(kind_real),intent(in) :: lon(n)  !< Points longitudes
real(kind_real),intent(in) :: lat(n)  !< Points latitudes
logical,intent(in) :: meshmask(n)     !< Masked points
type(esmf_field),intent(out) :: field !< ESMF field

! Local variables
integer :: lnew,info,nt,it,i,j,k
integer :: order(n),list(6*(n-2)),lptr(6*(n-2)),lend(n),near(n),next(n),ltri(6,2*(n-2))
integer :: nodeids(n),nodeowners(n),nodemask(n)
integer,allocatable :: elementids(:),elementtypes(:),elementconn(:,:),elementmask(:)
real(kind_real) :: xc,yc,zc,normc,latc,lonc,rad,v1(3),v2(3),v3(3),areas
real(kind_real) :: x(n),y(n),z(n),dist(n)
real(kind_real) :: nodecoords(2,n)
real(kind_real),allocatable :: elementcoords(:,:),elementarea(:)
type(esmf_mesh) :: mesh
type(esmf_arrayspec) :: arrayspec

! Shuffle arrays (more efficient to compute the Delaunay triangulation)
do i=1,n
   order(i) = i
end do
do i=n,2,-1
   call rand_integer(ndata%rng,1,n,j)
   k = order(j)
   order(j) = order(i)
   order(i) = k
end do

! Transform to cartesian coordinates
call trans(n,lat(order),lon(order),x,y,z)

! Create mesh
list = 0
call trmesh(n,x,y,z,list,lptr,lend,lnew,near,next,dist,info)

! Create triangles list
call trlist(n,list,lptr,lend,6,nt,ltri,info)

! Allocation
allocate(elementcoords(2,nt))
allocate(elementarea(nt))
allocate(elementids(nt))
allocate(elementtypes(nt))
allocate(elementconn(3,nt))
allocate(elementmask(nt))

! Nodes
nodecoords(1,:) = lon
nodecoords(2,:) = lat    
do i=1,n
   nodeids(i) = i
   if (meshmask(i)) then
      nodemask(i) = 1
   else
      nodemask(i) = 0
   end if
end do
nodeowners = 0

! Elements
do it=1,nt
   elementconn(:,it) = order(ltri(1:3,it))
   xc = sum(x(ltri(1:3,it)))
   yc = sum(y(ltri(1:3,it)))
   zc = sum(z(ltri(1:3,it)))
   normc = sqrt(xc**2+yc**2+zc**2)
   xc = xc/normc
   yc = yc/normc
   zc = zc/normc
   call scoord(xc,yc,zc,latc,lonc,rad)
   elementcoords(1,it) = lonc
   elementcoords(2,it) = latc
   v1 = (/x(ltri(1,it)),y(ltri(1,it)),z(ltri(1,it))/)
   v2 = (/x(ltri(2,it)),y(ltri(2,it)),z(ltri(2,it))/)
   v3 = (/x(ltri(3,it)),y(ltri(3,it)),z(ltri(3,it))/)
   elementarea(it) = areas(v1,v2,v3)
   elementids(it) = it
   elementtypes(it) = esmf_meshelemtype_tri
   if (all(meshmask(ltri(1:3,it)))) then
      elementMask(it) = 1
   else
      elementMask(it) = 0
   end if
end do

! Create ESMF mesh
mesh = esmf_meshcreate(parametricdim=2,spatialdim=2,nodeids=nodeids,nodecoords=pack(nodeCoords,mask=.true.), &
 & nodeowners=nodeowners,nodemask=nodemask,elementids=elementids,elementtypes=elementtypes, &
 & elementconn=pack(elementconn,mask=.true.),elementmask=elementmask,elementarea=elementarea, &
 & elementcoords=pack(elementcoords,mask=.true.),coordsys=esmf_coordsys_sph_rad,rc=info)
call esmf_check(info)

! Create arrayspec
call esmf_arrayspecset(arrayspec,1,esmffloat,rc=info)
call esmf_check(info)

! Create field
field = esmf_fieldcreate(mesh,arrayspec,meshloc=esmf_meshloc_node,rc=info)
call esmf_check(info)

! Validate field
call esmf_fieldvalidate(field,rc=info)
call esmf_check(info)

end subroutine esmf_create_field_from_mesh

!----------------------------------------------------------------------
! Subroutine: esmf_create_field_from_grid
!> Purpose: create an ESMF field
!----------------------------------------------------------------------
subroutine esmf_create_field_from_grid(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: imask2d(ndata%nlon,ndata%nlat),ilon,ilat,info
real(kind_real) :: lon2d(ndata%nlon,ndata%nlat)
real(kind_real) :: lat2d(ndata%nlon,ndata%nlat)
logical :: mask_unpack(ndata%nlon,ndata%nlat),lmask2d(ndata%nlon,ndata%nlat)
type(esmf_grid) :: grid
type(esmf_array) :: array
type(esmf_arrayspec) :: arrayspec

! Initialization
mask_unpack = .true.
lmask2d = .false.
call msr(lon2d)
call msr(lat2d)

! Unpack
lmask2d = unpack(any(ndata%mask,dim=2),mask=mask_unpack,field=lmask2d)
lon2d = unpack(ndata%lon,mask=mask_unpack,field=lon2d)
lat2d = unpack(ndata%lat,mask=mask_unpack,field=lat2d)

! Logical to array
imask2d = 0
do ilon=1,ndata%nlon
   do ilat=1,ndata%nlat
      if (lmask2d(ilon,ilat)) imask2d(ilon,ilat) = 1
   end do
end do

! Create grid
grid = esmf_gridcreate1peridim(minindex=(/1,1/),maxindex=(/ndata%nlon,ndata%nlat/),regdecomp=(/1,1/), &
 & coordsys=esmf_coordsys_sph_rad,gridedgelwidth=(/0,0/),gridedgeuwidth=(/0,1/), &
 & indexflag=esmf_index_global,rc=info)
call esmf_check(info)

! Add coordinate
call esmf_gridaddcoord(grid,staggerloc=esmf_staggerloc_center,rc=info)
call esmf_check(info)

! Add mask
call esmf_gridadditem(grid,staggerloc=esmf_staggerloc_center,itemflag=esmf_griditem_mask,rc=info)
call esmf_check(info)

! Set longitude coordinate
call esmf_gridgetcoord(grid,staggerloc=esmf_staggerloc_center,coorddim=1,array=array,rc=info)
call esmf_check(info)
call esmf_arrayscatter(array,lon2d,rootpet=0,rc=info)
call esmf_check(info)

! Set latitude coordinate
call esmf_gridgetcoord(grid,staggerloc=esmf_staggerloc_center,coordDim=2,array=array,rc=info)
call esmf_check(info)
call esmf_arrayscatter(array,lat2d,rootpet=0,rc=info)
call esmf_check(info)

! Set mask
call esmf_gridgetitem(grid,staggerloc=esmf_staggerloc_center,itemflag=esmf_griditem_mask,array=array,rc=info)
call esmf_check(info)
call esmf_arrayscatter(array,imask2d,rootpet=0,rc=info)
call esmf_check(info)

! Create arrayspec
call esmf_arrayspecset(arrayspec,2,esmffloat,rc=info)
call esmf_check(info)

! Create field
ndata%c0field = esmf_fieldcreate(grid,arrayspec,staggerloc=esmf_staggerloc_center,rc=info)
call esmf_check(info)

! Validate field
call esmf_fieldvalidate(ndata%c0field,rc=info)
call esmf_check(info)

end subroutine esmf_create_field_from_grid

!----------------------------------------------------------------------
! Subroutine: esmf_create_interp
!> Purpose: create an ESMF interpolation
!----------------------------------------------------------------------
subroutine esmf_create_interp(srcfield,dstfield,regionalflag,interp)

implicit none

! Passed variables
type(esmf_field),intent(in)   :: srcfield    !< Source field
type(esmf_field),intent(inout)   :: dstfield !< Destination field
logical,intent(in) :: regionalflag           !< Regional flag
type(linoptype),intent(inout) :: interp      !< Interpolation

! Local variables
integer :: regridscheme,petno,petcnt,total,localcount(1),maxcount,start,i,j,info
integer(esmf_kind_i4) :: maskvals(1)
integer(esmf_kind_i4),pointer :: factorindexlist(:,:)
integer(esmf_kind_i4),pointer :: allcounts(:) 
integer(esmf_kind_i4),pointer :: indexbuf(:)
real(esmf_kind_r8),pointer :: factorlist(:)
real(esmf_kind_r8),pointer :: weightbuf(:)
type(esmf_polemethod_flag) :: polemethod
type(esmf_vm) :: vm

! Get global vm information
call esmf_vmgetcurrent(vm,rc=info)
call esmf_check(info)

! Set up local pet info
call esmf_vmget(vm,localpet=petno,petcount=petcnt,rc=info)
call esmf_check(info)

! Regoinal options
if (regionalflag) then
   polemethod = esmf_polemethod_none
   regridscheme = esmf_regrid_scheme_region3d
else
   polemethod = esmf_polemethod_allavg
   regridscheme = esmf_regrid_scheme_full3d
end if

! Mask value
maskvals(1) = 0

! Genereate interpolation weights
call esmf_fieldregridstore(srcfield=srcfield,dstfield=dstfield,srcmaskvalues=maskvals, &
 & dstmaskvalues=maskvals,unmappedaction=esmf_unmappedaction_ignore, & 
 & factorindexList=factorindexList,factorlist=factorlist, &
 & regridmethod=esmf_regridmethod_bilinear,polemethod=polemethod,regridpolenpnts=1,rc=info)
call esmf_check(info)

! Gather data size
localcount(1) = size(factorlist,1)
allocate(allcounts(petcnt))
call esmf_vmallgather(vm,localcount,allcounts,1,rc=info)
call esmf_check(info)
total = 0
do i=1,petcnt
   total = allcounts(i)+total
end do
maxcount = 0
do i=1,petcnt
   if (allcounts(i) > maxcount) maxcount = allcounts(i)
end do

! Allocation
interp%n_s = total
call linop_alloc(interp)

if (petno==0) then 
   ! Allocation
   allocate(indexbuf(maxcount*2), weightbuf(maxcount))

   ! Loop over PETs
   start = 1
   do i=1,petcnt
      localcount(1) = allcounts(i)
      if (i==1) then 
         interp%col(start:start+localcount(1)-1) = factorindexlist(1,start:start+localcount(1)-1)
         interp%row(start:start+localcount(1)-1) = factorindexlist(2,start:start+localcount(1)-1)
         interp%S(start:start+localcount(1)-1) = factorlist(start:start+localcount(1)-1)
      else
         if (localcount(1)>0) then
            ! Receive the factorlist and factorindexlist 
            call esmf_vmrecv(vm,indexbuf,2*localcount(1),i-1,rc=info)
            call esmf_check(info)
            call esmf_vmrecv(vm,weightbuf,localcount(1),i-1,rc=info)
            call esmf_check(info)
            interp%col(start:start+localcount(1)-1) = indexbuf(1:localcount(1))
            interp%row(start:start+localcount(1)-1) = indexbuf(localcount(1)+1:2*localcount(1))
            interp%S(start:start+localcount(1)-1) = weightbuf(1:localcount(1))
         end if
      end if
      start = start + localcount(1)
   end do
else
   ! Allocation
   allocate(indexbuf(2*localcount(1)))
   
   if (localcount(1)>0) then
      ! Fill index buffer
      do j=1,localcount(1)
         indexbuf(j) = factorindexlist(1,j)
         indexbuf(j+localcount(1)) = factorindexlist(2,j)
      end do

      ! A non-root PET, send the results to PET 0
      call esmf_vmsend(vm,indexbuf,2*localcount(1),0,rc=info)
      call esmf_check(info)
      call esmf_vmsend(vm,factorList,localcount(1),0,rc=info)
      call esmf_check(info)
   end if
end if

! Broadcast
call esmf_vmbroadcast(vm,interp%row,interp%n_s,0,rc=info)
call esmf_check(info)
call esmf_vmbroadcast(vm,interp%col,interp%n_s,0,rc=info)
call esmf_check(info)
call esmf_vmbroadcast(vm,interp%S,interp%n_s,0,rc=info)
call esmf_check(info)

end subroutine esmf_create_interp

end module type_esmf
