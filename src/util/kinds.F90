! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module defining the kind of reals for the l95 and qg model

module kinds
  use, intrinsic :: iso_c_binding
  implicit none

  private
  public kind_real
  
  integer, parameter :: kind_real=c_double
end module kinds
