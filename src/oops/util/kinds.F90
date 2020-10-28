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
  public kind_float, kind_real, kind_single, kind_quad, kind_int, kind_long
  
  integer, parameter :: kind_single = c_float
  integer, parameter :: kind_float = kind_single
  integer, parameter :: kind_real   = c_double
  integer, parameter :: kind_quad   = c_long_double

  integer, parameter :: kind_int    = c_int
  integer, parameter :: kind_long   = c_long
end module kinds
