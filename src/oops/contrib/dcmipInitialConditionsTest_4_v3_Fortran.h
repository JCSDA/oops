/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#ifndef OOPS_CONTRIB_DCMIPINITIALCONDITIONSTEST_4_V3_FORTRAN_H_
#define OOPS_CONTRIB_DCMIPINITIALCONDITIONSTEST_4_V3_FORTRAN_H_

namespace oops {

extern "C" {

void test4_baroclinic_wave_f90(const int & moist, const double & X,     // IN
                               const double & lon, const double & lat,  // IN
                               double & p,                              // INOUT
                               const double & z, const int & zcoords,   // IN
                               double & u, double & v, double & w,      // OUT
                               double & t, double & phis, double & ps,  // OUT
                               double & rho, double & q,                // OUT
                               double & q1, double & q2);               // OUT

}

}  // namespace oops

#endif  // OOPS_CONTRIB_DCMIPINITIALCONDITIONSTEST_4_V3_FORTRAN_H_

