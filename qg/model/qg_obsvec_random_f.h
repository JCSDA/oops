/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_QG_OBSVEC_RANDOM_F_H_
#define QG_MODEL_QG_OBSVEC_RANDOM_F_H_

namespace qg {
  class ObsSpaceQG;

extern "C" {
  void qg_obsvec_random_f(const ObsSpaceQG &, const int &, double *);
}

}  // namespace qg

#endif  // QG_MODEL_QG_OBSVEC_RANDOM_F_H_
