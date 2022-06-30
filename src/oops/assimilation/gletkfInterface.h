/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_GLETKFINTERFACE_H_
#define OOPS_ASSIMILATION_GLETKFINTERFACE_H_

namespace oops {

  extern "C" {
    void letkf_core_f90(const int &, const float *, const float *, const float *,
                        float *, float *,
                        const float *, const int &, const int &,
                        const int &, const int &, const int &, const float &);
  }

}  // namespace oops
#endif  // OOPS_ASSIMILATION_GLETKFINTERFACE_H_
