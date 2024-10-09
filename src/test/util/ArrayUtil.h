/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <netcdf.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <set>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/utils/Hash.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/base/Variable.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#define ERR2(e) {ABORT(nc_strerror(e));}

namespace test {

using atlas::idx_t;

/// \detail This gather scatter tests checks the gatherSum
///         and scatter by using the two routines as part of an
///         adjoint test. For the gatherSum to be an adjoint of the scatter
///         one needs to set the array values on pes that are not the "root" to
///         to zero after applying the gatherSum.
///         (The gatherSum is doing an mpi gather on the root PE and then summing the
///         contributions from all PEs on this PE.)
///         The symmetric adjoint test is being used here.
///           scatter(fld) . scatter (fld) =
///              fld .(zero contributions on non root pes (gatherSum (scatter(fld))))
void test_gatherscatter_double_1D() {
  std::array<atlas::idx_t, 1> shape{10};

  // allocate fields for test
  auto Fld = atlas::Field(std::string("1D_double"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(shape[0]));

  auto Fld2 = atlas::Field(std::string("1D_double_2"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(shape[0]));

  // first testing that gatherSum works as intended.
  auto fview = atlas::array::make_view<double, 1>(Fld);
  fview.assign(1.0);

  std::size_t root(0);
  util::gatherSum<double>(oops::mpi::world(), root, fview);
  if (oops::mpi::world().rank() == root) {
    EXPECT(oops::is_close_absolute(fview(0),
                                   static_cast<double>(oops::mpi::world().size()),
                                   1.0e-15));
  }

  // now starting adjoint test
  auto f2view = atlas::array::make_view<double, 1>(Fld2);
  f2view.assign(fview);

  util::scatter<double>(oops::mpi::world(), root, fview);
  EXPECT(oops::is_close_absolute(fview(0),
                                 static_cast<double>(oops::mpi::world().size()),
                                 1.0e-15));
  double val(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    val += fview(t) * fview(t);
  }
  util::gatherSum<double>(oops::mpi::world(), root, fview);

  double val2(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    fview(t) = oops::mpi::world().rank() == root ? fview(t) : 0.0;
      // above line is to make the gather the adjoint of the scatter.
    val2 += f2view(t) * fview(t);
  }

  // adjoint test is across at PEs - so need to sum across all PEs
  // for adjoint test to work
  oops::mpi::world().allReduceInPlace(val, eckit::mpi::sum());
  oops::mpi::world().allReduceInPlace(val2, eckit::mpi::sum());

  EXPECT(oops::is_close_absolute(val, val2, 1.0e-15));
}

void test_gatherscatter_double_2D() {
  std::array<atlas::idx_t, 2> shape{10, 20};

  // allocate fields
  auto Fld = atlas::Field(std::string("2D_double"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(shape[0], shape[1]));

  auto Fld2 = atlas::Field(std::string("2D_double_2"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(shape[0], shape[1]));

  // check gatherSum is correct
  auto fview = atlas::array::make_view<double, 2>(Fld);
  fview.assign(1.0);

  std::size_t root(0);
  util::gatherSum<double>(oops::mpi::world(), root, fview);
  if (oops::mpi::world().rank() == root) {
    EXPECT(oops::is_close_absolute(fview(0, 0),
                                   static_cast<double>(oops::mpi::world().size()),
                                   1.0e-15));
  }

  // apply adjoint test
  auto f2view = atlas::array::make_view<double, 2>(Fld2);
  f2view.assign(fview);

  util::scatter<double>(oops::mpi::world(), root, fview);
  EXPECT(oops::is_close_absolute(fview(0, 0),
                                 static_cast<double>(oops::mpi::world().size()),
                                 1.0e-15));
  double val(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < fview.shape()[1]; ++t1) {
      val += fview(t, t1) * fview(t, t1);
    }
  }
  util::gatherSum<double>(oops::mpi::world(), root, fview);

  double val2(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < fview.shape()[1]; ++t1) {
      fview(t, t1) = oops::mpi::world().rank() == root ? fview(t, t1) : 0.0;
      // above line is to make the gather the adjoint of the scatter.
      val2 += f2view(t, t1) * fview(t, t1);
    }
  }

  oops::mpi::world().allReduceInPlace(val, eckit::mpi::sum());
  oops::mpi::world().allReduceInPlace(val2, eckit::mpi::sum());

  EXPECT(oops::is_close_absolute(val, val2, 1.0e-15));
}

void test_gatherscatter_double_3D() {
  std::array<atlas::idx_t, 3> shape{10, 20, 30};

  // allocate fields
  auto Fld = atlas::Field(std::string("3D_double"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(shape[0], shape[1], shape[2]));

  auto Fld2 = atlas::Field(std::string("3D_double_2"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(shape[0], shape[1], shape[2]));

  // check gatherSum is correct
  auto fview = atlas::array::make_view<double, 3>(Fld);
  fview.assign(1.0);

  std::size_t root(0);
  util::gatherSum<double>(oops::mpi::world(), root, fview);
  if (oops::mpi::world().rank() == root) {
    EXPECT(oops::is_close_absolute(fview(0, 0, 0),
                                   static_cast<double>(oops::mpi::world().size()),
                                   1.0e-15));
  }

  // apply adjoint test
  auto f2view = atlas::array::make_view<double, 3>(Fld2);
  f2view.assign(fview);

  util::scatter<double>(oops::mpi::world(), root, fview);
  EXPECT(oops::is_close_absolute(fview(0, 0, 0),
                                 static_cast<double>(oops::mpi::world().size()),
                                 1.0e-15));

  double val(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < fview.shape()[1]; ++t1) {
      for (idx_t t2 = 0; t2 < fview.shape()[2]; ++t2) {
        val += fview(t, t1, t2) * fview(t, t1, t2);
      }
    }
  }
  util::gatherSum<double>(oops::mpi::world(), root, fview);

  double val2(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < fview.shape()[1]; ++t1) {
      for (idx_t t2 = 0; t2 < fview.shape()[2]; ++t2) {
        fview(t, t1, t2) = oops::mpi::world().rank() == root ? fview(t, t1, t2) : 0.0;
        // above line is to make the gather the adjoint of the scatter.
        val2 += f2view(t, t1, t2) * fview(t, t1, t2);
      }
    }
  }

  oops::mpi::world().allReduceInPlace(val, eckit::mpi::sum());
  oops::mpi::world().allReduceInPlace(val2, eckit::mpi::sum());

  EXPECT(oops::is_close_absolute(val, val2, 1.0e-15));
}

void test_gatherscatter_int_2D() {
  std::array<atlas::idx_t, 2> shape{10, 20};

  // allocate fields
  auto Fld = atlas::Field(std::string("2D_int"),
    atlas::array::make_datatype<int>(),
    atlas::array::make_shape(shape[0], shape[1]));

  auto Fld2 = atlas::Field(std::string("2D_int_2"),
    atlas::array::make_datatype<int>(),
    atlas::array::make_shape(shape[0], shape[1]));

  // check gatherSum is correct
  auto fview = atlas::array::make_view<int, 2>(Fld);
  fview.assign(1.0);

  std::size_t root(0);
  util::gatherSum<int>(oops::mpi::world(), root, fview);
  if (oops::mpi::world().rank() == root) {
    EXPECT(fview(0, 0) == static_cast<int>(oops::mpi::world().size()));
  }

  // apply adjoint test
  auto f2view = atlas::array::make_view<int, 2>(Fld2);
  f2view.assign(fview);

  util::scatter<int>(oops::mpi::world(), root, fview);
  EXPECT(fview(0, 0) == static_cast<int>(oops::mpi::world().size()));

  int val(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < fview.shape()[1]; ++t1) {
      val += fview(t, t1) * fview(t, t1);
    }
  }
  util::gatherSum<int>(oops::mpi::world(), root, fview);

  int val2(0.0);
  for (idx_t t = 0; t < fview.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < fview.shape()[1]; ++t1) {
      fview(t, t1) = oops::mpi::world().rank() == root ? fview(t, t1) : 0;
      // above line is to make the gather the adjoint of the scatter.
      val2 += f2view(t, t1) * fview(t, t1);
    }
  }

  oops::mpi::world().allReduceInPlace(val, eckit::mpi::sum());
  oops::mpi::world().allReduceInPlace(val2, eckit::mpi::sum());

  EXPECT(val == val2);
}

void test_write_read_double_1D() {
  // 0. initialise data
  // 1. create a 1D field
  // 2. gather
  // 3. write to file on root PE.
  // 4. copy original field and reassign values to zero.
  // 5. read from file on root PE.
  // 6. compare with original fieldset

  // 0. Setup
  std::array<atlas::idx_t, 1> shape{10};
  std::size_t root(0);
  std::string ncfilepath{"./out_double_1D.nc"};
  std::vector<std::string> dimNames;
  std::vector<atlas::idx_t> dimSizes;
  std::vector<std::vector<std::string>> dimNamesForEveryVar;
  oops::Variables variables;

  //  the netcdf metadata in this test has the following structure
  //  root:    //   ie. NetCDF metadata root of local configuration
  //    global metadata:  // special key name for global meta data
  //    - string attribute:
  //        value: "my value"
  //        type: "string"
  //    - int32 attribute:
  //        value: 32
  //        type: "int32"
  //    - int64 attribute:
  //        value:  64
  //        type: "int64"
  //    - real32 attribute:
  //        value: 10.0f
  //        type: "real32"
  //    - real32 attribute 2:
  //        value: 100.0f
  //        type: "real32
  //    - real64 attribute:
  //        value: 20.0f
  //        type: "real64"
  //    1Ddouble:
  //    - string variable attribute:
  //        value: "my variable value"
  //        type: "string"
  //    - int32 variable attribute:
  //        value: 320
  //        type: "int32"
  //    - int64 variable attribute:
  //        value: 9147483647
  //        type: "int64"
  //    - real32 variable attribute:
  //        value: 10.1f
  //        type: "real32"
  //    - real32 variable attribute 2:
  //        value: 100.0f
  //        type: "real32"
  //    - real64 variable attribute:
  //        value: 20.1
  //        type: "real64"
  eckit::LocalConfiguration conf;
  util::setAttribute<std::string>(conf, "global metadata", "string attribute",
                                  "string", "my value");
  util::setAttribute<std::int32_t>(conf, "global metadata", "int32 attribute",
                                   "int32", static_cast<std::int32_t>(32));
  util::setAttribute<std::int64_t>(conf, "global metadata", "int64 attribute",
                                   "int64", static_cast<std::int64_t>(64));
  util::setAttribute<float>(conf, "global metadata", "real32 attribute",
                            "real32", float{10.0});
  util::setAttribute<float>(conf, "global metadata", "real32 attribute 2",
                            "real32", float{100.0});
  util::setAttribute<double>(conf, "global metadata", "real64 attribute",
                            "real64", double{20.0});
  util::setAttribute<std::string>(conf, "1Ddouble", "string variable attribute",
                                  "string", "my variable value");
  util::setAttribute<std::int32_t>(conf, "1Ddouble", "int32 variable attribute",
                                   "int32", static_cast<std::int32_t>(320));
  util::setAttribute<std::int64_t>(conf, "1Ddouble", "int64 variable attribute",
                                   "int64", static_cast<std::int64_t>(9147483647));
  util::setAttribute<float>(conf, "1Ddouble", "real32 variable attribute",
                            "real32", 10.1f);
  util::setAttribute<float>(conf, "1Ddouble", "real32 variable attribute 2",
                            "real32", 100.0f);
  util::setAttribute<double>(conf, "1Ddouble", "real64 variable attribute",
                            "real64", 20.1);
  oops::Log::info() << "conf input = " << conf << std::endl;

  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;
  std::string fieldname{"1Ddouble"};
  dimNames.push_back("index1");
  dimSizes.push_back(shape[0]);
  variables.push_back(fieldname);

  for (auto var : variables.variables()) {
    dimNamesForEveryVar.push_back(dimNames);
  }

  // 1. create a 1D field
  auto Fld = atlas::Field(std::string(fieldname),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape({shape[0]}));
  auto fview = atlas::array::make_view<double, 1>(Fld);
  fview.assign(1.0);


  // 2. gather
  util::gatherSum<double>(oops::mpi::world(), root, fview);

  // 3. write to file on root PE.
  if (oops::mpi::world().rank() == root) {
    auto fview_const = atlas::array::make_view<const double, 1>(Fld);

    oops::Log::info() << "write_read_double_1D::confToFile = " << conf << std::endl;
    util::atlasArrayWriteHeader(ncfilepath,
                                dimNames,
                                dimSizes,
                                variables,
                                dimNamesForEveryVar,
                                conf,
                                netcdfGeneralIDs,
                                netcdfDimIDs,
                                netcdfVarIDs,
                                netcdfDimVarIDs);

    util::atlasArrayWriteData(netcdfGeneralIDs,
                              netcdfVarIDs[0],
                              fview_const);

    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ERR2(retval);
  }

  // 4. copy original field and reassign values to zero.
  auto Fld2 = atlas::Field(std::string("1Ddouble2"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape({shape[0]}));
  auto f2view = atlas::array::make_view<double, 1>(Fld2);
  f2view.assign(0.0);

  // 5. read from file on root PE.
  if (oops::mpi::world().rank() == root) {
    eckit::LocalConfiguration confFromFile;
    oops::Variables variablesFromFile;
    util::atlasArrayInquire(ncfilepath,
                            dimNames,
                            dimSizes,
                            variablesFromFile,
                            dimNamesForEveryVar,
                            confFromFile,
                            netcdfGeneralIDs,
                            netcdfDimIDs,
                            netcdfVarIDs,
                            netcdfDimVarIDs);

    oops::Log::info() << "write_read_double_1D::confFromFile = " << confFromFile << std::endl;
    for (const auto & var : variablesFromFile) {
      oops::Log::info() << "write_read_double_1D::"  << var.name() << " "
                        << var.metaData() << std::endl;
    }

    util::atlasArrayReadData(netcdfGeneralIDs,
                             dimSizes,
                             netcdfVarIDs[0],
                             f2view);
    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ERR2(retval);

    // 6. compare with original field
    for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
      EXPECT(fview(t) == f2view(t));
    }
  }
}

void test_write_read_multiple() {
  // 0. setup
  // 1. create a fieldset with fields supporting existing interface.
  // 2. gather
  // 3. write to file on root PE.
  // 4. copy original fieldset and reassign values to zero.
  // 5. read from file on root PE.
  // 6. compare with original fieldset

  // 0. setup
  std::size_t root(1);
  std::string ncfilepath{"./out_multi.nc"};
  std::vector<std::string> dimNames{"index1", "index2", "index3", "index4"};
  std::vector<atlas::idx_t> dimSizes{10, 20, 30, 40, 50};
  oops::Variables variables;
  for (const std::string & name :
       std::vector<std::string>{"double_1", "double_2", "double_3", "int_1"}) {
    oops::VariableMetaData varMeta;
    if ( name.compare(0, 3, "int") == 0 ) {
      varMeta = oops::VariableMetaData(oops::VerticalStagger::CENTER, oops::ModelDataType::Int32);
    } else {
      varMeta = oops::VariableMetaData(oops::VerticalStagger::CENTER, oops::ModelDataType::Real64);
    }
    variables.push_back(oops::Variable(name, varMeta));
  }

  std::vector<std::vector<std::string>>
    dimNamesForEveryVar{{"index1"},
                        {"index1", "index2"},
                        {"index1", "index3", "index4"},
                        {"index2"}};
  std::vector<std::vector<atlas::idx_t>>
    dimSizesForEveryVar{{10},
                        {10, 20},
                        {10, 30, 40},
                        {20}};

  //  The netcdfMetaData for this test has the following structure
  //  root:    //   ie. NetCDF metadata root of local configuration
  //    global metadata:  // special key name for global meta data
  //    - int64 attribute:
  //        value: 64
  //        type: "int64"
  //    double_1:
  //    - int64 attribute:
  //        value: 128
  //        type: "int64"
  //    double_2:
  //    - int64 attribute:
  //        value: 256
  //        type: "int64"
  //    double_3:
  //    - int64 attribute:
  //        value: 512
  //        type: "int64"
  //    int1_1:
  //    - int64 attribute:
  //        value: 1024
  //        type: "int64"
  eckit::LocalConfiguration netcdfMetaData;
  util::setAttribute<std::int64_t>(netcdfMetaData, "global metadata", "int64 attribute",
                                  "int64", static_cast<std::int64_t>(64));
  util::setAttribute<std::int64_t>(netcdfMetaData, "double_1", "int64 attribute",
                                  "int64", static_cast<std::int64_t>(128));
  util::setAttribute<std::int64_t>(netcdfMetaData, "double_2", "int64 attribute",
                                  "int64", static_cast<std::int64_t>(256));
  util::setAttribute<std::int64_t>(netcdfMetaData, "double_3", "int64 attribute",
                                  "int64", static_cast<std::int64_t>(512));
  util::setAttribute<std::int64_t>(netcdfMetaData, "int1_1", "int64 attribute",
                                  "int64", static_cast<std::int64_t>(1024));

  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;

  // 1. create a fieldset with fields supporting existing interface.
  atlas::FieldSet fset;
  std::size_t t(0);

  for (auto & var : variables) {
    if (var.name().compare(0, 6, "double") == 0) {
      if (dimSizesForEveryVar[t].size() == 1) {
        auto Fld = atlas::Field(var.name(),
          atlas::array::make_datatype<double>(),
          atlas::array::make_shape({dimSizesForEveryVar[t][0]}));
        auto fview = atlas::array::make_view<double, 1>(Fld);
        util::UniformDistribution<double> x(dimSizesForEveryVar[t][0],
                                            1.0, 2.0);
        for (idx_t t0 = 0; t0 < static_cast<idx_t>(dimSizesForEveryVar[t][0]); ++t0) {
          fview(t0) = x[static_cast<std::size_t>(t0)];
        }
        fset.add(Fld);
      } else if (dimSizesForEveryVar[t].size() == 2) {
        auto Fld = atlas::Field(var.name(),
          atlas::array::make_datatype<double>(),
          atlas::array::make_shape({dimSizesForEveryVar[t][0],
                                    dimSizesForEveryVar[t][1]}));
        auto fview = atlas::array::make_view<double, 2>(Fld);
        util::UniformDistribution<double> x(
          (dimSizesForEveryVar[t][0] * dimSizesForEveryVar[t][1]), 1.0, 2.0);
        std::size_t n(0);
        for (idx_t t0 = 0; t0 < static_cast<idx_t>(dimSizesForEveryVar[t][0]); ++t0) {
          for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizesForEveryVar[t][1]); ++t1, ++n) {
            fview(t0, t1) = x[n];
          }
        }
        fset.add(Fld);
      } else if (dimSizesForEveryVar[t].size() == 3) {
        auto Fld = atlas::Field(var.name(),
          atlas::array::make_datatype<double>(),
          atlas::array::make_shape({dimSizesForEveryVar[t][0],
                                    dimSizesForEveryVar[t][1],
                                    dimSizesForEveryVar[t][2]}));
        auto fview = atlas::array::make_view<double, 3>(Fld);
        util::UniformDistribution<double> x(
          (dimSizesForEveryVar[t][0] * dimSizesForEveryVar[t][1] *
           dimSizesForEveryVar[t][2]), 1.0, 2.0, 1.0);
        std::size_t n(0);
        for (idx_t t0 = 0; t0 < static_cast<idx_t>(dimSizesForEveryVar[t][0]); ++t0) {
          for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizesForEveryVar[t][1]); ++t1) {
            for (idx_t t2 = 0; t2 < static_cast<idx_t>(dimSizesForEveryVar[t][2]); ++t2, ++n) {
              fview(t0, t1, t2) = x[n];
            }
          }
        }
        fset.add(Fld);
      }
    } else if (var.name().compare(0, 3, "int") == 0) {
      auto fld = atlas::Field(var.name(),
        atlas::array::make_datatype<int32_t>(),
        atlas::array::make_shape({dimSizesForEveryVar[t][0]}));
      auto fview = atlas::array::make_view<int32_t, 1>(fld);
      util::UniformIntDistribution<int32_t> x(dimSizesForEveryVar[t][0],
                                              1, 10);
      for (idx_t t0 = 0; t0 < static_cast<idx_t>(dimSizesForEveryVar[t][0]); ++t0) {
        fview(t0) = x[static_cast<std::size_t>(t0)];
      }
      fset.add(fld);
    }
    ++t;
  }

  // 2. gather
  for (atlas::Field & fld : fset) {
    if (fld.datatype().str() == "real64") {
      if (fld.shape().size() == 1) {
        auto fview = atlas::array::make_view<double, 1>(fld);
        util::gatherSum<double>(oops::mpi::world(), root, fview);
      } else if (fld.shape().size() == 2) {
        auto fview = atlas::array::make_view<double, 2>(fld);
        util::gatherSum<double>(oops::mpi::world(), root, fview);
      } else if (fld.shape().size() == 3) {
        auto fview = atlas::array::make_view<double, 3>(fld);
        util::gatherSum<double>(oops::mpi::world(), root, fview);
      }
    } else if (fld.datatype().str() == "int32") {
      auto fview = atlas::array::make_view<int32_t, 1>(fld);
      util::gatherSum<int32_t>(oops::mpi::world(), root, fview);
    }
  }

  // 3. write to file on root PE.
  if (oops::mpi::world().rank() == root) {
    oops::Log::info() << "write_read_multiple::confToFile = "
                      << netcdfMetaData << std::endl;
    util::atlasArrayWriteHeader(ncfilepath,
                                dimNames,
                                dimSizes,
                                variables,
                                dimNamesForEveryVar,
                                netcdfMetaData,
                                netcdfGeneralIDs,
                                netcdfDimIDs,
                                netcdfVarIDs,
                                netcdfDimVarIDs);


    t = 0;
    for (atlas::Field & fld : fset) {
      if (fld.datatype().str() == "real64") {
        if (fld.shape().size() == 1) {
          auto fview = atlas::array::make_view<const double, 1>(fld);
          util::atlasArrayWriteData(netcdfGeneralIDs,
                                    netcdfVarIDs[t],
                                    fview);
        } else if (fld.shape().size() == 2) {
          auto fview = atlas::array::make_view<const double, 2>(fld);
          util::atlasArrayWriteData(netcdfGeneralIDs,
                                    netcdfVarIDs[t],
                                    fview);
        } else if (fld.shape().size() == 3) {
          auto fview = atlas::array::make_view<const double, 3>(fld);
          util::atlasArrayWriteData(netcdfGeneralIDs,
                                    netcdfVarIDs[t],
                                    fview);
        }
      } else if (fld.datatype().str() == "int32") {
        auto fview = atlas::array::make_view<const int32_t, 1>(fld);
        util::atlasArrayWriteData(netcdfGeneralIDs,
                                  netcdfVarIDs[t],
                                  fview);
      }
      ++t;
    }

    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ERR2(retval);
  }

  // 4. copy original fieldset and reassign values to zero.
  t = 0;
  atlas::FieldSet fset2;
  for (const oops::Variable & var : variables) {
    const std::string varname = var.name();
    if (varname.compare(0, 6, "double") == 0) {
      if (dimSizesForEveryVar[t].size() == 1) {
        auto Fld = atlas::Field(varname,
          atlas::array::make_datatype<double>(),
          atlas::array::make_shape({dimSizesForEveryVar[t][0]}));
        atlas::array::make_view<double, 1>(Fld).assign(0.0);
        fset2.add(Fld);
      } else if (dimSizesForEveryVar[t].size() == 2) {
        auto Fld = atlas::Field(varname,
          atlas::array::make_datatype<double>(),
          atlas::array::make_shape({dimSizesForEveryVar[t][0],
                                    dimSizesForEveryVar[t][1]}));
        atlas::array::make_view<double, 2>(Fld).assign(0.0);
        fset2.add(Fld);
      } else if (dimSizesForEveryVar[t].size() == 3) {
        auto Fld = atlas::Field(varname,
          atlas::array::make_datatype<double>(),
          atlas::array::make_shape({dimSizesForEveryVar[t][0],
                                    dimSizesForEveryVar[t][1],
                                    dimSizesForEveryVar[t][2]}));
        atlas::array::make_view<double, 3>(Fld).assign(0.0);
        fset2.add(Fld);
      }
    } else if (varname.compare(0, 3, "int") == 0) {
      auto fld = atlas::Field(varname,
        atlas::array::make_datatype<int>(),
        atlas::array::make_shape({dimSizesForEveryVar[t][0]}));
      atlas::array::make_view<int, 1>(fld).assign(0);
      fset2.add(fld);
    }
    ++t;
  }

  // 5. read from file on root PE.
  eckit::LocalConfiguration confFromFile;
  if (oops::mpi::world().rank() == root) {
    util::atlasArrayInquire(ncfilepath,
                            dimNames,
                            dimSizes,
                            variables,
                            dimNamesForEveryVar,
                            confFromFile,
                            netcdfGeneralIDs,
                            netcdfDimIDs,
                            netcdfVarIDs,
                            netcdfDimVarIDs);
    oops::Log::info() << "write_read_multiple::confFromFile = "
                      << confFromFile << std::endl;

    t = 0;
    for (atlas::Field & fld : fset2) {
      if (fld.datatype().str() == "real64") {
        if (fld.shape().size() == 1) {
          auto fview = atlas::array::make_view<double, 1>(fld);
          util::atlasArrayReadData(netcdfGeneralIDs,
                                   dimSizesForEveryVar[t],
                                   netcdfVarIDs[t],
                                   fview);
        } else if (fld.shape().size() == 2) {
          auto fview = atlas::array::make_view<double, 2>(fld);
          util::atlasArrayReadData(netcdfGeneralIDs,
                                   dimSizesForEveryVar[t],
                                   netcdfVarIDs[t],
                                   fview);
        } else if (fld.shape().size() == 3) {
          auto fview = atlas::array::make_view<double, 3>(fld);
          util::atlasArrayReadData(netcdfGeneralIDs,
                                   dimSizesForEveryVar[t],
                                   netcdfVarIDs[t],
                                   fview);
        }
      } else if (fld.datatype().str() == "int32") {
        auto fview = atlas::array::make_view<int, 1>(fld);
        util::atlasArrayReadData(netcdfGeneralIDs,
                                 dimSizesForEveryVar[t],
                                 netcdfVarIDs[t],
                                 fview);
      }
      ++t;
    }

    // 6. compare with original fieldset
    atlas::idx_t i(0);
    for (const std::string & var : variables.variables()) {
      std::size_t ii = static_cast<std::size_t>(i);
      if (var.compare(0, 6, "double") == 0) {
        if (dimSizesForEveryVar[ii].size() == 1) {
          auto fview = atlas::array::make_view<double, 1>(fset[i]);
          auto f2view = atlas::array::make_view<double, 1>(fset2[i]);
          for (idx_t t0 = 0;
               t0 < static_cast<idx_t>(dimSizesForEveryVar[ii][0]);
               ++t0) {
            EXPECT(oops::is_close_absolute(fview(t0),
                                           f2view(t0),
                                           1.0e-15));
          }
        } else if (dimSizesForEveryVar[ii].size() == 2) {
          auto fview = atlas::array::make_view<double, 2>(fset[i]);
          auto f2view = atlas::array::make_view<double, 2>(fset2[i]);
          for (idx_t t0 = 0; t0 < static_cast<idx_t>(dimSizesForEveryVar[ii][0]); ++t0) {
            for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizesForEveryVar[ii][1]); ++t1) {
              EXPECT(oops::is_close_absolute(fview(t0, t1),
                                             f2view(t0, t1),
                                             1.0e-15));
            }
          }
        } else if (dimSizesForEveryVar[ii].size() == 3) {
          auto fview = atlas::array::make_view<double, 3>(fset[i]);
          auto f2view = atlas::array::make_view<double, 3>(fset2[i]);
          for (idx_t t0 = 0; t0 < static_cast<idx_t>(dimSizesForEveryVar[ii][0]); ++t0) {
            for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizesForEveryVar[ii][1]); ++t1) {
              for (idx_t t2 = 0; t2 < static_cast<idx_t>(dimSizesForEveryVar[ii][2]); ++t2) {
              EXPECT(oops::is_close_absolute(fview(t0, t1, t2),
                                             f2view(t0, t1, t2),
                                             1.0e-15));
              }
            }
          }
        }
      } else if (var.compare(0, 3, "int") == 0) {
        auto fview = atlas::array::make_view<int, 1>(fset[i]);
        auto f2view = atlas::array::make_view<int, 1>(fset2[i]);
        for (idx_t t0 = 0; t0 < static_cast<idx_t>(dimSizesForEveryVar[ii][0]); ++t0) {
          EXPECT(fview(t0) == f2view(t0));
        }
      }
      ++i;
    }
  }
}


class ArrayUtil : public oops::Test{
 public:
  using oops::Test::Test;
  virtual ~ArrayUtil() {}

 private:
  std::string testid() const override {
    return "oops::test::ArrayUtil";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("oops/test_gatherscatter_double_1D")
      { test_gatherscatter_double_1D(); });

    ts.emplace_back(CASE("oops/test_gatherscatter_double_2D")
      { test_gatherscatter_double_2D(); });

    ts.emplace_back(CASE("oops/test_gatherscatter_double_3D")
      { test_gatherscatter_double_3D(); });

    ts.emplace_back(CASE("oops/test_gatherscatter_int_2D")
      { test_gatherscatter_int_2D(); });

    ts.emplace_back(CASE("oops/test_write_read_double_1D")
      { test_write_read_double_1D(); });

    ts.emplace_back(CASE("oops/test_write_read_multiple")
      { test_write_read_multiple(); });
  }

  void clear() const override {}
};


}  // namespace test
