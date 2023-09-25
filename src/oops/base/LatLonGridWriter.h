/*
 * (C) Copyright 2022-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <netcdf.h>

#include <memory>
#include <optional>  // NOLINT(build/include_order): linter mis-identifies C++ header as C
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/GlobalInterpolator.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

// -----------------------------------------------------------------------------

namespace oops {
namespace detail {

static inline void check_nc_code(const int return_code) {
  if (return_code) {
    ABORT(nc_strerror(return_code));
  }
}

void writerForPressures(const atlas::FieldSet & fset,
                        const Variables & vars,
                        const std::vector<double> & pressureLevels,
                        const std::string & filepath,
                        const std::vector<double> & lats,
                        const std::vector<double> & lons) {
  // Get sizes
  const size_t nlats = lats.size();
  const size_t nlons = lons.size();
  const size_t nz = pressureLevels.size();

  // NetCDF IDs
  int ncid;  // file ID
  int lat_did, lon_did, lev_did;  // dim IDs
  int lat_vid, lon_vid, lev_vid, field_vid[vars.size()];  // var IDs

  // NetCDF file path
  std::string ncfilepath = filepath;
  ncfilepath.append(".pressureLevels.nc");
  oops::Log::info() << "Writing file: " << ncfilepath << std::endl;

  // Create NetCDF file
  check_nc_code(nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid));
  // Warning: the output here does not in general conform to the CF conventions, but we fake it
  //          via this netcdf attribute to allow METplus to handle certain regridding operations.
  check_nc_code(nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen("CF-1"), "CF-1"));

  // Create dimensions
  check_nc_code(nc_def_dim(ncid, "latitude", nlats, &lat_did));
  check_nc_code(nc_def_dim(ncid, "longitude", nlons, &lon_did));
  check_nc_code(nc_def_dim(ncid, "levels", nz, &lev_did));

  // Define coordinates
  check_nc_code(nc_def_var(ncid, "latitude", NC_DOUBLE, 1, &lat_did, &lat_vid));
  check_nc_code(nc_put_att_text(ncid, lat_vid, "units", strlen("degrees_north"), "degrees_north"));
  check_nc_code(nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_did, &lon_vid));
  check_nc_code(nc_put_att_text(ncid, lon_vid, "units", strlen("degrees_east"), "degrees_east"));
  check_nc_code(nc_def_var(ncid, "pressure_levels", NC_DOUBLE, 1, &lev_did, &lev_vid));
  check_nc_code(nc_put_att_text(ncid, lev_vid, "units", strlen("hPa"), "hPa"));

  // Define variables
  int field_dids[3] = {lat_did, lon_did, lev_did};
  const double missing = util::missingValue<double>();
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    check_nc_code(nc_def_var(ncid, vars[jvar].c_str(), NC_DOUBLE, 3, field_dids, &field_vid[jvar]));
    check_nc_code(
        nc_put_att_double(ncid, field_vid[jvar], "missing_value", NC_DOUBLE, 1, &missing));
  }

  // End definition mode
  check_nc_code(nc_enddef(ncid));

  // Write coordinates
  check_nc_code(nc_put_var_double(ncid, lat_vid, lats.data()));
  check_nc_code(nc_put_var_double(ncid, lon_vid, lons.data()));
  check_nc_code(nc_put_var_double(ncid, lev_vid, pressureLevels.data()));

  // Write fields
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    ASSERT(nz == static_cast<size_t>(fset.field(vars[jvar]).levels()));
    auto varView = atlas::array::make_view<double, 2>(fset[vars[jvar]]);
    double values[nlats][nlons][nz];
    for (size_t k = 0; k < nz; ++k) {
      for (size_t j = 0; j < nlats; ++j) {
        for (size_t i = 0; i < nlons; ++i) {
          values[j][i][k] = varView(j*nlons+i, k);
        }
      }
    }
    check_nc_code(nc_put_var_double(ncid, field_vid[jvar], &values[0][0][0]));
  }

  // Close file
  check_nc_code(nc_close(ncid));
}

void writerForLevels(const atlas::FieldSet & fset,
                     const Variables & vars,
                     const std::optional<std::vector<size_t>> & modelLevels,
                     const std::unordered_map<std::string, bool> isSurfaceVar,
                     const std::string & filepath,
                     const std::vector<double> & lats,
                     const std::vector<double> & lons) {
  // Check if FieldSet has surface or upper-air fields.
  // Note we can't use the number of levels from the FieldSet's Fields, because those have been
  // trimmed to the output levels, which could be 1 level even for an upper-air field.
  const bool haveSurfaceFields = std::any_of(isSurfaceVar.begin(), isSurfaceVar.end(),
      [&](const auto & kv) {return kv.second && fset.has(kv.first);});
  const bool haveUpperAirFields = std::any_of(isSurfaceVar.begin(), isSurfaceVar.end(),
      [&](const auto & kv) {return !kv.second && fset.has(kv.first);});

  // Check levels are provided if there are upper-air fields
  if (haveUpperAirFields) {
    ASSERT(modelLevels.has_value());
  }

  // Get sizes
  const size_t nlats = lats.size();
  const size_t nlons = lons.size();
  const size_t nz = (haveUpperAirFields ? modelLevels->size() : 1);

  // More checks on levels
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    const size_t var_levs = fset.field(vars[jvar]).levels();
    if (!isSurfaceVar.at(vars[jvar])) {
      ASSERT(var_levs == static_cast<size_t>(nz));
    } else {
      ASSERT(var_levs == 1);
    }
  }

  // NetCDF IDs
  int ncid;  // file ID
  int lat_did, lon_did, lev_did, sfc_did;  // dim IDs
  int lat_vid, lon_vid, lev_vid, field_vid[vars.size()];  // var IDs

  // NetCDF file path
  std::string ncfilepath = filepath;
  ncfilepath.append(".modelLevels.nc");
  oops::Log::info() << "Writing file: " << ncfilepath << std::endl;

  // Create NetCDF file
  check_nc_code(nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid));
  // Warning: the output here does not in general conform to the CF conventions, but we fake it
  //          via this netcdf attribute to allow METplus to handle certain regridding operations.
  check_nc_code(nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen("CF-1"), "CF-1"));

  // Create dimensions
  check_nc_code(nc_def_dim(ncid, "latitude", nlats, &lat_did));
  check_nc_code(nc_def_dim(ncid, "longitude", nlons, &lon_did));
  if (haveUpperAirFields) {
    check_nc_code(nc_def_dim(ncid, "levels", nz, &lev_did));
  }
  if (haveSurfaceFields) {
    check_nc_code(nc_def_dim(ncid, "dummy_surface_levels", 1, &sfc_did));
  }

  // Define coordinates
  check_nc_code(nc_def_var(ncid, "latitude", NC_DOUBLE, 1, &lat_did, &lat_vid));
  check_nc_code(nc_put_att_text(ncid, lat_vid, "units", strlen("degrees_north"), "degrees_north"));
  check_nc_code(nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_did, &lon_vid));
  check_nc_code(nc_put_att_text(ncid, lon_vid, "units", strlen("degrees_east"), "degrees_east"));
  if (haveUpperAirFields) {
    check_nc_code(nc_def_var(ncid, "model_levels", NC_INT, 1, &lev_did, &lev_vid));
  }

  // Define variables
  int field_dids[3] = {lat_did, lon_did, -1};  // note 3rd element is initialized below
  const double missing = util::missingValue<double>();
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    field_dids[2] = (isSurfaceVar.at(vars[jvar]) ? sfc_did : lev_did);
    check_nc_code(nc_def_var(ncid, vars[jvar].c_str(), NC_DOUBLE, 3, field_dids, &field_vid[jvar]));
    check_nc_code(
        nc_put_att_double(ncid, field_vid[jvar], "missing_value", NC_DOUBLE, 1, &missing));
  }

  // End definition mode
  check_nc_code(nc_enddef(ncid));

  // Write coordinates
  check_nc_code(nc_put_var_double(ncid, lat_vid, lats.data()));
  check_nc_code(nc_put_var_double(ncid, lon_vid, lons.data()));
  if (haveUpperAirFields) {
    // Convert size_t to int
    std::vector<int> values(modelLevels->size());
    for (size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<int>((*modelLevels)[i]);
    }
    check_nc_code(nc_put_var_int(ncid, lev_vid, values.data()));
  }

  // Write fields
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    const size_t var_levs = fset.field(vars[jvar]).levels();
    auto varView = atlas::array::make_view<double, 2>(fset[vars[jvar]]);
    double values[nlats][nlons][var_levs];
    for (size_t k = 0; k < var_levs; ++k) {
      for (size_t j = 0; j < nlats; ++j) {
        for (size_t i = 0; i < nlons; ++i) {
          values[j][i][k] = varView(j*nlons+i, k);
        }
      }
    }
    check_nc_code(nc_put_var_double(ncid, field_vid[jvar], &values[0][0][0]));
  }

  // Close file
  check_nc_code(nc_close(ncid));
}

std::vector<size_t> mergeLevels(const std::optional<std::vector<size_t>> modelLevels,
                                const bool bottomLevel, const bool modelLevelsAreTopDown,
                                const size_t nModelLevels) {
  const size_t levelBottomIndex = (modelLevelsAreTopDown ? nModelLevels - 1 : 0);
  std::vector<size_t> result;
  if (!modelLevels) {
    // If no model levels, just have bottom level
    result.push_back(levelBottomIndex);
  } else {
    result = *modelLevels;
    if (bottomLevel) {
      if (modelLevelsAreTopDown) {
        ASSERT(result.front() <= levelBottomIndex);
        if (result.front() < levelBottomIndex) {
          result.insert(result.begin(), levelBottomIndex);
        }
      } else {
        ASSERT(result.front() >= levelBottomIndex);
        if (result.front() > levelBottomIndex) {
          result.insert(result.begin(), levelBottomIndex);
        }
      }
    }
  }
  return result;
}

}  // namespace detail

// -----------------------------------------------------------------------------

class LatLonGridWriterParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(LatLonGridWriterParameters, Parameters)

 public:
  // Local interpolator type: currently has to be either "atlas interpolator" or
  // "oops unstructured grid interpolator". (see GlobalInterpolator for supported
  // interpolators)
  RequiredParameter<std::string> interpType{"local interpolator type", this};
  // Grid spacing in degrees. The output grid will be an atlas "LN" latlon grid with (2N+1) latitude
  // (incl both poles) and 4N longitude samples. N will be determined so that the resolution at the
  // equator is equal to, or slightly finer than, the requested resolution.
  RequiredParameter<double> gridRes{"resolution in degrees", this,
                                    {oops::exclusiveMinConstraint(0.0), oops::maxConstraint(90.0)}};

  RequiredParameter<Variables> variables{"variables to output", this};
  OptionalParameter<std::vector<double>> pressureLevels{"pressure levels in hPa", this};
  OptionalParameter<std::vector<size_t>> modelLevels{"model levels", this};
  Parameter<bool> bottomLevel{"bottom model level", false, this};

  // output file will be: <datapath>/<filename prefix>.YYYYMMDDThhmmssZ.nc
  // as with other writers, assumes datapath exists
  Parameter<std::string> path{"datapath", ".", this};
  RequiredParameter<std::string> prefix{"filename prefix", this};
};

// -----------------------------------------------------------------------------

/// Outputs model fields interpolated to a uniform lat-lon grid, optionally also interpolated to
/// pressure levels. This is useful for plotting model fields and for diagnostic computations.
///
/// The LatLonGridWriter can output a list of requested fields onto a list of requested levels,
/// either pressure levels or model levels. Usage of levels is as follows:
///
/// - if only surface fields are requested, then there is no need to request output levels; if any
///   are specified they will be ignored.
/// - if any upper-air fields are requested, then the user must request specific pressure and/or
///   model levels for the output. This can be done via any combination of options
///   - "pressure levels in hPa" to give a list of pressure levels; levels should be ordered
///     bottom-to-top, i.e., in decreasing pressure order
///   - "model levels" to give a list of model levels; levels should be ordered bottom-to-top
///   - "bottom model level" to request output at the bottom model level. This option is partly
///     redundant with "model levels", but allows the user to request the commonly-used lowest
///     level without needing to know the details of the indexing of model levels. If this is used
///     in conjuction with "model levels", the bottom level requested by this option will be
///     prepended to the model level list from the "model levels" option. If both options are
///     requesting the bottom level, it will only be printed once into the output file.
///   If more than one of these options is requested, then the union of the levels (both pressure
///   and model levels) will be output.
/// - a mix of surface and upper-air fields can be requested.
///
/// The LatLonGridWriter applies a VariableChange or LinearVariableChange to the provided State or
/// Increment before interpolations. This is currently not user-specifiable -- the MODEL's default
/// variable change is called. In most cases this variable change can compute fields of interest.
/// Control over this variable change can be added in the future.
///
/// If output is requested on pressure levels, then the vertical interpolation to these levels is
/// done by linear interpolation in log-pressure coordinate. If output is requested on model levels,
/// the requested level numbers are read off.
///
/// Results are written to one or two NetCDF files: one file for pressure-level output, and another
/// for model-level output. Surface fields are written into the model-level file.
template <typename MODEL>
class LatLonGridWriter : public util::Printable  {
 public:
  explicit LatLonGridWriter(
      const LatLonGridWriterParameters & parameters,
      const Geometry<MODEL> & sourceGeometry);
  ~LatLonGridWriter() = default;


  /// Writes data from xx, optionally interpolated to pressure levels
  /// Performs a variable change from xx to the requested output variables
  void interpolateAndWrite(const State<MODEL> & xx) const;

  /// Writes data from dx, optionally interpolated to pressure levels from xx
  /// Performs a linear variable change from dx to the requested output variables
  void interpolateAndWrite(const Increment<MODEL> & dx, const State<MODEL> & xx) const;

  /// Simple interface for contexts with no available background; therefore, this only works for
  /// writing fields on model levels
  /// Performs a linear variable change from dx to the requested output variables
  void interpolateAndWrite(const Increment<MODEL> & xx) const;

 private:
  void interpolateAndWrite(const atlas::FieldSet & fset, const util::DateTime & t) const;
  void interpolateToPressureLevels(const atlas::FieldSet & fsetLatLon,
                                   atlas::FieldSet & fsetPressureLevels) const;
  void trimToModelLevels(const atlas::FieldSet & fsetLatLon, atlas::FieldSet & fsetModelLevels,
                         const std::vector<size_t> targetLevels,
                         const size_t nModelLevels) const;

  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> targetFunctionSpace_;
  std::unique_ptr<oops::GlobalInterpolator> interp_;

  const Geometry<MODEL> & sourceGeometry_;
  const Variables vars_;
  const std::string path_;
  const std::string prefix_;
  const bool bottomLevel_;
  const bool modelLevelsAreTopDown_;

  size_t gridRes_;
  std::optional<std::vector<double>> pressureLevels_;
  std::optional<std::vector<size_t>> modelLevels_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LatLonGridWriter<MODEL>::LatLonGridWriter(
    const LatLonGridWriterParameters & parameters,
    const Geometry<MODEL> & sourceGeometry
    )
: comm_(sourceGeometry.getComm()),
  sourceGeometry_(sourceGeometry),
  vars_(parameters.variables),
  path_(parameters.path),
  prefix_(parameters.prefix),
  bottomLevel_(parameters.bottomLevel),
  modelLevelsAreTopDown_(sourceGeometry.levelsAreTopDown())
{
  // Sanity check any provided level options. It's fine to request none, but in this case only
  // surface fields may be output; the LatLonGridWriter will check this in the interpolate method
  // and may error there.
  if (parameters.pressureLevels.value() != boost::none) {
    pressureLevels_ = parameters.pressureLevels.value().value();
    ASSERT(pressureLevels_->size() > 0);
    // Check pressures requested (in hPa) are physically plausible and in bottom-to-top order
    for (size_t lev = 0; lev < pressureLevels_->size(); ++lev) {
      const double pressureLevel = (*pressureLevels_)[lev];
      ASSERT(pressureLevel >= 0.0 && pressureLevel <= 1200.0);
      if (lev > 0) {
        const double pressureLevelPrevious = (*pressureLevels_)[lev - 1];
        ASSERT(pressureLevel < pressureLevelPrevious);
      }
    }
  }
  if (parameters.modelLevels.value() != boost::none) {
    modelLevels_ = parameters.modelLevels.value().value();
    ASSERT(modelLevels_->size() > 0);
    // Check levels requested are in bottom-to-top order
    for (size_t lev = 1; lev < modelLevels_->size(); ++lev) {
      const size_t modelLevel = (*modelLevels_)[lev];
      const size_t modelLevelPrevious = (*modelLevels_)[lev - 1];
      if (modelLevelsAreTopDown_) {
        ASSERT(modelLevel < modelLevelPrevious);
      } else {
        ASSERT(modelLevel > modelLevelPrevious);
      }
    }
  }

  // Goal: find number of grid points on a 90deg arc giving the requested resolution or higher.
  //       This number N will be construct an atlas LN latlon grid with (2N+1) lats and 4N lons.
  const double ratio = 90.0 / parameters.gridRes;
  const auto isNearlyInt = [](const double x) { return fabs(round(x) - x) < 1.e-6; };
  gridRes_ = static_cast<size_t>(isNearlyInt(ratio) ? round(ratio) : ceil(ratio));

  // Want final results on proc0 only, so only proc0 gets points in the PointCloud.
  // This requires passing an empty but valid FunctionSpaces to the interpolator, something that we
  // can only do with PointCloud. (This is a bit hacky and could probably be cleaned up.)
  if (comm_.rank() == 0) {
    const std::string atlasGridName = "L" + std::to_string(gridRes_);
    const atlas::RegularLonLatGrid grid(atlasGridName);
    targetFunctionSpace_.reset(new atlas::functionspace::PointCloud(grid));
  } else {
    const std::vector<atlas::PointXY> empty{{}};
    targetFunctionSpace_.reset(new atlas::functionspace::PointCloud(empty));
  }

  interp_.reset(new oops::GlobalInterpolator(parameters.toConfiguration(),
        sourceGeometry.generic(), *targetFunctionSpace_,
        sourceGeometry.getComm()));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::interpolateAndWrite(const State<MODEL> & xx) const {
  const eckit::LocalConfiguration conf{};  // if needed, could read this from yaml
  const oops::VariableChange<MODEL> varchange(conf, sourceGeometry_);

  oops::Variables vars_for_latlon_interp = vars_;
  if (pressureLevels_) {
    ASSERT(!vars_.has("air_pressure"));
    vars_for_latlon_interp.push_back("air_pressure");
  }

  State<MODEL> tmp_xx = xx;
  varchange.changeVar(tmp_xx, vars_for_latlon_interp);

  interpolateAndWrite(tmp_xx.fieldSet(), xx.validTime());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::interpolateAndWrite(const Increment<MODEL> & dx,
                                                  const State<MODEL> & xx) const {
  const eckit::LocalConfiguration conf{};  // if needed, could read this from yaml
  const oops::VariableChange<MODEL> varchange(conf, sourceGeometry_);
  oops::LinearVariableChange<MODEL> linvarchange(sourceGeometry_, conf);

  if (pressureLevels_) {
    ASSERT(!vars_.has("air_pressure"));
  }
  Increment<MODEL> tmp_dx = dx;
  linvarchange.changeVarTraj(xx, vars_);
  linvarchange.changeVarTL(tmp_dx, vars_);

  State<MODEL> tmp_xx = xx;
  varchange.changeVar(tmp_xx, oops::Variables(std::vector<std::string>({"air_pressure"})));

  auto fset = tmp_dx.fieldSet();
  fset.add(tmp_xx.fieldSet().field("air_pressure"));

  interpolateAndWrite(fset, dx.validTime());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::interpolateAndWrite(const Increment<MODEL> & dx) const {
  // Sanity check: can't call this interface if output on pressure levels was requested
  if (pressureLevels_) {
    throw eckit::Exception("Writing an Increment to latlon without an available background State "
                           "is not compatible with request for output on pressure levels.");
  }
  // Sanity check: no background, so can't perform variable change
  for (const auto & v : vars_.variables()) {
    if (!dx.variables().has(v)) {
      throw eckit::Exception("Writing an Increment to latlon without an available background State "
                             "is not compatible with request for variable not in Increment.");
    }
  }

  interpolateAndWrite(dx.fieldSet(), dx.validTime());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::interpolateAndWrite(const atlas::FieldSet & fsetInput,
                                                  const util::DateTime & t) const {
  // First, some prep work with variables and levels
  size_t nModelLevels = 0;
  // Data structure to keep track of original levels of each variable. This is useful because we
  // may still want to distinguish surface variables from upper-air variables even if we've
  // interpolated upper-air variables to a single level, resulting in every variable having one
  // level.
  std::unordered_map<std::string, bool> isSurfaceVar = {};
  oops::Variables upperAirVars = vars_;
  for (const auto & var : vars_.variables()) {
    ASSERT(fsetInput.has(var));
    const size_t levels = fsetInput.field(var).levels();
    // Assume any var that comes from the model with a single level is a surface variable
    isSurfaceVar.insert({var, (levels == 1)});
    if (levels > 1) {
      // Sanity-check that requested upper-air fields all have the same number of levels, because
      // otherwise neither pressure interpolation nor selection by model level will be meaningful:
      if (nModelLevels == 0) { nModelLevels = levels; }
      ASSERT(levels == nModelLevels);
    } else {
      upperAirVars -= var;
    }
  }

  const bool haveUpperAirVars = (upperAirVars.size() > 0);

  // (The sanity check below could be moved to the constructor if we did not use the incoming
  //  FieldSet data to identify variables as surface or upper-air. If oops::Variables gain metadata
  //  indicating surface vs upper-air fields, we could simplify this logic a bit.)
  // To output upper-air variables, at least one kind of level must be requested
  if (haveUpperAirVars) {
    ASSERT_MSG(pressureLevels_ || modelLevels_ || bottomLevel_,
               "LatLonGridWriter requires pressure and/or model levels to be requested when "
               "upper-air fields are specified");
  }

  // Prepare interpolation source and target FieldSets
  // - source should contain just variables requested for output
  //   (with the addition of air pressure if output is requested on pressure levels...)
  // - target should be initialized to match source
  atlas::FieldSet fset;
  atlas::FieldSet fsetLatLon;

  oops::Variables vars_for_latlon_interp = vars_;
  if (pressureLevels_ && haveUpperAirVars) {
    ASSERT(fsetInput.has("air_pressure"));
    ASSERT(static_cast<size_t>(fsetInput.field("air_pressure").levels()) == nModelLevels);
    ASSERT(!vars_.has("air_pressure"));
    vars_for_latlon_interp.push_back("air_pressure");
  }

  for (const auto & var : vars_for_latlon_interp.variables()) {
    const auto & field = fsetInput.field(var);
    fset.add(field);

    const std::string name = field.name();
    const size_t levels = field.levels();
    const atlas::Field fieldLatLon = targetFunctionSpace_->createField<double>(
        atlas::option::name(name) | atlas::option::levels(levels));
    fsetLatLon.add(fieldLatLon);
  }

  // Interpolate input data to lat-lon grid
  interp_->apply(fset, fsetLatLon);

  // Note that fsetLatLon lives on rank-0 only, so from here we can do the vertical processing
  // and the writing exclusively on rank 0.
  if (comm_.rank() != 0) return;

  // Before handling the vertical direction and writing, prepare some writer-related data
  const std::string filepath = path_ + "/" + prefix_ + "." + t.toStringIO();

  // Because we're using a PointCloud FunctionSpace (for reasons discussed above), we have to
  // reconstruct the grid here to obtain the lats/lons for writing to file. (The PointCloud is
  // so generic it doesn't rely on or give access to a grid.)
  const std::string atlasGridName = "L" + std::to_string(gridRes_);
  const atlas::RegularLonLatGrid grid(atlasGridName);
  const size_t nlats = grid.ny();
  const size_t nlons = grid.nx();
  std::vector<double> lats(nlats);
  std::vector<double> lons(nlons);
  for (size_t i = 0; i < nlats; ++i) {
    lats[i] = grid.lat(i);
  }
  for (size_t i = 0; i < nlons; ++i) {
    lons[i] = grid.lon(i);
  }

  // Handle output on pressure levels: vertical interpolation then writing
  if (pressureLevels_ && haveUpperAirVars) {
    atlas::FieldSet fsetPressureLevels;
    interpolateToPressureLevels(fsetLatLon, fsetPressureLevels);

    detail::writerForPressures(fsetPressureLevels, upperAirVars, *pressureLevels_,
                               filepath, lats, lons);
  }

  // Handle output on model levels: trimming to requested level then writing
  if (modelLevels_ || bottomLevel_) {
    // If needed, add bottom level to model levels
    const std::vector<size_t> targetLevels = detail::mergeLevels(modelLevels_, bottomLevel_,
                                                                 modelLevelsAreTopDown_,
                                                                 nModelLevels);

    atlas::FieldSet fsetModelLevels;
    trimToModelLevels(fsetLatLon, fsetModelLevels, targetLevels, nModelLevels);

    detail::writerForLevels(fsetModelLevels, vars_, targetLevels, isSurfaceVar,
                            filepath, lats, lons);
  } else {
    // Once we are here, we know that
    // - any upper-air fields were written on either pressure or model levels in branches above
    // - surface fields were written above if any model levels were requested
    // But we still need to write surface fields in case no levels were requested, or only
    // pressure levels were requested. This falls outside the two cases above.
    const bool anyVarsForSurface = std::any_of(isSurfaceVar.begin(), isSurfaceVar.end(),
                                               [](const auto & kv) {return kv.second;});
    if (anyVarsForSurface) {
      atlas::FieldSet fsetSurface;
      oops::Variables varsForSurface{};
      for (const auto & kv : isSurfaceVar) {
        if (kv.second) {
          fsetSurface.add(fsetLatLon.field(kv.first));
          varsForSurface.push_back(kv.first);
        }
      }
      ASSERT(varsForSurface.size() >= 1);

      detail::writerForLevels(fsetSurface, varsForSurface, modelLevels_, isSurfaceVar,
                              filepath, lats, lons);
    }
    // No else-branch needed: the absence of surface fields implies upper-air fields, and we
    // checked above that if upper-air fields are requested than at least some level is requested.
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::interpolateToPressureLevels(const atlas::FieldSet & fsetLatLon,
    atlas::FieldSet & fsetPressureLevels) const
{
  ASSERT(pressureLevels_);
  ASSERT(fsetPressureLevels.empty());

  const size_t nPressureLevels = pressureLevels_->size();  // number of target levels

  // For linear interpolation, stencil size = 2
  std::vector<size_t> interp_is(2 * nPressureLevels);
  std::vector<double> interp_ws(2 * nPressureLevels);
  std::vector<bool> pressure_level_above_ground(nPressureLevels);

  const auto & p = atlas::array::make_view<double, 2>(fsetLatLon.field("air_pressure"));
  const size_t lastlev = p.shape(1) - 1;  // index for last source level

  // Sanity-check ordering of levels: top-down has pressure increasing along 2nd dimension:
  ASSERT(modelLevelsAreTopDown_ == (p(0, 0) < p(0, lastlev)));

  // Vector to store initial guesses for pressure brackets
  std::vector<int> guesses(nPressureLevels, 0);

  const auto computeVertInterpMatrix = [&](const int i) {
    for (size_t lev = 0; lev < nPressureLevels; ++lev) {
      const double targetPressure = 100.0 * (*pressureLevels_)[lev];  // convert hPa to Pa

      // Check whether this pressure can be bracketed
      // - if targetPressure is higher than model max, level is below terrain => return missing
      // - if targetPressure is lower than model min, level is above model top => error
      pressure_level_above_ground[lev] = true;  // default assumption
      const double maxModelPressure = (modelLevelsAreTopDown_ ? p(i, lastlev) : p(i, 0));
      const double minModelPressure = (modelLevelsAreTopDown_ ? p(i, 0) : p(i, lastlev));
      if (targetPressure > maxModelPressure) {
        pressure_level_above_ground[lev] = false;
        continue;
      } else if (targetPressure < minModelPressure) {
        throw eckit::Exception("Requested interpolation to pressure level above model top");
      }

      // Find source levels that bracket the desired pressure
      size_t l = guesses[lev];
      while ((p(i, l) - targetPressure) * (p(i, l+1) - targetPressure) > 0.0) {
        // Not yet a bracket, so update l:
        if ((modelLevelsAreTopDown_ && (p(i, l) < targetPressure))
            || (!modelLevelsAreTopDown_ && (p(i, l) > targetPressure))) {
          ++l;
        } else {
          --l;
        }
      }
      ASSERT((p(i, l) - targetPressure) * (p(i, l+1) - targetPressure) <= 0.0);
      // Save guess for next horizontal grid point
      guesses[lev] = l;

      // Solve for interp weights
      interp_is[2*lev] = l;
      interp_is[2*lev + 1] = l+1;

      const double logp0 = std::log(p(i, l));
      const double logp1 = std::log(p(i, l+1));
      const double logp = std::log(targetPressure);
      interp_ws[2*lev] = (logp1 - logp) / (logp1 - logp0);
      interp_ws[2*lev + 1] = (logp - logp0) / (logp1 - logp0);
    }
  };

  // Allocate multi-level (upper-air) fields in target FieldSet for vertical interpolation
  // (Do nothing for single-level (surface) fields, those are handled in model levels
  for (auto & f : fsetLatLon) {
    if (f.levels() > 1) {
      const std::string & name = f.name();
      const atlas::Field field = targetFunctionSpace_->createField<double>(
          atlas::option::name(name) | atlas::option::levels(nPressureLevels));
      fsetPressureLevels.add(field);
    }
  }

  // Outer loop over atmospheric columns, because can reuse the same vertical interpolation
  // matrix for every field in the column
  for (int i = 0; i < p.shape(0); ++i) {
    computeVertInterpMatrix(i);

    // Apply matrix to each multi-level field
    for (auto & f : fsetLatLon) {
      if (f.levels() > 1) {
        // Interpolate data from full-levels field to new subset-of-levels field
        const auto fullview = atlas::array::make_view<double, 2>(f);
        auto view = atlas::array::make_view<double, 2>(fsetPressureLevels.field(f.name()));
        for (size_t lev = 0; lev < nPressureLevels; ++lev) {
          if (pressure_level_above_ground[lev]) {
            view(i, lev) = interp_ws[2*lev] * fullview(i, interp_is[2*lev])
              + interp_ws[2*lev + 1] * fullview(i, interp_is[2*lev + 1]);
          } else {
            // Pressure level below model terrain, return missing
            view(i, lev) = util::missingValue<double>();
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::trimToModelLevels(const atlas::FieldSet & fsetLatLon,
                                                atlas::FieldSet & fsetModelLevels,
                                                const std::vector<size_t> targetLevels,
                                                const size_t nModelLevels) const
{
  const size_t nTargetLevels = targetLevels.size();

  for (auto & f : fsetLatLon) {
    // For single-level (typically = surface) fields, copy directly
    if (f.levels() == 1) {
      fsetModelLevels.add(f);
    } else {
      // Make new Field with fewer levels
      const std::string name = f.name();
      atlas::Field field = targetFunctionSpace_->createField<double>(
          atlas::option::name(name) | atlas::option::levels(nTargetLevels));

      // Copy data from full-levels field to new subset-of-levels field
      const auto fullview = atlas::array::make_view<double, 2>(f);
      auto view = atlas::array::make_view<double, 2>(field);
      for (size_t lev = 0; lev < nTargetLevels; ++lev) {
        ASSERT(targetLevels[lev] >= 0 && targetLevels[lev] < static_cast<size_t>(f.shape(1)));
        for (int i = 0; i < f.shape(0); ++i) {
          view(i, lev) = fullview(i, targetLevels[lev]);
        }
      }
      fsetModelLevels.add(field);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::print(std::ostream & os) const
{
  os << "LatLonGridWriter::<" << MODEL::name() << ">" << std::endl;
  os << "Grid Resolution: " << gridRes_ << std::endl;
  os << "Variables: " << vars_ << std::endl;
  os << "Output on Levels:";
  if (pressureLevels_) {
    os << " Pressure Levels (hPa) = (";
    for (size_t lev = 0; lev < pressureLevels_->size(); ++lev) {
      os << (*pressureLevels_)[lev] << ",";
    }
    os << ")";
  }
  if (modelLevels_) {
    if (pressureLevels_) os << " and";
    os << " Model Levels = (";
    for (size_t lev = 0; lev < modelLevels_->size(); ++lev) {
      os << (*modelLevels_)[lev] << ",";
    }
    os << ")";
  }
  if (bottomLevel_) {
    if (pressureLevels_ || modelLevels_) os << " and";
    os << " Bottom Model Level";
  }
  os << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
