/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <netcdf.h>

#include <memory>
#include <string>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/generic/GlobalInterpolator.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

// From SABER quench/Fields.cc file
#define ERR(e) {ABORT(nc_strerror(e));}

// -----------------------------------------------------------------------------

namespace oops {
namespace detail {

// Originally from SABER quench/Fields.cc file (PointCloud case), with minor changes to...
// - add dimensions for the number of lats/lons in the regulat latlon grid
// - write a subset of levels for each upper-air field (and only 1 level for surface fields)
void writer(const atlas::FunctionSpace & fs,
            const atlas::FieldSet & fset,
            const Variables & vars,
            const std::vector<size_t> & levels,
            const std::string & filepath,
            const size_t nlats,
            const size_t nlons) {
  // Get sizes
  atlas::idx_t npts = fs.size();
  atlas::idx_t nz = levels.size();
  ASSERT(nlats * nlons == static_cast<size_t>(npts));

  // NetCDF IDs
  int ncid, retval;
  int npts_id, nlats_id, nlons_id, nz_id, nzsfc_id;
  int d1D_id[1], d2D_id[2];
  int levels_id, lon_id, lat_id, var_id[vars.size()];

  // NetCDF file path
  std::string ncfilepath = filepath;
  ncfilepath.append(".nc");
  oops::Log::info() << "Writing file: " << ncfilepath << std::endl;

  // Create NetCDF file
  if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

  // Create dimensions
  if ((retval = nc_def_dim(ncid, "npoints", npts, &npts_id))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "nlats", nlats, &nlats_id))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "nlons", nlons, &nlons_id))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "nlevels", nz, &nz_id))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "nlevels_surface", 1, &nzsfc_id))) ERR(retval);

  // Define list of level indices
  d1D_id[0] = nz_id;
  if ((retval = nc_def_var(ncid, "model_levels_used", NC_INT, 1, d1D_id, &levels_id))) ERR(retval);

  // Define coordinates
  d1D_id[0] = npts_id;
  if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 1, d1D_id, &lon_id))) ERR(retval);
  if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 1, d1D_id, &lat_id))) ERR(retval);

  // Define variables
  d2D_id[0] = npts_id;
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    const size_t var_levs = fset.field(vars[jvar]).levels();
    if (var_levs > 1) {
      // Upper-air data case
      // Sanity check we didn't request more levels than the model has
      ASSERT(static_cast<size_t>(nz) <= var_levs);
      d2D_id[1] = nz_id;
    } else {
      // Surface data case
      d2D_id[1] = nzsfc_id;
    }
    if ((retval = nc_def_var(ncid, vars[jvar].c_str(), NC_DOUBLE, 2, d2D_id, &var_id[jvar]))) {
      ERR(retval);
    }
  }

  // End definition mode
  if ((retval = nc_enddef(ncid))) ERR(retval);

  // Copy levels
  int zlevels[nz][1];
  for (atlas::idx_t i = 0; i < nz; ++i) {
    zlevels[i][0] = static_cast<int>(levels[i]);
  }
  // Write levels
  if ((retval = nc_put_var_int(ncid, levels_id, &zlevels[0][0]))) ERR(retval);

  // Copy coordinates
  double zlon[npts][1];
  double zlat[npts][1];
  auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
  for (atlas::idx_t i = 0; i < npts; ++i) {
    zlon[i][0] = lonlatView(i, 0);
    zlat[i][0] = lonlatView(i, 1);
  }

  // Write coordinates
  if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
  if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    const size_t var_levs = fset.field(vars[jvar]).levels();
    auto varView = atlas::array::make_view<double, 2>(fset[vars[jvar]]);
    if (var_levs > 1) {
      // Upper-air data case
      // Copy data
      double zvar[npts][nz];
      for (atlas::idx_t k = 0; k < nz; ++k) {
        for (atlas::idx_t i = 0; i < npts; ++i) {
          zvar[i][k] = varView(i, levels[k]);
        }
      }
      // Write data
      if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);
    } else {
      // Surface data case
      // Copy data
      double zvar[npts][1];
      for (atlas::idx_t i = 0; i < npts; ++i) {
        zvar[i][0] = varView(i, 0);
      }
      // Write data
      if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);
    }
  }

  // Close file
  if ((retval = nc_close(ncid))) ERR(retval);
}

}  // namespace detail

// -----------------------------------------------------------------------------

class LatLonGridWriterParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(LatLonGridWriterParameters, Parameters)

 public:
  // Grid spacing in degrees. The output grid will be an atlas "LN" latlon grid with (2N+1) latitude
  // (incl both poles) and 4N longitude samples. N will be determined so that the resolution at the
  // equator is equal to, or slightly finer than, the requested resolution.
  RequiredParameter<double> gridRes{"resolution in degrees", this,
                                    {oops::exclusiveMinConstraint(0.0), oops::maxConstraint(90.0)}};
  RequiredParameter<Variables> variables{"variables to output", this};
  RequiredParameter<std::vector<std::size_t>> levels{"levels to output", this};
  // output file will be: <datapath>/<filename prefix>.YYYY-MM-DDThh:mm:ssZ.nc
  // as with other writers, assumes datapath exists
  Parameter<std::string> path{"datapath", ".", this};
  RequiredParameter<std::string> prefix{"filename prefix", this};
};

// -----------------------------------------------------------------------------

/// Interpolates requested model fields to a uniform lat-lon grid, then writes to a NetCDF file.
/// This is useful for plotting model fields and for diagnostic computations.
template <typename MODEL>
class LatLonGridWriter : public util::Printable  {
 public:
  explicit LatLonGridWriter(
      const LatLonGridWriterParameters & parameters,
      const Geometry<MODEL> & sourceGeometry);
  ~LatLonGridWriter() = default;

  template <typename FLDS>
  void interpolateAndWrite(const FLDS & xx) const;

 private:
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> targetFunctionSpace_;
  std::unique_ptr<oops::GlobalInterpolator<MODEL>> interp_;

  size_t gridRes_;
  Variables vars_;
  std::vector<std::size_t> levels_;
  std::string path_;
  std::string prefix_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LatLonGridWriter<MODEL>::LatLonGridWriter(
    const LatLonGridWriterParameters & parameters,
    const Geometry<MODEL> & sourceGeometry
    )
: comm_(sourceGeometry.getComm()),
  vars_(parameters.variables),
  levels_(parameters.levels),
  path_(parameters.path),
  prefix_(parameters.prefix)
{
  // Sanity check levels:
  ASSERT(levels_.size() > 0);

  // Goal: find number of grid points on a 90deg arc giving the requested resolution or higher.
  //       This number N will be construct an atlas LN latlon grid with (2N+1) lats and 4N lons.
  const double ratio = 90.0 / parameters.gridRes;
  const auto isNearlyInt = [](const double x) { return fabs(round(x) - x) < 1.e-6; };
  gridRes_ = static_cast<size_t>(isNearlyInt(ratio) ? round(ratio) : ceil(ratio));

  // Want final results on proc0 only, so only proc0 gets points in the PointCloud
  if (comm_.rank() == 0) {
    const std::string atlasGridName = "L" + std::to_string(gridRes_);
    const atlas::Grid grid(atlasGridName);
    targetFunctionSpace_.reset(new atlas::functionspace::PointCloud(grid));
  } else {
    const std::vector<atlas::PointXY> empty{{}};
    targetFunctionSpace_.reset(new atlas::functionspace::PointCloud(empty));
  }

  const eckit::LocalConfiguration emptyConf{};
  interp_.reset(new oops::GlobalInterpolator<MODEL>(emptyConf,
        sourceGeometry, *targetFunctionSpace_, sourceGeometry.getComm()));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
template <typename FLDS>
void LatLonGridWriter<MODEL>::interpolateAndWrite(const FLDS & xx) const {
  const auto & fset = xx.fieldSet();

  // Currently, must initialize target FieldSet before passing to GlobalInterp
  atlas::FieldSet fsetLatLon;
  for (const auto & f : fset) {
    const std::string name = f.name();
    const size_t nlev = f.levels();
    const atlas::Field field = targetFunctionSpace_->createField<double>(
        atlas::option::name(name) | atlas::option::levels(nlev));
    fsetLatLon.add(field);
  }

  // Interpolate input data to lat-lon grid
  interp_->apply(fset, fsetLatLon);

  // Write data on proc0 only
  if (comm_.rank() == 0) {
    const std::string filepath = path_ + "/" + prefix_ + "." + xx.validTime().toString();
    const size_t nlats = 2*gridRes_ + 1;
    const size_t nlons = 4*gridRes_;
    detail::writer(*targetFunctionSpace_, fsetLatLon, vars_, levels_, filepath, nlats, nlons);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LatLonGridWriter<MODEL>::print(std::ostream & os) const
{
  os << "LatLonGridWriter::<" << MODEL::name() << ">" << std::endl;
  os << "Grid Resolution: " << gridRes_ << std::endl;
  os << "Variables: " << vars_ << std::endl;
  os << "Levels: ";
  for (size_t i=0; i < levels_.size(); i++) {
    os << levels_[i] << " ";
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
