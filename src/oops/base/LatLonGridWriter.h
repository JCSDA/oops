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
// - add dimensions for the number of lats/lons in the regular latlon grid
// - write a subset of levels for each upper-air field (and only 1 level for surface fields)
// - write the variable in 3d (lat x lon x alt)
void writer(const atlas::FunctionSpace & fs,
            const atlas::FieldSet & fset,
            const Variables & vars,
            const std::vector<size_t> & levels,
            const std::string & filepath,
            const std::vector<double> & lats,
            const std::vector<double> & lons) {
  // Get sizes
  const size_t nlats = lats.size();
  const size_t nlons = lons.size();
  const atlas::idx_t npts = fs.size();
  const atlas::idx_t nz = levels.size();
  ASSERT(nlats * nlons == static_cast<size_t>(npts));

  // NetCDF IDs
  int ncid, retval;
  int npts_id, nz_id, nzsfc_id, lat_dimid, lon_dimid;
  int d1D_id[1], d3D_id[3];
  int levels_id, var_id[vars.size()];
  int lat_id, lon_id;


  // NetCDF file path
  std::string ncfilepath = filepath;
  ncfilepath.append(".nc");
  oops::Log::info() << "Writing file: " << ncfilepath << std::endl;

  // Create NetCDF file
  if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

  // Create dimensions
  if ((retval = nc_def_dim(ncid, "npoints", npts, &npts_id))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "latitude", nlats, &lat_dimid))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "longitude", nlons, &lon_dimid))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "nlevels", nz, &nz_id))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "nlevels_surface", 1, &nzsfc_id))) ERR(retval);

  // Define list of level indices
  d1D_id[0] = nz_id;
  if ((retval = nc_def_var(ncid, "model_levels_used", NC_INT, 1, d1D_id, &levels_id))) ERR(retval);

  // Define coordinates
  if ((retval = nc_def_var(ncid, "latitude", NC_DOUBLE, 1, &lat_dimid, &lat_id))) ERR(retval);
  if ((retval = nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_dimid, &lon_id))) ERR(retval);

  // add unit attributes
  if ((retval = nc_put_att_text(ncid, lat_id, "units",
    strlen("degrees"), "degrees")))
    ERR(retval);

  if ((retval = nc_put_att_text(ncid, lon_id, "units",
    strlen("degrees"), "degrees")))
    ERR(retval);

  // Define variables
  d3D_id[0] = lat_id;
  d3D_id[1] = lon_id;

  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    const size_t var_levs = fset.field(vars[jvar]).levels();
    if (var_levs > 1) {
      // Upper-air data case
      // Sanity check we didn't request more levels than the model has
      ASSERT(static_cast<size_t>(nz) <= var_levs);
      d3D_id[2] = nz_id;
    } else {
      // Surface data case
      d3D_id[2] = nzsfc_id;
    }

    if ((retval = nc_def_var(ncid, vars[jvar].c_str(), NC_DOUBLE, 3, d3D_id, &var_id[jvar]))) {
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

  // Write coordinates
  if ((retval = nc_put_var_double(ncid, lat_id, lats.data()))) ERR(retval);
  if ((retval = nc_put_var_double(ncid, lon_id, lons.data()))) ERR(retval);

  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    const size_t var_levs = fset.field(vars[jvar]).levels();
    auto varView = atlas::array::make_view<double, 2>(fset[vars[jvar]]);
    if (var_levs > 1) {
      double zvar2d[nlats][nlons][nz];
      for (atlas::idx_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < nlats; ++j) {
          for (size_t i = 0; i < nlons; ++i) {
            zvar2d[j][i][k] = varView(j*nlons+i, levels[k]);
          }
        }
      }
      if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar2d[0][0][0]))) ERR(retval);
    } else {
      double zvar2d[nlats][nlons][nz];
      for (size_t j = 0; j < nlats; ++j) {
        for (size_t i = 0; i < nlons; ++i) {
          zvar2d[j][i][0] = varView(j*nlons+i, 0);
        }
      }
      if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar2d[0][0][0]))) ERR(retval);
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
  RequiredParameter<std::vector<std::size_t>> levels{"levels to output", this};
  // output file will be: <datapath>/<filename prefix>.YYYYMMDDThhmmssZ.nc
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
  std::unique_ptr<oops::GlobalInterpolator> interp_;

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
    const std::string filepath = path_ + "/" + prefix_ + "." + xx.validTime().toStringIO();;

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

    detail::writer(*targetFunctionSpace_, fsetLatLon, vars_, levels_, filepath, lats, lons);
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
