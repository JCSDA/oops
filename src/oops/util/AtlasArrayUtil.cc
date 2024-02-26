/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "oops/util/AtlasArrayUtil.h"

#include <netcdf.h>
#include <string>

#include "eckit/mpi/Comm.h"

#include "atlas/array.h"
#include "atlas/array/LocalView.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/missingValues.h"


namespace util {

void atlasArrayWriteHeader(
    const std::string & ncfilepath,
    const std::vector<std::string> & dimNames,
    const std::vector<atlas::idx_t> & dimSizes,
    const std::vector<std::string> & variableNames,
    const std::vector<std::vector<std::string>> & dimNamesForEveryVar,
    std::vector<int> & netcdfGeneralIDs,
    std::vector<int> & netcdfDimIDs,
    std::vector<int> & netcdfVarIDs,
    std::vector<std::vector<int>> & netcdfDimVarIDs) {
  netcdfGeneralIDs.clear();
  netcdfDimIDs.clear();
  netcdfVarIDs.clear();
  netcdfDimVarIDs.clear();

  int ncid;
  int retval;

  // Missing value
  const double msv(::util::missingValue<double>());

  if ((retval = nc_create(ncfilepath.c_str(), NC_NETCDF4, &ncid))) ERR1(retval);
  netcdfGeneralIDs.push_back(ncid);

  for (std::size_t i = 0; i < dimNames.size(); ++i) {
    int dim_id;
    if ((retval = nc_def_dim(netcdfGeneralIDs[0],
                             dimNames[i].c_str(),
                             dimSizes[i], &dim_id))) ERR1(retval);
    netcdfDimIDs.push_back(dim_id);
  }

  // dim names and dim index are in the same order -> maybe use pair?
  std::vector<int> dimVarIDs;

  for (std::size_t i = 0; i < variableNames.size(); ++i) {
    std::vector<std::string> dimNamesForEachVar = dimNamesForEveryVar[i];

    dimVarIDs.clear();
    // find indices on dimNames of each
    for (std::size_t d2 = 0; d2 < dimNamesForEachVar.size(); ++d2) {
      std::size_t j(dimNamesForEachVar.size());
      for (std::size_t d = 0; d < dimNames.size(); ++d) {
        if (dimNames[d].compare(dimNamesForEachVar[d2]) == 0) {
          j = d;
        }
      }
      dimVarIDs.push_back(netcdfDimIDs.at(j));
    }
    netcdfDimVarIDs.push_back(dimVarIDs);

    int varID;
    if ((retval = nc_def_var(netcdfGeneralIDs[0],
                             variableNames[i].c_str(),
                             NC_DOUBLE,
                             dimVarIDs.size(),
                             dimVarIDs.data(),
                             &varID))) ERR1(retval);

    if ((retval = nc_put_att_double(netcdfGeneralIDs[0],
                                    varID, "_FillValue",
                                    NC_DOUBLE, 1,
                                    &msv))) ERR1(retval);
    netcdfVarIDs.push_back(varID);
  }

  // End definition mode
  if ((retval = nc_enddef(ncid))) ERR1(retval);
}

// Needs header info about all the fields in the file.
void atlasArrayInquire(
    const std::string & ncfilepath,
    std::vector<std::string> & dimNames,
    std::vector<atlas::idx_t> & dimSizes,
    std::vector<std::string> & variableNames,
    std::vector<std::vector<std::string>> & dimNamesForEveryVar,
    std::vector<int> & netcdfGeneralIDs,
    std::vector<int> & netcdfDimIDs,
    std::vector<int> & netcdfVarIDs,
    std::vector<std::vector<int>> & netcdfDimVarIDs) {
  dimNames.clear();
  dimSizes.clear();
  variableNames.clear();
  dimNamesForEveryVar.clear();
  netcdfGeneralIDs.clear();
  netcdfDimIDs.clear();
  netcdfVarIDs.clear();
  netcdfDimVarIDs.clear();

  int ncid, retval, ndims, nvars, ngatts, unlimdimid;

  // Open NetCDF file
  if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR1(retval);
  netcdfGeneralIDs.push_back(ncid);

  // Look at dataset
  if ((retval = nc_inq(netcdfGeneralIDs[0], &ndims, &nvars, &ngatts, &unlimdimid))) ERR1(retval);

  for (int dimid = 0; dimid < ndims; ++dimid) {
    char cdim_name[NC_MAX_NAME+1];
    size_t lenp;
    if ((retval = nc_inq_dim(netcdfGeneralIDs[0], dimid, cdim_name, &lenp))) ERR1(retval);
    dimNames.push_back(std::string(cdim_name));
    dimSizes.push_back(atlas::idx_t(lenp));
    netcdfDimIDs.push_back(dimid);
  }

  for (int varid = 0; varid < nvars; ++varid) {
    char cvar_name[NC_MAX_NAME+1];
    nc_type nc_var_type;
    int var_ndims;
    int var_dimids[NC_MAX_VAR_DIMS];
    int var_natts;
    if ((retval = nc_inq_var(netcdfGeneralIDs[0], varid, cvar_name, &nc_var_type,
                             &var_ndims, var_dimids, &var_natts))) ERR1(retval);
    variableNames.push_back(std::string(cvar_name));
    netcdfVarIDs.push_back(varid);

    std::vector<int> netcdf_dim_varID;
    std::vector<std::string> dimNamesForEachVar;
    for (int var_dim = 0; var_dim < var_ndims; ++var_dim) {
      netcdf_dim_varID.push_back(var_dimids[var_dim]);
      char cdim_name[NC_MAX_NAME+1];
      size_t lenp;
      if ((retval = nc_inq_dim(netcdfGeneralIDs[0],
                               var_dimids[var_dim],
                               cdim_name,
                               &lenp))) ERR1(retval);
      dimNamesForEachVar.push_back(std::string(cdim_name));
    }
    netcdfDimVarIDs.push_back(netcdf_dim_varID);
    dimNamesForEveryVar.push_back(dimNamesForEachVar);
  }
}


void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 1> & arrayInOut) {
  int retval;
  double zvar[dimSizes[0]];

  if ((retval = nc_get_var_double(netcdfGeneralIDs[0],
                                  varID, &zvar[0]))) ERR1(retval);

  for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
    arrayInOut(t) = zvar[t];
  }
}


void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 2> & arrayInOut) {
  int retval;
  double zvar[dimSizes[0]][dimSizes[1]];

  if ((retval = nc_get_var_double(netcdfGeneralIDs[0],
                                  varID, &zvar[0][0]))) ERR1(retval);

  for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
    for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizes[1]); ++t1) {
      arrayInOut(t, t1) = zvar[t][t1];
    }
  }
}

void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 3> & arrayInOut) {
  int retval;
  double zvar[dimSizes[0]][dimSizes[1]][dimSizes[2]];

  if ((retval = nc_get_var_double(netcdfGeneralIDs[0],
                                  varID, &zvar[0][0][0]))) ERR1(retval);

  for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
    for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizes[1]); ++t1) {
      for (idx_t t2 = 0; t2 < static_cast<idx_t>(dimSizes[2]); ++t2) {
        arrayInOut(t, t1, t2) = zvar[t][t1][t2];
      }
    }
  }
}


void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<int, 1> & arrayInOut) {
  int retval;
  int zvar[dimSizes[0]];

  if ((retval = nc_get_var_int(netcdfGeneralIDs[0],
                               varID, &zvar[0]))) ERR1(retval);

  for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
    arrayInOut(t) = zvar[t];
  }
}


void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 1> & arrayIn) {
  int retval;
  double zvar[arrayIn.shape()[0]];
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
    zvar[t] = arrayIn(t);
  }

  if ((retval = nc_put_var_double(netcdfGeneralIDs[0], varID, &zvar[0]))) ERR1(retval);
}


void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 2> & arrayIn) {
  int retval;
  double zvar[arrayIn.shape()[0]][arrayIn.shape()[1]];
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayIn.shape()[1]; ++t1) {
      zvar[t][t1] = arrayIn(t, t1);
    }
  }

  if ((retval = nc_put_var_double(netcdfGeneralIDs[0], varID, &zvar[0][0]))) ERR1(retval);
}

void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 3> & arrayIn) {
  int retval;
  double zvar[arrayIn.shape()[0]][arrayIn.shape()[1]][arrayIn.shape()[2]];
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayIn.shape()[1]; ++t1) {
      for (idx_t t2 = 0; t2 < arrayIn.shape()[2]; ++t2) {
        zvar[t][t1][t2] = arrayIn(t, t1, t2);
      }
    }
  }

  if ((retval = nc_put_var_double(netcdfGeneralIDs[0], varID, &zvar[0][0][0]))) ERR1(retval);
}

void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const int, 1> & arrayIn) {
  int retval;
  int zvar[arrayIn.shape()[0]];
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
      zvar[t] = arrayIn(t);
  }

  if ((retval = nc_put_var_int(netcdfGeneralIDs[0], varID, &zvar[0]))) ERR1(retval);
}


}  // namespace util
