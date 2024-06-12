/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "oops/util/AtlasArrayUtil.h"

#include <netcdf.h>
#include <string>


#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "atlas/array.h"
#include "atlas/array/LocalView.h"

#include "oops/base/Variable.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"


namespace util {

std::vector<std::string> getAttributeNames(const eckit::LocalConfiguration & conf,
                                           const std::string & varname) {
  std::vector<eckit::LocalConfiguration> aconfs = conf.has(varname) ?
    conf.getSubConfigurations(varname) : std::vector<eckit::LocalConfiguration>();

  std::vector<std::string> names;

  for (const auto & aconf : aconfs) {
    eckit::LocalConfiguration avconf = aconf.getSubConfiguration(aconf.keys()[0]);
    if (avconf.keys().size() < 2) {
      throw eckit::BadValue("we expect each 'attribute' to have at"
                            " least a 'value' and a 'type'. ", Here());
    }
    names.push_back(aconf.keys()[0]);
  }

  return names;
}

std::string getAttributeType(const eckit::LocalConfiguration & conf,
                             const std::string & varname,
                             const std::string & attname) {
  std::vector<eckit::LocalConfiguration> aconfs = conf.has(varname) ?
    conf.getSubConfigurations(varname) : std::vector<eckit::LocalConfiguration>();

  eckit::LocalConfiguration atconf;
  std::string typestr("");
  for (std::size_t i = 0; i < aconfs.size(); ++i) {
    if (aconfs[i].has(attname)) {
      aconfs[i].get(attname, atconf);
      return atconf.getString("type");
    }
  }
  return typestr;
}

// we use variable.datatype() to define the intent of what datatype
// needs to be stored in the NetCDF file and what missing value indicator
// should be used.
// for now we will support real64 (default) and int32
// the metadata for attributes come from the LocalConfiguration object.
// TO DO (Marek) - to refactor all this code to work from Field / FieldSets
void atlasArrayWriteHeader(
    const std::string & ncfilepath,
    const std::vector<std::string> & dimNames,
    const std::vector<atlas::idx_t> & dimSizes,
    const oops::Variables & variables,
    const std::vector<std::vector<std::string>> & dimNamesForEveryVar,
    const eckit::LocalConfiguration & conf,
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
  const double msvd(::util::missingValue<double>());
  const int32_t msvi(::util::missingValue<int32_t>());

  if ((retval = nc_create(ncfilepath.c_str(), NC_NETCDF4, &ncid))) ERR1(retval);

  const std::vector<std::string> attnames =
    util::getAttributeNames(conf, "global metadata");
  netcdfGeneralIDs.push_back(ncid);

  for (const std::string & attname : attnames) {
    const std::string type =
      util::getAttributeType(conf, "global metadata", attname);
    if (type == "string") {
      const std::string val =
        util::getAttributeValue<std::string>(conf, "global metadata", attname);
      if ((retval = nc_put_att_text(netcdfGeneralIDs[0], NC_GLOBAL, attname.data(),
                                    val.length() + 1, val.data() ))) ERR1(retval);
    } else if (type == "int32") {
      const std::int32_t val =
        util::getAttributeValue<std::int32_t>(conf, "global metadata", attname);
      if ((retval = nc_put_att_int(netcdfGeneralIDs[0], NC_GLOBAL, attname.data(),
                                   NC_INT, 1, &val))) ERR1(retval);
    } else if (type == "int64") {
      const long long int val =  //NOLINT
        util::getAttributeValue<std::int64_t>(conf, "global metadata", attname);
      if ((retval = nc_put_att_longlong(netcdfGeneralIDs[0], NC_GLOBAL, attname.data(),
                                        NC_INT64, 1, &val))) ERR1(retval);
    } else if (type == "real32") {
      const float val =
        util::getAttributeValue<float>(conf, "global metadata", attname);
      if ((retval = nc_put_att_float(netcdfGeneralIDs[0], NC_GLOBAL, attname.data(),
                                     NC_FLOAT, 1, &val))) ERR1(retval);
    } else if (type == "real64") {
      const double val =
        util::getAttributeValue<double>(conf, "global metadata", attname);
      if ((retval = nc_put_att_double(netcdfGeneralIDs[0], NC_GLOBAL, attname.data(),
                                      NC_DOUBLE, 1, &val))) ERR1(retval);
    }
  }

  for (std::size_t i = 0; i < dimNames.size(); ++i) {
    int dim_id;
    if ((retval = nc_def_dim(netcdfGeneralIDs[0],
                             dimNames[i].c_str(),
                             dimSizes[i], &dim_id))) ERR1(retval);
    netcdfDimIDs.push_back(dim_id);
  }

  // dim names and dim index are in the same order -> maybe use pair?
  std::vector<int> dimVarIDs;

  std::vector<std::vector<std::string>>   vakeys;
  std::vector<std::vector<std::int32_t>>  vavalsint32;
  std::vector<std::vector<long long int>> vavalsint64;  //NOLINT
  std::vector<std::vector<float>>         vavalsfloat;
  std::vector<std::vector<double>>        vavalsdouble;
  std::vector<std::vector<std::string>>   vavalsstring;
  std::vector<std::vector<std::string>>   varMetDatAttTypes;

  for (std::size_t i = 0; i < variables.size(); ++i) {
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

    if (variables[i].dataType() == oops::ModelDataType::Real64) {
      if ((retval = nc_def_var(netcdfGeneralIDs[0],
                               variables[i].name().c_str(),
                               NC_DOUBLE,
                               dimVarIDs.size(),
                               dimVarIDs.data(),
                               &varID))) ERR1(retval);

      if ((retval = nc_put_att_double(netcdfGeneralIDs[0],
                                      varID, "_FillValue",
                                      NC_DOUBLE, 1,
                                      &msvd))) ERR1(retval);

    } else if (variables[i].dataType() == oops::ModelDataType::Int32) {
      if ((retval = nc_def_var(netcdfGeneralIDs[0],
                               variables[i].name().c_str(),
                               NC_INT,
                               dimVarIDs.size(),
                               dimVarIDs.data(),
                               &varID))) ERR1(retval);

      if ((retval = nc_put_att_int(netcdfGeneralIDs[0],
                                   varID, "_FillValue",
                                   NC_INT, 1,
                                   &msvi))) ERR1(retval);
    } else {
      throw eckit::UserError("datatype currently not supported", Here());
    }

    const std::vector<std::string> varAttnames = util::getAttributeNames(conf, variables[i].name());
    for (const std::string & attname : varAttnames) {
      const std::string type =
        util::getAttributeType(conf, variables[i].name(), attname);
      if (type == "string") {
        const std::string val =
          util::getAttributeValue<std::string>(conf, variables[i].name(), attname);
        if ((retval = nc_put_att_text(netcdfGeneralIDs[0], varID, attname.data(),
                                      val.length() + 1, val.data() ))) ERR1(retval);
      } else if (type == "int32") {
        const std::int32_t val =
          util::getAttributeValue<std::int32_t>(conf, variables[i].name(), attname);
        if ((retval = nc_put_att_int(netcdfGeneralIDs[0], varID, attname.data(),
                                     NC_INT, 1, &val))) ERR1(retval);
      } else if (type == "int64") {
        const long long int val =  //NOLINT
          util::getAttributeValue<std::int64_t>(conf, variables[i].name(), attname);
        if ((retval = nc_put_att_longlong(netcdfGeneralIDs[0], varID, attname.data(),
                                          NC_INT64, 1, &val))) ERR1(retval);
      } else if (type == "real32") {
        const float val =
          util::getAttributeValue<float>(conf, variables[i].name(), attname);
        if ((retval = nc_put_att_float(netcdfGeneralIDs[0], varID, attname.data(),
                                       NC_FLOAT, 1, &val))) ERR1(retval);
      } else if (type == "real64") {
        const double val =
          util::getAttributeValue<double>(conf, variables[i].name(), attname);
        if ((retval = nc_put_att_double(netcdfGeneralIDs[0], varID, attname.data(),
                                        NC_DOUBLE, 1, &val))) ERR1(retval);
      }
    }

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
    oops::Variables & variables,
    std::vector<std::vector<std::string>> & dimNamesForEveryVar,
    eckit::LocalConfiguration & conf,
    std::vector<int> & netcdfGeneralIDs,
    std::vector<int> & netcdfDimIDs,
    std::vector<int> & netcdfVarIDs,
    std::vector<std::vector<int>> & netcdfDimVarIDs) {
  dimNames.clear();
  dimSizes.clear();
  oops::Variables variablesWork;
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

  std::vector<nc_type> nctypes(ngatts, 0);


  // Populate the global Configuration header
  // Assume all metadata is in string format.
  for (int gattsid = 0; gattsid < ngatts; ++gattsid) {
    char attname[NC_MAX_NAME + 1];
    std::fill(attname, attname + NC_MAX_NAME + 1, 0);
    if ((retval =
         nc_inq_attname(netcdfGeneralIDs[0], NC_GLOBAL,
                        gattsid, attname))) ERR1(retval);
    if ((retval =
         nc_inq_atttype(netcdfGeneralIDs[0], NC_GLOBAL,
                        attname, &nctypes[gattsid]))) ERR1(retval);

    // We ignore the "history" global attribute as it has
    // special internal structure different to standard attributes
    if (std::string {attname} == "history") {
      continue;
    }

    if (nctypes[gattsid] == NC_CHAR) {
      char attvalue[NC_MAX_NAME + 1];
      std::fill(attvalue, attvalue + NC_MAX_NAME + 1, 0);
      if ((retval = nc_get_att(netcdfGeneralIDs[0], NC_GLOBAL,
                               attname, attvalue))) ERR1(retval);
      util::setAttribute<std::string>(conf, "global metadata",
        std::string{attname}, "string", std::string{attvalue});
    } else if (nctypes[gattsid] == NC_INT) {
      std::int32_t attvalue_int32;
      if ((retval = nc_get_att(netcdfGeneralIDs[0], NC_GLOBAL,
                               attname, &attvalue_int32))) ERR1(retval);
      util::setAttribute<std::int32_t>(conf, "global metadata",
        std::string{attname}, "int32", attvalue_int32);
    } else if (nctypes[gattsid] == NC_INT64) {
      std::int64_t attvalue_int64;
      if ((retval = nc_get_att(netcdfGeneralIDs[0], NC_GLOBAL,
                               attname, &attvalue_int64))) ERR1(retval);
      util::setAttribute<std::int64_t>(conf, "global metadata",
        std::string{attname}, "int64", attvalue_int64);
    } else if (nctypes[gattsid] == NC_FLOAT) {
      float attvalue_float;
      if ((retval = nc_get_att(netcdfGeneralIDs[0], NC_GLOBAL,
                               attname, &attvalue_float))) ERR1(retval);
      util::setAttribute<float>(conf, "global metadata",
        std::string{attname}, "real32", attvalue_float);
    } else if (nctypes[gattsid] == NC_DOUBLE) {
      double attvalue_double;
      if ((retval = nc_get_att(netcdfGeneralIDs[0], NC_GLOBAL,
                               attname, &attvalue_double))) ERR1(retval);
      util::setAttribute<double>(conf, "global metadata",
        std::string{attname}, "real64", attvalue_double);
    } else {
      oops::Log::warning() << "attribute type not accounted for nctype number " +
                              std::to_string(nctypes[gattsid]) << std::endl;
    }
  }

  for (int dimid = 0; dimid < ndims; ++dimid) {
    char cdim_name[NC_MAX_NAME + 1];
    std::fill(cdim_name, cdim_name + NC_MAX_NAME + 1, 0);
    size_t lenp;
    if ((retval = nc_inq_dim(netcdfGeneralIDs[0], dimid, cdim_name, &lenp))) ERR1(retval);
    dimNames.push_back(std::string(cdim_name));
    dimSizes.push_back(atlas::idx_t(lenp));
    netcdfDimIDs.push_back(dimid);
  }

  for (int varid = 0; varid < nvars; ++varid) {
    char cvar_name[NC_MAX_NAME + 1];
    std::fill(cvar_name, cvar_name + NC_MAX_NAME + 1, 0);
    nc_type nc_var_type;
    int var_ndims;
    int var_dimids[NC_MAX_VAR_DIMS];
    int var_natts;
    if ((retval = nc_inq_var(netcdfGeneralIDs[0], varid, cvar_name, &nc_var_type,
                             &var_ndims, var_dimids, &var_natts))) ERR1(retval);
    variablesWork.push_back(std::string(cvar_name));
    netcdfVarIDs.push_back(varid);

    // getting all at the string based attributes.
    // The first attribute is always the _FillValue which we do not store in the variable metadata.
    std::vector<nc_type> varnctypes(var_natts, 0);

    for (int vasid = 0; vasid < var_natts; ++vasid) {
      char attname[NC_MAX_NAME + 1];
      std::fill(attname, attname + NC_MAX_NAME + 1, 0);
      if ((retval = nc_inq_attname(netcdfGeneralIDs[0], varid, vasid, attname))) ERR1(retval);

      if ((std::string {attname} == "_FillValue") ||
          (std::string {attname} == "type")) {
        continue;
      }

      if ((retval =
           nc_inq_atttype(netcdfGeneralIDs[0], varid,
                          attname, &varnctypes[vasid]))) ERR1(retval);

      if (varnctypes[vasid] == NC_CHAR) {
        char vattvalue[NC_MAX_NAME + 1];
        std::fill(vattvalue, vattvalue + NC_MAX_NAME + 1, 0);
        if ((retval = nc_get_att(netcdfGeneralIDs[0], varid,
                                 attname, vattvalue))) ERR1(retval);
        util::setAttribute<std::string>(conf, std::string{cvar_name},
          std::string{attname}, "string", std::string{vattvalue});
      } else if (varnctypes[vasid] == NC_INT) {
        std::int32_t attvalue_int32;
        if ((retval = nc_get_att(netcdfGeneralIDs[0], varid,
                                 attname, &attvalue_int32))) ERR1(retval);
        util::setAttribute<std::int32_t>(conf, std::string{cvar_name},
          std::string{attname}, "int32", attvalue_int32);
      } else if (varnctypes[vasid] == NC_INT64) {
        std::int64_t attvalue_int64;
        if ((retval = nc_get_att(netcdfGeneralIDs[0], varid,
                                 attname, &attvalue_int64))) ERR1(retval);
        util::setAttribute<std::int64_t>(conf, std::string{cvar_name},
          std::string{attname}, "int64", attvalue_int64);
      } else if (varnctypes[vasid] == NC_FLOAT) {
        float attvalue_float;
        if ((retval = nc_get_att(netcdfGeneralIDs[0], varid,
                                 attname, &attvalue_float))) ERR1(retval);
        util::setAttribute<float>(conf, std::string{cvar_name},
          std::string{attname}, "real32", attvalue_float);
      } else if (varnctypes[vasid] == NC_DOUBLE) {
        double attvalue_double;
        if ((retval = nc_get_att(netcdfGeneralIDs[0], varid,
                                 attname, &attvalue_double))) ERR1(retval);
        util::setAttribute<double>(conf, std::string{cvar_name},
          std::string{attname}, "real64", attvalue_double);
      } else {
        oops::Log::warning() << "attribute type not accounted for nctype number " +
                                std::to_string(nctypes[vasid]) << std::endl;
      }
    }

    std::vector<int> netcdf_dim_varID;
    std::vector<std::string> dimNamesForEachVar;
    for (int var_dim = 0; var_dim < var_ndims; ++var_dim) {
      netcdf_dim_varID.push_back(var_dimids[var_dim]);
      char cdimv_name[NC_MAX_NAME + 1];
      std::fill(cdimv_name, cdimv_name + NC_MAX_NAME + 1, 0);
      size_t lenp;
      if ((retval = nc_inq_dim(netcdfGeneralIDs[0],
                               var_dimids[var_dim],
                               cdimv_name,
                               &lenp))) ERR1(retval);
      dimNamesForEachVar.push_back(std::string{cdimv_name});
    }
    netcdfDimVarIDs.push_back(netcdf_dim_varID);
    dimNamesForEveryVar.push_back(dimNamesForEachVar);
  }

  variables = variablesWork;
}


void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 1> & arrayInOut) {
  int retval;
  std::vector<double> zvar(dimSizes[0]);

  if ((retval = nc_get_var_double(netcdfGeneralIDs[0],
                                  varID, zvar.data()))) ERR1(retval);

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
  std::vector<double> zvar(dimSizes[0] * dimSizes[1]);

  if ((retval = nc_get_var_double(netcdfGeneralIDs[0],
                                  varID, zvar.data()))) ERR1(retval);

  for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
    for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizes[1]); ++t1) {
      arrayInOut(t, t1) = zvar[t*dimSizes[1] + t1];
    }
  }
}

void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 3> & arrayInOut) {
  int retval;
  std::vector<double> zvar(dimSizes[0] * dimSizes[1] * dimSizes[2]);

  if ((retval = nc_get_var_double(netcdfGeneralIDs[0],
                                  varID, zvar.data()))) ERR1(retval);

  for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
    for (idx_t t1 = 0; t1 < static_cast<idx_t>(dimSizes[1]); ++t1) {
      for (idx_t t2 = 0; t2 < static_cast<idx_t>(dimSizes[2]); ++t2) {
        arrayInOut(t, t1, t2) = zvar[t*dimSizes[1]*dimSizes[2] + t1*dimSizes[2] + t2];
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
  std::vector<int> zvar(dimSizes[0]);

  if ((retval = nc_get_var_int(netcdfGeneralIDs[0],
                               varID, zvar.data()))) ERR1(retval);

  for (idx_t t = 0; t < static_cast<idx_t>(dimSizes[0]); ++t) {
    arrayInOut(t) = zvar[t];
  }
}


void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 1> & arrayIn) {
  int retval;
  std::vector<double> zvar(arrayIn.shape()[0]);
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
    zvar[t] = arrayIn(t);
  }

  if ((retval = nc_put_var_double(netcdfGeneralIDs[0], varID, zvar.data()))) ERR1(retval);
}


void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 2> & arrayIn) {
  int retval;
  std::vector<double> zvar(arrayIn.shape()[0] * arrayIn.shape()[1]);
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayIn.shape()[1]; ++t1) {
      zvar[t*arrayIn.shape()[1] + t1] = arrayIn(t, t1);
    }
  }

  if ((retval = nc_put_var_double(netcdfGeneralIDs[0], varID, zvar.data()))) ERR1(retval);
}

void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 3> & arrayIn) {
  int retval;
  std::vector<double> zvar(arrayIn.shape()[0] * arrayIn.shape()[1] * arrayIn.shape()[2]);
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayIn.shape()[1]; ++t1) {
      for (idx_t t2 = 0; t2 < arrayIn.shape()[2]; ++t2) {
        zvar[t*arrayIn.shape()[1]*arrayIn.shape()[2] + t1*arrayIn.shape()[2] + t2] =
          arrayIn(t, t1, t2);
      }
    }
  }

  if ((retval = nc_put_var_double(netcdfGeneralIDs[0], varID, zvar.data()))) ERR1(retval);
}

void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const int, 1> & arrayIn) {
  int retval;
  std::vector<int> zvar(arrayIn.shape()[0]);
  for (idx_t t = 0; t < arrayIn.shape()[0]; ++t) {
      zvar[t] = arrayIn(t);
  }

  if ((retval = nc_put_var_int(netcdfGeneralIDs[0], varID, zvar.data()))) ERR1(retval);
}


}  // namespace util
