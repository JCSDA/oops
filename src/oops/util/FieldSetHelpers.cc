/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FieldSetHelpers.h"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <memory>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/util/function/VortexRollup.h"

#include "eckit/utils/Hash.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/RandomField.h"

#define ERR(e) {ABORT(nc_strerror(e));}

namespace util {

// -----------------------------------------------------------------------------
atlas::FieldSet createRandomFieldSet(const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspace,
                                     const std::vector<size_t> & variableSizes,
                                     const std::vector<std::string> & vars) {
  oops::Log::trace() << "createRandomFieldSet starting" << std::endl;

  // Local ghost points
  atlas::Field localGhost;
  if (fspace.type() != "Spectral") {
    localGhost = fspace.ghost();
  }

  // Global ghost points
  atlas::Field globalGhost;
  if (fspace.type() != "PointCloud") {
    globalGhost = fspace.createField<int>(atlas::option::name("ghost")
     | atlas::option::global());

    // Gather masks on main processor
    if (fspace.type() == "StructuredColumns") {
      // StructuredColumns
      const atlas::functionspace::StructuredColumns fs(fspace);
      fs.gather(localGhost, globalGhost);
    } else if (fspace.type() == "NodeColumns") {
      // NodeColumns
      const atlas::functionspace::CubedSphereNodeColumns fs(fspace);
      fs.gather(localGhost, globalGhost);
      /* TODO(??): have to assume the NodeColumns is CubedSphere here, cannot differentiate with
         other NodeColumns from what is in the fspace.
      const atlas::functionspace::NodeColumns fs(fspace);
      fs.gather(localGhost, globalGhost);
      */
    } else if (fspace.type() == "Spectral") {
      // Not needed
    } else {
      ABORT(fspace.type() + " function space not supported yet");
    }
  }

  // Create FieldSet
  atlas::FieldSet fset;

  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    // Create field
    atlas::Field field = fspace.createField<double>(
      atlas::option::name(vars[jvar]) | atlas::option::levels(variableSizes[jvar]));

    // Get field owned size
    size_t n = 0;
    if (field.rank() == 2) {
      if (fspace.type() == "Spectral") {
        const atlas::functionspace::Spectral fs(fspace);
        const atlas::idx_t N = fs.truncation();
        const auto zonal_wavenumbers = fs.zonal_wavenumbers();
        const atlas::idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();
        for (int jm=0; jm < nb_zonal_wavenumbers; ++jm) {
          const atlas::idx_t m1 = zonal_wavenumbers(jm);
          for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
            if (m1 == 0) {
              n += field.shape(1);
            } else {
              n += 2*field.shape(1);
            }
          }
        }
      } else {
        auto localGhostView = atlas::array::make_view<int, 1>(localGhost);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          if (localGhostView(jnode) == 0) n += field.shape(1);
        }
      }
    }

    // Gather local sizes
    std::vector<size_t> nlocs(comm.size());
    comm.gather(n, nlocs, 0);

    // Get global size
    size_t nglb = 0;
    if (comm.rank() == 0) {
      for (const size_t & nloc : nlocs)
        nglb += nloc;
    }

    // Global field
    atlas::Field globalField = fspace.createField<double>(
      atlas::option::name(vars[jvar]) | atlas::option::levels(variableSizes[jvar])
      | atlas::option::global());

    std::vector<double> rand_vec_glb(nglb);
    if (comm.rank() == 0) {
      // Generate global random vector
      util::NormalDistributionField dist(nglb, 0.0, 1.0);
      for (size_t i = 0; i < nglb; ++i) {
        rand_vec_glb[i] = dist[i];
      }

      if (fspace.type() != "PointCloud") {
        // Copy to field
        n = 0;
        if (globalField.rank() == 2) {
          auto view = atlas::array::make_view<double, 2>(globalField);
          if (fspace.type() == "Spectral") {
            const atlas::functionspace::Spectral fs(fspace);
            const atlas::idx_t N = fs.truncation();
            int jnode = 0;
            for (int jm=0; jm <= N; ++jm) {
              const atlas::idx_t m1 = jm;
              for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
                if (m1 == 0) {
                  // Real part only
                  for (atlas::idx_t jlevel = 0; jlevel < globalField.shape(1); ++jlevel) {
                    view(jnode, jlevel) = dist[n];
                    ++n;
                  }
                  ++jnode;

                  // No imaginary part
                  for (atlas::idx_t jlevel = 0; jlevel < globalField.shape(1); ++jlevel) {
                    view(jnode, jlevel) = 0.0;
                  }
                  ++jnode;
                } else {
                  // Real part
                  for (atlas::idx_t jlevel = 0; jlevel < globalField.shape(1); ++jlevel) {
                    view(jnode, jlevel) = dist[n] * M_SQRT1_2;
                    ++n;
                  }
                  ++jnode;

                  // Imaginary part
                  for (atlas::idx_t jlevel = 0; jlevel < globalField.shape(1); ++jlevel) {
                    view(jnode, jlevel) = dist[n] * M_SQRT1_2;
                    ++n;
                  }
                  ++jnode;
                }
              }
            }
          } else {
            auto globalGhostView = atlas::array::make_view<int, 1>(globalGhost);
            for (atlas::idx_t jnode = 0; jnode < globalField.shape(0); ++jnode) {
              if (globalGhostView(jnode) == 0) {
                for (atlas::idx_t jlevel = 0; jlevel < globalField.shape(1); ++jlevel) {
                  view(jnode, jlevel) = dist[n];
                  ++n;
                }
              }
            }
          }
        }
      }
    }

    if (fspace.type() == "PointCloud") {
      // Scatter random vector
      std::vector<int> sendcounts(comm.size());
      std::vector<int> displs(comm.size());
      if (comm.rank() == 0) {
        int sum = 0;
        for (size_t i = 0; i < comm.size(); ++i) {
          sendcounts[i] = nlocs[i];
          displs[i] = sum;
          sum += sendcounts[i];
        }
      }
      std::vector<double> rand_vec_loc(n);
      comm.scatterv(rand_vec_glb.cbegin(), rand_vec_glb.cend(), sendcounts, displs,
      rand_vec_loc.begin(), rand_vec_loc.end(), 0);

      // Populate with random numbers
      n = 0;
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) = rand_vec_loc[n];
            ++n;
          }
        }
      }
    } else {
      // Scatter global field
      if (fspace.type() == "StructuredColumns") {
        // StructuredColumns
        const atlas::functionspace::StructuredColumns fs(fspace);
        fs.scatter(globalField, field);
      } else if (fspace.type() == "NodeColumns") {
        // CubedSphere
        const atlas::functionspace::CubedSphereNodeColumns fs(fspace);
        fs.scatter(globalField, field);
        /* TODO(??): have to assume the NodeColumns is CubedSphere here, cannot differentiate with
           other NodeColumns from what is in fspace.
        const atlas::functionspace::NodeColumns fs(fspace);
        fs.scatter(globalField, field);
        */
      } else if (fspace.type() == "Spectral") {
        const atlas::functionspace::Spectral fs(fspace);
        fs.scatter(globalField, field);
      } else {
        ABORT(fspace.type() + " function space not supported yet");
      }
    }

    // Set metadata for interpolation type
    if (fspace.type() != "Spectral") {
      field.metadata().set("interp_type", "default");
    }

    // Add field
    fset.add(field);
  }

  if (fspace.type() != "Spectral") {
    // Halo exchange
    fset.haloExchange();
  }

  // Return FieldSet
  return fset;
}

// -----------------------------------------------------------------------------
atlas::FieldSet createSmoothFieldSet(const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspace,
                                     const std::vector<size_t> & variableSizes,
                                     const std::vector<std::string> & vars) {
  if (fspace.type() == "Spectral") {
    // Create random spectral FieldSet
    atlas::FieldSet fset = createRandomFieldSet(comm, fspace, variableSizes, vars);

    // Convolve with Gaussian
    auto gaussian = [](double dist, double sigma){
      return std::exp(- dist * dist / (2 * sigma * sigma));};

    const atlas::functionspace::Spectral fs(fspace);
    const atlas::idx_t N = fs.truncation();
    const double sigma = N/4;
    const auto zonal_wavenumbers = fs.zonal_wavenumbers();
    const atlas::idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();
    for (auto & field : fset) {
      auto view = atlas::array::make_view<double, 2>(field);
      atlas::idx_t jnode = 0;
      for (int jm=0; jm < nb_zonal_wavenumbers; ++jm) {
        const atlas::idx_t m1 = zonal_wavenumbers(jm);
        for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
          for (const auto & part : {"real", "imaginary"}) {
            for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              if (n1 == 0 && (field.name() == "streamfunction" ||
                              field.name() == "velocity_potential")) {
                view(jnode, jlevel) = 0.0;  // Force zero mean
              } else {
                view(jnode, jlevel) *= gaussian(n1, sigma);
              }
            }
            (void)part;  // unused
            jnode++;
          }
        }
      }
    }
    return fset;
  }

  // Create FieldSet
  atlas::FieldSet fset;
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    // Create field
    atlas::Field field = fspace.createField<double>(
      atlas::option::name(vars[jvar]) | atlas::option::levels(variableSizes[jvar]));
    const auto lonlat = atlas::array::make_view<double, 2>(field.functionspace().lonlat());
    if (field.rank() == 2) {
      size_t nlev = field.shape(1);
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (size_t jlev = 0; jlev < nlev; ++jlev) {
          const double zz = static_cast<double>(jlev) / static_cast<double>(nlev);
          view(jnode, jlev) = atlas::util::function::vortex_rollup(lonlat(jnode, 0),
                                lonlat(jnode, 1), zz) * 2.0;
        }
      }
    }
    fset.add(field);
  }

  return fset;
}

// -----------------------------------------------------------------------------

void copyFieldSet(const atlas::FieldSet & otherFset, atlas::FieldSet & fset) {
  oops::Log::trace() << "copyFieldSet starting" << std::endl;
  fset.clear();
  for (const auto & otherField : otherFset) {
    // Create Field
    atlas::Field field = otherField.functionspace().createField<double>(
      atlas::option::name(otherField.name()) | atlas::option::levels(otherField.levels()));

    // Copy data
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto otherView = atlas::array::make_view<double, 2>(otherField);
      for (atlas::idx_t jnode = 0; jnode < otherField.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < otherField.shape(1); ++jlevel) {
            view(jnode, jlevel) = otherView(jnode, jlevel);
        }
      }
    } else {
      ABORT("copyFieldSet: wrong rank");
    }

    // Copy metadata
    field.metadata() = otherField.metadata();

    // Add field
    fset.add(field);
  }
}

// -----------------------------------------------------------------------------

atlas::FieldSet copyFieldSet(const atlas::FieldSet & otherFset) {
  oops::Log::trace() << "copyFieldSet starting" << std::endl;
  // Create FieldSet
  atlas::FieldSet fset;

  // Copy FieldSet
  copyFieldSet(otherFset, fset);

  // Return FieldSet
  return fset;
}

// -----------------------------------------------------------------------------

void shareFields(const atlas::FieldSet & otherFset, atlas::FieldSet & fset) {
  oops::Log::trace() << "shareFields starting" << std::endl;
  // Add other FieldSet fields
  for (const auto & field : otherFset) {
    fset.add(field);
  }
}

// -----------------------------------------------------------------------------

atlas::FieldSet shareFields(const atlas::FieldSet & otherFset) {
  oops::Log::trace() << "shareFields starting" << std::endl;
  // Create FieldSet
  atlas::FieldSet fset;

  // Share fields
  shareFields(otherFset, fset);

  // Return FieldSet
  return fset;
}

// -----------------------------------------------------------------------------

void removeFieldsFromFieldSet(atlas::FieldSet & fset,
                              const std::vector<std::string> & vars) {
  oops::Log::trace() << "removeFieldsFromFieldSet starting" << std::endl;

  // Create FieldSet
  atlas::FieldSet fsetTmp;

  for (const auto & field : fset) {
    if (std::find(vars.begin(), vars.end(), field.name()) == vars.end()) {
      // Add field
      fsetTmp.add(field);
    }
  }

  // Replace FieldSet
  fset = fsetTmp;

  oops::Log::trace() << "removeFieldsFromFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

bool compareFieldSets(const atlas::FieldSet & fset1,
                      const atlas::FieldSet & fset2,
                      const double & tol,
                      const bool & absolute) {
  oops::Log::trace() << "compareFieldSets starting" << std::endl;

  // Initialize flag
  bool sameFieldSets = true;

  // Compare FieldSets content
  std::vector<std::string> field_names1 = fset1.field_names();
  std::vector<std::string> field_names2 = fset2.field_names();
  std::sort(field_names1.begin(), field_names1.end());
  std::sort(field_names2.begin(), field_names2.end());
  sameFieldSets = (field_names1 == field_names2);
  if (!sameFieldSets) return sameFieldSets;

  for (const auto & field1 : fset1) {
    // Get second field
    atlas::Field field2 = fset2.field(field1.name());

    // Compare function spaces
    sameFieldSets = (field1.functionspace().type() == field2.functionspace().type());
    if (!sameFieldSets) return sameFieldSets;

    // Compare fields shapes
    sameFieldSets = (field1.shape() == field2.shape());
    if (!sameFieldSets) return sameFieldSets;

    // Compare data
    if (field1.rank() == 2) {
      auto view1 = atlas::array::make_view<double, 2>(field1);
      auto view2 = atlas::array::make_view<double, 2>(field2);
      for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          if (absolute) {
            sameFieldSets = oops::is_close_absolute(view1(jnode, jlevel), view2(jnode, jlevel),
              tol);
          } else {
            sameFieldSets = oops::is_close_relative(view1(jnode, jlevel), view2(jnode, jlevel),
              tol);
          }
          if (!sameFieldSets) return sameFieldSets;
        }
      }
    } else {
      ABORT("compareFieldSets: wrong rank");
    }
  }

  // Comparison successuful!
  return sameFieldSets;
}

// -----------------------------------------------------------------------------

std::string getGridUid(const atlas::FunctionSpace & fspace) {
  oops::Log::trace() << "getGridUid starting" << std::endl;

  // WARNING: this is a local ID per MPI task!
  // There is an unlikely failure mode where two FunctionSpaces have equal lonlats on some MPI
  // tasks but different lonlats on other MPI tasks. In this case the local IDs would compare equal
  // on some tasks and non-equal on others, which could lead to serious bugs.
  // The fix would be to allGather the strings and re-hash, ensuring a globally-consistent ID.
  auto customUidFromLonLat = [&](const atlas::FunctionSpace & fspace) {
    std::unique_ptr<eckit::Hash> hash(eckit::HashFactory::instance().build("md5"));
    // Add function space size to hash
    hash->add(fspace.size());
    // Add function space lon lat to hash
    const auto lonlatView = atlas::array::make_view<double, 2>(fspace.lonlat());
    for (size_t i = 0; i < lonlatView.shape(0) ; ++i) {
      hash->add(lonlatView(i, 0));
      hash->add(lonlatView(i, 1));
    }
    return hash->digest();
  };

  if (fspace.type() == "StructuredColumns") {
    // StructuredColumns
    const atlas::functionspace::StructuredColumns fs(fspace);
    return fs.grid().uid();
  } else if (fspace.type() == "NodeColumns") {
    // CubedSphere
    const atlas::functionspace::CubedSphereNodeColumns fs(fspace);
    return fs.mesh().grid().uid();
    /* TODO(??): have to assume the NodeColumns is CubedSphere here, cannot differentiate with
       other NodeColumns from what is in fspace.
    const atlas::functionspace::NodeColumns fs(fspace);
    return fs.grid().uid();
    */
  } else if (fspace.type() == "PointCloud") {
    return customUidFromLonLat(fspace);
  } else if (fspace.type() == "Spectral") {
    const atlas::functionspace::Spectral fs(fspace);
    return "Spectral" + std::to_string(fs.truncation());
  } else {
    ABORT(fspace.type() + " function space not supported yet");
    return "";
  }
}

// -----------------------------------------------------------------------------

std::string getGridUid(const atlas::FieldSet & fset) {
  oops::Log::trace() << "getGridUid starting" << std::endl;

  if (fset.size() > 0) {
    // Get grid UID of the first field
    std::string uid = getGridUid(fset[0].functionspace());

    // Check that other fields have the same UID
    for (const auto & field : fset) {
      if (getGridUid(field.functionspace()) != uid) {
        ABORT("All fields should have the same grid");
      }
    }

    // Return UID
    return uid;
  } else {
    // No field in the fieldset
    return "";
  }
}

// -----------------------------------------------------------------------------

void printDiagValues(const eckit::mpi::Comm & timeComm,
                     const eckit::mpi::Comm & comm,
                     const atlas::FunctionSpace & fspace,
                     const atlas::FieldSet & dataFset,
                     const atlas::FieldSet & diagFset) {
  oops::Log::trace() << "printDiagValues starting" << std::endl;

  // Global lon/lat field
  atlas::FieldSet globalCoords;
  atlas::Field lonGlobal = fspace.createField<double>(
    atlas::option::name("lon") | atlas::option::global() | atlas::option::levels(0));
  globalCoords.add(lonGlobal);
  atlas::Field latGlobal = fspace.createField<double>(
    atlas::option::name("lat") | atlas::option::global() | atlas::option::levels(0));
  globalCoords.add(latGlobal);
  auto lonViewGlobal = atlas::array::make_view<double, 1>(lonGlobal);
  auto latViewGlobal = atlas::array::make_view<double, 1>(latGlobal);
  const atlas::Field localLonlat = fspace.lonlat();
  const auto lonlatView = atlas::array::make_view<double, 2>(localLonlat);
  if (fspace.type() == "PointCloud") {
    // Copy local
    for (atlas::idx_t jnode = 0; jnode < localLonlat.shape(0); ++jnode) {
       lonViewGlobal(jnode) = lonlatView(jnode, 0);
       latViewGlobal(jnode) = lonlatView(jnode, 1);
    }
  } else {
    // Gather local
    atlas::FieldSet localCoords;
    atlas::Field lonLocal = fspace.createField<double>(
      atlas::option::name("lon") | atlas::option::levels(0));
    localCoords.add(lonLocal);
    atlas::Field latLocal = fspace.createField<double>(
      atlas::option::name("lat") | atlas::option::levels(0));
    localCoords.add(latLocal);
    auto lonViewLocal = atlas::array::make_view<double, 1>(lonLocal);
    auto latViewLocal = atlas::array::make_view<double, 1>(latLocal);
    for (atlas::idx_t jnode = 0; jnode < localLonlat.shape(0); ++jnode) {
       lonViewLocal(jnode) = lonlatView(jnode, 0);
       latViewLocal(jnode) = lonlatView(jnode, 1);
    }
    fspace.gather(localCoords, globalCoords);
  }

  for (const auto & diagField : diagFset) {
    // Get data field with the same name
    const atlas::Field dataField = dataFset.field(diagField.name());

    // Global fields
    atlas::Field globalDataField;
    atlas::Field globalDiagField;
    if (fspace.type() == "PointCloud") {
      // Copy local
      globalDataField = dataField;
      globalDiagField = diagField;
    } else {
      // Gather local
      globalDataField = fspace.createField<double>(
        atlas::option::name(dataField.name()) | atlas::option::levels(dataField.levels())
        | atlas::option::global());
      globalDiagField = fspace.createField<double>(
        atlas::option::name(diagField.name()) | atlas::option::levels(diagField.levels())
        | atlas::option::global());
      fspace.gather(dataField, globalDataField);
      fspace.gather(diagField, globalDiagField);
    }

    if (comm.rank() == 0) {
      // Print data values at diag points
      std::vector<double> lons;
      std::vector<double> lats;
      std::vector<size_t> levs;
      std::vector<size_t> subWindows;
      std::vector<double> values;

      if (globalDiagField.rank() == 2) {
        auto dataView = atlas::array::make_view<double, 2>(globalDataField);
        auto diagView = atlas::array::make_view<double, 2>(globalDiagField);
        for (int jnode = 0; jnode < globalDiagField.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < globalDiagField.shape(1); ++jlevel) {
            if (std::abs(diagView(jnode, jlevel) - 1.0) < 1.0e-12) {
              // Diagnostic point found
              lons.push_back(lonViewGlobal(jnode));
              lats.push_back(latViewGlobal(jnode));
              levs.push_back(jlevel+1);
              subWindows.push_back(timeComm.rank());
              values.push_back(dataView(jnode, jlevel));
            }
          }
        }
      } else {
        ABORT("getDiagValues: wrong rank");
      }

      // Gather sizes
      int size = lons.size();
      std::vector<int> sizes(timeComm.size());
      timeComm.gather(size, sizes, 0);

      // Gather data
      std::vector<int> displs;
      std::vector<int> recvcounts;
      if (timeComm.rank() == 0) {
        displs.resize(timeComm.size());
        recvcounts.resize(timeComm.size());
        for (size_t i = 0; i < timeComm.size(); ++i) {
          recvcounts[i] = sizes[i];
          displs[i] = static_cast<int>(i ? displs[i - 1] + recvcounts[i - 1] : 0);
        }
      }
      size_t recvsize = size_t(std::accumulate(recvcounts.begin(), recvcounts.end(), 0));
      std::vector<double> lonsOnRoot(recvsize);
      std::vector<double> latsOnRoot(recvsize);
      std::vector<size_t> levsOnRoot(recvsize);
      std::vector<size_t> subWindowsOnRoot(recvsize);
      std::vector<double> valuesOnRoot(recvsize);
      timeComm.gatherv(lons, lonsOnRoot, recvcounts, displs, 0);
      timeComm.gatherv(lats, latsOnRoot, recvcounts, displs, 0);
      timeComm.gatherv(levs, levsOnRoot, recvcounts, displs, 0);
      timeComm.gatherv(subWindows, subWindowsOnRoot, recvcounts, displs, 0);
      timeComm.gatherv(values, valuesOnRoot, recvcounts, displs, 0);

      // Print results
      if (timeComm.rank() == 0) {
        for (size_t i = 0; i < lonsOnRoot.size(); ++i) {
          oops::Log::test() << "  + Value for variable " << diagField.name()
                            <<", subwindow " << subWindowsOnRoot[i]
                            << std::fixed << std::setprecision(5)
                            << ", at (longitude, latitude, vertical index) point ("
                            << lonsOnRoot[i] << ", " << latsOnRoot[i]
                            << ", " << levsOnRoot[i] << "): "
                            << std::scientific << std::setprecision(16)
                            << valuesOnRoot[i] << std::endl;
        }
      }
    }
  }

  oops::Log::trace() << "printDiagValues done" << std::endl;
}

// -----------------------------------------------------------------------------

void readFieldSet(const eckit::mpi::Comm & comm,
                  const atlas::FunctionSpace & fspace,
                  const std::vector<size_t> & variableSizes,
                  const std::vector<std::string> & vars,
                  const eckit::Configuration & config,
                  atlas::FieldSet & fset) {
  // Options with one file per MPI task
  const bool oneFilePerTask = config.getBool("one file per task", false);
  ASSERT(oneFilePerTask || (fspace.type() != "PointCloud"));

  // Build filepath
  std::string filepath = config.getString("filepath");
  if (config.has("member")) {
    std::ostringstream out;
    out << std::setfill('0') << std::setw(6) << config.getInt("member");
    filepath.append("_");
    filepath.append(out.str());
  }
  if (oneFilePerTask) {
    size_t zpad = 6;
    std::string commSize = std::to_string(comm.size());
    commSize = std::string(zpad - std::min(zpad, commSize.length()), '0') + commSize;
    std::string commRank = std::to_string(comm.rank()+1);
    commRank = std::string(zpad - std::min(zpad, commRank.length()), '0') + commRank;
    filepath.append("_");
    filepath.append(commSize);
    filepath.append("-");
    filepath.append(commRank);
  }

  // Clear local fieldset
  fset.clear();

  // Create local fieldset
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    atlas::Field field = fspace.createField<double>(
      atlas::option::name(vars[jvar]) | atlas::option::levels(variableSizes[jvar]));
    fset.add(field);
  }

  // Initialize local fieldset
  for (auto & field : fset) {
    auto view = atlas::array::make_view<double, 2>(field);
    view.assign(0.0);
  }

  // NetCDF IDs
  int ncid, retval, var_id[vars.size()];

  if (oneFilePerTask) {
    // Case 1: one file per MPI task

    // NetCDF file path
    std::string ncfilepath = filepath;
    ncfilepath.append(".nc");
    oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

    // Open NetCDF file
    if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

    // Get variables
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      if ((retval = nc_inq_varid(ncid, vars[jvar].c_str(), &var_id[jvar]))) ERR(retval);
    }

    // Get number of nodes
    size_t nb_nodes = 0;
    const auto ghostView = atlas::array::make_view<int, 1>(fspace.ghost());
    for (atlas::idx_t jnode = 0; jnode < fset.field(vars[0]).shape(0); ++jnode) {
      if (ghostView(jnode) == 0) ++nb_nodes;
    }

    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      // Read data
      double zvar[nb_nodes][variableSizes[jvar]];
      if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);

      // Copy data
      auto varView = atlas::array::make_view<double, 2>(fset[vars[jvar]]);
      size_t inode = 0;
      for (atlas::idx_t jnode = 0; jnode < fset.field(vars[jvar]).shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          for (size_t k = 0; k < variableSizes[jvar]; ++k) {
            varView(jnode, k) = zvar[inode][k];
          }
          ++inode;
        }
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval);

    if (fspace.type() != "PointCloud") {
      // Exchange halo
      fset.haloExchange();
    }
  } else {
    // Case 2: one file for all MPI tasks

    // Global data
    atlas::FieldSet globalData;
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      atlas::Field field = fspace.createField<double>(
        atlas::option::name(vars[jvar])
        | atlas::option::levels(variableSizes[jvar]) | atlas::option::global());
      globalData.add(field);
    }

    // NetCDF input
    if (fspace.type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(fspace);

      if (comm.rank() == 0) {
        // Get grid
        atlas::StructuredGrid grid = fs.grid();

        // Get sizes
        atlas::idx_t nx = grid.nxmax();
        atlas::idx_t ny = grid.ny();

        // NetCDF IDs
        int ncid, retval, var_id[vars.size()];

        // NetCDF file path
        std::string ncfilepath = filepath;
        ncfilepath.append(".");
        ncfilepath.append(config.getString("netcdf extension", "nc"));
        oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

        // Open NetCDF file
        if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

        // Get variables
        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          if ((retval = nc_inq_varid(ncid, vars[jvar].c_str(), &var_id[jvar]))) ERR(retval);
        }

        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          // Read data
          double zvar[variableSizes[jvar]][ny][nx];
          if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0][0]))) ERR(retval);

          // Copy data
          auto varView = atlas::array::make_view<double, 2>(globalData[vars[jvar]]);
          for (size_t k = 0; k < variableSizes[jvar]; ++k) {
            for (atlas::idx_t j = 0; j < ny; ++j) {
              for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
                atlas::gidx_t gidx = grid.index(i, ny-1-j);
                varView(gidx, k) = zvar[k][j][i];
              }
            }
          }
        }

        // Close file
        if ((retval = nc_close(ncid))) ERR(retval);
      }
    } else if (fspace.type() == "NodeColumns") {
      // NodeColumns
      atlas::idx_t nb_nodes;

      // CubedSphere
      atlas::functionspace::CubedSphereNodeColumns fs(fspace);

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
      /* TODO(??): have to assume the NodeColumns is CubedSphere here, cannot differentiate with
         other NodeColumns from what is in fspace.
      // Other NodeColumns
      atlas::functionspace::NodeColumns fs(fspace);

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
      */

      if (comm.rank() == 0) {
        // NetCDF IDs
        int ncid, retval, var_id[vars.size()];

        // NetCDF file path
        std::string ncfilepath = filepath;
        ncfilepath.append(".nc");
        oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

        // Open NetCDF file
        if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

        // Get variables
        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          if ((retval = nc_inq_varid(ncid, vars[jvar].c_str(), &var_id[jvar]))) ERR(retval);
        }

        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          // Read data
          double zvar[nb_nodes][variableSizes[jvar]];
          if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);

          // Copy data
          auto varView = atlas::array::make_view<double, 2>(globalData[vars[jvar]]);
          for (size_t k = 0; k < variableSizes[jvar]; ++k) {
            for (atlas::idx_t i = 0; i < nb_nodes; ++i) {
              varView(i, k) = zvar[i][k];
            }
          }
        }

        // Close file
        if ((retval = nc_close(ncid))) ERR(retval);
      }
    } else {
      ABORT(fspace.type() + " function space not supported yet");
    }

    // Scatter data from main processor
    if (fspace.type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(fspace);
      fs.scatter(globalData, fset);
    } else if (fspace.type() == "NodeColumns") {
      // NodeColumns

      // CubedSphere
      atlas::functionspace::CubedSphereNodeColumns fs(fspace);
      fs.scatter(globalData, fset);
    }

    // Exchange halo
    fset.haloExchange();
  }
}

// -----------------------------------------------------------------------------

void writeFieldSet(const eckit::mpi::Comm & comm,
                   const eckit::Configuration & config,
                   const atlas::FieldSet & fset) {
  // Define variables vector from fset
  std::vector<std::string> vars = fset.field_names();

  // Get function space
  atlas::FunctionSpace fspace = fset.field(vars[0]).functionspace();

  // Options with one file per MPI task
  const bool oneFilePerTask = config.getBool("one file per task", false);
  ASSERT(oneFilePerTask || (fspace.type() != "PointCloud"));

  // Build filepath
  std::string filepath = config.getString("filepath");
  if (config.has("member")) {
    std::ostringstream out;
    out << std::setfill('0') << std::setw(6) << config.getInt("member");
    filepath.append("_");
    filepath.append(out.str());
  }
  if (oneFilePerTask) {
    size_t zpad = 6;
    std::string commSize = std::to_string(comm.size());
    commSize = std::string(zpad - std::min(zpad, commSize.length()), '0') + commSize;
    std::string commRank = std::to_string(comm.rank()+1);
    commRank = std::string(zpad - std::min(zpad, commRank.length()), '0') + commRank;
    filepath.append("_");
    filepath.append(commSize);
    filepath.append("-");
    filepath.append(commRank);
  }

  // Missing value
  const double msvalr = util::missingValue<double>();

  // NetCDF IDs
  int retval, ncid, nx_id, ny_id, nb_nodes_id, nz_id[vars.size()],
    d1D_id[1], d2D_id[2], d3D_id[3], lon_id, lat_id, var_id[vars.size()];

  if (oneFilePerTask) {
    // Case 1: one file per MPI task

    // NetCDF file path
    std::string ncfilepath = filepath;
    ncfilepath.append(".nc");
    oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

    // Get number of nodes
    size_t nb_nodes = 0;
    const auto ghostView = atlas::array::make_view<int, 1>(fspace.ghost());
    for (atlas::idx_t jnode = 0; jnode < fset.field(vars[0]).shape(0); ++jnode) {
      if (ghostView(jnode) == 0) ++nb_nodes;
    }

    // Definition mode

    // Create NetCDF file
    if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

    // Create horizontal dimension
    if ((retval = nc_def_dim(ncid, "nb_nodes", nb_nodes, &nb_nodes_id))) ERR(retval);

    // Dimensions arrays, horizontal part
    d1D_id[0] = nb_nodes_id;
    d2D_id[0] = nb_nodes_id;

    // Define coordinates
    if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 1, d1D_id, &lon_id))) ERR(retval);
    if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 1, d1D_id, &lat_id))) ERR(retval);

    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      // Create vertical dimension
      std::string nzName = "nz_" + vars[jvar];
      if ((retval = nc_def_dim(ncid, nzName.c_str(), fset.field(vars[jvar]).levels(),
        &nz_id[jvar]))) ERR(retval);

      // Dimensions array, vertical part
      d2D_id[1] = nz_id[jvar];

      // Define variable
      if ((retval = nc_def_var(ncid, vars[jvar].c_str(), NC_DOUBLE, 2, d2D_id,
        &var_id[jvar]))) ERR(retval);

      // Add missing value metadata
      if ((retval = nc_put_att_double(ncid, var_id[jvar], "_FillValue", NC_DOUBLE, 1,
        &msvalr))) ERR(retval);
    }

    // End definition mode
    if ((retval = nc_enddef(ncid))) ERR(retval);

    // Data mode

    // Copy coordinates
    const auto lonlatView = atlas::array::make_view<double, 2>(fspace.lonlat());
    double zlon[nb_nodes][1];
    double zlat[nb_nodes][1];
    size_t inode = 0;
    for (atlas::idx_t jnode = 0; jnode < fset.field(vars[0]).shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        zlon[inode][0] = lonlatView(jnode, 0);
        zlat[inode][0] = lonlatView(jnode, 1);
        ++inode;
      }
    }

    // Write coordinates
    if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
    if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      // Copy data
      const auto varView = atlas::array::make_view<double, 2>(fset.field(vars[jvar]));
      double zvar[nb_nodes][fset.field(vars[jvar]).levels()];
      inode = 0;
      for (atlas::idx_t jnode = 0; jnode < fset.field(vars[0]).shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          for (atlas::idx_t k = 0; k < fset.field(vars[jvar]).levels(); ++k) {
            zvar[inode][k] = varView(jnode, k);
          }
          ++inode;
        }
      }

      // Write data
      if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval);
  } else {
    // Case 2: one file for all MPI tasks

    // Prepare local coordinates and data
    atlas::FieldSet localData;
    atlas::Field lonLocal = fspace.createField<double>(atlas::option::name("lon"));
    localData.add(lonLocal);
    atlas::Field latLocal = fspace.createField<double>(atlas::option::name("lat"));
    localData.add(latLocal);
    const auto lonlatView = atlas::array::make_view<double, 2>(fspace.lonlat());
    auto lonViewLocal = atlas::array::make_view<double, 1>(localData.field("lon"));
    auto latViewLocal = atlas::array::make_view<double, 1>(localData.field("lat"));
    for (atlas::idx_t jnode = 0; jnode < fspace.lonlat().shape(0); ++jnode) {
       lonViewLocal(jnode) = lonlatView(jnode, 0);
       latViewLocal(jnode) = lonlatView(jnode, 1);
    }
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      localData.add(fset.field(vars[jvar]));
    }

    // Prepare global coordinates and data
    atlas::FieldSet globalData;
    atlas::Field lonGlobal = fspace.createField<double>(
      atlas::option::name("lon") | atlas::option::global());
    globalData.add(lonGlobal);
    atlas::Field latGlobal = fspace.createField<double>(
      atlas::option::name("lat") | atlas::option::global());
    globalData.add(latGlobal);
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      atlas::Field globalField = fspace.createField<double>(atlas::option::name(vars[jvar]) |
        atlas::option::levels(fset.field(vars[jvar]).levels()) | atlas::option::global());
      globalData.add(globalField);
    }

    if (fspace.type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(fspace);

      // Gather coordinates and data on main processor
      fs.gather(localData, globalData);

      if (comm.rank() == 0) {
        // NetCDF file path
        std::string ncfilepath = filepath;
        ncfilepath.append(".nc");
        oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

        // Get grid
        atlas::StructuredGrid grid = fs.grid();

        // Get sizes
        atlas::idx_t nx = grid.nxmax();
        atlas::idx_t ny = grid.ny();

        // Definition mode

        // Create NetCDF file
        if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

        // Create dimensions
        if ((retval = nc_def_dim(ncid, "nx", nx, &nx_id))) ERR(retval);
        if ((retval = nc_def_dim(ncid, "ny", ny, &ny_id))) ERR(retval);

        // Dimensions arrays, horizontal part
        d2D_id[0] = ny_id;
        d2D_id[1] = nx_id;
        d3D_id[1] = ny_id;
        d3D_id[2] = nx_id;

        // Define coordinates
        if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 2, d2D_id, &lon_id))) ERR(retval);
        if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 2, d2D_id, &lat_id))) ERR(retval);
        if ((retval = nc_put_att_double(ncid, lon_id, "_FillValue", NC_DOUBLE, 1, &msvalr)))
          ERR(retval);
        if ((retval = nc_put_att_double(ncid, lat_id, "_FillValue", NC_DOUBLE, 1, &msvalr)))
          ERR(retval);

        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          // Create vertical dimension
          std::string nzName = "nz_" + vars[jvar];
          if ((retval = nc_def_dim(ncid, nzName.c_str(), fset.field(vars[jvar]).levels(),
            &nz_id[jvar]))) ERR(retval);

          // Dimensions array, vertical part
          d3D_id[0] = nz_id[jvar];

          // Define variable
          if ((retval = nc_def_var(ncid, vars[jvar].c_str(), NC_DOUBLE, 3, d3D_id,
            &var_id[jvar]))) ERR(retval);

          // Add missing value metadata
          if ((retval = nc_put_att_double(ncid, var_id[jvar], "_FillValue", NC_DOUBLE, 1,
            &msvalr))) ERR(retval);
        }

        // End definition mode
        if ((retval = nc_enddef(ncid))) ERR(retval);

        // Data mode

        // Copy coordinates
        auto lonViewGlobal = atlas::array::make_view<double, 1>(globalData.field("lon"));
        auto latViewGlobal = atlas::array::make_view<double, 1>(globalData.field("lat"));
        double zlon[ny][nx];
        double zlat[ny][nx];
        for (atlas::idx_t j = 0; j < ny; ++j) {
          for (atlas::idx_t i = 0; i < nx; ++i) {
            zlon[j][i] = msvalr;
            zlat[j][i] = msvalr;
          }
          for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
            atlas::gidx_t gidx = grid.index(i, ny-1-j);
            zlon[j][i] = lonViewGlobal(gidx);
            zlat[j][i] = latViewGlobal(gidx);
          }
        }

        // Write coordinates
        if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
        if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          // Copy data
          auto varView = atlas::array::make_view<double, 2>(globalData[vars[jvar]]);
          double zvar[fset.field(vars[jvar]).levels()][ny][nx];
          for (atlas::idx_t k = 0; k < fset.field(vars[jvar]).levels(); ++k) {
            for (atlas::idx_t j = 0; j < ny; ++j) {
              for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
                atlas::gidx_t gidx = grid.index(i, ny-1-j);
                zvar[k][j][i] = varView(gidx, k);
              }
            }
          }

          // Write data
          if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0][0]))) ERR(retval);
        }

        // Close file
        if ((retval = nc_close(ncid))) ERR(retval);
      }
    } else if (fspace.type() == "NodeColumns") {
      /* TODO(??): have to assume the NodeColumns is CubedSphere here, cannot differentiate with
         other NodeColumns from what is in fspace.
      */

      // CubedSphere
      atlas::functionspace::CubedSphereNodeColumns fs(fspace);

      // Gather coordinates and data on main processor
      fs.gather(localData, globalData);

      if (comm.rank() == 0) {
        // NetCDF file path
        std::string ncfilepath = filepath;
        ncfilepath.append(".nc");
        oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

        // Get global number of nodes
        atlas::idx_t nb_nodes = fs.nb_nodes_global();

        // Definition mode

        // Create NetCDF file
        if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

        // Create dimensions
        if ((retval = nc_def_dim(ncid, "nb_nodes", nb_nodes, &nb_nodes_id))) ERR(retval);

        // Dimensions arrays, horizontal part
        d1D_id[0] = nb_nodes_id;
        d2D_id[0] = nb_nodes_id;

        // Define coordinates
        if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 1, d1D_id, &lon_id))) ERR(retval);
        if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 1, d1D_id, &lat_id))) ERR(retval);

        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          // Create vertical dimension
          std::string nzName = "nz_" + vars[jvar];
          if ((retval = nc_def_dim(ncid, nzName.c_str(), fset.field(vars[jvar]).levels(),
            &nz_id[jvar]))) ERR(retval);

          // Dimensions array, vertical part
          d2D_id[1] = nz_id[jvar];

          // Define variable
          if ((retval = nc_def_var(ncid, vars[jvar].c_str(), NC_DOUBLE, 2, d2D_id,
            &var_id[jvar]))) ERR(retval);

          // Add missing value metadata
          if ((retval = nc_put_att_double(ncid, var_id[jvar], "_FillValue", NC_DOUBLE, 1,
            &msvalr))) ERR(retval);
        }

        // End definition mode
        if ((retval = nc_enddef(ncid))) ERR(retval);

        // Data mode

        // Copy coordinates
        auto lonViewGlobal = atlas::array::make_view<double, 1>(globalData.field("lon"));
        auto latViewGlobal = atlas::array::make_view<double, 1>(globalData.field("lat"));
        double zlon[nb_nodes][1];
        double zlat[nb_nodes][1];
        for (atlas::idx_t i = 0; i < nb_nodes; ++i) {
          zlon[i][0] = lonViewGlobal(i);
          zlat[i][0] = latViewGlobal(i);
        }

        // Write coordinates
        if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
        if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

        for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
          // Copy data
          auto varView = atlas::array::make_view<double, 2>(globalData[vars[jvar]]);
          double zvar[nb_nodes][fset.field(vars[jvar]).levels()];
          for (atlas::idx_t k = 0; k < fset.field(vars[jvar]).levels(); ++k) {
            for (atlas::idx_t i = 0; i < nb_nodes; ++i) {
              zvar[i][k] = varView(i, k);
            }
          }

          // Write data
          if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);
        }

        // Close file
        if ((retval = nc_close(ncid))) ERR(retval);
      }
    } else {
      ABORT(fspace.type() + " function space not supported yet");
    }
  }
}

// -----------------------------------------------------------------------------

atlas::FieldSet createRandomFieldSet(const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspace,
                                     const oops::Variables & vars) {
  std::vector<size_t> variableSizes;
  for (const std::string & var : vars.variables()) {
    variableSizes.push_back(vars.getLevels(var));
  }
  return createRandomFieldSet(comm, fspace, variableSizes, vars.variables());
}

// -----------------------------------------------------------------------------

atlas::FieldSet createSmoothFieldSet(const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspace,
                                     const oops::Variables & vars) {
  std::vector<size_t> variableSizes;
  for (const std::string & var : vars.variables()) {
    variableSizes.push_back(vars.getLevels(var));
  }
  return createSmoothFieldSet(comm, fspace, variableSizes, vars.variables());
}

// -----------------------------------------------------------------------------

void readFieldSet(const eckit::mpi::Comm & comm,
                  const atlas::FunctionSpace & fspace ,
                  const oops::Variables & vars,
                  const eckit::Configuration & config,
                  atlas::FieldSet & fset) {
  std::vector<size_t> variableSizes;
  for (const std::string & var : vars.variables()) {
    variableSizes.push_back(vars.getLevels(var));
  }
  readFieldSet(comm, fspace, variableSizes, vars.variables(), config, fset);
}



// -----------------------------------------------------------------------------

}  // namespace util
