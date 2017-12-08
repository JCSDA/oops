//---------------------------------------------------------------------
/// Purpose: cover tree class API
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "type_ctree.h"
#include "type_ctree.hpp"

using namespace std;

// Constructor
CTREE* create_ctree(int n, double lon[], double lat[], int mask[]){
    return new cTree(n, lon, lat, mask);
}

// Destructor
void delete_ctree(CTREE* ctree){
    delete ctree;
}

// Find nearest neighbors
void find_nearest_neighbors(const CTREE* ctree, double lon, double lat, int nn, int nn_index[], double nn_dist[]) {
    ctree->find_nearest_neighbors(lon, lat, nn, nn_index, nn_dist);
}
