//---------------------------------------------------------------------
/// Purpose: cover tree class implementation
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "type_ctree.hpp"
#include <ostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Constructor
cTree::cTree(int n, double lon[], double lat[], int mask[]) {
    // Initialize tree: maximum distance is half the sphere circonference (sphere of radius 1, easier)
    tree = new CoverTree<CoverTreePoint>(M_PI);

    for(int i=0;i<n;i++) {
        // Check mask
        if(mask[i]==1) {
            // Insert point
            tree->insert(CoverTreePoint(i,lon[i],lat[i]));
        }
    }
}

// Destructor
cTree::~cTree(){}

// Find nearest neighbors
void cTree::find_nearest_neighbors(double lon, double lat, int nn, int nn_index[], double nn_dist[]) const{
    // Find the n nearest neighbors
    vector<CoverTreePoint> neighbors(tree->kNearestNeighbors(CoverTreePoint(-999,lon,lat),nn));
    for(int i=0;i<nn;i++) {
        // Copy neighbors index and disance
        nn_index[i]=neighbors[i].getIndex()+1;
        nn_dist[i]=neighbors[i].getDist();
    }
    return;
}
