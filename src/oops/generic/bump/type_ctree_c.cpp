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
    // Insert points
    int j=0;
    for(int i=0;i<n;i++) {
        // Check mask
        if(mask[i]==1) {
            if(j==0) {
                // Create tree with one point
                vector<CoverTreePoint> root;
                root.push_back(CoverTreePoint(i,lon[i],lat[i]));
                tree = new CoverTree<CoverTreePoint>(root);
                j++;
            } else {
                // Find the nearest neighbor
                vector<CoverTreePoint> neighbors(tree->kNearestNeighbors(CoverTreePoint(-999,lon[i],lat[i]),1));
                if(neighbors[0].getDist()>0.0) {
                   // Insert point
                   tree->insert(CoverTreePoint(i,lon[i],lat[i]));
                   j++;
                }
            }
        }
    }
    return;
}

// Destructor
cTree::~cTree(){}

// Find redundant
void cTree::find_redundant(int n, double lon[], double lat[], int redundant[]) {
    // Insert points
    int j=1;
    for(int i=1;i<n;i++) {
       // Find the nearest neighbor
       vector<CoverTreePoint> neighbors(tree->kNearestNeighbors(CoverTreePoint(-999,lon[i],lat[i]),1));
       if(neighbors[0].getDist()>0.0) {
          // Insert point
          tree->insert(CoverTreePoint(j,lon[i],lat[i]));
          j++;
       } else {
          // Mark as redundant
          redundant[i] = neighbors[0].getIndex()+1;
       }
    }
    return;
}

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
