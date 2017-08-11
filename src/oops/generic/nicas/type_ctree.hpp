//---------------------------------------------------------------------
/// Purpose: cover tree class header
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "external/Cover_Tree.h"
#include "external/Cover_Tree_Point.h"

class cTree {
    public:
        // Constructor
        cTree(int n, double lon[], double lat[], int mask[]);

        // Destructor
        ~cTree();

        // Find nearest neighbors
        void find_nearest_neighbors(double lon, double lat, int nn, int nn_index[], double nn_dist[]) const;
    private:
        // Cover tree
        CoverTree<CoverTreePoint> *tree;
};
