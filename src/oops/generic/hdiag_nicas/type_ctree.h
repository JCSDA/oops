//---------------------------------------------------------------------
/// Purpose: cover tree class header
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
extern "C"
{
    class cTree;
    typedef cTree CTREE;

    // Constructor
    CTREE* create_ctree(int n, double lon[], double lat[], int mask[]);

    // Destructor
    void delete_ctree(CTREE* ctree);

    // Find nearest neighbors
    void find_nearest_neighbors(const CTREE* ctree, double lon, double lat, int nn, int nn_index[], double nn_dist[]);
}
