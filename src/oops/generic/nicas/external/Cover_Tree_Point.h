//----------------------------------------------------------------------
// CoverTreePoint
// Modified by Benjamin Menetrier for nicas
//----------------------------------------------------------------------
#ifndef _CTREE_POINT_H
#define _CTREE_POINT_H
#include <cmath>

/**
 * A simple point class containing a vector of doubles and a single char name.
 */
class CoverTreePoint
{
private:
    int _index;
    double _lon;
    double _sinlat;
    double _coslat;
    double _dist;
public:
    CoverTreePoint(int index, double lon, double lat) : _index(index), _lon(lon), _sinlat(sin(lat)), _coslat(cos(lat)) {}
    double distance(const CoverTreePoint& p) const;
    const int& getIndex() const;
    const double& getLon() const;
    const double& getSinlat() const;
    const double& getCoslat() const;
    const double& getDist() const;
    bool operator==(const CoverTreePoint&) const;
    void setDist(double dist);
};

#endif // _CTREE_POINT_H

