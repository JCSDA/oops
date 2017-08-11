//----------------------------------------------------------------------
// CoverTreePoint
// Modified by Benjamin Menetrier for nicas
//----------------------------------------------------------------------
#include "Cover_Tree_Point.h"
#include <iostream>
#include <cmath>

using namespace std;

double CoverTreePoint::distance(const CoverTreePoint& p) const {
    const double& lon=p.getLon();
    const double& sinlat=p.getSinlat();
    const double& coslat=p.getCoslat();

    // Great-circle distance using Vincenty formula on the sphere
    double dist = atan2(sqrt(pow(coslat*sin(lon-_lon),2)+pow(_coslat*sinlat-_sinlat*coslat*cos(lon-_lon),2)),_sinlat*sinlat+_coslat*coslat*cos(lon-_lon));
    return dist;
}

const int& CoverTreePoint::getIndex() const {
    return _index;
}

const double& CoverTreePoint::getLon() const {
    return _lon;
}

const double& CoverTreePoint::getSinlat() const {
    return _sinlat;
}

const double& CoverTreePoint::getCoslat() const {
    return _coslat;
}

const double& CoverTreePoint::getDist() const {
    return _dist;
}

bool CoverTreePoint::operator==(const CoverTreePoint& p) const {
    return (_lon==p.getLon() && _sinlat==p.getSinlat() && _coslat==p.getCoslat());
}

void CoverTreePoint::setDist(double dist) {
    _dist = dist;
    return;
}
