#ifndef _POINT3_H
#define _POINT3_H

#include "../util/Vec3d.h"
using namespace util;

namespace geometry {
  class Point3 {
    public:
    double _x, _y, _z;
    double _nx, _ny, _nz;
    Point3() {}
    Point3(double x, double y, double z) {
      _x = x; _y = y; _z = z;
      _nx = 0; _ny = 0; _nz = 0;
    }
    Point3(double x, double y, double z, double nx, double ny, double nz) {
      _x = x; _y = y; _z = z;
      _nx = nx; _ny = ny; _nz = nz;
    }
    ~Point3() {}

    void set(double x, double y, double z) {
      _x = x; _y = y; _z = z;
    }
    void set(double x, double y, double z, double nx, double ny, double nz) {
      _x = x; _y = y; _z = z;
      _nx = nx; _ny = ny; _nz = nz;
    }
    void copy(Point3& p) {
      this->_x = p._x;
      this->_y = p._y;
      this->_z = p._z;
      this->_nx = p._nx;
      this->_ny = p._ny;
      this->_nz = p._nz;
    }
    Vec3d toVec() {
      return Vec3d(_x,_y,_z);
    }
    Vec3d normal() {
      return Vec3d(_nx,_ny,_nz);
    }
    double dis(Point3 p) {
      double d = (_x-p._x)*(_x-p._x)+(_y-p._y)*(_y-p._y)+(_z-p._z)*(_z-p._z);
      return sqrt(d);
    } 
    void print() {
      cout<<"("<<_x<<","<<_y<<","<<_z<<")  ";
    }
  };
  inline Point3 operator-(const Point3& p, const Vec3d &v) { 
    return Point3(p._x-v.x,p._y-v.y,p._z-v.z); 
  };
}

#endif
