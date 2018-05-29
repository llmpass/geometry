#ifndef _POINT3F_H
#define _POINT3F_H

#include "../util/Vec3f.h"
using namespace util;

namespace geometry {
  class Point3f {
    public:
    float _x, _y, _z;
    float _nx, _ny, _nz;
    Point3f() {}
    Point3f(float x, float y, float z) {
      _x = x; _y = y; _z = z;
      _nx = 0; _ny = 0; _nz = 0;
    }
    Point3f(float x, float y, float z, float nx, float ny, float nz) {
      _x = x; _y = y; _z = z;
      _nx = nx; _ny = ny; _nz = nz;
    }
    ~Point3f() {}

    void set(float x, float y, float z) {
      _x = x; _y = y; _z = z;
    }
    void set(float x, float y, float z, float nx, float ny, float nz) {
      _x = x; _y = y; _z = z;
      _nx = nx; _ny = ny; _nz = nz;
    }
    void copy(Point3f& p) {
      this->_x = p._x;
      this->_y = p._y;
      this->_z = p._z;
      this->_nx = p._nx;
      this->_ny = p._ny;
      this->_nz = p._nz;
    }
    Vec3f toVec() {
      return Vec3f(_x,_y,_z);
    }
    Vec3f normal() {
      return Vec3f(_nx,_ny,_nz);
    }
    float dis(Point3f p) {
      float d = (_x-p._x)*(_x-p._x)+(_y-p._y)*(_y-p._y)+(_z-p._z)*(_z-p._z);
      return sqrt(d);
    }
    void print() {
      cout<<"("<<_x<<","<<_y<<","<<_z<<")  ";
    }
  };

  inline Point3f operator+(const Point3f &a, const Vec3f &b)
  { return Point3f(a._x+b.x, a._y+b.y, a._z+b.z); };
  inline Point3f operator-(const Point3f &a, const Vec3f &b)
  { return Point3f(a._x-b.x, a._y-b.y, a._z-b.z); };
  inline Point3f operator-(const Point3f &v)
  { return Point3f(-v._x,-v._y,-v._z); };
}

#endif

