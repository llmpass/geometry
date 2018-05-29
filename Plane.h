#ifndef _PLANE_H
#define _PLANE_H

#include "Point3.h"
#include "../util/Vec3d.h"
using namespace util;

namespace geometry {
  /**
   * A plane is defined by a normal vector n and a 3d point p.
   * n /dot (x-p) = 0.
   */
  class Plane {
    public:
    Vec3d n;
    Point3 p;
    Plane() {}
    Plane(Vec3d n, Point3 p) {
      this->n = n;
      this->p = p;
    }
    Plane(Point3 p, Vec3d n) {
      this->n = n;
      this->p = p;
    }
    ~Plane() {}
    void print() {
      cout<<"Plane pass p=";
      p.print();
      cout<<"with normal "<<n.x<<","<<n.y<<","<<n.z;
    }
  };
}
#endif
