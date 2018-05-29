#ifndef _SEGEMENT_H
#define _SEGEMENT_H

#include "Point3.h"
using namespace util;

namespace geometry {
  class Segment {
    public:
    Point3 s, e;
    float t; // used in the parametric form: x(t) = (1-t)s+te, t\in[0,1]
    Vec3d v;
    Segment() {}
    Segment(Point3 p1, Point3 p2) {
      s.copy(p1); e.copy(p2);
      v = Vec3d(e._x-s._x,e._y-s._y,e._z-s._z);
    }
    Segment(const Segment& seg):s(seg.s),e(seg.e) {
      v = Vec3d(e._x-s._x,e._y-s._y,e._z-s._z);
    }
    ~Segment() {}

    void set(Point3 p1, Point3 p2) {
      s = p1; e = p2;
      v = Vec3d(e._x-s._x,e._y-s._y,e._z-s._z);
    }
    void print() {
      cout<<"Segment ";
      s.print(); cout<<"--  "; e.print(); cout<<endl;
    }

    double length() {
      v = Vec3d(e._x-s._x,e._y-s._y,e._z-s._z);
      return Length(v);
    }
    Point3 center() {
      float x = s._x+e._x;
      float y = s._y+e._y;
      float z = s._z+e._z;
      return Point3(x*0.5,y*0.5,z*0.5);
    }

    static Segment nullSeg() {
      return Segment(Point3(0,0,0), Point3(0,0,0));
    }
    void setNull() {
      set(Point3(0,0,0), Point3(0,0,0));
    }
    bool isNull() {
      if (length()==0) return true;
      else return false;
    }
  };
}

#endif
