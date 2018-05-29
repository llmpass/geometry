#ifndef _CIRCLE_H
#define _CIRCLE_H

#include "Point3.h"
#include "Plane.h"
#include <iostream>
using namespace std;

namespace geometry {
  class Circle {
    public:
    Point3 o;
    float r;
    Plane p;
    Circle() {}
    Circle(Point3 o, float r, Plane p) {
      this->o = o;
      this->r = r;
      this->p = p;
    }
    ~Circle() {}
    /*
     * true:  incircle or oncircle
     * false: outcircle
     */
    bool inCircle(Point3 p) {
      if (p.dis(o)<=r) return true;
      return false;
    }
    void print() {
      cout<<"Circle o=";
      o.print();
      cout<<" r="<<r<<endl<<"on ";
      p.print();
    }
  };
}
#endif
