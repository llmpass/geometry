#ifndef _POLYGON_H
#define _POLYGON_H

namespace geometry {
  /**
   * A polygon is expressed as a counterclockwise sequence of vertices (point3).
   * The edges (segements) are formed by p0p1, p1p2, ... pn-1pn, pnp0.
   **/
  class Polygon {
    public:
    int n; // number of vertices (segments)
    Point3* p; // point sequence
    Segment* s;
    Polygon() {n=0; p=0; s=0;}
    Polygon(Point3* pp, int n) {
      this->n = n;
      p = new Point3[n];
      s = new Segment[n];
      for (int i=0; i<n; ++i) 
        p[i].set(pp[i]._x,pp[i]._y,pp[i]._z);
      for (int i=0; i<n; ++i) 
        s[i].set(p[i],p[(i+1)%n]);
    }
    Polygon(const Polygon& p1):n(p1.n) {
      p = new Point3[n];
      s = new Segment[n];
      for (int i=0; i<n; ++i) p[i].set(p1.p[i]._x,p1.p[i]._y,p1.p[i]._z);
      for (int i=0; i<n; ++i) s[i].set(p[i],p[(i+1)%n]);
    }
    ~Polygon() {
      if (n>0) {delete [] p; delete [] s;}
    }
    void set(Point3* pp, int n1) {
      if (n>0) {delete [] p; delete [] s;}
      n = n1;
      p = new Point3[n];
      s = new Segment[n];
      for (int i=0; i<n; ++i) 
        p[i].set(pp[i]._x,pp[i]._y,pp[i]._z);
      for (int i=0; i<n; ++i) 
        s[i].set(p[i],p[(i+1)%n]);
    }
    void copy(Polygon pl) {
      if (n>0) {delete [] p; delete [] s;}
      n = pl.n;
      p = new Point3[n];
      s = new Segment[n];
      for (int i=0; i<n; ++i) 
        p[i].set(pl.p[i]._x,pl.p[i]._y,pl.p[i]._z);
      for (int i=0; i<n; ++i) s[i].set(p[i],p[(i+1)%n]);
    }
    Polygon operator=(const Polygon& p1) const {
      Polygon p2;
      p2.copy(p1);
      return p2;
    }
    void setNull() {
      n=0;
    }
    static Polygon nullPoly() {
      Polygon p;
      p.n=0;
      return p;
    }
    bool isNull() {
      if (n==0) return true;
      else return false;
    }
    void print() {
      cout<<"Polygon with "<<n<<" vertices:"<<endl;
      for (int i=0; i<n; ++i) p[i].print();
    }
    // mass center
    Point3 center() {
      float x=0, y=0, z=0, n1=1.0f/n;
      for (int i=0; i<n; ++i) {
        x+=p[i]._x; 
        y+=p[i]._y; 
        z+=p[i]._z;
      }
      return Point3(x*n1,y*n1,z*n1);
    }
    // reverse the order of the polygon
    // 0,1,2...,n-1 ---> n-1,n-2,...,0
    void reverse() {
      Point3* pp = new Point3[n];
      for (int i=0; i<n; ++i) pp[i] = p[n-1-i];
      set(pp,n);
      delete [] pp;
    }
  };
}

#endif
