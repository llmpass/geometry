#include "GeoAlg.h"
using namespace geometry;

  double geometry::leftToSegment(Point3& p, Segment& s, Vec3d& n) {
    // first define a plane using the input segment and input point
    Point3 p0 = s.s, p1 = s.e;
    Vec3d v01 = s.v;
    // a vector pointing from p0 to p.
    Vec3d v02(p._x-p0._x,p._y-p0._y,p._z-p0._z);
    // local plane
    Plane pl(n,p0);
    // project p1 and p2 to the local plane
    Vec3d v0u = Cross(v01,n); 
    // obtain 2d coords in the system defined by v01 x v0u, where
    // v0u is perpendicular to both v01 and n.
    double x0=0, y0=0;
    double x1=0, y1=Length(v01);
    double x2=Dot(v02,v0u), y2=Dot(v02,v01);
    double acx = x0 - x2;
    double bcx = x1 - x2;
    double acy = y0 - y2;
    double bcy = y1 - y2;
    return acx*bcy-acy*bcx;
  }

  bool geometry::intersect(Plane& p, Segment& s, double& t, Point3& interP) {
    double ndv = Dot(p.n,s.v);
    if (ndv==0) return false;
    Vec3d pms(p.p._x-s.s._x,p.p._y-s.s._y,p.p._z-s.s._z);
    double ndpms = Dot(p.n,pms);
    t = ndpms/ndv;
    if (t<0 || t>1) return false;
    interP.set(s.s._x+t*s.v.x,s.s._y+t*s.v.y,s.s._z+t*s.v.z);
    return true;
  }

  /**
   * Intersection between a plane and a triangle (stored in a polygon).
   * If intersect, return true, resulting segment stored in s.
   * Otherwise, return false.
   */
  bool geometry::intersect(Plane& pl, Polygon& poly, Segment& s) {
    double ta, tb, t;
    Point3 pa, pb, p;
    int i1, i2, i, n = poly.n;
    bool find1 = false, find2 = false;
    bool asInterPnt[n];
    for (i=0; i<n; ++i) asInterPnt[i] = false;
    // loop over vertices of the polygon and find out which two segments are
    // intersected with the plane.
    for (i=0; i<n; ++i) {
      if (intersect(pl,poly.s[i],t,p)) {
        // degenerate
        if (t==0) {
          if (asInterPnt[i]) continue;
          t=0.01f; 
          p.set(poly.s[i].s._x+t*poly.s[i].v.x, 
            poly.s[i].s._y+t*poly.s[i].v.y,
            poly.s[i].s._z+t*poly.s[i].v.z);
          asInterPnt[i] = true;
        }
        if (t==1) {
          if (asInterPnt[(i+1)%n]) continue;
          t=0.99f; 
          p.set(poly.s[i].s._x+t*poly.s[i].v.x,
            poly.s[i].s._y+t*poly.s[i].v.y,
            poly.s[i].s._z+t*poly.s[i].v.z);
          asInterPnt[(i+1)%n] = true;
        }
        if (find1) {pb.copy(p); i2 = i; find2 = true; break;}
        else {find1 = true; i1 = i; pa.copy(p);}
      }
    }
    if (find1 && find2) {
      s.set(pa,pb);
      return true;
    }
    else return false;
  }

  /**
   * Intersection between a segment and a circle.
   * The segment and the circle should be in the same plane.
   * If there is only one solution, we store this solution in p1 and t1.
   * Return the number of intersections.
   */
  int geometry::intersect(Circle& cir, Segment& s, double& t1, double& t2, 
    Point3& p1, Point3& p2) {
    Vec3d q = s.s.toVec()-cir.o.toVec();
    double a = Dot(s.v,s.v);
    double b = 2*Dot(q,s.v);
    double c = Dot(q,q)-cir.r*cir.r;
    double delta = b*b-4*a*c;
    if (delta<0) return 0;
    if (delta==0) {
      t1 = -b/(2*a);
      if (t1<1 || t1>0) return 0;
      p1.set(s.s._x+t1*s.v.x,s.s._y+t1*s.v.y,s.s._z+t1*s.v.z);
      return 1;
    }
    // two intersections: t1<t2
    t1 = (-b-sqrt(delta))/(2*a);
    t2 = (-b+sqrt(delta))/(2*a);
    if (t1<0 || t1>1) {
      t1 = t2;
      if (t1<0 || t1>1) return 0;
      p1.set(s.s._x+t1*s.v.x,s.s._y+t1*s.v.y,s.s._z+t1*s.v.z);
      return 1;
    }
    p1.set(s.s._x+t1*s.v.x,s.s._y+t1*s.v.y,s.s._z+t1*s.v.z);
    if (t2<0 || t2>1) return 1;
    p2.set(s.s._x+t2*s.v.x,s.s._y+t2*s.v.y,s.s._z+t2*s.v.z);
    return 2;
  }

  void geometry::clip(Segment s, Plane pl, Segment& s1, Segment& s2) {
    double t;
    Point3 p;
    if (intersect(pl,s,t,p)) {
      s1.set(s.s,p);
      s2.set(p,s.e);       
    }
    else {
      s1.set(s.s,s.e);
      s2.setNull();
    } 
  }

  /** 
   * Clip a polygon with a plane.
   * The vertex number of the polygon can only be increased by at most 1.
   * p1 and p2 are two output polygons.
   * If no intersection between the polygon and the plane, return false,
   * otherwise return true.
   */
  bool geometry::clip(Polygon poly, Plane pl, Polygon& p1, Polygon& p2) {
    double ta, tb, t;
    Point3 pa, pb, p;
    int i1, i2, i, n = poly.n;
    bool find1 = false, find2 = false;
    bool asInterPnt[n];
    for (i=0; i<n; ++i) asInterPnt[i] = false;
    // loop over vertices of the polygon and find out which two segments are
    // intersected with the plane.
    for (i=0; i<n; ++i) {
      if (intersect(pl,poly.s[i],t,p)) {
        // degenerate
        if (t==0) {
          if (asInterPnt[i]) continue;
          t=0.01f; 
          p.set(poly.s[i].s._x+t*poly.s[i].v.x, 
            poly.s[i].s._y+t*poly.s[i].v.y,
            poly.s[i].s._z+t*poly.s[i].v.z);
          asInterPnt[i] = true;
        }
        if (t==1) {
          if (asInterPnt[(i+1)%n]) continue;
          t=0.99f; 
          p.set(poly.s[i].s._x+t*poly.s[i].v.x,
            poly.s[i].s._y+t*poly.s[i].v.y,
            poly.s[i].s._z+t*poly.s[i].v.z);
          asInterPnt[(i+1)%n] = true;
        }
        if (find1) {pb.copy(p); i2 = i; find2 = true; break;}
        else {find1 = true; i1 = i; pa.copy(p);}
      }
    }
    if (find1 && find2) {
      // make sure i1 ia always smaller than i2
      if (i1>i2) {
        i = i1; i1 = i2; i2 = i;
        p.copy(pa); pa.copy(pb); pb.copy(p);
      }
      // form two new polygons p1 and p2
      // p1: pa pi1+1 ... pi2 pb (pa) i2-i1+2 vertices
      Point3* pp1 = new Point3[i2-i1+2];
      pp1[0].copy(pa);
      for (i=i1+1; i<=i2; ++i) pp1[i-i1].copy(poly.p[i]);
      pp1[i2-i1+1].copy(pb);
      p1.set(pp1,i2-i1+2);
      // p2: pb pi2+1 ... pi1 pa (pb) n-(i2-i1)+2 vertices
      Point3* pp2 = new Point3[n-(i2-i1)+2];
      pp2[0].copy(pb);
      for (i=1; i<n-(i2-i1)+1; ++i) pp2[i].copy(poly.p[(i+i2)%n]);
      pp2[n-(i2-i1)+1].copy(pa);
      p2.set(pp2,n-(i2-i1)+2);
      delete [] pp1; delete [] pp2;
      return true; 
    }
    else {
      // if there are no intersections, just set p1 as the input polygon and
      // p2 as a null polygon
      p1.copy(poly);
      p2.setNull();
      return false;
    }
  }

  /*
   * Check the position of a point respect to a plane.
   *  1: in the half space with the normal pointing
   *  0: on plane
   * -1: in the half space back to the normal pointing
   */
  double geometry::halfSpace(Point3 p, Plane pl) {
    Vec3d n  = pl.n;
    Vec3d px = p.toVec()-pl.p.toVec();
    Normalize(n); Normalize(px);
    return Dot(n,px);
  }

  /**
   * Check if a point appears in between two planes.
   * ==> this point should appear in different sides of the two planes.
   */
  bool geometry::inBetween(Point3 p, Plane p1, Plane p2) {
    double h1 = halfSpace(p,p1);
    double h2 = halfSpace(p,p2);
    if (h1*h2<0) return true;
    else return false;
  }

  Segment geometry::parallelClip(Segment s, Plane p1, Plane p2) {
    // check the normal directions of two planes
    Vec3d n1 = p1.n, n2 = p2.n;
    if (Dot(n1,n2)<0) n1 = -n1;
    // clip the segment by the first plane
    Segment s1, s2;
    clip(s,p1,s1,s2);
    Segment sa, sb, sc, sd;
    if (!s1.isNull()) clip(s1,p2,sa,sb);
    // check if sa or sb lie in between two planes
    if (!sa.isNull() && inBetween(sa.center(),p1,p2)) return sa; 
    if (!sb.isNull() && inBetween(sb.center(),p1,p2)) return sb; 
    if (!s2.isNull()) clip(s2, p2, sc, sd);
    else return Segment::nullSeg();
    // check if sc or sd lie in between two planes
    if (!sc.isNull() && inBetween(sc.center(),p1,p2)) return sc; 
    if (!sd.isNull() && inBetween(sd.center(),p1,p2)) return sd; 
    return Segment::nullSeg();
  }

  /**
   * Clip a polygon with 2 parallel planes.
   * Clip the polygon by these two planes sequentially. 
   * Check the resulting parts, return the polygon that lies between two
   * planes.
   */
  Polygon geometry::parallelClip(Polygon in, Plane p1, Plane p2) {
    // check the normal directions of two planes
    Vec3d n1 = p1.n, n2 = p2.n;
    if (Dot(n1,n2)<0) n1 = -n1;
    // clip the polygon by the first plane
    Polygon pl1, pl2;
    clip(in, p1, pl1, pl2);
    Polygon pla, plb, plc, pld;
    if (!pl1.isNull()) clip(pl1, p2, pla, plb);
    // check if pla and plb lie in between two planes
    if (!pla.isNull() && inBetween(pla.center(),p1,p2)) return pla; 
    if (!plb.isNull() && inBetween(plb.center(),p1,p2)) return plb; 
    if (!pl2.isNull()) clip(pl2, p2, plc, pld);
    else return Polygon::nullPoly();
    if (!plc.isNull() && inBetween(plc.center(),p1,p2)) return plc; 
    if (!pld.isNull() && inBetween(pld.center(),p1,p2)) return pld; 
    return Polygon::nullPoly();
  }

  void geometry::parallelClip(Polygon in, Plane p1, Plane p2, Polygon& out) {
    // check the normal directions of two planes
    Vec3d n1 = p1.n, n2 = p2.n;
    if (Dot(n1,n2)<0) n1 = -n1;
    out.setNull();
    // clip the polygon by the first plane
    Polygon pl1, pl2;
    clip(in, p1, pl1, pl2);
    Polygon pla, plb, plc, pld;
    if (!pl1.isNull()) {
      clip(pl1, p2, pla, plb);
      // check if pla and plb lie in between two planes
      if (!pla.isNull() && inBetween(pla.center(),p1,p2)) {
        out.copy(pla);
        return;
      }
      if (!plb.isNull() && inBetween(plb.center(),p1,p2)) { 
        out.copy(plb);
        return;
      }
    }
    if (!pl2.isNull()) {
      clip(pl2, p2, plc, pld);
      if (!plc.isNull() && inBetween(plc.center(),p1,p2)) {
        out.copy(plc);
        return;
      } 
      if (!pld.isNull() && inBetween(pld.center(),p1,p2)) {
        out.copy(pld); 
        return;
      }
    }
  }

  bool geometry::sameSide(Point3 a, Point3 b, Point3 c, Point3 d) {
    Vec3d ab = a.toVec()-b.toVec();
    Vec3d ac = a.toVec()-c.toVec();
    Vec3d ad = a.toVec()-d.toVec();
    Vec3d nabc = Cross(ab,ac);
    Vec3d nabd = Cross(ab,ad);
    double isSame = Dot(nabc,nabd);
    if (isSame>=0) return true;
    else return false;
  }

  /**
   * Projection of a point to a segment.
   */
  Point3 geometry::project(Point3& p, Segment& s) {
    Vec3d x0 = p.toVec();  
    Vec3d x1 = s.s.toVec();  
    Vec3d x2 = s.e.toVec();
    Vec3d x01 = x1-x0;
    Vec3d x21 = x1-x2;
    double t = Dot(x01,x21)/Dot(x21,x21); 
    Vec3d xProj = x1-t*x21;
    return Point3(xProj.x, xProj.y, xProj.z);
  }

  /**
   * distance between a point and a segment.
   * http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
   * x0: point vector of the given point
   * x1: point vector of the starting point of the segment
   * x2: point vector of the   ending point of the segment
   */
  double geometry::dis(Point3& p, Segment& s) {
    Vec3d x0 = p.toVec();  
    Vec3d x1 = s.s.toVec();  
    Vec3d x2 = s.e.toVec();
    Vec3d x01 = x1-x0;
    Vec3d x21 = x1-x2;
    double t = Dot(x01,x21)/Dot(x21,x21);
    if (t<1 && t>0) {
      // the projection point is inside the segment, just return the distance
      Vec3d xProj = x1-t*x21;
      Vec3d vD = x0-xProj;
      return Length(vD);
    }
    else {
      // outside, just return the smaller distance to one of the endpoints
      double d1 = Length(x01);
      double d2 = Length(x0-x2);
      return d1<d2?d1:d2;
    }
    return 0;
  }

  /**
   * Project the point to the plane.
   */
  Point3 geometry::project(Point3& p, Plane& pl) {
    // first project the point to the normal vector of the plane
    float k = Dot(p.toVec()-pl.p.toVec(),pl.n); // the length on normal
    Vec3d pp = pl.p.toVec()+pl.n*k;
    // second get the vector pointing from the origin of the plane to the
    // projection point, which is as same as pp to p
    Vec3d pp2p = p.toVec()-pp;
    return Point3(pl.p._x+pp2p.x,pl.p._y+pp2p.y,pl.p._z+pp2p.z);
  }

  Segment geometry::project(Segment& s, Plane& pl) {
    Point3 stP = project(s.s,pl);
    Point3 enP = project(s.e,pl);
    return Segment(stP,enP);
  }

  /**
   * Project the polygon to the plane. 
   */
  Polygon geometry::project(Polygon& in, Plane& pl) {
    Point3* pp = new Point3[in.n];
    for (int i=0; i<in.n; ++i) pp[i] = project(in.p[i],pl);
    Polygon poly(pp,in.n);
    delete [] pp;
    return poly;
  }

  /**
   * Signed area of a triangle.
   * area = ab \cross ac (normal)
   * The actual area is the magnitude of the normal vector. 
   */
  Vec3d geometry::area(Point3& a, Point3& b, Point3& c) {
    Vec3d ab = b.toVec()-a.toVec();
    Vec3d ac = c.toVec()-a.toVec();
    return Cross(ab,ac)*0.5f;
  }
  
  /**
   * Barycentric coordinates (l1,l2,l3) of a point in a triangle.
   * Point: p, Triangle: p1, p2, p3
   * p = l1*p1+l2*p2+l3*p3, where l1+l2+l3=1.
   */
  void geometry::baryCoord(Point3& p, Point3& p1, Point3& p2, Point3& p3, 
    double& l1, double& l2, double& l3) {
    // compute areas
    double ar  = Length(area(p1,p2,p3));
    double ar1 = Length(area(p ,p2,p3));
    double ar2 = Length(area(p ,p1,p3));
    double ar3 = Length(area(p ,p1,p2));
    if (ar!=0) {
      double ar11 = 1.0f/ar;
      l1 = ar1*ar11;
      l2 = ar2*ar11;
      l3 = ar3*ar11;
    }
  }

  bool geometry::inside(Point3& p, Polygon& poly) {
    // loop over edges of the polygon, and form triangles with the query
    // point. if the areas of these triangles have same signs, the point is
    // inside, otherwise it is outside.
    Vec3d area1 = area(p,poly.p[0],poly.p[1]);
    for (int i=1; i<poly.n; ++i) {
      Vec3d area2 = area(p,poly.p[i],poly.p[(i+1)%poly.n]);
      if (Dot(area1,area2)<0) return false;
    }
    return true;
  }

  /**
   * The center of an arc ab is computed by:
   * (+-)Normalize(oa+ob)*r;
   * space = +-1, since two points cannot detemine one single arc.
   */
  Point3 geometry::centerArc(Circle c, Point3 a, Point3 b, int space) {
    Vec3d oa = a.toVec()-c.o.toVec();
    Vec3d ob = b.toVec()-c.o.toVec();
    Vec3d od = oa+ob;
    Normalize(od);
    Vec3d d = (od*c.r+c.o.toVec())*space;
    return Point3(d.x,d.y,d.z);
  }

  /**
   * Area bounded by an arc.
   * mid is the midpoint of the arc.
   */
  double geometry::arcArea(Circle c, Point3 a, Point3 b, Point3 mid) {
    double r = c.r;
    Point3 o = c.o;
    Vec3d oa = a.toVec()-o.toVec();
    Vec3d ob = b.toVec()-o.toVec();
    // get the angle between oa and ob
    Normalize(oa); Normalize(ob);
    double abCor = max(min(Dot(oa,ob),1.0),0.0);
    double theta = acos(abCor); //theta\in[0,PI]
    // >PI or <PI
    Vec3d od = oa+ob;
    Normalize(od);
    Point3 d = centerArc(c,a,b,1);
    // <PI
    if (d._x==mid._x && d._y==mid._y && d._z==mid._z) return 0.5f*r*r*theta;
    else return 0.5f*r*r*(PI*2-theta);
  }

  double geometry::lengthOverCircleSegment(Circle c, Segment s) {
    double t1, t2, len = 0; 
    Point3 p1, p2;
    bool sIn = c.inCircle(s.s);
    bool eIn = c.inCircle(s.e);
    int inters = intersect(c,s,t1,t2,p1,p2);
    if (inters==0) {
      if (sIn && eIn) return s.length();
      else return 0;
    }
    if (inters==2) return p1.dis(p2);
    // only one intersection, need to check which endpoint is inside the circle
    if (sIn) return p1.dis(s.s);
    if (eIn) return p1.dis(s.e);
    return 0;
  }

  /**
   * Compute the overlapping area of a circle and a polygon.
   */
  double geometry::areaOverCirclePolygon(Circle c, Polygon p) {
    double t1, t2, ar = 0;
    vector<Point3> vl; // vertex list of the intersecting polygon (with arcs).
    vector<int> sl; // status list: the vertex is on:1 or inside:0 the circle.
    vector<int> from; // which segment does the intersection come from
    bool in[p.n];
    int i, next, interNum;
    Point3 p1, p2;
    if (p.n<3) return 0;
    // orient the input polygon to follow the orientation of the circle
    /*Vec3d ab = p.p[1].toVec()-p.p[0].toVec();
    Vec3d ac = p.p[2].toVec()-p.p[0].toVec();
    Vec3d polyNormal = Cross(ab,ac);
    if (Length(polyNormal)>1e-4) 
      if (Dot(polyNormal,c.p.n)<0) p.reverse(); 
    else {
      if (p.n<4) return 0;
      Vec3d ad = p.p[3].toVec()-p.p[0].toVec();
      polyNormal = Cross(ab,ad);
      if (Dot(polyNormal,c.p.n)<0) p.reverse();
    }*/
    for (i=0; i<p.n; ++i) {
      in[i] = c.inCircle(p.p[i]);
    }
    // loop over polygon vertices to form overlapping shape
    for (i=0; i<p.n; ++i) {
      if (in[i]) {vl.push_back(p.p[i]); sl.push_back(0); from.push_back(0);}
      // check each edge and compute the intersections
      next = (i+1)%p.n;
      interNum = intersect(c,p.s[i],t1,t2,p1,p2);
      if (interNum>0) {vl.push_back(p1); sl.push_back(1); from.push_back(i);
        if (interNum>1) {vl.push_back(p2); sl.push_back(1); from.push_back(i);}
      }
    }
    Point3 pCenter = p.center();
    // loop over vertex list of the overlapping polygon, determine the edge
    // type: segement or arc
    vector<bool> isArc;
    isArc.resize(vl.size());
    // special case: only two vertices!
    // in this case, one edge is an arc, another must be a segment.
    if (vl.size()==2) {
      isArc[0] = true; isArc[1] = false;
    }
    else {
      for (i=0; i<vl.size(); ++i) {
        next = (i+1)%vl.size();
        // if one of these two vertices is inside, must be a segment
        if (sl[i]==0 || sl[next]==0) isArc[i] = false;
        else {
          // only two vertices cannot determine one single arc, we must
          // determine which half of the circle is used by looking up the center
          // of the input polygon p.
          // more than 2 inter points, always choose the shorter one
          Point3 d1 = centerArc(c,vl[i],vl[next], 1);
          Point3 d2 = centerArc(c,vl[i],vl[next],-1);
          Point3 d = d1; 
          //if (sameSide(vl[i],vl[next],pCenter,d1)) d = d1;
          //else d = d2;
          // the parents of these two intersections are the same==>segment 
          if (from[i]==from[next]) isArc[i] = false;
          else isArc[i] = true;
        }
      }
    }
    // compute area
    Vec3d areav(0,0,0);
    // normal vector defined by the orientation of the arcs
    Vec3d normal(0,0,0);
    for (i=0; i<vl.size(); ++i) {
      //vl[i].print(); cout<<sl[i]<<endl;
      //cout<<vl[i].dis(c.o)<<endl;
      next = (i+1)%vl.size();
      // check edge i-->next, with vertex i and vertex next
      if (!isArc[i]) areav = areav+area(c.o,vl[i],vl[next]);
      else {
        // only two vertices cannot determine one single arc, we must
        // determine which half of the circle is used by looking up the center
        // of the input polygon p.
        Point3 d1 = centerArc(c,vl[i],vl[next], 1);
        Point3 d2 = centerArc(c,vl[i],vl[next],-1);
        Point3 d; 
        int side = 1;
        if (vl.size()==2) {
          if (sameSide(vl[i],vl[next],pCenter,d1)) d = d1;
          else {d = d2; side = -1;}
          ar += arcArea(c,vl[i],vl[next],d);
        }
        else {
          // always choose the smaller area (<PI)
          // side 1 is always smaller
          float ar1 = arcArea(c,vl[i],vl[next],d1);
          float ar2 = arcArea(c,vl[i],vl[next],d2);
          ar += min(ar1,ar2);
          if (ar1>ar2) {side = -1; cout<<"-1 side"<<endl;}
        }
        Vec3d oa = vl[i].toVec()-c.o.toVec();
        Vec3d ob = vl[next].toVec()-c.o.toVec();
        normal = Cross(oa,ob)*side;
        Normalize(normal);
        //if (Dot(normal,c.p.n)>0) cout<<Dot(normal,c.p.n)<<endl;
        //vl[i].print(); cout<<endl;
        //vl[next].print(); cout<<endl<<"--------"<<endl;
        //cout<<"side="<<side<<"  ";
        //cout<<normal.x<<"  "<<normal.y<<"  "<<normal.z<<endl;
      }
    }
    double sumArea = ar*Dot(c.p.n,normal); // sign of arc area
    //if (Length(normal)!=0) 
    sumArea += Dot(areav,c.p.n);
    //else sumArea += Length(areav);
    return fabs(sumArea);
  }
