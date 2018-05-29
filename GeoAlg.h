#ifndef _GEOALG_H
#define _GEOALG_H

#include "Plane.h"
#include "Circle.h"
#include "Segment.h"
#include "Polygon.h"
#include "../util/ArrayMath.h"
#include <vector>

using namespace util;
using namespace std;

namespace geometry {
  
  double leftToSegment(Point3& p, Segment& s, Vec3d& n);
  /**
   * Intersection between a plane and a segment.
   * Plane: n \dot (x-p) = 0 ==> n \dot x - n \dot p = 0.
   * Segement: x = (1-t)s+te = s+t(e-s) = s+tv, t\in[0,1].
   * Solution: t = [n \dot (p-s)]/(n \dot v).
   * if n \dot v = 0: parallel, degenerate (not considered)
   * if t<0 or t>1, return false.
   * interP: the intersection point
   */
  bool intersect(Plane& p, Segment& s, double& t, Point3& interP);
  /**
   * Intersection between a plane and a triangle (stored in a polygon).
   * If intersect, return true, resulting segment stored in s.
   * Otherwise, return false.
   */
  bool intersect(Plane& pl, Polygon& poly, Segment& s);
  /**
   * Intersection between a segment and a circle.
   * The segment and the circle should be in the same plane.
   * If there is only one solution, we store this solution in p1 and t1.
   * Return the number of intersections.
   */
  int intersect(Circle& cir, Segment& s, double& t1, double& t2, 
    Point3& p1, Point3& p2);
  void clip(Segment s, Plane pl, Segment& s1, Segment& s2);
  /** 
   * Clip a polygon with a plane.
   * The vertex number of the polygon can only be increased by at most 1.
   * p1 and p2 are two output polygons.
   * If no intersection between the polygon and the plane, return false,
   * otherwise return true.
   */
  bool clip(Polygon poly, Plane pl, Polygon& p1, Polygon& p2);
  /*
   * Check the position of a point respect to a plane.
   *  1: in the half space with the normal pointing
   *  0: on plane
   * -1: in the half space back to the normal pointing
   */
  double halfSpace(Point3 p, Plane pl);
  /**
   * Check if a point appears in between two planes.
   * ==> this point should appear in different sides of the two planes.
   */
  bool inBetween(Point3 p, Plane p1, Plane p2);
  Segment parallelClip(Segment s, Plane p1, Plane p2);
  /**
   * Clip a polygon with 2 parallel planes.
   * Clip the polygon by these two planes sequentially. 
   * Check the resulting parts, return the polygon that lies between two
   * planes.
   */
  Polygon parallelClip(Polygon in, Plane p1, Plane p2);
  void parallelClip(Polygon in, Plane p1, Plane p2, Polygon& out);
  bool sameSide(Point3 a, Point3 b, Point3 c, Point3 d);
  /**
   * Projection of a point to a segment.
   */
  Point3 project(Point3& p, Segment& s);
  /**
   * distance between a point and a segment.
   * http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
   * x0: point vector of the given point
   * x1: point vector of the starting point of the segment
   * x2: point vector of the   ending point of the segment
   */
  double dis(Point3& p, Segment& s);
  /**
   * Project the point to the plane.
   */
  Point3 project(Point3& p, Plane& pl);
  Segment project(Segment& s, Plane& pl);
  /**
   * Project the polygon to the plane. 
   */
  Polygon project(Polygon& in, Plane& pl);
  /**
   * Signed area of a triangle.
   * area = ab \cross ac (normal)
   * The actual area is the magnitude of the normal vector. 
   */
  Vec3d area(Point3& a, Point3& b, Point3& c);
  /**
   * Barycentric coordinates (l1,l2,l3) of a point in a triangle.
   * Point: p, Triangle: p1, p2, p3
   * p = l1*p1+l2*p2+l3*p3, where l1+l2+l3=1.
   */
  void baryCoord(Point3& p, Point3& p1, Point3& p2, Point3& p3, 
    double& l1, double& l2, double& l3);
  bool inside(Point3& p, Polygon& poly);
  /**
   * The center of an arc ab is computed by:
   * (+-)Normalize(oa+ob)*r;
   * space = +-1, since two points cannot detemine one single arc.
   */
  Point3 centerArc(Circle c, Point3 a, Point3 b, int space);
  /**
   * Area bounded by an arc.
   * mid is the midpoint of the arc.
   */
  double arcArea(Circle c, Point3 a, Point3 b, Point3 mid);
  double lengthOverCircleSegment(Circle c, Segment s);
  /**
   * Compute the overlapping area of a circle and a polygon.
   */
  double areaOverCirclePolygon(Circle c, Polygon p);
}

#endif
