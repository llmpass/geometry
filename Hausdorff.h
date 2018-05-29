#ifndef _HAUSDORFF_H
#define _HAUSDORFF_H

#include "PointSet.h"

/**
 * Hausdorff distance square between two point sets.
 */
namespace geometry {
  class Hausdorff {
    public:
    static double dis2(PointSet& s, PointSet& t) {
      if (s.pVertTree3==NULL) s.buildAnnKDTree();
      if (t.pVertTree3==NULL) t.buildAnnKDTree();
      // compute the min distance between a point in s and the t set.
      int ind; 
      double dis;
      // max of the min distances
      double maxD = -9999;
      for (int i=0; i<s.n; ++i) {
        t.nearestNabor(s.getP(i),ind,dis);
        if (dis>maxD) maxD = dis;
      }
      // compute the min distance between a point in t and the s set.
      for (int i=0; i<t.n; ++i) {
        s.nearestNabor(t.getP(i),ind,dis);
        if (dis>maxD) maxD = dis;      
      }
      return maxD;
    }
  };
}

#endif
