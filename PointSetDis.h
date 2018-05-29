#ifndef _POINTSETDIS_H
#define _POINTSETDIS_H

#include "PointSet.h"
#include "PointSetf.h"

namespace geometry {

  /**
   * L1, L2, Lmax (Hausdorff) distance between two point sets.
   */
  double Lmax(PointSet& s, PointSet& t) {
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
    return sqrt(maxD);
  }

  double Lmax(PointSetf& s, PointSetf& t) {
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
    return sqrt(maxD);
  }

  double L1(PointSet& s, PointSet& t) {
    //if (s.pVertTree3==NULL) s.buildAnnKDTree();
    //if (t.pVertTree3==NULL) t.buildAnnKDTree();
    // compute the min distance between a point in s and the t set.
    int ind; 
    double dis;
    // summation of min distances
    double sum = 0;
    for (int i=0; i<s.n; ++i) {
      t.nearestNabor(s.getP(i),ind,dis);
      sum += sqrt(dis);
    }
    // compute the min distance between a point in t and the s set.
    for (int i=0; i<t.n; ++i) {
      s.nearestNabor(t.getP(i),ind,dis);
      sum += sqrt(dis);     
    }
    return sum/(s.n+t.n);
  }

  double L1(PointSetf& s, PointSetf& t) {
    //if (s.pVertTree3==NULL) s.buildAnnKDTree();
    //if (t.pVertTree3==NULL) t.buildAnnKDTree();
    // compute the min distance between a point in s and the t set.
    int ind; 
    double dis;
    // summation of min distances
    double sum = 0;
    for (int i=0; i<s.n; ++i) {
      t.nearestNabor(s.getP(i),ind,dis);
      sum += sqrt(dis);
    }
    // compute the min distance between a point in t and the s set.
    for (int i=0; i<t.n; ++i) {
      s.nearestNabor(t.getP(i),ind,dis);
      sum += sqrt(dis);     
    }
    return sum/(s.n+t.n);
  }
  
  double L2(PointSet& s, PointSet& t) {
    if (s.pVertTree3==NULL) s.buildAnnKDTree();
    if (t.pVertTree3==NULL) t.buildAnnKDTree();
    // compute the min distance between a point in s and the t set.
    int ind; 
    double dis;
    // summation of min distances
    double sum2 = 0;
    for (int i=0; i<s.n; ++i) {
      t.nearestNabor(s.getP(i),ind,dis);
      sum2 += dis;
    }
    // compute the min distance between a point in t and the s set.
    for (int i=0; i<t.n; ++i) {
      s.nearestNabor(t.getP(i),ind,dis);
      sum2 += dis;     
    }
    return sqrt(sum2/(s.n+t.n));
  }

  double L2(PointSetf& s, PointSetf& t) {
    if (s.pVertTree3==NULL) s.buildAnnKDTree();
    if (t.pVertTree3==NULL) t.buildAnnKDTree();
    // compute the min distance between a point in s and the t set.
    int ind; 
    double dis;
    // summation of min distances
    double sum2 = 0;
    for (int i=0; i<s.n; ++i) {
      t.nearestNabor(s.getP(i),ind,dis);
      sum2 += dis;
    }
    // compute the min distance between a point in t and the s set.
    for (int i=0; i<t.n; ++i) {
      s.nearestNabor(t.getP(i),ind,dis);
      sum2 += dis;     
    }
    return sqrt(sum2/(s.n+t.n));
  }
}

#endif
