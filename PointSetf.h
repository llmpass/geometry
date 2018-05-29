#ifndef _POINTSETF_H
#define _POINTSETF_H

#include "Point3f.h"
#include "../ann/include/ANN/ANN.h"
#include <vector>

using namespace std;

namespace geometry {
  
class PointSetf {
  public:
    int n;
    Point3f* points;
    ANNkd_tree* pVertTree3;
    ANNdistArray dists;
    ANNidxArray idx;
    ANNpoint p;
    PointSetf() {
      n=0; points=NULL; pVertTree3=NULL;
      dists = new ANNdist[1];	
      idx = new ANNidx[1];
      p = annAllocPt(3);
    }
    PointSetf(float *x, float *y, float *z, int n) {
      this->n=n;
      // create Pnt arrays
      points = new Point3f[n];
      for (int i=0; i<n; ++i) points[i].set(x[i],y[i],z[i]);
      pVertTree3 = NULL;
      dists = new ANNdist[1];	
      idx = new ANNidx[1];
      p = annAllocPt(3);
    }
    PointSetf(float *x, float *y, float *z, 
      float* nx, float* ny, float* nz, int n) {
      this->n=n;
      // create Pnt arrays
      points = new Point3f[n];
      for (int i=0; i<n; ++i) points[i].set(x[i],y[i],z[i],nx[i],ny[i],nz[i]);
      pVertTree3 = NULL;
      dists = new ANNdist[1];	
      idx = new ANNidx[1];
      p = annAllocPt(3);
    }
    PointSetf(Point3f* pnts, int n) {
      this->n=n;
      // create Pnt arrays
      points = new Point3f[n];
      for (int i=0; i<n; ++i) points[i].set(pnts[i]._x,pnts[i]._y,pnts[i]._z);
      pVertTree3 = NULL;
      dists = new ANNdist[1];	
      idx = new ANNidx[1];
      p = annAllocPt(3);
    }
    ~PointSetf() {
      if (points!=NULL) delete [] points;
      if (pVertTree3!=NULL) delete pVertTree3;
    }

    Point3f& getP(int i) {
      return points[i];
    }

    void set(Point3f* pnts, int n) {
      this->n=n;
      // create Pnt arrays
      if (points!=NULL) delete [] points;
      points = new Point3f[n];
      for (int i=0; i<n; ++i) points[i].copy(pnts[i]);
      if (pVertTree3 != NULL) delete pVertTree3;
      pVertTree3 = NULL;
    }

    void buildAnnKDTree() {
      ANNpointArray vArray = annAllocPts(n,3);
      for (int i=0;  i<n; ++i) {
        vArray[i][0] = points[i]._x;
        vArray[i][1] = points[i]._y;
        vArray[i][2] = points[i]._z;
      }
      //Build kd Tree for nearest distance searching
      pVertTree3 = new ANNkd_tree(vArray,n,3);
    }
    // return the index and the distance of the nearest nabor to a given point 
    void nearestNabor(Point3f& pnt, int& ind, double& dis) {
      if (pVertTree3==NULL) buildAnnKDTree();
      // this is the square of distance array
      p[0] = pnt._x; p[1] = pnt._y; p[2] = pnt._z;
      pVertTree3->annkSearch(p,1,idx,dists,0);
      ind = idx[0]; 
      dis = dists[0];
    }

    void scale(float min1, float max1, float min2, float max2) {
      float x, y;
      float l2 = max2-min2;
      float l1 = max1-min1;
      if (l1<=0 || l2<=0) {
        cout<<"Invalid normalization arguments"<<endl;
        return;
      }
      for (int i=0; i<n; ++i) {
        points[i]._x = (points[i]._x-min2)/l2;
        points[i]._y = (points[i]._y-min1)/l1;
      }
    }
};
  
  /**
   * form a 2d point set from an image (ignore intensities)
   * get locations of nonzero pixels
   */
  void Img2PointSet(float** img, int n1, int n2, PointSetf& ps) {
    int np = 0;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        if (img[i2][i1]!=0) np++;
    Point3f* pnts = new Point3f[np];
    int c = 0;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        if (img[i2][i1]!=0) {
          pnts[c] = Point3f(i2,i1,0);
          c++;
        }
    ps.set(pnts,np);
    delete [] pnts;
  }

  void Img2PointSet(double** img, int n1, int n2, PointSetf& ps) {
    int np = 0;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        if (img[i2][i1]!=0) np++;
    Point3f* pnts = new Point3f[np];
    int c = 0;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        if (img[i2][i1]!=0) {
          pnts[c] = Point3f(i2,i1,0);
          c++;
        }
    ps.set(pnts,np);
    delete [] pnts;
  }

  void PointSet2Img(double**& img, int n1, int n2, double bin, 
    PointSetf& ps) {
    int i2, i1;
    double bin1 = 1.0/bin;
    for (int i=0; i<ps.n; ++i) {
      i2 = (int)(ps.points[i]._x*bin1+0.5);
      i1 = (int)(ps.points[i]._y*bin1+0.5); 
      if (i2>=0 && i2<n2 && i1>=0 && i1<n1) img[i2][i1] = 1;
    }
  }

  void PointSet2Img(double** img, int n1, int n2, 
    double bin1, double bin2, double min1, double min2, PointSetf& ps) {
    int i2, i1;
    double bin21 = 1.0/bin2;
    double bin11 = 1.0/bin1;
    for (int i=0; i<ps.n; ++i) {
      i2 = (int)((ps.points[i]._x-min2)*bin21+0.5);
      i1 = (int)((ps.points[i]._y-min1)*bin11+0.5); 
      if (i2>=0 && i2<n2 && i1>=0 && i1<n1) img[i2][i1] = 1;
    }
    /*for (i2=1; i2<n2-1; ++i2)
      for (i1=1; i1<n1-1; ++i1) {
        if (img[i2][i1]==1 && img[i2-1][i1]==0 && img[i2][i1-1]==0 && 
          img[i2+1][i1]==0 && img[i2][i1+1]==0) cout<<i2<<"  "<<i1<<endl;
        if (i2==260 && i1==166) {int iii; cin>>iii;}
      }*/
    cout<<"img from ps finished"<<endl;
  }

}

#endif

