/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import vec.*;
import util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.09.16
 */
public class ScreenFaultSurferOpen {

  /**
   * Returns fault blocks for specified fault images.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   * @return array of fault blocks.
   */
  public float[][][] findBlocks(int n1, int n2, int n3, FaultCell[] fc) {
    return blocks(n1,n2,n3,fc);
  }

  public FaultCell[] removeOutliers(
    int n1, int n2, int n3, FaultCell[] fci) 
  {
    int nc = fci.length;
    float[][] xf = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xf[0][ic] = fci[ic].i1;
      xf[1][ic] = fci[ic].i2;
      xf[2][ic] = fci[ic].i3;
    }
    KdTree kt = new KdTree(xf);
    int dh = 6;
    int dv = 12;
    int nk = 0;
    int[] mk = new int[nc];
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    for (int ic=0; ic<nc; ++ic) {
      int i1 = fci[ic].i1;
      int i2 = fci[ic].i2;
      int i3 = fci[ic].i3;
      int i1m = i1-dv; if(i1m<0){i1m=0;}
      int i2m = i2-dh; if(i2m<0){i2m=0;}
      int i3m = i3-dh; if(i3m<0){i3m=0;}
      int i1p = i1+dv; if(i1p>=n1){i1p=n1-1;}
      int i2p = i2+dh; if(i2p>=n2){i2p=n2-1;}
      int i3p = i3+dh; if(i3p>=n3){i3p=n3-1;}
      xmin[0] = i1m;
      xmin[1] = i2m;
      xmin[2] = i3m;
      xmax[0] = i1p;
      xmax[1] = i2p;
      xmax[2] = i3p;
      int[] id = kt.findInRange(xmin,xmax);
      int nd = id.length;
      MedianFinder mf = new MedianFinder(nd);
      float u1c = fci[ic].w1;
      float u2c = fci[ic].w2;
      float u3c = fci[ic].w3;
      float[] u1 = new float[nd];
      float[] u2 = new float[nd];
      float[] u3 = new float[nd];
      for (int i=0; i<nd; ++i) {
        int ip = id[i];
        u1[i] = fci[ip].w1;
        u2[i] = fci[ip].w2;
        u3[i] = fci[ip].w3;
      }
      float u1m = mf.findMedian(u1);
      float u2m = mf.findMedian(u2);
      float u3m = mf.findMedian(u3);
      float ns = dv*dv;
      if(bound(i1,i2,i3,n1,n2,n3)){ns=dv*dv*0.8f;}
      if(
         nd < ns ||
         u1m*u1c+u2m*u2c+u3m*u3c<0.95f||
         abs(u1c-u1m)>0.13f || 
         abs(u2c-u2m)>0.13f || 
         abs(u3c-u3m)>0.13f) 
      {
        nk ++;
        mk[ic] = 1;
      }
    }
    FaultCell[] fcr = new FaultCell[nc-nk];
    nk = -1;
    for (int ic=0; ic<nc; ++ic) {
      if(mk[ic]==0) {
        nk ++; fcr[nk] = fci[ic];
      }
    }
    return fcr; 
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  // Returns fault blocks.
  private static float[][][] blocks(int n1, int n2, int n3, FaultCell[] fc) {
    // Compute right-hand-side.
    float[][][] r  = new float[n3][n2][n1];
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] wp = new float[n3][n2][n1];
    float[][][] b =  new float[n3][n2][n1];
    //initialization(b);
    normalsFromCells(fc,u1,u2,u3,wp);
    for (int i3=0; i3<n3; ++i3) {
      int m3 = (i3>0)?i3-1:0;
      for (int i2=0; i2<n2; ++i2) {
        int m2 = (i2>0)?i2-1:0;
        for (int i1=0; i1<n1; ++i1) {
          int m1 = (i1>0)?i1-1:0;
          float fn1 = u1[i3][i2][i1];
          float fn2 = u2[i3][i2][i1];
          float fn3 = u3[i3][i2][i1];
          float fl1 = 0.5f*(wp[i3][i2][i1]+wp[i3][i2][m1]);
          float fl2 = 0.5f*(wp[i3][i2][i1]+wp[i3][m2][i1]);
          float fl3 = 0.5f*(wp[i3][i2][i1]+wp[m3][i2][i1]);
          float f1 = fn1*fl1;
          float f2 = fn2*fl2;
          float f3 = fn3*fl3;
          r[i3][i2][i1] += f1+f2+f3;
          r[i3][i2][m1] -= f1;
          r[i3][m2][i1] -= f2;
          r[m3][i2][i1] -= f3;
        }
      }
    }

    // Solve for fault blocks.
    A3 a = new A3(fc,wp);
    CgSolver cs = new CgSolver(0.01,100);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    cs.solve(a,vr,vb);
    return b;
  }

  private static void initialization(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float vi = 1.0f;
    for(int i2=0; i2<n2; ++i2) {
      for(int i1=0; i1<n1; ++i1) {
        x[0   ][i2][i1] = vi;
        x[1   ][i2][i1] = vi;
        x[n3-1][i2][i1] = vi;
        x[n3-2][i2][i1] = vi;
      }
    }
    for(int i3=0; i3<n3; ++i3) {
      for(int i1=0; i1<n1; ++i1) {
        x[i3][0   ][i1] = vi;
        x[i3][1   ][i1] = vi;
        x[i3][n2-1][i1] = vi;
        x[i3][n2-2][i1] = vi;
      }
    }
  }

  private static class A3 implements CgSolver.A {
    A3(FaultCell[] fc, float[][][] wp) {
      _fc = fc;
      _wp = wp;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      zero(y);
      applyLhs(_fc,_wp,x,y);
      //_ldk.apply(x,y);
    }
    private FaultCell[] _fc;
    private float[][][] _wp;
    /*
    private LocalDiffusionKernel _ldk =
      new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21);
    */
  }

  private static void applyLhs(
    FaultCell[] fc, float[][][] wp, float[][][] x, float[][][] y) 
  {
    LocalDiffusionKernel ldk =
      new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21);
    ldk.apply(x,y);
    int nc = fc.length;
    for (int ic=0; ic<nc; ++ic) {
      FaultCell fci = fc[ic];
      int i1 = fci.i1;
      int i2 = fci.i2;
      int i3 = fci.i3;
      y[i3][i2][i1] += x[i3][i2][i1]; //screen points with 0.0
    }
  }

  private static void normalsFromCells(FaultCell[] fc, 
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] wp) {
    int nc = fc.length;
    for (int ic=0; ic<nc; ic+=1) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      u1[i3][i2][i1] = fc[ic].w1;
      u2[i3][i2][i1] = fc[ic].w2;
      u3[i3][i2][i1] = fc[ic].w3;
      wp[i3][i2][i1] = fc[ic].fl;
    }
  }


  private static boolean bound(int i1, int i2, int i3, int n1, int n2, int n3) {
    if(i1<=2)   return true;
    if(i2<=2)   return true;
    if(i3<=2)   return true;
    if(i3>n3-3) return true;
    if(i2>n2-3) return true;
    if(i1>n1-3) return true;
    return false;
  }

}
