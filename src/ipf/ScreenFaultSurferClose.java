/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.09.16
 */
public class ScreenFaultSurferClose {

  /**
   * Returns fault blocks for specified fault images.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   * @return array of fault blocks.
   */
  public float[][][] findBlocks(int n1, int n2, int n3, FaultCell[] fc) {
    return blocks(n1,n2,n3,fc);
  }

  public float[][][] innerProductField(int n1, int n2, int n3, FaultCell[] fc) {
    int nc = fc.length;
    float[][] xf = new float[3][nc];
    float[][] uf = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xf[0][ic] = fc[ic].i1;
      xf[1][ic] = fc[ic].i2;
      xf[2][ic] = fc[ic].i3;
      uf[0][ic] = fc[ic].w1;
      uf[1][ic] = fc[ic].w2;
      uf[2][ic] = fc[ic].w3;

    }
    KdTree kt = new KdTree(xf);
    float[][][] f = fillfloat(-1.0f,n1,n2,n3);
    float[] x = new float[3];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          x[0] = i1;
          x[1] = i2;
          x[2] = i3;
          int id = kt.findNearest(x);
          float x1 = fc[id].i1;
          float x2 = fc[id].i2;
          float x3 = fc[id].i3;
          float u1 = fc[id].w1;
          float u2 = fc[id].w2;
          float u3 = fc[id].w3;
          float d1 = i1-x1;
          float d2 = i2-x2;
          float d3 = i3-x3;
          float ds = sqrt(d1*d1+d2*d2+d3*d3);
          //if(ds>5.0f){continue;}
          if (ds==0.0f) {
            f[i3][i2][i1]=0.0f;
          } else {
            float dot = d1*u1+d2*u2+d3*u3;
            f[i3][i2][i1] = dot/ds; 
          }
        }
      }
    }
    return f;
  }

  public FaultCell[] removeOutliers(
    int n1, int n2, int n3, float theta, float du, FaultCell[] fci) 
  {
    //FaultCell[] fc1 = removeOutliers1(n1,n2,n3,theta,du,fci);
    FaultCell[] fc2 = removeOutliers2(n1,n2,n3,fci);
    return fc2;
  }

  private FaultCell[] removeOutliers1( 
    int n1, int n2, int n3, float theta, float du, FaultCell[] fci) {
    int nc = fci.length;
    int[] d = new int[]{12,6,6};
    float[][] xf = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xf[0][ic] = fci[ic].i1;
      xf[1][ic] = fci[ic].i2;
      xf[2][ic] = fci[ic].i3;
    }
    KdTree kt = new KdTree(xf);
    int nk = 0;
    int[] mk = new int[nc];
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    for (int ic=0; ic<nc; ++ic) {
      int i1 = fci[ic].i1;
      int i2 = fci[ic].i2;
      int i3 = fci[ic].i3;
      int[] ix = new int[]{i1,i2,i3};
      defineRange(d,ix,n1,n2,n3,xmin,xmax);
      int[] id = kt.findInRange(xmin,xmax);
      int nd = id.length;
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
      MedianFinder mf = new MedianFinder(nd);
      float u1m = mf.findMedian(u1);
      float u2m = mf.findMedian(u2);
      float u3m = mf.findMedian(u3);
      if(u1m*u1c+u2m*u2c+u3m*u3c<theta||
         abs(u1c-u1m)>du || abs(u2c-u2m)>du || abs(u3c-u3m)>du) 
      {nk ++; mk[ic] = 1;}
    }
    FaultCell[] fcr = new FaultCell[nc-nk];
    updateFaultCells(mk,fci,fcr);
    return fcr;
  }

  private FaultCell[] removeOutliers2(
    final int n1, final int n2, final int n3, final FaultCell[] fci) 
  {
    final int nc = fci.length;
    float[][] xf = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xf[0][ic] = fci[ic].i1;
      xf[1][ic] = fci[ic].i2;
      xf[2][ic] = fci[ic].i3;
    }
    final KdTree kt = new KdTree(xf);
    final int[] mk = new int[nc];
    Parallel.loop(nc,new Parallel.LoopInt() {
    public void compute(int ic) {
      int[] d = new int[3];
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      int i1 = fci[ic].i1;
      int i2 = fci[ic].i2;
      int i3 = fci[ic].i3;
      float uc1 = fci[ic].w1;
      float uc2 = fci[ic].w2;
      float uc3 = fci[ic].w3;
      float dx = 8.0f;
      d[0] = round(abs(dx*uc1));
      d[1] = round(abs(dx*uc3));
      d[2] = round(abs(dx*uc2));
      if(d[1]<2){d[1]=2;}
      if(d[2]<2){d[2]=2;}
      int[] ix = new int[]{i1,i2,i3};
      defineRange(d,ix,n1,n2,n3,xmin,xmax);
      float dxv = xmax[0]-xmin[0];
      float dxh = max(xmax[1]-xmin[1],xmax[2]-xmin[2]);
      int th = round(dxv*dxh*0.8f);
      int[] id = kt.findInRange(xmin,xmax);
      int cot = 0;
      int nd = id.length;
      for (int k=0; k<nd; ++k) {
        int ip = id[k];
        int p1 = fci[ip].i1;
        int p2 = fci[ip].i2;
        int p3 = fci[ip].i3;
        float ui1 = fci[ip].w1;
        float ui2 = fci[ip].w2;
        float ui3 = fci[ip].w3;
        float dot = uc1*ui1+uc2*ui2+uc3*ui3;
        if(dot<0.90f) {continue;}
        float dp = (p1-i1)*uc1+(p2-i2)*uc2+(p3-i3)*uc3;
        if(abs(dp)<2.0f) {cot++;}
      }
      if (cot<th){mk[ic] = 1;}
    }});
    FaultCell[] fcr = new FaultCell[nc-sum(mk)];
    updateFaultCells(mk,fci,fcr);
    return fcr;
  }

  private void updateFaultCells(int[] mk, FaultCell[] fc1, FaultCell[] fc2) {
    int ik = -1;
    int nc = fc1.length;
    for (int ic=0; ic<nc; ++ic) {
      if(mk[ic]==0) {
        ik ++; 
        fc2[ik] = fc1[ic];
      }
    }
  }

  private static void defineRange(int[] d, int[] i, 
    int n1, int n2, int n3, float[] xmin, float[] xmax) 
  {
    int i1 = i[0];
    int i2 = i[1];
    int i3 = i[2];
    int d1 = d[0];
    int d2 = d[1];
    int d3 = d[2];
    int i1m = i1-d1; if(i1m<0){i1m=0;}
    int i2m = i2-d2; if(i2m<0){i2m=0;}
    int i3m = i3-d3; if(i3m<0){i3m=0;}
    int i1p = i1+d1; if(i1p>=n1){i1p=n1-1;}
    int i2p = i2+d2; if(i2p>=n2){i2p=n2-1;}
    int i3p = i3+d3; if(i3p>=n3){i3p=n3-1;}
    xmin[0] = i1m; xmin[1] = i2m; xmin[2] = i3m;
    xmax[0] = i1p; xmax[1] = i2p; xmax[2] = i3p;
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
    float[][][] mk = new float[n3][n2][n1];
    float[][][] b  = new float[n3][n2][n1];
    normalsFromCellsM(fc,u1,u2,u3,wp,mk);
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
    A3 a = new A3(fc,wp,mk);
    CgSolver cs = new CgSolver(0.01,100);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    cs.solve(a,vr,vb);
    float v = 6.0f;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float mki = mk[i3][i2][i1];
          if(mki==0.0f) {b[i3][i2][i1]=v;}
        }
      }
    }

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
    A3(FaultCell[] fc, float[][][] wp,float[][][] mk) {
      _fc = fc;
      _wp = wp;
      _mk = mk;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      zero(y);
      applyLhs(_wp,_mk,x,y);
      //_ldk.apply(x,y);
    }
    private FaultCell[] _fc;
    private float[][][] _wp;
    private float[][][] _mk;
    /*
    private LocalDiffusionKernel _ldk =
      new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21);
    */
  }

  private static void applyLhs(
    float[][][] wp, float[][][] mk, float[][][] x, float[][][] y) 
  {
    LocalDiffusionKernel ldk =
      new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21);
    ldk.apply(x,y);
    float[][][] z = mul(y,0.1f); 
    screen(wp,mk,x,y);
    //add(z,y,y);
  }

  private static void screen(
    float[][][] wp, float[][][] mk, float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float mki = mk[i3][i2][i1];
          float wpi = wp[i3][i2][i1];
          if(mki==1.0f)
            y[i3][i2][i1] += x[i3][i2][i1]; 
          //y[i3][i2][i1] += wpi*wpi*x[i3][i2][i1]; 
        }
      }
    }
  }

  public static void normalsFromCells(FaultCell[] fc, 
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] wp) {
    int nc = fc.length;
    int n3 = u1.length;
    int n2 = u1[0].length;
    int of = 1;
    for (int ic=0; ic<nc; ic+=1) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;

      int i2m = fc[ic].i2m;
      int i3m = fc[ic].i3m;
      int i2p = fc[ic].i2p;
      int i3p = fc[ic].i3p;
      if(i2m<=i2p) {i2m-=of; i2p+=of;}
      else         {i2m+=of; i2p-=of;} 
      if(i3m<=i3p) {i3m-=of; i3p+=of;}
      else         {i3m+=of; i3p-=of;} 
      if(i2m<0||i3m<0){continue;}
      if(i2p<0||i3p<0){continue;}
      if(i2p>n2-1||i3p>n3-1){continue;}
      if(i2m>n2-1||i3m>n3-1){continue;}

      u1[i3m][i2m][i1] = -fc[ic].w1;
      u2[i3m][i2m][i1] = -fc[ic].w2;
      u3[i3m][i2m][i1] = -fc[ic].w3;

      u1[i3p][i2p][i1] =  fc[ic].w1;
      u2[i3p][i2p][i1] =  fc[ic].w2;
      u3[i3p][i2p][i1] =  fc[ic].w3;

      wp[i3 ][i2 ][i1] =  fc[ic].fl;
      wp[i3m][i2m][i1] =  fc[ic].fl;
      wp[i3p][i2p][i1] =  fc[ic].fl;
    }
  }

  public static void normalsFromCellsM(
    FaultCell[] fc, float[][][] u1, float[][][] u2, 
    float[][][] u3, float[][][] wp, float[][][] mk) 
  {
    int nc = fc.length; 
    int n3 = u1.length;
    int n2 = u1[0].length;
    float[][] xf = new float[3][nc*3];
    float[][] uf = new float[3][nc*3];
    float[][] wf = new float[1][nc*3];
    float[][] sg = new float[1][nc*3];
    int of = 1;
    int kk = 0;
    for (int ic=0; ic<nc; ic++) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;

      float u1i = fc[ic].w1;
      float u2i = fc[ic].w2;
      float u3i = fc[ic].w3;
      float wpi = fc[ic].fl;

      mk[i3][i2][i1] = 1.0f;
      xf[0][kk] = i1; xf[1][kk] = i2;xf[2][kk] = i3;
      uf[0][kk] =u1i; uf[1][kk] =u2i;uf[2][kk] = u3i;
      wf[0][kk] = wpi;sg[0][kk] =1.0f; kk++;
      /*

      int i2m = fc[ic].i2m;
      int i3m = fc[ic].i3m;
      int i2p = fc[ic].i2p;
      int i3p = fc[ic].i3p;
      if(i2m<=i2p) {i2m-=of; i2p+=of;}
      else         {i2m+=of; i2p-=of;} 
      if(i3m<=i3p) {i3m-=of; i3p+=of;}
      else         {i3m+=of; i3p-=of;} 
      if(i2m<0||i3m<0){continue;}
      if(i2p<0||i3p<0){continue;}
      if(i2p>n2-1||i3p>n3-1){continue;}
      if(i2m>n2-1||i3m>n3-1){continue;}
      u1[i3m][i2m][i1] = -u1i;
      u2[i3m][i2m][i1] = -u2i;
      u3[i3m][i2m][i1] = -u3i;
      u1[i3p][i2p][i1] =  u1i;
      u2[i3p][i2p][i1] =  u2i;
      u3[i3p][i2p][i1] =  u3i;
      wp[i3m][i2m][i1] =  wpi;
      wp[i3p][i2p][i1] =  wpi;
      mk[i3m][i2m][i1] = 1.0f;
      mk[i3p][i2p][i1] = 1.0f;

      xf[0][kk] = i1; xf[1][kk] = i2m;xf[2][kk] = i3m;
      uf[0][kk] =-u1i;uf[1][kk] =-u2i;uf[2][kk] = -u3i;
      wf[0][kk] = wpi;sg[0][kk] =1.0f; kk++;

      xf[0][kk] = i1; xf[1][kk] = i2p;xf[2][kk] = i3p;
      uf[0][kk] = u1i;uf[1][kk] = u2i;uf[2][kk] = u3i;
      wf[0][kk] = wpi;sg[0][kk] =1.0f; kk++;
      */
    }
    copy(kk,3,xf,xf);
    copy(kk,3,uf,uf);
    copy(kk,1,wf,wf);
    nearestInterp(xf,uf,wf,sg,u1,u2,u3,wp,mk);

    for (int ic=0; ic<nc; ic++) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      float u1i = fc[ic].w1;
      float u2i = fc[ic].w2;
      float u3i = fc[ic].w3;
      float wpi = fc[ic].fl;
      u1[i3][i2][i1] = u1i;
      u2[i3][i2][i1] = u2i;
      u3[i3][i2][i1] = u3i;
      wp[i3][i2][i1] = wpi;
      /*
      u1[i3][i2][i1] = 0.0f;
      u2[i3][i2][i1] = 0.0f;
      u3[i3][i2][i1] = 0.0f;
      wp[i3][i2][i1] = 0.0f;

      mk[i3][i2][i1] = 0.0f;

      int i2m = fc[ic].i2m;
      int i3m = fc[ic].i3m;
      int i2p = fc[ic].i2p;
      int i3p = fc[ic].i3p;
      if(i2m<=i2p) {i2m-=of; i2p+=of;}
      else         {i2m+=of; i2p-=of;} 
      if(i3m<=i3p) {i3m-=of; i3p+=of;}
      else         {i3m+=of; i3p-=of;} 
      if(i2m<0||i3m<0){continue;}
      if(i2p<0||i3p<0){continue;}
      if(i2p>n2-1||i3p>n3-1){continue;}
      if(i2m>n2-1||i3m>n3-1){continue;}
      u1[i3m][i2m][i1] = -u1i;
      u2[i3m][i2m][i1] = -u2i;
      u3[i3m][i2m][i1] = -u3i;
      u1[i3p][i2p][i1] =  u1i;
      u2[i3p][i2p][i1] =  u2i;
      u3[i3p][i2p][i1] =  u3i;
      wp[i3m][i2m][i1] =  wpi;
      wp[i3p][i2p][i1] =  wpi;
      mk[i3m][i2m][i1] = 1.0f;
      mk[i3p][i2p][i1] = 1.0f;
      */
    }
  }

  public static void nearestInterp(
    final float[][] xf, final float[][] uf, final float[][] wf,final float[][] sg,
    final float[][][] u1,final float[][][] u2, final float[][][] u3, 
    final float[][][] wp,final float[][][] mk) 
  {
    final int n3 = u1.length;
    final int n2 = u1[0].length;
    final int n1 = u1[0][0].length;
    final int[] d = new int[]{8,4,4};
    final KdTree kt = new KdTree(xf);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] xmin = new float[3];
          float[] xmax = new float[3];
          int[] ix = new int[]{i1,i2,i3};
          defineRange(d,ix,n1,n2,n3,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd==0){continue;}
          int ic = id[0];
          float dd = 800.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float x1i = xf[0][ip];
            float x2i = xf[1][ip];
            float x3i = xf[2][ip];
            float u1i = uf[0][ip];
            float u2i = uf[1][ip];
            float u3i = uf[2][ip];
            float d1i = i1-x1i;
            float d2i = i2-x2i;
            float d3i = i3-x3i;
            float dsi = abs(d1i*u1i+d2i*u2i+d3i*u3i);
            if(dsi<dd) {ic = ip;dd = dsi;}
          }
          if(dd==800.0f){continue;}
          if(sg[0][ic]==0.0f) {continue;}
          float u1i = uf[0][ic];
          float u2i = uf[1][ic];
          float u3i = uf[2][ic];
          float wpi = wf[0][ic];
          u1[i3][i2][i1] = u1i;
          u2[i3][i2][i1] = u2i;
          u3[i3][i2][i1] = u3i;
          wp[i3][i2][i1] = wpi;
          mk[i3][i2][i1] = 2.0f;
          //if(dd>=0.0f&&dd<0.5f){mk[i3][i2][i1]=1.0f;}
        }
      }
    }});
  }


  public static void nearestInterpM(
    final float[][] xf, final float[][] uf, final float[][] wf,final float[][] sg,
    final float[][][] u1,final float[][][] u2, final float[][][] u3, 
    final float[][][] wp,final float[][][] mk) 
  {
    final int n3 = u1.length;
    final int n2 = u1[0].length;
    final int n1 = u1[0][0].length;
    final KdTree kt = new KdTree(xf);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] ix = new float[]{i1,i2,i3};
          int id = kt.findNearest(ix);
          float x1i = xf[0][id];
          float x2i = xf[1][id];
          float x3i = xf[2][id];
          float d1i = i1-x1i;
          float d2i = i2-x2i;
          float d3i = i3-x3i;
          float dsi = sqrt(d1i*d1i+d2i*d2i+d3i*d3i);
          if(dsi>5.0f){continue;}
          float u1i = uf[0][id];
          float u2i = uf[1][id];
          float u3i = uf[2][id];
          float wpi = wf[0][id];
          u1[i3][i2][i1] = u1i;
          u2[i3][i2][i1] = u2i;
          u3[i3][i2][i1] = u3i;
          wp[i3][i2][i1] = wpi;
        }
      }
    }});
  }
}
