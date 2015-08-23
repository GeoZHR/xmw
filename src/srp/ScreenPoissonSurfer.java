/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package srp;

import ipf.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Compute a volume whoes zero isosurfaces are fault surfaces.
 * <p>
 * Assume we can obtain an image with fault normal vectors, then we can 
 * solve a Poisson problem to compute a volume whose zero isosurfaces are fault 
 * surfaces.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.02.09
 */
public class ScreenPoissonSurfer {

  /**
   * Sets half-widths for smoothings in 1st and 2nd dimensions.
   * These smoothings serve as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute mappings.
   * @param sigma1 half-width for smoothing in 1st dimension, in samples.
   * @param sigma2 half-width for smoothing in 2nd dimension, in samples.
   */
  public void setSmoothings(double sigma1, double sigma2, double sigma3) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _sigma3 = (float)sigma3;
  }

  /**
   * Sets parameters that control the number of solver iterations.
   * @param small stop iterations when error norm is reduced by this fraction.
   * @param niter maximum number of solver iterations.
   */
  public void setIterations(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  public float[][][][] setNormalsFromCells(
    int n1, int n2, int n3, FaultCell[] fcs) {
    float[][][][] us = new float[5][n3][n2][n1];
    for (FaultCell fc:fcs) {
      int i1 = round(fc.getX1());
      int i2 = round(fc.getX2());
      int i3 = round(fc.getX3());
      us[0][i3][i2][i1] = fc.getFl();
      us[1][i3][i2][i1] = fc.getW1();
      us[2][i3][i2][i1] = fc.getW2();
      us[3][i3][i2][i1] = fc.getW3();
      us[4][i3][i2][i1] = 1.00f;
    }
    interpNormals(fcs,us);
    return us;
  }

  private void interpNormals(final FaultCell[] fc, final float[][][][] us) {
    final int di = 20;
    final int n3 = us[0].length;
    final int n2 = us[0][0].length;
    final int n1 = us[0][0][0].length;
    final int[][] bb2 = new int[n1][2];
    final int[][] bb3 = new int[n1][2];
    final float[][] xp = setKdTreeNodes(n1,n2,n3,fc,bb2,bb3);
    final int[] bs1 = setBounds(n1,xp[0]);
    final int[] bs2 = setBounds(n2,xp[1]);
    final int[] bs3 = setBounds(n3,xp[2]);
    final KdTree kt = new KdTree(xp);
    final float sth = (float)sin(Math.PI/18.0);
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      System.out.println("i3="+i3);
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          if(us[4][i3][i2][i1]==1.0f){continue;}
          if((i2<bb2[i1][0]||i2>bb2[i1][1])){continue;}
          if((i3<bb3[i1][0]||i3>bb3[i1][1])){continue;}
          int[] id = null;
          xmin[0]  = i1-di;
          xmin[1]  = i2-di;
          xmin[2]  = i3-di;
          xmax[0]  = i1+di;
          xmax[1]  = i2+di;
          xmax[2]  = i3+di;
          id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<1){continue;}
          int ic = id[0]; 
          float st = sth;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float x1 = fc[ip].getX1();
            float x2 = fc[ip].getX2();
            float x3 = fc[ip].getX3();
            float w1 = fc[ip].getW1();
            float w2 = fc[ip].getW2();
            float w3 = fc[ip].getW3();
            float d1 = x1-i1;
            float d2 = x2-i2;
            float d3 = x3-i3;
            float d11 = d1*d1;
            float d22 = d2*d2;
            float d33 = d3*d3;
            float dsi = d11+d22+d33;
            float wdi = w1*d1+w2*d2+w3*d3;
            if(dsi==0f){ic=ip;st=0f;break;}
            float sti = abs(wdi/sqrt(dsi));
            if(sti<st) {ic=ip;st=sti;}
          }
          if(st==sth){continue;}
          us[1][i3][i2][i1] = fc[ic].getW1();
          us[2][i3][i2][i1] = fc[ic].getW2();
          us[3][i3][i2][i1] = fc[ic].getW3();
          us[0][i3][i2][i1] = fc[ic].getFl()*(1.f-st);
        }
      }
    }});

  }

  private int[] setBounds(int n, float[] x) {
    int[] bs = new int[2];
    int n1m = (int)min(x)-10; 
    int n1p = (int)max(x)+10; 
    if(n1m<0){n1m=0;} if(n1p>n){n1p=n;}
    bs[0] = n1m; bs[1] = n1p;
    return bs;
  }


  private float[][] setKdTreeNodes(
    int n1, int n2, int n3, 
    FaultCell[] fc, int[][] bb2, int[][] bb3)
  {
    int nc = fc.length;
    float[][] xp = new float[3][nc];
    for (int i1=0; i1<n1; ++i1) {
      bb2[i1][0] = n2; bb2[i1][1] = -n2;
      bb3[i1][0] = n3; bb3[i1][1] = -n3;
    }
    for (int ic=0; ic<nc; ic++) {
      int i1 = round(fc[ic].getX1());
      int i2 = round(fc[ic].getX2());
      int i3 = round(fc[ic].getX3());

      xp[0][ic] = i1;
      xp[1][ic] = i2;
      xp[2][ic] = i3;

      int b2l = bb2[i1][0];
      int b2r = bb2[i1][1];
      if(i2<b2l) bb2[i1][0] = i2;
      if(i2>b2r) bb2[i1][1] = i2;
      int b3l = bb3[i1][0];
      int b3r = bb3[i1][1];
      if(i3<b3l) bb3[i1][0] = i3;
      if(i3>b3r) bb3[i1][1] = i3;
    }

    for (int i1=0; i1<n1; ++i1) {
      bb2[i1][0] -= 10; bb3[i1][0] -= 10;
      bb2[i1][1] += 10; bb3[i1][1] += 10;
    }
    return xp;
  }

  /**
   * @param us[0] array of fault likelihoods.
   * @param us[1] array of 1st component of fault normal vectors.
   * @param us[2] array of 2nd component of fault normal vectors.
   * @param us[3] array of 3rd component of fault normal vectors.
   */
  public float[][][] faultIndicator(float[][][][] us){
    int n3 = us[0].length;
    int n2 = us[0][0].length;
    int n1 = us[0][0][0].length;
    // Compute shifts r(x1,x2,x3), in samples.
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] f = new float[n3][n2][n1]; // fault isosurface volume, in samples
    //initialization(f);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vf = new VecArrayFloat3(f);
    Smoother3 smoother3 = new Smoother3(n1,n2,n3,_sigma1,_sigma2,_sigma3,us[0]);
    A3 a3 = new A3(smoother3,us[0],us[4]);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(us[0],us[1],us[2],us[3],b);
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vf);
    smoother3.apply(f);
    mark(us[0],f);
    return f;
  }

  private void mark(float[][][] w, float[][][] f) {
    float v = -30.0f;
    int n3 = w.length;
    int n2 = w[0].length;
    int n1 = w[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wi = w[i3][i2][i1];
      if(wi==0.0f) {f[i3][i2][i1] = v;}
    }}}
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 6.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 100; // maximum number of CG iterations

  private void initialization(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          f[i3][i2][i1] = (float)i1;
        }
      }
    }
  }

  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, float[][][] wp, float[][][] mk){
      _s3 = s3;
      _wp = wp;
      _mk = mk;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      _s3.apply(z);
      zero(y);
      applyLhs(_wp,_mk,z,y);
      _s3.applyTranspose(y);
    }
    private Smoother3 _s3;
    private float[][][] _wp;
    private float[][][] _mk;
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother3 {
    public Smoother3(
      int n1, int n2, int n3, 
      float sigma1, float sigma2, float sigma3, float[][][] ep) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      _ep = ep;
      //testSpd();
    }
    public void apply(float[][][] x) {
      smooth3(_sigma3,_ep,x);
      smooth2(_sigma2,_ep,x);
      smooth1(_sigma1,x);
      //zero1(x);
    }
    public void applyTranspose(float[][][] x) {
      //zero1(x);
      smooth1(_sigma1,x);
      smooth2(_sigma2,_ep,x);
      smooth3(_sigma3,_ep,x);
    }
    private float _sigma1,_sigma2,_sigma3;
    private float[][][] _ep;
    private void zero1(float[][][] x) {
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          x[i3][i2][   0] = 0.0f;
          x[i3][i2][n1-1] = 0.0f;
        }
      }
    }
    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n1 = _ep[0][0].length;
      int n2 = _ep[0].length;
      int n3 = _ep.length;
      float[][][] x = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] y = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] ax = copy(x);
      float[][][] ay = copy(y);
      VecArrayFloat3 vx = new VecArrayFloat3(x);
      VecArrayFloat3 vy = new VecArrayFloat3(y);
      VecArrayFloat3 vax = new VecArrayFloat3(ax);
      VecArrayFloat3 vay = new VecArrayFloat3(ay);
      apply(ax);
      apply(ay);
      applyTranspose(ax);
      applyTranspose(ay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("S3: yax="+yax+" xay="+xay);
      System.out.println("S3: xax="+xax+" yay="+yay);
    }
  }

  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 1.
  private static void smooth2(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 1.
  private static void smooth3(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply3(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] st = fillfloat(1.0f,n2);
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      if (s!=null) {
        for (int i2=0; i2<n2; ++i2)
          st[i2] = s[i2][i1];
      }
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }
  private static void smooth2(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] s3 = (s!=null)?s[i3]:null;
      float[][] x3 = x[i3];
      smooth2(sigma,s3,x3);
    }});
  }

  // Smoothing for dimension 3.
  private static void smooth3(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] s2 = (s!=null)?new float[n3][]:null;
      float[][] x2 = new float[n3][];
      for (int i3=0; i3<n3; ++i3) {
        if (s!=null)
          s2[i3] = s[i3][i2];
        x2[i3] = x[i3][i2];
      }
      smooth2(sigma,s2,x2);
    }});
  }

  private static void makeRhs(
    float[][][] wp, float[][][] u1, float[][][] u2, float[][][] u3, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1,i3m=0; i3<n3; ++i3,++i3m) {
      for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
        for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
          float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
          float wps = wpi*wpi;
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float y1 = wps*u1i;
          float y2 = wps*u2i;
          float y3 = wps*u3i;
          float ya = 0.25f*(y1+y2+y3);
          float yb = 0.25f*(y1-y2+y3);
          float yc = 0.25f*(y1+y2-y3);
          float yd = 0.25f*(y1-y2-y3);
          y[i3 ][i2 ][i1 ] += ya;
          y[i3 ][i2 ][i1m] -= yd;
          y[i3 ][i2m][i1 ] += yb;
          y[i3 ][i2m][i1m] -= yc;
          y[i3m][i2 ][i1 ] += yc;
          y[i3m][i2 ][i1m] -= yb;
          y[i3m][i2m][i1 ] += yd;
          y[i3m][i2m][i1m] -= ya;
        }
      }
    }
  }
  private static void applyLhs(
    final float[][][] wp, final float[][][] mk, 
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,mk,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,mk,x,y);
    }});
  }
  private static void applyLhsSlice3(
    int i3, float[][][] wp, float[][][] mk, float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
        float wps = wpi*wpi;
        float xa = 0.0f;
        float xb = 0.0f;
        float xc = 0.0f;
        float xd = 0.0f;
        xa += x00[i1 ];
        xd -= x00[i1m];
        xb += x01[i1 ];
        xc -= x01[i1m];
        xc += x10[i1 ];
        xb -= x10[i1m];
        xd += x11[i1 ];
        xa -= x11[i1m];
        float x1 = 0.25f*(xa+xb+xc+xd);
        float x2 = 0.25f*(xa-xb+xc-xd);
        float x3 = 0.25f*(xa+xb-xc-xd);
        float y1 = wps*x1;
        float y2 = wps*x2;
        float y3 = wps*x3;
        float ya = 0.25f*(y1+y2+y3);
        float yb = 0.25f*(y1-y2+y3);
        float yc = 0.25f*(y1+y2-y3);
        float yd = 0.25f*(y1-y2-y3);
        y00[i1 ] += ya;
        y00[i1m] -= yd;
        y01[i1 ] += yb;
        y01[i1m] -= yc;
        y10[i1 ] += yc;
        y10[i1m] -= yb;
        y11[i1 ] += yd;
        y11[i1m] -= ya;
        if(mk[i3][i2][i1]==1.0f){
          y00[i1] += wps*x00[i1];
        }
      }
    }
  }
}
