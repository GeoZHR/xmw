/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ssr;

import vec.*;
import util.*;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/** 
 * Estimate interval velocity from picked migrated velocity
 * @Author: Xinming Wu, University of Texas at Austin
 * @Version: 2016.07.21
 */

public class VelocityEstimator {

  /**
   * Constructs an impedance inverter.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   */
  public VelocityEstimator(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
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

  public void setSmoothness(float sc) {
    _sc = sc;
  }

  public void setTensors(Tensors2 d) {
    _d = d;
  }

  public float[][] applyForVelocity(
    float[][][] sp, float[][] w1, float[][] w2, float[][] v) 
  {
    int n2 = w1.length;
    int n1 = w2[0].length;
    float[][] b = new float[n2][n1];
    float[][] r = new float[n2][n1];
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(n1,n2,_sigma1,_sigma2,w2);
    A2 a2 = new A2(smoother2,_d,sp,w1,w2);
    CgSolver cg = new CgSolver(_small,_niter);
    makeRhs(w1,v,b);
    smoother2.applyTranspose(b);
    cg.solve(a2,vb,vr);
    smoother2.apply(r);
    return r;
  }

  public float[][] predictVelocity(float[][] v) {
    int n2 = v.length;
    int n1 = v[0].length;
    float[][] p = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float vt = 0f;
    for (int i1=0; i1<n1; ++i1) {
      float vi = v[i2][i1];
      vt += vi*vi;
      p[i2][i1] = sqrt(vt/(i1+1));
    }}
    return p;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private Tensors2 _d = null;
  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _small = 0.010f; // stop CG iterations if residuals are small
  private int _niter = 200; // maximum number of inner CG iterations
  private static float _sc = 0.5f;

  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, Tensors2 et,float[][][] sp, float[][] w1, float[][] w2) { 
      _s2=s2;
      _et=et;
      _sp=sp;
      _w1=w1;
      _w2=w2;
      testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      VecArrayFloat2 v2z = v2x.clone();
      float[][] x = v2z.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      _s2.apply(z);
      applyLhs(_et,_w1,_w2,z,y);
      if(_sp!=null) {
        screenLhs(_sp[0],_sp[1],z,y);
      }
      _s2.applyTranspose(y);
    }
    Smoother2 _s2;
    private Tensors2 _et = null;
    private float[][] _w1=null;
    private float[][] _w2=null;
    private float[][][] _sp=null;

    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n2 = _w1.length;
      int n1 = _w1[0].length;
      float[][] x = sub(randfloat(n1,n2),0.5f);
      float[][] y = sub(randfloat(n1,n2),0.5f);
      float[][] ax = zerofloat(n1,n2);
      float[][] ay = zerofloat(n1,n2);
      VecArrayFloat2 vx = new VecArrayFloat2(x);
      VecArrayFloat2 vy = new VecArrayFloat2(y);
      VecArrayFloat2 vax = new VecArrayFloat2(ax);
      VecArrayFloat2 vay = new VecArrayFloat2(ay);
      apply(vx,vax);
      apply(vy,vay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("A3: yax="+yax+" xay="+xay);
      System.out.println("A3: xax="+xax+" yay="+yay);
    }

  }

  private static void applyLhs(
    Tensors2 et, float[][] w1, float[][] w2, 
    float[][] x, float[][] y) 
  {
    applyLhs1(w1,x,y);
    applyLhs2(et,w2,x,y);
  }

  private static void applyLhs1(float[][] w, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] t = new float[n2][n1];
    applyIntegrate(x,t);
    mul(w,t,t);
    applyIntegrateTranspose(t,y);
  }

  private static void makeRhs(
    float[][] w, float[][] x, float[][] y) 
  {
    zero(y);
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] t = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      t[i2][i1] = x[i2][i1]*(i1+1);
    }}
    mul(w,t,t);
    applyIntegrateTranspose(t,y);
  }

  private static void applyLhs2(
    Tensors2 d, float[][] wp, float[][] x, float[][] y)
  {
    int n2 = x.length;
    int n1 = x[0].length;
    float[] ds = fillfloat(1.0f,3);
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if(d!=null){d.getTensor(i1,i2,ds);}
        float wpi = wp[i2][i1]*_sc;
        float wps = wpi*wpi;
        float d11 = ds[0]*wps;
        float d12 = ds[1]*wps;
        float d22 = ds[2]*wps;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }


  private static void screenLhs(
    float[][] cp, float[][] cm, float[][] x, float[][] y) 
  {
    int nc = cp[0].length;
    for (int ic=0; ic<nc; ++ic) {
      int i1p = (int)cp[0][ic];
      int i2p = (int)cp[1][ic];
      int i1m = (int)cm[0][ic];
      int i2m = (int)cm[1][ic];
      //float fls = fl[ic];//*fl[ic];

      float dx = 0.0f;

      dx += x[i2p][i1p];
      dx -= x[i2m][i1m];

      dx *= 5000;

      y[i2m][i1m] -= dx;
      y[i2p][i1p] += dx;
    }
  }


  private static void applyIntegrate(float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
      y[i2][0] = x[i2][0];
    }
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      y[i2][i1] = x[i2][i1]+y[i2][i1-1];
    }}
  }

  private static void applyIntegrateTranspose(float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
      y[i2][n1-1] = x[i2][n1-1];
    }
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=n1-2; i1>=0; --i1) {
      y[i2][i1] = x[i2][i1]+y[i2][i1+1];
    }}
  }


    // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(
      int n1, int n2, float sigma1, float sigma2, float[][] wp) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,_wp,x);
      smooth1(_sigma1,_wp,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,_wp,x);
      smooth2(_sigma2,_wp,x);
    }
    private float _sigma1,_sigma2;
    private float[][] _wp;
  }


  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  private static void smooth2(float sigma, float[][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }


    // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    int n2 = x.length;
    int n1 = x[0].length;
    float c = 0.5f*sigma*sigma;
    float[] st = fillfloat(1.0f,n1);
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      if (s!=null) {
        for (int i1=0; i1<n1; ++i1)
          st[i1] = s[i2][i1];
      }
      for (int i1=0; i1<n1; ++i1)
        xt[i1] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
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

}
