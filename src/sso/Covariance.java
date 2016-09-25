/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package sso;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;


import vec.*;
import util.*;

/**
 * Covariance-based semblance. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.15
 */

public class Covariance {

  // Computes fault semblance numerators and denominators.
  public float[][][] covarianceEigen(
    final int m1, final float[][] p, final float[][] f) 
  {
    final int d1 = (m1-1)/2;
    final int n2 = f.length;
    final int n1 = f[0].length;
    System.out.println("d1="+d1);
    final float[][] em = new float[n2][n1];
    final float[][] es = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    loop(n2,new LoopInt() {
    public void compute(int i2) {
      float[] xm = new float[n1];
      float[] xp = new float[n1];
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float[][] h = new float[3][n1];
      float[] fm  = f[i2m];
      float[] fp  = f[i2p];
      float[] f0  = f[i2 ];
      float[] pm  = p[i2m];
      float[] pp  = p[i2p];
      float[] em2 = em[i2];
      float[] es2 = es[i2];
      for (int i1=0; i1<n1; ++i1) {
        xm[i1] = i1-pm[i1];
        xp[i1] = i1+pp[i1];
      }
      h[0] = f0;
      si.interpolate(n1,1.0,0.0,fm,n1,xm,h[1]);
      si.interpolate(n1,1.0,0.0,fp,n1,xp,h[2]);
      if (            i2==0   ) h[1] = h[0];
      if (            i2==n2-1) h[2] = h[0];
      double[][] c = new double[2][2];
      for (int i1=d1; i1<n1-d1; ++i1) {
        for (int k2=0; k2<2; ++k2) {
          float[] h2 = h[k2];
          double[] c2 = c[k2];
        for (int k1=0; k1<2; ++k1) {
          float ci = 0.0f;
          float[] h1 = h[k1];
          for (int i=i1-d1; i<=i1+d1; ++i)
            ci += h1[i]*h2[i];
          c2[k1] = ci;
        }}
        DMatrix dm = new DMatrix(c);
        DMatrixEvd ed = new DMatrixEvd(dm);
        double[] es = ed.getRealEigenvalues();
        em2[i1] = (float)max(es);
        es2[i1] = (float)sum(es);
      }
    }});
    // pad the top and bottom boundaries
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<d1; ++i1) {
      em[i2][i1] = em[i2][d1];
      es[i2][i1] = es[i2][d1];
      em[i2][n1-1-i1] = em[i2][n1-d1-1];
      es[i2][n1-1-i1] = es[i2][n1-d1-1];
    }}
    return new float[][][]{em,es};
  }


  // Computes fault semblance numerators and denominators.
  public float[][][][] covarianceEigen(
    final int m1,
    final float[][][] p2, final float[][][] p3, final float[][][] f) 
  {
    final int d1 = (m1-1)/2;
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n1 = f[0][0].length;
    final float[] c = new float[1];
    final float[][][] em = new float[n3][n2][n1];
    final float[][][] es = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      float[] xmm = new float[n1];
      float[] xm0 = new float[n1];
      float[] xmp = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xpm = new float[n1];
      float[] xp0 = new float[n1];
      float[] xpp = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      float[][] h = new float[9][n1];
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm  = f[i3m][i2m];
        float[] fm0  = f[i3m][i2 ];
        float[] fmp  = f[i3m][i2p];
        float[] f0m  = f[i3 ][i2m];
        float[] f00  = f[i3 ][i2 ];
        float[] f0p  = f[i3 ][i2p];
        float[] fpm  = f[i3p][i2m];
        float[] fp0  = f[i3p][i2 ];
        float[] fpp  = f[i3p][i2p];
        float[] p2mm = p2[i3m][i2m];
        float[] p2mp = p2[i3m][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p20p = p2[i3 ][i2p];
        float[] p2pm = p2[i3p][i2m];
        float[] p2pp = p2[i3p][i2p];
        float[] p3mm = p3[i3m][i2m];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3mp = p3[i3m][i2p];
        float[] p3pm = p3[i3p][i2m];
        float[] p3p0 = p3[i3p][i2 ];
        float[] p3pp = p3[i3p][i2p];
        float[] em32 = em[i3][i2];
        float[] es32 = es[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          xmm[i1] = i1-p3mm[i1]-p2mm[i1];
          xm0[i1] = i1-p3m0[i1]         ;
          xmp[i1] = i1-p3mp[i1]+p2mp[i1];
          x0m[i1] = i1         -p20m[i1];
          x0p[i1] = i1         +p20p[i1];
          xpm[i1] = i1+p3pm[i1]-p2pm[i1];
          xp0[i1] = i1+p3p0[i1]         ;
          xpp[i1] = i1+p3pp[i1]+p2pp[i1];
        }
        h[0] = f00;
        si.interpolate(n1,1.0,0.0,f0m,n1,x0m,h[1]);
        si.interpolate(n1,1.0,0.0,f0p,n1,x0p,h[2]);
        si.interpolate(n1,1.0,0.0,fm0,n1,xm0,h[3]);
        si.interpolate(n1,1.0,0.0,fp0,n1,xp0,h[4]);
        si.interpolate(n1,1.0,0.0,fmm,n1,xmm,h[5]);
        si.interpolate(n1,1.0,0.0,fmp,n1,xmp,h[6]);
        si.interpolate(n1,1.0,0.0,fpm,n1,xpm,h[7]);
        si.interpolate(n1,1.0,0.0,fpp,n1,xpp,h[8]);
        if (            i2==0   ) h[1] = h[0];
        if (            i2==n2-1) h[2] = h[0];
        if (i3==0               ) h[3] = h[0];
        if (i3==n3-1            ) h[4] = h[0];
        if (i3==0    && i2==0   ) h[5] = h[0];
        if (i3==0    && i2==n2-1) h[6] = h[0];
        if (i3==n3-1 && i2==0   ) h[7] = h[0];
        if (i3==n3-1 && i2==n2-1) h[8] = h[0];
        double[][] cm = new double[9][9];
        for (int i1=d1; i1<n1-d1; ++i1) {
          for (int k2=0; k2<9; ++k2) {
            float[] h2 = h[k2];
            double[] cm2 = cm[k2];
          for (int k1=0; k1<9; ++k1) {
            float cmi = 0.0f;
            float[] h1 = h[k1];
            for (int i=i1-d1; i<=i1+d1; ++i)
              cmi += h1[i]*h2[i];
            cm2[k1] = cmi;
          }}
          DMatrix dm = new DMatrix(cm);
          DMatrixEvd ed = new DMatrixEvd(dm);
          double[] es = ed.getRealEigenvalues();
          em32[i1] = (float)max(es);
          es32[i1] = (float)sum(es);
        }
      }
      c[0]=c[0]+1f;
      System.out.println(100f*c[0]/n3+"% "+"done...");
    }});
    System.out.println("done");
    // pad the top and bottom boundaries
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<d1; ++i1) {
      em[i3][i2][i1] = em[i3][i2][d1];
      es[i3][i2][i1] = es[i3][i2][d1];
      em[i3][i2][n1-1-i1] = em[i3][i2][n1-d1-1];
      es[i3][i2][n1-1-i1] = es[i3][i2][n1-d1-1];
    }}}
    return new float[][][][]{em,es};
  }

  public float[][][] thin(final EigenTensors3 et, final float[][][] fx) {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final float[][][] gx = fillfloat(1f,n1,n2,n3);
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float[] u = et.getEigenvectorV(i1,i2,i3);
        float x1p = i1+u[0];
        float x2p = i2+u[1];
        float x3p = i3+u[2];
        float x1m = i1-u[0];
        float x2m = i2-u[1];
        float x3m = i3-u[2];
        float fxi = fx[i3][i2][i1];
        float fxp = si.interpolate(s1,s2,s3,fx,x1p,x2p,x3p);
        float fxm = si.interpolate(s1,s2,s3,fx,x1m,x2m,x3m);
        if(fxi<=fxp&&fxi<=fxm) {
          gx[i3][i2][i1] = fxi;
          if(abs(u[1])>abs(u[2])) {
            int i2m = max(i2-1,0);
            int i2p = min(i2+1,n2-1);
            gx[i3][i2m][i1] = fx[i3][i2m][i1];
            gx[i3][i2p][i1] = fx[i3][i2p][i1];
          } else {
            int i3m = max(i3-1,0);
            int i3p = min(i3+1,n3-1);
            gx[i3m][i2][i1] = fx[i3m][i2][i1];
            gx[i3p][i2][i1] = fx[i3p][i2][i1];

          }
        }
      }}
    }});
    return gx;
  }

  public float[][][][] covarianceEigenX(
    final int m1, final float[][][] fx, final float[][][] gx) 
  {
    final int d1 = (m1-1)/2;
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] em = new float[n3][n2][n1];
    final float[][][] es = new float[n3][n2][n1];
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      System.out.println("i3="+i3);
      float[] e = new float[2];
      float[][] a = new float[2][2];
      float[][] z = new float[2][2];
      for (int i2=0; i2<n2; ++i2) {
        float[] fx32 = fx[i3][i2];
        float[] gx32 = gx[i3][i2];
        float[] em32 = em[i3][i2];
        float[] es32 = es[i3][i2];
        for (int i1=d1; i1<n1-d1; ++i1) {
          float a11 = 0.0f;
          float a12 = 0.0f;
          float a22 = 0.0f;
          for (int i=i1-d1; i<=i1+d1; ++i) {
            float fxi = fx32[i];
            float gxi = gx32[i];
            a11 += fxi*fxi;
            a12 += fxi*gxi;
            a22 += gxi*gxi;
          }
          a[0][0] = a11;
          a[0][1] = a12;
          a[1][0] = a12;
          a[1][1] = a22;
          Eigen.solveSymmetric22(a,z,e);
          em32[i1] = max(e);
          es32[i1] = sum(e);
        }
      }
    }});
    // pad the top and bottom boundaries
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<d1; ++i1) {
      em[i3][i2][i1] = em[i3][i2][d1];
      es[i3][i2][i1] = es[i3][i2][d1];
      em[i3][i2][n1-1-i1] = em[i3][i2][n1-d1-1];
      es[i3][i2][n1-1-i1] = es[i3][i2][n1-d1-1];
    }}}
    return new float[][][][]{em,es};
  }

  public static float[][][] smooth(
    final double sigma, final float au, final float[][][] p2, 
    final float[][][] p3, final float[][][] g) {
    final int n3 = g.length;
    final int n2 = g[0].length;
    final int n1 = g[0][0].length;
    final EigenTensors3 d = new EigenTensors3(n1,n2,n3,true);
    d.setEigenvalues(au,0.001f,1.00f);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float u1i = 1.0f/sqrt(1.0f+p2i*p2i+p3i*p3i);
          float u2i = -p2i*u1i;
          float u3i = -p3i*u1i;
          float usi = 1.0f/sqrt(u1i*u1i+u2i*u2i);
          float w1i = -u2i*usi;
          float w2i =  u1i*usi;
          float w3i = 0.0f;
          d.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          d.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
        }
      }
    }});
    float c = (float)(0.5*sigma*sigma);
    float[][][] h = new float[n3][n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,c,g,h);
    return h;
  }

  public float[][][] smooth(
    float sigma, EigenTensors3 et3, float[][][] num, float[][][] den) 
  {
    int n3 = den.length;
    int n2 = den[0].length;
    int n1 = den[0][0].length;
    float[][][] w = sub(1f,div(num,den));
    float[][][] b = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    Smoother3 smoother3 = new Smoother3(sigma, et3);
    A3 a3 = new A3(smoother3,mul(w,den));
    CgSolver cs = new CgSolver(0.001,20);
    mul(w,num,b);
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vr);
    return r;
  }

  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, float[][][] wp) 
    {
      _s3 = s3;
      _wp = wp;
      float n3 = wp.length;
      float n2 = wp[0].length;
      float n1 = wp[0][0].length;
      _sc = sum(wp)/(n1*n2*n3);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      v3y.zero();
      _s3.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s3.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother3 _s3;
    private float[][][] _wp;
  }

  private static void applyLhs(float[][][] wp, float[][][] x, float[][][] y) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    for (int i3=0; i3<n3; ++i3)
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      y[i3][i2][i1] += wp[i3][i2][i1]*x[i3][i2][i1];
  }

  private static void addAndScale(float sc, float[][][] x, float[][][] y) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1) 
      y[i3][i2][i1] += sc*x[i3][i2][i1];
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother3 {
    public Smoother3(float sigma, EigenTensors3 et) {
      _et = et;
      _sigma = sigma;
    }
    public void apply(float[][][] x) {
      smooth(_sigma,_et,x);
      //smooth1(_sigma,x);
      //smooth2(_sigma,x);
      //smooth3(_sigma,x);
    }
    public void applyTranspose(float[][][] x) {
      smooth(_sigma,_et,x);
      //smooth3(_sigma,x);
      //smooth2(_sigma,x);
      //smooth1(_sigma,x);

    }
    private float _sigma;
    private EigenTensors3 _et;
  }

  // Smoothing for dimension 2.
  private static void smooth(float sigma, EigenTensors3 et, float[][][] x) {
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(et,sigma,copy(x),x);
  }

  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth3(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply3(x,x);
  }



}
