/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package slt;

import he.*;
import ipfx.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Compute a volume whoes zero isosurfaces are fault surfaces.
 * <p>
 * Assume we can obtain an image with fault normal vectors, then we can 
 * solve a Poisson problem to compute a volume whose zero isosurfaces are fault 
 * surfaces.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.14
 */
public class ScreenPoissonSurferC {

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

  public float[][][] getScreenMark(
    int n1, int n2, int n3, FaultCell[] fcs) {
    float[][][] mk = new float[n3][n2][n1];
    for (FaultCell fc:fcs) {
      if (fc==null){continue;}
      int i1 = round(fc.getX1());
      int i2 = round(fc.getX2());
      int i3 = round(fc.getX3());
      mk[i3][i2][i1] = fc.getFl();
    }
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k1m = 0;
      float mkm = 0.0f;
      for (int i1=0; i1<n1; ++i1) {
        float mki = mk[i3][i2][i1];
        if(mki>mkm) {
          k1m = i1;
          mkm = mki; 
        }
      }
      zero(mk[i3][i2]);
      mk[i3][i2][k1m] = mkm;
      if(k1m-1>=0) mk[i3][i2][k1m-1] = mkm;
      if(k1m+1<n1) mk[i3][i2][k1m+1] = mkm;
    }}
    return mk;
  }

  public float[][] updateScreenMark(
    float[][] ks, float[][] mk) {
    float[] k1 = ks[0];
    float[] k2 = ks[1];
    int n1 = mk[0].length;
    float[][] mp = copy(mk);
    CubicInterpolator ci = 
      new CubicInterpolator(CubicInterpolator.Method.LINEAR,k2,k1);
    int b2 = round(min(k2));
    int e2 = round(max(k2));
    float sigma = 10;
    for (int i2=b2; i2<=e2; ++i2) {
      int i1 = (int)ci.interpolate(i2);
      float d = min(abs(sub(k2,i2)));
      float s = exp(-0.5f*d*d/(sigma*sigma));
      System.out.println("s="+s);
      zero(mp[i2]);
      mp[i2][i1] = s;
      if(i1-1>=0) mp[i2][i1-1] = s;
      if(i1+1<n1) mp[i2][i1+1] = s;
    }
    /*
    int np = k1.length;
    for (int ip=0; ip<np; ++ip) {
      int i1 = (int)k1[ip];
      int i2 = (int)k2[ip];
      zero(mp[i2]);
      mp[i2][i1] = 1f;
      if(i1-1>=0) mp[i2][i1-1] = 1f;
      if(i1+1<n1) mp[i2][i1+1] = 1f;
    }
    */
    return mp;
  }

  public void updateScreenMark(
    float[][] ks, float[][][] mk) {
    float[] k1 = ks[0];
    float[] k2 = ks[1];
    float[] k3 = ks[2];
    int np = k1.length;
    int n1 = mk[0][0].length;
    for (int ip=0; ip<np; ++ip) {
      int i1 = (int)k1[ip];
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      zero(mk[i3][i2]);
      mk[i3][i2][i1] = 1f;
      if(i1-1>=0) mk[i3][i2][i1-1] = 1f;
      if(i1+1<n1) mk[i3][i2][i1+1] = 1f;
    }
  }


  /**
   * @param st array of thinned salt likelihoods.
   * @param u1 array of 1st component of salt normal vectors.
   * @param u2 array of 2nd component of salt normal vectors.
   */
  public float[][] saltIndicator(
    float[] k1, float[] k2, float[][] st, float[][] u1, float[][] u2)
  {
    int n2 = u1.length;
    int n1 = u1[0].length;
    float[][] b = new float[n2][n1]; // right-hand side
    float[][] f = new float[n2][n1]; // fault isosurface volume, in samples
    st = updateScreenMark(new float[][]{k1,k2},st);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vf = new VecArrayFloat2(f);
    A2 a2 = new A2(st);
    M2 m2 = new M2(_sigma1,_sigma2,k1,k2);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(u1,u2,b);
    cs.solve(a2,m2,vb,vf);
    return f;
  }

  /**
   * @param u1 array of 1st component of fault normal vectors.
   * @param u2 array of 2nd component of fault normal vectors.
   * @param u3 array of 3rd component of fault normal vectors.
   */
  public float[][][] saltIndicator(
    float[] k1, float[] k2, float[] k3, float[][][] mk, 
    float[][][] u1, float[][][] u2, float[][][] u3)
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] f = new float[n3][n2][n1]; // fault isosurface volume, in samples
    float[][] ks = extendControlPoints(k1,k2,k3,u1,u2,u3);
    updateScreenMark(ks,mk);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vf = new VecArrayFloat3(f);
    A3 a3 = new A3(mk);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,ks[0],ks[1],ks[2]);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(u1,u2,u3,b);
    cs.solve(a3,m3,vb,vf);
    return f;
  }


  /**
   * @param us[0] array of fault likelihoods.
   * @param us[1] array of 1st component of fault normal vectors.
   * @param us[2] array of 2nd component of fault normal vectors.
   * @param us[3] array of 3rd component of fault normal vectors.
   */
  public float[][][] saltIndicator(
    FaultCell[] cells, float[] k1, float[] k2, float[] k3, 
    float[][][] u1, float[][][] u2, float[][][] u3)
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float[][][] mk = getScreenMark(n1,n2,n3,cells);
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] f = new float[n3][n2][n1]; // fault isosurface volume, in samples
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vf = new VecArrayFloat3(f);
    A3 a3 = new A3(mk);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,k1,k2,k3);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(u1,u2,u3,b);
    cs.solve(a3,m3,vb,vf);
    return f;
  }

  public float[][][] findScreenPoints (
    final float smin,
    final float[][][] ss, final float[][][] u1, 
    final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = ss.length;
    final int n2 = ss[0].length;
    final int n1 = ss[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final float[][][] mk = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(1,n3-1,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=1; i2<n2-1 ;++i2) {
    for (int i1=1; i1<n1-1 ;++i1) {
      float sxi = ss[i3][i2][i1];
      float u1i = u1[i3][i2][i1]*2f;
      float u2i = u2[i3][i2][i1]*2f;
      float u3i = u3[i3][i2][i1]*2f;
      float x1m = i1-u1i;
      float x2m = i2-u2i;
      float x3m = i3-u3i;
      float x1p = i1+u1i;
      float x2p = i2+u2i;
      float x3p = i3+u3i;
      float sxm = si.interpolate(s1,s2,s3,ss,x1m,x2m,x3m);
      float sxp = si.interpolate(s1,s2,s3,ss,x1p,x2p,x3p);
      if (sxi>sxm && sxi>sxp && sxi>smin) 
        mk[i3][i2][i1] = sxi;
    }}}});
    return mk;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 6.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 400; // maximum number of CG iterations

  private float[][] extendControlPoints(float[] k1, float[] k2, float[] k3, 
    float[][][] u1, float[][][] u2, float[][][] u3) {
    int w2 = 3;
    int w3 = 3;
    int np = k1.length;
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    int ns = np*(w2+w2+1)*(w3+w3+1);
    float pmin = -8f;
    float pmax =  8f;
    float[] k1s = new float[ns];
    float[] k2s = new float[ns];
    float[] k3s = new float[ns];
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float p2i = -u2[i3][i2][i1]/u1[i3][i2][i1];
      float p3i = -u3[i3][i2][i1]/u1[i3][i2][i1];
      if(p2i<pmin) {p2i=pmin;}
      if(p2i>pmax) {p2i=pmax;}
      if(p3i<pmin) {p3i=pmin;}
      if(p3i>pmax) {p3i=pmax;}
      p2[i3][i2][i1] = p2i;
      p3[i3][i2][i1] = p3i;
    }}}
    int nk = 0;
    for (int ip=0; ip<np; ++ip) {
      int k1i = (int)k1[ip];
      int k2i = (int)k2[ip];
      int k3i = (int)k3[ip];
      int ib2 = k2i-w2; if(ib2<0   ) {ib2=0;   }
      int ib3 = k3i-w3; if(ib3<0   ) {ib3=0;   }
      int ie2 = k2i+w2; if(ie2>n2-1) {ie2=n2-1;} 
      int ie3 = k3i+w3; if(ie3>n3-1) {ie3=n3-1;} 
      int n2m = ie2-ib2+1;
      int n3m = ie3-ib3+1;
      float[][][] epm = fillfloat(1.0f,n1,n2m,n3m);
      float[][][] p2m = copy(n1,n2m,n3m,0,ib2,ib3,p2);
      float[][][] p3m = copy(n1,n2m,n3m,0,ib2,ib3,p3);
      int k2p = k2i-ib2;
      int k3p = k3i-ib3;
      float[] k1t = new float[]{k1i};
      float[] k2t = new float[]{k2p};
      float[] k3t = new float[]{k3p};
      SurfaceExtractorC se = new SurfaceExtractorC();
      se.setCG(0.01f,200);
      se.setExternalIterations(10);
      se.setSmoothings(10,10);
      float[][] surf = se.surfaceInitialization(n2m,n3m,n1-1,k1t,k2t,k3t);
      se.surfaceUpdateFromSlopes(epm,p2m,p3m,k1t,k2t,k3t,surf);
      for (int i3=0; i3<n3m; ++i3) {
      for (int i2=0; i2<n2m; ++i2) {
        int i1 = round(surf[i3][i2]);
        if(i1<0) {i1=0;} if(i1>=n1) {i1=n1-1;}
        k1s[nk] = i1;
        k2s[nk] = i2+ib2;
        k3s[nk] = i3+ib3;
        nk++;
      }}
    }
    k1s = copy(nk,0,k1s);
    k2s = copy(nk,0,k2s);
    k3s = copy(nk,0,k3s);
    return new float[][]{k1s,k2s,k3s};

  }

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(float[][] mk){
      _mk = mk;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      float[][] p = copy(x);
      VecArrayFloat2 v2p = new VecArrayFloat2(p);
      zero(p);
      zero(y);
      applyLhs(z,y);       //laplacian operator
      applyLhs(copy(y),p); //biharmonic operator
      screenLhs(_mk,z,y);  //screen points
      v2y.add(1f,v2p,20f);
    }
    private float[][] _mk;
  }

  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[] k1, float[] k2) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      if (k1!=null && k2!=null) {
        _k1 = copy(k1);
        _k2 = copy(k2);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k1,_k2,y);
      smooth2(_sigma2,y);
      smooth1(2.f*_sigma1,y);
      smooth2(_sigma2,y);
      constrain(_k1,_k2,y);
    }
    private float[] _k1,_k2;
    private float _sigma1,_sigma2;
  }

  private static class M3 implements CgSolver.A {
    M3(float sigma1, float sigma2, float sigma3, 
       float[] k1, float[] k2, float[] k3) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      if (k1!=null && k2!=null && k3!=null) {
        _k1 = copy(k1);
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      copy(x,y);
      constrain(_k1,_k2,_k3,y);
      smooth3(_sigma3,y);
      smooth2(_sigma2,y);
      smooth1(2.f*_sigma1,y);
      smooth2(_sigma2,y);
      smooth3(_sigma3,y);
      constrain(_k1,_k2,_k3,y);
    }
    private float[] _k1, _k2,_k3;
    private float _sigma1,_sigma2,_sigma3;
  }



  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(float[][][] mk){
      _mk = mk;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      float[][][] p = copy(x);
      VecArrayFloat3 v3p = new VecArrayFloat3(p);
      zero(p);
      zero(y);
      applyLhs(z,y);       //laplacian operator
      applyLhs(copy(y),p); //biharmonic operator
      screenLhs(_mk,z,y);  //screen points
      v3y.add(1f,v3p,20f);
    }
    private float[][][] _mk;
  }



  private static void constrain(float[] k1, float[] k2, float[][] x) {
    if (k1!=null && k2!=null) {
      int np = k1.length;
      for (int ip=0; ip<np; ++ip) {
        int i1 = (int)k1[ip]; 
        int i2 = (int)k2[ip]; 
        x[i2][i1] = 0.f;
      }
    }
  }

  private static void constrain(
    float[] k1, float[] k2, float[] k3, float[][][] x) {
    if (k1!=null && k2!=null && k3!=null) {
      int np = k1.length;
      for (int ip=0; ip<np; ++ip) {
        int i1 = (int)k1[ip]; 
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2][i1] = 0.f;
      }
    }
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

  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 3.
  private static void smooth3(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply3(x,x);
  }

  // right-hand side for 2D
  private static void makeRhs(float[][] u1, float[][] u2, float[][] y){
    int n2 = u1.length;
    int n1 = u1[0].length;
    for (int i2=1; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      float x1 = u1[i2][i1];
      float x2 = u2[i2][i1];
      float ya = 0.5f*(x1+x2);
      float yb = 0.5f*(x1-x2);
      y[i2  ][i1  ] += ya;
      y[i2  ][i1-1] -= yb;
      y[i2-1][i1  ] += yb;
      y[i2-1][i1-1] -= ya;
    }}
  }


  // right-hand side for 3D
  private static void makeRhs(
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1,i3m=0; i3<n3; ++i3,++i3m) {
    for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
    for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
      float y1 = u1[i3][i2][i1];
      float y2 = u2[i3][i2][i1];
      float y3 = u3[i3][i2][i1];
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
    }}}
  }



  // left-hand side for 2D
  private static void applyLhs(float[][] x, float[][] y){
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float ya = 0.5f*(x1+x2);
        float yb = 0.5f*(x1-x2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }


  // left-hand side for 3D
  private static void applyLhs(
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,x,y);
    }});
  }

  private static void applyLhsSlice3(
    int i3, float[][][] x, float[][][] y) 
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

        float y1 = 0.25f*(xa+xb+xc+xd);
        float y2 = 0.25f*(xa-xb+xc-xd);
        float y3 = 0.25f*(xa+xb-xc-xd);

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
      }
    }
  }

  private static void screenLhs(
    float[][] mk, float[][] x, float[][] y) 
  {
    int n2 = mk.length;
    int n1 = mk[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float mki = mk[i2][i1];
      if(mki>0f){
        y[i2][i1] += mki*x[i2][i1];
      }
    }}
  }


  private static void screenLhs(
    float[][][] mk, float[][][] x, float[][][] y) 
  {
    int n3 = mk.length;
    int n2 = mk[0].length;
    int n1 = mk[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float mki = mk[i3][i2][i1];
      if(mki>0f){
        mki = pow(mki,4);
        y[i3][i2][i1] += mki*x[i3][i2][i1];
      }
    }}}
  }
}
