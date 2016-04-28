/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package aii;

import vec.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/** 
 * Image-guided 2D acoustic impedance inversion with constraints from well logs
 * @Author: Xinming Wu and Dave Hale, Colorado School of Mines
 * @Version: 2015.10.10
 */

public class AcousticImpedanceInv3 {

  /**
   * Constructs an unfaulter.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   */
  public AcousticImpedanceInv3(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _sigma3 = (float)sigma2;
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

  public void setTensors(Tensors3 d) {
    _d = d;
  }

  public float[][][] applyForImpedance(float[][][] p,
    float[][][] r, float[][][] wp, float[] k1, float[] k2, float[] k3, float[] f) 
  {
    int n3 = r.length;
    int n2 = r[0].length;
    int n1 = r[0][0].length;
    float[][][] b = new float[n3][n2][n1];
    setInitial(p,k1,k2,k3,f);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vp = new VecArrayFloat3(p);
    CgSolver cs = new CgSolver(_small,_niter);
    A3 a3 = new A3(_d,wp);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,wp,k1,k2,k3);
    vb.zero();
    makeRhs(wp,r,b);
    cs.solve(a3,m3,vb,vp);
    return p;
  }

  public float[][][] initialInterp(
    float[][][] wp, float[] k1, float[] k2, float[] k3, float[] f) 
  {
    float sct = _sc;
    _sc = 1.0f;
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][] b = new float[n3][n2][n1];
    float[][][] p = new float[n3][n2][n1];
    setInitial(p,k1,k2,k3,f);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vp = new VecArrayFloat3(p);
    CgSolver cs = new CgSolver(_small,_niter);
    A3 a3 = new A3(_d,wp);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,wp,k1,k2,k3);
    vb.zero();
    cs.solve(a3,m3,vb,vp);
    _sc = sct;
    return p;
  }

  public void setInitial(
    float[][][] p, float[] k1, float[] k2, float[] k3, float[] f) 
  {
    int np = k1.length;
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    float fa = sum(f)/(float)np;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      p[i3][i2][i1] = fa;
    }}}
    for (int ip=0; ip<np; ++ip) {
      int i1 = round(k1[ip]);
      int i2 = round(k2[ip]);
      int i3 = round(k3[ip]);
      if(i1<0) i1=0; if(i1>n1-1) i1=n1-1;
      if(i2<0) i2=0; if(i2>n2-1) i2=n2-1;
      if(i3<0) i3=0; if(i3>n3-1) i3=n3-1;
      p[i3][i2][i1] = f[ip];
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private Tensors3 _d = null;
  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _sigma3 = 6.0f; // half-width of smoother in 3rd dimension
  private float _small = 0.010f; // stop CG iterations if residuals are small
  private int _niter = 200; // maximum number of inner CG iterations
  private static float _sc = 0.5f;

  private static class A3 implements CgSolver.A {
    A3(Tensors3 et, float[][][] w) { 
      _et=et;
      _w=w;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      VecArrayFloat3 v3z = v3x.clone();
      float[][][] x = v3z.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      v3y.zero();
      applyLhs(_et,_w,z,y);
    }
    private Tensors3 _et = null;
    private float[][][] _w=null;
  }

  private static void applyLhs(
    Tensors3 et, float[][][] w, float[][][] x, float[][][] y) 
  {
    if(_sc<1.0f) {applyLhs1(w,x,y);}
    if(_sc>0.0f) {applyLhs2(et,w,x,y);}
  }

  private static void applyLhs1(
    final float[][][] w, final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    final float sc = (1f-_sc);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
      float[] w32 = w[i3][i2];
      float[] x32 = x[i3][i2];
      float[] y32 = y[i3][i2];
      for (int i1=1; i1<n1; ++i1) {
        float xa=0.0f;
        float wi=w32[i1]*sc;
        float ws=wi*wi; 
        xa  = x32[i1  ];
        xa -= x32[i1-1];
        xa *= ws;
        y32[i1-1] -= xa;
        y32[i1  ]  = xa;
      }
    }}});
  }

  private static void makeRhs(
    final float[][][] w, final float[][][] r, final float[][][] y) 
  {
    zero(y);
    final int n3 = y.length;
    final int n2 = y[0].length;
    final int n1 = y[0][0].length;
    final float sc = (1f-_sc);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
      float[] w32 = w[i3][i2];
      float[] r32 = r[i3][i2];
      float[] y32 = y[i3][i2];
    for (int i1=1; i1<n1; ++i1) {
      float wi = w32[i1]*sc;
      float ws = wi*wi;
      float ri = ws*r32[i1];
      y32[i1  ] += ri;
      y32[i1-1] -= ri;
    }}}});
  }

  private static void applyLhs2(
    final Tensors3 d, final float[][][] wp, 
    final float[][][] x, final float[][][] y)
  { 
    final int n3 = y.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,wp,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,wp,x,y);
    }});
  }

  // 3D LHS
  private static void applyLhsSlice3(
    int i3, Tensors3 d, float[][][] wp, float[][][] x, float[][][] y)
  {
    int n2 = y[0].length;
    int n1 = y[0][0].length;
    float[] di = fillfloat(1.0f,6);
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
        if(d!=null){d.getTensor(i1,i2,i3,di);}
        float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
        float wps = wpi*wpi*_sc;
        float d11 = di[0];
        float d12 = di[1];
        float d13 = di[2];
        float d22 = di[3];
        float d23 = di[4];
        float d33 = di[5];
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

        float x1 = 0.25f*(xa+xb+xc+xd)*wps;
        float x2 = 0.25f*(xa-xb+xc-xd)*wps;
        float x3 = 0.25f*(xa+xb-xc-xd)*wps;

        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2+d23*x3;
        float y3 = d13*x1+d23*x2+d33*x3;

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


  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M3 implements CgSolver.A {
    M3(float sigma1, float sigma2, float sigma3, 
       float[][][] wp, float[] k1, float[] k2, float[] k3)
    {
      _wp = wp;
      _k1 = copy(k1);
      _k2 = copy(k2);
      _k3 = copy(k3);
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      copy(x,y);
      constrain(_k1,_k2,_k3,y);
      smooth3(_sigma3,_wp,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2f*_sigma1,y);
      smooth2(_sigma2,_wp,y);
      smooth3(_sigma3,_wp,y);
      constrain(_k1,_k2,_k3,y);
    }
    private float _sigma1;
    private float _sigma2;
    private float _sigma3;
    private float[][][] _wp;
    private float[] _k1,_k2,_k3;
  }

  private static void constrain(
    float[] k1, float[] k2, float[] k3, float[][][] x) 
  {
    int np = k2.length;
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for (int ip=0; ip<np; ++ip) {
      int i1 = round(k1[ip]);
      int i2 = round(k2[ip]);
      int i3 = round(k3[ip]);
      if(i1<0) i1=0; if(i1>n1-1) i1=n1-1;
      if(i2<0) i2=0; if(i2>n2-1) i2=n2-1;
      if(i3<0) i3=0; if(i3>n3-1) i3=n3-1;
      x[i3][i2][i1] = 0.f;
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

  private static void smooth1(
    final float sigma, final float[][][] x) 
  {
    final int n3 = x.length;
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] x3 = x[i3];
      smooth1(sigma,x3);
    }});
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



}
