package aii;

import vec.*;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class AcousticImpedanceInv2 {

  /**
   * Constructs an impedance inverter.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   */
  public AcousticImpedanceInv2(double sigma1, double sigma2) {
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

  public float[][] applyForImpedance(
    float[][] r, float[][] wp, float[] k1, float[] k2, float[] f) 
  {
    int n2 = r.length;
    int n1 = r[0].length;
    float[][] w = copy(wp);
    float[][] b = new float[n2][n1];
    float[][] p = new float[n2][n1];
    setInitial(w,p,k1,k2,f);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vp = new VecArrayFloat2(p);
    CgSolver cs = new CgSolver(_small,_niter);
    A2 a2 = new A2(_d,wp,wp);
    M2 m2 = new M2(_sigma1,_sigma2,wp,k1,k2);
    vb.zero();
    makeRhs(wp,r,b);
    cs.solve(a2,m2,vb,vp);
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
    A2(Tensors2 et,float[][] w1, float[][] w2) { 
      _et=et;
      _w1=w1;
      _w2=w2;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      VecArrayFloat2 v2z = v2x.clone();
      float[][] x = v2z.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      applyLhs(_et,_w1,_w2,z,y);
    }
    private Tensors2 _et = null;
    private float[][] _w1=null;
    private float[][] _w2=null;
  }

  private static void applyLhs(
    Tensors2 et, float[][] w1, float[][] w2, float[][] x, float[][] y) 
  {
    applyLhs1(w1,x,y);
    applyLhs2(et,w2,x,y);
  }

  private static void applyLhs1(float[][] w, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
      float[] w2 = w[i2];
      float[] x2 = x[i2];
      float[] y2 = y[i2];
    for (int i1=1; i1<n1; ++i1) {
      float xa=0.0f;
      float wi=w2[i1]*(1f-_sc);
      float ws=wi*wi; 
      xa  = x2[i1  ];
      xa -= x2[i1-1];
      xa *= ws;
      y2[i1-1] -= xa;
      y2[i1  ]  = xa;
    }}
  }

  private static void makeRhs(
    float[][] w, float[][] r, float[][] y) 
  {
    zero(y);
    int n2 = y.length;
    int n1 = y[0].length;
    for (int i2=0; i2<n2; ++i2) {
      float[] w2 = w[i2];
      float[] r2 = r[i2];
      float[] y2 = y[i2];
    for (int i1=1; i1<n1; ++i1) {
      float wi = w2[i1]*(1f-_sc);
      float ws = wi*wi;
      float ri = ws*r2[i1];
      y2[i1  ] += ri;
      y2[i1-1] -= ri;
    }}
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


  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[][] wp, float[] k1, float[] k2)
    {
      _wp = wp;
      _k1 = copy(k1);
      _k2 = copy(k2);
      _sigma1 = sigma1;
      _sigma2 = sigma2;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k1,_k2,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2f*_sigma1,y);
      smooth2(_sigma2,_wp,y);
      constrain(_k1,_k2,y);
    }
    private float _sigma1;
    private float _sigma2;
    private float[][] _wp;
    private float[] _k1,_k2;
  }

  private static void constrain(float[] k1, float[] k2, float[][] x) {
    int np = k2.length;
    for (int ip=0; ip<np; ++ip) {
      int i1 = (int)k1[ip]; 
      int i2 = (int)k2[ip]; 
      x[i2][i1] = 0.f;
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

  public void setInitial(
    float[][] w, float[][] p, float[] k1, float[] k2, float[] f) 
  {
    int np = k1.length;
    int n2 = p.length;
    int n1 = p[0].length;
    for (int ip=0; ip<np; ++ip) {
      int i1 = round(k1[ip]);
      int i2 = round(k2[ip]);
      if(i1<0) i1=0; if(i1>n1-1) i1=n1-1;
      if(i2<0) i2=0; if(i2>n2-1) i2=n2-1;
      w[i2][i1] = 0.0f;
      p[i2][i1] = f[ip];
    }
  }


}
