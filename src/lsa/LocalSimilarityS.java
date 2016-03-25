package lsa;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Local similarity attribute
 * @author Xinming Wu 
 * @version 2016.03.23
 */

public class LocalSimilarityS {
 
  public LocalSimilarityS(double smin, double smax, double ds) {
    int ns = (int)((smax-smin)/ds)+1;
    _ss = new Sampling(ns,ds,smin);
  }

  // scale the curvature term 
  public void setWeights(float w){
    _weight = w;
  }
  
  public void setSmoothings(float sigma1, float sigma2){
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }

  public float[][] apply(float[] x, float[] y) {
    int n1 = x.length;
    int ns = _ss.getCount();
    float[][] sa = new float[ns][n1];
    double[] s = _ss.getValues();
    for (int is=0; is<ns; ++is) {
      double si = s[is];
      float[] p = apply( 1,si,x,y);
      float[] q = apply(-1,si,y,x);
      mul(p,q,sa[is]);
    }
    return sa;
  }

  public float[] apply(int dir, double s, float[] x, float[] y) {
    int n1 = x.length;
    float[] r = new float[n1];
    float[] b = new float[n1];
    float[][] w = fillfloat(1f,n1,n1);
    VecArrayFloat1 vb = new VecArrayFloat1(b);
    VecArrayFloat1 vr = new VecArrayFloat1(r);
    Smoother2 smoother2 = new Smoother2(_sigma1,_sigma2,w);
    A2 a2 = new A2(smoother2,dir,_weight,w,s,x);
    CgSolver cs = new CgSolver(_small,_niter);
    if(dir>0) {
      makeRhsF(s,x,y,b);
    }else{ 
      makeRhsR(s,x,y,b);
    }
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    smoother2.apply(r);
    return r;
  }

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, int dir, float w1, float[][] wp, double s, float[] a) {
      _s2 = s2;
      _w1 = w1;
      _wp = wp;
      _s  = s;
      _a  = a;
      _dir = dir;
      _sc = 0.001f;
      //testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      zero(y);
      _s2.apply(z);
      scaleAndAdd(-_sc,z,y);
      if(_dir>0) {
        applyLhsF(_s,_a,z,y);
      }else{ 
        applyLhsR(_a,z,y);
      }
      _s2.applyTranspose(y);
      scaleAndAdd(_sc,x,y);
    }
    private Smoother2 _s2;
    private float _w1;
    private float[] _a;
    private double _s;
    private float[][] _wp;
    private int _dir;
    private float _sc;

  }

  private static void scaleAndAdd(float sc, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1) {
      y[i1] += sc*x[i1];
    }
  }


  // Smoother used as a preconditioner. 
  private static class Smoother2 {
    public Smoother2(
      float sigma1, float sigma2, float[][] el) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _el = el;
    }
    public void apply(float[] x) {
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[] x) {
      smooth1(_sigma1,x);
    }
    private float[][] _el;
    private float _sigma1,_sigma2;
  }


  private static void checkControlPoints(float[] k2, float[] k3, float[][] f) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
      }
    }
  }

  private static void constrain(float[] k2, float[] k3, float[][] x) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2] = 0.f;
      }
    }
  }

    // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
  }

  private static void smooth2(float sigma, float[][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }


  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    float[] st = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        xt[i1] = x[i2][i1];
        st[i1] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    float[] st = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        xt[i2] = x[i2][i1];
        st[i2] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }

  private static void applyLhsF(
    double s, float[] a, float[] x, float[] y) {
    zero(y);
    int n1 = x.length;
    float[] as = shift(s,a);
    for (int i1=0; i1<n1; ++i1) {
      y[i1] = as[i1]*as[i1]*x[i1];
    }

    /*
    float w1 = 0.1f;
    float ws = w1*w1;
    float[][] g = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float xa = 0.0f;
      float xb = 0.0f;
      xa += x[i2  ][i1  ];
      xb -= x[i2  ][i1-1];
      xb += x[i2-1][i1  ];
      xa -= x[i2-1][i1-1];
      float x1 = 0.5f*(xa+xb);
      float x2 = 0.5f*(xa-xb);
      float g1 = x1*ws;
      float g2 = x2*ws;
      float ga = 0.5f*(g1+g2);
      float gb = 0.5f*(g1-g2);
      g[i2  ][i1  ] += ga;
      g[i2  ][i1-1] -= gb;
      g[i2-1][i1  ] += gb;
      g[i2-1][i1-1] -= ga;
    }}
    add(g,y,y);
    */
  }

  private static void applyLhsR(
    float[] b, float[] x, float[] y) {
    zero(y);
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = b[i1]*b[i1]*x[i1];

    /*
    float w1 = 0.1f;
    float ws = w1*w1;
    float[][] g = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float xa = 0.0f;
      float xb = 0.0f;
      xa += x[i2  ][i1  ];
      xb -= x[i2  ][i1-1];
      xb += x[i2-1][i1  ];
      xa -= x[i2-1][i1-1];
      float x1 = 0.5f*(xa+xb);
      float x2 = 0.5f*(xa-xb);
      float g1 = x1*ws;
      float g2 = x2*ws;
      float ga = 0.5f*(g1+g2);
      float gb = 0.5f*(g1-g2);
      g[i2  ][i1  ] += ga;
      g[i2  ][i1-1] -= gb;
      g[i2-1][i1  ] += gb;
      g[i2-1][i1-1] -= ga;
    }}
    add(g,y,y);
    */

  }


  private static float[] shift(double s, float[] a) {
    int n1 = a.length;
    float[] as = new float[n1];
    Sampling s1 = new Sampling(n1);
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int i1=0; i1<n1; ++i1) {
      double x1 = i1+s;
      if(x1<0||x1>=n1) {x1=i1-s;}
      as[i1] = si.interpolate(s1,a,x1);
    }
    return as;
  }
  
  //private static void makeRhs
  private static void makeRhsF(
    double s, float[] x, float[] y, float[] b) 
  {
    int n1 = b.length;
    float[] xs = shift(s,x);
    for (int i1=0; i1<n1; ++i1) 
      b[i1] = xs[i1]*y[i1];
  }

  //private static void makeRhs
  private static void makeRhsR(
    double s, float[] x, float[] y, float[] b) 
  {
    int n1 = b.length;
    float[] ys = shift(s,y);
    for (int i1=0; i1<n1; ++i1) 
      b[i1] = ys[i1]*x[i1];
  }

 
  ///////////////////////////////////////////////////////////////////////////
  // private
  private Sampling _ss;
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 10.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 2.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.001f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
}
