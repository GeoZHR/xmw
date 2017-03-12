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

public class LocalSimilarity {
 
  public LocalSimilarity(double smin, double smax, double ds) {
    int ns = (int)((smax-smin)/ds)+1;
    _ss = new Sampling(ns,ds,smin);
  }

  public void setStrain(double strainMax) {
    _bstrain1 = (int)ceil(1.0/strainMax);
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
    float[][] p = apply( 1,x,y);
    float[][] q = apply(-1,y,x);
    return mul(p,q);
  }

  public float[][] accumulateForward(float[][] e) {
    float[][] d = like(e);
    accumulateForward(e,d);
    return d;
  }

  public float[][] accumulateReverse(float[][] e) {
    float[][] d = like(e);
    accumulateReverse(e,d);
    return d;
  }

  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }

  /**
   * Accumulates alignment errors in reverse direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse(float[][] e, float[][] d) {
    accumulate(-1,_bstrain1,e,d);
  }

  private static void accumulate(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = max3(dm,di,dp)+e[ii][il];
      }
    }
  }


  public float[][] apply(int dir, float[] x, float[] y) {
    int n1 = x.length;
    int ns = _ss.getCount();
    float[][] r = new float[ns][n1];
    float[][] b = new float[ns][n1];
    float[][] w = fillfloat(1f,n1,ns);
    double[] s = _ss.getValues();
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
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
    A2(Smoother2 s2, int dir, float w1, float[][] wp, double[] s, float[] a) {
      _s2 = s2;
      _w1 = w1;
      _wp = wp;
      _s  = s;
      _a  = a;
      _dir = dir;
      //_sc = sum(pow(a,2))/a.length;
      _sc = 0.1f;
      //testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
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
    private double[] _s;
    private float[][] _wp;
    private int _dir;
    private float _sc;
    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n1 = _a.length;
      int ns = _s.length;
      float[][] x = sub(randfloat(n1,ns),0.5f);
      float[][] y = sub(randfloat(n1,ns),0.5f);
      float[][] ax = zerofloat(n1,ns);
      float[][] ay = zerofloat(n1,ns);
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
      System.out.println("A2: yax="+yax+" xay="+xay);
      System.out.println("A2: xax="+xax+" yay="+yay);
    }

  }

  private static void scaleAndAdd(float sc, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      y[i2][i1] += sc*x[i2][i1];
    }}
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
    public void apply(float[][] x) {
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,x);
    }
    private float[][] _el;
    private float _sigma1,_sigma2;
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
    double[] s, float[] a, float[][] x, float[][] y) {
    zero(y);
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
      float[] as = shift(s[i2],a);
      for (int i1=0; i1<n1; ++i1) {
        y[i2][i1] = as[i1]*as[i1]*x[i2][i1];
      }
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
    float[] b, float[][] x, float[][] y) {
    zero(y);
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      y[i2][i1] = b[i1]*b[i1]*x[i2][i1];
    }}
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
      //if(x1<0||x1>=n1) {x1=i1-s;}
      as[i1] = si.interpolate(s1,a,x1);
    }
    return as;
  }
  
  //private static void makeRhs
  private static void makeRhsF(
    double[] s, float[] x, float[] y, float[][] b) 
  {
    int n2 = b.length;
    int n1 = b[0].length;
    for (int i2=0; i2<n2; ++i2) {
      float[] xs = shift(s[i2],x);
      for (int i1=0; i1<n1; ++i1) 
          b[i2][i1] = xs[i1]*y[i1];
    }
  }

  //private static void makeRhs
  private static void makeRhsR(
    double[] s, float[] x, float[] y, float[][] b) 
  {
    int n2 = b.length;
    int n1 = b[0].length;
    for (int i2=0; i2<n2; ++i2) {
      float[] ys = shift(s[i2],y);
      for (int i1=0; i1<n1; ++i1) 
          b[i2][i1] = ys[i1]*x[i1];
    }
  }

  private static float max3(float a, float b, float c) {
    return b>=a?(b>=c?b:c):(a>=c?a:c); // if equal, choose b
  }

  private static float[] like(float[] a) {
    return new float[a.length];
  }
  private static float[][] like(float[][] a) {
    return new float[a.length][a[0].length];
  }
  private static float[][][] like(float[][][] a) {
    return new float[a.length][a[0].length][a[0][0].length];
  }


 
  ///////////////////////////////////////////////////////////////////////////
  // private
  private Sampling _ss;
  private int _bstrain1;
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 10.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 10.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.001f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
}
