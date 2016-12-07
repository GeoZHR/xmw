/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package acm;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Compute gradient vector flow from a given image.
 * <p>
 * Reference: Snakes, Shapes, and Gradient Vector Flow by Chenyang Xu, 1998
 * The gradient vector flow in this class is efficiently solved using a CG method.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.22
 */

public class GradientVectorFlow {

   /**
   * Sets half-widths for smoothings in all dimensions.
   * These smoothings serve as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute mappings.
   * @param sigma half-width for smoothing in all dimensions, in samples.
   */
  public void setSmoothing(double sigma) {
    _sigma = (float)sigma;
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

  public void setScale(double scale) {
    _scale = (float)scale;
  }

  public float[][] applyForEdge(
    float sigma, float sigma1, float sigma2, float[][] f) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply1X(f,g1);
    rgf.applyX1(f,g2);
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    lof.applyForNormal(f,u1,u2);
    EigenTensors2 et = lof.applyForTensors(f);
    u1 = mul(g1,u1);
    u2 = mul(g2,u2);
    float[][] dg = abs(add(u1,u2));
    dg = smooth(8.0,et,dg);
    dg = div(dg,max(abs(dg)));
    return dg;
  }

  public EigenTensors2 applyForGVX(float sigma, float sigma1, float sigma2, 
    float[][] f, float[][][] up, float[][] wp) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] eu = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply1X(f,g1);
    rgf.applyX1(f,g2);
    up[0] = g1;
    up[1] = g2;

    LocalOrientFilterK lof = new LocalOrientFilterK(sigma1,sigma2);
    lof.applyForNormalLinear(f,u1,u2,eu);
    float[][] gu1 = mul(g1,u1);
    float[][] gu2 = mul(g2,u2);
    float[][] gus = abs(add(gu1,gu2));
    gus = div(gus,max(gus));
    for(int i2=0; i2<n2; ++i2) {
    for(int i1=0; i1<n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      float gui = gus[i2][i1];
      float gsi = sqrt(g1i*g1i+g2i*g2i);
      if(gsi>0.0f) {
        up[0][i2][i1] = gui*g1i/gsi;
        up[1][i2][i1] = gui*g2i/gsi;
      }
    }}
    copy(eu,wp);
    EigenTensors2 et = lof.applyForTensors(f);
    float[][] ev = sub(1.0f,eu);
    et.setEigenvalues(eu,ev);
    return et;
  }


  public EigenTensors2 applyForGV(
    float sigma, float sigma1, float sigma2, float[][] f, float[][][] up) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] eu = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply1X(f,g1);
    rgf.applyX1(f,g2);
    up[0] = g1;
    up[1] = g2;
    /*
    float[][] gs = sqrt(add(pow(g1,2f),pow(g2,2f)));
    float gsm = max(gs);
    for(int i2=0; i2<n2; ++i2) {
    for(int i1=0; i1<n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      up[0][i2][i1] = g1i/gsm;
      up[1][i2][i1] = g2i/gsm;
    }}
    */



    LocalOrientFilterK lof = new LocalOrientFilterK(sigma1,sigma2);
    lof.applyForNormalLinear(f,u1,u2,eu);
    EigenTensors2 et = lof.applyForTensors(f);
    float[][] ev = sub(1.0f,eu);
    et.setEigenvalues(eu,ev);
    return et;
  }

  public float[][][] applyForGradient(float sigma, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] gs = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply10(fx,g1);
    rgf.apply01(fx,g2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      gs[i2][i1] = g1i*g1i+g2i*g2i;
    }}
    return null;
  }

  public float[][][] applyForGVF(
    EigenTensors2 et, float[][] u1, float[][] u2, float[][] wp, float[][] ws) 
  {
    int n2 = u1.length;
    int n1 = u1[0].length;
    float[][][] v = new float[2][n2][n1];
    v[0] = linearSolver(et,u1,wp,ws);
    v[1] = linearSolver(et,u2,wp,ws);
    return v;
  }

  /*
  public float[][][][] applyForGVF(
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] wp) 
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float[][][][] v = new float[3][n3][n2][n1];
    v[0] = linearSolver(u1,wp);
    v[1] = linearSolver(u2,wp);
    v[2] = linearSolver(u3,wp);
    return v;
  }
  */

  public float[][] linearSolver(
    EigenTensors2 et, float[][] u, float[][] wp, float[][] ws) 
  {
    int n2 = u.length;
    int n1 = u[0].length;
    float[][] r = copy(u);
    float[][] b = new float[n2][n1];
    Smoother smoother = new Smoother(_sigma,ws);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    CgSolver cs = new CgSolver(_small,_niter);
    A2 a2 = new A2(_scale,et,smoother,wp);
    makeRhs(wp,u,b);
    smoother.applyTranspose(b);
    cs.solve(a2,vb,vr);
    smoother.apply(r);
    return r;
  }

  /*
  public float[][][] linearSolver(float[][][] u, float[][][] wp) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    float[][][] r = new float[n3][n2][n1];
    float[][][] b = new float[n3][n2][n1];
    Smoother smoother = new Smoother(_sigma);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    CgSolver cs = new CgSolver(_small,_niter);
    A3 a3 = new A3(smoother,wp);
    makeRhs(wp,u,b);
    smoother.applyTranspose(b);
    cs.solve(a3,vb,vr);
    smoother.apply(r);
    return r;
  }
  */

  public float[][] smooth(
    double sigma, EigenTensors2 d, float[][] g) {
    int n1 = g[0].length;
    int n2 = g.length;
    d.setEigenvalues(0.05f,1.00f);
    float c = (float)(0.5*sigma*sigma);
    float[][] h = new float[n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,c,g,h);
    return h;
  }



  ////////////////////////////////////////////////////////////////
  //private

  private float _scale = 1.0f;
  private float _sigma = 8.0f;
  private float _small = 0.01f;
  private int _niter = 200;

  private static class A2 implements CgSolver.A {
    A2(float scale, EigenTensors2 et, Smoother smoother, float[][] wp) {
      _et = et;
      _wp = wp;
      _scale = scale;
      _smoother = smoother;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      applyLhs(_scale,_et,_wp,z,y);
      _smoother.applyTranspose(y);
    }

    private float _scale;
    private float[][] _wp=null;
    private Smoother _smoother;
    private EigenTensors2 _et = null;
  }

  private static class A3 implements CgSolver.A {
    A3(Smoother smoother, float[][][] wp) {
      _wp = wp;
      _smoother = smoother;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      v3y.zero();
      applyLhs(_wp,z,y);
      //_smoother.applyTranspose(y);
    }

    private float[][][] _wp=null;
    private Smoother _smoother;
  }


  // Smoother used as a preconditioner. 
  private static class Smoother {
    public Smoother(float sigma, float[][] wp) {
      _wp = wp;
      _sigma = sigma;
    }

    public void apply(float[][] x) {
      smooth2(_sigma,_wp,x);
      smooth1(_sigma,_wp,x);
    }

    public void apply(float[][][] x) {
      smooth(_sigma,x);
    }

    public void applyTranspose(float[][] x) {
      smooth1(_sigma,_wp,x);
      smooth2(_sigma,_wp,x);
    }

    public void applyTranspose(float[][][] x) {
      smooth(_sigma,x);
    }
    private float _sigma;
    private float[][] _wp;
  }

  private static void smooth(float sigma, float[][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
  }

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



  private static void smooth(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
  }


  private static void applyLaplace(
    final Tensors2 et, final float[][] x, final float[][] y) {
    zero(y);
    int n2 = x.length;
    int n1 = x[0].length;
    float[] ds = new float[3];
    ds[0] = 1.0f;
    ds[1] = 0.0f;
    ds[2] = 1.0f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if(et!=null){et.getTensor(i1,i2,ds);}
        float d11 = ds[0];
        float d12 = ds[1];
        float d22 = ds[2];
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
  private static void applyLhs(final float scale,final Tensors2 et,
    final float[][] wp, final float[][] x, final float[][] y)
  {
    zero(y);
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] t = new float[n2][n1];
    applyLaplace(et,x,t);
    applyLaplace(et,t,y);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wpi = (wp!=null)?wp[i2][i1]:1.0f;
      float wps = wpi*wpi;
      y[i2][i1] = scale*y[i2][i1]-wps*x[i2][i1];
    }}
  }

  private static void applyLhs(
    final float[][][] wp, final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,x,y);
    }});
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
      float wps = wpi*wpi;
      y[i3][i2][i1] += wps*x[i3][i2][i1];
    }}}

  }
  private static void applyLhsSlice3(
    int i3, float[][][] wp, float[][][] x, float[][][] y) 
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
        float wpm = (1.0f-wpi)*(1.0f-wpi);
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
        float y1 = wpm*x1;
        float y2 = wpm*x2;
        float y3 = wpm*x3;
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

  private static void makeRhs(float[][] wp, float[][] u, float[][] y) {
    int n2 = u.length;
    int n1 = u[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wpi = (wp!=null)?wp[i2][i1]:1.0f;
      float wps = wpi*wpi;
      y[i2][i1] = -wps*u[i2][i1];
    }}
  }

  private static void makeRhs(float[][][] wp, float[][][] u, float[][][] y) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
      float wps = wpi*wpi;
      y[i3][i2][i1] = wps*u[i3][i2][i1];
    }}}
  }


}
