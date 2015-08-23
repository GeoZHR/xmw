/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package func;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Flattens and unflattens locally linear features in a 2D image.
 * <p>
 * We assume that a 2D image f(x1,x2) has coherent features that are locally
 * linear with finite (non-vertical) slope. This class computes and applies
 * coordinate mappings for flattening and unflattening such features.
 * In this version, the mappings change only the vertical coordinate.
 * <p> 
 * Let x1 denote the vertical coordinate for the image f(x1,x2). This class
 * computes a mapping function x1(u1,u2) such that features in a transformed
 * image g(u1,u2) = f(x1(u1,u2),u2) are approximately horizontal. This process
 * is often called "flattening", and the transformed image g(u1,u2) is "the
 * flattened image." Likewise, for any constant u1, the curve defined by the
 * function x1(u1,x2) is called a "horizon".
 * <p>
 * If the coordinates u1 and u2 are sampled finely enough, then the mapping
 * function x1(u1,u2) is invertible. Specifically, there exists an inverse
 * mapping u1(x1,x2) such that x1 = x1(u1(x1,x2),x2) for all x2. Currently,
 * samplings of u1 and u2 are set to be identical to those for x1 and x2. If
 * used directly, this sampling of the mapping x1(u1,u2) may cause aliasing of
 * the flattened image g(u1,u2), so that f(x1,x2) cannot be recovered by
 * unflattening.
 * <p>
 * Without additional constraints, the description above of the mapping
 * function x1(u1,u2) is ambiguous. For example, one could add any constant to
 * this mapping function, and features in the transformed image g(u1,u2) would
 * still be horizontal. The mapping function x1(u1,u2) is here constrained so
 * that x1(u1,u2) = u1 for the first and last sampled values of u1. In other
 * words, features at the top or bottom of an image are not shifted by
 * flattening or unflattening.
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.01.29
 */
public class Rgt2 {

  /**
   * Sets the relative weight of the PDE dr(x1,x2)/dx1 = 0.
   * Increasing this weight will cause shifts r(x1,x2) and s(u1,u2) to vary
   * more slowly with vertical coordinates x1 and u1, respectively. A weight
   * of 1.0 will cause this equation to get as much weight as other PDEs that
   * cause contours of constant u1 = u1(x1,x2) to be aligned with coherent
   * linear image features.
   * @param w1 the weight.
   */
  public void setWeight1(double w1) {
    _weight1 = (float)w1;
  }

  /**
   * Sets half-widths for smoothings in 1st and 2nd dimensions.
   * These smoothings serve as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute mappings.
   * @param sigma1 half-width for smoothing in 1st dimension, in samples.
   * @param sigma2 half-width for smoothing in 2nd dimension, in samples.
   */
  public void setSmoothings(double sigma1, double sigma2) {
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

  /**
   * Gets mappings computed from specified slopes and linearities.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param p2 array of slopes of image features.
   * @param el array of linearities of image features.
   */
  public void updateRgtFromNormals(
    float[][] u1, float[][] u2, float[][] el, float[][] r) 
  {
    // Sampling parameters.
    int n2 = u1.length;
    int n1 = u1[0].length;
    // Compute shifts r(x1,x2), in samples.
    float[][] b = new float[n2][n1]; // right-hand side
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(n1,n2,_sigma1,_sigma2,el);
    A2 a2 = new A2(smoother2,_weight1,el);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(el,u1,u2,b);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    smoother2.apply(r);

    float[][] s = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      s[i2][i1] = r[i2][i1]-i1;
    }}
    cleanShifts(s);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      r[i2][i1] = s[i2][i1]+i1;
    }}
  }

  public void setInitial(float[][] u) {
    int n2 = u.length;
    int n1 = u[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      u[i2][i1] = i1;
    }}
  }

  public float[][] flatten(float[][] u, float[][] f) {
    int n2 = u.length;
    int n1 = u[0].length;
    Sampling s1 = new Sampling(n1);
    float[][] x = new float[n2][n1];
    float[][] g = new float[n2][n1];
    InverseInterpolator ii = new InverseInterpolator(s1,s1);
    for (int i2=0; i2<n2; ++i2) 
      ii.invert(u[i2],x[i2]);
    SincInterpolator si = new SincInterpolator();
    for (int i2=0; i2<n2; ++i2)
      si.interpolate(n1,1.0,0.0,f[i2],n1,x[i2],g[i2]);
    return g;
  }

  public float[][] applyForG1(float[][] u) {
    int n2 = u.length;
    int n1 = u[0].length;
    float[][] g1 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      g1[i2][i1] = u[i2][i1]-u[i2][i1-1];
    }}
    for (int i2=0; i2<n2; ++i2) {
      g1[i2][0] = g1[i2][1];
    }
    return g1;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _weight1 = 0.010f; // weight for dr/d1 = 0 equation
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 1000; // maximum number of CG iterations

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float w1, float[][] wp) {
      _s2 = s2;
      _w1 = w1;
      _wp = wp;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      _s2.apply(z);
      zero(y);
      applyLhs(_w1,_wp,z,y);
      _s2.applyTranspose(y);
    }
    private Smoother2 _s2;
    private float _w1;
    private float[][] _wp;
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(
      int n1, int n2, float sigma1, float sigma2, float[][] el) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _el = el;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,_el,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,_el,x);
    }
    private float _sigma1,_sigma2;
    private float[][] _el;
    private void zero1(float[][] x) {
      int n1 = x[0].length;
      int n2 = x.length;
      for (int i2=0; i2<n2; ++i2) {
        x[i2][   0] = 0.0f;
        x[i2][n1-1] = 0.0f;
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

  private static void makeRhs(
    float[][] wp, float[][] u1, float[][] u2, float[][] y) 
  {
    int n2 = y.length;
    int n1 = y[0].length;
    for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float wpi = (wp!=null)?wp[i2][i1]:1.0f;
        float wps = wpi*wpi;
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float y1 = wps*u1i;
        float y2 = wps*u2i;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2 ][i1 ] += ya;
        y[i2 ][i1m] -= yb;
        y[i2m][i1 ] += yb;
        y[i2m][i1m] -= ya;
      }
    }
  }

  private static void applyLhs(
    float w1, float[][] wp, 
    float[][] x, float[][] y)
  {
    int n2 = x.length;
    int n1 = x[0].length;
    float w1s = w1*w1;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.0f;
        float wps = wpi*wpi;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = x1*(wps+w1s);
        float y2 = x2*wps;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Post-processing of computed shifts to ensure monotonic u1.
  private static void cleanShifts(float[][] r) {
    int n1 = r[0].length;
    int n2 = r.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if (r[i2][i1]<=r[i2][i1-1]-0.99f)
          r[i2][i1] = r[i2][i1-1]-0.99f;
      }
    }
  }

  private static void printStats(String s, int i1, float[][] a) {
    int n2 = a.length;
    float amin = a[0][i1];
    float amax = a[0][i1];
    for (int i2=1; i2<n2; ++i2) {
      if (a[i2][i1]<amin) amin = a[i2][i1];
      if (a[i2][i1]>amax) amax = a[i2][i1];
    }
    System.out.println(s+": i1="+i1+" min="+amin+" max="+amax);
  }
}
