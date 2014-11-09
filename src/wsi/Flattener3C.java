/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wsi;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Flattens and unflattens locally planar features in a 3D image.
 * <p>
 * We assume that a 3D image f(x1,x2,x3) has coherent features that are
 * locally planar with finite (non-vertical) slope. This class computes and
 * applies coordinate mappings for flattening and unflattening such features.
 * In this version, the mappings change only the vertical coordinate.
 * <p> 
 * Let x1 denote the vertical coordinate for the image f(x1,x2,x3). This class
 * computes a mapping function x1(u1,u2,u3) such that features in a
 * transformed image g(u1,u2,u3) = f(x1(u1,u2),u2,u3) are approximately
 * horizontal. This process is often called "flattening", and the transformed
 * image g(u1,u2,u3) is "the flattened image." Likewise, for any constant u1,
 * the curve defined by the function x1(u1,x2,x3) is called a "horizon".
 * <p>
 * If the coordinates u1, u2 and u3 are sampled finely enough, then the
 * mapping function x1(u1,u2,u3) is invertible. Specifically, there exists an
 * inverse mapping u1(x1,x2,x3) such that x1 = x1(u1(x1,x2,x3),x2,u3) for all
 * x2 and x3. Currently, samplings of u1, u2 and u3 are set to be identical to
 * those for x1, x2 and x3. If used directly, this sampling of the mapping
 * x1(u1,u2,u3) may cause aliasing of the flattened image g(u1,u2,u3), so that
 * f(x1,x2,x3) cannot be recovered by unflattening.
 * <p>
 * Without additional constraints, the description above of the mapping
 * function x1(u1,u2) is ambiguous. For example, one could add any constant to
 * this mapping function, and features in the transformed image g(u1,u2) would
 * still be horizontal. The mapping function x1(u1,u2) is here constrained so
 * that x1(u1,u2) = u1 for the first and last sampled values of u1. In other
 * words, features at the top or bottom of an image are not shifted by
 * flattening or unflattening.
 * @author Dave Hale and Xinming Wu, Colorado School of Mines
 * @version 2013.01.29
 */
public class Flattener3C {

  /** Coordinate mappings u1(x1,x2,x3) and x1(u1,u2,u3). */
  public static class Mappings {
    
    /** Sampling for the 1st dimension (the vertical coordinate). */
    public Sampling s1;
    
    /** Sampling for the 2nd dimension. */
    public Sampling s2;
    
    /** Sampling for the 3rd dimension. */
    public Sampling s3;

    /** Array of sampled u1(x1,x2,x3). */
    public float[][][] u1;
    
    /** Array of sampled x1(u1,u2,u3). */
    public float[][][] x1;

    /** Array of sampled sf(x1,x2,x3). */
    public float[][][] sf;

    /**
     * Uses these mappings to flatten the specified image.
     * @param f the image to flatten.
     * @return the flattened image.
     */
    public float[][][] flatten(float[][][] f) {
      return apply(x1,f);
    }

    /**
     * Uses these mappings to unflatten the specified image.
     * @param f the image to unflatten.
     * @return the unflattened image.
     */
    public float[][][] unflatten(float[][][] f) {
      return apply(u1,f);
    }

    /**
     * Gets the flattening shifts s(u1,u2,u3) = u1 - x1(u1,u2,u3).
     * @return the flattening shifts.
     */
    public float[][][] getShiftsS() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n3 = s3.getCount();
      float[][][] s = new float[n3][n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float u1i = f1+i1*d1;
            s[i3][i2][i1] = u1i-x1[i3][i2][i1];
          }
        }
      }
      return s;
    }

    /**
     * Gets the unflattening shifts r(x1,x2,x3) = u1(x1,x2,x3) - x1.
     * @return the unflattening shifts.
     */
    public float[][][] getShiftsR() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int n3 = s3.getCount();
      float[][][] r = new float[n3][n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            float x1i = f1+i1*d1;
            r[i3][i2][i1] = u1[i3][i2][i1]-x1i;
          }
        }
      }
      return r;
    }

    private Mappings(
      Sampling s1, Sampling s2, Sampling s3, 
      float[][][] u1, float[][][] x1, float[][][] r) 
    {
      this.s1 = s1;
      this.s2 = s2;
      this.s3 = s3;
      this.u1 = u1;
      this.x1 = x1;
      this.sf = r ;
    }

    private float[][][] apply(final float[][][] ux, final float[][][] f) {
      final int n1 = s1.getCount();
      final int n2 = s2.getCount();
      final int n3 = s3.getCount();
      final double d1 = s1.getDelta();
      final double f1 = s1.getFirst();
      final SincInterpolator si = new SincInterpolator();
      final float[][][] g = new float[n3][n2][n1];
      Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2)
          si.interpolate(n1,d1,f1,f[i3][i2],n1,ux[i3][i2],g[i3][i2]);
      }});
      return g;
    }
  }

  /**
   * Sets the relative weight of the PDE dr(x1,x2)/dx1 = 0.
   * Increasing this weight will cause shifts r(x1,x2) and s(u1,u2) to vary
   * more slowly with vertical coordinates x1 and u1, respectively. A weight
   * of 1.0 will cause this equation to get as much weight as other PDEs that
   * cause contours of constant u1 = u1(x1,x2) to be aligned with coherent
   * planar image features.
   * @param w1 the weight.
   */
  public void setWeight1(double w1) {
    _weight1 = (float)w1;
  }
  
  public void setScale(double sc) {
    _scale = (float)sc;
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
   * Gets mappings computed from specified slopes and planarities.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param p2 array of slopes of image features.
   * @param wp array of weights for slopes.
   * @param ws array of weights for smoothings.
   */
  public Mappings getMappingsFromSlopes(
    Sampling s1, Sampling s2, Sampling s3,
    float[][][] p2, float[][][] p3, float[][][] wp, float[][][] u) 
  {
    // Sampling parameters.
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    float d3 = (float)s3.getDelta();
    float f1 = (float)s1.getFirst();
    // If necessary, convert units for slopes to samples per sample.
    if (d1!=d2)
      p2 = mul(d2/d1,p2);
    if (d1!=d3)
      p3 = mul(d3/d1,p3);
    // Compute shifts r(x1,x2,x3), in samples.
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] r = new float[n3][n2][n1]; // shifts, in samples
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    //setWeightsForSmoothing(wp);
    //int[][] c = referenceTrace(mul(u,wp));
    int c = 0;
    //int c = referenceTrace1(mul(u,wp));
    A3 a3 = new A3(_scale,_weight1,wp,p2,p3);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,wp,c);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(wp,p2,p3,b);
    cs.solve(a3,m3,vb,vr);
    cleanShifts(r);
    final float[][][] u1 = r;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = f1+i1*d1;
          u1[i3][i2][i1] = x1i+r[i3][i2][i1]*d1;
        }
      }
    }

    // Compute x1(u1,u2).
    final float[][][] x1 = b;
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(u1[i3][i2],x1[i3][i2]);
    }});
    System.out.println("minX1"+min(x1));
    System.out.println("maxX1"+max(x1));

    return new Mappings(s1,s2,s3,u1,x1,r);
  }

  private int[][] referenceTrace(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    int[][] c = new int[2][n2];
    for (int i2=0; i2<n2; ++ i2) {
      float[][] x2 = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) 
        x2[i3] = x[i3][i2];
      c[1][i2] = i2;
      c[0][i2] = referenceTrace(x2);
    }
    return c;
  }

  private int referenceTrace(float[][] x) {
    int c = 0;
    int n3 = x.length;
    float sum = 0.0f;
    for (int i3=0; i3<n3; ++i3) {
      float smi = sum(abs(x[i3]));
      if (smi>sum) {c = i3;sum = smi;}
    }
    return c;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _scale = 0.0f;    // scale the curvature term
  private float _weight1 = 0.01f; // weight for dr/d1 = 0 equation
  private float _sigma1 = 4.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 4.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 4.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 1000; // maximum number of CG iterations

  private float[][][] setWeightsForSmoothing(float[][][] wp) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][] ws = copy(wp);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float wpi = wp[i3][i2][i1];
          if (wpi<0.00005f) {
            ws[i3][i2][i1] = 0.00005f;
          }
        }
      }
    }
    return ws;
  }


  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(float sc, float w1, float[][][] wp, 
       float[][][] p2, float[][][] p3) 
    {
      _wp = wp;
      _p2 = p2;
      _p3 = p3;
      _sc = sc;
      _w1 = w1;
      //testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      float[][][] y1 = copy(x);
      float[][][] y2 = copy(x);
      VecArrayFloat3 v3yy = new VecArrayFloat3(y2);
      v3yy.zero();
      zero(y);
      applyLhs(_w1,_wp,_p2,_p3,z,y);
      if (_sc>0.0f) {
        applyLaplace(_wp,z ,y1);
        applyLaplace(_wp,y1,y2);
        v3y.add(1.f,v3yy,_sc);
      }
    }
    private float _sc;
    private float _w1;
    private float[][][] _wp;
    private float[][][] _p2;
    private float[][][] _p3;
    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n1 = _wp[0][0].length;
      int n2 = _wp[0].length;
      int n3 = _wp.length;
      float[][][] x = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] y = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] ax = zerofloat(n1,n2,n3);
      float[][][] ay = zerofloat(n1,n2,n3);
      VecArrayFloat3 vx = new VecArrayFloat3(x);
      VecArrayFloat3 vy = new VecArrayFloat3(y);
      VecArrayFloat3 vax = new VecArrayFloat3(ax);
      VecArrayFloat3 vay = new VecArrayFloat3(ay);
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

  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M3 implements CgSolver.A {
    M3(float sigma1, float sigma2, float sigma3, float[][][] wp, 
       int c) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      _wp = wp;
      _c = c;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      copy(x,y);
      constrain(_c,y);
      //removeAverage(y);
      smooth3(_sigma3,_wp,y);
      smooth2(_sigma2,_wp,y);
      //smooth1(2.0f*_sigma1,_wp,y);
      smooth1(2.0f*_sigma1,y);
      smooth2(_sigma2,_wp,y);
      smooth3(_sigma3,_wp,y);
      //removeAverage(y);
      constrain(_c,y);
    }
    private float _sigma1,_sigma2,_sigma3;
    private float[][][] _wp;
    private int _c;
  }

  public static void constrain(int[][] c, float[][][] x) {
    int nc = c[0].length;
    int n1 = x[0][0].length;
    for (int ic=0; ic<nc; ++ic) {
      int i3 = c[0][ic];
      int i2 = c[1][ic];
      for (int i1=0; i1<n1; ++i1) {
        x[i3][i2][i1] = 0.0f;
      }
    }
  }

  public static void constrain(int c, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        x[i3][i2][c] = 0.0f;
      }
    }
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
  private static void smooth1(final float sigma, final float[][][] s, final float[][][] x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] x3 = x[i3];
      float[][] s3 = (s!=null)?s[i3]:null;
      smooth1(sigma,s3,x3);
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

  private static void removeAverage(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float nh = (float)(n2*n3);
    for (int i1=0; i1<n1; ++i1) {
      float sumx = 0.0f;
      for (int i3=0; i3<n3; ++i3)  
        for (int i2=0; i2<n2; ++i2)  
          sumx += x[i3][i2][i1];
      float avgx = sumx/nh;
      for (int i3=0; i3<n3; ++i3) 
        for (int i2=0; i2<n2; ++i2) 
          x[i3][i2][i1] -= avgx; 
    }
  }

  private static void makeRhs(
    float[][][] wp, float[][][] p2, float[][][] p3, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1,i3m=0; i3<n3; ++i3,++i3m) {
      for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
        for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
          float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
          if(wpi<0.01f) {wpi=0.01f;}
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float b12 = wpi*p2i;
          float b13 = wpi*p3i;
          float b22 = wpi;
          float b33 = wpi;
          float x2 = -wpi*p2i;
          float x3 = -wpi*p3i;
          float y1 = b12*x2+b13*x3;
          float y2 = b22*x2;
          float y3 = b33*x3;
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

  private static void applyLaplace(
    final float[][][] wp, final float[][][]x, final float[][][] y)
  {   
    final int n3 = x.length; 
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLaplaceSlice3(i3,wp,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLaplaceSlice3(i3,wp,x,y);
    }});
  }
  private static void applyLaplaceSlice3(
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
        float d11 = 1.0f;
        float d22 = 1.0f;
        float d33 = 1.0f;
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
        float y1 = d11*x1;
        float y2 = d22*x2;
        float y3 = d33*x3;
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

  private static void applyLhs(
    final float w1, final float[][][] wp, 
    final float[][][] p2, final float[][][] p3,
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,w1,wp,p2,p3,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,w1,wp,p2,p3,x,y);
    }});
  }
  private static void applyLhsSlice3(
    int i3, float w1,
    float[][][] wp, float[][][] p2, float[][][] p3,
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float w1s = w1*w1;
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
        if(wpi<0.01f) {wpi=0.01f;}
        float p2i = p2[i3][i2][i1];
        float p3i = p3[i3][i2][i1];
        float wps = wpi*wpi;
        float p2s = p2i*p2i;
        float p3s = p3i*p3i;
        float d11 = w1s+wps*(p2s+p3s);
        float d12 = wps*p2i;
        float d13 = wps*p3i;
        float d22 = wps;
        float d33 = wps;
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
        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2       ;
        float y3 = d13*x1       +d33*x3;
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
  private static void trace(String s) {
    System.out.println(s);
  }
  // Post-processing of computed shifts to ensure monotonic u1.
  private static void cleanShifts(float[][][] r) {
    int n1 = r[0][0].length;
    int n2 = r[0].length;
    int n3 = r.length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          if (r[i3][i2][i1]<=r[i3][i2][i1-1]-0.99f)
            r[i3][i2][i1] = r[i3][i2][i1-1]-0.99f;
        }
      }
    }
  }
}
