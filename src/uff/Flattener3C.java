/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package uff;

import vec.*;
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
      float[][][] u1, float[][][] x1) 
    {
      this.s1 = s1;
      this.s2 = s2;
      this.s3 = s3;
      this.u1 = u1;
      this.x1 = x1;
    }

    private float[][][] apply(final float[][][] ux, final float[][][] f) {
      final int n1 = s1.getCount();
      final int n2 = s2.getCount();
      final int n3 = s3.getCount();
      final double d1 = s1.getDelta();
      final double f1 = s1.getFirst();
      final SincInterp si = new SincInterp();
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
  public Mappings getMappingsFromShifts(
    Sampling s1, Sampling s2, Sampling s3, float[][][] r)
  {
    // Sampling parameters.
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    float d1 = (float)s1.getDelta();
    float f1 = (float)s1.getFirst();
    
    // Compute u1(x1,x2,x3).
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
    final float[][][] x1 = new float[n3][n2][n1];
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(u1[i3][i2],x1[i3][i2]);
    }});

    return new Mappings(s1,s2,s3,u1,x1);
  }

  public void computeShifts(
    float[][][] p2, float[][][] p3, float[][][] wp, 
    float[][][] cs, float[][][] r) 
  {
    int n3 = p2.length;
    int n2 = p2[0].length;
    int n1 = p2[0][0].length;
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    initializeShifts(cs,r); // initial shifts to satisfy constraints
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    //float[][][] w1 = setWeightsFromUnconformities(wp,uc);
    float[][][] ws = setWeightsForSmoothing(wp);
    A3 a3 = new A3(ws,p2,p3);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,ws,cs);
    CgSolver cg = new CgSolver(_small,_niter);
    makeRhs(wp,p2,p3,b);
    cg.solve(a3,m3,vb,vr);
    //checkShifts(cs[0],cs[1],cs[2],r);
    cleanShifts(r);
  }

  private float[][][] setWeightsFromUnconformities(float[][][] wp, int[][][] uc) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    if (uc==null) {
      return fillfloat(0.01f,n1,n2,n3);
    } else {
      int nc = uc[0].length;
      float[][][] w1 = fillfloat(0.06f,n1,n2,n3);
      float[][][] mk = fillfloat(0.00f,n1,n2,n3);
      for (int ic=0; ic<nc; ++ic) {
        int np = uc[0][ic].length;
        for (int ip=0; ip<np; ++ip) {
          int i1 = uc[0][ic][ip];
          int i2 = uc[1][ic][ip];
          int i3 = uc[2][ic][ip];
          wp[i3][i2][i1] = 0.0005f;
          w1[i3][i2][i1] = 0.0001f;
          mk[i3][i2][i1] = 1.0000f;
          int i1m = i1-1;
          int i1p = i1+1;
          if(i1m>=0) {
            wp[i3][i2][i1m] = 0.0025f;
            w1[i3][i2][i1m] = 0.0005f;
            mk[i3][i2][i1m] = 1.0000f;
          }
          if(i1p<n1) {
            wp[i3][i2][i1p] = 0.0025f;
            w1[i3][i2][i1p] = 0.0005f;
            mk[i3][i2][i1p] = 1.0000f;
          }
        }
      }
      extend(wp,w1,mk);
      return w1;
    }
  }

  private void extend(float[][] wp, float[][] w1, float[][] mk) {
    int n2 = wp.length;
    int n1 = wp[0].length;
    for (int i2=2; i2<n2-2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float mki = mk[i2][i1];
        float wpi = wp[i2][i1];
        float w1i = w1[i2][i1];
        if(mki==1.0f){
          wp[i2-1][i1] = wpi; 
          wp[i2-2][i1] = wpi; 
          wp[i2+1][i1] = wpi; 
          wp[i2+2][i1] = wpi; 
          w1[i2-1][i1] = w1i;
          w1[i2-2][i1] = w1i;
          w1[i2+1][i1] = w1i;
          w1[i2+2][i1] = w1i;
        }
      }
    }
  }

  private void extend(float[][][] wp, float[][][] w1, float[][][] mk) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    for (int i3=0; i3<n3; ++i3) 
      extend(wp[i3],w1[i3],mk[i3]);
    float[][] wpi2  = new float[n3][n1];
    float[][] w1i2  = new float[n3][n1];
    float[][] mki2  = new float[n3][n1];
    for (int i2=0; i2<n2; ++i2) { 
      for (int i3=0; i3<n3; ++i3) {
        wpi2[i3]  = wp[i3][i2];
        w1i2[i3]  = w1[i3][i2];
        mki2[i3]  = mk[i3][i2];
      }
      extend(wpi2,w1i2,mki2);
      for (int i3=0; i3<n3; ++i3) { 
        wp[i3][i2] = wpi2[i3];
        w1[i3][i2] = w1i2[i3];
      }
    }
  }


  private float[][][] setWeightsForSmoothing(float[][][] wp) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][] ws = copy(wp);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float wpi = wp[i3][i2][i1];
          if (wpi<0.0005f) {
            ws[i3][i2][i1] = 0.0005f;
          }
        }
      }
    }
    return ws;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1 = 4.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 4.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 4.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 1000; // maximum number of CG iterations

  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(float[][][] wp, float[][][] p2, float[][][] p3) 
    {
      _wp = wp;
      _p2 = p2;
      _p3 = p3;
      //testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      zero(y);
      applyLhs(_wp,_p2,_p3,z,y);
    }
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
    M3(float sigma1, float sigma2, float sigma3, 
       float[][][] wp, float[][][] cs) 
    {
      _wp = wp;
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      if (cs!=null) {
        _cs = copy(cs);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      copy(x,y);
      constrain(_cs,y);
      removeAverage(y);
      smooth3(_sigma3,_wp,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.0f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      smooth3(_sigma3,_wp,y);
      removeAverage(y);
      constrain(_cs,y);
    }
    private float _sigma1,_sigma2,_sigma3;
    private float[][][] _wp;
    private float[][][] _cs;
  }

  public static void initializeShifts(float[][][] cs, float[][][] r) {
    if (cs!=null) {
      int nc = cs[0].length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = cs[0][ic].length;
        float sum = 0.0f;
        for (int ik=0; ik<nk; ++ik) {
          sum += cs[0][ic][ik];
        }
        float avg = sum/(float)nk;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = (int)cs[0][ic][ik];
          int i2 = (int)cs[1][ic][ik];
          int i3 = (int)cs[2][ic][ik];
          r[i3][i2][i1] = avg-(float)i1;
        }

      }
    }
  }


  /*
  public static void initializeShifts(float[][][] cs, float[][][] r) {
    if (cs!=null) {
      int nc = cs[0].length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = cs[0][ic].length;
        int ik = 0;
        int i1 = (int)cs[0][ic][ik];
        int i2 = (int)cs[1][ic][ik];
        int i3 = (int)cs[2][ic][ik];
        float i1f = (float)i1 + cs[3][ic][ik]; 
        r[i3][i2][i1] = i1-i1f;
        for (ik=1; ik<nk; ++ik) {
          float ip = i1f;
          float rp = r[i3][i2][i1];
          i1  = (int)cs[0][ic][ik];
          i2  = (int)cs[1][ic][ik];
          i3  = (int)cs[2][ic][ik];
          i1f = (float)i1+cs[3][ic][ik];
          r[i3][i2][i1] = rp+ip-i1f;
        }
      }
    }
  }
  */

  public static void checkShifts(float[][] k1, float[][] k2, float[][] k3, float[][][] r) {
    if (k1!=null && k2!=null &&k3!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        trace("ic="+ic);
        int nk = k1[ic].length;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = (int)k1[ic][ik];
          int i2 = (int)k2[ic][ik];
          int i3 = (int)k3[ic][ik];
          trace("  i1="+i1+" i2="+i2+" i3="+i3+" r="+r[i3][i2][i1]+" u="+(i1+r[i3][i2][i1]));
          //assert r[i2][i1]==rp+ip-i1:"shifts r satisfy constraints";
        }
      }
    }
  }

  public static void constrain(float[][][] cs, float[][][] x) 
  {
    if (cs!=null) {
      int nc = cs[0].length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = cs[0][ic].length;
        float sum = 0.0f;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = (int)cs[0][ic][ik];
          int i2 = (int)cs[1][ic][ik];
          int i3 = (int)cs[2][ic][ik];
          sum += x[i3][i2][i1];
        }
        float avg = sum/(float)nk;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = (int)cs[0][ic][ik];
          int i2 = (int)cs[1][ic][ik];
          int i3 = (int)cs[2][ic][ik];
          x[i3][i2][i1] = avg;
        }
      }
    }
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
  private static void smooth1(final float sigma, final float[][][] s, final float[][][] x) {
    final int n3 = x.length;
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
          if(wpi<0.05f) {wpi=0.05f;}
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

  private static void applyLhs(
    final float[][][] wp, 
    final float[][][] p2, final float[][][] p3,
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,p2,p3,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,p2,p3,x,y);
    }});
  }
  private static void applyLhsSlice3(
    int i3, 
    float[][][] wp, float[][][] p2, float[][][] p3,
    float[][][] x, float[][][] y) 
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
        if(wpi<0.05f) {wpi=0.05f;}
        float w1i = 0.02f;
        //float w1i = w1[i3][i2][i1];
        float p2i = p2[i3][i2][i1];
        float p3i = p3[i3][i2][i1];
        float w1s = w1i*w1i;
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
