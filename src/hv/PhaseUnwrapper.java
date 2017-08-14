/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hv;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Weighted phase unwrapping. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2017.08.09
 */
public class PhaseUnwrapper{

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

  public void phaseGradient(
    float[][][] ph, float[][][] u1, float[][][] u2, float[][][] u3) {
    gradients(ph,u1,u2,u3);
  }

  public float[][][][] phaseSlope(
    float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3){
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0; i1<n1; ++i1){
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float usi = sqrt(u1i*u1i+u2i*u2i+u3i*u3i);
      if(u1i==0f) continue;
      u1i /= usi;
      u2i /= usi;
      u3i /= usi;
      float p2i = -u2i/u1i;
      float p3i = -u3i/u1i;
      p2i = max(p2i,-2);
      p3i = max(p3i,-2);
      p2i = min(p2i, 2);
      p3i = min(p3i, 2);

      p2[i3][i2][i1] = p2i;
      p3[i3][i2][i1] = p3i;
    }}}
    return new float[][][][]{p2,p3};

  }


  public void checkNans(
    float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    for (int i3=0; i3<n3; ++i3){
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0; i1<n1; ++i1){
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      if (Float.isNaN(u1i)||Float.isNaN(u2i)||Float.isNaN(u3i)) {
        u1[i3][i2][i1] = 0;
        u2[i3][i2][i1] = 0;
        u3[i3][i2][i1] = 0;
      }
    }}}

  }

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


  public float[][][] horizonVolumeFromRgt(float[][][] u1) {
    cleanRGT(u1);
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    Sampling s1 = new Sampling(n1);
    final float[][][] x1 = new float[n3][n2][n1];
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(u1[i3][i2],x1[i3][i2]);
    }});
    return x1;

  }

  /**
   * @param wp array of weights.
   * @param u1 array of 1st component of phase gradient.
   * @param u2 array of 2nd component of phase gradient.
   * @param u3 array of 3rd component of phase gradient.
   */
  public float[][][] phaseUnwrapping(
    float[][][] wp,
    float[][][] u1, float[][][] u2, float[][][] u3)
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] f = new float[n3][n2][n1]; // fault isosurface volume, in samples
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vf = new VecArrayFloat3(f);
    Smoother3 smoother3 = new Smoother3(_sigma1,_sigma2,_sigma3);
    A3 a3 = new A3(_weight1, wp, smoother3);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(wp,u1,u2,u3,b);
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vf);
    smoother3.apply(f);
    return f;
  }

    // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(float w1, float[][][] wp, Smoother3 s3){
      _w1 = w1;
      _wp = wp;
      _s3 = s3;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      _s3.apply(z);
      zero(y);
      applyLhs(_w1,_wp,z,y);       //laplacian operator
      _s3.applyTranspose(y);
    }
    private float[][][] _wp;
    private float _w1;
    private Smoother3 _s3;
  }

    // right-hand side for 3D
  private static void makeRhs(
    float[][][] wp,
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1,i3m=0; i3<n3; ++i3,++i3m) {
    for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
    for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
      //float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
      float wpi = 1f;
      float wps = wpi*wpi;
      float y1 = u1[i3][i2][i1];
      float y2 = u2[i3][i2][i1];
      float y3 = u3[i3][i2][i1];
      y1 *= wps;
      y2 *= wps;
      y3 *= wps;
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


  // left-hand side for 3D
  private static void applyLhs(
    final float w1, final float[][][] wp,
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,w1,wp,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,w1,wp,x,y);
    }});
  }

  private static void applyLhsSlice3(
    int i3, float w1, float[][][] wp, float[][][] x, float[][][] y) 
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
        //float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
        float wpi = 1f;
        float wps = wpi*wpi;
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
        y1 *= wps;
        y2 *= wps;
        y3 *= wps;
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

  // Smoother used as a preconditioner.
  private static class Smoother3 {
    public Smoother3(
      float sigma1, float sigma2, float sigma3) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
    }
    public void apply(float[][][] x) {
      smooth3(_sigma3,x);
      smooth2(_sigma2,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,x);
      smooth3(_sigma3,x);
    }
    private float _sigma1,_sigma2,_sigma3;
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

  public void toRadius(float[][][] ph) {
    int n1 = ph[0][0].length;
    int n2 = ph[0].length;
    int n3 = ph.length;
    float pi  = (float)(Math.PI/180.0);
    for (int i3=0; i3<n3; i3++){
    for (int i2=0; i2<n2; i2++){
    for (int i1=0; i1<n1; i1++){
      ph[i3][i2][i1] *= pi;
    }}}
  }


  private void gradients(
    float[][][] ph,  
    float[][][] ph1, float[][][] ph2, float[][][] ph3){
    float pi  = (float)Math.PI;
    float pi2 = 2.f*pi;
    int n1 = ph[0][0].length;
    int n2 = ph[0].length;
    int n3 = ph.length;
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        for (int i1=1; i1<n1-1; i1++){
          float phi1 = ph[i3][i2][i1-1];
          float phi2 = ph[i3][i2][i1  ];
          float d1 = phi2 - phi1;
          if (d1> pi){d1 = -pi2+d1;}
          if (d1<-pi){d1 =  pi2+d1;}
          ph1[i3][i2][i1] = d1;//(d1+d2)/2.f;
        }
      }
    }
   
    for (int i3=0; i3<n3; i3++){
      for (int i2=1; i2<n2-1; i2++){
        for (int i1=0; i1<n1; i1++){
          float phi1 = ph[i3][i2-1][i1];
          float phi2 = ph[i3][i2  ][i1];
          float d1 = phi2 - phi1;
          if (d1> pi){d1 = -pi2+d1;}
          if (d1<-pi){d1 =  pi2+d1;}
          ph2[i3][i2][i1] = d1;//(d1+d2)/2.f;
        }
      }
    }
    
    for (int i3=1; i3<n3-1; i3++){
      for (int i2=0; i2<n2; i2++){
        for (int i1=0; i1<n1; i1++){
          float phi1 = ph[i3-1][i2][i1];
          float phi2 = ph[i3  ][i2][i1];
          float d1 = phi2-phi1;
          if (d1> pi){d1 = -pi2+d1;}
          if (d1<-pi){d1 =  pi2+d1;}
          ph3[i3][i2][i1] = d1;//(d1+d2)/2.f;
        }
      }
    }
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        ph1[i3][i2][n1-1] = ph1[i3][i2][n1-2];
        ph1[i3][i2][0   ] = ph1[i3][i2][1   ];
      }
    }
    for (int i3=0; i3<n3; i3++){
      for (int i1=0; i1<n1; i1++){
        ph2[i3][n2-1][i1] = ph2[i3][n2-2][i1];
        ph2[i3][0   ][i1] = ph2[i3][1   ][i1];
      }
    }
    for (int i1=0; i1<n1; i1++){
      for (int i2=0; i2<n2; i2++){
        ph3[n3-1][i2][i1] = ph3[n3-2][i2][i1];
        ph3[0   ][i2][i1] = ph3[1   ][i2][i1];
      }
    }
    for (int i1=0; i1<n1;i1++){
      for (int i2=0; i2<n2;i2++){
        for (int i3=0; i3<n3;i3++){
          float ph1i = ph1[i3][i2][i1];
          float ph2i = ph2[i3][i2][i1];
          float ph3i = ph3[i3][i2][i1];
          float phs = sqrt(ph1i*ph1i + ph2i*ph2i + ph3i*ph3i);
          if(phs==0f) {
            ph1i = 1f;
            ph2i = 0f;
            ph3i = 0f;
          }
          ph1i = ph1i/phs;
          ph2i = ph2i/phs;
          ph3i = ph3i/phs;
          if (ph1i<0) {
            ph1i = -ph1i;
            ph2i = -ph2i;
            ph3i = -ph3i;
          }
          ph1[i3][i2][i1]  = ph1i;
          ph2[i3][i2][i1]  = ph2i;
          ph3[i3][i2][i1]  = ph3i;
        }
      }
    }
  }
  private void cleanRGT(float[][][] u) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    float dMin = 0.001f;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float ui = u[i3][i2][i1  ];
          float um = u[i3][i2][i1-1];
          float ut = um+dMin;
          if (ui<ut) {u[i3][i2][i1]=ut;}
        }
      }
    }
  }



  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 6.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private float _weight1 = 0.01f;
  private int _niter = 100; // maximum number of CG iterations

}
