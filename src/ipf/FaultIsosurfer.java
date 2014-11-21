/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Compute a volume whoes isosurfaces are fault surfaces.
 * <p>
 * Assume we can obtain an image with fault normal vectors, then we can 
 * solve a Poisson problem to compute a volume whose isosurfaces are fault 
 * surfaces.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.07.26
 */
public class FaultIsosurfer {

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

  /**
   * Gets mappings computed from specified slopes and planarities.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param u1 array of 1st component of fault normal vectors.
   * @param u2 array of 2nd component of fault normal vectors.
   * @param u3 array of 3rd component of fault normal vectors.
   * @param ep array of fault likelihoods.
   */
  public float[][][] findBlocks(float[][][][] flpt) {
    // Sampling parameters.
    final int n3 = flpt[0].length;
    final int n2 = flpt[0][0].length;
    final int n1 = flpt[0][0][0].length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    //float[][][] wp = new float[n3][n2][n1];
    float[][][] wp = fillfloat(0.0f,n1,n2,n3);
    FaultCell[] fc = faultCellsAndNormals(flpt,u1,u2,u3,wp);


    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] f = new float[n3][n2][n1]; // fault isosurface volume, in samples
    //initialization(f);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vf = new VecArrayFloat3(f);
    Smoother3 smoother3 = new Smoother3(_sigma1,_sigma2,_sigma3,wp,u1,u2,u3);
    A3 a3 = new A3(smoother3,wp,fc);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(wp,u1,u2,u3,b);
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vf);
    smoother3.apply(f);
    return f;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1 = 0.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 6.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations

  private FaultCell[] faultCellsAndNormals(float[][][][] flpt, 
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] wp) {
    FaultSkinner fs = new FaultSkinner();
    //fs.setMinSkinSize(3000);
    fs.setGrowLikelihoods(0.2f,0.5f);
    FaultCell[] fc = fs.findCells(flpt);
    int nc = fc.length;
    for (int ic=0; ic<nc; ++ic) {
      int i1 = fc[ic].i1;
      /*
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      u1[i3][i2][i1] = fc[ic].w1;
      u2[i3][i2][i1] = fc[ic].w2;
      u3[i3][i2][i1] = fc[ic].w3;
      wp[i3][i2][i1] = fc[ic].fl;
      */
      int i2m = fc[ic].i2m;
      int i3m = fc[ic].i3m;
      int i2p = fc[ic].i2p;
      int i3p = fc[ic].i3p;
      u1[i3m][i2m][i1] = -fc[ic].w1;
      u2[i3m][i2m][i1] = -fc[ic].w2;
      u3[i3m][i2m][i1] = -fc[ic].w3;

      u1[i3p][i2p][i1] =  fc[ic].w1;
      u2[i3p][i2p][i1] =  fc[ic].w2;
      u3[i3p][i2p][i1] =  fc[ic].w3;

      wp[i3m][i2m][i1] =  fc[ic].fl;
      wp[i3p][i2p][i1] =  fc[ic].fl;
    }
    return fc;
  }

  private void initialization(float[][][] f) {
    float vi = 1.0f;
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        f[i3][0   ][i1] = vi;
        f[i3][1   ][i1] = vi;
        f[i3][n2-1][i1] = vi;
        f[i3][n2-2][i1] = vi;
      }
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f[0   ][i2][i1] = vi;
        f[1   ][i2][i1] = vi;
        f[n3-1][i2][i1] = vi;
        f[n3-2][i2][i1] = vi;
      }
    }
  }

  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, float[][][] wp, FaultCell[] fc){
      _s3 = s3;
      _wp = wp;
      _fc = fc;
      //testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      _s3.apply(z);
      zero(y);
      applyLhs(_wp,_fc,z,y);
      _s3.applyTranspose(y);
    }
    private Smoother3 _s3;
    private float[][][] _wp;
    private FaultCell[] _fc;
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

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother3 {
    public Smoother3(
      float sigma1, float sigma2, float sigma3, 
      float[][][] ep,float[][][] u1,float[][][] u2,float[][][] u3) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      _ep = ep;
      _u1 = u1;
      _u2 = u2;
      _u3 = u3;
      //testSpd();
    }
    public void apply(float[][][] x) {
      /*
      smooth3(_sigma3,_ep,x);
      smooth2(_sigma2,_ep,x);
      smooth1(_sigma1,x);
      */
      //zero1(x);
      smooth(_sigma1,_u1,_u2,_u3,x);
    }
    public void applyTranspose(float[][][] x) {
      //zero1(x);
      /*
      smooth1(_sigma1,x);
      smooth2(_sigma2,_ep,x);
      smooth3(_sigma3,_ep,x);
      */
      smooth(_sigma1,_u1,_u2,_u3,x);
    }
    private float _sigma1,_sigma2,_sigma3;
    private float[][][] _ep;
    private float[][][] _u1;
    private float[][][] _u2;
    private float[][][] _u3;

    private void zero1(float[][][] x) {
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          x[i3][i2][   0] = 0.0f;
          x[i3][i2][n1-1] = 0.0f;
        }
      }
    }
    public void testSpd() {
      // symmetric: y'Ax = x'(A'y) = x'Ay
      // positive-semidefinite: x'Ax >= 0
      int n1 = _ep[0][0].length;
      int n2 = _ep[0].length;
      int n3 = _ep.length;
      float[][][] x = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] y = sub(randfloat(n1,n2,n3),0.5f);
      float[][][] ax = copy(x);
      float[][][] ay = copy(y);
      VecArrayFloat3 vx = new VecArrayFloat3(x);
      VecArrayFloat3 vy = new VecArrayFloat3(y);
      VecArrayFloat3 vax = new VecArrayFloat3(ax);
      VecArrayFloat3 vay = new VecArrayFloat3(ay);
      apply(ax);
      apply(ay);
      applyTranspose(ax);
      applyTranspose(ay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("S3: yax="+yax+" xay="+xay);
      System.out.println("S3: xax="+xax+" yay="+yay);
    }
  }

  public static void smooth(double sigma, 
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] g) 
  {
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    EigenTensors3 d = new EigenTensors3(n1,n2,n3,false);
    d.setEigenvalues(0.001f,1.00f,1.00f);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float usi = sqrt(u1i*u1i+u2i*u2i);
          float w1i = 0.0f;
          float w2i = 0.0f;
          float w3i = 0.0f;
          if(usi!=0.0f) {
            w1i = -u2i/usi;
            w2i =  u1i/usi;
            w3i = 0.0f;
          }
          d.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          d.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
        }
      }
    }
    float c = (float)(0.5*sigma*sigma);
    float[][][] h = new float[n3][n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,c,g,h);
    copy(h,g);
  }


  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
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

  private static void makeRhs(
    float[][][] wp, float[][][] u1, float[][][] u2, float[][][] u3, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1,i3m=0; i3<n3; ++i3,++i3m) {
      for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
        for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
          float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
          float wps = wpi*wpi;
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float y1 = wps*u1i;
          float y2 = wps*u2i;
          float y3 = wps*u3i;
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
  private static void applyLhs( final float[][][] wp, final FaultCell[] fc, 
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,x,y);
    }});

    final int nc = fc.length;
    Parallel.loop(nc,new Parallel.LoopInt() { 
    public void compute(int ic) {
      applyScreen(ic,fc,x,y);
    }});
  }

  private static void applyScreen(
    int ic, FaultCell[] fc, float[][][] x, float[][][] y)
  {
    int i1 = fc[ic].i1;
    int i2 = fc[ic].i2;
    int i3 = fc[ic].i3;
    float alpha = 1.0f;
    y[i3][i2][i1] += alpha*x[i3][i2][i1];
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
        float x1 = 0.25f*(xa+xb+xc+xd);
        float x2 = 0.25f*(xa-xb+xc-xd);
        float x3 = 0.25f*(xa+xb-xc-xd);
        float y1 = wps*x1;
        float y2 = wps*x2;
        float y3 = wps*x3;
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
}
