/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import vec.*;
import util.*;

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

  public float[][][][] normalsFromCellsOpen(int n1, int n2, int n3, FaultCell[] fc) {
    int nc = fc.length;
    float[] wf = new float[nc];
    float[][] xf = new float[3][nc];
    float[][] uf = new float[3][nc];
    float[][][][] us = new float[5][n3][n2][n1];
    for (int ic=0; ic<nc; ++ic) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      float u1 = fc[ic].w1;
      float u2 = fc[ic].w2;
      float u3 = fc[ic].w3;
      float wi = fc[ic].fl;
      //float wi = 1.0f;
      us[0][i3][i2][i1] = wi;
      us[1][i3][i2][i1] = u1;
      us[2][i3][i2][i1] = u2;
      us[3][i3][i2][i1] = u3;
      us[4][i3][i2][i1] = 1.0f;
      wf[ic] = wi;
      xf[0][ic] = i1;
      xf[1][ic] = i2;
      xf[2][ic] = i3;
      uf[0][ic] = u1;
      uf[1][ic] = u2;
      uf[2][ic] = u3;
    }
    //nearestInterp(xf,uf,wf,us);
    return us;
  }

  public float[][][][] normalsFromCellsClose(int n1, int n2, int n3, FaultCell[] fc) {
    int nc = fc.length;
    float[] wf = new float[nc*3];
    float[] sg = new float[nc*3];
    float[][] xf = new float[3][nc*3];
    float[][] uf = new float[3][nc*3];
    float[][][][] us = new float[5][n3][n2][n1];
    int ik = 0;
    for (int ic=0; ic<nc; ++ic) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      //float wi = fc[ic].fl;
      float wi = 1.0f;
      int i2m = fc[ic].i2m;
      int i3m = fc[ic].i3m;
      int i2p = fc[ic].i2p;
      int i3p = fc[ic].i3p;
      float u1 = fc[ic].w1;
      float u2 = fc[ic].w2;
      float u3 = fc[ic].w3;
      wf[ik] = wi;
      xf[0][ik] = i1;
      xf[1][ik] = i2;
      xf[2][ik] = i3;
      uf[0][ik] = u1;
      uf[1][ik] = u2;
      uf[2][ik] = u3;
      ik ++;
      wf[ik] = wi;
      sg[ik] = 1.0f;
      xf[0][ik] = i1;
      xf[1][ik] = i2m;
      xf[2][ik] = i3m;
      uf[0][ik] = -u1;
      uf[1][ik] = -u2;
      uf[2][ik] = -u3;
      ik ++;
      wf[ik] = wi;
      sg[ik] = 1.0f;
      xf[0][ik] = i1;
      xf[1][ik] = i2p;
      xf[2][ik] = i3p;
      uf[0][ik] = u1;
      uf[1][ik] = u2;
      uf[2][ik] = u3;
      ik ++;
    }
    nearestInterp(xf,uf,wf,sg,us);
    for (int ic=0; ic<nc; ++ic) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      //float wi = fc[ic].fl;
      float wi = 1.0f;
      int i2m = fc[ic].i2m;
      int i3m = fc[ic].i3m;
      int i2p = fc[ic].i2p;
      int i3p = fc[ic].i3p;
      float u1 = fc[ic].w1;
      float u2 = fc[ic].w2;
      float u3 = fc[ic].w3;
      us[0][i3][i2][i1] = wi;
      us[1][i3][i2][i1] = u1;
      us[2][i3][i2][i1] = u2;
      us[3][i3][i2][i1] = u3;

      us[0][i3m][i2m][i1] = wi;
      us[1][i3m][i2m][i1] = -u1;
      us[2][i3m][i2m][i1] = -u2;
      us[3][i3m][i2m][i1] = -u3;
      us[4][i3m][i2m][i1] = 1.0f;

      us[0][i3p][i2p][i1] = wi;
      us[1][i3p][i2p][i1] = u1;
      us[2][i3p][i2p][i1] = u2;
      us[3][i3p][i2p][i1] = u3;
      us[4][i3p][i2p][i1] = 1.0f;
    }
    return us;
  }


  public static void nearestInterp(
    final float[][]xf, final float[][] uf, 
    final float[] wf, final float[] sg, final float[][][][] us)
  {
    final int n3 = us[0].length;
    final int n2 = us[0][0].length;
    final int n1 = us[0][0][0].length;
    final int[] d = new int[]{16,4,4};
    final KdTree kt = new KdTree(xf);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] xmin = new float[3];
          float[] xmax = new float[3];
          int[] ix = new int[]{i1,i2,i3};
          defineRange(d,ix,n1,n2,n3,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd==0){continue;}
          int ic = id[0];
          float dd = 800.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float x1i = xf[0][ip];
            float x2i = xf[1][ip];
            float x3i = xf[2][ip];
            float u1i = uf[0][ip];
            float u2i = uf[1][ip];
            float u3i = uf[2][ip];
            float d1i = i1-x1i;
            float d2i = i2-x2i;
            float d3i = i3-x3i;
            float dsi = d1i*u1i+d2i*u2i+d3i*u3i;
            if(dsi<0.0f && sg[ic]==1.0f){continue;}
            if(abs(dsi)<dd) {ic = ip;dd = abs(dsi);}
          }
          if(dd==800.0f){continue;}
          us[0][i3][i2][i1] = wf[ic];
          float u1i = uf[0][ic];
          float u2i = uf[1][ic];
          float u3i = uf[2][ic];
          us[1][i3][i2][i1] = u1i;
          us[2][i3][i2][i1] = u2i;
          us[3][i3][i2][i1] = u3i;
          if(dd<0.1f&&sg[ic]==1.0f){us[4][i3][i2][i1]=1.0f;}
        }
      }
    }});
  }
  
  public static void nearestInterp(
    final float[][]xf, final float[][] uf, 
    final float[] wf, final float[][][][] us)
  {
    final int n3 = us[0].length;
    final int n2 = us[0][0].length;
    final int n1 = us[0][0][0].length;
    final int[] d = new int[]{15,6,6};
    final KdTree kt = new KdTree(xf);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float[] xmin = new float[3];
          float[] xmax = new float[3];
          int[] ix = new int[]{i1,i2,i3};
          defineRange(d,ix,n1,n2,n3,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd==0){continue;}
          int ic = id[0];
          float dd = 800.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float x1i = xf[0][ip];
            float x2i = xf[1][ip];
            float x3i = xf[2][ip];
            float u1i = uf[0][ip];
            float u2i = uf[1][ip];
            float u3i = uf[2][ip];
            float d1i = i1-x1i;
            float d2i = i2-x2i;
            float d3i = i3-x3i;
            float dsi = abs(d1i*u1i+d2i*u2i+d3i*u3i);
            if(dsi<dd) {ic = ip;dd = dsi;}
          }
          if(dd==800.0f){continue;}
          float u1i = uf[0][ic];
          float u2i = uf[1][ic];
          float u3i = uf[2][ic];
          us[1][i3][i2][i1] = u1i;
          us[2][i3][i2][i1] = u2i;
          us[3][i3][i2][i1] = u3i;
          us[0][i3][i2][i1] = wf[ic];
          //if(dd<0.2f){us[4][i3][i2][i1]=1.0f;}
        }
      }
    }});
  }

  private static void defineRange(int[] d, int[] i, 
    int n1, int n2, int n3, float[] xmin, float[] xmax) 
  {
    int i1 = i[0];
    int i2 = i[1];
    int i3 = i[2];
    int d1 = d[0];
    int d2 = d[1];
    int d3 = d[2];
    int i1m = i1-d1; if(i1m<0){i1m=0;}
    int i2m = i2-d2; if(i2m<0){i2m=0;}
    int i3m = i3-d3; if(i3m<0){i3m=0;}
    int i1p = i1+d1; if(i1p>=n1){i1p=n1-1;}
    int i2p = i2+d2; if(i2p>=n2){i2p=n2-1;}
    int i3p = i3+d3; if(i3p>=n3){i3p=n3-1;}
    xmin[0] = i1m; xmin[1] = i2m; xmin[2] = i3m;
    xmax[0] = i1p; xmax[1] = i2p; xmax[2] = i3p;
  }


  /**
   * @param us[0] array of fault likelihoods.
   * @param us[1] array of 1st component of fault normal vectors.
   * @param us[2] array of 2nd component of fault normal vectors.
   * @param us[3] array of 3rd component of fault normal vectors.
   */
  public float[][][] faultIndicator(float[][][][] us){
    int n3 = us[0].length;
    int n2 = us[0][0].length;
    int n1 = us[0][0][0].length;
    // Compute shifts r(x1,x2,x3), in samples.
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] f = new float[n3][n2][n1]; // fault isosurface volume, in samples
    //initialization(f);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vf = new VecArrayFloat3(f);
    Smoother3 smoother3 = new Smoother3(n1,n2,n3,_sigma1,_sigma2,_sigma3,us[0]);
    A3 a3 = new A3(smoother3,us[0],us[4]);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(us[0],us[1],us[2],us[3],b);
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vf);
    smoother3.apply(f);
    mark(us[0],f);
    return f;
  }

  private void mark(float[][][] w, float[][][] f) {
    float v = -30.0f;
    int n3 = w.length;
    int n2 = w[0].length;
    int n1 = w[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float wi = w[i3][i2][i1];
          if(wi==0.0f) {
            f[i3][i2][i1] = v;
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 6.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 100; // maximum number of CG iterations

  private void initialization(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          f[i3][i2][i1] = (float)i1;
        }
      }
    }
  }

  // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, float[][][] wp, float[][][] mk){
      _s3 = s3;
      _wp = wp;
      _mk = mk;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      _s3.apply(z);
      zero(y);
      applyLhs(_wp,_mk,z,y);
      _s3.applyTranspose(y);
    }
    private Smoother3 _s3;
    private float[][][] _wp;
    private float[][][] _mk;
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother3 {
    public Smoother3(
      int n1, int n2, int n3, 
      float sigma1, float sigma2, float sigma3, float[][][] ep) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      _ep = ep;
      //testSpd();
    }
    public void apply(float[][][] x) {
      smooth3(_sigma3,_ep,x);
      smooth2(_sigma2,_ep,x);
      smooth1(_sigma1,x);
      //zero1(x);
    }
    public void applyTranspose(float[][][] x) {
      //zero1(x);
      smooth1(_sigma1,x);
      smooth2(_sigma2,_ep,x);
      smooth3(_sigma3,_ep,x);
    }
    private float _sigma1,_sigma2,_sigma3;
    private float[][][] _ep;
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

  // Smoothing for dimension 1.
  private static void smooth2(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 1.
  private static void smooth3(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply3(x,x);
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
  private static void applyLhs(
    final float[][][] wp, final float[][][] mk, final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,mk,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,wp,mk,x,y);
    }});
  }
  private static void applyLhsSlice3(
    int i3, float[][][] wp, float[][][] mk, float[][][] x, float[][][] y) 
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
        if(mk[i3][i2][i1]==1.0f){y00[i1] += wps*x00[i1];}
      }
    }
  }
}
