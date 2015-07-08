/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package rgi;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Structure-oriented seismic impedance inversion.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.07.03
 */
public class ImpedanceInversion2 {

  public void setWavelet(int nh, float fpeak) {
    _nh = nh;
    _fpeak = fpeak;
  }

  public void setTensors(EigenTensors2 et) {
    _et = et;
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
   * Compute seismic impedance
   * @param f array of seismic amplitude image
   */
  public float[][] applyForImpedance(float[][]f, float[][] el) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[] h = getWavelet();
    int kh = -(h.length-1)/2;
    float[][] b = new float[n2][n1]; // right-hand side
    float[][] r = new float[n2][n1];
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(n1,n2,_sigma1,_sigma2,el);
    A2 a2 = new A2(kh,h,el,_et,smoother2);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(kh,h,f,b);
    //smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    //smoother2.apply(r);
    return r;
  }

  public float[][] applyForSyn(float[][] p) {
    int n2 = p.length;
    float[] h = getWavelet();
    int kh = -(h.length-1)/2;
    for (int i2=0; i2<n2; ++i2) {
      applyConv(kh,h,p[i2]);
    }
    return p;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nh = 181;
  private float _fpeak = 35.0f;
  private EigenTensors2 _et = null;
  private float _sigma1 = 24.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 24.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  

  private static class RickerWavelet {
    public RickerWavelet(double fpeak) {
      _a = PI*fpeak;
    }
    public float getValue(double t) {
      double b = _a*t;
      double c = b*b;
      return (float)((1.0-2.0*c)*exp(-c));
    }
    public float getWidth() {
      return (float)(6.0/_a);
    }
    private double _a;
  }

  private float[] getWavelet(){
    float dt = 0.002f;
    int kh = (_nh-1)/2;
    float[] w = new float[_nh];
    RickerWavelet rw = new RickerWavelet(_fpeak);
    for (int ik=-kh; ik<=kh; ++ik) {
      w[ik+kh] = rw.getValue(dt*ik);
    }
    return w;
  }

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(int kh, float[] h, float[][] wp, EigenTensors2 et,Smoother2 s2) {
      _kh = kh;
      _h  = h;
      _s2 = s2;
      _wp = wp;
      _et = et;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      zero(y);
      //_s2.apply(z);
      applyLhs1(_kh,_h,z,y);
      //_s2.applyTranspose(y);
    }
    private int _kh;
    private float[] _h;
    private Smoother2 _s2;
    private float[][] _wp;
    private EigenTensors2 _et;
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
    private float[][] _el;
    private float _sigma1,_sigma2;
  }
  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][] x) {
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
    lsf.setPreconditioner(true);
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

  private static void makeRhs(int kh, float[] h, float[][] f, float[][] y) {
    int n2 = f.length;
    copy(f,y);
    for (int i2=0; i2<n2; ++i2) {
      applyConvT(kh,h,y[i2]);
      //applyConvT(h,y[i2]);
    }
  }

  private static void applyLhs1(
    int kh, float[] h, float[][] x, float[][] y) 
  {
    copy(x,y);
    int n2 = y.length;
    for (int i2=0; i2<n2; ++i2) {
      applyConv(kh,h,y[i2]);
      applyConvT(kh,h,y[i2]);
      //applyConv(h,y[i2]);
      //applyConvT(h,y[i2]);

    }
  }


  private static void applyConv(float[] h, float[] x) {
    int n1 = x.length;
    int nh = h.length;
    float[] t = copy(x);
    float[] z = new float[nh+n1];
    float[] g = new float[nh+n1];
    copy(n1,0,t,nh,z);
    for (int i1=0; i1<n1+nh; ++i1) {
      int ik = i1;
      int i1k = 0;
      if(i1>=nh) {ik=nh-1;i1k=i1-nh+1;}
      float sum=0.0f;
      for (int k1=ik; k1>=0; --k1,++i1k) {
        sum += h[k1]*z[i1k];
      }
      g[i1] = sum;
    }
    copy(n1,0,g,0,x);
  }

  private static void applyConvT(float[] h, float[] x) {
    int n1 = x.length;
    int nh = h.length;
    float[] t = copy(x);
    float[] z = new float[nh+n1];
    float[] g = new float[nh+n1];
    copy(n1,0,t,0,z);
    for (int i1=0; i1<n1+nh; ++i1) {
      int ik = n1-i1;
      if(ik>=nh){ik=nh;}
      float sum = 0.0f;
      for (int k=0, i1k=i1; k<ik; ++k, ++i1k) {
        sum += h[k]*z[i1k];
      }
      g[i1] = sum;
    }
    copy(n1,0,g,0,x);
  }

  private static void applyConv(int kh, float[] h, float[] x) {
    int nt = x.length;
    int nh = h.length;
    float[] g = new float[nt];
    conv(nh,kh,h,nt,0,x,nt,0,g);
    copy(g,x);
  }

  private static void applyConvT(int kh, float[] h, float[] x) {
    reverse(x);
    int nt = x.length;
    int nh = h.length;
    float[] g = new float[nt];
    conv(nh,kh,h,nt,0,x,nt,0,g);
    copy(g,x);
    reverse(x);
  }



  private static void reverse(float[] z) {
    int n1 = z.length;
    for (int i1=0,j1=n1-1; i1<j1; ++i1,--j1) {
      float zt = z[i1];
      z[i1] = z[j1];
      z[j1] = zt;
    }
  }


  private static void applyLhs2(
    float[][] wp, EigenTensors2 d, float[][] x, float[][] y)
  {
    int n2 = x.length;
    int n1 = x[0].length;
    float[] ds = fillfloat(1.0f,3);
    ds[0] = 1.0f;
    ds[1] = 0.0f;
    ds[2] = 1.0f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if(d!=null){d.getTensor(i1,i2,ds);}
        float wpi = (wp!=null)?wp[i2][i1]:1.0f;
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
