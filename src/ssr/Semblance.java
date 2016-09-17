package ssr;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Structure- and stratigraphy-oriented semblance. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.15
 */

public class Semblance {

  public float[][][] smoothVW(float sigma, EigenTensors3 ets, float[][][] fx) {
    ets.setEigenvalues(0.0001f,1.0f,1.0f);
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = new float[n3][n2][n1];
    _lsf.apply(ets,sigma,fx,gx);
    return gx;
  }

  public float[][][][] applyForSemblanceNumDen(
    float[][][] p2, float[][][] p3, float[][][] f) 
  { 
    return semblanceNumDen(p2,p3,f);
  }


  public float[][][] shapeSemblance(
    EigenTensors3 et, float[][][] wp, float[][][] sn, float[][][] sd) 
  {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][] b = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    Smoother3 smoother3 = new Smoother3(10,et);
    A3 a3 = new A3(smoother3,et,wp,sd);
    CgSolver cs = new CgSolver(0.01,100);
    makeRhs(wp,sn,b);
    cs.solve(a3,vb,vr);
    return r;
  }


    // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, EigenTensors3 et, float[][][] wp, float[][][] sd) 
    {
      _s3 = s3;
      _et = et;
      _wp = wp;
      _sd = sd;
      /*
      float n3 = wp.length;
      float n2 = wp[0].length;
      float n1 = wp[0][0].length;
      _sc = sum(wp)/(n1*n2*n3);
      */
      _sc = 0.5f;
      testSpd();
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      v3y.zero();
      applyLhs(_sc,_et,_wp,_sd,z,y);
    }
    private float _sc;
    private Smoother3 _s3;
    private float[][][] _wp, _sd;
    private EigenTensors3 _et;
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
    public Smoother3(float sigma, EigenTensors3 ets) {
      _sigma = sigma;
      _ets = ets;
    }
    public void apply(float[][][] x) {
      //_lsf.apply(_ets,_sigma,x,x);
      smooth3(_sigma,x);
      smooth2(_sigma,x);
      smooth1(_sigma,x);
    }
    public void applyTranspose(float[][][] x) {
      //_lsf.apply(_ets,_sigma,x,x);
      smooth1(_sigma,x);
      smooth2(_sigma,x);
      smooth3(_sigma,x);

    }
    private float _sigma;
    private EigenTensors3 _ets;
  }

  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth3(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply3(x,x);
  }



  private static void applyLhs(
    float c, EigenTensors3 d, float[][][] wp, float[][][] sd, 
    float[][][] x, float[][][] y) {
    applyLhs1(wp,sd,x,y);
    applyLhs2(d,c,x,y);
  }

  private static void applyLhs1(
    float[][][] wp, float[][][] sd, float[][][] x, float[][][] y) 
  {
    zero(y);
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float sdi = sd[i3][i2][i1];
      float wpi = wp[i3][i2][i1];
      float  xi  = x[i3][i2][i1];
      y[i3][i2][i1] += sdi*wpi*xi;
    }}}
  }

  private static void makeRhs(
    float[][][] wp, float[][][] sn, float[][][] b) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float sni = sn[i3][i2][i1];
      float wpi = wp[i3][i2][i1];
      b[i3][i2][i1] = wpi*sni;
    }}}
  }

  private static void applyLhs2(
    final Tensors3 d, final float c, 
    final float[][][] x, final float[][][] y) 
  {
    int i3start = 0; // for apply71
    final int i3step  = 7; //for apply71
    //int i3start = 1; // for apply22
    //final int i3step  = 2; //for apply22
    final int i3stop = x.length;
    for (int i3pass=0; i3pass<i3step; ++i3pass,++i3start) {
      Parallel.loop(i3start,i3stop,i3step,new Parallel.LoopInt() {
        public void compute(int i3) {
          apply71(i3,d,c,x,y);
          //apply22(i3,d,c,x,y);
        }
      });
    }
  }

  private static void apply22(
    int i3, Tensors3 d, float c, float[][][] x, float[][][] y) 
  {
    c *= 0.0625f;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float[] di = new float[6];
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x0m = x[i3  ][i2-1];
      float[] xm0 = x[i3-1][i2  ];
      float[] xmm = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y0m = y[i3  ][i2-1];
      float[] ym0 = y[i3-1][i2  ];
      float[] ymm = y[i3-1][i2-1];
      for (int i1=1,m1=0; i1<n1; ++i1,++m1) {
        d.getTensor(i1,i2,i3,di);
        float csi = c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d13 = di[2]*csi;
        float d22 = di[3]*csi;
        float d23 = di[4]*csi;
        float d33 = di[5]*csi;
        float xa = x00[i1]-xmm[m1];
        float xb = x00[m1]-xmm[i1];
        float xc = x0m[i1]-xm0[m1];
        float xd = xm0[i1]-x0m[m1];
        float x1 = xa-xb+xc+xd;
        float x2 = xa+xb-xc+xd;
        float x3 = xa+xb+xc-xd;
        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2+d23*x3;
        float y3 = d13*x1+d23*x2+d33*x3;
        float ya = y1+y2+y3; y00[i1] += ya; ymm[m1] -= ya;
        float yb = y1-y2+y3; y0m[i1] += yb; ym0[m1] -= yb;
        float yc = y1+y2-y3; ym0[i1] += yc; y0m[m1] -= yc;
        float yd = y1-y2-y3; ymm[i1] += yd; y00[m1] -= yd;
      }
    }
  }

  private static void apply71(
    int i3, Tensors3 d, float c, float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[] di = new float[6];
    int i3m3 = i3-3; if (i3m3<0) i3m3 = 0;
    int i3m2 = i3-2; if (i3m2<0) i3m2 = 0;
    int i3m1 = i3-1; if (i3m1<0) i3m1 = 0;
    int i3p0 = i3;
    int i3p1 = i3+1; if (i3p1>=n3) i3p1 = n3-1;
    int i3p2 = i3+2; if (i3p2>=n3) i3p2 = n3-1;
    int i3p3 = i3+3; if (i3p3>=n3) i3p3 = n3-1;
    int i2m3,i2m2=0,i2m1=0,i2p0=0,i2p1=0,i2p2=1,i2p3=2;
    final float c1 =  C71[1], c2 = C71[2], c3 = C71[3];
    for (int i2=0; i2<n2; ++i2) {
      i2m3 = i2m2; i2m2 = i2m1; i2m1 = i2p0;
      i2p0 = i2p1; i2p1 = i2p2; i2p2 = i2p3; ++i2p3;
      if (i2p1>=n2) i2p1 = n2-1;
      if (i2p2>=n2) i2p2 = n2-1;
      if (i2p3>=n2) i2p3 = n2-1;
      float[] xp0p0 = x[i3p0][i2p0], yp0p0 = y[i3p0][i2p0];
      float[] xp0m3 = x[i3p0][i2m3], yp0m3 = y[i3p0][i2m3];
      float[] xp0m2 = x[i3p0][i2m2], yp0m2 = y[i3p0][i2m2];
      float[] xp0m1 = x[i3p0][i2m1], yp0m1 = y[i3p0][i2m1];
      float[] xp0p1 = x[i3p0][i2p1], yp0p1 = y[i3p0][i2p1];
      float[] xp0p2 = x[i3p0][i2p2], yp0p2 = y[i3p0][i2p2];
      float[] xp0p3 = x[i3p0][i2p3], yp0p3 = y[i3p0][i2p3];
      float[] xm3p0 = x[i3m3][i2p0], ym3p0 = y[i3m3][i2p0];
      float[] xm2p0 = x[i3m2][i2p0], ym2p0 = y[i3m2][i2p0];
      float[] xm1p0 = x[i3m1][i2p0], ym1p0 = y[i3m1][i2p0];
      float[] xp1p0 = x[i3p1][i2p0], yp1p0 = y[i3p1][i2p0];
      float[] xp2p0 = x[i3p2][i2p0], yp2p0 = y[i3p2][i2p0];
      float[] xp3p0 = x[i3p3][i2p0], yp3p0 = y[i3p3][i2p0];
      int m3,m2=0,m1=0,p0=0,p1=0,p2=1,p3=2;
      for (int i1=0; i1<n1; ++i1) {
        m3 = m2; m2 = m1; m1 = p0;
        p0 = p1; p1 = p2; p2 = p3; ++p3;
        if (p1>=n1) p1 = n1-1;
        if (p2>=n1) p2 = n1-1;
        if (p3>=n1) p3 = n1-1;
        d.getTensor(i1,i2,i3,di);
        float csi = c;
        float d11 = di[0]*csi;
        float d12 = di[1]*csi;
        float d13 = di[2]*csi;
        float d22 = di[3]*csi;
        float d23 = di[4]*csi;
        float d33 = di[5]*csi;
        float x1  = c1*(xp0p0[p1]-xp0p0[m1]) +
                    c2*(xp0p0[p2]-xp0p0[m2]) +
                    c3*(xp0p0[p3]-xp0p0[m3]);
        float x2  = c1*(xp0p1[p0]-xp0m1[p0]) +
                    c2*(xp0p2[p0]-xp0m2[p0]) +
                    c3*(xp0p3[p0]-xp0m3[p0]);
        float x3  = c1*(xp1p0[p0]-xm1p0[p0]) +
                    c2*(xp2p0[p0]-xm2p0[p0]) +
                    c3*(xp3p0[p0]-xm3p0[p0]);
        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2+d23*x3;
        float y3 = d13*x1+d23*x2+d33*x3;
        float c1y1 = c1*y1; yp0p0[p1] += c1y1; yp0p0[m1] -= c1y1;
        float c2y1 = c2*y1; yp0p0[p2] += c2y1; yp0p0[m2] -= c2y1;
        float c3y1 = c3*y1; yp0p0[p3] += c3y1; yp0p0[m3] -= c3y1;
        float c1y2 = c1*y2; yp0p1[p0] += c1y2; yp0m1[p0] -= c1y2;
        float c2y2 = c2*y2; yp0p2[p0] += c2y2; yp0m2[p0] -= c2y2;
        float c3y2 = c3*y2; yp0p3[p0] += c3y2; yp0m3[p0] -= c3y2;
        float c1y3 = c1*y3; yp1p0[p0] += c1y3; ym1p0[p0] -= c1y3;
        float c2y3 = c2*y3; yp2p0[p0] += c2y3; ym2p0[p0] -= c2y3;
        float c3y3 = c3*y3; yp3p0[p0] += c3y3; ym3p0[p0] -= c3y3;
      }
    }
  }

  private static final float[] C71 = {
    0.0f, 0.830893f, -0.227266f, 0.042877f
  };

  // Computes fault semblance numerators and denominators.
  private float[][][][] semblanceNumDen(
    final float[][][] p2, final float[][][] p3, final float[][][] f) 
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] sn = new float[n3][n2][n1];
    final float[][][] sd = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmm = new float[n1];
      float[] xm0 = new float[n1];
      float[] xmp = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xpm = new float[n1];
      float[] xp0 = new float[n1];
      float[] xpp = new float[n1];
      float[] gmm = new float[n1];
      float[] gm0 = new float[n1];
      float[] gmp = new float[n1];
      float[] g0m = new float[n1];
      float[] g0p = new float[n1];
      float[] gpm = new float[n1];
      float[] gp0 = new float[n1];
      float[] gpp = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p2mm = p2[i3m][i2m];
        float[] p2mp = p2[i3m][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p20p = p2[i3 ][i2p];
        float[] p2pm = p2[i3p][i2m];
        float[] p2pp = p2[i3p][i2p];
        float[] p3mm = p3[i3m][i2m];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3mp = p3[i3m][i2p];
        float[] p3pm = p3[i3p][i2m];
        float[] p3p0 = p3[i3p][i2 ];
        float[] p3pp = p3[i3p][i2p];
        float[] sn32 = sn[i3][i2];
        float[] sd32 = sd[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          xmm[i1] = i1-p3mm[i1]-p2mm[i1];
          xm0[i1] = i1-p3m0[i1]         ;
          xmp[i1] = i1-p3mp[i1]+p2mp[i1];
          x0m[i1] = i1         -p20m[i1];
          x0p[i1] = i1         +p20p[i1];
          xpm[i1] = i1+p3pm[i1]-p2pm[i1];
          xp0[i1] = i1+p3p0[i1]         ;
          xpp[i1] = i1+p3pp[i1]+p2pp[i1];
        }
        si.interpolate(n1,1.0,0.0,fmm,n1,xmm,gmm);
        si.interpolate(n1,1.0,0.0,fm0,n1,xm0,gm0);
        si.interpolate(n1,1.0,0.0,fmp,n1,xmp,gmp);
        si.interpolate(n1,1.0,0.0,f0m,n1,x0m,g0m);
        si.interpolate(n1,1.0,0.0,f0p,n1,x0p,g0p);
        si.interpolate(n1,1.0,0.0,fpm,n1,xpm,gpm);
        si.interpolate(n1,1.0,0.0,fp0,n1,xp0,gp0);
        si.interpolate(n1,1.0,0.0,fpp,n1,xpp,gpp);
        float[] hmm = gmm, hm0 = gm0, hmp = gmp;
        float[] h0m = g0m, h00 = f00, h0p = g0p;
        float[] hpm = gpm, hp0 = gp0, hpp = gpp;
        if (            i2==0   ) h0m = h00;
        if (            i2==n2-1) h0p = h00;
        if (i3==0               ) hm0 = h00;
        if (i3==n3-1            ) hp0 = h00;
        if (i3==0    && i2==0   ) hmm = h00;
        if (i3==0    && i2==n2-1) hmp = h00;
        if (i3==n3-1 && i2==0   ) hpm = h00;
        if (i3==n3-1 && i2==n2-1) hpp = h00;
        for (int i1=0; i1<n1; ++i1) {
          float hmmi = hmm[i1];
          float hm0i = hm0[i1];
          float hmpi = hmp[i1];
          float h0mi = h0m[i1];
          float h00i = h00[i1];
          float h0pi = h0p[i1];
          float hpmi = hpm[i1];
          float hp0i = hp0[i1];
          float hppi = hpp[i1];
          float sumn = hmmi+hm0i+hmpi+
                       h0mi+h00i+h0pi+
                       hpmi+hp0i+hppi;
          float sumd = hmmi*hmmi+hm0i*hm0i+hmpi*hmpi+
                       h0mi*h0mi+h00i*h00i+h0pi*h0pi+
                       hpmi*hpmi+hp0i*hp0i+hppi*hppi;
          sn32[i1] = sumn*sumn;
          sd32[i1] = 9.0f*sumd;
        }
      }
    }});
    return new float[][][][]{sn,sd};
  }





  private static final double _small = 0.001;
  private static final int _niter = 1000;
  private static final LocalDiffusionKernel _ldk = 
    //new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D22);
    new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D71);
  private static final LocalSmoothingFilter _lsf = 
    new LocalSmoothingFilter(_small,_niter,_ldk);


}

