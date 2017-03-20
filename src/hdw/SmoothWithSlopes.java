package hdw;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import vec.*;
import util.*;

public class SmoothWithSlopes {

  public float[][][] getHorizons(int nh, int dh, int fh, float[][][] ux) {
    int n3 = ux.length;
    int n2 = ux[0].length;
    int n1 = ux[0][0].length;
    float[][][] hs = new float[nh][n3][n2];
    int ik = 0;
    for (int ih=fh; ih<n1; ih+=dh) {
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        hs[ik][i3][i2] = ux[i3][i2][ih]+ih;
      }}
      ik++;
    }
    return hs;
  }

  public float[][][] flatten(float[][][] ux, float[][][] fx) {
    int n3 = ux.length;
    int n2 = ux[0].length;
    int n1 = ux[0][0].length;
    float[][][] gx = new float[n3][n2][n1];
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float[] u = new float[n1];
      for (int i1=0; i1<n1; ++i1)
        u[i1] = i1+ux[i3][i2][i1];
      si.interpolate(n1,1.0,0.0,fx[i3][i2],n1,u,gx[i3][i2]);
    }}
    return gx;
  }

  public float[][][] shiftsToRgt(float[][][] ux) {
    final int n3 = ux.length;
    final int n2 = ux[0].length;
    final int n1 = ux[0][0].length;
    final float[][][] tx = new float[n3][n2][n1];
    final float[] y = rampfloat(0,1,n1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] x = new float[n1];
      for (int i2=0; i2<n2; ++i2) {
        float[] tx32 = tx[i3][i2];
        CubicInterpolator ci = new CubicInterpolator(x,y);
        for (int i1=0; i1<n1; ++i1)
          tx32[i1] = ci.interpolate(i1);
      }
    }});
    return tx;
  }

  public float[][][] smooth(
    final float[][][] wx, final float[][][] p2, 
    final float[][][] p3, final float[][][] ux) {
    final int n3 = wx.length;
    final int n2 = wx[0].length;
    final int n1 = wx[0][0].length;
    final float[][][] us = new float[n3][n2][n1];
    //Parallel.loop(n1,new Parallel.LoopInt() {
    //public void compute(int i1) {
    for (int i1=0; i1<n1; ++i1) {
      System.out.println("i1="+i1);
      float[][] hx1 = new float[n3][n2];
      float[][] p21 = new float[n3][n2];
      float[][] p31 = new float[n3][n2];
      float[][] wx1 = new float[n3][n2];
      float[][] wh1 = new float[n3][n2];
      for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        hx1[i3][i2] = ux[i3][i2][i1]+i1;
      for (int iter=0; iter<5; iter++) {
        zero(p21);
        zero(p31);
        zero(wx1);
        zero(wh1);
        updateSlopesAndWeights(p2,p3,wx,hx1,p21,p31,wx1,wh1);
        hx1 = smooth(wx1,wh1,p21,p31,hx1);
      }
      for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        us[i3][i2][i1] = hx1[i3][i2]-i1;
    //}});
    }
    return us;
  }

  private static void updateSlopesAndWeights (
    float[][][] p, float[][][] q, float[][][] ep,
    float[][] surf, float[][] pi1, float[][] qi1,float[][] wi1, float[][] wh1)
  {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    SincInterpolator.Extrapolation extrap = SincInterpolator.Extrapolation.CONSTANT;
    SincInterpolator psi = new SincInterpolator();
    SincInterpolator qsi = new SincInterpolator();
    SincInterpolator wsi = new SincInterpolator();
    psi.setExtrapolation(extrap);
    qsi.setExtrapolation(extrap);
    wsi.setExtrapolation(extrap);
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        double x2i = (double)i2;
        double x3i = (double)i3;
        double x1i = (double)surf[i3][i2];
        wi1[i3][i2] = 
          wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,ep,x1i,x2i,x3i);
        pi1[i3][i2] = 
          psi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, p,x1i,x2i,x3i);
	      qi1[i3][i2] = 
          qsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, q,x1i,x2i,x3i);
        if (x1i<0||x1i>=n1) {
          wh1[i3][i2] = 0f;
        } else {
          wh1[i3][i2] = 1f;
        }
        if(x1i<0f)  surf[i3][i2] = 0f;
        if(x1i>=n1) surf[i3][i2] = n1-1f;
      }
    }
  }

  public float[][] smooth(
    float[][] wx, float[][] wh, float[][] p2, float[][] p3, float[][] hx) {
    int n2 = hx.length;
    int n1 = hx[0].length;
    float[][] b = new float[n2][n1];
    float[][] r = copy(hx);//new float[n2][n1];
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(_sigma1,_sigma2,wx);
    A2 a2 = new A2(smoother2,wx,wh);
    CgSolver cs = new CgSolver(0.001,100);
    makeRhs(hx,wx,wh,p2,p3,b);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    smoother2.apply(r);
    return r;
  }


  private float _sigma1 = 10;
  private float _sigma2 = 10;

    // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(float sigma1, float sigma2, float[][] wp) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,_wp,x);
      smooth1(_sigma1,_wp,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,_wp,x);
      smooth2(_sigma2,_wp,x);
    }
    private float[][] _wp;
    private float _sigma1,_sigma2;
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


    // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n2 = x.length;
    int n1 = x[0].length;
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      float[] xt = new float[n1];
      lsf.apply(c,s[i2],x[i2],xt);
      x[i2] = xt;
    }
  }



  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }



  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] wp, float[][] wh) 
    {
      _s2 = s2;
      _wp = wp;
      _wh = wh;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      _s2.apply(z);
      v2y.zero();
      applyLhs(_wp,_wh,z,y);
      _s2.applyTranspose(y);
    }
    private Smoother2 _s2;
    private float[][] _wp;
    private float[][] _wh;
  }

  private static void applyLhs( float[][] wp, float[][] wh, 
    float[][] x, float[][] y) {
    zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        //if(wpi<0.05f) {wpi=0.05f;}
        float wps = wpi*wpi*0.25f;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = (xa+xb);
        float x2 = (xa-xb);
        float y1 = wps*x1;
        float y2 = wps*x2;
        float ya = (y1+y2);
        float yb = (y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
        y[i2][i1] += wh[i2][i1]*x[i2][i1]*0.005f;
  }


  private static void makeRhs(
    float[][] x, float[][] wp, float[][] wh, 
    float[][] p2, float[][] p3, float[][] y) 
  {
    zero(y);
    int n2 = y.length;
    int n1 = y[0].length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        //if(wpi<0.05f) {wpi=0.05f;}
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float wps = wpi*wpi*0.5f;
        float y1 = wps*p2i;
        float y2 = wps*p3i;
        float ya = (y1+y2);
        float yb = (y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
    for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1)
      y[i2][i1] += wh[i2][i1]*x[i2][i1]*0.005f;
 
  }



}
