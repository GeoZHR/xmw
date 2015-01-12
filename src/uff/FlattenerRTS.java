package uff;

import vec.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates shift vectors to flatten features in 2D and 3D images.
 * In 2D, the shift vector field has two components r1(x1,x3) and r3(x1,x3).
 * In 3D, the vector field has three components r1(x1,x2,x3), r2(x1,x2,x3),
 * and r3(x1,x2,x3).
 * @author Xinming, Simon Luo and Dave Hale
 * @version 2012.02.01
 */
public class FlattenerRTS {

  /**
   * Constructs a flattener.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   */
  public FlattenerRTS(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  public void setIters(int inner, int outer) {
    _inner=inner;
    _outer=outer;
  }

  /**
   * Estimates shift vectors for a 2D image.
   * @param p array of parameters {u1,u2,el,a}.
   * @return array of shifts {r1,r2}.
   */
  public float[][][] findShifts(float[][][] p) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][][] p0 = copy(p);
    float[][][] b = new float[2][n2][n1];
    float[][][] r = new float[2][n2][n1];
    float[][][] rc = new float[2][n2][n1];
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    Smoother2 s2 = new Smoother2(_sigma1,_sigma2,p[2]);
    CgSolver cg = new CgSolver(_small,_inner);
    A2 ma = new A2(_epsilon,s2,p);
    for (int outer=0; outer<_outer; ++outer) {
      if (outer>0) {
        copy(r,rc);
        s2.apply(rc);
        updateParameters(rc,p0,p);
        vb.zero();
      }
      makeRhs(p,b);
      s2.applyTranspose(b);
      int inner = cg.solve(ma,vb,vr).niter;
      if (inner==0) break;
    }
    s2.apply(r);
    return r;
  }

  /**
   * Estimates shift vectors for a 3D image.
   * @param p array of parameters {u1,u2,u3,ep,a}.
   * @return array of shifts {r1,r2,r3}.
   */
  public float[][][][] findShifts(
    float[][] cp, float[][] cm, float[][] dr, 
    int[][][] fm, float[] fl, float[][][][] p) {
    int n1 = p[0][0][0].length;
    int n2 = p[0][0].length;
    int n3 = p[0].length;
    float[][][][] p0 = scopy(p);
    float[][][][] b = new float[3][n3][n2][n1];
    float[][][][] r = new float[3][n3][n2][n1];
    float[][][][] rc = new float[3][n3][n2][n1];
    VecArrayFloat4 vr = new VecArrayFloat4(r);
    VecArrayFloat4 vb = new VecArrayFloat4(b);
    Smoother3 s3 = new Smoother3(_sigma1,_sigma2,_sigma2,p[3]);
    CgSolver cg = new CgSolver(_small,_inner);
    A3 ma = new A3(_epsilon,s3,cp,cm,fl,p);
    for (int outer=0; outer<_outer; ++outer) {
      System.out.println("iter="+outer);
      if (outer>0) {
        scopy(r,rc);
        s3.applyOriginal(rc);
        updateParameters(rc,p0,p);
        vb.zero();
      }
      makeRhs(cp,cm,dr,fl,p,b);
      s3.applyTranspose(b);
      int inner = cg.solve(ma,vb,vr).niter;
      if (inner==0) break;
    }
    s3.applyOriginal(r);
    cleanShifts(fm,r);
    cleanShifts(fm,r);
    cleanShifts(fm,r);
    return r;
  }

  private void cleanShifts(int[][][] fm, float[][][][] r) {
    if(fm==null) {return;}
    int n4 = r.length;
    for (int i4=0; i4<n4; ++i4) {
      cleanShifts(fm,r[i4]);
    }
  }

  private void cleanShifts(int[][][] fm, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int m2 = n2-1;
    int m3 = n3-1;
    int fn = fm[0].length;
    for (int fi=0; fi<fn; ++fi) {
      int np = fm[0][fi].length;
      for (int ip=0; ip<np; ++ip) {
        int i1 = fm[0][fi][ip];
        int i2 = fm[1][fi][ip];
        int i3 = fm[2][fi][ip];
        int i2m1 = i2-1; if(i2m1<0) i2m1=0;
        int i3m1 = i3-1; if(i3m1<0) i3m1=0;
        int i2m2 = i2-2; if(i2m2<0) i2m2=0;
        int i3m2 = i3-2; if(i3m2<0) i3m2=0;
        int i2p1 = i2+1; if(i2p1>m2) i2p1=m2;
        int i3p1 = i3+1; if(i3p1>m3) i3p1=m3;
        int i2p2 = i2+2; if(i2p2>m2) i2p2=m2;
        int i3p2 = i3+2; if(i3p2>m3) i3p2=m3;
        x[i3][i2  ][i1] = x[i3][i2m2][i1];
        x[i3][i2m1][i1] = x[i3][i2m2][i1];
        x[i3][i2p1][i1] = x[i3][i2p2][i1];
        x[i3m1][i2][i1] = x[i3m2][i2][i1];
        x[i3p1][i2][i1] = x[i3p2][i2][i1];
      }
    }
  }


  /**
   * Applies shifts using sinc interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param r input array {r1,r2} of shifts.
   * @param f input image.
   * @param g output shifted image.
   */
  public void applyShifts(
    float[][][] r, float[][] f, float[][] g)
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] r1 = r[0], r2 = r[1];
    SincInterpolator si = new SincInterpolator();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = si.interpolate(
          n1,1.0,0.0,n2,1.0,0.0,f,i1-r1[i2][i1],i2-r2[i2][i1]);
      }
    }
  }

  /**
   * Applies shifts using sinc interpolation.
   * @param r input array {r1,r2,r3} of shifts.
   * @param f input image.
   * @param g output shifted image.
   */
  public void applyShifts(
    float[][][][] r, float[][][] f, float[][][] g)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] ff = f;
    final float[][][] gf = g;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          gf[i3][i2][i1] = si.interpolate(
            n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,
            ff,i1-r1[i3][i2][i1],i2-r2[i3][i2][i1],i3-r3[i3][i2][i1]);
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float W0 = 1.000f; // flatness
  private static final float W1 = 0.011f; // distance
  private static final float W2 = 0.001f; // thickness

  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _epsilon = 0.000f; // damping for stability?
  private float _small = 0.010f; // stop CG iterations if residuals are small
  //private float _small = 0.001f; // stop CG iterations if residuals are small
  private int _inner = 100; // maximum number of inner CG iterations
  private int _outer = 1; // maximum number of outer iterations

  private static class A2 implements CgSolver.A {
    A2(float epsilon, Smoother2 smoother, float[][][] p) {
      _epsilon = epsilon;
      _smoother = smoother;
      _p = p;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      VecArrayFloat3 v3z = v3x.clone();
      v3y.zero();
      float[][][] x = v3z.getArray();
      float[][][] y = v3y.getArray();
      _smoother.apply(x);
      applyLhs(_p,x,y);
      _smoother.applyTranspose(y);
      if (_epsilon>0.0f)
        v3y.add(1.0,v3x,_epsilon*_epsilon);
    }
    private float _epsilon = 0.0f;
    private Smoother2 _smoother;
    private float[][][] _p;
  }
  private static class A3 implements CgSolver.A {
    A3(float epsilon, Smoother3 smoother, 
       float[][] cp, float[][] cm, float[] fl, float[][][][] p) {
      _epsilon = epsilon;
      _smoother = smoother;
      _cp = cp;
      _cm = cm;
      _fl = fl;
      _p = p;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat4 v4x = (VecArrayFloat4)vx;
      VecArrayFloat4 v4y = (VecArrayFloat4)vy;
      VecArrayFloat4 v4z = v4x.clone();
      v4y.zero();
      float[][][][] x = v4z.getArray();
      float[][][][] y = v4y.getArray();
      _smoother.applyOriginal(x);
      applyLhs(_cp,_cm,_fl,_p,x,y);
      _smoother.applyTranspose(y);
      if (_epsilon>0.0f)
        v4y.add(1.0,v4x,_epsilon*_epsilon);
    }
    private float[] _fl;
    private float[][] _cp;
    private float[][] _cm;
    private float _epsilon = 0.0f;
    private Smoother3 _smoother;
    private float[][][][] _p;
  }

  // 2D RHS
  private static void makeRhs(float[][][] p, float[][][] y) {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    float[][] y1 = y[0]; float[][] y2 = y[1];
    float[][] u1 = p[0]; float[][] u2 = p[1]; float[][] el = p[2];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float eli = el[i2][i1]*0.5f;
        float u1m = 1.0f-u1i;
        float b1 = W2*eli*u1m;
        float b2 = W2*eli*u2i;
        float b3 = W0*eli*u2i;
        float b4 = W1*eli*u1m;
        float y11 = b1*u1i+b2*u2i;
        float y12 = b3*u1i+b4*u2i;
        float y21 = b1*u2i-b2*u1i;
        float y22 = b3*u2i-b4*u1i;
        float y1a = y11+y12;
        float y1b = y11-y12;
        float y2a = y21+y22;
        float y2b = y21-y22;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }

  // 2D LHS
  private static void applyLhs(
    float[][][] p, float[][][] x, float[][][] y)
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    float[][] x1 = x[0]; float[][] x2 = x[1];
    float[][] y1 = y[0]; float[][] y2 = y[1];
    float[][] u1 = p[0]; float[][] u2 = p[1]; float[][] el = p[2];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float x100 = x1[i2  ][i1  ];
        float x101 = x1[i2  ][i1-1];
        float x110 = x1[i2-1][i1  ];
        float x111 = x1[i2-1][i1-1];
        float x200 = x2[i2  ][i1  ];
        float x201 = x2[i2  ][i1-1];
        float x210 = x2[i2-1][i1  ];
        float x211 = x2[i2-1][i1-1];
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float eli = el[i2][i1]*0.25f;
        float x1a = x100-x111;
        float x1b = x101-x110;
        float x2a = x200-x211;
        float x2b = x201-x210;
        float x11 = x1a-x1b;
        float x12 = x1a+x1b;
        float x21 = x2a-x2b;
        float x22 = x2a+x2b;

        float b1 = W2*eli*( x11*u1i+x21*u2i);
        float b2 = W2*eli*(-x21*u1i+x11*u2i);
        float b3 = W0*eli*( x12*u1i+x22*u2i);
        float b4 = W1*eli*(-x22*u1i+x12*u2i);

        float y11 = b1*u1i+b2*u2i;
        float y12 = b3*u1i+b4*u2i;
        float y21 = b1*u2i-b2*u1i;
        float y22 = b3*u2i-b4*u1i;
        float y1a = y11+y12;
        float y1b = y11-y12;
        float y2a = y21+y22;
        float y2b = y21-y22;
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }

  private static void makeRhs(
    float[][] cp, float[][] cm, float[][] dr, float[] fl, 
    final float[][][][] p, final float[][][][] y)
  { 
    final int n3 = y[0].length;
    // i3 = 1,3,5,...
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      makeRhsSlice3(i3,p,y);
    }});
    // i3 = 2,4,6,...
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      makeRhsSlice3(i3,p,y);
    }});
    screenRhs(cp,cm,dr,fl,y);
  }

  private static void screenRhs(
    float[][] cp, float[][] cm, float[][] dr, float[] fl, float[][][][] y) {
    int nc = cp[0].length;
    for (int ic=0; ic<nc; ++ic) {
      float dr1 = dr[0][ic];
      float dr2 = dr[1][ic];
      float dr3 = dr[2][ic];
      int i1p = (int)cp[0][ic];
      int i2p = (int)cp[1][ic];
      int i3p = (int)cp[2][ic];
      int i1m = (int)cm[0][ic];
      int i2m = (int)cm[1][ic];
      int i3m = (int)cm[2][ic];
      float fls = fl[ic]*fl[ic];
      y[0][i3p][i2p][i1p] += dr1*fls;
      y[1][i3p][i2p][i1p] += dr2*fls;
      y[2][i3p][i2p][i1p] += dr3*fls;
      y[0][i3m][i2m][i1m] -= dr1*fls;
      y[1][i3m][i2m][i1m] -= dr2*fls;
      y[2][i3m][i2m][i1m] -= dr3*fls;
    }
  }

  private static void makeRhsSlice3(
    int i3, float[][][][] p, float[][][][] y)
  { 
    int n1 = y[0][0][0].length;
    int n2 = y[0][0].length;
    float[][][] y1 = y[0]; float[][][] y2 = y[1]; float[][][] y3 = y[2];
    float[][][] u1 = p[0]; float[][][] u2 = p[1]; float[][][] u3 = p[2];
    float[][][] ep = p[3];
    for (int i2=1; i2<n2; ++i2) {
      float[] y100 = y1[i3  ][i2  ];
      float[] y101 = y1[i3  ][i2-1];
      float[] y110 = y1[i3-1][i2  ];
      float[] y111 = y1[i3-1][i2-1];
      float[] y200 = y2[i3  ][i2  ];
      float[] y201 = y2[i3  ][i2-1];
      float[] y210 = y2[i3-1][i2  ];
      float[] y211 = y2[i3-1][i2-1];
      float[] y300 = y3[i3  ][i2  ];
      float[] y301 = y3[i3  ][i2-1];
      float[] y310 = y3[i3-1][i2  ];
      float[] y311 = y3[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float epi = ep[i3][i2][i1]*0.25f;
        float u1p = u1i+1.0f;
        float u1u1p = u1i*u1p;
        float u2u1p = u2i*u1p;
        float u3u1p = u3i*u1p;
        float u2u3 = u2i*u3i;
        float u2smu1p = u2i*u2i-u1p;
        float u3smu1p = u3i*u3i-u1p;

        float b1 = W2*epi*(u1i*u1i-1.0f);
        float b2 = W2*epi*(u2u1p);
        float b3 = W2*epi*(u3u1p);
        float b4 = W0*epi*(u2u1p);
        float b5 = W1*epi*(u2i*u2i);
        float b6 = W1*epi*(u2u3);
        float b7 = W0*epi*(u3u1p);
        float b8 = W1*epi*(u2u3);
        float b9 = W1*epi*(u3i*u3i);

        float y11 = b1*u1u1p + b2*u2u1p   + b3*u3u1p;
        float y12 = b4*u1u1p + b5*u2u1p   + b6*u3u1p;
        float y13 = b7*u1u1p + b8*u2u1p   + b9*u3u1p;
        float y21 = b1*u2u1p + b2*u2smu1p + b3*u2u3;
        float y22 = b4*u2u1p + b5*u2smu1p + b6*u2u3;
        float y23 = b7*u2u1p + b8*u2smu1p + b9*u2u3;
        float y31 = b1*u3u1p + b2*u2u3    + b3*u3smu1p;
        float y32 = b4*u3u1p + b5*u2u3    + b6*u3smu1p;
        float y33 = b7*u3u1p + b8*u2u3    + b9*u3smu1p;

        float y1a = y11+y12+y13;
        float y1b = y11-y12+y13;
        float y1c = y11+y12-y13;
        float y1d = y11-y12-y13;
        float y2a = y21+y22+y23;
        float y2b = y21-y22+y23;
        float y2c = y21+y22-y23;
        float y2d = y21-y22-y23;
        float y3a = y31+y32+y33;
        float y3b = y31-y32+y33;
        float y3c = y31+y32-y33;
        float y3d = y31-y32-y33;
        y100[i1 ] += y1a;
        y100[i1m] -= y1d;
        y101[i1 ] += y1b;
        y101[i1m] -= y1c;
        y110[i1 ] += y1c;
        y110[i1m] -= y1b;
        y111[i1 ] += y1d;
        y111[i1m] -= y1a;
        y200[i1 ] += y2a;
        y200[i1m] -= y2d;
        y201[i1 ] += y2b;
        y201[i1m] -= y2c;
        y210[i1 ] += y2c;
        y210[i1m] -= y2b;
        y211[i1 ] += y2d;
        y211[i1m] -= y2a;
        y300[i1 ] += y3a;
        y300[i1m] -= y3d;
        y301[i1 ] += y3b;
        y301[i1m] -= y3c;
        y310[i1 ] += y3c;
        y310[i1m] -= y3b;
        y311[i1 ] += y3d;
        y311[i1m] -= y3a;
      }
    }
  }

  private static void applyLhs(float[][] cp, float[][] cm, float[] fl, 
    final float[][][][] p, final float[][][][] x, final float[][][][] y)
  { 
    final int n3 = y[0].length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,p,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,p,x,y);
    }});
    screenLhs(cp,cm,fl,x,y);
  }

  private static void screenLhs(
    float[][] cp, float[][] cm, float[]fl, float[][][][] x, float[][][][] y) 
  {
    int nc = cp[0].length;
    for (int ic=0; ic<nc; ++ic) {
      int i1p = (int)cp[0][ic];
      int i2p = (int)cp[1][ic];
      int i3p = (int)cp[2][ic];
      int i1m = (int)cm[0][ic];
      int i2m = (int)cm[1][ic];
      int i3m = (int)cm[2][ic];
      float fls = fl[ic]*fl[ic];
      float x1p = x[0][i3p][i2p][i1p];
      float x2p = x[1][i3p][i2p][i1p];
      float x3p = x[2][i3p][i2p][i1p];
      float x1m = x[0][i3m][i2m][i1m];
      float x2m = x[1][i3m][i2m][i1m];
      float x3m = x[2][i3m][i2m][i1m];
      float dx1 = (x1p-x1m)*fls;
      float dx2 = (x2p-x2m)*fls;
      float dx3 = (x3p-x3m)*fls;
      y[0][i3p][i2p][i1p] += dx1;
      y[1][i3p][i2p][i1p] += dx2;
      y[2][i3p][i2p][i1p] += dx3;
      y[0][i3m][i2m][i1m] -= dx1;
      y[1][i3m][i2m][i1m] -= dx2;
      y[2][i3m][i2m][i1m] -= dx3;
    }
  }

  private static void applyLhsSlice3(
    int i3, float[][][] wp, float[][][][] x, float[][][][] y) 
  {
    int n2 = x[0][0].length;
    int n1 = x[0][0][0].length;
    float[][][] x1 = x[0]; float[][][] x2 = x[1]; float[][][] x3 = x[2];
    float[][][] y1 = y[0]; float[][][] y2 = y[1]; float[][][] y3 = y[2];

    for (int i2=1; i2<n2; ++i2) {
      float[] x100 = x1[i3  ][i2  ];
      float[] x101 = x1[i3  ][i2-1];
      float[] x110 = x1[i3-1][i2  ];
      float[] x111 = x1[i3-1][i2-1];

      float[] x200 = x2[i3  ][i2  ];
      float[] x201 = x2[i3  ][i2-1];
      float[] x210 = x2[i3-1][i2  ];
      float[] x211 = x2[i3-1][i2-1];

      float[] x300 = x3[i3  ][i2  ];
      float[] x301 = x3[i3  ][i2-1];
      float[] x310 = x3[i3-1][i2  ];
      float[] x311 = x3[i3-1][i2-1];

      float[] y100 = y1[i3  ][i2  ];
      float[] y101 = y1[i3  ][i2-1];
      float[] y110 = y1[i3-1][i2  ];
      float[] y111 = y1[i3-1][i2-1];

      float[] y200 = y2[i3  ][i2  ];
      float[] y201 = y2[i3  ][i2-1];
      float[] y210 = y2[i3-1][i2  ];
      float[] y211 = y2[i3-1][i2-1];

      float[] y300 = y3[i3  ][i2  ];
      float[] y301 = y3[i3  ][i2-1];
      float[] y310 = y3[i3-1][i2  ];
      float[] y311 = y3[i3-1][i2-1];

      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
        float wps = wpi*wpi;
        float x1a = 0.0f;
        float x1b = 0.0f;
        float x1c = 0.0f;
        float x1d = 0.0f;

        float x2a = 0.0f;
        float x2b = 0.0f;
        float x2c = 0.0f;
        float x2d = 0.0f;

        float x3a = 0.0f;
        float x3b = 0.0f;
        float x3c = 0.0f;
        float x3d = 0.0f;

        x1a += x100[i1 ];
        x1d -= x100[i1m];
        x1b += x101[i1 ];
        x1c -= x101[i1m];
        x1c += x110[i1 ];
        x1b -= x110[i1m];
        x1d += x111[i1 ];
        x1a -= x111[i1m];

        x2a += x200[i1 ];
        x2d -= x200[i1m];
        x2b += x201[i1 ];
        x2c -= x201[i1m];
        x2c += x210[i1 ];
        x2b -= x210[i1m];
        x2d += x211[i1 ];
        x2a -= x211[i1m];

        x3a += x300[i1 ];
        x3d -= x300[i1m];
        x3b += x301[i1 ];
        x3c -= x301[i1m];
        x3c += x310[i1 ];
        x3b -= x310[i1m];
        x3d += x311[i1 ];
        x3a -= x311[i1m];

        float y11 = 0.25f*(x1a+x1b+x1c+x1d)*wps;
        float y12 = 0.25f*(x1a-x1b+x1c-x1d)*wps;
        float y13 = 0.25f*(x1a+x1b-x1c-x1d)*wps;

        float y21 = 0.25f*(x2a+x2b+x2c+x2d)*wps;
        float y22 = 0.25f*(x2a-x2b+x2c-x2d)*wps;
        float y23 = 0.25f*(x2a+x2b-x2c-x2d)*wps;

        float y31 = 0.25f*(x3a+x3b+x3c+x3d)*wps;
        float y32 = 0.25f*(x3a-x3b+x3c-x3d)*wps;
        float y33 = 0.25f*(x3a+x3b-x3c-x3d)*wps;

        float y1a = 0.25f*(y11+y12+y13);
        float y1b = 0.25f*(y11-y12+y13);
        float y1c = 0.25f*(y11+y12-y13);
        float y1d = 0.25f*(y11-y12-y13);

        float y2a = 0.25f*(y21+y22+y23);
        float y2b = 0.25f*(y21-y22+y23);
        float y2c = 0.25f*(y21+y22-y23);
        float y2d = 0.25f*(y21-y22-y23);

        float y3a = 0.25f*(y31+y32+y33);
        float y3b = 0.25f*(y31-y32+y33);
        float y3c = 0.25f*(y31+y32-y33);
        float y3d = 0.25f*(y31-y32-y33);

        y100[i1 ] += y1a;
        y100[i1m] -= y1d;
        y101[i1 ] += y1b;
        y101[i1m] -= y1c;
        y110[i1 ] += y1c;
        y110[i1m] -= y1b;
        y111[i1 ] += y1d;
        y111[i1m] -= y1a;

        y200[i1 ] += y2a;
        y200[i1m] -= y2d;
        y201[i1 ] += y2b;
        y201[i1m] -= y2c;
        y210[i1 ] += y2c;
        y210[i1m] -= y2b;
        y211[i1 ] += y2d;
        y211[i1m] -= y2a;

        y300[i1 ] += y3a;
        y300[i1m] -= y3d;
        y301[i1 ] += y3b;
        y301[i1m] -= y3c;
        y310[i1 ] += y3c;
        y310[i1m] -= y3b;
        y311[i1 ] += y3d;
        y311[i1m] -= y3a;
      }
    }
  }


  // 3D LHS
  private static void applyLhsSlice3(
    int i3, float[][][][] p, float[][][][] x, float[][][][] y)
  {
    int n1 = y[0][0][0].length;
    int n2 = y[0][0].length;
    float[][][] x1 = x[0]; float[][][] x2 = x[1]; float[][][] x3 = x[2];
    float[][][] y1 = y[0]; float[][][] y2 = y[1]; float[][][] y3 = y[2];
    float[][][] u1 = p[0]; float[][][] u2 = p[1]; float[][][] u3 = p[2];
    float[][][] ep = p[3];
    for (int i2=1; i2<n2; ++i2) {
      float[] x100 = x1[i3  ][i2  ];
      float[] x101 = x1[i3  ][i2-1];
      float[] x110 = x1[i3-1][i2  ];
      float[] x111 = x1[i3-1][i2-1];
      float[] x200 = x2[i3  ][i2  ];
      float[] x201 = x2[i3  ][i2-1];
      float[] x210 = x2[i3-1][i2  ];
      float[] x211 = x2[i3-1][i2-1];
      float[] x300 = x3[i3  ][i2  ];
      float[] x301 = x3[i3  ][i2-1];
      float[] x310 = x3[i3-1][i2  ];
      float[] x311 = x3[i3-1][i2-1];
      float[] y100 = y1[i3  ][i2  ];
      float[] y101 = y1[i3  ][i2-1];
      float[] y110 = y1[i3-1][i2  ];
      float[] y111 = y1[i3-1][i2-1];
      float[] y200 = y2[i3  ][i2  ];
      float[] y201 = y2[i3  ][i2-1];
      float[] y210 = y2[i3-1][i2  ];
      float[] y211 = y2[i3-1][i2-1];
      float[] y300 = y3[i3  ][i2  ];
      float[] y301 = y3[i3  ][i2-1];
      float[] y310 = y3[i3-1][i2  ];
      float[] y311 = y3[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float x1000 = x100[i1 ];
        float x1001 = x100[i1m];
        float x1010 = x101[i1 ];
        float x1011 = x101[i1m];
        float x1100 = x110[i1 ];
        float x1101 = x110[i1m];
        float x1110 = x111[i1 ];
        float x1111 = x111[i1m];
        float x2000 = x200[i1 ];
        float x2001 = x200[i1m];
        float x2010 = x201[i1 ];
        float x2011 = x201[i1m];
        float x2100 = x210[i1 ];
        float x2101 = x210[i1m];
        float x2110 = x211[i1 ];
        float x2111 = x211[i1m];
        float x3000 = x300[i1 ];
        float x3001 = x300[i1m];
        float x3010 = x301[i1 ];
        float x3011 = x301[i1m];
        float x3100 = x310[i1 ];
        float x3101 = x310[i1m];
        float x3110 = x311[i1 ];
        float x3111 = x311[i1m];
        float x1a = x1000-x1111;
        float x1b = x1001-x1110;
        float x1c = x1010-x1101;
        float x1d = x1100-x1011;
        float x2a = x2000-x2111;
        float x2b = x2001-x2110;
        float x2c = x2010-x2101;
        float x2d = x2100-x2011;
        float x3a = x3000-x3111;
        float x3b = x3001-x3110;
        float x3c = x3010-x3101;
        float x3d = x3100-x3011;
        float x11 = x1a-x1b+x1c+x1d;
        float x12 = x1a+x1b-x1c+x1d;
        float x13 = x1a+x1b+x1c-x1d;
        float x21 = x2a-x2b+x2c+x2d;
        float x22 = x2a+x2b-x2c+x2d;
        float x23 = x2a+x2b+x2c-x2d;
        float x31 = x3a-x3b+x3c+x3d;
        float x32 = x3a+x3b-x3c+x3d;
        float x33 = x3a+x3b+x3c-x3d;

        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float epi = ep[i3][i2][i1]*0.0625f;
        float u1p = u1i+1.0f;
        float u1u1p = u1i*u1p;
        float u2u1p = u2i*u1p;
        float u3u1p = u3i*u1p;
        float u2u3 = u2i*u3i;
        float u2smu1p = u2i*u2i-u1p;
        float u3smu1p = u3i*u3i-u1p;

        float b1 = W2*epi* (x11*u1u1p + x21*u2u1p   + x31*u3u1p);
        float b2 = W2*epi* (x11*u2u1p + x21*u2smu1p + x31*u2u3);
        float b3 = W2*epi* (x11*u3u1p + x21*u2u3    + x31*u3smu1p);
        float b4 = W0*epi* (x12*u1u1p + x22*u2u1p   + x32*u3u1p);
        float b5 = W1*epi* (x12*u2u1p + x22*u2smu1p + x32*u2u3);
        float b6 = W1*epi* (x12*u3u1p + x22*u2u3    + x32*u3smu1p);
        float b7 = W0*epi* (x13*u1u1p + x23*u2u1p   + x33*u3u1p);
        float b8 = W1*epi* (x13*u2u1p + x23*u2smu1p + x33*u2u3);
        float b9 = W1*epi* (x13*u3u1p + x23*u2u3    + x33*u3smu1p);

        float y11 = b1*u1u1p + b2*u2u1p   + b3*u3u1p;
        float y12 = b4*u1u1p + b5*u2u1p   + b6*u3u1p;
        float y13 = b7*u1u1p + b8*u2u1p   + b9*u3u1p;
        float y21 = b1*u2u1p + b2*u2smu1p + b3*u2u3;
        float y22 = b4*u2u1p + b5*u2smu1p + b6*u2u3;
        float y23 = b7*u2u1p + b8*u2smu1p + b9*u2u3;
        float y31 = b1*u3u1p + b2*u2u3    + b3*u3smu1p;
        float y32 = b4*u3u1p + b5*u2u3    + b6*u3smu1p;
        float y33 = b7*u3u1p + b8*u2u3    + b9*u3smu1p;

        float y1a = y11+y12+y13;
        float y1b = y11-y12+y13;
        float y1c = y11+y12-y13;
        float y1d = y11-y12-y13;
        float y2a = y21+y22+y23;
        float y2b = y21-y22+y23;
        float y2c = y21+y22-y23;
        float y2d = y21-y22-y23;
        float y3a = y31+y32+y33;
        float y3b = y31-y32+y33;
        float y3c = y31+y32-y33;
        float y3d = y31-y32-y33;
        y100[i1 ] += y1a;
        y100[i1m] -= y1d;
        y101[i1 ] += y1b;
        y101[i1m] -= y1c;
        y110[i1 ] += y1c;
        y110[i1m] -= y1b;
        y111[i1 ] += y1d;
        y111[i1m] -= y1a;
        y200[i1 ] += y2a;
        y200[i1m] -= y2d;
        y201[i1 ] += y2b;
        y201[i1m] -= y2c;
        y210[i1 ] += y2c;
        y210[i1m] -= y2b;
        y211[i1 ] += y2d;
        y211[i1m] -= y2a;
        y300[i1 ] += y3a;
        y300[i1m] -= y3d;
        y301[i1 ] += y3b;
        y301[i1m] -= y3c;
        y310[i1 ] += y3c;
        y310[i1m] -= y3b;
        y311[i1 ] += y3d;
        y311[i1m] -= y3a;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  
  private static float[][][][] scopy(float[][][][] x) {
    int n1 = x[0][0][0].length;
    int n2 = x[0][0].length;
    int n3 = x[0].length;
    int n4 = x.length;
    float[][][][] y = new float[n4][n3][n2][n1];
    scopy(x,y);
    return y;
  }

  private static void scopy(float[][][][] x, float[][][][] y) {
    int n4 = x.length;
    for (int i4=0; i4<n4; i4++)
      copy(x[i4],y[i4]);
  }

  private static void updateParameters(
    float[][][] r, float[][][] p, float[][][] q)
  {
    for (int i=0; i<3; ++i) // shift normal vectors and ep
      applyShiftsLinear(r,p[i],q[i]);
    normalize(q[0],q[1]);
  }

  private static void updateParameters(
    float[][][][] r, float[][][][] p, float[][][][] q) {
    for (int i=0; i<4; ++i) // shift normal vectors and ep
      applyShiftsLinear(r,p[i],q[i]);
    normalize(q[0],q[1],q[2]);
  }

  // Normalize vectors
  private static void normalize(float[][] a1, float[][] a2) {
    int n2 = a1.length;
    int n1 = a1[0].length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float a1i = a1[i2][i1];
        float a2i = a2[i2][i1];
        float sca = 1.0f/sqrt(a1i*a1i+a2i*a2i);
        a1[i2][i1] *= sca;
        a2[i2][i1] *= sca;
      }
    }
  }
  private static void normalize(
    final float[][][] a1, final float[][][] a2, final float[][][] a3)
  {
    final int n1 = a1[0][0].length;
    final int n2 = a1[0].length;
    final int n3 = a1.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float a1i = a1[i3][i2][i1];
          float a2i = a2[i3][i2][i1];
          float a3i = a3[i3][i2][i1];
          float sca = 1.0f/sqrt(a1i*a1i+a2i*a2i+a3i*a3i);
          a1[i3][i2][i1] *= sca;
          a2[i3][i2][i1] *= sca;
          a3[i3][i2][i1] *= sca;
        }
      }
    }});
  }

  /**
   * Applies the specified shifts using linear interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param f input array to which shifts are to be applied.
   * @param r array {r1,r2} of shifts.
   * @return array with shifts applied.
   */
  private static void applyShiftsLinear(
    float[][][] r, float[][] f, float[][] g)
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] r1 = r[0], r2 = r[1];
    LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,f);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = li.interpolate(i1-r1[i2][i1],i2-r2[i2][i1]);
      }
    }
  }

  /**
   * Applies the specified shifts using linear interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param f input array to which shifts are to be applied.
   * @param r array {r1,r2,r3} of shifts.
   * @return array with shifts applied.
   */
  private static void applyShiftsLinear(
    float[][][][] r, float[][][] f, float[][][] g)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] gf = g;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          gf[i3][i2][i1] = li.interpolate(
            i1-r1[i3][i2][i1],i2-r2[i3][i2][i1],i3-r3[i3][i2][i1]);
    }});
  }


}
