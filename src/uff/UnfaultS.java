package uff;

import vec.*;
//import ifs.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates shift vectors to undo faulting of images.
 * @author Xinming
 * @version 2015.02.25
 */
public class UnfaultS {

  /**
   * Constructs an unfaulter.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   */
  public UnfaultS(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  public void setIters(int inner) {
    _inner=inner;
  }

  public void setTensors(Tensors3 d) {
    _d = d;
  }

  /**
   * Estimates unfault shift vectors in current coordinates for a 3D image.
   * @param sp screen points on faults.
   * @param wp weights, zeros on faults, ones elsewhere.
   * @return array of shifts {r1,r2,r3}.
   */
  public float[][][][] findShifts(
    float[][][] sp, float[][][] wp) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][][] b = new float[3][n3][n2][n1];
    float[][][][] r = new float[3][n3][n2][n1];
    VecArrayFloat4 vr = new VecArrayFloat4(r);
    VecArrayFloat4 vb = new VecArrayFloat4(b);
    Smoother3 s3 = new Smoother3(_sigma1,_sigma2,_sigma2,wp);
    CgSolver cg = new CgSolver(_small,_inner);
    A3 ma = new A3(s3,_d,sp,wp);
    vb.zero();
    makeRhs(sp,b);
    s3.applyTranspose(b);
    cg.solve(ma,vb,vr);
    s3.applyOriginal(r);
    //cleanShifts(wp,r);
    return r;
    //return convertShifts(40,r);
  }

  public float[][][][] convertShifts(int niter, float[][][][] r) {
    float[][][][] u = scopy(r);
    for (int iter=0; iter<niter; ++iter) {
      float[][][][] t = scopy(u);
      for (int i=0; i<3; ++i)
        nearestInterp(t,r[i],u[i]);
    }
    return u;
  }

  private void nearestInterp(float[][][][] r, float[][][] ri, float[][][] ui) {
    int n3 = ri.length;
    int n2 = ri[0].length;
    int n1 = ri[0][0].length;
    float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    for (int i3=0; i3<n3; ++i3) { 
    for (int i2=0; i2<n2; ++i2) { 
    for (int i1=0; i1<n1; ++i1) { 
      float x1 = i1+r1[i3][i2][i1];
      float x2 = i2+r2[i3][i2][i1];
      float x3 = i3+r3[i3][i2][i1];
      int k1 = round(x1);
      int k2 = round(x2);
      int k3 = round(x3);
      if(k1<0){k1=0;}if(k1>=n1){k1=n1-1;}
      if(k2<0){k2=0;}if(k2>=n2){k2=n2-1;}
      if(k3<0){k3=0;}if(k3>=n3){k3=n3-1;}
      ui[i3][i2][i1] = ri[k3][k2][k1];
    }}}
  }

  private void cleanShifts(float[][][] wp, float[][][][] r) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][] ds = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(0.0f,wp,ds,k1,k2,k3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wpi = wp[i3][i2][i1];
      if(wpi==0.0f) {
        int j1 = k1[i3][i2][i1];
        int j2 = k2[i3][i2][i1];
        int j3 = k3[i3][i2][i1];
        r[0][i3][i2][i1] = r[0][j3][j2][j1];
        r[1][i3][i2][i1] = r[1][j3][j2][j1];
        r[2][i3][i2][i1] = r[2][j3][j2][j1];
      }
    }}}
  }



  /**
   * Compute an unfaulted image
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
            ff,i1+r1[i3][i2][i1],i2+r2[i3][i2][i1],i3+r3[i3][i2][i1]);
    }});
  }

  /**
   * Compute a faulted image.
   * @param r input array {r1,r2,r3} of shifts.
   * @param f input image.
   * @param g output shifted image.
   */
  public void applyShiftsX(
    float[][][][] r, float[][][] f, float[][][] g)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] ff = f;
    final float[][][] gf = g;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          gf[i3][i2][i1] = si.interpolate(
            n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,
            ff,i1-r1[i3][i2][i1],i2-r2[i3][i2][i1],i3-r3[i3][i2][i1]);
    }});
  }

  public void applyShiftsR(
    float[][][][] r, float[][] f, float[][] g)
  {
    int n3 = r[0].length;
    int n2 = r[0][0].length;
    int n1 = r[0][0][0].length;
    float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final SincInterpolator si1 = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    float[] k1 = new float[n2*n3];
    float[] k2 = new float[n2*n3];
    float[] k3 = new float[n2*n3];
    int k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float i1  = f[i3][i2];
      float r1i = si1.interpolate(s1,s2,s3,r1,i1,i2,i3);
      float r2i = si1.interpolate(s1,s2,s3,r2,i1,i2,i3);
      float r3i = si1.interpolate(s1,s2,s3,r3,i1,i2,i3);
      k1[k] = i1+r1i;
      k2[k] = i2+r2i;
      k3[k] = i3+r3i;
      k++;
    }}
    NearestGridder2 ng = new NearestGridder2(k1,k2,k3);
    float[][] fg = ng.grid(s2,s3);
    copy(fg,g);
  }



  ///////////////////////////////////////////////////////////////////////////
  // private
  private Tensors3 _d;

  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _small = 0.010f; // stop CG iterations if residuals are small
  private int _inner = 100; // maximum number of inner CG iterations



  private static class A3 implements CgSolver.A {
    A3(Smoother3 smoother, Tensors3 et, 
       float[][][] sp, float[][][] wp) 
    {
      _et = et;
      _sp = sp;
      _wp = wp;
      _smoother = smoother;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat4 v4x = (VecArrayFloat4)vx;
      VecArrayFloat4 v4y = (VecArrayFloat4)vy;
      VecArrayFloat4 v4z = v4x.clone();
      v4y.zero();
      float[][][][] x = v4z.getArray();
      float[][][][] y = v4y.getArray();
      _smoother.applyOriginal(x);
      applyLhs(_et,_wp,x,y);
      if(_sp!=null) {
        screenLhs(_sp[0],_sp[1],_sp[3][0],x,y);
      }
      _smoother.applyTranspose(y);
    }

    private Tensors3 _et = null;
    private float[][][] _wp=null;
    private float[][][] _sp=null;
    private Smoother3 _smoother;
  }


  private static void makeRhs(float[][][] sp, final float[][][][] y) { 
    screenRhs(sp[0],sp[1],sp[2],sp[3][0],y);
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
      float fls = fl[ic];//*fl[ic];
      float dx1 = dr1*fls;
      float dx2 = dr2*fls;
      float dx3 = dr3*fls;
      y[0][i3p][i2p][i1p] += dx1;
      y[1][i3p][i2p][i1p] += dx2;
      y[2][i3p][i2p][i1p] += dx3;
      y[0][i3m][i2m][i1m] -= dx1;
      y[1][i3m][i2m][i1m] -= dx2;
      y[2][i3m][i2m][i1m] -= dx3;
    }
  }

  private static void applyLhs(
    final Tensors3 d, final float[][][] wp, 
    final float[][][][] x, final float[][][][] y)
  { 
    final int n3 = y[0].length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,wp,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,wp,x,y);
    }});
  }

  private static void screenLhs(
    float[][] cp, float[][] cm, float[] fl, float[][][][] x, float[][][][] y) 
  {
    int nc = cp[0].length;
    for (int ic=0; ic<nc; ++ic) {
      int i1p = (int)cp[0][ic];
      int i2p = (int)cp[1][ic];
      int i3p = (int)cp[2][ic];
      int i1m = (int)cm[0][ic];
      int i2m = (int)cm[1][ic];
      int i3m = (int)cm[2][ic];
      float fls = fl[ic];//*fl[ic];

      float dx1 = 0.0f;
      float dx2 = 0.0f;
      float dx3 = 0.0f;

      dx1 += x[0][i3p][i2p][i1p];
      dx2 += x[1][i3p][i2p][i1p];
      dx3 += x[2][i3p][i2p][i1p];

      dx1 -= x[0][i3m][i2m][i1m];
      dx2 -= x[1][i3m][i2m][i1m];
      dx3 -= x[2][i3m][i2m][i1m];

      dx1 *= fls;
      dx2 *= fls;
      dx3 *= fls;

      y[0][i3m][i2m][i1m] -= dx1;
      y[1][i3m][i2m][i1m] -= dx2;
      y[2][i3m][i2m][i1m] -= dx3;

      y[0][i3p][i2p][i1p] += dx1;
      y[1][i3p][i2p][i1p] += dx2;
      y[2][i3p][i2p][i1p] += dx3;
    }
  }

  // 3D LHS
  private static void applyLhsSlice3(
    int i3, Tensors3 d, float[][][] wp, float[][][][] x, float[][][][] y)
  {
    int n2 = y[0][0].length;
    int n1 = y[0][0][0].length;
    float[] di = fillfloat(1.0f,6);
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
        if(d!=null){d.getTensor(i1,i2,i3,di);}
        float wpi = (wp!=null)?wp[i3][i2][i1]:1.0f;
        float wps = wpi*wpi;
        float d11 = di[0];
        float d12 = di[1];
        float d13 = di[2];
        float d22 = di[3];
        float d23 = di[4];
        float d33 = di[5];
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

        float x11 = 0.25f*(x1a+x1b+x1c+x1d)*wps;
        float x12 = 0.25f*(x1a-x1b+x1c-x1d)*wps;
        float x13 = 0.25f*(x1a+x1b-x1c-x1d)*wps;

        float x21 = 0.25f*(x2a+x2b+x2c+x2d)*wps;
        float x22 = 0.25f*(x2a-x2b+x2c-x2d)*wps;
        float x23 = 0.25f*(x2a+x2b-x2c-x2d)*wps;

        float x31 = 0.25f*(x3a+x3b+x3c+x3d)*wps;
        float x32 = 0.25f*(x3a-x3b+x3c-x3d)*wps;
        float x33 = 0.25f*(x3a+x3b-x3c-x3d)*wps;

        float y11 = d11*x11+d12*x12+d13*x13;
        float y12 = d12*x11+d22*x12+d23*x13;
        float y13 = d13*x11+d23*x12+d33*x13;

        float y21 = d11*x21+d12*x22+d13*x23;
        float y22 = d12*x21+d22*x22+d23*x23;
        float y23 = d13*x21+d23*x22+d33*x23;

        float y31 = d11*x31+d12*x32+d13*x33;
        float y32 = d12*x31+d22*x32+d23*x33;
        float y33 = d13*x31+d23*x32+d33*x33;

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


}
