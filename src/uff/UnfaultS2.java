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
public class UnfaultS2 {

  /**
   * Constructs an unfaulter.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   */
  public UnfaultS2(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  public void setIters(int inner) {
    _inner=inner;
  }

  public void setTensors(Tensors2 d) {
    _d = d;
  }

  /**
   * Estimates unfault shift vectors in current coordinates for a 3D image.
   * @param sp screen points on faults.
   * @param wp weights, zeros on faults, ones elsewhere.
   * @return array of shifts {r1,r2,r3}.
   */
  public float[][][] findShifts(
    float[][][] sp, float[][] wp) {
    int n2 = wp.length;
    int n1 = wp[0].length;
    float[][][] b = new float[2][n2][n1];
    float[][][] r = new float[2][n2][n1];
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    Smoother2 s2 = new Smoother2(_sigma1,_sigma2,wp);
    CgSolver cg = new CgSolver(_small,_inner);
    A2 ma = new A2(s2,_d,sp,wp);
    vb.zero();
    makeRhs(sp,b);
    s2.applyTranspose(b);
    cg.solve(ma,vb,vr);
    s2.applyOriginal(r);
    //cleanShifts(wp,r);
    return r;
    //return convertShifts(40,r);
  }

  public float[][][] convertShifts(int niter, float[][][] r) {
    float[][][] u = copy(r);
    for (int iter=0; iter<niter; ++iter) {
      float[][][] t = copy(u);
      for (int i=0; i<2; ++i)
        nearestInterp(t,r[i],u[i]);
    }
    return u;
  }

  private void nearestInterp(float[][][] r, float[][] ri, float[][] ui) {
    int n2 = ri.length;
    int n1 = ri[0].length;
    float[][] r1 = r[0], r2 = r[1];
    for (int i2=0; i2<n2; ++i2) { 
    for (int i1=0; i1<n1; ++i1) { 
      float x1 = i1+r1[i2][i1];
      float x2 = i2+r2[i2][i1];
      int k1 = round(x1);
      int k2 = round(x2);
      if(k1<0){k1=0;}if(k1>=n1){k1=n1-1;}
      if(k2<0){k2=0;}if(k2>=n2){k2=n2-1;}
      ui[i2][i1] = ri[k2][k1];
    }}
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
    float[][][] r, float[][] f, float[][] g)
  {
    final int n2 = f.length;
    final int n1 = f[0].length;
    final float[][] ff = f;
    final float[][] gf = g;
    final float[][] r1 = r[0], r2 = r[1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1)
        gf[i2][i1] = si.interpolate(
          n1,1.0,0.0,n2,1.0,0.0,ff,i1+r1[i2][i1],i2+r2[i2][i1]);
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
  private Tensors2 _d;

  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _small = 0.010f; // stop CG iterations if residuals are small
  private int _inner = 100; // maximum number of inner CG iterations



  private static class A2 implements CgSolver.A {
    A2(Smoother2 smoother, Tensors2 et, 
       float[][][] sp, float[][] wp) 
    {
      _et = et;
      _sp = sp;
      _wp = wp;
      _sc = 0.001f;
      _smoother = smoother;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      VecArrayFloat3 v3z = v3x.clone();
      v3y.zero();
      //float[][][][] x = v4x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = v3z.getArray();
      _smoother.applyOriginal(z);
      applyLhs(_et,_wp,z,y);
      if(_sp!=null) {
        screenLhs(_sp[0],_sp[1],_sp[3][0],z,y);
      }
      //add(-_sc,z,y);
      _smoother.applyTranspose(y);
      //add(_sc,x,y);
    }

    private float _sc;
    private Tensors2 _et = null;
    private float[][] _wp=null;
    private float[][][] _sp=null;
    private Smoother2 _smoother;
  }

  private static void add(float sc, float[][][][] x, float[][][][] y) {
    int n4 = x.length;
    int n3 = x[0].length;
    int n2 = x[0][0].length;
    int n1 = x[0][0][0].length;
    for (int i4=0; i4<n4; ++i4) {
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      y[i4][i3][i2][i1] += sc*x[i4][i3][i2][i1];
    }}}}
  }


  private static void makeRhs(float[][][] sp, final float[][][] y) { 
    screenRhs(sp[0],sp[1],sp[2],sp[3][0],y);
  }

  private static void screenRhs(
    float[][] cp, float[][] cm, float[][] dr, float[] fl, float[][][] y) {
    int nc = cp[0].length;
    for (int ic=0; ic<nc; ++ic) {
      float dr1 = dr[0][ic];
      float dr2 = dr[1][ic];
      int i1p = (int)cp[0][ic];
      int i2p = (int)cp[1][ic];
      int i1m = (int)cm[0][ic];
      int i2m = (int)cm[1][ic];
      float fls = fl[ic];//*fl[ic];
      float dx1 = dr1*fls;
      float dx2 = dr2*fls;
      y[0][i2p][i1p] += dx1;
      y[1][i2p][i1p] += dx2;
      y[0][i2m][i1m] -= dx1;
      y[1][i2m][i1m] -= dx2;
    }
  }

  private static void screenLhs(
    float[][] cp, float[][] cm, float[] fl, 
    float[][][] x, float[][][] y) 
  {
    int nc = cp[0].length;
    for (int ic=0; ic<nc; ++ic) {
      int i1p = (int)cp[0][ic];
      int i2p = (int)cp[1][ic];
      int i1m = (int)cm[0][ic];
      int i2m = (int)cm[1][ic];
      float fls = fl[ic];//*fl[ic];

      float dx1 = 0.0f;
      float dx2 = 0.0f;

      dx1 += x[0][i2p][i1p];
      dx2 += x[1][i2p][i1p];

      dx1 -= x[0][i2m][i1m];
      dx2 -= x[1][i2m][i1m];

      dx1 *= fls;
      dx2 *= fls;

      y[0][i2m][i1m] -= dx1;
      y[1][i2m][i1m] -= dx2;

      y[0][i2p][i1p] += dx1;
      y[1][i2p][i1p] += dx2;
    }
  }

  // 3D LHS
  private static void applyLhs(
    Tensors2 d, float[][] wp, float[][][] x, float[][][] y)
  {
    int n2 = y[0].length;
    int n1 = y[0][0].length;
    float[] di = fillfloat(1.0f,3);
    float[][] x1 = x[0]; float[][] x2 = x[1];
    float[][] y1 = y[0]; float[][] y2 = y[1];
    for (int i2=1; i2<n2; ++i2) {
      float[] x10 = x1[i2  ];
      float[] x11 = x1[i2-1];
      float[] x20 = x2[i2  ];
      float[] x21 = x2[i2-1];
      float[] y10 = y1[i2  ];
      float[] y11 = y1[i2-1];
      float[] y20 = y2[i2  ];
      float[] y21 = y2[i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        if(d!=null){d.getTensor(i1,i2,di);}
        float wpi = (wp!=null)?wp[i2][i1]:1.0f;
        float wps = wpi*wpi;
        float d11 = di[0];
        float d12 = di[1];
        float d22 = di[2];
        float x1a = 0.0f;
        float x1b = 0.0f;

        float x2a = 0.0f;
        float x2b = 0.0f;

        x1a += x10[i1 ];
        x1b -= x10[i1m];
        x1b += x11[i1 ];
        x1a -= x11[i1m];

        x2a += x20[i1 ];
        x2b -= x20[i1m];
        x2b += x21[i1 ];
        x2a -= x21[i1m];

        float z11 = 0.5f*(x1a+x1b)*wps;
        float z12 = 0.5f*(x1a-x1b)*wps;

        float z21 = 0.5f*(x2a+x2b)*wps;
        float z22 = 0.5f*(x2a-x2b)*wps;

        float w11 = d11*z11+d12*z12;
        float w12 = d12*z11+d22*z12;

        float w21 = d11*z21+d12*z22;
        float w22 = d12*z21+d22*z22;

        float y1a = 0.5f*(w11+w12);
        float y1b = 0.5f*(w11-w12);

        float y2a = 0.5f*(w21+w22);
        float y2b = 0.5f*(w21-w22);

        y10[i1 ] += y1a;
        y10[i1m] -= y1b;
        y11[i1 ] += y1b;
        y11[i1m] -= y1a;

        y20[i1 ] += y2a;
        y20[i1m] -= y2b;
        y21[i1 ] += y2b;
        y21[i1m] -= y2a;
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
