package slt;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

import ad.*;

/**
 * Compute salt likelihoods. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.10
 */

public class SaltScanner {

  public SaltScanner (float sigma1, float sigma2) {
    _h1 = round(sigma1);
    _h2 = round(sigma2);
    //_h2 = round(sigma2*sigma2/2f);
    setScales(_h1);
  }

  public float[][][] scan(float[][] fx, float[][] u1, float[][] u2) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] cx = new float[n2][n1];
    float[][] ax = new float[n2][n1];
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] au = new float[n2][n1];
    float[][] av = new float[n2][n1];
    float[][] g11 = new float[n2][n1];
    float[][] g12 = new float[n2][n1];
    float[][] g22 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply10(fx,g1);
    rgf.apply01(fx,g2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      au[i2][i1] = 0.01f;
      av[i2][i1] = 1.00f;
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      g11[i2][i1] = g1i*g1i;
      g12[i2][i1] = g1i*g2i;
      g22[i2][i1] = g2i*g2i;
    }}
    EigenTensors2 ets = new EigenTensors2(u1,u2,au,av);
    g11=smooth2(ets,g11);
    g12=smooth2(ets,g12);
    g22=smooth2(ets,g22);
    float[][][] g11s = smooth1(g11,u1,u2);
    float[][][] g12s = smooth1(g12,u1,u2);
    float[][][] g22s = smooth1(g22,u1,u2);
    float[][] ca = new float[2][2];
    float[][] aa = new float[2][2];
    float[][] cz = new float[2][2];
    float[][] az = new float[2][2];
    float[] ce = new float[2];
    float[] ae = new float[2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      ca[0][0] = g11s[0][i2][i1];
      ca[0][1] = g12s[0][i2][i1];
      ca[1][0] = g12s[0][i2][i1];
      ca[1][1] = g22s[0][i2][i1];
      aa[0][0] = g11s[1][i2][i1];
      aa[0][1] = g12s[1][i2][i1];
      aa[1][0] = g12s[1][i2][i1];
      aa[1][1] = g22s[1][i2][i1];
      Eigen.solveSymmetric22(ca,cz,ce);
      Eigen.solveSymmetric22(aa,az,ae);
      div(ce,ce[0],ce);
      div(ae,ae[0],ae);
      float ceui = ce[0];
      float cevi = ce[1];
      if (cevi<0.0f) cevi = 0.0f;
      if (ceui<cevi) ceui = cevi;
      float aeui = ae[0];
      float aevi = ae[1];
      if (aevi<0.0f) aevi = 0.0f;
      if (aeui<aevi) aeui = aevi;
      float cei = (ceui-cevi)/ceui;
      float aei = (aeui-aevi)/aeui;
      cx[i2][i1] = cei;//abs(cei-aei);
      ax[i2][i1] = aei;//abs(cei-aei);
    }}
    return new float[][][]{cx,ax};
  }

  public float[][][] scanX(
    EigenTensors3 ets, float[][][] fx, 
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;

    float[][][] g11c = new float[n3][n2][n1];
    float[][][] g12c = new float[n3][n2][n1];
    float[][][] g13c = new float[n3][n2][n1];
    float[][][] g22c = new float[n3][n2][n1];
    float[][][] g23c = new float[n3][n2][n1];
    float[][][] g33c = new float[n3][n2][n1];

    float[][][] g11s = new float[n3][n2][n1];
    float[][][] g12s = new float[n3][n2][n1];
    float[][][] g13s = new float[n3][n2][n1];
    float[][][] g22s = new float[n3][n2][n1];
    float[][][] g23s = new float[n3][n2][n1];
    float[][][] g33s = new float[n3][n2][n1];


    computeGradientProducts(fx,g11c,g12c,g13c,g22c,g23c,g33c);
    trace("structure tensors done...");

    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,_h2,g11c,g11s);
    trace("1st smooth parallel to structures done...");
    lsf.apply(ets,_h2,g12c,g12s);
    trace("2nd smooth parallel to structures done...");
    lsf.apply(ets,_h2,g13c,g13s);
    trace("3rd smooth parallel to structures done...");
    lsf.apply(ets,_h2,g22c,g22s);
    trace("4th smooth parallel to structures done...");
    lsf.apply(ets,_h2,g23c,g23s);
    trace("5th smooth parallel to structures done...");
    lsf.apply(ets,_h2,g33c,g33s);
    trace("6th smooth parallel to structures done...");

    smooth1X(g11s,u1,u2,u3,g11c);
    trace("1st smooth normal to structures done...");
    smooth1X(g12s,u1,u2,u3,g12c);
    trace("2nd smooth normal to structures done...");
    smooth1X(g13s,u1,u2,u3,g13c);
    trace("3rd smooth normal to structures done...");
    smooth1X(g22s,u1,u2,u3,g22c);
    trace("4th smooth normal to structures done...");
    smooth1X(g23s,u1,u2,u3,g23c);
    trace("5th smooth normal to structures done...");
    smooth1X(g33s,u1,u2,u3,g33c);
    trace("6th smooth normal to structures done...");

    solveEigenproblemsX(g11c,g12c,g13c,g22c,g23c,g33c,g11s);
    trace("planarities done...");

    return g11s;
  }


  public float[][][][] scan(
    EigenTensors3 ets, float[][][] fx, 
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;

    float[][][] g11c = new float[n3][n2][n1];
    float[][][] g12c = new float[n3][n2][n1];
    float[][][] g13c = new float[n3][n2][n1];
    float[][][] g22c = new float[n3][n2][n1];
    float[][][] g23c = new float[n3][n2][n1];
    float[][][] g33c = new float[n3][n2][n1];

    float[][][] g11a = new float[n3][n2][n1];
    float[][][] g12a = new float[n3][n2][n1];
    float[][][] g13a = new float[n3][n2][n1];
    float[][][] g22a = new float[n3][n2][n1];
    float[][][] g23a = new float[n3][n2][n1];
    float[][][] g33a = new float[n3][n2][n1];

    float[][][] g11s = new float[n3][n2][n1];
    float[][][] g12s = new float[n3][n2][n1];
    float[][][] g13s = new float[n3][n2][n1];
    float[][][] g22s = new float[n3][n2][n1];
    float[][][] g23s = new float[n3][n2][n1];
    float[][][] g33s = new float[n3][n2][n1];


    computeGradientProducts(fx,g11c,g12c,g13c,g22c,g23c,g33c);
    trace("structure tensors done...");

    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,_h2,g11c,g11s);
    trace("1st smooth parallel to structures done...");
    lsf.apply(ets,_h2,g12c,g12s);
    trace("2nd smooth parallel to structures done...");
    lsf.apply(ets,_h2,g13c,g13s);
    trace("3rd smooth parallel to structures done...");
    lsf.apply(ets,_h2,g22c,g22s);
    trace("4th smooth parallel to structures done...");
    lsf.apply(ets,_h2,g23c,g23s);
    trace("5th smooth parallel to structures done...");
    lsf.apply(ets,_h2,g33c,g33s);
    trace("6th smooth parallel to structures done...");
    /*

    float[][][] g11s = smooth2(ets,g11c);
    trace("1st smooth parallel to structures done...");
    float[][][] g12s = smooth2(ets,g12c);
    trace("2nd smooth parallel to structures done...");
    float[][][] g13s = smooth2(ets,g13c);
    trace("3rd smooth parallel to structures done...");
    float[][][] g22s = smooth2(ets,g22c);
    trace("4th smooth parallel to structures done...");
    float[][][] g23s = smooth2(ets,g23c);
    trace("5th smooth parallel to structures done...");
    float[][][] g33s = smooth2(ets,g33c);
    trace("6th smooth parallel to structures done...");
    */

    smooth1(g11s,u1,u2,u3,g11c,g11a);
    trace("1st smooth normal to structures done...");
    smooth1(g12s,u1,u2,u3,g12c,g12a);
    trace("2nd smooth normal to structures done...");
    smooth1(g13s,u1,u2,u3,g13c,g13a);
    trace("3rd smooth normal to structures done...");
    smooth1(g22s,u1,u2,u3,g22c,g22a);
    trace("4th smooth normal to structures done...");
    smooth1(g23s,u1,u2,u3,g23c,g23a);
    trace("5th smooth normal to structures done...");
    smooth1(g33s,u1,u2,u3,g33c,g33a);
    trace("6th smooth normal to structures done...");

    solveEigenproblems(g11c,g12c,g13c,g22c,g23c,g33c,g11s);
    solveEigenproblems(g11a,g12a,g13a,g22a,g23a,g33a,g22s);
    trace("planarities done...");

    return new float[][][][]{g11s,g22s};
  }

  public void smooth1X(
    final float[][][] fx, final float[][][] u1, final float[][][] u2, 
    final float[][][] u3, final float[][][] gx) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] g1 = new float[n3][n2][n1];
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(_h1);
    rgf.apply100(fx,g1);
    rgf.apply010(fx,g2);
    rgf.apply001(fx,g3);
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float g1i = g1[i3][i2][i1];
        float g2i = g2[i3][i2][i1];
        float g3i = g3[i3][i2][i1];
        gx[i3][i2][i1] = g1i*u1i+g2i*u2i+g3i*u3i;
      }}
    }});
  }


  public void smooth1(
    final float[][][] fx, final float[][][] u1, final float[][][] u2, 
    final float[][][] u3, final float[][][] cx, final float[][][] ax) 
  {
    zero(cx);
    zero(ax);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply000(fx,fx);
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        for (int k=1; k<=_h1; k++) {
          float u1k = k*u1i;
          float u2k = k*u2i;
          float u3k = k*u3i;
          float x1p = i1+u1k;
          float x2p = i2+u2k;
          float x3p = i3+u3k;
          float x1m = i1-u1k;
          float x2m = i2-u2k;
          float x3m = i3-u3k;
          float sci = _c1[k-1];
          float fxm = si.interpolate(s1,s2,s3,fx,x1m,x2m,x3m);
          float fxp = si.interpolate(s1,s2,s3,fx,x1p,x2p,x3p);
          cx[i3][i2][i1] += sci*fxm;
          ax[i3][i2][i1] += sci*fxp;
        }
      }}
    }});
  }


  public float[][][] smooth1(
    final float[][] fx, final float[][] u1, final float[][] u2) 
  {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final float[][] cx = new float[n2][n1];
    final float[][] ax = new float[n2][n1];
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n2, new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        for (int k=1; k<=_h1; k++) {
          float u1k = k*u1i;
          float u2k = k*u2i;
          float x1p = i1+u1k;
          float x2p = i2+u2k;
          float x1m = i1-u1k;
          float x2m = i2-u2k;
          float sci = _c1[k-1];
          float fxm = si.interpolate(s1,s2,fx,x1m,x2m);
          float fxp = si.interpolate(s1,s2,fx,x1p,x2p);
          cx[i2][i1] += sci*fxm;
          ax[i2][i1] += sci*fxp;
        }
      }
    }});
    return new float[][][]{cx,ax};
  }

  public float[][] smooth2(EigenTensors2 et, float[][] fx) {
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setParameters(_h2,5,0.5f);
    return fed.apply(et,fx);
  }

  public float[][][] smooth2(EigenTensors3 et, float[][][] fx) {
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setParameters(_h2,5,0.5f);
    return fed.apply(et,fx);
  }


  public float[][] apply(boolean causal, float a, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[n2][n1];
    if (causal) {
      causal1(a,fx,gx);
      causal2(a,gx,gx);
      causal1(a,gx,gx);
    } else {
      anticausal1(a,fx,gx);
      anticausal2(a,gx,gx);
      anticausal1(a,gx,gx);
    }
    return gx;
  }

  private void computeGradientProducts(
    final float[][][] fx,
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33)
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] g1 = new float[n3][n2][n1];
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply100(fx,g1);
    rgf.apply010(fx,g2);
    rgf.apply001(fx,g3);
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] g1i = g1[i3][i2];
          float[] g2i = g2[i3][i2];
          float[] g3i = g3[i3][i2];
          float[] g11i = g11[i3][i2];
          float[] g12i = g12[i3][i2];
          float[] g13i = g13[i3][i2];
          float[] g22i = g22[i3][i2];
          float[] g23i = g23[i3][i2];
          float[] g33i = g33[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float g1ii = g1i[i1];
            float g2ii = g2i[i1];
            float g3ii = g3i[i1];
            g11i[i1] = g1ii*g1ii;
            g22i[i1] = g2ii*g2ii;
            g33i[i1] = g3ii*g3ii;
            g12i[i1] = g1ii*g2ii;
            g13i[i1] = g1ii*g3ii;
            g23i[i1] = g2ii*g3ii;
          }
        }
      }
    });
  }

  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] ep)
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        double[][] a = new double[3][3];
        double[][] z = new double[3][3];
        double[] e = new double[3];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            a[0][0] = g11[i3][i2][i1];
            a[0][1] = g12[i3][i2][i1];
            a[0][2] = g13[i3][i2][i1];
            a[1][0] = g12[i3][i2][i1];
            a[1][1] = g22[i3][i2][i1];
            a[1][2] = g23[i3][i2][i1];
            a[2][0] = g13[i3][i2][i1];
            a[2][1] = g23[i3][i2][i1];
            a[2][2] = g33[i3][i2][i1];
            Eigen.solveSymmetric33(a,z,e);
            float eui = (float)e[0];
            float evi = (float)e[1];
            float ewi = (float)e[2];
            if (ewi<0.0f) ewi = 0.0f;
            if (evi<ewi) evi = ewi;
            if (eui<evi) eui = evi;
            float esi = (eui>0.0f)?1.0f/eui:1.0f;
            ep[i3][i2][i1] = (eui-evi)*esi;
          }
        }
      }
    });
  }

  private void solveEigenproblemsX(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] ep)
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        double[][] a = new double[3][3];
        double[][] z = new double[3][3];
        double[] e = new double[3];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            a[0][0] = g11[i3][i2][i1];
            a[0][1] = g12[i3][i2][i1];
            a[0][2] = g13[i3][i2][i1];
            a[1][0] = g12[i3][i2][i1];
            a[1][1] = g22[i3][i2][i1];
            a[1][2] = g23[i3][i2][i1];
            a[2][0] = g13[i3][i2][i1];
            a[2][1] = g23[i3][i2][i1];
            a[2][2] = g33[i3][i2][i1];
            Eigen.solveSymmetric33(a,z,e);
            float eui = (float)e[0];
            float evi = (float)e[1];
            float ewi = (float)e[2];
            if (ewi<0.0f) ewi = 0.0f;
            if (evi<ewi) evi = ewi;
            if (eui<evi) eui = evi;
            ep[i3][i2][i1] = eui-evi;
          }
        }
      }
    });
  }



  // Vertically causal and anti-causal filters are implemented 
  // by one-side recursive exponential filters.  
  private void causal2(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float b = 1.0f - a;
    for (int i1=0; i1<n1; ++i1){ 
      float yi = y[0][i1] = x[0][i1];
      for (int i2=1; i2<n2; ++i2) 
        y[i2][i1] = yi = a*yi + b*x[i2][i1];
    }
  }

  private void causal1(float a, float[] x, float[] y) {
    int n1 = x.length;
    float b = 1.0f - a;
    float yi = y[0] = x[0];
    for (int i1=1; i1<n1; ++i1) 
      y[i1] = yi = a*yi + b*x[i1];
  }

  private void causal1(final float a, final float[][] x, final float[][] y) {
    final int n2 = x.length;
    Parallel.loop(n2, new Parallel.LoopInt() {
      public void compute(int i2) {
        causal1(a,x[i2],y[i2]);
      }
    });
  }

  private void anticausal1(float a, float[] x, float[] y) {
    int n1 = x.length;
    float b = 1.0f - a;
    float yi = y[n1-1] = x[n1-1];
    for(int i1=n1-2; i1>=0; --i1)
      y[i1] = yi = a*yi + b*x[i1];
  }
  private void anticausal2(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float b = 1.0f - a;
    for(int i1=0; i1<n1; ++i1) {
      float yi = y[n2-1][i1] = x[n2-1][i1];
      for(int i2=n2-2; i2>=0; --i2)
        y[i2][i1] = yi = a*yi + b*x[i2][i1];
    }
  }

  private void anticausal1(final float a, final float[][] x, final float[][] y) {
    final int n2 = x.length;
    Parallel.loop(n2, new Parallel.LoopInt() {
      public void compute(int i2) {
        anticausal1(a,x[i2],y[i2]);
      }
    });
  }

  private void setScales(int h1) {
    _c1 = new float[h1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(h1/2);
    float[] x = new float[20*h1];
    float[] y = new float[20*h1];
    x[10*h1] = 1;
    rgf.apply0(x,y);
    copy(h1,10*h1,y,0,_c1);
    mul(_c1,10f,_c1);
  }

  private static void trace(String s) {
    System.out.println(s);
  }


  private int _h1 = 5;
  private int _h2 = 10;
  private float[] _c1 = new float[_h1];
}

