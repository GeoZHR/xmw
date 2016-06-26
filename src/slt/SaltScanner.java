package slt;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Stopwatch;
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

  /*
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
  */
  public void applyForPhase(float[][] u, float[][] ph){
    int n2 = u.length;
    int n1 = u[0].length; 
    float[] ui = new float[n1];
    HilbertTransformFilter hbt = new HilbertTransformFilter();
    for (int i2=0; i2<n2; i2++){
      hbt.apply(n1,u[i2],ui);
      for (int i1=0; i1<n1; i1++){
        float uri =   u[i2][i1];
        float uii =  ui[i1];
        ph[i2][i1] = -atan2(uii,uri);
      }
    } 
  }

  public float[][] applyForLinear(
    float sigma, EigenTensors2 ets, float[][] fx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    sigma = sigma*sigma*0.5f;

    float[][] g11c = new float[n2][n1];
    float[][] g12c = new float[n2][n1];
    float[][] g22c = new float[n2][n1];

    float[][] g11s = new float[n2][n1];
    float[][] g12s = new float[n2][n1];
    float[][] g22s = new float[n2][n1];


    computeGradientProducts(fx,g11c,g12c,g22c);
    trace("structure tensors done...");
    Stopwatch sw = new Stopwatch();
    sw.start();
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,sigma,g11c,g11s);
    trace("1st smooth parallel to structures done...");
    lsf.apply(ets,sigma,g12c,g12s);
    trace("2nd smooth parallel to structures done...");
    lsf.apply(ets,sigma,g22c,g22s);
    trace("3rd smooth parallel to structures done...");
    /*
    g11s = smooth2(sigma,ets,g11c);
    g12s = smooth2(sigma,ets,g12c);
    g13s = smooth2(sigma,ets,g13c);
    g22s = smooth2(sigma,ets,g22c);
    g23s = smooth2(sigma,ets,g23c);
    g33s = smooth2(sigma,ets,g33c);
    */
    trace("lsf smooth: done in "+sw.time()+" seconds");
    sw.stop();
    return solveEigenproblems(g11s,g12s,g22s);

  }


  public float[][][] applyForPlanar(
    float sigma, EigenTensors3 ets, float[][][] fx) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    sigma = sigma*sigma*0.5f;

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
    Stopwatch sw = new Stopwatch();
    sw.start();
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,sigma,g11c,g11s);
    trace("1st smooth parallel to structures done...");
    lsf.apply(ets,sigma,g12c,g12s);
    trace("2nd smooth parallel to structures done...");
    lsf.apply(ets,sigma,g13c,g13s);
    trace("3rd smooth parallel to structures done...");
    lsf.apply(ets,sigma,g22c,g22s);
    trace("4th smooth parallel to structures done...");
    lsf.apply(ets,sigma,g23c,g23s);
    trace("5th smooth parallel to structures done...");
    lsf.apply(ets,sigma,g33c,g33s);
    trace("6th smooth parallel to structures done...");
    trace("lsf smooth: done in "+sw.time()+" seconds");
    /*
    g11s = smooth2(sigma,ets,g11c);
    g12s = smooth2(sigma,ets,g12c);
    g13s = smooth2(sigma,ets,g13c);
    g22s = smooth2(sigma,ets,g22c);
    g23s = smooth2(sigma,ets,g23c);
    g33s = smooth2(sigma,ets,g33c);
    */
    trace("lsf smooth: done in "+sw.time()+" seconds");
    sw.stop();
    return solveEigenproblems(g11s,g12s,g13s,g22s,g23s,g33s);

  }

  public float[][] saltLikelihood(
    float sigma, float[][] ep, float[][] u1, float[][] u2) 
  {
    final int n2 = ep.length;
    final int n1 = ep[0].length;
    final float[][] g1 = new float[n2][n1];
    final float[][] g2 = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    rgf.apply10(ep,g1);
    rgf.apply01(ep,g2);
    final float[][] sl = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      sl[i2][i1] = abs(g1i*u1i+g2i*u2i); 
    }}
    sub(sl,min(sl),sl);
    div(sl,max(sl),sl);
    return sl;
  }

  public float[][][] saltLikelihood(
    float sigma, float[][][] ep, 
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = ep.length;
    final int n2 = ep[0].length;
    final int n1 = ep[0][0].length;
    final float[][][] g1 = new float[n3][n2][n1];
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    rgf.apply100(ep,g1);
    rgf.apply010(ep,g2);
    rgf.apply001(ep,g3);
    final float[][][] sl = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i3][i2][i1];
        float g2i = g2[i3][i2][i1];
        float g3i = g3[i3][i2][i1];
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        sl[i3][i2][i1] = abs(g1i*u1i+g2i*u2i+g3i*u3i); 
      }}
    }});
    normalize(sl);
    return sl;
  }

  public float[][][] saltLikelihoodX(
    float sigma, float[][][] ep, float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n3 = ep.length;
    int n2 = ep[0].length;
    int n1 = ep[0][0].length;
    float[][][] sl = new float[n3][n2][n1];
    float[][][] ap1 = new float[n3][n2][n1];
    float[][][] cp1 = new float[n3][n2][n1];
    float[][][] ap2 = new float[n3][n2][n1];
    float[][][] cp2 = new float[n3][n2][n1];
    float[][][] ap3 = new float[n3][n2][n1];
    float[][][] cp3 = new float[n3][n2][n1];
    causal1(sigma,ep,cp1);
    causal2(sigma,ep,cp2);
    causal3(sigma,ep,cp3);
    anticausal1(sigma,ep,ap1);
    anticausal2(sigma,ep,ap2);
    anticausal3(sigma,ep,ap3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
        float g1i = cp1[i3][i2][i1]-ap1[i3][i2][i1];
        float g2i = cp2[i3][i2][i1]-ap2[i3][i2][i1];
        float g3i = cp3[i3][i2][i1]-ap3[i3][i2][i1];
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        //sl[i3][i2][i1] = abs(g1i*u1i+g2i*u2i+g3i*u3i); 
        sl[i3][i2][i1] = abs(g1i*g1i+g2i*g2i+g3i*g3i); 
    }}}
    normalize(sl);
    return sl;
  }


  private void normalize(float[][][] x) {
    sub(x,min(x),x);
    div(x,max(x),x);
  }

  public float[][] smooth2(float sigma, EigenTensors2 et, float[][] fx) {
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setCycles(5,0.5f);
    return fed.applyLinearDiffusion(sigma,et,fx);
  }

  public float[][][] smooth2(float sigma, EigenTensors3 et, float[][][] fx) {
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setCycles(5,0.5f);
    return fed.applyLinearDiffusion(sigma,et,fx);
  }


  private void computeGradientProducts(
    final float[][] fx,
    final float[][] g11, final float[][] g12, final float[][] g22)
  {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final float[][] g1 = new float[n2][n1];
    final float[][] g2 = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
    rgf.apply10(fx,g1);
    rgf.apply01(fx,g2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      g11[i2][i1] = g1[i2][i1]*g1[i2][i1];
      g12[i2][i1] = g1[i2][i1]*g2[i2][i1];
      g22[i2][i1] = g2[i2][i1]*g2[i2][i1];
    }}
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
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
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

  private float[][] solveEigenproblems(
    final float[][] g11, final float[][] g12, final float[][] g22)
  {
    final int n2 = g11.length;
    final int n1 = g11[0].length;
    final float[][] ed = new float[n2][n1];
    final float[][] el = new float[n2][n1];
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
        a[0][0] = g11[i2][i1];
        a[0][1] = g12[i2][i1];
        a[1][0] = g12[i2][i1];
        a[1][1] = g22[i2][i1];
        Eigen.solveSymmetric22(a,z,e);
        float u1i = z[0][0];
        float u2i = z[0][1];
        if (u1i<0.0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        float eui = e[0];
        float evi = e[1];
        if (evi<0.0f) evi = 0.0f;
        if (eui<evi) eui = evi;
        float esi = (eui>0.0f)?1.0f/eui:1.0f;
        ed[i2][i1] = (eui-evi)/(eui+evi);///eui;
        el[i2][i1] = 1-evi/eui;//(eui-evi)*esi;
    }}
    return el;//new float[][][]{ed,el};
  }

  private float[][][] solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33)
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    final float[][][] ep = new float[n3][n2][n1];
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
            //ep[i3][i2][i1] = 0.5f*(evi+ewi)/eui;
          }
        }
      }
    });
    return ep;
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

  private void causal(float a, float[] x, float[] y) {
    int n1 = x.length;
    float b = 1.0f - a;
    float yi = y[0] = x[0];
    for (int i1=1; i1<n1; ++i1) 
      y[i1] = yi = a*yi + b*x[i1];
  }

  private void causal(final float a, final float[][] x, final float[][] y) {
    final int n2 = x.length;
    Parallel.loop(n2, new Parallel.LoopInt() {
      public void compute(int i2) {
        causal(a,x[i2],y[i2]);
      }
    });
  }

  private void causal1(
    final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        causal(a,x[i3],y[i3]);
      }
    });
  }

  private void causal2(
    final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        causal2(a,x[i3],y[i3]);
      }
    });
  }

  private void causal3(
    final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[][] x2 = new float[n3][n1];
        float[][] y2 = new float[n3][n1];
        for(int i3=0; i3<n3; ++i3) {
        for(int i1=0; i1<n1; ++i1) {
          x2[i3][i1] = x[i3][i2][i1];
        }}
        causal2(a,x2,y2);
        for(int i3=0; i3<n3; ++i3) {
        for(int i1=0; i1<n1; ++i1) {
          y[i3][i2][i1] = y2[i3][i1];
        }}
      }
    });
  }



  private void anticausal(float a, float[] x, float[] y) {
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

  private void anticausal(final float a, final float[][] x, final float[][] y) {
    final int n2 = x.length;
    Parallel.loop(n2, new Parallel.LoopInt() {
      public void compute(int i2) {
        anticausal(a,x[i2],y[i2]);
      }
    });
  }

  private void anticausal1(
    final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3, new Parallel.LoopInt() {
      public void compute(int i3) {
        anticausal(a,x[i3],y[i3]); 
      }
    });
  }

  private void anticausal2(
    final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3, new Parallel.LoopInt() {
      public void compute(int i3) {
        anticausal2(a,x[i3],y[i3]); 
      }
    });
  }

  private void anticausal3(
    final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[][] x2 = new float[n3][n1];
        float[][] y2 = new float[n3][n1];
        for(int i3=0; i3<n3; ++i3) {
        for(int i1=0; i1<n1; ++i1) {
          x2[i3][i1] = x[i3][i2][i1];
        }}
        anticausal2(a,x2,y2);
        for(int i3=0; i3<n3; ++i3) {
        for(int i1=0; i1<n1; ++i1) {
          y[i3][i2][i1] = y2[i3][i1];
        }}
      }
    });
  }




  private static void trace(String s) {
    System.out.println(s);
  }


}

