package uff;

import java.util.concurrent.atomic.AtomicInteger;
import edu.mines.jtk.util.AtomicFloat;

import edu.mines.jtk.util.Threads;
import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.dsp.LocalSmoothingFilter;
import edu.mines.jtk.dsp.RecursiveExponentialFilter;
import static edu.mines.jtk.util.ArrayMath.*;

public class Smoother3 {

  public Smoother3(float sigma1, float sigma2, float sigma3, float[][][] wp) {
    _wp = wp;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
    _sigma3 = sigma3;
  }

  private void weightsForHorizontalSmooth(float[][][] wp) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][] _wh = copy(wp);
    for (int i3=0; i3<n3; ++i3) 
      for (int i2=0; i2<n2; ++i2) 
        for (int i1=0; i1<n1; ++i1) 
          _wh[i3][i2][i1]=(_wh[i3][i2][i1]<0.05f)?0.05f:_wh[i3][i2][i1];
  }

  public void apply(float[][][][] x) {
    int n4 = x.length;
    for (int i4=0; i4<n4; ++i4) {
      smooth1(_sigma1,_wp,x[i4]);
      smooth2(_sigma2,_wp,x[i4]);
      smooth3(_sigma3,_wp,x[i4]);
      smooth3(_sigma3,_wp,x[i4]);
      smooth2(_sigma2,_wp,x[i4]);
      smooth1(_sigma1,_wp,x[i4]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float[][][] _wp;

  private float _sigma1,_sigma2,_sigma3;
  private static final boolean PARALLEL = true;

  private static void subtract(final float[][][] x) {
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    final AtomicFloat af = new AtomicFloat();
    final AtomicInteger ai = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=ai.getAndIncrement(); i3<n3; i3=ai.getAndIncrement())
            af.getAndAdd(sum(x[i3])/(float)n1/(float)n2/(float)n3);
        }
      });
    }
    Threads.startAndJoin(threads);
    sub(x,af.get(),x);
  }

  private static void subtract(float[] e0, float[][][] x) {
    if (PARALLEL) { // TODO: check parallel implementation
      subtractS(e0,x);
    } else {
      subtractS(e0,x);
    }
  }
  private static void subtractS(float[] e0, float[][][] x) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    double d0 = 0.0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          d0 += e0[i1]*x[i3][i2][i1];
        }
      }
    }
    float f0 = (float)d0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          x[i3][i2][i1] -= f0*e0[i1];
        }
      }
    }
  }


    // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    int n2 = x.length;
    int n1 = x[0].length;
    float c = 0.5f*sigma*sigma;
    float[] st = fillfloat(1.0f,n1);
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      if (s!=null) {
        for (int i1=0; i1<n1; ++i1)
          st[i1] = s[i2][i1];
      }
      for (int i1=0; i1<n1; ++i1)
        xt[i1] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  private static void smooth1(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] x3 = x[i3];
      float[][] s3 = (s!=null)?s[i3]:null;
      smooth1(sigma,s3,x3);
    }});
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



  // Smoothing for dimension 1.
  private void smooth1(float sigma, float[][][] x) {
    new RecursiveExponentialFilter(sigma).apply1(x,x);
  }
  /*

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
  private static void smooth2(float sigma, float[][][] s, float[][][] x) {
    if (PARALLEL) {
      smooth2P(sigma,s,x);
    } else {
      smooth2S(sigma,s,x);
    }
  }

  private static void smooth2S(float sigma, float[][][] s, float[][][] x) {
    int n3 = x.length;
    for (int i3=0; i3<n3; ++i3) {
      float[][] s3 = (s!=null)?s[i3]:null;
      float[][] x3 = x[i3];
      smooth2(sigma,s3,x3);
    }
  }
  private static void smooth2P(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    final AtomicInteger ai = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          for (int i3=ai.getAndIncrement(); i3<n3; i3=ai.getAndIncrement()) {
            float[][] s3 = (s!=null)?s[i3]:null;
            float[][] x3 = x[i3];
            smooth2(sigma,s3,x3);
          }
        }
      });
    }
    Threads.startAndJoin(threads);
  }

  // Smoothing for dimension 3.
  private static void smooth3(float sigma, float[][][] s, float[][][] x) {
    if (PARALLEL) {
      smooth3P(sigma,s,x);
    } else {
      smooth3S(sigma,s,x);
    }
  }
  private static void smooth3S(float sigma, float[][][] s, float[][][] x) {
    int n2 = x[0].length;
    int n3 = x.length;
    float[][] s2 = (s!=null)?new float[n3][]:null;
    float[][] x2 = new float[n3][];
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
        if (s!=null)
          s2[i3] = s[i3][i2];
        x2[i3] = x[i3][i2];
      }
      smooth2(sigma,s2,x2);
    }
  }
  private static void smooth3P(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n2 = x[0].length;
    final int n3 = x.length;
    final AtomicInteger ai = new AtomicInteger(0);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
          float[][] s2 = (s!=null)?new float[n3][]:null;
          float[][] x2 = new float[n3][];
          for (int i2=ai.getAndIncrement(); i2<n2; i2=ai.getAndIncrement()) {
            for (int i3=0; i3<n3; ++i3) {
              if (s!=null)
                s2[i3] = s[i3][i2];
              x2[i3] = x[i3][i2];
            }
            smooth2(sigma,s2,x2);
          }
        }
      });
    }
    Threads.startAndJoin(threads);
  }
  */

}
