/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package flc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Refinement of a flattened image using dynamic warping
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.08.17
 */
public class FlattenerR {


  public float[][] alignTraces(int smax, float epow, float[][] ts) {
    WellLogWarpingD wlw = new WellLogWarpingD();
    wlw.setPowError(epow);
    wlw.setMaxShift(smax);
    float[][] sx = wlw.findShifts(ts);
    float[][] gs = wlw.applyShiftsX(ts,sx);
    return gs;
  }

  public float[][][] flattenImage(int smax, 
    float r1min, float r1max, float r2min, 
    float r2max, float r3min, float r3max, 
    float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    DynamicWarpingK dwk = new DynamicWarpingK(8,-smax,smax,s1,s2,s3);
    dwk.setStrainLimits(r1min,r1max,r2min,r2max,r3min,r3max);
    dwk.setSmoothness(4,2,2);
    float[][][] gx = getReferImageX(n2/2,n3/2,fx);
    float[][][] fs = dwk.findShifts(s1,gx,s1,fx);
    return dwk.applyShifts(s1,fx,fs);
  }

  public float[][] flattenTraces(
    int smax, float r1min, float r1max, float[][] fx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    Sampling s1 = new Sampling(n1);
    float[][] gx = new float[n2][n1];
    DynamicWarpingK dwk = new DynamicWarpingK(8,-smax,smax,s1);
    dwk.setStrainLimits(r1min,r1max);
    dwk.setSmoothness(2);
    /*
    float sc = 1.0f/n2;
    float[] fr = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
    for (int i2=0; i2<n2; ++i2) {
      fr[i1] += fx[i2][i1]*sc;
    }}
    */
    float[] fr = getMedianTrace(fx);
    for (int i2=0; i2<n2; ++i2) {
      float[] fs = dwk.findShifts(s1,fr,s1,fx[i2]);
      gx[i2] = dwk.applyShifts(s1,fx[i2],fs);
    }
    return gx;
  }

  public float[][][] sincInterp(
    final Sampling si2, final Sampling si3, 
    final Sampling so2, final Sampling so3, final float[][][] fx) {
    int n3 = so3.getCount();
    final int n2 = so2.getCount();
    final int n1 = fx[0][0].length;
    final double f2 = so2.getFirst();
    final double f3 = so3.getFirst();
    final double d2 = so2.getDelta();
    final double d3 = so3.getDelta();
    final Sampling s1 = new Sampling(n1);
    final float[][][] gx = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      double x1 = i1;
      double x2 = f2+i2*d2;
      double x3 = f3+i3*d3;
      gx[i3][i2][i1] = si.interpolate(s1,si2,si3,fx,x1,x2,x3);
    }}}});
    return gx;
  }

  public float[][] imageToTraces(float[][][] fx) {
    int k = 0;
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][] gx = new float[n2*n3][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      gx[k] = copy(fx[i3][i2]);
      k++;
    }}
    return gx;
  }

  public float[][][] tracesToImage(int n2, int n3, float[][] gx) {
    int k = 0;
    int n1 = gx[0].length;
    float[][][] fx = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      fx[i3][i2] = copy(gx[k]);
      k++;
    }}
    return fx;
  }

  public float[][][] resample(
    Sampling si2, Sampling si3, 
    Sampling so2, Sampling so3, float[][][] fx) 
  {
    int n1 = fx.length;
    int n2 = so2.getCount();
    int n3 = so3.getCount();
    double f2 = si2.getFirst();
    double f3 = si3.getFirst();
    double d2 = si2.getDelta();
    double d3 = si3.getDelta();
    float[][][] gx = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k2 = (int)((so2.getValue(i2)-f2)/d2); 
      int k3 = (int)((so3.getValue(i3)-f3)/d3); 
      gx[i3][i2] = copy(fx[k3][k2]);
    }}
    return gx;
  }

  public float[][] getReferImage(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[] fa = new float[n1];
    float[][] g = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      fa = add(fa,f[i2]);
    }
    fa = div(fa,n2);
    for (int i2=0; i2<n2; ++i2) {
      g[i2] = fa;
    }
    return g;
  }

  public float[][][] getReferImage(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float dm = 0.0f;
    float[] fr = new float[n1];
    float[][][] g  = new float[n3][n2][n1];
    float[][][] g1 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply1XX(f,g1);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float[] g1i = g1[i3][i2];
      float g1s = sum(abs(g1i));
      if(g1s>dm) {fr=f[i3][i2];dm=g1s;}
    }}
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      g[i3][i2] = fr;
    }}

    return g;
  }

  public float[] getMedianTrace(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[] fm = new float[n1];
    float[] f1 = new float[n2];
    MedianFinder mf = new MedianFinder(n2);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        f1[i2] = f[i2][i1];
      }
      fm[i1] = mf.findMedian(f1);
    }
    return fm;
  }


  public float[][] getReferImageM(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[] fa = new float[n1];
    float[] f1 = new float[n2];
    float[][] g = new float[n2][n1];
    MedianFinder mf = new MedianFinder(n2);
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        f1[i2] = f[i2][i1];
      }
      fa[i1] = mf.findMedian(f1);
    }
    for (int i2=0; i2<n2; ++i2) {
      g[i2] = fa;
    }
    return g;
  }

  public float[][][] getReferImageM(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[] fa = new float[n1];
    float[] f1 = new float[n3*n2];
    float[][][] g = new float[n3][n2][n1];
    MedianFinder mf = new MedianFinder(n2*n3);
    for (int i1=0; i1<n1; ++i1) {
      int k = 0;
      for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
        f1[k] = f[i3][i2][i1];
        k++;
      }}
      fa[i1] = mf.findMedian(f1);
    }
    float dm = Float.MAX_VALUE;
    float[] fr = new float[n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float[] fi = f[i3][i2];
      float di = sum(abs(sub(fi,fa)));
      if(di<dm) {fr=copy(fi);dm=di;}
    }}
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      g[i3][i2] = fr;
    }}
    return g;
  }

  public float[][] getReferImageX(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      g[i2] = f[43];
    }
    return g;
  }


  public float[][][] getReferImageX(int k2, int k3, float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] g = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      g[i3][i2] = f[k3][k2];
    }}
    return g;
  }

}
