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
  public void setMaxShift(int smax) {
    _smax = smax;
  }
  public void setPowError(float epow) {
    _epow = epow;
  }

  public float[][] alignTraces(float[][] ts) {
    float[] ws  = new float[1];
    double[] ps = new double[1];
    ws[0] = 1.0f; ps[0] = _epow;
    float[][][] tss = new float[][][]{ts};
    WellLogWarping wlw = new WellLogWarping();
    wlw.setPowError(ps);
    wlw.setMaxShift(_smax);
    float[][] sx = wlw.findShifts(ws,tss);
    float[][] gs = wlw.applyShiftsX(ts,sx);
    return gs;
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


  public float[][][] getReferImageX(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] g = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      g[i3][i2] = f[98][175];
    }}
    return g;
  }


  private int _smax=100;
  private float _epow = 1.0f;

}
