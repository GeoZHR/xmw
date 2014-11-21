/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wsi;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Migrated shot gather image flattening before stacking.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.11.08
 */
public class WarpAndStack {

  public void setForWarp(
    int minlagSc, int maxlagSc, 
    int esmooth, float usmooth, float strainMax1, float strainMax2) 
  {
    _minlagSc = minlagSc;
    _maxlagSc = maxlagSc;
    _esmooth = esmooth;
    _usmooth = usmooth;
    _strainMax1 = strainMax1;
    _strainMax2 = strainMax2;
  }

  public void applyWarp(final float[][][] x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int[] sid = new int[n2];
    smooth3(32.0f,x);
    final float[][] f = nearestShot(sid,x);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      System.out.println("i2="+i2);
      int r3 = sid[i2];
      for (int i3=0; i3<n3; ++i3) {
        int d3 = i3-r3;
        warp(d3,f[i2],x[i3][i2]);
       }
    }});
  }

  public void applyWarp2(final float[][][] x) {
    smooth3(32.0f,x);
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int[] sid = new int[n2];
    final float[][] f = nearestShot(sid,x);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      System.out.println("i3="+i3);
      warp2(f,x[i3]);
    }});
  }
  private void warp2(float[][] f, float[][] x) {
    DynamicWarping dw  = new DynamicWarping(-55,20);
    dw.setErrorSmoothing(_esmooth);
    dw.setShiftSmoothing(_usmooth);
    dw.setStrainMax(_strainMax1,_strainMax2);
    float[][] u = dw.findShifts(f,x);
    float[][] y = dw.applyShifts(u,x);
    copy(y,x);
  }


  public void smooth3(float sigma, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] y = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8.0,1.0,1.0);
    EigenTensors3 et = lof.applyForTensors(x);
    et.setEigenvalues(0.001f,0.001f,1.0f);
    float c = 0.5f*sigma*sigma;
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(et,c,x,y);
    copy(y,x);
  }

  public void smooth(final float sigma, final float[][][]x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      System.out.println("i2="+i2);
      float[][] x2 = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) 
        x2[i3] = x[i3][i2];
      smooth(sigma,x2);
      for (int i3=0; i3<n3; ++i3) 
        x[i3][i2] = x2[i3];
    }});
  }

  public void smooth(float sigma, float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] y = new float[n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8.0,1.0);
    EigenTensors2 et = lof.applyForTensors(x);
    et.setEigenvalues(0.001f,1.0f);
    float c = 0.5f*sigma*sigma;
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(et,c,x,y);
    copy(y,x);
  }
  /*
  public void applyWarp2(final float[][][] x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    final float[][][] p3 = new float[n3][n2][n1];
    final float[][][] ep = new float[n3][n2][n1];
    computeSlope(x,p3,ep);
    final float[][] f = stack(x);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      System.out.println("i3="+i3);
      warp2(f,x[i3]);
    }});
  }
  */


  private void cleanImage(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] x3 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(3.0f);
    rgf.apply001(x,x3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x3i = x3[i3][i2][i1];
          if (abs(x3i)>0.05f) {
            x[i3][i2][i1]  = 0.0f;
          }
        }
      }
    }
  }


  private int maxSlope(float[] p, float[] w) {
    int n1 = p.length;
    float maxP = 0.0f;
    for (int i1=0; i1<n1; ++i1) {
      float wi = w[i1];
      if (wi>0.2f) {
        float pi = abs(p[i1]);
        if (pi>maxP) {maxP=pi;}
      }
    }
    return round(maxP*10.0f);
  }

  private void warp(int offset, float[] f, float[] x) {
    if(offset==0){return;}
    int minlag=-10;
    int maxlag= 10;
    if(offset<0){
      minlag =  offset*_minlagSc; 
      maxlag = -offset*_maxlagSc;
    } else {
      minlag = -offset*_minlagSc; 
      maxlag = offset*_maxlagSc;
    }
    DynamicWarping dw  = new DynamicWarping(minlag,maxlag);
    dw.setStrainMax(_strainMax1);
    dw.setErrorSmoothing(_esmooth);
    dw.setShiftSmoothing(_usmooth);
    float[] u = dw.findShifts(f,x);
    float[] y = dw.applyShifts(u,x);
    copy(y,x);
  }

  /*
  private void warp2(float[][] f, float[][] x) {
    DynamicWarping dw  = new DynamicWarping(_minlag,_maxlag);
    dw.setErrorSmoothing(_esmooth);
    dw.setShiftSmoothing(_usmooth);
    dw.setStrainMax(_strainMax1,_strainMax2);
    float[][] u = dw.findShifts(f,x);
    float[][] y = dw.applyShifts(u,x);
    copy(y,x);
  }

  */

  private void computeSlope(final float[][][] x, final float[][][] p3, float[][][] ep) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] p2 = new float[n3][n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,4.0f,20.0f);
    lsf.findSlopes(x,p2,p3,ep);
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) 
    if(abs(p3[i3][i2][i1])>16.0f) 
       x[i3][i2][i1] = 0.0f;
  }

  private int referenceTrace(float[] xsi, float[][] x2i) {
    int c = 0;
    int n3 = x2i.length;
    float sum = 50000000000.0f;
    for (int i3=0; i3<n3; ++i3) {
      float smi = sum(abs(sub(x2i[i3],xsi)));
      if (smi<sum) {
        c = i3;
        sum = smi;
      }
    }
    return c;
  }


  private int   _maxlagSc =  10;
  private int   _minlagSc = -10;
  private int   _esmooth = 4;
  private float _usmooth = 2.0f;
  private float _strainMax1 = 0.25f;
  private float _strainMax2 = 0.25f;

  private float[][] stack(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][] xs = new float[n2][n1];
    for (int i3=0; i3<n3; ++i3) 
      add(x[i3],xs,xs);
    return div(xs,n3);
  }

  private float[][] nearestShot(int[] sid, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] y = zerofloat(n1,n2,n3);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(3.0f);
    ref.apply3(x,y);
    //float[][][] y = copy(x);
    float[][] r = new float[n2][n1];
    float os = 3.33756f;
    float ox = 3.04800f;
    float dx = 0.02286f;
    float ds = 0.04572f*5.0f*5.0f;
    float[] xs = new float[n3];
    for (int i3=0; i3<n3; ++i3)
      xs[i3] = os+i3*ds;
    for (int i2=0; i2<n2; ++i2) {
      float xi = ox+dx*i2;
      int[] id = new int[1];
      float min = min(abs(sub(xs,xi)),id);
      int i3 = id[0];
      sid[i2] = i3;
      System.out.println("i3="+i3);
      r[i2] = y[i3][i2];
    }
    return r;
  }

  private void stack(float[] g, float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        g[i1] += f[i2][i1];
  }

}
