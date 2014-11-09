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
    int mlag, int esmooth, float usmooth, float strainMax1, float strainMax2) {
    _mlag = mlag;
    _esmooth = esmooth;
    _usmooth = usmooth;
    _strainMax1 = strainMax1;
    _strainMax2 = strainMax2;
  }

  public void applyWarp(final float[][][] x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    final float[][][] p3 = new float[n3][n2][n1];
    final float[][][] ep = new float[n3][n2][n1];
    computeSlope(x,p3,ep);
    final float[][] f = stack(x);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      System.out.println("i2="+i2);
      for (int i3=0; i3<n3; ++i3) 
        warp(f[i2],x[i3][i2]);
    }});
  }

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

  private void warp(float[] f, float[] x) {
    DynamicWarping dw  = new DynamicWarping(-_mlag,_mlag);
    dw.setStrainMax(_strainMax1);
    dw.setErrorSmoothing(_esmooth);
    dw.setShiftSmoothing(_usmooth);
    float[] u = dw.findShifts(f,x);
    float[] y = dw.applyShifts(u,x);
    copy(y,x);
  }

  private void warp2(float[][] f, float[][] x) {
    DynamicWarping dw  = new DynamicWarping(-_mlag,_mlag);
    dw.setErrorSmoothing(_esmooth);
    dw.setShiftSmoothing(_usmooth);
    dw.setStrainMax(_strainMax1,_strainMax2);
    float[][] u = dw.findShifts(f,x);
    float[][] y = dw.applyShifts(u,x);
    copy(y,x);
  }


  private void computeSlope(final float[][][] x, final float[][][] p3, float[][][] ep) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] p2 = new float[n3][n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,4.0f,10.0f);
    lsf.findSlopes(x,p2,p3,ep);
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) 
    if(abs(p3[i3][i2][i1])>10.0f) 
       x[i3][i2][i1] = 0.0f;
  }



  private int _mlag = 20;
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

  private void stack(float[] g, float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        g[i1] += f[i2][i1];
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


  private int[] referenceTrace(float[][] f) {
    int n2 = f.length;
    float sum = 0.0f;
    int[] c = new int[1];
    for (int i2=0; i2<n2; ++i2) {
      float smi = sum(abs(f[i2]));
      if (smi>sum) {
        c[0] = i2;
        sum = smi;
      }
    }
    System.out.println("ci="+c[0]);
    return c;
  }

}
