/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.List;
import java.util.Iterator;
import java.util.LinkedList;
import edu.mines.jtk.mesh.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Run multiple FastLevelSet2 simultanously in multiple threads.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class FastLevelSets2 {
  public FastLevelSets2(int n1, int n2, int[] c1, int[] c2, int[] r) {
    _c1 = c1;
    _c2 = c2;
    _r = r;
  }

  /**
   * Sets the number of update iterations.
   * @param outerIters  total outer iterations.
   * @param speedIters  speed evolution iterations.
   * @param smoothIters smooth iterations.
   */
  public void setIterations(int outerIters, int speedIters, int smoothIters) {
    _outerIters = outerIters;
    _speedIters = speedIters;
    _smoothIters = smoothIters;
  }

  static public float[][] downSample(int d1, int d2, float[][] fx) {
    int m1 = 0;
    int m2 = 0;
    int n2 = fx.length;
    int n1 = fx[0].length;
    for (int i1=0; i1<n1; i1+=d1) m1++;
    for (int i2=0; i2<n2; i2+=d2) m2++;
    float[][] fxs = new float[m2][m1];
    for (int i1=0,k1=0; i1<n1; i1+=d1,k1++) 
    for (int i2=0,k2=0; i2<n2; i2+=d2,k2++)
      fxs[k2][k1] = fx[i2][i1];
    return fxs;
  }

  public float[][] upSample(int d1, int d2, int n1, int n2, float[][] fx) {
    int m1 = fx.length;
    int m2 = fx[0].length;
    int w1 = d1/2;
    int w2 = d2/2;
    float[][] fxu = new float[n2][n1];
    for (int k1=0; k1<m1; ++k1) {
    for (int k2=0; k2<m2; ++k2) {
      int c1 = k1*d1;
      int c2 = k2*d2;
      int b1 = c1-w1;
      int b2 = c2-w2;
      int e1 = c1+w1;
      int e2 = c2+w2;
      b1 = max(b1,0);
      b2 = max(b2,0);
      e1 = min(e1,n1-1);
      e2 = min(e2,n2-1);
      for (int i1=b1;i1<=e1;++i1)
      for (int i2=b2;i2<=e2;++i2)
        fxu[i2][i1] = fx[k2][k1];
    }}
    return fxu;
  }

  public float[][][][] applySegments(
    int gw, float sigma, final float[][] el, 
    final float[][] fx, final float[][] p2, final float[][] ph) {
    final int np = _c1.length;
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    zero(ph);
    add(ph,3,ph);
    final float[][][][] xs = new float[2][np][2][];
    Parallel.loop(np,new Parallel.LoopInt() {
    public void compute(int ip) {
      System.out.println("ip="+ip);
      FastLevelSet2 ls = new FastLevelSet2(n1,n2,_c1[ip],_c2[ip],_r[ip]);
      xs[0][ip] = ls.getLout();
      ls.setIterations(_outerIters,_speedIters,_smoothIters);
      ls.updateLevelSet(gw,sigma,el,fx,p2);
      xs[1][ip] = ls.getLout();
      float[][] phi = ls.getPhi();
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if(phi[i2][i1]<=1) {ph[i2][i1] = phi[i2][i1];}
      }}
    }});
    return xs;
  }

  public float[][][][] applySegments(
    int gw, float sigma, final float[][] fx, final float[][] ph) {
    final int np = _c1.length;
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    zero(ph);
    add(ph,3,ph);
    final float[][][][] xs = new float[2][np][2][];
    Parallel.loop(np,new Parallel.LoopInt() {
    public void compute(int ip) {
      System.out.println("ip="+ip);
      FastLevelSet2 ls = new FastLevelSet2(n1,n2,_c1[ip],_c2[ip],_r[ip]);
      xs[0][ip] = ls.getLout();
      ls.setIterations(_outerIters,_speedIters,_smoothIters);
      ls.updateLevelSet(gw,sigma,fx);
      xs[1][ip] = ls.getLout();
      float[][] phi = ls.getPhi();
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if(phi[i2][i1]<=1) {ph[i2][i1] = phi[i2][i1];}
      }}
    }});
    return xs;
  }

  public float[][][] applySegmentsX(
    int gw, float sigma, final float[][] fx, final float[][] ph) {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final float[][][] xs = new float[2][2][];
    FastLevelSet2 ls = new FastLevelSet2(ph);
    xs[0] = ls.getLout();
    ls.setIterations(_outerIters,_speedIters,_smoothIters);
    ls.updateLevelSet(gw,sigma,fx);
    xs[1] = ls.getLout();
    zero(ph);
    add(ph,3,ph);
    float[][] phi = ls.getPhi();
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(phi[i2][i1]<=1) {ph[i2][i1] = phi[i2][i1];}
    }}
    return xs;
  }



  private int[] _c1;
  private int[] _c2;
  private int[] _r;
  private int _outerIters = 100;
  private int _speedIters = 10;
  private int _smoothIters = 5;

}

