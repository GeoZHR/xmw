/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hvc;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates local slopes of features in 2D and 3D images.
 * 
 * For a 2D image f(x1,x2), slope p2(x1,x2) is the ratio dx1/dx2 of 
 * linear features nearest the image sample at point (x1,x2). An 
 * estimate of the linearity (a number in [0,1]) of features nearest 
 * that sample may also be computed.
 * 
 * Likewise, for a 3D image f(x1,x2,x3), slopes p2(x1,x2,x3) and
 * p3(x1,x2,x3) are the ratios dx1/dx2 and dx1/dx3, respectively, 
 * of planar features nearest the image sample at point (x1,x2,x3). 
 * An estimate of the planarity (a number in [0,1]) of features 
 * nearest that sample may also be computed.
 *
 * All slopes are measured in unitless samples (dx1) per sample (dx2).
 * Minimum and maximum slopes may be specified, and default min and 
 * max bounds are -100 and 100 samples per sample, respectively.
 * Estimated slopes will be within the default or specified bounds.
 *
 * @author Xinming Wu, University of Texas at Austin.
 * @version 2018.08.03
 */
public class GlobalCorrelationFinder {

  /**
   * Constructs a global correlation finder with shift bounds.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public GlobalCorrelationFinder(int shiftMin, int shiftMax) {
    _shiftMin = shiftMin;
    _shiftMax = shiftMax;
  }

  public void setStrainMax(double strainMax) {
    _strainMax = strainMax;
  }

  /**
   * Finds 2D indexes of traces to be correlated.
   */
  public int[][] getTraceIndexes(int d2, int dm, int n2) { 
    ArrayList<Integer> ks = new ArrayList<Integer>();
    for (int i2=0; i2<n2; i2+=d2)
      ks.add(i2);
    int ns = ks.size();
    ArrayList<int[]> kc = new ArrayList<int[]>();
    for (int is=0; is<ns; ++is) {
    for (int js=is+1; js<ns; ++js) {
      int p2 = ks.get(is);
      int m2 = ks.get(js);
      int b2 = abs(p2-m2);
      if(b2<=dm)
        kc.add(new int[]{p2,m2,0});
    }}
    return kc.toArray(new int[0][]);
  }

  /**
   * Finds 3D indexes of traces to be correlated.
   */
  public int[][] getTraceIndexes(int d2, int d3, int dm, int n2, int n3) { 
    ArrayList<int[]> ks = new ArrayList<int[]>();
    for (int i3=0; i3<n3; i3+=d3) {
    for (int i2=0; i2<n2; i2+=d2) {
      ks.add(new int[]{i2,i3});
    }}
    int ns = ks.size();
    ArrayList<int[]> kc = new ArrayList<int[]>();
    for (int is=0; is<ns; ++is) {
    for (int js=is+1; js<ns; ++js) {
      int p2 = ks.get(is)[0];
      int p3 = ks.get(is)[1];
      int m2 = ks.get(js)[0];
      int m3 = ks.get(js)[1];
      int b2 = p2-m2;
      int b3 = p3-m3;
      float ds = sqrt(b2*b2+b3*b3);
      if(ds<=dm)
        kc.add(new int[]{p2,p3,m2,m3,0});
    }}
    return kc.toArray(new int[0][]);
  }

  public int[][] getTraceIndexes(int d2, int d3, int dm, float wm, float[][] wp) { 
    int n3 = wp.length;
    int n2 = wp[0].length;
    ArrayList<int[]> ks = new ArrayList<int[]>();
    for (int i3=0; i3<n3; i3+=d3) {
    for (int i2=0; i2<n2; i2+=d2) {
      if(wp[i3][i2]>wm)
        ks.add(new int[]{i2,i3});
    }}
    int ns = ks.size();
    ArrayList<int[]> kc = new ArrayList<int[]>();
    for (int is=0; is<ns; ++is) {
    for (int js=is+1; js<ns; ++js) {
      int p2 = ks.get(is)[0];
      int p3 = ks.get(is)[1];
      int m2 = ks.get(js)[0];
      int m3 = ks.get(js)[1];
      int b2 = p2-m2;
      int b3 = p3-m3;
      float ds = sqrt(b2*b2+b3*b3);
      if(ds<=dm)
        kc.add(new int[]{p2,p3,m2,m3,0});
    }}
    return kc.toArray(new int[0][]);
  }


  /**
   * Finds correlations of specified traces in a 2D image.
   * @param dm only close traces (distacen smaller than dm) are correlated.
   * @param ps array[np] of trace indexes.
   * @param fx array[n2][n1] of input 2D image samples.
   */
  public float[][] findCorrelations(final int[][] pc, final float[][] fx) {
    final int nc = pc.length;
    final int n1 = fx[0].length;
    final Sampling s1 = new Sampling(n1);
    final DynamicWarping dw = new DynamicWarping(_shiftMin,_shiftMax);
    dw.setStrainMax(_strainMax);
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    final float[][] t = new float[nc][n1];
    Parallel.loop(nc,new Parallel.LoopInt() {
      public void compute(int ic) {
      int p2 = pc[ic][0];
      int m2 = pc[ic][1];
      float[] fp = fx[p2];
      float[] fm = fx[m2];
      float[] u = dw.findShifts(fp,fm);
      for (int i1=0; i1<n1; ++i1)
        u[i1] += i1;
      pc[ic][2] = round(max(u[0],0));
      ii.invert(u,t[ic]);
    }});
    return t;
  }

  /**
   * Finds correlations of specified traces in a 3D image.
   * @param dm only close (distance smaller than dm) traces are correlated.
   * @param ps array[nc][5] of trace indexes.
   * @param fx array[n3][n2][n1] of input 3D image samples.
   */
  public float[][] findCorrelations(int[][] pc, float[][][] fx) {
    int nc = pc.length;
    int n1 = fx[0][0].length;
    Sampling s1 = new Sampling(n1);
    DynamicWarping dw = new DynamicWarping(_shiftMin,_shiftMax);
    dw.setStrainMax(_strainMax);
    InverseInterpolator ii = new InverseInterpolator(s1,s1);
    final float[][] t = new float[nc][n1];
    Parallel.loop(nc,new Parallel.LoopInt() {
      public void compute(int ic) {
        int p2 = pc[ic][0];
        int p3 = pc[ic][1];
        int m2 = pc[ic][2];
        int m3 = pc[ic][3];
        float[] fp = fx[p3][p2];
        float[] fm = fx[m3][m2];
        float[] u = dw.findShifts(fp,fm);
        for (int i1=0; i1<n1; ++i1)
          u[i1] += i1;
        pc[ic][4] = round(max(u[0],0));
        ii.invert(u,t[ic]);
    }});
    return t;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _shiftMin,_shiftMax; 
  private double _strainMax;
}
