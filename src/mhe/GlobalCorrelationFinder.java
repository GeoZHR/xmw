/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mhe;

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
  public float[][] getTraceIndexes(
    int d2, int dm, int dc, int n2, float[] k2, float scale) { 
    ArrayList<Integer> ks = new ArrayList<Integer>();
    for (int i2=0; i2<n2; i2+=d2)
      ks.add(i2);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(dm);
    float[] wm = new float[n2*2+1]; wm[n2] = 1f;
    rgf.apply0(wm,wm);
    wm = div(wm,max(wm));
    int ns = ks.size();
    ArrayList<float[]> kc = new ArrayList<float[]>();
    for (int is=0; is<ns; ++is) {
    for (int js=is+1; js<ns; ++js) {
      int p2 = ks.get(is);
      int m2 = ks.get(js);
      int b2 = abs(p2-m2);
      if(b2<=dm)
        kc.add(new float[]{p2,m2,scale*wm[n2+b2]});
    }}
    RecursiveGaussianFilter rgfc = new RecursiveGaussianFilter(dc);
    float[] wc = new float[n2*2+1]; wc[n2] = 1f;
    rgfc.apply0(wc,wc);
    wc = div(wc,max(wc));
    int np = k2.length;
    for (int ip=0; ip<np; ++ip) {
      int k2i = (int)k2[ip];
      for (int i2=0; i2<n2; i2+=d2) {
        int s2 = abs(i2-k2i);
        kc.add(new float[]{k2i,i2,scale*wc[n2+s2]});
      }
    }
    return kc.toArray(new float[0][]);
  }

  /**
   * Finds 3D indexes of traces to be correlated.
   */
  public float[][] getTraceIndexes(
    int d2, int d3, int dm, int dc, 
    float[] k2, float[] k3, int n2, int n3, float scale) { 
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(dm);
    float[][] wm = new float[n3*2+1][n2*2+1]; 
    wm[n3][n2] = 1f;
    rgf.apply00(wm,wm);
    wm = div(wm,max(wm));
    ArrayList<int[]> ks = new ArrayList<int[]>();
    for (int i3=0; i3<n3; i3+=d3) {
    for (int i2=0; i2<n2; i2+=d2) {
      ks.add(new int[]{i2,i3});
    }}
    int ns = ks.size();
    ArrayList<float[]> kc = new ArrayList<float[]>();
    for (int is=0; is<ns; ++is) {
    for (int js=is+1; js<ns; ++js) {
      int p2 = ks.get(is)[0];
      int p3 = ks.get(is)[1];
      int m2 = ks.get(js)[0];
      int m3 = ks.get(js)[1];
      int b2 = abs(p2-m2);
      int b3 = abs(p3-m3);
      float ds = sqrt(b2*b2+b3*b3);
      if(ds<=dm)
        kc.add(new float[]{p2,p3,m2,m3,scale*wm[n3+b3][n2+b2]});
    }}

    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(1);
    RecursiveGaussianFilter rgfc = new RecursiveGaussianFilter(dc);
    float[][] wc = new float[n3*2+1][n2*2+1]; 
    wc[n3  ][n2  ] = 10f;
    rgf1.apply00(wc,wc);
    rgfc.apply00(wc,wc);
    wc = sub(wc,min(wc));
    wc = div(wc,max(wc));
    int np = k2.length;
    for (int ip=0; ip<np; ++ip) {
      int k2i = (int)k2[ip];
      int k3i = (int)k3[ip];
      int b2 = max(k2i-dc,0);
      int e2 = min(k2i+dc,n2-1);
      int b3 = max(k3i-dc,0);
      int e3 = min(k3i+dc,n3-1);
      for (int i3=b3; i3<=e3; i3+=d3) {
      for (int i2=b2; i2<=e2; i2+=d2) {
        int s2 = abs(i2-k2i);
        int s3 = abs(i3-k3i);
        kc.add(new float[]{k2i,k3i,i2,i3,5*scale*wc[n3+s3][n2+s2]});
      }}
    }
    return kc.toArray(new float[0][]);
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
  public float[][] findCorrelations(final float[][] pc, final float[][] fx) {
    final int nc = pc.length;
    final int n1 = fx[0].length;
    final Sampling s1 = new Sampling(n1);
    final DynamicWarping dw = new DynamicWarping(_shiftMin,_shiftMax);
    dw.setStrainMax(_strainMax);
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    final float[][] t = new float[nc][n1];
    Parallel.loop(nc,new Parallel.LoopInt() {
      public void compute(int ic) {
      int p2 = (int)pc[ic][0];
      int m2 = (int)pc[ic][1];
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

  public float[][] findSlopes(final float[][] pc, final float[][] fx) {
    final int nc = pc.length;
    final int n1 = fx[0].length;
    final DynamicWarping dw = new DynamicWarping(_shiftMin,_shiftMax);
    dw.setStrainMax(_strainMax);
    final float[][] u = new float[nc][n1];
    Parallel.loop(nc,new Parallel.LoopInt() {
      public void compute(int ic) {
      int p2 = (int)pc[ic][0];
      int m2 = (int)pc[ic][1];
      float[] fp = fx[p2];
      float[] fm = fx[m2];
      u[ic] = dw.findShifts(fp,fm);
    }});
    return u;
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
