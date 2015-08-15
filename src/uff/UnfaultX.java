/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package uff;

import ipfx.*;
import java.util.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.dsp.LocalDiffusionKernel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * First compute dip slips for samples alongside faults,
 * then extend dip slips away from faults
 * @author Xinming Wu
 * @version 2015.08.14
 */

public class UnfaultX {

  public float[][][][] getDipSlips(
    int n1, int n2, int n3, FaultSkin[] skins, float smark){
    computeUnfaultShifts(n1,n2,n3,skins);
    float[][][] ss = new float[n3][n2][n1];
    float[][][] s1 = fillfloat(smark,n1,n2,n3);
    float[][][] s2 = fillfloat(smark,n1,n2,n3);
    float[][][] s3 = fillfloat(smark,n1,n2,n3);
    for (FaultSkin skin:skins) {
      for (FaultCell fc:skin) {
        int i1i  = fc.getI1();
        int i2m = bound(fc.getM2(),n2);
        int i3m = bound(fc.getM3(),n3);
        int i2p = bound(fc.getP2(),n2);
        int i3p = bound(fc.getP3(),n3);

        float t1 = 0.5f*fc.getT1();
        float t2 = 0.5f*fc.getT2();
        float t3 = 0.5f*fc.getT3();
        // Set or accumulate slip on the plus side.
        if (s1[i3p][i2p][i1i]==smark) {
          s1[i3p][i2p][i1i]  = t1;
          s2[i3p][i2p][i1i]  = t2;
          s3[i3p][i2p][i1i]  = t3;
          ss[i3p][i2p][i1i]  = 1f;
        } else {
          s1[i3p][i2p][i1i] += t1;
          s2[i3p][i2p][i1i] += t2;
          s3[i3p][i2p][i1i] += t3;
          ss[i3p][i2p][i1i] += 1f;
        }
        // Set or accumulate slip on the plus side.
        if (s1[i3m][i2m][i1i]==smark) {
          s1[i3m][i2m][i1i]  = -t1;
          s2[i3m][i2m][i1i]  = -t2;
          s3[i3m][i2m][i1i]  = -t3;
          ss[i3m][i2m][i1i]  = 1.f;
        } else {
          s1[i3m][i2m][i1i] -= t1;
          s2[i3m][i2m][i1i] -= t2;
          s3[i3m][i2m][i1i] -= t3;
          ss[i3m][i2m][i1i] += 1f;
        }
      }
    }
    // Where more than one slip was accumulated, compute the average.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (ss[i3][i2][i1]>1.0f) {
            float si = 1.0f/ss[i3][i2][i1];
            s1[i3][i2][i1] *= si;
            s2[i3][i2][i1] *= si;
            s3[i3][i2][i1] *= si;
          }
        }
      }
    }
    return new float[][][][]{s1,s2,s3};
  }

  /**
   * Interpolates specified dip-slip vectors.
   * @param s array {s1,s2,s3} of dip-slip vectors.
   * @param smark the mark for slips not adjacent to a fault.
   * @return interpolated dip-slip vectors.
   */
  public static float[][][][] interpolateDipSlips(
      float[][][][] s, float smark) {
    int n3 = s[0].length;
    int n2 = s[0][0].length;
    int n1 = s[0][0][0].length;
    float[][][] sp = new float[n3][n2][n1];
    float[][][] st = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    float[][][][] sq = new float[3][n3][n2][n1];
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(smark,s[0],st,k1,k2,k3);
    clip(0.0f,100.0f,st,st);
    LocalDiffusionKernel.Stencil stencil = LocalDiffusionKernel.Stencil.D21;
    LocalDiffusionKernel ldk = new LocalDiffusionKernel(stencil);
    BlendedGridder3 bg = new BlendedGridder3();
    bg.setBlendingKernel(ldk);
    bg.setSmoothness(0.5);
    for (int is=0; is<3; ++is) {
      float[][][] si = s[is];
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            int j1 = k1[i3][i2][i1];
            int j2 = k2[i3][i2][i1];
            int j3 = k3[i3][i2][i1];
            sp[i3][i2][i1] = si[j3][j2][j1];
          }
        }
      }
      bg.gridBlended(st,sp,sq[is]);
    }
    return sq;
  }

    /**
   * Unfaults an image using interpolated dip-slip vectors.
   * @param s array {s1,s2,s3} of interpolated dip-slip vectors.
   * @param g image to be unfaulted.
   * @return unfaulted image.
   */
  public static float[][][] unfault(float[][][][] s, final float[][][] g) {
    final int n1 = g[0][0].length;
    final int n2 = g[0].length;
    final int n3 = g.length;
    final float[][][] s1 = s[0];
    final float[][][] s2 = s[1];
    final float[][][] s3 = s[2];
    final float[][][] gs = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1 = i1+s1[i3][i2][i1];
          float x2 = i2+s2[i3][i2][i1];
          float x3 = i3+s3[i3][i2][i1];
          gs[i3][i2][i1] = si.interpolate(
              n1,1.0,0.0,
              n2,1.0,0.0,
              n3,1.0,0.0,
              g,x1,x2,x3);
        }
      }
    }});
    return gs;
  }

  private void computeUnfaultShifts(
    final int n1, final int n2, final int n3, final FaultSkin[] skins) {
    final int nk = skins.length;
    Parallel.loop(nk,new Parallel.LoopInt() {
    public void compute(int ik) {
      System.out.println("skin="+ik);
      FaultSkin skin = skins[ik];
      FloatList x1l = new FloatList();
      FloatList x2l = new FloatList();
      FloatList x3l = new FloatList();
      FloatList s1l = new FloatList();
      for (FaultCell cell:skin) {
        x1l.add(cell.getX1());
        x2l.add(cell.getX2());
        x3l.add(cell.getX3());
        s1l.add(cell.getS1());
      }
      float[] x1a = x1l.trim();
      float[] x2a = x2l.trim();
      float[] x3a = x3l.trim();
      float[] s1a = s1l.trim();
      float x2min = 0;
      float x3min = 0;
      float x2max = n2;
      float x3max = n3;
      float x1min = min(x1a);
      float x1max = max(x1a);
      SibsonInterp s1i = new SibsonInterp(s1a,x1a,x2a,x3a);
      s1i.setBounds(x1min,x1max,x2min,x2max,x3min,x3max);
      s1i.setGradientPower(2.0);
      for (FaultCell cell:skin) {
        float s1 = cell.getS1();
        float s2 = cell.getS2();
        float s3 = cell.getS3();
        float p1 = cell.getX1()+s1;
        float p2 = cell.getX2()+s2;
        float p3 = cell.getX3()+s3;
        float m1 = cell.getX1()-s1;
        float m2 = cell.getX2()-s2;
        float m3 = cell.getX3()-s3;

        if(p1>x1max){p1=x1max;}
        if(p2>x2max){p2=x2max;}
        if(p3>x3max){p3=x3max;}
        if(p1<x1min){p1=x1min;}
        if(p2<x2min){p2=x2min;}
        if(p3<x3min){p3=x3min;}
        if(m1>x1max){m1=x1max;}
        if(m2>x2max){m2=x2max;}
        if(m3>x3max){m3=x3max;}
        if(m1<x1min){m1=x1min;}
        if(m2<x2min){m2=x2min;}
        if(m3<x3min){m3=x3min;}

        float dp = s1i.interpolate(p1,p2,p3)-s1;
        float dm = s1-s1i.interpolate(m1,m2,m3);
        s1 -= (dp+dm)*0.5f;
        cell.setUnfaultShifts(new float[]{s1,s2,s3});
      }
    }});
  }

  private int bound(int i, int n) {
    if(i<0) {i=0;}
    if(i>=n){i=n-1;}
    return i;
  }

  private class FloatList {
    public int n = 0;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public void add(float[] f) {
      int m = f.length;
      for (int i=0; i<m; ++i)
        add(f[i]);
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

}

