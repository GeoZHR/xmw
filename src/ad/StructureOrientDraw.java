/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ad;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import ipfx.*;
import ipfx.FaultCell;
import static ipfx.FaultGeometry.*;

/**
 * Structure oriented line drawing. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.04
 */
public class StructureOrientDraw {

  public void setSmoothings(float sigmac, float sigmam) {
    _sigmac = sigmac;
    _sigmam = sigmam;
  }

  public void setThreshold(float tau) {
    _tau = tau;
  }

  public float[][][] applyDraw(float[][] fx, float[][] u1, float[][] u2) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] hx = new float[n2][n1];
    float[][] au = fillfloat(1.0f,n1,n2);
    float[][] av = fillfloat(0.0001f,n1,n2);
    EigenTensors2 et1 = new EigenTensors2(u1,u2,au,av);
    EigenTensors2 et2 = new EigenTensors2(u1,u2,av,au);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(et1,_sigmac,fx,g1);
    lsf.apply(et1,_sigmac*1.6f,fx,g2);
    float[][] gx = sub(g1,g2);
    lsf.apply(et2,_sigmam,gx,hx);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float hxi = hx[i2][i1];
      float hxt = 1f+tanh(hxi);
      if (hxi<0f && hxt<_tau)
        hx[i2][i1] = 0f;
      else
        hx[i2][i1] = 1f;
    }}
    return new float[][][]{hx,u1,u2};
  }

  public float[][][] applyDraw(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] hx = new float[n2][n1];
    LocalOrientFilter lof1 = new LocalOrientFilter(1.0,1.0);
    LocalOrientFilter lof2 = new LocalOrientFilter(2.0,2.0);
    EigenTensors2 et1 = lof1.applyForTensors(fx);
    EigenTensors2 et2 = lof2.applyForTensors(fx);
    et1.setEigenvalues(1.0f,0.0001f);
    et2.setEigenvalues(0.0001f,1.0f);
    lof2.applyForNormal(fx,u1,u2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(et1,_sigmac,fx,g1);
    lsf.apply(et1,_sigmac*1.6f,fx,g2);
    float[][] gx = sub(g1,g2);
    lsf.apply(et2,_sigmam,gx,hx);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float hxi = hx[i2][i1];
      float hxt = 1f+tanh(hxi);
      if (hxi<0f && hxt<_tau)
        hx[i2][i1] = 0f;
      else
        hx[i2][i1] = 1f;
    }}
    return new float[][][]{hx,u1,u2};
  }

  public float[][] applyForEdge(
    float[][] fx, float[][] u1, float[][] u2) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ed = new float[n2][n1];
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      float x1p1 = i1+u1i;
      float x2p1 = i2+u2i;
      float x1m1 = i1-u1i;
      float x2m1 = i2-u2i;
      float x1p2 = i1+u1i*2;
      float x2p2 = i2+u2i*2;
      float x1m2 = i1-u1i*2;
      float x2m2 = i2-u2i*2;
      float fxp1 = si.interpolate(s1,s2,fx,x1p1,x2p1);
      float fxm1 = si.interpolate(s1,s2,fx,x1m1,x2m1);
      float fxp2 = si.interpolate(s1,s2,fx,x1p2,x2p2);
      float fxm2 = si.interpolate(s1,s2,fx,x1m2,x2m2);
      ed[i2][i1] = 0.5f*(fxp1+fxp2-fxm1-fxm2);
    }}
    return ed;
  }

  public float[][][] randSample(
    int np, int seed, float[][] fx, float[][] u1, float[][] u2) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] rx = fillfloat(0,n1,n2);
    float[][] r1 = fillfloat(1,n1,n2);
    float[][] r2 = fillfloat(0,n1,n2);
    Random r = new Random(seed);
    int [][] mark = zeroint(n1,n2);

    /*
    float[][] rt = randfloat(n1,n2);
    float[][] t1 = randfloat(n1,n2);
    float[][] t2 = randfloat(n1,n2);
    mul(rt,2f,rt);
    sub(1f,rt,rt);
    t1 = sub(t1,0.5f);
    t2 = sub(t2,0.5f);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float t1i = t1[i2][i1];
      float t2i = t2[i2][i1];
      float tsi = sqrt(t1i*t1i+t2i*t2i);
      t1[i2][i1] /= tsi;
      t2[i2][i1] /= tsi;
    }}
    for (int ip=0; ip<n1*n2/100; ++ip) {
      boolean marked = false;
      while (!marked) {
        int i1 = r.nextInt(n1);
        int i2 = r.nextInt(n2);
        boolean ok = true;
        int m = 1;
        for (int j2=max(0,i2-m); j2<min(n2,i2+m+1); ++j2){
        for (int j1=max(0,i1-m); j1<min(n1,i1+m+1); ++j1){
          if (mark[j2][j1]>0) ok = false;
        }}
        if (ok) {
          marked = true;
          mark[i2][i1] = 1;
          float rxi = rt[i2][i1];
          float u1i = t1[i2][i1];
          float u2i = t2[i2][i1];
          rx[i2 ][i1 ] = rxi;
          r1[i2 ][i1 ] = u1i;
          r2[i2 ][i1 ] = u2i;
        }
      }
    }
    */
    mark = zeroint(n1,n2);

    for (int ip=0; ip<np; ++ip) {
      boolean marked = false;
      while (!marked) {
        int i1 = r.nextInt(n1);
        int i2 = r.nextInt(n2);
        boolean ok = true;
        int m = 1;
        for (int j2=max(0,i2-m); j2<min(n2,i2+m+1); ++j2){
        //for (int j1=max(0,i1-m); j1<min(n1,i1+m+1); ++j1){
          if (mark[j2][i1]>0) ok = false;
        }
        //}
        if (ok) {
          marked = true;
          mark[i2][i1] = 1;
          float rxi = 1-fx[i2][i1];
          float u1i =   u1[i2][i1];
          float u2i =   u2[i2][i1];
          int i1m = i1-1; if(i1m<0) i1m=0;
          int i1p = i1+1; if(i1p>=n1) i1p=n1-1;
          int i2m = i2-1; if(i2m<0) i2m=0;
          int i2p = i2+1; if(i2p>=n2) i2p=n2-1;
          rx[i2 ][i1 ] = rxi;
          r1[i2 ][i1 ] = u1i;
          r2[i2 ][i1 ] = u2i;
          /*
          rx[i2m][i1 ] = rxi;
          rx[i2p][i1 ] = rxi;
          rx[i2m][i1m] = rxi;
          rx[i2 ][i1m] = rxi;
          rx[i2p][i1m] = rxi;
          rx[i2 ][i1p] = rxi;
          rx[i2m][i1p] = rxi;
          rx[i2p][i1p] = rxi;

          r1[i2m][i1 ] = u1i;
          r1[i2p][i1 ] = u1i;
          r1[i2m][i1m] = u1i;
          r1[i2 ][i1m] = u1i;
          r1[i2p][i1m] = u1i;
          r1[i2m][i1p] = u1i;
          r1[i2 ][i1p] = u1i;
          r1[i2p][i1p] = u1i;

          r2[i2m][i1 ] = u2i;
          r2[i2p][i1 ] = u2i;
          r2[i2m][i1m] = u2i;
          r2[i2 ][i1m] = u2i;
          r2[i2p][i1m] = u2i;
          r2[i2m][i1p] = u2i;
          r2[i2 ][i1p] = u2i;
          r2[i2p][i1p] = u2i;
          */

        }
      }
    }
    return new float[][][]{rx,r1,r2};
  }

  private float _tau = 0.5f;
  private float _sigmac = 1.0f;
  private float _sigmam = 3.0f;

}
