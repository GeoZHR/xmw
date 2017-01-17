/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import pik.*;

/**
 * 2D fast level set method.
 * <p>
 * Based on the works by Yonggang Shi and William Clem Karl, 2008, 
 * A Real-Time Algorithm for the Approximation of Level-Set-Based 
 * Curve Evolution.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class SaltPicker2 {

  public float[][] initialBoundary(
    float d, float[] c1, float[] c2, float[] u1, float[][] pa) 
  {
    int nc = c1.length;
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    ArrayList<Float> u1a = new ArrayList<Float>();
    ArrayList<Float> u2a = new ArrayList<Float>();
    for (int ic=1; ic<nc; ++ic) {
      float x1p = c1[ic-1];
      float x1c = c1[ic  ];
      float x2p = c2[ic-1];
      float x2c = c2[ic  ];
      float dx1 = x1c-x1p;
      float dx2 = x2c-x2p;
      float dxc = sqrt(dx1*dx1+dx2*dx2);
      float u1i =  dx2/dxc;
      float u2i = -dx1/dxc;
      if (u1[ic]*u1i<0) {
        u1i = -u1i;
        u2i = -u2i;
      }
      u1a.add(u1i); 
      u2a.add(u2i);
      x1a.add(c1[ic-1]);
      x2a.add(c2[ic-1]);
      for (float di=d; di<dxc; di+=d) {
        float x1i = x1p+dx1*di/dxc;
        float x2i = x2p+dx2*di/dxc;
        x1a.add(x1i); x2a.add(x2i);
        u1a.add(u1i); u2a.add(u2i);
      }
    }
    int np = x1a.size();
    float[][] xus = new float[4][np];
    for (int ip=0; ip<np; ++ip) {
      xus[0][ip] = x1a.get(ip);
      xus[1][ip] = x2a.get(ip);
      xus[2][ip] = u1a.get(ip);
      xus[3][ip] = u2a.get(ip);
    }
    smooth(8,xus);
    return xus;
  }

  private void smooth(float sig, float[][] xu) {
    int np = xu[0].length;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sig);
    ref.setEdges(edges);
    ref.apply1(xu,xu);
    for (int ip=0; ip<np; ++ip) {
      float u1i = xu[2][ip];
      float u2i = xu[3][ip];
      float usa = 1f/sqrt(u1i*u1i+u2i*u2i);
      xu[2][ip] = u1i*usa;
      xu[3][ip] = u2i*usa;
    }
  }

  public float[][] refine(int r, float d, float[][] xu, float[][] fx) {
    float[][] bs = bandSample(r,d,xu,fx);
    int m2 = bs.length;
    int m1 = bs[0].length;
    OptimalPathPicker opp = new OptimalPathPicker(40,2.5f);
    float[][] ft = opp.applyTransform(bs);
    float[][] wht = opp.applyForWeight(ft);
    float[][] tms1 = zerofloat(m2,m1);
    float[][] tms2 = zerofloat(m2,m1);
    float[] pik1 = opp.forwardPick(r,wht,tms1);
    float[] pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2);
    int np = xu[0].length;
    for (int ip=0; ip<np; ++ip) {
      float u1i = xu[2][ip];
      float u2i = xu[3][ip];
      xu[0][ip] += u1i*(pik2[ip]-r)*d;
      xu[1][ip] += u2i*(pik2[ip]-r)*d;
    }
    return bs;
  }


  public float[][] bandSample(
    int r, float d, float[][] xu, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int np = xu[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    float[][] fbs = new float[np][2*r+1];
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    float sig=50f;
    float pi = (float)Math.PI;
    float sigs = sig*sig;
    float gaus = 1f/sqrt(2f*sigs*pi);
    for (int ip=0; ip<np; ++ip) {
      float x1i = xu[0][ip];
      float x2i = xu[1][ip];
      float u1i = xu[2][ip];
      float u2i = xu[3][ip];
      fbs[ip][r] = si.interpolate(s1,s2,fx,x1i,x2i)*gaus;
      for (int ir=1; ir<=r; ++ir) {
        float ird = ir*d;
        float y1i = x1i+u1i*ird;
        float y2i = x2i+u2i*ird;
        float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
        float gui = exp(-0.5f*ird*ird/sigs);
        fbs[ip][ir+r]=fxi*gui*gaus;
      }
      for (int ir=-r; ir<=-1; ++ir) {
        float ird = ir*d;
        float y1i = x1i+u1i*ird;
        float y2i = x2i+u2i*ird;
        float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
        float gui = exp(-0.5f*ird*ird/sigs);
        fbs[ip][ir+r]=fxi*gui*gaus;
      }
    }
    sub(fbs,min(fbs),fbs);
    div(fbs,max(fbs),fbs);
    return fbs;
  }


  public float[][] applyForInsAmp(final float[][] fx) {
    final int n2 = fx.length;
    final int n1 = fx[0].length; 
    final float[][] pa = new float[n2][n1];
    final HilbertTransformFilter hbt = new HilbertTransformFilter();
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
      float[] fi = new float[n1];
      hbt.apply(n1,fx[i2],fi);
      for (int i1=0; i1<n1; i1++){
        float fxi = fi[i1];
        float fxr = fx[i2][i1];
        float pai = sqrt(fxr*fxr+fxi*fxi);
        if(Float.isInfinite(pai)||Float.isNaN(pai)){
          pa[i2][i1] = 0f;
        } else { pa[i2][i1] = pai; }
      }
    }}); 
    return pa;
  }



}
