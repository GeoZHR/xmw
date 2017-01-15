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
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mesh.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import pik.*;
import util.*;

/**
 * 2D fast level set method.
 * <p>
 * Based on the works by Yonggang Shi and William Clem Karl, 2008, 
 * A Real-Time Algorithm for the Approximation of Level-Set-Based 
 * Curve Evolution.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class LevelSet2 {

  public LevelSet2(float mu, float lamda, float alpha, 
    int r, float dt, float niter) {
    _mu = mu;
    _lamda = lamda;
    _alpha = alpha;
    _r = r;
    _dt = dt;
    _niter = niter;
  }

  public float[][] getMark() {
    int n2 = _mark.length;
    int n1 = _mark[0].length;
    float[][] mk = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      mk[i2][i1] = _mark[i2][i1];
    }}
    return mk;
  }

  public float[][] initialLevelSet(
    int n1, int n2, int[] c1, int[] c2, int[] r1, int[] r2, float v) {
    float[][] phi = fillfloat(v,n1,n2);
    int np = c1.length;
    for (int ip=0; ip<np; ++ip) {
      int b1 = c1[ip]-r1[ip];
      int e1 = c1[ip]+r1[ip];
      int b2 = c2[ip]-r2[ip];
      int e2 = c2[ip]+r2[ip];
      for (int i2=b2; i2<=e2; ++i2) 
      for (int i1=b1; i1<=e1; ++i1) 
        phi[i2][i1] = -v;
    }
    updateBand(_r,phi);
    //_mark=fillbyte((byte)1,n1,n2);
    return phi;
  }

  public int[][] toGrayIntegers(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int[][] gx = new int[n2][n1];
    fx = sub(fx,min(fx));
    float fmax = 255f/max(fx);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i2][i1] = round(fx[i2][i1]*fmax);
    }}
    return gx;
  }

  public float[][] applyForEdge(float sigma, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int[][] ft = toGrayIntegers(fx);
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] ge = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    //rgf.apply1X(ft,g1);
    //rgf.applyX1(ft,g2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      ge[i2][i1] = 1f/(1f+g1i*g1i+g2i*g2i);
    }}
    return ge;
  }

  public void updateLevelSet(float sigma, float[][] g, float[][] phi) {
    int n2 = phi.length;
    int n1 = phi[0].length;
    for (int iter=0; iter<_niter; ++iter) {
      if(iter%100==0) System.out.println("iter="+iter);
      float[][] g1 = new float[n2][n1];
      float[][] g2 = new float[n2][n1];
      float[][] ph1 = new float[n2][n1];
      float[][] ph2 = new float[n2][n1];
      float[][] dp1 = new float[n2][n1];
      float[][] dp2 = new float[n2][n1];
      float[][] dps = new float[n2][n1];
      float[][] cvs = new float[n2][n1];
      gradient(g,g1,g2);
      gradient(phi,ph1,ph2);
      distance(ph1,ph2,dp1,dp2);
      divergence(dp1,dp2,dps);
      curevature(ph1,ph2,cvs);
      float[][] ph = copy(phi);
      for (int i2=1; i2<n2-1; ++i2) {
        int i2m = i2-1;
        int i2p = i2+1;
        for (int i1=1; i1<n1-1; ++i1) {
          if(_mark[i2][i1]==0){continue;}
          int i1m = i1-1;
          int i1p = i1+1;
          float phii = ph[i2 ][i1 ]; 
          float phpi = ph[i2p][i1 ]; 
          float phmi = ph[i2m][i1 ]; 
          float phip = ph[i2 ][i1p]; 
          float phim = ph[i2 ][i1m]; 
          float distTerm = dps[i2][i1];
          distTerm += (phip+phim+phpi+phmi-4*phii);
          float gi  = g[i2][i1];
          float g1i = g1[i2][i1];
          float g2i = g2[i2][i1];
          float ph1i = ph1[i2][i1];
          float ph2i = ph2[i2][i1];
          float deltaPhi = delta(phii,sigma);
          float areaTerm = deltaPhi*gi;
          float edgeTerm = deltaPhi*(g1i*ph1i+g2i*ph2i+gi*cvs[i2][i1]);
          phi[i2][i1] += _dt*(_mu*distTerm+_lamda*edgeTerm+_alpha*areaTerm);
        }
      }
      updateBand(_r,phi);
    }
  }

  public void updateLevelSetK(float sigma, float[][] g, float[][] phi) {
    int n2 = phi.length;
    int n1 = phi[0].length;
    for (int iter=0; iter<_niter; ++iter) {
      if(iter%100==0) System.out.println("iter="+iter);
      float[][] g1 = new float[n2][n1];
      float[][] g2 = new float[n2][n1];
      float[][] ph1 = new float[n2][n1];
      float[][] ph2 = new float[n2][n1];
      float[][] dp1 = new float[n2][n1];
      float[][] dp2 = new float[n2][n1];
      float[][] dps = new float[n2][n1];
      float[][] cvs = new float[n2][n1];
      gradient(g,g1,g2);
      gradient(phi,ph1,ph2);
      distance(ph1,ph2,dp1,dp2);
      divergence(dp1,dp2,dps);
      curevature(ph1,ph2,cvs);
      float[][] ph = copy(phi);
      for (int i2=1; i2<n2-1; ++i2) {
        int i2m = i2-1;
        int i2p = i2+1;
        for (int i1=1; i1<n1-1; ++i1) {
          if(_mark[i2][i1]==0){continue;}
          int i1m = i1-1;
          int i1p = i1+1;
          float phii = ph[i2 ][i1 ]; 
          float phpi = ph[i2p][i1 ]; 
          float phmi = ph[i2m][i1 ]; 
          float phip = ph[i2 ][i1p]; 
          float phim = ph[i2 ][i1m]; 
          float distTerm = dps[i2][i1];
          distTerm += (phip+phim+phpi+phmi-4*phii);
          float gi  = g[i2][i1];
          float deltaPhi = delta(phii,sigma);
          float p1 = getGaussianDensity(gi,0f,0.5f);
          float p2 = getGaussianDensity(gi,0f,1.0f);
          phi[i2][i1] += _mu*distTerm+(log(p2)-log(p1)+4f*cvs[i2][i1])*deltaPhi;
        }
      }
      updateBand(_r,phi);
    }
  }

  public void updateLevelSetPK(
    final float sigma, final float[][] el, final float[][] gx, 
    final float[][] p2, final float[][] phi) 
  {
    final int n2 = phi.length;
    final int n1 = phi[0].length;
    final float[][] damp = density(0.6f,el,gx);
    final float[][] ddip = density(0.6f,el,p2);
    balanceDensity(new float[][][]{damp,ddip});
    for (int iter=0; iter<_niter; ++iter) {
      updateBand(_r,phi);
      if(iter%100==0) {
        System.out.println("iter="+iter);
      }
      final float[][] ph1 = new float[n2][n1];
      final float[][] ph2 = new float[n2][n1];
      final float[][] phs = new float[n2][n1];
      final float[][] dps = new float[n2][n1];
      final float[][] cvs = new float[n2][n1];
      compute(phi,ph1,ph2,phs,cvs,dps);
      Parallel.loop(1,n2-1,1,new Parallel.LoopInt() { 
      public void compute(int i2) {
        for (int i1=1; i1<n1-1; ++i1) {
          if(_mark[i2][i1]==0){continue;}
          float cvsi = cvs[i2][i1];
          float distTerm = dps[i2][i1]+phs[i2][i1];
          float deltaPhi = delta(phi[i2][i1],sigma)*_alpha;
          float ext = cvsi*_lamda;
          ext += damp[i2][i1];
          ext += ddip[i2][i1];
          phi[i2][i1] += _mu*distTerm+ext*deltaPhi;
        }
      }});
    }
  }

  public float[][][][] refine(int r, float d, float[][] ph, float[][] fx) {
    int n2 = ph.length;
    int n1 = ph[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(8.0);
    rgf1.apply00(ph,ph);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply1X(ph,g1);
    rgf.applyX1(ph,g2);
    ContourFinder cf = new ContourFinder();
    ContourFinder.Contour ct = cf.findContour(0,s1,s2,ph);
    float[][][] ps = ct.getCoords();
    int nc = ps.length;
    float[][][] bs = bandSample(r,d,ps,ph,g1,g2,fx);
    SincInterpolator si = new SincInterpolator();
    for (int ic=0; ic<nc; ++ic) {
      int m2 = bs[ic].length;
      int m1 = bs[ic][0].length;
      sub(bs[ic],min(bs[ic]),bs[ic]);
      div(bs[ic],max(bs[ic]),bs[ic]);
      OptimalPathPicker opp = new OptimalPathPicker(20,1.0f);
      float[][] ft = opp.applyTransform(bs[ic]);
      float[][] wht = opp.applyForWeight(ft);
      float[][] tms1 = zerofloat(m2,m1);
      float[][] tms2 = zerofloat(m2,m1);
      float[] pik1 = opp.forwardPick(r,wht,tms1);
      float[] pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2);
      int np = ps[ic][0].length;
      float[] x1s = ps[ic][0];
      float[] x2s = ps[ic][1];
      for (int ip=0; ip<np; ++ip) {
        float x1i = x1s[ip];
        float x2i = x2s[ip];
        float g1i = si.interpolate(s1,s2,g1,x1i,x2i);
        float g2i = si.interpolate(s1,s2,g2,x1i,x2i);
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        g1i /= gsi;
        g2i /= gsi;
        x1s[ip] = x1i+g1i*(pik2[ip]-r)*d;
        x2s[ip] = x2i+g2i*(pik2[ip]-r)*d;
      }
    }

    /*
    bs = bandSample(r,d,ps,ph,g1,g2,fx);
    for (int ic=0; ic<nc; ++ic) {
      int m2 = bs[ic].length;
      int m1 = bs[ic][0].length;
      sub(bs[ic],min(bs[ic]),bs[ic]);
      div(bs[ic],max(bs[ic]),bs[ic]);
      OptimalPathPicker opp = new OptimalPathPicker(20,10f);
      float[][] ft = opp.applyTransform(bs[ic]);
      float[][] wht = opp.applyForWeight(ft);
      float[][] tms1 = zerofloat(m2,m1);
      float[][] tms2 = zerofloat(m2,m1);
      float[] pik1 = opp.forwardPick(r,wht,tms1);
      float[] pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2);
      int np = ps[ic][0].length;
      float[] x1s = ps[ic][0];
      float[] x2s = ps[ic][1];
      for (int ip=0; ip<np; ++ip) {
        float x1i = x1s[ip];
        float x2i = x2s[ip];
        float g1i = si.interpolate(s1,s2,g1,x1i,x2i);
        float g2i = si.interpolate(s1,s2,g2,x1i,x2i);
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        g1i /= gsi;
        g2i /= gsi;
        x1s[ip] = x1i+g1i*(pik2[ip]-r)*d;
        x2s[ip] = x2i+g2i*(pik2[ip]-r)*d;
      }
    }
    */


    return new float[][][][]{ps,bs};
  }

  public float[][][] bandSample(
    int r, float d, float[][][] ps, float[][] ph, 
    float[][] g1, float[][] g2, float[][] fx) {
    int n2 = ph.length;
    int n1 = ph[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    int nc = ps.length;
    float[][][] fbs = new float[nc][][];
    SincInterpolator si = new SincInterpolator();
    for (int ic=0; ic<nc; ++ic) {
      float[] x1s = ps[ic][0];
      float[] x2s = ps[ic][1];
      int np = x1s.length;
      fbs[ic] = new float[np][2*r+1];
      for (int ip=0; ip<np; ++ip) {
        float x1i = x1s[ip];
        float x2i = x2s[ip];
        float g1i = si.interpolate(s1,s2,g1,x1i,x2i);
        float g2i = si.interpolate(s1,s2,g2,x1i,x2i);
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        g1i /= gsi;
        g2i /= gsi;
        fbs[ic][ip][r] = si.interpolate(s1,s2,fx,x1i,x2i);
        for (int ir=1; ir<=r; ++ir) {
          float y1i = x1i+g1i*ir*d;
          float y2i = x2i+g2i*ir*d;
          float phi = si.interpolate(s1,s2,ph,y1i,y2i);
          float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
          fbs[ic][ip][ir+r]=fxi;
          //if (phi>0f){fbs[ic][ip][ir+r]=fxi;}
          //else {fbs[ic][ip][ir+r]=-10f;}
        }
        for (int ir=-r; ir<=-1; ++ir) {
          float y1i = x1i+g1i*ir*d;
          float y2i = x2i+g2i*ir*d;
          float phi = si.interpolate(s1,s2,ph,y1i,y2i);
          float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
          fbs[ic][ip][ir+r]=fxi;
          //if (phi<0f){fbs[ic][ip][ir+r]=fxi;}
          //else {fbs[ic][ip][ir+r]=-10f;}
        }
      }
    }
    return fbs;
  }

  private void balanceDensity(float[][][] ds) {
    int n3 = ds.length;
    int n2 = ds[0].length;
    int n1 = ds[0][0].length;
    float[] sc = new float[n3];
    float[] sd = new float[n3];
    for (int i3=0; i3<n3; ++i3) {
      float[][] d3 = ds[i3];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if(d3[i2][i1]<0f) {
          sc[i3] += 1f;
          sd[i3] += abs(d3[i2][i1]);
        }
      }}
    }
    div(sd,sc,sd);
    float sdmax = max(sd);
    for (int i3=0; i3<n3; ++i3) {
      float sci = sdmax/sd[i3];
      mul(ds[i3],sci,ds[i3]);
    }

  }

  public float[][] density(
    float f, float[][] fx, float[][] gx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[] p1 = new float[256];
    float[] p2 = new float[256];
    float[][] dp = new float[n2][n1];
    int[][] gg = toGrayIntegers(gx);
    density(f,fx,gg,p1,p2);
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      int ggi = gg[i2][i1];
      float di = log(p2[ggi])-log(p1[ggi]);
      if(Float.isNaN(di)) {continue;}
      if(Float.isInfinite(di)) {continue;}
      dp[i2][i1] = di;
    }}
    return dp;
  }

  public void density(float f, float[][] fx, 
    int[][] gx, float[] p1, float[] p2) {
    float c1 = 0f;
    float c2 = 0f;
    int n2 = fx.length;
    int n1 = fx[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int gxi = gx[i2][i1];
      float fxi = fx[i2][i1];
      if(fxi<f) {
        c1 += 1f;
        p1[gxi] += 1f;
      } else {
        c2 += 1f;
        p2[gxi] += 1f;
      }
    }}
    div(p1,c1,p1);
    div(p2,c2,p2);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2);
    rgf.apply0(p1,p1);
    rgf.apply0(p2,p2);
  }



  public float[][] density(
    float mu1, float mu2, float sig1, float sig2, float[][] fx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] dp = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float d1 = getGaussianDensity(fx[i2][i1],mu1,sig1);
      float d2 = getGaussianDensity(fx[i2][i1],mu2,sig2);
      float di = log(d2)-log(d1);
      if(Float.isNaN(di)) {continue;}
      if(Float.isInfinite(di)) {continue;}
      dp[i2][i1] = di;
    }}
    return dp;
  }

  public float[][] delta(float[][] phi) {
    int n2 = phi.length;
    int n1 = phi[0].length;
    float[][] dp = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      dp[i2][i1] = delta(phi[i2][i1],1.5f);
    }}
    return dp;

  }


  public void updateLevelSetP(
    final float sigma, final float[][] g, final float[][] phi) {
    final int n2 = phi.length;
    final int n1 = phi[0].length;
    final float[][] g1 = new float[n2][n1];
    final float[][] g2 = new float[n2][n1];
    gradient(g,g1,g2);
    for (int iter=0; iter<_niter; ++iter) {
      if(iter%100==0) System.out.println("iter="+iter);
      final float[][] ph1 = new float[n2][n1];
      final float[][] ph2 = new float[n2][n1];
      final float[][] phs = new float[n2][n1];
      final float[][] dps = new float[n2][n1];
      final float[][] cvs = new float[n2][n1];
      compute(phi,ph1,ph2,phs,cvs,dps);
      Parallel.loop(1,n2-1,1,new Parallel.LoopInt() { 
      public void compute(int i2) {
        for (int i1=1; i1<n1-1; ++i1) {
          if(_mark[i2][i1]==0){continue;}
          float distTerm = dps[i2][i1]+phs[i2][i1];
          float gi  = g[i2][i1];
          float g1i = g1[i2][i1];
          float g2i = g2[i2][i1];
          float ph1i = ph1[i2][i1];
          float ph2i = ph2[i2][i1];
          float cvsi = cvs[i2][i1];
          float deltaPhi = delta(phi[i2][i1],sigma);
          float areaTerm = deltaPhi*gi;
          float edgeTerm = deltaPhi*(g1i*ph1i+g2i*ph2i+gi*cvsi);
          phi[i2][i1] += _dt*(_mu*distTerm+_lamda*edgeTerm+_alpha*areaTerm);
        }
      }});
      updateBand(_r,phi);
    }
  }


  private float getGaussianDensity(float x, float mu, float sigma) {
    float dx = x-mu;
    dx = dx*dx;
    sigma = 1f/sigma;
    float pi = (float)(1f/sqrt(2*Math.PI));
    return (sigma*pi)*exp(-0.5f*dx*sigma*sigma);
  }


  private void compute(final float[][] ph, 
    final float[][] ph1, final float[][] ph2, final float[][] phs, 
    final float[][] cvs, final float[][] dps) {
    final int n2 = ph.length;
    final int n1 = ph[0].length;
    final float[][] dp1 = new float[n2][n1];
    final float[][] dp2 = new float[n2][n1];
    Parallel.loop(1,n2-1,3,new Parallel.LoopInt() { // i1 = 0, 3, 6, ...
      public void compute(int i2) {
        computeGSlice2(i2,ph,ph1,ph2,phs,dp1,dp2);
    }});
    Parallel.loop(2,n2-1,3,new Parallel.LoopInt() { // i1 = 1, 4, 7, ...
      public void compute(int i2) {
        computeGSlice2(i2,ph,ph1,ph2,phs,dp1,dp2);
    }});
    Parallel.loop(3,n2-1,3,new Parallel.LoopInt() { // i1 = 2, 5, 8, ...
      public void compute(int i2) {
        computeGSlice2(i2,ph,ph1,ph2,phs,dp1,dp2);
    }});

    Parallel.loop(1,n2-1,3,new Parallel.LoopInt() { // i1 = 0, 3, 6, ...
      public void compute(int i2) {
        computeCSlice2(i2,dp1,dp2,ph1,ph2,dps,cvs);
    }});
    Parallel.loop(2,n2-1,3,new Parallel.LoopInt() { // i1 = 1, 4, 7, ...
      public void compute(int i2) {
        computeCSlice2(i2,dp1,dp2,ph1,ph2,dps,cvs);
    }});
    Parallel.loop(3,n2-1,3,new Parallel.LoopInt() { // i1 = 2, 5, 8, ...
      public void compute(int i2) {
        computeCSlice2(i2,dp1,dp2,ph1,ph2,dps,cvs);
    }});
  }

  private void computeGSlice2(int i2, float[][] ph, 
    float[][] ph1, float[][] ph2, float[][] phs,
    float[][] dp1, float[][] dp2) {
    int n1 = ph[0].length;
    for (int i1=1; i1<n1-1; ++i1) {
      if(_mark[i2][i1]==0){continue;}
      int i1p = i1+1;
      int i1m = i1-1;
      int i2p = i2+1;
      int i2m = i2-1;
      float phip = ph[i2 ][i1p];
      float phim = ph[i2 ][i1m];
      float phpi = ph[i2p][i1 ];
      float phmi = ph[i2m][i1 ];
      float phii = ph[i2 ][i1 ];
      float ph1ii = (phip-phim)*0.5f; 
      float ph2ii = (phpi-phmi)*0.5f; 
      phs[i2][i1] = phip+phim+phmi+phpi-phii*4;
      float psi = sqrt(ph1ii*ph1ii+ph2ii*ph2ii);
      float dsi = distance(psi)-1;
      dp1[i2][i1] = dsi*ph1ii;
      dp2[i2][i1] = dsi*ph2ii;
      psi = 1f/(psi+0.0000000001f);
      ph1[i2][i1] = ph1ii*psi;
      ph2[i2][i1] = ph2ii*psi;
    }
  }

  private void computeCSlice2(int i2, float[][] dp1, float[][] dp2, 
    float[][] ph1, float[][] ph2, float[][] dps, float[][] cvs) { 
    int n1 = dp1[0].length;
    for (int i1=1; i1<n1-1; ++i1) {
      if(_mark[i2][i1]==0){continue;}
      int i2p = i2+1;
      int i2m = i2-1;
      int i1p = i1+1;
      int i1m = i1-1;
      float dp1ip = dp1[i2 ][i1p]; 
      float dp1im = dp1[i2 ][i1m]; 
      float dp2pi = dp2[i2p][i1 ]; 
      float dp2mi = dp2[i2m][i1 ]; 
      float ph1ip = ph1[i2 ][i1p]; 
      float ph1im = ph1[i2 ][i1m]; 
      float ph2pi = ph2[i2p][i1 ]; 
      float ph2mi = ph2[i2m][i1 ]; 
      dps[i2][i1] = 0.5f*(dp1ip-dp1im+dp2pi-dp2mi);
      cvs[i2][i1] = 0.5f*(ph1ip-ph1im+ph2pi-ph2mi);
    }
  }

  public void updateLevelSetX(float sigma, float[][] g, float[][] phi) {
    int n2 = phi.length;
    int n1 = phi[0].length;
    float d1,d2;
    float d11,d12,d22;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    gradient(g,g1,g2);
    for (int iter=0; iter<_niter; ++iter) {
      if(iter%100==0) System.out.println("iter="+iter);
      float[][] ph = copy(phi);
      float[][] dp1 = new float[n2][n1];
      float[][] dp2 = new float[n2][n1];
      distance(ph,dp1,dp2);
      for (int i2=1; i2<n2-1; ++i2) {
        int i2m = i2-1;
        int i2p = i2+1;
        for (int i1=1; i1<n1-1; ++i1) {
          int i1m = i1-1;
          int i1p = i1+1;
          float phii = ph[i2 ][i1 ]; 
          float phpi = ph[i2p][i1 ]; 
          float phmi = ph[i2m][i1 ]; 
          float phip = ph[i2 ][i1p]; 
          float phim = ph[i2 ][i1m]; 
          float phmm = ph[i2m][i1m]; 
          float phpp = ph[i2p][i1p]; 
          float phpm = ph[i2p][i1m]; 
          float phmp = ph[i2m][i1p]; 
          phii *= 2f;
          d11 = phip+phim-phii;
          d22 = phpi+phmi-phii;
          d1  = 0.5f*(phip-phim);
          d2  = 0.5f*(phpi-phmi);
          if (d1*d2<0f) {
            d12 = 0.5f*(phii+phmm+phpp-phmi-phpi-phim-phip);
          } else {
            d12 = 0.5f*(phmi+phpi+phim+phip-phii-phmp-phpm);
          }
          float ds = 1f/(pow(d1*d1+d2*d2,1.5f)+0.00000000001f);
          float cv = (d1*d1*d22-2*d2*d1*d12+d2*d2*d11)*ds;
          float distTerm = 0f;
          distTerm += 0.5f*(dp1[i2 ][i1p]-dp1[i2 ][i1m]);
          distTerm += 0.5f*(dp2[i2p][i1 ]-dp2[i2m][i1 ]);
          distTerm += (d11+d22);
          float deltaPhi = delta(ph[i2][i1],sigma);
          float gi  = g[i2][i1];
          float g1i = g1[i2][i1];
          float g2i = g2[i2][i1];
          float di = 1f/(sqrt(d1*d1+d2*d2)+0.0000000001f);
          d1 *= di;
          d2 *= di;
          float areaTerm = deltaPhi*gi;
          float edgeTerm = deltaPhi*(g1i*d1+g2i*d2+gi*cv);
          phi[i2][i1] += _dt*(_mu*distTerm+_lamda*edgeTerm+_alpha*areaTerm);
        }
      }
      updateBand(_r,phi);
    }
  }

  private int _r;
  private float _mu;
  private float _dt;
  private float _niter;
  private float _lamda;
  private float _alpha;
  private byte _mark[][]=null;

  private void updateBand(int d, float[][] phi) {
    int n2 = phi.length;
    int n1 = phi[0].length;
    byte[][] mkt = new byte[n2][n1];
    if (_mark!=null) {
      mkt = copy(_mark);
      zero(_mark);
    } else {
      _mark = new byte[n2][n1];
    }
    for (int i2=1; i2<n2-1; i2++) {
    for (int i1=1; i1<n1-1; i1++) {
      float phl = phi[i2-1][i1  ];
      float phr = phi[i2+1][i1  ];
      float pha = phi[i2  ][i1-1];
      float phb = phi[i2  ][i1+1];
      if (phl*phr<0||pha*phb<0) {
        int b1 = i1-d; b1=max(b1,0);
        int e1 = i1+d; e1=min(e1,n1-1);
        int b2 = i2-d; b2=max(b2,0);
        int e2 = i2+d; e2=min(e2,n2-1);
        for (int k2=b2; k2<=e2; ++k2) {
        for (int k1=b1; k1<=e1; ++k1) {
          if (mkt[k2][k1]==0) {
            if (phi[k2][k1]<0) phi[k2][k1] = -d;
            if (phi[k2][k1]>0) phi[k2][k1] =  d;
          }
        }}
        for (int k2=b2+1; k2<e2; ++k2) {
        for (int k1=b1+1; k1<e1; ++k1) {
          _mark[k2][k1] = (byte)1;
        }}
      }
    }}

  }

  private float delta(float x, float sigma) {
    double pi = Math.PI;
    if (abs(x)<=sigma) {
      sigma = 1f/sigma;
      return (float)(0.5*sigma*(1f+cos(pi*x*sigma)));
    } else {
      return 0f;
    }
  }

  private void curevature(float[][] g1, float[][] g2, float[][] cv) {
    int n2 = g1.length;
    int n1 = g2[0].length;
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=1; i1<n1-1; ++i1) {
        if(_mark[i2][i1]==0){continue;}
        float g2i = g2[i2][i1]; 
        float g1i = g1[i2][i1]; 
        float gsi = 1f/(sqrt(g1i*g1i+g2i*g2i)+0.0000000001f);
        g1[i2][i1] *= gsi;
        g2[i2][i1] *= gsi;
      }
    }
    divergence(g1,g2,cv);
  }

  private void gradient(float[][] f, float[][] g1, float[][] g2) {
    int n2 = f.length;
    int n1 = f[0].length;
    for (int i2=1; i2<n2-1; ++i2) {
      int i2m = i2-1;
      int i2p = i2+1;
      for (int i1=1; i1<n1-1; ++i1) {
        int i1m = i1-1;
        int i1p = i1+1;
        float fpi = f[i2p][i1 ]; 
        float fmi = f[i2m][i1 ]; 
        float fip = f[i2 ][i1p]; 
        float fim = f[i2 ][i1m]; 
        g1[i2][i1] = 0.5f*(fip-fim);
        g2[i2][i1] = 0.5f*(fpi-fmi);
      }
    }
  }

  private void divergence(float[][] g1, float[][] g2, float[][] cv) {
    int n2 = g1.length;
    int n1 = g2[0].length;
    for (int i2=1; i2<n2-1; ++i2) {
      int i2m = i2-1;
      int i2p = i2+1;
      for (int i1=1; i1<n1-1; ++i1) {
        if(_mark[i2][i1]==0){continue;}
        int i1m = i1-1;
        int i1p = i1+1;
        float g1ip = g1[i2 ][i1p]; 
        float g1im = g1[i2 ][i1m]; 
        float g2pi = g2[i2p][i1 ]; 
        float g2mi = g2[i2m][i1 ]; 
        cv[i2][i1] = 0.5f*(g1ip-g1im+g2pi-g2mi);
      }
    }
  }

  private void distance(
    float[][] g1, float[][] g2, float[][] dp1, float[][] dp2) {
    int n2 = g1.length;
    int n1 = g1[0].length;
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=1; i1<n1-1; ++i1) {
        if(_mark[i2][i1]==0){continue;}
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        float dsi = distance(gsi);
        dp1[i2][i1] = dsi*g1i-g1i;
        dp2[i2][i1] = dsi*g2i-g2i;
      }
    }
  }

  private void distance(
    float[][] g, float[][] dp1, float[][] dp2) {
    int n2 = g.length;
    int n1 = g[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    gradient(g,g1,g2);
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=1; i1<n1-1; ++i1) {
        if(_mark[i2][i1]==0){continue;}
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        float dsi = distance(gsi);
        dp1[i2][i1] = dsi*g1i-g1i;
        dp2[i2][i1] = dsi*g2i-g2i;
      }
    }
  }


  public float distance(float x) {
    float pt = (float)(2*x*Math.PI);
    if (x>0&&x<=1) {
      return sin(pt)/pt;
    } else if(x>1){
      return (x-1)/x;
    } else {
      return 1;
    }
  }

}
