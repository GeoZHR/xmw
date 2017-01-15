/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

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


public class LevelSet3 {

  public LevelSet3(float mu, float lamda, 
    float alpha, int r, float niter) {
    _mu = mu;
    _lamda = lamda;
    _alpha = alpha;
    _r = r;
    _niter = niter;
  }

  public float[][][] getMark() {
    int n3 = _mark.length;
    int n2 = _mark[0].length;
    int n1 = _mark[0][0].length;
    float[][][] mk = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      mk[i3][i2][i1] = _mark[i3][i2][i1];
    }}}
    return mk;
  }

  public float[][][] initialLevelSet(
    int n1, int n2, int n3, int[] c1, int[] c2, int[] c3, 
    int[] r, float v) {
    float[][][] phi = fillfloat(v,n1,n2,n3);
    int np = c1.length;
    for (int ip=0; ip<np; ++ip) {
      int b1 = c1[ip]-r[ip];
      int e1 = c1[ip]+r[ip];
      int b2 = c2[ip]-r[ip];
      int e2 = c2[ip]+r[ip];
      int b3 = c3[ip]-r[ip];
      int e3 = c3[ip]+r[ip];
      for (int i3=b3; i3<=e3; ++i3) 
      for (int i2=b2; i2<=e2; ++i2) 
      for (int i1=b1; i1<=e1; ++i1) 
        phi[i3][i2][i1] = -v;
    }
    updateBand(_r,phi);
    return phi;
  }

  public int[][][] toGrayIntegers(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    int[][][] gx = new int[n3][n2][n1];
    fx = sub(fx,min(fx));
    float fmax = 255f/max(fx);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i3][i2][i1] = round(fx[i3][i2][i1]*fmax);
    }}}
    return gx;
  }

  public void updateLevelSetPK(
    final float sigma, final float[][][] el, 
    final float[][][][] gs, final float[][][] phi) 
  {
    final int ns = gs.length;
    final int n3 = phi.length;
    final int n2 = phi[0].length;
    final int n1 = phi[0][0].length;
    final float[][][][] ds = new float[ns][n3][n2][n1];
    for (int is=0; is<ns; ++is)
      ds[is] = density(0.6f,el,gs[is]);
    balanceDensity(ds);
    for (int iter=0; iter<_niter; ++iter) {
      updateBand(_r,phi);
      System.out.println("iter="+iter);
      final float[][][] phs = new float[n3][n2][n1];
      final float[][][] dps = new float[n3][n2][n1];
      final float[][][] cvs = new float[n3][n2][n1];
      compute(phi,phs,cvs,dps);
      Parallel.loop(1,n3-1,1,new Parallel.LoopInt() { 
      public void compute(int i3) {
        for (int i2=1; i2<n2-1; ++i2) {
        for (int i1=1; i1<n1-1; ++i1) {
          if(_mark[i3][i2][i1]==0){continue;}
          float cvsi = cvs[i3][i2][i1];
          float distTerm = dps[i3][i2][i1]+phs[i3][i2][i1];
          float deltaPhi = delta(phi[i3][i2][i1],sigma)*_alpha;
          float ext = cvsi*_lamda;
          for (int is=0; is<ns; ++is)
            ext += ds[is][i3][i2][i1];
          phi[i3][i2][i1] += _mu*distTerm+ext*deltaPhi;
        }}
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
          float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
          fbs[ic][ip][ir+r]=fxi;
        }
        for (int ir=-r; ir<=-1; ++ir) {
          float y1i = x1i+g1i*ir*d;
          float y2i = x2i+g2i*ir*d;
          float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
          fbs[ic][ip][ir+r]=fxi;
        }
      }
    }
    return fbs;
  }

  private void balanceDensity(float[][][][] ds) {
    int n4 = ds.length;
    int n3 = ds[0].length;
    int n2 = ds[0][0].length;
    int n1 = ds[0][0][0].length;
    float[] sc = new float[n4];
    float[] sd = new float[n4];
    for (int i4=0; i4<n4; ++i4) {
      float[][][] d4 = ds[i4];
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if(d4[i3][i2][i1]<0f) {
          sc[i4] += 1f;
          sd[i4] += abs(d4[i3][i2][i1]);
        }
      }}}
    }
    div(sd,sc,sd);
    float sdmax = max(sd);
    for (int i4=0; i4<n4; ++i4) {
      float sci = sdmax/sd[i4];
      mul(ds[i4],sci,ds[i4]);
    }

  }

  public float[][][] density(
    float f, float[][][] fx, float[][][] gx) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] p1 = new float[256];
    float[] p2 = new float[256];
    float[][][] dp = new float[n3][n2][n1];
    int[][][] gg = toGrayIntegers(gx);
    density(f,fx,gg,p1,p2);
    for (int i3=1; i3<n3-1; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      int ggi = gg[i3][i2][i1];
      float di = log(p2[ggi])-log(p1[ggi]);
      if(Float.isNaN(di)) {continue;}
      if(Float.isInfinite(di)) {continue;}
      dp[i3][i2][i1] = di;
    }}}
    return dp;
  }

  public void density(float f, float[][][] fx, 
    int[][][] gx, float[] p1, float[] p2) {
    float c1 = 0f;
    float c2 = 0f;
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int gxi = gx[i3][i2][i1];
      float fxi = fx[i3][i2][i1];
      if(fxi<f) {
        c1 += 1f;
        p1[gxi] += 1f;
      } else {
        c2 += 1f;
        p2[gxi] += 1f;
      }
    }}}
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


  private float getGaussianDensity(float x, float mu, float sigma) {
    float dx = x-mu;
    dx = dx*dx;
    sigma = 1f/sigma;
    float pi = (float)(1f/sqrt(2*Math.PI));
    return (sigma*pi)*exp(-0.5f*dx*sigma*sigma);
  }


  private void compute(final float[][][] ph, 
    final float[][][] phs, final float[][][] cvs, final float[][][] dps) {
    final int n3 = ph.length;
    final int n2 = ph[0].length;
    final int n1 = ph[0][0].length;
    final float[][][] dp1 = new float[n3][n2][n1];
    final float[][][] dp2 = new float[n3][n2][n1];
    final float[][][] dp3 = new float[n3][n2][n1];
    final float[][][] ph1 = new float[n3][n2][n1];
    final float[][][] ph2 = new float[n3][n2][n1];
    final float[][][] ph3 = new float[n3][n2][n1];
    Parallel.loop(1,n3-1,3,new Parallel.LoopInt() { // i1 = 0, 3, 6, ...
      public void compute(int i3) {
        computeGSlice3(i3,ph,ph1,ph2,ph3,phs,dp1,dp2,dp3);
    }});
    Parallel.loop(2,n3-1,3,new Parallel.LoopInt() { // i1 = 1, 4, 7, ...
      public void compute(int i3) {
        computeGSlice3(i3,ph,ph1,ph2,ph3,phs,dp1,dp2,dp3);
    }});
    Parallel.loop(3,n3-1,3,new Parallel.LoopInt() { // i1 = 2, 5, 8, ...
      public void compute(int i3) {
        computeGSlice3(i3,ph,ph1,ph2,ph3,phs,dp1,dp2,dp3);
    }});

    Parallel.loop(1,n3-1,3,new Parallel.LoopInt() { // i1 = 0, 3, 6, ...
      public void compute(int i3) {
        computeCSlice3(i3,dp1,dp2,dp3,ph1,ph2,ph3,dps,cvs);
    }});
    Parallel.loop(2,n3-1,3,new Parallel.LoopInt() { // i1 = 1, 4, 7, ...
      public void compute(int i3) {
        computeCSlice3(i3,dp1,dp2,dp3,ph1,ph2,ph3,dps,cvs);
    }});
    Parallel.loop(3,n3-1,3,new Parallel.LoopInt() { // i1 = 2, 5, 8, ...
      public void compute(int i3) {
        computeCSlice3(i3,dp1,dp2,dp3,ph1,ph2,ph3,dps,cvs);
    }});
  }

  private void computeGSlice3(int i3, float[][][] ph, 
    float[][][] ph1, float[][][] ph2, float[][][] ph3, float[][][] phs,
    float[][][] dp1, float[][][] dp2, float[][][] dp3) {
    int n2 = ph[0].length;
    int n1 = ph[0][0].length;
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      if(_mark[i3][i2][i1]==0){continue;}
      int i1p = i1+1;
      int i1m = i1-1;
      int i2p = i2+1;
      int i2m = i2-1;
      int i3p = i3+1;
      int i3m = i3-1;
      float ph1p = ph[i3 ][i2 ][i1p];
      float ph1m = ph[i3 ][i2 ][i1m];
      float ph2p = ph[i3 ][i2p][i1 ];
      float ph2m = ph[i3 ][i2m][i1 ];
      float ph3p = ph[i3p][i2 ][i1 ];
      float ph3m = ph[i3m][i2 ][i1 ];
      float phii = ph[i3 ][i2 ][i1 ];
      float dph1 = (ph1p-ph1m)*0.5f; 
      float dph2 = (ph2p-ph2m)*0.5f; 
      float dph3 = (ph3p-ph3m)*0.5f; 
      phs[i3][i2][i1] = ph1p+ph1m+ph2p+ph2m+ph3p+ph3m-phii*6;
      float psi = sqrt(dph1*dph1+dph2*dph2+dph3*dph3);
      float dsi = distance(psi)-1;
      dp1[i3][i2][i1] = dsi*dph1;
      dp2[i3][i2][i1] = dsi*dph2;
      dp3[i3][i2][i1] = dsi*dph3;
      psi = 1f/(psi+0.0000000001f);
      ph1[i3][i2][i1] = dph1*psi;
      ph2[i3][i2][i1] = dph2*psi;
      ph3[i3][i2][i1] = dph3*psi;
    }}
  }

  private void computeCSlice3(int i3, 
    float[][][] ph1, float[][][] ph2, float[][][] ph3,  
    float[][][] dp1, float[][][] dp2, float[][][] dp3,  
    float[][][] cvs, float[][][] dps) { 
    int n2 = dp1[0].length;
    int n1 = dp1[0][0].length;
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      if(_mark[i3][i2][i1]==0){continue;}
      int i3p = i3+1;
      int i3m = i3-1;
      int i2p = i2+1;
      int i2m = i2-1;
      int i1p = i1+1;
      int i1m = i1-1;
      float ddp1 = dp1[i3][i2][i1p]-dp1[i3][i2][i1m];
      float ddp2 = dp2[i3][i2p][i1]-dp2[i3][i2m][i1];
      float ddp3 = dp3[i3p][i2][i1]-dp3[i3m][i2][i1];
      float dph1 = ph1[i3][i2][i1p]-ph1[i3][i2][i1m];
      float dph2 = ph2[i3][i2p][i1]-ph2[i3][i2m][i1];
      float dph3 = ph3[i3p][i2][i1]-ph3[i3m][i2][i1];
      dps[i3][i2][i1] = 0.5f*(ddp1+ddp2+ddp3);
      cvs[i3][i2][i1] = 0.5f*(dph1+dph2+dph3);
    }}
  }

  private int _r;
  private float _mu;
  private float _niter;
  private float _lamda;
  private float _alpha;
  private byte _mark[][][]=null;

  private void updateBand(int d, float[][][] phi) {
    int n3 = phi.length;
    int n2 = phi[0].length;
    int n1 = phi[0][0].length;
    byte[][][] mkt = new byte[n3][n2][n1];
    if (_mark!=null) {
      mkt = copy(_mark);
      zero(_mark);
    } else {
      _mark = new byte[n3][n2][n1];
    }
    for (int i3=1; i3<n3-1; i3++) {
    for (int i2=1; i2<n2-1; i2++) {
    for (int i1=1; i1<n1-1; i1++) {
      float ph3m = phi[i3-1][i2][i1  ];
      float ph3p = phi[i3+1][i2][i1  ];
      float ph2m = phi[i3][i2-1][i1  ];
      float ph2p = phi[i3][i2+1][i1  ];
      float ph1m = phi[i3][i2  ][i1-1];
      float ph1p = phi[i3][i2  ][i1+1];
      if (ph1m*ph1p<0||ph2m*ph2p<0||ph3m*ph3p<0) {
        int b1 = i1-d; b1=max(b1,0);
        int b2 = i2-d; b2=max(b2,0);
        int b3 = i3-d; b3=max(b3,0);
        int e1 = i1+d; e1=min(e1,n1-1);
        int e2 = i2+d; e2=min(e2,n2-1);
        int e3 = i3+d; e3=min(e3,n3-1);
        for (int k3=b2; k3<=e3; ++k3) {
        for (int k2=b2; k2<=e2; ++k2) {
        for (int k1=b1; k1<=e1; ++k1) {
          if (mkt[k3][k2][k1]==0) {
            if (phi[k3][k2][k1]<0) phi[k3][k2][k1] = -d;
            if (phi[k3][k2][k1]>0) phi[k3][k2][k1] =  d;
          }
        }}}
        for (int k3=b3+1; k3<e3; ++k3) {
        for (int k2=b2+1; k2<e2; ++k2) {
        for (int k1=b1+1; k1<e1; ++k1) {
          _mark[k3][k2][k1] = (byte)1;
        }}}
      }
    }}}
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
