/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package spv;

import java.util.*;
import javax.imageio.ImageIO;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;
import util.*;

/**
 * Scan for fault orientations.
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.07.15
 */

public class FaultOrientScanner2 {
  /**
   * Constructs a scanner with specified parameters.
   * @param sigmaTheta half-width for smoothing along dip of edges.
   */
  public FaultOrientScanner2(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  /**
   * Gets a sampling of edge dip theta appropriate for this scanner.
   * @param thetaMin minimum edge dip, in degrees.
   * @param thetaMax maximum edge dip, in degrees.
   */
  public Sampling getThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigma1,thetaMin,thetaMax);
  }

  /**
   * Scans a specified image for edge dips.
   * @param thetaMin minimum edge dip, in degrees.
   * @param thetaMax maximum edge dip, in degrees.
   * @param g the image to be scanned.
   * @return array {el,et} of edge gradients and dips.
   */
  public float[][][] scan(double thetaMin, double thetaMax, float[][] g) {
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    return scanTheta(st,g);
  }

  /**
   * Scans a specified color image for edge dips.
   * @param thetaMin minimum edge dip, in degrees.
   * @param thetaMax maximum edge dip, in degrees.
   * @param g the image to be scanned.
   * @return array {el,et} of edge gradients and dips.
   */
  public float[][][] scanColorImage(
    double thetaMin, double thetaMax, float[][][] g) {
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    return scanThetaColor(st,g);
  }


  /**
   * Thins fault images to include only ridges in image edges.
   * After thinning, may be only one voxel wide. Thinned fault strikes and
   * dips are set to zero where thinned fault likelihoods are zero.
   * @param fet array {el,et} of edge gradeitns and dips.
   * @return array {elt,ett} of thinned edge gradeints and dips.
   */
  public float[][][] thin(float[][][] fet) {
    int n1 = fet[0][0].length;
    int n2 = fet[0].length;
    float[][] f = fet[0];
    float[][] t = fet[1];
    float[][] ff = new float[n2][n1];
    float[][] tt = new float[n2][n1];
    float pi = (float)(Math.PI/180.0);
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float ti = t[i2][i1]+90f;
      if(ti>180) ti -= 180;
      ti *= pi;
      float d1 = 0f;
      float d2 = 0f;
      d1 = cos(ti);
      d2 = sin(ti);
      float x1p = i1+d1;
      float x2p = i2+d2;
      float x1m = i1-d1;
      float x2m = i2-d2;
      float fi = f[i2][i1];
      float fp = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,x1p,x2p);
      float fm = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,x1m,x2m);
      if(fp<fi&&fm<fi) {
        ff[i2][i1] = fi;
        tt[i2][i1] = t[i2][i1];
      }
    }}
    return new float[][][]{ff,tt};
  }

  public float[][] edgeLikeFit2(int r, float[][] fl) {
    int n2 = fl.length;
    int n1 = fl[0].length;
    float[][] flr = new float[n2][n1];
    for (int i2=r; i2<n2-r-1; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float[] ft = new float[r*2+1];
      for (int r2=-r;r2<=r;r2++)
        ft[r2+r] = fl[i2+r2][i1];
      float[] abc = parabolicFit(ft);
      float a = abc[0];
      float b = abc[1];
      float s = -abs(2*a*r+b)/(2*a);
      flr[i2][i1] = fl[i2][i1]/(s+0.0001f);
    }}
    return flr;
  }

  public float[][] edgeLikeFit(int r, float[][] el, float[][] et) {
    int n2 = el.length;
    int n1 = el[0].length;
    float[][] elr = new float[n2][n1];
    float pi = (float)(Math.PI/180.0);
    SincInterpolator si = new SincInterpolator();
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float[] ft = new float[r*2+1];
      float ti = et[i2][i1]+90f;
      if(ti>180) ti -= 180;
      ti *= pi;
      float t1 = cos(ti);
      float t2 = sin(ti);
      for (int ir=-r;ir<=r;ir++)
        ft[ir+r] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,el,i1+ir*t1,i2+ir*t2);
      float[] abc = parabolicFit(ft);
      float a = abc[0];
      float b = abc[1];
      float s = -abs(2*a*r+b)/(2*a);
      elr[i2][i1] = el[i2][i1]/(s+0.0001f);
    }}
    return elr;
  }


  public float[][] edgeLikeFitG(int r, float[][] el, float[][] et) {
    int n2 = el.length;
    int n1 = el[0].length;
    float[][] elr = new float[n2][n1];
    float pi = (float)(Math.PI/180.0);
    SincInterpolator si = new SincInterpolator();
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(r);
    float[] ft = new float[10*r+1];
    float[] fs = new float[10*r+1];
    float[] f1 = new float[10*r+1];
    float[] f2 = new float[10*r+1];
    int hr = 5*r;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float ti = et[i2][i1]+90f;
      if(ti>180) ti -= 180;
      ti *= pi;
      float t1 = cos(ti);
      float t2 = sin(ti);
      for (int ir=-hr;ir<=hr;ir++)
        ft[ir+hr] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,el,i1+ir*t1,i2+ir*t2);
      rgf.apply0(ft,fs);
      rgf.apply1(ft,f1);
      rgf.apply2(ft,f2);
      elr[i2][i1] = fs[hr]*(-f2[hr]/(abs(f1[hr])+0.0001f));
    }}
    return elr;
  }

  public float[] parabolicFit(float[] f) {
    int n1 = f.length;
    double[][] A = new double[3][3];
    double[][] B = new double[3][1];
    double s4=0.0;
    double s3=0.0;
    double s2=0.0;
    double s1=0.0;
    double s0=0.0;
    double b0=0.0;
    double b1=0.0;
    double b2=0.0;
    for (int i1=0; i1<n1; ++i1) {
      double x1 = i1;
      double x2 = i1*x1;
      double x3 = i1*x2;
      double x4 = i1*x3;
      s0 += 1.;
      s1 += x1;
      s2 += x2;
      s3 += x3;
      s4 += x4;
      b0 += x2*f[i1];
      b1 += x1*f[i1];
      b2 += f[i1];
    }
    A[0][0] = s4;
    A[1][0] = s3;
    A[2][0] = s2;
    A[0][1] = s3;
    A[1][1] = s2;
    A[2][1] = s1;
    A[0][2] = s2;
    A[1][2] = s1;
    A[2][2] = s0;
    B[0][0] = b0;
    B[1][0] = b1;
    B[2][0] = b2;
    DMatrix da = new DMatrix(A);
    DMatrix db = new DMatrix(B);
    DMatrix dx = da.solve(db);
    double[][] x = dx.get();
    float a = (float)x[0][0];
    float b = (float)x[1][0];
    float c = (float)x[2][0];
    return new float[]{a,b,c};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1;
  private float _sigma2;

  private static final float NO_DIP    = -0.00001f;

  private Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigma1,thetaMin,thetaMax);
  }
  private static Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  // Smoothing filter
  private static RecursiveExponentialFilter makeRef(double sigma) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    return ref;
  }

  // Scans over all edge dips theta. 
  private float[][][] scanTheta(Sampling thetaSampling, float[][] g) {
    final int n2 = g.length;
    final int n1 = g[0].length;
    final float[][] f = new float[n2][n1];
    final float[][] t = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    int nt = thetaSampling.getCount();
    float[][][] rsc = getAngleMap(n1,n2);
    RecursiveExponentialFilter ref = makeRef(_sigma1);
    for (int it=0; it<nt; ++it) {
      System.out.println(it+"/"+(nt-1)+" done...");
      float ti = (float)thetaSampling.getValue(it);
      float theta = toRadians(ti);
      float[][] gr = rotate(theta,rsc[0],rsc[1],rsc[2],g);
      ref.apply1(gr,gr);
      float[][] s2 = unrotate(n1,n2,theta,rsc[0],rsc[1],rsc[2],gr);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float st = s2[i2][i1];
        if (st>f[i2][i1]) {
          f[i2][i1] = st;
          t[i2][i1] = ti;
        }
      }}
    }
    return new float[][][]{f,t};
  }

  private float[][][] scanThetaColor(Sampling thetaSampling, float[][][] g) {
    final int n2 = g[0].length;
    final int n1 = g[0][0].length;
    final float[][] f = new float[n2][n1];
    final float[][] t = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    int nt = thetaSampling.getCount();
    float[][][] rsc = getAngleMap(n1,n2);
    float ss = _sigma2*_sigma2;
    float a = (1.0f+ss-sqrt(1.0f+2.0f*ss))/ss;
    RecursiveExponentialFilter ref = makeRef(_sigma1);
    for (int it=0; it<nt; ++it) {
      System.out.println(it+"/"+(nt-1)+" done...");
      float ti = (float)thetaSampling.getValue(it);
      float theta = toRadians(ti);
      float[][][] gr = rotate(theta,rsc[0],rsc[1],rsc[2],g);
      ref.apply1(gr,gr);
      float[][][] gg = computeGradient(a,gr);
      float[][][] gu = unrotate(n1,n2,theta,rsc[0],rsc[1],rsc[2],gg);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g0 = gu[0][i2][i1];
        float g1 = gu[1][i2][i1];
        float g2 = gu[2][i2][i1];
        float gm = max(g0,g1,g2);
        if (gm>f[i2][i1]) {
          f[i2][i1] = gm;
          t[i2][i1] = ti;
        }
      }}
    }
    return new float[][][]{f,t};
  }


  public float[][][] getAngleMap(int n1, int n2) {
    int h2 = (n2-1)/2;
    int h1 = (n1-1)/2;
    int nr = 2*round(sqrt(h1*h1+h2*h2))+1;
    int nh = (nr-1)/2;
    float[][] sr = new float[nr][nr];
    float[][] cr = new float[nr][nr];
    float[][] rr = new float[nr][nr];
    for (int k2=-nh; k2<=nh; ++k2) {
    for (int k1=-nh; k1<=nh; ++k1) {
      int i1 = k1+nh;
      int i2 = k2+nh;
      float ri = sqrt(k1*k1+k2*k2);
      float rc = 1.0f/ri;
      rr[i2][i1] = ri;
      sr[i2][i1] = k1*rc;
      cr[i2][i1] = k2*rc;
    }}
    return new float[][][]{rr,sr,cr};
  }

  public float[][][] rotate(float theta, 
    float[][] rr, float[][] sr, float[][] cr, float[][][] gx) {
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    int h2 = (n2-1)/2;
    int h1 = (n1-1)/2;
    int nr = 2*round(sqrt(h1*h1+h2*h2))+1;
    int nh = (nr-1)/2;
    float[][][] gr = new float[3][nr][nr];
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int k2=-nh; k2<=nh; ++k2) {
    for (int k1=-nh; k1<=nh; ++k1) {
      int i2 = k2+nh;
      int i1 = k1+nh;
      float st = sin(theta);
      float ct = cos(theta);
      float rk = rr[i2][i1]; 
      float sk = sr[i2][i1]; 
      float ck = cr[i2][i1]; 
      float sp = sk*ct-ck*st;
      float cp = ck*ct+sk*st;
      float x1 = rk*sp+h1;
      float x2 = rk*cp+h2;
      gr[0][i2][i1] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,gx[0],x1,x2);
      gr[1][i2][i1] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,gx[1],x1,x2);
      gr[2][i2][i1] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,gx[2],x1,x2);
    }}
    return gr;
  }


  public float[][] rotate(float theta, 
    float[][] rr, float[][] sr, float[][] cr, float[][] gx) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    int h2 = (n2-1)/2;
    int h1 = (n1-1)/2;
    int nr = 2*round(sqrt(h1*h1+h2*h2))+1;
    int nh = (nr-1)/2;
    float[][] gr = new float[nr][nr];
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int k2=-nh; k2<=nh; ++k2) {
    for (int k1=-nh; k1<=nh; ++k1) {
      int i2 = k2+nh;
      int i1 = k1+nh;
      float st = sin(theta);
      float ct = cos(theta);
      float rk = rr[i2][i1]; 
      float sk = sr[i2][i1]; 
      float ck = cr[i2][i1]; 
      float sp = sk*ct-ck*st;
      float cp = ck*ct+sk*st;
      float x1 = rk*sp+h1;
      float x2 = rk*cp+h2;
      gr[i2][i1] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,gx,x1,x2);
    }}
    return gr;
  }

  public float[][] unrotate(int n1, int n2, float theta, 
    float[][] rr, float[][] sr, float[][] cr, float[][] gr) {
    int h2 = (n2-1)/2;
    int h1 = (n1-1)/2;
    int nr = 2*round(sqrt(h1*h1+h2*h2))+1;
    int nh = (nr-1)/2;
    float[][] gu = new float[n2][n1];
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int k2=-h2; k2<=h2; ++k2) {
    for (int k1=-h1; k1<=h1; ++k1) {
      int i2 = k2+nh;
      int i1 = k1+nh;
      float st = sin(theta);
      float ct = cos(theta);
      float rk = rr[i2][i1]; 
      float sk = sr[i2][i1]; 
      float ck = cr[i2][i1]; 
      float sp = sk*ct+ck*st;
      float cp = ck*ct-sk*st;
      float x1 = rk*sp+nh;
      float x2 = rk*cp+nh;
      gu[k2+h2][k1+h1] = si.interpolate(nr,1.0,0.0,nr,1.0,0.0,gr,x1,x2);
    }}
    return gu;
  }

  public float[][][] unrotate(int n1, int n2, float theta, 
    float[][] rr, float[][] sr, float[][] cr, float[][][] gr) {
    int h2 = (n2-1)/2;
    int h1 = (n1-1)/2;
    int nr = 2*round(sqrt(h1*h1+h2*h2))+1;
    int nh = (nr-1)/2;
    float[][][] gu = new float[3][n2][n1];
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int k2=-h2; k2<=h2; ++k2) {
    for (int k1=-h1; k1<=h1; ++k1) {
      int i2 = k2+nh;
      int i1 = k1+nh;
      float st = sin(theta);
      float ct = cos(theta);
      float rk = rr[i2][i1]; 
      float sk = sr[i2][i1]; 
      float ck = cr[i2][i1]; 
      float sp = sk*ct+ck*st;
      float cp = ck*ct-sk*st;
      float x1 = rk*sp+nh;
      float x2 = rk*cp+nh;
      gu[0][k2+h2][k1+h1] = si.interpolate(nr,1.0,0.0,nr,1.0,0.0,gr[0],x1,x2);
      gu[1][k2+h2][k1+h1] = si.interpolate(nr,1.0,0.0,nr,1.0,0.0,gr[1],x1,x2);
      gu[2][k2+h2][k1+h1] = si.interpolate(nr,1.0,0.0,nr,1.0,0.0,gr[2],x1,x2);
    }}
    return gu;
  }


  private float[][] computeGradient(float a, float[][] gs) {
    int n2 = gs.length;
    int n1 = gs[0].length;
    float[][] gc = new float[n2][n1];
    float[][] ga = new float[n2][n1];
    float[][] gg = new float[n2][n1];
    causal2(a,gs,gc);
    anticausal2(a,gs,ga);
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      gg[i2][i1] = abs(gc[i2][i1]-ga[i2][i1]);
    return gg;
  }

  private float[][][] computeGradient(float a, float[][][] gs) {
    int n2 = gs[0].length;
    int n1 = gs[0][0].length;
    float[][][] gg = new float[3][n2][n1];
    for (int i3=0; i3< 3; ++i3) {
      float[][] gs3 = gs[i3];
      float[][] gg3 = gg[i3];
      float[][] gc = new float[n2][n1];
      float[][] ga = new float[n2][n1];
      causal2(a,gs3,gc);
      anticausal2(a,gs3,ga);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        gg3[i2][i1] = abs(gc[i2][i1]-ga[i2][i1]);
      }}
    }
    return gg;
  }



  // Horizontally causal and anti-causal filters are implemented 
  // by one-side recursive exponential filters.  
  private void causal2(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float b = 1.0f - a;
    for (int i1=0; i1<n1; ++i1){ 
      float yi = y[0][i1] = x[0][i1];
      for (int i2=1; i2<n2; ++i2) 
        y[i2][i1] = yi = a*yi + b*x[i2][i1];
    }
  }


 private void anticausal2(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float b = 1.0f - a;
    for(int i1=0; i1<n1; ++i1) {
      float yi = y[n2-1][i1] = x[n2-1][i1];
      for(int i2=n2-2; i2>=0; --i2)
        y[i2][i1] = yi = a*yi + b*x[i2][i1];
    }
  }


}
