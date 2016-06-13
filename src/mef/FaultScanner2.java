/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package mef;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

/**
 * Computes fault likelihoods, and dips, by scanning over fault
 * orientations. Fault likelihoods are in the range [0,1], where 0 and 1
 * denote lowest and highest likelihoods, respectively. Computed fault strike
 * and dip angles are those for which maximum fault likelihoods occurred. 
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2016.04.29
 */
public class FaultScanner2 {

  /**
   * Constructs a fault scanner with specified parameters.
   * @param sigmaTheta half-width for smoothing along dip of faults.
   */
  public FaultScanner2(double sigmaTheta) {
    _sigmaTheta = sigmaTheta;
  }

  /**
   * Gets a sampling of fault dip theta appropriate for this scanner.
   * @param thetaMin minimum fault dip, in degrees.
   * @param thetaMax maximum fault dip, in degrees.
   */
  public Sampling getThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigmaTheta,thetaMin,thetaMax);
  }

  /**
   * Returns slopes and planarities of features in a specified image.
   * Image features are assumed to be locally linear, with slopes that may
   * vary throughout the image. Smoothing parameters control the extents of
   * Gaussian windows within which slope is estimated for each image sample.
   * Typically, because slopes tend to vary most rapidly in the 2nd
   * dimension, the extent of smoothing for the 1st dimension should be
   * greater than that for the 2nd dimension.
   * @param sigma1 half-width for smoothing in 1st dimension.
   * @param sigma2 half-width for smoothing in 2nd dimension.
   * @param slopeMax upper bound on computed slopes.
   * @param f input image for which to compute slopes.
   * @return array {p2,el} of slopes and planarities.
   */
  public static float[][][] slopes(
      double sigma1, double sigma2, 
      double slopeMax, float[][] f) {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float p2min = (float)(-slopeMax);
    final float p2max = (float)( slopeMax);

    // Normal vectors.
    final float[][] u1 = new float[n2][n1];
    final float[][] u2 = new float[n2][n1];
    final float[][] el = new float[n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    lof.applyForNormalLinear(f,u1,u2,el);

    // Slopes from normal vectors.
    final float[][] p2 = u2;
    loop(n2,new LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        if (-u2i<p2min*u1i) u2i = -p2min*u1i;
        if (-u2i>p2max*u1i) u2i = -p2max*u1i;
        if (u1i==0.0f) {
            p2[i2][i1] = (u2i<0.0f)?p2max:p2min;
        } else {
          p2[i2][i1] = -u2i/u1i;
        }
      }
    }});
    return new float[][][]{p2,el};
  }

  /**
   * Returns an image with specified samples taper to zero at edges.
   * Tapering enables simple zero-value boundary conditions to be
   * used in fault scanning without artifacts caused by abrupt 
   * truncations of strong image features near edges.
   * @param m1 width of the tapered samples near edges in 1st dimension.
   * @param m2 width of the tapered samples near edges in 2nd dimension.
   * @param f input image.
   * @return the tapered image.
   */
  public static float[][] taper(int m1, int m2, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = copy(f);
    float[] t1 = new float[m1];
    float[] t2 = new float[m2];
    for (int i1=0; i1<m1; ++i1)
      t1[i1] = (float)(0.54+0.46*cos(PI*(m1-i1)/m1));
    for (int i2=0; i2<m2; ++i2)
      t2[i2] = (float)(0.54+0.46*cos(PI*(m2-i2)/m2));
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0,j1=n1-1; i1<m1; ++i1,--j1) {
        float ti = t1[i1];
        g[i2][i1] *= ti;
        g[i2][j1] *= ti;
      }
    }
    for (int i2=0,j2=n2-1; i2<m2; ++i2,--j2) {
      float ti = t2[i2];
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] *= ti;
        g[j2][i1] *= ti;
      }
    }
    return g;
  }

  /**
   * Scans a specified image for fault strikes and dips.
   * @param thetaMin minimum fault dip, in degrees.
   * @param thetaMax maximum fault dip, in degrees.
   * @param p2 slopes in the 2nd dimension.
   * @param g the image to be scanned.
   * @return array {fl,fp,ft} of fault likelihoods, strikes, and dips.
   */
  public float[][][] scan(
      double thetaMin, double thetaMax,
      float[][] p2, float[][] g) {
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    return scan(st,p2,g);
  }

  public float[][][] scan(
      double thetaMin, double thetaMax,
      float sig1, float sig2, float smooth, float[][] g) {
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    return scan(st,sig1,sig2,smooth,g);
  }


  /**
   * Scans with the specified sampling of fault strikes and dips.
   * @param thetaSampling sampling of fault dip angles, in degrees.
   * @param p2 slopes in the 2nd dimension.
   * @param g the image to be scanned.
   * @return array {fl,fp,ft} of fault likelihoods, strikes, and dips.
   */
  public float[][][] scan(
      Sampling thetaSampling, float[][] p2, float[][] g) {
    float[][][] snd = semblanceNumDen(p2,g);
    return scanTheta(thetaSampling,snd);
  }

  public float[][][] scan(
      Sampling thetaSampling, float sig1, float sig2, float smooth, float[][] g) {
    float[][][] snd = semblanceNumDen(sig1,sig2, smooth,g);
    return scanTheta(thetaSampling,snd);
  }


  /**
   * Thins fault images to include only ridges in fault likelihoods.
   * After thinning, may be only one voxel wide. Thinned fault strikes and
   * dips are set to zero where thinned fault likelihoods are zero.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes, and dips.
   * @return array {fl,fp,ft} of thinned fault likelihoods, strikes, and dips.
   */
  public float[][][] thin(float[][][] flpt) {
    int n1 = flpt[0][0].length;
    int n2 = flpt[0].length;
    float[][] f = flpt[0];
    float[][] t = flpt[1];
    f = copy(f);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.applyX0(f,f);
    float[][] ff = new float[n2][n1];
    float[][] tt = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float[] fm = f[i2m];
      float[] f0 = f[i2 ];
      float[] fp = f[i2p];
      float[] t0 = t[i2 ];
      for (int i1=0; i1<n1; ++i1) {
        float f00 = f0[i1];
        float t00 = t0[i1];
        if ((fm[i1]<f00 && fp[i1]<f00) ||
            (fp[i1]<f00 && fm[i1]<f00)){
          ff[i2][i1] = f00;
          tt[i2][i1] = t00;
        } else {
          tt[i2][i1] = NO_DIP;
        }
      }
    }
    return new float[][][]{ff,tt};
  }

  /**
   * Applies structure-oriented smoothing limited by fault likelihoods.
   * For this method, faults are assumed to exist and smoothing stops
   * at samples where fault likelihood exceeds a specified value.
   * This method is usually applied using thinned fault likelihoods.
   * @param flstop smoothing stops where fault likelihood &gt; this value.
   * @param sigma smoothing radius (except near faults).
   * @param p2 array of slopes in 2nd dimension.
   * @param fl array of fault likelihoods, typically thinned.
   * @param g image to be smoothed.
   */
  public static float[][][] smooth(
      double flstop, double sigma, float[][][] p2, float[][][] p3, 
      float[][][] fl, float[][][] g) {
    int n1 = g[0][0].length;
    int n2 = g[0].length;
    int n3 = g.length;
    EigenTensors3 d = new EigenTensors3(n1,n2,n3,true);
    d.setEigenvalues(0.001f,1.00f,1.00f);
    float[][][] s = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          s[i3][i2][i1] = (fl[i3][i2][i1]<flstop)?1.0f:0.0f;
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float u1i = 1.0f/sqrt(1.0f+p2i*p2i+p3i*p3i);
          float u2i = -p2i*u1i;
          float u3i = -p3i*u1i;
          float usi = 1.0f/sqrt(u1i*u1i+u2i*u2i);
          float w1i = -u2i*usi;
          float w2i =  u1i*usi;
          float w3i = 0.0f;
          d.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          d.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
        }
      }
    }
    float c = (float)(0.5*sigma*sigma);
    float[][][] h = new float[n3][n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,c,s,g,h);
    return h;
  }

  /**
   * Adjusts dips for a specified aspect ratio dz/dx. A fault scanner computes
   * fault dips in degrees measured in sample coordinates. Because vertical
   * image sampling intervals are often less than horizontal sampling
   * intervals, fault dips measured in sample coordinates tend to be greater
   * than those in physical coordinates. This method converts fault dips to
   * degrees measured in physical coordinates, using the specified ratio of
   * vertical-to-horizontal sampling intervals.
   * @param dzdx ratio of vertical to horizontal sampling intervals.
   * @param ft array of fault dips measured in sample coordinates.
   * @return array of fault dips measured in physical coordinates.
   */
  public static float[][] convertDips(double dzdx, float[][] ft) {
    float scale = (float)dzdx;
    int n1 = ft[0].length;
    int n2 = ft.length;
    float[][] gt = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fti = ft[i2][i1];
        if (fti!=NO_DIP)
          gt[i2][i1] = toDegrees(atan(scale*tan(toRadians(fti))));
      }
    }
    return gt;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _sigmaTheta;

  private static final float NO_DIP    = -0.00001f;

  private static void trace(String s) {
    System.out.println(s);
  }


  private Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigmaTheta,thetaMin,thetaMax);
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

  // Shear horizontally such that q(i1,i2) = p(i1,i2+s*i1).
  private static float[][] shear(
    SincInterpolator si, double s, float[][] p)
  {
    int n1 = p[0].length;
    int n2p = p.length;
    int n2q = n2p+(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] q = new float[n2q][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2p; ++i2)
        pp[i2] = p[i2][i1];
      double f2q = (s<0.0f)?s*i1:s*i1-dqp;
      si.interpolate(n2p,1.0,0.0,pp,n2q,1.0f,f2q,qq);
      for (int i2=0; i2<n2q; ++i2)
        q[i2][i1] = qq[i2];
    }
    return q;
  }

  // Unshear horizontally such that p(i1,i2) = q(i1,i2-s*i1).
  private static float[][] unshear(
    SincInterpolator si, double s, float[][] q)
  {
    int n1 = q[0].length;
    int n2q = q.length;
    int n2p = n2q-(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] p = new float[n2p][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2q; ++i2)
        qq[i2] = q[i2][i1];
      double f2p = (s<0.0f)?-s*i1:-s*i1+dqp;
      si.interpolate(n2q,1.0,0.0,qq,n2p,1.0f,f2p,pp);
      for (int i2=0; i2<n2p; ++i2)
        p[i2][i1] = pp[i2];
    }
    return p;
  }


  // Smoothing filter
  private static RecursiveExponentialFilter makeRef(double sigma) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    return ref;
  }

  // For one fault strike, scans over all fault dips theta. The fault strike
  // vector has already been aligned with image axis 2, and semblance
  // numerators and denominators have already been smoothed in that direction.
  // Therefore, this scan over fault dip can be performed independently for
  // each i2. For each fault dip theta, this method shears semblance num and
  // den to align any faults having that dip with image axis 1. Semblance
  // num and den are then smoothed vertically, with an extent sigma that is
  // dip-adjusted (shorter for smaller fault dips), so that after unshearing
  // the extent of smoothing is roughly the same for all fault dips.
  private float[][][] scanTheta(Sampling thetaSampling, float[][][] snd) {
    final int n2 = snd[0].length;
    final int n1 = snd[0][0].length;
    final float[][] sn = snd[0];
    final float[][] sd = snd[1];
    final float[][] f = new float[n2][n1];
    final float[][] t = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    int nt = thetaSampling.getCount();
    for (int it=0; it<nt; ++it) {
      System.out.println(it+"/"+(nt-1)+" done...");
      float ti = (float)thetaSampling.getValue(it);
      float theta = toRadians(ti);
      float shear = -1.0f/tan(theta);
      float[][] sns = shear(si,shear,sn);
      float[][] sds = shear(si,shear,sd);
      float sigma = (float)_sigmaTheta*sin(theta);
      RecursiveExponentialFilter ref = makeRef(sigma);
      ref.apply1(sns,sns);
      ref.apply1(sds,sds);
      float[][] ss = semblanceFromNumDen(sns,sds);
      float[][] s2 = unshear(si,shear,ss);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float st = s2[i2][i1]; // semblance
        st = st*st; // semblance^2
        st = st*st; // semblance^4
        st = st*st; // semblance^8
        float fi = 1.0f-st;
        if (fi>f[i2][i1]) {
          f[i2][i1] = fi;
          t[i2][i1] = ti;
        }
      }}
    }

    for (int it=0; it<nt; ++it) {
      System.out.println(it+"/"+(nt-1)+" done...");
      float ti = (float)thetaSampling.getValue(it);
      float theta = toRadians(ti);
      float shear = 1.0f/tan(theta);
      float[][] sns = shear(si,shear,sn);
      float[][] sds = shear(si,shear,sd);
      float sigma = (float)_sigmaTheta*sin(theta);
      RecursiveExponentialFilter ref = makeRef(sigma);
      ref.apply1(sns,sns);
      ref.apply1(sds,sds);
      float[][] ss = semblanceFromNumDen(sns,sds);
      float[][] s2 = unshear(si,shear,ss);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float st = s2[i2][i1]; // semblance
        st = st*st; // semblance^2
        st = st*st; // semblance^4
        st = st*st; // semblance^8
        float fi = 1.0f-st;
        if (fi>f[i2][i1]) {
          f[i2][i1] = fi;
          t[i2][i1] = -ti;
        }
      }}
    }
    return new float[][][]{f,t};
  }


  // Computes fault semblance numerators and denominators.
  private static float[][][] semblanceNumDen(
    float sig1, float sig2, float smooth, float[][] f) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] sn = new float[n2][n1];
    float[][] sd = new float[n2][n1];
    float[][] fs = mul(f,f);
    LocalOrientFilter lof = new LocalOrientFilter(sig1,sig2);
    EigenTensors2 ets = lof.applyForTensors(f);
    ets.setEigenvalues(0.0001f,1.0f);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,smooth,f, sn);
    lsf.apply(ets,smooth,fs,sd);
    return new float[][][]{mul(sn,sn),sd};
  }


  // Computes fault semblance numerators and denominators.
  private static float[][][] semblanceNumDen(
    final float[][] p2, final float[][] f) 
  {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] sn = new float[n2][n1];
    final float[][] sd = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    loop(n2,new LoopInt() {
    public void compute(int i2) {
      float[] xm = new float[n1];
      float[] xp = new float[n1];
      float[] gm = new float[n1];
      float[] gp = new float[n1];
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float[] fm = f[i2m];
      float[] f0 = f[i2 ];
      float[] fp = f[i2p];
      float[] p2m = p2[i2m];
      float[] p2p = p2[i2p];
      float[] sn32 = sn[i2];
      float[] sd32 = sd[i2];
      for (int i1=0; i1<n1; ++i1) {
          xm[i1] = i1-p2m[i1];
          xp[i1] = i1+p2p[i1];
      }
      si.interpolate(n1,1.0,0.0,fm,n1,xm,gm);
      si.interpolate(n1,1.0,0.0,fp,n1,xp,gp);
      float[] hm = gm, h0=f0, hp = gp;
      if (            i2==0   ) hm = h0;
      if (            i2==n2-1) hp = h0;
      for (int i1=0; i1<n1; ++i1) {
        float hmi = hm[i1];
        float h0i = h0[i1];
        float hpi = hp[i1];
        float sumn = hmi+h0i+hpi;
        float sumd = hmi*hmi+h0i*h0i+hpi*hpi;
        sn32[i1] = sumn*sumn;
        sd32[i1] = 3.0f*sumd;
      }
    }});
    return new float[][][]{sn,sd};
  }

  // Computes semblance ratios from numerators and denominators.
  // Takes care to ensure that semblances are in range [0,1].
  private static float[][] semblanceFromNumDen(float[][] sn, float[][] sd) {
    int n1 = sn[0].length;
    int n2 = sn.length;
    float[][] sr = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] sn2 = sn[i2];
      float[] sd2 = sd[i2];
      float[] sr2 = sr[i2];
      for (int i1=0; i1<n1; ++i1) {
        float sni = sn2[i1];
        float sdi = sd2[i1];
        if (sdi<=0.0f || sni<=0.0f) {
          sr2[i1] = 0.0f;
        } else if (sdi<sni) {
          sr2[i1] = 1.0f;
        } else {
          sr2[i1] = sni/sdi;
        }
      }
    }
    return sr;
  }

}
