/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;


import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Uses image samples alongside fault skins to estimate 2D fault dip-slips.
 * <p>
 * Before estimating slips, the seismic images should be smoothed along
 * reflectors, but not across faults. (Thinned fault likelihoods can be used
 * to stop smoothing at potential fault locations.) This smoothing does what
 * seismic interpreters do visually when estimating fault throws. In effect,
 * such smoothing enables us to use samples of a seismic image located away
 * from faults to estimate slips along those faults.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.21
 */
public class FaultCorrelater {

  /**
   * Constructs a fault slipper for the specified seismic image.
   * @param gs seismic image, smoothed up to (but not across) faults.
   */
  public FaultCorrelater(float[][] gs, float[][] p2) {
    _p2 = p2;
    _gs = gs;
  }

  public void setErrorPower(float power) {
    _power = power;
  }

  public void setOffset(float offset) {
    _offset = offset;
  }

    /**
   * Enables or disables a zero-slope reflector assumption. For educational
   * use, only. Ignoring reflector slopes can yield significant errors in
   * estimated dip slips. The default is false.
   * @param zeroSlope true, to assume zero slopes; false, otherwise.
   */
  public void setZeroSlope(boolean zeroSlope) {
    _zeroSlope = zeroSlope;
  }


  /**
   * Computes integer fault throws that correlate 
   * hanging-wall and foot-wall samples.
   * @param curves array of curves for which to compute throws.
   * @param smin an estimate for minimum fault throw, in samples.
   * @param smax an estimate for maximum fault throw, in samples.
   */
  public void computeThrow(FaultCurve[] curves, double smin, double smax) {
    for (FaultCurve curve:curves)
      computeThrow(curve,smin,smax);
  }

  /**
   * Computes integer fault throws that correlate 
   * hanging-wall and foot-wall samples.
   * @param curve the curve for which to compute throws.
   * @param smin an estimate for minimum fault throw, in samples.
   * @param smax an estimate for maximum fault throw, in samples.
   */

  public void computeThrow(FaultCurve curve, double smin, double smax) {
    Check.argument(smax>=0.0f,"smax not less than zero");
    FaultPoint[][] cab = curve.getCellsAB();
    int lmin = (int)smin;
    int lmax = (int)smax;
    DynamicWarping dw = new DynamicWarping(lmin,lmax);
    dw.setStrainMax(0.25); // TODO: always 0.25? goes with 4 below?
    computeAlignmentErrors(curve,lmin,lmax,_offset,_gs);
    extrapolateAlignmentErrors(lmin,lmax,cab);
    computeShifts(dw,cab);
    clearErrors(curve);
    for (int nsmooth=0; nsmooth<2; ++nsmooth) // TODO: 2?
      smoothShifts(curve);
    computeDipSlips(curve);
    findCorrespondPoint(curve);
  }

  public int[][][] getCorrelatedIndex(FaultCurve[] curves) {
    FaultPoint[] fps = FaultCurve.getPoints(curves);
    int np = fps.length;
    int[][][] mps = new int[2][2][np];
    for (int ip=0; ip<np; ++ip) {
      FaultPoint fpi = fps[ip];
      FaultPoint fci = fpi.cc;
      int i1i = fpi.i1;
      int i2m = fpi.i2m;
      int c1i = fci.i1 ;
      int c2p = fci.i2p;
      mps[0][0][ip] = i1i;
      mps[0][1][ip] = i2m;
      mps[1][0][ip] = c1i;
      mps[1][1][ip] = c2p;
    }
    return mps;
  }

  /**
   * Interpolates specified dip-slip vectors.
   * @param s array {s1,s2,s3} of dip-slip vectors.
   * @param smark the mark for slips not adjacent to a fault.
   * @return interpolated dip-slip vectors.
   */
  public static float[][][] interpolateDipSlips(
    float[][][] s, float smark) {
    int n2 = s[0].length;
    int n1 = s[0][0].length;
    float[][] p = new float[n2][n1];
    float[][] t = new float[n2][n1];
    short[][] k1 = new short[n2][n1];
    short[][] k2 = new short[n2][n1];
    float[][][] sq = new float[2][n2][n1];
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(smark,s[0],t,k1,k2);
    clip(0.0f,100.0f,t,t);
    LocalDiffusionKernel.Stencil stencil = LocalDiffusionKernel.Stencil.D21;
    LocalDiffusionKernel ldk = new LocalDiffusionKernel(stencil);
    BlendedGridder2 bg = new BlendedGridder2();
    bg.setBlendingKernel(ldk);
    bg.setSmoothness(0.5);
    for (int is=0; is<2; ++is) {
      float[][] si = s[is];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          int j1 = k1[i2][i1];
          int j2 = k2[i2][i1];
          p[i2][i1] = si[j2][j1];
        }
      }
      bg.gridBlended(t,p,sq[is]);
    }
    return sq;
  }

  /**
   * Unfaults an image using interpolated dip-slip vectors.
   * @param s array {s1,s2,s3} of interpolated dip-slip vectors.
   * @param g image to be unfaulted.
   * @return unfaulted image.
   */
  public static float[][] unfault(float[][][] s, final float[][] g) {
    final int n1 = g[0].length;
    final int n2 = g.length;
    final float[][] s1 = s[0];
    final float[][] s2 = s[1];
    final float[][] gs = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x1 = i1+s1[i2][i1];
        float x2 = i2+s2[i2][i1];
        gs[i2][i1] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,g,x1,x2);
      }
    }});
    return gs;
  }

  public void findCorrespondPoint(FaultCurve curve) {
    for (FaultPoint point:curve) {
      FaultPoint cellBegin = point;
      float smp = point.smp;

      // Begin at point location.
      float[] y = {point.x1,point.x2};

      // Walk down-dip (for normal fault) or up-dip (for reverse fault).
      if (smp>0.0f) {
        for (; smp>=1.0f; smp-=1.0f)
          point = point.walkDownDipFrom(y);
      } else {
        for (; smp<=-1.0f; smp+=1.0f)
          point = point.walkUpDipFrom(y);
      }
      cellBegin.fc = point;
    }
  }



  public float[][][] getDipSlips(
    int n1, int n2, FaultCurve[] curves, float smark) {
    float[][] s1 = new float[n2][n1];
    float[][] s2 = new float[n2][n1];

    // Initially set all slip vectors to the specified mark.
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        s1[i2][i1] = smark;
        s2[i2][i1] = smark;
      }
    }
    for (FaultCurve curve:curves) {
    for (FaultPoint point:curve) {
      int i1 = point.i1;
      int m2 = point.i2m;
      int p2 = point.i2p;
      s1[m2][i1] = 0f;
      s2[m2][i1] = 0f;

      s1[p2][i1] = point.s1;
      s2[p2][i1] = point.s2;
    }}
    return new float[][][]{s1,s2};
  }



  

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _offset = 2;
  private boolean _zeroSlope; // if true, assume reflectors have zero slope
  private float[][] _gs,_p2; // seismic image (smoothed) and slopes
  private static float _power = 4f;


  private static void trace(String s) {
    System.out.println(s);
  }



  /**
   * Computes dip-slip vectors from vertical shifts for specified curve.
   */
  private void computeDipSlips(FaultCurve curve) {

    // For all cells in the curve, ...
    for (FaultPoint point:curve) {
      FaultPoint cellBegin = point;

      // Cell coordinates, dip vector, and shift.
      float x1 = point.x1;
      float x2 = point.x2;
      float u1 = point.u1;
      float u2 = point.u2;
      float smp = point.smp;

      // Offset vector d.
      float d1 = 0.0f;
      float d2 = _offset;
      if(point.ft<0) d2=-d2;

      // Reflector slopes at point x-d.
      float p2 = imageValueAt(x1-d1,x2-d2,_p2);

      // Unit-vector a normal to reflector.
      float a1 = 1.0f/sqrt(1.0f+p2*p2);
      float a2 = -p2*a1;

      // Slip adjustment for reflector slope on minus side.
      float am = (a1*d1+a2*d2)/(a1*u1+a2*u2);
      float u1m = am*u1;
      float u2m = am*u2;

      // Begin at point location.
      float[] y = {point.x1,point.x2};

      // Walk down-dip (for normal fault) or up-dip (for reverse fault).
      if (smp>0.0f) {
        for (; smp>=1.0f; smp-=1.0f)
          point = point.walkDownDipFrom(y);
      } else {
        for (; smp<=-1.0f; smp+=1.0f)
          point = point.walkUpDipFrom(y);
      }

      // Account for any remaining fractional shift.
      float y1 = y[0], y2 = y[1];
      y1 += smp;
      y2 += smp*point.us*point.u2;

      // Unit dip vector.
      u1 = point.u1;
      u2 = point.u2;

      // Offset vector d.
      d1 = 0.0f;
      d2 =  _offset;
      if(point.ft<0) d2=-d2;

      // Reflector slopes at point y+d.
      p2 = imageValueAt(y1+d1,y2+d2,_p2);

      // Unit-vector a normal to reflector.
      a1 = 1.0f/sqrt(1.0f+p2*p2);
      a2 = -p2*a1;

      // Slip adjustment for reflector slope on plus side.
      float ap = (a1*d1+a2*d2)/(a1*u1+a2*u2);
      float u1p = ap*u1;
      float u2p = ap*u2;

      // Record total dip slip in point at which we began the walk.
      cellBegin.s1 = y1-x1;
      cellBegin.s2 = y2-x2;
      if (!_zeroSlope) {
        cellBegin.s1 += u1m+u1p;
        cellBegin.s2 += u2m+u2p;
      }
    }
  }



    /**
   * Smooths shifts in the specified curve.
   */
  private static void smoothShifts(FaultCurve curve) {
    FaultPoint.Get1 getter = new FaultPoint.Get1() { 
      public float get(FaultPoint point) { return point.smp; }
    };
    FaultPoint.Set1 setter = new FaultPoint.Set1() { 
      public void set(FaultPoint point, float smp) { point.smp = smp; }
    };
    curve.smooth1(getter,setter);
  }



    /**
   * Uses dynamic warping to compute minus-plus shifts. Assumes that
   * minus-plus errors have already been computed and stored in the cells
   * referenced in the specified above-below and left-right arrays.
   */
  private static void computeShifts(
      DynamicWarping dw, FaultPoint[][] cab) {
    
    // Arrays of arrays of errors, linked above and below.
    int nab = cab.length;
    float[][][] eab = new float[nab][][];
    for (int iab=0; iab<nab; ++iab) {
      int mab = cab[iab].length;
      eab[iab] = new float[mab][];
      for (int jab=0; jab<mab; ++jab) {
        FaultPoint c = cab[iab][jab];
        eab[iab][jab] = c.emp;
      }
    }


    // Smooth alignment errors in above-below and left-right directions.
    for (int ismooth=0; ismooth<2; ++ismooth) { // TODO: how many?
      dw.smoothErrors1(eab,eab);
      normalizeErrors(eab); // TODO: helpful?
    }

    // Find shifts by accumulating once more and then backtracking.
    for (int iab=0; iab<nab; ++iab) {
      float[][] dab = dw.accumulateForward(eab[iab]);
      float[] s = dw.backtrackReverse(dab,eab[iab]);
      int mab = s.length;
      for (int jab=0; jab<mab; ++jab) {
        FaultPoint c = cab[iab][jab];
        c.smp = s[jab];
      }
    }
  }

    /**
   * Normalizes errors to account for varying lengths of arrays of cells.
   * This normalization is different from that performed by dynamic warping
   * when smoothing alignment errors, because that normalization assumes that
   * all arrays in the array of arrays of alignment errors have the same
   * length. This assumption is false for the arrays cellsAB and cellsLR in
   * fault skins.
   */
  private static void normalizeErrors(float[][][] e) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2) {
      int n1 = e[i2].length;
      float scale = 1.0f/n1;
      for (int i1=0; i1<n1; ++i1) {
        float[] ei = e[i2][i1];
        int nlag = ei.length;
        for (int ilag=0; ilag<nlag; ++ilag)
          ei[ilag] *= scale;
      }
    }
  }

  private static void clearErrors(FaultCurve curve) {
    for (FaultPoint point:curve)
      point.emp = null;
  }



    /**
   * Computes alignment errors and initializes shifts for specified curve.
   */
  private static void computeAlignmentErrors(
      FaultCurve curve, int lmin, int lmax, float offset, float[][] f) {
    for (FaultPoint point:curve)
      computeAlignmentErrors(point,lmin,lmax,offset,f);
  }


    /**
   * Computes minus-plus alignment errors for one point. These errors
   * correspond to differences between the sample value on the minus side of
   * the point and those for the plus sides of cells up and down dip from the
   * point.
   * <p> 
   * For lags where image sample values are unavailable (e.g., near image
   * boundaries), errors are extrapolated from other lags, but are negated, so
   * that extrapolated errors can be detected and modified later, after errors
   * for all relevant cells have been computed.
   */
  private static void computeAlignmentErrors(
    FaultPoint point, int lmin, int lmax, float offset, float[][] f) {
    Check.argument(lmin<=0,"lmin<=0");
    Check.argument(lmax>=0,"lmax>=0");
    int n1 = f[0].length;
    float[] y = new float[2];

    // New arrays for alignment errors.
    int lag0 = -lmin;
    float[] emp = point.emp = new float[1+lmax-lmin];

    // Errors for lag zero.
    float d2 =  offset;
    if(point.ft<0) d2=-offset;
    float y1 = point.x1, y2 = point.x2;
    float fm = imageValueAt(y1,y2-d2,f);
    float gp = imageValueAt(y1,y2+d2,f);
    float empl = emp[lag0] = alignmentError(fm,gp);

    // Errors for samples above; make any extrapolated errors negative.
    FaultPoint ca = point;
    int nlaga = min(-lmin,ca.i1);
    y1 = point.x1; y2 = point.x2;
    for (int ilag=1; ilag<=-lmin; ++ilag) {
      if (ilag<=nlaga) {
        y[0] = y1; y[1] = y2;
        ca = ca.walkUpDipFrom(y);
        y1 = y[0]; y2 = y[1];
        d2 = offset;
        if(ca.ft<0) d2=-offset;
        gp = imageValueAt(y1,y2+d2,f);
        empl = emp[lag0-ilag] = alignmentError(fm,gp);
      } else {
        emp[lag0-ilag] = -empl;
      }
    }

    // Errors for samples below; make any extrapolated errors negative.
    FaultPoint cb = point;
    int nlagb = min(lmax,n1-1-cb.i1);
    y1 = point.x1; y2 = point.x2;
    for (int ilag=1; ilag<=lmax; ++ilag) {
      if (ilag<=nlagb) {
        y[0] = y1; y[1] = y2;
        cb = cb.walkDownDipFrom(y);
        y1 = y[0]; y2 = y[1];
        d2 = offset;
        if(cb.ft<0) d2=-offset;
        gp = imageValueAt(y1,y2+d2,f);
        empl = emp[lag0+ilag] = alignmentError(fm,gp);
      } else {
        emp[lag0+ilag] = -empl;
      }
    }
  }

  /**
   * Extrapolates alignment errors emp where not computed. Errors that could
   * not be computed are negative, and are copies of errors for smaller lags
   * that could be computed. (Errors for zero lag can always be computed.) 
   * <p>
   * For each lag with a negative error, this method first attempts to
   * extrapolate using other errors for the same lag stored in point nabors
   * above or below. This first extrapolation, where feasible, works well when
   * shifts vary slowly with depth.
   * <p> 
   * If this first extrapolation is infeasible, because the number of above
   * and below nabors for some lag is too small, then errors are extrapolated
   * using the errors already computed for other lags. Those errors are
   * already stored in the cells, but are negative, so in this second
   * extrapolation we simply change their sign.
   */
  private static void extrapolateAlignmentErrors(
      int lmin, int lmax, FaultPoint[][] cab) {
    int nab = cab.length;

    // For all arrays of cells linked above-below, ...
    for (int iab=0; iab<nab; ++iab) {
      int mab = cab[iab].length;
      float[][] emp = new float[mab][];

      // Make arrays of errors for all cells in one above-below array.
      for (int jab=0; jab<mab; ++jab) {
        FaultPoint c = cab[iab][jab];
        emp[jab] = c.emp;
      }

      // For each point (each array of errors), ...
      for (int jab=0; jab<mab; ++jab) {

        // For all lags, ...
        for (int lag=lmin,ilag=0; lag<=lmax; ++lag,++ilag) {

          // The error for one point and one lag.
          float empi = emp[jab][ilag];

          // If negative, this error was extrapolating using errors for
          // other lags, so search for an error with the same lag.
          if (empi<0.0f) {

            // If lag is negative, search below; otherwise, search above.
            if (lag<0) {
              for (int kab=jab; kab<mab && empi<0.0f; ++kab)
                empi = emp[kab][ilag];
            } else if (lag>0) {
              for (int kab=jab; kab>=0 && empi<0.0f; --kab)
                empi = emp[kab][ilag];
            }

            // If no good error found, use what we have (but made positive).
            if (empi<0.0f) 
              empi = -emp[jab][ilag];
          }

          // Update the error stored in the point for one lag.
          emp[jab][ilag] = empi;
        } 
      }
    }
  }

  private static float imageValueAt(
    float p1, float p2, float[][]f) {
    int n2 = f.length;
    int n1 = f[0].length;
    int i1 = max(0,min(n1-1,round(p1)));
    int i2 = max(0,min(n2-1,round(p2)));
    return f[i2][i1];
  }


  private static float alignmentError(float f, float g) {
    float fmg = abs(f-g);
    return pow(fmg,2f);
  }



}
