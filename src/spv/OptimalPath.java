/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package spv;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.*;

/**
 * Enhance fault attributes and estimates fault strikes, and dips, 
 * by optimal surface voting.
 *
 * @author Xinming Wu, Univerity of Texas at Austin.
 * @version 2017.07.18
 */
public class OptimalPath {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public OptimalPath(int ru, int rv) {
    _ru = ru;
    _rv = rv;
    _lmin = -ru;
    _lmax =  ru;
    _nl = 1+_lmax-_lmin;
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

  /**
   * Sets bound on fault surface slopes in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
  }

  /**
   * Sets the number of nonlinear smoothings of fault attributes.
   * The default number of smoothings is one.
   * @param esmooth number of nonlinear smoothings.
   */
  public void setAttributeSmoothing(int esmooth) {
    _esmooth = esmooth;
    System.out.println(_esmooth);
  }

    /**
   * Sets extents of smoothing filters used to smooth an extracted fault surface.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   */
  public void setPathSmoothing(double usmooth1) {
    _usmooth1 = usmooth1;
    updateSmoothingFilters();
  }

  public float[][] transpose(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      ft[i1][i2] = fx[i2][i1];
    }}
    return ft;
  }

  /**
   * Extract optimal fault surface from an input fault attribute image.
   * @param fx input array for the fault attribute image.
   */
  public float[] findPath(int c1, int c2, float[][] fx) {
    final int n1 = fx.length;
    final int nl = fx[0].length;
    final float[] u = new float[n1];
    for (int is=0; is<_esmooth; ++is) {
      smoothFaultAttributes(fx,fx);
      setControlPoint(c1,c2,fx);
    }
    float[][] d = new float[n1][nl];
    accumulateForward(fx,d);
    setControlPoint(c1,c2,d);
    backtrackReverse(d,fx,u);
    sub(u,_lmin,u);
    smoothPath(u,u);
    return u;
  }


  public void setControlPoint(int c1, int c2, float[][] fb) {
    int nv = fb.length;
    int nu = fb[0].length;
    for (int iu=0; iu<c1-1; ++iu) {
      int dm = round((c1-iu-2)*0.5f);
      int vb = max(c2-dm,0);
      int ve = min(c2+dm,nv-1);
      for (int iv=vb; iv<=ve; ++iv) {
        fb[iv][iu] = 2f;
      }
    }
    for (int iu=nu-1; iu>c1+1; --iu) {
      int dm = round((iu-c1-2)*0.5f);
      int vb = max(c2-dm,0);
      int ve = min(c2+dm,nv-1);
      for (int iv=vb; iv<=ve; ++iv) {
        fb[iv][iu] = 2f;
      }
    }

  }

 /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothPath(float[] u, float[] us) {
    if (_ref1!=null) {
      _ref1.apply(u,us);
    } 
  }

    /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothFaultAttributes(float[][] fx, float[][] fs) {
    smoothFaultAttributes1(_bstrain1,fx,fs);
    normalizeErrors(fs);
  }

  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }

  public void accumulateForwardX(float[][] e, float[][] d) {
    accumulateX( 1,_bstrain1,e,d);
    //normalizeErrors(d);
  }
  public void accumulateBackwardX(float[][] e, float[][] d) {
    accumulateX(-1,_bstrain1,e,d);
    //normalizeErrors(d);
  }



  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateBackward(float[][] e, float[][] d) {
    accumulate(-1,_bstrain1,e,d);
  }


  /**
   * Returns shifts found by backtracking in reverse.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[] backtrackReverse(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    backtrackReverse(d,e,u);
    return u;
  }

    /**
   * Computes shifts by backtracking in reverse direction.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse(float[][] d, float[][] e, float[] u) {
    backtrack(-1,_bstrain1,_lmin,d,e,u);
  }

  /**
   * Finds shifts by backtracking in accumulated alignment errors.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param lmin minimum lag corresponding to lag index zero.
   * @param d input array[ni][nl] of accumulated errors.
   * @param e input array[ni][nl] of alignment errors.
   * @param u output array[ni] of computed shifts.
   */
  private static void backtrack(
    int dir, int b, int lmin, float[][] d, float[][] e, float[] u) 
  {
    float ob = 1.0f/b;
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?nim1:0;
    int is = (dir>0)?1:-1;
    int ii = ib;
    int il = max(0,min(nlm1,-lmin));
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        il = jl;
        dl = d[ii][jl];
      }
    }
    u[ii] = il+lmin;
    while (ii!=ie) {
      int ji = max(0,min(nim1,ii+is));
      int jb = max(0,min(nim1,ii+is*b));
      int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[jb][ilm1];
      float di = d[ji][il  ];
      float dp = d[jb][ilp1];
      for (int kb=ji; kb!=jb; kb+=is) {
        dm += e[kb][ilm1];
        dp += e[kb][ilp1];
      }
      dl = min3(dm,di,dp);
      if (dl!=di) {
        if (dl==dm) {
          il = ilm1;
        } else {
          il = ilp1;
        }
      }
      ii += is;
      u[ii] = il+lmin;
      if (il==ilm1 || il==ilp1) {
        float du = (u[ii]-u[ii-is])*ob;
        u[ii] = u[ii-is]+du;
        for (int kb=ji; kb!=jb; kb+=is) {
          ii += is;
          u[ii] = u[ii-is]+du;
        }
      }
    }
  }

  /**
   * Smooths fault attributes in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothFaultAttributes1(
    int b, float[][] e, float[][] es) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] ef = new float[n1][nl];
    float[][] er = new float[n1][nl];
    accumulate( 1,b,e,ef);
    accumulate(-1,b,e,er);
    for (int i1=0; i1<n1; ++i1)
      for (int il=0; il<nl; ++il)
        es[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private static void accumulate(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }

  private static void accumulateX(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = max3(dm,di,dp)+e[ii][il];
      }
    }
  }


  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private static float max3(float a, float b, float c) {
    return max(max(a,b),c);
  }


  private void updateSmoothingFilters() {
    _ref1 = (_usmooth1<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth1*_bstrain1);
  }

    /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float emin = e[0][0];
    float emax = e[0][0];
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        float ei = e[i1][il];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float eshift = emin;
    float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        e[i1][il] = (e[i1][il]-eshift)*escale;
      }
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _nl; // number of lags
  private int _lmin,_lmax; // min,max lags
  private int _ru,_rv;
  private int _esmooth = 1; // number of nonlinear smoothings of attributes
  private int _bstrain1 = 4; // inverse of bound on slope in 1st dimension
  private double _usmooth1 = 2.0; // extent of smoothing shifts in 1st dim
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private static final float NO_STRIKE = -0.00001f;
  private static final float NO_DIP    = -0.00001f;
  private SincInterpolator _si;

}
