package he;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import mef.*;
import util.*;

/**
 * Dynamic programming for refining a picked seismic horizon.
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2016.06.15
 */
public class SurfaceRefinerDp {


  /**
   * Sets bound on strain for all dimensions. Must be in (0,1].
   * The actual bound on strain is 1.0/ceil(1.0/strainMax), which
   * is less than the specified strainMax when 1.0/strainMax is not
   * an integer. The default bound on strain is 1.0 (100%).
   * @param strainMax the bound, a value less than or equal to one.
   */
  public void setStrainMax(double strainMax) {
    //Check.argument(strainMax<=1.0,"strainMax<=1.0");
    //Check.argument(strainMax>0.0,"strainMax>0.0");
    setStrainMax(strainMax,strainMax);
  }

  /**
   * Sets bound on strains in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1, double strainMax2) {
    //Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    //Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    //Check.argument(strainMax1>0.0,"strainMax1>0.0");
    //Check.argument(strainMax2>0.0,"strainMax2>0.0");
    setStrainMax(strainMax1,strainMax2,strainMax2);
  }

  /**
   * Sets bound on strains in 1st, 2nd and 3rd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   * @param strainMax3 bound on strain in the 3rd dimension.
   */
  public void setStrainMax(
    double strainMax1, double strainMax2, double strainMax3) 
  {
    //Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    //Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    //Check.argument(strainMax3<=1.0,"strainMax3<=1.0");
    //Check.argument(strainMax1>0.0,"strainMax1>0.0");
    //Check.argument(strainMax2>0.0,"strainMax2>0.0");
    //Check.argument(strainMax3>0.0,"strainMax3>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
    _bstrain2 = (int)ceil(1.0/strainMax2);
    _bstrain3 = (int)ceil(1.0/strainMax3);
    updateSmoothingFilters();
  }

  public void setControlPoints(int[] cp, int[] cs) {
    _cp = cp;
    _cs = cs;
  }

  /**
   * Sets known points and shifts.
   * @param n1 length of signal f and g
   * @param cp indexes of known control points
   * @param cs shifts of known control points
   */
  public void setControlPoints(int n1, int n2, int strain, int[] cp, int[] cs) {
    _cp = cp;
    _cs = cs;
    _lmins = fillint(0,n1);
    _lmaxs = fillint(n2-1,n1);
    int np = cp.length;
    //int strain = round(1.0f/_bstrain1);
    for (int ip=0; ip<np; ++ip) {
      int ic = cp[ip];
      float lmin = cs[ip];
      float lmax = cs[ip];
      _lmins[ic] = round(lmin);
      _lmaxs[ic] = round(lmax);
      for (int i1=ic+1; i1<n1; ++i1) {
        lmin -= strain;
        lmax += strain;
        int rmin = round(lmin);
        int rmax = round(lmax);
        if (rmin>_lmins[i1]) {_lmins[i1]=rmin;}
        if (rmax<_lmaxs[i1]) {_lmaxs[i1]=rmax;}
      }
      lmin = cs[ip]; lmax = cs[ip];
      for (int i1=ic-1; i1>=0; --i1) {
        lmin -= strain;
        lmax += strain;
        int rmin = round(lmin);
        int rmax = round(lmax);
        if (rmin>_lmins[i1]) {_lmins[i1]=rmin;}
        if (rmax<_lmaxs[i1]) {_lmaxs[i1]=rmax;}
      }
    }
  }

  public void setErrors(float[][] e) {
    int n1 = e.length;
    int nl = e[0].length;
    int np = _cp.length;
    float em = max(e);
    for (int ip=0; ip<np; ++ip) {
      int ic = _cp[ip];
      for (int il=0; il<nl; ++il) {
        e[ic][il] = n1*em;
      }
    }

    for (int ip=0; ip<np; ++ip) {
      int ic = _cp[ip];
      int is = _cs[ip];
      e[ic][is] = min(e);
    }
  }

  /**
   * Sets the number of nonlinear smoothings of alignment errors.
   * In dynamic warping, alignment errors are smoothed the specified 
   * number of times, along all dimensions (in order 1, 2, ...), 
   * before estimating shifts by accumulating and backtracking along 
   * only the 1st dimension. 
   * <p> 
   * The default number of smoothings is zero, which is best for 1D
   * sequences. For 2D and 3D images, two smoothings are recommended.
   * @param esmooth number of nonlinear smoothings.
   */
  public void setErrorSmoothing(int esmooth) {
    _esmooth = esmooth;
  }

  /**
   * Sets extent of smoothing filters used to smooth shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factor. Default 
   * factor is zero, for no smoothing.
   * @param usmooth extent of smoothing filter in all dimensions.
   */
  public void setShiftSmoothing(double usmooth) {
    setShiftSmoothing(usmooth,usmooth);
  }

  /**
   * Sets extents of smoothing filters used to smooth shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   */
  public void setShiftSmoothing(double usmooth1, double usmooth2) {
    setShiftSmoothing(usmooth1,usmooth2,usmooth2);
  }

  /**
   * Sets extents of smoothing filters used to smooth shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   * @param usmooth3 extent of smoothing filter in 3rd dimension.
   */
  public void setShiftSmoothing(
    double usmooth1, double usmooth2, double usmooth3) 
  {
    _usmooth1 = usmooth1;
    _usmooth2 = usmooth2;
    _usmooth3 = usmooth3;
    updateSmoothingFilters();
  }

  public float[][][] getErrorMatrix(
    int m1, float d1, float[][][] gx, float[][] sf) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    float[][][] em = new float[n3][n2][m1];
    int c1 = (m1-1)/2;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int k1=-c1; k1<=c1; ++k1) {
      float x1 = k1*d1+sf[i3][i2];
      em[i3][i2][k1+c1] = si.interpolate(n1,1,0,n2,1,0,n3,1,0,gx,x1,i2,i3);
    }}}
    return em;
  }


  public float[][][][] getErrorMatrix(
    int m1, float d1, FaultSkin[] sks, float[][][] gx, float[][] sf) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    final float[][][] fl = new float[n3][n2][n1];
    for (FaultSkin ski:sks) {
    for (FaultCell fci:ski) {
      int i1 = fci.getI1();
      int i2 = fci.getI2();
      int i3 = fci.getI3();
      fl[i3][i2][i1] = 1f;
    }}
    float[][][] em = new float[n3][n2][m1];
    float[][][] fm = new float[n3][n2][m1];
    int c1 = (m1-1)/2;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int k1=-c1; k1<=c1; ++k1) {
      float x1 = k1*d1+sf[i3][i2];
      em[i3][i2][k1+c1] = si.interpolate(n1,1,0,n2,1,0,n3,1,0,gx,x1,i2,i3);
      fm[i3][i2][k1+c1] = fl[i3][i2][round(x1)];
    }}}

    return new float[][][][] {fm,em};
  }


  /**
   * Finds a path for a given error matrix.
   * @param f input array for the sequence f.
   * @param g input array for the sequence g.
   * @param u output array of shifts u.
   */
  public void findPath(float[][] e, float[] u) {
    int n2 = e.length;
    int n1 = e[0].length;
    if(_cp==null) _lmins = fillint(0,n2);
    if(_cp==null) _lmaxs = fillint(n1-1,n2);
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    float[][] d = accumulateForward(e);
    backtrackReverse(d,e,u);
    smoothShifts(u,u);
  }

    /**
   * Computes shifts for specified images.
   * @param f input array for the image f.
   * @param g input array for the image g.
   * @param u output array of shifts u.
   */
  public void findSurface(float[][][] e, float[][] u) {
    final int n2 = e.length;
    final int n1 = e[0].length;
    final int nl = e[0][0].length;
    final float[][] uf = u;
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] d = du.get();
      if (d==null) du.set(d=new float[n1][nl]);
      accumulateForward(e[i2],d);
      backtrackReverse(d,e[i2],uf[i2]);
    }});
    smoothShifts(u,u);
  }

  public void findSurface(float[][][] fl, float[][][] e, float[][] u) {
    final int n2 = e.length;
    final int n1 = e[0].length;
    final int nl = e[0][0].length;
    final float[][] uf = u;
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(fl,e,e);
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] d = du.get();
      if (d==null) du.set(d=new float[n1][nl]);
      accumulateForward(e[i2],d);
      backtrackReverse(d,e[i2],uf[i2]);
    }});
    smoothShifts(u,u);
  }



  /**
   * Returns smoothed (and normalized) alignment errors.
   * @param e array[n1][nl] of alignment errors.
   * @return array[n1][nl] of smoothed errors.
   */
  public float[][] smoothErrors(float[][] e) {
    float[][] es = like(e);
    smoothErrors(e,es);
    return es;
  }

  /**
   * Returns smoothed (and normalized) alignment errors.
   * @param e array[n2][n1][nl] of alignment errors.
   * @return array[n2][n1][nl] of smoothed errors.
   */
  public float[][][] smoothErrors(float[][][] e) {
    float[][][] es = like(e);
    smoothErrors(e,es);
    return es;
  }

  /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n1][nl] of alignment errors.
   * @param es output array[n1][nl] of smoothed errors.
   */
  public void smoothErrors(float[][] e, float[][] es) {
    smoothErrors1(_bstrain1,e,es);
    normalizeErrors(es);
  }

  /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothErrors(float[][][] e, float[][][] es) {
    smoothErrors1(_bstrain1,e,es);
    normalizeErrors(es);
    smoothErrors2(_bstrain2,es,es);
    normalizeErrors(es);
  }

  public void smoothErrors(float[][][] fl, float[][][] e, float[][][] es) {
    smoothErrors1(_bstrain1,fl,e,es);
    normalizeErrors(es);
    smoothErrors2(_bstrain2,fl,es,es);
    normalizeErrors(es);
  }


  /**
   * Smooths (and normalizes) alignment errors in only the 1st dimension.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothErrors1(float[][][] e, float[][][] es) {
    smoothErrors1(_bstrain1,e,es);
    normalizeErrors(es);
  }

  /**
   * Returns smoothed shifts.
   * @param u array of shifts to be smoothed.
   * @return array of smoothed shifts
   */
  public float[] smoothShifts(float[] u) {
    float[] us = like(u);
    smoothShifts(u,us);
    return us;
  }

  /**
   * Returns smoothed shifts.
   * @param u array of shifts to be smoothed.
   * @return array of smoothed shifts
   */
  public float[][] smoothShifts(float[][] u) {
    float[][] us = like(u);
    smoothShifts(u,us);
    return us;
  }

  /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothShifts(float[] u, float[] us) {
    if (_ref1!=null) {
      _ref1.apply(u,us); 
    } else if (u!=us) {
      copy(u,us);
    }
  }

  /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothShifts(float[][] u, float[][] us) {
    if (_ref1!=null) {
      _ref1.apply1(u,us);
    } else {
      copy(u,us);
    }
    if (_ref2!=null)
      _ref2.apply2(us,us);
  }

  /**
   * Returns errors accumulated in forward direction.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][] accumulateForward(float[][] e) {
    float[][] d = like(e);
    accumulateForward(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][] accumulateReverse(float[][] e) {
    float[][] d = like(e);
    accumulateReverse(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in forward direction in 1st dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateForward1(float[][][] e) {
    float[][][] d = like(e);
    accumulateForward1(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction in 1st dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateReverse1(float[][][] e) {
    float[][][] d = like(e);
    accumulateReverse1(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in forward direction in 2nd dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateForward2(float[][][] e) {
    float[][][] d = like(e);
    accumulateForward2(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction in 2nd dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateReverse2(float[][][] e) {
    float[][][] d = like(e);
    accumulateReverse2(e,d);
    return d;
  }

  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }

  public void accumulateForward(float[][] fl, float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,fl,e,d);
  }


  /**
   * Accumulates alignment errors in reverse direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse(float[][] e, float[][] d) {
    accumulate(-1,_bstrain1,e,d);
  }

  /**
   * Accumulates alignment errors in forward direction in 1st dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward1(float[][][] e, float[][][] d) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateForward(e[i2],d[i2]);
  }

  /**
   * Accumulates alignment errors in reverse direction in 1st dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse1(float[][][] e, float[][][] d) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateReverse(e[i2],d[i2]);
  }

  /**
   * Accumulates alignment errors in forward direction in 2nd dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward2(float[][][] e, float[][][] d) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] di1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        ei1[i2] = e[i2][i1];
        di1[i2] = d[i2][i1];
      }
      accumulate( 1,_bstrain2,ei1,di1);
    }
  }

  /**
   * Accumulates alignment errors in reverse direction in 2nd dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse2(float[][][] e, float[][][] d) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] di1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        ei1[i2] = e[i2][i1];
        di1[i2] = d[i2][i1];
      }
      accumulate(-1,_bstrain2,ei1,di1);
    }
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
   * Returns shifts found by backtracking in reverse in 1st dimension.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[][] backtrackReverse1(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    backtrackReverse1(d,e,u);
    return u;
  }

  /**
   * Returns shifts found by backtracking in reverse in 2nd dimension.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[][] backtrackReverse2(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    backtrackReverse2(d,e,u);
    return u;
  }

  /**
   * Computes shifts by backtracking in reverse direction.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse(float[][] d, float[][] e, float[] u) {
    backtrack(-1,_bstrain1,d,e,u);
  }

  public void backtrackReverse(
    float[][] fl, float[][] d, float[][] e, float[] u) 
  {
    backtrack(-1,_bstrain1,fl, d,e,u);
  }


  /**
   * Computes shifts by backtracking in reverse direction in 1st dimension.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse1(float[][][] d, float[][][] e, float[][] u) {
    int n2 = d.length;
    for (int i2=0; i2<n2; ++i2)
      backtrackReverse(d[i2],e[i2],u[i2]);
  }

  /**
   * Computes shifts by backtracking in reverse direction in 2nd dimension.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse2(float[][][] d, float[][][] e, float[][] u) {
    int n1 = d[0].length;
    int n2 = d.length;
    float[][] di1 = new float[n2][];
    float[][] ei1 = new float[n2][];
    float[] ui1 = new float[n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        di1[i2] = d[i2][i1];
        ei1[i2] = e[i2][i1];
      }
      backtrack(-1,_bstrain2,di1,ei1,ui1);
      for (int i2=0; i2<n2; ++i2)
        u[i2][i1] = ui1[i2];
    }
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public void normalizeErrors(float[][] e) {
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
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public void normalizeErrors(float[][][] e) {
    final float[][][] ef = e;
    int n2 = e.length;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      int nl = ef[i2][0].length;
      int n1 = ef[i2].length;
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          float ei = ef[i2][i1][il];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n1][nl] of errors.
   * @return transposed array[nl][n1] of errors.
   */
  public float[][] transposeLag(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] t = new float[nl][n1];
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        t[il][i1] = e[i1][il];
      }
    }
    return t;
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n2][n1][nl] of errors.
   * @return transposed array[nl][n2][n1] of errors.
   */
  public float[][][] transposeLag(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[nl][n2][n1];
    for (int il=0; il<nl; ++il) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[il][i2][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }

  public float[][] getError(float[][][] e, float[][] u) {
    int n3 = e.length;
    int n2 = e[0].length;
    int n1 = e[0][0].length;
    float[][] eu = new float[n3][n2];
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3)
    for (int i2=0; i2<n2; ++i2)
        eu[i3][i2]=si.interpolate(n1,1,0,e[i3][i2],u[i3][i2]);
    return eu;
  }

  public float[][] smooth(float sig1, float sig2, float[][] w, float[][] u) {
    int n2 = w.length;
    int n1 = w[0].length;
    float[][] b = new float[n2][n1];
    float[][] r = new float[n2][n1];
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(sig1,sig2);
    A2 a2 = new A2(smoother2,w);
    CgSolver cs = new CgSolver(0.001,200);
    mul(w,u,b);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    return r;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _esmooth = 0; // number of nonlinear smoothings of errors
  private double _usmooth1 = 0.0; // extent of smoothing shifts in 1st dim
  private double _usmooth2 = 0.0; // extent of smoothing shifts in 2nd dim
  private double _usmooth3 = 0.0; // extent of smoothing shifts in 3rd dim
  private int _bstrain1 = 1; // inverse of bound on strain in 1st dimension
  private int _bstrain2 = 1; // inverse of bound on strain in 2nd dimension
  private int _bstrain3 = 1; // inverse of bound on strain in 3rd dimension
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private RecursiveExponentialFilter _ref2; // for smoothing shifts
  private RecursiveExponentialFilter _ref3; // for smoothing shifts
  private int[] _cp = null;
  private int[] _cs = null;
  private int[] _lmins = null;
  private int[] _lmaxs = null;
  private int[][] _lmins1 = null;
  private int[][] _lmins2 = null;
  private int[][] _lmaxs1 = null;
  private int[][] _lmaxs2 = null;

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] wp) 
    {
      _s2 = s2;
      _wp = wp;
      float n2 = wp.length;
      float n1 = wp[0].length;
      _sc = sum(wp)/(n1*n2);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      _s2.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s2.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother2 _s2;
    private float[][] _wp;
  }

  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    int n2 = wp.length;
    int n1 = wp[0].length;
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      y[i2][i1] += wp[i2][i1]*x[i2][i1];
  }

  private static void addAndScale(float sc, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      y[i2][i1] += sc*x[i2][i1];
    }}
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(float sigma1, float sigma2) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,x);
    }
    private float _sigma1,_sigma2;
  }

  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  private void updateSmoothingFilters() {
    _ref1 = (_usmooth1<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth1*_bstrain1);
    _ref2 = (_usmooth2<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth2*_bstrain2);
    _ref3 = (_usmooth3<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth3*_bstrain3);
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private void accumulate(int dir, int b, float[][] e, float[][] d) {
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

  private void accumulate(
    int dir, int b, float[][] fl, float[][] e, float[][] d) 
  {
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
          if(fl[kb][ilm1]==1f||fl[kb][ilp1]==1f) {break;}
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }


    /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private void accumulateX(int dir, int b, float[][] e, float[][] d) {
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    int ilb = _lmins[ib];
    int ile = _lmaxs[ib];
    for (int il=ilb; il<=ile; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      int lbi = _lmins[ii];
      int lei = _lmaxs[ii];

      int lbj = _lmins[ji];
      int lej = _lmaxs[ji];

      int lbb = _lmins[jb];
      int leb = _lmaxs[jb];
      for (int il=lbi; il<=lei; ++il) {
        int jl = il;
        if (jl<lbj) {jl=lbj;}
        if (jl>lej) {jl=lej;}
        int ilm1 = il-1; if (ilm1<lbb) ilm1 = lbb;
        int ilp1 = il+1; if (ilp1>leb) ilp1 = leb;
        float dm = d[jb][ilm1];
        float di = d[ji][jl  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          int ilmk = il-1;
          int ilpk = il+1;
          int lbk = _lmins[kb];
          int lek = _lmaxs[kb];
          if (ilmk<lbk) ilmk = lbk;
          if (ilpk>lek) ilpk = lek;
          dm += e[kb][ilmk];
          dp += e[kb][ilpk];
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }

  private void accumulateXX(int dir, int b, float[][] e, float[][] d) {
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    int ilb = _lmins[ib];
    int ile = _lmaxs[ib];
    for (int il=ilb; il<=ile; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int lbi = _lmins[ii];
      int lei = _lmaxs[ii];

      int lbj = _lmins[ji];
      int lej = _lmaxs[ji];

      for (int il=lbi; il<=lei; ++il) {
        int jl = il;
        if (jl<lbj) {jl=lbj;}
        if (jl>lej) {jl=lej;}

        int dn = 40;
        float ds[] = new float[dn*2+1];
        int k = 0;
        for (int dk=-dn; dk<=dn; ++dk) {
          int ilp1 = il+dk; 
          if (ilp1>lej) ilp1 = lej;
          if (ilp1<lbj) ilp1 = lbj;
          ds[k] = d[ji][ilp1]*gauss(40,dk);
          if(dk==0) ds[k] = d[ji][jl];
          k++;
        }
        d[ii][il] = min(ds)+e[ii][il];
      }
    }
  }

  private float gauss(float sig, float dx) {
    return exp(-dx*dx*0.5f/(sig*sig));
  }



  /**
   * Finds shifts by backtracking in accumulated alignment errors.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param d input array[ni][nl] of accumulated errors.
   * @param e input array[ni][nl] of alignment errors.
   * @param u output array[ni] of computed shifts.
   */
  private void backtrack(
    int dir, int b, float[][] d, float[][] e, float[] u) 
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
    int il = 0;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il;
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
      u[ii] = il;
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

  private void backtrack(
    int dir, int b, float[][] fl, float[][] d, float[][] e, float[] u) 
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
    int il = 0;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il;
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
        if(fl[kb][ilm1]==1f||fl[kb][ilp1]==1f) {break;}
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
      u[ii] = il;
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


  private void backtrackX(
    int dir, int b, float[][] d, float[][] e, float[] u) 
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
    //int il = max(0,min(nlm1,-lmin));
    int il = 0;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il;
    while (ii!=ie) {

      int ji = max(0,min(nim1,ii+is));
      int jb = max(0,min(nim1,ii+is*b));
      int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      int lbi = _lmins[ji];
      int lei = _lmaxs[ji];
      int lbb = _lmins[jb];
      int leb = _lmaxs[jb];
      if (il<lbi)  il = lbi;
      if (il>lei)  il = lei;
      if (ilm1<lbb)  ilm1 = lbb;
      if (ilp1>leb)  ilp1 = leb;
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
      u[ii] = il;
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


  private void backtrackXX(
    int dir, int b, float[][] d, float[][] e, float[] u) 
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
    //int il = max(0,min(nlm1,-lmin));
    int il = 0;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il;
    while (ii!=ie) {
      int ji = max(0,min(nim1,ii+is));
      int lbj = _lmins[ji];
      int lej = _lmaxs[ji];
      int lm = il;
      float dm = d[ji][il];
      int dn = 40;
      for (int dk=-dn; dk<=dn; ++dk) {
        int ilp = il+dk; 
        if (ilp>lej) ilp = lej;
        if (ilp<lbj) ilp = lbj;
        if(d[ji][ilp]<dm) {lm = ilp;}
      }

      /*
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
      */
      ii += is;
      u[ii] = lm;
      il = lm;
    }
  }



  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private void shiftAndScale(float emin, float emax, float[][] e) {
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

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private void shiftAndScale(float emin, float emax, float[][][] e) {
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      int nl = ef[i2][0].length;
      int n1 = ef[i2].length;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          ef[i2][i1][il] = (ef[i2][i1][il]-eshift)*escale;
        }
      }
    }});
  }

  /**
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private void smoothErrors1(int b, float[][] e, float[][] es) {
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

  private void smoothErrors1(int b, float[][] fl, float[][] e, float[][] es) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] ef = new float[n1][nl];
    float[][] er = new float[n1][nl];
    accumulate( 1,b,fl,e,ef);
    accumulate(-1,b,fl,e,er);
    for (int i1=0; i1<n1; ++i1)
      for (int il=0; il<nl; ++il)
        es[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
  }


  /**
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private void smoothErrors1(int b, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int bf = b;
    final float[][][] ef = e;
    final float[][][] esf = es;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      smoothErrors1(bf,ef[i2],esf[i2]);
    }});
  }

  private void smoothErrors1(
    int b, float[][][] fl, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int bf = b;
    final float[][][] ef = e;
    final float[][][] esf = es;
    final float[][][] flf = fl;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      smoothErrors1(bf,flf[i2],ef[i2],esf[i2]);
    }});
  }


  /**
   * Smooths alignment errors in 2nd dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 2nd dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private void smoothErrors2(int b, float[][][] e, float[][][] es) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final int bf = b;
    final float[][][]  ef = e;
    final float[][][] esf = es;
    final Parallel.Unsafe<float[][][]> eeu = 
      new Parallel.Unsafe<float[][][]>();
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][][] ee = eeu.get();
      if (ee==null) eeu.set(ee=new float[4][n2][nl]);
      float[][]  e1 = ee[0];
      float[][] es1 = ee[1];
      float[][] ef1 = ee[2];
      float[][] er1 = ee[3];
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  ef[i2][i1];
        es1[i2] = esf[i2][i1];
        for (int il=0; il<nl; ++il) {
          ef1[i2][il] = 0.0f;
          er1[i2][il] = 0.0f;
        }
      }
      accumulate( 1,bf,e1,ef1);
      accumulate(-1,bf,e1,er1);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es1[i2][il] = ef1[i2][il]+er1[i2][il]-e1[i2][il];
        }
      }
    }});
  }

  private void smoothErrors2(
    int b, float[][][] fl, float[][][] e, float[][][] es) 
  {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final int bf = b;
    final float[][][]  ef = e;
    final float[][][] esf = es;
    final float[][][] flf = fl;
    final Parallel.Unsafe<float[][][]> eeu = 
      new Parallel.Unsafe<float[][][]>();
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][][] ee = eeu.get();
      if (ee==null) eeu.set(ee=new float[5][n2][nl]);
      float[][]  e1 = ee[0];
      float[][] es1 = ee[1];
      float[][] ef1 = ee[2];
      float[][] er1 = ee[3];
      float[][] fl1 = ee[4];
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  ef[i2][i1];
        es1[i2] = esf[i2][i1];
        fl1[i2] = flf[i2][i1];
        for (int il=0; il<nl; ++il) {
          ef1[i2][il] = 0.0f;
          er1[i2][il] = 0.0f;
        }
      }
      accumulate( 1,bf,fl1,e1,ef1);
      accumulate(-1,bf,fl1,e1,er1);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es1[i2][il] = ef1[i2][il]+er1[i2][il]-e1[i2][il];
        }
      }
    }});
  }


  private float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private float[] like(float[] a) {
    return new float[a.length];
  }
  private float[][] like(float[][] a) {
    return new float[a.length][a[0].length];
  }
  private float[][][] like(float[][][] a) {
    return new float[a.length][a[0].length][a[0][0].length];
  }

  ///////////////////////////////////////////////////////////////////////////
  // for 3D image warping

  private void normalizeErrors(float[][][][] e) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    MinMax mm = Parallel.reduce(n3,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i3) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int il=0; il<nl; ++il) {
            float ei = ef[i3][i2][i1][il];
            if (ei<emin) emin = ei;
            if (ei>emax) emax = ei;
          }
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }
  private void shiftAndScale(float emin, float emax, float[][][][] e) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int il=0; il<nl; ++il) {
            ef[i3][i2][i1][il] = (ef[i3][i2][i1][il]-eshift)*escale;
          }
        }
      }
    }});
  }
  private void smoothErrors(float[][][][] e) {
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      smoothErrors1(_bstrain1,ef[i3],ef[i3]);
    }});
    normalizeErrors(e);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      smoothErrors2(_bstrain2,ef[i3],ef[i3]);
    }});
    normalizeErrors(e);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] ei2 = new float[n3][][];
      for (int i3=0; i3<n3; ++i3)
        ei2[i3] = ef[i3][i2];
      smoothErrors2(_bstrain3,ei2,ei2);
    }});
    normalizeErrors(e);
  }
  private void computeShifts(float[][][][] e, float[][][] u) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    final float[][][] uf = u;
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] d = du.get();
      if (d==null) du.set(d=new float[n1][nl]);
      for (int i2=0; i2<n2; ++i2) {
        accumulateForward(ef[i3][i2],d);
        backtrackReverse(d,ef[i3][i2],uf[i3][i2]);
      }
    }});
  }
  private void smoothShifts(float[][][] u) {
    if (_ref1!=null) _ref1.apply1(u,u);
    if (_ref2!=null) _ref2.apply2(u,u);
    if (_ref3!=null) _ref3.apply3(u,u);
  }
  private class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }
  private class OverlappingWindows2 {
    public OverlappingWindows2(
      int n1, int n2, int l1, int l2, double f1, double f2) 
    {
      Check.argument(0.0<=f1 && f1<1.0,"0 <= f1 < 1");
      Check.argument(0.0<=f2 && f2<1.0,"0 <= f2 < 1");
      _n1 = n1;
      _n2 = n2;
      _l1 = min(l1,n1);
      _l2 = min(l2,n2);
      _m1 = 1+(int)ceil((_n1-_l1)/(_l1*(1.0-f1)));
      _m2 = 1+(int)ceil((_n2-_l2)/(_l2*(1.0-f2)));
      _s1 = (double)(_n1-_l1)/max(1,_m1-1);
      _s2 = (double)(_n2-_l2)/max(1,_m2-1);
      makeWeights();
      makeScalars();
    }
    public int getN1() { return _n1; }
    public int getN2() { return _n2; }
    public int getL1() { return _l1; }
    public int getL2() { return _l2; }
    public int getM1() { return _m1; }
    public int getM2() { return _m2; }
    public int getI1(int k1) { return (int)(k1*_s1+0.5); }
    public int getI2(int k2) { return (int)(k2*_s2+0.5); }
    public float getWeight(int i1, int i2, int j1, int j2) {
      return _w[j2][j1]*_s[i2+j2][i1+j1];
    }
    public float[][] getWeights() { return _w; }
    public float[][] getScalars() { return _s; }
    private int _n1,_n2; // numbers of samples
    private int _l1,_l2; // window lengths
    private int _m1,_m2; // numbers of windows
    private double _s1,_s2; // nominal window spacings
    private float[][] _w; // weights[l2][l1] for windowing
    private float[][] _s; // scalars[n2][n1] for normalization
    private void makeWeights() {
      _w = new float[_l2][_l1];
      for (int i2=0; i2<_l2; ++i2) {
        for (int i1=0; i1<_l1; ++i1) {
          double s1 = sin((i1+1.0)*PI/(_l1+1.0));
          double s2 = sin((i2+1.0)*PI/(_l2+1.0));
          _w[i2][i1] = (float)(s1*s1*s2*s2);
        }
      }
    }
    private void makeScalars() {
      _s = new float[_n2][_n1];
      for (int k2=0; k2<_m2; ++k2) {
        int i2 = getI2(k2);
        for (int k1=0; k1<_m1; ++k1) {
          int i1 = getI1(k1);
          for (int j2=0; j2<_l2; ++j2) {
            for (int j1=0; j1<_l1; ++j1) {
              _s[i2+j2][i1+j1] += _w[j2][j1];
            }
          }
        }
      }
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _s[i2][i1] = 1.0f/_s[i2][i1];
        }
      }
    }
  }

}
