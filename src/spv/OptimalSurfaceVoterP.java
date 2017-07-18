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
 * @version 2017.07.16
 */
public class OptimalSurfaceVoterP {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public OptimalSurfaceVoterP(int shiftMin, int shiftMax, int rv, int rw) {
    _rv = rv;
    _rw = rw;
    _lmin = shiftMin;
    _lmax = shiftMax;
    _nl = 1+_lmax-_lmin;
    updateShiftRanges();
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

    /**
   * Sets bound on fault surface slopes in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1, double strainMax2) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    Check.argument(strainMax2>0.0,"strainMax2>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
    _bstrain2 = (int)ceil(1.0/strainMax2);
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
   * Sets extents of smoothing filters used to smooth shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   */
  public void setShiftSmoothing(double usmooth1, double usmooth2) {
    _usmooth1 = usmooth1;
    _usmooth2 = usmooth2;
    updateSmoothingFilters();
  }

  public FaultCell[] pickSeeds(
    int d, float fm, float[][][] ft, float[][][] pt, float[][][] tt) {
    int n3 = ft.length;
    int n2 = ft[0].length;
    int n1 = ft[0][0].length;
    ArrayList<FaultCell> cs = new ArrayList<FaultCell>();
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float fti = ft[i3][i2][i1];
      float pti = pt[i3][i2][i1];
      float tti = tt[i3][i2][i1];
      if(fti>fm) {
        FaultCell cell = new FaultCell(i1,i2,i3,fti,pti,tti);
        cs.add(cell);
      }
    }}}
    int np = cs.size();
    int[] is = new int[np];
    float[] fs = new float[np];
    for (int ip=0; ip<np; ++ip) {
      is[ip] = ip;
      fs[ip] = cs.get(ip).getFl();
    }
    quickIndexSort(fs,is);
    int[][][] mark = new int[n3][n2][n1];
    ArrayList<FaultCell> seeds = new ArrayList<FaultCell>();
    for (int ip=np-1; ip>=0; --ip) {
      FaultCell cell = cs.get(is[ip]);
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      int b1 = i1-d; b1=max(b1,0);
      int b2 = i2-d; b2=max(b2,0);
      int b3 = i3-d; b3=max(b3,0);
      int e1 = i1+d; e1=min(e1,n1-1);
      int e2 = i2+d; e2=min(e2,n2-1);
      int e3 = i3+d; e3=min(e3,n3-1);
      boolean ok = true;
      for (int k3=b3;k3<=e3;k3++) {
      for (int k2=b2;k2<=e2;k2++) {
      for (int k1=b1;k1<=e1;k1++) {
        if(mark[k3][k2][k1]==1) {
          ok=false;
          break;
        }
      }}}
      if(ok) {
        seeds.add(cell);
        mark[i3][i2][i1] = 1;
      }
    }
    return seeds.toArray(new FaultCell[0]);
  }

  public float[][][] applyVoting(int d, float fm,
    float[][][] ft, float[][][] pt, float[][][] tt) {
    final int n3 = ft.length;
    final int n2 = ft[0].length;
    final int n1 = ft[0][0].length;
    final FaultCell[] seeds = pickSeeds(d,fm,ft,pt,tt);
    final int ns = seeds.length;
    final int nu = _nl;
    final int ru = -_lmin;
    final int nv = _rv*2+1;
    final int nw = _rw*2+1;
    final int[] ct = new int[1];
    final float[][][] fs = smooth(ft);
    final float[][][] fe = new float[n3][n2][n1];
    Stopwatch sw = new Stopwatch();
    sw.start();
    final Parallel.Unsafe<float[][][]> fsu = 
      new Parallel.Unsafe<float[][][]>();
    final Parallel.Unsafe<float[][][]> feu = 
      new Parallel.Unsafe<float[][][]>();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      float[][][] fst = fsu.get();
      if(fst==null) fsu.set(fst=new float[n3][n2][n1]);
      fst = fs;
      float[][][] fet = feu.get();
      if(fet==null) feu.set(fet=new float[n3][n2][n1]);
      fet = fe;
      float[][] dws = new float[3][nw];
      float[][] dvs = new float[3][nv];
      float[][] dus = new float[3][nu];
      ct[0] += 1;
      if(ct[0]%1000==0)
        System.out.println("done: "+ct[0]+"/"+ns);
      FaultCell cell = seeds[is];
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      float tti = cell.getFt();
      float pti = cell.getFp();
      float[] u = faultNormalVectorFromStrikeAndDip(pti,tti);
      float[] v = faultDipVectorFromStrikeAndDip(pti,tti);
      float[] w = faultStrikeVectorFromStrikeAndDip(pti,tti);
      for (int iw=-_rw; iw<=_rw; ++iw) {
        int kw = iw+_rw;
        dws[0][kw] = iw*w[0];
        dws[1][kw] = iw*w[1];
        dws[2][kw] = iw*w[2];
      }
      for (int iv=-_rv; iv<=_rv; ++iv) {
        int kv = iv+_rv;
        dvs[0][kv] = iv*v[0];
        dvs[1][kv] = iv*v[1];
        dvs[2][kv] = iv*v[2];
      }
      for (int iu=-ru; iu<=ru; ++iu) {
        int ku = iu+ru;
        dus[0][ku] = iu*u[0];
        dus[1][ku] = iu*u[1];
        dus[2][ku] = iu*u[2];
      }
      findSurface(i1,i2,i3,u,dws,dvs,dus,fst,fet);
    }});
    double timeUsed = sw.time();
    System.out.println("time used: "+timeUsed+" seconds");
    return fe;
  }

  public void findSurface(
    int c1, int c2, int c3, float[] u, 
    float[][] dws, float[][] dvs, float[][] dus, 
    float[][][] fx, float[][][] fe) {
    int nu = _nl;
    int ru = -_lmin;
    int nv = dvs[0].length;
    int nw = dws[0].length;
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] fs = fillfloat(1f,nu,nv,nw);
    for (int kw=0; kw<nw; kw++) {
      float dw1 = dws[0][kw]+c1;
      float dw2 = dws[1][kw]+c2;
      float dw3 = dws[2][kw]+c3;
      for (int kv=0; kv<nv; kv++) {
        float dv1 = dw1+dvs[0][kv];
        float dv2 = dw2+dvs[1][kv];
        float dv3 = dw3+dvs[2][kv];
        int um = _lmins[kw][kv];
        int up = _lmaxs[kw][kv];
        for (int ku=um+ru; ku<=up+ru; ku++) {
          int i1 = round(dv1+dus[0][ku]);
          int i2 = round(dv2+dus[1][ku]);
          int i3 = round(dv3+dus[2][ku]);
          i1 = min(max(i1,0),n1-1);
          i2 = min(max(i2,0),n2-1);
          i3 = min(max(i3,0),n3-1);
         fs[kw][kv][ku] = 1-fx[i3][i2][i1];
        }
      }
    }
    float fa = 0.0f;
    float[][] sf = findSurface(fs);
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    ArrayList<Integer> k3s = new ArrayList<Integer>();
    for (int kw=0; kw<nw; ++kw) {
      float dw1 = dws[0][kw]+c1;
      float dw2 = dws[1][kw]+c2;
      float dw3 = dws[2][kw]+c3;
      for (int kv=0; kv<nv; ++kv) {
        float iu = sf[kw][kv];
        int i1 = round(iu*u[0]+dvs[0][kv]+dw1);
        int i2 = round(iu*u[1]+dvs[1][kv]+dw2);
        int i3 = round(iu*u[2]+dvs[2][kv]+dw3);
        boolean inbox = true;
        if(i1<=0||i1>=n1-1) inbox = false;
        if(i2<=0||i2>=n2-1) inbox = false;
        if(i3<=0||i3>=n3-1) inbox = false;
        if(inbox) {
          k1s.add(i1);
          k2s.add(i2);
          k3s.add(i3);
          fa += fx[i3][i2][i1];
        }
      }
    }
    int np = k1s.size();
    fa /= np;
    boolean alignX2 = false;
    if(abs(u[2])>abs(u[1])) alignX2 = true;
    for (int ip=0; ip<np; ++ip) {
      int i1 = k1s.get(ip);
      int i2 = k2s.get(ip);
      int i3 = k3s.get(ip);
      fe[i3][i2][i1] += fa;
      if (alignX2) {
        fe[i3-1][i2][i1] += fa;
        fe[i3+1][i2][i1] += fa;
      } else {
        fe[i3][i2-1][i1] += fa;
        fe[i3][i2+1][i1] += fa;
      }
    }
  }

  public static float[][][][] thin(float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    f = copy(f);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.applyX0X(f,f);
    rgf.applyXX0(f,f);
    float[][][] ff = new float[n3][n2][n1];
    float[][][] pp = new float[n3][n2][n1];
    float[][][] tt = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p00 = p[i3 ][i2 ];
        float[] t00 = t[i3 ][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          float f000 = f00[i1];
          float p000 = p00[i1];
          float t000 = t00[i1];
          if ((                p000<= 22.5f && f0m[i1]<f000 && f0p[i1]<f000) ||
              ( 22.5f<=p000 && p000<= 67.5f && fpm[i1]<f000 && fmp[i1]<f000) ||
              ( 67.5f<=p000 && p000<=112.5f && fp0[i1]<f000 && fm0[i1]<f000) ||
              (112.5f<=p000 && p000<=157.5f && fpp[i1]<f000 && fmm[i1]<f000) ||
              (157.5f<=p000 && p000<=202.5f && f0p[i1]<f000 && f0m[i1]<f000) ||
              (202.5f<=p000 && p000<=247.5f && fmp[i1]<f000 && fpm[i1]<f000) ||
              (247.5f<=p000 && p000<=292.5f && fm0[i1]<f000 && fp0[i1]<f000) ||
              (292.5f<=p000 && p000<=337.5f && fmm[i1]<f000 && fpp[i1]<f000) ||
              (337.5f<=p000                 && f0m[i1]<f000 && f0p[i1]<f000)) {
            ff[i3][i2][i1] = f000;
            pp[i3][i2][i1] = p000;
            tt[i3][i2][i1] = t000;
          } else {
            pp[i3][i2][i1] = NO_STRIKE;
            tt[i3][i2][i1] = NO_DIP;
          }
        }
      }
    }
    float[][][][] flptn = new float[][][][]{ff,pp,tt};
    return flptn;
  }


    /**
   * Returns fault dip vector for specified strike and dip angles.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {u1,u2,u3} of components for dip vector.
   */
  public static float[] faultDipVectorFromStrikeAndDip(
    double phi, double theta) {
    double p = toRadians(phi);
    double t = toRadians(theta);
    double cp = cos(p);
    double sp = sin(p);
    double ct = cos(t);
    double st = sin(t);
    float u1 = (float)( st);
    float u2 = (float)( ct*cp);
    float u3 = (float)(-ct*sp);
    return new float[]{u1,u2,u3};
  }

  /**
   * Returns fault strike vector for specified strike and dip angles.
   * The dip angle theta is not used, but is provided for consistency.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {v1,v2,v3} of components for strike vector.
   */
  public static float[] faultStrikeVectorFromStrikeAndDip(
      double phi, double theta) {
    double p = toRadians(phi);
    double cp = cos(p);
    double sp = sin(p);
    float v1 = 0.0f;
    float v2 = (float)sp;
    float v3 = (float)cp;
    return new float[]{v1,v2,v3};
  }

  /**
   * Returns fault normal vector for specified strike and dip angles.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {w1,w2,w3} of components for normal vector.
   */
  public static float[] faultNormalVectorFromStrikeAndDip(
      double phi, double theta) {
    double p = toRadians(phi);
    double t = toRadians(theta);
    double cp = cos(p);
    double sp = sin(p);
    double ct = cos(t);
    double st = sin(t);
    float w1 = (float)(-ct);
    float w2 = (float)( st*cp);
    float w3 = (float)(-st*sp);
    return new float[]{w1,w2,w3};
  }

  /**
   * Computes shifts for specified images.
   * @param f input array for the image f.
   * @param g input array for the image g.
   * @param u output array of shifts u.
   */
  public float[][] findSurface(float[][][] fx) {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final int nl = fx[0][0].length;
    final float[][] u = new float[n2][n1];
    final float[][] uf = u;
    for (int is=0; is<_esmooth; ++is)
      smoothFaultAttributes(fx,fx);
    float[][] d = new float[n1][nl];
    for (int i2=0; i2<n2; ++i2) {
      accumulateForward(fx[i2],d);
      backtrackReverse(d,fx[i2],uf[i2]);
    }
    smoothShifts(u,u);
    return u;
  }

  private float[][][] smooth(float[][][] ft) {
    int n3 = ft.length;
    int n2 = ft[0].length;
    int n1 = ft[0][0].length;
    float[][][] fs = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply000(ft,fs);
    fs = sub(fs,min(fs));
    fs = mul(fs,1f/max(fs));
    return fs;
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
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothFaultAttributes(float[][][] fx, float[][][] fs) {
    smoothErrors1(_bstrain1,fx,fs);
    smoothErrors2(_bstrain2,fs,fs);
  }

  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
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
        dl = d[ii][jl];
        il = jl;
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
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors1(int b, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int bf = b;
    final float[][][] ef = e;
    final float[][][] esf = es;
    for (int i2=0; i2<n2; ++i2)
      smoothErrors1(bf,ef[i2],esf[i2]);
  }

    /**
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors1(int b, float[][] e, float[][] es) {
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

 /**
   * Smooths alignment errors in 2nd dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 2nd dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors2(int b, float[][][] e, float[][][] es) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    int bf = b;
    float[][] e1  = new float[n2][nl];
    float[][] es1 = new float[n2][nl];
    float[][] ef1 = new float[n2][nl];
    float[][] er1 = new float[n2][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  e[i2][i1];
        es1[i2] = es[i2][i1];
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
    }
  }


  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private void updateSmoothingFilters() {
    _ref1 = (_usmooth1<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth1*_bstrain1);
    _ref2 = (_usmooth2<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth2*_bstrain2);
  }

  private void updateShiftRanges() {
    int nw = _rw*2+1;
    int nv = _rv*2+1;
    _lmins = new int[nw][nv];
    _lmaxs = new int[nw][nv];
    for (int iw=-_rw; iw<=_rw; ++iw) {
    for (int iv=-_rv; iv<=_rv; ++iv) {
      float wv = sqrt(iw*iw+iv*iv);
      if(wv>2) {
        _lmins[iw+_rw][iv+_rv] = max(-round(wv),_lmin);
        _lmaxs[iw+_rw][iv+_rv] = min( round(wv),_lmax);
      }
    }}
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _nl; // number of lags
  private int _lmin,_lmax; // min,max lags
  private int[][] _lmins,_lmaxs;
  private int _rv,_rw;
  private SincInterpolator _si; // for warping with non-integer shifts
  private int _bstrain1 = 1; // inverse of bound on slope in 1st dimension
  private int _bstrain2 = 1; // inverse of bound on slope in 2nd dimension
  private int _esmooth = 1; // number of nonlinear smoothings of errors
  private double _usmooth1 = 0.0; // extent of smoothing shifts in 1st dim
  private double _usmooth2 = 0.0; // extent of smoothing shifts in 2nd dim
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private RecursiveExponentialFilter _ref2; // for smoothing shifts
  private static final float NO_STRIKE = -0.00001f;
  private static final float NO_DIP    = -0.00001f;
}
