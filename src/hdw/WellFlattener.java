/****************************************************************************
Copyright 2012, Colorado School of Mines and others.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****************************************************************************/
package hdw;

import java.util.*;
import java.util.Random;
import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;

import util.SmoothWithShaping;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Seismic flattening with dynamic warping.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2017.02.16
 */
public class WellFlattener {
    /**
   * The method used to extrapolate alignment errors.
   * Alignment errors |f[i]-g[i+l]| cannot be computed for indices
   * i and lags l for which the sum i+l is out of bounds. For such
   * indices and lags, errors are missing and must be extrapolated.
   * <p>
   * The extrapolation methods provided are designed to work best 
   * in the case where errors are low for one particular lag l, that 
   * is, when the sequences f and g are related by a constant shift.
   */
  public enum ErrorExtrapolation {
    /**
     * For each lag, extrapolate alignment errors using the nearest
     * error not missing for that lag.
     * <p>
     * This is the default extrapolation method.
     */
    NEAREST,
    /**
     * For each lag, extrapolate alignment errors using the average
     * of all errors not missing for that lag.
     */
    AVERAGE,
    /**
     * For each lag, extrapolate alignment errors using a reflection
     * of nearby errors not missing for that lag.
     */
    REFLECT
  }

    /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public WellFlattener(int shiftMin, int shiftMax) {
    Check.argument(shiftMax-shiftMin>1,"shiftMax-shiftMin>1");
    _lmin = shiftMin;
    _lmax = shiftMax;
    _nl = 1+_lmax-_lmin;
    _extrap = ErrorExtrapolation.NEAREST;
    //_extrap = ErrorExtrapolation.AVERAGE;
    //_extrap = ErrorExtrapolation.REFLECT;
  }

    /**
   * Sets bound on strains in 1st dimension.
   * @param strainMax bound on strain in the 1st dimension.
   */
  public void setStrainMax(double strainMax1) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
  }

  /**
   * Sets the method used to extrapolate alignment errors.
   * Extrapolation is necessary when the sum i+l of sample index
   * i and lag l is out of bounds. The default method is to use 
   * the error computed for the nearest index i and the same lag l.
   * @param ee the error extrapolation method.
   */
  public void setErrorExtrapolation(ErrorExtrapolation ee) {
    _extrap = ee;
  }


    /**
   * Sets the exponent used to compute alignment errors |f-g|^e.
   * The default exponent is 2.
   * @param e the exponent.
   */
  public void setErrorExponent(double e) {
    _epower = (float)e;
  }


  public void setGate(int w, float smin) {
    _w = w;
    _smin = smin;
  }

  public void replaceNulls(float fv, float[][] fx) {
    int nw = fx.length;
    int nz = fx[0].length;
    for (int iw=0; iw<nw; ++iw) {
    for (int iz=0; iz<nz; ++iz) {
      if(fx[iw][iz]==_nullValue)
        fx[iw][iz] = fv;
    }}
  }

  public float[][][] normalizeWells(float[][][] fx) {
    int nc = fx.length;
    int nw = fx[0].length;
    int nz = fx[0][0].length;
    float[][][] gx = new float[nc][nw][nz];
    for (int ic=0; ic<nc; ++ic) {
    for (int iw=0; iw<nw; ++iw) {
    }}
    return gx;
  }

  public float[][] pickTops(int it, Sampling sz, Sampling sw, 
    float[][] us, float[][] ws) {
    int nw = us.length;
    ArrayList<Float> wl = new ArrayList<Float>();
    ArrayList<Float> zl = new ArrayList<Float>();
    wl.add((float)sw.getValue(0));
    zl.add((float)(it*sz.getDelta()+sz.getFirst()));
    for (int iw=1; iw<nw; ++iw) {
      float iz = it+us[iw][it];
      int kz = round(iz);
      if(ws[iw][kz]!=_nullValue) {
        iz = (float)(iz*sz.getDelta()+sz.getFirst());
        zl.add(iz);
        wl.add((float)sw.getValue(iw));
      }
    }
    int np = zl.size();
    float[] wls = new float[np];
    float[] zls = new float[np];
    for (int ip=0; ip<np; ++ip) {
      wls[ip] = wl.get(ip);
      zls[ip] = zl.get(ip);
    }
    return new float[][]{zls,wls};
  }

  public float[][][] flattenX(float[][] gx) {
    int nw = gx.length;
    int nz = gx[0].length;
    int[][] gs = findBounds(gx);
    int[][] mk = zeroint(nz,nw);
    float[][] gp = padNullValues(gs,gx);
    float[][] fx = zerofloat(nz,nw);
    float[][] ft = zerofloat(nz,nw);
    float[][] us = zerofloat(nz,nw);
    ft[0] = gp[0];
    fx[0] = gx[0];
    markNullValues(fx,mk);
    float[] sm = new float[nw];
    sm[0] = 1f;
    for (int iw=1; iw<nw; ++iw) {
      System.out.println("iw="+iw);;
      float[] gw = gp[iw];
      float[][] e = computeErrors(iw,sm,mk,ft,gw);
      int b1 = gs[iw][0];
      int e1 = gs[iw][1];
      setNulls(b1,e1,e);
      normalizeErrors(e);
      float[] u = findShifts(_bstrain1,e);
      us[iw] = u;
      applyShifts(u,gw,ft[iw]);
      applyShifts(u,gx[iw],fx[iw]);
      markNullValues(fx,mk);
      /*
      float sms = 0f;
      float scs = 0f;
      for (int kw=0; kw<iw; kw++){
        if(sm[kw]>_smin) {
          scs += 1f;
          sms += correlate(ft[kw],ft[iw]);
        }
      }
      sm[iw] = sms/scs;
      */
      sm[iw] = correlate(ft[0],ft[iw]);
      System.out.println("sm="+sm[iw]);;
    }
    return new float[][][]{fx,us};
  }

  public float[][] flatten(float[][] gx, float[] cm) {
    int nw = gx.length;
    int nz = gx[0].length;
    int[][] gs = findBounds(gx);
    int[][] mk = zeroint(nz,nw);
    float[][] gp = padNullValues(gs,gx);
    float[][] fx = zerofloat(nz,nw);
    float[][] ft = zerofloat(nz,nw);
    ft[0] = gp[0];
    fx[0] = gx[0];
    markNullValues(fx,mk);
    float[] sm = new float[nw];
    sm[0] = 1f;
    cm[0] = 1f;
    for (int iw=1; iw<nw; ++iw) {
      System.out.println("iw="+iw);;
      float[] gw = gp[iw];
      float[][] e = computeErrors(iw,sm,mk,ft,gw);
      int b1 = gs[iw][0];
      int e1 = gs[iw][1];
      setNulls(b1,e1,e);
      normalizeErrors(e);
      float[] u = findShifts(_bstrain1,e);
      applyShifts(u,gw,ft[iw]);
      applyShifts(u,gx[iw],fx[iw]);
      markNullValues(fx,mk);
      float sms = 0f;
      float scs = 0f;
      for (int kw=0; kw<iw; kw++){
        float cx = correlateX(fx[kw],fx[iw]);
        if (cx>0f) {
          scs += 1f;
          sms += cx;
        }
      }
      cm[iw] = sms/scs;
      sm[iw] = correlate(ft[0],ft[iw]);
      System.out.println("sm="+sm[iw]);;
    }
    return fx;
  }

  public float[] confidence(float[][] fw) {
    int nw = fw.length;
    float[] cw = new float[nw];
    cw[0] = 1;
    for (int iw=1; iw<nw; ++iw) {
      float sms = 0f;
      float scs = 0f;
      for (int kw=0; kw<iw; kw++){
        float cx = correlateX(fw[kw],fw[iw]);
        if (cx>0f) {
          scs += 1f;
          sms += cx;
        }
      }
      cw[iw] = sms/scs;
    }
    return cw;
  }


  public float[][] flatten(float[][] gx) {
    int nw = gx.length;
    int nz = gx[0].length;
    int[][] gs = findBounds(gx);
    int[][] mk = zeroint(nz,nw);
    float[][] gp = padNullValues(gs,gx);
    float[][] fx = zerofloat(nz,nw);
    float[][] ft = zerofloat(nz,nw);
    ft[0] = gp[0];
    fx[0] = gx[0];
    markNullValues(fx,mk);
    float[] sm = new float[nw];
    sm[0] = 1f;
    for (int iw=1; iw<nw; ++iw) {
      System.out.println("iw="+iw);;
      float[] gw = gp[iw];
      float[][] e = computeErrors(iw,sm,mk,ft,gw);
      int b1 = gs[iw][0];
      int e1 = gs[iw][1];
      setNulls(b1,e1,e);
      normalizeErrors(e);
      float[] u = findShifts(_bstrain1,e);
      applyShifts(u,gw,ft[iw]);
      applyShifts(u,gx[iw],fx[iw]);
      markNullValues(fx,mk);
      /*
      float sms = 0f;
      float scs = 0f;
      for (int kw=0; kw<iw; kw++){
        if(sm[kw]>_smin) {
          scs += 1f;
          sms += correlate(ft[kw],ft[iw]);
        }
      }
      sm[iw] = sms/scs;
      */
      sm[iw] = correlate(ft[0],ft[iw]);
      System.out.println("sm="+sm[iw]);;
    }
    return fx;
  }


  // for plots only
  public float[][] flatten(float[][] gx, float[][][] et) {
    int nw = gx.length;
    int nz = gx[0].length;
    int[][] gs = findBounds(gx);
    int[][] mk = zeroint(nz,nw);
    float[][] gp = padNullValues(gs,gx);
    float[][] fx = zerofloat(nz,nw);
    float[][] ft = zerofloat(nz,nw);
    ft[0] = gp[0];
    fx[0] = gx[0];
    markNullValues(fx,mk);
    float[] sm = new float[nw];
    sm[0] = 1f;
    for (int iw=1; iw<nw; ++iw) {
      System.out.println("iw="+iw);;
      float[] gw = gp[iw];
      float[][] e = computeErrors(iw,sm,mk,ft,gw);
      int b1 = gs[iw][0];
      int e1 = gs[iw][1];
      setNulls(b1,e1,e);
      normalizeErrors(e);
      if(iw==15) {
        copy(e,et[0]);
        float[][] d = fillfloat(_nullValue,_nl,nz);
        accumulateForward(_bstrain1,e,d);
        copy(d,et[1]);
      }
      float[] u = findShifts(_bstrain1,e);
      applyShifts(u,gw,ft[iw]);
      applyShifts(u,gx[iw],fx[iw]);
      markNullValues(fx,mk);
      /*
      float sms = 0f;
      float scs = 0f;
      for (int kw=0; kw<iw; kw++){
        if(sm[kw]>_smin) {
          scs += 1f;
          sms += correlate(ft[kw],ft[iw]);
        }
      }
      sm[iw] = sms/scs;
      */
      sm[iw] = correlate(ft[0],ft[iw]);
      System.out.println("sm="+sm[iw]);;
    }
    return fx;
  }


  public float[][][] flatten(float[][][] gx) {
    int nc = gx.length;
    int nw = gx[0].length;
    int nz = gx[0][0].length;
    int[][][] gs = findBounds(gx);
    int[][][] mk = zeroint(nz,nw,nc);
    float[][][] gp = padNullValues(gs,gx);
    float[][][] fx = zerofloat(nz,nw,nc);
    float[][][] ft = zerofloat(nz,nw,nc);
    float[] sm = new float[nw];
    sm[0] = 1f;
    for (int ic=0; ic<nc; ++ic) {
      ft[ic][0] = gp[ic][0];
      fx[ic][0] = gx[ic][0];
    }
    markNullValues(fx,mk);
    for (int iw=1; iw<nw; ++iw) {
      System.out.println("iw="+iw);;
      float[][] gw = new float[nc][nz];
      for (int ic=0; ic<nc; ++ic) {
        gw[ic]  = gp[ic][iw];
      }
      float[][] e = computeErrors(iw,sm,gs,mk,ft,gw);
      float[] u = findShifts(_bstrain1,e);
      for (int ic=0; ic<nc; ++ic) {
        applyShifts(u,gw[ic],ft[ic][iw]);
        applyShifts(u,gx[ic][iw],fx[ic][iw]);
      }
      markNullValues(fx,mk);
    }
    return fx;
  }

  private void applyShifts(float[] u, float[] g, float[] h) {
    int n1 = u.length;
    for (int i1=0; i1<n1; ++i1) {
      float u1 = u[i1];
      if(u1!=_nullValue) {
        int s1 = round(u1);
        int k1 = i1+s1;
        if(k1>=0&&k1<n1)  {
          h[i1] = g[k1];
        }
      }
    }
  }

  private float[] findShifts(int b, float[][] e) {
    int n1 = e.length;
    int nl = e[0].length;
    float[] u = new float[n1];
    float[][] d = fillfloat(_nullValue,nl,n1);
    accumulateForward(b,e,d);
    backtrackReverse(b,_lmin,d,e,u);
    return smooth(u,e);
  }

  private float[] smooth(float[] u, float[][] e) {
    int n1 = u.length;
    float[] w = new float[n1];
    SmoothWithShaping sws = new SmoothWithShaping(20);
    for (int i1=0; i1<n1; ++i1) {
      int kl = round(u[i1]-_lmin);
      float ei = e[i1][kl];
      //TODO assume shifts are non-zero?
      if(ei==_nullValue){
        ei=0f;
      } else {
        ei = 1f-ei;
      }
      w[i1] = ei;
    }
    return sws.smooth(w,u);
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
  public void backtrackReverse(
    int b, int lmin, float[][] d, float[][] e, float[] u) 
  {
    float ob = 1.0f/b;
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = kmaxNotNull(d);
    int ie = kminNotNull(d);
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
      int ji = max(ie,min(nim1,ii-1));
      int jb = max(ie,min(nim1,ii-b));
      int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[jb][ilm1];
      float di = d[ji][il  ];
      float dp = d[jb][ilp1];
      for (int kb=ji; kb!=jb; kb+=-1) {
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
      ii+=-1;
      u[ii] = il+lmin;
      /*
      if (il==ib || il==ie) {
        float du = (u[ii]-u[ii+1])*ob;
        u[ii] = u[ii+1]+du;
        for (int kb=ji; kb!=jb; kb+=-1) {
          ii+=-1;
          u[ii] = u[ii+1]+du;
        }
      }
      */
    }
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private void accumulateForward(int b, float[][] e, float[][] d) {
    int ni = e.length;
    int nl = e[0].length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = kminNotNull(e);
    int ie = kmaxNotNull(e);
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0f;
    for (int ii=ib; ii!=ie; ii++) {
      int ji = max(ib,min(nim1,ii-1));
      int jb = max(ib,min(nim1,ii-b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb--) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private int kminNotNull(float[][] e) {
    int nk = e.length;
    for (int k=0; k<nk; ++k) {
      if (e[k][0]!=_nullValue || e[k][1]!=_nullValue)
        return k;
    }
    return nk;
  }

  private int kmaxNotNull(float[][] e) {
    int nk = e.length;
    for (int k=nk-1; k>=0; --k) {
      if (e[k][0]!=_nullValue || e[k][1]!=_nullValue)
        return k;
    }
    return -1;
  }


  private float[][] computeErrors(int iw, 
    float[] sm, int[][][] gs, int[][][] mk, float[][][] fx, float[][] gw) 
  {
    int nc = fx.length;
    int n1 = fx[0][0].length;
    float[][] sc = zerofloat(_nl,n1);
    float[][] es = fillfloat(_nullValue,_nl,n1);
    for (int ic=0; ic<nc; ++ic) {
      int b1 = gs[ic][iw][0];
      int e1 = gs[ic][iw][1];
      if(b1>=e1) {continue;}
      float[][] ec = computeErrors(iw,sm,mk[ic],fx[ic],gw[ic]);
      setNulls(b1,e1,ec);
      normalizeErrors(ec);
      for (int i1=0; i1< n1; ++i1) {
      for (int il=0; il<_nl; ++il) {
        float esi = es[i1][il];
        float eci = ec[i1][il];
        if(eci!=_nullValue) {
          //eci = 1-eci;
          sc[i1][il] +=1f;
          if(esi==_nullValue) es[i1][il]  = eci;
          else                es[i1][il] += eci;
        }
      }}
    }
    for (int i1=0; i1<n1; ++i1) {
    for (int il=0; il<_nl; ++il) {
      if(es[i1][il]!=_nullValue) 
        es[i1][il]  /= sc[i1][il];
    }}
    normalizeErrors(es);
    return es;
  }

  private void setNulls(int b1, int e1, float[][] e) {
    int n1 = e.length;
    int nl = e[0].length;
    for (int i1=0; i1<b1; i1++)
    for (int il=0; il<nl; il++)
      e[i1][il] = _nullValue;
    for (int i1=e1+1; i1<n1; i1++)
    for (int il=0; il<nl; il++)
      e[i1][il] = _nullValue;
    float emax = max(e);
    for (int i1=b1;i1<=e1; i1++) {
    for (int il=0; il<nl; il++) {
      if(e[i1][il]==_nullValue)
        e[i1][il] = emax;
    }}
  }


  /**
   * Computes alignment errors, not normalized.
   * @param f input array[ni] for sequence f.
   * @param g input array[ni] for sequence g.
   * @param e output array[ni][nl] of alignment errors.
   */
  private float[][] computeErrors(
    final int m2, final float[] sm, final int[][] mk, 
    final float[][] f, final float[] g) {
    final int nl = _nl;
    final int n1 = g.length;
    final int n1m = n1-1;
    final float[][] e = fillfloat(_nullValue,_nl,n1);
    final boolean average = _extrap==ErrorExtrapolation.AVERAGE;
    final boolean nearest = _extrap==ErrorExtrapolation.NEAREST;
    final boolean reflect = _extrap==ErrorExtrapolation.REFLECT;
    final float[] eavg = average?new float[nl]:null; 
    final int[] navg = average?new int[nl]:null;
    float emax = 0.0f;
    // Notes for indexing:
    // 0 <= il < nl, where il is index for lag
    // 0 <= i1 < n1, where i1 is index for sequence f
    // 0 <= j1 < n1, where j1 is index for sequence g
    // j1 = i1+il+lmin, where il+lmin = lag
    // 0 <= i1+il+lmin < n1, so that j1 is in bounds
    // max(0,-lmin-i1) <= il < min(nl,n1-lmin-i1)
    // max(0,-lmin-il) <= i1 < min(n1,n1-lmin-il)
    // j1 = 0    => i1 =     -lmin-il
    // j1 = n1-1 => i1 = n1-1-lmin-il

    // Compute errors where indices are in bounds for both f and g.
    int w = min(_w,m2);
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      int illo = max(0,   -_lmin-i1); // see notes
      int ilhi = min(nl,n1-_lmin-i1); // above
      for (int il=illo,j1=i1+il+_lmin; il<ilhi; ++il,++j1) {
        float es = 0f;
        float ct = 0f;
        for (int i2=0; i2<w; i2++) {
          float gj = g[j1];
          float fi = f[i2][i1];
          if(sm[i2]<_smin) {continue;}
          //if(mk[i2][i1]==1&&i2>10) {continue;}
          //if(mk[i2][i1]==1) {continue;}
          es += error(fi,gj);
          ct += 1f;
        }
        e[i1][il] = es/ct;
      }
    }});
    if (average) {
      for (int i1=0; i1<n1; ++i1) {
        int illo = max(0,   -_lmin-i1); // see notes
        int ilhi = min(nl,n1-_lmin-i1); // above
        for (int il=illo; il<ilhi; ++il) {
          float ei = e[i1][il];
          if(ei!=_nullValue) {
          eavg[il] += ei;
          navg[il] += 1;
          if (ei>emax) 
            emax = ei;
          }
        }
      }
    }


    /*
    for (int i1=0; i1<n1; ++i1) {
      int illo = max(0,   -_lmin-i1); // see notes
      int ilhi = min(nl,n1-_lmin-i1); // above
      for (int il=illo,j1=i1+il+_lmin; il<ilhi; ++il,++j1) {
        float es = 0f;
        float ct = 0f;
        for (int i2=0; i2<w; i2++) {
          float gj = g[j1];
          float fi = f[i2][i1];
          if(sm[i2]<_smin) {continue;}
          if(mk[i2][i1]==1&&i2>10) {continue;}
          es += error(fi,gj);
          ct +=1f;
        }
        if(ct>0f) {es/=ct;e[i1][il] = es;}
        if (average) {
          eavg[il] += es;
          navg[il] += 1;
        }
        if (es>emax) 
          emax = es;
      }
    }
    */

    // If necessary, complete computation of average errors for each lag.
    if (average) {
      for (int il=0; il<nl; ++il) {
        if (navg[il]>0)
          eavg[il] /= navg[il];
      }
    }

    // For indices where errors have not yet been computed, extrapolate.
    for (int i1=0; i1<n1; ++i1) {
      int illo = max(0,   -_lmin-i1); // same as
      int ilhi = min(nl,n1-_lmin-i1); // above
      for (int il=0; il<nl; ++il) {
        if (il<illo || il>=ilhi) {
          if (average) {
            if (navg[il]>0) {
              e[i1][il] = eavg[il];
            } else {
              e[i1][il] = emax;
            }
          } else if (nearest || reflect) {
            int k1 = (il<illo)?-_lmin-il:n1m-_lmin-il;
            if (reflect)
              k1 += k1-i1;
            if (0<=k1 && k1<n1) {
              e[i1][il] = e[k1][il];
            } else {
              e[i1][il] = emax;
            }
          } else {
            e[i1][il] = emax;
          }
        }
      }
    }
    return e;
  }

    /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public void normalizeErrors(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float emin = FLT_MAX;
    float emax = -emin;
    for (int i1=0; i1<n1; ++i1) {
    for (int il=0; il<nl; ++il) {
      float ei = e[i1][il];
      if (ei!=_nullValue) {
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }}
    shiftAndScale(emin,emax,e);
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
      float ei = e[i1][il];
      if (ei!=_nullValue) {
        e[i1][il] = (e[i1][il]-eshift)*escale;
      }
    }}
  }

  private float error(float f, float g) {
    float d = f-g;
    if(_epower==1f) return abs(d);
    return pow(abs(d),_epower);
  }

  private float correlateX(float[] f, float[] g) {
    int n1 = f.length;
    float ff = 0f;
    float gg = 0f;
    float fg = 0f;
    float fa = 0f;
    float ga = 0f;
    float ct = 0f;
    float vt = 0f;
    float fm = FLT_MAX;
    float gm = FLT_MAX;
    for (int i1=0; i1<n1; ++i1) {
      float fi = f[i1];
      float gi = g[i1];
      if (gi!=_nullValue) vt+=1f;
      if(gi!=_nullValue&&fi!=_nullValue) {
        fa += fi;
        ga += gi;
        ct += 1f;
        if(fi<fm) fm = fi;
        if(gi<gm) gm = gi;
      }
    }
    if (ct/vt<0.5f) return -1f;
    fa /= ct;
    ga /= ct;
    for (int i1=0; i1<n1; ++i1) {
      float fi = f[i1];
      float gi = g[i1];
      if(gi!=_nullValue&&fi!=_nullValue) {
        fi -= fa;
        gi -= ga;
        fi -= (fm-fa);
        gi -= (gm-ga);
        fg += fi*gi;
        ff += fi*fi;
        gg += gi*gi;
      }
    }
    return (fg*fg)/(ff*gg);
  }


  private float correlate(float[] f, float[] g) {
    int n1 = f.length;
    float ff = 0f;
    float gg = 0f;
    float fg = 0f;
    for (int i1=0; i1<n1; ++i1) {
      float fi = f[i1];
      float gi = g[i1];
      fg += fi*gi;
      ff += fi*fi;
      gg += gi*gi;
    }
    return fg*fg/(ff*gg);
  }

  private float correlate(Random random, int d, int f1, int g1, 
    float[] f, float[] g) 
  {
    int n1 = f.length;
    float ff = 0f;
    float gg = 0f;
    float fg = 0f;
    for (int k=-d; k<=d; ++k) {
      int fk = f1+k;
      int gk = g1+k;
      if(fk<0||fk>=n1)
        fk = random.nextInt(n1);
      if(gk<0||gk>=n1) 
        gk = random.nextInt(n1);
      float fi = f[fk];
      float gi = g[gk];
      fg += fi*gi;
      ff += fi*fi;
      gg += gi*gi;
    }
    float cx =(fg*fg)/(ff*gg);
    if(Float.isNaN(cx)) cx=0;
    if(Float.isInfinite(cx)) cx=0;
    return cx;
  }

  private float correlate(int d, int f1, int g1, float[] f, float[] g) {
    int n1 = f.length;
    int dfm = f1-max(0,f1-d);
    int dgm = g1-max(0,g1-d);
    int dfp = min(f1+d,n1-1)-f1;
    int dgp = min(g1+d,n1-1)-g1;
    int dm = min(dfm,dgm);
    int dp = min(dfp,dgp);
    float ff = 0f;
    float gg = 0f;
    float fg = 0f;
    for (int k=-dm; k<=dp; ++k) {
      float fi = f[f1+k];
      float gi = g[g1+k];
      fg += fi*gi;
      ff += fi*fi;
      gg += gi*gi;
    }
    float cx =(fg*fg)/(ff*gg);
    if(Float.isNaN(cx)) cx=0;
    if(Float.isInfinite(cx)) cx=0;
    return cx;
  }


  /**
   * Find indexes of the first and last non-null 
   * value for each well log.
   * @param fx array of well logs.
   * @return array of top and bottom indexes.
   */
  private int[] findBounds(float[] fx) {
    int nz = fx.length;
    int[] id = new int[2];
    // Index of first non-null value.
    int iz = 0;
    while(iz<nz&&fx[iz]==_nullValue)
      ++iz;
    // Index of first non-null value.
    int lz = nz-1;
    while(lz>=0&&fx[lz]==_nullValue)
      --lz;
    id[0] = iz;
    id[1] = lz;
    return id;
  }

  private int[][] findBounds(float[][] fx) {
    int nw = fx.length;
    int[][] id = new int[nw][2];
    for (int iw=0; iw<nw; ++iw)
      id[iw] = findBounds(fx[iw]);
    return id;
  }

  private int[][][] findBounds(float[][][] fx) {
    int nc = fx.length;
    int nw = fx[0].length;
    int[][][] id = new int[nc][nw][2];
    for (int ic=0; ic<nc; ++ic)
      id[ic] = findBounds(fx[ic]);
    return id;
  }



  /**
   * Replace null values with a non-null value randomly 
   * selected from the log sequence.
   * @param fx array of well logs.
   * @param id array of top and bottem indexes.
   * @return array of logs with padded non-null values.
   */
  private float[][][] padNullValues(int[][][] id, float[][][] fx) {
    int nc = fx.length;
    int nw = fx[0].length;
    int nz = fx[0][0].length;
    Random random = new Random(314159);
    float[][][] gx = new float[nc][nw][nz];
    for (int ic=0; ic<nc; ++ic) {
    for (int iw=0; iw<nw; ++iw) {
      int bz = id[ic][iw][0];
      int ez = id[ic][iw][1];
      if(bz>=ez) {continue;}
      int[] igood = findGood(fx[ic][iw]);
      for (int iz=0; iz<nz; ++iz) {
        if (fx[ic][iw][iz]==_nullValue) {
          int jz = random.nextInt(igood.length);
          int kz = igood[jz];
          gx[ic][iw][iz] = fx[ic][iw][kz];
        } else {
          gx[ic][iw][iz] = fx[ic][iw][iz];
        }
      }
    }}
    return gx;
  }

  private float[][] padNullValues(int[][] id, float[][] fx) {
    int nw = fx.length;
    int nz = fx[0].length;
    Random random = new Random(314159);
    float[][] gx = new float[nw][nz];
    for (int iw=0; iw<nw; ++iw) {
      int bz = id[iw][0];
      int ez = id[iw][1];
      if(bz>=ez) {continue;}
      int[] igood = findGood(fx[iw]);
      for (int iz=0; iz<nz; ++iz) {
        if (fx[iw][iz]==_nullValue) {
          int jz = random.nextInt(igood.length);
          int kz = igood[jz];
          gx[iw][iz] = fx[iw][kz];
        } else {
          gx[iw][iz] = fx[iw][iz];
        }
      }
    }
    return gx;
  }


  private void markNullValues(float[][] fx, int[][] mk) {
    zero(mk);
    int nw = fx.length;
    int nz = fx[0].length;
    for (int iw=0; iw<nw; ++iw) {
    for (int iz=0; iz<nz; ++iz) {
      if (fx[iw][iz]==_nullValue) {
        mk[iw][iz] = 1;
      } 
    }}
  }

  private void markNullValues(float[][][] fx, int[][][] mk) {
    zero(mk);
    int nc = fx.length;
    int nw = fx[0].length;
    int nz = fx[0][0].length;
    for (int ic=0; ic<nc; ++ic) {
    for (int iw=0; iw<nw; ++iw) {
    for (int iz=0; iz<nz; ++iz) {
      if (fx[ic][iw][iz]==_nullValue) {
        mk[ic][iw][iz] = 1;
      } 
    }}}
  }



  public int[] findGood(float[] f) {
    int ngood = 0;
    int n = f.length;
    int[] igood = new int[n];
    for (int i=0; i<n; ++i) {
      if (f[i]!=_nullValue) {
        igood[ngood] = i;
        ++ngood;
      }
    }
    return copy(ngood,igood);
  }


  private int _nl;
  private int _lmin;
  private int _lmax;
  private int _w=100;
  private float _smin=0.5f;
  private float _nullValue=-999.25f;
  private int _bstrain1;
  private float _epower=1f;
  private ErrorExtrapolation _extrap; // method for error extrapolation

}

    
