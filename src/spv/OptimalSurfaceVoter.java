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
public class OptimalSurfaceVoter {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public OptimalSurfaceVoter(int shiftMin, int shiftMax, int rv, int rw) {
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
    int n3 = ft.length;
    int n2 = ft[0].length;
    int n1 = ft[0][0].length;
    float[][][] fe = new float[n3][n2][n1];
    float[][][] fs = new float[n3][n2][n1];
    FaultCell[] seeds = pickSeeds(d,fm,ft,pt,tt);
    int ns = seeds.length;
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply000(ft,fs);
    fs = sub(fs,min(fs));
    fs = mul(fs,1f/max(fs));
    for (int is=0; is<ns; ++is) {
      if(is%100==0)
        System.out.println("is="+is+"/"+ns);
      FaultCell cell = seeds[is];
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      float tti = cell.getFt();
      float pti = cell.getFp();
      findSurface(i1,i2,i3,_rv,_rw,tti,pti,fs,fe);
    }
    return fe;
  }

  public void findSurface(
    int c1, int c2, int c3, int rv, int rw, 
    float ft, float fp, float[][][] fx, float[][][] fe) {
    int nu = _nl;
    int ru = -_lmin;
    int nv = rv*2+1;
    int nw = rw*2+1;
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    float[][][] fs = new float[nw][nv][nu];
    float[] u = faultNormalVectorFromStrikeAndDip(fp,ft);
    float[] v = faultDipVectorFromStrikeAndDip(fp,ft);
    float[] w = faultStrikeVectorFromStrikeAndDip(fp,ft);
    for (int iw=-rw; iw<=rw; iw++) {
    for (int iv=-rv; iv<=rv; iv++) {
      int um = _lmins[iw+rw][iv+rv];
      int up = _lmaxs[iw+rw][iv+rv];
      for (int iu=um; iu<=up; iu++) {
        float x1 = c1+iu*u[0]+iv*v[0]+iw*w[0];
        float x2 = c2+iu*u[1]+iv*v[1]+iw*w[1];
        float x3 = c3+iu*u[2]+iv*v[2]+iw*w[2];
        fs[iw+rw][iv+rv][iu+ru] = 1-_si.interpolate(s1,s2,s3,fx,x1,x2,x3);
      }
    }}
    float[][] sf = findSurface(fs);
    float fa = 0.0f;
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    ArrayList<Integer> k3s = new ArrayList<Integer>();
    for (int iw=-rw; iw<=rw; ++iw) {
    for (int iv=-rv; iv<=rv; ++iv) {
      float iu = sf[iw+rw][iv+rv];
      float x3 = iu*u[2]+iv*v[2]+iw*w[2]+c3;
      float x2 = iu*u[1]+iv*v[1]+iw*w[1]+c2;
      float x1 = iu*u[0]+iv*v[0]+iw*w[0]+c1;
      int i1 = round(x1);
      int i2 = round(x2);
      int i3 = round(x3);
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
    }}
    int np = k1s.size();
    fa /= np;
    for (int ip=0; ip<np; ++ip) {
      int i1 = k1s.get(ip);
      int i2 = k2s.get(ip);
      int i3 = k3s.get(ip);
      for (int d3=-1;d3<=1;d3++) {
      for (int d2=-1;d2<=1;d2++) {
        int p2 = i2+d2;
        int p3 = i3+d3;
        p2 = max(p2,0);
        p3 = max(p3,0);
        p2 = min(p2,n2-1);
        p3 = min(p3,n3-1);
        fe[p3][p2][i1] += fa;
      }}
    }
    /*
    Sampling sv = new Sampling(nv);
    Sampling sw = new Sampling(nw);
    float[] wvu = buildTrigs(sw,sv,sf);
    int nc = wvu.length;
    float[] xyz = new float[nc];
    for (int ic=0; ic<nc; ic+=3) {
      float iu = wvu[ic+2];
      float iv = wvu[ic+1]-rv;
      float iw = wvu[ic  ]-rw;
      xyz[ic  ] = iu*u[2]+iv*v[2]+iw*w[2]+c3;
      xyz[ic+1] = iu*u[1]+iv*v[1]+iw*w[1]+c2;
      xyz[ic+2] = iu*u[0]+iv*v[0]+iw*w[0]+c1;
    }
    return xyz;
    */
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


  private float[] buildTrigs(Sampling sx, Sampling sy, float[][] z){
    int i = 0;
    int nx = z.length;
    int ny = z[0].length;
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=0;ix<nx-1; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      //if(x0>180f) {continue;}
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny-1; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        if(y0<0f||y0>360f) {continue;}
        float y1 = (float)sy.getValue(iy+1);
        xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
      }
    }
    return copy(i,0,xyz);
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
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] d = du.get();
      if (d==null) du.set(d=new float[n1][nl]);
      accumulateForward(fx[i2],d);
      backtrackReverse(d,fx[i2],uf[i2]);
    }});
    smoothShifts(u,u);
    return u;
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
    //normalizeErrors(fs);
    smoothErrors2(_bstrain2,fs,fs);
    //normalizeErrors(fs);
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
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][] e) {
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
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
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


  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
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
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      smoothErrors1(bf,ef[i2],esf[i2]);
    }});
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
      if(wv>1) {
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
