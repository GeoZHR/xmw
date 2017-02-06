package hdw;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Horizon picking with dynamic programming.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.11.24
 */
import edu.mines.jtk.dsp.*;

public class DynamicFlattening {

  public DynamicFlattening(int k, int dm, int dp) {
    _k = k;
    _dm = dm;
    _dp = dp;
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

  public void setWeights(float c1, float c2, float c3) {
    _c1 = c1;
    _c2 = c2;
    _c3 = c3;
  }

  public float[][] findPeaks(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] fp = new float[n2][n1];
    float[][] f1 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply10(fx,f1);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float fxi = fx[i2][i1  ];
      float f1m = f1[i2][i1-1];
      float f1p = f1[i2][i1+1];
      if(fxi>0 && f1m*f1p<0f) {
        fp[i2][i1] = 1f;
      }
    }}
    return fp;
  }

  public float[][] findThroughs(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] fp = new float[n2][n1];
    float[][] f1 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply10(fx,f1);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float fxi = fx[i2][i1  ];
      float f1m = f1[i2][i1-1];
      float f1p = f1[i2][i1+1];
      if(fxi<0 && f1m*f1p<0f) {
        fp[i2][i1] = 1f;
      }
    }}
    return fp;
  }

  public float[][] findZeroCrossings(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] fp = new float[n2][n1];
    float[][] f1 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply10(fx,f1);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float fxi = fx[i2][i1  ];
      float fxm = fx[i2][i1-1];
      float fxp = fx[i2][i1+1];
      float f1m = f1[i2][i1-1];
      float f1p = f1[i2][i1+1];
      if(fxm*fxp<0 && f1m*f1p>0f&&fxi<0.1f) {
        fp[i2][i1] = 1f;
      }
    }}
    return fp;
  }

  public float[][][] findShifts(
    int d1, int dl, float[][] p, float[][] el, float[][]fx) 
  {
    float[][][] e = computeErrors(d1,dl,fx);
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    float[][][] es = new float[n2][n1][nl];
    smoothErrors2S(dl,p,el,e,es);
    return es;
  }

  public float[][][] computeErrors(int dl, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int nl = dl*2+1;
    float[][][] e = new float[n2][n1][nl];
    for (int i2=0; i2<n2; ++i2)
      computeErrors(fx[0],fx[i2],e[i2]);
    return e;
  }



    /**
   * Computes alignment errors, not normalized.
   * @param f input array[ni] for sequence f.
   * @param g input array[ni] for sequence g.
   * @param e output array[ni][nl] of alignment errors.
   */
  private void computeErrors(float[] f, float[] g, float[][] e) {
    int n1 = f.length;
    int nl = e[0].length;
    int n1m = n1-1;
    int lmin = (nl-1)/2;
    boolean average = false;
    boolean nearest = true;
    boolean reflect = false;
    float[] eavg = average?new float[nl]:null; 
    int[] navg = average?new int[nl]:null;
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
    for (int i1=0; i1<n1; ++i1) {
      int illo = max(0,   -lmin-i1); // see notes
      int ilhi = min(nl,n1-lmin-i1); // above
      for (int il=illo,j1=i1+il+lmin; il<ilhi; ++il,++j1) {
        float ei = abs(f[i1]-g[j1]);
        e[i1][il] = ei;
        if (average) {
          eavg[il] += ei;
          navg[il] += 1;
        }
        if (ei>emax) 
          emax = ei;
      }
    }

    // If necessary, complete computation of average errors for each lag.
    if (average) {
      for (int il=0; il<nl; ++il) {
        if (navg[il]>0)
          eavg[il] /= navg[il];
      }
    }

    // For indices where errors have not yet been computed, extrapolate.
    for (int i1=0; i1<n1; ++i1) {
      int illo = max(0,   -lmin-i1); // same as
      int ilhi = min(nl,n1-lmin-i1); // above
      for (int il=0; il<nl; ++il) {
        if (il<illo || il>=ilhi) {
          if (average) {
            if (navg[il]>0) {
              e[i1][il] = eavg[il];
            } else {
              e[i1][il] = emax;
            }
          } else if (nearest || reflect) {
            int k1 = (il<illo)?-lmin-il:n1m-lmin-il;
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
  }



  public float[][][] computeErrors(float sigma, int dl, float[][] fx) {
    int nl = dl*2+1;
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][][] e = new float[n2][n1][nl];
    LocalCorrelationFilter.Type type = LocalCorrelationFilter.Type.SIMPLE;
    LocalCorrelationFilter.Window wind = LocalCorrelationFilter.Window.GAUSSIAN;
    LocalCorrelationFilter lcf = new LocalCorrelationFilter(type,wind,sigma);
    float[] c = new float[n1];
    float[] fx0 = fx[0];
    lcf.unbias(fx0);
    for (int i2=0; i2<n2; ++i2) {
      lcf.setInputs(fx0,lcf.unbias(fx[i2]));
      for (int il=-dl; il<=dl; il++) {
        lcf.correlate(il,c);
        lcf.normalize(il,c);
        for (int i1=0; i1<n1; ++i1) {
          float ci = c[i1];
          if(Float.isNaN(ci)) {ci=-1;}
          e[i2][i1][il+dl] = ci;
        }
      }
    }
    return e;
  }

  public float[][][] computeRelativeErrors(float sigma, int dl, float[][] fx) {
    int nl = dl*2+1;
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][][] e = new float[n2][n1][nl];
    LocalCorrelationFilter.Type type = LocalCorrelationFilter.Type.SIMPLE;
    LocalCorrelationFilter.Window wind = LocalCorrelationFilter.Window.GAUSSIAN;
    LocalCorrelationFilter lcf = new LocalCorrelationFilter(type,wind,sigma);
    float[] c = new float[n1];
    float[] fx0 = fx[0];
    lcf.unbias(fx0);
    for (int i2=0; i2<n2; ++i2) {
      lcf.setInputs(fx0,lcf.unbias(fx[i2]));
      for (int il=-dl; il<=dl; il++) {
        lcf.correlate(il,c);
        lcf.normalize(il,c);
        for (int i1=0; i1<n1; ++i1) {
          float ci = c[i1];
          if(Float.isNaN(ci)) {ci=-1;}
          e[i2][i1][il+dl] = ci;
        }
      }
      fx0 = lcf.unbias(fx[i2]);
    }
    return e;
  }


  public float[][][] computeErrorsX(
    int d1, int dl, float[][] fx) 
  {
    int nl = dl*2+1;
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[] fr = new float[d1*2+1];
    float[] fs = new float[d1*2+1];
    float[][][] e = new float[n2][n1][nl];
    float[] fx0 = fx[0];
    for (int i2=0; i2<n2; ++i2) {
      float[] fx2 = fx[i2];
      float[][] e2 = e[i2];
    for (int i1=0; i1<n1; ++i1) {
      int i1b = i1-d1; i1b = max(i1b,0);
      int i1e = i1+d1; i1e = min(i1e,n1-1);
      float frs = 0f;
      for (int k1=i1b; k1<=i1e; ++k1) {
        frs += fx0[k1];
        fr[d1+k1-i1] = fx0[k1];
      }
      float fsig = 0.0f;
      float fra = frs/(i1e-i1b+1f);
      for (int k1=0; k1<d1*2+1; k1++) {
        float fri = fr[k1]-fra;
        fsig += fri*fri;
        fr[k1] = fri;
      }
      fsig = sqrt(fsig);
      int ilb = -dl; if(i1<dl) ilb=-i1;
      int ile =  dl; if((n1-1-i1)<dl) i1e=n1-1-i1;
      for (int il=ilb; il<=ile; ++il) {
        i1b = il+i1-d1; i1b = max(i1b,0);
        i1e = il+i1+d1; i1e = min(i1e,n1-1);
        for (int k1=i1b; k1<=i1e; ++k1)
          fs[d1+k1-i1-il] = fx2[k1];
        e2[i1][il+dl] = crossCorrelation(fsig,fr,fs);
      }
    }}
    return e;
  }

  public void smooth1(float strain, float[][][] e, float[][][] es) {
    smoothErrors1(strain,e,es);
    normalizeErrors(es);
  }


  public float[][][] smooth1(float strain, float[][][] e) {
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    float[][][] es = new float[n2][n1][nl];
    smoothErrors1(strain,e,es);
    normalizeErrors(es);
    return es;
  }

  public float[][] pickX(final int lmax, 
    float[][] p, float[][] el, float[][][] e, float[][][] d) {
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    float[][] u = new float[n1][n2];
    for (int i1=0; i1<n1; ++i1) {
      float[][] e1 = new float[n2][nl]; 
      float[][] ef1 = fillfloat(_emin,nl,n2);//new float[n2][nl]; 
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  e[i2][i1];
      }
      accumulateForward2(i1,lmax,p,el,e1,ef1);
      for (int i2=0; i2<n2; ++i2) {
      for (int il=0; il<nl; ++il) {
        d[i2][i1][il] = ef1[i2][il];
      }}
      backtrackK(i1,lmax,p,ef1,u[i1]);
    }
    return u;
  }


  public float[][] pick(float strain, final int lmax, 
    float[][] p, float[][] el, float[][][] e, float[][][] d) {
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    float[][] u = new float[n1][n2];
    /*
    float[][][] es2 = new float[n2][n1][nl];
    //float[][][] es1 = smooth1(strain,e);
    for (int i1=0; i1<n1; ++i1) {
      float[][] e1 = new float[n2][nl]; 
      float[][] ef1 = fillfloat(_emin,nl,n2);//new float[n2][nl]; 
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  e[i2][i1];
      }
      accumulateForward2X(i1,lmax,p,el,e1,ef1);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es2[i2][i1][il] = ef1[i2][il];
        }
      }
    }
    float[][][] es1 = smooth1(strain,es2);
    System.out.println("esmin="+min(es1));
    System.out.println("esmax="+max(es1));
    */
    float[][][] es1 = e;
    for (int i1=0; i1<n1; ++i1) {
      float[][] e1 = new float[n2][nl]; 
      float[][] ef1 = fillfloat(_emin,nl,n2);//new float[n2][nl]; 
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  es1[i2][i1];
      }
      accumulateForward2X(i1,lmax,p,el,e1,ef1);
      backtrackK(i1,lmax,p,ef1,u[i1]);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          d[i2][i1][il] = ef1[i2][il];
        }
      }
    }
    return u;
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

  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
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




  private void smoothErrors1(float strain1, float[][][] e, float[][][] d) {
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    int bstrain1 = (int)ceil(1.0/strain1);
    for (int i2=0; i2<n2; ++i2) {
      float[][] ef = new float[n1][nl];
      float[][] er = new float[n1][nl];
      accumulate( 1,bstrain1,e[i2],ef);
      accumulate(-1,bstrain1,e[i2],er);
      for (int i1=0; i1<n1; ++i1)
        for (int il=0; il<nl; ++il)
          d[i2][i1][il] = ef[i1][il]+er[i1][il]-e[i2][i1][il];
    }
  }

    /**
   * Computes shifts for specified sequences.
   * @param f input array for the sequence f.
   * @param g input array for the sequence g.
   * @param u output array of shifts u.
   */
  public void findShifts(float strain, int lmax, float[][][]e, float[][] u) {
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    float[][][] d = new float[n2][n1][nl];
    float[][][] es = new float[n2][n1][nl];
    for (int is=0; is<2; ++is)
      smooth1(strain,e,es);
    int b = (int)ceil(1.0/strain);
    for (int i2=0; i2<n2; ++i2) {
      accumulate(1,b,es[i2],d[i2]);
      backtrackX(-1,b,-lmax,d[i2],es[i2],u[i2]);
    }
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
        d[ii][il] = max3(dm,di,dp)+e[ii][il];
      }
    }
  }

  private static float max3(float a, float b, float c) {
    return b>=a?(b>=c?b:c):(a>=c?a:c); // if equal, choose b
  }



  public void smoothErrors2S(final int lmax, 
    float[][] p, float[][] el, float[][][] e, float[][][] es) {
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    for (int i1=0; i1<n1; ++i1) {
      float[][] e1 = new float[n2][nl]; 
      float[][] es1 = new float[n2][nl]; 
      float[][] ef1 = fillfloat(_emin,nl,n2);//new float[n2][nl]; 
      //float[][] er1 = new float[n2][nl]; 
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  e[i2][i1];
        es1[i2] =  es[i2][i1];
      }
      accumulateForward2(i1,lmax,p,el,e1,ef1);
      //accumulateBackward2(i1,lmax,p,e1,er1);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es1[i2][il] = ef1[i2][il];//+er1[i2][il]-e1[i2][il];
        }
      }
    }
  }

  private void backtrackK(
    int i1, int lmax, float[][] p, float[][] d, float[] u) 
  {
    int n2 = p.length;
    int n1 = p[0].length;
    int nl = d[0].length;
    int nlm = nl-1;
    int ilp = 0;
    float dmax = d[n2-1][0];
    for (int il=0; il<nl; ++il) {
      float di = d[n2-1][il];
      if (di>dmax) {ilp = il; dmax = di;}
    }
    u[n2-1] = ilp;
    for (int i2=n2-2; i2>=0; i2--) {
      float[] d2 = d[i2];
      int ih = i1+ilp-lmax;
      if(ih<0||ih>n1-1){continue;}
      int ilm = ilp-round(p[i2+1][ih]);
      int klm = ilm+_dm;
      int klp = ilm+_dp;
      if(klm<0){klm=0;} if(klm>nlm) {klm=nlm;}
      if(klp<0){klp=0;} if(klp>nlm) {klp=nlm;}
      ilp = klm; dmax = d2[klm];
      for (int kl=klm+1; kl<=klp; ++kl) {
        float dk = d[i2][kl];
        if (dk>dmax) {dmax=dk;ilp = kl;}
      }
      u[i2] = ilp;
    }
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
  private static void backtrackX(
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
      dl = max3(dm,di,dp);
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


  public void smoothErrors2(final int lmax, 
    float[][] p, float[][] el, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int n1 = e[0].length;
    final int nl = e[0][0].length;
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
      accumulateForward2(i1,lmax,p,el,e1,ef1);
      //accumulateBackward2(i1,lmax,p,el,e1,er1);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es1[i2][il] = ef1[i2][il]+er1[i2][il]-e1[i2][il];
        }
      }
    }});
  }

  public void accumulateForward2(
    int i1, int lmax, float[][] p, float[][] el, float[][] e, float[][] d) {
    int n2 = p.length;
    int n1 = p[0].length;
    int ilb = 0;
    int ile = 2*lmax;
    int nl = 2*lmax+1;
    d[0][lmax] = 1f;
    int[][] lsk = new int[n2][nl];
    for (int il=0; il<nl; il++)
      lsk[0][il] = il;
    for (int i2=1; i2<n2; i2++) {
      int i2m = i2-1;
      float[] p2 = p[i2];
      float[] d2 = d[i2m];
      for (int il=ilb; il<=ile; ++il) {
        int ih = i1+il-lmax;
        if(ih<0) {continue;}
        if(ih>=n1) {continue;}
        float eli = el[i2][ih];
        //if(e[i2][il]<0f&&eli>0.9f) {continue;}
        float pi = p2[ih];
        int ilm = il-round(pi);
        int klm = ilm+_dm; 
        int klp = ilm+_dp; 
        if(klm<ilb){klm=ilb;} if(klm>ile) {klm=ile;}
        if(klp<ilb){klp=ilb;} if(klp>ile) {klp=ile;}
        int klt = klm; 
        float dmax = d2[klm];
        for (int kl = klm; kl<=klp; ++kl) {
          float dk = d[i2m][kl];
          if (dk>dmax) {dmax=dk;klt=kl;}
        }
        lsk[i2][il] = klt;
        int lm = il;
        int ke = min(i2,5);
        float e1 = 0;
        for (int k=0; k<=ke; k++) {
          lm = lsk[i2-k][lm];
          if(k<1) {continue;}
          int hm = i1+lm-lmax;
          hm = max(hm,0);
          hm = min(hm,n1-1);
          float plm = p[i2-k][hm];
          float dlm = il-lm;
          float spm = sqrt(1+plm*plm);
          float slm = sqrt((1+k)*(1+k)+dlm*dlm);
          e1 += (1+k+dlm*pi)/(spm*slm);
        }
        e1 /= ke;
        eli = pow(eli,8);
        float e3 = e[i2][il];
        float dt = dmax+_c1*e1*eli+_c3*e3;
        if(dmax!=_emin)
          d[i2][il] = dt;
      }
    }
  }

  public void accumulateForward2X(
    int i1, int lmax, float[][] p, float[][] el, float[][] e, float[][] d) {
    int n2 = p.length;
    int n1 = p[0].length;
    int ilb = 0;
    int ile = 2*lmax;
    int nl = 2*lmax+1;
    d[0][lmax] = 1f;
    int[][] lsk = new int[n2][nl];
    for (int il=0; il<nl; il++)
      lsk[0][il] = il;
    for (int i2=1; i2<n2; i2++) {
      int i2m = i2-1;
      float[] p2 = p[i2];
      float[] d2 = d[i2m];
      for (int il=ilb; il<=ile; ++il) {
        int ih = i1+il-lmax;
        if(ih<0) {continue;}
        if(ih>=n1) {continue;}
        float eli = el[i2][ih];
        //if(e[i2][il]<0f&&eli>0.9f) {continue;}
        float pi = p2[ih];
        int ilm = il-round(pi);
        int klm = ilm+_dm; 
        int klp = ilm+_dp; 
        if(klm<ilb){klm=ilb;} if(klm>ile) {klm=ile;}
        if(klp<ilb){klp=ilb;} if(klp>ile) {klp=ile;}
        int klt = klm; 
        float dmax = d2[klm];
        for (int kl = klm; kl<=klp; ++kl) {
          float dk = d[i2m][kl];
          if (dk>dmax) {dmax=dk;klt=kl;}
        }
        lsk[i2][il] = klt;
        int lm = il;
        int ke = min(i2,5);
        float e1 = 0;
        float pc = 0;
        for (int k=0; k<=ke; k++) {
          lm = lsk[i2-k][lm];
          if(k<1) {continue;}
          int hm = i1+lm-lmax;
          hm = max(hm,0);
          hm = min(hm,n1-1);
          pc -= p[i2-k][hm];
          float dlm = il-lm;
          float x2s = (1+k)*(1+k);
          float spm = sqrt(x2s+pc*pc);
          float slm = sqrt(x2s+dlm*dlm);
          e1 += (x2s+dlm*pc)/(spm*slm);
        }
        eli = pow(eli,8);
        float e3 = e[i2][il];
        float dt = dmax+_c1*e1*eli+_c3*e3;
        if(dmax!=_emin)
          d[i2][il] = dt;
      }
    }
  }



  public void accumulateBackward2(
    int i1, int lmax, float[][] p, float[][] e, float[][] d) 
  {
    int n2 = p.length;
    int n1 = p[0].length;
    int ilb = 0;
    int ile = 2*lmax;
    for (int i2=n2-2; i2>=0; i2--) {
      int i2p = i2+1;
      float[] p2 = p[i2 ];
      float[] d2 = d[i2p];
      for (int il=ilb; il<=ile; ++il) {
        int ih = i1+il-lmax;
        if(ih<0) {continue;}
        if(ih>=n1) {continue;}
        float pi = p2[ih];
        int ilm = il+round(pi);
        int klm = ilm+_dm;
        int klp = ilm+_dp; 
        if(klm<ilb){klm=ilb;} if(klm>ile) {klm=ile;}
        if(klp<ilb){klp=ilb;} if(klp>ile) {klp=ile;}
        int klt = klm; float dmax = d2[klm];
        for (int kl = klm; kl<=klp; ++kl) {
          float dk = d[i2p][kl];
          if (dmax>dk) {dmax=dk;klt=kl;}
        }
        float dl = klt-il;
        float sl = sqrt(1+dl*dl);
        float sp = sqrt(1+pi*pi);
        float e1 = (1f+dl*pi)/(sp*sl);
        float e3 = e[i2][il];
        float dt = dmax+_c1*e1+_c3*e3;
        d[i2][il] += dt;
      }
    }
  }

  public float[] copyInWind(int c1, float[] g) {
    int n1 = g.length;
    int i1b = max(c1-_wind,0);
    int i1e = min(c1+_wind,n1-1);
    float[] gc = new float[2*_wind+1];
    for (int i1=i1b; i1<=i1e; ++i1)
      gc[i1+c1] = g[i1];
    return gc;
  }

  public float crossCorrelation(float fsig, float[] f, float[] g) {
    int n1 = g.length;
    float gs = 0.0f;
    for (int i1=0; i1<n1; ++i1)
      gs += g[i1];
    float ga = gs/n1;
    float gsig = 0.0f;
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i1]-ga;
      gsig += gi*gi;
      g[i1] = gi;
    }
    gsig = sqrt(gsig);
    float fgc = 0.0f;
    for (int i1=0; i1<n1; ++i1)
      fgc += f[i1]*g[i1];
    return  fgc/(fsig*gsig);
  }



  public float[] pickX(
    int[] k1, int[] k2, float[][] p, float[][] e) {
    _emax = sum(abs(e));
    int n2 = p.length;
    int n1 = p[0].length;
    float[] u = new float[n2];
    float[][] d = fillfloat(100,n1,n2);
    //setControlPoints(k1,k2,di);
    int ip2 = k2[0];
    int ip1 = k1[0];
    accumulateForward(ip1,ip2,p,e,d);
    track(d,u);
    //backtrackReverse(pi,di,u);
    //accumulateBackward(i1p,i2p,pi,eli,ei,di);
    //backtrackReverse(pi,di,u);
    return u;
  }


  public float[] pick(
    float[] k1, float[] k2, float[][] p, float[][] e) {
    _emax = sum(abs(e));
    final float[][] ei = mapResample(e);
    final float[][] pi = slopeResample(p);
    int n2 = pi.length;
    int n1 = pi[0].length;
    final float[][] di = fillfloat(_emax,n1,n2);
    final int np = k2.length;
    final int[] id = rampint(0,1,np);
    quickIndexSort(k2,id);
    final float[][] u = new float[2][n2];
    for (int ip=0; ip<np; ++ip) {
      int k = id[ip];
      int i2p = round(k2[k]);
      int i1p = round(k1[k]*_k);
      System.out.println("k="+k);
      System.out.println("i1p="+i1p);
      System.out.println("i2p="+i2p);
      int i1l = -1;
      int i1r = -1;
      int i2l = 0;
      int i2r = n2-1;
      if(ip>0&&ip<np-1) {
        int kl = id[ip-1];
        int kr = id[ip+1];
        i2l=round(k2[kl]);
        i2r=round(k2[kr]);
        i1l=round(k1[kl]*_k);
        i1r=round(k1[kr]*_k);
      }else if(ip==np-1 && np>1) {
        int kl = id[ip-1];
        i2l=round(k2[kl]);
        i1l=round(k1[kl]*_k);
      }else if(ip==0 && np>1) {
        int kr = id[ip+1];
        i2r=round(k2[kr]);
        i1r=round(k1[kr]*_k);
      }
      accumulateForward(i1p,i2p,i2r,pi,ei,di);
      backtrack(i1r,i2r,i2p,pi,di,u);
      accumulateBackward(i1p,i2p,i2l,pi,ei,di);
      backtrackReverse(i1l,i2l,i2p,pi,di,u);
    }
    return mul(u[0],1.f/_k);
  }

  public void pickForward(int i1b, int i1e, int i2b, int i2e, 
    float[][] p, float[][] e, float[][] u) {


  }

  public float[] pickForward(
    int[] k1, int[] k2, float[][] p, float[][] el, float[][] e) {
    _emax = sum(abs(e));
    float[][] ei = mapResample(e);
    float[][] eli = mapResample(el);
    float[][] pi = slopeResample(p);
    int n2 = pi.length;
    int n1 = pi[0].length;
    float[] u = new float[n2];
    float[][] di = new float[n2][n1];
    int i2p = k2[0];
    int i1p = k1[0]*_k;
    accumulateForward(i1p,i2p,pi,eli,ei,di);
    backtrack(pi,di,u);
    //accumulateBackward(i1p,i2p,pi,eli,ei,di);
    //backtrackReverse(pi,di,u);
    return mul(u,1.f/_k);
  }

  public float[] pickBackward(
    int[] k1, int[] k2, float[][] p, float[][] el, float[][] e) {
    _emax = sum(abs(e));
    float[][] ei = mapResample(e);
    float[][] eli = mapResample(el);
    float[][] pi = slopeResample(p);
    int n2 = pi.length;
    int n1 = pi[0].length;
    float[] u = new float[n2];
    float[][] di = new float[n2][n1];
    int i2p = k2[0];
    int i1p = k1[0]*_k;
    accumulateBackward(i1p,i2p,pi,eli,ei,di);
    backtrackReverse(pi,di,u);
    return mul(u,1.f/_k);
  }


  public float[] trackForward(int k1, int k2, float[][] p) {
    float[][] pi = slopeResample(p);
    int n2 = pi.length;
    int n1 = pi[0].length;
    float[] u = new float[n2];
    u[0] = k1*_k;
    for (int i2=1; i2<n2; ++i2) {
      int i1 = round(u[i2-1]);
      i1 = max(i1,0);
      i1 = min(i1,n1-1);
      u[i2] = u[i2-1]+pi[i2-1][i1];
    }
    return mul(u,1.0f/_k);
  }

  public float[] trackBackward(int k1, int k2, float[][] p) {
    float[][] pi = slopeResample(p);
    int n2 = pi.length;
    int n1 = pi[0].length;
    float[] u = new float[n2];
    u[n2-1] = k1*_k;
    for (int i2=n2-2; i2>=0; i2--) {
      int i2p = i2+1;
      int i1 = round(u[i2p]);
      i1 = max(i1,0);
      i1 = min(i1,n1-1);
      u[i2] = u[i2p]-pi[i2p][i1];
    }
    return mul(u,1.0f/_k);

  }


  public float[][] mapResample(float[][] mp) {
    int n2 = mp.length;
    int n1 = mp[0].length;
    int m1 = (n1-1)*_k+1;
    double dm = 1.0/_k;
    float[][] mpi = new float[n2][m1];
    for (int i2=0; i2<n2; ++i2) {
      _si.interpolate(n1,1,0,mp[i2],m1,dm,0,mpi[i2]);
    }
    return mpi;
  }

  public void setControlPoints(int[] k1, int[] k2, float[][] d) {
    zero(d);
    int np = k2.length;
    int n1 = d[0].length;
    for (int ip=0; ip<np; ++ip) {
      int p2 = k2[ip];
      int p1 = k1[ip]*_k;
      for (int i1=0; i1<n1; ++i1) {
        d[p2][i1] = _emax;
      }
      d[p2][p1] = 0f;
    }
  }

  public float[][] slopeResample(float[][] p) {
    int n2 = p.length;
    int n1 = p[0].length;
    int m1 = (n1-1)*_k+1;
    double dm = 1.0/_k;
    float[][] pi = new float[n2][m1];
    float[] pt = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] pi2 = pi[i2];
      float[] p2 = p[i2];
      for (int i1=0; i1<n1; ++i1) 
        pt[i1] = _k*p2[i1];
      _si.interpolate(n1,1,0,pt,m1,dm,0,pi2);
    }
    return pi;
  }

  public void accumulateForward(
    float[][] p, float[][] el, float[][] e, float[][] d) 
  {
    int n2 = e.length;
    int n1 = e[0].length;
    int n1m = n1-1;
    int[] i1t = new int[n1];
    for (int i2=1; i2<n2; i2++) {
      int i2m = i2-1;
      float[] p2 = p[i2];
      float[] d2 = d[i2m];
      int[] kts = new int[n1];
      for (int i1=0; i1<n1; ++i1) {
        float p2i = p2[i1];
        int i1m = i1-round(p2i);
        int k1m = i1m+_dm; 
        int k1p = i1m+_dp; 
        if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
        if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
        int k1t = k1m; float dmin = d2[k1m];
        for (int k1 = k1m+1; k1<=k1p; ++k1) {
          float dk = d[i2m][k1];
          if (dk<dmin) {dmin=dk;k1t=k1;}
        }
        kts[i1] = k1t;
        float e1 = i1-k1t-p2i;
        float e2 = i1-2f*k1t+i1t[k1t];
        float e3 = e[i2][i1];
        e3 = e[i2][i1];
        e1 *=e1;
        e2 *=e2;
        d[i2][i1] += dmin+_c1*e1+_c2*e2+_c3*e3;
      }
      copy(kts,i1t);
    }
  }

  public void accumulateForward(int i1b, int i2b, int i2e, 
    float[][] p, float[][] e, float[][] d) 
  {
    int n1 = e[0].length;
    int n1m = n1-1;
    int[] i1t = new int[n1];
    d[i2b][i1b] = 0f;
    int i1p = i1b;
    for (int i2=i2b+1; i2<=i2e; i2++) {
      int i2m = i2-1;
      float[] p2 = p[i2];
      float[] d2 = d[i2m];
      int[] kts = new int[n1];
      int k1b = i1p+round(p2[i1p])-200;
      int k1e = k1b+400;
      k1b = max(k1b,0);
      k1e = min(k1e,n1-1);
      float dmk = FLT_MAX;
      //for (int i1=k1b; i1<=k1e; ++i1) {
      for (int i1=0; i1<=n1m; ++i1) {
        float p2i = p2[i1];
        int i1m = i1-round(p2i);
        int k1m = i1m+_dm; 
        int k1p = i1m+_dp; 
        if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
        if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
        int k1t = k1m; float dmin = d2[k1m];
        for (int k1 = k1m+1; k1<=k1p; ++k1) {
          float dk = d[i2m][k1];
          if (dk<dmin) {dmin=dk;k1t=k1;}
        }
        kts[i1] = k1t;
        float e1 = i1-k1t-p2i;
        float e2 = i1-2f*k1t+i1t[k1t];
        float e3 = e[i2][i1];
        e3 = e[i2][i1];
        e1 =abs(e1);
        e2 *=e2;
        float dt = dmin+_c1*e1+_c2*e2+_c3*e3;
        float di = d[i2][i1]; if(di!=_emax) {dt+=di;}
        d[i2][i1] = dt;
        if(dt<dmk) {dmk=dt;i1p=i1;}
      }
      copy(kts,i1t);
    }
  }

  public void accumulateBackward(int i1b, int i2b, int i2e,
    float[][] p, float[][] e, float[][] d) {
    int n1 = e[0].length;
    int n1m = n1-1;
    int[] i1t = new int[n1];
    int i1p = i1b;
    d[i2b][i1b] = 0f;
    for (int i2=i2b-1; i2>=i2e; i2--) {
      int i2p = i2+1;
      float[] p2 = p[i2];
      float[] d2 = d[i2p];
      int[] kts = new int[n1];
      int k1b = i1p-round(p[i2+1][i1p])-200;
      int k1e = k1b+400;
      k1b = max(k1b,0);
      k1e = min(k1e,n1-1);
      float dmk = FLT_MAX;
      //for (int i1=k1b; i1<=k1e; ++i1) {
      for (int i1=0; i1<=n1m; ++i1) {
        d[i2][i1] = 0f;
        float p2i = p2[i1];
        int i1m = i1+round(p2i);
        int k1m = i1m+_dm;
        int k1p = i1m+_dp; 
        if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
        if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
        int k1t = k1m; float dmin = d2[k1m];
        for (int k1 = k1m+1; k1<=k1p; ++k1) {
          float dk = d[i2p][k1];
          if (dmin>dk) {dmin=dk;k1t=k1;}
        }
        kts[i1] = k1t;
        float e1 = k1t-i1-p2i;
        float e2 = i1-2f*k1t+i1t[k1t];
        float e3 = e[i2][i1];
        e1 =abs(e1);
        e2 *=e2;
        float dt = dmin+_c1*e1+_c2*e2+_c3*e3;
        float di = d[i2][i1]; if(di!=_emax) {dt+=di;}
        d[i2][i1] = dt;
        if(dt<dmk) {dmk=dt;i1p=i1;}
      }
      copy(kts,i1t);
    }
  }


  private void backtrack(
    int i1b, int i2b, int i2e, float[][] p, float[][] d, float[][] u) 
  {
    int n1 = d[0].length;
    int n1m = n1-1;
    int i1p = 0;
    float dmin = d[i2b][0];
    if (i1b>=0) {
      i1p = i1b;
      dmin = d[i2b][i1b];
    } else {
      for (int i1=1; i1<n1; ++i1) {
        float di = d[i2b][i1];
        if (di<dmin) {i1p = i1; dmin = di;}
      }
    }
    u[0][i2b] = i1p;
    u[1][i2b] = dmin;
    for (int i2=i2b-1; i2>=i2e; i2--) {
      float[] d2 = d[i2];
      int i1m = i1p-round(p[i2+1][i1p]);
      int k1m = i1m+_dm;
      int k1p = i1m+_dp;
      if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
      if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
      i1p = k1m; dmin = d2[k1m];
      for (int k1=k1m+1; k1<=k1p; ++k1) {
        float dk = d[i2][k1];
        if (dmin>dk) {dmin=dk;i1p = k1;}
      }
      u[0][i2] = i1p;
      u[1][i2] = dmin;
    }
  }

  private void backtrackReverse(
    int i1b, int i2b, int i2e, float[][] p, float[][] d, float[][] u) 
  {
    int n1 = d[0].length;
    int n1m = n1-1;
    int i1p = 0;
    float dmin = d[i2b][0];
    if (i1b>=0) {
      i1p = i1b;
      dmin = d[i2b][i1b];
    } else {
      for (int i1=1; i1<n1; ++i1) {
        float di = d[i2b][i1];
        if (di<dmin) {i1p = i1; dmin = d[i2b][i1];}
      }
    }
    u[0][i2b] = i1p;
    u[1][i2b] = dmin;
    for (int i2=i2b+1; i2<=i2e; i2++) {
      float[] d2 = d[i2];
      int i1m = i1p+round(p[i2-1][i1p]);
      int k1m = i1m+_dm;
      int k1p = i1m+_dp;
      if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
      if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
      i1p = k1m; dmin = d2[k1m];
      for (int k1=k1m+1; k1<=k1p; ++k1) {
        float dk = d[i2][k1];
        if (dmin>dk) {dmin=dk;i1p = k1;}
      }
      u[0][i2] = i1p;
      u[1][i2] = dmin;
    }
  }


  public void accumulateForward(
    int i1p, int i2b, float[][] p, float[][] el, float[][] e, float[][] d) 
  {
    int n2 = e.length;
    int n1 = e[0].length;
    int n1m = n1-1;
    int[] i1t = new int[n1];
    i1p = nearestIndex(i1p,p[1]);
    d[i2b][i1p] = 0f;
    for (int i2=i2b+1; i2<n2; i2++) {
      int i2m = i2-1;
      float[] p2 = p[i2];
      float[] d2 = d[i2m];
      int[] kts = new int[n1];
      int i1b = i1p+round(p2[i1p])-200;
      int i1e = i1b+400;
      i1b = max(i1b,0);
      i1e = min(i1e,n1-1);
      float dmk = FLT_MAX;
      for (int i1=i1b; i1<=i1e; ++i1) {
        float p2i = p2[i1];
        int i1m = i1-round(p2i);
        int k1m = i1m+_dm; 
        int k1p = i1m+_dp; 
        if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
        if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
        int k1t = k1m; float dmin = d2[k1m];
        for (int k1 = k1m+1; k1<=k1p; ++k1) {
          float dk = d[i2m][k1];
          if (dk<dmin) {dmin=dk;k1t=k1;}
        }
        kts[i1] = k1t;
        float e1 = i1-k1t-p2i;
        float e2 = i1-2f*k1t+i1t[k1t];
        float e3 = e[i2][i1];
        e3 = e[i2][i1];
        e1 =abs(e1);
        e2 *=e2;
        float dt = dmin+_c1*e1+_c2*e2+_c3*e3;
        d[i2][i1] = dt;
        if(dt<dmk) {dmk=dt;i1p=i1;}
      }
      copy(kts,i1t);
    }
  }

  public void accumulateForward(
    int i1p, int i2b, float[][] p, float[][] e, float[][] d) 
  {
    int n2 = e.length;
    int n1 = e[0].length;
    d[i2b][i1p  ] = 0f;
    d[i2b][i1p-1] = 0f;
    d[i2b][i1p-2] = 0f;
    d[i2b][i1p+1] = 0f;
    d[i2b][i1p+2] = 0f;
    for (int i2=i2b+1; i2<n2; i2++) {
      int i2m = i2-1;
      float[] p2 = p[i2];
      float[] e2 = e[i2];
      float[] d2 = d[i2m];
      for (int i1=0; i1<n1; ++i1) {
        float p2i = p2[i1];
        float i1m = i1-p2i;
        float k1t = i1m; 
        float dmin = FLT_MAX;
        for (float di=-0.8f; di<=0.8f; di+=0.1f) {
          float x1i = i1m+di;
          float dmi = _si.interpolate(n1,1.0,0.0,d2,x1i);
          if (dmi<dmin) {dmin=dmi;k1t=x1i;}
        }
        float e1i = i1-k1t-p2i;
        float e3i = e2[i1];
        e1i = abs(e1i);
        d[i2][i1] = dmin+_c1*e1i+_c3*e3i;
      }
    }
  }


  private int nearestIndex(int i1p, float[] p) {
    int k1 = i1p;
    int n1 = p.length;
    int d1 = n1;
    for (int i1=0; i1<n1; ++i1) {
      int t1 = i1-round(p[i1]);
      int di = abs(t1-i1p);
      if(di<d1) {k1=t1;d1=di;}
    }
    return k1;
  }

  public void accumulateBackwardX(
    int i2b,float[][] p, float[][] el, float[][] e, float[][] d) 
  {
    int n1 = e[0].length;
    int n1m = n1-1;
    int[] i1t = new int[n1];
    for (int i2=i2b-1; i2>=0; i2--) {
      int i2p = i2+1;
      float[] p2 = p[i2];
      float[] d2 = d[i2p];
      int[] kts = new int[n1];
      for (int i1=0; i1<n1; ++i1) {
        float p2i = p2[i1];
        int i1m = i1+round(p2i);
        int k1m = i1m+_dm;
        int k1p = i1m+_dp; 
        if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
        if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
        int k1t = k1m; float dmin = d2[k1m];
        for (int k1 = k1m+1; k1<=k1p; ++k1) {
          float dk = d[i2p][k1];
          if (dmin>dk) {dmin=dk;k1t=k1;}
        }
        float e1 = i1-k1t-p2i;
        float e2 = i1-2f*k1t+i1t[k1t];
        e1 *=e1;
        e2 *=e2;
        d[i2][i1] = dmin+e2;
        kts[i1] = k1t;
      }
      copy(kts,i1t);
    }
  }


  public void accumulateBackward(
    int i1p, int i2b,float[][] p, float[][] el, float[][] e, float[][] d) 
  {
    int n1 = e[0].length;
    int n1m = n1-1;
    int[] i1t = new int[n1];
    d[i2b][i1p] = 0f;
    for (int i2=i2b-1; i2>=0; i2--) {
      int i2p = i2+1;
      float[] p2 = p[i2];
      float[] d2 = d[i2p];
      int[] kts = new int[n1];
      int i1b = i1p-round(p[i2+1][i1p])-200;
      int i1e = i1b+400;
      i1b = max(i1b,0);
      i1e = min(i1e,n1-1);
      float dmk = FLT_MAX;
      for (int i1=i1b; i1<=i1e; ++i1) {
        d[i2][i1] = 0f;
        float p2i = p2[i1];
        int i1m = i1+round(p2i);
        int k1m = i1m+_dm;
        int k1p = i1m+_dp; 
        if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
        if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
        int k1t = k1m; float dmin = d2[k1m];
        for (int k1 = k1m+1; k1<=k1p; ++k1) {
          float dk = d[i2p][k1];
          if (dmin>dk) {dmin=dk;k1t=k1;}
        }
        kts[i1] = k1t;
        float e1 = k1t-i1-p2i;
        float e2 = i1-2f*k1t+i1t[k1t];
        float e3 = e[i2][i1];
        e1 =abs(e1);
        e2 *=e2;
        float dt = dmin+_c1*e1+_c2*e2+_c3*e3;
        d[i2][i1] = dt;
        if(dt<dmk) {dmk=dt;i1p=i1;}
      }
      copy(kts,i1t);
    }
  }

  public void accumulateBackward(
    float[][] p, float[][] el, float[][] e, float[][] d) 
  {
    int n2 = e.length;
    int n1 = e[0].length;
    int n1m = n1-1;
    int[] i1t = new int[n1];
    for (int i2=n2-2; i2>=0; i2--) {
      int i2p = i2+1;
      float[] p2 = p[i2];
      float[] d2 = d[i2p];
      int[] kts = new int[n1];
      for (int i1=0; i1<n1; ++i1) {
        float p2i = p2[i1];
        int i1m = i1+round(p2i);
        int k1m = i1m+_dm;
        int k1p = i1m+_dp; 
        if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
        if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
        int k1t = k1m; float dmin = d2[k1m];
        for (int k1 = k1m+1; k1<=k1p; ++k1) {
          float dk = d[i2p][k1];
          if (dmin>dk) {dmin=dk;k1t=k1;}
        }
        kts[i1] = k1t;
        float e1 = k1t-i1-p2i;
        float e2 = i1-2f*k1t+i1t[k1t];
        float e3 = e[i2][i1];
        e1 *=e1;
        e2 *=e2;
        d[i2][i1] += dmin+_c1*e1+_c2*e2+_c3*e3;
      }
      copy(kts,i1t);
    }
  }


  private void backtrack(
    float[][] p, float[][] d, float[] u) 
  {
    int n2 = d.length;
    int n1 = d[0].length;
    int n2m = n2-1;
    int n1m = n1-1;
    int i1p = 0;
    float dmin = d[n2m][0];
    for (int i1=1; i1<n1; ++i1) {
      if (d[n2m][i1]<dmin) {
        i1p = i1;
        dmin = d[n2m][i1];
      }
    }
    u[n2m] = i1p;
    for (int i2=n2m-1; i2>=0; i2--) {
      float[] d2 = d[i2];
      int i1m = i1p-round(p[i2+1][i1p]);
      int k1m = i1m+_dm;
      int k1p = i1m+_dp;
      if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
      if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
      i1p = k1m; dmin = d2[k1m];
      for (int k1=k1m+1; k1<=k1p; ++k1) {
        float dk = d[i2][k1];
        if (dmin>dk) {dmin=dk;i1p = k1;}
      }
      u[i2] = i1p;
    }
  }

  public void track(float[][] d, float[] u) {
    int n2 = u.length;
    int[] id = new int[1];
    for (int i2=0; i2<n2; ++i2) {
      float dmin = min(d[i2],id);
      u[i2] = id[0];
    }
  }

  private void backtrackReverse(
    float[][] p, float[][] d, float[] u) 
  {
    int n2 = d.length;
    int n1 = d[0].length;
    int n1m = n1-1;
    int i1p = 0;
    float dmin = d[0][0];
    for (int i1=1; i1<n1; ++i1) {
      if (d[0][i1]<dmin) {
        i1p = i1;
        dmin = d[0][i1];
      }
    }
    u[0] = i1p;
    for (int i2=1; i2<n2; i2++) {
      float[] d2 = d[i2];
      int i1m = i1p+round(p[i2-1][i1p]);
      int k1m = i1m+_dm;
      int k1p = i1m+_dp;
      if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
      if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
      i1p = k1m; dmin = d2[k1m];
      for (int k1=k1m+1; k1<=k1p; ++k1) {
        float dk = d[i2][k1];
        if (dmin>dk) {dmin=dk;i1p = k1;}
      }
      u[i2] = i1p;
    }
  }

  private void backtrackReverseX(
    int i1p, float[][] p, float[][] d, float[] u) 
  {
    int n2 = d.length;
    int n1 = d[0].length;
    int n1m = n1-1;
    float dmin = d[0][0];
    u[0] = i1p;
    for (int i2=1; i2<n2; i2++) {
      float[] d2 = d[i2];
      int i1m = i1p+round(p[i2-1][i1p]);
      int k1m = i1m+_dm;
      int k1p = i1m+_dp;
      if(k1m<0){k1m=0;} if(k1m>n1m) {k1m=n1m;}
      if(k1p<0){k1p=0;} if(k1p>n1m) {k1p=n1m;}
      i1p = k1m; dmin = d2[k1m];
      for (int k1=k1m+1; k1<=k1p; ++k1) {
        float dk = d[i2][k1];
        if (dmin>dk) {dmin=dk;i1p = k1;}
      }
      u[i2] = i1p;
    }
  }


  private int _k = 2;
  private int _dm = -1;
  private int _dp =  1;
  private int[] _lmin, _lmax;
  private float _emax = 1000000f;
  private float _emin = -100f;
  private SincInterpolator _si;
  private float _c1=1f;
  private float _c2=0.001f;
  private float _c3=1f;
  float _mark = 555f;
  int _wind = 10;
}

