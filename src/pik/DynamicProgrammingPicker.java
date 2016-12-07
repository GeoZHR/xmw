package pik;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping of sequences and images.
 * <p>
 * For sequences f and g, dynamic warping finds a sequence of 
 * shifts u such that f[i1] ~ g[i1+u[i1]], subject to a bound b1 
 * on strain, the rate at which the shifts u[i1] vary with sample 
 * index i1.
 * <p>
 * An increasing u[i1] = u[i1-1] + 1 implies that, between indices
 * i1-1 and i1, g[i1] is a stretched version of f[i1] ~ g[i1+u[i1]].
 * For, in this case, values in f for indices i1 and i1-1 are one 
 * sample apart, but corresponding values in g are two samples 
 * apart, which implies stretching by 100%. Likewise, a decreasing 
 * u[i1] = u[i1-1] - 1 implies squeezing by 100%.
 * <p>
 * In practice, 100% strain (stretching or squeezing) may be extreme.
 * Therefore, the upper bound on strain may be smaller than one. For 
 * example, if the bound b1 = 0.5, then |u[i1]-u[i1-1]| &le; 0.5.
 * <p>
 * For 2D images f and g, dynamic warping finds a 2D array of shifts
 * u[i2][i1] such that f[i2][i1] ~ g[i2][i1+u[i2][i1]], subject to 
 * bounds b1 and b2 on strains, the rates at which shifts u[i2][i1] 
 * vary with samples indices i1 and i2, respectively.
 * <p>
 * For 3D images f and g, dynamic warping finds a 3D array of shifts
 * u[i3][i2][i1] in a similar way. However, finding shifts for 3D 
 * images may require an excessive amount of memory. Dynamic image 
 * warping requires a temporary array of nlag*nsample floats, where 
 * the number of lags nlag = 1+shiftMax-shiftMin and nsample is the 
 * number of image samples. For 3D images, the product nlag*nsample 
 * is likely to be too large for the temporary array to fit in random-
 * access memory (RAM). In this case, shifts u are obtained by blending 
 * together shifts computed from overlapping subsets of the 3D image.
 * <p>
 * Estimated shifts u can be smoothed, and the extent of smoothing 
 * along each dimension is inversely proportional to the strain limit 
 * for that dimension. These extents can be scaled by specified factors 
 * for more or less smoothing. The default scale factors are zero, for 
 * no smoothing.
 * <p>
 * This class provides numerous methods, but typical applications
 * require only several of these, usually only the methods that find
 * and apply shifts. The many other methods are provided only for 
 * atypical applications and research.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.11.24
 */
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.dsp.LocalCorrelationFilter;

public class DynamicProgrammingPicker {

  public DynamicProgrammingPicker(int k, int dm, int dp) {
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
  private SincInterpolator _si;
  private float _c1=1f;
  private float _c2=0.001f;
  private float _c3=1f;
}

