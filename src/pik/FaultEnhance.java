package pik;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;



/**
 * Fault enhancing with ray tracing
 * @author Xinming Wu and Sergey Fomel
 * @version 2016.09.29
 */

public class FaultEnhance {
 
  public FaultEnhance(int gate, float an) {
    _gate = gate;
    _an = an;
  }

  public float[][][][] thin1(float[][][] f, float[][][] t) {
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n1 = f[0][0].length;
    final float[][][] ft = new float[n3][n2][n1];
    final float[][][] tt = new float[n3][n2][n1];
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][] f1 = new float[n3][n2];
      float[][] t1 = new float[n3][n2];
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        f1[i3][i2] = f[i3][i2][i1];
        t1[i3][i2] = t[i3][i2][i1];
      }}
      float[][][] ftt = thin1(f1,t1);
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        ft[i3][i2][i1] = ftt[0][i3][i2];
        tt[i3][i2][i1] = ftt[1][i3][i2];
      }}
    }});
    return new float[][][][]{ft,tt};
  }


  public float[][][] thin1(float[][] f, float[][] t) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] ff = new float[n2][n1];
    float[][] tt = new float[n2][n1];
    float pi = (float)(Math.PI/180.0);
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    int r = 5;
    float[] fs = new float[r*2+1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float ti = t[i2][i1]*pi;
      float d1 = -cos(ti);
      float d2 =  sin(ti);
      for (int i=-r; i<=r; ++i) 
        fs[i+r] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,i1+d1*i,i2+d2*i);
      rgf.apply0(fs,fs);
      float fi = fs[r];
      float fp = fs[r+1];
      float fm = fs[r-1];
      if(fp<fi&&fm<fi) {
        ff[i2][i1] = f[i2][i1];
        tt[i2][i1] = t[i2][i1];
      }
    }}
    return new float[][][]{ff,tt};
  }

  public float[][][] applyTracing2X(int rx, 
    int minTheta, int maxTheta, float[][][] fx) {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] fe = new float[n3][n2][n1];
    final RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(2);
    Parallel.loop(650,900,1,new Parallel.LoopInt() {
    public void compute(int i2) {
      System.out.println("i2="+i2);
      float[][] fx2 = new float[n3][n1];
      float[][] fs2 = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3)
        fx2[i3] = fx[i3][i2];
      float[][] ft2 = findRidges(0.001f,fx2);
      int[][] seeds = pickSeeds(4,0.1f,ft2);
      rgf.apply00(ft2,fs2);
      fs2 = sub(fs2,min(fs2));
      fs2 = div(fs2,max(fs2));
      float[][][] fph = enhanceInPolarSpace1(rx,minTheta,maxTheta,seeds,fs2);
      float[][] fe2 = fph[0];
      fe2 = sub(fe2,min(fe2));
      fe2 = div(fe2,max(fe2));
      fe2 = pow(fe2,0.5f);
      fe2 = sub(fe2,min(fe2));
      fe2 = div(fe2,max(fe2));
      for (int i3=0; i3<n3; ++i3) {
        fe[i3][i2] = fe2[i3];
      }
    }});
    return fe;
  }

  public float[][][] applyTracing3(int rx, 
    int minTheta, int maxTheta, float[][][] ft, float[][][] tt) {
    final int n3 = ft.length;
    final int n2 = ft[0].length;
    final int n1 = ft[0][0].length;
    final float[][][] fe = new float[n3][n2][n1];
    final RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(2);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      System.out.println("i3="+i3);
      float[][] fs3 = new float[n2][n1];
      float[][] ft3 = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        boolean ok = false;
        float tti = tt[i3][i2][i1];
        if(tti>=0f   && tti<45f  ) ok=true;
        if(tti>=135f && tti<=180f) ok=true;
        if(ok) ft3[i2][i1] = ft[i3][i2][i1];
      }}
      int[][] seeds = pickSeeds(4,0.3f,ft3);
      rgf.apply00(ft3,fs3);
      fs3 = sub(fs3,min(fs3));
      fs3 = div(fs3,max(fs3));
      float[][][] fph = enhanceInPolarSpace1(rx,minTheta,maxTheta,seeds,fs3);
      float[][] fe3 = fph[0];
      fe3 = sub(fe3,min(fe3));
      fe3 = div(fe3,max(fe3));
      fe3 = pow(fe3,0.5f);
      fe3 = sub(fe3,min(fe3));
      fe3 = div(fe3,max(fe3));
      fe[i3] = fe3;
    }});
    return fe;
  }


  public float[][][] applyTracing2(int rx, 
    int minTheta, int maxTheta, float[][][] ft, float[][][] tt) {
    final int n3 = ft.length;
    final int n2 = ft[0].length;
    final int n1 = ft[0][0].length;
    final float[][][] fe = new float[n3][n2][n1];
    final RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(2);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      System.out.println("i2="+i2);
      float[][] ft2 = new float[n3][n1];
      float[][] fs2 = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        float tti = tt[i3][i2][i1];
        if(tti>=45f && tti<=135f)
          ft2[i3][i1] = ft[i3][i2][i1];
      }}
      int[][] seeds = pickSeeds(4,0.3f,ft2);
      rgf.apply00(ft2,fs2);
      fs2 = sub(fs2,min(fs2));
      fs2 = div(fs2,max(fs2));
      float[][][] fph = enhanceInPolarSpace1(rx,minTheta,maxTheta,seeds,fs2);
      float[][] fe2 = fph[0];
      fe2 = sub(fe2,min(fe2));
      fe2 = div(fe2,max(fe2));
      fe2 = pow(fe2,0.5f);
      fe2 = sub(fe2,min(fe2));
      fe2 = div(fe2,max(fe2));
      for (int i3=0; i3<n3; ++i3)
        fe[i3][i2] = fe2[i3];
    }});
    return fe;
  }

  public float[][][] applyTracingV(
    int rx, int minTheta, int maxTheta, float[][][] fx, float[][][] ft) {
    float[][][][] ftt = thin1(fx,ft);
    float[][][] fe2 = applyTracing2(rx,minTheta,maxTheta,ftt[0],ftt[1]);
    //float[][][] fe3 = applyTracing3(rx,minTheta,maxTheta,ftt[0],ftt[1]);
    float[][][] fes = fe2;//add(fe2,fe3);
    fes = sub(fes,min(fes));
    fes = div(fes,max(fes));
    fes = pow(fes,0.5f);
    fes = sub(fes,min(fes));
    fes = div(fes,max(fes));
    return fes;
  }

  public float[][][][] applyTracing1(int rx, float sigma, float[][][] fx) {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] fe = new float[n3][n2][n1];
    final float[][][] ph = new float[n3][n2][n1];
    final RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(2);
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      System.out.println("i1="+i1);
      float[][] fx1 = new float[n3][n2];
      float[][] fs1 = new float[n3][n2];
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
          fx1[i3][i2] = fx[i3][i2][i1];
      float[][] ft1 = findRidges(0.001f,fx1);
      int[][] seeds = pickSeeds(4,0.1f,ft1);
      rgf.apply00(ft1,fs1);
      fs1 = sub(fs1,min(fs1));
      fs1 = div(fs1,max(fs1));
      float[][][] fph = enhanceInPolarSpace(rx,seeds,fs1);
      float[][] fe1 = fph[0];
      fe1 = sub(fe1,min(fe1));
      fe1 = div(fe1,max(fe1));
      fe1 = pow(fe1,0.5f);
      fe1 = sub(fe1,min(fe1));
      fe1 = div(fe1,max(fe1));
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        fe[i3][i2][i1] = fe1[i3][i2];
        ph[i3][i2][i1] = fph[1][i3][i2];
      }}
    }});
    return new float[][][][] {fe,ph};
  }

  public float[][][] enhanceInPolarSpace1(int rx, 
    int minTheta, int maxTheta, int[][] seeds, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int ns = seeds[0].length;
    int nr = round(rx);
    int np = nr*2;
    int[] k1s = new int[np];
    int[] k2s = new int[np];
    float[][] gx = new float[n2][n1];
    float[][] gm = new float[n2][n1];
    float[][] ph = new float[n2][n1];
    float[][] wx = exp(mul(-1,fx));
    for (int is=0; is<ns; ++is) {
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      float[][] gs = applyPolarTransform(minTheta,k1,k2,rx,wx);
      int i2p = 90;
      int na = gs.length;
      float gmax = FLT_MAX;
      for (int i2=0; i2<=maxTheta-minTheta; i2++) {
        float gsum1 = sum(gs[i2]);
        float gsum2 = sum(gs[na-i2-1]);
        if(gsum1<gmax) {
          i2p = i2;
          gmax = gsum1;
        }
        if(gsum2<gmax) {
          i2p = na-i2-1;
          gmax = gsum2;
        }
      }
      float[][] ws = matrixTransform(gs);
      float[] p1 = forwardPick(i2p,ws);
      float[] p2 = backwardPick(round(p1[np-1]),ws);
      float[] rs = new float[np];
      float[] as = new float[np];
      for (int ip=0; ip<nr; ip++) {
        rs[ip   ] = nr-ip;
        as[ip   ] = p2[ip]+180+minTheta;
        rs[ip+nr] = ip;
        as[ip+nr] = p2[ip+nr]+minTheta;
      }
      float[][] xs = reverseTransform(k1,k2,rs,as);
      float fpa = 0.0f;
      for (int ip=0; ip<np; ip++) {
        int i1 = round(xs[0][ip]);
        int i2 = round(xs[1][ip]);
        k1s[ip] = i1;
        k2s[ip] = i2;
        i1 = max(i1,0);
        i2 = max(i2,0);
        i1 = min(i1,n1-1);
        i2 = min(i2,n2-1);
        fpa += fx[i2][i1];
      }
      fpa /= np;
      for (int ip=0; ip<np; ++ip) {
        int i1 = k1s[ip];
        int i2 = k2s[ip];
        for (int d2=-1; d2<=1; d2++) {
        for (int d1=-1; d1<=1; d1++) {
          int c1 = i1+d1;
          int c2 = i2+d2;
          if(c1>=0&&c1<n1&&c2>=0&&c2<n2) {
            gx[c2][c1] += fpa;
            float asi = as[ip];
            if(asi>180) asi -=180;
            if (fpa>gm[c2][c1]) {
              ph[c2][c1] = asi;
              gm[c2][c1] = fpa;
            }
          }
        }}
      }
    }
    return new float[][][]{gx,ph};
  }

  public float[][][] enhanceInPolarSpace(
    int rx, int[][] seeds, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int ns = seeds[0].length;
    int nr = rx;
    int np = nr*2;
    int[] k1s = new int[np];
    int[] k2s = new int[np];
    float[][] gx = new float[n2][n1];
    float[][] gm = new float[n2][n1];
    float[][] ph = new float[n2][n1];
    float[][] wx = exp(mul(-1,fx));
    for (int is=0; is<ns; ++is) {
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      float[][] gs = applyPolarTransform(k1,k2,rx,wx);
      int i2p = 90;
      float gmax = FLT_MAX;
      for (int i2=0; i2<180; i2++) {
        float gsum = sum(gs[i2]);
        if(gsum<gmax) {
          i2p = i2;
          gmax = gsum;
        }
      }
      float[][] ws = matrixTransform(gs);
      float[] p1 = forwardPick(i2p,ws);
      float[] p2 = backwardPick(round(p1[np-1]),ws);
      float[] rs = new float[np];
      float[] as = new float[np];
      for (int ip=0; ip<nr; ip++) {
        rs[ip   ] = nr-ip;
        as[ip   ] = p2[ip]+180;
        rs[ip+nr] = ip;
        as[ip+nr] = p2[ip+nr];
      }
      float[][] xs = reverseTransform(k1,k2,rs,as);
      float fpa = 0.0f;
      for (int ip=0; ip<np; ip++) {
        int i1 = round(xs[0][ip]);
        int i2 = round(xs[1][ip]);
        k1s[ip] = i1;
        k2s[ip] = i2;
        i1 = max(i1,0);
        i2 = max(i2,0);
        i1 = min(i1,n1-1);
        i2 = min(i2,n2-1);
        fpa += fx[i2][i1];
      }
      fpa /= np;
      for (int ip=0; ip<np; ++ip) {
        int i1 = k1s[ip];
        int i2 = k2s[ip];
        for (int d2=-1; d2<=1; d2++) {
        for (int d1=-1; d1<=1; d1++) {
          int c1 = i1+d1;
          int c2 = i2+d2;
          if(c1>=0&&c1<n1&&c2>=0&&c2<n2) {
            gx[c2][c1] += fpa;
            float asi = as[ip];
            if(asi>180) asi-=180;
            if(fpa>gm[c2][c1]) {
              ph[c2][c1] = asi;
              gm[c2][c1] = fpa;
            }
          }
        }}
      }
    }
    return new float[][][]{gx,ph};
  }
  public float[][] applyPolarTransform(
    int minTheta, float x1, float x2, float rx, float[][] fx) {
    int nr = round(rx);
    int np = 180-2*minTheta+1;
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[np][nr*2];
    for (int ip=0; ip<np; ++ip) {
      float ph1 = (float)toRadians(ip+minTheta);
      float sph = sin(ph1);
      float cph = cos(ph1);
      for (int ir=0; ir<nr; ++ir) {
        float r11 = ir*sph;
        float r12 = ir*cph;
        int i11 = round(x1-r11);
        int i12 = round(x2-r12);
        int i21 = round(x1+r11);
        int i22 = round(x2+r12);
        i11 = min(i11,n1-1);
        i12 = min(i12,n2-1);
        i21 = min(i21,n1-1);
        i22 = min(i22,n2-1);
        i11 = max(i11,0);
        i12 = max(i12,0);
        i21 = max(i21,0);
        i22 = max(i22,0);
        gx[ip][ir+nr  ] = fx[i12][i11];
        gx[ip][nr-ir-1] = fx[i22][i21];
      }
    }
    return gx;
  }

  public float[][] applyPolarTransform(
    float x1, float x2, float rx, float[][] fx) {
    int np = 180;
    int nr = round(rx);
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[np][nr*2];
    for (int ip=0; ip<np; ++ip) {
      float ph1 = (float)toRadians(ip);
      float sph = sin(ph1);
      float cph = cos(ph1);
      for (int ir=0; ir<nr; ++ir) {
        float r11 = ir*sph;
        float r12 = ir*cph;
        int i11 = round(x1-r11);
        int i12 = round(x2-r12);
        int i21 = round(x1+r11);
        int i22 = round(x2+r12);
        i11 = min(i11,n1-1);
        i12 = min(i12,n2-1);
        i21 = min(i21,n1-1);
        i22 = min(i22,n2-1);
        i11 = max(i11,0);
        i12 = max(i12,0);
        i21 = max(i21,0);
        i22 = max(i22,0);
        gx[ip][ir+nr  ] = fx[i12][i11];
        gx[ip][nr-ir-1] = fx[i22][i21];
      }
    }
    return gx;
  }
  public float[][] reverseTransform(
    float c1, float c2, float[] rs, float[] as) {
    int np = rs.length;
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    for (int ip=0; ip<np; ++ip) {
      float ri = rs[ip];
      float phi = (float)toRadians(as[ip]);
      float rx1 = ri*sin(phi);
      float rx2 = ri*cos(phi);
      x1[ip] = c1-rx1;
      x2[ip] = c2-rx2;
    }
    return new float[][]{x1,x2};
  }

  public float[][] applyForEnhanceX(
    int d1, int d2, float sigma, int[][] seeds, float[][] fx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int ns = seeds[0].length;
    float[][] gx = new float[n2][n1];
    float[][] wx = applyForWeight(fx);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    for (int is=0; is<ns; ++is) {
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      int b1 = max(k2-d2,0);
      int e1 = min(k2+d2+1,n2);
      int bb2 = max(k1-d1,0);
      int eb2 = k1+1;
      int bf2 = eb2;
      int ef2 = min(k1+d1+1,n1);
      int m1 = ef2-bb2; 
      float[] fp = new float[m1];
      float[] pb = backwardPick(k2, b1, bb2, e1, eb2, wx);
      float[] pf =  forwardPick(k2, b1, bf2, e1, ef2, wx);
      for (int i1=bb2; i1<ef2; ++i1) {
        int i2=0;
        if (i1<eb2) {i2 = round(pb[i1-bb2]+b1);}
        else        {i2 = round(pf[i1-bf2]+b1);}
        i2 = max(i2,0);
        i2 = min(i2,n2-1);
        fp[i1-bb2] = fx[i2][i1];
      }
      rgf.apply0(fp,fp);
      for (int i1=bb2; i1<ef2; ++i1) {
        int i2=0;
        if (i1<eb2) {i2 = round(pb[i1-bb2]+b1);}
        else        {i2 = round(pf[i1-bf2]+b1);}
        i2 = max(i2,0);
        i2 = min(i2,n2-1);
        if (abs(i1-eb2)>10) {
          gx[i2][i1] += fp[i1-bb2];
        }
      }
    }
    return gx;
  }


  public float[][] applyForEnhance(float sigma, int[][] seeds, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int ns = seeds[0].length;
    float[][] gx = new float[n2][n1];
    float[][] wx = applyForWeight(fx);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    for (int is=0; is<ns; ++is) {
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      if (k1==0) {
        float[] pt = forwardPick(k2,wx);
        float[] fp = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          fp[i1] = fx[i2][i1];
        }
        rgf.apply0(fp,fp);
        for (int i1=0; i1<n1; ++i1) {
          if (abs(i1-k1)>20) {continue;}
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          float fpi = fp[i1];
          float gxi = gx[i2][i1];
          if (fpi>gxi) {gx[i2][i1] = fpi;}
        }
      }else if (k1==n1-1) {
        float[] pt = backwardPick(k2,wx);
        float[] fp = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          fp[i1] = fx[i2][i1];
        }
        rgf.apply0(fp,fp);
        for (int i1=0; i1<n1; ++i1) {
          if (abs(i1-k1)>20) {continue;}
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          float fpi = fp[i1];
          float gxi = gx[i2][i1];
          if (fpi>gxi) {gx[i2][i1] = fpi;}
        }
      } else {
        int d12 = round((k1+1)/2);
        int b12 = k2-d12; b12 = max(0,b12);
        int e12 = k2+d12; e12 = min(n2-1,e12);
        int d22 = round((n1-k1-1)/2);
        int b22 = k2-d22; b22 = max(0,b22);
        int e22 = k2+d22; e22 = min(n2-1,e22);
        float[][] wx1 = copy(e12-b12+1,k1+1,b12,0,wx); 
        float[][] wx2 = copy(e22-b22+1,n1-k1-1,b22,k1,wx); 
        float[] p1 = backwardPick(k2-b12,wx1);
        float[] p2 = forwardPick(k2-b22,wx2);
        float[] fp = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          if (abs(i1-k1)>20) {continue;}
          int i2=0;
          if (i1<=k1) {
            i2 = round(p1[i1]+b12);
          } else {
            i2 = round(p2[i1-k1-1]+b22);
          }
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          fp[i1] = fx[i2][i1];
        }
        rgf.apply0(fp,fp);
        for (int i1=0; i1<n1; ++i1) {
          int i2 = 0;
          if (i1<=k1) {
            i2 = round(p1[i1]+b12);
          } else {
            i2 = round(p2[i1-k1-1]+b22);
          }
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          float fpi = fp[i1];
          float gxi = gx[i2][i1];
          if (fpi>gxi) {gx[i2][i1] = fpi;}
        }
      }
    }
    return gx;
  }

  public float[][] thin(float fm, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n2][n1];
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=1;i2<n2-1;i2++) {
    for (int i1=0;i1<n1  ;i1++) {
      float fsi = fs[i2][i1];
      float fsm = fs[i2-1][i1];
      float fsp = fs[i2+1][i1];
      if (fsm<fsi&&fsp<fsi&&fsi>fm) {
        ft[i2  ][i1] = fsi;
        ft[i2-1][i1] = fsm;
        ft[i2+1][i1] = fsp;
      }
    }}
    return ft;
  }

  public float[][] thinX(float fm, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n2][n1];
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=1;i2<n2-1;i2++) {
    for (int i1=1;i1<n1-1;i1++) {
      float fsi = fs[i2][i1];
      float fsm2 = fs[i2-1][i1];
      float fsp2 = fs[i2+1][i1];
      float fsm1 = fs[i2][i1-1];
      float fsp1 = fs[i2][i1+1];
      if (fsm2<fsi&&fsp2<fsi&&fsi>fm) {
        ft[i2  ][i1] = fsi;
        //ft[i2-1][i1] = fsm;
        //ft[i2+1][i1] = fsp;
      }
      if (fsm1<fsi&&fsp1<fsi&&fsi>fm) {
        ft[i2][i1] = fsi;
        //ft[i2][i1-1] = fsm;
        //ft[i2][i1+1] = fsp;
      }
    }}

    return ft;
  }

  public float[][] thin2(float fm, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n2][n1];
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=1;i2<n2-1;i2++) {
    for (int i1=1;i1<n1-1;i1++) {
      float fsi  = fs[i2][i1];
      float fxi  = fx[i2][i1];
      float fsm2 = fs[i2-1][i1];
      float fsp2 = fs[i2+1][i1];
      if (fsm2<fsi&&fsp2<fsi&&fxi>fm)
        ft[i2][i1] = fxi;
    }}
    return ft;
  }


  public float[][] findRidges(float fm, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n2][n1];
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=1;i2<n2-1;i2++) {
    for (int i1=1;i1<n1-1;i1++) {
      float fsi  = fs[i2][i1];
      float fxi  = fx[i2][i1];
      float fsm2 = fs[i2-1][i1];
      float fsp2 = fs[i2+1][i1];
      float fsm1 = fs[i2][i1-1];
      float fsp1 = fs[i2][i1+1];
      boolean onRidge = false;
      if (fsm2<fsi&&fsp2<fsi&&fxi>fm)
        onRidge = true;
      if (fsm1<fsi&&fsp1<fsi&&fxi>fm)
        onRidge = true;
      if (onRidge) {
        ft[i2][i1] = fxi;
      }
    }}
    return ft;
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

  public float[][] faultLikeFit1(int r, float[][] fl) {
    int n2 = fl.length;
    int n1 = fl[0].length;
    float[][] flr = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=r; i1<n1-r-1; ++i1) {
      float[] ft = new float[r*2+1];
      for (int r1=-r;r1<=r;r1++)
        ft[r1+r] = fl[i2][i1+r1];
      float[] abc = parabolicFit(ft);
      float a = abc[0];
      float b = abc[1];
      float s = -abs(2*a*r+b)/(2*a);
      flr[i2][i1] = fl[i2][i1]/(s+0.0001f);
    }}
    return flr;
  }

  public float[][] faultLikeFit(int r, float[][] fl) {
    int n2 = fl.length;
    int n1 = fl[0].length;
    float[][] f1 = faultLikeFit1(r,fl);
    float[][] f2 = faultLikeFit2(r,fl);
    float[][] ft = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float f1i = f1[i2][i1];
      float f2i = f2[i2][i1];
      ft[i2][i1] = max(f1i,f2i);
    }}
    return ft;
  }


  public float[][] faultLikeFit2(int r, float[][] fl) {
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

  public float[][] seedsToImage(int[][] seeds, float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    int ns = seeds[0].length;
    float[][] fs = new float[n2][n1];
    for (int is=0; is<ns; ++is) {
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      fs[k2][k1] = ft[k2][k1];
    }
    return fs;
  }

  public int[][] pickSeeds3(int d, float fm, float[][] ft, float[][] tt) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    ArrayList<Float> fs = new ArrayList<Float>();
    ArrayList<Integer> ks = new ArrayList<Integer>();
    for (int i2=0;i2<n2;i2+=1) {
    for (int i1=0;i1<n1;i1+=1) {
      float fti = ft[i2][i1];
      float tti = tt[i2][i1];
      boolean ok = false;
      if(tti>=0f   && tti<45f  ) ok=true;
      if(tti>=135f && tti<=180f) ok=true;
      if(fti>fm && ok) {
        fs.add(fti);
        ks.add(n1*i2+i1);
      }
    }}
    int np = fs.size();
    int[] ka = new int[np];
    int[] ia = new int[np];
    float[] fa = new float[np];
    for (int ip=0; ip<np; ++ip) {
      ia[ip] = ip;
      ka[ip] = ks.get(ip);
      fa[ip] = fs.get(ip);
    }
    quickIndexSort(fa,ia);
    int[][] mark = new int[n2][n1];
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    for (int ip=np-1; ip>=0; --ip) {
      int ki = ka[ia[ip]];
      int i1 = ki%n1;
      int i2 = round(ki/n1);
      int b1 = i1-d; b1=max(b1,0);
      int b2 = i2-d; b2=max(b2,0);
      int e1 = i1+d; e1=min(e1,n1-1);
      int e2 = i2+d; e2=min(e2,n2-1);
      boolean ok = true;
      for (int k2=b2;k2<=e2;k2++) {
      for (int k1=b1;k1<=e1;k1++) {
        if(mark[k2][k1]==1) {
          ok=false;
          break;
        }
      }}
      if(ok) {
        k1s.add(i1);
        k2s.add(i2);
        mark[i2][i1] = 1;
      }
    }
    int ns = k1s.size();
    int[][] seeds = new int[2][ns];
    for (int is=0; is<ns; ++is) {
      seeds[0][is] = k1s.get(is); 
      seeds[1][is] = k2s.get(is); 
    }
    return seeds;
  }

  public int[][] pickSeeds2(int d, float fm, float[][] ft, float[][] tt) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    ArrayList<Float> fs = new ArrayList<Float>();
    ArrayList<Integer> ks = new ArrayList<Integer>();
    for (int i2=0;i2<n2;i2+=1) {
    for (int i1=0;i1<n1;i1+=1) {
      float fti = ft[i2][i1];
      float tti = tt[i2][i1];
      if(fti>fm && tti>=45f && tti<=135f) {
        fs.add(fti);
        ks.add(n1*i2+i1);
      }
    }}
    int np = fs.size();
    int[] ka = new int[np];
    int[] ia = new int[np];
    float[] fa = new float[np];
    for (int ip=0; ip<np; ++ip) {
      ia[ip] = ip;
      ka[ip] = ks.get(ip);
      fa[ip] = fs.get(ip);
    }
    quickIndexSort(fa,ia);
    int[][] mark = new int[n2][n1];
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    for (int ip=np-1; ip>=0; --ip) {
      int ki = ka[ia[ip]];
      int i1 = ki%n1;
      int i2 = round(ki/n1);
      int b1 = i1-d; b1=max(b1,0);
      int b2 = i2-d; b2=max(b2,0);
      int e1 = i1+d; e1=min(e1,n1-1);
      int e2 = i2+d; e2=min(e2,n2-1);
      boolean ok = true;
      for (int k2=b2;k2<=e2;k2++) {
      for (int k1=b1;k1<=e1;k1++) {
        if(mark[k2][k1]==1) {
          ok=false;
          break;
        }
      }}
      if(ok) {
        k1s.add(i1);
        k2s.add(i2);
        mark[i2][i1] = 1;
      }
    }
    int ns = k1s.size();
    int[][] seeds = new int[2][ns];
    for (int is=0; is<ns; ++is) {
      seeds[0][is] = k1s.get(is); 
      seeds[1][is] = k2s.get(is); 
    }
    return seeds;
  }


  public int[][] pickSeeds(int d, float fm, float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    ArrayList<Float> fs = new ArrayList<Float>();
    ArrayList<Integer> ks = new ArrayList<Integer>();
    for (int i2=0;i2<n2;i2+=1) {
    for (int i1=0;i1<n1;i1+=1) {
      if(ft[i2][i1]>fm) {
        ks.add(n1*i2+i1);
        fs.add(ft[i2][i1]);
      }
    }}
    int np = fs.size();
    int[] ka = new int[np];
    int[] ia = new int[np];
    float[] fa = new float[np];
    for (int ip=0; ip<np; ++ip) {
      ia[ip] = ip;
      ka[ip] = ks.get(ip);
      fa[ip] = fs.get(ip);
    }
    quickIndexSort(fa,ia);
    int[][] mark = new int[n2][n1];
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    for (int ip=np-1; ip>=0; --ip) {
      int ki = ka[ia[ip]];
      int i1 = ki%n1;
      int i2 = round(ki/n1);
      int b1 = i1-d; b1=max(b1,0);
      int b2 = i2-d; b2=max(b2,0);
      int e1 = i1+d; e1=min(e1,n1-1);
      int e2 = i2+d; e2=min(e2,n2-1);
      boolean ok = true;
      for (int k2=b2;k2<=e2;k2++) {
      for (int k1=b1;k1<=e1;k1++) {
        if(mark[k2][k1]==1) {
          ok=false;
          break;
        }
      }}
      if(ok) {
        k1s.add(i1);
        k2s.add(i2);
        mark[i2][i1] = 1;
      }
    }
    int ns = k1s.size();
    int[][] seeds = new int[2][ns];
    for (int is=0; is<ns; ++is) {
      seeds[0][is] = k1s.get(is); 
      seeds[1][is] = k2s.get(is); 
    }
    return seeds;
  }



  public int[][] findSeedsX(int d1, int d2, float fm, 
    float[][] fx, float[][] ft) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] fs = new float[n2][n1];
    ArrayList<Integer> k1 = new ArrayList<Integer>();
    ArrayList<Integer> k2 = new ArrayList<Integer>();
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=5;i2<n2-5;i2+=d2) {
    for (int i1=5;i1<n1-5;i1+=d1) {
      float fsi = fs[i2][i1];
      float fsm2 = fs[i2-1][i1];
      float fsp2 = fs[i2+1][i1];
      float fsm1 = fs[i2][i1-1];
      float fsp1 = fs[i2][i1+1];
      boolean onRidge = false;
      if (fsm2<fsi&&fsp2<fsi&&fsi>fm)
        onRidge = true;
      if (fsm1<fsi&&fsp1<fsi&&fsi>fm)
        onRidge = true;
      if (onRidge) {
        k1.add(i1);
        k2.add(i2);
        ft[i2][i1] = fsi;
      }
    }}
    int ns = k1.size();
    int[][] seeds = new int[2][ns];
    for (int is=0; is<ns; ++is) {
      seeds[0][is] = k1.get(is); 
      seeds[1][is] = k2.get(is); 
    }
    return seeds;
  }


  public int[][] findSeeds(int d1, int d2, float fm, 
    float[][] fx, float[][] ft) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] fs = new float[n2][n1];
    ArrayList<Integer> k1 = new ArrayList<Integer>();
    ArrayList<Integer> k2 = new ArrayList<Integer>();
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=5;i2<n2-5;i2+=d2) {
    for (int i1=5;i1<n1-5;i1+=d1) {
      float fsi = fs[i2][i1];
      float fsm = fs[i2-1][i1];
      float fsp = fs[i2+1][i1];
      if (fsm<fsi&&fsp<fsi&&fsi>fm) {
        k1.add(i1);
        k2.add(i2);
        ft[i2][i1] = fsi;
      }
    }}
    int ns = k1.size();
    int[][] seeds = new int[2][ns];
    for (int is=0; is<ns; ++is) {
      seeds[0][is] = k1.get(is); 
      seeds[1][is] = k2.get(is); 
    }
    return seeds;
  }

  public float[][] matrixTransform(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i1][i2]  = fx[i2][i1];
    }}
    return gx;
  }

  public float[][] applyForWeight(float[][] vel) {
    int n2 = vel.length;
    int n1 = vel[0].length;
    float[][] w = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      w[i1][i2]  = exp(-vel[i2][i1]);
    }}
    return w;
  }

  public float[][] applyForWeightX(float[][] vel) {
    int n2 = vel.length;
    int n1 = vel[0].length;
    float[][] w = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      w[i2][i1]  = exp(-vel[i2][i1]);
    }}
    return w;
  }


  public float[] applyForPath(
    int i0, float sig, float[][] vel) {
    int n1 = vel[0].length;
    float[][] wx = applyForWeight(vel);
    float[] p1 = forwardPick(i0,wx);
    float[] p2 = backwardPick(round(p1[n1-1]),wx);
    float[] b = new float[n1];
    float[] r = new float[n1];
    float[] ws = new float[n1];
    makeRhsWeights(p1,p2,vel,b,ws);
    VecArrayFloat1 vb = new VecArrayFloat1(b);
    VecArrayFloat1 vr = new VecArrayFloat1(r);
    Smoother1 smoother1 = new Smoother1(sig);
    A1 a1 = new A1(smoother1,ws);
    CgSolver cs = new CgSolver(0.001,200);
    smoother1.applyTranspose(b);
    cs.solve(a1,vb,vr);
    return r;
  }

  public float[][][] applyForSurfaceInline(
    float sig1, float sig2, final float[][][] vel) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    float[][] p1 = new float[n2][n3];
    float[][] w1 = new float[n3][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        w1[i3][i1] = exp(-vel[i3][i2][i1]);
      }}
      p1[i2] = forwardPick(n1-1,w1);
    }
    float[][] b = new float[n3][n2];
    float[][] r = new float[n3][n2];
    float[][] ws = new float[n3][n2];
    makeRhsWeightsInline(p1,vel,b,ws);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(sig1,sig2);
    A2 a2 = new A2(smoother2,ws);
    CgSolver cs = new CgSolver(0.001,200);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    float[][] p1t = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      p1t[i3][i2] = p1[i2][i3];
    }}
    return new float[][][]{p1t,r};
  }


  public float[][][] applyForSurface(
    float sig1, float sig2, final float[][][] vel) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    float[][] p1 = new float[n2][n3];
    float[][] w1 = new float[n3][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        w1[i3][i1] = exp(-vel[i3][i2][i1]);
      }}
      p1[i2] = forwardPick(n1-1,w1);
    }
    float[][] p2 = new float[n3][n2];
    float[][] w2 = new float[n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        w2[i2][i1] = exp(-vel[i3][i2][i1]);
      }}
      int i0 = round(p1[0][i3]);
      i0 = min(i0,n1-1);
      i0 = max(i0,0);
      p2[i3] = forwardPick(i0,w2);
    }
    float[][] b = new float[n3][n2];
    float[][] r = new float[n3][n2];
    float[][] ws = new float[n3][n2];
    makeRhsWeights(p1,p2,vel,b,ws);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(sig1,sig2);
    A2 a2 = new A2(smoother2,ws);
    CgSolver cs = new CgSolver(0.001,200);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    float[][] p1t = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      p1t[i3][i2] = p1[i2][i3];
    }}
    return new float[][][]{p1t,p2,r};
  }

  public float[] backwardPick(int i0, float[][] wx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    float[] prob = new float[_gate*2-1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-1][i1]+wx[n2-1][i0]);
	    tt[n2-1][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-2][i1]+wx[n2-1][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][i1] = i0;
	    tt[n2-2][i1]=prev[i1];
    }

    for (int i2=n2-3; i2>=0; i2--) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2+1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    forwardTrack(p, next, what);
    return p;
  }

  public float[] backwardPick(int i0, int b1, int b2, 
    int e1, int e2, float[][] wx) {
    int n2 = e2-b2;
    int n1 = e1-b1;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    float[] prob = new float[_gate*2-1];
    for (int k1=0; k1<n1; ++k1) {
      dist[k1] = sqrt(k1*k1+_an*_an); 
    }
	  for (int i1=b1; i1<e1; i1++) {
      int k1 = i1-b1;
	    float wi = 0.5f*(wx[e2-1][i1]+wx[e2-1][i0]);
	    tt[n2-1][k1] = abs(i1-i0)*wi;
	  }
    for (int i1=b1; i1<e1; i1++) {
      int k1 = i1-b1;
	    float wi = 0.5f*(wx[e2-2][i1]+wx[e2-1][i0]);
	    prev[k1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][k1] = i0-b1;
	    tt[n2-2][k1]=prev[k1];
    }
    for (int i2=e2-3; i2>=b2; i2--) {
      int k2 = i2-b2;
	    for (int i1=b1; i1<e1; i1++) {
        int k1 = i1-b1;
	      float wi = wx[i2][i1];
	      int kb = max(k1-_gate,-1);
	      int ke = min(k1+_gate,n1);
	      float c = FLT_MAX;
	      int kc = -1;
	      for (int k=kb+1; k<ke; k++) {
		      float w2 = 0.5f*(wi+wx[i2+1][k+b1]);
		      float d = dist[abs(k1-k)]*w2+prev[k];
		      int kt = k-kb-1;
		      if (d < c) {
		        c =	d;
		        kc = kt;
		      }
		      prob[kt]=d;
	      }
        float[] vs = find_minimum(kc,ke-kb-1,kb+1,c,what[k2][k1],prob);
	      next[k1]= vs[0];
        what[k2][k1] = vs[1];
	    }
	    for (int k1=0; k1<n1; k1++) {
	      prev[k1]=next[k1];
	      tt[k2][k1]=prev[k1];
	    }
    }
    forwardTrack(p, next, what);
    return p;
  }



  public float[] backwardPick(int i0, float[][] wx, float[][] tx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    float[] prob = new float[_gate*2-1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-1][i1]+wx[n2-1][i0]);
	    tt[n2-1][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-2][i1]+wx[n2-1][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][i1] = i0;
	    tt[n2-2][i1]=prev[i1];
    }

    for (int i2=n2-3; i2>=0; i2--) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2+1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    forwardTrack(p, next, what);
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      tx[i1][i2]  = tt[i2][i1];
    }}
    return p;
  }

  public float[] forwardPick(int i0, float[][] wx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    float[] prob = new float[_gate*2-1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[0][i1]+wx[0][i0]);
	    tt[0][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[1][i1]+wx[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[1][i1] = i0;
	    tt[1][i1]=prev[i1];
    }

    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      int ic = -1;
	      float c = FLT_MAX;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    backwardTrack(p, next, what);
    return p;
  }

  public float[] forwardPick(int i0, int b1, int b2, 
    int e1, int e2, float[][] wx) {
    int n2 = e2-b2;
    int n1 = e1-b1;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    float[] prob = new float[_gate*2-1];
    for (int k1=0; k1<n1; ++k1) {
      dist[k1] = sqrt(k1*k1+_an*_an); 
    }
	  for (int i1=b1; i1<e1; i1++) {
	    float wi = 0.5f*(wx[b2][i1]+wx[b2][i0]);
	    tt[0][i1-b1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=b1; i1<e1; i1++) {
      int k1 = i1-b1;
	    float wi = 0.5f*(wx[b2+1][i1]+wx[b2][i0]);
	    prev[k1] = dist[abs(i1-i0)]*wi;
	    what[1][k1] = i0-b1;
	    tt[1][k1]=prev[k1];
    }

    for (int i2=b2+2; i2<e2; i2++) {
      int k2 = i2-b2; 
	    for (int i1=b1; i1<e1; i1++) {
        int k1 = i1-b1; 
	      float wi = wx[i2][i1];
	      int kb = max(k1-_gate,-1);
	      int ke = min(k1+_gate,n1);
	      int kc = -1;
	      float c = FLT_MAX;
	      for (int k=kb+1; k<ke; k++) {
		      float w2 = 0.5f*(wi+wx[i2-1][k+b1]);
		      float d = dist[abs(k1-k)]*w2+prev[k];
		      int kt = k-kb-1;
		      if (d < c) {
		        c =	d;
		        kc = kt;
		      }
		      prob[kt]=d;
	      }
        float[] vs = find_minimum(kc,ke-kb-1,kb+1,c,what[k2][k1],prob);
	      next[k1]= vs[0];
        what[k2][k1] = vs[1];
	    }
	    for (int k1=0; k1<n1; k1++) {
	      prev[k1]=next[k1];
	      tt[k2][k1]=prev[k1];
	    }
    }
    backwardTrack(p, next, what);
    return p;
  }


  public float[] forwardPick(int i0, float[][] wx, float[][] tx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    float[] prob = new float[_gate*2-1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[0][i1]+wx[0][i0]);
	    tt[0][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[1][i1]+wx[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[1][i1] = i0;
	    tt[1][i1]=prev[i1];
    }

    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    backwardTrack(p, next, what);
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      tx[i1][i2]  = tt[i2][i1];
    }}
    return p;
  }



  private float[] find_minimum(
    int ic, int nc, int jc, float c, float pick, float[] prob)
  {
    float fm, f0, fp, a, b;
    if (0==ic) {
	    ic++;
	    fm=c;
	    f0=prob[ic];
	    fp=prob[ic+1];
    } else if (nc-1==ic) {
	    ic--;
	    fm=prob[ic-1];
	    f0=prob[ic];
	    fp=c;
    } else {
	    fm=prob[ic-1];
	    f0=c;
	    fp=prob[ic+1];
    }

    ic += jc;
    a = fm+fp-2f*f0;
    if (a <= 0.) { /* no minimum */
	    if (fm < f0 && fm < fp) {
	      pick = ic-1;
	      return new float[]{fm,pick};
	    } 
	    if (fp < f0 && fp < fm) {
	      pick = ic+1;
	      return new float[]{fp,pick};
	    } 
	    pick = ic;
	    return new float[]{f0,pick};
    }

    b = 0.5f*(fm-fp);
    a = b/a;
    if (a > 1.) {
	    pick = ic+1;
	    return new float[]{fp,pick};
    }

    if (a < -1.) {
	    pick = ic-1;
	    return new float[]{fm,pick};
    }

    if (f0 < 0.5*b*a) {
	    pick = ic;
	    return new float[]{f0,pick};
    }

    f0 -= 0.5*b*a;
    pick=ic+a;
    return new float[]{f0,pick};
  }

  private void forwardTrack(
    float[] path, float[] next, float[][] what)
  {
    float c, d, fc;
    int n1 = next.length;
    int n2 = path.length;
    c = FLT_MAX;
    fc = 0;
    /* minimum at the bottom */
    for (int i1=0; i1 < n1; i1++) {
	    d = next[i1];
	    if (d < c) {
	      c = d;
	      fc = i1;
	    }
    }
    /* coming up */
    for (int i2=0; i2<n2; i2++) {
	    path[i2]=fc;
	    fc = interpolate(fc,i2,what);
    }
  }


  private void backwardTrack(
    float[] path, float[] next, float[][] what)
  {
    float c, d, fc;
    int n1 = next.length;
    int n2 = path.length;
    c = FLT_MAX;
    fc = 0;
    /* minimum at the bottom */
    for (int i1=0; i1 < n1; i1++) {
	    d = next[i1];
	    if (d < c) {
	      c = d;
	      fc = i1;
	    }
    }
    /* coming up */
    for (int i2=n2-1; i2 >= 0; i2--) {
	    path[i2]=fc;
	    fc = interpolate(fc,i2,what);
    }
  }

  private float interpolate(float fc, int i2, float[][] what) {
    int n1 = what[0].length;
    int ic = round(fc-0.5f);
    fc -= ic;
    if (n1-1 <= ic) return what[i2][n1-1];
    if (0 > ic) return what[i2][0];
    fc = what[i2][ic]*(1f-fc)+what[i2][ic+1]*fc;
    return fc;
  }

  private void makeRhsWeightsInline(
    float[][] p23, float[][][] vel, float[][] b, float[][] ws) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k23 = round(p23[i2][i3]);
      k23 = max(0,k23);
      k23 = min(n1-1,k23);
      float w23i = vel[i3][i2][k23];
      w23i *= w23i;
      ws[i3][i2] = w23i;
      b[i3][i2] = p23[i2][i3]*w23i;
    }}
  }


  private void makeRhsWeights(
    float[][] p23, float[][] p32, float[][][] vel, float[][] b, float[][] ws) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k23 = round(p23[i2][i3]);
      int k32 = round(p32[i3][i2]);
      k23 = max(0,k23);
      k23 = min(n1-1,k23);
      k32 = max(0,k32);
      k32 = min(n1-1,k32);
      float w23i = vel[i3][i2][k23];
      float w32i = vel[i3][i2][k32];
      w23i *= w23i;
      w32i *= w32i;
      ws[i3][i2] = w23i+w32i;
      b[i3][i2] = p23[i2][i3]*w23i+p32[i3][i2]*w32i;
    }}
  }


  private void makeRhsWeights(
    float[] p1, float[] p2, float[][] vel, float[] b, float[] ws) 
  {
    int n2 = vel.length;
    int n1 = vel[0].length;
    for (int i1=0; i1<n1; ++i1) {
      int k12 = round(p1[i1]);
      int k22 = round(p2[i1]);
      k12 = max(0,k12);
      k12 = min(n2-1,k12);
      k22 = max(0,k22);
      k22 = min(n2-1,k22);
      float w1i = vel[k12][i1];
      float w2i = vel[k22][i1];
      w1i *= w1i;
      w2i *= w2i;
      ws[i1] = w1i+w2i;
      b[i1] = p1[i1]*w1i+p2[i1]*w2i;
    }
  }

  // Conjugate-gradient operators.
  private static class A1 implements CgSolver.A {
    A1(Smoother1 s1, float[] wp) 
    {
      _s1 = s1;
      _wp = wp;
      float n1 = wp.length;
      _sc = sum(wp)/(n1);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      v1y.zero();
      _s1.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s1.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother1 _s1;
    private float[] _wp;
  }


  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] wp) 
    {
      _s2 = s2;
      _wp = wp;
      float n2 = wp.length;
      float n1 = wp[0].length;
      _sc = 4f*sum(wp)/(n1*n2);
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

  private static void applyLhs(float[] wp, float[] x, float[] y) {
    int n1 = wp.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += wp[i1]*x[i1];
  }


  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    int n2 = wp.length;
    for (int i2=0; i2<n2; ++i2)
      applyLhs(wp[i2],x[i2],y[i2]);
  }

  private static void addAndScale(float sc, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      addAndScale(sc,x[i2],y[i2]);
    }
  }

  private static void addAndScale(float sc, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1) {
      y[i1] += sc*x[i1];
    }
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


  private static class Smoother1 {
    public Smoother1(float sigma) {
      _sigma = sigma;
    }
    public void apply(float[] x) {
      smooth1(_sigma,x);
    }
    public void applyTranspose(float[] x) {
      smooth1(_sigma,x);
    }
    private float _sigma;
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
  private static void smooth1(float sigma, float[] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
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



  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _gate;
  private float _an;


}
