package spv;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;



/**
 * Edge enhancing with shortest path voting
 * @author Xinming Wu
 * @version 2017.07.16
 */

public class ShortestPathVoting2 {
 
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

  public void causalFilter(float a, float[] x, float[] y) {
    int n1 = x.length;
    float b = 1.0f - a;
    float yi = y[0] = x[0];
    for (int i1=1; i1<n1; ++i1) 
      y[i1] = yi = a*yi + b*x[i1];
  }

  public float[][] votingMap(int r, int sig, int[][] seeds, 
    float[][] u1, float[][] u2, float[][] ft, float[][] tt) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    int ns = seeds[0].length;
    float[][] fe = new float[n2][n1];
    float[][] fs = copy(ft);
    RecursiveGaussianFilterP rgfs = new RecursiveGaussianFilterP(2);
    rgfs.apply00(fs,fs);
    for (int is=0; is<ns; ++is) {
      System.out.println("is="+is);
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      int b1 = max(k1-r,0);
      int b2 = max(k2-r,0);
      int e1 = min(k1+r,n1-1);
      int e2 = min(k2+r,n2-1);
      int m1 = min(k1-b1,e1-k1);
      int m2 = min(k2-b2,e2-k2);
      int ri = min(m1,m2);
      float pi = (float)(Math.PI/180.0);
      for (int i2=-ri; i2<=ri; ++i2) {
      for (int i1=-ri; i1<=ri; ++i1) {
        float ti = tt[k2][k1]+90f;
        if(ti>180) ti -= 180;
        ti *= pi;
        float t1 = cos(ti);
        float t2 = sin(ti);
        float rs = sqrt(i1*i1+i2*i2);
        if (rs>0f) {
          float r1 = i1/rs;
          float r2 = i2/rs;
          float st = t1*r1+t2*r2;
          float ct = 1-st*st;
          ct = 1-ct;
          fe[k2+i2][k1+i1] = ct;
        }
      }}
    }
    return fe;
  }


  public float[][] applyEnhance(int r, int sig, int[][] seeds, 
    float[][] u1, float[][] u2, float[][] ft, float[][] tt) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    int ns = seeds[0].length;
    float[][] fe = new float[n2][n1];
    float[][] fs = copy(ft);
    RecursiveGaussianFilterP rgfs = new RecursiveGaussianFilterP(2);
    rgfs.apply00(fs,fs);
    int nr = 2*r+1;
    float[][] u1r = new float[nr][nr];
    float[][] u2r = new float[nr][nr];
    for (int i2=0; i2<nr; ++i2) {
      for (int i1=0; i1<nr; ++i1) {
        float d1 = i1-r;
        float d2 = i2-r;
        float ds = sqrt(d1*d1+d2*d2);
        if(ds>0f) {
          float u1i = -d2/ds;
          float u2i =  d1/ds;
          if(u1i<0f) {
            u1i = -u1i;
            u2i = -u2i;
          }
          u1r[i2][i1] = u1i;
          u2r[i2][i1] = u2i;
        }
    }}
    for (int is=0; is<ns; ++is) {
      System.out.println("is="+is);
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      int b1 = max(k1-r,0);
      int b2 = max(k2-r,0);
      int e1 = min(k1+r,n1-1);
      int e2 = min(k2+r,n2-1);
      int m1 = min(k1-b1,e1-k1);
      int m2 = min(k2-b2,e2-k2);
      int ri = min(m1,m2);
      int c1 = ri*2+1;
      int c2 = ri*2+1;

      float[][] mp = new float[c2][c1];
      float[][] u1s = new float[c2][c1];
      float[][] u2s = new float[c2][c1];

      float[][] fc  = copy(c1,c2,k1-ri,k2-ri,ft);
      float pi = (float)(Math.PI/180.0);
      for (int i2=-ri; i2<=ri; ++i2) {
      for (int i1=-ri; i1<=ri; ++i1) {
        float ti = tt[k2][k1]+90f;
        if(ti>180) ti -= 180;
        ti *= pi;
        float t1 = cos(ti);
        float t2 = sin(ti);
        u1s[i2+ri][i1+ri] = t1;
        u2s[i2+ri][i1+ri] = t2;
        /*
        float rs = sqrt(i1*i1+i2*i2);
        if (rs>0f) {
          float r1 = i1/rs;
          float r2 = i2/rs;
          float st = t1*r1+t2*r2;
          float ct = 1-st*st;
          fc[i2+ri][i1+ri] *= ct;
        }
        */
      }}
      fc  = sub(fc,min(fc));
      fc  = mul(fc,1f/max(fc));
      float[][][] ps = applyEnhance(ri,ri,ri-2,u1s,u2s,fc);
      if(ps!=null) {
      int np = ps.length;
      for (int ip=0; ip<np; ++ip)
        collect(ri,k1,k2,ps[ip],ft,mp);
      }
      for (int i2=0; i2<c2; ++i2) {
      for (int i1=0; i1<c1; ++i1) {
        if(mp[i2][i1]>0f) {
          fe[i2+k2-ri][i1+k1-ri] += mp[i2][i1];
        }
      }}
    }
    return fe;
  }

  public float[][] applyEnhance(int r, int sig, int[][] seeds, 
    float[][] u1, float[][] u2, float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    int ns = seeds[0].length;
    float[][] fe = new float[n2][n1];
    float[][] fs = copy(ft);
    RecursiveGaussianFilterP rgfs = new RecursiveGaussianFilterP(2);
    rgfs.apply00(fs,fs);
    int nr = 2*r+1;
    float[][] u1r = new float[nr][nr];
    float[][] u2r = new float[nr][nr];
    for (int i2=0; i2<nr; ++i2) {
      for (int i1=0; i1<nr; ++i1) {
        float d1 = i1-r;
        float d2 = i2-r;
        float ds = sqrt(d1*d1+d2*d2);
        if(ds>0f) {
          float u1i = -d2/ds;
          float u2i =  d1/ds;
          if(u1i<0f) {
            u1i = -u1i;
            u2i = -u2i;
          }
          u1r[i2][i1] = u1i;
          u2r[i2][i1] = u2i;
        }
    }}
    for (int is=0; is<ns; ++is) {
      //System.out.println("is="+is);
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      int b1 = max(k1-r,0);
      int b2 = max(k2-r,0);
      int e1 = min(k1+r,n1-1);
      int e2 = min(k2+r,n2-1);
      int m1 = min(k1-b1,e1-k1);
      int m2 = min(k2-b2,e2-k2);
      int ri = min(m1,m2);
      int c1 = ri*2+1;
      int c2 = ri*2+1;
      float[][] mp = new float[c2][c1];
      float[][] fc  = copy(c1,c2,k1-ri,k2-ri,ft);
      float[][] u1s = copy(c1,c2,r-ri,r-ri,u1r);
      float[][] u2s = copy(c1,c2,r-ri,r-ri,u2r);
      //float[][] u1s = copy(c1,c2,k1-ri,k2-ri,u1);
      //float[][] u2s = copy(c1,c2,k1-ri,k2-ri,u2);

      float[][][] ps = applyEnhance(ri,ri,ri-2,u1s,u2s,fc);
      if(ps!=null) {
      int np = ps.length;
      for (int ip=0; ip<np; ++ip)
        collect(ri,k1,k2,ps[ip],fs,mp);
      }
      for (int i2=0; i2<c2; ++i2) {
      for (int i1=0; i1<c1; ++i1) {
        if(mp[i2][i1]>0f) {
          fe[i2+k2-ri][i1+k1-ri] += mp[i2][i1];
        }
      }}
    }
    return fe;
  }

  public void collect(int ri, int k1, int k2, 
    float[][] ps, float[][] fs, float[][] mp) {
    int m2 = mp.length;
    int m1 = mp[0].length;
    int np = ps[0].length;
    int n2 = fs.length;
    int n1 = fs[0].length;
    float fa = 0.0f;
    for (int ip=0; ip<np; ++ip) {
      int i1 = round(ps[0][ip])+k1-ri;
      int i2 = round(ps[1][ip])+k2-ri;
      i1 = min(i1,n1-1);
      i1 = max(i1,0);
      i2 = max(i2,0);
      i2 = min(i2,n2-1);
      fa += fs[i2][i1];
    }
    fa /= np;
    for (int ip=0; ip<np; ++ip) {
      int c1 = round(ps[0][ip]);
      int c2 = round(ps[1][ip]);
      for (int d2=-0;d2<=0;d2++){
      for (int d1=-0;d1<=0;d1++){
        int i1 = c1+d1;
        int i2 = c2+d2;
        if(i1>=0&&i1<m1&&i2>=0&&i2<m2)
          mp[i2][i1] = fa;;
      }}
    }
  }
  public float[][][] applyEnhance(
    int k1, int k2, int r, float[][] u1, float[][] u2, float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    int[][] mk = new int[n2][n1];
    float[][] ts = fillfloat(1f,n1,n2);
    ts[k2][k1] = 0.0f;
    float[][] av = clip(0.05f,1.0f,ft);
    float[][] au = mul(av,0.1f);//fillfloat(0.01f,n1,n2);
    EigenTensors2 et = new EigenTensors2(u1,u2,au,av);
    TimeMarker2 tm = new TimeMarker2(n1,n2,et);
    tm.apply(ts,mk);
    int n = 360; // number of sampling angles
    int[] ir = new int[n];
    float[] fr = new float[n];
    float[] r1 = new float[n];
    float[] r2 = new float[n];
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    SincInterpolator si = new SincInterpolator();
    float dt = (float)(2.0*Math.PI/n);
    for (int i=0; i<n; i++) {
      float theta = i*dt; 
      float r1i = k1+r*sin(theta);
      float r2i = k2+r*cos(theta);
      r1[i] = r1i;
      r2[i] = r2i;
      fr[i] = si.interpolate(s1,s2,ts,r1i,r2i);
      ir[i] = i;
    }
    quickIndexSort(fr,ir);
    float e11 = r1[ir[0]];
    float e12 = r2[ir[0]];
    float a0 =  ir[0]*360f/n;
    float[][][] pss = new float[2][2][];
    pss[0] = backTrack(e11,e12,k1,k2,ts);
    for (int i=1; i<n; i++) {
      int ii = ir[i];
      float ai = ii*360f/n;
      float da = abs(ai-a0);
      if(da>180f) da = 360f-da;
      if (da>30f) {
        float e21 = r1[ii];
        float e22 = r2[ii];
        pss[1] = backTrack(e21,e22,k1,k2,ts);
        break;
      }
    }
    return pss;
  }

  public float[][] backTrack(
    float e1, float e2, float b1, float b2, float[][] ts) {
    int n2 = ts.length;
    int n1 = ts[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    float[][] t1 = new float[n2][n1];
    float[][] t2 = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply10(ts,t1);
    rgf.apply01(ts,t2);
    ArrayList<Float> k1 = new ArrayList<Float>();
    ArrayList<Float> k2 = new ArrayList<Float>();
    SincInterpolator si = new SincInterpolator();
    float ti = si.interpolate(s1,s2,ts,e1,e2);
    k1.add(e1);
    k2.add(e2);
    float t1i = 0.0f;
    float t2i = 0.0f;
    float tsi = 100.0f;
    int ct = 0;
    while (ti>=10f&&ct<200) {
      t1i = si.interpolate(s1,s2,t1,e1,e2);
      t2i = si.interpolate(s1,s2,t2,e1,e2);
      tsi = sqrt(t1i*t1i+t2i*t2i);
      t1i /= tsi;
      t2i /= tsi;
      e1 -= t1i;
      e2 -= t2i;
      k1.add(e1);
      k2.add(e2);
      ti = si.interpolate(s1,s2,ts,e1,e2);
      ct++;
    }
    k1.add(b1);
    k2.add(b2);
    int np = k1.size();
    float[][] ps = new float[2][np];
    for (int ip=0; ip<np; ++ip) {
      ps[0][ip] = k1.get(ip);
      ps[1][ip] = k2.get(ip);
    }
    //rgf.apply0X(ps,ps);
    return ps;
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
        ft[i2][i1] = fx[i2][i1];
      }
    }}
    return ft;
  }

  public int[][] pickSeeds(int d, float fm, float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    ArrayList<Float> fs = new ArrayList<Float>();
    ArrayList<Integer> ks = new ArrayList<Integer>();
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
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
    //quickIndexSort(fa,ia);
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

}
