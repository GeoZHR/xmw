package swt;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Seismic well ties.
 * 
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.09
 */

public class SeismicWellTie {


  public float[][][] logsToArray(WellLog[] logs, float[] zs) {
    int nc = 2;
    float dz = 0.0f;
    int nl = logs.length;
    float fk = 0.0003048f; // 1 ft = 0.0003048 km
    float[] fzs = new float[nl];
    float[] lzs = new float[nl];
    for (int il=0; il<nl; ++il){
      float[] z = logs[il].z;
      float[] d = logs[il].d;
      float[] v = logs[il].v;
      int nz = z.length;
      int jz = 0;
      int lz = nz-1;
      while (jz<nz && d[jz]==_vnull && v[jz]==_vnull)
        ++jz;
      while (lz>=0 && d[lz]==_vnull && v[lz]==_vnull)
        --lz;
      fzs[il] = z[jz];
      lzs[il] = z[lz];
      dz = z[1]-z[0];
    }
    float fz = min(fzs);
    float lz = max(lzs);
    int nz = round((lz-fz)/dz)+1;
    float[][][] wa = fillfloat(_vnull,nz,nl,nc);
    for (int il=0; il<nl; ++il) {
      float[] z = logs[il].z;
      float[] d = logs[il].d;
      float[] v = logs[il].v;
      int jz = 0;
      while (jz<z.length && d[jz]==_vnull && v[jz]==_vnull)
        ++jz;
      int jw = round((z[jz]-fz)/dz);
      //for( ; jz<z.length; jz++) {
      for( ; jw<nz && jz<z.length; jw++) {
        wa[0][il][jw] = d[jz];
        wa[1][il][jw] = v[jz];
        jz ++;
      }
    }
    dz *= fk;
    fz *= fk;
    zs[0] = nz;
    zs[1] = dz;
    zs[2] = fz;
    return wa;
  }

  public float[][] computeSyns(
    boolean simple, WellLog[] logs, float[][] ndft, float[][] ndfz) 
  {
    float fpeak = 35f;
    float q = 100.0f;
    int il = 0;
    float ds = 0.002f;
    int nl = logs.length;
    float[][] sa = new float[nl][];
    for (WellLog log:logs) {
      SynSeis.Model md = SynSeis.getModel(log);
      float fsi = md.tmin();
      float lsi = md.tmax();
      fsi = round(fsi/ds)*ds;
      lsi = round(lsi/ds)*ds;
      int nsi = round((lsi-fsi)/ds)+1;
      Sampling ssi = new Sampling(nsi,ds,fsi);
      Sampling szi = md.sz;
      ndft[il][0] = nsi;
      ndft[il][1] = ds;
      ndft[il][2] = fsi;
      ndfz[il][0] = szi.getCount();
      ndfz[il][1] = (float)szi.getDelta();
      ndfz[il][2] = (float)szi.getFirst();
      if (simple) {
        sa[il] = SynSeis.makeSimpleSeismogram(md,fpeak,ssi);
      } else {
        sa[il] = SynSeis.makeBetterSeismogram(md,q,fpeak,ssi);
      }
      il++;
    }
    return sa;
  }

  public float[][] computeSyns(boolean simple, WellLog[] logs, float[] ndf) {
    float ds = 0.002f;
    float fpeak = 35f;
    float q = 100.0f;
    int nl = logs.length;
    float fs =  FLT_MAX;
    float ls = -FLT_MAX;
    for (WellLog log:logs) {
      SynSeis.Model md = SynSeis.getModel(log);
      if (fs>md.tmin()){fs=md.tmin();}
      if (ls<md.tmax()){ls=md.tmax();}
    }
    fs = round(fs/ds)*ds;
    ls = round(ls/ds)*ds;
    int ns = round((ls-fs)/ds)+1;
    ndf[0] = ns;
    ndf[1] = ds;
    ndf[2] = fs;
    int il = 0;
    float[][] sa = fillfloat(_vnull,ns,nl);
    for (WellLog log:logs) {
      SynSeis.Model md = SynSeis.getModel(log);
      float fsi = round(md.tmin()/ds)*ds;
      float lsi = round(md.tmax()/ds)*ds;
      int nsi = round((lsi-fsi)/ds)+1;
      Sampling ssi = new Sampling(nsi,ds,fsi);
      float[] fsyn = new float[nsi];
      if (simple) {
        fsyn = SynSeis.makeSimpleSeismogram(md,fpeak,ssi);
      } else {
        fsyn = SynSeis.makeBetterSeismogram(md,q,fpeak,ssi);
      }
      fsyn = normalize(fsyn);
      int j1 = round((fsi-fs)/ds);
      copy(nsi,0,fsyn,j1,sa[il]);
      il++;
    }
    return sa;
  }

  public float[][] getValidSamples(Sampling sz, float[][] sa, float[][] ndf) {
    int nl = sa.length;
    int nz = sa[0].length;
    float[][] va = new float[nl][];
    float dz = (float)sz.getDelta();
    float fz = (float)sz.getFirst();
    for (int il=0; il<nl; ++il){
      int jz = 0;
      int lz = nz-1;
      while (jz<nz && sa[il][jz]==_vnull)
        ++jz;
      while (lz>=0 && sa[il][lz]==_vnull)
        --lz;
      int nzi = lz-jz+1;
      ndf[il][0] = nzi;
      ndf[il][1] = dz;
      ndf[il][2] = jz*dz+fz;
      va[il] = new float[nzi];
      copy(nzi,jz,sa[il],0,va[il]);
    }
    System.out.println("mva="+min(va));
    return va;
  }

  public float[][] applyShifts(float vnull, Sampling[] sfs,
    Sampling sfm, float[] u, float[][] f, float[][] ndfs) {
    int n2 = f.length;
    float[][] us = invertShifts(vnull,sfs,sfm,u,f,ndfs);
    SincInterpolator si = new SincInterpolator();
    float[][] fs = new float[n2][];
    for (int i2=0; i2<n2; ++i2) {
      int nt = (int)ndfs[i2][0];
      float dt = ndfs[i2][1];
      float ft = ndfs[i2][2];
      float[] s = new float[nt];
      for (int i1=0; i1<nt; ++i1) {
        s[i1] = ft+i1*dt-us[i2][i1];
      }
      fs[i2] = new float[nt];
      int ns2    = sfs[i2].getCount();
      double ds2 = sfs[i2].getDelta();
      double fs2 = sfs[i2].getFirst();
      si.interpolate(ns2,ds2,fs2,f[i2],nt,s,fs[i2]);
    }
    return fs;
  }

  public float[][] invertShifts(float vnull, Sampling[] sfs,
    Sampling sfm, float[] u, float[][] f, float[][] ndfs) 
  {
    int n2 = f.length;
    float dt = (float)sfm.getDelta();
    float[][] us = new float[n2][];
    for (int i2=0; i2<n2; ++i2) {
      SincInterpolator si = new SincInterpolator();
      float fs = (float)sfs[i2].getFirst();
      float ft = fs+si.interpolate(sfm,u,fs);
      ft = round(ft/dt)*dt;
      int nt = sfs[i2].getCount();
      ndfs[i2][0] = nt;
      ndfs[i2][1] = dt;
      ndfs[i2][2] = ft;
      float[] t = new float[nt];
      for (int i1=0; i1<nt; ++i1) {
        float x1 = (float)sfs[i2].getValue(i1);
        t[i1] = x1+si.interpolate(sfm,u,x1);
      }
      us[i2] = new float[nt];
      float[] s = new float[nt];
      Sampling ss = new Sampling(nt,dt,fs);
      Sampling st = new Sampling(nt,dt,ft);
      InverseInterpolator ii = new InverseInterpolator(ss,st);
      ii.invert(t,s);
      for (int i1=0; i1<nt; ++i1) {
        us[i2][i1] = (float)st.getValue(i1)-s[i1];
      }
    }
    return us;
  }

  public float[] stack(Sampling sz, float[][] sa, float[] ndf) {
    int nl = sa.length;
    int nz = sa[0].length;
    float[] ss = new float[nz];
    float[] sc = new float[nz];
    for (int il=0; il<nl; ++il) {
    for (int iz=0; iz<nz; ++iz) {
      float si = sa[il][iz];
      if (si!=_vnull) {
        ss[iz] += si;
        sc[iz] += 1f;
      }
    }}
    for (int iz=0; iz<nz; ++iz) {
      float ssi = ss[iz];
      float sci = sc[iz];
      if (sci!=0.0f) {
        ss[iz] = ssi/sci;
      } else {
        ss[iz] = _vnull;
      }
    }
    float dz = (float)sz.getDelta();
    float fz = (float)sz.getFirst();
    int jz = 0;
    int lz = nz-1;
    while (jz<nz && ss[jz]==_vnull)
      ++jz;
    while (lz>=0 && ss[lz]==_vnull)
      --lz;
    int nzi = lz-jz+1;
    ndf[0] = nzi;
    ndf[1] = dz;
    ndf[2] = jz*dz+fz;
    float[] sv = new float[nzi];
    for (int ik=0; ik<nzi; ++ik) {
      sv[ik] = ss[ik+jz];
    }
    return sv;
  }

  public float[][] getTimeDepth(WellLog[] logs) {
    int il = 0;
    int nl = logs.length;
    float[][] td = new float[nl][];
    for (WellLog log:logs) {
      SynSeis.Model md = SynSeis.getModel(log);
      td[il] = md.t;
      il++;
    }
    return td;
  }

  public float[][] applySynsFlat(
    Sampling sz, Sampling sm[], Sampling[] st, float[][] ndfw, 
    float[][] ss, float[][] tz, float[][] gt) 
  {
    int nl = sm.length;
    float[][] gw = new float[nl][];
    for (int il=0; il<nl; ++il) {
      int nz    = sz.getCount();
      double fz = sz.getFirst();
      double dz = sz.getDelta();
      double lm = sm[il].getLast();
      double fm = sm[il].getFirst();
      double dm = sm[il].getDelta();
      int nw = 0;
      double fw = 0.0;
      ArrayList<Float> gl = new ArrayList<Float>();
      SincInterpolator sit = new SincInterpolator();
      SincInterpolator sig = new SincInterpolator();
      for (int ik=0; ik<nz; ++ik) {
        float s = ss[il][ik]; 
        double w = ik*dz+fz;
        double z = w-s*dz;
        if(s!=_vnull && z>fm+dm && z<lm-2*dm) {
          if(nw==0){fw=w;}
          float t = sit.interpolate(sm[il],tz[il],z);
          float g = sig.interpolate(st[il],gt[il],t);
          gl.add(g);
          nw++;
        }
      }
      gw[il] = new float[nw];
      for (int iw=0; iw<nw; ++iw) {
        gw[il][iw] = gl.get(iw);
      }
      ndfw[il][0] = nw;
      ndfw[il][1] = (float)dz;
      ndfw[il][2] = (float)fw;
    }
    return gw;
  }

  private float[] normalize(float[] x) {
    float sigma=100.0f;
    float[] y = mul(x,x);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.apply(y,y);
    return div(x,sqrt(y));
  }

  private float _vnull = -999.2500f;

}

