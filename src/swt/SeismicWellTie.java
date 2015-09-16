package swt;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Seismic well ties.
 * 
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.09
 */

public class SeismicWellTie {


  public Object[] updateTimeDepth(Sampling s1, Sampling s2, Sampling s3, 
    WellLog[] logs, float[][][] gx) 
  {
    // flatten real seismic traces
    int smax = 40;
    float epow = 0.125f;
    int nl = logs.length;
    int n1 = s1.getCount();
    double f1 = s1.getFirst();
    double d1 = s1.getDelta();
    float[][] fx = new float[nl][n1];
    SynSeis.Model[] models = new SynSeis.Model[nl];
    for (int il=0; il<nl; ++il) {
      WellLog log = logs[il];
      SynSeis.Model md = SynSeis.getModel(log);
      int i2 = s2.indexOfNearest(md.x2);
      int i3 = s3.indexOfNearest(md.x3);
      fx[il] = gx[i3][i2];
      models[il] = md;
    }
    float[][][] frs = seisFlatten(smax,epow,fx);

    float[][] r1 = null;
    float[][] r2 = null;
    float[][] r3 = null;
    float[][] r6 = null;
    Sampling[] r4 = null;
    Sampling[] r5 = null;
    Sampling[] r7 = null;

    float umin = -0.30f;
    float umax =  0.30f;
    float rmin = -0.15f;
    float rmax =  0.15f;
    float dmin =  (float)(10*d1);
    int n = 2;
    for (int k=0; k<n; ++k) {
      umin += k*abs(umin)/(2f*n);
      umax -= k*abs(umax)/(2f*n);
    // flatten synthetic seismograms
      boolean simple = true; //use simple model
      double[][] ndfwx = new double[nl][3];
      double[][] ndfwu = new double[nl][3];
      float[][][] wrs = synsFlatten(models,simple,smax,epow,logs,ndfwx,ndfwu);
      Sampling[] swx = new Sampling[nl];
      Sampling[] swu = new Sampling[nl];
      for(int il=0; il<nl; ++il) {
        swx[il] = new Sampling((int)ndfwx[il][0],ndfwx[il][1],ndfwx[il][2]);
        swu[il] = new Sampling((int)ndfwu[il][0],ndfwu[il][1],ndfwu[il][2]);
      }

    // align flattened syns to flattened seis
      float[][] wu = wrs[2];
      float[][] fu = frs[0];
      Sampling[] sw = new Sampling[nl];
      for (int il=0; il<nl; ++il) {
        double[] ndfi = ndfwu[il];
        sw[il] = new Sampling((int)ndfi[0],ndfi[1],ndfi[2]);
      }
      float[] ndfwm = new float[3];
      float[] ws = synsToSeisShift(sw,s1,umin,umax,rmin,rmax,dmin,wu,fu,ndfwm); 
      Sampling sws = new Sampling((int)ndfwm[0],ndfwm[1],ndfwm[2]);
      System.out.println("wgs="+sum(abs(ws)));

    // update time depth
      float[][] wfs = wrs[1]; // shifts for syns flattening
      float[][] ffs = frs[1]; // shifts for seis flattening
      double[][] wts = new double[nl][];
      SincInterpolator si = new SincInterpolator();
      si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
      for (int il=0; il<nl; ++il) {
        float[] fxt = new float[n1];
        float[] fxu = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          double fxti = f1+i1*d1;
          double fxui = fxti+ffs[il][i1]*d1;
          fxt[i1] = (float)fxti;
          fxu[i1] = (float)fxui;
          if(i1>0 && fxu[i1]<=fxu[i1-1]) {
            fxu[i1] = fxu[i1-1]+0.000001f;
          }
        }

        double fwt = ndfwx[il][2];
        double dwt = ndfwx[il][1];
        int nwt = (int)ndfwx[il][0];
        wts[il] = new double[nwt];
        CubicInterpolator.Method method = CubicInterpolator.Method.LINEAR;
        CubicInterpolator ci = new CubicInterpolator(method,fxu,fxt);
        for (int iwt=0; iwt<nwt; ++iwt) {
          double wti = fwt+iwt*dwt;
          double wui = wti+wfs[il][iwt]*dwt;   // shift to RGT of syns
          double fui = wui+si.interpolate(sws,ws,wui); // shift to RGT of seis
          double fti = ci.interpolate((float)fui);    // convert to time of seis
          wts[il][iwt] = fti;
        }
      }

    // compute syns in time of seis 
      float[][] wft = new float[nl][];
      Sampling[] sft = new Sampling[nl];
      for (int il=0; il<nl; ++il) {
        int nt = wts[il].length;
        float[] x = new float[nt];
        float[] y = new float[nt];
        float[] z = new float[nt];
        for (int it=0; it<nt; ++it) {
          y[it] = wrs[0][il][it];
          x[it] = (float)wts[il][it];
          z[it] = (float)(wts[il][0]+it*d1);
          if(it>0 && x[it]<=x[it-1]) {
            x[it] = x[it-1]+0.000001f;
          }
        }
        CubicInterpolator.Method method = CubicInterpolator.Method.LINEAR;
        CubicInterpolator ci = new CubicInterpolator(method,x,y);
        wft[il] = ci.interpolate(z);
        sft[il] = new Sampling(nt,d1,wts[il][0]);
      }
      r1 = fx;
      r2 = wrs[0];
      r3 = wft;
      r4 = swx;
      r5 = sft;
      r6 = wrs[2];
      r7 = swu;
      updateModels(models,swx,wts);
    }
    return new Object[]{r1,r2,r3,r4,r5,r6,r7};
  }

  private void updateModels(
    SynSeis.Model[] models, Sampling[] swx, double[][] wts) 
  {
    int nl = models.length;
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int il=0; il<nl; ++il) {
      SynSeis.Model model = models[il];
      float[] t1 = model.t;
      float[] t2 = copy(t1);
      /*
      int nt = t1.length;
      float[] t2 = new float[nt];
      int nx = swx[il].getCount();
      double dx = swx[il].getDelta();
      double fx = swx[il].getFirst();
      float[] fa = new float[nx];
      for (int ix=0; ix<nx; ++ix) {
        fa[ix] = (float)wts[il][ix];
      }
      si.interpolate(nx,dx,fx,fa,nt,t1,t2);
      */
      double dt = wts[il][0]-swx[il].getFirst();
      t2 = add(t1,(float)dt);
      model.updateTime(t2);
      models[il] = model;
    }
  }

  public float[][] synsToSeis(Sampling[] sw, Sampling sg, 
    float smin, float smax, float rmin, float rmax, 
    float dmin, float[][] wu, float[][] gu, float[][] ndf) 
  {
    int nl = sw.length;
    double[] fs = new double[nl];
    double[] ls = new double[nl];
    for (int il=0; il<nl; ++il) {
      ls[il] = sw[il].getLast();
      fs[il] = sw[il].getFirst();
    }
    double fw = min(fs);
    double lw = max(ls);
    double dw = sw[0].getDelta();
    int nw = (int)((lw-fw)/dw)+1; 
    Sampling swm = new Sampling(nw,dw,fw);
    DynamicWarpingX dwx = new DynamicWarpingX(smin,smax,swm);
    dwx.setStrainLimits(rmin,rmax);
    dwx.setSmoothness(dmin);
    float[][] es = dwx.computeErrorsX(sw,wu,sg,gu);
    float[]   ws = dwx.findShifts(es);
    float[][] wsu = applyShifts(_vnull,sw,swm,ws,wu,ndf);
    System.out.println("nw="+nw);
    System.out.println("ns="+ws.length);
    return wsu;
  }


  public float[] synsToSeisShift(Sampling[] sw, Sampling sg, 
    float smin, float smax, float rmin, float rmax, 
    float dmin, float[][] wu, float[][] gu, float[] ndfwm) 
  {
    int nl = sw.length;
    double[] fs = new double[nl];
    double[] ls = new double[nl];
    for (int il=0; il<nl; ++il) {
      ls[il] = sw[il].getLast();
      fs[il] = sw[il].getFirst();
    }
    double fw = min(fs);
    double lw = max(ls);
    double dw = sw[0].getDelta();
    int nw = (int)((lw-fw)/dw)+1; 
    Sampling swm = new Sampling(nw,dw,fw);
    ndfwm[0] = nw;
    ndfwm[1] = (float)dw;
    ndfwm[2] = (float)fw;
    DynamicWarpingX dwx = new DynamicWarpingX(smin,smax,swm);
    dwx.setStrainLimits(rmin,rmax);
    dwx.setSmoothness(dmin);
    float[][] es = dwx.computeErrorsX(sw,wu,sg,gu);
    return dwx.findShifts(es);
  }

  public float[][][] seisFlatten(int smax, float epow, float[][] gx) {
    MultiTraceWarping mtw = new MultiTraceWarping();
    mtw.setMaxShift(smax);
    mtw.setPowError(epow);
    mtw.setNullValue(max(gx)+10f);
    float[][] gs = mtw.findShifts(gx);
    float[][] gu = mtw.applyShiftsX(gx,gs);
    System.out.println("gs="+sum(abs(gs)));
    return new float[][][]{gu,gs};
  }

  public float[][][] synsFlatten(SynSeis.Model[] models, boolean simple,  
    int smax, float epow, WellLog[] logs,double[][] ndfx,double[][] ndfu) 
  {
    double[] ndf = new double[3];
    float[][] wx = computeSyns(models,simple,logs,ndf);
    Sampling st = new Sampling((int)ndf[0],ndf[1],ndf[2]);
    MultiTraceWarping mtw = new MultiTraceWarping();
    mtw.setMaxShift(smax);
    mtw.setPowError(epow);
    float[][] ws = mtw.findShifts(wx);
    float[][] wu = mtw.applyShiftsW(wx,ws);
    float[][] vu = getValidSamples(st,wu,ndfu);
    float[][][] xs = getValidSamples(st,wx,ws,ndfx);
    System.out.println("ws="+sum(abs(ws)));
    return new float[][][]{xs[0],xs[1],vu};
  }

  public float[][] synsFlatten(
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

  public float[][] computeSyns(
    boolean simple, WellLog[] logs, float[][] ndft, float[][] ndfz) 
  {
    float fpeak = 35f;
    float q = 100.0f;
    float ds = 0.002f;
    int nl = logs.length;
    float[][] sa = new float[nl][];
    for (int il=0; il<nl; ++il) {
      WellLog log = logs[il];
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
    }
    return sa;
  }

  public float[][] computeSyns(SynSeis.Model[] models,
    boolean simple, WellLog[] logs, double[] ndf) {
    float fpeak = 35f;
    float q = 100.00f;
    double ds = 0.002f;
    int nl = logs.length;
    double fs =  FLT_MAX;
    double ls = -FLT_MAX;
    for (int il=0; il<nl; ++il) {
      WellLog log=logs[il];
      SynSeis.Model md = (models!=null)?models[il]:SynSeis.getModel(log);
      if (fs>md.tmin()){fs=md.tmin();}
      if (ls<md.tmax()){ls=md.tmax();}
    }
    fs = round(fs/ds)*ds;
    ls = round(ls/ds)*ds;
    int ns = (int)((ls-fs)/ds)+1;
    ndf[0] = ns;
    ndf[1] = ds;
    ndf[2] = fs;
    float[][] sa = fillfloat(_vnull,ns,nl);
    for (int il=0; il<nl; ++il) {
      WellLog log=logs[il];
      SynSeis.Model md = (models!=null)?models[il]:SynSeis.getModel(log);
      double fsi = round(md.tmin()/ds)*ds;
      double lsi = round(md.tmax()/ds)*ds;
      int nsi = (int)((lsi-fsi)/ds)+1;
      Sampling ssi = new Sampling(nsi,ds,fsi);
      float[] fsyn = new float[nsi];
      if (simple) {
        fsyn = SynSeis.makeSimpleSeismogram(md,fpeak,ssi);
      } else {
        fsyn = SynSeis.makeBetterSeismogram(md,q,fpeak,ssi);
      }
      fsyn = normalize(fsyn);
      int j1 = (int)((fsi-fs)/ds);
      copy(nsi,0,fsyn,j1,sa[il]);
    }
    return sa;
  }


  public float[][] getValidSamples(
    Sampling sz, float[][] sa, double[][] ndf) {
    int nl = sa.length;
    int nz = sa[0].length;
    float[][] va = new float[nl][];
    double dz = sz.getDelta();
    double fz = sz.getFirst();
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
    return va;
  }

  public float[][][] getValidSamples(
    Sampling sz, float[][] sa, float[][] sb, double[][] ndf) {
    int nl = sa.length;
    int nz = sa[0].length;
    float[][] va = new float[nl][];
    float[][] vb = new float[nl][];
    double dz = sz.getDelta();
    double fz = sz.getFirst();
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
      vb[il] = new float[nzi];
      copy(nzi,jz,sa[il],0,va[il]);
      copy(nzi,jz,sb[il],0,vb[il]);
    }
    return new float[][][]{va,vb};
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

  /*
  Given a uniformly sampled u(s) such that s+u(s) increases monotonically, 
  computes a uniformly sampled r(t) = u(t-r(t)) such that t-r(t) increases 
  monotonically. Uses the following 3-step algorithm:
  (1) computes t(s) = s+u(s)
  (2) computes s(t) = by inverse linear interpolation of t(s)
  (3) computes r(t) = t-s(t)
  Returns the sampling of time st and the sequence of sampled r(t).
  */

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



  private float[] normalize(float[] x) {
    float sigma=100.0f;
    float[] y = mul(x,x);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.apply(y,y);
    return div(x,sqrt(y));
  }

  private float _vnull = -999.2500f;

}

