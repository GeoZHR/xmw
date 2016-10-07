package swt;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Seismic well ties.
 * 
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.09
 */

public class SeismicWellTie {

  public float[][] getSyns(
    boolean simple, WellLog[] logs, double[][] ndfs) 
  {
    double[] ndf = new double[3];
    float[][] wx = computeSyns(simple,logs,ndf);
    Sampling sx = new Sampling((int)ndf[0],ndf[1],ndf[2]);
    float[][] vx = getValidSamples(sx,wx,ndfs);
    return vx;
  }

  public float[][] getSyns(
    boolean simple, SynSeis.Model[] mds, double[][] ndfs) 
  {
    double[] ndf = new double[3];
    float[][] wx = computeSyns(simple,mds,ndf);
    Sampling sx = new Sampling((int)ndf[0],ndf[1],ndf[2]);
    float[][] vx = getValidSamples(sx,wx,ndfs);
    return vx;
  }

  public SynSeis.Model[] updateTimeDepthS(
    Sampling s1, Sampling s2, Sampling s3, WellLog[] lgs, float[][][] gx) 
  {
    int nl = lgs.length;
    SynSeis.Model[] mds = new SynSeis.Model[nl];
    for (int il=0; il<nl; ++il) 
      mds[il]=updateTimeDepthS(s1,s2,s3,lgs[il],gx);
    return mds;
  }

  public SynSeis.Model updateTimeDepthS(
    Sampling s1, Sampling s2, Sampling s3, WellLog lg, float[][][] gx) 
  {
    int k = 0;
    float wgs = 1.0f;
    double d1 = s1.getDelta();
    float dmin= (float)(20*d1);
    float umin=-0.10f, umax=0.20f;
    float rmin=-0.05f, rmax=0.05f;
    SynSeis.Model md = SynSeis.getModel(lg);
    int i2 = s2.indexOfNearest(md.x2);
    int i3 = s3.indexOfNearest(md.x3);
    while(k<5 && wgs>0.0f) {
      double[] ndf = new double[3];
      float[] wx = computeSyn(true,d1,md,ndf); 
      Sampling sw = new Sampling((int)ndf[0],ndf[1],ndf[2]);
      float[] ws = synsToSeisShift(sw,s1,umin,umax,rmin,rmax,dmin,wx,gx[i3][i2]); 
      wgs = sum(abs(ws));
      System.out.println("wgs="+wgs);
      updateModel(md,sw,ws);
      k++;
    }
    return md;
  }

  /**
   * Simultaneously update time-depth functions of multiple well logs.
   * @param s1 vertical sampling of the real seismic image.
   * @param logs array of logs to be tied to seismic 
   * @param fx array of seismic traces close to wells.
   * @param fs array of shifts that flatten seismic traces.
   * @param fu array of flattened seismic traces.
   * @return array for the updated models of the multiple well logs.
   */
  public Object[] updateTimeDepthM(
    int niter, Sampling s1, WellLog[] logs, 
    float[][] fx, float[][] fs, float[][] fu) 
  {
    int nl = logs.length;
    double d1 = s1.getDelta();
    SynSeis.Model[] mds = new SynSeis.Model[nl];
    for (int il=0; il<nl; ++il) 
      mds[il] = SynSeis.getModel(logs[il]);

    int k = 0;
    float wgs = 1.0f;
    float dmin= (float)(20*d1);
    float umin=-0.10f, umax=0.30f;
    float rmin=-0.05f, rmax=0.05f;
    float[][] wws = new float[nl][];
    float[][] wwu = new float[nl][];
    float[][] wsu = new float[nl][];
    Sampling[] swwx = new Sampling[nl];
    Sampling[] swwu = new Sampling[nl];
    Sampling[] swsu = new Sampling[nl];
    while(k<niter && wgs>0.0f) {
      // flatten synthetic seismograms
      int smax = 40;
      float epow = 0.0625f;
      boolean simple = true; //use simple model
      double[][] ndfwx = new double[nl][3];
      double[][] ndfwu = new double[nl][3];
      float[][][] wrs = null;
      if(k==0)
        wrs=synsFlatten(simple,smax,epow,mds,ndfwx,ndfwu);
      else
        wrs=synsFlatten(simple,smax,epow,s1,mds,fx,ndfwx,ndfwu);
      wws = wrs[1]; // shifts for syns flattening
      wwu = wrs[2]; // flattened syns
      for(int il=0; il<nl; ++il) {
        swwx[il] = new Sampling((int)ndfwx[il][0],ndfwx[il][1],ndfwx[il][2]);
        swwu[il] = new Sampling((int)ndfwu[il][0],ndfwu[il][1],ndfwu[il][2]);
      }

      // align flattened syns to flattened seis
      float[] ndfwm = new float[3];
      float[] wus = synsToSeisShift(swwu,s1,umin,umax,rmin,rmax,dmin,wwu,fu,ndfwm); 
      Sampling sws = new Sampling((int)ndfwm[0],ndfwm[1],ndfwm[2]);
      wgs = sum(abs(wus));
      System.out.println("wgs="+wgs);
      float[][] ndf = new float[nl][3];
      wsu = applyShifts(_vnull,swwu,sws,wus,wwu,ndf);
      for(int il=0; il<nl; ++il)
        swsu[il] = new Sampling((int)ndf[il][0],ndf[il][1],ndf[il][2]);

      // update time depth
      double[][] wts = compositeShifts(s1,sws,wus,ndfwx,wws,fs);
      updateModels(mds,swwx,wts);
      k++;
    }
    return new Object[]{swsu,wsu,mds};
  }

  public double[][] compositeShifts(
    Sampling s1, Sampling swu, float[] wus, 
    double[][] ndfw, float[][] wxs, float[][] fxs) 
  {
    int nl = ndfw.length;
    int n1 = s1.getCount();
    double d1 = s1.getDelta();
    double f1 = s1.getFirst();
    double[][] wts = new double[nl][];
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int il=0; il<nl; ++il) {
      float[] fxt = new float[n1];
      float[] fxu = new float[n1];
      for (int i1=0; i1<n1; ++i1) {
        double fxti = f1+i1*d1;
        double fxui = fxti+fxs[il][i1]*d1;
        fxt[i1] = (float)fxti;
        fxu[i1] = (float)fxui;
        if(i1>0 && fxu[i1]<=fxu[i1-1]) {
          fxu[i1] = fxu[i1-1]+0.000001f;
        }
      }
      double fwt = ndfw[il][2];
      double dwt = ndfw[il][1];
      int nwt = (int)ndfw[il][0];
      wts[il] = new double[nwt];
      CubicInterpolator.Method method = CubicInterpolator.Method.LINEAR;
      CubicInterpolator ci = new CubicInterpolator(method,fxu,fxt);
      for (int iwt=0; iwt<nwt; ++iwt) {
        double wti = fwt+iwt*dwt;
        double wui = wti+wxs[il][iwt]*dwt;   // shift to RGT of syns
        double fui = wui+si.interpolate(swu,wus,wui); // shift to RGT of seis
        double fti = ci.interpolate((float)fui);    // convert to time of seis
        wts[il][iwt] = fti;
      }
    } 
    return wts;
  }

  public float[][][][] getSamples(Sampling s1, WellLog[] logs) {
    int nl = logs.length;
    SynSeis.Model[] models = new SynSeis.Model[nl];
    for (int il=0; il<nl; ++il) {
      models[il] = SynSeis.getModel(logs[il]);
      System.out.println("id="+logs[il].id);
      System.out.println("x2="+models[il].x2);
      System.out.println("x3="+models[il].x3);
    }
    return getSamples(s1,models);
  }

  public float[][][][] getSamples(Sampling s1, SynSeis.Model[] models) {
    int nl = models.length;
    float f1 = (float)s1.getFirst();
    float d1 = (float)s1.getDelta();
    float[][][] sds = new float[4][nl][];
    float[][][] svs = new float[4][nl][];
    float[][][] sps = new float[5][nl][];
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = models[il];
      Sampling sz = md.sz;
      float ft = f1+round((min(md.t)-f1)/d1)*d1;
      int nt = round((max(md.t)-min(md.t))/d1)+1;
      Sampling st = new Sampling(nt,d1,ft);
      InverseInterpolator ii = new InverseInterpolator(sz,st);
      float[] ts = md.t;
      float[] zs = new float[nt];
      for (int i=0; i<md.nz; ++i) {
        if(i>0 && ts[i]<=ts[i-1]) {
          ts[i] = ts[i-1]+0.000001f;
        }
      }
      ii.invert(ts,zs);
      sds[0][il] = new float[nt];
      sds[1][il] = new float[nt];
      sds[2][il] = new float[nt];
      sds[3][il] = new float[nt];
      svs[0][il] = new float[nt];
      svs[1][il] = new float[nt];
      svs[2][il] = new float[nt];
      svs[3][il] = new float[nt];
      sps[0][il] = new float[nt];
      sps[1][il] = new float[nt];
      sps[2][il] = new float[nt];
      sps[3][il] = new float[nt];
      sps[4][il] = new float[nt];
      float[] fvs = md.v;
      float[] fds = md.d;
      for (int it=0; it<nt; ++it) {
        int itm = it-1;
        int itp = it+1;
        itm = max(itm,0);
        itp = min(itp,nt-1);
        float zm = zs[itm];
        float zp = zs[itp];
        int im = sz.indexOfNearest(zm);
        int ip = sz.indexOfNearest(zp);
        int np = ip-im+1;
        float[] fv = copy(np,im,fvs); 
        float[] fd = copy(np,im,fds); 
        MedianFinder mf = new MedianFinder(np);
        float fvm = mf.findMedian(fv);
        float fdm = mf.findMedian(fd);
        svs[0][il][it] = fvm;
        svs[1][il][it] = (float)st.getValue(it);
        svs[2][il][it] = (float)md.x2;
        svs[3][il][it] = (float)md.x3;
        sds[0][il][it] = fdm;
        sds[1][il][it] = (float)st.getValue(it);
        sds[2][il][it] = (float)md.x2;
        sds[3][il][it] = (float)md.x3;
        sps[1][il][it] = (float)st.getValue(it);
        sps[2][il][it] = (float)md.x2;
        sps[3][il][it] = (float)md.x3;
        float pi = fvm*fdm;
        sps[0][il][it] = 0.5f*log(pi);
        if (it>0)
          sps[4][il][it] = sps[0][il][it]-sps[0][il][it-1];
      }
    }
    return new float[][][][]{svs,sds,sps};
  }

  public float[][][][] getSamples(
    Sampling s1, WellLog[] logs, SynSeis.Model[] models) 
  {
    int nl = logs.length;
    float f1 = (float)s1.getFirst();
    float d1 = (float)s1.getDelta();
    float[][][] spw = new float[4][nl][];
    float[][][] sps = new float[4][nl][];
    for (int il=0; il<nl; ++il) {
      SynSeis.Model mds = models[il];
      SynSeis.Model mdw = SynSeis.getModel(logs[il]);
      Sampling szw = mdw.sz;
      Sampling szs = mds.sz;
      float ftw = f1+round((min(mdw.t)-f1)/d1)*d1;
      int ntw = round((max(mdw.t)-min(mdw.t))/d1)+1;
      Sampling stw = new Sampling(ntw,d1,ftw);
      float fts = f1+round((min(mds.t)-f1)/d1)*d1;
      int nts = round((max(mds.t)-min(mds.t))/d1)+1;
      Sampling sts = new Sampling(nts,d1,fts);
      InverseInterpolator wi = new InverseInterpolator(szw,stw);
      InverseInterpolator si = new InverseInterpolator(szs,sts);
      float[] zw = new float[ntw];
      float[] zs = new float[nts];
      float[] tw = mdw.t;
      float[] ts = mds.t;
      for (int i=0; i<mds.nz; ++i) {
        if(i>0 && ts[i]<=ts[i-1]) {
          ts[i] = ts[i-1]+0.000001f;
        }
      }
      wi.invert(tw,zw);
      si.invert(ts,zs);
      sps[0][il] = new float[nts];
      sps[1][il] = new float[nts];
      sps[2][il] = new float[nts];
      sps[3][il] = new float[nts];
      spw[0][il] = new float[ntw];
      spw[1][il] = new float[ntw];
      spw[2][il] = new float[ntw];
      spw[3][il] = new float[ntw];
      float[] fvs = mds.v;
      float[] fvw = mdw.v;
      for (int its=0; its<nts; ++its) {
        int itm = its-1;
        int itp = its+1;
        itm = max(itm,0);
        itp = min(itp,nts-1);
        float zm = zs[itm];
        float zp = zs[itp];
        float zi = zs[its];
        int im = szs.indexOfNearest(zm);
        int ii = szs.indexOfNearest(zi);
        int ip = szs.indexOfNearest(zp);
        int ib = min(im+1,ii);
        int ie = max(ip-1,ii);
        int np = ie-ib+1;
        float[] ws = new float[np];
        ws[round(np/2f)-1] =1.0f; 
        RecursiveExponentialFilter ref = new RecursiveExponentialFilter(np/2f);
        ref.apply(ws,ws);
        float[] fv = copy(np,ib,fvs); 
        MedianFinder mf = new MedianFinder(np);
        float fvm = mf.findMedian(ws,fv);
        sps[0][il][its] = fvm;
        sps[1][il][its] = (float)sts.getValue(its);
        sps[2][il][its] = (float)mds.x2;
        sps[3][il][its] = (float)mds.x3;
      }
      for (int itw=0; itw<ntw; ++itw) {
        int itm = itw-1;
        int itp = itw+1;
        itm = max(itm,0);
        itp = min(itp,ntw-1);
        float zm = zw[itm];
        float zp = zw[itp];
        float zi = zw[itw];
        int im = szw.indexOfNearest(zm);
        int ii = szw.indexOfNearest(zi);
        int ip = szw.indexOfNearest(zp);
        int ib = min(im+1,ii);
        int ie = max(ip-1,ii);
        int np = ie-ib+1;
        float[] ws = new float[np];
        ws[round(np/2f)-1] =1.0f; 
        RecursiveExponentialFilter ref = new RecursiveExponentialFilter(np/2f);
        ref.apply(ws,ws);
        float[] fv = copy(np,ib,fvw); 
        MedianFinder mf = new MedianFinder(np);
        float fvm = mf.findMedian(ws,fv);
        spw[0][il][itw] = fvm;
        spw[1][il][itw] = (float)stw.getValue(itw);
        spw[2][il][itw] = (float)mdw.x2;
        spw[3][il][itw] = (float)mdw.x3;
      }

    }
    return new float[][][][]{spw,sps};
  }

  public float[][] convertPoints(float[][][] sps) {
    int nl = sps[0].length;
    FloatList fv = new FloatList();
    FloatList x1 = new FloatList();
    FloatList x2 = new FloatList();
    FloatList x3 = new FloatList();
    for (int il=0; il<nl; ++il) {
      int n1 = sps[0][il].length;
      for (int i1=0; i1<n1; ++i1) {
        float fvi = sps[0][il][i1];
        float x1i = sps[1][il][i1];
        float x2i = sps[2][il][i1];
        float x3i = sps[3][il][i1];
        fv.add(fvi);
        x1.add(x1i);
        x2.add(x2i);
        x3.add(x3i);
      }
    }
    float[] fva = fv.trim();
    float[] x1a = x1.trim();
    float[] x2a = x2.trim();
    float[] x3a = x3.trim();
    return new float[][]{fva,x1a,x2a,x3a};
  }

  private void updateModel(
    SynSeis.Model model, Sampling swx, float[] ws) 
  {
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    float[] t1 = model.t;
    int nt = t1.length;
    float[] t2 = new float[nt];
    int nx = swx.getCount();
    double dx = swx.getDelta();
    double fx = swx.getFirst();
    float[] fa = new float[nx];
    for (int ix=0; ix<nx; ++ix) 
      fa[ix] = (float)(swx.getValue(ix)+ws[ix]);
    si.interpolate(nx,dx,fx,fa,nt,t1,t2);
    model.updateTime(t2);
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
      int nt = t1.length;
      float[] t2 = new float[nt];
      int nx = swx[il].getCount();
      double dx = swx[il].getDelta();
      double fx = swx[il].getFirst();
      float[] fa = new float[nx];
      for (int ix=0; ix<nx; ++ix) 
        fa[ix] = (float)wts[il][ix];
      si.interpolate(nx,dx,fx,fa,nt,t1,t2);
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
    return wsu;
  }

  public float[] synsToSeisShift(Sampling[] sw, Sampling sg, 
    float smin, float smax, float rmin, float rmax, 
    float dmin, float[][] wu, float[][] gu, float[] ndfwm) 
  {
    int nl = sw.length;
    double[] fs = new double[nl];
    double[] ls = new double[nl];
    float[][] wuc = copy(wu);
    float[][] guc = copy(gu);
    for (int il=0; il<nl; ++il) {
      ls[il] = sw[il].getLast();
      fs[il] = sw[il].getFirst();
      float sig = (wu[il].length/20f);
      wuc[il] = gain(sig,wu[il]);
      guc[il] = gain(sig,gu[il]);
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
    float[][] es = dwx.computeErrorsX(sw,wuc,sg,guc);
    return dwx.findShifts(es);
  }

  public float[] synsToSeisShift(
    Sampling sw, Sampling sg, 
    float smin, float smax, float rmin, float rmax, 
    float dmin, float[] wx, float[] gx) 
  {
    DynamicWarpingR dwr = new DynamicWarpingR(smin,smax,sw);
    dwr.setStrainLimits(rmin,rmax);
    dwr.setSmoothness(dmin);
    float sig = (wx.length/20f);
    wx = gain(sig,wx);
    gx = gain(sig,gx);
    float[][] e = dwr.computeErrors(sw,wx,sg,gx);
    float[] ws = dwr.findShifts(e);
    return ws;
  }

  private float[] gain(float sig, float[] fx) {
    float[] gx = mul(fx,fx);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sig);
    ref.apply(gx,gx);
    return div(fx,sqrt(gx));
  }

  public float[][][] seisFlatten(int smax, float epow, float[][] gx) {
    MultiTraceWarping mtw = new MultiTraceWarping();
    mtw.setMaxShift(smax);
    mtw.setPowError(epow);
    mtw.setNullValue(max(gx)+10f);
    float[][] gs = mtw.findShifts(gx);
    float[][] gu = mtw.applyShiftsX(gx,gs);
    return new float[][][]{gs,gu};
  }

  public float[][][] synsFlatten(boolean simple,  
    int smax, float epow, WellLog[] logs, double[][] ndfx,double[][] ndfu) 
  {
    int nl = logs.length;
    SynSeis.Model[] mds = new SynSeis.Model[nl];
    for (int il=0; il<nl; ++il)
      mds[il] = SynSeis.getModel(logs[il]);
    return synsFlatten(simple,smax,epow,mds,ndfx,ndfu);
  }


  public float[][][] synsFlatten(boolean simple,  
    int smax, float epow, SynSeis.Model[] models, 
    double[][] ndfx,double[][] ndfu) 
  {
    double[] ndf = new double[3];
    float[][] wx = computeSyns(simple,models,ndf);
    Sampling st = new Sampling((int)ndf[0],ndf[1],ndf[2]);
    MultiTraceWarping mtw = new MultiTraceWarping();
    mtw.setMaxShift(smax);
    mtw.setPowError(epow);
    float[][] ws = mtw.findShifts(wx);
    //float[][] ws = mtw.findShifts(x2,x3,wx);
    float[][] wu = mtw.applyShiftsW(wx,ws);
    float[][] vu = getValidSamples(st,wu,ndfu);
    float[][][] xs = getValidSamples(st,wx,ws,ndfx);
    System.out.println("ws="+sum(abs(ws)));
    return new float[][][]{xs[0],xs[1],vu};
  }

  public float[][][] synsFlatten(
    boolean simple, int smax, float epow, 
    Sampling s1, SynSeis.Model[] models,float[][] fx, 
    double[][] ndfx,double[][] ndfu) 
  {
    double[] ndf = new double[3];
    float[][] wx = computeSyns(simple,models,ndf);
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[][] wxt = copy(wx);
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x1 = (float)(ndf[2]+i1*ndf[1]);
        float fi = si.interpolate(s1,fx[i2],x1);
        if(wxt[i2][i1]==_vnull){wxt[i2][i1]=fi;}
    }}
    Sampling st = new Sampling((int)ndf[0],ndf[1],ndf[2]);
    MultiTraceWarping mtw = new MultiTraceWarping();
    mtw.setMaxShift(smax);
    mtw.setPowError(epow);
    float[][] ws = mtw.findShifts(wxt);
    float[][] wu = mtw.applyShiftsW(wx,ws);
    float[][] vu = getValidSamples(st,wu,ndfu);
    float[][][] xs = getValidSamples(st,wx,ws,ndfx);
    return new float[][][]{xs[0],xs[1],vu};
  }

  public float[][] computeSyns(
    boolean simple, WellLog[] logs, double[] ndf) 
  {
    int nl = logs.length;
    SynSeis.Model[] mds = new SynSeis.Model[nl];
    for (int il=0; il<nl; ++il)
      mds[il] = SynSeis.getModel(logs[il]);
    return computeSyns(simple,mds,ndf);
  }

  public float[] computeSyn(
    boolean simple, double ds, SynSeis.Model md, double[] ndf) 
  {
    float fpeak = 35f;
    float q = 100.00f;
    double fs = md.tmin();
    double ls = md.tmax();
    fs = round(fs/ds)*ds;
    ls = round(ls/ds)*ds;
    int ns = (int)((ls-fs)/ds)+1;
    ndf[0] = ns; ndf[1] = ds; ndf[2] = fs;
    Sampling sw = new Sampling(ns,ds,fs);
    float[]fa = fillfloat(0.0f,ns);
    if (simple)
      fa = SynSeis.makeSimpleSeismogram(md,fpeak,sw);
    else
      fa = SynSeis.makeBetterSeismogram(md,q,fpeak,sw);
    return fa;
  }

  public float[][] computeSyns(
    boolean simple, SynSeis.Model[] models, double[] ndf) 
  {
    float fpeak = 35f;
    float q = 100.00f;
    double ds = 0.002f;
    double fs =  FLT_MAX;
    double ls = -FLT_MAX;
    int nl = models.length;
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = models[il];
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
      SynSeis.Model md = models[il];
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

