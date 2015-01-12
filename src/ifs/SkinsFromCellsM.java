/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;


import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import static ifs.FaultGeometry.*;

/**
 * Fault cells are decomposed by strikes to construct fault skins.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.12.29
 */
public class SkinsFromCellsM {

  public SkinsFromCellsM(
    Sampling s1d, Sampling s2d, Sampling s3d,
    Sampling s1u, Sampling s2u, Sampling s3u, FaultCell[] fc) 
  {
    _fc = fc;
    _s1d = s1d;
    _s2d = s2d;
    _s3d = s3d;
    _s1u = s1u;
    _s2u = s2u;
    _s3u = s3u;
    _nc = fc.length;
    _n1 = s1d.getCount();
    _n2 = s2d.getCount();
    _n3 = s3d.getCount();
  }

  public FaultSkin[] applyForSkins(FaultSkin[] fs) {
    HashSet<Integer> hsc = new HashSet<Integer>();
    HashSet<FaultSkin> hss = new HashSet<FaultSkin>();
    for (int ic= 0; ic<_nc; ++ic) hsc.add(ic);
    /*
    int minSize = 200;
    FaultSkin[] fss = sortSkins(fs);
    int ns = fss.length;
    System.out.println("ns="+ns);
    for (int is=0; is<ns; ++is) {
      System.out.println("is="+is);
      FaultSkin fsi = fss[is];
      if(fsi.size()<minSize) break;
      FaultSkin fsk = reskinBigSkin(fsi,hsc);
      if(fsk!=null) {
        hss.add(fsk);
        System.out.println("cellsInBigSkin="+fsk.size());
        removeUsedCells(new FaultSkin[]{fsk},hsc);
      }
    }
    */
    int nct = hsc.size();
    while(nct>600) {
      FaultCell[] fc = findStrike(hsc);
      System.out.println("cells="+fc.length);
      if(fc.length<500){break;}
      FaultSkin[] sks = reskin(fc);
      int nk = sks.length;
      if(nk<1) {break;}
      for (int ik=0; ik<nk; ++ik)
        hss.add(sks[ik]);
      removeUsedCells(sks,hsc);
      if(nct-hsc.size()==0){break;}
      nct = hsc.size();
    }
    System.out.println("Skin with remaining cells...");
    int nc = hsc.size();
    System.out.println("nc="+nc);
    if(nc<500){return getSkins(hss);}
    FaultCell[] fcr = new FaultCell[nc];
    int ik=-1;
    for (int ic:hsc) {
      ik++;
      fcr[ik] = _fc[ic];
    }
    FaultSkin[] sks = reskin(fcr);
    int nk = sks.length;
    if(nk<1) {return getSkins(hss);}
    for (ik=0; ik<nk; ++ik) 
      hss.add(sks[ik]);
    return getSkins(hss);
  }

  private FaultSkin reskinBigSkin(FaultSkin fs, HashSet<Integer> hsc) {
    int nc = fs.size();
    FaultCell[] fc = fs.getCells();
    int p2 = fc[0].i2, m2=p2;
    int p3 = fc[0].i3, m3=p3;
    float fpt = fc[0].fp, fps = fpt;
    for (int ic=1; ic<nc; ++ic) {
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      float fpi = fc[ic].fp;
      fps += fpi;
      if(abs(fpi-fpt)>180) {
        System.out.println("test test dfp="+(fpi-fpt));
      }
      fpt = fpi;
      if(i2>p2) p2=i2; if(i2<m2) m2=i2;
      if(i3>p3) p3=i3; if(i3<m3) m3=i3;
    }
    float dp = 20.0f;
    float fp = fps/nc;
    m2 -= 20; m3 -= 20; p2 += 20; p3 += 20;
    FaultCell[] fcs = findInStrikeWindow(m2,p2,m3,p3,dp,fp,hsc); 
    if(fcs.length<500){return null;}
    FaultSkin[] sks = reskin(fc);
    int nk = sks.length;
    if(nk<1) {return null;}
    int maxSize = 0;
    FaultSkin fsm = null;
    for (int ik=0; ik<nk; ++ik) {
      int sz = sks[ik].size();
      if(sz>maxSize) {
        maxSize = sz;
        fsm = sks[ik];
      }
    }
    return fsm;
  }

  private FaultSkin[] sortSkins(FaultSkin[] fsr) {
    int ns = fsr.length;
    FaultSkin[] fss = new FaultSkin[ns];
    int[] ids = new int[ns];
    int[] szs = new int[ns];
    for (int is=0; is<ns; ++is) {
      ids[is] = is;
      szs[is] = fsr[is].size();
    }
    quickIndexSort(szs,ids);
    for (int is=0; is<ns; ++is) {
      int id = ids[ns-is-1];
      fss[is] = fsr[id];
    }
    return fss;
  }

  private FaultSkin[] reskin(FaultCell[] fc) {
    int d = 20;
    float[][][][] flpt = faultImagesFromCells(d,fc);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    fs.setMinSkinSize(round(d*d*4));
    FaultCell[] cells = fs.findCells(flpt);
    return fs.findSkins(cells);
  }


  private void removeUsedCells(FaultSkin[] sk, HashSet<Integer> hsc) {
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    int[] ds = new int[]{3,3,3};
    int[] ns = new int[]{_n1,_n2,_n3};
    FaultCell[] fck = FaultSkin.getCells(sk);
    float[][] xc = setKdTreeCoords(fck);
    KdTree kt = new KdTree(xc);
    HashSet<Integer> hst = new HashSet<Integer>();
    for (int ic:hsc) hst.add(ic);
    for (int ic:hst) {
      int i1 = _fc[ic].i1;
      int i2 = _fc[ic].i2;
      int i3 = _fc[ic].i3;
      int[] is = new int[]{i1,i2,i3};
      getRange(ds,is,ns,xmin,xmax);
      int[] id = kt.findInRange(xmin,xmax);
      int nd = id.length;
      if(nd>ds[0]*ds[0]*2){hsc.remove(ic);}
    }
  }

  private FaultCell[] findStrike(HashSet<Integer> hsc) {
    int dp = 8;
    int maxN = 0;
    int nc = hsc.size();
    int[] dc = new int[nc];
    float[] pm1 = new float[1];
    float[] pp1 = new float[1];
    float[] pm2 = new float[1];
    float[] pp2 = new float[1];
    float[] pm3 = new float[1];
    float[] pp3 = new float[1];
    float[][] xc = setKdTreeStrike(hsc,dc);
    KdTree kt = new KdTree(xc);
    HashSet<Integer> hst = new HashSet<Integer>();
    for (int pt=0; pt<=360; ++pt) {
      float pmi = pt-dp;
      float ppi = pt+dp;
      int nd1=0, nd2=0, nd3=0;
      pm3[0] = pmi; pp3[0] = ppi;
      int[] id1=null, id2=null, id3=null;
      if (pmi<0.0f) {
        pm3[0] = 0.0f;
        pp1[0] = 360f;
        pm1[0] = 360f+pmi;
        id1 = kt.findInRange(pm1,pp1);
      }
      if(ppi>360f) {
        pp3[0] = 360f;
        pm2[0] = 0.0f;
        pp2[0] = ppi-360f;
        id2 = kt.findInRange(pm2,pp2);
      }
      id3 = kt.findInRange(pm3,pp3);
      if(id1!=null){nd1=id1.length;}
      if(id2!=null){nd2=id2.length;}
      if(id3!=null){nd3=id3.length;}
      int nd = nd1+nd2+nd3;
      if(nd>maxN) {
        maxN = nd;
        hst.clear();
        if(id1!=null) {
          for (int ik=0; ik<nd1; ++ik){
            int ip = id1[ik];
            hst.add(dc[ip]);
          }
        }
        if(id2!=null) {
          for (int ik=0; ik<nd2; ++ik){
            int ip = id2[ik];
            hst.add(dc[ip]);
          }
        }
        if(id3!=null) {
          for (int ik=0; ik<nd3; ++ik){
            int ip = id3[ik];
            hst.add(dc[ip]);
          }
        }
      }
    }
    int ik=-1;
    FaultCell[] fc = new FaultCell[maxN];
    for (int ic:hst) {
      ik++;
      fc[ik] = _fc[ic];
    }
    return fc;
  }

  private FaultCell[] findInStrikeWindow(
    int m2, int p2, int m3, int p3,
    float dp, float fp, HashSet<Integer> hsc) 
  {
    int nc = hsc.size();
    int[] dc = new int[nc];
    float[] pm1 = new float[1];
    float[] pp1 = new float[1];
    float[] pm2 = new float[1];
    float[] pp2 = new float[1];
    float[] pm3 = new float[1];
    float[] pp3 = new float[1];
    float[][] xc = setKdTreeStrike(hsc,dc);
    KdTree kt = new KdTree(xc);
    HashSet<Integer> hst = new HashSet<Integer>();
    float pmi = fp-dp;
    float ppi = fp+dp;
    int nd1=0, nd2=0, nd3=0;
    pm3[0] = pmi; pp3[0] = ppi;
    int[] id1=null, id2=null, id3=null;
    if (pmi<0.0f) {
      pm3[0] = 0.0f;
      pp1[0] = 360f;
      pm1[0] = 360f+pmi;
      id1 = kt.findInRange(pm1,pp1);
    }
    if(ppi>360f) {
      pp3[0] = 360f;
      pm2[0] = 0.0f;
      pp2[0] = ppi-360f;
      id2 = kt.findInRange(pm2,pp2);
    }
    id3 = kt.findInRange(pm3,pp3);
    if(id1!=null){nd1=id1.length;}
    if(id2!=null){nd2=id2.length;}
    if(id3!=null){nd3=id3.length;}
    hst.clear();
    int nd = 0;
    if(id1!=null) {
      for (int ik=0; ik<nd1; ++ik){
        int ip = id1[ik];
        int ic = dc[ip];
        nd++; hst.add(ic);
      }
    } 
    if(id2!=null) {
      for (int ik=0; ik<nd2; ++ik){
        int ip = id2[ik];
        int ic = dc[ip];
        nd++; hst.add(ic);
      }
    }
    if(id3!=null) {
      for (int ik=0; ik<nd3; ++ik){
        int ip = id3[ik];
        int ic = dc[ip];
        nd++; hst.add(ic);
      }
    }
    int ik=-1;
    FaultCell[] fc = new FaultCell[nd];
    for (int ic:hst) {
      ik++;
      fc[ik] = _fc[ic];
    }
    return fc;
  }

  private float[][][][] faultImagesFromCells( 
    final int dt, final FaultCell[] fc) {
    int nc = fc.length;
    float sigmaNor = 2.0f;
    final int[][] bb2 = new int[_n1][2];
    final int[][] bb3 = new int[_n1][2];
    final float[][] xc = new float[3][nc];
    final float[][] ws = new float[nc][6];
    final float[][] us = new float[nc][6];
    final float[][] vs = new float[nc][6];
    setKdTreeNodes(fc,xc,ws,us,vs,bb2,bb3);
    final int[] bs1 = setBounds(_n1,xc[0]);
    final int[] bs2 = setBounds(_n2,xc[1]);
    final int[] bs3 = setBounds(_n3,xc[2]);
    final KdTree kt = new KdTree(xc);
    final int[] ns = new int[]{_n1,_n2,_n3};
    final float st = (float)sin(Math.PI/8.0);
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    final float[][][] fl = new float[_n3][_n2][_n1];
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      System.out.println("i3="+i3);
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          if((i2<bb2[i1][0]||i2>bb2[i1][1])){continue;}
          if((i3<bb3[i1][0]||i3>bb3[i1][1])){continue;}
          int[] id = null;
          int di = 10,nd = 0;
          int[] is = new int[]{i1,i2,i3};
          while(nd<20 && di<dt) {
            int[] ds=new int[]{di,di,di};
            getRange(ds,is,ns,xmin,xmax);
            id = kt.findInRange(xmin,xmax);
            nd = id.length;
            di += 2;
          }
          if(nd<10){continue;}
          float wps = 0.0f;
          float sv = 0.25f/(di*di); 
          float su = 0.25f/(di*di); 
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float w1i = fc[ip].w1;
            float w2i = fc[ip].w2;
            float w3i = fc[ip].w3;
            float x1i = fc[ip].i1;
            float x2i = fc[ip].i2;
            float x3i = fc[ip].i3;

            float dx1 = x1i-i1;
            float dx2 = x2i-i2;
            float dx3 = x3i-i3;
            float d11 = dx1*dx1;
            float d22 = dx2*dx2;
            float d33 = dx3*dx3;
            float dsi = d11+d22+d33;
            float wdi = w1i*dx1+w2i*dx2+w3i*dx3;
            if(dsi!=0.0f&&abs(wdi/sqrt(dsi))>st){continue;}
            float d12 = dx1*dx2;
            float d13 = dx1*dx3;
            float d23 = dx2*dx3;

            float w11 = ws[ip][0];
            float w22 = ws[ip][1];
            float w33 = ws[ip][2];
            float w12 = ws[ip][3];
            float w13 = ws[ip][4];
            float w23 = ws[ip][5];

            float u11 = us[ip][0];
            float u22 = us[ip][1];
            float u33 = us[ip][2];
            float u12 = us[ip][3];
            float u13 = us[ip][4];
            float u23 = us[ip][5];

            float v11 = vs[ip][0];
            float v22 = vs[ip][1];
            float v33 = vs[ip][2];
            float v12 = vs[ip][3];
            float v13 = vs[ip][4];
            float v23 = vs[ip][5];

            float wd1 = w12*d12*2.0f;
            float wd2 = w13*d13*2.0f;
            float wd3 = w23*d23*2.0f;

            float ud1 = u12*d12*2.0f;
            float ud2 = u13*d13*2.0f;
            float ud3 = u23*d23*2.0f;

            float vd1 = v12*d12*2.0f;
            float vd2 = v13*d13*2.0f;
            float vd3 = v23*d23*2.0f;

            float wds = w11*d11+w22*d22+w33*d33;
            float uds = u11*d11+u22*d22+u33*d33;
            float vds = v11*d11+v22*d22+v33*d33;

            float gss = 0.0f;
            float wpi = fc[ip].fl;
            wpi = pow(wpi,10.f);
            gss += (wd1+wd2+wd3+wds)*sw;
            gss += (ud1+ud2+ud3+uds)*su;
            gss += (vd1+vd2+vd3+vds)*sv;
            float fli = exp(-gss)*wpi;
            fl[i3][i2][i1] += fli;
            wps += wpi;
          }
          if(wps!=0f) fl[i3][i2][i1] /= wps;
        }
      }
    }});
    return interpFaultImages(bs2,bs3,fl);
  }

  private float[][][][] interpFaultImages(
    int[] bs2, int[] bs3, float[][][] fld) 
  {
    int n1u = _s1u.getCount();
    int n2u = _s2u.getCount();
    int n3u = _s3u.getCount();
    double d1 = _s1d.getDelta();
    double d2 = _s2d.getDelta();
    double d3 = _s3d.getDelta();
    int i2b = bs2[0], i2e = bs2[1];
    int i3b = bs3[0], i3e = bs3[1];
    int n1s = fld[0][0].length;
    int n2s=i2e-i2b,  n3s=i3e-i3b;
    float[][][] flds = new float[n3s][n2s][n1s];
    for (int i3s=0,i3=i3b; i3s<n3s; ++i3,++i3s) 
      for (int i2s=0,i2=i2b; i2s<n2s; ++i2,++i2s) 
        flds[i3s][i2s] = fld[i3][i2];
    int n2us = round(n2s*(float)d2); 
    int n3us = round(n3s*(float)d3);
    float[][][] flus = new float[n3us][n2us][n1u];
    float[][][] fpus = new float[n3us][n2us][n1u];
    float[][][] ftus = new float[n3us][n2us][n1u];
    Sampling s1s = new Sampling(n1s,d1,0.0);
    Sampling s2s = new Sampling(n2s,d2,0.0);
    Sampling s3s = new Sampling(n3s,d3,0.0);
    upgridLikelihood(s1s,s2s,s3s,flds,flus);
    strikeAndDipFromLikelihood(flus,fpus,ftus);
    float[][][] flu = new float[n3u][n2u][n1u];
    float[][][] fpu = new float[n3u][n2u][n1u];
    float[][][] ftu = new float[n3u][n2u][n1u];
    for (int i3s=0,i3=i3b*(int)d3; i3s<n3us; ++i3,++i3s) 
      for (int i2s=0,i2=i2b*(int)d2; i2s<n2us; ++i2,++i2s) { 
        if(i3<n3u&&i2<n2u) {
          flu[i3][i2] = flus[i3s][i2s];
          fpu[i3][i2] = fpus[i3s][i2s];
          ftu[i3][i2] = ftus[i3s][i2s];
        }
      }
    return new float[][][][] {flu,fpu,ftu};
  }

  private void upgridLikelihood(
    final Sampling s1s, final Sampling s2s, final Sampling s3s,
    final float[][][] fld, final float[][][] flu) 
  {
    final int n3 = flu.length;
    final int n2 = flu[0].length;
    final int n1 = flu[0][0].length;
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1 = (float) i1;
          float x2 = (float) i2;
          float x3 = (float) i3;
          flu[i3][i2][i1] = si.interpolate(
          s1s,s2s,s3s,fld,x1,x2,x3);
        }
      }
    }});
  }

  private void strikeAndDipFromLikelihood(final float[][][] fl, 
    final float[][][] fp, final float[][][] ft) {
    div(fl,max(fl),fl);
    final int n3 = fl.length;
    final int n2 = fl[0].length;
    final int n1 = fl[0][0].length;
    final float[][][] u1 = new float[n3][n2][n1];
    final float[][][] u2 = new float[n3][n2][n1];
    final float[][][] u3 = new float[n3][n2][n1];
    final float[][][] sc = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8,4);
    lof.applyForNormal(fl,u1,u2,u3);
    final EigenTensors3 d = new EigenTensors3(n1,n2,n3,false);
    final float db = 15f;
    final float b1 = n1-db;
    final float b2 = n2-db;
    final float b3 = n3-db;
    final float dd = 1.0f/20f;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fli = fl[i3][i2][i1];
          if(fli<0.2f){continue;}
          fli = 1.2f-fli;
          float sci = dd*fli;
          sc[i3][i2][i1]=fli;
          float u1i = -u1[i3][i2][i1];
          float u2i = -u2[i3][i2][i1];
          float u3i = -u3[i3][i2][i1];
          float usi = 1.0f/sqrt(u1i*u1i+u2i*u2i);
          float w1i = -u2i*usi;
          float w2i =  u1i*usi;
          d.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          d.setEigenvectorW(i1,i2,i3,w1i,w2i,0.f);
          if(u2i!=0.0f && u3i!=0.0f) {
            ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
            fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
          }
          if(i1<db){sc[i3][i2][i1]=i1*sci;}
          if(i2<db){sc[i3][i2][i1]=i2*sci;}
          if(i3<db){sc[i3][i2][i1]=i3*sci;}
          if(i1>b1){sc[i3][i2][i1]=(db-i1+b1)*sci;}
          if(i2>b2){sc[i3][i2][i1]=(db-i2+b2)*sci;}
          if(i3>b3){sc[i3][i2][i1]=(db-i3+b3)*sci;}
        }
      }
    }});
    float sigma = 20f;
    float c = sigma*sigma*0.5f;
    d.setEigenvalues(0.001f,1.0f,1.0f);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.1,20);
    lsf.apply(d,c,sc,fl,fl);
    div(fl,max(fl),fl);
  }

  private float[][] setKdTreeCoords(FaultCell[] fc) {
    int nc = fc.length;
    float d1 = (float)_s1d.getDelta();
    float d2 = (float)_s2d.getDelta();
    float d3 = (float)_s3d.getDelta();
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ic++) {
      xc[0][ic] = fc[ic].x1/d1;
      xc[1][ic] = fc[ic].x2/d2;
      xc[2][ic] = fc[ic].x3/d3;
    }
    return xc;
  }

  private void setKdTreeNodes(FaultCell[] fc, float[][] xc, 
    float[][] ws, float[][] us, float[][] vs, int[][] bb2, int[][] bb3) {
    int nc = fc.length;
    for (int i1=0; i1<_n1; ++i1) {
      bb2[i1][0] = _n2; bb2[i1][1] = -_n2;
      bb3[i1][0] = _n3; bb3[i1][1] = -_n3;
    }
    for (int ic=0; ic<nc; ic++) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;

      float w1 = fc[ic].w1;
      float w2 = fc[ic].w2;
      float w3 = fc[ic].w3;

      float u1 = fc[ic].u1;
      float u2 = fc[ic].u2;
      float u3 = fc[ic].u3;

      float v1 = fc[ic].v1;
      float v2 = fc[ic].v2;
      float v3 = fc[ic].v3;

      xc[0][ic] = i1;
      xc[1][ic] = i2;
      xc[2][ic] = i3;

      ws[ic][0] = w1*w1;
      ws[ic][1] = w2*w2;
      ws[ic][2] = w3*w3;
      ws[ic][3] = w1*w2;
      ws[ic][4] = w1*w3;
      ws[ic][5] = w2*w3;

      us[ic][0] = u1*u1;
      us[ic][1] = u2*u2;
      us[ic][2] = u3*u3;
      us[ic][3] = u1*u2;
      us[ic][4] = u1*u3;
      us[ic][5] = u2*u3;

      vs[ic][0] = v1*v1;
      vs[ic][1] = v2*v2;
      vs[ic][2] = v3*v3;
      vs[ic][3] = v1*v2;
      vs[ic][4] = v1*v3;
      vs[ic][5] = v2*v3;

      int b2l = bb2[i1][0];
      int b2r = bb2[i1][1];
      if(i2<b2l) bb2[i1][0] = i2;
      if(i2>b2r) bb2[i1][1] = i2;
      int b3l = bb3[i1][0];
      int b3r = bb3[i1][1];
      if(i3<b3l) bb3[i1][0] = i3;
      if(i3>b3r) bb3[i1][1] = i3;
    }

    for (int i1=0; i1<_n1; ++i1) {
      bb2[i1][0] -= 3; bb3[i1][0] -= 3;
      bb2[i1][1] += 3; bb3[i1][1] += 3;
    }

  }

  private float[][] setKdTreeStrike(HashSet<Integer> fcs, int[] dc) {
    int ik = -1;
    int nc = fcs.size();
    float[][] xc = new float[1][nc];
    for (int ic:fcs) {
      ik++;
      dc[ik] = ic;
      xc[0][ik] = _fc[ic].fp;
    }
    return xc;
  }

  private FaultSkin[] getSkins(HashSet<FaultSkin> hsk) {
    int ik = 0;
    int nk = hsk.size();
    int[] np = new int[nk];
    int[] ip = new int[nk];
    FaultSkin[] sk  = new FaultSkin[nk];
    FaultSkin[] sks = new FaultSkin[nk];
    for (FaultSkin ski:hsk) {
      ip[ik] = ik;
      sk[ik] = ski;
      np[ik] = ski.size();
      ik++;
    }
    quickIndexSort(np,ip);
    for (ik=0; ik<nk; ++ik) {
      int id = ip[nk-ik-1];
      sks[ik] = sk[id];
    }
    return sks;
  }


  private int[] setBounds(int n, float[] x) {
    int[] bs = new int[2];
    int n1m = (int)min(x)-5; 
    int n1p = (int)max(x)+5; 
    if(n1m<0){n1m=0;}
    if(n1p>n){n1p=n;}
    bs[0] = n1m;
    bs[1] = n1p;
    return bs;
  }


  private static void getRange(int[] ds, int[] is, int[] ns, 
    float[] xmin, float[] xmax) {
    int i1m = is[0]-ds[0]; if(i1m<0){i1m=0;}
    int i2m = is[1]-ds[1]; if(i2m<0){i2m=0;}
    int i3m = is[2]-ds[2]; if(i3m<0){i3m=0;}
    int i1p = is[0]+ds[0]; if(i1p>=ns[0]){i1p=ns[0]-1;}
    int i2p = is[1]+ds[1]; if(i2p>=ns[1]){i2p=ns[1]-1;}
    int i3p = is[2]+ds[2]; if(i3p>=ns[2]){i3p=ns[2]-1;}
    xmin[0] = i1m; xmin[1] = i2m; xmin[2] = i3m;
    xmax[0] = i1p; xmax[1] = i2p; xmax[2] = i3p;
  }


  @SuppressWarnings("unchecked")  
  <T> HashSet<T>[] newListArray(int size) {  
    return (HashSet<T>[])new HashSet<?>[size];  
  }  

  private int _n1;
  private int _n2;
  private int _n3;
  private int _nc;
  private Sampling _s1d;
  private Sampling _s2d;
  private Sampling _s3d;
  private Sampling _s1u;
  private Sampling _s2u;
  private Sampling _s3u;
  private FaultCell[] _fc;
}
