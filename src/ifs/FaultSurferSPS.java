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
 * Reconstruct fault surfaces from oriented fault cells. 
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.12.17
 */
public class FaultSurferSPS {

  public FaultSurferSPS(int n1, int n2, int n3, FaultCell[] fc) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fc = fc;
    _nc = fc.length;
  }

  public void setStrikePartitions(int ds) {
    _ds = ds;
  }


  public float getStrike(FaultSkin sk) {
    float fpa = 0.0f;
    FaultCell[] fc = FaultSkin.getCells(new FaultSkin[] {sk});
    int nc = fc.length;
    for (int ic=0; ic<nc; ++ic)
      fpa += fc[ic].fp;
    return fpa/nc;
  }

  public FaultSkin[] applySurferM(int minSkinSize) {
    HashSet<Integer> hsc = new HashSet<Integer>();
    HashSet<FaultSkin> hss = new HashSet<FaultSkin>();
    for (int ic=0; ic<_nc; ++ic) hsc.add(ic);
    int nct = hsc.size();
    //int minCellSize = round(max(_n2,_n3)*_n1*0.5f);
    //int minSkinSize = round(max(_n2,_n3)*_n1*0.4f);
    while(nct>minSkinSize) {
      FaultCell[] fc = findStrike(hsc);
      System.out.println("cells="+fc.length);
      if(fc.length<minSkinSize){break;}
      FaultSkin[] sks = reskin(minSkinSize,fc);
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
    if(nc<minSkinSize){return getSkins(hss);}
    FaultCell[] fcr = new FaultCell[nc];
    int ik=-1;
    for (int ic:hsc) {
      ik++;
      fcr[ik] = _fc[ic];
    }
    FaultSkin[] sks = reskin(minSkinSize,fcr);
    int nk = sks.length;
    if(nk<1) {return getSkins(hss);}
    for (ik=0; ik<nk; ++ik) 
      hss.add(sks[ik]);
    return getSkins(hss);
  }

  public float[][][] test() {
    HashSet<Integer> hsc = new HashSet<Integer>();
    for (int ic=0; ic<_nc; ++ic) hsc.add(ic);
    FaultCell[] fc = findStrike(hsc);
    System.out.println("cells="+fc.length);
    float[][][][] flpt = faultImagesFromCells(fc);
    return flpt[0];
  }


  private FaultSkin[] reskin(int minSkinSize,FaultCell[] fc) {
    float[][][][] flpt = faultImagesFromCells(fc);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    //fs.setMinSkinSize(round(d*d*4));
    fs.setMinSkinSize(minSkinSize);
    //fs.setMaxDeltaLikelihood(0.2);
    FaultCell[] cells = fs.findCells(flpt);
    return fs.findSkins(cells);
  }

  private float[][][][] faultImagesFromCells(FaultCell[] fc) {
    ScreenPoissonSurfer sps = new ScreenPoissonSurfer();
    float[][][][] us = sps.setNormalsFromCells(_n1,_n2,_n3,fc);
    float[][][] fb = sps.faultIndicator(us);
    float[][][] fl = convertToFaultLikelihood(fb);
    float[][][] fp = new float[_n3][_n2][_n1];
    float[][][] ft = new float[_n3][_n2][_n1];
    computeStrikeAndDip(fl,fp,ft);
    return new float[][][][]{fl,fp,ft};
  }

  private float[][][] convertToFaultLikelihood(float[][][] fb) {
    int n3 = fb.length;
    int n2 = fb[0].length;
    int n1 = fb[0][0].length;
    float hpi = (float)Math.PI*0.5f;
    float[][][] fl = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) {
        float fbi = fb[i3][i2][i1];
        if(abs(fbi)>hpi){continue;}
        fl[i3][i2][i1] = cos(fbi);
      }
    return fl;
  }

  private void computeStrikeAndDip(final float[][][] fl, 
    final float[][][] fp, final float[][][] ft) {
    final int n3 = fl.length;
    final int n2 = fl[0].length;
    final int n1 = fl[0][0].length;
    final float[][][] u1 = new float[n3][n2][n1];
    final float[][][] u2 = new float[n3][n2][n1];
    final float[][][] u3 = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8,4);
    lof.applyForNormal(fl,u1,u2,u3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = -u1[i3][i2][i1];
          float u2i = -u2[i3][i2][i1];
          float u3i = -u3[i3][i2][i1];
          if(u2i!=0.0f && u3i!=0.0f) {
            ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
            fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
          }
        }
      }
    }});
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



  private void smooth(final float[][][] sf, 
    final float[][][] fp, final float[][][] ft) {
    div(sf,max(sf),sf);
    final int n3 = sf.length;
    final int n2 = sf[0].length;
    final int n1 = sf[0][0].length;
    final float[][][] u1 = new float[n3][n2][n1];
    final float[][][] u2 = new float[n3][n2][n1];
    final float[][][] u3 = new float[n3][n2][n1];
    final float[][][] sc = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8,4);
    lof.applyForNormal(sf,u1,u2,u3);
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
          float sfi = sf[i3][i2][i1];
          if(sfi<0.2f){continue;}
          sfi = 1.2f-sfi;
          float sci = dd*sfi;
          sc[i3][i2][i1]=sfi;
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
    lsf.apply(d,c,sc,sf,sf);
    div(sf,max(sf),sf);
  }

  private float[][] setKdTreeCoords(FaultCell[] fc) {
    int nc = fc.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ic++) {
      xc[0][ic] = fc[ic].x1;
      xc[1][ic] = fc[ic].x2;
      xc[2][ic] = fc[ic].x3;
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
      bb2[i1][0] -= 3;
      bb3[i1][0] -= 3;
      bb2[i1][1] += 3;
      bb3[i1][1] += 3;
    }

  }


  private void setKdTreePoints(FaultCell[] fc, 
    float[][] xc, float[][] ws, float[][] us, float[][] vs, int[][][] mk) {
    int nc = fc.length;
    float d = _dx;
    for (int ic=0; ic<nc; ic++) {
      float x1 = fc[ic].x1;
      float x2 = fc[ic].x2;
      float x3 = fc[ic].x3;

      float w1 = fc[ic].w1;
      float w2 = fc[ic].w2;
      float w3 = fc[ic].w3;

      float u1 = fc[ic].u1;
      float u2 = fc[ic].u2;
      float u3 = fc[ic].u3;

      float v1 = fc[ic].v1;
      float v2 = fc[ic].v2;
      float v3 = fc[ic].v3;

      xc[0][ic] = x1;
      xc[1][ic] = x2;
      xc[2][ic] = x3;

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
      int i1m = round(x1-d);if(i1m<0){i1m=0;}
      int i2m = round(x2-d);if(i2m<0){i2m=0;}
      int i3m = round(x3-d);if(i3m<0){i3m=0;}
      int i1p = round(x1+d);if(i1p>=_n1){i1p=_n1;}
      int i2p = round(x2+d);if(i2p>=_n2){i2p=_n2;}
      int i3p = round(x3+d);if(i3p>=_n3){i3p=_n3;}
      for (int i3=i3m; i3<i3p; ++i3) 
      for (int i2=i2m; i2<i2p; ++i2) 
      for (int i1=i1m; i1<i1p; ++i1) 
        mk[i3][i2][i1]=1;
    }
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

  private static float distance(float[] x, float[] y) {
    float d1 = y[0]-x[0];
    float d2 = y[1]-x[1];
    float d3 = y[2]-x[2];
    return sqrt(d1*d1+d2*d2+d3*d3);
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

  private FaultCell[][] strikePartition() {
    int ns = (int)(360/_ds)+1;
    FaultCell[][] fcs = new FaultCell[ns][];
    HashSet<Integer>[] fls = newListArray(ns);
    for (int is=0; is<ns; ++is) 
      fls[is] = new HashSet<Integer>();

    for (int ic=0; ic<_nc; ++ic) {
      int is = (int)(_fc[ic].fp/_ds);
      fls[is].add(ic);
    }
    for (int is=0; is<ns; ++is) {
      int k = 0;
      int nc = fls[is].size();
      fcs[is] = new FaultCell[nc];
      for (int ic:fls[is]) {
        fcs[is][k] = _fc[ic];
        k++;
      }
    }
    return fcs;
  }

  @SuppressWarnings("unchecked")  
  <T> HashSet<T>[] newListArray(int size) {  
    return (HashSet<T>[])new HashSet<?>[size];  
  }  

  private int _n1;
  private int _n2;
  private int _n3;
  private int _nc;
  private int _ds=15;
  private int _dx=20;
  private FaultCell[] _fc;
}
