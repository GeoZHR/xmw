/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import java.util.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import static ifs.FaultGeometry.*;

/**
 * Reconstruct fault surfaces from oriented fault cells. 
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.12.17
 */
public class FaultSurfer {

  public FaultSurfer(int n1, int n2, int n3, FaultCell[] fc) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fc = fc;
    _nc = fc.length;
  }

  public void setStrikePartitions(int ds) {
    _ds = ds;
  }

  public FaultSkin[] applySurfer() {
    HashSet<FaultSkin> hfs = new HashSet<FaultSkin>();
    FaultCell[][] fcs = strikePartition();
    int ns = fcs.length;
    //build fault skins in each strike window
    System.out.println("Skin with cells in each strike window...");
    for (int is=0; is<ns; ++is) {
      if(fcs[is].length<2000){continue;}
      FaultSkin[] sks = reskin(fcs[is]);
      int nk = sks.length;
      if(nk<1) {continue;}
      for (int ik=0; ik<nk; ++ik) {
        hfs.add(sks[ik]);
      }
    }
    //build fault skins with multiple strikes
    System.out.println("Skin with remaining cells...");
    FaultCell[] fcr = checkUnusedCells(getSkins(hfs));
    if(fcr.length<2000){return getSkins(hfs);}
    FaultSkin[] sks = reskin(fcr);
    int nk = sks.length;
    if(nk<1) {return getSkins(hfs);}
    for (int ik=0; ik<nk; ++ik) {
      hfs.add(sks[ik]);
    }
    return getSkins(hfs);
  }

  private FaultSkin[] getSkins(HashSet<FaultSkin> sk) {
    int nk = sk.size();
    FaultSkin[] sks = new FaultSkin[nk];
    int ik = 0;
    for (FaultSkin ski:sk) {
      sks[ik] = ski;
      ik++;
    }
    return sks;
  }

  private FaultSkin[] reskin(FaultCell[] fc) {
    float[][][] fl = new float[_n3][_n2][_n1];
    float[][][] fp = new float[_n3][_n2][_n1];
    float[][][] ft = new float[_n3][_n2][_n1];
    fl = reconstructFaultImagesFromCells(fc,fp,ft);
    //fl = recomputeFaultImagesFromCells(fc,fp,ft);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.3f,0.6f);
    fs.setMinSkinSize(2000);
    //fs.setMaxDeltaLikelihood(0.2);
    FaultCell[] cells = fs.findCells(new float[][][][] {fl,fp,ft});
    return fs.findSkins(cells);
  }

  private FaultCell[] checkUnusedCells(FaultSkin[] sk) {
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    int[] ds = new int[]{3,3,3};
    int[] ns = new int[]{_n1,_n2,_n3};
    FaultCell[] fck = FaultSkin.getCells(sk);
    float[][] xc = setKdTreeCoords(fck);
    KdTree kt = new KdTree(xc);
    HashSet<Integer> hfc = new HashSet<Integer>();
    for (int ic=0; ic<_nc; ++ic) {
      hfc.add(ic);
      int i1 = _fc[ic].i1;
      int i2 = _fc[ic].i2;
      int i3 = _fc[ic].i3;
      int[] is = new int[]{i1,i2,i3};
      getRange(ds,is,ns,xmin,xmax);
      int[] id = kt.findInRange(xmin,xmax);
      int nd = id.length;
      if(nd>ds[0]*ds[0]*2){hfc.remove(ic);}
    }
    int k = 0;
    int np = hfc.size();
    FaultCell[] fcr = new FaultCell[np];
    for (int ic:hfc) {
      fcr[k] = _fc[ic];
      k++;
    }
    return fcr;
  }

  private float[][][] recomputeFaultImagesFromCells(
    final FaultCell[] fc, final float[][][] fp, final float[][][] ft) 
  {
    float sigmaNor = 4.0f;
    float sigmaPhi = 10.0f;
    float sigmaTheta = 40.0f;
    final int[] ds = new int[]{10,10,10};
    final int[] ns = new int[]{_n1,_n2,_n3};
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    final float sv = 1.0f/(sigmaPhi*sigmaPhi); 
    final float su = 1.0f/(sigmaTheta*sigmaTheta); 
    int nc = fc.length;
    final float[][] xc = new float[3][nc];
    final float[][] ws = new float[nc][6];
    final float[][] us = new float[nc][6];
    final float[][] vs = new float[nc][6];
    setKdTreePoints(fc,xc,ws,us,vs);
    final int[] bs1 = setBounds(_n1,xc[0]);
    final int[] bs2 = setBounds(_n2,xc[1]);
    final int[] bs3 = setBounds(_n3,xc[2]);
    final KdTree kt = new KdTree(xc);
    final float[][][] sf = zerofloat(_n1,_n2,_n3);
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      System.out.println("i3="+i3);
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          float[] y = new float[]{i1,i2,i3};
          int ne = kt.findNearest(y);
          float x1 = fc[ne].x1;
          float x2 = fc[ne].x2;
          float x3 = fc[ne].x3;
          float dd = distance(new float[]{x1,x2,x3},y);
          if(dd>10.0f){continue;}
          int[] is = new int[]{i1,i2,i3};
          getRange(ds,is,ns,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<20){continue;}
          float wps = 0.0f;
          float w1a = 0.0f;
          float w2a = 0.0f;
          float w3a = 0.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float w1i = fc[ip].w1;
            float w2i = fc[ip].w2;
            float w3i = fc[ip].w3;
            /*
            float dx1 = i1-fc[ip].x1;
            float dx2 = i2-fc[ip].x2;
            float dx3 = i3-fc[ip].x3;
            float dxs = dx1*dx1+dx2*dx2+dx3*dx3;
            float wpi = exp(-dxs/100.0f);
            */
            w1a += w1i;
            w2a += w2i;
            w3a += w3i;
          }
          float wsa = 1.0f/sqrt(w1a*w1a+w2a*w2a+w3a*w3a);
          w1a *= wsa;
          w2a *= wsa;
          w3a *= wsa;
          float fti = faultDipFromNormalVector(w1a,w2a,w3a);
          float fpi = faultStrikeFromNormalVector(w1a,w2a,w3a);
          float[] u = faultDipVectorFromStrikeAndDip(fpi,fti);
          float[] v = faultStrikeVectorFromStrikeAndDip(fpi,fti);
          fp[i3][i2][i1] = fpi;
          ft[i3][i2][i1] = fti;

          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float fli = fc[ip].fl;

            float dx1 = i1-fc[ip].x1;
            float dx2 = i2-fc[ip].x2;
            float dx3 = i3-fc[ip].x3;

            float w1s = w1a*w1a;
            float w2s = w2a*w2a;
            float w3s = w3a*w3a;
            float w12 = w1a*w2a;
            float w13 = w1a*w3a;
            float w23 = w2a*w3a;

            float u1s = u[0]*u[0];
            float u2s = u[1]*u[1];
            float u3s = u[2]*u[2];
            float u12 = u[0]*u[1];
            float u13 = u[0]*u[2];
            float u23 = u[1]*u[2];

            float v1s = v[0]*v[0];
            float v2s = v[1]*v[1];
            float v3s = v[2]*v[2];
            float v12 = v[0]*v[1];
            float v13 = v[0]*v[2];
            float v23 = v[1]*v[2];

            float d11 = dx1*dx1;
            float d12 = dx1*dx2;
            float d13 = dx1*dx3;
            float d22 = dx2*dx2;
            float d23 = dx2*dx3;
            float d33 = dx3*dx3;

            float g11 = w1s*d11+w12*d12+w13*d13;
            float g12 = u1s*d11+u12*d12+u13*d13;
            float g13 = v1s*d11+v12*d12+v13*d13;
            float g21 = w12*d12+w2s*d22+w23*d23;
            float g22 = u12*d12+u2s*d22+u23*d23;
            float g23 = v12*d12+v2s*d22+v23*d23;
            float g31 = w13*d13+w23*d23+w3s*d33;
            float g32 = u13*d13+u23*d23+u3s*d33;
            float g33 = v13*d13+v23*d23+v3s*d33;

            float gss = 0.0f;
            gss += (g11+g21+g31)*sw;
            gss += (g12+g22+g32)*su;
            gss += (g13+g23+g33)*sv;
            float wpi = pow(fli,8.0f);
            float sfi = exp(-gss)*wpi;
            sf[i3][i2][i1] += sfi;

            wps += wpi;
          }
          sf[i3][i2][i1] /= wps;
        }
      }
    }});
    return div(sf,max(sf));
  }


  private float[][][] reconstructFaultImagesFromCells(
    final FaultCell[] fc, final float[][][] fp, final float[][][] ft) 
  {
    float sigmaNor = 4.0f;
    float sigmaPhi = 10.0f;
    float sigmaTheta = 40.0f;
    final int[] ds = new int[]{10,10,10};
    final int[] ns = new int[]{_n1,_n2,_n3};
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    final float sv = 1.0f/(sigmaPhi*sigmaPhi); 
    final float su = 1.0f/(sigmaTheta*sigmaTheta); 
    int nc = fc.length;
    final float[][] xc = new float[3][nc];
    final float[][] ws = new float[nc][6];
    final float[][] us = new float[nc][6];
    final float[][] vs = new float[nc][6];
    setKdTreePoints(fc,xc,ws,us,vs);
    final int[] bs1 = setBounds(_n1,xc[0]);
    final int[] bs2 = setBounds(_n2,xc[1]);
    final int[] bs3 = setBounds(_n3,xc[2]);
    final KdTree kt = new KdTree(xc);
    final float[][][] sf = zerofloat(_n1,_n2,_n3);
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      System.out.println("i3="+i3);
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          float[] y = new float[]{i1,i2,i3};
          int ne = kt.findNearest(y);
          float x1 = fc[ne].x1;
          float x2 = fc[ne].x2;
          float x3 = fc[ne].x3;
          float dd = distance(new float[]{x1,x2,x3},y);
          if(dd>10.0f){continue;}
          int[] is = new int[]{i1,i2,i3};
          getRange(ds,is,ns,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<20){continue;}
          float wps = 0.0f;
          float w1a = 0.0f;
          float w2a = 0.0f;
          float w3a = 0.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float fli = fc[ip].fl;
            float w1i = fc[ip].w1;
            float w2i = fc[ip].w2;
            float w3i = fc[ip].w3;

            float dx1 = i1-fc[ip].x1;
            float dx2 = i2-fc[ip].x2;
            float dx3 = i3-fc[ip].x3;

            float w1s = ws[ip][0];
            float w2s = ws[ip][1];
            float w3s = ws[ip][2];
            float w12 = ws[ip][3];
            float w13 = ws[ip][4];
            float w23 = ws[ip][5];

            float u1s = us[ip][0];
            float u2s = us[ip][1];
            float u3s = us[ip][2];
            float u12 = us[ip][3];
            float u13 = us[ip][4];
            float u23 = us[ip][5];

            float v1s = vs[ip][0];
            float v2s = vs[ip][1];
            float v3s = vs[ip][2];
            float v12 = vs[ip][3];
            float v13 = vs[ip][4];
            float v23 = vs[ip][5];

            float d11 = dx1*dx1;
            float d12 = dx1*dx2;
            float d13 = dx1*dx3;
            float d22 = dx2*dx2;
            float d23 = dx2*dx3;
            float d33 = dx3*dx3;

            float g11 = w1s*d11+w12*d12+w13*d13;
            float g12 = u1s*d11+u12*d12+u13*d13;
            float g13 = v1s*d11+v12*d12+v13*d13;
            float g21 = w12*d12+w2s*d22+w23*d23;
            float g22 = u12*d12+u2s*d22+u23*d23;
            float g23 = v12*d12+v2s*d22+v23*d23;
            float g31 = w13*d13+w23*d23+w3s*d33;
            float g32 = u13*d13+u23*d23+u3s*d33;
            float g33 = v13*d13+v23*d23+v3s*d33;

            float gss = 0.0f;
            gss += (g11+g21+g31)*sw;
            gss += (g12+g22+g32)*su;
            gss += (g13+g23+g33)*sv;
            float wpi = pow(fli,8.0f);
            float sfi = exp(-gss)*wpi;
            sf[i3][i2][i1] += sfi;

            w1a += w1i*sfi;
            w2a += w2i*sfi;
            w3a += w3i*sfi;
            wps += wpi;
          }
          float wsa = 1.0f/sqrt(w1a*w1a+w2a*w2a+w3a*w3a);
          w1a *= wsa;
          w2a *= wsa;
          w3a *= wsa;
          ft[i3][i2][i1] = faultDipFromNormalVector(w1a,w2a,w3a);
          fp[i3][i2][i1] = faultStrikeFromNormalVector(w1a,w2a,w3a);
          sf[i3][i2][i1] /= wps;
        }
      }
    }});
    return div(sf,max(sf));
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

  private void setKdTreePoints(
    FaultCell[] fc,float[][] xc, float[][] ws, float[][] us, float[][] vs) {
    int nc = fc.length;
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
    xmin[0] = i1m;
    xmin[1] = i2m;
    xmin[2] = i3m;
    xmax[0] = i1p;
    xmax[1] = i2p;
    xmax[2] = i3p;
  }

  private FaultCell[][] strikePartition() {
    int ns = (int)(360/_ds);
    FaultCell[][] fcs = new FaultCell[ns][];
    HashSet<Integer>[] fls = new HashSet[ns];
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

  private int _n1;
  private int _n2;
  private int _n3;
  private int _nc;
  private int _ds=15;
  private FaultCell[] _fc;
}
