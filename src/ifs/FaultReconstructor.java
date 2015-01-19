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
 * Reconstruct fault images from fault cells.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.15
 */

public class FaultReconstructor {

  public FaultReconstructor(int n1, int n2, int n3, FaultCell[] fcs){
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fcs = fcs;
  }

  public FaultSkin[] recomputeSkins(int minSkinSize) {
    HashSet<Integer> hsc = new HashSet<Integer>();
    for (int ic=0; ic<_fcs.length; ++ic) hsc.add(ic);
    FaultCell[][][] cg = new FaultCell[_n3][_n2][_n1];
    while(hsc.size()>minSkinSize/10) {
      System.out.println("cells remaining:"+hsc.size());
      FaultCell[] fci = findStrike(hsc);
      findCells(20,fci,cg,hsc);
    }
    //return getCells(cg);
    FaultCell[] cells = getCells(cg);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    fs.setMinSkinSize(minSkinSize);
    return fs.findSkins(cells);
  }

  public FaultCell[] applyForFaultImages() {
    int np = 360;
    int owl = 20;
    double owf = 0.5;
    final int dt = 30;
    final float sigmaNor = 2.0f;
    KdTree kt1 = setStrikeKdTree();
    final int[] ns = new int[]{_n1,_n2,_n3};
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    final OverlapWinds ow = new OverlapWinds(np,owl,owf);
    final int l1 = ow.getL1(), m1 = ow.getM1();
    final FaultCell[][][] fcg = new FaultCell[_n3][_n2][_n1];
    final float[][][] fl  = new float[_n3][_n2][_n1];
    final float[][][] g11 = new float[_n3][_n2][_n1];
    final float[][][] g12 = new float[_n3][_n2][_n1];
    final float[][][] g13 = new float[_n3][_n2][_n1];
    final float[][][] g22 = new float[_n3][_n2][_n1];
    final float[][][] g23 = new float[_n3][_n2][_n1];
    final float[][][] g33 = new float[_n3][_n2][_n1];
    float fpt = -1f;
    for (int k1=0; k1<m1; ++k1) {
      System.out.println("k1="+k1);
      final int p1 = ow.getI1(k1);
      if(k1>2){fpt=ow.getI1(k1-3);}
      float[] pmin = new float[]{p1};
      float[] pmax = new float[]{p1+l1-1};
      int[] id = kt1.findInRange(pmin,pmax);
      int nd = id.length;
      if(nd<1){continue;}
      final float[][] xc = new float[3][nd];
      final float[][] ws = new float[nd][6];
      final float[][] us = new float[nd][6];
      final float[][] vs = new float[nd][6];
      final int[][] bb2 = new int[_n1][2];
      final int[][] bb3 = new int[_n1][2];
      final FaultCell[] fc = setForInterp(id,xc,ws,us,vs,bb2,bb3);
      final KdTree kt2 = new KdTree(xc);
      final int[] bs1 = setBounds(_n1,xc[0]);
      final int[] bs2 = setBounds(_n2,xc[1]);
      final int[] bs3 = setBounds(_n3,xc[2]);
      Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
      public void compute(int i3) {
        float[] xmin = new float[3];
        float[] xmax = new float[3];
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
              id = kt2.findInRange(xmin,xmax);
              nd = id.length;
              di += 2;
            }
            if(nd<10){continue;}
            float fss = 0.0f;
            float g11s = 0.0f;
            float g12s = 0.0f;
            float g13s = 0.0f;
            float g22s = 0.0f;
            float g23s = 0.0f;
            float g33s = 0.0f;
            float sv = 0.25f/(di*di); 
            float su = 0.25f/(di*di); 
            float scs = 0.0f;
            for (int ik=0; ik<nd; ++ik) {
              int ip = id[ik];
              float x1i = fc[ip].i1;
              float x2i = fc[ip].i2;
              float x3i = fc[ip].i3;

              float dx1 = x1i-i1;
              float dx2 = x2i-i2;
              float dx3 = x3i-i3;
              float d11 = dx1*dx1;
              float d22 = dx2*dx2;
              float d33 = dx3*dx3;
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
              int fpi = round(fc[ip].fp);
              int fli = round(fc[ip].fl);
              float wpi = 1.0f;//fli*ow.getWeight(p1,fpi-p1);
              gss += (wd1+wd2+wd3+wds)*sw;
              gss += (ud1+ud2+ud3+uds)*su;
              gss += (vd1+vd2+vd3+vds)*sv;
              float sfi = exp(-gss)*wpi;
              scs  += wpi;
              fss  += sfi;
              g11s += sfi*w11;
              g12s += sfi*w12;
              g13s += sfi*w13;
              g22s += sfi*w22;
              g23s += sfi*w23;
              g33s += sfi*w33;
            }
            g11[i3][i2][i1] += g11s;
            g12[i3][i2][i1] += g12s;
            g13[i3][i2][i1] += g13s;
            g22[i3][i2][i1] += g22s;
            g23[i3][i2][i1] += g23s;
            g33[i3][i2][i1] += g33s;
            fl[i3][i2][i1]  += fss/scs;
          }
        }
      }});
      findCells(fpt,fcg,fl,g11,g12,g13,g22,g23,g33);
    }
    findCells(361,fcg,fl,g11,g12,g13,g22,g23,g33);
    return getCells(fcg);
  }

  private void findCells(
    float fpt, FaultCell[][][] fcg, float[][][] fl,
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33) 
  {
    float[][][] flt = copy(fl);
    float[][][] u1 = new float[_n3][_n2][_n1];
    float[][][] u2 = new float[_n3][_n2][_n1];
    float[][][] u3 = new float[_n3][_n2][_n1];
    float[][][] fp = new float[_n3][_n2][_n1];
    float[][][] ft = new float[_n3][_n2][_n1];
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    smooth(flt,fp,ft,u1,u2,u3);
    div(flt,max(flt),flt);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    FaultCell[] fcs = fs.findCells(new float[][][][]{flt,fp,ft});
    for (FaultCell fci:fcs){
      int i1 = fci.i1;
      int i2 = fci.i2;
      int i3 = fci.i3;
      int i1m = i1-1; if(i1m<0){i1m=0;}
      int i2m = i2-1; if(i2m<0){i2m=0;}
      int i3m = i3-1; if(i3m<0){i3m=0;}
      int i1p = i1+1; if(i1p>=_n1){i1p=_n1-1;}
      int i2p = i2+1; if(i2p>=_n2){i2p=_n2-1;}
      int i3p = i3+1; if(i3p>=_n3){i3p=_n3-1;}
      if(fci.fp<fpt) {
        if(fcg[i3][i2][i1 ]==null) {fcg[i3][i2][i1 ]=fci; continue;}
        if(fcg[i3m][i2][i1]==null) {fcg[i3m][i2][i1]=fci; continue;}
        if(fcg[i3p][i2][i1]==null) {fcg[i3p][i2][i1]=fci; continue;}
        if(fcg[i3][i2m][i1]==null) {fcg[i3][i2m][i1]=fci; continue;}
        if(fcg[i3][i2p][i1]==null) {fcg[i3][i2p][i1]=fci; continue;}
        if(fcg[i3][i2][i1m]==null) {fcg[i3][i2][i1m]=fci; continue;}
        if(fcg[i3][i2][i1p]==null) {fcg[i3][i2][i1p]=fci; continue;}
      }
    }
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      if(fp[i3][i2][i1]<fpt-5) {
        fl[i3][i2][i1]  = 0f;
        g11[i3][i2][i1] = 0f;
        g12[i3][i2][i1] = 0f;
        g13[i3][i2][i1] = 0f;
        g22[i3][i2][i1] = 0f;
        g23[i3][i2][i1] = 0f;
        g33[i3][i2][i1] = 0f;
      }
    }}}
  }

  private float computeStrike(
    float g11, float g12, float g13,
    float g22, float g23, float g33) 
  {
    double[][] a = new double[3][3];
    double[][] z = new double[3][3];
    double[] e = new double[3];
    a[0][0] = g11;
    a[0][1] = g12;
    a[0][2] = g13;
    a[1][0] = g12;
    a[1][1] = g22;
    a[1][2] = g23;
    a[2][0] = g13;
    a[2][1] = g23;
    a[2][2] = g33;
    Eigen.solveSymmetric33(a,z,e);
    float u1i = (float)z[0][0];
    float u2i = (float)z[0][1];
    float u3i = (float)z[0][2];
    if (u1i>0.0f) {
      u1i = -u1i;
      u2i = -u2i;
      u3i = -u3i;
    }
    if(u2i!=0.0f && u3i!=0.0f) {
      return faultStrikeFromNormalVector(u1i,u2i,u3i);
    } else {return -1;}
  }

  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n1 = g11[0][0].length;
    final int n2 = g11[0].length;
    final int n3 = g11.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        double[][] a = new double[3][3];
        double[][] z = new double[3][3];
        double[] e = new double[3];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            a[0][0] = g11[i3][i2][i1];
            a[0][1] = g12[i3][i2][i1];
            a[0][2] = g13[i3][i2][i1];
            a[1][0] = g12[i3][i2][i1];
            a[1][1] = g22[i3][i2][i1];
            a[1][2] = g23[i3][i2][i1];
            a[2][0] = g13[i3][i2][i1];
            a[2][1] = g23[i3][i2][i1];
            a[2][2] = g33[i3][i2][i1];
            Eigen.solveSymmetric33(a,z,e);
            float u1i = (float)z[0][0];
            float u2i = (float)z[0][1];
            float u3i = (float)z[0][2];
            if (u1i>0.0f) {
              u1i = -u1i;
              u2i = -u2i;
              u3i = -u3i;
            }
            u1[i3][i2][i1] = u1i;
            u2[i3][i2][i1] = u2i;
            u3[i3][i2][i1] = u3i;
          }
        }
      }
    });
  }


  private void smooth(
    final float[][][] sf, final float[][][] fp, final float[][][] ft,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    div(sf,max(sf),sf);
    final int n3 = sf.length;
    final int n2 = sf[0].length;
    final int n1 = sf[0][0].length;
    /*
    final float[][][] sc = new float[n3][n2][n1];
    final EigenTensors3 d = new EigenTensors3(n1,n2,n3,false);
    final float db = 15f;
    final float b1 = n1-db;
    final float b2 = n2-db;
    final float b3 = n3-db;
    final float dd = 1.0f/20f;
    */
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          /*
          float sfi = sf[i3][i2][i1];
          if(sfi<0.2f){continue;}
          sfi = 1.2f-sfi;
          float sci = dd*sfi;
          sc[i3][i2][i1]=sfi;
          */
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          /*
          float usi = 1.0f/sqrt(u1i*u1i+u2i*u2i);
          float w1i = -u2i*usi;
          float w2i =  u1i*usi;
          d.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          d.setEigenvectorW(i1,i2,i3,w1i,w2i,0.f);
          */
          if(u2i!=0.0f && u3i!=0.0f) {
            ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
            fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
          }
          /*
          if(i1<db){sc[i3][i2][i1]=i1*sci;}
          if(i2<db){sc[i3][i2][i1]=i2*sci;}
          if(i3<db){sc[i3][i2][i1]=i3*sci;}
          if(i1>b1){sc[i3][i2][i1]=(db-i1+b1)*sci;}
          if(i2>b2){sc[i3][i2][i1]=(db-i2+b2)*sci;}
          if(i3>b3){sc[i3][i2][i1]=(db-i3+b3)*sci;}
          */
        }
      }
    }});
    /*
    float sigma = 20f;
    float c = sigma*sigma*0.5f;
    d.setEigenvalues(0.001f,1.0f,1.0f);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.1,20);
    lsf.apply(d,c,sc,sf,sf);
    div(sf,max(sf),sf);
    */
  }


  private FaultCell[] setForInterp(final int[] id, 
    final float[][] xc, final float[][] ws,  
    final float[][] us, final float[][] vs, final int[][] bb2, final int[][] bb3) 
  { 
    final int nd = id.length;
    final FaultCell[] fc = new FaultCell[nd];
    for (int i1=0; i1<_n1; ++i1) {
      bb2[i1][0] = _n2; bb2[i1][1] = -_n2;
      bb3[i1][0] = _n3; bb3[i1][1] = -_n3;
    }

    for (int ip=0; ip<nd; ++ip) {
      int ic = id[ip];
      fc[ip] = _fcs[ic];
      int i1 = _fcs[ic].i1;
      int i2 = _fcs[ic].i2;
      int i3 = _fcs[ic].i3;

      float w1 = _fcs[ic].w1;
      float w2 = _fcs[ic].w2;
      float w3 = _fcs[ic].w3;

      float u1 = _fcs[ic].u1;
      float u2 = _fcs[ic].u2;
      float u3 = _fcs[ic].u3;

      float v1 = _fcs[ic].v1;
      float v2 = _fcs[ic].v2;
      float v3 = _fcs[ic].v3;

      xc[0][ip] = i1;
      xc[1][ip] = i2;
      xc[2][ip] = i3;

      ws[ip][0] = w1*w1;
      ws[ip][1] = w2*w2;
      ws[ip][2] = w3*w3;
      ws[ip][3] = w1*w2;
      ws[ip][4] = w1*w3;
      ws[ip][5] = w2*w3;

      us[ip][0] = u1*u1;
      us[ip][1] = u2*u2;
      us[ip][2] = u3*u3;
      us[ip][3] = u1*u2;
      us[ip][4] = u1*u3;
      us[ip][5] = u2*u3;

      vs[ip][0] = v1*v1;
      vs[ip][1] = v2*v2;
      vs[ip][2] = v3*v3;
      vs[ip][3] = v1*v2;
      vs[ip][4] = v1*v3;
      vs[ip][5] = v2*v3;
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

    return fc;
  }



  private FaultCell[] getCells(FaultCell[][][] cg) {
    HashSet<FaultCell> hsf = new HashSet<FaultCell>();
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          FaultCell fc = cg[i3][i2][i1];
          if(fc!=null){hsf.add(fc);}
        }
      }
    }
    int nc = hsf.size();
    FaultCell[] fcs = new FaultCell[nc];
    int ic = 0;
    for (FaultCell fc:hsf) {
      fcs[ic] = fc;
      ic++;
    }
    return fcs;
  }

  private void findCells(final int dt,
    final FaultCell[] fc, FaultCell[][][] cg, HashSet<Integer> hsc) 
  {
    int nc = fc.length;
    float sigNor = 2.0f;
    final int[][] bb2 = new int[_n1][2];
    final int[][] bb3 = new int[_n1][2];

    final float[][] xc = new float[3][nc];
    final float[][] ws = new float[nc][6];
    final float[][] us = new float[nc][6];
    final float[][] vs = new float[nc][6];
    final int[] ns = new int[]{_n1,_n2,_n3};
    checkForInterp(fc,xc,ws,us,vs,bb2,bb3);
    final KdTree kt = new KdTree(xc);
    final int[] bs1 = setBounds(_n1,xc[0]);
    final int[] bs2 = setBounds(_n2,xc[1]);
    final int[] bs3 = setBounds(_n3,xc[2]);
    final float st = (float)sin(Math.PI/4.0);
    final float sw = 1.0f/(sigNor*sigNor); 
    final float[][][] fl = new float[_n3][_n2][_n1];
    final float[][][] fp = new float[_n3][_n2][_n1];
    final float[][][] ft = new float[_n3][_n2][_n1];
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
            float fli = fc[ip].fl;
            float wpi = pow(fli,10.f);
            gss += (wd1+wd2+wd3+wds)*sw;
            gss += (ud1+ud2+ud3+uds)*su;
            gss += (vd1+vd2+vd3+vds)*sv;
            float sfi = exp(-gss)*wpi;
            fl[i3][i2][i1] += sfi;
            wps += wpi;
          }
          if(wps!=0f) fl[i3][i2][i1] /= wps;
        }
      }
    }});
    int i1b = bs1[0], i1e = bs1[1];
    int i2b = bs2[0], i2e = bs2[1];
    int i3b = bs3[0], i3e = bs3[1];
    int n1s = i1e-i1b,n2s=i2e-i2b,n3s=i3e-i3b;
    float[][][] sfs = new float[n3s][n2s][n1s];
    float[][][] fps = new float[n3s][n2s][n1s];
    float[][][] fts = new float[n3s][n2s][n1s];
    for (int i3s=0,i3=i3b; i3s<n3s; ++i3,++i3s) 
      for (int i2s=0,i2=i2b; i2s<n2s; ++i2,++i2s) 
        for (int i1s=0,i1=i1b; i1s<n1s; ++i1,++i1s) 
          sfs[i3s][i2s][i1s] = fl[i3][i2][i1];
    smooth(sfs,fps,fts);
    for (int i3s=0,i3=i3b; i3s<n3s; ++i3,++i3s) 
      for (int i2s=0,i2=i2b; i2s<n2s; ++i2,++i2s) 
        for (int i1s=0,i1=i1b; i1s<n1s; ++i1,++i1s){ 
          fl[i3][i2][i1] = sfs[i3s][i2s][i1s];
          fp[i3][i2][i1] = fps[i3s][i2s][i1s];
          ft[i3][i2][i1] = fts[i3s][i2s][i1s];
        }
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    fs.setMinSkinSize(10);
    FaultCell[] fcs = fs.findCells(new float[][][][] {fl,fp,ft});
    removeUsedCells(fcs,hsc);
    for (FaultCell fci:fcs) {
      int i1 = fci.i1;
      int i2 = fci.i2;
      int i3 = fci.i3;

      if (cg[i3][i2][i1]==null) {
        cg[i3][i2][i1] = fci;
      } else {
        float fli = fci.fl;
        float flp = cg[i3][i2][i1].fl;
        if(fli>flp) {
          cg[i3][i2][i1] = fci;
        }
      }
    }
  }

  private void removeUsedCells(FaultCell[] fck, HashSet<Integer> hsc) {
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    int[] ds = new int[]{3,3,3};
    int[] ns = new int[]{_n1,_n2,_n3};
    float[][] xc = setKdTreeCoords(fck);
    KdTree kt = new KdTree(xc);
    HashSet<Integer> hst = new HashSet<Integer>();
    for (int ic:hsc) hst.add(ic);
    for (int ic:hst) {
      int i1 = _fcs[ic].i1;
      int i2 = _fcs[ic].i2;
      int i3 = _fcs[ic].i3;
      int[] is = new int[]{i1,i2,i3};
      getRange(ds,is,ns,xmin,xmax);
      int[] id = kt.findInRange(xmin,xmax);
      int nd = id.length;
      if(nd>ds[0]*ds[0]*2){hsc.remove(ic);}
    }
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


  private void checkForInterp(FaultCell[] fc, float[][] xc, 
    float[][] ws, float[][] us, float[][] vs, int[][] bb2, int[][] bb3) {
    int dd = 30;
    int nc = fc.length;
    for (int i1=0; i1<_n1; ++i1) {
      bb2[i1][0] = _n2; bb2[i1][1] = -_n2;
      bb3[i1][0] = _n3; bb3[i1][1] = -_n3;
    }

    for (int ic=0; ic<nc; ++ic) {
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
      bb2[i1][0] -= 5;
      bb3[i1][0] -= 5;
      bb2[i1][1] += 5;
      bb3[i1][1] += 5;
    }

    /*
    for (int ic=0; ic<nc; ++ic) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      float w1 = fc[ic].w1;
      float w2 = fc[ic].w2;
      float w3 = fc[ic].w3;

      int i1m = i1-dd; if(i1m<0){i1m=0;}
      int i2m = i2-dd; if(i2m<0){i2m=0;}
      int i3m = i3-dd; if(i3m<0){i3m=0;}
      int i1p = i1+dd; if(i1p>=_n1){i1p=_n1-1;}
      int i2p = i2+dd; if(i2p>=_n2){i2p=_n2-1;}
      int i3p = i3+dd; if(i3p>=_n3){i3p=_n3-1;}
      for (int k3=i3m; k3<=i3p; ++k3) {
        for (int k2=i2m; k2<=i2p; ++k2) {
          for (int k1=i1m; k1<=i1p; ++k1) {
            if((k2<bb2[k1][0]||k2>bb2[k1][1])){continue;}
            if((k3<bb3[k1][0]||k3>bb3[k1][1])){continue;}
            float d1 = k1-i1;
            float d2 = k2-i2;
            float d3 = k3-i3;
            float ds = abs(d1*w1+d2*w2+d3*w3);
            if(ds<5f) {mk[k3][k2][k1] = 1;}
          }
        }
      }
    }
    return mk;
    */
  }

  private FaultCell[] findStrike(HashSet<Integer> hsc) {
    int dp = 179;
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
      hsc.remove(ic);
      fc[ik] = _fcs[ic];
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
      xc[0][ik] = _fcs[ic].fp;
    }
    return xc;
  }


  private int[][] getBounds(int[] id) {
    int np = id.length;
    int[][] bs = new int[2][3];
    bs[0][0] = _n1; bs[1][0] = 0;
    bs[0][1] = _n2; bs[1][1] = 0;
    bs[0][2] = _n3; bs[1][2] = 0;
    for (int ip=0; ip<np; ++ip) {
      int ic = id[ip];
      int i1 = _fcs[ic].i1;
      int i2 = _fcs[ic].i2;
      int i3 = _fcs[ic].i3;
      if(i1<bs[0][0]) bs[0][0]=i1;
      if(i2<bs[0][1]) bs[0][1]=i2;
      if(i3<bs[0][2]) bs[0][2]=i3;
      if(i1>bs[1][0]) bs[1][0]=i1;
      if(i2>bs[1][1]) bs[1][1]=i2;
      if(i3>bs[1][2]) bs[1][2]=i3;
    }
    int db = 5;
    bs[0][0] -= db; bs[1][0] += db;
    bs[0][1] -= db; bs[1][1] += db;
    bs[0][2] -= db; bs[1][2] += db;
    if(bs[0][0]<0) bs[0][0] = 0;
    if(bs[0][1]<0) bs[0][1] = 0;
    if(bs[0][2]<0) bs[0][2] = 0;
    if(bs[1][0]>=_n1) bs[1][0] = _n1-1;
    if(bs[1][1]>=_n2) bs[1][1] = _n2-1;
    if(bs[1][2]>=_n3) bs[1][2] = _n3-1;
    return bs;
  }



  private float strikeDif(float fp1, float fp2) {
    float dfp = abs(fp1-fp2);
    if(dfp>(360-30)){dfp=abs(360-dfp);}
    return dfp;
  }

  private float[][][][] smooth(float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float sigma1 = 8.0f;
    float sigma2 = 2.0f;
    float sigma  = 40.f;
    float c = sigma*sigma*0.5f;
    float[][][] u1  = new float[n3][n2][n1];
    float[][][] u2  = new float[n3][n2][n1];
    float[][][] u3  = new float[n3][n2][n1];
    float[][][] fls = new float[n3][n2][n1];
    float[][][] fps = new float[n3][n2][n1];
    float[][][] fts = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    EigenTensors3 d = lof.applyForTensors(fl);
    d.setEigenvalues(0.001f,1.0f,1.0f);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.1,20);
    lsf.apply(d,c,fl,fls);
    div(fls,max(fls),fls); 
    lof.applyForNormal(fls,u1,u2,u3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float fli = fls[i3][i2][i1];
          if(u2i!=0f && u3i!=0f && fli>0.1f) {
            fts[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
            fps[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
          }
        }
      }
    }
    return new float[][][][]{fls,fps,fts};
  }




  /*

  private HashSet<FaultCell[]> sortFaultCells() {
    HashSet<FaultCell> hfc = new HashSet<FaultCell>();
    for (FaultCell fc:_fcs){hfc.add(fc);}
    int nc = hfc.size();
    while(nc>100) {

    }

    final KdTree kt = setStrikeKdTree(hfc,id);
    final HashSet<Inte> hsf = new HashSet<FaultCell[]>();
    Parallel.loop(0,360,1,new Parallel.LoopInt() {
    public void compute(int fp) {
      FaultCell[] fcs = findStrikeInWindow(fp,kt);
      if(fcs!=null){hsf.add(fcs);}
    }});
    return hsf;
    //return hsf.toArray(new FaultCell[hsf.size()][]);
  }

  private KdTree setStrikeKdTree() {
    int nc = _fcs.length;
    float[][] xc = new float[1][nc];
    for (int ic=0; ic<nc; ++ic) {
      xc[0][ic] = _fcs[ic].fp;
    }
    return new KdTree(xc);
  }

  private FaultCell[] cellsInStrikeWd(HashSet<Integer> hsc) {
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
  */

  private KdTree setStrikeKdTree() {
    int nc = _fcs.length;
    float[][] pc = new float[1][nc];
    for (int ic=0; ic<nc; ++ic) {
      pc[0][ic] = _fcs[ic].fp; 
    }
    return new KdTree(pc);
  }

  private int _n1, _n2, _n3;
  private FaultCell[] _fcs;

  private static class OverlapWinds {
    public OverlapWinds(int n1, int l1, double f1) 
    {
      Check.argument(0.0<=f1 && f1<1.0,"0 <= f1 < 1");
      _n1 = n1;
      _l1 = min(l1,n1);
      _m1 = 1+(int)ceil((_n1-_l1)/(_l1*(1.0-f1)));
      _s1 = (double)(_n1-_l1)/max(1,_m1-1);
      makeWeights();
      makeScalars();
    }
    public int getN1() { return _n1; }
    public int getL1() { return _l1; }
    public int getM1() { return _m1; }
    public int getI1(int k1) { return (int)(k1*_s1+0.5); }
    public float getWeight(int i1, int j1) {
      return _w[j1]*_s[i1+j1];
    }
    public float[] getWeights() { return _w; }
    public float[] getScalars() { return _s; }

    private void makeWeights() {
      _w = new float[_l1];
      for (int i1=0; i1<_l1; ++i1) {
        double s1 = sin((i1+1.0)*PI/(_l1+1.0));
        _w[i1] = (float)(s1*s1);
      }
    }

    private void makeScalars() {
      _s = new float[_n1];
      for (int k1=0; k1<_m1; ++k1) {
        int i1 = getI1(k1);
        for (int j1=0; j1<_l1; ++j1) {
          _s[i1+j1] += _w[j1];
        }
      }
      for (int i1=0; i1<_n1; ++i1) {
        _s[i1] = 1.0f/_s[i1];
      }
    }

    private int _n1; // numbers of samples
    private int _l1; // window lengths
    private int _m1; // numbers of windows
    private double _s1; // nominal window spacings
    private float[] _w; // weights[l1] for windowing
    private float[] _s; // scalars[n1] for normalization
    private int _sigNor;
    private int _sigPhi;
    private int _sigTheta;
  }

}


