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

  public float[][][][] applyForFaultImages() {
    int np = 360;
    int owl = 30;
    double owf = 0.5;
    float[][][] fw = new float[_n3][_n2][_n1];
    float[][][] fl = new float[_n3][_n2][_n1];
    float[][][] fp = new float[_n3][_n2][_n1];
    float[][][] ft = new float[_n3][_n2][_n1];
    OverlapWinds ow = new OverlapWinds(np,owl,owf);
    int m1 = ow.getM1();
    int l1 = ow.getL1();
    KdTree kt = setStrikeKdTree();
    for (int k1=0; k1<m1; ++k1) {
      System.out.println("k1="+k1);
      int p1 = ow.getI1(k1);
      float[] pmin = new float[]{p1};
      float[] pmax = new float[]{p1+l1-1};
      int[] id = kt.findInRange(pmin,pmax);
      int ns = id.length;
      if(ns<10){continue;}
      int[][] bs = getBounds(id);
      int n1s = bs[1][0]-bs[0][0]+1;
      int n2s = bs[1][1]-bs[0][1]+1;
      int n3s = bs[1][2]-bs[0][2]+1;
      float[][][] g = new float[n3s][n2s][n1s];
      for (int is=0; is<ns; ++is) {
        int ic = id[is];
        float fpi = _fcs[ic].fp;
        float fli = _fcs[ic].fl;
        int i1 = _fcs[ic].i1-bs[0][0];
        int i2 = _fcs[ic].i2-bs[0][1];
        int i3 = _fcs[ic].i3-bs[0][2];
        g[i3][i2][i1] = fli*ow.getWeight(p1,round(fpi-p1));
      }
      RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
      rgf.apply0XX(g,g);
      rgf.applyX0X(g,g);
      rgf.applyXX0(g,g);
      div(g,max(g),g);
      float[][][][] flpt = smooth(g);
      for (int j3=0; j3<n3s; ++j3) {
        for (int j2=0; j2<n2s; ++j2) {
          for (int j1=0; j1<n1s; ++j1) {
            int c1 = j1+bs[0][0];
            int c2 = j2+bs[0][1];
            int c3 = j3+bs[0][2];
            float flj = flpt[0][j3][j2][j1];
            float fpj = flpt[1][j3][j2][j1];
            float ftj = flpt[2][j3][j2][j1];
            fl[c3][c2][c1] += flj;
            if(flj>fw[c3][c2][c1]) {
              fw[c3][c2][c1] = flj;
              fp[c3][c2][c1] = fpj;
              ft[c3][c2][c1] = ftj;
            }
          }
        }
      }
    }
    return new float[][][][]{fl,fp,ft}; 
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
    float sigma  = 20.f;
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
  }

}


