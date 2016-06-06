/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import util.*;
import static mef.FaultGeometry.*;

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

  public FaultSkin[] reskin(int minSkinSize, float dpi) {
    HashSet<Integer> hsc = new HashSet<Integer>();
    for (int ic=0; ic<_fcs.length; ++ic) hsc.add(ic);
    HashSet<FaultSkin> hsk = new HashSet<FaultSkin>();
    while(hsc.size()>minSkinSize/10) {
      System.out.println("cells remaining:"+hsc.size());
      FaultCell[] fcs = findStrike(20,hsc,dpi);
      FaultSkin[] sks = reskin(minSkinSize,fcs);
      int ns = sks.length;
      if (ns<1) {continue;}
      for (FaultSkin ski:sks) {
        if(ski.size()>=minSkinSize)
          hsk.add(ski);
      }
    }
    return hsk.toArray(new FaultSkin[0]);
  }

  private FaultSkin[] reskin(int minSkinSize, FaultCell[] cells) {
    float[][][][] flpt = applyVote(_n1,_n2,_n3,cells);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.2f,0.5f);
    fs.setMaxDeltaStrike(10);
    fs.setMaxPlanarDistance(0.2f);
    fs.setMinSkinSize(minSkinSize);
    FaultCell[] fcs = fs.findCells(flpt);
    return fs.findSkins(fcs);
  }

  public float[][][][] applyVote(
    final int n1, final int n2, final int n3, FaultCell[] cells) {
    int nc = cells.length;
    final float[] fs = new float[nc]; 
    final float[][] xs = new float[3][nc];
    final float[][] us = new float[3][nc];
    setKdTreeNodes(cells,xs,us,fs);
    final KdTree kt = new KdTree(xs);
    final float[][][] ws = gaussianWeights();
    final float wm = max(ws);
    final float[][][] ss = new float[n3][n2][n1];
    final float[][][] fp = new float[n3][n2][n1];
    final float[][][] ft = new float[n3][n2][n1];
    final int[] c = new int[1];
    final int b3 = (int)min(xs[2]);
    final int e3 = (int)max(xs[2]);
    final int b2 = (int)min(xs[1]);
    final int e2 = (int)max(xs[1]);
    final int b1 = (int)min(xs[0]);
    final int e1 = (int)max(xs[0]);
    Parallel.loop(b3,e3,1,new Parallel.LoopInt() {
    public void compute(int i3) {
      c[0] = c[0]+1;
      System.out.println("c="+c[0]+"/"+n3);
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=b2; i2<e2; ++i2) {
      for (int i1=b1; i1<e1; ++i1) {
        int i1min = max(i1-_w1,0);
        int i2min = max(i2-_w2,0);
        int i3min = max(i3-_w3,0);
        int i1max = min(i1+_w1,n1-1);
        int i2max = min(i2+_w2,n2-1);
        int i3max = min(i3+_w3,n3-1);
        xmin[0] = i1min; xmax[0] = i1max;
        xmin[1] = i2min; xmax[1] = i2max;
        xmin[2] = i3min; xmax[2] = i3max;
        int[] id = kt.findInRange(xmin,xmax);
        int nd = id.length;
        float g11 = 0.0f;
        float g22 = 0.0f;
        float g33 = 0.0f;
        float g12 = 0.0f;
        float g13 = 0.0f;
        float g23 = 0.0f;
        if(nd<1) {
          g11 = wm;
          g22 = wm;
          g33 = wm;
          continue;
        }
        double[][] a = new double[3][3];
        for (int ik=0; ik<nd; ++ik) {
          int ip = id[ik];
          float fl = fs[ip];
          float x1 = xs[0][ip];
          float x2 = xs[1][ip];
          float x3 = xs[2][ip];
          float u1 = us[0][ip];
          float u2 = us[1][ip];
          float u3 = us[2][ip];
          float r1 = i1-x1;
          float r2 = i2-x2;
          float r3 = i3-x3;
          float rs = sqrt(r1*r1+r2*r2+r3*r3);
          int k1 = round(abs(r1));
          int k2 = round(abs(r2));
          int k3 = round(abs(r3));
          float wsi = ws[k3][k2][k1]*fl;
          if (rs==0) {
            g11 += u1*u1*wsi;
            g12 += u1*u2*wsi;
            g13 += u1*u3*wsi;
            g22 += u2*u2*wsi;
            g23 += u2*u3*wsi;
            g33 += u3*u3*wsi;
          } else {
            r1 /= rs; r2 /= rs; r3 /= rs;
            float ur = u1*r1+u2*r2+u3*r3;
            if (abs(ur)>0.5f){continue;}
            float v1 = u1;
            float v2 = u2;
            float v3 = u3;
            if(abs(ur)>0.0001f) {
              float cx = 0.5f*rs/ur; // find a better way?
              float c1 = x1+u1*cx;
              float c2 = x2+u2*cx;
              float c3 = x3+u3*cx;
              v1 = c1-i1; v2 = c2-i2; v3 = c3-i3;
              float vs = 1.0f/sqrt(v1*v1+v2*v2+v3*v3);
              v1 *= vs; v2 *= vs; v3 *= vs; 
            }
            ur = 1f-ur*ur;
            ur *= ur;
            ur *= ur;
            float sc = wsi*ur;
            g11 += sc*v1*v1;
            g12 += sc*v1*v2;
            g13 += sc*v1*v3;
            g22 += sc*v2*v2;
            g23 += sc*v2*v3;
            g33 += sc*v3*v3;
            ss[i3][i2][i1] += sc;
          }
        }
        a[0][0] = g11; a[0][1] = g12; a[0][2] = g13;
        a[1][0] = g12; a[1][1] = g22; a[1][2] = g23;
        a[2][0] = g13; a[2][1] = g23; a[2][2] = g33;
        float[][] ue = solveEigenproblems(a);
        float u1i = ue[0][0];
        float u2i = ue[0][1];
        float u3i = ue[0][2];
        if (u2i==0.0f&&u3i==0f){continue;}
        if (u1i>0.0f) {
          u1i = -u1i;
          u2i = -u2i;
          u3i = -u3i;
        }
        ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
        fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
      }}
    }});
    sub(ss,min(ss),ss);
    div(ss,max(ss),ss);
    return new float[][][][]{ss,fp,ft};
  }



  private void setKdTreeNodes(
    FaultCell[] cells, float[][] xs, float[][] us, float[] fls) {
    int nc = cells.length;
    for (int ic=0; ic<nc; ic++) {
      FaultCell fc = cells[ic];
      fls[ic] = fc.getFl();
      xs[0][ic] = fc.getI1();
      xs[1][ic] = fc.getI2();
      xs[2][ic] = fc.getI3();
      us[0][ic] = fc.getW1();
      us[1][ic] = fc.getW2();
      us[2][ic] = fc.getW3();
    }
  }

  private float[][] solveEigenproblems(double[][] a) {
    double[] e = new double[3];
    double[][] z = new double[3][3];
    Eigen.solveSymmetric33(a,z,e);
    float eui = (float)e[0];
    float evi = (float)e[1];
    float ewi = (float)e[2];
    float u1i = (float)z[0][0];
    float u2i = (float)z[0][1];
    float u3i = (float)z[0][2];
    if (ewi<0.0f) ewi = 0.0f;
    if (evi<ewi) evi = ewi;
    if (eui<evi) eui = evi;
    float[] es = new float[]{eui,evi,ewi};
    float[] us = new float[]{u1i,u2i,u3i};
    return new float[][]{us,es};
  }

  private float[][][] gaussianWeights() {
    int n1 = _w1*2+1;
    int n2 = _w2*2+1;
    int n3 = _w3*2+1;
    float[][][] fx = new float[n3][n2][n1];
    fx[_w3][_w2][_w1] = 1f;
    float[][][] fs = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(_sigma);
    rgf.apply000(fx,fs);
    return copy(_w1+1,_w2+1,_w3+1,_w1,_w2,_w3,fs);
  }

  private FaultCell[] findStrike(int dp, HashSet<Integer> hsc, float dpi) {
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
    for (float pt=0f; pt<=360f; pt+=dpi) {
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

  private int _w1 = 20;
  private int _w2 = 20;
  private int _w3 = 20;
  private float _sigma=20f;
  private int _n1, _n2, _n3;
  private FaultCell[] _fcs;

}


