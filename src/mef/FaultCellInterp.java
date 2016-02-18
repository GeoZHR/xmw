/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

import util.*;
import static mef.FaultGeometry.*;

/**
 * Interpolate new fault cells from known cells, not finished.......
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.02.17
 */

public class FaultCellInterp {

  public void setParameters(float dfp, float dft, float dnp) {
    _dfpmax = dfp;
    _dftmax = dft;
    _dnpmax = dnp;
  }

  public FaultCell[] interp(int n1, int n2, int n3, FaultSkin[] skins) {

    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fc = FaultSkin.getCells(skins);
    setKdTree();
    // Grid of cells used to quickly find cell nabors.
    FaultCellGrid cellGrid = new FaultCellGrid(n1,n2,n3);

    // Make a list of cells with missing nabors.
    ArrayList<FaultCell> seedList = new ArrayList<FaultCell>();
    for (FaultSkin skin:skins) {
    for (FaultCell cell:skin) {
      cellGrid.set(cell);
      if (cell.ca==null) {seedList.add(cell);continue;}
      if (cell.cb==null) {seedList.add(cell);continue;}
      if (cell.cr==null) {seedList.add(cell);continue;}
      if (cell.cl==null) {seedList.add(cell);continue;}
    }}

    // While potential seeds remain, ...
    while (seedList.size()>0) {
      FaultCell[] seeds = seedList.toArray(new FaultCell[0]);
      // clear the seed list for adding new seeds with missing nabors
      seedList.clear(); 
      int nseed = seeds.length;
      for (int kseed=0; kseed<nseed; ++kseed) {
        FaultCell seed = seeds[kseed];
      }
    }

    return null;
  }

  public FaultCell[] interpNabors(float sigma, FaultCell cell) {
    FaultCell[] cells = nabors(cell);
    if(cells==null) {return null;}
    int nc = cells.length;
    int[] bs1 = new int[2];
    int[] bs2 = new int[2];
    int[] bs3 = new int[2];
    defineBox(cell,bs1,bs2,bs3);
    int n1 = bs1[1]-bs1[0]+1;
    int n2 = bs2[1]-bs2[0]+1;
    int n3 = bs3[1]-bs3[0]+1;
    float[][][] g11 = new float[n3][n2][n1];
    float[][][] g12 = new float[n3][n2][n1];
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g22 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    float[][][] g33 = new float[n3][n2][n1];
    float[][][] gw = gaussWeight(_w1,_w2,_w3,sigma);
    for (int k3=0; k3<n3; ++k3) {
    int i3 = k3+bs3[0];
    for (int k2=0; k2<n2; ++k2) {
    int i2 = k2+bs2[0];
    for (int k1=0; k1<n1; ++k1) {
    int i1 = k1+bs1[0];
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fci = cells[ic];
        float w1 = fci.w1;
        float w2 = fci.w2;
        float w3 = fci.w3;
        int r1 = fci.i1-i1;
        int r2 = fci.i2-i2;
        int r3 = fci.i3-i3;
        float rs = sqrt(r1*r1+r2*r2+r3*r3);
        float sc = gw[abs(r1)][abs(r2)][abs(r3)];
        if (rs>0f) {
          r1 /= rs; r2 /= rs; r3 /= rs;
          float wr = w1*r1+w2*r2+w3*r3;
          wr = 1-wr*wr;
          wr *= wr;
          wr *= wr;
          wr *= wr;
          wr *= wr;
          wr *= wr;
          sc *= wr;
        }
        g11[k3][k2][k1] += w1*w1*sc;
        g12[k3][k2][k1] += w1*w2*sc;
        g13[k3][k2][k1] += w1*w3*sc;
        g22[k3][k2][k1] += w2*w2*sc;
        g23[k3][k2][k1] += w2*w3*sc;
        g33[k3][k2][k1] += w3*w3*sc;
      }
    }}}
    float[][][][] flpt = solveEigenproblems(g11,g12,g13,g22,g23,g33);
    return null;
    //return findcells(0.2f,flpt);
    //return findCells(bs1,bs2,bs3,fls,g11,g12,g13,g22,g23,g33);
  }


  private float[][][][] solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33) {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    final float[][][] ss = new float[n3][n2][n1];
    final float[][][] fp = new float[n3][n2][n1];
    final float[][][] ft = new float[n3][n2][n1];
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        double[] e = new double[3];
        double[][] z = new double[3][3];
        double[][] a = new double[3][3];
        float g11i = g11[i3][i2][i1];
        float g12i = g12[i3][i2][i1];
        float g13i = g13[i3][i2][i1];
        float g22i = g22[i3][i2][i1];
        float g23i = g23[i3][i2][i1];
        float g33i = g33[i3][i2][i1];
        a[0][0]=g11i; a[0][1]=g12i; a[0][2]=g13i;
        a[1][0]=g12i; a[1][1]=g22i; a[1][2]=g23i;
        a[2][0]=g13i; a[2][1]=g23i; a[2][2]=g33i;
        Eigen.solveSymmetric33(a,z,e);
        float eui = (float)e[0];
        float evi = (float)e[1];
        float u1i = (float)z[0][0];
        float u2i = (float)z[0][1];
        float u3i = (float)z[0][2];
        if (u2i==0.0f&&u3i==0f){continue;}
        if (u1i>0.0f) {
          u1i = -u1i;
          u2i = -u2i;
          u3i = -u3i;
        }
        ss[i3][i2][i1] = (eui-evi);
        ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
        fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
      }}
    }});
    sub(ss,min(ss),ss);
    div(ss,max(ss),ss);
    return new float[][][][]{ss,fp,ft};
  }

  private void defineBox(FaultCell fc, int[] bs1, int[] bs2, int[] bs3) {
    float d = 5f;
    int i1m = round(fc.x1-d);if(i1m<0){i1m=0;}
    int i2m = round(fc.x2-d);if(i2m<0){i2m=0;}
    int i3m = round(fc.x3-d);if(i3m<0){i3m=0;}
    int i1p = round(fc.x1+d);if(i1p>=_n1){i1p=_n1-1;}
    int i2p = round(fc.x2+d);if(i2p>=_n2){i2p=_n2-1;}
    int i3p = round(fc.x3+d);if(i3p>=_n3){i3p=_n3-1;}
    bs1[0] = i1m; bs1[1] = i1p;
    bs2[0] = i2m; bs2[1] = i2p;
    bs3[0] = i3m; bs3[1] = i3p;
  }

  public FaultCell[] nabors(FaultCell cell) {
    float dv = 30f;
    float dh = 30f;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float v2 = cell.v2;
    float v3 = cell.v3;
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    xmin[0] = x1-dv; xmax[0] = x1+dv;
    xmin[1] = x2-dh; xmax[1] = x2+dh;
    xmin[2] = x3-dh; xmax[2] = x3+dh;
    int[] id = _kt.findInRange(xmin,xmax);
    int nd = id.length;
    if(nd<1) {return null;}
    int nbc=0, nac=0;
    ArrayList<Float> dsL = new ArrayList<Float>();
    ArrayList<Float> dsR = new ArrayList<Float>();
    ArrayList<FaultCell> fcL = new ArrayList<FaultCell>();
    ArrayList<FaultCell> fcR = new ArrayList<FaultCell>();
    for (int ik=0; ik<nd; ++ik) {
      int ic = id[ik];
      FaultCell fci = _fc[ic];
      if(canBeNabors(fci,cell)){
        float d1 = fci.x1-x1;
        float d2 = fci.x2-x2;
        float d3 = fci.x3-x3;
        float dd = d2*v2+d3*v3;
        float ds = d1*d1+d2*d2+d3*d3;
        if(abs(d1)>1f) {ds *= abs(d1);}
        if(dd<=0f) {dsL.add(ds);fcL.add(fci);} 
        else       {dsR.add(ds);fcR.add(fci);}
        if(fci.x1>x1){nac++;}
        if(fci.x1<x1){nbc++;}
      }
    }
    int nb = 30;
    int nl = dsL.size();
    int nr = dsR.size();
    if(nbc<5||nac<5) {return null;}
    if(nl<10||nr<10) {return null;}
    ArrayList<FaultCell> nbs = new ArrayList<FaultCell>();
    if(nl<nb) {
      for (FaultCell fc:fcL) {nbs.add(fc);}
    } else {
      int il = 0;
      int[] idl = new int[nl];
      float[] dsl = new float[nl];
      FaultCell[] fcl = fcL.toArray(new FaultCell[0]);
      for (float dsi:dsL) {dsl[il]=dsi; idl[il]=il; il++;}
      quickIndexSort(dsl,idl);
      for (int ik=0; ik<nb; ++ik) {
        int ic = idl[ik];
        nbs.add(fcl[ic]);
      }
    }
    if(nr<nb) {
      for (FaultCell fc:fcR) {nbs.add(fc);}
    } else {
      int ir = 0;
      int[] idr = new int[nr];
      float[] dsr = new float[nr];
      FaultCell[] fcr = fcR.toArray(new FaultCell[0]);
      for (float dsi:dsR) {dsr[ir]=dsi; idr[ir]=ir; ir++;}
      quickIndexSort(dsr,idr);
      for (int ik=0; ik<nb; ++ik) {
        int ic = idr[ik];
        nbs.add(fcr[ic]);
      }
    }
    return nbs.toArray(new FaultCell[0]);
  }

  private boolean canBeNabors(FaultCell ci, FaultCell cn) {
    boolean can = true;
    if (absDeltaFp(ci,cn)>_dfpmax) {
      can = false;
    } else if (absDeltaFt(ci,cn)>_dftmax) {
      can = false;
    } else if (maxDistanceToPlane(ci,cn)>_dnpmax) {
      can = false;
    }
    return can;
  }
  private static float absDeltaFp(FaultCell ca, FaultCell cb) {
    float del = ca.fp-cb.fp;
    return min(abs(del),abs(del+360.0f),abs(del-360.0f));
  }
  private static float absDeltaFt(FaultCell ca, FaultCell cb) {
    return abs(ca.ft-cb.ft);
  }
  private static float maxDistanceToPlane(FaultCell ca, FaultCell cb) {
    float aw1 = ca.w1, aw2 = ca.w2, aw3 = ca.w3;
    float ax1 = ca.x1, ax2 = ca.x2, ax3 = ca.x3;
    float bw1 = cb.w1, bw2 = cb.w2, bw3 = cb.w3;
    float bx1 = cb.x1, bx2 = cb.x2, bx3 = cb.x3;
    float dx1 = ax1-bx1;
    float dx2 = ax2-bx2;
    float dx3 = ax3-bx3;
    float dab = abs(aw1*dx1+aw2*dx2+aw3*dx3);
    float dba = abs(bw1*dx1+bw2*dx2+bw3*dx3);
    return max(dab,dba);
  }


  public float[][][] gaussWeight(int w1, int w2, int w3, float sigma) {
    int n3 = w3*2+1;
    int n2 = w2*2+1;
    int n1 = w1*2+1;
    float[][][] gx = new float[n3][n2][n1];
    float[][][] gs = new float[n3][n2][n1];
    gx[w3][w2][w1] = 1.0f;
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    rgf.apply000(gx,gs);
    return copy(w1+1,w2+1,w3+1,w1,w2,w3,gs);
  }

  private void setKdTree() {
    int nc = _fc.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xc[0][ic] = _fc[ic].x1;
      xc[1][ic] = _fc[ic].x2;
      xc[2][ic] = _fc[ic].x3;
    }
    _kt =  new KdTree(xc);
  }

  private KdTree _kt;
  private FaultCell[] _fc;
  private int _n1, _n2, _n3;
  private int _w1, _w2, _w3; // search box
  private float _dfpmax=20f; // max difference between strikes of nabors
  private float _dftmax=20f; // max difference between dips of nabors
  private float _dnpmax=2.f; // max distance to planes of nabors

}


