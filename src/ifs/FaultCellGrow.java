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
 * Regrid fault cells from known cells.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.24
 */

public class FaultCellGrow {

  public FaultCellGrow(FaultCell[] fc, float[][][] fl) {
    _fc = fc; 
    _fl = fl;
    _kt = setKdTree();
    _n3 = fl.length;
    _n2 = fl[0].length;
    _n1 = fl[0][0].length;
  }

  public FaultCell[] combineCells(FaultCell[] fc1, FaultCell[] fc2) {
    int ic = 0;
    int nc = fc1.length+fc2.length;
    FaultCell[] cells = new FaultCell[nc];
    for (FaultCell fc:fc1) {
      cells[ic] = fc;
      ic++;
    }
    for (FaultCell fc:fc2) {
      cells[ic] = fc;
      ic++;
    }
    return cells;

  }

  public FaultCell[] applyForCells(FaultCell cell) {
    float sigNor = 2.0f;
    final float[] da = new float[1];
    final FaultCell[] cells = findNabors(da,cell);
    if(cells==null) {return null;}
    final int nc = cells.length;
    final float su = 0.25f/da[0];
    final float sv = 0.25f/da[0];
    final float sw = 1.0f/(sigNor*sigNor);
    final int[] bs1 = new int[2];
    final int[] bs2 = new int[2];
    final int[] bs3 = new int[2];
    defineBox(cell,bs1,bs2,bs3);
    final int n1 = bs1[1]-bs1[0]+1;
    final int n2 = bs2[1]-bs2[0]+1;
    final int n3 = bs3[1]-bs3[0]+1;
    final float[][][] fls = new float[n3][n2][n1];
    final float[][][] g11 = new float[n3][n2][n1];
    final float[][][] g12 = new float[n3][n2][n1];
    final float[][][] g13 = new float[n3][n2][n1];
    final float[][][] g22 = new float[n3][n2][n1];
    final float[][][] g23 = new float[n3][n2][n1];
    final float[][][] g33 = new float[n3][n2][n1];
    final float x1 = cell.x1;
    final float x2 = cell.x2;
    final float x3 = cell.x3;
    final float v1 = cell.v1;
    final float v2 = cell.v2;
    final float v3 = cell.v3;
    final float u1 = cell.u1;
    final float u2 = cell.u2;
    final float u3 = cell.u3;
    final float w1 = cell.w1;
    final float w2 = cell.w2;
    final float w3 = cell.w3;
    Parallel.loop(n3,new Parallel.LoopInt(){
    public void compute(int k3) {
    int i3 = k3+bs3[0];
    for (int k2=0; k2<n2; ++k2) {
    int i2 = k2+bs2[0];
    for (int k1=0; k1<n1; ++k1) {
    int i1 = k1+bs1[0];
      float d1 = i1-x1; 
      float d2 = i2-x2; 
      float d3 = i3-x3; 
      float dv = abs(d1*v1+d2*v2+d3*v3);
      float du = abs(d1*u1+d2*u2+d3*u3);
      float dw = abs(d1*w1+d2*w2+d3*w3);
      if(dv>3){continue;}
      if(du>3){continue;}
      if(dw>5){continue;}
      float wps = 0.0f;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fci = cells[ic];
        float dx1 = fci.x1-i1;
        float dx2 = fci.x2-i2;
        float dx3 = fci.x3-i3;
        float d11 = dx1*dx1;
        float d22 = dx2*dx2;
        float d33 = dx3*dx3;
        float d12 = dx1*dx2;
        float d13 = dx1*dx3;
        float d23 = dx2*dx3;
        float w11 = fci.w11;
        float w12 = fci.w12;
        float w13 = fci.w13;
        float w22 = fci.w22;
        float w23 = fci.w23;
        float w33 = fci.w33;
        float u11 = fci.u11;
        float u12 = fci.u12;
        float u13 = fci.u13;
        float u22 = fci.u22;
        float u23 = fci.u23;
        float u33 = fci.u33;
        float v11 = fci.v11;
        float v12 = fci.v12;
        float v13 = fci.v13;
        float v22 = fci.v22;
        float v23 = fci.v23;
        float v33 = fci.v33;
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
        float wpi = fci.fl;//pow(fci.fl,10f);
        wps += wpi;
        gss += (wd1+wd2+wd3+wds)*sw;
        gss += (ud1+ud2+ud3+uds)*su;
        gss += (vd1+vd2+vd3+vds)*sv;
        float fli = exp(-gss)*wpi;
        fls[k3][k2][k1] += fli;
        g11[k3][k2][k1] += fli*w11;
        g12[k3][k2][k1] += fli*w12;
        g13[k3][k2][k1] += fli*w13;
        g22[k3][k2][k1] += fli*w22;
        g23[k3][k2][k1] += fli*w23;
        g33[k3][k2][k1] += fli*w33;
      }
      fls[k3][k2][k1] /= wps;
    }}}});
    return findCells(bs1,bs2,bs3,fls,g11,g12,g13,g22,g23,g33);
  }

  /*
  public FaultCell[] applyForCellsM(FaultCell cell) {
    float sigNor = 1.0f;
    final float[] da = new float[1];
    final FaultCell[] cells = findNabors(da,cell);
    if(cells==null) {return null;}
    final int nc = cells.length;
    final float su = 0.25f/da[0];
    final float sv = 0.25f/da[0];
    final float sw = 1.0f/(sigNor*sigNor);
    final int[] bs1 = new int[2];
    final int[] bs2 = new int[2];
    final int[] bs3 = new int[2];
    defineBox(cell,bs1,bs2,bs3);
    final int n1 = bs1[1]-bs1[0]+1;
    final int n2 = bs2[1]-bs2[0]+1;
    final int n3 = bs3[1]-bs3[0]+1;
    final float x1 = cell.x1;
    final float x2 = cell.x2;
    final float x3 = cell.x3;
    final float v1 = cell.v1;
    final float v2 = cell.v2;
    final float v3 = cell.v3;
    final float u1 = cell.u1;
    final float u2 = cell.u2;
    final float u3 = cell.u3;
    final float w1 = cell.w1;
    final float w2 = cell.w2;
    final float w3 = cell.w3;
    final float[][][] fls = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt(){
    public void compute(int k3) {
    int i3 = k3+bs3[0];
    for (int k2=0; k2<n2; ++k2) {
    int i2 = k2+bs2[0];
    for (int k1=0; k1<n1; ++k1) {
    int i1 = k1+bs1[0];
      float d1 = i1-x1; 
      float d2 = i2-x2; 
      float d3 = i3-x3; 
      float dv = abs(d1*v1+d2*v2+d3*v3);
      float du = abs(d1*u1+d2*u2+d3*u3);
      float dw = abs(d1*w1+d2*w2+d3*w3);
      if(dv>3){continue;}
      if(du>3){continue;}
      if(dw>8){continue;}
      float wps = 0.0f;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fci = cells[ic];
        float dx1 = fci.x1-i1;
        float dx2 = fci.x2-i2;
        float dx3 = fci.x3-i3;
        float d11 = dx1*dx1;
        float d22 = dx2*dx2;
        float d33 = dx3*dx3;
        float d12 = dx1*dx2;
        float d13 = dx1*dx3;
        float d23 = dx2*dx3;
        float w11 = fci.w11;
        float w12 = fci.w12;
        float w13 = fci.w13;
        float w22 = fci.w22;
        float w23 = fci.w23;
        float w33 = fci.w33;
        float u11 = fci.u11;
        float u12 = fci.u12;
        float u13 = fci.u13;
        float u22 = fci.u22;
        float u23 = fci.u23;
        float u33 = fci.u33;
        float v11 = fci.v11;
        float v12 = fci.v12;
        float v13 = fci.v13;
        float v22 = fci.v22;
        float v23 = fci.v23;
        float v33 = fci.v33;
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
        float wpi = pow(fci.fl,10f);
        wps += wpi;
        gss += (wd1+wd2+wd3+wds)*sw;
        gss += (ud1+ud2+ud3+uds)*su;
        gss += (vd1+vd2+vd3+vds)*sv;
        float fli = exp(-gss)*wpi;
        fls[k3][k2][k1] += fli;
      }
      fls[k3][k2][k1] /= wps;
    }}}});
    return findCells(bs1,bs2,bs3,fls);
  }

  private FaultCell[] findCells(
    int[] bs1, int[] bs2, int[] bs3, float[][][] fl) 
  {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] fp = new float[n3][n2][n1];
    float[][][] ft = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8,4);
    lof.applyForNormal(fl,u1,u2,u3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = -u1[i3][i2][i1];
      float u2i = -u2[i3][i2][i1];
      float u3i = -u3[i3][i2][i1];
      if(u2i!=0.0f && u3i!=0.0f) {
        ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
        fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
      }
    }}}
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    FaultCell[] fcs = fs.findCells(new float[][][][]{fl,fp,ft});
    int nc = fcs.length; if(nc<1){return null;}
    for (int ic=0; ic<nc; ++ic) {
      float x1i = fcs[ic].x1+bs1[0];
      float x2i = fcs[ic].x2+bs2[0];
      float x3i = fcs[ic].x3+bs3[0];
      int i1 = round(x1i);
      int i2 = round(x2i);
      int i3 = round(x3i);
      if(i1<0) {continue;}
      if(i2<0) {continue;}
      if(i3<0) {continue;}
      if(i1>=_n1) {continue;}
      if(i2>=_n2) {continue;}
      if(i3>=_n3) {continue;}
      float fpi = fcs[ic].fp;
      float fti = fcs[ic].ft;
      float fli = _fl[i3][i2][i1];
      fcs[ic] = new FaultCell(x1i,x2i,x3i,fli,fpi,fti);
    }
    return fcs;
  }
  */

  private FaultCell[] findCells(
    int[] bs1, int[] bs2, int[] bs3, float[][][] fl, 
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33) 
  {
    div(fl,max(fl),fl);
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] fp = new float[n3][n2][n1];
    float[][][] ft = new float[n3][n2][n1];
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    faultImages(fp,ft,u1,u2,u3);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    FaultCell[] fcs = fs.findCells(new float[][][][]{fl,fp,ft});
    int nc = fcs.length; if(nc<1){return null;}
    for (int ic=0; ic<nc; ++ic) {
      float x1i = fcs[ic].x1+bs1[0];
      float x2i = fcs[ic].x2+bs2[0];
      float x3i = fcs[ic].x3+bs3[0];
      int i1 = round(x1i);
      int i2 = round(x2i);
      int i3 = round(x3i);
      if(i1<0) {continue;}
      if(i2<0) {continue;}
      if(i3<0) {continue;}
      if(i1>=_n1) {continue;}
      if(i2>=_n2) {continue;}
      if(i3>=_n3) {continue;}
      float fpi = fcs[ic].fp;
      float fti = fcs[ic].ft;
      float fli = _fl[i3][i2][i1];
      //float flt = _fl[i3][i2][i1];
      //float fli = exp(flt)/exp(1f);//_fl[i3][i2][i1];
      fcs[ic] = new FaultCell(x1i,x2i,x3i,fli,fpi,fti);
    }
    return fcs;
  }

  private void solveEigenproblems(
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33,
    float[][][] u1,  float[][][] u2,  float[][][] u3) 
  {
    int n3 = g11.length;
    int n2 = g11[0].length;
    int n1 = g11[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      double[] e = new double[3];
      double[][] a = new double[3][3];
      double[][] z = new double[3][3];
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
  }

  private void faultImages(
    final float[][][] fp, final float[][][] ft,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = fp.length;
    final int n2 = fp[0].length;
    final int n1 = fp[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      if(u2i!=0.0f && u3i!=0.0f) {
        ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
        fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
      }
    }}}
  }

  private void defineBox(FaultCell fc, int[] bs1, int[] bs2, int[] bs3) {
    float d = 15f;
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

  private FaultCell[] findNabors(FaultCell cell) {
    int dd = 10;
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    xmin[0] = cell.x1-dd;
    xmin[1] = cell.x2-dd;
    xmin[2] = cell.x3-dd;
    xmax[0] = cell.x1+dd;
    xmax[1] = cell.x2+dd;
    xmax[2] = cell.x3+dd;
    int[] id = _kt.findInRange(xmin,xmax);
    int nd = id.length;
    if(nd<1) {return null;}
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (int ik=0; ik<nd; ++ik) {
      int ic = id[ik];
      FaultCell fc = _fc[ic];
      if(canBeNabors(cell,fc)) {
        hsc.add(fc);
      }
    }
    return getCells(hsc);
  }

  private FaultCell[] findNabors(float[] da, FaultCell cell) {
    float[] dl = new float[1];
    float[] dr = new float[1];
    FaultCell[] fcL = findNaborsL(dl,cell);
    FaultCell[] fcR = findNaborsR(dr,cell);
    if(fcL==null||fcR==null){return null;}
    float nl = fcL.length, nr = fcR.length;
    //if(nl<2||nr<2){return null;}
    da[0] = (dl[0]+dr[0])*0.5f;
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (int il=0; il<nl; ++il){
      FaultCell fc = fcL[il];
      hsc.add(fc);
    }
    for (int ir=0; ir<nr; ++ir){
      FaultCell fc = fcR[ir];
      if(hsc.contains(fc)){continue;}
      hsc.add(fc);
    }
    return getCells(hsc);
  }

  public FaultCell[] findNaborsL(float[] da, FaultCell cell) {
    int dd = 10;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float w2 = abs(cell.w2);
    float w3 = abs(cell.w3);
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    int nd = 0;
    int[] id = null;
    xmax[2] = x3+dd;
    xmax[1] = x2+dd;
    xmax[0] = x1+dd;
    xmin[2] = x3-dd;
    xmin[1] = x2-dd;
    xmin[0] = x1-dd;
    while(nd<30&&dd<20) { 
      if(w2>=w3) {
        xmax[2]  = x3;
        xmin[2] -= 2f;
      } else {
        xmax[1]  = x2;
        xmin[1] -= 2f;
      }
      dd +=2;
      id = _kt.findInRange(xmin,xmax);
      nd = id.length;
    }
    if(nd<1) {return null;}
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (int ik=0; ik<nd; ++ik) {
      int ic = id[ik];
      FaultCell fc = _fc[ic];
      if(canBeNabors(cell,fc)) {
        hsc.add(fc);
      }
    }
    if(hsc.size()<1){
      return null;
    }
    //return getCells(hsc);
    return distance(da,cell,getCells(hsc));
  }

  public FaultCell[] findNaborsR(float[] da, FaultCell cell) {
    int dd = 10;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float w2 = abs(cell.w2);
    float w3 = abs(cell.w3);
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    xmin[2] = x3-dd;
    xmin[1] = x2-dd;
    xmin[0] = x1-dd;
    xmax[2] = x3+dd;
    xmax[1] = x2+dd;
    xmax[0] = x1+dd;
    int nd = 0;
    int[] id = null;
    while(nd<30&&dd<20) { 
      if(w2>=w3) {
        xmin[2]  = x3;
        xmax[2] += 2f;
      } else {
        xmin[1]  = x2;
        xmax[1] += 2f;
      }
      dd +=2;
      id = _kt.findInRange(xmin,xmax);
      nd = id.length;
    }
    if(nd<1) {return null;}
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (int ik=0; ik<nd; ++ik) {
      int ic = id[ik];
      FaultCell fc = _fc[ic];
      if(canBeNabors(cell,fc)) {
        hsc.add(fc);
      }
    }
    if(hsc.size()<1){
      return null;
    }
    //return getCells(hsc);
    return distance(da,cell,getCells(hsc));
  }


  private FaultCell[] distance(
    float[] da, FaultCell fcc, FaultCell[] fcs) 
  {
    int ns = 50;
    float x1 = fcc.x1;
    float x2 = fcc.x2;
    float x3 = fcc.x3;
    int nc = fcs.length;
    int[] id = new int[nc];
    float[] ds = new float[nc];
    float[] dc = new float[nc];
    for (int ic=0; ic<nc; ++ic) {
      id[ic] = ic;
      float d1 = fcs[ic].x1-x1;
      float d2 = fcs[ic].x2-x2;
      float d3 = fcs[ic].x3-x3;
      ds[ic] = d1*d1+d2*d2+d3*d3;
      dc[ic] = ds[ic]*abs(d1);
    }
    if(nc<ns){da[0]=sum(ds)/nc;return fcs;}
    else {
      float dd = 0.0f;
      quickIndexSort(dc,id);
      FaultCell[] cells = new FaultCell[ns];
      for (int ik=0; ik<ns; ++ik) {
        int ic = id[ik];
        dd += ds[ic];
        cells[ik] = fcs[ic];
      }
      da[0] = dd/ns;
      return cells;
    }
  }


  private FaultCell[] getCells(HashSet<FaultCell> hsc) {
    int nc = hsc.size();
    if(nc==0){return null;}
    int ic = 0;
    FaultCell[] fcs = new FaultCell[nc];
    for (FaultCell fc:hsc) {
      fcs[ic] = fc;
      ic++;
    }
    return fcs;
  }

  private boolean canBeNabors(FaultCell ci, FaultCell cn) {
    boolean can = true;
    if (cn.fl<_fllo) {
      can = false;
    } else if (absDeltaFp(ci,cn)>_dfpmax) {
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

  private KdTree setKdTree() {
    int nc = _fc.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xc[0][ic] = _fc[ic].x1;
      xc[1][ic] = _fc[ic].x2;
      xc[2][ic] = _fc[ic].x3;
    }
    return new KdTree(xc);
  }

  
  private KdTree _kt;
  private FaultCell[] _fc;
  private float[][][] _fl;
  private int _n1, _n2, _n3;
  private float _fllo=0.2f;
  //private float _dfpmax=15f; // max difference between strikes of nabors
  private float _dfpmax=40f; // max difference between strikes of nabors
  private float _dftmax=10f; // max difference between dips of nabors
  private float _dnpmax=5.f; // max distance to planes of nabors
}

