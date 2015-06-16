/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package figi;

import ipfx.*;
import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Set weights and constraints for flattening: 
 * set zero weights at fault, 
 * set hard constraints using fault slip vectors.
 * @author Xinming Wu
 * @version 2014.02.09
 */

public class FaultSlipConstraints {

  public FaultSlipConstraints(FaultSkin[] sks) {
    _sks = sks;
  }

  public void setValuesNearFaults(float v, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] mk = new float[n3][n2][n1];
    float[][][] ds = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    for (FaultSkin skin:_sks) {
    for (FaultCell cell:skin) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      mk[i3][i2][i1] = 1.0f;
    }}
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(0f,mk,ds,k1,k2,k3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float dsi = ds[i3][i2][i1];
      if(dsi<=6f&&dsi>0f){fl[i3][i2][i1]=v;}
    }}}
 
  }


  public float[][][] controlPoints(
    float[][][] ws, float[][][] wp, float[][][] cp) {
    int n3 = ws.length;
    int n2 = ws[0].length;
    int n1 = ws[0][0].length;
    computeUnfaultShifts(n1,n2,n3);
    setCellsC(n1,n2,n3);
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin fs:_sks) {
    for (FaultCell fc:fs) {
      float fl = fc.getFl();
      if(fl==0.0f) {continue;}
      float[] hx = new float[3];
      float[] fx = new float[3];
      int i1i = fc.getI1();
      int i2i = fc.getI2();
      int i3i = fc.getI3();
      int i2p = bound(fc.getP2(),n2);
      int i3p = bound(fc.getP3(),n3);
      int i2m = bound(fc.getM2(),n2);
      int i3m = bound(fc.getM3(),n3);
      int d2m = i2m-i2i;
      int d3m = i3m-i3i;
      int d2p = i2p-i2i;
      int d3p = i3p-i3i;
      if(abs(d2m)>0){i2p -= d2m;}
      else          {i2m -= d2p;}
      if(abs(d3m)>0){i3p -= d3m;}
      else          {i3m -= d3p;}
      if(i2m<0||i2m>=n2){continue;}
      if(i2p<0||i2p>=n2){continue;}
      if(i3m<0||i3m>=n3){continue;}
      if(i3p<0||i3p>=n3){continue;}
      ws[i3i][i2i][i1i] = 0.000f;
      wp[i3i][i2i][i1i] = 0.005f;
      wp[i3m][i2m][i1i] = 0.005f;
      wp[i3p][i2p][i1i] = 0.005f;
      cp[i3m][i2m][i1i] = 1.000f;
      cp[i3p][i2p][i1i] = 4.000f;
      hx[0] = i1i; fx[0] = i1i;
      hx[1] = i2p; fx[1] = i2m;
      hx[2] = i3p; fx[2] = i3m;
      float t1 = fc.getT1();
      float t2 = fc.getT2();
      float t3 = fc.getT3();
      float[] ts = new float[]{t1,t2,t3};
      cl.add(new float[][]{hx,fx,mul(ts,0.5f)});
    }}
    int ns = cl.size();
    System.out.println("sets of control points:"+ns);
    float[][][] cs = new float[3][ns][3];
    for (int is=0; is<ns; ++is) {
      float[][] ps = cl.get(is);
      for (int ip=0; ip<3; ++ip) {
        cs[0][is][ip] = ps[ip][0];
        cs[1][is][ip] = ps[ip][1];
        cs[2][is][ip] = ps[ip][2];
      }
    }
    return cs;
  }

  public float[][][] screenPoints(float[][][] wp){
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    computeUnfaultShifts(n1,n2,n3);
    setCells(n1,n2,n3);
    flNormalization();
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin skin:_sks) {
      for (FaultCell fc:skin) {
        float fl = fc.getFl();
        float[] h1 = new float[3];
        float[] f1 = new float[3];
        float[] c1 = new float[3];
        float[] h2 = new float[3];
        float[] f2 = new float[3];
        int i1i  = fc.getI1();
        int i2i  = fc.getI2();
        int i3i  = fc.getI3();
        int i2m1 = bound(fc.getM2(),n2);
        int i3m1 = bound(fc.getM3(),n3);
        int i2p1 = bound(fc.getP2(),n2);
        int i3p1 = bound(fc.getP3(),n3);

        h1[0] = i1i;  f1[0] = i1i;
        h1[1] = i2p1; f1[1] = i2m1;
        h1[2] = i3p1; f1[2] = i3m1;
        wp[i3m1][i2m1][i1i] = 0.0f;
        wp[i3p1][i2p1][i1i] = 0.0f;

        int d2m = i2m1-i2i;
        int d3m = i3m1-i3i;
        int d2p = i2p1-i2i;
        int d3p = i3p1-i3i;
        int i2m2 = i2m1, i2p2=i2p1; 
        int i3m2 = i3m1, i3p2=i3p1; 
        if(abs(d2m)>0){i2p2-=d2m;i2m2+=d2m;}
        else          {i2m2-=d2p;i2p2+=d2p;}
        if(abs(d3m)>0){i3p2-=d3m;i3m2+=d3m;}
        else          {i3m2-=d3p;i3p2+=d3p;}
        i2m2=bound(i2m2,n2); 
        i2p2=bound(i2p2,n2);
        i3m2=bound(i3m2,n3); 
        i3p2=bound(i3p2,n3);
        h2[0] = i1i;  f2[0] = i1i;
        h2[1] = i2p2; f2[1] = i2m2;
        h2[2] = i3p2; f2[2] = i3m2;

        float[] ch = new float[3];
        float[] cf = new float[3];
        float t1 = fc.getT1();
        float t2 = fc.getT2();
        float t3 = fc.getT3();
        float[] ts = new float[]{t1,t2,t3};
        c1[0] = fl; c1[1] = fl; c1[2] = fl;
        cl.add(new float[][]{h1,f1,ts,c1});
        cl.add(new float[][]{h2,h1,ch,c1});
        cl.add(new float[][]{f2,f1,cf,c1});
      }
    }
    int ns = cl.size();
    System.out.println("sets of control points:"+ns);
    float[][][] cs = new float[4][3][ns];
    for (int is=0; is<ns; ++is) {
      float[][] ps = cl.get(is);
      cs[0][0][is] = ps[0][0];
      cs[0][1][is] = ps[0][1];
      cs[0][2][is] = ps[0][2];

      cs[1][0][is] = ps[1][0];
      cs[1][1][is] = ps[1][1];
      cs[1][2][is] = ps[1][2];

      cs[2][0][is] = ps[2][0];
      cs[2][1][is] = ps[2][1];
      cs[2][2][is] = ps[2][2];

      cs[3][0][is] = ps[3][0];
      cs[3][1][is] = ps[3][1];
      cs[3][2][is] = ps[3][2];
    }
    return cs;
  }

  private void computeUnfaultShifts(int n1, int n2, int n3) {
    FaultCell[][][] fcg = new FaultCell[n3][n2][n1];
    for (FaultSkin skin:_sks) {
    for (FaultCell cell:skin) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      fcg[i3][i2][i1] = cell;
    }}
    for (FaultSkin skin:_sks) {
    for (FaultCell cell:skin) {
      float[] sa = cell.getS();
      float x1 = cell.getX1()+cell.getR1();
      float x2 = cell.getX2()+cell.getR2();
      float x3 = cell.getX3()+cell.getR3();
      int k1 = round(x1);
      int k2 = round(x2);
      int k3 = round(x3);
      FaultCell cb = fcg[k3][k2][k1];
      if(cb!=null) {
        float ds1 = cb.getS1()-sa[0];
        float ds2 = cb.getS2()-sa[1];
        float ds3 = cb.getS3()-sa[2];
        sa[0] -= ds1;  //not sure
        sa[1] -= ds2;  //not sure
        sa[2] -= ds3;  //not sure
      }
      cell.setUnfaultShifts(sa);
    }}
    checkUnfaultShifts();
  }

  private void checkUnfaultShifts() {
    for (FaultSkin skin:_sks) {
      checkUnfaultShiftsAB(skin);
      checkUnfaultShiftsLR(skin);
    }
  }

  private void checkUnfaultShiftsAB(FaultSkin skin) {
    FaultCell[][] cells = skin.getCellsAB();
    int na = cells.length;
    for (int ia=0; ia<na; ++ia) {
      FaultCell[] celli = cells[ia];
      int nc = celli.length;
      for (int ic=1; ic<nc; ++ic) {
        FaultCell cb = celli[ic  ];
        FaultCell ca = celli[ic-1];
        float t1b = cb.getT1();
        float t1a = ca.getT1();
        float dba = t1b-t1a+1;
        if(dba<0f){cb.setT1(t1b-dba);}
      }
    }
  }

  private void checkUnfaultShiftsLR(FaultSkin skin) {
    FaultCell[][] cells = skin.getCellsLR();
    int na = cells.length;
    for (int ia=0; ia<na; ++ia) {
      FaultCell[] celli = cells[ia];
      int nc = celli.length;
      for (int ic=1; ic<nc; ++ic) {
        FaultCell cr = celli[ic  ];
        FaultCell cl = celli[ic-1];
        float x2r = cr.getI2();
        float x3r = cr.getI3();
        float t2r = cr.getT2();
        float t3r = cr.getT3();
        float x2l = cl.getI2();
        float x3l = cl.getI3();
        float t2l = cl.getT2();
        float t3l = cl.getT3();
        float dx2 = x2r-x2l;
        float dx3 = x3r-x3l;
        float ds2 = dx2+t2r-t2l;
        float ds3 = dx3+t3r-t3l;
        if(dx2*ds2<0f) cr.setT2(t2r-ds2);
        if(dx3*ds3<0f) cr.setT3(t3r-ds3);
      }
    }
  }

  private int bound(int i, int n) {
    if(i<0) {i=0;}
    if(i>=n){i=n-1;}
    return i;
  }



  private void flNormalization() {
    for (FaultSkin skin:_sks) {
      FaultCell[] fcs = FaultSkin.getCells(new FaultSkin[]{skin});
      int nc = fcs.length;
      float[] fls = new float[nc];
      for (int ic=0; ic<nc; ++ic) {
        fls[ic] = fcs[ic].getFl();
      }
      div(fls,max(fls),fls);
      for (int ic=0; ic<nc; ++ic) {
        fcs[ic].setFl(fls[ic]);
      }
    }
  }

  private void setCells(int n1, int n2, int n3) {
    FaultCell[][][] cells = new FaultCell[n3][n2][n1];
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
        int i1i = fc.getI1();
        int i2i = fc.getI2();
        int i3i = fc.getI3();
        int i2m1 = bound(fc.getM2(),n2);
        int i3m1 = bound(fc.getM3(),n3);
        int i2p1 = bound(fc.getP2(),n2);
        int i3p1 = bound(fc.getP3(),n3);

        int d2m = i2m1-i2i;
        int d3m = i3m1-i3i;
        int d2p = i2p1-i2i;
        int d3p = i3p1-i3i;
        int i2m2 = i2m1, i2p2=i2p1; 
        int i3m2 = i3m1, i3p2=i3p1; 
        if(abs(d2m)>0){i2p2-=d2m;i2m2+=d2m;}
        else          {i2m2-=d2p;i2p2+=d2p;}
        if(abs(d3m)>0){i3p2-=d3m;i3m2+=d3m;}
        else          {i3m2-=d3p;i3p2+=d3p;}
        i2m2 = bound(i2m2,n2); i3m2 = bound(i3m2,n3);
        i2p2 = bound(i2p2,n2); i3p2 = bound(i3p2,n3);
        FaultCell cm1 = cells[i3m1][i2m1][i1i];
        FaultCell cp1 = cells[i3p1][i2p1][i1i];
        FaultCell cm2 = cells[i3m2][i2m2][i1i];
        FaultCell cp2 = cells[i3p2][i2p2][i1i];
        FaultCell[] cls = new FaultCell[]{cm1,cp1,cm2,cp2};
        int[] i2s = new int[]{i2m1,i2p1,i2m2,i2p2};
        int[] i3s = new int[]{i3m1,i3p1,i3m2,i3p2};
        int ik = -1;
        for (FaultCell cell:cls) {
          ik ++;
          int i1 = i1i;
          int i2 = i2s[ik];
          int i3 = i3s[ik];
          if(cell==null) {
            cells[i3][i2][i1] = fc;
          } else {
            float s1i = abs(fc.getT1());
            float s1c = abs(cell.getT1());
            if(s1i>s1c) {
              cell.setFl(0.0f);
              cells[i3][i2][i1] = fc;
            }else{fc.setFl(0.0f);}
          }
        }
      }
    }
  }


  private void setCellsC(int n1, int n2, int n3) {
    FaultCell[][][] cells = new FaultCell[n3][n2][n1];
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
        int i1i = fc.getI1();
        int i2i = fc.getI2();
        int i3i = fc.getI3();
        int i2m = bound(fc.getM2(),n2);
        int i3m = bound(fc.getM3(),n3);
        int i2p = bound(fc.getP2(),n2);
        int i3p = bound(fc.getP3(),n3);

        int d2m = i2m-i2i;
        int d3m = i3m-i3i;
        int d2p = i2p-i2i;
        int d3p = i3p-i3i;
        if(abs(d2m)>0){i2p-=d2m;}
        else          {i2m-=d2p;}
        if(abs(d3m)>0){i3p-=d3m;}
        else          {i3m-=d3p;}
        i2m = bound(i2m,n2); i3m = bound(i3m,n3);
        i2p = bound(i2p,n2); i3p = bound(i3p,n3);
        FaultCell cm = cells[i3m][i2m][i1i];
        FaultCell cp = cells[i3p][i2p][i1i];
        FaultCell[] cls = new FaultCell[]{cm,cp};
        int[] i2s = new int[]{i2m,i2p};
        int[] i3s = new int[]{i3m,i3p};
        int ik = -1;
        for (FaultCell cell:cls) {
          ik ++;
          int i1 = i1i;
          int i2 = i2s[ik];
          int i3 = i3s[ik];
          if(cell==null) {
            cells[i3][i2][i1] = fc;
          } else {
            float s1i = abs(fc.getT1());
            float s1c = abs(cell.getT1());
            if(s1i>s1c) {
              cell.setFl(0.0f);
              cells[i3][i2][i1] = fc;
            }else{fc.setFl(0.0f);}
          }
        }
      }
    }
  }

  private FaultSkin[] _sks;
}

