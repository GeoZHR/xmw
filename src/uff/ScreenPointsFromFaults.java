/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package uff;

import ifs.*;
import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Set weights and constraints for flattening: 
 * set zero weights at fault, 
 * set hard constraints using fault slip vectors.
 * @author Xinming Wu
 * @version 2014.09.18
 */

public class ScreenPointsFromFaults {

  public float[][][] getScreenPoints(FaultSkin[] fss, float[][][] w) {
    _n3 = w.length;
    _n2 = w[0].length;
    _n1 = w[0][0].length;
    FaultCell[] fcs = FaultSkin.getCells(fss);
    int[][][] fm = new int[_n3][_n2][_n1];
    faultMap(fcs,fm,w);
    setWeightsOnFault(fss,w);
    FaultCellGrid fcg = new FaultCellGrid(fcs);
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    int nc = fcs.length;
    for (int ic=0; ic<nc; ++ic) {
      float[] cx = fcs[ic].getX();
      float[] cs = fcs[ic].getS();
      float[] px = add(cx,cs);
      if(round(px[0])>=_n1){continue;}
      if(nearestFaultCell(px,fm)) {
        FaultCell fc = fcg.get((int)px[0],(int)px[1],(int)px[2]);
        int[] xm = fc.getIm();
        int[] xp = fc.getIp();
        //cl.add(new float[][]{xm,xp,cs});
      }
    }
    int ns = cl.size();
    System.out.println("sets of control points:"+ns);
    float[][][] cs = new float[3][3][ns];
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
    }
    return cs;

  }


  private boolean nearestFaultCell(float[] xi, int[][][] fm) {
    int d = 4;
    float ds = 500;
    int i1 = round(xi[0]);
    int i2 = round(xi[1]);
    int i3 = round(xi[2]);
    for (int d3=-d; d3<=d; ++d3) {
      for (int d2=-d; d2<=d; ++d2) {
        int x2 = bound2(i2+d2);
        int x3 = bound3(i3+d3);
        if(fm[x3][x2][i1]==1) {
          float dsi = d2*d2+d3*d3;
          if(dsi<ds) {
            ds = dsi;
            xi[0] = i1;
            xi[1] = x2;
            xi[2] = x3;
          }
        }
      }
    }
    if(ds<500.0f) {return true;}
    else          {return false;}
  }

  private void setWeightsOnFault(FaultSkin[] fss, float[][][] ws) {
    int n3 = ws.length;
    int n2 = ws[0].length;
    int m2 = n2-1;
    int m3 = n3-1;
    for (FaultSkin fs:fss) {
      for (FaultCell fc:fs) {
        float[] pi = fc.getX();
        int pi1 = round(pi[0]);
        int pi2 = round(pi[1]);
        int pi3 = round(pi[2]);
        int p2m = pi2-1;if(p2m<0) {p2m=0;}
        int p2p = pi2+1;if(p2p>m2){p2p=m2;}
        int p3m = pi3-1;if(p3m<0) {p3m=0;}
        int p3p = pi3+1;if(p3p>m3){p3p=m3;}
        ws[pi3][pi2][pi1] = 0.0f;
        ws[p3m][pi2][pi1] = 0.0f;
        ws[pi3][p2m][pi1] = 0.0f;
        ws[p3p][pi2][pi1] = 0.0f;
        ws[pi3][p2p][pi1] = 0.0f;
      }
    }
  }


  private int bound1(int x) {
    if(x<0)   {return 0;}
    if(x>=_n1){return _n1-1;}
    return x;
  }
  private int bound2(int x) {
    if(x<0)   {return 0;}
    if(x>=_n2){return _n2-1;}
    return x;
  }
  private int bound3(int x) {
    if(x<0)   {return 0;}
    if(x>=_n3){return _n3-1;}
    return x;
  }

  public void faultMap(FaultCell[] fc, int[][][] fm, float[][][] w) {
    int nc = fc.length;
    for (int ic=0; ic<nc; ++ic) {
      int[] id = fc[ic].getI();
      int i1 = id[0];
      int i2 = id[1];
      int i3 = id[2];
      fm[i3][i2][i1] = 1;
       w[i3][i2][i1] = 0.0f;
    }
  }

  private void addPoints(float[] c, float[] k, float[][][] cp) {
    int c1 = round(c[0]);
    int c2 = round(c[1]);
    int c3 = round(c[2]);
    int k1 = round(k[0]);
    int k2 = round(k[1]);
    int k3 = round(k[2]);
    cp[c3][c2][c1] = 1.0f;
    cp[k3][k2][k1] = 1.0f;
  }


  private void onFault(float[] p, float[][][] w) {
    int i1 = round(p[0]);
    int i2 = round(p[1]);
    int i3 = round(p[2]);
    float wi = w[i3][i2][i1];
    if (wi==0.0f){w[i3][i2][i1]=0.1f;} 
  }

  /*
  private boolean onFault(float[] p, float[][][] w) {
    int i1 = round(p[0]);
    int i2 = round(p[1]);
    int i3 = round(p[2]);
    float wi = w[i3][i2][i1];
    if (wi==0.0f){return true;} 
    else {return false;}
  }

 */
  private boolean onBound(int p1, int p2, int p3) {
    if(p1<0||p1>=_n1){return true;}
    if(p2<0||p2>=_n2){return true;}
    if(p3<0||p3>=_n3){return true;}
    return false;
  }


  private int _n1,_n2,_n3;
}

