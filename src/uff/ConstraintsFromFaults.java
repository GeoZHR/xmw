/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package uff;

import ipf.*;
import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Set weights and constraints for flattening: 
 * set zero weights at fault, 
 * set hard constraints using fault slip vectors.
 * @author Xinming Wu
 * @version 2014.09.18
 */

public class ConstraintsFromFaults {

  public ConstraintsFromFaults(
    FaultSkin[] fss, float[][][] w) {
    _w = w;
    _fss = fss;
    _n3 = w.length;
    _n2 = w[0].length;
    _n1 = w[0][0].length;
    _mk = new int[_n3][_n2][_n1];
    _fm = new int[_n3][_n2][_n1];
    faultMap(_fm);
  }

  public float[][][] getWeightsAndConstraints(float[][][] ws, float[][][] cp) {
    setWeightsOnFault(ws);
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin fs:_fss) {
      for (FaultCell fc:fs) {
        float[] cx = fc.getX();
        float[] cs = fc.getS();
        float[] cw = fc.getW();
        float[] fx = new float[3];
        float[] hx = new float[3];
        hx[0] = bound1(round(cx[0]));//+cs[0]));
        hx[1] = bound2(round(cx[1]));//+cs[1]));
        hx[2] = bound3(round(cx[2]));//+cs[2]));
        //if(!nearestFaultCell(hx)) {continue;}
        fx = copy(hx);
        boolean valid = false;
        float w2 = abs(cw[1]);
        float w3 = abs(cw[2]);
        if (w2>w3) {valid = shift2(cw[1],fx,hx);} 
        else       {valid = shift3(cw[2],fx,hx);}
        if(valid) {
          onFault(fx,ws);
          onFault(hx,ws);
          //if (onFault(fx,ws)) {continue;}
          //if (onFault(hx,ws)) {continue;}
          cl.add(new float[][]{fx,hx,mul(cs,0.5f)});
          addPoints(fx,hx,cp);
        }
      }
    }
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

  private boolean badQuality(float[] cx, float[] cs) {
    float qx = 0.0f;
    float qy = 0.0f;
    float sc = 1.f/4.f;
    int x1 = round(cx[0]);
    int x2 = round(cx[1]);
    int x3 = round(cx[2]);
    int x2m = bound2(x2-2);
    int x2p = bound2(x2+2);
    int x3m = bound2(x3-2);
    int x3p = bound2(x3+2);
    qx += _w[x3m][x2 ][x1];
    qx += _w[x3p][x2 ][x1];
    qx += _w[x3 ][x2m][x1];
    qx += _w[x3 ][x2p][x1];
    if(qx*sc<.3f){return true;}
    float[] cy = add(cx,cs);
    cy[0] = bound1(round(cy[0]));
    if(nearestFaultCell(cy)) {
      int y1 = round(cy[0]);
      int y2 = round(cy[1]);
      int y3 = round(cy[2]);
      int y2m = bound2(y2-2);
      int y2p = bound2(y2+2);
      int y3m = bound2(y3-2);
      int y3p = bound2(y3+2);
      qy += _w[y3m][y2 ][y1];
      qy += _w[y3p][y2 ][y1];
      qy += _w[y3 ][y2m][y1];
      qy += _w[y3 ][y2p][y1];
      if(qy*sc<.3f){return true;}
    }
    return false;
  }

  private boolean nearestFaultCell(float[] xi) {
    int d = 4;
    float ds = 500;
    int i1 = round(xi[0]);
    int i2 = round(xi[1]);
    int i3 = round(xi[2]);
    for (int d3=-d; d3<=d; ++d3) {
      for (int d2=-d; d2<=d; ++d2) {
        int x2 = bound2(i2+d2);
        int x3 = bound3(i3+d3);
        if(_fm[x3][x2][i1]==1) {
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

  public void faultMap(int[][][] fm) {
    for (FaultSkin fs:_fss) {
      for (FaultCell fc:fs) {
        float[] x = fc.getX();
        int x1 = round(x[0]);
        int x2 = round(x[1]);
        int x3 = round(x[2]);
        fm[x3][x2][x1] = 1;
      }
    }
  }

  public int[][][] getFaultMap() {
    int fn = _fss.length;
    int[][][] fm = new int[3][fn][]; 
    for (int fi=0; fi<fn; ++fi) { 
      FaultCell[] fcs = _fss[fi].getCells();
      int nc = fcs.length;
      fm[0][fi] = new int[nc];
      fm[1][fi] = new int[nc];
      fm[2][fi] = new int[nc];
      for (int ic=0; ic<nc; ++ic) { 
        float[] xc = fcs[ic].getX();
        fm[0][fi][ic] = round(xc[0]);
        fm[1][fi][ic] = round(xc[1]);
        fm[2][fi][ic] = round(xc[2]);
      }
    }
    return fm;
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

  private boolean shift2(float w2, float[] c, float[] k) {
    float sn2 = (w2<0.f)?-1.f:1.f;
    float ds2 = sn2*2.0f;
    float ds1 = sn2*1.0f;

    c[1] -= ds2;
    c[2] -= ds1;
    int c1 = round(c[0]);
    int c3 = round(c[2]);
    int c2 = round(c[1]);
    if(onBound(c1,c2,c3)){return false;}

    k[1] +=ds2;
    k[2] +=ds1;
    int k1 = round(k[0]);
    int k2 = round(k[1]);
    int k3 = round(k[2]);
    if(onBound(k1,k2,k3)){return false;}

    _mk[c3][c2][c1] += 1;
    if(_mk[c3][c2][c1]>1) {return false;}

    _mk[k3][k2][k1] += 1;
    if(_mk[k3][k2][k1]>1) {return false;}
 
    return true;
  }

  private boolean shift3(float w3, float[] c, float[] k) {
    float sn3 = (w3<0.f)?-1.f:1.f;
    float ds3 = sn3*2.0f;
    float ds1 = sn3*1.0f;

    c[1] -= ds1;
    c[2] -= ds3;
    int c1 = round(c[0]);
    int c2 = round(c[1]);
    int c3 = round(c[2]);
    if(onBound(c1,c2,c3)){return false;}

    k[1] += ds1;
    k[2] += ds3;
    int k1 = round(k[0]);
    int k2 = round(k[1]);
    int k3 = round(k[2]);
    if(onBound(k1,k2,k3)){return false;}

    _mk[c3][c2][c1] += 1;
    if(_mk[c3][c2][c1]>1) {return false;}
    _mk[k3][k2][k1] += 1;
    if(_mk[k3][k2][k1]>1) {return false;}

    return true;
  }


  private void setWeightsOnFault(float[][][] ws) {
    int n3 = ws.length;
    int n2 = ws[0].length;
    int m2 = n2-1;
    int m3 = n3-1;
    for (FaultSkin fs:_fss) {
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
  private FaultSkin[] _fss;
  private float[][][] _w = null;
  private float[][][] _p2 = null;
  private float[][][] _p3 = null;
  private int[][][] _mk = null;
  private int[][][] _fm = null;
}

