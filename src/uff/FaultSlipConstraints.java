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
 * @version 2014.02.01
 */

public class FaultSlipConstraints {

  public FaultSlipConstraints(FaultSkin[] sks) {
    _sks = sks;
  }

  public int[][][] getFaultMap() {
    int fn = _sks.length;
    int[][][] fm = new int[3][fn][]; 
    for (int fi=0; fi<fn; ++fi) { 
      FaultCell[] fcs = _sks[fi].getCells();
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


  public void setNormals(float[][][][] p) {
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
        int[] ip = fc.getIp();
        int[] im = fc.getIm();
        float[] cs = fc.getS();
        int i1p = ip[0];
        int i2p = ip[1];
        int i3p = ip[2];
        int i1m = im[0];
        int i2m = im[1];
        int i3m = im[2];
        float s1 = cs[0];
        float s2 = cs[1];
        float s3 = cs[2];
        float ss = 1.0f/sqrt(s1*s1+s2*s2+s3*s3);
        if(s1<0.0f) {ss = -ss;}
        p[0][i3p][i2p][i1p] = s1*ss;
        p[1][i3p][i2p][i1p] = s2*ss;
        p[2][i3p][i2p][i1p] = s3*ss;
        p[0][i3m][i2m][i1m] = s1*ss;
        p[1][i3m][i2m][i1m] = s2*ss;
        p[2][i3m][i2m][i1m] = s3*ss;
      }
    }
  }

  /*
  public float[][][] getWeightsAndConstraints(
    float[][][] ws, float[][][] fm, float[][][] cp) {
    int ct = 2;
    int n3 = ws.length;
    int n2 = ws[0].length;
    int n1 = ws[0][0].length;
    short[][][] mk = new short[n3][n2][n1];
    float[][][] ds = setWeightsNearFaults(ws,fm);
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin fs:_sks) {
      ct++;
    for (FaultCell fc:fs) {
      float[] cx = fc.getX();
      float[] cw = fc.getW();
      float[] cs = fc.getS();
      float[] hx = add(cx,cs);
      float[] fx = add(cx,cs);
      boolean sh = shiftToHangWall(hx,cw,ds);
      boolean sf = shiftToFootWall(fx,cw,ds);
      if(sh && sf) {
        int h1 = round(hx[0]);
        int h2 = round(hx[1]);
        int h3 = round(hx[2]);
        int f1 = round(fx[0]);
        int f2 = round(fx[1]);
        int f3 = round(fx[2]);
        int mh = mk[h3][h2][h1];
        int mf = mk[f3][f2][f1];
        if(mh==1||mf==1){continue;}
        mk[h3][h2][h1] = 1;
        mk[f3][f2][f1] = 1;
        cp[h3][h2][h1] = ct;
        cp[f3][f2][f1] = ct;
        hx[0] = h1; fx[0] = f1;
        hx[1] = h2; fx[1] = f2;
        hx[2] = h3; fx[2] = f3;
        cl.add(new float[][]{fx,hx,mul(cs,0.5f)});
      }
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

  */
  public float[][][] screenPoints(
    float[][][] ws, float[][][] wp, float[][][] cp) 
  {
    int ct = 2;
    int nk = _sks.length;
    setWeightsNearFaults(ws,wp);
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (int ik=0; ik<nk; ++ik) {
      ct++;
      FaultSkin[] ski = new FaultSkin[]{_sks[ik]};
      FaultCell[] fcs = FaultSkin.getCells(ski);
      int nc = fcs.length;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fc = fcs[ic];
        float fl = fc.getFl();
        int[] im = fc.getIm();
        int[] ip = fc.getIp();
        float[] cs = fc.getS();
        float[] hx = new float[3];
        float[] fx = new float[3];
        float[] xf = new float[3];
        int f1 = im[0];
        int f2 = im[1];
        int f3 = im[2];
        int h1 = ip[0];
        int h2 = ip[1];
        int h3 = ip[2];
        cp[h3][h2][h1] += ct;
        cp[f3][f2][f1] += ct;
        hx[0] = h1; fx[0] = f1;
        hx[1] = h2; fx[1] = f2;
        hx[2] = h3; fx[2] = f3;
        xf[0] = fl; xf[1] = fl; xf[2] = fl;
        cl.add(new float[][]{hx,fx,cs,xf});
        ws[f3][f2][f1] = 0.0f;
        ws[h3][h2][h1] = 0.0f;
        wp[f3][f2][f1] = 0.0f;
        wp[h3][h2][h1] = 0.0f;
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

  private boolean shiftToHangWall(float[] hx, float[] cw, float[][][] ds) {
    float w2 = cw[1];
    float w3 = cw[2];
    float x1, x2, x3;
    int n3 = ds.length;
    int n2 = ds[0].length;
    int n1 = ds[0][0].length;
    for (float di=0.5f; di<3.0f; di+=0.5f) {
      x1 = hx[0];
      x2 = hx[1]+di*w2;
      x3 = hx[2]+di*w3;
      int i1 = round(x1);
      int i2 = round(x2);
      int i3 = round(x3);
      if(i1<0||i1>=n1){return false;}
      if(i2<0||i2>=n2){return false;}
      if(i3<0||i3>=n3){return false;}
      if(ds[i3][i2][i1]>=1.0f) {
        hx[0] = x1;
        hx[1] = x2;
        hx[2] = x3;
        return true;
      }
    }
    return false;
  }

  private boolean shiftToFootWall(float[] fx, float[] cw, float[][][] ds) {
    float w2 = cw[1];
    float w3 = cw[2];
    float x1, x2, x3;
    int n3 = ds.length;
    int n2 = ds[0].length;
    int n1 = ds[0][0].length;
    for (float di=0.5f; di<3.0f; di+=0.5f) {
      x1 = fx[0];
      x2 = fx[1]-di*w2;
      x3 = fx[2]-di*w3;
      int i1 = round(x1);
      int i2 = round(x2);
      int i3 = round(x3);
      if(i1<0||i1>=n1){return false;}
      if(i2<0||i2>=n2){return false;}
      if(i3<0||i3>=n3){return false;}
      if(ds[i3][i2][i1]>=1.0f) {
        fx[0] = x1;
        fx[1] = x2;
        fx[2] = x3;
        return true;
      }
    }
    return false;
  }


  public void setWeightsNearFaults(float[][][] ws, float[][][] wp) {
    float fnull = -99f;
    int n3 = ws.length;
    int n2 = ws[0].length;
    int n1 = ws[0][0].length;
    float[][][] ds = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    float[][][] fm = fillfloat(fnull,n1,n2,n3);
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
        int[] is = fc.getI();
        int i1 = is[0];
        int i2 = is[1];
        int i3 = is[2];
        fm[i3][i2][i1] = 0.0f;
      }
    }
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(fnull,fm,ds,k1,k2,k3);
    for (int i3=0; i3<n3; ++i3){
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0; i1<n1; ++i1){
      float di = ds[i3][i2][i1];
      float wi = wp[i3][i2][i1];
      if(di<=2f) { 
        ws[i3][i2][i1] = 1.0f;
        if(wi>0.6f) {wp[i3][i2][i1] = 0.6f;}
      }
    }}}
  }

  FaultSkin[] _sks;
}

