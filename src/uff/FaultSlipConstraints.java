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
 * @version 2014.02.09
 */

public class FaultSlipConstraints {

  public FaultSlipConstraints(FaultSkin[] sks) {
    _sks = sks;
  }

  public float[][][] controlPoints(
    float[][][] ws, float[][][] wp, float[][][] cp) {
    int ct = 2;
    int n3 = ws.length;
    int n2 = ws[0].length;
    int n1 = ws[0][0].length;
    short[][][] mk = new short[n3][n2][n1];
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin fs:_sks) {
      ct++;
    for (FaultCell fc:fs) {
      int[] im = fc.getIm();
      int[] ip = fc.getIp();
      float[] cs = fc.getS();
      float[] hx = new float[3];
      float[] fx = new float[3];
      int h1 = ip[0];
      int h2 = ip[1];
      int h3 = ip[2];
      int f1 = im[0];
      int f2 = im[1];
      int f3 = im[2];
      int mh = mk[h3][h2][h1];
      int mf = mk[f3][f2][f1];
      ws[h3][h2][h1] = 0.000f;
      ws[f3][f2][f1] = 0.000f;
      wp[h3][h2][h1] = 0.005f;
      wp[f3][f2][f1] = 0.005f;
      if(mh==1||mf==1){continue;}
      mk[h3][h2][h1] = 1;
      mk[f3][f2][f1] = 1;
      cp[h3][h2][h1] = ct;
      cp[f3][f2][f1] = ct;
      hx[0] = h1; fx[0] = f1;
      hx[1] = h2; fx[1] = f2;
      hx[2] = h3; fx[2] = f3;
      cl.add(new float[][]{hx,fx,mul(cs,0.5f)});
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
    int nk = _sks.length;
    float[][] fln = normalizeFl();
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (int ik=0; ik<nk; ++ik) {
      FaultSkin[] ski = new FaultSkin[]{_sks[ik]};
      FaultCell[] fcs = FaultSkin.getCells(ski);
      int nc = fcs.length;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fc = fcs[ic];
        int[] im = fc.getIm();
        int[] ip = fc.getIp();
        float[] cs = fc.getS();
        float fl = fln[ik][ic];
        float[] hx = new float[3];
        float[] fx = new float[3];
        float[] xf = new float[3];
        int f1 = im[0];
        int f2 = im[1];
        int f3 = im[2];
        int h1 = ip[0];
        int h2 = ip[1];
        int h3 = ip[2];
        wp[f3][f2][f1] = 0.0f;
        wp[h3][h2][h1] = 0.0f;
        hx[0] = h1; fx[0] = f1;
        hx[1] = h2; fx[1] = f2;
        hx[2] = h3; fx[2] = f3;
        xf[0] = fl; xf[1] = fl; xf[2] = fl;
        //cs[1] *= 2f; cs[2] *= 2f;
        cl.add(new float[][]{hx,fx,cs,xf});
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

  public float[][][] screenPointsNearFaults(int n2, int n3) {
    int nk = _sks.length;
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (int ik=0; ik<nk; ++ik) {
      FaultSkin[] ski = new FaultSkin[]{_sks[ik]};
      FaultCell[] fcs = FaultSkin.getCells(ski);
      int nc = fcs.length;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fc = fcs[ic];
        int[] im = fc.getIm();
        int[] ip = fc.getIp();
        int i1m = im[0], i2m = im[1], i3m = im[2];
        int i1p = ip[0], i2p = ip[1], i3p = ip[2];
        if(i2p>i2m) {
          i2m -= 1; if(i2m<0)   {i2m=0;}
          i2p += 1; if(i2p>=n2) {i2p=n2-1;}
        } else {
          i2p -= 1; if(i2p<0)   {i2p=0;}
          i2m += 1; if(i2m>=n2) {i2m=n2-1;}
        }
        if(i3p>i3m) {
          i3m -= 1; if(i3m<0)   {i3m=0;}
          i3p += 1; if(i3p>=n3) {i3p=n3-1;}
        } else {
          i3p -= 1; if(i3p<0)   {i3p=0;}
          i3m += 1; if(i3m>=n3) {i3m=n3-1;}
        }
        float[] hx = new float[]{i1p,i2p,i3p};
        float[] fx = new float[]{i1m,i2m,i3m};
        cl.add(new float[][]{hx,fx});
      }
    }
    int ns = cl.size();
    float[][][] cs = new float[2][3][ns];
    for (int is=0; is<ns; ++is) {
      float[][] ps = cl.get(is);
      cs[0][0][is] = ps[0][0];
      cs[0][1][is] = ps[0][1];
      cs[0][2][is] = ps[0][2];

      cs[1][0][is] = ps[1][0];
      cs[1][1][is] = ps[1][1];
      cs[1][2][is] = ps[1][2];
    }
    return cs;
  }


  public void setNormals(float[][][][] p) {
    //setNormalsNearFaults(p);
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

  private void setNormalsNearFaults(float[][][][] ps) {
    float fnull = -99f;
    int n3 = ps[0].length;
    int n2 = ps[0][0].length;
    int n1 = ps[0][0][0].length;
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
      if(di<=1f) { 
        ps[0][i3][i2][i1] = 1.0f;
        ps[1][i3][i2][i1] = 0.0f;
        ps[2][i3][i2][i1] = 0.0f;
      }
    }}}
  }


  private float[][] normalizeFl() {
    int nk = _sks.length;
    float[][] fl = new float[nk][];
    for (int ik=0; ik<nk; ++ik) {
      FaultSkin[] ski = new FaultSkin[]{_sks[ik]};
      FaultCell[] fcs = FaultSkin.getCells(ski);
      int nc = fcs.length;
      fl[ik] = new float[nc];
      float fmin = fcs[0].getFl();
      float fmax = fcs[0].getFl();
      for (int ic=0; ic<nc; ++ic) {
        float flc = fcs[ic].getFl();
        fl[ik][ic] = flc;
        if(flc<fmin) {fmin=flc;}
        if(flc>fmax) {fmax=flc;}
      }
      //sub(fl[ik],fmin,fl[ik]);
      div(fl[ik],fmax,fl[ik]);
    }
    return fl;
  }

  FaultSkin[] _sks;
}

