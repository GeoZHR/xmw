/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package uff;

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

  public float[][][] controlPoints(
    float[][][] ws, float[][][] wp, float[][][] cp) {
    int ct = 2;
    int n3 = ws.length;
    int n2 = ws[0].length;
    int n1 = ws[0][0].length;
    setCellsC(n1,n2,n3);
    short[][][] mk = new short[n3][n2][n1];
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin fs:_sks) {
      ct++;
    for (FaultCell fc:fs) {
      int[] is = fc.getI();
      int[] im = fc.getIm();
      int[] ip = fc.getIp();
      float[] cs = fc.getS();
      float[] hx = new float[3];
      float[] fx = new float[3];
      int i1i = is[0];
      int i2i = is[1];
      int i3i = is[2];
      int i2p = ip[1];
      int i3p = ip[2];
      int i2m = im[1];
      int i3m = im[2];
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
      int mh = mk[i3m][i2m][i1i];
      int mf = mk[i3p][i2p][i1i];
      ws[i3i][i2i][i1i] = 0.000f;
      wp[i3i][i2i][i1i] = 0.005f;
      wp[i3m][i2m][i1i] = 0.005f;
      wp[i3p][i2p][i1i] = 0.005f;
      if(mh==1||mf==1){continue;}
      if(fc.getFl()==0.0f){continue;}
      mk[i3m][i2m][i1i] = 1;
      mk[i3p][i2p][i1i] = 1;
      cp[i3m][i2m][i1i] = ct;
      cp[i3p][i2p][i1i] = ct;
      hx[0] = i1i; fx[0] = i1i;
      hx[1] = i2p; fx[1] = i2m;
      hx[2] = i3p; fx[2] = i3m;
      cs[1] *= 2f; cs[2] *= 2f;
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

  public float[][][] screenPointsX(
    float[][][] ws, float[][][][] p){
    int nk = _sks.length;
    int n3 = p[0].length;
    int n2 = p[0][0].length;
    int n1 = p[0][0][0].length;
    setCells(n1,n2,n3);
    float[][] fln = normalizeFl();
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (int ik=0; ik<nk; ++ik) {
      FaultSkin[] ski = new FaultSkin[]{_sks[ik]};
      FaultCell[] fcs = FaultSkin.getCells(ski);
      int nc = fcs.length;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fc = fcs[ic];
        int[] is = fc.getI();
        int[] im = fc.getIm();
        int[] ip = fc.getIp();
        float[] cs = fc.getS();
        float[] hx = new float[3];
        float[] fx = new float[3];
        float[] sc = new float[3];
        int i1i = is[0];
        int i2i = is[1];
        int i3i = is[2];
        int i2p = ip[1];
        int i3p = ip[2];
        int i2m = im[1];
        int i3m = im[2];
        int d2m = i2m-i2i;
        int d3m = i3m-i3i;
        int d2p = i2p-i2i;
        int d3p = i3p-i3i;
        if(abs(d2m)>0){i2p -= d2m;}
        else          {i2m -= d2p;}
        if(abs(d3m)>0){i3p -= d3m;}
        else          {i3m -= d3p;}
        if(i2m<0){i2m=0;}if(i2m>=n2){i2m=n2-1;}
        if(i2p<0){i2p=0;}if(i2p>=n2){i2p=n2-1;}
        if(i3m<0){i3m=0;}if(i3m>=n3){i3m=n3-1;}
        if(i3p<0){i3p=0;}if(i3p>=n3){i3p=n3-1;}
        hx[0] = i1i; fx[0] = i1i;
        hx[1] = i2p; fx[1] = i2m;
        hx[2] = i3p; fx[2] = i3m;
        ws[i3i][i2i][i1i] = 0.00f;
        ws[i3m][i2m][i1i] = 0.01f;
        ws[i3p][i2p][i1i] = 0.01f;
        p[3][i3i][i2i][i1i] = 0.f;
        p[3][i3m][i2m][i1i] = 0.f;
        p[3][i3p][i2p][i1i] = 0.f;
        float fl = fln[ik][ic];
        sc[0] = 1f; sc[1] = 1f; sc[2] = 1f;
        //sc[0] = fl; sc[1] = fl; sc[2] = fl;
        cs[1] *= 2f; cs[2] *= 2f;
        cl.add(new float[][]{hx,fx,cs,sc});
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

  public float[][][] screenPoints(
    float[][][] ws, float[][][][] p){
    int nk = _sks.length;
    int n3 = p[0].length;
    int n2 = p[0][0].length;
    int n1 = p[0][0][0].length;
    setCellsC(n1,n2,n3);
    float[][] fln = normalizeFl();
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (int ik=0; ik<nk; ++ik) {
      FaultSkin[] ski = new FaultSkin[]{_sks[ik]};
      FaultCell[] fcs = FaultSkin.getCells(ski);
      int nc = fcs.length;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fc = fcs[ic];
        int[] is = fc.getI();
        int[] im = fc.getIm();
        int[] ip = fc.getIp();
        float[] cs = fc.getS();
        float[] h1 = new float[3];
        float[] f1 = new float[3];
        float[] c1 = new float[3];
        float[] h2 = new float[3];
        float[] f2 = new float[3];
        float[] c2 = new float[3];
        int i1i = is[0];
        int i2i = is[1];
        int i3i = is[2];
        int i2m1 = im[1];
        int i3m1 = im[2];
        int i2p1 = ip[1];
        int i3p1 = ip[2];
        if(i2m1<0)   {i2m1=0;}
        if(i2p1<0)   {i2p1=0;}
        if(i3m1<0)   {i3m1=0;}
        if(i3p1<0)   {i3p1=0;}
        if(i2m1>=n2) {i2m1=n2-1;}
        if(i2p1>=n2) {i2p1=n2-1;}
        if(i3m1>=n3) {i3m1=n3-1;}
        if(i3p1>=n3) {i3p1=n3-1;}
        h1[0] = i1i;  f1[0] = i1i;
        h1[1] = i2p1; f1[1] = i2m1;
        h1[2] = i3p1; f1[2] = i3m1;
        ws[i3m1][i2m1][i1i] = 0.00f;
        ws[i3p1][i2p1][i1i] = 0.00f;
        p[3][i3m1][i2m1][i1i] = 0.f;
        p[3][i3p1][i2p1][i1i] = 0.f;

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
        if(i2m2<0)   {i2m2=0;}
        if(i2p2<0)   {i2p2=0;}
        if(i3m2<0)   {i3m2=0;}
        if(i3p2<0)   {i3p2=0;}
        if(i2m2>=n2) {i2m2=n2-1;}
        if(i2p2>=n2) {i2p2=n2-1;}
        if(i3m2>=n3) {i3m2=n3-1;}
        if(i3p2>=n3) {i3p2=n3-1;}
        h2[0] = i1i;  f2[0] = i1i;
        h2[1] = i2p2; f2[1] = i2m2;
        h2[2] = i3p2; f2[2] = i3m2;
        p[3][i3m2][i2m2][i1i] = 0.f;
        p[3][i3p2][i2p2][i1i] = 0.f;
        float fl = fln[ik][ic];
        //c1[0] = 1f; c2[1] = 1f; c2[2] = 1f;
        c1[0] = fl; c1[1] = fl; c1[2] = fl;
        cs[1] *= 2f; cs[2] *= 2f;
        cl.add(new float[][]{h1,f1,cs,c1});
        cl.add(new float[][]{h2,h1,new float[3],c1});
        cl.add(new float[][]{f2,f1,new float[3],c1});
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

  public void setValuesOnFaults(float v, FaultSkin[] sks, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        int[] ms = fc.getIm();
        int[] ps = fc.getIm();
        int i1m = ms[0];
        int i2m = ms[1];
        int i3m = ms[2];
        int i1p = ps[0];
        int i2p = ps[1];
        int i3p = ps[2];
        if(i2m>=n2){continue;}
        if(i3m>=n3){continue;}
        fl[i3m][i2m][i1m] = v;
        fl[i3p][i2p][i1p] = v;
      }
    }
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

  /*
  private void setCells(int n1, int n2, int n3) {
    FaultCell[][][] cells = new FaultCell[n3][n2][n1];
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
        int[] ic = fc.getI();
        int i1 = ic[0];
        int i2 = ic[1];
        int i3 = ic[2];
        FaultCell cl = cells[i3][i2][i1];
        if(cl==null){cells[i3][i2][i1]=cl;}
        else {
          float fli = fc.getFl();
          float flc = cl.getFl();
          if(fli>flc) {
            cl.setFl(0.0f);
            cells[i3][i2][i1] = fc;
          } else {
            fc.setFl(0.0f);
          }
        }
      }
    }
  }
  */
  private void setCellsC(int n1, int n2, int n3) {
    FaultCell[][][] cells = new FaultCell[n3][n2][n1];
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
        int[] ic = fc.getI();
        int[] im = fc.getIm();
        int[] ip = fc.getIp();
        int i1i = ic[0];
        int i2i = ic[1];
        int i3i = ic[2];
        int i2m = im[1];
        int i3m = im[2];
        int i2p = ip[1];
        int i3p = ip[2];
        if(i2m<0)   {i2m=0;}
        if(i2p<0)   {i2p=0;}
        if(i3m<0)   {i3m=0;}
        if(i3p<0)   {i3p=0;}
        if(i2m>=n2) {i2m=n2-1;}
        if(i2p>=n2) {i2p=n2-1;}
        if(i3m>=n3) {i3m=n3-1;}
        if(i3p>=n3) {i3p=n3-1;}
        int d2m = i2m-i2i;
        int d3m = i3m-i3i;
        int d2p = i2p-i2i;
        int d3p = i3p-i3i;
        if(abs(d2m)>0){i2p -= d2m;}
        else          {i2m -= d2p;}
        if(abs(d3m)>0){i3p -= d3m;}
        else          {i3m -= d3p;}
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
            float s1i = abs(fc.getS1());
            float s1c = abs(cell.getS1());
            if(s1i>s1c) {
              cell.setFl(0.0f);
              cells[i3][i2][i1] = fc;
            }else{fc.setFl(0.0f);}
          }
        }
      }
    }
  }


  private void setCells(int n1, int n2, int n3) {
    FaultCell[][][] cells = new FaultCell[n3][n2][n1];
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
        int[] ic = fc.getI();
        int[] im = fc.getIm();
        int[] ip = fc.getIp();
        int i1i = ic[0];
        int i2i = ic[1];
        int i3i = ic[2];
        int i2m1 = im[1];
        int i3m1 = im[2];
        int i2p1 = ip[1];
        int i3p1 = ip[2];
        if(i2m1<0)   {i2m1=0;}
        if(i2p1<0)   {i2p1=0;}
        if(i3m1<0)   {i3m1=0;}
        if(i3p1<0)   {i3p1=0;}
        if(i2m1>=n2) {i2m1=n2-1;}
        if(i2p1>=n2) {i2p1=n2-1;}
        if(i3m1>=n3) {i3m1=n3-1;}
        if(i3p1>=n3) {i3p1=n3-1;}

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
        if(i2m2<0)   {i2m2=0;}
        if(i2p2<0)   {i2p2=0;}
        if(i3m2<0)   {i3m2=0;}
        if(i3p2<0)   {i3p2=0;}
        if(i2m2>=n2) {i2m2=n2-1;}
        if(i2p2>=n2) {i2p2=n2-1;}
        if(i3m2>=n3) {i3m2=n3-1;}
        if(i3p2>=n3) {i3p2=n3-1;}
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
            float s1i = abs(fc.getS1());
            float s1c = abs(cell.getS1());
            if(s1i>s1c) {
              cell.setFl(0.0f);
              cells[i3][i2][i1] = fc;
            }else{fc.setFl(0.0f);}
          }
        }
      }
    }
  }
  FaultSkin[] _sks;
}

