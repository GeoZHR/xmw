/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.util.*;
import edu.mines.jtk.util.*;
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
    for (FaultSkin sk:_sks) {
      for (FaultCell cell:sk) {
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
    computeUnfaultShifts(n1,n2,n3,_sks);
    setCellsC(n1,n2,n3);
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
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
    computeUnfaultShifts(_sks);
    //computeUnfaultShifts(n1,n2,n3,_sks);
    setCells(n1,n2,n3);
    flNormalization();
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin sk:_sks) {
      for (FaultCell fc:sk) {
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

        float t1 = fc.getT1();
        float t2 = fc.getT2();
        float t3 = fc.getT3();
        float[] ch = new float[3];
        float[] cf = new float[3];
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


  private void computeUnfaultShifts(FaultSkin[] skins) {
    int ik = 0;
    for (FaultSkin skin:skins) {
      System.out.println("ik="+ik);
      ik++;
      final FaultCell[] cells = skin.getCells();
      final int nc = cells.length;
      final float[][] sc = new float[3][nc];
      final float[][] xc = new float[3][nc];
      for (int ic=0; ic<nc; ++ic) {
        sc[0][ic] = cells[ic].getS1();
        sc[1][ic] = cells[ic].getS2();
        sc[2][ic] = cells[ic].getS3();
        xc[0][ic] = cells[ic].getX1();
        xc[1][ic] = cells[ic].getX2();
        xc[2][ic] = cells[ic].getX3();
      }
      final KdTree kt = new KdTree(xc);
      Parallel.loop(nc,new Parallel.LoopInt() {
      public void compute(int ic) {
      //for (int ic=0; ic<nc; ++ic) {
        FaultCell cell = cells[ic];
        float[] xminp = new float[3];
        float[] xmaxp = new float[3];
        float[] xminm = new float[3];
        float[] xmaxm = new float[3];
        float s1 = sc[0][ic];
        float s2 = sc[1][ic];
        float s3 = sc[2][ic];
        float x1 = xc[0][ic];
        float x2 = xc[1][ic];
        float x3 = xc[2][ic];
        float p1 = x1+s1;
        float p2 = x2+s2;
        float p3 = x3+s3;
        float m1 = x1-s1;
        float m2 = x2-s2;
        float m3 = x3-s3;
        xminp[0] = p1-6;
        xminp[1] = p2-6;
        xminp[2] = p3-6;
        xmaxp[0] = p1+6;
        xmaxp[1] = p2+6;
        xmaxp[2] = p3+6;
        xminm[0] = m1-6;
        xminm[1] = m2-6;
        xminm[2] = m3-6;
        xmaxm[0] = m1+6;
        xmaxm[1] = m2+6;
        xmaxm[2] = m3+6;
        int[] idp = kt.findInRange(xminp,xmaxp);
        int[] idm = kt.findInRange(xminm,xmaxm);
        int ndp = idp.length;
        int ndm = idm.length;
        if (ndp>1&&ndm>1){
          float[] s1p = new float[ndp];
          float[] x1p = new float[ndp];
          float[] x2p = new float[ndp];
          float[] x3p = new float[ndp];
          float[] s1m = new float[ndm];
          float[] x1m = new float[ndm];
          float[] x2m = new float[ndm];
          float[] x3m = new float[ndm];
          for (int ip=0; ip<ndp; ++ip) {
            int pc = idp[ip];
            s1p[ip] = sc[0][pc];
            x1p[ip] = xc[0][pc];
            x2p[ip] = xc[1][pc];
            x3p[ip] = xc[2][pc];
          }
          SibsonInterp sip = new SibsonInterp(s1p,x1p,x2p,x3p);
          sip.setBounds(min(x1p)-5,max(x1p)+5,min(x2p)-5,max(x2p)+5,
                        min(x3p)-5,max(x3p)+5);
          float dp = sip.interpolate(p1,p2,p3)-s1;
          for (int im=0; im<ndm; ++im) {
            int mc = idm[im];
            s1m[im] = sc[0][mc];
            x1m[im] = xc[0][mc];
            x2m[im] = xc[1][mc];
            x3m[im] = xc[2][mc];
          }
          SibsonInterp sim = new SibsonInterp(s1m,x1m,x2m,x3m);
          sim.setBounds(min(x1m)-5,max(x1m)+5,min(x2m)-5,max(x2m)+5,
                        min(x3m)-5,max(x3m)+5);
          float dm = s1-sim.interpolate(m1,m2,m3);
          s1 -= (dp+dm)*0.5f;
        }
        cell.setUnfaultShifts(new float[]{s1,s2,s3});
      }});

    }
  }
  private void computeUnfaultShifts(
    final int n1, final int n2, final int n3, final FaultSkin[] skins) {
    final int nk = skins.length;
    for (int ik=0; ik<nk; ++ik) {
    //Parallel.loop(nk,new Parallel.LoopInt() {
    //public void compute(int ik) {
      System.out.println("skin="+ik);
      final FaultCell[] cells = skins[ik].getCells();
      final int nc = cells.length;
      FloatList x1l = new FloatList();
      FloatList x2l = new FloatList();
      FloatList x3l = new FloatList();
      FloatList s1l = new FloatList();
      for (int ic=0; ic<nc; ++ic) {
        FaultCell cell = cells[ic];
        x1l.add(cell.getX1());
        x2l.add(cell.getX2());
        x3l.add(cell.getX3());
        s1l.add(cell.getS1());
      }
      float[] x1a = x1l.trim();
      float[] x2a = x2l.trim();
      float[] x3a = x3l.trim();
      float[] s1a = s1l.trim();
      float x1min = 0;
      float x2min = 0;
      float x3min = 0;
      float x1max = n1-1;
      float x2max = n2-1;
      float x3max = n3-1;
      SibsonInterp s1i = new SibsonInterp(s1a,x1a,x2a,x3a);
      s1i.setBounds(x1min,x1max,x2min,x2max,x3min,x3max);
      for (int ic=0; ic<nc; ++ic) {
        //System.out.println("ic="+ic);
        FaultCell cell = cells[ic];
        float s1 = cell.getS1();
        float s2 = cell.getS2();
        float s3 = cell.getS3();
        float p1 = cell.getX1()+s1;
        float p2 = cell.getX2()+s2;
        float p3 = cell.getX3()+s3;
        float m1 = cell.getX1()-s1;
        float m2 = cell.getX2()-s2;
        float m3 = cell.getX3()-s3;

        if(p1>x1max){p1=x1max;}
        if(p2>x2max){p2=x2max;}
        if(p3>x3max){p3=x3max;}
        if(p1<x1min){p1=x1min;}
        if(p2<x2min){p2=x2min;}
        if(p3<x3min){p3=x3min;}
        if(m1>x1max){m1=x1max;}
        if(m2>x2max){m2=x2max;}
        if(m3>x3max){m3=x3max;}
        if(m1<x1min){m1=x1min;}
        if(m2<x2min){m2=x2min;}
        if(m3<x3min){m3=x3min;}

        float dp = s1i.interpolate(p1,p2,p3)-s1;
        float dm = s1-s1i.interpolate(m1,m2,m3);
        s1 -= (dp+dm)*0.5f;
        cell.setUnfaultShifts(new float[]{s1,s2,s3});
      }
    }
    //}});
    //checkUnfaultShifts();
  }

  private void checkUnfaultShifts() {
    int nk = _sks.length;
    for (int ik=0; ik<nk; ++ik) {
      FaultSkin skin = _sks[ik];
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
    for (FaultSkin skin:_sks){
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
        int mc = cls.length;
        for (int im=0; im<mc; ++im) {
          FaultCell cell=cls[im];
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


  private class FloatList {
    public int n = 0;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public void add(float[] f) {
      int m = f.length;
      for (int i=0; i<m; ++i)
        add(f[i]);
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

  private FaultSkin[] _sks;
}

