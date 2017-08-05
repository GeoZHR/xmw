/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.*;
import pik.*;

/**
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class SaltPicker3 {

  public float[][][] pick3(
    int d2, float[][] c3, float[] c2, float[][] c1, float[][][] pa, 
    float[][][] fss) 
  {
    int nc = c2.length;
    int n3 = pa.length;
    int n2 = pa[0].length;
    int n1 = pa[0][0].length;
    SaltPicker2 sp2 = new SaltPicker2();
    float[][][] pks = new float[n2][4][];
    for (int ic=0; ic<nc; ++ic) {
      int k2 = round(c2[ic]);
      float[][] pa2 = slice2(k2,pa);
      float[][] xu = sp2.initialBoundary(1,c1[ic],c3[ic]);
      sp2.refine(50,1,4,1,xu,pa2);
      float[][] xp = copy(xu);
      float[][] xm = copy(xu);
      float[][] xt = sp2.regridBoundary(1.0f,xu[0],xu[1]);
      pks[k2] = xt;
      float[][] wx = new float[n3][n1];
      pointsToImage(xt[1],xt[0],pa2,wx);
      signAsignmentH(k2,wx,fss);
      int e2 = min(k2+d2,n2-1);
      int b2 = max(k2-d2,0);
      for (int i2=k2+1; i2<=e2; i2++) {
        pa2 = slice2(i2,pa);
        xp = sp2.pickNext(5,1,3,1,xp[0],xp[1],pa2);
        xt = sp2.regridBoundary(1.0f,xp[0],xp[1]);
        pks[i2] = xt;
        zero(wx);
        pointsToImage(xt[1],xt[0],pa2,wx);
        signAsignmentH(i2,wx,fss);
      }
      for (int i2=k2-1; i2>=b2; i2--) {
        pa2 = slice2(i2,pa);
        xm = sp2.pickNext(5,1,3,1,xm[0],xm[1],pa2);
        xt = sp2.regridBoundary(1.0f,xm[0],xm[1]);
        pks[i2] = xt;
        zero(wx);
        pointsToImage(xt[1],xt[0],pa2,wx);
        signAsignmentH(i2,wx,fss);
      }
    }
    System.out.println("fssmin="+min(fss));
    System.out.println("fssmax="+max(fss));
    return pks;
  }

  private float[][] slice2(int k2, float[][][] fx) {
    int n3 = fx.length;
    int n1 = fx[0][0].length;
    float[][] fs = new float[n3][n1];
    for (int i3=0; i3<n3; ++i3)
      for (int i1=0; i1<n1; ++i1)
        fs[i3][i1] = fx[i3][k2][i1];
    return fs;
  }


  public float[][][] pick3(
    int d3, float[] c3, float[][] c2, float[][] c1, float[][][] pa, 
    float[][][] fss) 
  {
    int nc = c3.length;
    int n3 = pa.length;
    int n2 = pa[0].length;
    int n1 = pa[0][0].length;
    SaltPicker2 sp2 = new SaltPicker2();
    float[][][] pks = new float[n3][4][];
    for (int ic=0; ic<nc; ++ic) {
      int k3 = round(c3[ic]);
      float[][] xu = sp2.initialBoundary(1,c1[ic],c2[ic]);
      sp2.refine(50,1,10,1,xu,pa[k3]);
      float[][] xp = copy(xu);
      float[][] xm = copy(xu);
      float[][] xt = sp2.regridBoundary(1.0f,xu[0],xu[1]);
      pks[k3] = xt;
      float[][] wx = new float[n2][n1];
      pointsToImage(xt[1],xt[0],pa[k3],wx);
      signAsignmentH(wx,fss[k3]);
      int e3 = min(k3+d3,n3-1);
      int b3 = max(k3-d3,0);
      for (int i3=k3+1; i3<=e3; i3++) {
        xp = sp2.pickNext(5,1,3,1,xp[0],xp[1],pa[i3]);
        xt = sp2.regridBoundary(1.0f,xp[0],xp[1]);
        pks[i3] = xt;
        zero(wx);
        pointsToImage(xt[1],xt[0],pa[i3],wx);
        signAsignmentH(wx,fss[i3]);
      }
      for (int i3=k3-1; i3>=b3; i3--) {
        xm = sp2.pickNext(5,1,3,1,xm[0],xm[1],pa[i3]);
        xt = sp2.regridBoundary(1.0f,xm[0],xm[1]);
        pks[i3] = xt;
        zero(wx);
        pointsToImage(xt[1],xt[0],pa[k3],wx);
        signAsignmentH(wx,fss[i3]);
      }
    }
    return pks;
  }


  public float[][][][] signedPoints(int dp, float d, 
    float[][][] pks, float[][][] fss) {
    int n3 = fss.length;
    int n2 = fss[0].length;
    int n1 = fss[0][0].length;
    float[][][][] ps = new float[n3][][][];
    for (int i3=0; i3<n3; ++i3) {
      float[] x1 = pks[i3][0];
      float[] x2 = pks[i3][1];
      float[] u1 = pks[i3][2];
      float[] u2 = pks[i3][3];
      float[][] fs3 = fss[i3];
      int np = x1.length;
      ArrayList<float[]> fz = new ArrayList<float[]>();
      ArrayList<float[]> fp = new ArrayList<float[]>();
      ArrayList<float[]> fn = new ArrayList<float[]>();
      for (int ip=0; ip<np; ip+=dp) {
        float x1i = x1[ip];
        float x2i = x2[ip];
        float u1i = u1[ip];
        float u2i = u2[ip];
        float x1p = x1i+u1i*d;
        float x2p = x2i+u2i*d;
        float x1m = x1i-u1i*d;
        float x2m = x2i-u2i*d;
        int i1p = round(x1p);
        int i2p = round(x2p);
        int i1m = round(x1m);
        int i2m = round(x2m);
        if(i1p<0) continue;
        if(i2p<0) continue;
        if(i1m<0) continue;
        if(i2m<0) continue;
        if(i1p>n1-1) continue;
        if(i2p>n2-1) continue;
        if(i1m>n1-1) continue;
        if(i2m>n2-1) continue;
        float fsp = fs3[i2p][i1p];
        fz.add(new float[]{x1i,x2i,i3, 0f});
        if(fsp>0f) {
          fp.add(new float[]{x1p,x2p,i3, d});
          fn.add(new float[]{x1m,x2m,i3,-d});
        } else {
          fp.add(new float[]{x1m,x2m,i3, d});
          fn.add(new float[]{x1p,x2p,i3,-d});
        }
      }
      int ns = fz.size();
      ps[i3] = new float[3][4][ns];
      float[][][] ps3 = ps[i3];
      for (int is=0; is<ns; ++is) {
        float[] fzi = fz.get(is);
        float[] fpi = fp.get(is);
        float[] fni = fn.get(is);
        ps3[0][0][is] = fzi[0]; 
        ps3[0][1][is] = fzi[1]; 
        ps3[0][2][is] = fzi[2]; 
        ps3[0][3][is] = fzi[3]; 

        ps3[1][0][is] = fpi[0]; 
        ps3[1][1][is] = fpi[1]; 
        ps3[1][2][is] = fpi[2]; 
        ps3[1][3][is] = fpi[3]; 

        ps3[2][0][is] = fni[0]; 
        ps3[2][1][is] = fni[1]; 
        ps3[2][2][is] = fni[2]; 
        ps3[2][3][is] = fni[3]; 
      }
    }
    return ps;
  }


  public float[][] regridBoundary(float d, float[] x1, float[] x2) {
    int np = x1.length;
    x1[np-1] = x1[0];
    x2[np-1] = x2[0];
    int k = 0;
    float[] ds = new float[np];
    for (int ip=1; ip<np; ++ip) {
      float x1i = x1[ip];
      float x2i = x2[ip];
      float x1m = x1[ip-1];
      float x2m = x2[ip-1];
      float dx1 = x1i-x1m;
      float dx2 = x2i-x2m;
      float dsi = sqrt(dx1*dx1+dx2*dx2);
      if(dsi>0.0f) {
        k++;
        x1[k] = x1i;
        x2[k] = x2i;
        ds[k] = ds[ip-1]+dsi;
      }
    }
    x1 = copy(k+1,x1);
    x2 = copy(k+1,x2);
    ds = copy(k+1,ds);
    double l = ds[k];
    int n = (int)round(l/d);
    Sampling ss = new Sampling(n,d,0);
    double[] sv = ss.getValues();
    CubicInterpolator cx1 = new CubicInterpolator(ds,x1);
    CubicInterpolator cx2 = new CubicInterpolator(ds,x2);
    float[] x1n = new float[n];
    float[] x2n = new float[n];
    float[] u1n = new float[n];
    float[] u2n = new float[n];
    for (int i=0; i<n; ++i) {
      float si = (float)sv[i];
      x1n[i] = cx1.interpolate(si);
      x2n[i] = cx2.interpolate(si);
    }
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(2);
    float[] g1 = new float[n];
    float[] g2 = new float[n];
    rgf.apply1(x1n,g1);
    rgf.apply1(x2n,g2);
    for (int i=0; i<n; ++i) {
      float g1i = g1[i];
      float g2i = g2[i];
      float gsi = sqrt(g1i*g1i+g2i*g2i);
      if(gsi>0.0f){g1i/=gsi;g2i/=gsi;}
      u1n[i] = -g2i;
      u2n[i] =  g1i;
    }
    return new float[][]{x1n,x2n,u1n,u2n};
  }



  public void pointsToImage(
    float[] ys, float[] zs, float[][] pa, float[][] wx) 
  {
    int n2 = pa.length;
    int n1 = pa[0].length;
    int np = ys.length;
    for (int ip=0; ip<np; ++ip) {
      int i2 = round(ys[ip]);
      int i1 = round(zs[ip]);
      i2 = max(i2,0);
      i2 = min(i2,n2-1);
      i1 = max(i1,0);
      i1 = min(i1,n1-1);
      wx[i2][i1] = pa[i2][i1];
    }
  }

  public void pointsToImage(
    float[] ys, float[] zs, float[][] wx) 
  {
    int n2 = wx.length;
    int n1 = wx[0].length;
    int np = ys.length;
    for (int ip=0; ip<np; ++ip) {
      int i2 = round(ys[ip]);
      int i1 = round(zs[ip]);
      i2 = max(i2,0);
      i2 = min(i2,n2-1);
      i1 = max(i1,0);
      i1 = min(i1,n1-1);
      wx[i2][i1] = 1;
    }
  }


  public void signAsignmentH(
    float[][] wx, float[][] fx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    for (int i1=0; i1<n1; ++i1) {
      float mk = 1f;
    for (int i2=0; i2<n2; ++i2) {
      if(wx[i2][i1]!=0f) {
        fx[i2][i1] = 0f;
        if(i2==0) mk*=-1f;
        if(i2>0&&wx[i2-1][i1]==0) mk *= -1f;
      } else {
        fx[i2][i1] = mk;
      }
    }}
    //despike
    float[][] fs = zerofloat(n1,n2);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(10);
    ref.apply(fx,fs);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(fx[i2][i1]!=0) {
        if(fs[i2][i1]>0) fx[i2][i1] =  1f;
        if(fs[i2][i1]<0) fx[i2][i1] = -1f;
      }
    }}
  }

  public void signAsignmentH(
    int k2, float[][] wx, float[][][] fx) {
    int n3 = wx.length;
    int n1 = wx[0].length;
    float[][] fx2 = new float[n3][n1];
    for (int i1=0; i1<n1; ++i1) {
      float mk = 1f;
    for (int i3=0; i3<n3; ++i3) {
      if(wx[i3][i1]!=0f) {
        fx2[i3][i1] = 0f;
        if(i3==0) mk*=-1f;
        if(i3>0&&wx[i3-1][i1]==0) mk*=-1f;
      } else {
        fx2[i3][i1] = mk;
      }
    }}
    //despike
    float[][] fs2 = zerofloat(n1,n3);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(7);
    ref.apply(fx2,fs2);
    for (int i3=0; i3<n3; ++i3) {
    for (int i1=0; i1<n1; ++i1) {
      if(fx2[i3][i1]!=0) {
        if(fs2[i3][i1]>0) fx[i3][k2][i1] =  1f;
        if(fs2[i3][i1]<0) fx[i3][k2][i1] = -1f;
      } else {
        fx[i3][k2][i1] = 1f;
      }
    }}
  }



}
