/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.09.16
 */
public class FloatScaleSurface {

  public float[][][] findScalarField(int n1, int n2, int n3, FaultCell[] fc) {
    return scalarField(n1,n2,n3,fc);
    //return accumulateGaussian(n1,n2,n3,fc);
  }

  public FaultCell[] getCells(FaultSkin[] fs) {
    FaultSkin[] fsi = new FaultSkin[2];
    fsi[0] = fs[2];
    fsi[1] = fs[3];
    return FaultSkin.getCells(fs);
  }

  private void setKdTreePointsM(FaultCell[] fc,
    float[][] xf, float[][] xs, float[][] ws, float[][] us, float[][] vs, float[] w) {
    int nc = fc.length;
    for (int ic=0; ic<nc; ic++) {
      w[ic] = fc[ic].fl;
      float x1 = fc[ic].x1;
      float x2 = fc[ic].x2;
      float x3 = fc[ic].x3;

      float w1 = fc[ic].w1;
      float w2 = fc[ic].w2;
      float w3 = fc[ic].w3;

      float u1 = fc[ic].u1;
      float u2 = fc[ic].u2;
      float u3 = fc[ic].u3;

      float v1 = fc[ic].v1;
      float v2 = fc[ic].v2;
      float v3 = fc[ic].v3;

      xf[0][ic] = x1;
      xf[1][ic] = x2;
      xf[2][ic] = x3;
    

      xs[ic][0] = x1;
      xs[ic][1] = x2;
      xs[ic][2] = x3;

      ws[ic][0] = w1*w1;
      ws[ic][1] = w2*w2;
      ws[ic][2] = w3*w3;
      ws[ic][3] = w1*w2;
      ws[ic][4] = w1*w3;
      ws[ic][5] = w2*w3;

      us[ic][0] = u1*u1;
      us[ic][1] = u2*u2;
      us[ic][2] = u3*u3;
      us[ic][3] = u1*u2;
      us[ic][4] = u1*u3;
      us[ic][5] = u2*u3;

      vs[ic][0] = v1*v1;
      vs[ic][1] = v2*v2;
      vs[ic][2] = v3*v3;
      vs[ic][3] = v1*v2;
      vs[ic][4] = v1*v3;
      vs[ic][5] = v2*v3;
    }
  }

  public float[][][] scalarFieldM(final int n1, final int n2, final int n3, 
    final FaultCell[] fc, final float[][][] fp, final float[][][] ft) {
    final int d2=10;
    final int d3=10;
    final int d1=10;
    //final float v = -1.f;
    final float v = 0.0f;
    float sigmaNor = 4.0f;
    float sigmaPhi = 10.0f;
    float sigmaTheta = 20.0f;
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    final float sv = 1.0f/(sigmaPhi*sigmaPhi); 
    final float su = 1.0f/(sigmaTheta*sigmaTheta); 
    int nc = fc.length;
    final float[] wp = new float[nc];
    final float[][] xf = new float[3][nc];
    final float[][] xs = new float[nc][3];
    final float[][] ws = new float[nc][6];
    final float[][] us = new float[nc][6];
    final float[][] vs = new float[nc][6];
    setKdTreePointsM(fc,xf,xs,ws,us,vs,wp);
    final int[] bs1 = setBounds(n1,xf[0]);
    final int[] bs2 = setBounds(n2,xf[1]);
    final int[] bs3 = setBounds(n3,xf[2]);
    final KdTree kt = new KdTree(xf);
    final float[][][] sf = fillfloat(v,n1,n2,n3);
    final float[][][] fpp = copy(fp);
    final float[][][] ftp = copy(ft);
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      System.out.println("i3="+i3);
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          float[] y = new float[]{i1,i2,i3};
          int ne = kt.findNearest(y);
          float x1 = xf[0][ne];
          float x2 = xf[1][ne];
          float x3 = xf[2][ne];
          float dd = distance(new float[]{x1,x2,x3},y);
          if(dd>50.0f){continue;}
          getRange(d1,d2,d3,i1,i2,i3,n1,n2,n3,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<1){continue;}
          sf[i3][i2][i1] = 0.0f;
          float wps = 0.0f;
          float fpa = 0.0f;
          float fta = 0.0f;
          float fps = 0.0f;
          float fts = 0.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float wpi = pow(wp[ip],8.0f);
            float w1s = ws[ip][0];
            float w2s = ws[ip][1];
            float w3s = ws[ip][2];
            float w12 = ws[ip][3];
            float w13 = ws[ip][4];
            float w23 = ws[ip][5];

            float u1s = us[ip][0];
            float u2s = us[ip][1];
            float u3s = us[ip][2];
            float u12 = us[ip][3];
            float u13 = us[ip][4];
            float u23 = us[ip][5];

            float v1s = vs[ip][0];
            float v2s = vs[ip][1];
            float v3s = vs[ip][2];
            float v12 = vs[ip][3];
            float v13 = vs[ip][4];
            float v23 = vs[ip][5];

            float dx1 = i1-xs[ip][0];
            float dx2 = i2-xs[ip][1];
            float dx3 = i3-xs[ip][2];

            float d11 = dx1*dx1;
            float d12 = dx1*dx2;
            float d13 = dx1*dx3;
            float d22 = dx2*dx2;
            float d23 = dx2*dx3;
            float d33 = dx3*dx3;

            float g11 = w1s*d11+w12*d12+w13*d13;
            float g12 = u1s*d11+u12*d12+u13*d13;
            float g13 = v1s*d11+v12*d12+v13*d13;
            float g21 = w12*d12+w2s*d22+w23*d23;
            float g22 = u12*d12+u2s*d22+u23*d23;
            float g23 = v12*d12+v2s*d22+v23*d23;
            float g31 = w13*d13+w23*d23+w3s*d33;
            float g32 = u13*d13+u23*d23+u3s*d33;
            float g33 = v13*d13+v23*d23+v3s*d33;
            float gss = 0.0f;
            gss += (g11+g21+g31)*sw;
            gss += (g12+g22+g32)*su;
            gss += (g13+g23+g33)*sv;
            float sfi = exp(-gss)*wpi;
            sf[i3][i2][i1] += sfi;

            int j1 = fc[ip].i1;
            int j2 = fc[ip].i2;
            int j3 = fc[ip].i3;

            int k1 = fc[ne].i1;
            int k2 = fc[ne].i2;
            int k3 = fc[ne].i3;
            float fpr = fpp[k3][k2][k1];
            float fpi = fpp[j3][j2][j1];
            if(abs(fpr-fpi)<=30) {
              fps += sfi;
              fpa += fpp[j3][j2][j1]*sfi;
            }
            fts += sfi;
            fta += ftp[j3][j2][j1]*sfi;
            wps += wpi;
          }
          ft[i3][i2][i1] = fta/fts;
          fp[i3][i2][i1] = fpa/fps;
          sf[i3][i2][i1] /= wps;
        }
      }
    }});
    return sf;
  }



  private void setKdTreePoints(FaultCell[] fc,
    float[][] xf, float[][] xs, float[][] ws, float[][] us, float[][] vs, float[] w) {
    int nc = fc.length;
    for (int ic=0; ic<nc; ic++) {
      w[ic] = fc[ic].fl;
      float x1 = fc[ic].x1;
      float x2 = fc[ic].x2;
      float x3 = fc[ic].x3;

      float w1 = fc[ic].w1;
      float w2 = fc[ic].w2;
      float w3 = fc[ic].w3;

      float u1 = fc[ic].u1;
      float u2 = fc[ic].u2;
      float u3 = fc[ic].u3;

      float v1 = fc[ic].v1;
      float v2 = fc[ic].v2;
      float v3 = fc[ic].v3;

      xf[0][ic] = x1;
      xf[1][ic] = x2;
      xf[2][ic] = x3;
    

      xs[ic][0] = x1-w1;
      xs[ic][1] = x2-w2;
      xs[ic][2] = x3-w3;

      xs[ic][3] = x1+w1;
      xs[ic][4] = x2+w2;
      xs[ic][5] = x3+w3;


      ws[ic][0] = w1*w1;
      ws[ic][1] = w2*w2;
      ws[ic][2] = w3*w3;
      ws[ic][3] = w1*w2;
      ws[ic][4] = w1*w3;
      ws[ic][5] = w2*w3;

      us[ic][0] = u1*u1;
      us[ic][1] = u2*u2;
      us[ic][2] = u3*u3;
      us[ic][3] = u1*u2;
      us[ic][4] = u1*u3;
      us[ic][5] = u2*u3;

      vs[ic][0] = v1*v1;
      vs[ic][1] = v2*v2;
      vs[ic][2] = v3*v3;
      vs[ic][3] = v1*v2;
      vs[ic][4] = v1*v3;
      vs[ic][5] = v2*v3;
    }
  }

  public float[][][] scalarField(
    final int n1, final int n2, final int n3, FaultCell[] fc) 
  {
    final int d2=10;
    final int d3=10;
    final int d1=20;
    final float v = -1.f;
    float sigmaNor = 4.0f;
    float sigmaPhi = 6.0f;
    float sigmaTheta = 20.0f;
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    final float sv = 1.0f/(sigmaPhi*sigmaPhi); 
    final float su = 1.0f/(sigmaTheta*sigmaTheta); 
    int nc = fc.length;
    final float[] wp = new float[nc];
    final float[][] xf = new float[3][nc];
    final float[][] xs = new float[nc][6];
    final float[][] ws = new float[nc][6];
    final float[][] us = new float[nc][6];
    final float[][] vs = new float[nc][6];
    setKdTreePoints(fc,xf,xs,ws,us,vs,wp);
    final int[] bs1 = setBounds(n1,xf[0]);
    final int[] bs2 = setBounds(n2,xf[1]);
    final int[] bs3 = setBounds(n3,xf[2]);
    final KdTree kt = new KdTree(xf);
    final float[][][] sf = fillfloat(v,n1,n2,n3);
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      System.out.println("i3="+i3);
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          float[] y = new float[]{i1,i2,i3};
          int ne = kt.findNearest(y);
          float x1 = xf[0][ne];
          float x2 = xf[1][ne];
          float x3 = xf[2][ne];
          float dd = distance(new float[]{x1,x2,x3},y);
          if(dd>16.0f){continue;}
          getRange(d1,d2,d3,i1,i2,i3,n1,n2,n3,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<50){continue;}
          sf[i3][i2][i1] = 0.0f;
          float wps = 0.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            float wpi = pow(wp[ip],8.0f);
            float w1s = ws[ip][0];
            float w2s = ws[ip][1];
            float w3s = ws[ip][2];
            float w12 = ws[ip][3];
            float w13 = ws[ip][4];
            float w23 = ws[ip][5];

            float u1s = us[ip][0];
            float u2s = us[ip][1];
            float u3s = us[ip][2];
            float u12 = us[ip][3];
            float u13 = us[ip][4];
            float u23 = us[ip][5];

            float v1s = vs[ip][0];
            float v2s = vs[ip][1];
            float v3s = vs[ip][2];
            float v12 = vs[ip][3];
            float v13 = vs[ip][4];
            float v23 = vs[ip][5];

            float x1m = i1-xs[ip][0];
            float x2m = i2-xs[ip][1];
            float x3m = i3-xs[ip][2];

            float x1p = i1-xs[ip][3];
            float x2p = i2-xs[ip][4];
            float x3p = i3-xs[ip][5];

            float x1sm = x1m*x1m;
            float x2sm = x2m*x2m;
            float x3sm = x3m*x3m;
            float x12m = x1m*x2m;
            float x13m = x1m*x3m;
            float x23m = x2m*x3m;

            float x1sp = x1p*x1p;
            float x2sp = x2p*x2p;
            float x3sp = x3p*x3p;
            float x12p = x1p*x2p;
            float x13p = x1p*x3p;
            float x23p = x2p*x3p;

            float g11m = w1s*x1sm+w12*x12m+w13*x13m;
            float g12m = u1s*x1sm+u12*x12m+u13*x13m;
            float g13m = v1s*x1sm+v12*x12m+v13*x13m;
            float g21m = w12*x12m+w2s*x2sm+w23*x23m;
            float g22m = u12*x12m+u2s*x2sm+u23*x23m;
            float g23m = v12*x12m+v2s*x2sm+v23*x23m;
            float g31m = w13*x13m+w23*x23m+w3s*x3sm;
            float g32m = u13*x13m+u23*x23m+u3s*x3sm;
            float g33m = v13*x13m+v23*x23m+v3s*x3sm;

            float g11p = w1s*x1sp+w12*x12p+w13*x13p;
            float g12p = u1s*x1sp+u12*x12p+u13*x13p;
            float g13p = v1s*x1sp+v12*x12p+v13*x13p;
            float g21p = w12*x12p+w2s*x2sp+w23*x23p;
            float g22p = u12*x12p+u2s*x2sp+u23*x23p;
            float g23p = v12*x12p+v2s*x2sp+v23*x23p;
            float g31p = w13*x13p+w23*x23p+w3s*x3sp;
            float g32p = u13*x13p+u23*x23p+u3s*x3sp;
            float g33p = v13*x13p+v23*x23p+v3s*x3sp;

            float gssm = (g11m+g21m+g31m)*sw+(g12m+g22m+g32m)*su+(g13m+g23m+g33m)*sv;
            float gssp = (g11p+g21p+g31p)*sw+(g12p+g22p+g32p)*su+(g13p+g23p+g33p)*sv;

            sf[i3][i2][i1] += (exp(-gssp)-exp(-gssm))*wpi;
            wps += wpi;
          }
          sf[i3][i2][i1] /= wps;
        }
      }
    }});
    return sf;
  }

  private float _sigma1 = 1.0f/15.0f; 
  private float _sigma2 = 1.0f/15.0f; 
  private float _sigma3 = 1.0f/15.0f; 
  private float _sigmas1 = _sigma1*_sigma1;
  private float _sigmas2 = _sigma2*_sigma2;
  private float _sigmas3 = _sigma3*_sigma3;
  private float _twoPi = (float)(0.5/Math.PI);
  private float _pi = (float)Math.PI;
  private float _sigma = 4.0f;
  private float _sigmas = 1.0f/(_sigma*_sigma);
  private float _scale = 1.0f/_sigma/sqrt(2.0f*_pi);
  private float basicFunction(float d1, float d2, float d3) {
    float f2 = _sigma2*exp(-d2*d2*_sigmas3*0.5f);
    float f3 = _sigma3*exp(-d3*d3*_sigmas2*0.5f);
    float f1 = _sigmas1*d1*exp(-d1*d1*_sigmas1*0.5f);
    return f1*f2*f3*_twoPi;
  }

  private float weight(float d) {
    float w = _scale*exp(-d*d*_sigmas*0.5f);
    return w;
  }


  private int[] setBounds(int n, float[] x) {
    int[] bs = new int[2];
    int n1m = (int)min(x)-5; 
    int n1p = (int)max(x)+5; 
    if(n1m<0){n1m=0;}
    if(n1p>n){n1p=n;}
    bs[0] = n1m;
    bs[1] = n1p;
    return bs;
  }

  private static float distance(float[] x, float[] y) {
    float d1 = y[0]-x[0];
    float d2 = y[1]-x[1];
    float d3 = y[2]-x[2];
    return sqrt(d1*d1+d2*d2+d3*d3);
  }


  private static float distance(float[] x, float[] u, float[] y) {
    float d1 = y[0]-x[0];
    float d2 = y[1]-x[1];
    float d3 = y[2]-x[2];
    float du = d1*u[0]+d2*u[1]+d3*u[2];
    return du;
  }


  private static void getRange(int d1, int d2, int d3, int i1, int i2, int i3, 
    int n1, int n2, int n3, float[] xmin, float[] xmax) 
  {
    int i1m = i1-d1; if(i1m<0){i1m=0;}
    int i2m = i2-d2; if(i2m<0){i2m=0;}
    int i3m = i3-d3; if(i3m<0){i3m=0;}
    int i1p = i1+d1; if(i1p>=n1){i1p=n1-1;}
    int i2p = i2+d2; if(i2p>=n2){i2p=n2-1;}
    int i3p = i3+d3; if(i3p>=n3){i3p=n3-1;}
    xmin[0] = i1m;
    xmin[1] = i2m;
    xmin[2] = i3m;
    xmax[0] = i1p;
    xmax[1] = i2p;
    xmax[2] = i3p;
  }


}
