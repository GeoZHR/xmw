/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package stv;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static stv.FaultGeometry.*;

/**
 * 3D tensor voting. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.11.09
 */
public class HermiteSurface {

  public void setSigma (float sigma) {
    _sigma = sigma;
    _d3 = round(sigma);
    _d2 = round(sigma);
    _d1 = round(1.5f*sigma);
  }

  public void setWindow (int d1, int d2, int d3) {
    _d1 = d1;
    _d2 = d2;
    _d3 = d3;
  }

  public float[][][] surfer(
    final int n1, final int n2, final int n3, FaultCell[] cells) {
    int nc = cells.length;
    final float[] fs = new float[nc]; 
    final float[][] xs = new float[3][nc];
    final float[][] us = new float[3][nc];
    setKdTreeNodes(cells,xs,us,fs);
    final KdTree kt = new KdTree(xs);
    final float sigmas = _sigma*_sigma;
    final float sigmad = -20f/_sigma*_sigma;
    final float[][][] ss = new float[n3][n2][n1];
    final float scale = sigmas/120f;
    loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      System.out.println("i3="+i3);
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        int i1min = max(i1-_d1,0);
        int i2min = max(i2-_d2,0);
        int i3min = max(i3-_d3,0);
        int i1max = min(i1+_d1,n1-1);
        int i2max = min(i2+_d2,n2-1);
        int i3max = min(i3+_d3,n3-1);
        xmin[0] = i1min; xmax[0] = i1max;
        xmin[1] = i2min; xmax[1] = i2max;
        xmin[2] = i3min; xmax[2] = i3max;
        int[] id = kt.findInRange(xmin,xmax);
        int nd = id.length;
        for (int ik=0; ik<nd; ++ik) {
          int ip = id[ik];
          float fl = fs[ip];
          float x1 = xs[0][ip];
          float x2 = xs[1][ip];
          float x3 = xs[2][ip];
          float u1 = us[0][ip];
          float u2 = us[1][ip];
          float u3 = us[2][ip];
          float r1 = i1-x1;
          float r2 = i2-x2;
          float r3 = i3-x3;
          float rs = sqrt(r1*r1+r2*r2+r3*r3);
          float cs = 1-rs/_sigma;
          cs = cs*cs*cs*sigmad;
          float p1 = cs*r1; 
          float p2 = cs*r2;
          float p3 = cs*r3;
          ss[i3][i2][i1] += (u1*p1+u2*p2+u3*p3)*scale*fl;
        }
      }}
    }});
    return ss;
  }


   // Uses fault images to find cells, oriented points located on ridges.
  public FaultCell[] findCells( float fmin,
    float[][][] f, float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    // Vertical image boundaries are discontinuities that may look like
    // faults. If a fault appears to be near and nearly parallel to image
    // boundaries, then assume it is a boundary artifact and not truly a
    // fault.
    int imax = 5; // max number of samples considered to be near boundary
    float wwmax = 0.75f; // cosine of 30 degrees, squared

    // Loop over all samples. Construct cells for samples nearest to ridges.
    ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmi = f[i3m][i2 ];
        float[] fim = f[i3 ][i2m];
        float[] fip = f[i3 ][i2p];
        float[] fpi = f[i3p][i2 ];
        float[] fmm = f[i3m][i2m];
        float[] fpp = f[i3p][i2p];
        float[] fmp = f[i3m][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fii = f[i3 ][i2 ];
        float[] u1i = u1[i3][i2];
        float[] u2i = u2[i3][i2];
        float[] u3i = u3[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          float fmii = fmi[i1 ];
          float fimi = fim[i1 ];
          float fipi = fip[i1 ];
          float fpii = fpi[i1 ];
          float fmmi = fmm[i1 ];
          float fppi = fpp[i1 ];
          float fmpi = fmp[i1 ];
          float fpmi = fpm[i1 ];
          float fiii = fii[i1 ];
          float u1ii = u1i[i1 ];
          float u2ii = u2i[i1 ];
          float u3ii = u3i[i1 ];
          if (u2ii==0f&&u3ii==0f){continue;}
          if (u1ii>0.0f) {
            u1ii = -u1ii;
            u2ii = -u2ii;
            u3ii = -u3ii;
          }
          float piii = faultStrikeFromNormalVector(u1ii,u2ii,u3ii);
          // Most image samples will not have a fault cell.
          FaultCell cell = null;

          // Accumulators for ridge likelihoods and locations. Depending on
          // the limits on fault strike used below, we may find more than one
          // ridge.
          float nr = 0;
          float fl = 0.0f;
          float d2 = 0.0f;
          float d3 = 0.0f;

          // If S-N ridge, ...
          if ((fipi<fiii && fimi<fiii) &&
              ((337.5f<=piii || piii<= 22.5f) || 
               (157.5f<=piii && piii<=202.5f))) {
            float f1 = 0.5f*(fipi-fimi); // 1st derivative
            float f2 = fipi-2.0f*fiii+fimi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if (imax<=i2 && i2<n2-imax || u2ii*u2ii<=wwmax) {
                fl += fr;
                d2 += dr;
                nr += 1;
              }
            }
          }

          // If SW-NE ridge, ...
          if ((fmpi<fiii && fpmi<fiii) &&
              (( 22.5f<=piii && piii<= 67.5f) || 
               (202.5f<=piii && piii<=247.5f))) {
            float f1 = 0.5f*(fmpi-fpmi); // 1st derivative
            float f2 = fmpi-2.0f*fiii+fpmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if ((imax<=i2 && i2<n2-imax || u2ii*u2ii<=wwmax) &&
                  (imax<=i3 && i3<n3-imax || u3ii*u3ii<=wwmax)) {
                fl += fr;
                d2 += dr;
                d3 -= dr;
                nr += 1;
              }
            }
          }

          // If W-E ridge, ...
          if ((fpii<fiii && fmii<fiii) &&
              (( 67.5f<=piii && piii<=112.5f) ||
               (247.5f<=piii && piii<=292.5f))) {
            float f1 = 0.5f*(fpii-fmii); // 1st derivative
            float f2 = fmii-2.0f*fiii+fpii; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if (imax<=i3 && i3<n3-imax || u3ii*u3ii<=wwmax) {
                fl += fr;
                d3 += dr;
                nr += 1;
              }
            }
          }

          // If NW-SE ridge, ...
          if ((fppi<fiii && fmmi<fiii) &&
              ((112.5f<=piii && piii<=157.5f) || 
               (292.5f<=piii && piii<=337.5f))) {
            float f1 = 0.5f*(fppi-fmmi); // 1st derivative
            float f2 = fppi-2.0f*fiii+fmmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=fmin) {
              if ((imax<=i2 && i2<n2-imax || u2ii*u2ii<=wwmax) &&
                  (imax<=i3 && i3<n3-imax || u3ii*u3ii<=wwmax)) {
                fl += fr;
                d2 += dr;
                d3 += dr;
                nr += 1;
              }
            }
          }

          // If at least one ridge, construct a cell and add to list.
          if (nr>0) {
            fl /= nr;
            d2 /= nr;
            d3 /= nr;
            float tiii = faultDipFromNormalVector(u1ii,u2ii,u3ii);
            cell = new FaultCell(i1,i2+d2,i3+d3,fl,piii,tiii);
            cellList.add(cell);
          }
        }
      }
    }
    return cellList.toArray(new FaultCell[0]);
  }

  public FaultCell[] findCells (
    float[][][] ss, float[][][] u1, float[][][] u2, float[][][] u3) {
    ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    int n3 = ss.length;
    int n2 = ss[0].length;
    int n1 = ss[0][0].length;
    SincInterpolator si = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    for (int i3=0; i3<n3 ;++i3) {
    for (int i2=0; i2<n2 ;++i2) {
    for (int i1=0; i1<n1 ;++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float x1m = i1-u1i;
      float x2m = i2-u2i;
      float x3m = i3-u3i;
      float x1p = i1+u1i;
      float x2p = i2+u2i;
      float x3p = i3+u3i;
      float sxi = ss[i3][i2][i1];
      float sxm = si.interpolate(s1,s2,s3,ss,x1m,x2m,x3m);
      float sxp = si.interpolate(s1,s2,s3,ss,x1p,x2p,x3p);
      if (sxi>sxm && sxi>sxp && sxi>0.1f) {
        if (u1i>0.0f) {
          u1i = -u1i;
          u2i = -u2i;
          u3i = -u3i;
        }
        float ft = faultDipFromNormalVector(u1i,u2i,u3i);
        float fp = faultStrikeFromNormalVector(u1i,u2i,u3i);
        FaultCell cell = new FaultCell(i1,i2,i3,sxi,fp,ft);
        cellList.add(cell);
      }
    }}}
    return cellList.toArray(new FaultCell[0]);
  }

  public float[][][] direvative(
    float[][][] sm, float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = sm.length;
    int n2 = sm[0].length;
    int n1 = sm[0][0].length;
    float[][][] sp = new float[n3][n2][n1];
    sm = sub(sm,min(sm));
    sm = div(sm,max(sm));
    SincInterpolator si = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    for (int i3=0; i3<n3 ;++i3) {
    for (int i2=0; i2<n2 ;++i2) {
    for (int i1=0; i1<n1 ;++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float x1m = i1-u1i;
      float x2m = i2-u2i;
      float x3m = i3-u3i;
      float x1p = i1+u1i;
      float x2p = i2+u2i;
      float x3p = i3+u3i;
      float sxi = sm[i3][i2][i1];
      float sxm = si.interpolate(s1,s2,s3,sm,x1m,x2m,x3m);
      float sxp = si.interpolate(s1,s2,s3,sm,x1p,x2p,x3p);
      sp[i3][i2][i1] = (sxp-sxm)*0.5f;
      if(sxi<0.2f) {sp[i3][i2][i1] = -30f;}
    }}}
    return sp;
  }

  public float[][][] findRidges (
    float[][][] sm, float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = sm.length;
    int n2 = sm[0].length;
    int n1 = sm[0][0].length;
    float[][][] sp = new float[n3][n2][n1];
    SincInterpolator si = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    for (int i3=0; i3<n3 ;++i3) {
    for (int i2=0; i2<n2 ;++i2) {
    for (int i1=0; i1<n1 ;++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float x1m = i1-u1i;
      float x2m = i2-u2i;
      float x3m = i3-u3i;
      float x1p = i1+u1i;
      float x2p = i2+u2i;
      float x3p = i3+u3i;
      float sxi = sm[i3][i2][i1];
      float sxm = si.interpolate(s1,s2,s3,sm,x1m,x2m,x3m);
      float sxp = si.interpolate(s1,s2,s3,sm,x1p,x2p,x3p);
      if (sxi>sxm && sxi>sxp) {
        sp[i3][i2][i1] = sxi;
      }
    }}}
    return sp;
  }

  public EigenTensors3 mark (
    float[][][] g11, float[][][] g12, float[][][] g13, 
    float[][][] g22, float[][][] g23, float[][][] g33, FaultCell[] cells) {
    int n3 = g11.length;
    int n2 = g11[0].length;
    int n1 = g11[0][0].length;
    float[][][] dx = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    for (FaultCell fc:cells) {
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      float w1 = fc.w1;
      float w2 = fc.w2;
      float w3 = fc.w3;
      float fl = fc.fl;
      g11[i3][i2][i1] = w1*w1*fl;
      g12[i3][i2][i1] = w1*w2*fl;
      g13[i3][i2][i1] = w1*w3*fl;
      g22[i3][i2][i1] = w2*w2*fl;
      g23[i3][i2][i1] = w2*w3*fl;
      g33[i3][i2][i1] = w3*w3*fl;
    }
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(0.0f,g11,dx,k1,k2,k3);
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] eu = fillfloat(0.01f,n1,n2,n3);
    float[][][] ev = fillfloat(1.00f,n1,n2,n3);
    float[][][] ew = fillfloat(1.00f,n1,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int c1 = k1[i3][i2][i1];
      int c2 = k2[i3][i2][i1];
      int c3 = k3[i3][i2][i1];
      float g11i = g11[c3][c2][c1];
      float g12i = g12[c3][c2][c1];
      float g13i = g13[c3][c2][c1];
      float g22i = g22[c3][c2][c1];
      float g23i = g23[c3][c2][c1];
      float g33i = g33[c3][c2][c1];
      double[][] a = new double[3][3];
      a[0][0] = g11i; a[0][1] = g12i; a[0][2] = g13i;
      a[1][0] = g12i; a[1][1] = g22i; a[1][2] = g23i;
      a[2][0] = g13i; a[2][1] = g23i; a[2][2] = g33i;
      double[] e = new double[3];
      double[][] z = new double[3][3];
      Eigen.solveSymmetric33(a,z,e);
      float u1i = (float)z[0][0];
      float u2i = (float)z[0][1];
      float w1i = (float)z[2][0];
      float w2i = (float)z[2][1];
      u1[i3][i2][i1] = u1i;
      u2[i3][i2][i1] = u2i;
      w1[i3][i2][i1] = w1i;
      w2[i3][i2][i1] = w2i;
    }}}
    return new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,false);
  }

  public float[][][][] applyVoteX(float sigma,
    final int n1, final int n2, final int n3, final FaultCell[] cells) {
    int nc = cells.length;
    final float[] w11s = new float[nc];
    final float[] w12s = new float[nc];
    final float[] w13s = new float[nc];
    final float[] w22s = new float[nc];
    final float[] w23s = new float[nc];
    final float[] w33s = new float[nc];
    final float[][][] g11 = new float[n3][n2][n1];
    final float[][][] g12 = new float[n3][n2][n1];
    final float[][][] g13 = new float[n3][n2][n1];
    final float[][][] g22 = new float[n3][n2][n1];
    final float[][][] g23 = new float[n3][n2][n1];
    final float[][][] g33 = new float[n3][n2][n1];
    final float sigmas = 0.5f/(_sigma*_sigma);
    precompute(cells,w11s,w12s,w13s,w22s,w23s,w33s);
    for (int ic=0; ic<nc; ++ic) {
      FaultCell cell = cells[ic];
      final int c1 = cell.i1;
      final int c2 = cell.i2;
      final int c3 = cell.i3;
      final float w1 = cell.w1;
      final float w2 = cell.w2;
      final float w3 = cell.w3;
      final float w11 = w11s[ic];
      final float w12 = w12s[ic];
      final float w13 = w13s[ic];
      final float w22 = w22s[ic];
      final float w23 = w23s[ic];
      final float w33 = w33s[ic];
      loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          float d1 = i1-c1;
          float d2 = i2-c2;
          float d3 = i3-c3;
          float ds = sqrt(d1*d1+d2*d2+d3*d3);
          if(ds==0f) {
            g11[i3][i2][i1] += w11;
            g12[i3][i2][i1] += w12;
            g13[i3][i2][i1] += w13;
            g22[i3][i2][i1] += w22;
            g23[i3][i2][i1] += w23;
            g33[i3][i2][i1] += w33;
          } else {
            d1 /= ds; d2 /= ds; d3 /= ds;
            float wd = w1*d1+w2*d2+w3*d3;
            if(abs(wd)>0.4f){continue;}
            float sc = exp(-ds*ds*sigmas)*pow((1f-wd*wd),4);
            g11[i3][i2][i1] += w11*sc;
            g12[i3][i2][i1] += w12*sc;
            g13[i3][i2][i1] += w13*sc;
            g22[i3][i2][i1] += w22*sc;
            g23[i3][i2][i1] += w23*sc;
            g33[i3][i2][i1] += w33*sc;
          }
        }}
      }});
    }
    return solveEigenproblems(g11,g12,g13,g22,g23,g33);
  }

  public float[][][][] solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33) {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    final float[][][] ss = new float[n3][n2][n1];
    final float[][][] cs = new float[n3][n2][n1];
    final float[][][] js = new float[n3][n2][n1];
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        double[] e = new double[3];
        double[][] z = new double[3][3];
        double[][] a = new double[3][3];
        float g11i = g11[i3][i2][i1];
        float g12i = g12[i3][i2][i1];
        float g13i = g13[i3][i2][i1];
        float g22i = g22[i3][i2][i1];
        float g23i = g23[i3][i2][i1];
        float g33i = g33[i3][i2][i1];
        a[0][0]=g11i; a[0][1]=g12i; a[0][2]=g13i;
        a[1][0]=g12i; a[1][1]=g22i; a[1][2]=g23i;
        a[2][0]=g13i; a[2][1]=g23i; a[2][2]=g33i;
        Eigen.solveSymmetric33(a,z,e);
        float eui = (float)e[0];
        float evi = (float)e[1];
        float ewi = (float)e[2];
        if (ewi<0.0f) ewi = 0.0f;
        if (evi<ewi) evi = ewi;
        if (eui<evi) eui = evi;
        ss[i3][i2][i1] = eui-evi;
        cs[i3][i2][i1] = evi-ewi;
        js[i3][i2][i1] = ewi;
      }}
    }});
    return new float[][][][]{ss,cs,js};

  }


  private void precompute (
    final FaultCell[] cells,
    final float[] w11s, final float[] w12s, final float[] w13s, 
    final float[] w22s, final float[] w23s, final float[] w33s) {
    final int nc = cells.length;
    loop(nc,new Parallel.LoopInt() {
    public void compute(int ic) {
      FaultCell cell = cells[ic];
      float fl = cell.fl;
      float w1 = cell.w1;
      float w2 = cell.w2;
      float w3 = cell.w3;
      w11s[ic] = fl*w1*w1;
      w12s[ic] = fl*w1*w2;
      w13s[ic] = fl*w1*w3;
      w22s[ic] = fl*w2*w2;
      w23s[ic] = fl*w2*w3;
      w33s[ic] = fl*w3*w3;
    }});
  }

  public FaultCell[] getSparseCells (
    int n1, int n2, int n3, int d1, int d2, int d3, FaultCell[] cells) {
    FaultCell[][][] fcg = new FaultCell[n3][n2][n1];
    for (FaultCell fc:cells) {
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      fcg[i3][i2][i1] = fc;
    }
    ArrayList<FaultCell> fcs = new ArrayList<FaultCell>();
    for (int i3=0; i3<n3; i3+=d3) {
    for (int i2=0; i2<n2; i2+=d2) {
    for (int i1=0; i1<n1; i1+=d1) {
      if (fcg[i3][i2][i1]!=null) {
        fcs.add(fcg[i3][i2][i1]);
      }
    }}}
    return fcs.toArray(new FaultCell[0]);
  }


  
  private float[][] solveEigenproblems(double[][] a) {
    double[] e = new double[3];
    double[][] z = new double[3][3];
    Eigen.solveSymmetric33(a,z,e);
    float eui = (float)e[0];
    float evi = (float)e[1];
    float ewi = (float)e[2];
    float u1i = (float)z[0][0];
    float u2i = (float)z[0][1];
    float u3i = (float)z[0][2];
    if (ewi<0.0f) ewi = 0.0f;
    if (evi<ewi) evi = ewi;
    if (eui<evi) eui = evi;
    float[] es = new float[]{eui,evi,ewi};
    float[] us = new float[]{u1i,u2i,u3i};
    return new float[][]{us,es};
  }

  private void setKdTreeNodes(
    FaultCell[] cells, float[][] xs, float[][] us, float[] fls) {
    int nc = cells.length;
    for (int ic=0; ic<nc; ic++) {
      FaultCell fc = cells[ic];
      float w1 = fc.w1;
      float w2 = fc.w2;
      float w3 = fc.w3;
      float fl = fc.fl;
      fls[ic] = fl;
      xs[0][ic] = fc.i1;
      xs[1][ic] = fc.i2;
      xs[2][ic] = fc.i3;
      us[0][ic] = w1;
      us[1][ic] = w2;
      us[2][ic] = w3;
    }
  }


  private float _sigma=20f;
  private int _d1 = 10;
  private int _d2 = 10;
  private int _d3 = 10;
}
