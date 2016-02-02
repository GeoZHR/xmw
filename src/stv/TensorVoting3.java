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
import edu.mines.jtk.util.Stopwatch;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.*;
import ipfx.FaultCell;
import static ipfx.FaultGeometry.*;

/**
 * 3D tensor voting. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.11.09
 */
public class TensorVoting3 {

  public void setSigma (float sigma) {
    _sigma = sigma;
    _w3 = round(sigma);
    _w2 = round(sigma);
    _w1 = round(1.5f*sigma);
  }

  public void setVoteWindow(int w1, int w2, int w3) {
    _w1 = w1;
    _w2 = w2;
    _w3 = w3;
  }

  public float[][][][] applyVote(
    final int n1, final int n2, final int n3, FaultCell[] cells) {
    int nc = cells.length;
    final float[] fs = new float[nc]; 
    final float[][] xs = new float[3][nc];
    final float[][] us = new float[3][nc];
    setKdTreeNodes(cells,xs,us,fs);
    final KdTree kt = new KdTree(xs);
    final float[][][] ws = gaussianWeights();
    final float[][][] ss = new float[n3][n2][n1];
    final float[][][] cs = new float[n3][n2][n1];
    final float[][][] fp = new float[n3][n2][n1];
    final float[][][] ft = new float[n3][n2][n1];
    final int[] c = new int[1];
    loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      c[0] = c[0]+1;
      System.out.println("c="+c[0]);
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        int i1min = max(i1-_w1,0);
        int i2min = max(i2-_w2,0);
        int i3min = max(i3-_w3,0);
        int i1max = min(i1+_w1,n1-1);
        int i2max = min(i2+_w2,n2-1);
        int i3max = min(i3+_w3,n3-1);
        xmin[0] = i1min; xmax[0] = i1max;
        xmin[1] = i2min; xmax[1] = i2max;
        xmin[2] = i3min; xmax[2] = i3max;
        int[] id = kt.findInRange(xmin,xmax);
        int nd = id.length;
        float g11 = 0.0f;
        float g22 = 0.0f;
        float g33 = 0.0f;
        float g12 = 0.0f;
        float g13 = 0.0f;
        float g23 = 0.0f;
        if(nd<1) continue;
        double[][] a = new double[3][3];
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
          int k1 = round(abs(r1));
          int k2 = round(abs(r2));
          int k3 = round(abs(r3));
          float wsi = ws[k3][k2][k1]*fl*fl;
          if (rs==0) {
            g11 += u1*u1*wsi;
            g12 += u1*u2*wsi;
            g13 += u1*u3*wsi;
            g22 += u2*u2*wsi;
            g23 += u2*u3*wsi;
            g33 += u3*u3*wsi;
          } else {
            r1 /= rs; r2 /= rs; r3 /= rs;
            float ur = u1*r1+u2*r2+u3*r3;
            if (abs(ur)>0.5f){continue;}
            float v1 = u1;
            float v2 = u2;
            float v3 = u3;
            if(abs(ur)>0.0001f) {
              float cx = 0.5f*rs/ur; // find a better way?
              float c1 = x1+u1*cx;
              float c2 = x2+u2*cx;
              float c3 = x3+u3*cx;
              v1 = c1-i1; v2 = c2-i2; v3 = c3-i3;
              float vs = 1.0f/sqrt(v1*v1+v2*v2+v3*v3);
              v1 *= vs; v2 *= vs; v3 *= vs; 
            }
            ur = 1f-ur*ur;
            ur *= ur;
            ur *= ur;
            ur *= ur;
            ur *= ur;
            ur *= ur;
            ur *= ur;
            float sc = wsi*ur;
            g11 += sc*v1*v1;
            g12 += sc*v1*v2;
            g13 += sc*v1*v3;
            g22 += sc*v2*v2;
            g23 += sc*v2*v3;
            g33 += sc*v3*v3;
          }
        }
        a[0][0] = g11; a[0][1] = g12; a[0][2] = g13;
        a[1][0] = g12; a[1][1] = g22; a[1][2] = g23;
        a[2][0] = g13; a[2][1] = g23; a[2][2] = g33;
        float[][] ue = solveEigenproblems(a);
        float eui = ue[1][0];
        float evi = ue[1][1];
        float ewi = ue[1][2];
        float u1i = ue[0][0];
        float u2i = ue[0][1];
        float u3i = ue[0][2];
        if (u2i==0.0f&&u3i==0f){continue;}
        if (u1i>0.0f) {
          u1i = -u1i;
          u2i = -u2i;
          u3i = -u3i;
        }
        ss[i3][i2][i1] = (eui-evi);
        cs[i3][i2][i1] = (evi-ewi)*(eui-evi);
        ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
        fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
      }}
    }});
    sub(ss,min(ss),ss);
    div(ss,max(ss),ss);
    sub(cs,min(cs),cs);
    div(cs,max(cs),cs);
    return new float[][][][]{ss,cs,fp,ft};
  }

  private float[][][] gaussianWeights() {
    int n1 = _w1*2+1;
    int n2 = _w2*2+1;
    int n3 = _w3*2+1;
    float[][][] fx = new float[n3][n2][n1];
    fx[_w3][_w2][_w1] = 1f;
    float[][][] fs = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(_sigma);
    rgf.apply000(fx,fs);
    return copy(_w1+1,_w2+1,_w3+1,_w1,_w2,_w3,fs);
  }

  public float[][][][] applyVoteFast(
    final int n1, final int n2, final int n3, 
    Sampling sp, Sampling st, FaultCell[] fcs) {
      int np = sp.getCount();
      int nt = st.getCount();
      float[][][][][] scs = voteScale(sp,st,_sigma);
      FaultCell[][][] fca = getCellArray(n1,n2,n3,sp,st,fcs);
      final float[][][] g11 = new float[n3][n2][n1];
      final float[][][] g12 = new float[n3][n2][n1];
      final float[][][] g13 = new float[n3][n2][n1];
      final float[][][] g22 = new float[n3][n2][n1];
      final float[][][] g23 = new float[n3][n2][n1];
      final float[][][] g33 = new float[n3][n2][n1];
      Stopwatch sw = new Stopwatch();
      sw.start();
      for (int ip=0; ip<np-2; ++ip) {
        final float phi = (float)sp.getValue(ip);
        if (ip>0) {
          double timeUsed = sw.time();
          double timeLeft = ((double)np/(double)ip-1.0)*timeUsed;
          int timeLeftSec = 1+(int)timeLeft;
          trace("FaultScanner.scan: done in "+timeLeftSec+" seconds");
        }
        for (int it=0; it<nt; ++it) {
          FaultCell[] fc = fca[ip][it];
          int nc = fc.length;
          if(nc==0) {continue;}
          final float[][][] sci = scs[ip][it];
          float theta = (float)st.getValue(it);
          float[] us = faultNormalVectorFromStrikeAndDip(phi,theta);
          final float u11 = us[0]*us[0];
          final float u12 = us[0]*us[1];
          final float u13 = us[0]*us[2];
          final float u22 = us[1]*us[1];
          final float u23 = us[1]*us[2];
          final float u33 = us[2]*us[2];
          for (int ic=0; ic<nc; ++ic) {
            FaultCell cell = fc[ic];
            float fli = cell.getFl();
            final int c1 = cell.getI1();
            final int c2 = cell.getI2();
            final int c3 = cell.getI3();
            final int d3 = c3-_w3;
            final int d2 = c2-_w2;
            final int d1 = c1-_w1;
            final int b3 = max(d3,0);
            final int b2 = max(d2,0);
            final int b1 = max(d1,0);
            final int e3 = min(c3+_w3,n3-1);
            final int e2 = min(c2+_w2,n2-1);
            final int e1 = min(c1+_w1,n1-1);
            final int m3 = e3-b3+1;
            final int m2 = e2-b2+1;
            final int m1 = e1-b1+1;
            final float g11i = u11*fli;
            final float g12i = u12*fli;
            final float g13i = u13*fli;
            final float g22i = u22*fli;
            final float g23i = u23*fli;
            final float g33i = u33*fli;
            loop(m3,new LoopInt() {
            public void compute(int i3) {
              int k3 = i3+b3;
              for (int i2=0; i2<m2; ++i2) {
                int k2 = i2+b2;
                for (int i1=0; i1<m1; ++i1) {
                  int k1 = i1+b1;
                  int p1 = i1; if(d1<0) p1-=d1;
                  int p2 = i2; if(d2<0) p2-=d2;
                  int p3 = i3; if(d3<0) p3-=d3;
                  float sct = sci[p3][p2][p1];
                  g11[k3][k2][k1] += sct*g11i;
                  g12[k3][k2][k1] += sct*g12i;
                  g13[k3][k2][k1] += sct*g13i;
                  g22[k3][k2][k1] += sct*g22i;
                  g23[k3][k2][k1] += sct*g23i;
                  g33[k3][k2][k1] += sct*g33i;
              }}
            }});
          }
        }
      }
      sw.stop();
      trace("FaultScanner.scan: done");
      return solveEigenproblems(g11,g12,g13,g22,g23,g33);
  }


  public float[][][][][] voteScale(Sampling sp, Sampling st, float sigma) {
    int np = sp.getCount();
    int nt = st.getCount();
    int n3 = _w3*2+1;
    int n2 = _w2*2+1;
    int n1 = _w1*2+1;
    float[][][][][] scs = new float[np][nt][][][];
    float[][][] gx = new float[n3][n2][n1];
    float[][][] gs = new float[n3][n2][n1];
    gx[_w3][_w2][_w1] = 1.0f;
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    rgf.apply000(gx,gs);
    for (int ip=0; ip<np; ++ip) {
    for (int it=0; it<nt; ++it) {
      float fp = (float)sp.getValue(ip);
      float ft = (float)st.getValue(it);
      float[] us = faultNormalVectorFromStrikeAndDip(fp,ft);
      scs[ip][it] = voteScale(sigma,us[0],us[1],us[2],gs);
    }}
    return scs;
  }

  private float[][][] voteScale(
    float sigma, final float u1, final float u2, final float u3,
    final float[][][] gs) 
  {
    final float[][][] scs = new float[_w3*2+1][_w2*2+1][_w1*2+1];
    loop(-_w3,_w3+1,1,new LoopInt() {
    public void compute(int i3) {
      int k3 = i3+_w3;
      for (int i2=-_w2; i2<=_w2; ++i2) {
        int k2 = i2+_w2;
      for (int i1=-_w1; i1<=_w1; ++i1) {
        int k1 = i1+_w1;
        float r1 = i1;
        float r2 = i2;
        float r3 = i3;
        float rs = sqrt(r1*r1+r2*r2+r3*r3);
        if (rs==0) {
          scs[k3][k2][k1] = gs[k3][k2][k1];
        } else {
          r1 /= rs; r2 /= rs; r3 /= rs;
          float ur = u1*r1+u2*r2+u3*r3;
          ur = 1-ur*ur;
          ur *= ur;
          ur *= ur;
          ur *= ur;
          ur *= ur;
          ur *= ur;
          scs[k3][k2][k1] = gs[k3][k2][k1]*ur;
        }
      }}
    }});
    return scs;
  }


  private FaultCell[][][] getCellArray(
    int n1, int n2, int n3, Sampling sp, Sampling st, FaultCell[] fcs) 
  {
    int np = sp.getCount();
    int nt = st.getCount();
    int[][] ct = new int[np][nt];
    FaultCell[][][] fca = new FaultCell[np][nt][];
    FaultCell[][][] fct = new FaultCell[np][nt][n1*max(n2,n3)*5];
    for (FaultCell fci:fcs) {
      int ip = sp.indexOfNearest(fci.getFp());
      int it = st.indexOfNearest(fci.getFt());
      int ic = ct[ip][it];
      fct[ip][it][ic] = fci;
      ct[ip][it] += 1;
    }
    for (int ip=0; ip<np; ++ip) {
    for (int it=0; it<nt; ++it) {
      int nc = ct[ip][it];
      fca[ip][it] = new FaultCell[nc];
      for (int ic=0; ic<nc; ++ic) {
        fca[ip][it][ic] = fct[ip][it][ic];
      }
    }}
    fct = null;
    return fca;
  }

  public FaultCell[] randCells (
    int np, int sd, int n1, int n2, int n3, FaultCell[] fcs) {
    FaultCell[][][] fcg = new FaultCell[n3][n2][n1];
    for (FaultCell fc:fcs) {
      int i1 = fc.getI1();
      int i2 = fc.getI2();
      int i3 = fc.getI3();
      fcg[i3][i2][i1] = fc;
    }
    Random r = new Random(sd);
    int[][][] mark = zeroint(n1,n2,n3);
    FaultCell[] cells = new FaultCell[np];
    trace("np="+np);
    for (int ip=0; ip<np; ++ip) {
      boolean marked = false;
      while (!marked) {
        int i1 = r.nextInt(n1);
        int i2 = r.nextInt(n2);
        int i3 = r.nextInt(n3);
        boolean ok = true;
        int m = 1;
        for (int j3=max(0,i3-m); j3<min(n3,i3+m+1); ++j3){
        for (int j2=max(0,i2-m); j2<min(n2,i2+m+1); ++j2){
        for (int j1=max(0,i1-m); j1<min(n1,i1+m+1); ++j1){
          if (mark[j3][j2][j1]>0) ok = false;
        }}}
        FaultCell cell = fcg[i3][i2][i1];
        if (ok && cell!=null) {
          marked = true;
          mark[i3][i2][i1] = 1;
          cells[ip] = cell;
        }
      }
    }
    return cells;
  }

  public FaultCell[] getFaultCells(int n1, int n2, int n3, FaultCell[] fcs) {
    FaultCell[][][] fcg = new FaultCell[n3][n2][n1];
    ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    for (FaultCell fc:fcs) {
      int i1 = fc.getI1();
      int i2 = fc.getI2();
      int i3 = fc.getI3();
      fcg[i3][i2][i1] = fc;
    }

    for (int i3=0; i3<n3; i3+=1) {
    for (int i2=0; i2<n2; i2+=2) {
    for (int i1=0; i1<n1; i1+=2) {
      FaultCell fc = fcg[i3][i2][i1];
      if(fc!=null) cellList.add(fc);
    }}}
    return cellList.toArray(new FaultCell[0]);
  }
  public void getFlImage(FaultCell[] fcs, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    for (FaultCell fc:fcs) {
      int i1 = fc.getI1();
      int i2 = fc.getI2();
      int i3 = fc.getI3();
      i1 = min(i1,n1-1); i1 = max(i1,0);
      i2 = min(i2,n2-1); i2 = max(i2,0);
      i3 = min(i3,n3-1); i3 = max(i3,0);
      /*
      int m1 = fc.getM1();
      int m2 = fc.getM2();
      int m3 = fc.getM3();
      int p1 = fc.getP1();
      int p2 = fc.getP2();
      int p3 = fc.getP3();
      m1 = min(m1,n1-1); m1 = max(m1,0);
      m2 = min(m2,n2-1); m2 = max(m2,0);
      m3 = min(m3,n3-1); m3 = max(m3,0);
      p1 = min(p1,n1-1); p1 = max(p1,0);
      p2 = min(p2,n2-1); p2 = max(p2,0);
      p3 = min(p3,n3-1); p3 = max(p3,0);
      fl[m3][m2][m1] = fc.getFl();
      fl[p3][p2][p1] = fc.getFl();
      */
      fl[i3][i2][i1] = fc.getFl();
    }
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
            float x2 = i2+d2;
            float x3 = i3+d3;
            x2 = min(x2,n2-1); x2 = max(x2,0);
            x3 = min(x3,n3-1); x3 = max(x3,0);
            float tiii = faultDipFromNormalVector(u1ii,u2ii,u3ii);
            cell = new FaultCell(i1,i2+d2,i3+d3,fl,piii,tiii);
            cellList.add(cell);
          }
        }
      }
    }
    return cellList.toArray(new FaultCell[0]);
  }

  public FaultCell[] findCellsX (
    final float smin,
    final float[][][] ss, final float[][][] u1, 
    final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = ss.length;
    final int n2 = ss[0].length;
    final int n1 = ss[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final float[][][] sm = new float[n3][n2][n1];
    final float[][][] sp = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    final ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    loop(1,n3-1,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=1; i2<n2-1 ;++i2) {
    for (int i1=1; i1<n1-1 ;++i1) {
      float u1i = u1[i3][i2][i1]*1f;
      float u2i = u2[i3][i2][i1]*1f;
      float u3i = u3[i3][i2][i1]*1f;
      float x1m = i1-u1i;
      float x2m = i2-u2i;
      float x3m = i3-u3i;
      float x1p = i1+u1i;
      float x2p = i2+u2i;
      float x3p = i3+u3i;
      sm[i3][i2][i1] = si.interpolate(s1,s2,s3,ss,x1m,x2m,x3m);
      sp[i3][i2][i1] = si.interpolate(s1,s2,s3,ss,x1p,x2p,x3p);
    }}}});
    for (int i3=1; i3<n3-1; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float sxi = ss[i3][i2][i1];
      float sxm = sm[i3][i2][i1];
      float sxp = sp[i3][i2][i1];
      if (sxi>sxm && sxi>sxp && sxi>smin && abs(u1i)<0.9f) {
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

  public float[][][] thin(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] fl = new float[n3][n2][n1];
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    LocalOrientFilterP lof = new LocalOrientFilterP(8.0,2.0,2.0);
    lof.applyForNormal(fx,u1,u2,u3);
    FaultCell[] cells = findCells(0.1f,fx,u1,u2,u3);
    getFlImage(cells,fl);
    //return findRidges(fs,u1,u2,u3);
    return fl;
  }
  public float[][][] findRidges (
    final float[][][] sm, final float[][][] u1, 
    final float[][][] u2, final float[][][] u3) {
    final int n3 = sm.length;
    final int n2 = sm[0].length;
    final int n1 = sm[0][0].length;
    final float[][][] sp = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2 ;++i2) {
      for (int i1=0; i1<n1 ;++i1) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        if(abs(u1i)>0.5f) {continue;}
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
      }}
    }});
    return sp;
  }

  public float[][][][] solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33) {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    final float[][][] ss = new float[n3][n2][n1];
    final float[][][] cs = new float[n3][n2][n1];
    final float[][][] fp = new float[n3][n2][n1];
    final float[][][] ft = new float[n3][n2][n1];
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
        float u1i = (float)z[0][0];
        float u2i = (float)z[0][1];
        float u3i = (float)z[0][2];
        if (u2i==0.0f&&u3i==0f){continue;}
        if (u1i>0.0f) {
          u1i = -u1i;
          u2i = -u2i;
          u3i = -u3i;
        }
        ss[i3][i2][i1] = (eui-evi);
        cs[i3][i2][i1] = (evi-ewi)*(eui-evi);
        ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
        fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
      }}
    }});
    sub(ss,min(ss),ss);
    div(ss,max(ss),ss);
    sub(cs,min(cs),cs);
    div(cs,max(cs),cs);
    return new float[][][][]{ss,cs,fp,ft};
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
      fls[ic] = fc.getFl();
      xs[0][ic] = fc.getI1();
      xs[1][ic] = fc.getI2();
      xs[2][ic] = fc.getI3();
      us[0][ic] = fc.getW1();
      us[1][ic] = fc.getW2();
      us[2][ic] = fc.getW3();
    }
  }

  private float[][] setKdTreeNodes(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    ArrayList<Float> x1s = new ArrayList<Float>();
    ArrayList<Float> x2s = new ArrayList<Float>();
    ArrayList<Float> x3s = new ArrayList<Float>();
    ArrayList<Float> fxs = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fxi = fx[i3][i2][i1];
      if(fxi>0.0f) {
        x1s.add((float)i1);
        x2s.add((float)i2);
        x3s.add((float)i3);
        fxs.add((float)fxi);
      }
    }}}
    int np = x1s.size();
    float[][] xs = new float[4][np];
    for (int ip=0; ip<np; ++ip) {
      xs[0][ip] = x1s.get(ip);
      xs[1][ip] = x2s.get(ip);
      xs[2][ip] = x3s.get(ip);
      xs[3][ip] = fxs.get(ip);
    }
    return xs;
  }

  private int _w1 = 20;
  private int _w2 = 20;
  private int _w3 = 20;
  private float _sigma=10f;
}
