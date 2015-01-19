/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import static ifs.FaultGeometry.*;

/**
 * Regrid fault cells from known cells.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.19
 */

public class FaultCellRegrid {
  public FaultCellRegrid(int n1, int n2, int n3, FaultCell[] fcs) {
    _n1  = n1;
    _n2  = n2;
    _n3  = n3;
    _fcs = fcs;
    _nc  = fcs.length;
  }


  public FaultCell[] regridding() {
    int dd = 4;
    int dk = 20;
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for(FaultCell fc:_fcs){hsc.add(fc);}
    KdTree kt = setKdTreeCells(_fcs);
    int ct = 0;
    for (FaultCell fc:_fcs) {
      System.out.println("ct="+ct);
      ct++;
      float x1 = fc.x1;
      float x2 = fc.x2;
      float x3 = fc.x3;
      float fpc = fc.fp;
      int[] bs1 = bound(x1,dd,_n1);
      int[] bs2 = bound(x2,dd,_n2);
      int[] bs3 = bound(x3,dd,_n3);
      xmin[0] = x1-dk; xmax[0] = x1+dk;
      xmin[1] = x2-dk; xmax[1] = x2+dk;
      xmin[2] = x3-dk; xmax[2] = x3+dk;
      int[] id = kt.findInRange(xmin,xmax);
      int nd = id.length; if(nd<10){continue;}
      int n1 = bs1[1]-bs1[0]+1;
      int n2 = bs2[1]-bs2[0]+1;
      int n3 = bs3[1]-bs3[0]+1;
      float[][][] flt = new float[n3][n2][n1];
      float[][][] g11 = new float[n3][n2][n1];
      float[][][] g12 = new float[n3][n2][n1];
      float[][][] g13 = new float[n3][n2][n1];
      float[][][] g22 = new float[n3][n2][n1];
      float[][][] g23 = new float[n3][n2][n1];
      float[][][] g33 = new float[n3][n2][n1];
      for (int i3=0; i3<n3; ++i3) {
        int k3 = i3+bs3[0];
        for (int i2=0; i2<n2; ++i2) {
          int k2 = i2+bs2[0];
          for (int i1=0; i1<n1; ++i1) {
            int k1 = i1+bs1[0];
            float sv = 0.25f/(dk*dk); 
            float su = 0.25f/(dk*dk); 
            float sw = 0.25f/(2f*2f); 
            for (int ik=0; ik<nd; ++ik) {
              FaultCell fci = _fcs[id[ik]];
              if(abs(fci.fp-fpc)>8f){continue;}
              float x1i = fci.x1; 
              float x2i = fci.x2; 
              float x3i = fci.x3; 
              float dx1 = x1i-k1;
              float dx2 = x2i-k2;
              float dx3 = x3i-k3;
              float d11 = dx1*dx1;
              float d22 = dx2*dx2;
              float d33 = dx3*dx3;
              float d12 = dx1*dx2;
              float d13 = dx1*dx3;
              float d23 = dx2*dx3;
              float w11 = fci.w11;
              float w12 = fci.w12;
              float w13 = fci.w13;
              float w22 = fci.w22;
              float w23 = fci.w23;
              float w33 = fci.w33;
              float u11 = fci.u11;
              float u12 = fci.u12;
              float u13 = fci.u13;
              float u22 = fci.u22;
              float u23 = fci.u23;
              float u33 = fci.u33;
              float v11 = fci.v11;
              float v12 = fci.v12;
              float v13 = fci.v13;
              float v22 = fci.v22;
              float v23 = fci.v23;
              float v33 = fci.v33;
              float wd1 = w12*d12*2.0f;
              float wd2 = w13*d13*2.0f;
              float wd3 = w23*d23*2.0f;

              float ud1 = u12*d12*2.0f;
              float ud2 = u13*d13*2.0f;
              float ud3 = u23*d23*2.0f;

              float vd1 = v12*d12*2.0f;
              float vd2 = v13*d13*2.0f;
              float vd3 = v23*d23*2.0f;

              float wds = w11*d11+w22*d22+w33*d33;
              float uds = u11*d11+u22*d22+u33*d33;
              float vds = v11*d11+v22*d22+v33*d33;
              float gss = 0.0f;
              float wpi = pow(fci.fl,10f);
              gss += (wd1+wd2+wd3+wds)*sw;
              gss += (ud1+ud2+ud3+uds)*su;
              gss += (vd1+vd2+vd3+vds)*sv;
              float fli = exp(-gss)*wpi;
              flt[i3][i2][i1] += fli;
              g11[i3][i2][i1] += fli*w11;
              g12[i3][i2][i1] += fli*w12;
              g13[i3][i2][i1] += fli*w13;
              g22[i3][i2][i1] += fli*w22;
              g23[i3][i2][i1] += fli*w23;
              g33[i3][i2][i1] += fli*w33;
            }
          }
        }
      }
      FaultCell[] fsn = findCells(bs1,bs2,bs3,flt,g11,g12,g13,g22,g23,g33);
      KdTree ktt = setKdTreeCells(fsn);
      xmin[0] = x1-1; xmax[0] = x1+1;
      xmin[1] = x2-1; xmax[1] = x2+1;
      xmin[2] = x3-1; xmax[2] = x3+1;
      int[] idt = ktt.findInRange(xmin,xmax);
      int ndt = idt.length; if(ndt<1){continue;}
      for (int ik=0; ik<ndt; ++ik) {
        hsc.add(fsn[idt[ik]]);
      }
    }
    return getCells(hsc);
  }

  private FaultCell[] getCells(HashSet<FaultCell> hsc) {
    int ic = 0;
    int nc = hsc.size();
    FaultCell[] fcs = new FaultCell[nc];
    for (FaultCell fc:hsc) {
      fcs[ic] = fc;
      ic++;
    }
    return fcs;
  }

  private FaultCell[] findCells(
    int[] bs1, int[] bs2, int[] bs3, float[][][] fl, 
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33) 
  {
    div(fl,max(fl),fl);
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] fp = new float[n3][n2][n1];
    float[][][] ft = new float[n3][n2][n1];
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    faultImages(fp,ft,u1,u2,u3);
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1f,0.3f);
    float[][][] flp = new float[_n3][_n2][_n1];
    float[][][] fpp = new float[_n3][_n2][_n1];
    float[][][] ftp = new float[_n3][_n2][_n1];
    for (int i3=0; i3<n3; ++i3) {
      int k3 = i3+bs3[0];
      for (int i2=0; i2<n2; ++i2) {
        int k2 = i2+bs2[0];
        for (int i1=0; i1<n1; ++i1) {
          int k1 = i1+bs1[0];
          flp[k3][k2][k1] = fl[i3][i2][i1];
          fpp[k3][k2][k1] = fp[i3][i2][i1];
          ftp[k3][k2][k1] = ft[i3][i2][i1];
        }
      }
    }
    FaultCell[] fcs = fs.findCells(new float[][][][]{flp,fpp,ftp});
    return fcs;
  }


  private void faultImages(
    final float[][][] fp, final float[][][] ft,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = fp.length;
    final int n2 = fp[0].length;
    final int n1 = fp[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          if(u2i!=0.0f && u3i!=0.0f) {
            ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
            fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
          }
        }
      }
    }});
  }



  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n1 = g11[0][0].length;
    final int n2 = g11[0].length;
    final int n3 = g11.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        double[][] a = new double[3][3];
        double[][] z = new double[3][3];
        double[] e = new double[3];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            a[0][0] = g11[i3][i2][i1];
            a[0][1] = g12[i3][i2][i1];
            a[0][2] = g13[i3][i2][i1];
            a[1][0] = g12[i3][i2][i1];
            a[1][1] = g22[i3][i2][i1];
            a[1][2] = g23[i3][i2][i1];
            a[2][0] = g13[i3][i2][i1];
            a[2][1] = g23[i3][i2][i1];
            a[2][2] = g33[i3][i2][i1];
            Eigen.solveSymmetric33(a,z,e);
            float u1i = (float)z[0][0];
            float u2i = (float)z[0][1];
            float u3i = (float)z[0][2];
            if (u1i>0.0f) {
              u1i = -u1i;
              u2i = -u2i;
              u3i = -u3i;
            }
            u1[i3][i2][i1] = u1i;
            u2[i3][i2][i1] = u2i;
            u3[i3][i2][i1] = u3i;
          }
        }
      }
    });
  }



  private int[] bound(float x, float d, int n) {
    int im = round(x-d); if(im<0){im=0;}
    int ip = round(x+d); if(ip>=n){ip=n-1;}
    return new int[]{im,ip};
  }

  private KdTree setKdTreeCells(FaultCell[] fcs) {
    int nc = fcs.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xc[0][ic] = fcs[ic].x1;
      xc[1][ic] = fcs[ic].x2;
      xc[2][ic] = fcs[ic].x3;
      ic++;
    }
    return new KdTree(xc);
  }

  private FaultCell[] _fcs;
  private int _n1, _n2, _n3, _nc;
  private FaultCell[][][] fcg;
}

