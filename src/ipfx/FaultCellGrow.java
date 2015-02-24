/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipfx;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import static ipfx.FaultGeometry.*;

/**
 * Regrid fault cells from known cells.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.24
 */

public class FaultCellGrow {

  public FaultCellGrow(FaultCell[] fc, float[][][] fl) {
    _fc = fc; 
    _fl = fl;
    setKdTree();
    _n3 = fl.length;
    _n2 = fl[0].length;
    _n1 = fl[0][0].length;
  }

  public void setParameters(float dfp, float dft, float dnp) {
    _dfpmax = dfp;
    _dftmax = dft;
    _dnpmax = dnp;
  }

  public FaultCell[] combineCells(FaultCell[] fc1, FaultCell[] fc2) {
    int ic = 0;
    int nc = fc1.length+fc2.length;
    FaultCell[] cells = new FaultCell[nc];
    for (FaultCell fc:fc1) {
      cells[ic] = fc;
      ic++;
    }
    for (FaultCell fc:fc2) {
      cells[ic] = fc;
      ic++;
    }
    return cells;

  }

  public FaultCell[] applyForCells(FaultCell cell) {
    float ds = 20.f;
    float sigNor = 2.0f;
    final FaultCell[] cells = nabors(cell);
    if(cells==null) {return null;}
    final int nc = cells.length;
    final float su = 0.25f/(ds*ds);
    final float sv = 0.25f/(ds*ds);
    final float sw = 1.0f/(sigNor*sigNor);
    final int[] bs1 = new int[2];
    final int[] bs2 = new int[2];
    final int[] bs3 = new int[2];
    defineBox(cell,bs1,bs2,bs3);
    final int n1 = bs1[1]-bs1[0]+1;
    final int n2 = bs2[1]-bs2[0]+1;
    final int n3 = bs3[1]-bs3[0]+1;
    final float[][][] fls = new float[n3][n2][n1];
    final float[][][] g11 = new float[n3][n2][n1];
    final float[][][] g12 = new float[n3][n2][n1];
    final float[][][] g13 = new float[n3][n2][n1];
    final float[][][] g22 = new float[n3][n2][n1];
    final float[][][] g23 = new float[n3][n2][n1];
    final float[][][] g33 = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt(){
    public void compute(int k3) {
    int i3 = k3+bs3[0];
    for (int k2=0; k2<n2; ++k2) {
    int i2 = k2+bs2[0];
    for (int k1=0; k1<n1; ++k1) {
    int i1 = k1+bs1[0];
      float wps = 0.0f;
      for (int ic=0; ic<nc; ++ic) {
        FaultCell fci = cells[ic];
        float dx1 = fci.x1-i1;
        float dx2 = fci.x2-i2;
        float dx3 = fci.x3-i3;
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
        float wpi = fci.fl;//pow(fci.fl,10f);
        wps += wpi;
        gss += (wd1+wd2+wd3+wds)*sw;
        gss += (ud1+ud2+ud3+uds)*su;
        gss += (vd1+vd2+vd3+vds)*sv;
        float fli = exp(-gss)*wpi;
        fls[k3][k2][k1] += fli;
        g11[k3][k2][k1] += fli*w11;
        g12[k3][k2][k1] += fli*w12;
        g13[k3][k2][k1] += fli*w13;
        g22[k3][k2][k1] += fli*w22;
        g23[k3][k2][k1] += fli*w23;
        g33[k3][k2][k1] += fli*w33;
      }
      fls[k3][k2][k1] /= wps;
    }}}});
    return findCells(bs1,bs2,bs3,fls,g11,g12,g13,g22,g23,g33);
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
    FaultCell[] fcs = cells(new float[][][][]{fl,fp,ft});
    int nc = fcs.length; if(nc<1){return null;}
    for (int ic=0; ic<nc; ++ic) {
      float x1i = fcs[ic].x1+bs1[0];
      float x2i = fcs[ic].x2+bs2[0];
      float x3i = fcs[ic].x3+bs3[0];
      int i1 = round(x1i);
      int i2 = round(x2i);
      int i3 = round(x3i);
      if(i1<0) {continue;}
      if(i2<0) {continue;}
      if(i3<0) {continue;}
      if(i1>=_n1) {continue;}
      if(i2>=_n2) {continue;}
      if(i3>=_n3) {continue;}
      float fpi = fcs[ic].fp;
      float fti = fcs[ic].ft;
      float fli = bilinearInterp(x1i,x2i,x3i);
      fcs[ic] = new FaultCell(x1i,x2i,x3i,fli,fpi,fti);
    }
    return fcs;
  }

  private float bilinearInterp(float x1, float x2, float x3) {
    int i1i = round(x1);
    int i2i = round(x2);
    int i3i = round(x3);
    int i2p=i2i,i2m=i2i;
    int i3p=i3i,i3m=i3i;
    if(i2i>x2) {i2m -= 1;}
    else       {i2p += 1;}
    if(i3i>x3) {i3m -= 1;}
    else       {i3p += 1;}
    if(i2m<0||i2p>=_n2) {return _fl[i3i][i2i][i1i];}
    if(i3m<0||i3p>=_n3) {return _fl[i3i][i2i][i1i];}
    x2 -= i2i; x3 -= i3i;
    float fla = _fl[i3m][i2m][i1i];
    float flb = _fl[i3m][i2p][i1i];
    float flc = _fl[i3p][i2m][i1i];
    float fld = _fl[i3p][i2p][i1i];
    float fli = fla*(1f-x2)*(1f-x3)+flb*(1f-x3)*x2+flc*(1f-x2)*x3+fld*x2*x3;
    return fli;
  }

  // Uses fault images to find cells, oriented points located on ridges.
  private FaultCell[] cells(float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

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
        float[] pii = p[i3 ][i2 ];
        float[] tii = t[i3 ][i2 ];
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
          float piii = pii[i1 ];
          float tiii = tii[i1 ];

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
            if (fr>=_fllo) {
              fl += fr;
              d2 += dr;
              nr += 1;
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
            if (fr>=_fllo) {
              fl += fr;
              d2 += dr;
              d3 -= dr;
              nr += 1;
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
            if (fr>=_fllo) {
              fl += fr;
              d3 += dr;
              nr += 1;
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
            if (fr>=_fllo) {
              fl += fr;
              d2 += dr;
              d3 += dr;
              nr += 1;
            }
          }

          // If at least one ridge, construct a cell and add to list.
          if (nr>0) {
            fl /= nr;
            d2 /= nr;
            d3 /= nr;
            cell = new FaultCell(i1,i2+d2,i3+d3,fl,piii,tiii);
            cellList.add(cell);
          }
        }
      }
    }
    return cellList.toArray(new FaultCell[0]);
  }


  private void solveEigenproblems(
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33,
    float[][][] u1,  float[][][] u2,  float[][][] u3) 
  {
    int n3 = g11.length;
    int n2 = g11[0].length;
    int n1 = g11[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      double[] e = new double[3];
      double[][] a = new double[3][3];
      double[][] z = new double[3][3];
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
  }

  private void faultImages(
    final float[][][] fp, final float[][][] ft,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = fp.length;
    final int n2 = fp[0].length;
    final int n1 = fp[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      if(u2i!=0.0f && u3i!=0.0f) {
        ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
        fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
      }
    }}}
  }

  private void defineBox(FaultCell fc, int[] bs1, int[] bs2, int[] bs3) {
    float d = 5f;
    int i1m = round(fc.x1-d);if(i1m<0){i1m=0;}
    int i2m = round(fc.x2-d);if(i2m<0){i2m=0;}
    int i3m = round(fc.x3-d);if(i3m<0){i3m=0;}
    int i1p = round(fc.x1+d);if(i1p>=_n1){i1p=_n1-1;}
    int i2p = round(fc.x2+d);if(i2p>=_n2){i2p=_n2-1;}
    int i3p = round(fc.x3+d);if(i3p>=_n3){i3p=_n3-1;}
    bs1[0] = i1m; bs1[1] = i1p;
    bs2[0] = i2m; bs2[1] = i2p;
    bs3[0] = i3m; bs3[1] = i3p;
  }

  public FaultCell[] nabors(FaultCell cell) {
    float dv = 10f;
    //float dh = 65f;
    float dh = 30f;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float v2 = cell.v2;
    float v3 = cell.v3;
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    xmin[0] = x1-dv; xmax[0] = x1+dv;
    xmin[1] = x2-dh; xmax[1] = x2+dh;
    xmin[2] = x3-dh; xmax[2] = x3+dh;
    int[] id = _kt.findInRange(xmin,xmax);
    int nd = id.length;
    if(nd<1) {return null;}
    ArrayList<Float> dsL = new ArrayList<Float>();
    ArrayList<Float> dsR = new ArrayList<Float>();
    ArrayList<FaultCell> fcL = new ArrayList<FaultCell>();
    ArrayList<FaultCell> fcR = new ArrayList<FaultCell>();
    for (int ik=0; ik<nd; ++ik) {
      int ic = id[ik];
      FaultCell fci = _fc[ic];
      if(canBeNabors(fci,cell)){
        float d1 = fci.x1-x1;
        float d2 = fci.x2-x2;
        float d3 = fci.x3-x3;
        float dd = d2*v2+d3*v3;
        float ds = d1*d1+d2*d2+d3*d3;
        if(abs(d1)>1f) {ds *= abs(d1);}
        if(dd<=0f) {dsL.add(ds);fcL.add(fci);} 
        else       {dsR.add(ds);fcR.add(fci);}
      }
    }
    int nb = 20;
    int nl = dsL.size();
    int nr = dsR.size();
    if(nl<10||nr<10) {return null;}
    ArrayList<FaultCell> nbs = new ArrayList<FaultCell>();
    if(nl<nb) {
      for (FaultCell fc:fcL) {nbs.add(fc);}
    } else {
      int il = 0;
      int[] idl = new int[nl];
      float[] dsl = new float[nl];
      FaultCell[] fcl = fcL.toArray(new FaultCell[0]);
      for (float dsi:dsL) {dsl[il]=dsi; idl[il]=il; il++;}
      quickIndexSort(dsl,idl);
      for (int ik=0; ik<nb; ++ik) {
        int ic = idl[ik];
        nbs.add(fcl[ic]);
      }
    }
    if(nr<nb) {
      for (FaultCell fc:fcR) {nbs.add(fc);}
    } else {
      int ir = 0;
      int[] idr = new int[nr];
      float[] dsr = new float[nr];
      FaultCell[] fcr = fcR.toArray(new FaultCell[0]);
      for (float dsi:dsR) {dsr[ir]=dsi; idr[ir]=ir; ir++;}
      quickIndexSort(dsr,idr);
      for (int ik=0; ik<nb; ++ik) {
        int ic = idr[ik];
        nbs.add(fcr[ic]);
      }
    }
    return nbs.toArray(new FaultCell[0]);
  }


  private boolean canBeNabors(FaultCell ci, FaultCell cn) {
    boolean can = true;
    if (absDeltaFp(ci,cn)>_dfpmax) {
      can = false;
    } else if (absDeltaFt(ci,cn)>_dftmax) {
      can = false;
    } else if (maxDistanceToPlane(ci,cn)>_dnpmax) {
      can = false;
    }
    return can;
  }
  private static float absDeltaFp(FaultCell ca, FaultCell cb) {
    float del = ca.fp-cb.fp;
    return min(abs(del),abs(del+360.0f),abs(del-360.0f));
  }
  private static float absDeltaFt(FaultCell ca, FaultCell cb) {
    return abs(ca.ft-cb.ft);
  }
  private static float maxDistanceToPlane(FaultCell ca, FaultCell cb) {
    float aw1 = ca.w1, aw2 = ca.w2, aw3 = ca.w3;
    float ax1 = ca.x1, ax2 = ca.x2, ax3 = ca.x3;
    float bw1 = cb.w1, bw2 = cb.w2, bw3 = cb.w3;
    float bx1 = cb.x1, bx2 = cb.x2, bx3 = cb.x3;
    float dx1 = ax1-bx1;
    float dx2 = ax2-bx2;
    float dx3 = ax3-bx3;
    float dab = abs(aw1*dx1+aw2*dx2+aw3*dx3);
    float dba = abs(bw1*dx1+bw2*dx2+bw3*dx3);
    return max(dab,dba);
  }

  private void setKdTree() {
    int nc = _fc.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xc[0][ic] = _fc[ic].x1;
      xc[1][ic] = _fc[ic].x2;
      xc[2][ic] = _fc[ic].x3;
    }
    _kt =  new KdTree(xc);
  }

  
  private KdTree _kt;
  private FaultCell[] _fc;
  private float[][][] _fl;
  private int _n1, _n2, _n3;
  private float _fllo = 0.2f;
  //private float _dfpmax=15f; // max difference between strikes of nabors
  //private float _dfpmax=13f; // max difference between strikes of nabors
  private float _dfpmax=20f; // max difference between strikes of nabors
  private float _dftmax=10f; // max difference between dips of nabors
  //private float _dnpmax=4.f; // max distance to planes of nabors
  private float _dnpmax=5.f; // max distance to planes of nabors
}

