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
  public void setCells(FaultCell[] fcs) {
    _fcs = fcs;
    _nc  = fcs.length;
  }

  public void resetCells(FaultCell[] fcs) {
    int nc = fcs.length;
    for (int ic=0; ic<nc; ++ic) {
      FaultCell fc = fcs[ic];
      float x1 = fc.x1;
      float x2 = fc.x2;
      float x3 = fc.x3;
      float fl = fc.fl;
      float fp = fc.fp;
      float ft = fc.ft;
      fcs[ic] = new FaultCell(x1,x2,x3,fl,fp,ft);
    }
  }

  public FaultCell[] addNabors(final FaultCell[] fcs, final float[][][] fl) {
    _fcs1 = fcs;
    _kt = setKdTreeCells(_fcs);
    _kt1 = setKdTreeCells(_fcs1);
    Stopwatch sw = new Stopwatch();
    sw.start();
    final FaultReskin fr = new FaultReskin();
    fr.checkNabors(fcs);
    double timeUsed = sw.time();
    System.out.println("time for check cells: "+timeUsed);
    final FaultCell[][][] cg = new FaultCell[_n3][_n2][_n1];
    final int nc = fcs.length;
    Parallel.loop(nc,new Parallel.LoopInt(){
    public void compute(int ic) {
      FaultCell fc = fcs[ic];
      int c1 = fc.i1;
      int c2 = fc.i2;
      int c3 = fc.i3;
      cg[c3][c2][c1] = fc;
      FaultCell[] fcsn = null;
      if(fc.needInterp){fcsn=interpNabors(fc);}
      if(fcsn!=null){
        //FaultCell[] fcsn = nearestNabors(fc,cells);
        if(fc.ca==null){
          fr.addNaborAbove(fc,fcsn); 
          FaultCell ca = fc.ca;
          if(ca!=null){
          int i1 = ca.i1;
          int i2 = ca.i2;
          int i3 = ca.i3;
          ca.fl = fl[i3][i2][i1];
          cg[i3][i2][i1] = ca;
          }
        }
        if(fc.cb==null){
          fr.addNaborBelow(fc,fcsn); 
          FaultCell cb = fc.cb;
          if(cb!=null){
            int i1 = cb.i1;
            int i2 = cb.i2;
            int i3 = cb.i3;
            cb.fl = fl[i3][i2][i1];
            cg[i3][i2][i1] = cb;
          }
        }
        if(fc.cl==null){
          fr.addNaborLeft(fc,fcsn);  
          FaultCell cl = fc.cl;
          if(cl!=null){
          int i1 = cl.i1;
          int i2 = cl.i2;
          int i3 = cl.i3;
          cl.fl = fl[i3][i2][i1];
          cg[i3][i2][i1] = cl;
          }
        }
        if(fc.cr==null){
          fr.addNaborRight(fc,fcsn); 
          FaultCell cr = fc.cr;
          if(cr!=null){
          int i1 = cr.i1;
          int i2 = cr.i2;
          int i3 = cr.i3;
          cr.fl = fl[i3][i2][i1];
          cg[i3][i2][i1] = cr;
          }
        }
      }
      fc.setInterp(false);
      cg[c3][c2][c1] = fc;
    }});
    timeUsed = sw.time();
    System.out.println("time for interp cells: "+timeUsed);
    sw.stop();
    return getCells(cg);
  }

  private boolean canBeNabor(int dk, FaultCell fc) {
    float x1 = fc.x1;
    float x2 = fc.x2;
    float x3 = fc.x3;
    float fp = fc.fp;
    float[] xmin = new float[]{x1-dk,x2-dk,x3-dk};
    float[] xmax = new float[]{x1+dk,x2+dk,x3+dk};
    int[] id = _kt.findInRange(xmin,xmax);
    int nd = id.length; if(nd<10){return false;}
    float ns = 0f, ds=0f;
    for (int ik=0; ik<nd; ++ik) {
      int ic = id[ik];
      FaultCell fci = _fcs[ic];
      float fpc = fci.fp;
      if(abs(fpc-fp)>8f){continue;}
      ns += 1f;
      float w1c = fci.w1;
      float w2c = fci.w2;
      float w3c = fci.w3;
      float d1c = fci.x1-x1;
      float d2c = fci.x2-x2;
      float d3c = fci.x3-x3;
      ds += abs(d1c*w1c+d2c*w2c+d3c*w3c);
    }
    if(ns/nd<0.2f){return false;}
    if(ds/ns>0.5f){return false;}
    return true;
  }

  private void checkNabors(final FaultCell[] fcs) {
    final int nc = fcs.length;
    final FaultReskin fr = new FaultReskin();
    Parallel.loop(nc, new Parallel.LoopInt() {
    public void compute(int ic) {
      if(fcs[ic].needInterp){
        if(fcs[ic].ca==null){fr.addNaborAbove(fcs[ic],fcs);}
        if(fcs[ic].cb==null){fr.addNaborBelow(fcs[ic],fcs);}
        if(fcs[ic].cl==null){fr.addNaborLeft(fcs[ic],fcs);}
        if(fcs[ic].cr==null){fr.addNaborRight(fcs[ic],fcs);}
        if(fcs[ic].ca!=null&&fcs[ic].cb!=null&&
         fcs[ic].cl!=null&&fcs[ic].cr!=null) {
         fcs[ic].setInterp(false);
        }
      }
    }});
  }

  public FaultCell[] smoothCells(final float[][][] fl) {
    _kt = setKdTreeCells(_fcs);
    final FaultCell[][][] cg = new FaultCell[_n3][_n2][_n1];
    Parallel.loop(_nc,new Parallel.LoopInt(){
    public void compute(int ic) {
      FaultCell fc = _fcs[ic];
      FaultCell[] fcs = interpNabors(fc);
      if(fcs!=null){
        FaultCell fn = nearestCell(fc,fcs);
        if(fn!=null) {
          int i1 = fn.i1;
          int i2 = fn.i2;
          int i3 = fn.i3;
          fn.fl = fl[i3][i2][i1];
          cg[i3][i2][i1] = fn;
        }
      }
    }});
    return getCells(cg);
  }

  private FaultCell[] getCells(FaultCell[][][] cg) {
    int n3 = cg.length;
    int n2 = cg[0].length;
    int n1 = cg[0][0].length;
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (int i3=0; i3<n3; ++i3){
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0; i1<n1; ++i1){
      FaultCell fc = cg[i3][i2][i1];
      if(fc!=null){hsc.add(fc);}
    }}}
    return getCells(hsc);
  }

  private FaultCell nearestCell(FaultCell fc, FaultCell[] fcs) {
    float x1 = fc.x1;
    float x2 = fc.x2;
    float x3 = fc.x3;
    FaultCell fa = null;
    float dm = Float.MAX_VALUE;
    for (FaultCell fci:fcs) {
      float dx1 = fci.x1-x1;
      float dx2 = fci.x2-x2;
      float dx3 = fci.x3-x3;
      float dxs = dx1*dx1+dx2*dx2+dx3*dx3;
      if(dxs<dm) {dm=dxs; fa=fci;}
    }
    return fa;
  }

  private FaultCell[] nearestNabors(FaultCell fc, FaultCell[] fcs) {
    float x1 = fc.x1;
    float x2 = fc.x2;
    float x3 = fc.x3;
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (int i3=-1; i3<=1; ++i3) {
    for (int i2=-1; i2<=1; ++i2) {
    for (int i1=-1; i1<=1; ++i1) {
      float x1i = x1+i1;
      float x2i = x2+i2;
      float x3i = x3+i3;
      FaultCell fn = closeCell(x1i,x2i,x3i,fcs);
      if(fn!=null){hsc.add(fn);}
    }}}
    return getCells(hsc);
  }

  private FaultCell closeCell(float x1, float x2, float x3, FaultCell[] fcs) {
    FaultCell fn = null;
    float dm = Float.MAX_VALUE;
    for (FaultCell fci:fcs) {
      float dx1 = fci.x1-x1;
      float dx2 = fci.x2-x2;
      float dx3 = fci.x3-x3;
      float dxs = dx1*dx1+dx2*dx2+dx3*dx3;
      if(dxs<dm) {dm=dxs; fn=fci;}
    }
    return fn;
  }


  private FaultCell[] interpNabors(FaultCell fc) {
    int dd = 5;
    int dk = 10;
    final float x1 = fc.x1;
    final float x2 = fc.x2;
    final float x3 = fc.x3;
    final float fp = fc.fp;
    final float sv = 0.25f/(dk*dk); 
    final float su = 0.25f/(dk*dk); 
    final float sw = 0.25f/(2f*2f); 
    final int[] bs1 = bound(x1,dd,_n1);
    final int[] bs2 = bound(x2,dd,_n2);
    final int[] bs3 = bound(x3,dd,_n3);
    final int n1 = bs1[1]-bs1[0]+1;
    final int n2 = bs2[1]-bs2[0]+1;
    final int n3 = bs3[1]-bs3[0]+1;
    final int[] ds = new int[]{dk,dk,dk};
    final int[] ns = new int[]{_n1,_n2,_n3};
    final float[][][] flt = new float[n3][n2][n1];
    final float[][][] g11 = new float[n3][n2][n1];
    final float[][][] g12 = new float[n3][n2][n1];
    final float[][][] g13 = new float[n3][n2][n1];
    final float[][][] g22 = new float[n3][n2][n1];
    final float[][][] g23 = new float[n3][n2][n1];
    final float[][][] g33 = new float[n3][n2][n1];
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    int[] is = new int[]{round(x1),round(x2),round(x3)};
    getRange(new int[]{12,12,12},is,ns,xmin,xmax);
    int[] cd = _kt.findInRange(xmin,xmax);
    int nd = cd.length; if(nd<10){return null;}
    float nc=0f, dc=0f;
    for (int ik=0; ik<nd; ++ik) {
      FaultCell fci = _fcs[cd[ik]];
      float dp = fci.fp-fp;
      dp=min(abs(dp),abs(dp+360f),abs(dp-360f));
      if(dp>10f){continue;}
      nc += 1f;
      float w1i = fci.w1; 
      float w2i = fci.w2; 
      float w3i = fci.w3; 
      float dx1 = fci.x1-x1;
      float dx2 = fci.x2-x2;
      float dx3 = fci.x3-x3;
      dc += abs(w1i*dx1+w2i*dx2+w3i*dx3);
    }
    if(nc/nd<0.1f){return null;}
    if(dc/nc>1.0f){return null;}

    for (int k3=0; k3<n3; ++k3) {
      int i3 = k3+bs3[0];
      for (int k2=0; k2<n2; ++k2) {
        int i2 = k2+bs2[0];
        for (int k1=0; k1<n1; ++k1) {
          float wps = 0f;
          int i1 = k1+bs1[0];
          int[] xs = new int[]{i1,i2,i3};
          getRange(ds,xs,ns,xmin,xmax);
          int[]  id = _kt1.findInRange(xmin,xmax);
          nd = id.length; if(nd<10){continue;}
          for (int ik=0; ik<nd; ++ik) {
            FaultCell fci = _fcs1[id[ik]];
            float dp = fci.fp-fp;
            dp=min(abs(dp),abs(dp+360f),abs(dp-360f));
            if(dp>10f){continue;}
            float x1i = fci.x1; 
            float x2i = fci.x2; 
            float x3i = fci.x3; 
            float dx1 = x1i-i1;
            float dx2 = x2i-i2;
            float dx3 = x3i-i3;
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
            wps += wpi;
            gss += (wd1+wd2+wd3+wds)*sw;
            gss += (ud1+ud2+ud3+uds)*su;
            gss += (vd1+vd2+vd3+vds)*sv;
            float fli = exp(-gss)*wpi;
            flt[k3][k2][k1] += fli;
            g11[k3][k2][k1] += fli*w11;
            g12[k3][k2][k1] += fli*w12;
            g13[k3][k2][k1] += fli*w13;
            g22[k3][k2][k1] += fli*w22;
            g23[k3][k2][k1] += fli*w23;
            g33[k3][k2][k1] += fli*w33;
          }
          if(wps!=0f) flt[k3][k2][k1] /= wps;
        }
      }
    }
    //}});
    return findCells(bs1,bs2,bs3,flt,g11,g12,g13,g22,g23,g33);
  }

  private static void getRange(int[] ds, int[] is, int[] ns, 
    float[] xmin, float[] xmax) {
    int i1m = is[0]-ds[0]; if(i1m<0){i1m=0;}
    int i2m = is[1]-ds[1]; if(i2m<0){i2m=0;}
    int i3m = is[2]-ds[2]; if(i3m<0){i3m=0;}
    int i1p = is[0]+ds[0]; if(i1p>=ns[0]){i1p=ns[0]-1;}
    int i2p = is[1]+ds[1]; if(i2p>=ns[1]){i2p=ns[1]-1;}
    int i3p = is[2]+ds[2]; if(i3p>=ns[2]){i3p=ns[2]-1;}
    xmin[0] = i1m; xmin[1] = i2m; xmin[2] = i3m;
    xmax[0] = i1p; xmax[1] = i2p; xmax[2] = i3p;
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
    int dt = 10;
    float[][][] flt = new float[n3+20][n2+20][n1+20];
    float[][][] fpt = new float[n3+20][n2+20][n1+20];
    float[][][] ftt = new float[n3+20][n2+20][n1+20];
    for (int i3=0; i3<n3; ++i3) {
      int k3 = i3+dt;
    for (int i2=0; i2<n2; ++i2) {
      int k2 = i2+dt;
    for (int i1=0; i1<n1; ++i1) {
      int k1 = i1+dt;
      flt[k3][k2][k1] = fl[i3][i2][i1];
      fpt[k3][k2][k1] = fp[i3][i2][i1];
      ftt[k3][k2][k1] = ft[i3][i2][i1];
    }}}
    FaultCell[] fcs = fs.findCells(new float[][][][]{flt,fpt,ftt});
    int nc = fcs.length; if(nc<1){return null;}
    for (int ic=0; ic<nc; ++ic) {
      float fli = fcs[ic].fl;
      float fpi = fcs[ic].fp;
      float fti = fcs[ic].ft;
      float x1i = fcs[ic].x1+bs1[0]-dt;
      float x2i = fcs[ic].x2+bs2[0]-dt;
      float x3i = fcs[ic].x3+bs3[0]-dt;
      int i1 = round(x1i);
      int i2 = round(x2i);
      int i3 = round(x3i);
      if(i1<0) x1i+=1;
      if(i2<0) x2i+=1;
      if(i3<0) x3i+=1;
      if(i1>=_n1) x1i-=1;
      if(i2>=_n2) x2i-=1;
      if(i3>=_n3) x3i-=1;
      fcs[ic] = new FaultCell(x1i,x2i,x3i,fli,fpi,fti);
    }
    return fcs;
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



  private int[] bound(float x, float d, int n) {
    int im = round(x-d); if(im<0) {im=0;}
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

  private KdTree _kt,_kt1;
  private FaultCell[] _fcs,_fcs1;
  private int _n1, _n2, _n3, _nc;
  private FaultCell[][][] fcg;
}

