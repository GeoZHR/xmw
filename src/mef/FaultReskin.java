package mef;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

import slt.*;
import static mef.FaultGeometry.*;


/**
 * Construct smooth and single-valued fault surface using 
 * the screen poisson method.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.02.29
 */

public class FaultReskin {

 
 public float[][] faultSurfer(int n1, int n2, int n3, FaultSkin skin) {
   float[][] sf = new float[n3][n2];
   return sf;
 }

 public FaultCell[] getCells(FaultSkin sk) {
   ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
   for (FaultCell cell:sk) {
     int i1 = cell.i1;
     int i3 = cell.i3;
     cell.ca = null;
     cell.cb = null;
     cell.cl = null;
     cell.cr = null;
     cell.skin = null;
     if (i1>290&&i3>260 &&i3<294) {continue;}
     fcl.add(cell);
   }
   return fcl.toArray(new FaultCell[0]);
 }


 public FaultSkin[] reskin(int m1, int m2, FaultSkin skin) {
   ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
   for (FaultCell cell:skin) {
     if(cell.i1<m1) {continue;}
     if(cell.x2<m2) {continue;}
     cell.skin=null;
     fcl.add(cell);
   }
   FaultSkinner fs = new FaultSkinner();
   fs.setGrowLikelihoods(0.1f,0.6f);
   fs.setMinSkinSize(5000);
   FaultCell[] cells = fcl.toArray(new FaultCell[0]);
   return fs.findSkins(cells);
 }

 public FaultSkin[] reskinJake(int n1, int n2, int n3, FaultSkin[] skins) {
   FaultCell[][][] fca = new FaultCell[n3][n2][n1];
   ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
   for (FaultCell cell:skins[1]) {
     cell.ca = null;
     cell.cb = null;
     cell.cl = null;
     cell.cr = null;
     cell.skin = null;
     int i1 = cell.i1;
     int i2 = cell.i2;
     int i3 = cell.i3;
     int b2 = max(0,i2-2);
     int b3 = max(0,i3-2);
     int e2 = min(n2-1,i2+2);
     int e3 = min(n3-1,i3+2);
     for (int k3=b3; k3<=e3; ++k3) {
     for (int k2=b2; k2<=e2; ++k2) {
       fca[k3][k2][i1] = cell;
     }}
     fcl.add(cell);
   }
   int nk = skins.length;
   for (int k=0; k<nk; ++k) {
     if(k==1) {continue;}
     for (FaultCell cell:skins[k]) {
       cell.ca = null;
       cell.cb = null;
       cell.cl = null;
       cell.cr = null;
       cell.skin = null;
       int i1 = cell.i1;
       int i2 = cell.i2;
       int i3 = cell.i3;
       if (fca[i3][i2][i1]==null) {
         fcl.add(cell);
       }
     }
   }
   FaultSkinner  fs = new FaultSkinner();
   fs.setGrowLikelihoods(0.005f,0.6f);
   fs.setMaxDeltaStrike(10);
   fs.setMaxPlanarDistance(0.2f);
   fs.setMinSkinSize(10000);
   return fs.findSkins(fcl.toArray(new FaultCell[0]));
 }

 public float[][][][] rescan(
    int n1, int n2, int n3, 
    Sampling sp, Sampling st, FaultCell[] cells) 
 {
   _fcs = cells;
   float[][][][][] gws = gaussWeights(10.f,2.f,40,40,40,sp,st);
   KdTree kt = setStrikeKdTree();
   float[][][] fl = new float[n3][n2][n1];
   float[][][] fp = new float[n3][n2][n1];
   float[][][] ft = new float[n3][n2][n1];
   float dp1 = 20f;
   float dp2 = 20f;
   for (float fpi=0; fpi<=360; fpi+=dp1) {
     float fpm = fpi-dp2;
     float fpp = fpi+dp2;
     FaultCell[] fcs = cellsInStrikeWd(fpm,fpp,kt);
     float[][][][] flpt = new float[3][n3][n2][n1];
     faultImagesFromCells(sp,st,fcs,gws,flpt);
     float[][][] flt = flpt[0];
     float[][][] fpt = flpt[1];
     float[][][] ftt = flpt[2];
     for (int i3=0; i3<n3; ++i3) {
     for (int i2=0; i2<n2; ++i2) {
     for (int i1=0; i1<n1; ++i1) {
       float fli = flt[i3][i2][i1];
       if(fli>fl[i3][i2][i1]) {
         fl[i3][i2][i1] = fli;
         fp[i3][i2][i1] = fpt[i3][i2][i1];
         ft[i3][i2][i1] = ftt[i3][i2][i1];
       }
     }}}
   }
   return new float[][][][]{fl,fp,ft};
 }

 public FaultSkin[] faultSkinsFromCellsJake(
   int n1, int n2, int n3, FaultCell[] cells) {
   setCells(cells);
   float[][][] fls = new float[_n3][_n2][_n1];
   float[][][] fps = new float[_n3][_n2][_n1];
   float[][][] fts = new float[_n3][_n2][_n1];
   float[][][] g11 = new float[_n3][_n2][_n1];
   float[][][] g12 = new float[_n3][_n2][_n1];
   float[][][] g13 = new float[_n3][_n2][_n1];
   float[][][] g22 = new float[_n3][_n2][_n1];
   float[][][] g23 = new float[_n3][_n2][_n1];
   float[][][] g33 = new float[_n3][_n2][_n1];
   for (FaultCell cell:cells) {
     float fli = cell.fl;
     int i1 = cell.i1-_j1;
     int i2 = cell.i2-_j2;
     int i3 = cell.i3-_j3;
     fls[i3][i2][i1] = fli;
     g11[i3][i2][i1] = cell.w11*fli;
     g12[i3][i2][i1] = cell.w12*fli;
     g13[i3][i2][i1] = cell.w13*fli;
     g22[i3][i2][i1] = cell.w22*fli;
     g23[i3][i2][i1] = cell.w23*fli;
     g33[i3][i2][i1] = cell.w33*fli;
   }
   System.out.println("assignments done...");
   RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(2);
   rgf1.apply000(fls,fls);
   System.out.println("fl smoothing done...");
   RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(10);
   float[][][][] gs = {g11,g22,g33,g12,g13,g23};
   for (float[][][] g:gs) {
     rgf2.apply000(g,g);
   }
   System.out.println("gaussian smoothing done...");
   float[][][] w1 = new float[_n3][_n2][_n1];
   float[][][] w2 = new float[_n3][_n2][_n1];
   float[][][] u1 = new float[_n3][_n2][_n1];
   float[][][] u2 = new float[_n3][_n2][_n1];
   float[][][] u3 = new float[_n3][_n2][_n1];
   solveEigenproblems(g11,g12,g13,g22,g23,g33,w1,w2,u1,u2,u3);
   // Compute u1 such that u3 > 0.
   for (int i3=0; i3<_n3; ++i3) {
     for (int i2=0; i2<_n2; ++i2) {
       for (int i1=0; i1<_n1; ++i1) {
         float u1i = u2[i3][i2][i1];
         float u2i = u2[i3][i2][i1];
         float u3i = u3[i3][i2][i1];
         float u1s = 1.0f-u2i*u2i-u3i*u3i;
         u1i = (u1s>0.0f)?sqrt(u1s):0.0f;
         if (u3i<0.0f) {
           u1i = -u1i;
           u2i = -u2i;
         }
         u1[i3][i2][i1] = u1i;
         u2[i3][i2][i1] = u2i;
       }
     }
   }
   System.out.println("eigentensors done...");
   float[][][] eu = fillfloat(0.01f,_n1,_n2,_n3);
   float[][][] ev = fillfloat(1.00f,_n1,_n2,_n3);
   float[][][] ew = fillfloat(1.00f,_n1,_n2,_n3);
   EigenTensors3 et = new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
   LocalSmoothingFilter lsf = new LocalSmoothingFilter();
   lsf.apply(et,400,fls,fls);
   computeStrikeDip(fls,fps,fts);
   System.out.println("structure-oriented smoothing done...");
   FaultSkinner  fs = new FaultSkinner();
   fs.setGrowLikelihoods(0.005f,0.6f);
   fs.setMaxDeltaStrike(10);
   fs.setMaxPlanarDistance(0.2f);
   fs.setMinSkinSize(10000);
   div(fls,max(fls),fls);
   FaultCell[] fcs = fs.findCells(new float[][][][]{fls,fps,fts});
   int nc = fcs.length;
   for (int ic=0; ic<nc; ++ic) {
     FaultCell fci = fcs[ic];
     float x1 = fci.x1+_j1;
     float x2 = fci.x2+_j2;
     float x3 = fci.x3+_j3;
     float fl = fci.fl;
     float fp = fci.fp;
     float ft = fci.ft;
     fcs[ic] = new FaultCell(x1,x2,x3,fl,fp,ft);
   }
   return fs.findSkins(fcs);
 }


 public float[][][][] faultImagesFromCellsJake(
   int n1, int n2, int n3, FaultCell[] cells) {
   setCells(cells);
   float[][][] fls = new float[_n3][_n2][_n1];
   float[][][] fps = new float[_n3][_n2][_n1];
   float[][][] fts = new float[_n3][_n2][_n1];
   float[][][] g11 = new float[_n3][_n2][_n1];
   float[][][] g12 = new float[_n3][_n2][_n1];
   float[][][] g13 = new float[_n3][_n2][_n1];
   float[][][] g22 = new float[_n3][_n2][_n1];
   float[][][] g23 = new float[_n3][_n2][_n1];
   float[][][] g33 = new float[_n3][_n2][_n1];
   for (FaultCell cell:cells) {
     float fli = cell.fl;
     int i1 = cell.i1-_j1;
     int i2 = cell.i2-_j2;
     int i3 = cell.i3-_j3;
     fls[i3][i2][i1] = fli;
     g11[i3][i2][i1] = cell.w11*fli;
     g12[i3][i2][i1] = cell.w12*fli;
     g13[i3][i2][i1] = cell.w13*fli;
     g22[i3][i2][i1] = cell.w22*fli;
     g23[i3][i2][i1] = cell.w23*fli;
     g33[i3][i2][i1] = cell.w33*fli;
   }
   System.out.println("assignments done...");
   RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(2);
   rgf1.apply000(fls,fls);
   System.out.println("fl smoothing done...");
   RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(10);
   float[][][][] gs = {g11,g22,g33,g12,g13,g23};
   for (float[][][] g:gs) {
     rgf2.apply000(g,g);
   }
   System.out.println("gaussian smoothing done...");
   float[][][] w1 = new float[_n3][_n2][_n1];
   float[][][] w2 = new float[_n3][_n2][_n1];
   float[][][] u1 = new float[_n3][_n2][_n1];
   float[][][] u2 = new float[_n3][_n2][_n1];
   float[][][] u3 = new float[_n3][_n2][_n1];
   solveEigenproblems(g11,g12,g13,g22,g23,g33,w1,w2,u1,u2,u3);
   // Compute u1 such that u3 > 0.
   for (int i3=0; i3<_n3; ++i3) {
     for (int i2=0; i2<_n2; ++i2) {
       for (int i1=0; i1<_n1; ++i1) {
         float u1i = u2[i3][i2][i1];
         float u2i = u2[i3][i2][i1];
         float u3i = u3[i3][i2][i1];
         float u1s = 1.0f-u2i*u2i-u3i*u3i;
         u1i = (u1s>0.0f)?sqrt(u1s):0.0f;
         if (u3i<0.0f) {
           u1i = -u1i;
           u2i = -u2i;
         }
         u1[i3][i2][i1] = u1i;
         u2[i3][i2][i1] = u2i;
       }
     }
   }
   System.out.println("eigentensors done...");
   float[][][] eu = fillfloat(0.01f,_n1,_n2,_n3);
   float[][][] ev = fillfloat(1.00f,_n1,_n2,_n3);
   float[][][] ew = fillfloat(1.00f,_n1,_n2,_n3);
   EigenTensors3 et = new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
   LocalSmoothingFilter lsf = new LocalSmoothingFilter();
   lsf.apply(et,400,fls,fls);
   computeStrikeDip(fls,fps,fts);
   System.out.println("structure-oriented smoothing done...");
   float[][][] fl = new float[n3][n2][n1];
   float[][][] fp = new float[n3][n2][n1];
   float[][][] ft = new float[n3][n2][n1];
   for (int i3=0; i3<_n3; ++i3) {
   for (int i2=0; i2<_n2; ++i2) {
   for (int i1=0; i1<_n1; ++i1) {
     int k1 = i1+_j1;
     int k2 = i2+_j2;
     int k3 = i3+_j3;
     fl[k3][k2][k1] = fls[i3][i2][i1];
     fp[k3][k2][k1] = fps[i3][i2][i1];
     ft[k3][k2][k1] = fts[i3][i2][i1];
   }}}
   return new float[][][][]{fl,fp,ft};
 }


 public float[][][][] faultImagesFromCells(
   int n1, int n2, int n3, FaultCell[] cells) {
   setCells(cells);
   float[][][] fls = new float[_n3][_n2][_n1];
   float[][][] fps = new float[_n3][_n2][_n1];
   float[][][] fts = new float[_n3][_n2][_n1];
   float[][][] g11 = new float[_n3][_n2][_n1];
   float[][][] g12 = new float[_n3][_n2][_n1];
   float[][][] g13 = new float[_n3][_n2][_n1];
   float[][][] g22 = new float[_n3][_n2][_n1];
   float[][][] g23 = new float[_n3][_n2][_n1];
   float[][][] g33 = new float[_n3][_n2][_n1];
   for (FaultCell cell:cells) {
     float fli = cell.fl;
     int i1 = cell.i1-_j1;
     int i2 = cell.i2-_j2;
     int i3 = cell.i3-_j3;
     fls[i3][i2][i1] = 1f;
     g11[i3][i2][i1] = cell.w11*fli;
     g12[i3][i2][i1] = cell.w12*fli;
     g13[i3][i2][i1] = cell.w13*fli;
     g22[i3][i2][i1] = cell.w22*fli;
     g23[i3][i2][i1] = cell.w23*fli;
     g33[i3][i2][i1] = cell.w33*fli;
   }
   System.out.println("assignments done...");
   RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(1);
   rgf1.apply000(fls,fls);
   System.out.println("fl smoothing done...");
   RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(10);
   float[][][][] gs = {g11,g22,g33,g12,g13,g23};
   for (float[][][] g:gs) {
     rgf2.apply000(g,g);
   }
   System.out.println("gaussian smoothing done...");
   float[][][] w1 = new float[_n3][_n2][_n1];
   float[][][] w2 = new float[_n3][_n2][_n1];
   float[][][] u1 = new float[_n3][_n2][_n1];
   float[][][] u2 = new float[_n3][_n2][_n1];
   float[][][] u3 = new float[_n3][_n2][_n1];
   solveEigenproblems(g11,g12,g13,g22,g23,g33,w1,w2,u1,u2,u3);
   // Compute u1 such that u3 > 0.
   for (int i3=0; i3<_n3; ++i3) {
     for (int i2=0; i2<_n2; ++i2) {
       for (int i1=0; i1<_n1; ++i1) {
         float u1i = u2[i3][i2][i1];
         float u2i = u2[i3][i2][i1];
         float u3i = u3[i3][i2][i1];
         float u1s = 1.0f-u2i*u2i-u3i*u3i;
         u1i = (u1s>0.0f)?sqrt(u1s):0.0f;
         if (u3i<0.0f) {
           u1i = -u1i;
           u2i = -u2i;
         }
         u1[i3][i2][i1] = u1i;
         u2[i3][i2][i1] = u2i;
       }
     }
   }
   System.out.println("eigentensors done...");
   float[][][] eu = fillfloat(0.01f,_n1,_n2,_n3);
   float[][][] ev = fillfloat(1.00f,_n1,_n2,_n3);
   float[][][] ew = fillfloat(1.00f,_n1,_n2,_n3);
   EigenTensors3 et = new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
   LocalSmoothingFilter lsf = new LocalSmoothingFilter();
   lsf.apply(et,2000,fls,fls);
   computeStrikeDip(fls,fps,fts);
   System.out.println("structure-oriented smoothing done...");
   float[][][] fl = new float[n3][n2][n1];
   float[][][] fp = new float[n3][n2][n1];
   float[][][] ft = new float[n3][n2][n1];
   for (int i3=0; i3<_n3; ++i3) {
   for (int i2=0; i2<_n2; ++i2) {
   for (int i1=0; i1<_n1; ++i1) {
     int k1 = i1+_j1;
     int k2 = i2+_j2;
     int k3 = i3+_j3;
     fl[k3][k2][k1] = fls[i3][i2][i1];
     fp[k3][k2][k1] = fps[i3][i2][i1];
     ft[k3][k2][k1] = fts[i3][i2][i1];
   }}}
   return new float[][][][]{fl,fp,ft};
 }


 public float[][][][] faultImagesFromCells(
   int n1, int n2, int n3, Sampling sp, Sampling st, FaultCell[] cells) {
   _fcs = cells;
   float[][][][] flpt = new float[3][n3][n2][n1];
   float[][][][][] gws = gaussWeights(10.f,2.f,40,40,40,sp,st);
   faultImagesFromCells(sp,st,cells,gws,flpt);
   return flpt;
 }

 private void faultImagesFromCells(
   final Sampling sp, final Sampling st,
   final FaultCell[] cells, final float[][][][][] gws, 
   final float[][][][] flpt) 
 {
   final int nc = cells.length;
   final float[][][] fl = flpt[0];
   final float[][][] fp = flpt[1];
   final float[][][] ft = flpt[2];
   final int n3 = fl.length;
   final int n2 = fl[0].length;
   final int n1 = fl[0][0].length;
   final int m3 = gws[0][0].length;
   final int m2 = gws[0][0][0].length;
   final int m1 = gws[0][0][0][0].length;
   final int d1 = (m1-1)/2;
   final int d2 = (m2-1)/2;
   final int d3 = (m3-1)/2;
   final float[][][] g11 = new float[n3][n2][n1];
   final float[][][] g12 = new float[n3][n2][n1];
   final float[][][] g13 = new float[n3][n2][n1];
   final float[][][] g22 = new float[n3][n2][n1];
   final float[][][] g23 = new float[n3][n2][n1];
   final float[][][] g33 = new float[n3][n2][n1];
   for (int ic=0; ic<nc; ++ic) {
     if(ic%1000==0) {
       System.out.println("ic="+(float)ic/(float)nc);
     }
     FaultCell fc = cells[ic];
     float fpi = fc.fp;
     float fti = fc.ft;
     final int c1 = fc.i1-d1;
     final int c2 = fc.i2-d2;
     final int c3 = fc.i3-d3;
     final float w11 = fc.w11;
     final float w12 = fc.w12;
     final float w13 = fc.w13;
     final float w22 = fc.w22;
     final float w23 = fc.w23;
     final float w33 = fc.w33;
     int it = st.indexOfNearest(fti);
     int ip = sp.indexOfNearest(fpi);
     final float[][][] gw = gws[ip][it];
     loop(m3,new Parallel.LoopInt(){
     public void compute(int i3) {
       int k3 = i3+c3; 
       if(k3>=0 && k3<n3) {
       for (int i2=0; i2<m2; ++i2) {
         int k2 = i2+c2;
         if(k2>=0 && k2<n2) {
         for (int i1=0; i1<m1; ++i1) {
           int k1 = i1+c1;
           if(k1>=0 && k1<n1) {
             float gwi = gw[i3][i2][i1];
             fl[k3][k2][k1]  += gwi;
             g11[k3][k2][k1] += w11*gwi;
             g12[k3][k2][k1] += w12*gwi;
             g13[k3][k2][k1] += w13*gwi;
             g22[k3][k2][k1] += w22*gwi;
             g23[k3][k2][k1] += w23*gwi;
             g33[k3][k2][k1] += w33*gwi;
           }
         }}
       }}
     }});
   }
   System.out.println("accumulation done...");
   solveEigenproblems(g11,g12,g13,g22,g23,g33,fp,ft);
   System.out.println("eigen-decomposition done...");
   //computeStrikeDip(fl,fp,ft);
 }

 public float[][][][][] gaussWeights(
    float sigu, float sigw, int d1, int d2, int d3, 
    Sampling sp, Sampling st) 
  {
    final int c1 = d1;
    final int c2 = d2;
    final int c3 = d3;
    final int n1 = d1*2+1;
    final int n2 = d2*2+1;
    final int n3 = d3*2+1;
    final float sw = 0.25f/(sigw*sigw);
    final float su = 0.25f/(sigu*sigu);
    final float sv = su;
    int np = sp.getCount();
    int nt = st.getCount();
    final float[][][][][] gws = new float[np][nt][n3][n2][n1];
    for (int ip=0; ip<np; ++ip) {
    for (int it=0; it<nt; ++it) {
      float fpi = (float)sp.getValue(ip);
      float fti = (float)st.getValue(it);
      final float[][][] gwi = gws[ip][it];
      float[] u = faultDipVectorFromStrikeAndDip(fpi,fti);
      float[] v = faultStrikeVectorFromStrikeAndDip(fpi,fti);
      float[] w = faultNormalVectorFromStrikeAndDip(fpi,fti);
      float u1 = u[0];
      float u2 = u[1];
      float u3 = u[2];
      float v1 = v[0];
      float v2 = v[1];
      float v3 = v[2];
      float w1 = w[0];
      float w2 = w[1];
      float w3 = w[2];
      final float w11 = w1*w1;
      final float w12 = w1*w2;
      final float w13 = w1*w3;
      final float w22 = w2*w2;
      final float w23 = w2*w3;
      final float w33 = w3*w3;
      final float v11 = v1*v1;
      final float v12 = v1*v2;
      final float v13 = v1*v3;
      final float v22 = v2*v2;
      final float v23 = v2*v3;
      final float v33 = v3*v3;
      final float u11 = u1*u1;
      final float u12 = u1*u2;
      final float u13 = u1*u3;
      final float u22 = u2*u2;
      final float u23 = u2*u3;
      final float u33 = u3*u3;
      loop(n3,new LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float dx1 = i1-c1;
          float dx2 = i2-c2;
          float dx3 = i3-c3;
          float d11 = dx1*dx1;
          float d22 = dx2*dx2;
          float d33 = dx3*dx3;
          float d12 = dx1*dx2;
          float d13 = dx1*dx3;
          float d23 = dx2*dx3;

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
          gss += (wd1+wd2+wd3+wds)*sw;
          gss += (ud1+ud2+ud3+uds)*su;
          gss += (vd1+vd2+vd3+vds)*sv;
          gwi[i3][i2][i1] = exp(-gss);
        }}
      }});
    }}
    return gws;
  }


 private void faultImagesFromCells(
    final int dt, final FaultCell[] fc, final float[][][][] flpt)
  {
    final float[][][] fl = flpt[0];
    final float[][][] fp = flpt[1];
    final float[][][] ft = flpt[2];
    int nc = fc.length;
    final int n3 = fl.length;
    final int n2 = fl[0].length;
    final int n1 = fl[0][0].length;
    float sigmaNor = 4.0f;
    final float mark = -360f;
    final float[][][] fpt = fillfloat(mark,n1,n2,n3);
    final float[][] xc = new float[3][nc];
    setKdTreeNodes(fc,xc,fpt);
    final KdTree kt = new KdTree(xc);
    final float sv = 0.25f/(dt*dt); 
    final float su = 0.25f/(dt*dt); 
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          xmin[1] = i2-dt; xmax[1] = i2+dt;
          xmin[2] = i3-dt; xmax[2] = i3+dt;
          xmin[0] = i1-dt/4; xmax[0] = i1+dt/4;
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<10){continue;}
          ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            FaultCell fci = fc[ip];
            fcl.add(fci);
          }
          FaultCell[] cells = fcl.toArray(new FaultCell[0]);
          for (FaultCell fci:cells) {
            float x1i = fci.i1;
            float x2i = fci.i2;
            float x3i = fci.i3;

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
            float w22 = fci.w22;
            float w33 = fci.w33;
            float w12 = fci.w12;
            float w13 = fci.w13;
            float w23 = fci.w23;

            float u11 = fci.u11;
            float u22 = fci.u22;
            float u33 = fci.u33;
            float u12 = fci.u12;
            float u13 = fci.u13;
            float u23 = fci.u23;

            float v11 = fci.v11;
            float v22 = fci.v22;
            float v33 = fci.v33;
            float v12 = fci.v12;
            float v13 = fci.v13;
            float v23 = fci.v23;

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
            float flc = fci.fl;
            float wpi = pow(flc,2.f);
            gss += (wd1+wd2+wd3+wds)*sw;
            gss += (ud1+ud2+ud3+uds)*su;
            gss += (vd1+vd2+vd3+vds)*sv;
            float fli = exp(-gss)*wpi;
            fl[i3][i2][i1] += fli;
          }
        }
      }
    }});
    computeStrikeDip(fl,fp,ft);
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


  private void setKdTreeNodes(
    FaultCell[] fc, float[][] xc, float[][][] fp) {
    float mark = -360f;
    int nc = fc.length;
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    float[][][] pt = fillfloat(mark,n1,n2,n3);
    for (int ic=0; ic<nc; ic++) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      xc[0][ic] = i1;
      xc[1][ic] = i2;
      xc[2][ic] = i3;
      pt[i3][i2][i1] = fc[ic].fp;
      w1[i3][i2][i1] = fc[ic].w1;
      w2[i3][i2][i1] = fc[ic].w2;
      w3[i3][i2][i1] = fc[ic].w3;
    }
  }

  private void computeStrikeDip(
    float[][][] fl, float[][][] fp, float[][][] ft) 
  {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    LocalOrientFilterP lof = new LocalOrientFilterP(8,4);
    lof.applyForNormal(fl,u1,u2,u3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fli = fl[i3][i2][i1];
      if(fli>0f) {
        int k1 = i1;
        int k2 = i2;
        int k3 = i3;
        if(k1==0) {k1=1;}
        if(k2==0) {k2=1;}
        if(k3==0) {k3=1;}
        if(k1==n1-1) {k1=n1-2;}
        if(k2==n2-1) {k2=n2-2;}
        if(k3==n3-1) {k3=n3-2;}
        float u1i = -u1[k3][k2][k1];
        float u2i = -u2[k3][k2][k1];
        float u3i = -u3[k3][k2][k1];
        if(u2i!=0.0f && u3i!=0.0f) {
          ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
          fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
        }
      }
    }}}
  }
 


 public float[][][][] faultSlopes(int n1, int n2, int n3, FaultSkin skin) {
    float[][][] fls = new float[n3][n2][n1];
    float[][][] g11 = new float[n3][n2][n1];
    float[][][] g12 = new float[n3][n2][n1];
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g22 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    float[][][] g33 = new float[n3][n2][n1];
    float pmin = FLT_MAX;
    float pmax = FLT_MIN;
    for (FaultCell cell:skin) {
      int i1 = cell.i1;
      int i2 = cell.i2;
      int i3 = cell.i3;
      float w1 = cell.w1;
      float w2 = cell.w2;
      float w3 = cell.w3;
      float fl = cell.fl;
      fls[i3][i2][i1] = fl;
      g11[i3][i2][i1] = w1*w1*fl;
      g12[i3][i2][i1] = w1*w2*fl;
      g13[i3][i2][i1] = w1*w3*fl;
      g22[i3][i2][i1] = w2*w2*fl;
      g23[i3][i2][i1] = w2*w3*fl;
      g33[i3][i2][i1] = w3*w3*fl;
      float p2 = -w2/w1;
      float p3 = -w3/w1;
      if(p2<pmin) {pmin=p2;}
      if(p3<pmin) {pmin=p3;}
      if(p2>pmax) {pmax=p2;}
      if(p3>pmax) {pmax=p3;}
    }
    pmin -=5; pmax +=5;
    RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(8.0);
    RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(64.0);
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {fls,g11,g22,g33,g12,g13,g23};
    for (float[][][] g:gs) {
      rgf1.apply0XX(g,h); copy(g,h);
      rgf2.applyX0X(h,g); copy(h,g);
      rgf2.applyXX0(g,h); copy(h,g);
    }
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float p2i = pmin;
      float p3i = pmin;
      if (u1i!=0f) {
        p2i = -u2i/u1i;
        p3i = -u3i/u1i;
      } 
      if(p2i<=pmin||p3i<=pmin) {continue;}
      if(p2i>=pmax||p3i>=pmax) {continue;}
      p2[i3][i2][i1] = p2i;
      p3[i3][i2][i1] = p3i;
    }}}
    div(fls,max(fls),fls);
    return new float[][][][]{p2,p3,fls};
  }

  
  public float[][][] faultIndicator(int n1, int n2, int n3, FaultSkin skin) {
    FaultCell[] fcs = skin.getCells();
    setCells(n1,n2,n3,fcs);
    System.out.println("fault setting done...");
    float[][][] sfs = fillfloat(-30f,n1,n2,n3);
    float[][][] fl  = new float[_n3][_n2][_n1];
    float[][][] ws  = new float[_n3][_n2][_n1];
    float[][][] u1  = new float[_n3][_n2][_n1];
    float[][][] u2  = new float[_n3][_n2][_n1];
    float[][][] u3  = new float[_n3][_n2][_n1];
    float[][][] g11 = new float[_n3][_n2][_n1];
    float[][][] g12 = new float[_n3][_n2][_n1];
    float[][][] g13 = new float[_n3][_n2][_n1];
    float[][][] g22 = new float[_n3][_n2][_n1];
    float[][][] g23 = new float[_n3][_n2][_n1];
    float[][][] g33 = new float[_n3][_n2][_n1];
    initialTensors(fl,ws,g11,g12,g13,g22,g23,g33);
    System.out.println("tensors done...");
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    System.out.println("normals done...");
    ScreenPoissonSurfer sps = new ScreenPoissonSurfer();
    sps.setSmoothings(20,20,20);
    mul(ws,u1,u1);
    mul(ws,u2,u2);
    mul(ws,u3,u3);
    float[][][] sft = sps.saltIndicator(fl,u1,u2,u3);
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      sfs[i3+_j3][i2+_j2][i1+_j1] = sft[i3][i2][i1];
    }}}
    System.out.println("fault indicator done...");
    return sfs;
  }

  private void setCells(int n1, int n2, int n3, FaultCell[] cells) {
    _j1 = 0;
    _j2 = 0;
    _j3 = 0;
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _cells = new FaultCell[_n3][_n2][_n1];
    for (FaultCell cell:cells)
      set(cell);
  }


  private void setCells(FaultCell[] cells) {
    int i1min = Integer.MAX_VALUE;
    int i2min = Integer.MAX_VALUE;
    int i3min = Integer.MAX_VALUE;
    int i1max = -i1min;
    int i2max = -i2min;
    int i3max = -i3min;
    for (FaultCell cell:cells) {
      if (cell.i1<i1min) i1min = cell.i1;
      if (cell.i2<i2min) i2min = cell.i2;
      if (cell.i3<i3min) i3min = cell.i3;
      if (cell.i1>i1max) i1max = cell.i1;
      if (cell.i2>i2max) i2max = cell.i2;
      if (cell.i3>i3max) i3max = cell.i3;
    }
    _j1 = i1min;
    _j2 = i2min;
    _j3 = i3min;
    _n1 = 1+i1max-i1min;
    _n2 = 1+i2max-i2min;
    _n3 = 1+i3max-i3min;
    _cells = new FaultCell[_n3][_n2][_n1];
    for (FaultCell cell:cells)
      set(cell);
  }


  private void set(FaultCell cell) {
    int i1 = cell.i1-_j1;
    int i2 = cell.i2-_j2;
    int i3 = cell.i3-_j3;
    i1 = max(i1,0); i1 = min(i1,_n1-1);
    i2 = max(i2,0); i2 = min(i2,_n2-1);
    i3 = max(i3,0); i3 = min(i3,_n3-1);
    if(_cells[i3][i2][i1]==null){
      _cells[i3][i2][i1] = cell;
    }

  }

  private void setNull(FaultCell cell) {
    int i1 = cell.i1-_j1;
    int i2 = cell.i2-_j2;
    int i3 = cell.i3-_j3;
    i1 = max(i1,0); i1 = min(i1,_n1-1);
    i2 = max(i2,0); i2 = min(i2,_n2-1);
    i3 = max(i3,0); i3 = min(i3,_n3-1);
    _cells[i3][i2][i1]=null;
  }

  public float[][][] initialTensorsTest(
    int n1, int n2, int n3, FaultSkin skin) {
    setCells(n1,n2,n3,skin.getCells());
    int[][][] mk = new int[_n3][_n2][_n1];
    float[][][] fls = new float[n3][n2][n1];
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      FaultCell cell = _cells[i3][i2][i1];
      if (cell!=null&&mk[i3][i2][i1]!=1) {
        FaultCell cm = cell;
        FaultCell[] cells = findOverlapCells(i1,i2,i3,cell);
        int nc = cells.length;
        if (nc>1) {
        int nbm = nabors(cell);
        for (int ic=0; ic<nc; ++ic) {
          FaultCell fci = cells[ic];
          if(notNearbyCells(cell,fci)) {
            int nbi = nabors(fci);
            if(nbi>nbm) {cm = fci;nbm = nbi;} 
            else {setNull(fci);}
          }
        }}
        float fl = cm.fl;
        int k1 = cm.i1-_j1;
        int k2 = cm.i2-_j2;
        int k3 = cm.i3-_j3;
        mk[k3][k2][k1] = 1;
        fls[k3][k2][k1] = fl; 
      }
    }}}
    return fls;
  }


  private void initialTensors(
    float[][][] fls, float[][][] wss,
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33) 
  {
    int[][][] mk = new int[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      FaultCell cell = _cells[i3][i2][i1];
      if (cell!=null&&mk[i3][i2][i1]!=1) {
        FaultCell cm = cell;
        FaultCell[] cells = findOverlapCells(i1,i2,i3,cell);
        int nc = cells.length;
        if (nc>1) {
        int nbm = nabors(cell);
        for (int ic=0; ic<nc; ++ic) {
          FaultCell fci = cells[ic];
          if(notNearbyCells(cell,fci)) {
            int nbi = nabors(fci);
            if(nbi>nbm) {cm = fci;nbm = nbi;} 
            else {setNull(fci);}
          }
        }}
        float w1 = cm.w1;
        float w2 = cm.w2;
        float w3 = cm.w3;
        float fl = cm.fl;
        int k1 = cm.i1-_j1;
        int k2 = cm.i2-_j2;
        int k3 = cm.i3-_j3;
        mk[k3][k2][k1] = 1;
        fls[k3][k2][k1] = fl; 
        wss[k3][k2][k1] = fl; 
        g11[k3][k2][k1] = w1*w1*fl; 
        g12[k3][k2][k1] = w1*w2*fl; 
        g13[k3][k2][k1] = w1*w3*fl; 
        g22[k3][k2][k1] = w2*w2*fl; 
        g23[k3][k2][k1] = w2*w3*fl; 
        g33[k3][k2][k1] = w3*w3*fl; 
      }
    }}}
    RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(8.0);
    RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(64.0);
    float[][][] h = new float[_n3][_n2][_n1];
    float[][][][] gs = {wss,g11,g22,g33,g12,g13,g23};
    for (float[][][] g:gs) {
      rgf1.apply0XX(g,h); copy(g,h);
      rgf2.applyX0X(h,g); copy(h,g);
      rgf2.applyXX0(g,h); copy(h,g);
    }
  }

  private boolean notNearbyCells(FaultCell c1, FaultCell c2) {
    int d2 = c1.i2-c2.i2;
    int d3 = c1.i3-c2.i3;
    float ds = d2*d2+d3*d3;
    if(ds<=8f) {return false;}
    else {return true;}
  }

  public FaultCell[] findOverlapCells(int c1, int c2, int c3, FaultCell cell) {
    ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
    for (int i2=0; i2<_n2;  ++i2) {
      FaultCell fc = _cells[c3][i2][c1];
      if(fc!=null) fcl.add(fc);
    }
    for (int i3=0; i3<_n3;  ++i3) {
      FaultCell fc = _cells[i3][c2][c1];
      if(fc!=null) fcl.add(fc);
    }
    return fcl.toArray(new FaultCell[0]);
  }

  private int nabors(FaultCell cell) {
    final int c1 = cell.i1-_j1;
    final int c2 = cell.i2-_j2;
    final int c3 = cell.i3-_j3;
    final int[] nc = new int[1];
    final int b2 = max(c2-100,0);
    final int b3 = max(c3-100,0);
    final int e2 = min(c2+100,_n2-1);
    final int e3 = min(c3+100,_n3-1);
    final float fp = cell.fp;
    Parallel.loop(b3,e3+1,1,new Parallel.LoopInt() {
    public void compute(int k3) {
      for (int k2=b2; k2<=e2; ++k2) {
      FaultCell fc = _cells[k3][k2][c1];
      if (fc!=null) {
        float del = fc.fp-fp;
        float fp = min(abs(del),abs(del+360.0f),abs(del-360.0f));
        if(fp<=10f) {nc[0] +=1;}
      }
    }
    }});
    return nc[0];
  }

  private FaultCell[] findOverlapCells(FaultCell cell) {
    //search in direction of fault normal
    float d  = 3f;
    float dm = sqrt(_n1*_n1+_n2*_n2+_n3*_n3);
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    while (d<dm) {
    }
    //search in opposite direction of fault normal
    return null;
  }


  private int cellsLR(FaultCell cell) {
    int nc = 0;
    FaultCell c = cell;
    for (FaultCell cl=c.cl; cl!=null && cl!=cell; cl=c.cl)
      c = cl;
    FaultCell cLeft = c;
    for (c=c.cr; c!=null && c!=cLeft; c=c.cr) 
      nc ++;
    return nc;
  }


  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] fp, final float[][][] ft)
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      double[][] a = new double[3][3];
      double[][] z = new double[3][3];
      double[] e = new double[3];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          a[0][0] = g11[i3][i2][i1];
          if(a[0][0]==0f) {continue;}
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
    final float[][][] w1, final float[][][] w2, final float[][][] u1, 
    final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
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
            float w1i = (float)z[2][0];
            float w2i = (float)z[2][1];
            float w3i = (float)z[2][2];
            if (u1i<0.0f) {
              u1i = -u1i;
              u2i = -u2i;
              u3i = -u3i;
            }

            if (w3i<0.0f) {
              w1i = -w1i;
              w2i = -w2i;
            }
            u1[i3][i2][i1] = u1i;
            u2[i3][i2][i1] = u2i;
            u3[i3][i2][i1] = u3i;
            w1[i3][i2][i1] = w1i;
            w2[i3][i2][i1] = w2i;
          }
        }
      }
    });
  }



  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
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
            if (u1i<0.0f) {
              u1i = -u1i;
              u2i = -u2i;
              u3i = -u3i;
            }
            if (u1!=null) u1[i3][i2][i1] = u1i;
            if (u2!=null) u2[i3][i2][i1] = u2i;
            if (u3!=null) u3[i3][i2][i1] = u3i;
          }
        }
      }
    });
  }

  private KdTree setStrikeKdTree() {
    int nc = _fcs.length;
    float[][] pc = new float[1][nc];
    for (int ic=0; ic<nc; ++ic) {
      pc[0][ic] = _fcs[ic].fp; 
    }
    return new KdTree(pc);
  }

  // find cells in a strike window
  private FaultCell[] cellsInStrikeWd(float pmi, float ppi, KdTree kt) {
    float[] pm1 = new float[1];
    float[] pp1 = new float[1];
    float[] pm2 = new float[1];
    float[] pp2 = new float[1];
    float[] pm3 = new float[1];
    float[] pp3 = new float[1];
    pm3[0] = pmi; pp3[0] = ppi;
    int[] id1=null, id2=null, id3=null;
    if (pmi<0.0f) {
      pm3[0] = 0.0f;
      pp1[0] = 360f;
      pm1[0] = 360f+pmi;
      id1 = kt.findInRange(pm1,pp1);
    } 
    if(ppi>360f) {
      pp3[0] = 360f;
      pm2[0] = 0.0f;
      pp2[0] = ppi-360f;
      id2 = kt.findInRange(pm2,pp2);
    }
    id3 = kt.findInRange(pm3,pp3);
    int nd1=0, nd2=0, nd3=0;
    if(id1!=null){nd1=id1.length;}
    if(id2!=null){nd2=id2.length;}
    if(id3!=null){nd3=id3.length;}
    int nd = nd1+nd2+nd3;
    int ic = 0;
    FaultCell[] fcs = new FaultCell[nd];
    if(nd>0) {
      if(id1!=null) {
        for (int ik=0; ik<nd1; ++ik){
          fcs[ic] = _fcs[id1[ik]];
          ic++;
        }
      }
      if(id2!=null) {
        for (int ik=0; ik<nd2; ++ik){
          fcs[ic] = _fcs[id2[ik]];
          ic++;
        }
      }
      if(id3!=null) {
        for (int ik=0; ik<nd3; ++ik){
          fcs[ic] = _fcs[id3[ik]];
          ic++;
        }
      }
    }
    return fcs;
  }



  private int _j1,_j2,_j3; // min cell indices
  private int _n1,_n2,_n3; // numbers of cells
  private FaultCell[][][] _cells; // array of cells
  private FaultCell[] _fcs; // list of cells

}


