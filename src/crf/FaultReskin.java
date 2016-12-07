package crf;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

import mef.*;
import util.*;
import static mef.FaultGeometry.*;


/**
 * Construct smooth and single-valued fault surface using 
 * the screen poisson method.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.02.29
 */

public class FaultReskin {


 public FaultSkin[] applyForSkins(
   int n1, int n2, int n3, int size, FaultSkin[] skins) 
 {
   int nk = skins.length; 
   ArrayList<FaultSkin> skinList = new ArrayList<FaultSkin>();
   for (int ik=0; ik<nk; ++ik) {
     System.out.println("ik="+ik);
     FaultSkin skin = skins[ik];
     int nc = skin.size();
     int[] k1 = new int[nc];
     int[] k2 = new int[nc];
     int[] k3 = new int[nc];
     float[] fl = new float[nc];
     float[] w1 = new float[nc];
     float[] w2 = new float[nc];
     float[] w3 = new float[nc];
     FaultCell[] cells = skin.getCells();
     for (int ic=0; ic<nc; ++ic) {
       FaultCell cell = cells[ic];
       int i1 = cell.getI1();
       int i2 = cell.getI2();
       int i3 = cell.getI3();
       i1 = min(i1,n1-1);
       i2 = min(i2,n2-1);
       i3 = min(i3,n3-1);
       i1 = max(i1,0);
       i2 = max(i2,0);
       i3 = max(i3,0);
       k1[ic] = i1;
       k2[ic] = i2;
       k3[ic] = i3;
       fl[ic] = cell.getFl();
       w1[ic] = cell.getW1();
       w2[ic] = cell.getW2();
       w3[ic] = cell.getW3();
     }
     int b1 = min(k1);
     int b2 = min(k2);
     int b3 = min(k3);
     int e1 = max(k1);
     int e2 = max(k2);
     int e3 = max(k3);
     int m1 = e1-b1+1;
     int m2 = e2-b2+1;
     int m3 = e3-b3+1;
     float[][][] fls = getFaultImages(b1,b2,b3,m1,m2,m3,k1,k2,k3,fl,w1,w2,w3);
     FaultSkin[] sks = reskin(b1,b2,b3,size,fls);
     for (FaultSkin ski:sks) {
       skinList.add(ski);
     }
   }
   return skinList.toArray(new FaultSkin[0]);
 }

 public FaultSkin[] reskin(int b1, int b2, int b3, int size, float[][][] fls) {
   int n3 = fls.length;
   int n2 = fls[0].length;
   int n1 = fls[0][0].length;
   float[][][] fps = new float[n3][n2][n1];
   float[][][] fts = new float[n3][n2][n1];
   computeStrikeDip(fls,fps,fts);
   FaultSkinner fs = new FaultSkinner();
   fs.setGrowLikelihoods(0.1f,0.6f);
   fs.setMaxDeltaStrike(10);
   fs.setMaxPlanarDistance(0.2f);
   fs.setMinSkinSize(size);
   FaultCell[] fcs = fs.findCells(new float[][][][]{fps,fls,fts});
   int nc = fcs.length;
   for (int ic=0; ic<nc; ++ic) {
     FaultCell cell = fcs[ic];
     float x1i = cell.getX1();
     float x2i = cell.getX2();
     float x3i = cell.getX3();
     float fli = cell.getFl();
     float fpi = cell.getFp();
     float fti = cell.getFt();
     fcs[ic] = new FaultCell(x1i+b1,x2i+b2,x3i+b3,fli,fpi,fti);
   }
   return fs.findSkins(fcs);
 }

 public float[][][][] getFaultImagesX(
   int b1, int b2, int b3, int n1, int n2, int n3, 
   int[] k1, int[] k2, int[] k3, float[] fl, float[] v1, float[] v2, float[] v3) 
 {
   float[][][] fls = new float[n3][n2][n1];
   float[][][] fps = new float[n3][n2][n1];
   float[][][] fts = new float[n3][n2][n1];
   float[][][] g11 = new float[n3][n2][n1];
   float[][][] g12 = new float[n3][n2][n1];
   float[][][] g13 = new float[n3][n2][n1];
   float[][][] g22 = new float[n3][n2][n1];
   float[][][] g23 = new float[n3][n2][n1];
   float[][][] g33 = new float[n3][n2][n1];
   int nc = k1.length;
   for (int ic=0; ic<nc; ++ic) {
     int i1 = k1[ic]-b1;
     int i2 = k2[ic]-b2;
     int i3 = k3[ic]-b3;
     float fli = fl[ic];
     float w1i = v1[ic];
     float w2i = v2[ic];
     float w3i = v3[ic];
     float w11 = w1i*w1i;
     float w12 = w1i*w2i;
     float w13 = w1i*w3i;
     float w22 = w2i*w2i;
     float w23 = w2i*w3i;
     float w33 = w3i*w3i;
     fls[i3][i2][i1] = fli;
     g11[i3][i2][i1] = w11*fli;
     g12[i3][i2][i1] = w12*fli;
     g13[i3][i2][i1] = w13*fli;
     g22[i3][i2][i1] = w22*fli;
     g23[i3][i2][i1] = w23*fli;
     g33[i3][i2][i1] = w33*fli;
   }
   RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(4);
   rgf1.apply000(fls,fls);
   fls = pow(abs(fls),0.5f);
   fls = sub(fls,min(fls));
   fls = div(fls,max(fls));
   RecursiveGaussianFilter rgf2 = new RecursiveGaussianFilter(10);
   float[][][][] gs = {g11,g22,g33,g12,g13,g23};
   for (float[][][] g:gs) {
     rgf2.apply0XX(g,g);
     rgf1.applyX0X(g,g);
     rgf1.applyXX0(g,g);
   }
   System.out.println("gaussian smoothing done...");
   float[][][] w1 = new float[n3][n2][n1];
   float[][][] w2 = new float[n3][n2][n1];
   float[][][] u1 = new float[n3][n2][n1];
   float[][][] u2 = new float[n3][n2][n1];
   float[][][] u3 = new float[n3][n2][n1];
   solveEigenproblems(g11,g12,g13,g22,g23,g33,w1,w2,u1,u2,u3);
   // Compute u1 such that u3 > 0.
   for (int i3=0; i3<n3; ++i3) {
     for (int i2=0; i2<n2; ++i2) {
       for (int i1=0; i1<n1; ++i1) {
         float u1i = u1[i3][i2][i1];
         float u2i = u2[i3][i2][i1];
         float u3i = u3[i3][i2][i1];
         /*
         float u1s = 1.0f-u2i*u2i-u3i*u3i;
         u1i = (u1s>0.0f)?sqrt(u1s):0.0f;
         if (u3i<0.0f) {
           u1i = -u1i;
           u2i = -u2i;
         }
         u1[i3][i2][i1] = u1i;
         u2[i3][i2][i1] = u2i;
         */
         if(u2i!=0.0f && u3i!=0.0f) {
          fts[i3][i2][i1] = faultDipFromNormalVector(-u1i,-u2i,-u3i);
          fps[i3][i2][i1] = faultStrikeFromNormalVector(-u1i,-u2i,-u3i);
         }
       }
     }
   }
   System.out.println("fpmax="+max(fps));
   System.out.println("ftmax="+max(fts));
   /*
   System.out.println("eigentensors done...");
   float[][][] eu = fillfloat(0.1f,n1,n2,n3);
   float[][][] ev = fillfloat(1.0f,n1,n2,n3);
   float[][][] ew = fillfloat(1.0f,n1,n2,n3);
   EigenTensors3 et = new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
   LocalSmoothingFilter lsf = new LocalSmoothingFilter();
   lsf.apply(et,10,fls,fls);
   fls = sub(fls,min(fls));
   fls = div(fls,max(fls));
   */
   return new float[][][][]{fls,fps,fts};
 }


 public float[][][] getFaultImages(
   int b1, int b2, int b3, int n1, int n2, int n3, 
   int[] k1, int[] k2, int[] k3, float[] fl, float[] v1, float[] v2, float[] v3) 
 {
   float[][][] fls = new float[n3][n2][n1];
   float[][][] g11 = new float[n3][n2][n1];
   float[][][] g12 = new float[n3][n2][n1];
   float[][][] g13 = new float[n3][n2][n1];
   float[][][] g22 = new float[n3][n2][n1];
   float[][][] g23 = new float[n3][n2][n1];
   float[][][] g33 = new float[n3][n2][n1];
   int nc = k1.length;
   for (int ic=0; ic<nc; ++ic) {
     int i1 = k1[ic]-b1;
     int i2 = k2[ic]-b2;
     int i3 = k3[ic]-b3;
     float fli = fl[ic];
     float w1i = v1[ic];
     float w2i = v2[ic];
     float w3i = v3[ic];
     float w11 = w1i*w1i;
     float w12 = w1i*w2i;
     float w13 = w1i*w3i;
     float w22 = w2i*w2i;
     float w23 = w2i*w3i;
     float w33 = w3i*w3i;
     fls[i3][i2][i1] = fli;
     g11[i3][i2][i1] = w11*fli;
     g12[i3][i2][i1] = w12*fli;
     g13[i3][i2][i1] = w13*fli;
     g22[i3][i2][i1] = w22*fli;
     g23[i3][i2][i1] = w23*fli;
     g33[i3][i2][i1] = w33*fli;
   }
   RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(1);
   rgf1.apply000(fls,fls);
   fls = sub(fls,min(fls));
   fls = div(fls,max(fls));

   RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(10);
   float[][][][] gs = {g11,g22,g33,g12,g13,g23};
   for (float[][][] g:gs) {
     rgf2.apply000(g,g);
   }
   System.out.println("gaussian smoothing done...");
   float[][][] w1 = new float[n3][n2][n1];
   float[][][] w2 = new float[n3][n2][n1];
   float[][][] u1 = new float[n3][n2][n1];
   float[][][] u2 = new float[n3][n2][n1];
   float[][][] u3 = new float[n3][n2][n1];
   solveEigenproblems(g11,g12,g13,g22,g23,g33,w1,w2,u1,u2,u3);
   // Compute u1 such that u3 > 0.
   for (int i3=0; i3<n3; ++i3) {
     for (int i2=0; i2<n2; ++i2) {
       for (int i1=0; i1<n1; ++i1) {
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
   float[][][] eu = fillfloat(0.1f,n1,n2,n3);
   float[][][] ev = fillfloat(1.0f,n1,n2,n3);
   float[][][] ew = fillfloat(1.0f,n1,n2,n3);
   EigenTensors3 et = new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
   LocalSmoothingFilter lsf = new LocalSmoothingFilter();
   lsf.applySmoothS(fls,fls);
   lsf.apply(et,10,fls,fls);
   fls = sub(fls,min(fls));
   fls = div(fls,max(fls));
   return fls;
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
    LocalOrientFilter lof = new LocalOrientFilter(8,4);
    lof.applyForNormal(fl,u1,u2,u3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fli = fl[i3][i2][i1];
      if(fli>0.001f) {
        int k1 = i1;
        int k2 = i2;
        int k3 = i3;
        if(k1==0) {k1=1;}
        if(k2==0) {k2=1;}
        if(k3==0) {k3=1;}
        if(k1==n1-1) {k1=n1-2;}
        if(k2==n2-1) {k2=n2-2;}
        if(k3==n3-1) {k3=n3-2;}
        float u1i = u1[k3][k2][k1];
        float u2i = u2[k3][k2][k1];
        float u3i = u3[k3][k2][k1];
        if(u2i!=0.0f && u3i!=0.0f) {
          ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
          fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
        }
      }
    }}}
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

  private int _j1,_j2,_j3; // min cell indices
  private int _n1,_n2,_n3; // numbers of cells
  private FaultCell[][][] _cells; // array of cells
  private FaultCell[] _fcs; // list of cells

}


