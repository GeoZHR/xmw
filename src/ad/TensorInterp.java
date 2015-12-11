/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ad;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import ipfx.FaultCell;

/**
 * Complete surface reconstruction from noisy and incomplete point cloud. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.06
 */


public class TensorInterp {

  public void setParameters(float t, int m, float tm) {
    _t = t;
    _m = m;
    _tm = tm;
  }

  public float[][][] apply(
    float[][] fx, float[][] u1, float[][] u2) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g11 = new float[n2][n1];
    float[][] g12 = new float[n2][n1];
    float[][] g22 = new float[n2][n1];
    EigenTensors2 et = initialTensors(fx,u1,u2,g11,g12,g22);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.applySmoothS(g11,g11);
    lsf.applySmoothS(g12,g12);
    lsf.applySmoothS(g22,g22);
    FedStep fs = new FedStep(_t,_m,_tm);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    for (int im=0; im<_m; im++) {
    for (int ic=0; ic<nc; ++ic) {
      float tc = ts[ic];
      applyLaplacianX(et,-tc,copy(g11),g11);
      applyLaplacianX(et,-tc,copy(g12),g12);
      applyLaplacianX(et,-tc,copy(g22),g22);
      et = updateTensors(g11,g12,g22);
    }}
    return salientMap(g11,g12,g22);
  }

  public float[][] applyX(
    float[][] fx, float[][] u1, float[][] u2) {
    mul(fx,u1,u1);
    mul(fx,u2,u2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.applySmoothS(u1,u1);
    lsf.applySmoothS(u2,u2);
    correctVectors(u1,u2);
    EigenTensors2 et = updateTensorsX(u1,u2);
    FedStep fs = new FedStep(_t,_m,_tm);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    for (int im=0; im<_m; im++) {
    for (int ic=0; ic<nc; ++ic) {
      float tc = ts[ic];
      applyLaplacianX(et,-tc,copy(u1),u1);
      applyLaplacianX(et,-tc,copy(u2),u2);
      correctVectors(u1,u2);
      et = updateTensorsX(u1,u2);
    }}
    mul(u1,u1,u1);
    mul(u2,u2,u2);
    add(u1,u2,u2);
    return sqrt(u2);//)salientMap(g11,g12,g22);
  }

  private void correctVectors(float[][] u1, float[][] u2) {
    int n2 = u1.length;
    int n1 = u1[0].length;
    float th = 0.05f;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      if (u1i>0f&&abs(u2i/u1i)>th) {
        u1[i2][i1] = -u1i;
        u2[i2][i1] = -u2i;
      }
    }}
  }

  private EigenTensors2 updateTensorsX(float[][] u1, float[][] u2) {
    int n2 = u1.length;
    int n1 = u1[0].length;
    float[] es = new float[]{0.0001f,1.0f};
    EigenTensors2 et = new EigenTensors2(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      float usi = sqrt(u1i*u1i+u2i*u2i);
      if (usi!=0f) {
        u1i /= usi;
        u2i /= usi;
      }
      et.setEigenvalues(i1,i2,es);
      et.setEigenvectorU(i1,i2,u1i,u2i);
    }}
    return et;
  }


  // input are oriented points
  public float[][][][] apply(int n1, int n2, int n3, FaultCell[] cells) {
    float[][][] fx = new float[n3][n2][n1];
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    for (FaultCell cell:cells) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      fx[i3][i2][i1] = cell.getFl();
      u1[i3][i2][i1] = cell.getW1();
      u2[i3][i2][i1] = cell.getW2();
      u3[i3][i2][i1] = cell.getW3();
    }
    return apply(fx,u1,u2,u3);
  }

  // input are positions and orientations of points
  public float[][][][] apply(
    float[][][] fx, float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] g11 = new float[n3][n2][n1];
    float[][][] g12 = new float[n3][n2][n1];
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g22 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    float[][][] g33 = new float[n3][n2][n1];
    EigenTensors3 et = initialTensors(fx,u1,u2,u3,g11,g12,g13,g22,g23,g33);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.applySmoothS(g11,g11);
    lsf.applySmoothS(g12,g12);
    lsf.applySmoothS(g13,g13);
    lsf.applySmoothS(g22,g22);
    lsf.applySmoothS(g23,g23);
    lsf.applySmoothS(g33,g33);
    FedStep fs = new FedStep(_t,_m,_tm);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    trace("nc="+nc);
    Stopwatch sw = new Stopwatch();
    sw.start();
    int k=0;
    int nk = nc*_m;
    for (int im=0; im<_m; im++) {
    for (int ic=0; ic<nc; ++ic) {
      if (k>0) {
        double timeUsed = sw.time();
        double timeLeft = ((double)nk/(double)k-1.0)*timeUsed;
        int timeLeftSec = 1+(int)timeLeft;
        trace("Interpolation: done in "+timeLeftSec+" seconds");
      }
      k++;
      float tc = ts[ic];
      applyLaplacian(et,-tc,copy(g11),g11);
      applyLaplacian(et,-tc,copy(g12),g12);
      applyLaplacian(et,-tc,copy(g13),g13);
      applyLaplacian(et,-tc,copy(g22),g22);
      applyLaplacian(et,-tc,copy(g23),g23);
      applyLaplacian(et,-tc,copy(g33),g33);
      et = updateTensors(g11,g12,g13,g22,g23,g33);
    }}
    sw.stop();
    trace("Interpolation: done");
    return salientMap(g11,g12,g13,g22,g23,g33);
  }

  public float[][] applyForEdge(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] fy = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    LocalOrientFilter lof = new LocalOrientFilter(4,4);
    lof.applyForNormal(fx,u1,u2);
    rgf.apply10(fx,g1);
    rgf.apply01(fx,g2);
    mul(g1,u1,g1);
    mul(g2,u2,g2);
    add(g1,g2,fy);
    return abs(fy);
  }

  public EigenTensors2 initialTensorsX(
    float[][] fx, float[][] u1, float[][] u2) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g11 = new float[n2][n1];
    float[][] g12 = new float[n2][n1];
    float[][] g22 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fxi = fx[i2][i1];
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      g11[i2][i1] = fxi*u1i*u1i;
      g12[i2][i1] = fxi*u1i*u2i;
      g22[i2][i1] = fxi*u2i*u2i;
    }}
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(4);
    rgf.apply00(g11,g11);
    rgf.apply00(g12,g12);
    rgf.apply00(g22,g22);
    float[] es = new float[]{0.0001f,1f};
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    EigenTensors2 ets = new EigenTensors2(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g11i = g11[i2][i1];
      float g12i = g12[i2][i1];
      float g22i = g22[i2][i1];
      a[0][0] = g11i;
      a[0][1] = g12i;
      a[1][0] = g12i;
      a[1][1] = g22i;
      Eigen.solveSymmetric22(a,z,e);
      float u1i = z[0][0];
      float u2i = z[0][1];
      ets.setEigenvalues(i1,i2,es);
      ets.setEigenvectorU(i1,i2,u1i,u2i);
    }}
    return ets;
  }

  public EigenTensors2 initialTensors(
    float[][] fx, float[][] u1, float[][] u2,
    float[][] g11, float[][] g12, float[][] g22) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fxi = fx[i2][i1];
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      g11[i2][i1] = fxi*u1i*u1i;
      g12[i2][i1] = fxi*u1i*u2i;
      g22[i2][i1] = fxi*u2i*u2i;
    }}
    return updateTensors(g11,g12,g22);
  }


  public EigenTensors3 initialTensors(
    final float[][][] fx, 
    final float[][][] u1,  final float[][][] u2,  final float[][][] u3,
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fxi = fx[i3][i2][i1];
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        g11[i3][i2][i1] = fxi*u1i*u1i;
        g12[i3][i2][i1] = fxi*u1i*u2i;
        g13[i3][i2][i1] = fxi*u1i*u3i;
        g22[i3][i2][i1] = fxi*u2i*u2i;
        g23[i3][i2][i1] = fxi*u2i*u3i;
        g33[i3][i2][i1] = fxi*u3i*u3i;
      }}
    }});
    return updateTensors(g11,g12,g13,g22,g23,g33);
  }


  private EigenTensors2 updateTensors(
    float[][] g11, float[][] g12, float[][] g22)
  {
    int n2 = g11.length;
    int n1 = g11[0].length;
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    float[] es = new float[]{0.001f,1f};
    EigenTensors2 ets = new EigenTensors2(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g11i = g11[i2][i1];
      float g12i = g12[i2][i1];
      float g22i = g22[i2][i1];
      a[0][0] = g11i;
      a[0][1] = g12i;
      a[1][0] = g12i;
      a[1][1] = g22i;
      Eigen.solveSymmetric22(a,z,e);
      float u1i = z[0][0];
      float u2i = z[0][1];
      ets.setEigenvalues(i1,i2,es);
      ets.setEigenvectorU(i1,i2,u1i,u2i);
    }}
    return ets;
  }

  private EigenTensors3 updateTensors(
    final float[][][] g11, final float[][][] g12, final float[][][] g13, 
    final float[][][] g22, final float[][][] g23, final float[][][] g33) 
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    final float[] es = new float[]{0.001f,1f,1f};
    final EigenTensors3 ets = new EigenTensors3(n1,n2,n3,true);
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
            Eigen.solveSymmetric33Fast(a,z,e);
            float u1i = (float)z[0][0];
            float u2i = (float)z[0][1];
            float u3i = (float)z[0][2];
            float w1i = (float)z[2][0];
            float w2i = (float)z[2][1];
            float w3i = (float)z[2][2];
            ets.setEigenvalues(i1,i2,i3,es);
            ets.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
            ets.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
          }
        }
      }
    });
    return ets;
  }

  private float[][][] salientMap(
    float[][] g11, float[][] g12, float[][] g22) 
  {
    final int n2 = g11.length;
    final int n1 = g11[0].length;
    final float[][] sm = new float[n2][n1];
    final float[][] u1 = new float[n2][n1];
    final float[][] u2 = new float[n2][n1];
    double[][] a = new double[2][2];
    double[][] z = new double[2][2];
    double[] e = new double[2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g11i = g11[i2][i1];
      float g12i = g12[i2][i1];
      float g22i = g22[i2][i1];
      a[0][0] = g11i;
      a[0][1] = g12i;
      a[1][0] = g12i;
      a[1][1] = g22i;
      Eigen.solveSymmetric22(a,z,e);
      u1[i2][i1] = (float)z[0][0];
      u2[i2][i1] = (float)z[0][1];
      sm[i2][i1] = (float)(e[0]-e[1]);
    }}
    div(sm,max(sm),sm);
    return new float[][][]{sm,u1,u2}; 
  }


  private float[][][][] salientMap(
    final float[][][] g11, final float[][][] g12, final float[][][] g13, 
    final float[][][] g22, final float[][][] g23, final float[][][] g33) 
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    final float[][][] sm = new float[n3][n2][n1];
    final float[][][] u1 = new float[n3][n2][n1];
    final float[][][] u2 = new float[n3][n2][n1];
    final float[][][] u3 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply000(g11,g11);
    rgf.apply000(g12,g12);
    rgf.apply000(g13,g13);
    rgf.apply000(g22,g22);
    rgf.apply000(g23,g23);
    rgf.apply000(g33,g33);
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
        u1[i3][i2][i1] = u1i;
        u2[i3][i2][i1] = u2i;
        u3[i3][i2][i1] = u3i;
        sm[i3][i2][i1] = (float)(e[0]-e[1]);
      }}
    }});
    div(sm,max(sm),sm);
    return new float[][][][]{sm,u1,u2,u3}; 
  }

  private void applyLaplacian(
    EigenTensors2 et, float s,float[][] x, float[][] y){
    int n2 = x.length;
    int n1 = x[0].length;
    float[] ds = fillfloat(1.0f,3);
    ds[0] = 1.0f;
    ds[1] = 0.0f;
    ds[2] = 1.0f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if(et!=null){et.getTensor(i1,i2,ds);}
        float d11 = ds[0]*s;
        float d12 = ds[1]*s;
        float d22 = ds[2]*s;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private void applyLaplacianX(
    EigenTensors2 et, float s, float[][] f, float[][] g){
    int n2 = f.length;
    int n1 = f[0].length;
    float[] ds = new float[3];
    for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
      for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
      if(et!=null){et.getTensor(i1p,i2p,ds);}
      float d11 = ds[0]*s;
      float d12 = ds[1]*s;
      float d22 = ds[2]*s;
      float a = 0.5f*d11;
      float b = 0.5f*d12;
      float c = 0.5f*d22;
      float t = 2.0f*(a+c)/12f;
      float fpp = f[i2p][i1p];
      float fpm = f[i2p][i1m];
      float fmp = f[i2m][i1p];
      float fmm = f[i2m][i1m];
      float apppm = (a-t)*(fpp-fpm);
      float ampmm = (a-t)*(fmp-fmm);
      float bppmm = (b+t)*(fpp-fmm);
      float bpmmp = (b-t)*(fpm-fmp);
      float cppmp = (c-t)*(fpp-fmp);
      float cpmmm = (c-t)*(fpm-fmm);
      g[i2p][i1p] += apppm+bppmm+cppmp;
      g[i2p][i1m] -= apppm+bpmmp-cpmmm;
      g[i2m][i1p] += ampmm+bpmmp-cppmp;
      g[i2m][i1m] -= ampmm+bppmm+cpmmm;
    }}
  }

  private static void applyLaplacian(
    final EigenTensors3 d, final float c,
    final float[][][] x, final float[][][] y)
  { 
    final int n3 = y.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,c,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,c,x,y);
    }});
  }

  // 3D LHS
  private static void applyLhsSlice3(
    int i3, EigenTensors3 d, float c, float[][][] x, float[][][] y)
  {
    int n2 = y[0].length;
    int n1 = y[0][0].length;
    float[] di = new float[6];
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        d.getTensor(i1,i2,i3,di);
        float d11 = di[0]*c;
        float d12 = di[1]*c;
        float d13 = di[2]*c;
        float d22 = di[3]*c;
        float d23 = di[4]*c;
        float d33 = di[5]*c;
        float xa = 0.0f;
        float xb = 0.0f;
        float xc = 0.0f;
        float xd = 0.0f;

        xa += x00[i1 ];
        xd -= x00[i1m];
        xb += x01[i1 ];
        xc -= x01[i1m];
        xc += x10[i1 ];
        xb -= x10[i1m];
        xd += x11[i1 ];
        xa -= x11[i1m];

        float x1 = 0.25f*(xa+xb+xc+xd);
        float x2 = 0.25f*(xa-xb+xc-xd);
        float x3 = 0.25f*(xa+xb-xc-xd);

        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2+d23*x3;
        float y3 = d13*x1+d23*x2+d33*x3;

        float ya = 0.25f*(y1+y2+y3);
        float yb = 0.25f*(y1-y2+y3);
        float yc = 0.25f*(y1+y2-y3);
        float yd = 0.25f*(y1-y2-y3);

        y00[i1 ] += ya;
        y00[i1m] -= yd;
        y01[i1 ] += yb;
        y01[i1m] -= yc;
        y10[i1 ] += yc;
        y10[i1m] -= yb;
        y11[i1 ] += yd;
        y11[i1m] -= ya;  
      }
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private float _t = 10f; //stop time
  private int _m = 5; //number of cycles
  private float _tm = 0.5f; //stability limit
}
