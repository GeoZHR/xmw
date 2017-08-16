/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipfx;

import vec.*;
import util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;
import static ifs.FaultGeometry.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.09.16
 */
public class FaultExtension {

  /**
   * Returns fault blocks for specified fault images.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   * @return array of fault blocks.
   */

  /*
  public float[][][] findBlocks(float[][][][] flpt) {
    return blocks(flpt);
  }
  */

    /**
   * @param us[0] array of fault likelihoods.
   * @param us[1] array of 1st component of fault normal vectors.
   * @param us[2] array of 2nd component of fault normal vectors.
   * @param us[3] array of 3rd component of fault normal vectors.
   */
  public float[][][] faultExtension(
    float[][][] u1, float[][][] u2, float[][][] u3)
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float[][][] b = new float[n3][n2][n1]; // right-hand side
    float[][][] f = new float[n3][n2][n1]; // fault isosurface volume, in samples
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vf = new VecArrayFloat3(f);
    Smoother3 smoother3 = new Smoother3(n1,n2,n3,_sigma1,_sigma2,_sigma3);
    A3 a3 = new A3(smoother3);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(u1,u2,u3,b);
    System.out.println("bmax="+max(b));
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vf);
    smoother3.apply(f);
    return f;
  }

  public float[][][][] getSlipVectors(int n1, int n2, int n3, FaultSkin[] skins) {
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    for (FaultSkin skin: skins) {
    for (FaultCell cell: skin) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      float c2 = cell.getS2();
      float c3 = cell.getS3();
      if(c2<0f) {
        c2 = -c2;
        c3 = -c3;
      }
      w2[i3][i2][i1] = c2;
      w3[i3][i2][i1] = c3;
    }}
    return new float[][][][]{w1,w2,w3};
  }

    // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3){
      _s3 = s3;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      _s3.apply(z);
      zero(y);
      applyLhs(z,y);       //laplacian operator
      _s3.applyTranspose(y);
    }
    private Smoother3 _s3;
  }

    // right-hand side for 3D
  private static void makeRhs(
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1,i3m=0; i3<n3; ++i3,++i3m) {
    for (int i2=1,i2m=0; i2<n2; ++i2,++i2m) {
    for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
      float y1 = u1[i3][i2][i1];
      float y2 = u2[i3][i2][i1];
      float y3 = u3[i3][i2][i1];
      float ya = 0.25f*(y1+y2+y3);
      float yb = 0.25f*(y1-y2+y3);
      float yc = 0.25f*(y1+y2-y3);
      float yd = 0.25f*(y1-y2-y3);
      y[i3 ][i2 ][i1 ] += ya;
      y[i3 ][i2 ][i1m] -= yd;
      y[i3 ][i2m][i1 ] += yb;
      y[i3 ][i2m][i1m] -= yc;
      y[i3m][i2 ][i1 ] += yc;
      y[i3m][i2 ][i1m] -= yb;
      y[i3m][i2m][i1 ] += yd;
      y[i3m][i2m][i1m] -= ya;
    }}}
  }


  // left-hand side for 3D
  private static void applyLhs(
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() { // i3 = 1, 3, 5, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() { // i3 = 2, 4, 6, ...
    public void compute(int i3) {
      applyLhsSlice3(i3,x,y);
    }});
  }

  private static void applyLhsSlice3(
    int i3, float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
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

        float y1 = 0.25f*(xa+xb+xc+xd);
        float y2 = 0.25f*(xa-xb+xc-xd);
        float y3 = 0.25f*(xa+xb-xc-xd);

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

  // Smoother used as a preconditioner.
  private static class Smoother3 {
    public Smoother3(int n1, int n2, int n3, 
      float sigma1, float sigma2, float sigma3) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
    }
    public void apply(float[][][] x) {
      smooth3(_sigma3,x);
      smooth2(_sigma2,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,x);
      smooth3(_sigma3,x);
    }
    private float _sigma1,_sigma2,_sigma3;
  }


  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 3.
  private static void smooth3(float sigma, float[][][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply3(x,x);
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  //
  /*
  private static float[][][] extension(float[][][][] fls) {
    float[][][] fl = fls[0];
    float[][][] s1 = fls[0];
    float[][][] s2 = fls[1];
    float[][][] s3 = fls[2];
    int n3 = s1.length;
    int n2 = s1[0].length;
    int n1 = s1[0][0].length;

    // Compute right-hand-side.
    float[][][] r = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      int m3 = (i3>0)?i3-1:0;
      for (int i2=0; i2<n2; ++i2) {
        int m2 = (i2>0)?i2-1:0;
        for (int i1=0; i1<n1; ++i1) {
          int m1 = (i1>0)?i1-1:0;
          float fs1 = s1[i3][i2][i1];
          float fs2 = s2[i3][i2][i1];
          float fs3 = s3[i3][i2][i1];
          float fl1 = 1f;//0.5f*(fl[i3][i2][i1]+fl[i3][i2][m1]);
          float fl2 = 1f;//0.5f*(fl[i3][i2][i1]+fl[i3][m2][i1]);
          float fl3 = 1f;//0.5f*(fl[i3][i2][i1]+fl[m3][i2][i1]);
          float f1 = fl1*fs1;
          float f2 = fl2*fs2;
          float f3 = fl3*fs3;
          r[i3][i2][i1] += f1+f2+f3;
          r[i3][i2][m1] -= f1;
          r[i3][m2][i1] -= f2;
          r[m3][i2][i1] -= f3;
        }
      }
    }

    // Solve for fault blocks.
    A3 a = new A3();
    CgSolver cs = new CgSolver(0.01,100);
    float[][][] b = new float[n3][n2][n1];
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    cs.solve(a,vr,vb);
    return b;
  }


  // Returns fault blocks.
  private static float[][][] blocks(float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;

    // Compute right-hand-side.
    float[][][] r = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      int m3 = (i3>0)?i3-1:0;
      for (int i2=0; i2<n2; ++i2) {
        int m2 = (i2>0)?i2-1:0;
        for (int i1=0; i1<n1; ++i1) {
          int m1 = (i1>0)?i1-1:0;
          float fl = f[i3][i2][i1];
          float fp = p[i3][i2][i1];
          float ft = t[i3][i2][i1];
          float[] fn = faultNormalVectorFromStrikeAndDip(fp,ft);
          float fn1 = fn[0];
          float fn2 = fn[1];
          float fn3 = fn[2];
          float fl1 = 0.5f*(f[i3][i2][i1]+f[i3][i2][m1]);
          float fl2 = 0.5f*(f[i3][i2][i1]+f[i3][m2][i1]);
          float fl3 = 0.5f*(f[i3][i2][i1]+f[m3][i2][i1]);
          float f1 = fl1*fn1;
          float f2 = fl2*fn2;
          float f3 = fl3*fn3;
          r[i3][i2][i1] += f1+f2+f3;
          r[i3][i2][m1] -= f1;
          r[i3][m2][i1] -= f2;
          r[m3][i2][i1] -= f3;
        }
      }
    }

    // Solve for fault blocks.
    A3 a = new A3();
    CgSolver cs = new CgSolver(0.01,100);
    float[][][] b = new float[n3][n2][n1];
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    cs.solve(a,vr,vb);
    return b;
  }
  private static class A3 implements CgSolver.A {
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      zero(y);
      _ldk.apply(x,y);
    }
    private LocalDiffusionKernel _ldk =
      new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21);
  }
  */


  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _sigma3 = 6.0f; // precon smoothing extent for 3rd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 100; // maximum number of CG iterations

}
