/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Nonlinear structure tensors.
 * <p>
 * Based on the works by Yonggang Shi and William Clem Karl, 2008, 
 * A Real-Time Algorithm for the Approximation of Level-Set-Based 
 * Curve Evolution.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class NonlinearStructureTensor {
  /**
   * Construct with a 2D initial square-shape level-set function 
   * with integes -3,-1,1,3.
   * @param n1 the 1st dimension number.
   * @param n2 the 2nd dimension number.
   * @param c1 the x coordinate of the center.
   * @param c2 the y coordinate of the center.
   * @param r  the radius.
   */
  public NonlinearStructureTensor(float p, float t, float d) {
    _p = p;
    _t = t;
    _d = d;
  }

  public void applyForNormalLinear(
    float[][] fx, float[][] u1, float[][] u2, float[][] el) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g11 = new float[n2][n1];
    float[][] g22 = new float[n2][n1];
    float[][] g12 = new float[n2][n1];
    computeGradientProduct(fx,g11,g12,g22);
    nonlinearDiffusion(g11,g12,g22);
    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = g11[i2][i1];
        a[0][1] = g12[i2][i1];
        a[1][0] = g12[i2][i1];
        a[1][1] = g22[i2][i1];
        Eigen.solveSymmetric22(a,z,e);
        float u1i = z[0][0];
        float u2i = z[0][1];
        if (u1i<0.0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        float eui = e[0];
        float evi = e[1];
        if (evi<0.0f) evi = 0.0f;
        if (eui<evi) eui = evi;
        if (u1!=null) u1[i2][i1] = u1i;
        if (u2!=null) u2[i2][i1] = u2i;
        if (el!=null) el[i2][i1] = (eui-evi)/eui;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _p; // edge enhancing parameter
  private float _t; // diffustion time 
  private float _d; // diffustion time 

  private void computeGradientProduct(
    float[][] fx, float[][] g11, float[][] g12, float[][] g22) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply10(fx,g1);
    rgf.apply01(fx,g2);
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      g11[i2][i1] = g1i*g1i;
      g12[i2][i1] = g1i*g2i;
      g22[i2][i1] = g2i*g2i;
    }}
  }

  private void nonlinearDiffusion(
    final float[][] g11, final float[][] g12, final float[][] g22) {
    float nt = round(_t/_d);
    for (int i=0; i<=nt; ++i) {
      float[][] s = diffusivity(new float[][][]{g11,g12,g22});
      applyLaplacian(-_d,s,copy(g11),g11);
      applyLaplacian(-_d,s,copy(g12),g12);
      applyLaplacian(-_d,s,copy(g22),g22);
    }
  }

  public void diffusion(float[][] fx) {
    float nt = round(_t/_d);
    for (int i=0; i<=nt; ++i) {
      float[][] s = diffusivity(fx);
      applyLaplacian(_d,s,copy(fx),fx);
    }
  }


  private float[][] diffusivity(float[][] gs) {
    int n2 = gs.length;
    int n1 = gs[0].length;
    float[][] s = new float[n2][n1];
    float[][] ds = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    float[][] g1t = new float[n2][n1];
    float[][] g2t = new float[n2][n1];
    rgf.apply10(gs,g1t);
    rgf.apply01(gs,g2t);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g1i = g1t[i2][i1];
      float g2i = g2t[i2][i1];
      ds[i2][i1] = sqrt(g1i*g1i+g2i*g2i); 
    }}
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float dsi = ds[i2][i1]+0.00001f;
      s[i2][i1] = pow(dsi,-_p);
    }}
    return s;
  }


  private float[][] diffusivity(float[][][] gs) {
    int n3 = gs.length;
    int n2 = gs[0].length;
    int n1 = gs[0][0].length;
    float[][] s = new float[n2][n1];
    float[][] ds = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    for (int i3=0; i3<n3; ++i3) {
      float[][] g1t = new float[n2][n1];
      float[][] g2t = new float[n2][n1];
      rgf.apply10(gs[i3],g1t);
      rgf.apply01(gs[i3],g2t);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1t[i2][i1];
        float g2i = g2t[i2][i1];
        ds[i2][i1] += sqrt(g1i*g1i+g2i*g2i); 
      }}
    }
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float dsi = ds[i2][i1]+0.000001f;
      s[i2][i1] = pow(dsi,-_p);
    }}
    return s;
  }

  private void applyLaplacian(
    float c, float[][] s, float[][] x, float[][] y){
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float si = c*s[i2][i1]*0.25f;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = (xa+xb);
        float x2 = (xa-xb);
        float y1 = x1*si;
        float y2 = x2*si;
        float ya = (y1+y2);
        float yb = (y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

}
