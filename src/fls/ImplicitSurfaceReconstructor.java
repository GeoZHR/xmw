/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.*;
import vec.*;
import pik.*;

/**
 * 2D fast level set method.
 * <p>
 * Based on the works by Yonggang Shi and William Clem Karl, 2008, 
 * A Real-Time Algorithm for the Approximation of Level-Set-Based 
 * Curve Evolution.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class ImplicitSurfaceReconstructor {

  public void testSquare(float[][][] wx, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;

    int ri = 100;
    int c1 = round(n1/2f);
    int c2 = round(n2/2f);
    int c3 = round(n3/2f);
    int b1 = c1-ri;
    int e1 = c1+ri;

    int b2 = c2-ri;
    int e2 = c2+ri;
    int b3 = c3-ri;
    int e3 = c3+ri;
    int di = 5;
    for (int i3=b3; i3<e3; i3+=di) {
    for (int i2=b2; i2<e2; i2+=di) {
      fx[i3][i2][b1-1] = 1;
      fx[i3][i2][b1  ] = 0;
      fx[i3][i2][b1+1] =-1;
      fx[i3][i2][e1-1] =-1;
      fx[i3][i2][e1  ] = 0;
      fx[i3][i2][e1+1] = 1;

      wx[i3][i2][b1-1] = 1;
      wx[i3][i2][b1  ] = 1;
      wx[i3][i2][b1+1] = 1;
      wx[i3][i2][e1-1] = 1;
      wx[i3][i2][e1  ] = 1;
      wx[i3][i2][e1+1] = 1;

    }}
    for (int i3=b3; i3<e3; i3+=di) {
    for (int i1=b1; i1<e1; i1+=di) {
      fx[i3][b2-1][i1] = 1;
      fx[i3][b2  ][i1] = 0;
      fx[i3][b2+1][i1] =-1;
      fx[i3][e2-1][i1] =-1;
      fx[i3][e2  ][i1] = 0;
      fx[i3][e2+1][i1] = 1;
      wx[i3][b2-1][i1] = 1;
      wx[i3][b2  ][i1] = 1;
      wx[i3][b2+1][i1] = 1;
      wx[i3][e2-1][i1] = 1;
      wx[i3][e2  ][i1] = 1;
      wx[i3][e2+1][i1] = 1;

    }}
    for (int i2=b2; i2<e2; i2+=di) {
    for (int i1=b1; i1<e1; i1+=di) {
      fx[b3-1][i2][i1] = 1;
      fx[b3  ][i2][i1] = 0;
      fx[b3+1][i2][i1] =-1;
      fx[e3-1][i2][i1] =-1;
      fx[e3  ][i2][i1] = 0;
      fx[e3+1][i2][i1] = 1;

      wx[b3-1][i2][i1] = 1;
      wx[b3  ][i2][i1] = 1;
      wx[b3+1][i2][i1] = 1;
      wx[e3-1][i2][i1] = 1;
      wx[e3  ][i2][i1] = 1;
      wx[e3+1][i2][i1] = 1;
    }}
  }

  public void testSphere(float[][][] wx, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;

    float ri = 100f;
    float rm = ri-2;
    float rp = ri+2;
    float c1 = n1/2f;
    float c2 = n2/2f;
    float c3 = n3/2f;
    float pi = (float)(Math.PI);
    int da = 5;
    for (int i=0; i<180; i+=da) {
      float theta = i*pi/180f;
    for (int k=0; k<360; k+=da) {
      float phi = k*pi/180f;
      int i1 = round(c1+ri*cos(theta));
      int i2 = round(c2+ri*sin(theta)*sin(phi));
      int i3 = round(c3+ri*sin(theta)*cos(phi));
      int m1 = round(c1+rm*cos(theta));
      int m2 = round(c2+rm*sin(theta)*sin(phi));
      int m3 = round(c3+rm*sin(theta)*cos(phi));
      int p1 = round(c1+rp*cos(theta));
      int p2 = round(c2+rp*sin(theta)*sin(phi));
      int p3 = round(c3+rp*sin(theta)*cos(phi));
      fx[i3][i2][i1] = 0;
      fx[p3][p2][p1] = 1;
      fx[m3][m2][m1] =-1;
      wx[i3][i2][i1] = 1;
      wx[p3][p2][p1] = 1;
      wx[m3][m2][m1] = 1;
    }}
  }

  public void pointsToImage(
    float[] xs, float[][] ys, float[][] zs, float[][][] pa, float[][][] wx) 
  {
    int ns = xs.length;
    int n3 = pa.length;
    int n2 = pa[0].length;
    int n1 = pa[0][0].length;
    for (int is=0; is<ns; ++is) {
      int i3 = round(xs[is]);
      i3 = max(i3,0);
      i3 = min(i3,n3-1);
      int np = ys[is].length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = round(ys[is][ip]);
        int i1 = round(zs[is][ip]);
        i2 = max(i2,0);
        i2 = min(i2,n2-1);
        i1 = max(i1,0);
        i1 = min(i1,n1-1);
        wx[i3][i2][i1] = pa[i3][i2][i1];
      }
    }
  }

  public void signAsignment(
    float[][][] wx, float[][][] fx) {
    int n3 = wx.length;
    int n2 = wx[0].length;
    int n1 = wx[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float[] wx32 = wx[i3][i2];
      float[] fx32 = fx[i3][i2];
      float mk = 1f;
    for (int i1=0; i1<n1; ++i1) {
      if(wx32[i1]!=0f) {
        fx32[i1] = 0f;
        if(i1==0) mk*=-1f;
        if(i1>0&&wx32[i1-1]==0) mk *= -1f;
      } else {
        fx32[i1] = mk;
      }
    }}}
  }

  public float[][][] smooth(float sig1, float sig2, float sig3, 
    float[][][] wx, float[][][] fx) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] b = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    makeRhs(wx,fx,b);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    Smoother3 smoother3 = new Smoother3(sig1,sig2,sig3);
    A3 a3 = new A3(smoother3,wx);
    CgSolver cs = new CgSolver(0.001,200);
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vr);
    return r;
  }

    // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, float[][][] wp) 
    {
      _s3 = s3;
      _wp = wp;
      float n3 = wp.length;
      float n2 = wp[0].length;
      float n1 = wp[0][0].length;
      _sc = sum(wp)/(n1*n2*n3);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = copy(x);
      v3y.zero();
      _s3.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s3.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother3 _s3;
    private float[][][] _wp;
  }


  private static void addAndScale(float sc, float[][][] x, float[][][] y) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      y[i3][i2][i1] += sc*x[i3][i2][i1];
    }}}
  }


  private static void applyLhs(float[][][] wp, float[][][] x, float[][][] y) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wpi = wp[i3][i2][i1];
      y[i3][i2][i1] = x[i3][i2][i1]*wpi*wpi;
    }}}
  }



  private void makeRhs(
    float[][][] wx, float[][][] fx, float[][][] b) 
  {
    int n3 = wx.length;
    int n2 = wx[0].length;
    int n1 = wx[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wxi = wx[i3][i2][i1];
      b[i3][i2][i1] = fx[i3][i2][i1]*wxi*wxi;
    }}}
  }


    // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother3 {
    public Smoother3(float sigma1, float sigma2, float sigma3) {
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
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  // Smoothing for dimension 3.
  private static void smooth3(float sigma, float[][][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply3(x,x);
  }


}
