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
 * Compute an implicit function for 3D salt boundary reconstruction.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2017.03.06
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

  public void getImplicitFunc(
    float sig, float[][] ds, float[][][][] wxs, 
    float[][][][] ws, float[][][][] fs) {
    int n3 = wxs[0].length;
    for (int i3=0; i3<n3; ++i3) {
      float[][][] fw0 = getImplicitFuncAndWeights(wxs[0][i3]);
      float[][][] fw1 = getImplicitFuncAndWeights(wxs[0][i3]);
      fs[0][i3] = fw0[0];
      fs[1][i3] = fw1[0];
      ws[0][i3] = mul(fw0[1],gauss(sig,ds[0][i3]));
      ws[1][i3] = mul(fw1[1],gauss(sig,ds[1][i3]));
    }
  }

  private float gauss(float sig, float dx) {
    return exp(-dx*dx*0.5f/(sig*sig));
  }


  public float[][][] getImplicitFuncAndWeights(float[][] wx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[][] pf = new float[n2][n1];
    float[][] ws = new float[n2][n1];
    float[][] dx = new float[n2][n1];
    short[][] k1 = new short[n2][n1];
    short[][] k2 = new short[n2][n1];
    float[][] sg = new float[n2][n1];
    ClosestPointTransform cpt = new ClosestPointTransform(1,1);
    cpt.apply(0.0f,wx,dx,k1,k2);
    signAsignmentH(wx,sg);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int c1 = k1[i2][i1];
      int c2 = k2[i2][i1];
      ws[i2][i1] = wx[c2][c1];
      pf[i2][i1] = sg[i2][i1]*dx[i2][i1];
    }}
    return new float[][][] {pf,ws};
  }



  public float[][][] smoothFit(float sig1, float sig2, float sig3, 
    float[][][][] ws, float[][][][] fs) {
    int n3 = ws[0].length;
    int n2 = ws[0][0].length;
    int n1 = ws[0][0][0].length;
    float[][][] b = new float[n3][n2][n1];
    float[][][] r = new float[n3][n2][n1];
    makeRhs(ws,fs,b);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    Smoother3 smoother3 = new Smoother3(sig1,sig2,sig3);
    A3 a3 = new A3(smoother3,ws);
    CgSolver cs = new CgSolver(0.001,200);
    smoother3.applyTranspose(b);
    cs.solve(a3,vb,vr);
    smoother3.apply(r);
    return r;
  }

    // Conjugate-gradient operators.
  private static class A3 implements CgSolver.A {
    A3(Smoother3 s3, float[][][][] wp) 
    {
      _s3 = s3;
      _wp = wp;
      float n4 = wp.length;
      float n3 = wp[0].length;
      float n2 = wp[0][0].length;
      float n1 = wp[0][0][0].length;
      _sc = 4f*(sum(wp[0])+sum(wp[1]))/(n1*n2*n3*n4);
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
    private float[][][][] _wp;
  }


  private static void applyLhs(float[][][][] ws, float[][][] x, float[][][] y) {
    int n3 = ws[0].length;
    int n2 = ws[0][0].length;
    int n1 = ws[0][0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float ws1 = ws[0][i3][i2][i1];
      float ws2 = ws[1][i3][i2][i1];
      y[i3][i2][i1] += (ws1+ws2)*x[i3][i2][i1];
    }}}

  }


  private void makeRhs(float[][][][] ws, float[][][][] fs, float[][][] b) {
    int n3 = ws[0].length;
    int n2 = ws[0][0].length;
    int n1 = ws[0][0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fs1 = fs[0][i3][i2][i1];
      float ws1 = ws[0][i3][i2][i1];
      float fs2 = fs[1][i3][i2][i1];
      float ws2 = ws[1][i3][i2][i1];
      b[i3][i2][i1] = fs1*ws1+fs2*ws2;
    }}}

  }

  public float[][] pointsToArray(
    float[] ds, float[] xs, float[][] ys, float[][] zs) {
    int np = 0;
    int ns = xs.length;
    for (int is=0; is<ns; ++is)
      np += ys[is].length;
    float[][] ps = new float[4][np];
    int ip = 0;
    for (int is=0; is<ns; ++is) {
      int nc = ys[is].length;
      for (int ic=0; ic<nc; ++ic) {
        ps[3][ip] = ds[is];
        ps[2][ip] = xs[is];
        ps[1][ip] = ys[is][ic];
        ps[0][ip] = zs[is][ic];
      }
    }
    return ps;
  }

  public void pointsToImage(
    float[] ys, float[] zs, float[][] pa, float[][] wx) 
  {
    int n2 = pa.length;
    int n1 = pa[0].length;
    int np = ys.length;
    for (int ip=0; ip<np; ++ip) {
      int i2 = round(ys[ip]);
      int i1 = round(zs[ip]);
      i2 = max(i2,0);
      i2 = min(i2,n2-1);
      i1 = max(i1,0);
      i1 = min(i1,n1-1);
      wx[i2][i1] = pa[i2][i1];
    }
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

  public void signAsignmentH(
    float[][] wx, float[][] fx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    for (int i1=0; i1<n1; ++i1) {
      float mk = 1f;
    for (int i2=0; i2<n2; ++i2) {
      if(wx[i2][i1]!=0f) {
        fx[i2][i1] = 0f;
        if(i2==0) mk*=-1f;
        if(i2>0&&wx[i2-1][i1]==0) mk *= -1f;
      } else {
        fx[i2][i1] = mk;
      }
    }}
    //despike
    float[][] fs = zerofloat(n1,n2);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(10);
    ref.apply(fx,fs);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(fx[i2][i1]!=0) {
        if(fs[i2][i1]>0) fx[i2][i1] =  1f;
        if(fs[i2][i1]<0) fx[i2][i1] = -1f;
      }
    }}
  }

  public void signAsignmentV(
    float[][] wx, float[][] fx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    for (int i2=0; i2<n2; ++i2) {
      float[] wx2 = wx[i2];
      float[] fx2 = fx[i2];
      float mk = 1f;
    for (int i1=0; i1<n1; ++i1) {
      if(wx2[i1]!=0f) {
        fx2[i1] = 0f;
        if(i1==0) mk*=-1f;
        if(i1>0&&wx2[i1-1]==0) mk *= -1f;
      } else {
        fx2[i1] = mk;
      }
    }}

  }

  public void signAsignmentV(
    float[][][] wx, float[][][] fx) {
    int n3 = wx.length;
    for (int i3=0; i3<n3; ++i3)
      signAsignmentV(wx[i3],fx[i3]);
  }

  public void signAsignmentH(
    float[][][] wx, float[][][] fx) {
    int n3 = wx.length;
    for (int i3=0; i3<n3; ++i3)
      signAsignmentH(wx[i3],fx[i3]);
  }


  /*
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
    smoother3.apply(r);
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
  */


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
