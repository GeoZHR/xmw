/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ad;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import ipfx.*;
import ipfx.FaultCell;
import static ipfx.FaultGeometry.*;

/**
 * Automatic drawing. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.04
 */
public class CoherenceShockFilter {

  public void setIterations(int niter) {
    _niter = niter;
  }

  public void setSmoothings(float sigma1, float sigma2) {
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }

  public float[][] apply(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] fy = copy(fx);
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float dt = 0.5f;
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(1.0);
    //LocalOrientFilter lof = new LocalOrientFilter(_sigma1,_sigma2);
    for (int iter=0; iter<_niter; ++iter) {
      //EigenTensors2 et = lof.applyForTensors(fy);
      //et.setEigenvalues(1f,0.0001f);
      //applyLaplacian(et,fy,fd);
      float[][] fd = new float[n2][n1];
      applyLaplacian(fy,fd);
      rgf1.apply10(fy,g1);
      rgf1.apply01(fy,g2);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        float fdi = fd[i2][i1];
        float gsi = dt*sqrt(g1i*g1i+g2i*g2i);
        if(fdi<0.0f)
          fy[i2][i1] += gsi;
        else
          fy[i2][i1] -= gsi;
      }}
    }
    return fy;
  }

  public float[][][] apply(int nt, EigenTensors3 et, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = copy(fx);
    float[][][] g1 = new float[n3][n2][n1];
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] gd = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    for (int it=0; it<nt; ++it) {
      trace("it="+it);
      zero(gd);
      rgf.apply100(gx,g1);
      rgf.apply010(gx,g2);
      rgf.apply001(gx,g3);
      applyLaplacian(et,gx,gd);
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i3][i2][i1];
        float g2i = g2[i3][i2][i1];
        float g3i = g3[i3][i2][i1];
        float gsi = sqrt(g1i*g1i+g2i*g2i+g3i*g3i);
        float gdi = -gd[i3][i2][i1];
        if(gdi<0.0f)
          gx[i3][i2][i1] += gsi;
        else
          gx[i3][i2][i1] -= gsi;
      }}}
    }
    return gx;
  }

  public float[][][] apply(
    int nt, float[][][] fx, 
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = copy(fx);
    float[][][] g1 = new float[n3][n2][n1];
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] p1 = new float[n3][n2][n1];
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    float[][][] gd = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    for (int it=0; it<nt; ++it) {
      trace("it="+it);
      zero(gd);
      rgf.apply100(gx,g1);
      rgf.apply010(gx,g2);
      rgf.apply001(gx,g3);
      mul(g1,u1,p1);
      mul(g2,u2,p2);
      mul(g3,u3,p3);
      add(p1,p2,gd);
      add(gd,p3,gd);
      rgf.apply100(gd,p1);
      rgf.apply010(gd,p2);
      rgf.apply001(gd,p3);
      mul(p1,u1,p1);
      mul(p2,u2,p2);
      mul(p3,u3,p3);
      add(p1,p2,gd);
      add(gd,p3,gd);
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i3][i2][i1];
        float g2i = g2[i3][i2][i1];
        float g3i = g3[i3][i2][i1];
        float gsi = sqrt(g1i*g1i+g2i*g2i+g3i*g3i);
        float gdi = gd[i3][i2][i1];
        if(gdi<0.0f)
          gx[i3][i2][i1] += gsi;
        else
          gx[i3][i2][i1] -= gsi;
      }}}
    }
    return gx;
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
    LocalOrientFilter lof = new LocalOrientFilter(_sigma1,_sigma2);
    lof.applyForNormal(fx,u1,u2);
    rgf.apply10(fx,g1);
    rgf.apply01(fx,g2);
    mul(g1,u1,g1);
    mul(g2,u2,g2);
    add(g1,g2,fy);
    return abs(fy);
  }

  private void applyLaplacian(EigenTensors2 et, float[][] x, float[][] y){
    int n2 = x.length;
    int n1 = x[0].length;
    float[] ds = fillfloat(1.0f,3);
    ds[0] = 1.0f;
    ds[1] = 0.0f;
    ds[2] = 1.0f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if(et!=null){et.getTensor(i1,i2,ds);}
        float d11 = ds[0];
        float d12 = ds[1];
        float d22 = ds[2];
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
    mul(-1f,y,y);
  }

  private void applyLaplacian(float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    LocalOrientFilter lof = new LocalOrientFilter(_sigma1,_sigma2);
    lof.applyForNormal(x,u1,u2);
    rgf.apply10(x,g1);
    rgf.apply01(x,g2);
    mul(g1,u1,g1);
    mul(g2,u2,g2);
    add(g1,g2,y);
    rgf.apply10(y,g1);
    rgf.apply01(y,g2);
    mul(g1,u1,g1);
    mul(g2,u2,g2);
    add(g1,g2,y);
  }

  private static void applyLaplacian(
    final EigenTensors3 d,
    final float[][][] x, final float[][][] y)
  { 
    final int n3 = y.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,x,y);
    }});
  }

  // 3D LHS
  private static void applyLhsSlice3(
    int i3, EigenTensors3 d, float[][][] x, float[][][] y)
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
        float d11 = di[0];
        float d12 = di[1];
        float d13 = di[2];
        float d22 = di[3];
        float d23 = di[4];
        float d33 = di[5];
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


  private void updateDirections(
    float[][] fx, float[][] v1, float[][] v2)
  {
    LocalOrientFilter lof = new LocalOrientFilter(4.0,4.0);
    lof.apply(fx,null,null,null,v1,v2,null,null,null);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private int _niter = 5;
  private float _sigma1 = 5f;
  private float _sigma2 = 5f;
}
