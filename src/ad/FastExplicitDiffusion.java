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

/**
 * Fast explicit diffusion filter. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.05
 */

public class FastExplicitDiffusion {

  public void setParameters(float t, int m, float tm) {
    _t = t;
    _m = m;
    _tm = tm;
  }

  public float[][] apply(EigenTensors2 et, float[][] fx) {
    return apply(et,null,fx);
  }

  public float[][] apply(
    EigenTensors2 et, float[][] wp, float[][] fx) 
  {
    float[][] gx = copy(fx);
    FedStep fs = new FedStep(_t,_m,_tm);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      applyLaplacianX(et,-ts[ic],wp,copy(gx),gx);
    }}
    return gx;
  }

  public float[][][] apply( EigenTensors3 et, float[][][] fx) {
    return apply(et,null,fx);
  }

  public float[][][] applyX(
    final float lambda,
    final EigenTensors3 et, final float[][][] fx) 
  {
    final float cp = 3.31488f;
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] gx = copy(fx);
    final float[][][] wp = new float[n3][n2][n1];
    final float[][][] g1 = new float[n3][n2][n1];
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    FedStep fs = new FedStep(_t,_m,_tm);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    for (int m=0; m<_m; ++m) {
      trace("m="+m);
    for (int ic=0; ic<nc; ++ic) {
      rgf.apply100(gx,g1);
      rgf.apply010(gx,g2);
      rgf.apply001(gx,g3);
      Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float g1i = g1[i3][i2][i1];
          float g2i = g2[i3][i2][i1];
          float g3i = g3[i3][i2][i1];
          float gsi = g1i*g1i+g2i*g2i+g3i*g3i;
          if (gsi==0f)
            wp[i3][i2][i1] = 1f;
          else
            wp[i3][i2][i1] = 1-exp(-cp/pow(gsi/lambda,4));
        }}
      }});
      applyLaplacian(et,-ts[ic],wp,copy(gx),gx);
    }}
    return gx;
  }

  public float[][][] apply(
    EigenTensors3 et, float[][][] wp, float[][][] fx) 
  {
    float[][][] gx = copy(fx);
    FedStep fs = new FedStep(_t,_m,_tm);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    trace("nc="+nc);
    for (int m=0; m<_m; ++m) {
      trace("m="+m);
    for (int ic=0; ic<nc; ++ic) {
      applyLaplacian(et,-ts[ic],wp,copy(gx),gx);
    }}
    return gx;
  }

  private void applyLaplacianX(
    EigenTensors2 et, float s, float[][] w, float[][] f, float[][] g){
    int n2 = f.length;
    int n1 = f[0].length;
    float[] ds = new float[3];
    for (int i2m=0,i2p=1; i2p<n2; ++i2m,++i2p) {
      for (int i1m=0,i1p=1; i1p<n1; ++i1m,++i1p) {
      if(et!=null){et.getTensor(i1p,i2p,ds);}
      float wpi = (w!=null)?w[i2p][i1p]:s;
      float d11 = ds[0]*wpi;
      float d12 = ds[1]*wpi;
      float d22 = ds[2]*wpi;
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


  private void applyLaplacian(
    EigenTensors2 et, float[][] wp, float[][] x, float[][] y){
    int n2 = x.length;
    int n1 = x[0].length;
    float[] ds = fillfloat(1.0f,3);
    ds[0] = 1.0f;
    ds[1] = 0.0f;
    ds[2] = 1.0f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if(et!=null){et.getTensor(i1,i2,ds);}
        float wpi = wp[i2][i1];
        float d11 = ds[0]*wpi;
        float d12 = ds[1]*wpi;
        float d22 = ds[2]*wpi;
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

  private static void applyLaplacian(
    final EigenTensors3 d, final float c, 
    final float[][][] w, final float[][][] x, final float[][][] y)
  { 
    final int n3 = y.length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,c,w,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,d,c,w,x,y);
    }});
  }

  // 3D LHS
  private static void applyLhsSlice3(
    int i3, EigenTensors3 d, float c, float[][][] w, float[][][] x, float[][][] y)
  {
    int n2 = y[0].length;
    int n1 = y[0][0].length;
    float[] di = fillfloat(1.0f,6);
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
        if(d!=null){d.getTensor(i1,i2,i3,di);}
        float wpi = (w!=null)?w[i3][i2][i1]:c;
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

        float x1 = 0.25f*(xa+xb+xc+xd)*wpi;
        float x2 = 0.25f*(xa-xb+xc-xd)*wpi;
        float x3 = 0.25f*(xa+xb-xc-xd)*wpi;

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
