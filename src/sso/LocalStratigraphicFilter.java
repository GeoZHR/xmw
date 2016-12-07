/****************************************************************************
Copyright 2007, Colorado School of Mines and others.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****************************************************************************/
package sso;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

import ad.FastExplicitDiffusion;

/**
 * Estimating seismic stratigraphic orientations.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.08.05
 */
public class LocalStratigraphicFilter {

  public LocalStratigraphicFilter(
    float[][][] u1, float[][][] u2, float[][][] u3) { 
    _u1 = u1;
    _u2 = u2;
    _u3 = u3;
  }

  public void setGradientSmoothing(double scale) {
    _scaleG = (float)scale;
  }

  public void applyForStratigraphy(
    float scale, float au, float av, float aw,
    float[][][] x, float[][][] w2, float[][][] w3, float[][][] el) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    _au = au;
    _av = av;
    _aw = aw;
    _scale = scale;
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] xs = new float[n3][n2][n1];
    EigenTensors3 ets = getTensors();
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setCycles(3,0.1f);
    if (_scaleG>0.0f) {
      ets.setEigenvalues(0.001f,1.0f,1.0f);
      xs = fed.apply(_scaleG,ets,x);
      computeOrientGradient(xs,g2,g3);
    } else {
      computeOrientGradient(x,g2,g3);
    }
    float[][][] g22 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    float[][][] g33 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g2i = g2[i3][i2][i1];
      float g3i = g3[i3][i2][i1];
      g22[i3][i2][i1] = g2i*g2i;
      g23[i3][i2][i1] = g2i*g3i;
      g33[i3][i2][i1] = g3i*g3i;
    }}}

    // Smoothed gradient products comprise the structure tensor.
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g22,g33,g23};
    ets.setEigenvalues(_au,_av,_aw);
    for (float[][][] g:gs) {
      smoothS(g,h);
      h = fed.apply(_scale,ets,h);
      copy(h,g);
    }
    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = g22[i3][i2][i1];
        a[0][1] = g23[i3][i2][i1];
        a[1][0] = g23[i3][i2][i1];
        a[1][1] = g33[i3][i2][i1];
        Eigen.solveSymmetric22(a,z,e);
        float eui = e[0];
        float evi = e[1];
        if (evi<0.0f) evi = 0.0f;
        if (eui<evi) eui = evi;
        float x2i = z[1][0];
        float x3i = z[1][1];
        if (x2i<0.0f) {
          x2i = -x2i;
          x3i = -x3i;
        }
        float u1i = _u1[i3][i2][i1];
        float u2i = _u2[i3][i2][i1];
        float u3i = _u3[i3][i2][i1];
        if (u1i<0f) {
          u1i = -u1i; u2i = -u2i; u3i = -u3i;
        }
        float u12s = u1i*u1i+u2i*u2i;
        float u23i = u2i*u3i;
        float u13i = u1i*u3i;
        float u3r = sqrt(u12s*u12s+u23i*u23i+u13i*u13i);
        float w2i =-(u23i)/u3r;
        float w3i =   u12s/u3r;
        float wsi = 1f/sqrt(w2i*w2i+w3i*w3i);
        w2i *= wsi; w3i *= wsi;
        float a2i = x2i+w2i*x3i;
        float a3i = w3i*x3i;
        if (a2i<0.0f) {
          a2i = -a2i;
          a3i = -a3i;
        }
        w2[i3][i2][i1] = a2i;
        w3[i3][i2][i1] = a3i;
        el[i3][i2][i1] = (eui-evi)/eui;
      }
    }}
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  public void computeOrientGradient(
    final float[][][] fx, final float[][][] g2, final float[][][] g3) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] g2i = g2[i3][i2];
        float[] g3i = g3[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          float u1i = _u1[i3][i2][i1];
          float u2i = _u2[i3][i2][i1];
          float u3i = _u3[i3][i2][i1];
          if (u1i<0f) {
            u1i = -u1i;
            u2i = -u2i;
            u3i = -u3i;
          }
          float u12s = u1i*u1i+u2i*u2i;
          float u12r = sqrt(u12s);
          float u23i = u2i*u3i;
          float u13i = u1i*u3i;
          float u3r = sqrt(u12s*u12s+u23i*u23i+u13i*u13i);
          float v1i = -u2i/u12r;
          float v2i =  u1i/u12r;
          float w1i =-(u13i)/u3r;
          float w2i =-(u23i)/u3r;
          float w3i = u12s/u3r;

          float v1p = i1+v1i;
          float v2p = i2+v2i;
          float v1m = i1-v1i;
          float v2m = i2-v2i;

          float w1p = i1+w1i;
          float w2p = i2+w2i;
          float w3p = i3+w3i;
          float w1m = i1-w1i;
          float w2m = i2-w2i;
          float w3m = i3-w3i;

          float gvp = si.interpolate(s1,s2,s3,fx,v1p,v2p,i3);
          float gvm = si.interpolate(s1,s2,s3,fx,v1m,v2m,i3);

          float gwp = si.interpolate(s1,s2,s3,fx,w1p,w2p,w3p);
          float gwm = si.interpolate(s1,s2,s3,fx,w1m,w2m,w3m);

          g2i[i1] = gvp-gvm;
          g3i[i1] = gwp-gwm;
        }
      }
    }});
  }

  private EigenTensors3 getTensors() {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    EigenTensors3 ets = new EigenTensors3(n1,n2,n3,true);
    for (int i3=0; i3<n3;++i3) {
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      float u1i = _u1[i3][i2][i1];
      float u2i = _u2[i3][i2][i1];
      float u3i = _u3[i3][i2][i1];
      if (u1i<0f) {
        u1i = -u1i; u2i = -u2i; u3i = -u3i;
      }
      float u12s = u1i*u1i+u2i*u2i;
      float u23i = u2i*u3i;
      float u13i = u1i*u3i;
      float u3r = sqrt(u12s*u12s+u23i*u23i+u13i*u13i);
      float w1i =-(u13i)/u3r;
      float w2i =-(u23i)/u3r;
      float w3i = u12s/u3r;
      ets.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
      ets.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
    }}}
    return ets;
  }

  /*
   * Computes y = S'Sx. Arrays x and y may be the same array.
   */
  private static void smoothS(float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    float[][][] t = new float[3][n2][n1];
    copy(x[0],t[0]);
    copy(x[0],t[1]);
    for (int i3=0; i3<n3; ++i3) {
      int i3m = (i3>0)?i3-1:0;
      int i3p = (i3<n3m)?i3+1:n3m;
      int j3m = i3m%3;
      int j3  = i3%3;
      int j3p = i3p%3;
      copy(x[i3p],t[j3p]);
      float[][] x3m = t[j3m];
      float[][] x3p = t[j3p];
      float[][] x30 = t[j3];
      float[][] y30 = y[i3];
      for (int i2=0; i2<n2; ++i2) {
        int i2m = (i2>0)?i2-1:0;
        int i2p = (i2<n2m)?i2+1:n2m;
        float[] x3m2m = x3m[i2m];
        float[] x3m20 = x3m[i2 ];
        float[] x3m2p = x3m[i2p];
        float[] x302m = x30[i2m];
        float[] x3020 = x30[i2 ];
        float[] x302p = x30[i2p];
        float[] x3p2m = x3p[i2m];
        float[] x3p20 = x3p[i2 ];
        float[] x3p2p = x3p[i2p];
        float[] y3020 = y30[i2 ];
        for (int i1=0; i1<n1; ++i1) {
          int i1m = (i1>0)?i1-1:0;
          int i1p = (i1<n1m)?i1+1:n1m;
          y3020[i1] = 0.125000f*(x3020[i1 ]) +
                      0.062500f*(x3020[i1m]+x3020[i1p]+
                                 x302m[i1 ]+x302p[i1 ]+
                                 x3m20[i1 ]+x3p20[i1 ]) +
                      0.031250f*(x3m20[i1m]+x3m20[i1p]+
                                 x3m2m[i1 ]+x3m2p[i1 ]+
                                 x302m[i1m]+x302m[i1p]+
                                 x302p[i1m]+x302p[i1p]+
                                 x3p20[i1m]+x3p20[i1p]+
                                 x3p2m[i1 ]+x3p2p[i1 ]) +
                      0.015625f*(x3m2m[i1m]+x3m2m[i1p]+
                                 x3m2p[i1m]+x3m2p[i1p]+
                                 x3p2m[i1m]+x3p2m[i1p]+
                                 x3p2p[i1m]+x3p2p[i1p]);
        }
      }
    }
  }

  private float _au=1.00f;
  private float _av=0.05f;
  private float _aw=1.00f;
  private float _scaleG = 0.0f;
  private float _scale  = 0.0f;
  private float[][][] _u1 = null;
  private float[][][] _u2 = null;
  private float[][][] _u3 = null;
}
