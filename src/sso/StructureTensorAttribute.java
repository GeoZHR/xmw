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
 * Structure tensors for estimating local structural 
 * and stratigraphic orientations.
 * Methods of this class can compute for each image sample numerous
 * parameters related to orientation. All orientation information 
 * is derived from eigenvectors and eigenvalues of the structure tensor
 * (also called the "gradient-squared tensor"). This tensor is equivalent 
 * to a matrix of 2nd partial derivatives of an autocorrelation evaluated 
 * at zero lag. In other words, orientation is here determined by the 
 * (2-D) ellipse or (3-D) ellipsoid that best fits the peak of the 
 * autocorrelation of image samples in a local window.
 * <p>
 * The coordinate system for a 2-D image has two orthogonal axes 1 and 2, 
 * which correspond to the 1st and 2nd indices of the array containing 
 * image samples. For 2-D images, the eigenvectors are the unit vectors 
 * u = (u1,u2) and v = (v1,v2). The 1st eigenvector u is perpendicular 
 * to the best fitting line, and the 1st component u1 of u is always 
 * non-negative. The 2nd eigenvector v is perpendicular to u such that 
 * the cross product u1*v2-u2*v1 = 1; that is, v1 = -u2 and v2 = u1. 
 * The angle theta = asin(u2) is the angle measured counter-clockwise 
 * between the 1st eigenvector u and axis 1; -pi/2 &lt;= theta &lt;= pi/2.
 * <p>
 * The coordinate system for a 3-D image has three orthogonal axes 1, 2 
 * and 3, which correspond to the 1st, 2nd and 3rd indices of the array 
 * containing image samples. For 3-D images, the eigenvectors are unit 
 * vectors u = (u1,u2,u3), v = (v1,v2,v3), and w = (w1,w2,w3). The 1st 
 * eigenvector u is orthogonal to the best fitting plane, and the 1st 
 * component u1 of u is always non-negative. The 2nd eigenvector v is 
 * orthogonal to the best fitting line within the best fitting plane.
 * The 3rd eigenvector w is orthogonal to both u and v and is aligned
 * with the direction in which the images changes least. The dip angle 
 * theta = acos(u1) is the angle between the 1st eigenvector u and axis 1; 
 * 0 &lt;= theta &lt;= pi/2. The azimuthal angle phi = atan2(u3,u2)
 * is well-defined for only non-zero theta; -pi &lt;= phi &lt;= pi.
 * <p>
 * The local linearity or planarity of features is determined by the
 * eigenvalues. For 2-D images with eigenvalues eu and ev (corresponding 
 * to the eigenvectors u and v), linearity is (eu-ev)/eu. For 3-D
 * images with eigenvalues eu, ev, and ew, planarity is (eu-ev)/eu
 * and linearity is (ev-ew)/eu. Both linearity and planarity are
 * in the range [0,1].
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.08.18
 */
public class StructureTensorAttribute {

  /**
   * Constructs a filter with anisotropic structure-oriented smoothing.
   * @param scale factor for the smoothing
   */
  public StructureTensorAttribute(EigenTensors2 et, double scale) {
    _et2 = et;
    _scale = (float)scale;
    setEigenvalues(1.0f,0.05f);
  }

  /**
   * Constructs a filter with anisotropic structure-oriented smoothing.
   * @param scale factor for the smoothing
   */
  public StructureTensorAttribute(EigenTensors3 et, double scale) { 
    _et3 = et;
    _scale = (float)scale;
    setEigenvalues(1.0f,0.05f,1.0f);
    //resetTensors();
  }

  /**
   * Sets eigenvalues of structure tensors for anisotropic smoothing.
   * @param au eigenvalue corresponding to the vector normal to structures
   * @param av eigenvalue corresponding to the vector parallel to structures
   */
  public void setEigenvalues(double au, double av) {
    _au = (float)au;
    _av = (float)av;
  }

  /**
   * Sets eigenvalues of structure tensors for anisotropic smoothing.
   * @param au eigenvalue corresponding to the vector normal to structures
   * @param av eigenvalue corresponding to the vector parallel to structures 
   * but normal to stratigraphic features
   * @param aw eigenvalue corresponding to the vector parallel to structures 
   * and stratigraphic features
   */
  public void setEigenvalues(double au, double av, double aw) {
    _au = (float)au;
    _av = (float)av;
    _aw = (float)aw;
  }

  public void setGradientSmoothing(double scale) {
    _scaleG = (float)scale;
  }

  /**
   * Applies this filter for the specified image and outputs. All
   * outputs are optional and are computed for only non-null arrays.
   * @param x input array for 2-D image
   * @param el (eu-ev)/eu, a measure of linearity.
   */
  public void applyForLinear(float[][] x,float[][] el)
  {
    // Gradient.
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] g1 = el;
    float[][] g2 = new float[n2][n1]; 
    //_et2.setEigenvalues(0.001f,1.0f);
    //_lsf.apply(_et2,3,x,h);
    computeOrientGradient(x,g1,g2);

    // Gradient products.
    float[][] g11 = g1;
    float[][] g22 = g2;
    float[][] g12 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        g11[i2][i1] = g1i*g1i;
        g22[i2][i1] = g2i*g2i;
        g12[i2][i1] = g1i*g2i;
      }
    }
    
    //Smoothed gradient products comprise the structure tensor.
    float[][] h = new float[n2][n1];
    float[][][] gs = {g11,g22,g12};
    _et2.setEigenvalues(_au,_av);
    for (float[][] g:gs) {
      _lsf.applySmoothS(g,h);
      _lsf.apply(_et2,_scale,h,g);
    }

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
        float eui = e[0];
        float evi = e[1];
        if (evi<0.0f) evi = 0.0f;
        if (eui<evi) eui = evi;
        el[i2][i1] = (eui-evi)/eui;
      }
    }
  }


  /**
   * Applies this filter for the specified image and outputs. All
   * outputs are optional and are computed for only non-null arrays.
   * @param x input array for 3-D image.
   * @param ep (eu-ev)/eu, a measure of planarity.
   * @param el (ev-ew)/eu, a measure of linearity.
   */
  public void applyForPlanarLinear(
    float[][][] x,float[][][] ep, float[][][] el)
  {
    // Gradient.
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g1 = ep;
    float[][][] g2 = el;
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] xs = new float[n3][n2][n1];
    if (_scaleG>0.0f) {
      _et3.setEigenvalues(0.001f,1.0f,1.0f);
      _lsf.apply(_et3,_scaleG,x,xs);
      computeOrientGradientX(xs,g1,g2,g3);
    } else {
      computeOrientGradientX(x,g1,g2,g3);
    }

    // Gradient products.
    float[][][] g11 = g1;
    float[][][] g22 = g2;
    float[][][] g33 = g3;
    float[][][] g12 = xs;
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g11,g22,g33,g12,g13,g23};
    _et3.setEigenvalues(_au,_av,_aw);
    for (float[][][] g:gs) {
      _lsf.applySmoothS(g,h);
      _lsf.apply(_et3,_scale,h,g);
    }

    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(g11,g12,g13,g22,g23,g33,ep,el);
  }

  public void applyForPlanar( float[][][] x,float[][][] ep) {
    // Gradient.
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g1 = ep;
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] xs = new float[n3][n2][n1];
    if (_scaleG>0.0f) {
      _et3.setEigenvalues(0.001f,1.0f,1.0f);
      _lsf.apply(_et3,_scaleG,x,xs);
      computeOrientGradientX(xs,g1,g2,g3);
    } else {
      computeOrientGradientX(x,g1,g2,g3);
    }

    // Gradient products.
    float[][][] g11 = g1;
    float[][][] g22 = g2;
    float[][][] g33 = g3;
    float[][][] g12 = xs;
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g11,g22,g33,g12,g13,g23};
    _et3.setEigenvalues(_au,_av,_aw);
    for (float[][][] g:gs) {
      _lsf.applySmoothS(g,h);
      _lsf.apply(_et3,_scale,h,g);
    }
    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(g11,g12,g13,g22,g23,g33,ep);
  }


  public void applyForLinear(
    float[][][] x,float[][][] el)
  {
    // Gradient.
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g2 = el;
    float[][][] g3 = new float[n3][n2][n1];
    computeOrientGradientX(x,g2,g3);
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

    //Smoothed gradient products comprise the structure tensor.
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g22,g33,g23};
    _et3.setEigenvalues(_au,_av,_aw);
    for (float[][][] g:gs) {
      _lsf.applySmoothS(g,h);
      _lsf.apply(_et3,_scale,h,g);
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
        el[i3][i2][i1] = (eui-evi)/eui;
      }
    }}
  }



  public void updateTensors(float scale, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    applyForW(scale,x,w2,w3);
    updateTensors(_et3,w2,w3);
  }

  public void updateTensors(
    final EigenTensors3 et, final float[][][] w2, final float[][][] w3) {
    final int n3 = w2.length;
    final int n2 = w2[0].length;
    final int n1 = w2[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float[] u = et.getEigenvectorU(i1,i2,i3);
        float u1i = u[0];
        float u2i = u[1];
        float u3i = u[2];
        float w2i = w2[i3][i2][i1];
        float w3i = w3[i3][i2][i1];
        float w1i = 1f;
        if (u1i!=0f) w1i = -(w2i*u2i+w3i*u3i)/u1i;
        et.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
        et.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
      }}
    }});
  }

  public void applyForW(
    float scale, float[][][] x, float[][][] w2, float[][][] w3) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] xs = new float[n3][n2][n1];
    _et3.setEigenvalues(0.001f,1.0f,1.0f);
    _lsf.apply(_et3,3,x,xs);
    computeOrientGradient(xs,g2,g3);
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

    //Smoothed gradient products comprise the structure tensor.
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g22,g33,g23};
    _et3.setEigenvalues(0.1f,1.0f,1.0f);
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setCycles(5,0.1f);
    for (float[][][] g:gs) {
      _lsf.applySmoothS(g,h);
      h = fed.apply(scale,_et3,h);
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
        float x2i = z[0][0];
        float x3i = z[0][1];
        if (x2i<0.0f) {
          x2i = -x2i;
          x3i = -x3i;
        }
        w2[i3][i2][i1] = -x3i;
        w3[i3][i2][i1] =  x2i;
      }
    }}
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  public void computeOrientGradient(
    final float[][] fx, final float[][] g1, final float[][] g2)
  {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final float[][] p2 = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float[] u = _et2.getEigenvectorU(i1,i2);
        float u1i = u[0];
        float u2i = u[1];
        if (u1i<0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        float p2t = -u2i/u1i;
        p2t = max(p2t,-10f);
        p2t = min(p2t, 10f);
        p2[i2][i1] = p2t;
        float u1p = i1+u1i;
        float u2p = i2+u2i;
        float u1m = i1-u1i;
        float u2m = i2-u2i;
        float fup = si.interpolate(s1,s2,fx,u1p,u2p);
        float fum = si.interpolate(s1,s2,fx,u1m,u2m);
        g1[i2][i1] = fup-fum;
      }
    }});
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float[] g2i = g2[i2];
      float[] p2i = p2[i2];
      float[] fmi = fx[i2m];
      float[] fpi = fx[i2p];
      float[] xmi = new float[n2];
      float[] xpi = new float[n2];
      float[] g2m = new float[n2];
      float[] g2p = new float[n2];
      for (int i1=0; i1<n1; ++i1) {
        xmi[i1] = i1-p2i[i1];
        xpi[i1] = i1+p2i[i1];
      }
      si.interpolate(n1,1.0,0.0,fmi,n1,xmi,g2m);
      si.interpolate(n1,1.0,0.0,fpi,n1,xpi,g2p);
      for (int i1=0; i1<n1; ++i1)
        g2i[i1] = g2p[i1]-g2m[i1];
    }});
  }

  public void computeOrientGradient(
    final float[][][] fx, final float[][][] g1, 
    final float[][][] g2, final float[][][] g3) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final float[][][] p2 = new float[n3][n2][n1];
    final float[][][] p3 = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] p2i = p2[i3][i2];
        float[] p3i = p3[i3][i2];
        float[] g1i = g1[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          float[] u = _et3.getEigenvectorU(i1,i2,i3);
          float u1i = u[0];
          float u2i = u[1];
          float u3i = u[2];
          if (u1i<0f) {
            u1i = -u1i;
            u2i = -u2i;
            u3i = -u3i;
          }
          if (u1i==0f) {
            u1i = 1f;
            u2i = 0f;
            u3i = 0f;
          }
          float p2t = -u2i/u1i;
          float p3t = -u3i/u1i;
          p2t = max(p2t,-10f);
          p2t = min(p2t, 10f);
          p3t = max(p3t,-10f);
          p3t = min(p3t, 10f);
          p2i[i1] = p2t;
          p3i[i1] = p3t;
          float u1p = i1+u1i;
          float u2p = i2+u2i;
          float u3p = i3+u3i;
          float u1m = i1-u1i;
          float u2m = i2-u2i;
          float u3m = i3-u3i;
          float g1p = si.interpolate(s1,s2,s3,fx,u1p,u2p,u3p);
          float g1m = si.interpolate(s1,s2,s3,fx,u1m,u2m,u3m);
          g1i[i1] = g1p-g1m;
      }}
    }});
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] x2m = new float[n1];
      float[] x2p = new float[n1];
      float[] x3m = new float[n1];
      float[] x3p = new float[n1];
      float[] g2m = new float[n1];
      float[] g2p = new float[n1];
      float[] g3m = new float[n1];
      float[] g3p = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] g2i = g2[i3 ][i2 ];
        float[] g3i = g3[i3 ][i2 ];
        float[] p2i = p2[i3 ][i2 ];
        float[] p3i = p3[i3 ][i2 ];
        float[] f2m = fx[i3 ][i2m];
        float[] f2p = fx[i3 ][i2p];
        float[] f3m = fx[i3m][i2 ];
        float[] f3p = fx[i3p][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          x2m[i1] = i1-p2i[i1];
          x2p[i1] = i1+p2i[i1];
          x3m[i1] = i1-p3i[i1];
          x3p[i1] = i1+p3i[i1];
        }
        si.interpolate(n1,1.0,0.0,f2m,n1,x2m,g2m);
        si.interpolate(n1,1.0,0.0,f2p,n1,x2p,g2p);
        si.interpolate(n1,1.0,0.0,f3m,n1,x3m,g3m);
        si.interpolate(n1,1.0,0.0,f3p,n1,x3p,g3p);
        for (int i1=0; i1<n1; ++i1) {
          g2i[i1] = g2p[i1]-g2m[i1];
          g3i[i1] = g3p[i1]-g3m[i1];
        }
      }
    }});
  }

  public void computeOrientGradient(
    final float[][][] fx, final float[][][] g2, final float[][][] g3) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] p2 = new float[n3][n2][n1];
    final float[][][] p3 = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] p2i = p2[i3][i2];
        float[] p3i = p3[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          float[] u = _et3.getEigenvectorU(i1,i2,i3);
          float u1i = u[0];
          float u2i = u[1];
          float u3i = u[2];
          if (u1i<0f) {
            u1i = -u1i;
            u2i = -u2i;
            u3i = -u3i;
          }
          float p2t = -u2i/u1i;
          float p3t = -u3i/u1i;
          p2t = max(p2t,-10f);
          p2t = min(p2t, 10f);
          p3t = max(p3t,-10f);
          p3t = min(p3t, 10f);
          p2i[i1] = p2t;
          p3i[i1] = p3t;
        }
      }
    }});
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] x2m = new float[n1];
      float[] x2p = new float[n1];
      float[] x3m = new float[n1];
      float[] x3p = new float[n1];
      float[] g2m = new float[n1];
      float[] g2p = new float[n1];
      float[] g3m = new float[n1];
      float[] g3p = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] g2i = g2[i3 ][i2 ];
        float[] g3i = g3[i3 ][i2 ];
        float[] p2m = p2[i3 ][i2m];
        float[] p2p = p2[i3 ][i2p];
        float[] p3m = p3[i3m][i2 ];
        float[] p3p = p3[i3p][i2 ];
        float[] f2m = fx[i3 ][i2m];
        float[] f2p = fx[i3 ][i2p];
        float[] f3m = fx[i3m][i2 ];
        float[] f3p = fx[i3p][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          x2m[i1] = i1-p2m[i1];
          x2p[i1] = i1+p2p[i1];
          x3m[i1] = i1-p3m[i1];
          x3p[i1] = i1+p3p[i1];
        }
        si.interpolate(n1,1.0,0.0,f2m,n1,x2m,g2m);
        si.interpolate(n1,1.0,0.0,f2p,n1,x2p,g2p);
        si.interpolate(n1,1.0,0.0,f3m,n1,x3m,g3m);
        si.interpolate(n1,1.0,0.0,f3p,n1,x3p,g3p);
        for (int i1=0; i1<n1; ++i1) {
          g2i[i1] = g2p[i1]-g2m[i1];
          g3i[i1] = g3p[i1]-g3m[i1];
        }
      }
    }});
  }

  public void computeOrientGradientX(
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
          float[] v = _et3.getEigenvectorV(i1,i2,i3);
          float[] w = _et3.getEigenvectorW(i1,i2,i3);
          float v1i = v[0];
          float v2i = v[1];
          float v3i = v[2];
          float w1i = w[0];
          float w2i = w[1];
          float w3i = w[2];
          if (v2i<0f) {
            v1i = -v1i; v2i = -v2i; v3i = -v3i;
          }
          if (w3i<0f) {
            w1i = -w1i; w2i = -w2i; w3i = -w3i;
          }
          float v1p = i1+v1i;
          float v2p = i2+v2i;
          float v3p = i3+v3i;
          float v1m = i1-v1i;
          float v2m = i2-v2i;
          float v3m = i3-v3i;
          float w1p = i1+w1i;
          float w2p = i2+w2i;
          float w3p = i3+w3i;
          float w1m = i1-w1i;
          float w2m = i2-w2i;
          float w3m = i3-w3i;
          float g2p = si.interpolate(s1,s2,s3,fx,v1p,v2p,v3p);
          float g2m = si.interpolate(s1,s2,s3,fx,v1m,v2m,v3m);
          float g3p = si.interpolate(s1,s2,s3,fx,w1p,w2p,w3p);
          float g3m = si.interpolate(s1,s2,s3,fx,w1m,w2m,w3m);
          g2i[i1] = g2p-g2m;
          g3i[i1] = g3p-g3m;
      }}
    }});
  }



  public void computeOrientGradientX(
    final float[][][] fx, final float[][][] g1, 
    final float[][][] g2, final float[][][] g3) 
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
        float[] g1i = g1[i3][i2];
        float[] g2i = g2[i3][i2];
        float[] g3i = g3[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          float[] u = _et3.getEigenvectorU(i1,i2,i3);
          float[] v = _et3.getEigenvectorV(i1,i2,i3);
          float[] w = _et3.getEigenvectorW(i1,i2,i3);
          float u1i = u[0];
          float u2i = u[1];
          float u3i = u[2];
          float v1i = v[0];
          float v2i = v[1];
          float v3i = v[2];
          float w1i = w[0];
          float w2i = w[1];
          float w3i = w[2];
          if (u1i<0f) {
            u1i = -u1i; u2i = -u2i; u3i = -u3i;
          }
          if (v2i<0f) {
            v1i = -v1i; v2i = -v2i; v3i = -v3i;
          }
          if (w3i<0f) {
            w1i = -w1i; w2i = -w2i; w3i = -w3i;
          }
          float u1p = i1+u1i;
          float u2p = i2+u2i;
          float u3p = i3+u3i;
          float u1m = i1-u1i;
          float u2m = i2-u2i;
          float u3m = i3-u3i;
          float v1p = i1+v1i;
          float v2p = i2+v2i;
          float v3p = i3+v3i;
          float v1m = i1-v1i;
          float v2m = i2-v2i;
          float v3m = i3-v3i;
          float w1p = i1+w1i;
          float w2p = i2+w2i;
          float w3p = i3+w3i;
          float w1m = i1-w1i;
          float w2m = i2-w2i;
          float w3m = i3-w3i;
          float g1p = si.interpolate(s1,s2,s3,fx,u1p,u2p,u3p);
          float g1m = si.interpolate(s1,s2,s3,fx,u1m,u2m,u3m);
          float g2p = si.interpolate(s1,s2,s3,fx,v1p,v2p,v3p);
          float g2m = si.interpolate(s1,s2,s3,fx,v1m,v2m,v3m);
          float g3p = si.interpolate(s1,s2,s3,fx,w1p,w2p,w3p);
          float g3m = si.interpolate(s1,s2,s3,fx,w1m,w2m,w3m);
          g1i[i1] = g1p-g1m;
          g2i[i1] = g2p-g2m;
          g3i[i1] = g3p-g3m;
      }}
    }});
  }

  private void computeGradientProducts(
    final float[][][] g1, final float[][][] g2, final float[][][] g3,
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33)
  {
    final int n1 = g1[0][0].length;
    final int n2 = g1[0].length;
    final int n3 = g1.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] g1i = g1[i3][i2];
          float[] g2i = g2[i3][i2];
          float[] g3i = g3[i3][i2];
          float[] g11i = g11[i3][i2];
          float[] g12i = g12[i3][i2];
          float[] g13i = g13[i3][i2];
          float[] g22i = g22[i3][i2];
          float[] g23i = g23[i3][i2];
          float[] g33i = g33[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float g1ii = g1i[i1];
            float g2ii = g2i[i1];
            float g3ii = g3i[i1];
            g11i[i1] = g1ii*g1ii;
            g22i[i1] = g2ii*g2ii;
            g33i[i1] = g3ii*g3ii;
            g12i[i1] = g1ii*g2ii;
            g13i[i1] = g1ii*g3ii;
            g23i[i1] = g2ii*g3ii;
          }
        }
      }
    });
  }
  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] ep)
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
            float eai = (float)e[0];
            float ebi = (float)e[1];
            float eci = (float)e[2];
            if (eci<0.0f)eci = 0.0f;
            if (ebi<eci) ebi = eci;
            if (eai<ebi) eai = ebi;
            float esi = (eai>0.0f)?1.0f/eai:1.0f;
            ep[i3][i2][i1] = (eai-ebi)*esi;
          }
        }
      }
    });
  }


  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] ep, final float[][][] el)
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
            float eai = (float)e[0];
            float ebi = (float)e[1];
            float eci = (float)e[2];
            if (eci<0.0f)eci = 0.0f;
            if (ebi<eci) ebi = eci;
            if (eai<ebi) eai = ebi;
            float esi = (eai>0.0f)?1.0f/eai:1.0f;
            ep[i3][i2][i1] = (eai-ebi)*esi;
            el[i3][i2][i1] = (ebi-eci)*esi;
          }
        }
      }
    });
  }


  private float _au=1.00f;
  private float _av=0.05f;
  private float _aw=1.00f;
  private float _scale = 20f;
  private EigenTensors2 _et2=null;
  private EigenTensors3 _et3=null;
  private float _scaleG = 0.0f;
  private LocalSmoothingFilter _lsf = new LocalSmoothingFilter(0.0001f,1000);
}
