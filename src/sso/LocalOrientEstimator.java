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
 * @version 2016.08.05
 */
public class LocalOrientEstimator {

  /**
   * Constructs a filter with anisotropic structure-oriented smoothing.
   * @param scale factor for the smoothing
   */
  public LocalOrientEstimator(EigenTensors2 et, double scale) {
    _et2 = et;
    _scale = (float)scale;
    setEigenvalues(1.0f,0.05f);
  }

  /**
   * Constructs a filter with anisotropic structure-oriented smoothing.
   * @param scale factor for the smoothing
   */
  public LocalOrientEstimator(EigenTensors3 et, double scale) { 
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
   * Applies this filter to estimate orientation angles.
   * @param x input array for 2-D image.
   * @param theta orientation angle; -pi &lt;= theta &lt;= pi
   */
  public void applyForTheta(float[][] x, float[][] theta) {
    apply(x,
      theta,
      null,null,
      null,null,
      null,null,
      null);
  }

  /**
   * Applies this filter to estimate normal vectors (1st eigenvectors).
   * @param x input array for 2-D image.
   * @param u1 1st component of normal vector.
   * @param u2 2nd component of normal vector.
   */
  public void applyForNormal(float[][] x, float[][] u1, float[][] u2) {
    apply(x,
      null,
      u1,u2,
      null,null,
      null,null,
      null);
  }

  /**
   * Applies this filter to estimate normal vectors and linearities.
   * @param x input array for 2-D image.
   * @param u1 1st component of normal vector.
   * @param u2 2nd component of normal vector.
   * @param el linearity in range [0,1].
   */
  public void applyForNormalLinear(float[][] x, 
    float[][] u1, float[][] u2, float[][] el) 
  {
    apply(x,
      null,
      u1,u2,
      null,null,
      null,null,
      el);
  }

  /**
   * Applies this filter to estimate 2-D structure tensors.
   * @param x input array for 2-D image.
   * @return structure tensors.
   */
  public EigenTensors2 applyForTensors(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] eu = new float[n2][n1];
    float[][] ev = new float[n2][n1];
    apply(x,
      null,
      u1,u2,
      null,null,
      eu,ev,
      null);
    return new EigenTensors2(u1,u2,eu,ev);
  }

  /**
   * Applies this filter for the specified image and outputs. All
   * outputs are optional and are computed for only non-null arrays.
   * @param x input array for 2-D image
   * @param theta orientation angle = asin(u2); -pi &lt;= theta &lt;= pi
   * @param u1 1st component of 1st eigenvector.
   * @param u2 2nd component of 1st eigenvector.
   * @param v1 1st component of 2nd eigenvector.
   * @param v2 2nd component of 2nd eigenvector.
   * @param eu largest eigenvalue corresponding to the eigenvector u.
   * @param ev smallest eigenvalue corresponding to the eigenvector v.
   * @param el (eu-ev)/eu, a measure of linearity.
   */
  public void apply(float[][] x,
    float[][] theta,
    float[][] u1, float[][] u2, 
    float[][] v1, float[][] v2,
    float[][] eu, float[][] ev, 
    float[][] el)
  {
    // Where possible, use output arrays for workspace.
    float[][][] t = new float[8][][];
    int nt = 0;
    if (theta!=null) t[nt++] = theta;
    if (u1!=null) t[nt++] = u1;
    if (u2!=null) t[nt++] = u2;
    if (v1!=null) t[nt++] = v1;
    if (v2!=null) t[nt++] = v2;
    if (eu!=null) t[nt++] = eu;
    if (ev!=null) t[nt++] = ev;
    if (el!=null) t[nt++] = el;

    // Gradient.
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] g1 = (nt>0)?t[0]:new float[n2][n1];
    float[][] g2 = (nt>1)?t[1]:new float[n2][n1];
    float[][] h = (nt>2)?t[2]:new float[n2][n1];
    //_et2.setEigenvalues(0.001f,1.0f);
    //_lsf.apply(_et2,3,x,h);
    computeOrientGradient(x,g1,g2);

    // Gradient products.
    float[][] g11 = g1;
    float[][] g22 = g2;
    float[][] g12 = (nt>3)?t[3]:new float[n2][n1];
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
        float[] u = _et2.getEigenvectorU(i1,i2);
        float u1i = u[0];
        float u2i = u[1];
        if (u1i<0f) {
          u1i = -u1i; u2i = -u2i;
        }
        float v1i = -u2i/u1i;
        float v2i = 1;
        v2i  = 1f/(v1i*v1i+1f);
        v1i *= v2i;

        a[0][0] = g11[i2][i1];
        a[0][1] = g12[i2][i1];
        a[1][0] = g12[i2][i1];
        a[1][1] = g22[i2][i1];
        Eigen.solveSymmetric22(a,z,e);
        float eui = e[0];
        float evi = e[1];
        if (evi<0.0f) evi = 0.0f;
        if (eui<evi) eui = evi;
        float x1i = z[0][0];
        float x2i = z[0][1];
        if (x1i<0.0f) {
          x1i = -x1i;
          x2i = -x2i;
        }
        float a1i = u1i*x1i+v1i*x2i;
        float a2i = u2i*x1i+v2i*x2i;
        float asi = 1f/sqrt(a1i*a1i+a2i*a2i);
        a1i *= asi; a2i *= asi;
        float b1i = -a2i;
        float b2i =  a1i;

        if (theta!=null) theta[i2][i1] = asin(a2i);
        if (u1!=null) u1[i2][i1] = a1i;
        if (u2!=null) u2[i2][i1] = a2i;
        if (v1!=null) v1[i2][i1] = b1i;
        if (v2!=null) v2[i2][i1] = b2i;
        if (eu!=null) eu[i2][i1] = eui;
        if (ev!=null) ev[i2][i1] = evi;
        if (el!=null) el[i2][i1] = (eui-evi)/eui;
      }
    }
  }

  /**
   * Applies this filter to estimate orientation angles.
   * @param x input array for 3-D image.
   * @param theta orientation dip angle; 0 &lt;= theta &lt;= pi/2.
   * @param phi orientation azimuthal angle; -pi &lt;= phi &lt;= pi.
   */
  public void applyForThetaPhi(float[][][] x, 
    float[][][] theta, float[][][] phi) 
  {
    apply(x,
      theta,phi,
      null,null,null,
      null,null,null,
      null,null,null,
      null,null,null,
      null,null);
  }



  /**
   * Applies this filter to estimate normal vectors (1st eigenvectors).
   * @param x input array for 3-D image.
   * @param u1 1st component of normal vector.
   * @param u2 2nd component of normal vector.
   * @param u3 3rd component of normal vector.
   */
  public void applyForNormal(float[][][] x, 
    float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    apply(x,
      null,null,
      u1,u2,u3,
      null,null,null,
      null,null,null,
      null,null,null,
      null,null);
  }

  /**
   * Applies this filter to estimate normal vectors and planarities.
   * Normal vectors are 1st eigenvectors corresponding to largest eigenvalues.
   * @param x input array for 3-D image.
   * @param u1 1st component of normal vector.
   * @param u2 2nd component of normal vector.
   * @param u3 3rd component of normal vector.
   * @param ep planarity in range [0,1].
   */
  public void applyForNormalPlanar(float[][][] x, 
    float[][][] u1, float[][][] u2, float[][][] u3, float[][][] ep) 
  {
    apply(x,
      null,null,
      u1,u2,u3,
      null,null,null,
      null,null,null,
      null,null,null,
      ep,null);
  }

  /**
   * Applies this filter to estimate inline vectors (3rd eigenvectors).
   * @param x input array for 3-D image.
   * @param w1 1st component of inline vector.
   * @param w2 2nd component of inline vector.
   * @param w3 3rd component of inline vector.
   */
  public void applyForInline(float[][][] x, 
    float[][][] w1, float[][][] w2, float[][][] w3)
  {
    apply(x,
      null,null,
      null,null,null,
      null,null,null,
      w1,w2,w3,
      null,null,null,
      null,null);
  }

  /**
   * Applies this filter to estimate inline vectors and linearities.
   * Inline vectors are 3rd eigenvectors corresponding to smallest eigenvalues.
   * @param x input array for 3-D image.
   * @param w1 1st component of inline vector.
   * @param w2 2nd component of inline vector.
   * @param w3 3rd component of inline vector.
   * @param el linearity in range [0,1].
   */
  public void applyForInlineLinear(float[][][] x, 
    float[][][] w1, float[][][] w2, float[][][] w3,
    float[][][] el) 
  {
    apply(x,
      null,null,
      null,null,null,
      null,null,null,
      w1,w2,w3,
      null,null,null,
      null,el);
  }

  /**
   * Applies this filter to estimate seismic slopes.
   * @param x input array for 3-D image.
   * @param p2 inline slopes.
   * @param p3 crossline slopes.
   * @param ep planarity in range [0,1].
   */
  public void applyForSlopePlanar(float pmax, float[][][] x, 
    float[][][] p2, float[][][] p3, float[][][] ep) 
  {
    float pmin = -pmax;
    int n3 = p2.length;
    int n2 = p2[0].length;
    int n1 = p3[0][0].length;
    float[][][] u1 = p2;
    float[][][] u2 = p3;
    float[][][] u3 = new float[n3][n2][n1];
    apply(x,
      null,null,
      u1,u2,u3,
      null,null,null,
      null,null,null,
      null,null,null,
      ep,null);
    // Compute slopes from normal vectors.
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      if (-u2i<pmin*u1i) u2i = -pmin*u1i;
      if (-u2i>pmax*u1i) u2i = -pmax*u1i;
      if (-u3i<pmin*u1i) u3i = -pmin*u1i;
      if (-u3i>pmax*u1i) u3i = -pmax*u1i;
      if (u1i==0.0f) {
        p2[i3][i2][i1] = (u2i<0.0f)?pmax:pmin;
        p3[i3][i2][i1] = (u3i<0.0f)?pmax:pmin;
      } else {
        p2[i3][i2][i1] = -u2i/u1i;
        p3[i3][i2][i1] = -u3i/u1i;
      }
    }}}

  }

    /**
   * Applies this filter to estimate 3-D structure tensors.
   * @param x input array for 3-D image.
   * @param compressed true, for compressed tensors; false, otherwise.
   * @return structure tensors.
   */
  public EigenTensors3 applyForTensors(float[][][] x) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] eu = new float[n3][n2][n1];
    float[][][] ev = new float[n3][n2][n1];
    float[][][] ew = new float[n3][n2][n1];
    apply(x,
      null,null,
      null,u2,u3,
      null,null,null,
      w1,w2,null,
      eu,ev,ew,
      null,null);

    // Compute u1 such that u3 > 0.
    float[][][] u1 = u3;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float u1s = 1.0f-u2i*u2i-u3i*u3i;
          float u1i = (u1s>0.0f)?sqrt(u1s):0.0f;
          if (u3i<0.0f) {
            u1i = -u1i;
            u2i = -u2i;
          }
          u1[i3][i2][i1] = u1i;
          u2[i3][i2][i1] = u2i;
        }
      }
    }
    return new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
  }

  public EigenTensors3 applyForTensorsX(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] eu = new float[n3][n2][n1];
    float[][][] ev = new float[n3][n2][n1];
    float[][][] ew = new float[n3][n2][n1];
    applyX(x,
      null,null,
      null,u2,u3,
      null,null,null,
      w1,w2,null,
      eu,ev,ew,
      null,null);

    // Compute u1 such that u3 > 0.
    float[][][] u1 = u3;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float u1s = 1.0f-u2i*u2i-u3i*u3i;
          float u1i = (u1s>0.0f)?sqrt(u1s):0.0f;
          if (u3i<0.0f) {
            u1i = -u1i;
            u2i = -u2i;
          }
          u1[i3][i2][i1] = u1i;
          u2[i3][i2][i1] = u2i;
        }
      }
    }
    return new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
  }

  public float[][][][] slopesFromNormals(float pmax, 
    float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float p2min = -pmax;
    float p3min = -pmax;
    float p2max =  pmax;
    float p3max =  pmax;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      if (u1i<0f) {
        u1i = -u1i;
        u2i = -u2i;
        u3i = -u3i;
      }
      if (-u2i<p2min*u1i) u2i = -p2min*u1i;
      if (-u2i>p2max*u1i) u2i = -p2max*u1i;
      if (-u3i<p3min*u1i) u3i = -p3min*u1i;
      if (-u3i>p3max*u1i) u3i = -p3max*u1i;
      if (u1i==0.0f) {
        p2[i3][i2][i1] = (u2i<0.0f)?p2max:p2min;
        p3[i3][i2][i1] = (u3i<0.0f)?p3max:p3min;
      } else {
        p2[i3][i2][i1] = -u2i/u1i;
        p3[i3][i2][i1] = -u3i/u1i;
      }
    }}}
    return new float[][][][]{p2,p3};
  }

  public float[][][][] slopesFromTensors(float pmax, EigenTensors3 et) {
    int n3 = et.getN3();
    int n2 = et.getN2();
    int n1 = et.getN1();
    float p2min = -pmax;
    float p3min = -pmax;
    float p2max =  pmax;
    float p3max =  pmax;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float[] u = et.getEigenvectorU(i1,i2,i3);
      float u1i = u[0];
      float u2i = u[1];
      float u3i = u[2];
      if (u1i<0f) {
        u1i = -u1i;
        u2i = -u2i;
        u3i = -u3i;
      }
      if (-u2i<p2min*u1i) u2i = -p2min*u1i;
      if (-u2i>p2max*u1i) u2i = -p2max*u1i;
      if (-u3i<p3min*u1i) u3i = -p3min*u1i;
      if (-u3i>p3max*u1i) u3i = -p3max*u1i;
      if (u1i==0.0f) {
        p2[i3][i2][i1] = (u2i<0.0f)?p2max:p2min;
        p3[i3][i2][i1] = (u3i<0.0f)?p3max:p3min;
      } else {
        p2[i3][i2][i1] = -u2i/u1i;
        p3[i3][i2][i1] = -u3i/u1i;
      }
    }}}
    return new float[][][][]{p2,p3};
  }


  /**
   * Applies this filter for the specified image and outputs. All
   * outputs are optional and are computed for only non-null arrays.
   * @param x input array for 3-D image.
   * @param theta orientation dip angle; 0 &lt;= theta &lt;= pi/2.
   * @param phi orientation azimuthal angle; -pi &lt;= phi &lt;= pi.
   * @param u1 1st component of 1st eigenvector.
   * @param u2 2nd component of 1st eigenvector.
   * @param u3 3rd component of 1st eigenvector.
   * @param v1 1st component of 2nd eigenvector.
   * @param v2 2nd component of 2nd eigenvector.
   * @param v3 3rd component of 2nd eigenvector.
   * @param w1 1st component of 3rd eigenvector.
   * @param w2 2nd component of 3rd eigenvector.
   * @param w3 3rd component of 3rd eigenvector.
   * @param eu largest eigenvalue corresponding to the eigenvector u.
   * @param ev middle eigenvalue corresponding to the eigenvector v.
   * @param ew smallest eigenvalue corresponding to the eigenvector w.
   * @param ep (eu-ev)/eu, a measure of planarity.
   * @param el (ev-ew)/eu, a measure of linearity.
   */
  public void apply(float[][][] x,
    float[][][] theta, float[][][] phi,
    float[][][] a1, float[][][] a2, float[][][] a3, 
    float[][][] b1, float[][][] b2, float[][][] b3, 
    float[][][] c1, float[][][] c2, float[][][] c3, 
    float[][][] ea, float[][][] eb, float[][][] ec, 
    float[][][] ep, float[][][] el)
  {
    // Where possible, use output arrays for workspace.
    int nt = 0;
    float[][][][] t = new float[16][][][];
    if (theta!=null) t[nt++] = theta;
    if (phi!=null) t[nt++] = phi;
    if (a1!=null) t[nt++] = a1;
    if (a2!=null) t[nt++] = a2;
    if (a3!=null) t[nt++] = a3;
    if (b1!=null) t[nt++] = b1;
    if (b2!=null) t[nt++] = b2;
    if (b3!=null) t[nt++] = b3;
    if (c1!=null) t[nt++] = c1;
    if (c2!=null) t[nt++] = c2;
    if (c3!=null) t[nt++] = c3;
    if (ea!=null) t[nt++] = ea;
    if (eb!=null) t[nt++] = eb;
    if (ec!=null) t[nt++] = ec;
    if (ep!=null) t[nt++] = ep;
    if (el!=null) t[nt++] = el;

    // Gradient.
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g1 = (nt>0)?t[0]:new float[n3][n2][n1];
    float[][][] g2 = (nt>1)?t[1]:new float[n3][n2][n1];
    float[][][] g3 = (nt>2)?t[2]:new float[n3][n2][n1];
    float[][][] xs = (nt>3)?t[3]:new float[n3][n2][n1];
    if (_scaleG>0.0f) {
      _et3.setEigenvalues(0.001f,1.0f,1.0f);
      _lsf.apply(_et3,_scaleG,x,xs);
      computeOrientGradient(xs,g1,g2,g3);
    } else {
      computeOrientGradient(x,g1,g2,g3);
    }

    // Gradient products.
    float[][][] g11 = g1;
    float[][][] g22 = g2;
    float[][][] g33 = g3;
    float[][][] g12 = xs;
    float[][][] g13 = (nt>4)?t[4]:new float[n3][n2][n1];
    float[][][] g23 = (nt>5)?t[5]:new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    
    // Smoothed gradient products comprise the structure tensor.
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(_scale);
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setCycles(3,0.1f);
    float[][][] h = (nt>6)?t[6]:new float[n3][n2][n1];
    float[][][][] gs = {g11,g22,g33,g12,g13,g23};
    _et3.setEigenvalues(_au,_av,_aw);
    for (float[][][] g:gs) {
      rgf1.apply0XX(g,h);
      h = fed.apply(_scale,_et3,h);
      copy(h,g);
    }

    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(g11,g12,g13,g22,g23,g33,
      theta,phi,a1,a2,a3,b1,b2,b3,c1,c2,c3,ea,eb,ec,ep,el);
  }

  public void applyX(float[][][] x,
    float[][][] theta, float[][][] phi,
    float[][][] a1, float[][][] a2, float[][][] a3, 
    float[][][] b1, float[][][] b2, float[][][] b3, 
    float[][][] c1, float[][][] c2, float[][][] c3, 
    float[][][] ea, float[][][] eb, float[][][] ec, 
    float[][][] ep, float[][][] el)
  {
    // Where possible, use output arrays for workspace.
    int nt = 0;
    float[][][][] t = new float[16][][][];
    if (theta!=null) t[nt++] = theta;
    if (phi!=null) t[nt++] = phi;
    if (a1!=null) t[nt++] = a1;
    if (a2!=null) t[nt++] = a2;
    if (a3!=null) t[nt++] = a3;
    if (b1!=null) t[nt++] = b1;
    if (b2!=null) t[nt++] = b2;
    if (b3!=null) t[nt++] = b3;
    if (c1!=null) t[nt++] = c1;
    if (c2!=null) t[nt++] = c2;
    if (c3!=null) t[nt++] = c3;
    if (ea!=null) t[nt++] = ea;
    if (eb!=null) t[nt++] = eb;
    if (ec!=null) t[nt++] = ec;
    if (ep!=null) t[nt++] = ep;
    if (el!=null) t[nt++] = el;

    // Gradient.
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g1 = (nt>0)?t[0]:new float[n3][n2][n1];
    float[][][] g2 = (nt>1)?t[1]:new float[n3][n2][n1];
    float[][][] g3 = (nt>2)?t[2]:new float[n3][n2][n1];
    float[][][] xs = (nt>3)?t[3]:new float[n3][n2][n1];
    if (_scaleG>0.0f) {
      _et3.setEigenvalues(0.001f,1.0f,1.0f);
      _lsf.apply(_et3,_scaleG,x,xs);
      computeOrientGradient(xs,g1,g2,g3);
    } else {
      computeOrientGradient(x,g1,g2,g3);
    }

    // Gradient products.
    float[][][] g11 = g1;
    float[][][] g22 = g2;
    float[][][] g33 = g3;
    float[][][] g12 = xs;
    float[][][] g13 = (nt>4)?t[4]:new float[n3][n2][n1];
    float[][][] g23 = (nt>5)?t[5]:new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    
    // Smoothed gradient products comprise the structure tensor.
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(_scale);
    //FastExplicitDiffusion fed = new FastExplicitDiffusion();
    //fed.setCycles(3,0.1f);
    float[][][] h = (nt>6)?t[6]:new float[n3][n2][n1];
    float[][][][] gs = {g11,g22,g33,g12,g13,g23};
    _et3.setEigenvalues(_au,_av,_aw);
    for (float[][][] g:gs) {
      _lsf.applySmoothS(g,h);
      _lsf.apply(_et3,_scale,h,g);
      //rgf1.apply0XX(g,h);
      //h = fed.apply(_scale,_et3,h);
      //copy(h,g);
    }

    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(g11,g12,g13,g22,g23,g33,
      theta,phi,a1,a2,a3,b1,b2,b3,c1,c2,c3,ea,eb,ec,ep,el);
  }


  public void applyForStratigraphy(
    float[][][] x, float[][][] w2, float[][][] w3, float[][][] el) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] xs = new float[n3][n2][n1];
    FastExplicitDiffusion fed = new FastExplicitDiffusion();
    fed.setCycles(5,0.1f);
    if (_scaleG>0.0f) {
      _et3.setEigenvalues(0.001f,1.0f,1.0f);
      _lsf.apply(_et3,_scaleG,x,xs);
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
    _et3.setEigenvalues(_au,_av,_aw);
    for (float[][][] g:gs) {
      _lsf.applySmoothS(g,h);
      //_lsf.apply(_et3,_scale,h,g);
      h = fed.apply(_scale,_et3,h);
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
        float x2i = z[0][0];
        float x3i = z[0][1];
        if (x2i<0.0f) {
          x2i = -x2i;
          x3i = -x3i;
        }
        w2[i3][i2][i1] = -x3i;
        w3[i3][i2][i1] =  x2i;
        el[i3][i2][i1] = (eui-evi)/eui;
      }
    }}
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
    final float[][][] theta, final float[][][] phi,
    final float[][][] a1, final float[][][] a2, final float[][][] a3, 
    final float[][][] b1, final float[][][] b2, final float[][][] b3, 
    final float[][][] c1, final float[][][] c2, final float[][][] c3, 
    final float[][][] ea, final float[][][] eb, final float[][][] ec, 
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
            float x1i = (float)z[0][0];
            float x2i = (float)z[0][1];
            float x3i = (float)z[0][2];
            float y1i = (float)z[1][0];
            float y2i = (float)z[1][1];
            float y3i = (float)z[1][2];
            float z1i = (float)z[2][0];
            float z2i = (float)z[2][1];
            float z3i = (float)z[2][2];
            if (x1i<0.0f) {
              x1i = -x1i; x2i = -x2i; x3i = -x3i;
            }
            if (y2i<0.0f) {
              y1i = -y1i; y2i = -y2i; y3i = -y3i;
            }
            if (z3i<0.0f) {
              z1i = -z1i; z2i = -z2i; z3i = -z3i;
            }
            float[] u = _et3.getEigenvectorU(i1,i2,i3);
            float u1i = u[0];
            float u2i = u[1];
            float u3i = u[2];
            if (u1i<0f) {
              u1i = -u1i; u2i = -u2i; u3i = -u3i;
            }
            if (u1i==0f) {
              u1i = 1f; u2i = 0f; u3i = 0f;
            }
            float v1i = -u2i/u1i;
            float v2i = 1;
            float w1i = -u3i/u1i;
            float w3i = 1;
            v2i  = 1f/sqrt(v1i*v1i+1);
            v1i *= v2i;
            w3i  = 1f/sqrt(w1i*w1i+1);
            w1i *= w3i;
            if (a1!=null||a2!=null||a3!=null) {
              float a1i = u1i*x1i+v1i*x2i+w1i*x3i;
              float a2i = u2i*x1i+v2i*x2i;
              float a3i = u3i*x1i+w3i*x3i;
              float asi = 1f/sqrt(a1i*a1i+a2i*a2i+a3i*a3i);
              if (a1!=null) a1[i3][i2][i1] = a1i*asi;
              if (a2!=null) a2[i3][i2][i1] = a2i*asi;
              if (a3!=null) a3[i3][i2][i1] = a3i*asi;
            }

            if (b1!=null||b2!=null||b3!=null) {
              float b1i = u1i*y1i+v1i*y2i+w1i*y3i;
              float b2i = u2i*y1i+v2i*y2i;
              float b3i = u3i*y1i+w3i*y3i;
              float bsi = 1f/sqrt(b1i*b1i+b2i*b2i+b3i*b3i);
              if (b1!=null) b1[i3][i2][i1] = b1i*bsi;
              if (b2!=null) b2[i3][i2][i1] = b2i*bsi;
              if (b3!=null) b3[i3][i2][i1] = b3i*bsi;
            }

            if (c1!=null||c2!=null||c3!=null) {
              float c1i = u1i*z1i+v1i*z2i+w1i*z3i;
              float c2i = u2i*z1i+v2i*z2i;
              float c3i = u3i*z1i+w3i*z3i;
              float csi = 1f/sqrt(c1i*c1i+c2i*c2i+c3i*c3i);
              if (c1!=null) c1[i3][i2][i1] = c1i*csi;
              if (c2!=null) c2[i3][i2][i1] = c2i*csi;
              if (c3!=null) c3[i3][i2][i1] = c3i*csi;
            }
            float eai = (float)e[0];
            float ebi = (float)e[1];
            float eci = (float)e[2];
            if (eci<0.0f)eci = 0.0f;
            if (ebi<eci) ebi = eci;
            if (eai<ebi) eai = ebi;
            if (theta!=null) theta[i3][i2][i1] = acos(u1i);
            if (phi!=null) phi[i3][i2][i1] = atan2(u3i,u2i);
            if (ea!=null) ea[i3][i2][i1] = eai;
            if (eb!=null) eb[i3][i2][i1] = ebi;
            if (ec!=null) ec[i3][i2][i1] = eci;
            if (ep!=null || el!=null) {
              float esi = (eai>0.0f)?1.0f/eai:1.0f;
              if (ep!=null) ep[i3][i2][i1] = (eai-ebi)*esi;
              if (el!=null) el[i3][i2][i1] = (ebi-eci)*esi;
            }
          }
        }
      }
    });
  }


  private void solveEigenproblemsX(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] theta, final float[][][] phi,
    final float[][][] a1, final float[][][] a2, final float[][][] a3, 
    final float[][][] b1, final float[][][] b2, final float[][][] b3, 
    final float[][][] c1, final float[][][] c2, final float[][][] c3, 
    final float[][][] ea, final float[][][] eb, final float[][][] ec, 
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
            float x1i = (float)z[0][0];
            float x2i = (float)z[0][1];
            float x3i = (float)z[0][2];
            float y1i = (float)z[1][0];
            float y2i = (float)z[1][1];
            float y3i = (float)z[1][2];
            float z1i = (float)z[2][0];
            float z2i = (float)z[2][1];
            float z3i = (float)z[2][2];
            if (x1i<0.0f) {
              x1i = -x1i; x2i = -x2i; x3i = -x3i;
            }
            if (y2i<0.0f) {
              y1i = -y1i; y2i = -y2i; y3i = -y3i;
            }
            if (z3i<0.0f) {
              z1i = -z1i; z2i = -z2i; z3i = -z3i;
            }
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

            if (a1!=null) {
              float a1i = u1i*x1i+v1i*x2i+w1i*x3i;
              float a2i = u2i*x1i+v2i*x2i+w2i*x3i;
              float a3i = u3i*x1i+v3i*x2i+w3i*x3i;
              float asi = 1f/sqrt(a1i*a1i+a2i*a2i+a3i*a3i);
              if (a1i<0f) {
                a1[i3][i2][i1] = -a1i*asi;
                a2[i3][i2][i1] = -a2i*asi;
                a3[i3][i2][i1] = -a3i*asi;
              } else {
                a1[i3][i2][i1] = a1i*asi;
                a2[i3][i2][i1] = a2i*asi;
                a3[i3][i2][i1] = a3i*asi;
              }
            }

            if (b1!=null) {
              float b1i = u1i*y1i+v1i*y2i+w1i*y3i;
              float b2i = u2i*y1i+v2i*y2i+w2i*y3i;
              float b3i = u3i*y1i+v2i*y2i+w3i*y3i;
              float bsi = 1f/sqrt(b1i*b1i+b2i*b2i+b3i*b3i);
              b1[i3][i2][i1] = b1i*bsi;
              b2[i3][i2][i1] = b2i*bsi;
              b3[i3][i2][i1] = b3i*bsi;
            }

            if (c1!=null) {
              float c1i = u1i*z1i+v1i*z2i+w1i*z3i;
              float c2i = u2i*z1i+v2i*z2i+w2i*z3i;
              float c3i = u3i*z1i+v3i*z2i+w3i*z3i;
              float csi = 1f/sqrt(c1i*c1i+c2i*c2i+c3i*c3i);
              c1[i3][i2][i1] = c1i*csi;
              c2[i3][i2][i1] = c2i*csi;
              c3[i3][i2][i1] = c3i*csi;
            }

            float eai = (float)e[0];
            float ebi = (float)e[1];
            float eci = (float)e[2];
            if (eci<0.0f)eci = 0.0f;
            if (ebi<eci) ebi = eci;
            if (eai<ebi) eai = ebi;
            if (theta!=null) theta[i3][i2][i1] = acos(u1i);
            if (phi!=null) phi[i3][i2][i1] = atan2(u3i,u2i);
            if (ea!=null) ea[i3][i2][i1] = eai;
            if (eb!=null) eb[i3][i2][i1] = ebi;
            if (ec!=null) ec[i3][i2][i1] = eci;
            if (ep!=null || el!=null) {
              float esi = (eai>0.0f)?1.0f/eai:1.0f;
              if (ep!=null) ep[i3][i2][i1] = (eai-ebi)*esi;
              if (el!=null) el[i3][i2][i1] = (ebi-eci)*esi;
            }
          }
        }
      }
    });
  }


  private void resetTensors() {
    int n1 = _et3.getN1();
    int n2 = _et3.getN2();
    int n3 = _et3.getN3();
    for (int i3=0; i3<n3;++i3) {
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      float[] u = _et3.getEigenvectorU(i1,i2,i3);
      float u1i = u[0];
      float u2i = u[1];
      float u3i = u[2];
      if (u1i<0f) {
        u1i = -u1i; u2i = -u2i; u3i = -u3i;
      }
      float w1i = -u3i/u1i;
      float w2i = 0f;
      float w3i = 1;
      w3i  = 1f/(w1i*w1i+1);
      w1i *= w3i;
      _et3.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
      _et3.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
    }}}
    System.out.println("test!!!!");

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
