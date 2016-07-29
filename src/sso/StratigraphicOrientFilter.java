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

import vec.*;
import util.*;

/**
 * Shaping regulerization for estimating orientations of features in images.
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
 * @version 2016.07.01
 */
public class StratigraphicOrientFilter {

  /**
   * Constructs a filter with an isotropic Gaussian window.
   * @param sigma half-width of window; same for all dimensions.
   */
  public StratigraphicOrientFilter(double sigma) {
    this(sigma,sigma,sigma);
  }

  
  /**
   * Constructs a filter with a possibly anisotropic Gaussian window.
   * @param sigma1 half-width of window in 1st dimension.
   * @param sigma2 half-width of window in 2nd and higher dimensions.
   */
  public StratigraphicOrientFilter(double sigma1, double sigma2) {
    this(sigma1,sigma2,sigma2);
  }

  /**
   * Constructs a filter with a possibly anisotropic Gaussian window.
   * @param sigma1 half-width of window in 1st dimension.
   * @param sigma2 half-width of window in 2nd dimension.
   * @param sigma3 half-width of window in 3rd and higher dimensions.
   */
  public StratigraphicOrientFilter(double sigma1, double sigma2, double sigma3) {
    _rgfSmoother1 = (sigma1>=1.0)?new RecursiveGaussianFilter(sigma1):null;
    if (sigma2==sigma1) {
      _rgfSmoother2 = _rgfSmoother1;
    } else {
      _rgfSmoother2 = (sigma2>=1.0)?new RecursiveGaussianFilter(sigma2):null;
    }
    if (sigma3==sigma2) {
      _rgfSmoother3 = _rgfSmoother2;
    } else {
      _rgfSmoother3 = (sigma3>=1.0)?new RecursiveGaussianFilter(sigma3):null;
    }
    setGradientSmoothing(1.0);
  }

  /**
   * Sets half-width of Gaussian derivative filter used to compute gradients.
   * Typically, this half-width should not exceed one-fourth that of the
   * the corresponding Gaussian window used to compute local averages of 
   * gradient products.
   * The default half-width for Gaussian derivatives is 1.0.
   * @param sigma half-width of derivatives; same for all dimensions.
   */
  public void setGradientSmoothing(double sigma) {
    setGradientSmoothing(sigma,sigma,sigma);
  }

  /**
   * Sets half-widths of Gaussian derivative filters used to compute gradients.
   * Typically, these half-widths should not exceed one-fourth those of the
   * the corresponding Gaussian windows used to compute local averages of 
   * gradient-squared tensors.
   * The default half-widths for Gaussian derivatives is 1.0.
   * @param sigma1 half-width of derivative in 1st dimension.
   * @param sigma2 half-width of derivatives in 2nd and higher dimensions.
   */
  public void setGradientSmoothing(double sigma1, double sigma2) {
    setGradientSmoothing(sigma1,sigma2,sigma2);
  }

  /**
   * Sets half-widths of Gaussian derivative filters used to compute gradients.
   * Typically, these half-widths should not exceed one-fourth those of the
   * the corresponding Gaussian windows used to compute local averages of 
   * gradient-squared tensors.
   * The default half-widths for Gaussian derivatives is 1.0.
   * @param sigma1 half-width of derivative in 1st dimension.
   * @param sigma2 half-width of derivative in 2nd dimension.
   * @param sigma3 half-width of derivatives in 3rd and higher dimensions.
   */
  public void setGradientSmoothing(
    double sigma1, double sigma2, double sigma3) 
  {
    _rgfGradient1 = new RecursiveGaussianFilter(sigma1);
    if (sigma2==sigma1) {
      _rgfGradient2 = _rgfGradient1;
    } else {
      _rgfGradient2 = new RecursiveGaussianFilter(sigma2);
    }
    if (sigma3==sigma2) {
      _rgfGradient3 = _rgfGradient2;
    } else {
      _rgfGradient3 = new RecursiveGaussianFilter(sigma3);
    }
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
    _rgfGradient1.apply10(x,g1);
    _rgfGradient2.apply01(x,g2);

    // Gradient products.
    float[][] g11 = g1;
    float[][] g22 = g2;
    float[][] g12 = (nt>2)?t[2]:new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        g11[i2][i1] = g1i*g1i;
        g22[i2][i1] = g2i*g2i;
        g12[i2][i1] = g1i*g2i;
      }
    }
    
    /*
    // Smoothed gradient products comprise the structure tensor.
    if (_rgfSmoother1!=null || _rgfSmoother2!=null) {
      float[][] h = (nt>3)?t[3]:new float[n2][n1];
      float[][][] gs = {g11,g22,g12};
      for (float[][] g:gs) {
        if (_rgfSmoother1!=null) {
          _rgfSmoother1.apply0X(g,h);
        } else {
          copy(g,h);
        }
        if (_rgfSmoother2!=null) {
          _rgfSmoother2.applyX0(h,g);
        } else {
          copy(h,g);
        }
      }
    }
    */
    //Smoothed gradient products comprise the structure tensor.
    float[][] h = (nt>3)?t[3]:new float[n2][n1];
    float[][][] gs = {g11,g22,g12};
    for (float[][] g:gs) {
      h = smooth(8,8,x,g);
      copy(h,g);
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
        float u1i = z[0][0];
        float u2i = z[0][1];
        if (u1i<0.0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        float v1i = -u2i;
        float v2i = u1i;
        float eui = e[0];
        float evi = e[1];
        if (evi<0.0f) evi = 0.0f;
        if (eui<evi) eui = evi;
        if (theta!=null) theta[i2][i1] = asin(u2i);
        if (u1!=null) u1[i2][i1] = u1i;
        if (u2!=null) u2[i2][i1] = u2i;
        if (v1!=null) v1[i2][i1] = v1i;
        if (v2!=null) v2[i2][i1] = v2i;
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
   * Applies this filter to estimate 3-D structure tensors.
   * @param x input array for 3-D image.
   * @param compressed true, for compressed tensors; false, otherwise.
   * @return structure tensors.
   */
  public EigenTensors3 applyForTensors(
    float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] v1, float[][][] v2, float[][][] v3,
    float[][][] w1, float[][][] w2, float[][][] w3,
    float[][][] ep, float[][][] el,
    float[][][] fx) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] au = fillfloat(0.5f,n1,n2,n3);
    float[][][] av = fillfloat(1.0f,n1,n2,n3);
    float[][][] aw = fillfloat(1.0f,n1,n2,n3);
    EigenTensors3 d = new EigenTensors3(u1,u2,w1,w2,au,av,aw,true);
    
    float[][][] gu = new float[n3][n2][n1];
    float[][][] gv = new float[n3][n2][n1];
    float[][][] gw = new float[n3][n2][n1];

    float[][][] guu = new float[n3][n2][n1];
    float[][][] guv = new float[n3][n2][n1];
    float[][][] guw = new float[n3][n2][n1];
    float[][][] gvv = new float[n3][n2][n1];
    float[][][] gvw = new float[n3][n2][n1];
    float[][][] gww = new float[n3][n2][n1];
    computeOrientGradient(u1,u2,u3,v1,v2,v3,w1,w2,w3,fx,gu,gv,gw);
    computeGradientProducts(gu,gv,gw,guu,guv,guw,gvv,gvw,gww);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,20,guu,guu);
    lsf.apply(d,20,guv,guv);
    lsf.apply(d,20,guw,guw);
    lsf.apply(d,20,gvv,gvv);
    lsf.apply(d,20,gvw,gvw);
    lsf.apply(d,20,gww,gww);

    float[][][] x1 = new float[n3][n2][n1];
    float[][][] x2 = new float[n3][n2][n1];
    float[][][] x3 = new float[n3][n2][n1];
    float[][][] z1 = new float[n3][n2][n1];
    float[][][] z2 = new float[n3][n2][n1];
    float[][][] z3 = new float[n3][n2][n1];
    float[][][] a1 = new float[n3][n2][n1];
    float[][][] a2 = new float[n3][n2][n1];
    float[][][] a3 = new float[n3][n2][n1];
    float[][][] c1 = new float[n3][n2][n1];
    float[][][] c2 = new float[n3][n2][n1];
    float[][][] c3 = new float[n3][n2][n1];
    solveEigenproblems(guu,guv,guw,gvv,gvw,gww,
      null,null,x1,x2,x3,null,null,null,z1,z2,z3,au,av,aw,ep,el);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];

      float v1i = v1[i3][i2][i1];
      float v2i = v2[i3][i2][i1];
      float v3i = v3[i3][i2][i1];

      float w1i = w1[i3][i2][i1];
      float w2i = w2[i3][i2][i1];
      float w3i = w3[i3][i2][i1];

      float x1i = x1[i3][i2][i1];
      float x2i = x2[i3][i2][i1];
      float x3i = x3[i3][i2][i1];

      float z1i = z1[i3][i2][i1];
      float z2i = z2[i3][i2][i1];
      float z3i = z3[i3][i2][i1];

      if (u3i<0f) {
        u1i = -u1i;
        u2i = -u2i;
        u3i = -u3i;
      }
      if (w3i<0f) {
        w1i = -w1i;
        w2i = -w2i;
        w3i = -w3i;
      }
      if (x3i<0f) {
        x1i = -x1i;
        x2i = -x2i;
        x3i = -x3i;
      }
      if (z3i<0f) {
        z1i = -z1i;
        z2i = -z2i;
        z3i = -z3i;
      }


      float a1i = u1i*x1i+u2i*x2i+u3i*x3i;
      float a2i = v1i*x1i+v2i*x2i+v3i*x3i;
      float a3i = w1i*x1i+w2i*x2i+w3i*x3i;

      float c1i = u1i*z1i+u2i*z2i+u3i*z3i;
      float c2i = v1i*z1i+v2i*z2i+v3i*z3i;
      float c3i = w1i*z1i+w2i*z2i+w3i*z3i;
      if (a3i<0f) {
        a1i = -a1i;
        a2i = -a2i;
        a3i = -a3i;
      }
      if (c3i<0f) {
        c1i = -c1i;
        c2i = -c2i;
        c3i = -c3i;
      }
      a1[i3][i2][i1] = a1i;
      a2[i3][i2][i1] = a2i;
      c1[i3][i2][i1] = c1i;
      c2[i3][i2][i1] = c2i;
      u1[i3][i2][i1] = u1i;
      u2[i3][i2][i1] = u2i;
      u3[i3][i2][i1] = u3i;
      w1[i3][i2][i1] = w1i;
      w2[i3][i2][i1] = w2i;
      w3[i3][i2][i1] = w3i;
    }}}
    return new EigenTensors3(u1,u2,w1,w2,au,av,aw,true);
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
    float[][][] u1, float[][][] u2, float[][][] u3, 
    float[][][] v1, float[][][] v2, float[][][] v3, 
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] eu, float[][][] ev, float[][][] ew, 
    float[][][] ep, float[][][] el)
  {
    // Where possible, use output arrays for workspace.
    float[][][][] t = new float[16][][][];
    int nt = 0;
    if (theta!=null) t[nt++] = theta;
    if (phi!=null) t[nt++] = phi;
    if (u1!=null) t[nt++] = u1;
    if (u2!=null) t[nt++] = u2;
    if (u3!=null) t[nt++] = u3;
    if (v1!=null) t[nt++] = v1;
    if (v2!=null) t[nt++] = v2;
    if (v3!=null) t[nt++] = v3;
    if (w1!=null) t[nt++] = w1;
    if (w2!=null) t[nt++] = w2;
    if (w3!=null) t[nt++] = w3;
    if (eu!=null) t[nt++] = eu;
    if (ev!=null) t[nt++] = ev;
    if (ew!=null) t[nt++] = ew;
    if (ep!=null) t[nt++] = ep;
    if (el!=null) t[nt++] = el;

    // Gradient.
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    float[][][] g1 = (nt>0)?t[0]:new float[n3][n2][n1];
    float[][][] g2 = (nt>1)?t[1]:new float[n3][n2][n1];
    float[][][] g3 = (nt>2)?t[2]:new float[n3][n2][n1];
    _rgfGradient1.apply100(x,g1);
    _rgfGradient2.apply010(x,g2);
    _rgfGradient3.apply001(x,g3);

    // Gradient products.
    float[][][] g11 = g1;
    float[][][] g22 = g2;
    float[][][] g33 = g3;
    float[][][] g12 = (nt>3)?t[3]:new float[n3][n2][n1];
    float[][][] g13 = (nt>4)?t[4]:new float[n3][n2][n1];
    float[][][] g23 = (nt>5)?t[5]:new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    /*
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float g1i = g1[i3][i2][i1];
          float g2i = g2[i3][i2][i1];
          float g3i = g3[i3][i2][i1];
          g11[i3][i2][i1] = g1i*g1i;
          g22[i3][i2][i1] = g2i*g2i;
          g33[i3][i2][i1] = g3i*g3i;
          g12[i3][i2][i1] = g1i*g2i;
          g13[i3][i2][i1] = g1i*g3i;
          g23[i3][i2][i1] = g2i*g3i;
        }
      }
    }
    */
    
    // Smoothed gradient products comprise the structure tensor.
    if (_rgfSmoother1!=null || _rgfSmoother2!=null || _rgfSmoother3!=null) {
      float[][][] h = (nt>6)?t[6]:new float[n3][n2][n1];
      float[][][][] gs = {g11,g22,g33,g12,g13,g23};
      for (float[][][] g:gs) {
        if (_rgfSmoother1!=null) {
          _rgfSmoother1.apply0XX(g,h);
        } else {
          copy(g,h);
        }
        if (_rgfSmoother2!=null) {
          _rgfSmoother2.applyX0X(h,g);
        } else {
          copy(h,g);
        }
        if (_rgfSmoother3!=null) {
          _rgfSmoother3.applyXX0(g,h);
          copy(h,g);
        }
      }
    }

    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(g11,g12,g13,g22,g23,g33,
      theta,phi,u1,u2,u3,v1,v2,v3,w1,w2,w3,eu,ev,ew,ep,el);
    /*
    float[][] a = new float[3][3];
    float[][] z = new float[3][3];
    float[] e = new float[3];
    for (int i3=0; i3<n3; ++i3) {
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
          float u1i = z[0][0];
          float u2i = z[0][1];
          float u3i = z[0][2];
          float v1i = z[1][0];
          float v2i = z[1][1];
          float v3i = z[1][2];
          if (u1i<0.0f) {
            u1i = -u1i;
            u2i = -u2i;
            u3i = -u3i;
          }
          if (v2i<0.0f) {
            v1i = -v1i;
            v2i = -v2i;
            v3i = -v3i;
          }
          float w1i = u2i*v3i-u3i*v2i;
          float w2i = u3i*v1i-u1i*v3i;
          float w3i = u1i*v2i-u2i*v1i;
          float eui = e[0];
          float evi = e[1];
          float ewi = e[2];
          if (ewi<0.0f) ewi = 0.0f;
          if (evi<ewi) evi = ewi;
          if (eui<evi) eui = evi;
          if (theta!=null) theta[i3][i2][i1] = acos(u1i);
          if (phi!=null) phi[i3][i2][i1] = atan2(u3i,u2i);
          if (u1!=null) u1[i3][i2][i1] = u1i;
          if (u2!=null) u2[i3][i2][i1] = u2i;
          if (u3!=null) u3[i3][i2][i1] = u3i;
          if (v1!=null) v1[i3][i2][i1] = v1i;
          if (v2!=null) v2[i3][i2][i1] = v2i;
          if (v3!=null) v3[i3][i2][i1] = v3i;
          if (w1!=null) w1[i3][i2][i1] = w1i;
          if (w2!=null) w2[i3][i2][i1] = w2i;
          if (w3!=null) w3[i3][i2][i1] = w3i;
          if (eu!=null) eu[i3][i2][i1] = eui;
          if (ev!=null) ev[i3][i2][i1] = evi;
          if (ew!=null) ew[i3][i2][i1] = ewi;
          if (ep!=null || el!=null) {
            float esi = (eui>0.0f)?1.0f/eui:1.0f;
            if (ep!=null) ep[i3][i2][i1] = (eui-evi)*esi;
            if (el!=null) el[i3][i2][i1] = (evi-ewi)*esi;
          }
        }
      }
    }
    */
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private RecursiveGaussianFilter _rgfGradient1;
  private RecursiveGaussianFilter _rgfGradient2;
  private RecursiveGaussianFilter _rgfGradient3;
  private RecursiveGaussianFilter _rgfSmoother1;
  private RecursiveGaussianFilter _rgfSmoother2;
  private RecursiveGaussianFilter _rgfSmoother3;

  public void computeOrientGradient(
    final float[][][] u1, final float[][][] u2, final float[][][] u3, 
    final float[][][] v1, final float[][][] v2, final float[][][] v3, 
    final float[][][] w1, final float[][][] w2, final float[][][] w3, 
    final float[][][] fx, final float[][][] gu, final float[][][] gv, 
    final float[][][] gw) 
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
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];

        float v1i = v1[i3][i2][i1];
        float v2i = v2[i3][i2][i1];
        float v3i = v3[i3][i2][i1];

        float w1i = w1[i3][i2][i1];
        float w2i = w2[i3][i2][i1];
        float w3i = w3[i3][i2][i1];

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

        float fup = si.interpolate(s1,s2,s3,fx,u1p,u2p,u3p);
        float fum = si.interpolate(s1,s2,s3,fx,u1m,u2m,u3m);
        float fvp = si.interpolate(s1,s2,s3,fx,v1p,v2p,v3p);
        float fvm = si.interpolate(s1,s2,s3,fx,v1m,v2m,v3m);
        float fwp = si.interpolate(s1,s2,s3,fx,w1p,w2p,w3p);
        float fwm = si.interpolate(s1,s2,s3,fx,w1m,w2m,w3m);
        gu[i3][i2][i1] = fup-fum;
        gv[i3][i2][i1] = fvp-fvm;
        gw[i3][i2][i1] = fwp-fwm;
      }}
    }});
  }

  private void computeStructureOrientedGradientX(
    final EigenTensors3 et, 
    final float[][][] fx, final float[][][] p2, final float[][][] p3, 
    final float[][][] g2, final float[][][] g3)
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xm0 = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xp0 = new float[n1];
      float[] gm0 = new float[n1];
      float[] g0m = new float[n1];
      float[] g0p = new float[n1];
      float[] gp0 = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fxm0 = fx[i3m][i2 ];
        float[] fxp0 = fx[i3p][i2 ];
        float[] fx0m = fx[i3 ][i2m];
        float[] fx0p = fx[i3 ][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p20p = p2[i3 ][i2p];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3p0 = p3[i3p][i2 ];
        float[] g232 = g2[i3 ][i2 ];
        float[] g332 = g3[i3 ][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          x0m[i1] = i1-p20m[i1];
          x0p[i1] = i1+p20p[i1];
          xm0[i1] = i1-p3m0[i1];
          xp0[i1] = i1+p3p0[i1];
        }
        si.interpolate(n1,1.0,0.0,fx0m,n1,x0m,g0m);
        si.interpolate(n1,1.0,0.0,fx0p,n1,x0p,g0p);
        si.interpolate(n1,1.0,0.0,fxm0,n1,xm0,gm0);
        si.interpolate(n1,1.0,0.0,fxp0,n1,xp0,gp0);
        float[] f00 = fx[i3][i2];
        if (i2==0   ) g0m = fx[i3][i2];
        if (i2==n2-1) g0p = f00;
        if (i3==0   ) gm0 = f00;
        if (i3==n3-1) gp0 = f00;
        for (int i1=0; i1<n1; ++i1) {
          g232[i1] = g0p[i1]-g0m[i1];
          g332[i1] = gp0[i1]-gm0[i1];
        }
      }
    }});
  }


  private void computeStructureOrientedGradientX(
    final float[][][] fx, final float[][][] p2, final float[][][] p3, 
    final float[][][] g2, final float[][][] g3)
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xm0 = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xp0 = new float[n1];
      float[] gm0 = new float[n1];
      float[] g0m = new float[n1];
      float[] g0p = new float[n1];
      float[] gp0 = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fxm0 = fx[i3m][i2 ];
        float[] fxp0 = fx[i3p][i2 ];
        float[] fx0m = fx[i3 ][i2m];
        float[] fx0p = fx[i3 ][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p20p = p2[i3 ][i2p];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3p0 = p3[i3p][i2 ];
        float[] g232 = g2[i3 ][i2 ];
        float[] g332 = g3[i3 ][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          x0m[i1] = i1-p20m[i1];
          x0p[i1] = i1+p20p[i1];
          xm0[i1] = i1-p3m0[i1];
          xp0[i1] = i1+p3p0[i1];
        }
        si.interpolate(n1,1.0,0.0,fx0m,n1,x0m,g0m);
        si.interpolate(n1,1.0,0.0,fx0p,n1,x0p,g0p);
        si.interpolate(n1,1.0,0.0,fxm0,n1,xm0,gm0);
        si.interpolate(n1,1.0,0.0,fxp0,n1,xp0,gp0);
        float[] f00 = fx[i3][i2];
        if (i2==0   ) g0m = fx[i3][i2];
        if (i2==n2-1) g0p = f00;
        if (i3==0   ) gm0 = f00;
        if (i3==n3-1) gp0 = f00;
        for (int i1=0; i1<n1; ++i1) {
          g232[i1] = g0p[i1]-g0m[i1];
          g332[i1] = gp0[i1]-gm0[i1];
        }
      }
    }});
  }


  private void computStructureOrientedGradient(
    final float[][][] fx, final float[][][] p2, final float[][][] p3, 
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
        float[] p232 = p2[i3][i2];
        float[] p332 = p3[i3][i2];
        float[] g232 = g2[i3][i2];
        float[] g332 = g3[i3][i2];
        for (int i1=0; i1<n1;++i1) {
          float p2i = p232[i1];
          float p3i = p332[i1];
          float p2s = 1f/sqrt(1f+p2i*p2i);
          float p3s = 1f/sqrt(1f+p3i*p3i);
          float x2m = i2-p2s;
          float x2p = i2+p2s;
          float x3m = i3-p3s;
          float x3p = i3+p3s;
          float z2m = i1-p2i*p2s;
          float z2p = i1+p2i*p2s;
          float z3m = i1-p3i*p3s;
          float z3p = i1+p3i*p3s;
          float x2i = i2;
          float x3i = i3;
          float f2m = si.interpolate(s1,s2,s3,fx,z2m,x2m,x3i);
          float f2p = si.interpolate(s1,s2,s3,fx,z2p,x2p,x3i);
          float f3m = si.interpolate(s1,s2,s3,fx,z3m,x2i,x3m);
          float f3p = si.interpolate(s1,s2,s3,fx,z3p,x2i,x3p);
          g232[i1] = (f2p-f2m);
          g332[i1] = (f3p-f3m);
        }
      }
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
    final float[][][] u1, final float[][][] u2, final float[][][] u3, 
    final float[][][] v1, final float[][][] v2, final float[][][] v3, 
    final float[][][] w1, final float[][][] w2, final float[][][] w3, 
    final float[][][] eu, final float[][][] ev, final float[][][] ew, 
    final float[][][] ep, final float[][][] el)
  {
    final int n1 = g11[0][0].length;
    final int n2 = g11[0].length;
    final int n3 = g11.length;
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
            float v1i = (float)z[1][0];
            float v2i = (float)z[1][1];
            float v3i = (float)z[1][2];
            float w1i = (float)z[2][0];
            float w2i = (float)z[2][1];
            float w3i = (float)z[2][2];
            if (u1i<0.0f) {
              u1i = -u1i;
              u2i = -u2i;
              u3i = -u3i;
            }
            if (v2i<0.0f) {
              v1i = -v1i;
              v2i = -v2i;
              v3i = -v3i;
            }
            if (w3i<0.0f) {
              w1i = -w1i;
              w2i = -w2i;
              w3i = -w3i;
            }
            float eui = (float)e[0];
            float evi = (float)e[1];
            float ewi = (float)e[2];
            if (ewi<0.0f) ewi = 0.0f;
            if (evi<ewi) evi = ewi;
            if (eui<evi) eui = evi;
            if (theta!=null) theta[i3][i2][i1] = acos(u1i);
            if (phi!=null) phi[i3][i2][i1] = atan2(u3i,u2i);
            if (u1!=null) u1[i3][i2][i1] = u1i;
            if (u2!=null) u2[i3][i2][i1] = u2i;
            if (u3!=null) u3[i3][i2][i1] = u3i;
            if (v1!=null) v1[i3][i2][i1] = v1i;
            if (v2!=null) v2[i3][i2][i1] = v2i;
            if (v3!=null) v3[i3][i2][i1] = v3i;
            if (w1!=null) w1[i3][i2][i1] = w1i;
            if (w2!=null) w2[i3][i2][i1] = w2i;
            if (w3!=null) w3[i3][i2][i1] = w3i;
            if (eu!=null) eu[i3][i2][i1] = eui;
            if (ev!=null) ev[i3][i2][i1] = evi;
            if (ew!=null) ew[i3][i2][i1] = ewi;
            if (ep!=null || el!=null) {
              float esi = (eui>0.0f)?1.0f/eui:1.0f;
              if (ep!=null) ep[i3][i2][i1] = (eui-evi)*esi;
              if (el!=null) el[i3][i2][i1] = (evi-ewi)*esi;
            }
          }
        }
      }
    });
  }

  public float[][] smooth(float sig1, float sig2, float[][] w, float[][] u) {
    int n2 = w.length;
    int n1 = w[0].length;
    float[][] b = new float[n2][n1];
    float[][] r = new float[n2][n1];
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(sig1,sig2);
    A2 a2 = new A2(smoother2,w);
    CgSolver cs = new CgSolver(0.001,200);
    mul(w,u,b);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    return r;
  }

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] wp) 
    {
      _s2 = s2;
      _wp = wp;
      float n2 = wp.length;
      float n1 = wp[0].length;
      _sc = 2f*sum(wp)/(n1*n2);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      _s2.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s2.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother2 _s2;
    private float[][] _wp;
  }

  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    int n2 = wp.length;
    int n1 = wp[0].length;
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      y[i2][i1] += wp[i2][i1]*x[i2][i1];
  }

  private static void addAndScale(float sc, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      y[i2][i1] += sc*x[i2][i1];
    }}
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(float sigma1, float sigma2) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,x);
    }
    private float _sigma1,_sigma2;
  }

  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

}
