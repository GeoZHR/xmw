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
package mef;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import java.util.*;

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
public class FaultOrientEstimator {

    /**
   * Constructs a filter with a possibly anisotropic Gaussian window.
   * @param sigma1 half-width of window in 1st dimension.
   * @param sigma2 half-width of window in 2nd and higher dimensions.
   */
  public FaultOrientEstimator(double sigma1, double sigma2) {
    this(sigma1,sigma2,sigma2);
  }

  /**
   * Constructs a filter with a possibly anisotropic Gaussian window.
   * @param sigma1 half-width of window in 1st dimension.
   * @param sigma2 half-width of window in 2nd dimension.
   * @param sigma3 half-width of window in 3rd and higher dimensions.
   */
  public FaultOrientEstimator(double sigma1, double sigma2, double sigma3) {
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
  }

  public float[][] thin(float fm, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n2][n1];
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply00(fx,fs);
    for (int i2=1;i2<n2-1;i2++) {
    for (int i1=0;i1<n1  ;i1++) {
      float fsi = fs[i2][i1];
      float fsm = fs[i2-1][i1];
      float fsp = fs[i2+1][i1];
      if (fsm<fsi&&fsp<fsi&&fsi>fm) {
        ft[i2  ][i1] = fsi;
      }
    }}
    return ft;
  }

  public float[][][] thin(float fm, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] ft = new float[n3][n2][n1];
    float[][][] fs = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply000(fx,fs);
    for (int i3=0;i3<n3  ;i3++) {
    for (int i2=1;i2<n2-1;i2++) {
    for (int i1=0;i1<n1  ;i1++) {
      float fsi = fs[i3][i2  ][i1];
      float fsm = fs[i3][i2-1][i1];
      float fsp = fs[i3][i2+1][i1];
      if (fsm<fsi&&fsp<fsi&&fsi>fm) {
        ft[i3][i2][i1] = fsi;
      }
    }}}
    for (int i3=1;i3<n3-1;i3++) {
    for (int i2=0;i2<n2  ;i2++) {
    for (int i1=0;i1<n1  ;i1++) {
      float fsi = fs[i3  ][i2][i1];
      float fsm = fs[i3-1][i2][i1];
      float fsp = fs[i3+1][i2][i1];
      if (fsm<fsi&&fsp<fsi&&fsi>fm) {
        ft[i3][i2][i1] = fsi;
      }
    }}}
    return ft;
  }


  public void applyForNormal(
    float d1, float d2, float d3, float[][][] x,
    float[][][] u1, float[][][] u2, float[][][] u3) {
    apply(d1,d2,d3,x,u1,u2,u3,null,null,null,null,null,null);
  }

  public EigenTensors3 applyForTensors(
    float d1, float d2, float d3, float[][][] x, boolean compressed) {
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
    apply(d1,d2,d3,x,
      null,u2,u3,
      null,null,null,
      w1,w2,null);
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
    return new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,compressed);
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
  public void apply(float d1, float d2, float[][] x,
    float[][] u1, float[][] u2) {
    // average position.
    int n2 = x.length;
    int n1 = x[0].length;
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    ArrayList<Float> xva = new ArrayList<Float>();
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float xi = x[i2][i1];
      if (xi>0f) {
        xva.add(xi);
        x1a.add((float)i1);
        x2a.add((float)i2);
      }
    }}
    int np = x1a.size();
    float[] xv = new float[np];
    float[][] xs = new float[2][np];
    for (int ip=0; ip<np; ++ip) {
      xv[ip] = x1a.get(ip);
      xs[0][ip] = x1a.get(ip);
      xs[1][ip] = x2a.get(ip);
    }
    KdTree kd = new KdTree(xs);
    zero(u1);
    zero(u2);
    float[][] t11 = new float[n2][n1];
    float[][] t12 = new float[n2][n1];
    float[][] t22 = new float[n2][n1];
    for (int ip=0; ip<np; ++ip) {
      int i1 = (int)xs[0][ip];
      int i2 = (int)xs[1][ip];
      float[] xmin = new float[]{i1-d1,i2-d2};
      float[] xmax = new float[]{i1+d1,i2+d2};
      int[] id = kd.findInRange(xmin,xmax);
      int nc = id.length;
      if(nc<3) {continue;}
      float c1  = 0f;
      float c2  = 0f;
      float x11 = 0f;
      float x12 = 0f;
      float x22 = 0f;
      for (int ic=0; ic<nc; ++ic) {
        int kp = id[ic];
        float x1 = xs[0][kp];
        float x2 = xs[1][kp];
        x11 += x1*x1;
        x12 += x1*x2;
        x22 += x2*x2;
        c1 += x1;
        c2 += x2;
      }
      x11 -= c1*c1/nc;
      x12 -= c1*c2/nc;
      x22 -= c2*c2/nc;
      float xvi = xv[ip];
      t11[i2][i1] = xvi*x11;
      t12[i2][i1] = xvi*x12;
      t22[i2][i1] = xvi*x22;
    }
    // Smoothed gradient products comprise the structure tensor.
    if (_rgfSmoother1!=null || _rgfSmoother2!=null) {
      float[][] h = new float[n2][n1];
      float[][][] ts = {t11,t22,t12};
      for (float[][] t:ts) {
        if (_rgfSmoother1!=null) {
          _rgfSmoother1.apply0X(t,h);
        } else {
          copy(t,h);
        }
        if (_rgfSmoother2!=null) {
          _rgfSmoother2.applyX0(h,t);
        } else {
          copy(h,t);
        }
      }
    }

    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = t11[i2][i1];
        a[0][1] = t12[i2][i1];
        a[1][0] = t12[i2][i1];
        a[1][1] = t22[i2][i1];
        Eigen.solveSymmetric22(a,z,e);
        float u1i = z[1][0];
        float u2i = z[1][1];
        if (u1i<0.0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        u1[i2][i1] = u1i;
        u2[i2][i1] = u2i;
      }
    }
  }

  public void apply(
    float d1, float d2, float d3, float[][][] x,
    float[][][] u1, float[][][] u2, float[][][] u3, 
    float[][][] v1, float[][][] v2, float[][][] v3, 
    float[][][] w1, float[][][] w2, float[][][] w3)
  {
    // Where possible, use output arrays for workspace.
    float[][][][] t = new float[9][][][];
    int nt = 0;
    if (u1!=null) t[nt++] = u1;
    if (u2!=null) t[nt++] = u2;
    if (u3!=null) t[nt++] = u3;
    if (v1!=null) t[nt++] = v1;
    if (v2!=null) t[nt++] = v2;
    if (v3!=null) t[nt++] = v3;
    if (w1!=null) t[nt++] = w1;
    if (w2!=null) t[nt++] = w2;
    if (w3!=null) t[nt++] = w3;

    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    // Tensor components.
    float[][][] t11 = (nt>0)?t[0]:new float[n3][n2][n1];
    float[][][] t12 = (nt>1)?t[1]:new float[n3][n2][n1];
    float[][][] t13 = (nt>2)?t[2]:new float[n3][n2][n1];
    float[][][] t22 = (nt>3)?t[3]:new float[n3][n2][n1];
    float[][][] t23 = (nt>4)?t[4]:new float[n3][n2][n1];
    float[][][] t33 = (nt>5)?t[5]:new float[n3][n2][n1];
    computeTensors(d1,d2,d3,x,t11,t12,t13,t22,t23,t33);
    if (_rgfSmoother1!=null || _rgfSmoother2!=null || _rgfSmoother3!=null) {
      float[][][] h = (nt>6)?t[6]:new float[n3][n2][n1];
      float[][][][] ts = {t11,t12,t13,t22,t23,t33};
      for (float[][][] ti:ts) {
        if (_rgfSmoother1!=null) {
          _rgfSmoother1.apply0XX(ti,h);
        } else {
          copy(ti,h);
        }
        if (_rgfSmoother2!=null) {
          _rgfSmoother2.applyX0X(h,ti);
        } else {
          copy(h,ti);
        }
        if (_rgfSmoother3!=null) {
          _rgfSmoother3.applyXX0(ti,h);
          copy(h,ti);
        }
      }
    }
    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(t11,t12,t13,t22,t23,t33,u1,u2,u3,v1,v2,v3,w1,w2,w3);
  }


  public void computeTensors(
    final float d1, final float d2, final float d3, final float[][][] x,
    final float[][][] t11, final float[][][] t12, final float[][][] t13, 
    final float[][][] t22, final float[][][] t23, final float[][][] t33) {
    // average position.
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    ArrayList<Float> x3a = new ArrayList<Float>();
    ArrayList<Float> xva = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float xi = x[i3][i2][i1];
      if (xi>0f) {
        xva.add(xi);
        x1a.add((float)i1);
        x2a.add((float)i2);
        x3a.add((float)i3);
      }
    }}}
    final int np = x1a.size();
    final float[] xv = new float[np];
    final float[][] xs = new float[3][np];
    for (int ip=0; ip<np; ++ip) {
      xv[ip] = x1a.get(ip);
      xs[0][ip] = x1a.get(ip);
      xs[1][ip] = x2a.get(ip);
      xs[2][ip] = x3a.get(ip);
    }
    final KdTree kd = new KdTree(xs);
    loop(np,new LoopInt() {
    public void compute(int ip) {
      int i1 = (int)xs[0][ip];
      int i2 = (int)xs[1][ip];
      int i3 = (int)xs[2][ip];
      float[] xmin = new float[]{i1-d1,i2-d2,i3-d3};
      float[] xmax = new float[]{i1+d1,i2+d2,i3+d3};
      int[] id = kd.findInRange(xmin,xmax);
      int nc = id.length;
      if(nc>5) {
        float c1  = 0f;
        float c2  = 0f;
        float c3  = 0f;
        float x11 = 0f;
        float x12 = 0f;
        float x13 = 0f;
        float x22 = 0f;
        float x23 = 0f;
        float x33 = 0f;
        for (int ic=0; ic<nc; ++ic) {
          int kp = id[ic];
          float x1 = xs[0][kp];
          float x2 = xs[1][kp];
          float x3 = xs[2][kp];
          c1 += x1;
          c2 += x2;
          c3 += x3;
          x11 += x1*x1;
          x12 += x1*x2;
          x13 += x1*x3;
          x22 += x2*x2;
          x23 += x2*x3;
          x33 += x3*x3;
        }
        float sc = 1f/nc;
        x11 -= c1*c1*sc;
        x12 -= c1*c2*sc;
        x13 -= c1*c3*sc;
        x22 -= c2*c2*sc;
        x23 -= c2*c3*sc;
        x33 -= c3*c3*sc;
        float xvi = xv[ip];
        t11[i3][i2][i1] = xvi*x11;
        t12[i3][i2][i1] = xvi*x12;
        t13[i3][i2][i1] = xvi*x13;
        t22[i3][i2][i1] = xvi*x22;
        t23[i3][i2][i1] = xvi*x23;
        t33[i3][i2][i1] = xvi*x33;
      }
    }});
  }

  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] u1, final float[][][] u2, final float[][][] u3, 
    final float[][][] v1, final float[][][] v2, final float[][][] v3, 
    final float[][][] w1, final float[][][] w2, final float[][][] w3)
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    loop(n3,new LoopInt() {
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
            float u1i = (float)z[2][0];
            float u2i = (float)z[2][1];
            float u3i = (float)z[2][2];
            float v1i = (float)z[1][0];
            float v2i = (float)z[1][1];
            float v3i = (float)z[1][2];
            float w1i = (float)z[0][0];
            float w2i = (float)z[0][1];
            float w3i = (float)z[0][2];
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
            if (u1!=null) u1[i3][i2][i1] = u1i;
            if (u2!=null) u2[i3][i2][i1] = u2i;
            if (u3!=null) u3[i3][i2][i1] = u3i;
            if (v1!=null) v1[i3][i2][i1] = v1i;
            if (v2!=null) v2[i3][i2][i1] = v2i;
            if (v3!=null) v3[i3][i2][i1] = v3i;
            if (w1!=null) w1[i3][i2][i1] = w1i;
            if (w2!=null) w2[i3][i2][i1] = w2i;
            if (w3!=null) w3[i3][i2][i1] = w3i;
          }
        }
      }
    });
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private RecursiveGaussianFilter _rgfSmoother1;
  private RecursiveGaussianFilter _rgfSmoother2;
  private RecursiveGaussianFilter _rgfSmoother3;

}
