/****************************************************************************
Copyright 2008, Colorado School of Mines and others.
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
package ad;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.util.UnitSphereSampling;


/**
 * An array of eigen-decompositions of tensors for 3D image processing.
 * Each tensor is a symmetric positive-semidefinite 3-by-3 matrix:
 * <pre><code>
 *     |a11 a12 a13|
 * A = |a12 a22 a23|
 *     |a13 a23 a33|
 * </code></pre>
 * Such tensors can be used to parameterize anisotropic image processing.
 * <p>
 * The eigen-decomposition of the matrix A is
 * <pre><code>
 * A = au*u*u' + av*v*v' + aw*w*w' 
 *   = (au-av)*u*u' + (aw-av)*w*w' + av*I
 * </code></pre>
 * where u, v, and w are orthogonal unit eigenvectors of A. (The notation 
 * u' denotes the transpose of u.) The outer products of eigenvectors are
 * scaled by the non-negative eigenvalues au, av, and aw. The second
 * equation exploits the identity u*u' + v*v' + w*w' = I, and makes
 * apparent the redundancy of the vector v.
 * <p>
 * Only the 1st and 2nd components of the eigenvectors u and w are stored. 
 * Except for a sign, the 3rd components may be computed from the 1st and 
 * 2nd. Because the tensors are independent of the choice of sign, the 
 * eigenvectors u and w are stored with an implied non-negative 3rd 
 * component.
 * <p>
 * Storage may be further reduced by compression, whereby eigenvalues
 * and eigenvectors are quantized. Quantization errors for eigenvalues
 * (au,av,aw) are less than 0.001*(au+av+aw). Quantization errors for 
 * eigenvectors are less than one degree of arc on the unit sphere.
 * Memory required to store each tensor is 12 bytes if compressed, and
 * 28 bytes if not compressed.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.07
 */
public class EigenTensors {

  /**
   * Constructs tensors for specified array dimensions and eigenvalues.
   * The 3rd components of eigenvectors u and w are computed from the 1st 
   * and 2nd components and are assumed to be non-negative.
   * @param u1 array of 1st components of u.
   * @param u2 array of 2nd components of u.
   * @param w1 array of 1st components of w.
   * @param w2 array of 2nd components of w.
   * @param au array of eigenvalues au.
   * @param av array of eigenvalues av.
   * @param aw array of eigenvalues aw.
   * @param compressed true, for compressed tensors; false, otherwise.
   */
  public EigenTensors(
    float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] au, float[][][] av, float[][][] aw)
  {
    _au = au;
    _av = av;
    _aw = aw;
    _u1 = u1;
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    _n3 = n3;
    _n2 = n2;
    _n1 = n1;
    _u1 = u1;
    _u2 = u2;
    _u3 = u3;
    float[][][] v1 = new float[n3][n2][n1];
    float[][][] v2 = new float[n3][n2][n1];
    float[][][] v3 = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          if (u1i<0f) {
            u1i = -u1i;
            u2i = -u2i;
            u3i = -u3i;
          }
          _u1[i3][i2][i1] = u1i;
          _u2[i3][i2][i1] = u2i;
          _u3[i3][i2][i1] = u3i;
          float u12s = u1i*u1i+u2i*u2i;
          float u12r = 1f/sqrt(u12s);
          float u23i = u2i*u3i;
          float u13i = u1i*u3i;
          float u3r = 1f/sqrt(u12s*u12s+u23i*u23i+u13i*u13i);
          v1[i3][i2][i1] =-u2i*u12r;
          v2[i3][i2][i1] = u1i*u12r;
          v3[i3][i2][i1] = 0;
          w1[i3][i2][i1] =-u13i*u3r;
          w2[i3][i2][i1] =-u23i*u3r;
          w3[i3][i2][i1] = u12s*u3r;
        }
      }
    }
    _v1 = v1;
    _v2 = v2;
    _v3 = v3;
    _w1 = w1;
    _w2 = w2;
    _w3 = w3;
  }

  public EigenTensors(
    float[][][] u1, float[][][] u2, float[][][] u3,
    float au, float av, float aw)
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    _n3 = n3;
    _n2 = n2;
    _n1 = n1;
    _au = fillfloat(au,n1,n2,n3);
    _av = fillfloat(av,n1,n2,n3);
    _aw = fillfloat(aw,n1,n2,n3);
    _u1 = u1;
    _u2 = u2;
    _u3 = u3;
    float[][][] v1 = new float[n3][n2][n1];
    float[][][] v2 = new float[n3][n2][n1];
    float[][][] v3 = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          if (u1i<0f) {
            u1i = -u1i;
            u2i = -u2i;
            u3i = -u3i;
          }
          _u1[i3][i2][i1] = u1i;
          _u2[i3][i2][i1] = u2i;
          _u3[i3][i2][i1] = u3i;
          float u12s = u1i*u1i+u2i*u2i;
          float u12r = 1f/sqrt(u12s);
          float u23i = u2i*u3i;
          float u13i = u1i*u3i;
          float u3r = 1f/sqrt(u12s*u12s+u23i*u23i+u13i*u13i);
          v1[i3][i2][i1] =-u2i*u12r;
          v2[i3][i2][i1] = u1i*u12r;
          v3[i3][i2][i1] = 0;
          w1[i3][i2][i1] =-u13i*u3r;
          w2[i3][i2][i1] =-u23i*u3r;
          w3[i3][i2][i1] = u12s*u3r;
        }
      }
    }
    _v1 = v1;
    _v2 = v2;
    _v3 = v3;
    _w1 = w1;
    _w2 = w2;
    _w3 = w3;
  }


  /**
   * Gets tensor elements for specified indices.
   * @param i1 index for 1st dimension.
   * @param i2 index for 2nd dimension.
   * @param i3 index for 3rd dimension.
   * @param a array {a11,a12,a13,a22,a23,a33} of tensor elements.
   */
  public void getTensor(int i1, int i2, int i3, float[] a) {
    float au,av,aw,u1,u2,u3,v1,v2,v3,w1,w2,w3;
    au = _au[i3][i2][i1];
    av = _av[i3][i2][i1];
    aw = _aw[i3][i2][i1];
    u1 = _u1[i3][i2][i1];
    u2 = _u2[i3][i2][i1];
    u3 = _u3[i3][i2][i1];
    v1 = _v1[i3][i2][i1];
    v2 = _v2[i3][i2][i1];
    v3 = _v3[i3][i2][i1];
    w1 = _w1[i3][i2][i1];
    w2 = _w2[i3][i2][i1];
    w3 = _w3[i3][i2][i1];
    a[0] = au*u1*u1+av*v1*v1+aw*w1*w1; // a11
    a[1] = au*u1*u2+av*v1*v2+aw*w1*w2; // a12
    a[2] = au*u1*u3+av*v1*v3+aw*w1*w3; // a13
    a[3] = au*u2*u2+av*v2*v2+aw*w2*w2; // a22
    a[4] = au*u2*u3+av*v2*v3+aw*w2*w3; // a23
    a[5] = au*u3*u3+av*v3*v3+aw*w3*w3; // a33
  }

  public float[] getEigenvectorU(int i1, int i2, int i3) {
    return new float[]{_u1[i3][i2][i1],_u2[i3][i2][i1],_u3[i3][i2][i1]};
  }

  public float[] getEigenvalues(int i1, int i2, int i3) {
    return new float[]{_au[i3][i2][i1],_av[i3][i2][i1],_aw[i3][i2][i1]};
  }



  public void setEigenvalues(float au, float av, float aw) {
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      setEigenvalues(i1,i2,i3,au,av,aw);
    }}}
  }

  public void setEigenvalues(
    int i1, int i2, int i3, float au, float av, float aw) {
    _au[i3][i2][i1] = au;
    _av[i3][i2][i1] = av;
    _aw[i3][i2][i1] = aw;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private static UnitSphereSampling _uss; // for compressing unit vectors

  private boolean _compressed; // true if tensors compressed
  private int _n1,_n2,_n3; // array dimensions
  private float[][][] _au; // au not compressed
  private float[][][] _av; // au not compressed
  private float[][][] _aw; // aw not compressed
  private float[][][] _u1; // u1 not compressed
  private float[][][] _u2; // u2 not compressed
  private float[][][] _u3; // u2 not compressed
  private float[][][] _v1; // w1 not compressed
  private float[][][] _v2; // w2 not compressed
  private float[][][] _v3; // w2 not compressed
  private float[][][] _w1; // w1 not compressed
  private float[][][] _w2; // w1 not compressed
  private float[][][] _w3; // w2 not compressed

}
