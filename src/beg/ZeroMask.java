package beg;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import he.*;

/**
 * A mask for samples that are zero or near zero.
 * Values in the mask image are either true or false.
 * Samples for which the mask is false may be unreliable
 * in some applications.
 * For example, where the mask is false, we might set structure 
 * tensors to default tensors that represent horizontal layering.
 * <p> 
 * Algorithm:
 * (1) Compute global mean absolute amplitude (gabs) for the entire image.
 * (2) Compute local mean absolute amplitude (labs) in Gaussian windows.
 * Mask is zero for all samples where labs is less than small*gabs.
 */
public class ZeroMask {

  /**
   * Constructs a zero mask.
   * @param small small value; zeros in mask where labs &lt; small*gabs.
   * @param sigma1 Gaussian window half-width for 1st dimension.
   * @param sigma2 Gaussian window half-width for 2nd dimension.
   * @param sigma3 Gaussian window half-width for 3rd dimension.
   * @param x array of image values from which mask is derived.
   */
  public ZeroMask(
    double small, double sigma1, double sigma2, double sigma3,
    float[][][] x) 
  {
    _n1 = x[0][0].length;
    _n2 = x[0].length;
    _n3 = x.length;
    float[][][] t = abs(x);
    float a = ((sum(t)/_n1)/_n2)/_n3; // global mean absolute amplitude
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(sigma1);
    RecursiveGaussianFilter rgf2 = new RecursiveGaussianFilter(sigma2);
    RecursiveGaussianFilter rgf3 = new RecursiveGaussianFilter(sigma3);
    float[][][] b = zerofloat(_n1,_n2,_n3);
    rgf1.apply0XX(t,b);
    rgf2.applyX0X(b,t);
    rgf3.applyXX0(t,b); // local mean absolute amplitude
    _mask = new boolean[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (b[i3][i2][i1]<small*a) {
            _mask[i3][i2][i1] = false;
          } else {
            _mask[i3][i2][i1] = true;
          }
        }
      }
    }
  }

  /**
   * Constructs a zero mask from specified array of floats.
   * Mask is true for all non-zero samples in the array.
   * @param m array of values from which mask is derived.
   */
  public ZeroMask(float[][][] m) {
    _n1 = m[0][0].length;
    _n2 = m[0].length;
    _n3 = m.length;
    _mask = new boolean[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (m[i3][i2][i1]!=0.0f)
            _mask[i3][i2][i1] = true;
        }
      }
    }
  }

  /**
   * Constructs a zero mask from a water bottom surface.
   * Control points are specified to extract the water bottom surface.
   * Mask is true for all samples above the water bottom surface.
   * @param k1 array of 1st coordinate of control points.
   * @param k2 array of 2nd coordinate of control points.
   * @param k3 array of 3rd coordinate of control points.
   * @param p2 array of inline slopes.
   * @param p3 array of crossline slopes.
   * @param ep array of planarities.
   * @param gx array of input seismic image.
   */

  public ZeroMask(float[] k1, float[] k2, float[] k3,
    float[][][] p2, float[][][] p3, float[][][] ep,float[][][] gx) 
  {
    _n3 = p2.length;
    _n2 = p2[0].length;
    _n1 = p2[0][0].length;
    SurfaceExtractorC se = new SurfaceExtractorC();
    se.setCG(0.01f,200);
    se.setExternalIterations(20);
    se.setSmoothings(6.0f,6.0f);
    float lmt = _n1-1.f;
    k1 = se.refineConstraints(k1,k2,k3,gx);
    float[][] surf = se.surfaceInitialization(_n2,_n3,lmt,k1,k2,k3);
    se.surfaceUpdateFromSlopes(ep,p2,p3,k1,k2,k3,surf);
    _mask = new boolean[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
      int n1b = round(surf[i3][i2])-5;
      if(n1b<0){n1b=0;}
      for (int i1=0; i1<_n1; ++i1) {
        if(i1<n1b){_mask[i3][i2][i1]=false;}
        else      {_mask[i3][i2][i1]=true;}
      }
    }}
  }

  /**
   * Returns an array of floats representing this mask.
   * @return array of floats.
   */
  public float[][][] getAsFloats() {
    float[][][] f = new float[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3)
      for (int i2=0; i2<_n2; ++i2)
        for (int i1=0; i1<_n1; ++i1)
          f[i3][i2][i1] = (_mask[i3][i2][i1])?1.0f:0.0f;
    return f;
  }

  /**
   * Applies this mask to a specified array of values.
   * @param vfalse value to use where mask is false.
   * @param v array of values to be masked.
   */
  public void setValue(float vfalse, float[][][] v) {
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (!_mask[i3][i2][i1]) {
            v[i3][i2][i1] = vfalse;
          }
        }
      }
    }
  }

  /**
   * Applies this mask to a specified eigentensor field.
   * @param efalse eigentensor to use where mask is false.
   * @param e eigentensors to be masked.
   */
  public void setValue(float[] efalse, EigenTensors3 e) {
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (!_mask[i3][i2][i1])
            e.setTensor(i1,i2,i3,efalse);
        }
      }
    }
  }

  private int _n1,_n2,_n3;
  private boolean[][][] _mask;
}
