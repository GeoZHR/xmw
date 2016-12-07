/****************************************************************************
Copyright 2009, Colorado School of Mines and others.
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
package sis;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Tensor-guided blended neighbor gridding in 2D.
 * Gridding is interpolation of a set of known sample values on a 
 * uniformly sampled grid. Here the interpolation is performed by 
 * a two-step blended neighbor process described by Hale (2009).
 * <p>
 * The first step is to compute for all samples the distance to the 
 * nearest known sample and the value of that known sample. This first
 * step produces a distance map and a nearest-neighbor interpolant.
 * <p>
 * The second step is to blend (smooth) the nearest-neighbor interpolant,
 * where the extent of smoothing varies spatially and is proportional to 
 * distances in the distance map.
 * <p>
 * In tensor-guided gridding, we replace distance with time. Time is a
 * simple term for non-Euclidean distance measured in a metric-tensor
 * field. So "nearest" now means nearest in time. In the first step we 
 * compute a time map by solving an eikonal equation with coefficients 
 * that may be both anisotropic and spatially varying. In the second 
 * step, we blend the nearest-neighbor interpolant with an anisotropic 
 * and spatially varying smoothing filter.
 * <p>
 * The default tensor field is homogeneous and isotropic. In this
 * special case, time is equivalent to distance, and tensor-guided
 * gridding is similar to gridding with Sibson's natural neighbor 
 * interpolant.
 * <p>
 * Reference: 
 * <a href="http://www.mines.edu/~dhale/papers/Hale09ImageGuidedBlendedNeighborInterpolation.pdf">
 * Hale, D., 2009, Image-guided blended neighbor interpolation, CWP-634</a>
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.07.21
 */
public class GeodesicSegmentation {

  public void setTensors(float[][] w) {
    int n2 = w.length;
    int n1 = w[0].length;
    _tensors = new EigenTensors2(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wi = w[i2][i1];
      _tensors.setTensor(i1,i2,wi,0f,wi);
    }}
  }

  /**
   * Sets the maximum time computed by this gridder. The gridder has 
   * linear precision where times are less than the maximum time.
   * @param tmax the maximum time.
   */
  public void setTimeMax(double tmax) {
    _tmax = (float)tmax;
  }

  /**
   * Experimental use only.
   */
  public double getTimeMarkerS() {
    return _tms;
  }

  public float[][] applyForAmplitude(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[n2][n1];
    HilbertTransformFilter htf = new HilbertTransformFilter();
    for (int i2=0; i2<n2; ++i2) {
      float[] ft = new float[n1];
      htf.apply(n1,fx[i2],ft);
      for (int i1=0; i1<n1; ++i1) {
        float fti = ft[i1];
        float fxi = fx[i2][i1];
        gx[i2][i1] = sqrt(fti*fti+fxi*fxi);
      }
    }
    return gx;
  }

  public float[][] apply(
    float[] fx1, float[] fx2, 
    float[] bx1, float[] bx2, float[][] w) {
    int n2 = w.length;
    int n1 = w[0].length;
    int nf = fx1.length;
    float[][] pf = new float[n2][n1];
    for (int k=0; k<nf; ++k) {
      int i1 = round(fx1[k]);
      int i2 = round(fx2[k]);
      pf[i2][i1] = 1f;
    }
    setTensors(sub(1f,w));
    float[][] ft = gridNearest(0f,pf);
    int nb = bx1.length;
    float[][] pb = new float[n2][n1];
    for (int k=0; k<nb; ++k) {
      int i1 = round(bx1[k]);
      int i2 = round(bx2[k]);
      pb[i2][i1] = 1f;
    }
    setTensors(w);
    float[][] fb = gridNearest(0f,pb);
    float[][] mk = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fti = ft[i2][i1];
      float fbi = fb[i2][i1];
      if (fti<fbi) {
        mk[i2][i1] = 1f;
      }
    }}
    return mk;
  }

  /**
   * Computes gridded values using nearest neighbors.
   * Gridded values in the array p are computed for only unknown 
   * samples with value equal to the specified null value. Any
   * known (non-null) sample values in the array p are not changed.
   * <p>
   * This method also computes and returns an array of times to
   * nearest-neighbor samples. Times are zero for known samples 
   * and are positive for gridded samples.
   * @param pnull the null value representing unknown samples.
   * @param p array of sample values to be gridded.
   * @return array of times to nearest known samples.
   */
  public float[][] gridNearest(float pnull, float[][] p) {
    int n1 = p[0].length;
    int n2 = p.length;

    // Make array of times, zero for samples with non-null values.
    // Also, count the number of marks that we will need.
    int nmark = 0;
    float[][] t = new float[n2][n1];
    float tnull = -FLT_MAX;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (p[i2][i1]==pnull) {
          t[i2][i1] = tnull;
        } else {
          t[i2][i1] = 0.0f;
          ++nmark;
        }
      }
    }
    gridNearest(nmark,t,p);
    return t;
  }

  /**
   * Computes gridded values using nearest neighbors.
   * Gridded values in the array p are computed for only unknown 
   * samples denoted by corresponding non-zero times in the array t. 
   * This method does not change known values in p, which correspond
   * to zero times in t.
   * @param t array of times to nearest known samples.
   * @param p array of nearest-neighbor gridded values.
   */
  public void gridNearest(float[][] t, float[][] p) {
    int n1 = t[0].length;
    int n2 = t.length;

    // Count the known samples, the number of marks we need.
    int nmark = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (t[i2][i1]==0.0f)
          ++nmark;
      }
    }
    gridNearest(nmark,t,p);
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _tms; // time marker CPU time in seconds

  private EigenTensors2 _tensors;
  private float _tmax = FLT_MAX;

  private void gridNearest(int nmark, float[][] t, float[][] p) {
    int n1 = t[0].length;
    int n2 = t.length;

    // Make an array for marks, while storing values of known samples
    // in an array of values indexed by the mark.
    float[] pmark = new float[nmark];
    int[][] m = new int[n2][n1];
    int mark = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (t[i2][i1]==0.0f) {
          pmark[mark] = p[i2][i1];
          m[i2][i1] = mark;
          ++mark;
        }
      }
    }

    // Use the time marker to compute both times and marks.
    edu.mines.jtk.util.Stopwatch sw = new edu.mines.jtk.util.Stopwatch();
    TimeMarker2 tm = new TimeMarker2(n1,n2,_tensors);
    //tm.setConcurrency(TimeMarker2.Concurrency.SERIAL);
    sw.start();
    tm.apply(t,m);
    sw.stop();
    _tms = sw.time();

    // Adjust times to ensure interpolation of known samples.
    adjustTimes(nmark,m,t);

    // Use the marks to compute the nearest-neighbor interpolant.
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float ti = t[i2][i1];
        if (ti!=0.0f)
          p[i2][i1] = pmark[m[i2][i1]];
        if (ti>_tmax)
          t[i2][i1] = _tmax;
      }
    }
  }

  // Adjusts times to be nearly zero in neighborhood of known samples.
  private void adjustTimes(int nmark, int[][] m, float[][] t) {
    int n1 = m[0].length;
    int n2 = m.length;
    int n1m = n1-1;
    int n2m = n2-1;
    float[] s = new float[nmark];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = (i2>0  )?i2-1:i2;
      int i2p = (i2<n2m)?i2+1:i2;
      for (int i1=0; i1<n1; ++i1) {
        int i1m = (i1>0  )?i1-1:i1;
        int i1p = (i1<n1m)?i1+1:i1;
        if (t[i2][i1]==0.0f) {
          float tmax = 0.0f;
          tmax = max(tmax,t[i2m][i1m]);
          tmax = max(tmax,t[i2m][i1 ]);
          tmax = max(tmax,t[i2m][i1p]);
          tmax = max(tmax,t[i2 ][i1m]);
          tmax = max(tmax,t[i2 ][i1p]);
          tmax = max(tmax,t[i2p][i1m]);
          tmax = max(tmax,t[i2p][i1 ]);
          tmax = max(tmax,t[i2p][i1p]);
          s[m[i2][i1]] = tmax;
        }
      }
    }
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (t[i2][i1]>0.0f)
          t[i2][i1] = max(FLT_MIN,t[i2][i1]-s[m[i2][i1]]);
      }
    }
  }
}
