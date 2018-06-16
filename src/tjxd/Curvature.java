
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

package tjxd;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Methods of computing
 * 1) seismic horizon surface curvature, 
 * 2) seismic volumetric curvature and 
 * 3) horizon volumetric curvature.
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2018.05.10
 */

public class Curvature {
  /**
   * Computes horizon volumetric curvature
   * @param hv input seismic horizon volume
   * @return most positive and negative curvatures
   */
  public float[][][][] horizonVolumetricCurvature (float sigma, float[][][] hv) {
    int n3 = hv.length;
    int n2 = hv[0].length;
    int n1 = hv[0][0].length;
    float[][][] cn = new float[n3][n2][n1];
    float[][][] cp = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    float[][][] a = new float[n3][n2][n1];
    float[][][] b = new float[n3][n2][n1];
    float[][][] c = new float[n3][n2][n1];
    for (int i1=0; i1<n1; ++i1) {
      float s1i = 0.0f;
      for (int i3=0; i3<n3; ++i3) 
      for (int i2=0; i2<n2; ++i2) 
        s1i += hv[i3][i2][i1];
      s1i /= n2*n3;
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float hvi = hv[i3][i2][i1]-s1i;
        hv[i3][i2][i1] = hvi;
      }}
    }
    rgf.apply020(hv,a);
    rgf.apply002(hv,b);
    rgf.apply010(hv,c);
    rgf.apply001(c,c);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] cp3 = cp[i3];
      float[][] cn3 = cn[i3];
      float[][] a3 = a[i3];
      float[][] b3 = b[i3];
      float[][] c3 = c[i3];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float abp = a3[i2][i1]+b3[i2][i1];
        float abm = a3[i2][i1]-b3[i2][i1];
        float cis = c3[i2][i1]*c3[i2][i1];
        float abc = sqrt(abm*abm+cis);
        cp3[i2][i1] = abp+abc;
        cn3[i2][i1] = abp-abc;
      }}
    }});
    return new float[][][][]{cn,cp};
  }

}
