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
 * Set up eigen tensors.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.08.05
 */
public class EigenTensorSetup {

  /**
   * Constructs eigentensors with normal vectors.
   * @param u1 vertical components of reflection normal vectors
   * @param u2 inline components of reflection normal vectors
   * @param u3 crossline components of reflection normal vectors
   */
  public EigenTensorSetup(float[][][] u1, float[][][] u2, float[][][] u3) {
    _u1 = u1;
    _u2 = u2;
    _u3 = u3;
  }

  public EigenTensors3 applyForTensors() {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    return null;
  }

  private float[][][] _u1, _u2, _u3;
}


  
