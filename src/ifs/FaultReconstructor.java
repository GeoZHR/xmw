/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Reconstruct fault images from fault cells.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.15
 */

public class FaultReconstructor {

  public FaultReconstructor(int n1, int n2, int n3, FaultCell[] fcs){
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fcs = fcs;
  }

  private float[][][] accumulateGaussians(FaultCell[] fcs) {
    float[][][] fg = new float[_n1][_n2][_n3];
    return fg;
  }


  private KdTree setKdTree(FaultCell[] fcs) {
    int nc = fcs.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {

    }
    return new KdTree(xc);
  }

  private int _n1, _n2, _n3;
  private FaultCell[] _fcs;
}


