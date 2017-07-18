/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package spv;

import java.io.Serializable;

import edu.mines.jtk.io.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;
/**
 * A fault cell is an oriented point located on a fault. 
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.07.15
 */

public class FaultCell {

  public FaultCell(int i1,int i2, int i3, float fl, float fp, float ft) {
    _i1 = i1;
    _i2 = i2;
    _i3 = i3;
    _fl = fl;
    _fp = fp;
    _ft = ft;
  }

  public int getI1() {
    return _i1;
  }
  public int getI2() {
    return _i2;
  }
  public int getI3() {
    return _i3;
  }



  public int[] getIndex() {
    return new int[]{_i1,_i2,_i3};
  }

  public float getFl() {
    return _fl;
  }

  public float getFp() {
    return _fp;
  }

  public float getFt() {
    return _ft;
  }


  /////////////////////////////////////////////////////////////////////////
  // package

  int _i1,_i2,_i3; // cell indices
  float _fl,_fp,_ft; // likelihood, strike (phi) and dip (theta)

}
