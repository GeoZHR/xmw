/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package hv;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * From a precomputed horizon volume x1(u1,x2,x3) or 
 * a relative geologic time volume u1(x1,x2,x3), seismic horizons can be easily 
 * extracted by simply horizontal slicing in the horizon volume at each relative 
 * geologic time u1
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.05.27
 */

public class HorizonFromRgt {
  /**
   * Sets horizon volume and relative geologc time volume for horizon extraction.
   * @param x1 an array of horizon volume x1(u1,x2,x3).
   * @param u1 an array of relative geologic time volume u1(x1,x2,x3).
   */
  public HorizonFromRgt(float[][][] x1, float[][][] u1) {
    _x1 = x1;
    _u1 = u1;
  }

  public float[][] singleHorizon(float u1i) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    float top = 0;
    float down = n1-1;
    int i1 = round(u1i);
    float[][] sfi = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      sfi[i3][i2] = _x1[i3][i2][i1];
    }}
    float[] xyz = trianglesForSurface(sfi,top,down);
    int nt = xyz.length;
    float[] u1s = fillfloat(u1i,nt);
    float rgtMin = min(_u1);
    float rgtMax = max(_u1);
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);
    float[] rgb = cp.getRgbFloats(u1s);
    return new float[][]{xyz,rgb};
  }

  public float[] trianglesForSurface(float[][] surf, float top, float down) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    float[] xyz = new float[n1*n2*18];
    int i = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float surfi = surf[i2][i1];
        if (surfi<down && surfi>top && i1<n1-1 && i2<n2-1) {
          xyz[i   ] = (float) i2;
          xyz[i+1 ] = (float) i1;
          xyz[i+2 ] = surf[i2][i1];
          xyz[i+3 ] = (float) i2+1;
          xyz[i+4 ] = (float) i1;
          xyz[i+5 ] = surf[i2+1][i1];
          xyz[i+6 ] = (float) i2+1;
          xyz[i+7 ] = (float) i1+1;
          xyz[i+8 ] = surf[i2+1][i1+1];

          xyz[i+9 ] = (float) i2;
          xyz[i+10] = (float) i1;
          xyz[i+11] = surf[i2][i1];
          xyz[i+12] = (float) i2+1;
          xyz[i+13] = (float) i1+1;
          xyz[i+14] = surf[i2+1][i1+1];
          xyz[i+15] = (float) i2;
          xyz[i+16] = (float) i1+1;
          xyz[i+17] = surf[i2][i1+1];
          i = i+18;
        } 
      } 
    }
    return copy(i,0,xyz);
  }

  private float[][][] _x1;
  private float[][][] _u1;
}
