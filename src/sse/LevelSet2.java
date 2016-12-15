/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package sse;


import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * 2D level set method
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */

public class LevelSet2 {
  /**
   * Initialize a 2D level-set function with a center point and radias.
   * This method first computes some sparse points at a circle defined 
   * by the center (c1,c2) and radius, then places zero values at these 
   * points and pairs of negative and positive values, respectively, 
   * inside and outside the circle in directions normal to the circle.
   * A radial function interpolation method is applied to interpolate an 
   * initial level-set function using these points defined with zeros, 
   * negative and positive values.
   * @param n1 the 1st dimension number.
   * @param n2 the 2nd dimension number.
   * @param c1 the x coordinate of the center.
   * @param c2 the y coordinate of the center.
   * @param r  the radius.
   */
  public float[][] initialLeveSet(
    Sampling s1, Sampling s2, float c1, float c2, float r) {
    float dp = 20f;
    int np = round(360f/dp);
    float[] x1 = new float[np*3];
    float[] x2 = new float[np*3];
    float[] fx = new float[np*3];
    for (int ip=0; ip<np; ip++) {
      float phi = (float)toRadians(ip*dp);
      float sph = sin(phi);
      float cph = cos(phi);
      float x1c = c1+r*sph;
      float x2c = c1+r*cph;
      float x1p = c1+(r+5f)*sph;
      float x2p = c1+(r+5f)*cph;
      float x1m = c1+(r-5f)*sph;
      float x2m = c1+(r-5f)*cph;
      int kp = ip*3;
      x1[kp] = x1c; x1[kp+1] = x1p; x1[kp+2] = x1m;
      x2[kp] = x2c; x2[kp+1] = x2p; x1[kp+2] = x2m;
      fx[kp] = 0.f; fx[kp+1] = 10f; fx[kp+2] = -10;
    }
    RadialInterpolator2.Biharmonic basis = new RadialInterpolator2.Biharmonic();
    RadialInterpolator2 ri = new RadialInterpolator2(basis,fx,x1,x2);
    return ri.interpolate(s1,s2); 
  }

  public void updateLevelSet(float[][] u1, float[][] u2, float[][] fx) {
    int n2 = u1.length;
    int n1 = u1[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {

    }}

  }

}
