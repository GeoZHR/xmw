/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package spv;

import static edu.mines.jtk.util.ArrayMath.*;
/**
 * A fault cell is an oriented point located on a fault. 
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.07.15
 */

public class FaultPoint {

  public FaultPoint(int i1,int i2, int i3, float fl, float fp, float ft) {
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

  public float[] getFaultNormal() {
    return faultNormalVectorFromStrikeAndDip(_fp,_ft);
  }
  public float[] getFaultDipVector() {
    return faultDipVectorFromStrikeAndDip(_fp,_ft);
  }
  public float[] getFaultStrikeVector() {
    return faultStrikeVectorFromStrikeAndDip(_fp,_ft);
  }


      /**
   * Returns fault dip vector for specified strike and dip angles.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {u1,u2,u3} of components for dip vector.
   */
  public static float[] faultDipVectorFromStrikeAndDip(
    double phi, double theta) {
    double p = toRadians(phi);
    double t = toRadians(theta);
    double cp = cos(p);
    double sp = sin(p);
    double ct = cos(t);
    double st = sin(t);
    float u1 = (float)( st);
    float u2 = (float)( ct*cp);
    float u3 = (float)(-ct*sp);
    return new float[]{u1,u2,u3};
  }

  /**
   * Returns fault strike vector for specified strike and dip angles.
   * The dip angle theta is not used, but is provided for consistency.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {v1,v2,v3} of components for strike vector.
   */
  public static float[] faultStrikeVectorFromStrikeAndDip(
      double phi, double theta) {
    double p = toRadians(phi);
    double cp = cos(p);
    double sp = sin(p);
    float v1 = 0.0f;
    float v2 = (float)sp;
    float v3 = (float)cp;
    return new float[]{v1,v2,v3};
  }

  /**
   * Returns fault normal vector for specified strike and dip angles.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {w1,w2,w3} of components for normal vector.
   */
  public static float[] faultNormalVectorFromStrikeAndDip(
      double phi, double theta) {
    double p = toRadians(phi);
    double t = toRadians(theta);
    double cp = cos(p);
    double sp = sin(p);
    double ct = cos(t);
    double st = sin(t);
    float w1 = (float)(-ct);
    float w2 = (float)( st*cp);
    float w3 = (float)(-st*sp);
    return new float[]{w1,w2,w3};
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  int _i1,_i2,_i3; // cell indices
  float _fl,_fp,_ft; // likelihood, strike (phi) and dip (theta)
  boolean _marked = false;

}
