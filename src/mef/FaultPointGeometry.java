/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Methods for 2D fault geometry in seismic image coordinates.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.21
 */
public class FaultPointGeometry {

  /**
   * Returns fault dip vector for specified dip angle.
   * @param theta fault dip angle, in degrees.
   * @return array {u1,u2} of components for dip vector.
   */
  public static float[] faultDipVectorFromDip(double theta) {
    double t = toRadians(theta);
    double ct = cos(t);
    double st = sin(t);
    float u1 = (float)(st);
    float u2 = (float)(ct);
    if(u1<0f) {u1=-u1;}
    else      {u2=-u2;}
    return new float[]{u1,u2};
  }

  /**
   * Returns fault normal vector for specified dip angle.
   * @param theta fault dip angle, in degrees.
   * @return array {w1,w2} of components for normal vector.
   */
  public static float[] faultNormalVectorFromDip(double theta) {
    double t = toRadians(theta);
    double ct = cos(t);
    double st = sin(t);
    float w1 = (float)(-ct);
    float w2 = (float)(-st);
    return new float[]{w1,w2};
  }

  /**
   * Returns fault dip angle for specified fault dip vector.
   * @param u1 1st component of fault dip vector.
   * @param u2 2nd component of fault dip vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromDipVector(float u1, float u2) {
    float theta = acos(u2);
    if(u2>0f) {theta = -theta;}
    return toDegrees(theta);
  }

  /**
   * Returns fault dip angle for specified fault dip vector.
   * @param u array {u1,u2} of components of fault dip vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromDipVector(float[] u) {
    return faultDipFromDipVector(u[0],u[1]);
  }

  /**
   * Returns fault dip angle for specified fault normal vector.
   * @param w1 1st component of fault normal vector.
   * @param w2 2nd component of fault normal vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromNormalVector(float w1, float w2) {
    float theta = -asin(w2);
    return toDegrees(theta);
  }

  /**
   * Returns fault dip angle for specified fault normal vector.
   * @param w array {w1,w2} of components of fault normal vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromNormalVector(float[] w) {
    return faultDipFromNormalVector(w[0],w[1]);
  }


  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    float[] fp = {  0.0f, 90.0f,180.0f,270.0f,
                    0.0f, 90.0f,180.0f,270.0f};
    float[] ft = { 90.0f, 90.0f, 90.0f, 90.0f,
                   89.0f, 89.0f, 89.0f, 89.0f};
    for (int i=0; i<fp.length; ++i)
      test(ft[i]);
  }
  public static void test(float thetaa) {
    trace("theta="+thetaa);
    float[] ua = faultNormalVectorFromDip(thetaa);
    trace("u1="+ua[0]);
    trace("u2="+ua[1]);
    float thetab = faultDipFromNormalVector(ua);
    trace("thetb="+thetab);
    float[] ub = faultNormalVectorFromDip(thetab);
    assertEqual(ua,ub);
  }
  public static void assertEqual(float x, float y) {
    assert abs(x-y)<0.01f;
  }
  public static void assertEqual(float[] x, float[] y) {
    for (int i=0; i<x.length; ++i)
      assertEqual(x[i],y[i]);
  }
  public static void trace(String s) {
    System.out.println(s);
  }
}
