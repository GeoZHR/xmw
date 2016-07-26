package mef;

import java.io.Serializable;
import static edu.mines.jtk.util.ArrayMath.*;
import static mef.FaultPointGeometry.*;

/**
 * A 2D fault point is an oriented point located on a fault. Fault points can
 * be linked to form fault curves, which may be used to analyze faults.
 * <p>
 * Fault points are computed from images of fault likelihoods, and
 * dips. Each fault point is an oriented point located on a ridge in an image
 * of fault likelihood. Fault cells have indices (i1,i2) that indicate
 * which image sample is nearest to this ridge. In this way an image sample
 * is associated with either no point or one point.
 * <p>
 * A fault point has up to two neighbors ("nabors") that lie above and below.
 * Links to nabors enables cells to form a curve of connected points, 
 * which represents a fault.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.21
 */


public class FaultPoint implements Serializable {
  private static final long serialVersionUID = 1L;

  /**
   * Gets the fault likelihood for this point.
   * @return the fault likelihood.
   */
  public float getFl() {
    return fl;
  }

  public void setFl(float fl) {
    this.fl = fl;
  }



  public static void getFlImage(FaultPoint[] points, float[][] fl) {
    for (FaultPoint point:points) {
      int i1 = point.i1;
      int i2 = point.i2;
      fl[i2][i1] = point.fl;
    }
  }

  public static void getFtImage(FaultPoint[] points, float[][] ft) {
    for (FaultPoint point:points) {
      int i1 = point.i1;
      int i2 = point.i2;
      ft[i2][i1] = point.ft;
    }
  }

  public void setUnfaultShifts(float[] ts) {
    this.t1=ts[0];
    this.t2=ts[1];
  }


  /////////////////////////////////////////////////////////////////////////
  // package

  int i1,i2; // point indices
  FaultPoint fc; // corresponding point
  float x1,x2; // point coordinates
  float fl,ft; // likelihood, strike (phi) and dip (theta)
  float u1,u2,us; // dip vector and scale factor = 1/sin(theta)
  float w1,w2; // normal vector
  FaultPoint ca,cb; // nabors above, below
  FaultPoint cc; //corresponding point
  FaultCurve curve; // if not null, the curve to which this point belongs
  int i2m,i2p; // sample indices i2 for minus and plus sides of point
  float[] emp; // array of minus-plus alignment errors
  float smp; // shift from minus side to plus side of point
  float s1,s2; // fault dip-slip vector
  float t1,t2; // fault dip-slip vector

  interface Get1 { public float get(FaultPoint point); }
  interface GetN { public float[] get(FaultPoint point); }
  interface Set1 { public void set(FaultPoint point, float value); }
  interface SetN { public void set(FaultPoint point, float[] values); }

  public FaultPoint(float x1, float x2, float fl, float ft) {
    set(x1,x2,fl,ft);
  }


    /////////////////////////////////////////////////////////////////////////
  // private

  private void set(float x1, float x2, float fl, float ft) {
    this.x1 = x1; 
    this.x2 = x2; 
    this.fl = fl; 
    this.ft = ft;
    i1 = round(x1);
    i2 = round(x2);
    float[] u = faultDipVectorFromDip(ft);
    float[] w = faultNormalVectorFromDip(ft);
    u1 = u[0]; u2 = u[1]; us = 1.0f/u1;
    w1 = w[0]; w2 = w[1]; 

    // Indices (i2m,i2p) and (i3m,i3p) for minus-plus pairs of samples.
    // Cell normal vector w points from the minus side to the plus side.
    i2m = i2p = i2;
    if (x2>i2) {
      ++i2p;
    } else if (x2<i2) {
      --i2m;
    }
    if ((i2p-i2m)*w2<0.0f) {
      int i2t = i2m; 
      i2m = i2p; 
      i2p = i2t;
    }
  }

    /**
   * Returns the distance squared from this point to a sample.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return distance squared.
   */
  float distanceSquaredTo(float p1, float p2) {
    float d1 = p1-x1;
    float d2 = p2-x2;
    return d1*d1+d2*d2;
  }


  private static float distanceSquared(
      FaultPoint c, float p1, float p2) {
    return c!=null ? c.distanceSquaredTo(p1,p2) : Float.MAX_VALUE;
  }


  FaultPoint walkUpDipFrom(float[] p) {
    FaultPoint point = this;
    float p1 = p[0];
    float p2 = p[1];
    //assert abs(point.distanceFromPlaneTo(p1,p2,p3))<0.01f;

    // If a point below is found, project point horizontally onto its plane.
    FaultPoint ca = point.ca;
    if (ca!=null) {
      point = ca;
      p1 = point.x1;
      p2 = point.x2;
    }

    // Return updated point and point.
    //assert abs(point.distanceFromPlaneTo(p1,p2,p3))<0.01f;
    p[0] = p1;
    p[1] = p2;
    return point;
  }

  /**
   * Walks a point down the fault along a curve tangent to fault dip.
   * The input point is assumed to lie in the plane of this point,
   * and the output point will lie in the plane of the returned point.
   * That returned point will be one immediately below this point, if
   * sufficient point nabors exist, or this point, otherwise.
   * @param p input and output array {p1,p2,p3} of point coordinates.
   * @return the point with a plane that contains the output point.
   */
  FaultPoint walkDownDipFrom(float[] p) {
    FaultPoint point = this;
    float p1 = p[0];
    float p2 = p[1];
    //assert abs(point.distanceFromPlaneTo(p1,p2,p3))<0.01f;

    // If a point below is found, project point horizontally onto its plane.
    FaultPoint cb = point.cb;
    if (cb!=null) {
      point = cb;
      p1 = point.x1;
      p2 = point.x2;
    }

    // Return updated point and point.
    //assert abs(point.distanceFromPlaneTo(p1,p2,p3))<0.01f;
    p[0] = p1;
    p[1] = p2;
    return point;
  }






}
