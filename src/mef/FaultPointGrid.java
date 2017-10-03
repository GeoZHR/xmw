/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

/**
 * Fault points in a 2D sampling grid. Each grid sample indexed by (i1,i2)
 * contains either one fault point or null. The grid facilitates searches for
 * point nabors in skins and fast iterations along fault traces tangent to
 * fault strike and fault curves tangent to fault dip.
 * <p> 
 * Grid indices need not (and typically do not) begin at zero. Index bounds
 * for a fault point grid are determined by the minima and maxima of indices of
 * points used to construct the grid.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.21
 */

import static edu.mines.jtk.util.ArrayMath.*;

public class FaultPointGrid {

  /**
   * Constructs a fault grid for specified points. Grid index bounds are
   * determined by the minimum and maximum indices of the specified points.
   * @param points array of points to be included in the grid.
   */
  public FaultPointGrid(FaultPoint[] points) {
    int i1min = Integer.MAX_VALUE;
    int i2min = Integer.MAX_VALUE;
    int i1max = -i1min;
    int i2max = -i2min;
    for (FaultPoint point:points) {
      if (point.i1<i1min) i1min = point.i1;
      if (point.i2<i2min) i2min = point.i2;
      if (point.i1>i1max) i1max = point.i1;
      if (point.i2>i2max) i2max = point.i2;
    }
    _j1 = i1min;
    _j2 = i2min;
    _n1 = 1+i1max-i1min;
    _n2 = 1+i2max-i2min;
    _points = new FaultPoint[_n2][_n1];
    for (FaultPoint point:points)
      set(point);
  }

  public FaultPointGrid(int n1, int n2) {
    _j1 = 0;
    _j2 = 0;
    _n1 = n1;
    _n2 = n2;
    _points = new FaultPoint[n2][n1];
  }

  public FaultPointGrid(int n1, int n2, FaultPoint[] points) {
    _j1 = 0;
    _j2 = 0;
    _n1 = n1;
    _n2 = n2;
    _points = new FaultPoint[n2][n1];
    for (FaultPoint point:points)
      set(point);
  }


  /**
   * Gets the number of points in the 1st dimension.
   * @return the number of points.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of points in the 2nd dimension.
   * @return the number of points.
   */
  public int getN2() {
    return _n2;
  }

  /**
   * Gets the lower bound on grid indices in the 1st dimension.
   * @return the lower bound.
   */
  public int getI1Min() {
    return _j1;
  }

  /**
   * Gets the lower bound on grid indices in the 2nd dimension.
   * @return the lower bound.
   */
  public int getI2Min() {
    return _j2;
  }

  /**
   * Gets the upper bound on grid indices in the 1st dimension.
   * @return the upper bound.
   */
  public int getI1Max() {
    return _j1+_n1-1;
  }

  /**
   * Gets the upper bound on grid indices in the 2nd dimension.
   * @return the upper bound.
   */
  public int getI2Max() {
    return _j2+_n2-1;
  }


  /**
   * Gets the fault point with specified indices, if any.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @return the fault point; null, if none or if indices are out of bounds.
   */
  public FaultPoint get(int i1, int i2) {
    i1 -= _j1; 
    i2 -= _j2; 
    if (0<=i1 && i1<_n1 && 
        0<=i2 && i2<_n2) {
      return _points[i2][i1];
    } else {
      return null;
    }
  }


  /**
   * Sets the specified fault point. Uses the point's {x1,x2,x3} coordinates to
   * determine the indices of the point in this grid.
   * @param point the fault point.
   */
  public void set(FaultPoint point) {
    int i1 = point.i1-_j1;
    int i2 = point.i2-_j2;
    i1 = max(i1,0); i1 = min(i1,_n1-1);
    i2 = max(i2,0); i2 = min(i2,_n2-1);
    if(_points[i2][i1]==null){
      _points[i2][i1] = point;
    }

  }

  /**
   * Finds a fault point above the specified point. Searches for a point above
   * that lies nearest to the line containing the specified point and its dip
   * vector. If the specified point is already linked to a nabor point above,
   * this method skips the search and simply returns that nabor point. 
   * @param point the point for which to find a point above.
   * @return the point above; null, if none.
   */
  public FaultPoint findPointAbove(FaultPoint point) {
    if (point==null) return null;
    if (point.ca!=null) return point.ca;
    return findPointAboveBelow(true,point);
  }

  /**
   * Finds a fault point below the specified point. Searches for a point below
   * that lies nearest to the line containing the specified point and its dip
   * vector. If the specified point is already linked to a nabor point below,
   * this method skips the search and simply returns that nabor point. 
   * @param point the point for which to find a point below.
   * @return the point below; null, if none.
   */
  public FaultPoint findPointBelow(FaultPoint point) {
    if (point==null) return null;
    if (point.cb!=null) return point.cb;
    return findPointAboveBelow(false,point);
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _j1,_j2; // min point indices
  private int _n1,_n2; // numbers of points
  private FaultPoint[][] _points; // array of points

  private FaultPoint findPointAboveBelow(boolean above, FaultPoint point) {
    int i1 = point.i1;
    int i2 = point.i2;
    float x1 = point.x1;
    float x2 = point.x2;
    float u1 = point.u1;
    float u2 = point.u2;
    int k1 = 1;
    if (above) {
      k1 = -k1;
      u1 = -u1;
      u2 = -u2;
    }
    FaultPoint cmin = null;
    float dmin = Float.MAX_VALUE;
      for (int k2=-5; k2<=5; ++k2) {
        FaultPoint c = get(i1+k1,i2+k2);
        if (c!=null) {
          float d1 = c.x1-x1;
          float d2 = c.x2-x2;
          float du = d1*u1+d2*u2;
          if (du>0.0f) {
            d1 -= du*u1;
            d2 -= du*u2;
            float d = d1*d1+d2*d2; // squared distance to dip line
            if (d<dmin) {
              cmin = c;
              dmin = d;
            }
          }
        }
      }
    return cmin;
  }

}
