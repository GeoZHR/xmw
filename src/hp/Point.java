/****************************************************************************
Copyright (c) 2009, University of Texas at Austin and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hp;

import java.util.List;
import java.util.LinkedList;


/**
 * A vector.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.11.22
 */
public class Point {

  private int[] _x; // location
  private float[] _v; //features
  private float[][] _vs; //features
  private int _n; //feature dimension
  private int _c=-1; // cluster ID
  private Point _left =null;
  private Point _right=null;

  public Point(float[] features) {
    _v = features;
    _n = features.length;
  }

  public Point(float[][] features) {
    _vs = features;
    _n  = features.length;
  }



  public Point(int[] x, float[] features) {
    _x = x;
    _v = features;
    _n = features.length;
  }

  public Point(int[] x, float[][] features) {
    _x = x;
    _vs = features;
    _n  = features.length;
  }


  public void setClusterId(int cid) {
    _c = cid;
  }

  public void setLeftPoint(Point left) {
    _left = left;
  }

  public void setRightPoint(Point right) {
    _right = right;
  }

  public int getClusterId() {
    return _c;
  }

  public int[] getLocation() {
    return _x;
  }

  public int getI1() {
    return _x[0];
  }

  public int getI2() {
    return _x[1];
  }

  public int getI3() {
    return _x[2];
  }

  public Point getLeftPoint() {
    return _left;
  }

  public Point getRightPoint() {
    return _right;
  }

  public int getLeftCid() {
    if (_left==null) return -1;
    return _left.getClusterId();
  }

  public int getRightCid() {
    if (_right==null) return -1;
    return _right.getClusterId();
  }

  public float[] getFeatures() {
    return _v;
  }

  public float[][] getFeaturesx() {
    return _vs;
  }


  public int getFeatureDimension() {
    return _n;
  }

  public float getDistance(Point b) {
    float ds = 0;
    float[] bv = b.getFeatures();
    for (int i=0; i<_n; i++) {
      float di = _v[i]-bv[i];
      ds += di*di;
    }
    return (float) Math.sqrt(ds);
  }

  public float getDistancex(Point b) {
    float ds = 0;
    float[][] bvs = b.getFeaturesx();
    int n2 = bvs.length;
    int n1 = bvs[0].length;
    for (int i2=0; i2<n2; i2++) {
      float d1 = 0f;
      for (int i1=0; i1<n1; i1++) {
        float di = _vs[i2][i1]-bvs[i2][i1];
        d1 += di*di;
      }
      ds += d1;
    }
    return (float) Math.sqrt(ds);
  }


  public float getSemblance(Point b) {
    float d = 0f;
    float n = 0f;
    int n1 = _n;
    float[] bv = b.getFeatures();
    for (int i1=0; i1<n1; ++i1) {
      float fi = _v[i1];
      float gi = bv[i1];
      float si = (fi+gi)*0.5f;
      d += si*si;
      float fs = fi*fi;
      float gs = gi*gi;
      n += (fs+gs)*0.5f;
    }
    return d/n;
  }



}
