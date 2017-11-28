/****************************************************************************
Copyright (c) 2009, University of Texas at Austin and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package kmc;

import java.util.List;
import java.util.LinkedList;


/**
 * A vector.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.11.22
 */
public class Point {

  private int _pid; // point ID
  private float[] _x; // location
  private Long _cid = null; // cluster ID
  private List<Integer> _linkPoints = null; // must linked constraints
  private List<Integer> _unlinkPoints = null; // cannot linked constraints

  public Point(float[] x) {
    _x = x;
  }


  public Point(int pid, float[] x) {
    _x = x;
    _pid = pid;
  }

  public Point(int pid, Long cid, float[] x) {
    _pid = pid;
    _cid = cid;
    _x = x;
  }

  public void setClusterId(Long cid) {
    _cid = cid;
  }

  public void addMustLinkPoint(int pid) {
    if(_linkPoints==null)
      _linkPoints = new LinkedList<Integer>();
    _linkPoints.add(pid);
  }

  public void addCannotLinkPoint(int pid) {
    if(_unlinkPoints==null)
      _unlinkPoints = new LinkedList<Integer>();
    _unlinkPoints.add(pid);
  }

  public int getId() {
    return _pid;
  }

  public Long getClusterId() {
    return _cid;
  }

  public float[] getLocation() {
    return _x;
  }

  public int getDimension() {
    return _x.length;
  }

  public List<Integer> getLinkPoints() {
    return _linkPoints;
  }

  public List<Integer> getUnlinkPoints() {
    return _unlinkPoints;
  }

  public float getDistance(Point b) {
    float ds = 0;
    float[] bx = b.getLocation();
    int n = bx.length;
    for (int i=0; i<n; i++) {
      float di = _x[i]-bx[i];
      ds += di*di;
    }
    return (float) Math.sqrt(ds);
  }




}
