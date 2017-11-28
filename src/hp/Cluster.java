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
 * A cluster.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.11.22
 */

public class Cluster {

  private int _id;

  private Point _centroid;

  private List<Point> _points = new LinkedList<Point>();

  public Cluster(int id, Point centroid) {
    _id = id;
    _centroid = centroid;
  }

  public void setCentroid(Point centroid) {
    _centroid = centroid;
  }

  public void addPoint(Point point) {
    _points.add(point);
  }

  public int getId() {
    return _id;
  }

  public Point getCentroid() {
    return _centroid;
  }

  public List<Point> getPoints() {
    return _points;
  }

  public int getSize() {
    return _points.size();
  }

  public void clearPoints() {
    if(_points==null||_points.isEmpty())
      return;
    for (Point point:_points)
      point.setClusterId(-1);
    _points.clear();
  }

  public void removePoint(Point point) {
    _points.remove(point);
  }

}
