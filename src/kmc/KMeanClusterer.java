/****************************************************************************
Copyright (c) 2009, University of Texas at Austin and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package kmc;

import java.util.List;
import java.util.Collections;
import java.util.Random;
import java.util.Set;
import java.util.HashSet;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * A cluster.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.11.22
 */

public class KMeanClusterer {

  public void setConvergence(int niter, double convergeDistance) {
    _niter = niter;
    _convergeDistance = convergeDistance;
  }

  public Cluster[] applyClustering(Cluster[] clusters, List<Point> points) {
    int iter = 0;
    assignPointsToClusters(clusters,points);
    while (iter<_niter) {
      System.out.println("iter="+iter);
      Cluster[] updatedClusters = updateClusterCentroids(clusters);
      while (!assignPointsToClusters(updatedClusters,points));
      if (converge(clusters,updatedClusters))
        return updatedClusters;
      clusters = updatedClusters;
      iter++;
    }
    return clusters;
  }


  public Cluster[] applyClustering(int nc, List<Point> points) {
    Cluster[] clusters = findInitialClusters(nc, points);
 		while (!assignPointsToClusters(clusters, points)) {
			clusters = findInitialClusters(nc, points);
		}
    int niter = 0;
    while (niter<_niter) {
      Cluster[] updatedClusters = updateClusterCentroids(clusters);
      while (!assignPointsToClusters(updatedClusters,points));
      if (converge(clusters,updatedClusters))
        return updatedClusters;
      clusters = updatedClusters;
      niter++;
    }
    return clusters;
  }

  private boolean converge(Cluster[] previous, Cluster[] current) {
    int nc = previous.length;
    for (int ic=0; ic<nc; ++ic) {
      Point cc = current[ic].getCentroid();
      Point cp = previous[ic].getCentroid();
      if(cc.getDistance(cp)>_convergeDistance)
        return false;
    }
    return true;
  }


  private Cluster[] findInitialClusters (int nc, List<Point> points) {
    Cluster[] clusters = new Cluster[nc];
		Random random = new Random(System.currentTimeMillis());
		Set<Integer> clusterCenters = new HashSet<Integer>();
    int np = points.size();
		for (int ic=0; ic<nc; ic++) {
			int ip = random.nextInt(np);
			while (clusterCenters.contains(ip)) {
				ip = random.nextInt(np);
			}
			clusterCenters.add(ip);
			clusters[ic] = new Cluster((long)ic, points.get(ip));
		}
		return clusters;
  }

  private Cluster[] updateClusterCentroids(Cluster[] clusters) {
    int nc = clusters.length;
    Cluster[] newClusters = new Cluster[nc];
    for (int ic=0; ic<nc; ++ic) {
      Cluster ci = clusters[ic];
      Point cc = ci.getCentroid();
      int nx = cc.getDimension();
      List<Point> cp = ci.getPoints();
      if (cp!=null && !cp.isEmpty()) {
        //cc = getAverageCentroid(nx,cp);
        cc = getMeanCentroid(nx,cp);
      }
      newClusters[ic] = new Cluster((long)ic,cc);
    }
    return newClusters;
  }

  private Point getAverageCentroid(int nx, List<Point> points) {
    int np = points.size();
    float[] xc = new float[nx];
    for (Point point:points) {
      float[] xp = point.getLocation();
      for (int ix=0; ix<nx; ++ix)
        xc[ix] += xp[ix];
    }
    div(xc,(float)np,xc);
    return new Point(xc);
  }

  private Point getMeanCentroid(int nx, List<Point> points) {
    int np = points.size();
    float[] xc = new float[nx];
    float[][] xps = new float[nx][np];
    int ip = 0;
    for (Point point:points) {
      float[] xp = point.getLocation();
      for (int ix=0; ix<nx; ++ix)
        xps[ix][ip] = (float)xp[ix];
      ip++;
    }
    MedianFinder mf = new MedianFinder(np);
    for (int ix=0; ix<nx; ++ix)
      xc[ix] = mf.findMedian(xps[ix]);
    return new Point(xc);
  }


  private boolean assignPointsToClusters(Cluster[] clusters, List<Point> points) {
    for (Cluster cluster:clusters)
      cluster.clearPoints();
    Collections.shuffle(points);
    for (Point point:points) {
      Cluster nearestCluster = null;
      double minDistance = Float.MAX_VALUE;
      for (Cluster cluster:clusters) {
        Long cid = cluster.getId();
        Point centroid = cluster.getCentroid();
        double ds = centroid.getDistance(point);
        if (ds < minDistance) {
          //if(!violateConstraints(cid,point,points)) {
            minDistance = ds;
            nearestCluster = cluster;
          //}
        }
      }
      if (nearestCluster == null)
        return false;
      nearestCluster.addPoint(point);
      point.setClusterId(nearestCluster.getId());
    }
    return true;
  }

  private boolean violateConstraints(Long cid, Point point, List<Point> points) {
    List<Integer> unlinkPoints = point.getUnlinkPoints();
    if(unlinkPoints!=null) {
      int np = unlinkPoints.size();
      for (int ip=0; ip<np; ++ip) {
        int ulpid = unlinkPoints.get(ip);
        Long ulpCid = points.get(ulpid).getClusterId();
        if(ulpCid!=null&&ulpCid==cid)
          return true;
      }
    }
    return false;
  }


  private int _niter = 100;
  private double _convergeDistance=10;

}
