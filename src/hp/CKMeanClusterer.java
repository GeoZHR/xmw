/****************************************************************************
Copyright (c) 2009, University of Texas at Austin and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package hp;

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

public class CKMeanClusterer {

  public void setConvergence(int niter, double convergeDistance) {
    _niter = niter;
    _convergeDistance = convergeDistance;
  }

  public Cluster[] applyClustering(
    boolean dist, int[][] mks, Cluster[] clusters, List<Point> points) {
    int iter = 0;
    int n2 = mks.length;
    int nc = clusters.length;
    int[][] cmark = fillint(-1,nc,n2);
    assignPointsToClusters(dist,cmark,clusters,points);
    while (iter<_niter) {
      System.out.println("iter="+iter);
      Cluster[] updatedClusters = updateClusterCentroids(clusters);
      cmark = fillint(-1,nc,n2);
      //while (!assignPointsToClusters(cmark,updatedClusters,points));
      assignPointsToClusters(dist,cmark,updatedClusters,points);
      if (converge(clusters,updatedClusters))
        return updatedClusters;
      clusters = updatedClusters;
      iter++;
    }
    return clusters;
  }

  /*

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
  */

  private boolean converge(Cluster[] previous, Cluster[] current) {
    int nc = previous.length;
    for (int ic=0; ic<nc; ++ic) {
      Point cc = current[ic].getCentroid();
      Point cp = previous[ic].getCentroid();
      if(cc.getDistancex(cp)>_convergeDistance)
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
			clusters[ic] = new Cluster(ic, points.get(ip));
		}
		return clusters;
  }

  private Cluster[] updateClusterCentroids(Cluster[] clusters) {
    int nc = clusters.length;
    Cluster[] newClusters = new Cluster[nc];
    for (int ic=0; ic<nc; ++ic) {
      Cluster ci = clusters[ic];
      Point cc = ci.getCentroid();
      float[][] vs = cc.getFeaturesx();
      int nx = vs.length;
      int nd = vs[0].length;
      List<Point> cp = ci.getPoints();
      if (cp!=null && !cp.isEmpty()) {
        cc = getMeanCentroid(nx,nd,cp);
      }
      newClusters[ic] = new Cluster(ic,cc);
    }
    return newClusters;
  }

  private Point getAverageCentroid(int nx, List<Point> points) {
    int np = points.size();
    float[] vc = new float[nx];
    for (Point point:points) {
      float[] vp = point.getFeatures();
      for (int ix=0; ix<nx; ++ix)
        vc[ix] += vp[ix];
    }
    div(vc,(float)np,vc);
    return new Point(vc);
  }

  private Point getMeanCentroid(int nx, List<Point> points) {
    int np = points.size();
    float[] vc = new float[nx];
    float[][] vps = new float[nx][np];
    int ip = 0;
    for (Point point:points) {
      float[] vp = point.getFeatures();
      for (int ix=0; ix<nx; ++ix)
        vps[ix][ip] = (float)vp[ix];
      ip++;
    }
    MedianFinder mf = new MedianFinder(np);
    for (int ix=0; ix<nx; ++ix)
      vc[ix] = mf.findMedian(vps[ix]);
    return new Point(vc);
  }

  private Point getMeanCentroid(int nx, int nd, List<Point> points) {
    int np = points.size();
    float[][] vc = new float[nx][nd];
    float[][][] vps = new float[nx][nd][np];
    int ip = 0;
    for (Point point:points) {
      float[][] vp = point.getFeaturesx();
      for (int ix=0; ix<nx; ++ix)
      for (int id=0; id<nd; ++id)
        vps[ix][id][ip] = vp[ix][id];
      ip++;
    }
    MedianFinder mf = new MedianFinder(np);
    for (int ix=0; ix<nx; ++ix)
    for (int id=0; id<nd; ++id)
      vc[ix][id] = mf.findMedian(vps[ix][id]);
    return new Point(vc);
  }


  private boolean assignPointsToClustersx(
    boolean dist, int[][] cmark, Cluster[] clusters, List<Point> points) {
    int nc = clusters.length;
    for (Cluster cluster:clusters)
      cluster.clearPoints();
    Collections.shuffle(points);
    for (int ic=0; ic<nc; ++ic) {
      Cluster cluster = clusters[ic];
      Point centroid = cluster.getCentroid();
      for (Point point:points) {
      }
    }


    for (Point point:points) {
      int p2 = point.getI2();
      Cluster nearestCluster = null;
      float maxSemblance = Float.MIN_VALUE;
      float minDistance = Float.MAX_VALUE;
      int cn = 0;
      for (int ic=0; ic<nc; ++ic) {
        Cluster cluster = clusters[ic];
        Point centroid = cluster.getCentroid();
        //if (cmark[p2][ic]==-1) {
        if(dist) {
          float ds = centroid.getDistancex(point);
          if (ds<minDistance) {
            minDistance = ds;
            nearestCluster = cluster;
            cn = ic;
          }
        }else {
          float ds = centroid.getSemblance(point);
          if (ds>maxSemblance) {
            maxSemblance = ds;
            nearestCluster = cluster;
            cn = ic;
          }
        }
        //}
      }
      if (nearestCluster == null)
        return false;
      nearestCluster.addPoint(point);
      cmark[p2][cn] = 1;
      point.setClusterId(nearestCluster.getId());
    }
    return true;
  }


  private boolean assignPointsToClusters(
    boolean dist, int[][] cmark, Cluster[] clusters, List<Point> points) {
    int nc = clusters.length;
    for (Cluster cluster:clusters)
      cluster.clearPoints();
    Collections.shuffle(points);
    for (Point point:points) {
      int p2 = point.getI2();
      Cluster nearestCluster = null;
      float maxSemblance = Float.MIN_VALUE;
      float minDistance = Float.MAX_VALUE;
      int cn = 0;
      for (int ic=0; ic<nc; ++ic) {
        Cluster cluster = clusters[ic];
        Point centroid = cluster.getCentroid();
        //if (cmark[p2][ic]==-1) {
        if(dist) {
          float ds = centroid.getDistancex(point);
          if (ds<minDistance) {
            minDistance = ds;
            nearestCluster = cluster;
            cn = ic;
          }
        }else {
          float ds = centroid.getSemblance(point);
          if (ds>maxSemblance) {
            maxSemblance = ds;
            nearestCluster = cluster;
            cn = ic;
          }
        }
        //}
      }
      if (nearestCluster == null)
        return false;
      nearestCluster.addPoint(point);
      cmark[p2][cn] = 1;
      point.setClusterId(nearestCluster.getId());
    }
    return true;
  }

  private boolean assignPointsToClusters(
    int[][] cmark, int[][] mks, Cluster[] clusters, List<Point> points) {
    int nc = clusters.length;
    for (Point point:points)
      point.setClusterId(-1);

    for (Cluster cluster:clusters)
      cluster.clearPoints();
    Collections.shuffle(points);
    for (Point point:points) {
      int p2 = point.getI2();
      int p1 = point.getI1();
      if(point.getClusterId()>=0) continue;
      Cluster nearestCluster = null;
      double minDistance = Float.MAX_VALUE;
      int cn = 0;
      for (int ic=0; ic<nc; ++ic) {
        Cluster cluster = clusters[ic];
        Point centroid = cluster.getCentroid();
        if (cmark[p2][ic]==0) {
          //double ds = centroid.getDistance(point);
          double ds = centroid.getSemblance(point);
          if (ds < minDistance) {
            minDistance = ds;
            nearestCluster = cluster;
            cn = ic;
          }
        }
      }
      if (nearestCluster==null) return false;
      nearestCluster.addPoint(point);
      cmark[p2][cn] = 1;
      mks[p2][p1] = cn;
      point.setClusterId(nearestCluster.getId());
    }
    return true;
  }


  private int _niter = 100;
  private double _convergeDistance=10;

}
