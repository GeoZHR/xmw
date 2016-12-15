/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.List;
import java.util.Iterator;
import java.util.LinkedList;
import edu.mines.jtk.mesh.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * 2D fast level set method.
 * <p>
 * Based on the works by Yonggang Shi and William Clem Karl, 2008, 
 * A Real-Time Algorithm for the Approximation of Level-Set-Based 
 * Curve Evolution.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class FastLevelSet2 {
  /**
   * Construct with a 2D initial square-shape level-set function 
   * with integes -3,-1,1,3.
   * @param n1 the 1st dimension number.
   * @param n2 the 2nd dimension number.
   * @param c1 the x coordinate of the center.
   * @param c2 the y coordinate of the center.
   * @param r  the radius.
   */
  public FastLevelSet2(int n1, int n2, int c1, int c2, int r) {
    _n1 = n1;
    _n2 = n2;
    _phi = fillbyte((byte)3,_n1,_n2);
    int b1 = c1-r;
    int e1 = c1+r;
    int b2 = c2-r;
    int e2 = c2+r;
    e2 = min(e2,_n2-5);
    e1 = min(e1,_n1-5);
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      _phi[i2][i1] = (byte)-3;
    }}
    for (int i2=b2-1; i2<=e2+1; i2++) {
      int i1b = b1-1;
      int i1e = e1+1;
      _phi[i2][i1b] = (byte)1;
      _phi[i2][i1e] = (byte)1;
      Point p1b = new Point(i1b,i2);
      Point p1e = new Point(i1e,i2);
      _lout.add(p1b);
      _lout.add(p1e);
    }
    for (int i1=b1-1; i1<=e1+1; i1++) {
      int i2b = b2-1;
      int i2e = e2+1;
      _phi[i2b][i1] = (byte)1;
      _phi[i2e][i1] = (byte)1;
      Point p2b = new Point(i1,i2b);
      Point p2e = new Point(i1,i2e);
      _lout.add(p2b);
      _lout.add(p2e);
    }
    for (int i2=b2; i2<=e2; i2++) {
      _phi[i2][b1] = (byte)-1;
      _phi[i2][e1] = (byte)-1;
      Point p1b = new Point(b1,i2);
      Point p1e = new Point(e1,i2);
      _lin.add(p1b);
      _lin.add(p1e);
    }
    for (int i1=b1; i1<=e1; i1++) {
      _phi[b2][i1] = (byte)-1;
      _phi[e2][i1] = (byte)-1;
      Point p2b = new Point(i1,b2);
      Point p2e = new Point(i1,e2);
      _lin.add(p2b);
      _lin.add(p2e);
    }
  }

  /**
   * Sets the number of update iterations.
   * @param outerIters  total outer iterations.
   * @param speedIters  speed evolution iterations.
   * @param smoothIters smooth iterations.
   */
  public void setIterations(int outerIters, int speedIters, int smoothIters) {
    _outerIters = outerIters;
    _speedIters = speedIters;
    _smoothIters = smoothIters;
  }

  public void updateLevelSet(int gw, double sigma, float[][] fx) {
    _gWindow = gw;
    createGaussFilter(gw,sigma);
    updateEvolutionSpeedX(fx);
    for (int iter=0; iter<_outerIters; iter++) {
      boolean converged = cycleOne(fx);
      cycleTwo();
      if(converged){break;}
    }
  }

  public float[][] getPhi() {
    float[][] phi = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2)
    for (int i1=0; i1<_n1; ++i1)
      phi[i2][i1] = _phi[i2][i1];
    return phi;
  }

  public float[][] getLout() {
    int np = _lout.size();
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    for (int ip=0; ip<np; ++ip) {
      Point pi = _lout.get(ip);
      x1[ip] = pi.p1;
      x2[ip] = pi.p2;
    }
    return new float[][]{x1,x2};
  }

  public void checkNan(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(Float.isNaN(x[i2][i1])) {
        x[i2][i1] = 0;
      }
    }}
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _n1;
  private int _n2;
  private byte[][] _phi;
  private int _gWindow;
  private byte[][] _gauss;
  private int _gaussThreshold;
  private int _outerIters = 100;
  private int _speedIters = 10;
  private int _smoothIters = 5;
  private TriMesh _tmesh = new TriMesh();
  private List<Point> _lin  = new LinkedList<Point>(); //inside boundary points;
  private List<Point> _lout = new LinkedList<Point>(); //outside boundary points;
  private List<Point> _linAdd  = new LinkedList<Point>(); //inside boundary points;
  private List<Point> _loutAdd  = new LinkedList<Point>(); //inside boundary points;

  private boolean cycleOne(float[][] fx) {
    boolean converged = false;
    for (int iter=0; iter<_speedIters; iter++) {
      //outward evolution
      Iterator<Point> louti = _lout.iterator();
      while(louti.hasNext()) {
        Point pi = louti.next();
        if(pi.fe>0) {switchIn(louti,pi);}
      }
      _lout.addAll(0,_loutAdd);
      _loutAdd.clear();
      cleanLin();
      //inward evolution
      Iterator<Point> lini = _lin.iterator();
      while(lini.hasNext()) {
        Point pi = lini.next();
        if(pi.fe<0) {switchOut(lini,pi);}
      }
      _lin.addAll(0,_linAdd);
      _linAdd.clear();
      cleanLout();
      //update evolution speed and check convergence
      updateEvolutionSpeedX(fx); 
      converged = checkConvergence();
      if(converged) {break;}
    }
    return converged;
  }

  private void cycleTwo() {
    for (int iter=0; iter<_smoothIters; iter++) {
      updateSmoothSpeed();
      //outward evolution
      Iterator<Point> louti = _lout.iterator();
      while(louti.hasNext()) {
        Point pi = louti.next();
        if(pi.fs>_gaussThreshold) {switchIn(louti,pi);}
      }
      _lout.addAll(0,_loutAdd);
      _loutAdd.clear();
      cleanLin();
      //inward evolution
      Iterator<Point> lini = _lin.iterator();
      while(lini.hasNext()) {
        Point pi = lini.next();
        if(pi.fs<_gaussThreshold) {switchOut(lini,pi);}
      }
      _lin.addAll(0,_linAdd);
      _linAdd.clear();
      cleanLout();
    }
  }

  /**
	 * Creates the Gaussian filter matrix, scaled up to integers
	 */
	protected void createGaussFilter(int gw, double sigma) {
		int m2 = 2*gw+1;
		int m1 = 2*gw+1;
		int gfScale = 0;
    double sigmas = 1.0/(sigma*sigma);
		// Rough heuristic: scale by number of elements in filter
		double scale1 = m1*m2;
    _gauss = new byte[m2][m1];
		// In theory could just calculate 1/8th and duplicate instead
		for (int i2=0; i2<m2; ++i2) {
      double d2 = i2-gw;
			for (int i1=0; i1<m1; ++i1) {
        double d1 = i1-gw;
				double ds = d2*d2+d1*d1;
				double gf = sigmas*Math.exp(-0.5*sigmas*ds)*scale1;
        _gauss[i2][i1]=(byte)gf;
				gfScale += gf;
			}
		}
		_gaussThreshold = gfScale/2;
	}

	/**
	 * Update evolution speeed and check convergence.
	 * @return true if convergence has been reached
	 */
  private boolean checkConvergence() {
      //updateEvolutionSpeed();
      //updateEvolutionSpeedXX();
      for (Point pi:_lout) {
        if(pi.fe>0) return false;
      }
      for (Point pi:_lin) {
        if(pi.fe<0) return false;
      }
      return true;
  }

  /**
	 * Update evolution speed fv for pixels in Lout and Lin.
	 */
  private void updateEvolutionSpeed(float[][] u1, float[][] u2) {
    for (Point pi:_lin) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      float u1i = u1[p2i][p1i];
      float u2i = u2[p2i][p1i];
      int p1m = p1i-1; p1m=max(p1m,0);
      int p2m = p2i-1; p2m=max(p2m,0);
      int p1p = p1i+1; p1p=min(p1p,_n1-1);
      int p2p = p2i+1; p2p=min(p2p,_n2-1);
      float d1i = signum(_phi[p2i][p1p]-_phi[p2i][p1m]);
      float d2i = signum(_phi[p2p][p1i]-_phi[p2m][p1i]);
      pi.setEvolutionSpeed((int)signum(d1i*u1i+d2i*u2i));
    }
    for (Point pi:_lout) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      float u1i = u1[p2i][p1i];
      float u2i = u2[p2i][p1i];
      int p1m = p1i-1; p1m=max(p1m,0);
      int p2m = p2i-1; p2m=max(p2m,0);
      int p1p = p1i+1; p1p=min(p1p,_n1-1);
      int p2p = p2i+1; p2p=min(p2p,_n2-1);
      float d1i = signum(_phi[p2i][p1p]-_phi[p2i][p1m]);
      float d2i = signum(_phi[p2p][p1i]-_phi[p2m][p1i]);
      pi.setEvolutionSpeed((int)signum(d1i*u1i+d2i*u2i));
    }
  }


  /**
	 * Use distance map to update evolution speed for pixels 
   * in Lout and Lin.
	 */
  private void updateEvolutionSpeed(float[][] fx) {
    for (Point pi:_lin) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      int b1 = p1i-2;
      int b2 = p2i-2;
      int e1 = p1i+2;
      int e2 = p2i+2;
      b1 = max(b1,0);
      b2 = max(b2,0);
      e1 = min(e1,_n1-1);
      e2 = min(e2,_n2-1);
      float ns = 0f;
      float fs = 0f;
      for (int i2=b2; i2<=e2; ++i2) {
      for (int i1=b1; i1<=e1; ++i1) {
        int phi = _phi[i2][i1];
        if(phi==1) {fs+=fx[i2][i1];ns+=1f;}
      }}
      if(ns>0) {
        pi.setEvolutionSpeed((int)signum(fs/ns-fx[p2i][p1i]));
      } else {
        pi.setEvolutionSpeed(0);
      }
    }
    for (Point pi:_lout) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      int b1 = p1i-2;
      int b2 = p2i-2;
      int e1 = p1i+2;
      int e2 = p2i+2;
      b1 = max(b1,0);
      b2 = max(b2,0);
      e1 = min(e1,_n1-1);
      e2 = min(e2,_n2-1);
      float ns = 0f;
      float fs = 0.0f;
      for (int i2=b2; i2<=e2; ++i2) {
      for (int i1=b1; i1<=e1; ++i1) {
        int phi = _phi[i2][i1];
        if(phi==-1) {fs+=fx[i2][i1];ns+=1f;}
      }}
      if(ns>0) {
        pi.setEvolutionSpeed((int)signum(fx[p2i][p1i]-fs/ns));
      } else {
        pi.setEvolutionSpeed(0);
      }
    }

  }

  private void updateEvolutionSpeedX(float[][] fx) {
    for (Point pi:_lin) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      if (fx[p2i][p1i]>1.2f)
        pi.setEvolutionSpeed(-1);
      else
        pi.setEvolutionSpeed(0);
    }
    for (Point pi:_lout) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      if (fx[p2i][p1i]<0.2f)
        pi.setEvolutionSpeed(1);
      else
        pi.setEvolutionSpeed(0);
    }

  }



  /**
	 * Update smooth speed fs for pixels in Lout and Lin.
	 */
  private void updateSmoothSpeed() {
    for (Point pi:_lin)
      calculateSmoothSpeed(pi);
    for (Point pi:_lout)
      calculateSmoothSpeed(pi);
  }

  	/**
	 * Calculate the smoothing field at a point
	 */
	private void calculateSmoothSpeed(Point p) {
		// Convolve neighbourhood of a point with a gaussian
    int p2 = p.p2;
    int p1 = p.p1;
    int gw = _gWindow;
		int d2m = max(-gw,-p2);
		int d1m = max(-gw,-p1);
		int d2p = min(gw+1,_n2-p2);
		int d1p = min(gw+1,_n1-p1);
		int f = 0;
		for(int d2=d2m; d2<d2p; ++d2) {
		for(int d1=d1m; d1<d1p; ++d1) {
			if (_phi[p2+d2][p1+d1]< 0) 
			f+= _gauss[gw+d2][gw+d1];
		}}
		p.setSmoothSpeed(f);
	}



  /**
	 * Switch a point from Lout to the Lin.
	 * @param louti Iterator to the point to be moved
	 * @param p The point to be moved (needed to update phi)
	 */
	private void switchIn(Iterator<Point> louti, Point p) {
    //step 1: delete the point from lout and add it to lin; set phi=-1
    int p1i = p.p1;
    int p2i = p.p2;
    _lin.add(p);
    _phi[p2i][p1i]=(byte)-1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p2i][p1m]==3) {
      Point pa = new Point(p1m,p2i);
      _loutAdd.add(pa);
      _phi[p2i][p1m] = 1;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p2i][p1p]==3) {
      Point pb = new Point(p1p,p2i);
      _loutAdd.add(pb);
      _phi[p2i][p1p] = 1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p2m][p1i]==3) {
      Point pl = new Point(p1i,p2m);
      _loutAdd.add(pl);
      _phi[p2m][p1i] = 1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p2p][p1i]==3) {
      Point pr = new Point(p1i,p2p);
      _loutAdd.add(pr);
      _phi[p2p][p1i] = 1;
    }
    louti.remove();
	}

  /**
	 * Switch a point from Lin to the Lout.
	 * @param lini Iterator to the point to be moved
	 * @param p The point to be moved (needed to update phi)
	 */
	private void switchOut(Iterator<Point> lini, Point p) {
    //step 1: delete the point from lout and add it to lin; set phi=-1
    int p1i = p.p1;
    int p2i = p.p2;
    _lout.add(p);
    _phi[p2i][p1i]= 1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p2i][p1m]==-3) {
      Point pa = new Point(p1m,p2i);
      _linAdd.add(pa);
      _phi[p2i][p1m] = -1;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p2i][p1p]==-3) {
      Point pb = new Point(p1p,p2i);
      _linAdd.add(pb);
      _phi[p2i][p1p] = -1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p2m][p1i]==-3) {
      Point pl = new Point(p1i,p2m);
      _linAdd.add(pl);
      _phi[p2m][p1i] = -1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p2p][p1i]==-3) {
      Point pr = new Point(p1i,p2p);
      _linAdd.add(pr);
      _phi[p2p][p1i] = -1;
    }
    lini.remove();
	}

  /**
	 * Eliminate redundant point in Lin.
	 */
  private void cleanLin() {
    Iterator<Point> lini = _lin.iterator();
    while (lini.hasNext()) {
      Point pi = lini.next();
      int p1i = pi.p1;
      int p2i = pi.p2;
      int p1m = p1i-1; 
      boolean isRedundant=true;
      //check the neighbor above
      if(p1m>=0&&_phi[p2i][p1m]>0)
        isRedundant = false;
      //check the neighbor below
      int p1p = p1i+1; 
      if(p1p<_n1&&_phi[p2i][p1p]>0)
        isRedundant = false;
      //check the left neighbor
      int p2m = p2i-1; 
      if(p2m>=0&&_phi[p2m][p1i]>0)
        isRedundant = false;
      //check the right neighbor
      int p2p = p2i+1; 
      if(p2p<_n2&&_phi[p2p][p1i]>0) 
        isRedundant = false;
      if(isRedundant) {
        lini.remove();
        _phi[p2i][p1i] = -3;
      }
    }
  }

  /**
	 * Eliminate redundant point in Lout.
	 */
  private void cleanLout() {
    Iterator<Point> louti = _lout.iterator();
    while (louti.hasNext()) {
      Point pi = louti.next();
      int p1i = pi.p1;
      int p2i = pi.p2;
      int p1m = p1i-1; 
      boolean isRedundant=true;
      //check the neighbor above
      if(p1m>=0&&_phi[p2i][p1m]<0)
        isRedundant = false;
      //check the neighbor below
      int p1p = p1i+1; 
      if(p1p<_n1&&_phi[p2i][p1p]<0)
        isRedundant = false;
      //check the left neighbor
      int p2m = p2i-1; 
      if(p2m>=0&&_phi[p2m][p1i]<0)
        isRedundant = false;
      //check the right neighbor
      int p2p = p2i+1; 
      if(p2p<_n2&&_phi[p2p][p1i]<0) 
        isRedundant = false;
      if(isRedundant) {
        louti.remove();
        _phi[p2i][p1i] = 3;
      }
    }
  }

  /**
   * Build a triangle mesh from known points.
	 * @param fx a map with non-zero values only at know points
   */
  private void buildTmesh(float[][] fx) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      if(fx[i2][i1]>0f) {
        TriMesh.Node nd = new TriMesh.Node(i1,i2);
        _tmesh.addNode(nd);
      }
    }}
  }

  /**
   * Find nearest distance using the preconstructed triangle mesh
   */
  private float findDistance(float x1, float x2) {
    TriMesh.Node nd = _tmesh.findNodeNearest(x1,x2);
    float d1 = x1-nd.x();
    float d2 = x2-nd.y();
    return d1*d1+d2*d2;
  }



  /**
   * A point with 2D integer coordinates
   */
  public class Point {
	  public int p1; //1st coordinate
	  public int p2; //2nd coordinate
    public int fe; //evolution speed
    public int fs; //smooth speed
	  /**
	   * Constructor
	   * @param p1 The 1st coordinate
	   * @param p2 The 2nd coordinate
	   */
	  public Point(int p1, int p2) {
	  	this.p1 = p1;
	  	this.p2 = p2;
	  }

    public void setEvolutionSpeed(int fe) {
      this.fe = fe;
    }
    public void setSmoothSpeed(int fs) {
      this.fs = fs;
    }
  }
}
