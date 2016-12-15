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
import edu.mines.jtk.util.*;
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


public class FastLevelSet3 {
  /**
   * Construct with a 2D initial square-shape level-set function 
   * with integes -3,-1,1,3.
   * @param n1 the 1st dimension number.
   * @param n2 the 2nd dimension number.
   * @param c1 the x coordinate of the center.
   * @param c2 the y coordinate of the center.
   * @param r  the radius.
   */
  public FastLevelSet3(int n1, int n2, int n3, 
    int c1, int c2, int c3, int r) 
  {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _lins  = new LinkedList[n3];
    _louts = new LinkedList[n3];
    _linsAdd  = new LinkedList[n3];
    _loutsAdd = new LinkedList[n3];
    for (int i3=0; i3<n3; ++i3) {
      _lins[i3] = new LinkedList<Sample>();
      _louts[i3] = new LinkedList<Sample>();
      _linsAdd[i3] = new LinkedList<Sample>();
      _loutsAdd[i3] = new LinkedList<Sample>();
      _lins[i3].clear();
      _louts[i3].clear();
      _linsAdd[i3].clear();
      _loutsAdd[i3].clear();
    }
    _phi = fillbyte((byte)3,n1,n2,n3);
    int b1 = c1-r;
    int e1 = c1+r;
    int b2 = c2-r;
    int e2 = c2+r;
    int b3 = c3-r;
    int e3 = c3+r;
    b1 = max(b1,1);
    b2 = max(b2,1);
    b3 = max(b3,1);
    e1 = min(e1,_n1-2);
    e2 = min(e2,_n2-2);
    e3 = min(e3,_n3-2);
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      _phi[i3][i2][i1] = (byte)-3;
    }}}
    for (int i3=b3-1; i3<=e3+1; i3++) {
    for (int i2=b2-1; i2<=e2+1; i2++) {
      int i1b = b1-1;
      int i1e = e1+1;
      _phi[i3][i2][i1b] = (byte)1;
      _phi[i3][i2][i1e] = (byte)1;
      Point p1b = new Point(i1b,i2,i3);
      Point p1e = new Point(i1e,i2,i3);
      _lout.add(p1b);
      _lout.add(p1e);
      _louts[i3].add(new Sample(i1b,i2));
      _louts[i3].add(new Sample(i1e,i2));
    }}
    for (int i3=b3-1; i3<=e3+1; i3++) {
    for (int i1=b1-1; i1<=e1+1; i1++) {
      int i2b = b2-1;
      int i2e = e2+1;
      _phi[i3][i2b][i1] = (byte)1;
      _phi[i3][i2e][i1] = (byte)1;
      Point p2b = new Point(i1,i2b,i3);
      Point p2e = new Point(i1,i2e,i3);
      _lout.add(p2b);
      _lout.add(p2e);
      _louts[i3].add(new Sample(i1,i2b));
      _louts[i3].add(new Sample(i1,i2e));
    }}
    for (int i2=b2-1; i2<=e2+1; i2++) {
    for (int i1=b1-1; i1<=e1+1; i1++) {
      int i3b = b3-1;
      int i3e = e3+1;
      _phi[i3b][i2][i1] = (byte)1;
      _phi[i3e][i2][i1] = (byte)1;
      Point p3b = new Point(i1,i2,i3b);
      Point p3e = new Point(i1,i2,i3e);
      _lout.add(p3b);
      _lout.add(p3e);
      _louts[i3b].add(new Sample(i1,i2));
      _louts[i3e].add(new Sample(i1,i2));
    }}

    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
      _phi[i3][i2][b1] = (byte)-1;
      _phi[i3][i2][e1] = (byte)-1;
      Point p1b = new Point(b1,i2,i3);
      Point p1e = new Point(e1,i2,i3);
      _lin.add(p1b);
      _lin.add(p1e);
      _lins[i3].add(new Sample(b1,i2));
      _lins[i3].add(new Sample(e1,i2));
    }}
    for (int i3=b3; i3<=e3; i3++) {
    for (int i1=b1; i1<=e1; i1++) {
      _phi[i3][b2][i1] = (byte)-1;
      _phi[i3][e2][i1] = (byte)-1;
      Point p2b = new Point(i1,b2,i3);
      Point p2e = new Point(i1,e2,i3);
      _lin.add(p2b);
      _lin.add(p2e);
      _lins[i3].add(new Sample(i1,b2));
      _lins[i3].add(new Sample(i1,e2));
    }}
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      _phi[b3][i2][i1] = (byte)-1;
      _phi[b3][i2][i1] = (byte)-1;
      Point p3b = new Point(i1,i2,b3);
      Point p3e = new Point(i1,i2,e3);
      _lin.add(p3b);
      _lin.add(p3e);
      _lins[b3].add(new Sample(i1,i2));
      _lins[e3].add(new Sample(i1,i2));
    }}
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

  public void updateLevelSet(int gw, double sigma, float[][][] fx) {
    _gWindow = gw;
    createGaussFilter(gw,sigma);
    updateEvolutionSpeed(fx);
    for (int iter=0; iter<_outerIters; iter++) {
      System.out.println("iter="+iter);
      //System.out.println("nout="+_lout.size());
      //boolean converged = cycleOne(fx);
      //cycleTwo();
      cycleOneX(fx);
      cycleTwoX();
      //if(converged){break;}
    }
  }

  public float[][][] getPhi() {
    float[][][] phi = new float[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3)
    for (int i2=0; i2<_n2; ++i2)
    for (int i1=0; i1<_n1; ++i1)
      phi[i3][i2][i1] = _phi[i3][i2][i1];
    return phi;
  }

  public float[][] getLout() {
    int np = _lout.size();
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    float[] x3 = new float[np];
    for (int ip=0; ip<np; ++ip) {
      Point pi = _lout.get(ip);
      x1[ip] = pi.p1;
      x2[ip] = pi.p2;
      x3[ip] = pi.p3;
    }
    return new float[][]{x1,x2,x3};
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
  private int _n3;
  private byte[][][] _phi;
  private int _gWindow;
  private byte[][][] _gauss;
  private int _gaussThreshold;
  private int _outerIters = 100;
  private int _speedIters = 10;
  private int _smoothIters = 5;
  private List<Point> _lin  = new LinkedList<Point>(); //inside boundary points;
  private List<Point> _lout = new LinkedList<Point>(); //outside boundary points;
  private List<Point> _linAdd  = new LinkedList<Point>(); //inside boundary points;ko
  private List<Point> _loutAdd  = new LinkedList<Point>(); //inside boundary points;
  LinkedList<Sample>[] _lins;
  LinkedList<Sample>[] _louts;// = (ArrayList<Sample>[])new ArrayList[];
  LinkedList<Sample>[] _linsAdd;
  LinkedList<Sample>[] _loutsAdd;

  private void cycleOneX(final float[][][] fx) {
    for (int iter=0; iter<_smoothIters; iter++) {
      Parallel.loop(0,_n3,3,new Parallel.LoopInt() { // i1 = 0, 2, 4, ...
      public void compute(int i3) {
        switchInSlice3(i3,_louts[i3].iterator(),fx[i3]);
      }});
      Parallel.loop(1,_n3,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i3) {
        switchInSlice3(i3,_louts[i3].iterator(),fx[i3]);
      }});
      Parallel.loop(2,_n3,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i3) {
        switchInSlice3(i3,_louts[i3].iterator(),fx[i3]);
      }});
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        _louts[i3].addAll(0,_loutsAdd[i3]);
        _loutsAdd[i3].clear();
      }});
      //cleanLinX();
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        cleanLin(i3,_lins[i3].iterator());
      }});

      Parallel.loop(0,_n3,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i3) {
        switchOutSlice3(i3,_lins[i3].iterator(),fx[i3]);
      }});
      Parallel.loop(1,_n3,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        switchOutSlice3(i3,_lins[i3].iterator(),fx[i3]);
      }});
      Parallel.loop(2,_n3,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        switchOutSlice3(i3,_lins[i3].iterator(),fx[i3]);
      }});
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        _lins[i3].addAll(0,_linsAdd[i3]);
        _linsAdd[i3].clear();
      }});
      //cleanLoutX();
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        cleanLout(i3,_louts[i3].iterator());
      }});
    }
  }

  private void cycleTwoX() {
    for (int iter=0; iter<_smoothIters; iter++) {
      Parallel.loop(0,_n3,3,new Parallel.LoopInt() { // i1 = 0, 2, 4, ...
      public void compute(int i3) {
        switchInSlice3(i3,_louts[i3].iterator());
      }});
      Parallel.loop(1,_n3,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i3) {
        switchInSlice3(i3,_louts[i3].iterator());
      }});
      Parallel.loop(2,_n3,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i3) {
        switchInSlice3(i3,_louts[i3].iterator());
      }});
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        _louts[i3].addAll(0,_loutsAdd[i3]);
        _loutsAdd[i3].clear();
      }});

      //cleanLinX();
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        cleanLin(i3,_lins[i3].iterator());
      }});

      Parallel.loop(0,_n3,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i3) {
        switchOutSlice3(i3,_lins[i3].iterator());
      }});
      Parallel.loop(1,_n3,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        switchOutSlice3(i3,_lins[i3].iterator());
      }});
      Parallel.loop(2,_n3,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        switchOutSlice3(i3,_lins[i3].iterator());
      }});
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        _lins[i3].addAll(0,_linsAdd[i3]);
        _linsAdd[i3].clear();
      }});
      //cleanLoutX();
      Parallel.loop(_n3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i3) {
        cleanLout(i3,_louts[i3].iterator());
      }});
    }
  }

  private void cleanLin(int p3i, Iterator<Sample> lini) {
    while (lini.hasNext()) {
      Sample pi = lini.next();
      int p1i = pi.p1;
      int p2i = pi.p2;
      int p1m = p1i-1; 
      boolean isRedundant=true;
      //check the neighbor above
      if(p1m>=0&&_phi[p3i][p2i][p1m]>0)
        isRedundant = false;
      //check the neighbor below
      int p1p = p1i+1; 
      if(p1p<_n1&&_phi[p3i][p2i][p1p]>0)
        isRedundant = false;
      //check the left neighbor
      int p2m = p2i-1; 
      if(p2m>=0&&_phi[p3i][p2m][p1i]>0)
        isRedundant = false;
      //check the right neighbor
      int p2p = p2i+1; 
      if(p2p<_n2&&_phi[p3i][p2p][p1i]>0) 
        isRedundant = false;
      //check the left neighbor
      int p3m = p3i-1; 
      if(p3m>=0&&_phi[p3m][p2i][p1i]>0)
        isRedundant = false;
      //check the right neighbor
      int p3p = p3i+1; 
      if(p3p<_n3&&_phi[p3p][p2i][p1i]>0) 
        isRedundant = false;
      if(isRedundant) {
        lini.remove();
        _phi[p3i][p2i][p1i] = -3;
      }
    }
  }

  private void cleanLout(int p3i, Iterator<Sample> louti) {
    while (louti.hasNext()) {
      Sample pi = louti.next();
      int p1i = pi.p1;
      int p2i = pi.p2;
      int p1m = p1i-1; 
      boolean isRedundant=true;
      //check the neighbor above
      if(p1m>=0&&_phi[p3i][p2i][p1m]<0)
        isRedundant = false;
      //check the neighbor below
      int p1p = p1i+1; 
      if(p1p<_n1&&_phi[p3i][p2i][p1p]<0)
        isRedundant = false;
      //check the left neighbor
      int p2m = p2i-1; 
      if(p2m>=0&&_phi[p3i][p2m][p1i]<0)
        isRedundant = false;
      //check the right neighbor
      int p2p = p2i+1; 
      if(p2p<_n2&&_phi[p3i][p2p][p1i]<0) 
        isRedundant = false;
      //check the left neighbor
      int p3m = p3i-1; 
      if(p3m>=0&&_phi[p3m][p2i][p1i]<0)
        isRedundant = false;
      //check the right neighbor
      int p3p = p3i+1; 
      if(p3p<_n3&&_phi[p3p][p2i][p1i]<0) 
        isRedundant = false;
      if(isRedundant) {
        louti.remove();
        _phi[p3i][p2i][p1i] = 3;
      }
    }
  }


  private void switchInSlice3(int i3, Iterator<Sample>louti,float[][] fx3) {
    while(louti.hasNext()) {
      Sample sp = louti.next();
      int i1 = sp.p1;
      int i2 = sp.p2;
      if(fx3[i2][i1]<0.2f){switchIn(i3,louti,sp);}
    }
  }

  private void switchInSlice3(int i3, Iterator<Sample>louti) {
    while(louti.hasNext()) {
      Sample sp = louti.next();
      int i1 = sp.p1;
      int i2 = sp.p2;
      if(calculateSmoothSpeed(i1,i2,i3)>_gaussThreshold){switchIn(i3,louti,sp);}
    }
  }

  private void switchOutSlice3(int i3, Iterator<Sample>lini,float[][] fx3) {
    while(lini.hasNext()) {
      Sample sp = lini.next();
      int i1 = sp.p1;
      int i2 = sp.p2;
      if(fx3[i2][i1]>0.2f){switchOut(i3,lini,sp);}
    }
  }

  private void switchOutSlice3(int i3, Iterator<Sample>lini) {
    while(lini.hasNext()) {
      Sample sp = lini.next();
      int i1 = sp.p1;
      int i2 = sp.p2;
      if(calculateSmoothSpeed(i1,i2,i3)<_gaussThreshold){switchOut(i3,lini,sp);}
    }
  }


  private boolean cycleOne(float[][][] fx) {
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
      updateEvolutionSpeed(fx); 
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
		int m3 = 2*gw+1;
		int m2 = 2*gw+1;
		int m1 = 2*gw+1;
		int gfScale = 0;
    double sigmas = 1.0/(sigma*sigma);
    double sigmat = sigmas/sigma;
		// Rough heuristic: scale by number of elements in filter
		double scale1 = m1*m2*m3;
    _gauss = new byte[m3][m2][m1];
		// In theory could just calculate 1/8th and duplicate instead
		for (int i3=0; i3<m3; ++i3) {
      double d3  = i3-gw;
		for (int i2=0; i2<m2; ++i2) {
      double d2  = i2-gw;
			for (int i1=0; i1<m1; ++i1) {
        double d1  = i1-gw;
				double ds = d3*d3+d2*d2+d1*d1;
				double gf = sigmat*exp(-0.5*sigmas*ds)*scale1;
        _gauss[i3][i2][i1]=(byte)gf;
				gfScale += gf;
			}
		}}
		_gaussThreshold = gfScale/2;
	}

	/**
	 * Update evolution speeed and check convergence.
	 * @return true if convergence has been reached
	 */
  private boolean checkConvergence() {
      for (Point pi:_lout) {
        if(pi.fe>0) return false;
      }
      for (Point pi:_lin) {
        if(pi.fe<0) return false;
      }
      return true;
  }

  private void updateEvolutionSpeed(float[][][] fx) {
    for (Point pi:_lin) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      int p3i = pi.p3;
      if (fx[p3i][p2i][p1i]>1.2f)
        pi.setEvolutionSpeed(-1);
      else
        pi.setEvolutionSpeed(0);
    }
    for (Point pi:_lout) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      int p3i = pi.p3;
      if (fx[p3i][p2i][p1i]<0.2f)
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
    int p1 = p.p1;
    int p2 = p.p2;
    int p3 = p.p3;
    int gw = _gWindow;
		int d1m = max(-gw,-p1);
		int d2m = max(-gw,-p2);
		int d3m = max(-gw,-p3);
		int d1p = min(gw+1,_n1-p1);
		int d2p = min(gw+1,_n2-p2);
		int d3p = min(gw+1,_n3-p3);
		int f = 0;
		for(int d3=d3m; d3<d3p; ++d3) {
		for(int d2=d2m; d2<d2p; ++d2) {
		for(int d1=d1m; d1<d1p; ++d1) {
			if (_phi[p3+d3][p2+d2][p1+d1]<0) 
			  f += _gauss[gw+d3][gw+d2][gw+d1];
		}}}
		p.setSmoothSpeed(f);
	}

	private float calculateSmoothSpeed(int p1, int p2, int p3) {
		// Convolve neighbourhood of a point with a gaussian
    int gw = _gWindow;
		int d1m = max(-gw,-p1);
		int d2m = max(-gw,-p2);
		int d3m = max(-gw,-p3);
		int d1p = min(gw+1,_n1-p1);
		int d2p = min(gw+1,_n2-p2);
		int d3p = min(gw+1,_n3-p3);
		int f = 0;
		for(int d3=d3m; d3<d3p; ++d3) {
		for(int d2=d2m; d2<d2p; ++d2) {
		for(int d1=d1m; d1<d1p; ++d1) {
			if (_phi[p3+d3][p2+d2][p1+d1]<0) 
			  f += _gauss[gw+d3][gw+d2][gw+d1];
		}}}
    return f;
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
    int p3i = p.p3;
    _lin.add(p);
    _phi[p3i][p2i][p1i]=(byte)-1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p3i][p2i][p1m]==3) {
      Point pa = new Point(p1m,p2i,p3i);
      _loutAdd.add(pa);
      _phi[p3i][p2i][p1m] = 1;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p3i][p2i][p1p]==3) {
      Point pb = new Point(p1p,p2i,p3i);
      _loutAdd.add(pb);
      _phi[p3i][p2i][p1p] = 1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p3i][p2m][p1i]==3) {
      Point pl = new Point(p1i,p2m,p3i);
      _loutAdd.add(pl);
      _phi[p3i][p2m][p1i] = 1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p3i][p2p][p1i]==3) {
      Point pr = new Point(p1i,p2p,p3i);
      _loutAdd.add(pr);
      _phi[p3i][p2p][p1i] = 1;
    }

    // check the left neighbor
    int p3m = p3i-1; 
    if(p3m>=0&&_phi[p3m][p2i][p1i]==3) {
      Point pl = new Point(p1i,p2i,p3m);
      _loutAdd.add(pl);
      _phi[p3m][p2i][p1i] = 1;
    }
     // check the right neighbor
    int p3p = p3i+1; 
    if(p3p<_n3&&_phi[p3p][p2i][p1i]==3) {
      Point pr = new Point(p1i,p2i,p3p);
      _loutAdd.add(pr);
      _phi[p3p][p2i][p1i] = 1;
    }
    louti.remove();
	}

	private void switchIn(int p3i, Iterator<Sample> louti, Sample p) {
    //step 1: delete the point from lout and add it to lin; set phi=-1
    int p1i = p.p1;
    int p2i = p.p2;
    _lins[p3i].add(p);
    _phi[p3i][p2i][p1i]=(byte)-1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p3i][p2i][p1m]==3) {
      Sample pa = new Sample(p1m,p2i);
      _loutsAdd[p3i].add(pa);
      _phi[p3i][p2i][p1m] = 1;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p3i][p2i][p1p]==3) {
      Sample pb = new Sample(p1p,p2i);
      _loutsAdd[p3i].add(pb);
      _phi[p3i][p2i][p1p] = 1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p3i][p2m][p1i]==3) {
      Sample pl = new Sample(p1i,p2m);
      _loutsAdd[p3i].add(pl);
      _phi[p3i][p2m][p1i] = 1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p3i][p2p][p1i]==3) {
      Sample pr = new Sample(p1i,p2p);
      _loutsAdd[p3i].add(pr);
      _phi[p3i][p2p][p1i] = 1;
    }

    // check the left neighbor
    int p3m = p3i-1; 
    if(p3m>=0&&_phi[p3m][p2i][p1i]==3) {
      Sample pl = new Sample(p1i,p2i);
      _loutsAdd[p3m].add(pl);
      _phi[p3m][p2i][p1i] = 1;
    }
     // check the right neighbor
    int p3p = p3i+1; 
    if(p3p<_n3&&_phi[p3p][p2i][p1i]==3) {
      Sample pr = new Sample(p1i,p2i);
      _loutsAdd[p3p].add(pr);
      _phi[p3p][p2i][p1i] = 1;
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
    int p3i = p.p3;
    _lout.add(p);
    _phi[p3i][p2i][p1i]= 1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p3i][p2i][p1m]==-3) {
      Point pa = new Point(p1m,p2i,p3i);
      _linAdd.add(pa);
      _phi[p3i][p2i][p1m] = -1;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p3i][p2i][p1p]==-3) {
      Point pb = new Point(p1p,p2i,p3i);
      _linAdd.add(pb);
      _phi[p3i][p2i][p1p] = -1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p3i][p2m][p1i]==-3) {
      Point pl = new Point(p1i,p2m,p3i);
      _linAdd.add(pl);
      _phi[p3i][p2m][p1i] = -1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p3i][p2p][p1i]==-3) {
      Point pr = new Point(p1i,p2p,p3i);
      _linAdd.add(pr);
      _phi[p3i][p2p][p1i] = -1;
    }
    // check the left neighbor
    int p3m = p3i-1; 
    if(p3m>=0&&_phi[p3m][p2i][p1i]==-3) {
      Point pl = new Point(p1i,p2i,p3m);
      _linAdd.add(pl);
      _phi[p3m][p2i][p1i] = -1;
    }
     // check the right neighbor
    int p3p = p3i+1; 
    if(p3p<_n3&&_phi[p3p][p2i][p1i]==-3) {
      Point pr = new Point(p1i,p2i,p3p);
      _linAdd.add(pr);
      _phi[p3p][p2i][p1i] = -1;
    }
    lini.remove();
	}


	private void switchOut(int p3i, Iterator<Sample> lini, Sample p) {
    //step 1: delete the point from lout and add it to lin; set phi=-1
    int p1i = p.p1;
    int p2i = p.p2;
    _louts[p3i].add(p);
    _phi[p3i][p2i][p1i]= 1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p3i][p2i][p1m]==-3) {
      Sample pa = new Sample(p1m,p2i);
      _linsAdd[p3i].add(pa);
      _phi[p3i][p2i][p1m] = -1;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p3i][p2i][p1p]==-3) {
      Sample pb = new Sample(p1p,p2i);
      _linsAdd[p3i].add(pb);
      _phi[p3i][p2i][p1p] = -1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p3i][p2m][p1i]==-3) {
      Sample pl = new Sample(p1i,p2m);
      _linsAdd[p3i].add(pl);
      _phi[p3i][p2m][p1i] = -1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p3i][p2p][p1i]==-3) {
      Sample pr = new Sample(p1i,p2p);
      _linsAdd[p3i].add(pr);
      _phi[p3i][p2p][p1i] = -1;
    }
    // check the left neighbor
    int p3m = p3i-1; 
    if(p3m>=0&&_phi[p3m][p2i][p1i]==-3) {
      Sample pl = new Sample(p1i,p2i);
      _phi[p3m][p2i][p1i] = -1;
      _linsAdd[p3m].add(pl);
    }
     // check the right neighbor
    int p3p = p3i+1; 
    if(p3p<_n3&&_phi[p3p][p2i][p1i]==-3) {
      Sample pr = new Sample(p1i,p2i);
      _phi[p3p][p2i][p1i] = -1;
      _linsAdd[p3p].add(pr);
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
      int p3i = pi.p3;
      int p1m = p1i-1; 
      boolean isRedundant=true;
      //check the neighbor above
      if(p1m>=0&&_phi[p3i][p2i][p1m]>0)
        isRedundant = false;
      //check the neighbor below
      int p1p = p1i+1; 
      if(p1p<_n1&&_phi[p3i][p2i][p1p]>0)
        isRedundant = false;
      //check the left neighbor
      int p2m = p2i-1; 
      if(p2m>=0&&_phi[p3i][p2m][p1i]>0)
        isRedundant = false;
      //check the right neighbor
      int p2p = p2i+1; 
      if(p2p<_n2&&_phi[p3i][p2p][p1i]>0) 
        isRedundant = false;
      //check the left neighbor
      int p3m = p3i-1; 
      if(p3m>=0&&_phi[p3m][p2i][p1i]>0)
        isRedundant = false;
      //check the right neighbor
      int p3p = p3i+1; 
      if(p3p<_n3&&_phi[p3p][p2i][p1i]>0) 
        isRedundant = false;
      if(isRedundant) {
        lini.remove();
        _phi[p3i][p2i][p1i] = -3;
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
      int p3i = pi.p3;
      int p1m = p1i-1; 
      boolean isRedundant=true;
      //check the neighbor above
      if(p1m>=0&&_phi[p3i][p2i][p1m]<0)
        isRedundant = false;
      //check the neighbor below
      int p1p = p1i+1; 
      if(p1p<_n1&&_phi[p3i][p2i][p1p]<0)
        isRedundant = false;
      //check the left neighbor
      int p2m = p2i-1; 
      if(p2m>=0&&_phi[p3i][p2m][p1i]<0)
        isRedundant = false;
      //check the right neighbor
      int p2p = p2i+1; 
      if(p2p<_n2&&_phi[p3i][p2p][p1i]<0) 
        isRedundant = false;
      //check the left neighbor
      int p3m = p3i-1; 
      if(p3m>=0&&_phi[p3m][p2i][p1i]<0)
        isRedundant = false;
      //check the right neighbor
      int p3p = p3i+1; 
      if(p3p<_n3&&_phi[p3p][p2i][p1i]<0) 
        isRedundant = false;
      if(isRedundant) {
        louti.remove();
        _phi[p3i][p2i][p1i] = 3;
      }
    }
  }


  /**
   * A point with 2D integer coordinates
   */
  public class Point {
	  public int p1; //1st coordinate
	  public int p2; //2nd coordinate
	  public int p3; //2nd coordinate
    public int fe; //evolution speed
    public int fs; //smooth speed
	  /**
	   * Constructor
	   * @param p1 The 1st coordinate
	   * @param p2 The 2nd coordinate
	   */
	  public Point(int p1, int p2, int p3) {
	  	this.p1 = p1;
	  	this.p2 = p2;
	  	this.p3 = p3;
	  }

    public void setEvolutionSpeed(int fe) {
      this.fe = fe;
    }
    public void setSmoothSpeed(int fs) {
      this.fs = fs;
    }
  }

  public class Sample {
	  public int p1; //1st coordinate
	  public int p2; //2nd coordinate
	  /**
	   * Constructor
	   * @param p1 The 1st coordinate
	   * @param p2 The 2nd coordinate
	   */
	  public Sample(int p1, int p2) {
	  	this.p1 = p1;
	  	this.p2 = p2;
	  }

  }

}
