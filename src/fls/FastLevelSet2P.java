/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.Iterator;
import java.util.LinkedList;
import edu.mines.jtk.dsp.*;
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


public class FastLevelSet2P {
  /**
   * Construct with a 2D initial square-shape level-set function 
   * with integes -3,-1,1,3.
   * @param n1 the 1st dimension number.
   * @param n2 the 2nd dimension number.
   * @param c1 the x coordinate of the center.
   * @param c2 the y coordinate of the center.
   * @param r  the radius.
   */
  public FastLevelSet2P(int n1, int n2,int c1, int c2, int r) {
    _n1 = n1;
    _n2 = n2;
    _lins  = new LinkedList[n2];
    _louts = new LinkedList[n2];
    _linsAdd  = new LinkedList[n2];
    _loutsAdd = new LinkedList[n2];
    for (int i2=0; i2<n2; ++i2) {
      _lins[i2] = new LinkedList<Sample>();
      _louts[i2] = new LinkedList<Sample>();
      _linsAdd[i2] = new LinkedList<Sample>();
      _loutsAdd[i2] = new LinkedList<Sample>();
      _lins[i2].clear();
      _louts[i2].clear();
      _linsAdd[i2].clear();
      _loutsAdd[i2].clear();
    }
    _phi = fillbyte((byte)3,n1,n2);
    initialize(c1,c2,r);
  }

  public FastLevelSet2P(int n1, int n2,int[] c1, int[] c2, int[] r) {
    _n1 = n1;
    _n2 = n2;
    _lins  = new LinkedList[n2];
    _louts = new LinkedList[n2];
    _linsAdd  = new LinkedList[n2];
    _loutsAdd = new LinkedList[n2];
    for (int i2=0; i2<n2; ++i2) {
      _lins[i2] = new LinkedList<Sample>();
      _louts[i2] = new LinkedList<Sample>();
      _linsAdd[i2] = new LinkedList<Sample>();
      _loutsAdd[i2] = new LinkedList<Sample>();
      _lins[i2].clear();
      _louts[i2].clear();
      _linsAdd[i2].clear();
      _loutsAdd[i2].clear();
    }
    int np = c1.length;
    _phi = fillbyte((byte)3,n1,n2);
    for (int ip=0; ip<np; ++ip)
      initialize(c1[ip],c2[ip],r[ip]);
  }


  private void initialize(int c1, int c2, int r) {
    int b1 = c1-r;
    int e1 = c1+r;
    int b2 = c2-r;
    int e2 = c2+r;
    b1 = max(b1,1);
    b2 = max(b2,1);
    e1 = min(e1,_n1-2);
    e2 = min(e2,_n2-2);
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      _phi[i2][i1] = (byte)-3;
    }}
    for (int i2=b2-1; i2<=e2+1; i2++) {
      int i1b = b1-1;
      int i1e = e1+1;
      _phi[i2][i1b] = (byte)1;
      _phi[i2][i1e] = (byte)1;
      _louts[i2].add(new Sample(i1b));
      _louts[i2].add(new Sample(i1e));
    }
    for (int i1=b1-1; i1<=e1+1; i1++) {
      int i2b = b2-1;
      int i2e = e2+1;
      _phi[i2b][i1] = (byte)1;
      _phi[i2e][i1] = (byte)1;
      _louts[i2b].add(new Sample(i1));
      _louts[i2e].add(new Sample(i1));
    }

    for (int i2=b2; i2<=e2; i2++) {
      _phi[i2][b1] = (byte)-1;
      _phi[i2][e1] = (byte)-1;
      _lins[i2].add(new Sample(b1));
      _lins[i2].add(new Sample(e1));
    }

    for (int i1=b1; i1<=e1; i1++) {
      _phi[b2][i1] = (byte)-1;
      _phi[e2][i1] = (byte)-1;
      _lins[b2].add(new Sample(i1));
      _lins[e2].add(new Sample(i1));
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


  public float[][] applyForDensity(float ev, float[][] el, float[][][] gs) {
    int n3 = gs.length;
    int n2 = gs[0].length;
    int n1 = gs[0][0].length;
    float[][][] dps = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3)
      dps[i3] = density(ev,el,gs[i3]);
    balanceDensity(dps);
    float[][] dp = new float[n2][n1];
    for (int i3=0; i3<n3; ++i3)
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      dp[i2][i1] += dps[i3][i2][i1];
    return dp;
  }

  public float[][][] applySegments(
    int gw, float sigma, float[][] dp, float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    float[][] fc = copy(ft);
    zero(ft);
    add(ft,4,ft);
    float[][][] xs = new float[2][2][];
    xs[0] = getLout();
    updateLevelSet(gw,sigma,dp);
    xs[1] = getLout();
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(_phi[i2][i1]<=1) {ft[i2][i1] = fc[i2][i1];}
    }}
    return xs;
  }

  public float[][] getLout() {
    int k = 0;
    float[] x1 = new float[_n2*_n1];
    float[] x2 = new float[_n2*_n1];
    for (int i2=10; i2<_n2-10; ++i2) {
      for (Sample sp:_louts[i2]) {
        int i1 = sp.p1;
        if(i1<10||i1>=_n1-10){continue;}
        x1[k] = i1;
        x2[k] = i2;
        k++;
      }
    }
    return new float[][]{copy(k,0,x1),copy(k,0,x2)};
  }

  public void updateLevelSet(int gw, double sigma, float[][] dp) {
    _gWindow = gw;
    createGaussFilter(gw,sigma);
    for (int iter=0; iter<_outerIters; iter++) {
      if(iter%50==0) System.out.println("iter="+iter);
      cycleOne(dp);
      cycleTwo();
      //if(converged){break;}
    }
  }

  public float[][] getPhi() {
    float[][] phi = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2)
    for (int i1=0; i1<_n1; ++i1)
      phi[i2][i1] = _phi[i2][i1];
    return phi;
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

  private void balanceDensity(float[][][] ds) {
    int n3 = ds.length;
    int n2 = ds[0].length;
    int n1 = ds[0][0].length;
    float[] sc = new float[n3];
    float[] sd = new float[n3];
    for (int i3=0; i3<n3; ++i3) {
      float[][] d3 = ds[i3];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if(d3[i2][i1]<0f) {
          sc[i3] += 1f;
          sd[i3] += abs(d3[i2][i1]);
        }
      }}
    }
    div(sd,sc,sd);
    float sdmax = max(sd);
    for (int i3=0; i3<n3; ++i3) {
      float sci = sdmax/sd[i3];
      mul(ds[i3],sci,ds[i3]);
    }

  }


  public float[][] density(
    float f, float[][] fx, float[][] gx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[] p1 = new float[256];
    float[] p2 = new float[256];
    float[][] dp = new float[n2][n1];
    int[][] gg = toGrayIntegers(gx);
    density(f,fx,gg,p1,p2);
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      int ggi = gg[i2][i1];
      float di = log(p2[ggi])-log(p1[ggi]);
      if(Float.isNaN(di)) {continue;}
      if(Float.isInfinite(di)) {continue;}
      dp[i2][i1] = di;
    }}
    return dp;
  }

  public void density(float f, float[][] fx, 
    int[][] gx, float[] p1, float[] p2) {
    float c1 = 0f;
    float c2 = 0f;
    int n2 = fx.length;
    int n1 = fx[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int gxi = gx[i2][i1];
      float fxi = fx[i2][i1];
      if(fxi<f) {
        c1 += 1f;
        p1[gxi] += 1f;
      } else {
        c2 += 1f;
        p2[gxi] += 1f;
      }
    }}
    div(p1,c1,p1);
    div(p2,c2,p2);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2);
    rgf.apply0(p1,p1);
    rgf.apply0(p2,p2);
  }

  private int[][] toGrayIntegers(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int[][] gx = new int[n2][n1];
    fx = sub(fx,min(fx));
    float fmax = 255f/max(fx);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i2][i1] = round(fx[i2][i1]*fmax);
    }}
    return gx;
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
  LinkedList<Sample>[] _lins;
  LinkedList<Sample>[] _louts;// = (ArrayList<Sample>[])new ArrayList[];
  LinkedList<Sample>[] _linsAdd;
  LinkedList<Sample>[] _loutsAdd;

  private void cycleOne(final float[][] fx) {
    for (int iter=0; iter<_speedIters; iter++) {
      Parallel.loop(0,_n2,3,new Parallel.LoopInt() { // i1 = 0, 2, 4, ...
      public void compute(int i2) {
        switchInSlice2(i2,_louts[i2].iterator(),fx[i2]);
      }});
      Parallel.loop(1,_n2,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i2) {
        switchInSlice2(i2,_louts[i2].iterator(),fx[i2]);
      }});
      Parallel.loop(2,_n2,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i2) {
        switchInSlice2(i2,_louts[i2].iterator(),fx[i2]);
      }});
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        _louts[i2].addAll(0,_loutsAdd[i2]);
        _loutsAdd[i2].clear();
      }});

      //cleanLinX();
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        cleanLin(i2,_lins[i2].iterator());
      }});

      Parallel.loop(0,_n2,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i2) {
        switchOutSlice2(i2,_lins[i2].iterator(),fx[i2]);
      }});
      Parallel.loop(1,_n2,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        switchOutSlice2(i2,_lins[i2].iterator(),fx[i2]);
      }});
      Parallel.loop(2,_n2,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        switchOutSlice2(i2,_lins[i2].iterator(),fx[i2]);
      }});
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        _lins[i2].addAll(0,_linsAdd[i2]);
        _linsAdd[i2].clear();
      }});
      //cleanLoutX();
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        cleanLout(i2,_louts[i2].iterator());
      }});
    }
  }

  private void cycleTwo() {
    for (int iter=0; iter<_smoothIters; iter++) {
      Parallel.loop(0,_n2,3,new Parallel.LoopInt() { // i1 = 0, 2, 4, ...
      public void compute(int i2) {
        switchInSlice2(i2,_louts[i2].iterator());
      }});
      Parallel.loop(1,_n2,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i2) {
        switchInSlice2(i2,_louts[i2].iterator());
      }});
      Parallel.loop(2,_n2,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i2) {
        switchInSlice2(i2,_louts[i2].iterator());
      }});
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        _louts[i2].addAll(0,_loutsAdd[i2]);
        _loutsAdd[i2].clear();
      }});

      //cleanLinX();
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        cleanLin(i2,_lins[i2].iterator());
      }});

      Parallel.loop(0,_n2,3,new Parallel.LoopInt() { // i1 = 1, 3, 5, ...
      public void compute(int i2) {
        switchOutSlice2(i2,_lins[i2].iterator());
      }});
      Parallel.loop(1,_n2,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        switchOutSlice2(i2,_lins[i2].iterator());
      }});
      Parallel.loop(2,_n2,3,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        switchOutSlice2(i2,_lins[i2].iterator());
      }});
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        _lins[i2].addAll(0,_linsAdd[i2]);
        _linsAdd[i2].clear();
      }});
      //cleanLoutX();
      Parallel.loop(_n2,new Parallel.LoopInt() { // i1 = 2, 4, 6, ...
      public void compute(int i2) {
        cleanLout(i2,_louts[i2].iterator());
      }});
    }
  }

  private void cleanLin(int p2i, Iterator<Sample> lini) {
    while (lini.hasNext()) {
      Sample pi = lini.next();
      int p1i = pi.p1;
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

  private void cleanLout(int p2i, Iterator<Sample> louti) {
    while (louti.hasNext()) {
      Sample pi = louti.next();
      int p1i = pi.p1;
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


  private void switchInSlice2(int i2, Iterator<Sample>louti,float[] fx2) {
    while(louti.hasNext()) {
      Sample sp = louti.next();
      int i1 = sp.p1;
      if(fx2[i1]<0){switchIn(i2,louti,sp);}
    }
  }

  private void switchInSlice2(int i2, Iterator<Sample>louti) {
    while(louti.hasNext()) {
      Sample sp = louti.next();
      int i1 = sp.p1;
      if(smoothSpeed(i1,i2)>_gaussThreshold){switchIn(i2,louti,sp);}
    }
  }

  private void switchOutSlice2(int i2, Iterator<Sample>lini,float[] fx2) {
    while(lini.hasNext()) {
      Sample sp = lini.next();
      int i1 = sp.p1;
      if(fx2[i1]>0){switchOut(i2,lini,sp);}
    }
  }

  private void switchOutSlice2(int i2, Iterator<Sample>lini) {
    while(lini.hasNext()) {
      Sample sp = lini.next();
      int i1 = sp.p1;
      if(smoothSpeed(i1,i2)<_gaussThreshold){switchOut(i2,lini,sp);}
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
  /*
  private boolean checkConvergence() {
      for (Point pi:_lout) {
        if(pi.fe>0) return false;
      }
      for (Point pi:_lin) {
        if(pi.fe<0) return false;
      }
      return true;
  }
  */


  /**
	 * Calculate the smoothing field at a point
	 */
	private float smoothSpeed(int p1, int p2) {
		// Convolve neighbourhood of a point with a gaussian
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
    return f;
	}

	private void switchIn(int p2i, Iterator<Sample> louti, Sample p) {
    //step 1: delete the point from lout and add it to lin; set phi=-1
    int p1i = p.p1;
    p.ct=0;
    _lins[p2i].add(p);
    _phi[p2i][p1i]=(byte)-1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p2i][p1m]==3) {
      Sample pa = new Sample(p1m);
      _loutsAdd[p2i].add(pa);
      _phi[p2i][p1m] = 1;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p2i][p1p]==3) {
      Sample pb = new Sample(p1p);
      _loutsAdd[p2i].add(pb);
      _phi[p2i][p1p] = 1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p2m][p1i]==3) {
      Sample pl = new Sample(p1i);
      _loutsAdd[p2m].add(pl);
      _phi[p2m][p1i] = 1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p2p][p1i]==3) {
      Sample pr = new Sample(p1i);
      _loutsAdd[p2p].add(pr);
      _phi[p2p][p1i] = 1;
    }

    louti.remove();
	}


	private void switchOut(int p2i, Iterator<Sample> lini, Sample p) {
    //step 1: delete the point from lout and add it to lin; set phi=-1
    int p1i = p.p1;
    p.ct=0;
    _louts[p2i].add(p);
    _phi[p2i][p1i]= 1;
    //step 2: check neighbors
    // check the neighbor above
    int p1m = p1i-1; 
    if(p1m>=0&&_phi[p2i][p1m]==-3) {
      Sample pa = new Sample(p1m);
      _linsAdd[p2i].add(pa);
      _phi[p2i][p1m] = -1;
    }

    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p2i][p1p]==-3) {
      Sample pb = new Sample(p1p);
      _linsAdd[p2i].add(pb);
      _phi[p2i][p1p] = -1;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p2m][p1i]==-3) {
      Sample pl = new Sample(p1i);
      _linsAdd[p2m].add(pl);
      _phi[p2m][p1i] = -1;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p2p][p1i]==-3) {
      Sample pr = new Sample(p1i);
      _linsAdd[p2p].add(pr);
      _phi[p2p][p1i] = -1;
    }

    lini.remove();
	}

  public class Sample {
	  public int p1; //1st coordinate
    public int ct=0;
	  /**
	   * Constructor
	   * @param p1 The 1st coordinate
	   */
	  public Sample(int p1) {
	  	this.p1 = p1;
	  }

  }

}
