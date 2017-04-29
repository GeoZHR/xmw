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
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mesh.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import pik.*;
import util.*;

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
    int i1b = b1-1;
    int i1e = e1+1;
    for (int i2=b2-1; i2<=e2+1; i2++) {
      _phi[i2][i1b] = (byte)1;
      _phi[i2][i1e] = (byte)1;
      Point p1b = new Point(i1b,i2);
      Point p1e = new Point(i1e,i2);
      _lout.add(p1b);
      _lout.add(p1e);
    }
    int i2b = b2-1;
    int i2e = e2+1;
    for (int i1=b1-1; i1<=e1+1; i1++) {
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
    if(_sin==0f||_sout==0f) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        byte phi = _phi[i2][i1];
        if (phi==-3) _sin += 1f;
        if (phi== 3) _sout += 1f;
      }}
    }

  }

  /**
   * Construct a 2D level set function using a known phi map
   * with integers -3,-1,1,3.
   * @param phi  a map with integers -3, -1, 1, 3.
   */
  public FastLevelSet2(float[][] ph) {
    _lin.clear();
    _lout.clear();
    _linAdd.clear();
    _loutAdd.clear();
    _n2 = ph.length;
    _n1 = ph[0].length;
    _phi = fillbyte((byte)3,_n1,_n2);
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      float phi = ph[i2][i1];
      if(phi==-1f){_lin.add(new Point(i1,i2));}
      if(phi==1f){_lout.add(new Point(i1,i2));}
      _phi[i2][i1] = (byte)phi;
    }}
    if(_sin==0f||_sout==0f) {
      for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
        byte phi = _phi[i2][i1];
        if (phi==-3) _sin += 1f;
        if (phi== 3) _sout += 1f;
      }}
    }

  }

  public float[][] getPhi() {
    float[][] phi = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2)
    for (int i1=0; i1<_n1; ++i1)
      phi[i2][i1] = _phi[i2][i1];
    return phi;

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

  public void updateLevelSet(int gw, double sigma, 
    float[][] el, float[][] fx, float[][] p2) {
    _gWindow = gw;
    createGaussFilter(gw,sigma);
    float[][] damp = density(0.3f,el,fx);
    float[][] ddip = density(0.3f,el,p2);
    //float[][][] dps = new float[][][]{damp,ddip};
    float[][][] dps = new float[][][]{damp};
    balanceDensity(dps);
    updateEvolutionSpeed(dps); 
    for (int iter=0; iter<_outerIters; iter++) {
      if(iter%50==0) System.out.println("iter="+iter);
      boolean converged = cycleOne(dps);
      cycleTwo();
      if(converged){break;}
    }
  }

  public float[][][][] refine(int r, float d, float[][] ph, float[][] fx) {
    int n2 = ph.length;
    int n1 = ph[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(8.0);
    rgf1.apply00(ph,ph);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply1X(ph,g1);
    rgf.applyX1(ph,g2);
    ContourFinder cf = new ContourFinder();
    ContourFinder.Contour ct = cf.findContour(0,s1,s2,ph);
    float[][][] ps = ct.getCoords();
    int nc = ps.length;
    float[][][] bs = bandSample(r,d,ps,ph,g1,g2,fx);
    SincInterpolator si = new SincInterpolator();
    for (int ic=0; ic<nc; ++ic) {
      int m2 = bs[ic].length;
      int m1 = bs[ic][0].length;
      sub(bs[ic],min(bs[ic]),bs[ic]);
      div(bs[ic],max(bs[ic]),bs[ic]);
      OptimalPathPicker opp = new OptimalPathPicker(20,1.0f);
      float[][] ft = opp.applyTransform(bs[ic]);
      float[][] wht = opp.applyForWeight(ft);
      float[][] tms1 = zerofloat(m2,m1);
      float[][] tms2 = zerofloat(m2,m1);
      float[] pik1 = opp.forwardPick(r,wht,tms1);
      float[] pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2);
      int np = ps[ic][0].length;
      float[] x1s = ps[ic][0];
      float[] x2s = ps[ic][1];
      for (int ip=0; ip<np; ++ip) {
        float x1i = x1s[ip];
        float x2i = x2s[ip];
        float g1i = si.interpolate(s1,s2,g1,x1i,x2i);
        float g2i = si.interpolate(s1,s2,g2,x1i,x2i);
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        g1i /= gsi;
        g2i /= gsi;
        x1s[ip] = x1i+g1i*(pik2[ip]-r)*d;
        x2s[ip] = x2i+g2i*(pik2[ip]-r)*d;
      }
    }

    /*
    bs = bandSample(r,d,ps,ph,g1,g2,fx);
    for (int ic=0; ic<nc; ++ic) {
      int m2 = bs[ic].length;
      int m1 = bs[ic][0].length;
      sub(bs[ic],min(bs[ic]),bs[ic]);
      div(bs[ic],max(bs[ic]),bs[ic]);
      OptimalPathPicker opp = new OptimalPathPicker(20,10f);
      float[][] ft = opp.applyTransform(bs[ic]);
      float[][] wht = opp.applyForWeight(ft);
      float[][] tms1 = zerofloat(m2,m1);
      float[][] tms2 = zerofloat(m2,m1);
      float[] pik1 = opp.forwardPick(r,wht,tms1);
      float[] pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2);
      int np = ps[ic][0].length;
      float[] x1s = ps[ic][0];
      float[] x2s = ps[ic][1];
      for (int ip=0; ip<np; ++ip) {
        float x1i = x1s[ip];
        float x2i = x2s[ip];
        float g1i = si.interpolate(s1,s2,g1,x1i,x2i);
        float g2i = si.interpolate(s1,s2,g2,x1i,x2i);
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        g1i /= gsi;
        g2i /= gsi;
        x1s[ip] = x1i+g1i*(pik2[ip]-r)*d;
        x2s[ip] = x2i+g2i*(pik2[ip]-r)*d;
      }
    }
    */
    return new float[][][][]{ps,bs};
  }

    public float[][][] bandSample(
    int r, float d, float[][][] ps, float[][] ph, 
    float[][] g1, float[][] g2, float[][] fx) {
    int n2 = ph.length;
    int n1 = ph[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    int nc = ps.length;
    float[][][] fbs = new float[nc][][];
    SincInterpolator si = new SincInterpolator();
    for (int ic=0; ic<nc; ++ic) {
      float[] x1s = ps[ic][0];
      float[] x2s = ps[ic][1];
      int np = x1s.length;
      fbs[ic] = new float[np][2*r+1];
      for (int ip=0; ip<np; ++ip) {
        float x1i = x1s[ip];
        float x2i = x2s[ip];
        float g1i = si.interpolate(s1,s2,g1,x1i,x2i);
        float g2i = si.interpolate(s1,s2,g2,x1i,x2i);
        float gsi = sqrt(g1i*g1i+g2i*g2i);
        g1i /= gsi;
        g2i /= gsi;
        fbs[ic][ip][r] = si.interpolate(s1,s2,fx,x1i,x2i);
        for (int ir=1; ir<=r; ++ir) {
          float y1i = x1i+g1i*ir*d;
          float y2i = x2i+g2i*ir*d;
          float phi = si.interpolate(s1,s2,ph,y1i,y2i);
          float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
          fbs[ic][ip][ir+r]=fxi;
          //if (phi>0f){fbs[ic][ip][ir+r]=fxi;}
          //else {fbs[ic][ip][ir+r]=-10f;}
        }
        for (int ir=-r; ir<=-1; ++ir) {
          float y1i = x1i+g1i*ir*d;
          float y2i = x2i+g2i*ir*d;
          float phi = si.interpolate(s1,s2,ph,y1i,y2i);
          float fxi = si.interpolate(s1,s2,fx,y1i,y2i);
          fbs[ic][ip][ir+r]=fxi;
          //if (phi<0f){fbs[ic][ip][ir+r]=fxi;}
          //else {fbs[ic][ip][ir+r]=-10f;}
        }
      }
    }
    return fbs;
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




  public void updateLevelSet(int gw, double sigma, float[][] fx) {
    _gWindow = gw;
    createGaussFilter(gw,sigma);
    float[][] g1 = new float[_n2][_n1];
    float[][] g2 = new float[_n2][_n1];
    LocalOrientFilter lof = new LocalOrientFilter(1,1);
    float[][] u1 = zerofloat(_n1,_n2);
    float[][] u2 = zerofloat(_n1,_n2);
    float[][] el = zerofloat(_n1,_n2);
    lof.applyForNormalLinear(fx,u1,u2,el);
    LocalSlopeFinder lsf = new LocalSlopeFinder(1,1,5);
    lsf.findSlopes(fx,u1);

    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(1.0);
    RecursiveGaussianFilter rgf2 = new RecursiveGaussianFilter(1.0);
    rgf1.apply10(fx,g1);
    rgf1.apply01(fx,g2);
    float[][] g11 = new float[_n2][_n1];
    float[][] g12 = new float[_n2][_n1];
    float[][] g22 = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      float eli = el[i2][i1];
      g11[i2][i1] = g1i*g1i;
      g12[i2][i1] = g1i*g2i;
      g22[i2][i1] = g2i*g2i;
      if(Float.isNaN(eli))el[i2][i1] = 0f;
    }}
    //rgf2.apply00(g11,g11);
    //rgf2.apply00(g12,g12);
    //rgf2.apply00(g22,g22);
    int[][] gx = toGrayIntegers(u1);
    int[][] t1 = toGrayIntegers(fx);
    //int[][] t2 = toGrayIntegers(g12);
    //int[][] t3 = toGrayIntegers(g22);
    //float[][][] gs = new float[][][]{fx,g11,g12,g22};
    System.out.println("g11max="+max(fx));
    System.out.println("g11min="+min(g11));
    System.out.println("g11avg="+sum(fx)/(_n1*_n2));
    float[][][] gs = new float[][][]{fx};//,g11,g12,g22};
    _muin[0] = 0.0f;
    _muin[1] = 0.0f;
    _muin[2] = 0.0f;
    _muin[3] = 0.0f;
    _muout[0] = 0.0f;//2*sum(fx)/(_n1*_n2);//
    _muout[1] = 0.0f;//sum(g11)/(_n1*_n2);
    _muout[2] = 0f;//sum(g12)/(_n1*_n2);
    _muout[3] = 0f;//sum(g22)/(_n1*_n2);
    _sigin[0] = 0.5f;
    _sigin[1] = 1.5f;
    _sigout[0] = 1.0f;
    _sigout[1] = 0.3f;
    updateEvolutionSpeed(gs); 
    for (int iter=0; iter<_outerIters; iter++) {
      System.out.println("iter="+iter);
      boolean converged = cycleOne(gs);
      cycleTwo();
      if(converged){break;}
    }
  }


  public float[][] getLout() {
    int np = _lout.size();
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    int k = 0;
    for (int ip=0; ip<np; ++ip) {
      Point pi = _lout.get(ip);
      if(pi.p1<10||pi.p1>=_n1-2){continue;}
      if(pi.p2<10||pi.p2>=_n2-2){continue;}
      x1[k] = pi.p1;
      x2[k] = pi.p2;
      k++;
    }
    return new float[][]{copy(k,0,x1),copy(k,0,x2)};
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
  private float[][] _gauss;
  private float _gaussThreshold;
  private float _sout=0f,_sin=0f;
  private float[][] _pin = new float[4][256];
  private float[][] _pout = new float[4][256];
  private float[] _muin = new float[4];
  private float[] _muout = new float[4];
  private float[] _sigin = new float[4];
  private float[] _sigout = new float[4];
  private int _outerIters = 100;
  private int _speedIters = 10;
  private int _smoothIters = 5;
  private TriMesh _tmesh = new TriMesh();
  private List<Point> _lin  = new LinkedList<Point>(); //inside boundary points;
  private List<Point> _lout = new LinkedList<Point>(); //outside boundary points;
  private List<Point> _linAdd  = new LinkedList<Point>(); //inside boundary points;
  private List<Point> _loutAdd  = new LinkedList<Point>(); //inside boundary points;

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
      updateEvolutionSpeed(fx); 

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
		int m2 = 2*gw+1;
		int m1 = 2*gw+1;
		float gfScale = 0f;
    double sigmas = 1.0/(sigma*sigma);
		// Rough heuristic: scale by number of elements in filter
    _gauss = new float[m2][m1];
		// In theory could just calculate 1/8th and duplicate instead
		for (int i2=0; i2<m2; ++i2) {
      double d2 = i2-gw;
			for (int i1=0; i1<m1; ++i1) {
        double d1 = i1-gw;
				double ds = d2*d2+d1*d1;
				float gf = (float)(sigmas*Math.exp(-0.5*sigmas*ds));
        _gauss[i2][i1]=gf;
				gfScale += gf;
			}
		}
		_gaussThreshold = gfScale/2f;
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
  /*
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
  */

  /*
  private void updateEvolutionSpeedX(float[][] fx) {
    for (Point pi:_lin) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      //if (fx[p2i][p1i]>1.2f)
      if (fx[p2i][p1i]>0.6f)
        pi.setEvolutionSpeed(-1);
      else
        pi.setEvolutionSpeed(0);
    }
    for (Point pi:_lout) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      if (fx[p2i][p1i]<0.5f)
        pi.setEvolutionSpeed(1);
      else
        pi.setEvolutionSpeed(0);
    }

  }
  */

  private void updateEvolutionSpeed(int[][][] fx) {
    int n3 = fx.length;
    updateParzenDensity(fx);
    for (Point pi:_lin) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      float pin =0f;
      float pout =0f;
      for (int i3=0; i3<n3; ++i3) {
        int fxi = fx[i3][p2i][p1i];
        pin += _pin[i3][fxi];
        pout += _pout[i3][fxi];
      }
      if (pout-pin>0f)
        pi.setEvolutionSpeed(-1);
      else
        pi.setEvolutionSpeed(0);
    }
    for (Point pi:_lout) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      float pin =0f;
      float pout =0f;
      for (int i3=0; i3<n3; ++i3) {
        int fxi = fx[i3][p2i][p1i];
        pin += _pin[i3][fxi];
        pout += _pout[i3][fxi];
      }
      if (pin-pout>0f)
        pi.setEvolutionSpeed(1);
      else
        pi.setEvolutionSpeed(0);
    }

  }

  private void updateEvolutionSpeed(float[][][] fx) {
    int n3 = fx.length;
    for (Point pi:_lin) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      float ext = 0f;
      for (int i3=0; i3<n3; ++i3)
        ext += fx[i3][p2i][p1i];
      if (ext>0f)
        pi.setEvolutionSpeed(-1);
      else
        pi.setEvolutionSpeed(0);
    }
    for (Point pi:_lout) {
      int p1i = pi.p1;
      int p2i = pi.p2;
      float ext = 0f;
      for (int i3=0; i3<n3; ++i3)
        ext += fx[i3][p2i][p1i];
      if (ext<0f)
        pi.setEvolutionSpeed(1);
      else
        pi.setEvolutionSpeed(0);
    }

  }


  private float getGaussianDensity(float x, float mu, float sigma) {
    float dx = x-mu;
    dx = dx*dx;
    sigma = 1f/sigma;
    float pi = (float)(1f/sqrt(2*Math.PI));
    return (sigma*pi)*exp(-0.5f*dx*sigma*sigma);
  }


  private void updateParzenDensity(int[][][] fx) {
    int n3 = fx.length;
    for (int i3=0; i3< n3; i3++) {
      updateParzenDensity(fx[i3],_pin[i3],_pout[i3]);
    }
  }

  private void updateParzenDensity(
    final int[][] fx, final float[] pin, final float[] pout) {
    final float din = 1f/_sin;
    final float dout = 1f/_sout;
    Parallel.loop(255,new Parallel.LoopInt() { // i1 = 0, 2, 4, ...
    public void compute(int i) {
      for (int i2=0; i2<_n2; i2++) {
      for (int i1=0; i1<_n1; i1++) {
        int fxi = fx[i2][i1];
        int phi = _phi[i2][i1];
        if (fxi==i) {
          if (phi==-3) pin[i] += din;
          if (phi==3) pout[i] += dout;
        }
      }}
    }});
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(8);
    rgf.apply0(pin,pin);
    rgf.apply0(pout,pout);
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
		float f = 0f;
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
      _sout -= 1f;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p2i][p1p]==3) {
      Point pb = new Point(p1p,p2i);
      _loutAdd.add(pb);
      _phi[p2i][p1p] = 1;
      _sout -= 1f;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p2m][p1i]==3) {
      Point pl = new Point(p1i,p2m);
      _loutAdd.add(pl);
      _phi[p2m][p1i] = 1;
      _sout -= 1f;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p2p][p1i]==3) {
      Point pr = new Point(p1i,p2p);
      _loutAdd.add(pr);
      _phi[p2p][p1i] = 1;
      _sout -= 1f;
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
      _sin -= 1f;
    }
    // check the neighbor below
    int p1p = p1i+1; 
    if(p1p<_n1&&_phi[p2i][p1p]==-3) {
      Point pb = new Point(p1p,p2i);
      _linAdd.add(pb);
      _phi[p2i][p1p] = -1;
      _sin -= 1f;
    }
    // check the left neighbor
    int p2m = p2i-1; 
    if(p2m>=0&&_phi[p2m][p1i]==-3) {
      Point pl = new Point(p1i,p2m);
      _linAdd.add(pl);
      _phi[p2m][p1i] = -1;
      _sin -= 1f;
    }
     // check the right neighbor
    int p2p = p2i+1; 
    if(p2p<_n2&&_phi[p2p][p1i]==-3) {
      Point pr = new Point(p1i,p2p);
      _linAdd.add(pr);
      _phi[p2p][p1i] = -1;
      _sin -= 1f;
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
        _sin += 1f;
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
        _sout += 1f;
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



  /**
   * A point with 2D integer coordinates
   */
  public class Point {
	  public int p1; //1st coordinate
	  public int p2; //2nd coordinate
    public float fe; //evolution speed
    public float fs; //smooth speed
	  /**
	   * Constructor
	   * @param p1 The 1st coordinate
	   * @param p2 The 2nd coordinate
	   */
	  public Point(int p1, int p2) {
	  	this.p1 = p1;
	  	this.p2 = p2;
	  }

    public void setEvolutionSpeed(float fe) {
      this.fe = fe;
    }
    public void setSmoothSpeed(float fs) {
      this.fs = fs;
    }
  }
}
