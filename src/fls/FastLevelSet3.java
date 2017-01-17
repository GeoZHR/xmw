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
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import stv.*;
import sso.*;

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

  static public float[][][] downSample(int d1, int d2, int d3, float[][][] fx) {
    int m1 = 0;
    int m2 = 0;
    int m3 = 0;
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    for (int i1=0; i1<n1; i1+=d1) m1++;
    for (int i2=0; i2<n2; i2+=d2) m2++;
    for (int i3=0; i3<n3; i3+=d3) m3++;
    float[][][] fxs = new float[m3][m2][m1];
    for (int i1=0,k1=0; i1<n1; i1+=d1,k1++) 
    for (int i2=0,k2=0; i2<n2; i2+=d2,k2++)
    for (int i3=0,k3=0; i3<n3; i3+=d3,k3++)
      fxs[k3][k2][k1] = fx[i3][i2][i1];
    return fxs;

  }

  public TriangleGroup getTriangleGroup(int zm, int zp, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    Sampling s1 = new Sampling(zp-zm+1,1,zm);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    float[][][] fs = copy(zp-zm+1,n2,n3,zm,0,0,fx);
    MarchingCubes mc = new MarchingCubes(s1,s2,s3,fs);
    MarchingCubes.Contour ct = mc.getContour(0);
    return new TriangleGroup(ct.i,ct.x,ct.u);
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

  public void updateLevelSet(int gw, double sigma, float[][][] dp) {
    _gWindow = gw;
    createGaussFilter(gw,sigma);
    //updateEvolutionSpeed(fx);
    for (int iter=0; iter<_outerIters; iter++) {
      System.out.println("iter="+iter);
      //boolean converged = cycleOne(dp);
      //cycleTwo();
      cycleOneX(dp);
      cycleTwoX();
      //if(converged){break;}
    }
  }

  public void updateLevelSet(int gw, double sigma, float[][][] dp, float[][][] ph) {
    _gWindow = gw;
    createGaussFilter(gw,sigma);
    initialize(ph);
    //updateEvolutionSpeed(fx);
    for (int iter=0; iter<_outerIters; iter++) {
      System.out.println("iter="+iter);
      boolean converged = cycleOne(dp);
      cycleTwo();
      //cycleOneX(fx);
      //cycleTwoX();
      //if(converged){break;}
    }
  }

  public void initialize(float[][][] ph) {
    _n3 = ph.length;
    _n2 = ph[0].length;
    _n1 = ph[0][0].length;
    _phi = fillbyte((byte)3,_n1,_n2,_n3);
    for (int i3=0; i3<_n3; i3++) {
    for (int i2=0; i2<_n2; i2++) {
    for (int i1=0; i1<_n1; i1++) {
      float phi = ph[i3][i2][i1];
      _phi[i3][i2][i1] = (byte)phi;
      if(phi==-1) _lin.add(new Point(i1,i2,i3));
      if(phi== 1) _lout.add(new Point(i1,i2,i3));
    }}}
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

  public float[][][][] applyForDip(float sig1, float sig2, float[][][] fx) {
    final int n3 = fx.length;
    final int n2 = fx[0].length; 
    final int n1 = fx[0][0].length; 
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    float p2min = -5;
    float p3min = -5;
    float p2max =  5;
    float p3max =  5;
    LocalOrientFilter lof = new LocalOrientFilter(sig1,sig2,sig2);
    lof.applyForNormal(fx,u1,u2,u3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          if (-u2i<p2min*u1i) u2i = -p2min*u1i;
          if (-u2i>p2max*u1i) u2i = -p2max*u1i;
          if (-u3i<p3min*u1i) u3i = -p3min*u1i;
          if (-u3i>p3max*u1i) u3i = -p3max*u1i;
          if (u1i==0.0f) {
            p2[i3][i2][i1] = (u2i<0.0f)?p2max:p2min;
            p3[i3][i2][i1] = (u3i<0.0f)?p3max:p3min;
          } else {
            p2[i3][i2][i1] = -u2i/u1i;
            p3[i3][i2][i1] = -u3i/u1i;
          }
        }
      }
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2);
    float[][][] g1 = new float[n3][n2][n1];
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    rgf.apply100(p2,g1);
    rgf.apply010(p2,g2);
    rgf.apply001(p2,g3);
    g1 = mul(g1,u1);
    g2 = mul(g2,u2);
    g3 = mul(g3,u3);
    p2 = add(g1,g2);
    p2 = add(g3,p2);

    rgf.apply100(p3,g1);
    rgf.apply010(p3,g2);
    rgf.apply001(p3,g3);
    g1 = mul(g1,u1);
    g2 = mul(g2,u2);
    g3 = mul(g3,u3);
    p3 = add(g1,g2);
    p3 = add(g3,p3);
    return new float[][][][]{p2,p3};
  }

  public void applyForInsAttributes(
    final float[][][] fx, final float[][][] pa, 
    final float[][][] ph, final float[][][] pf){
    final int n3 = fx.length;
    final int n2 = fx[0].length; 
    final int n1 = fx[0][0].length; 
    final float[][][] f1 = pa;
    final float[][][] f2 = ph;
    final float[][][] f3 = pf;
    HilbertTransform ht = new HilbertTransform(5,5);
    ht.applyInFrequency(fx,f1,f2,f3);
    final float pi = (float)(Math.PI);
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
      for (int i2=0; i2<n2; i2++){
        float[] pf32 = pf[i3][i2];
        float[] ph32 = ph[i3][i2];
        float[] pa32 = pa[i3][i2];
        float[] fx32 = fx[i3][i2];
        float[] f132 = f1[i3][i2];
        float[] f232 = f1[i3][i2];
        float[] f332 = f1[i3][i2];
        for (int i1=0; i1<n1; i1++){
          float f1i = f132[i1];
          float f2i = f232[i1];
          float f3i = f332[i1];
          float fxr = fx32[i1];
          float phi = -atan2(f1i,fxr);
          float pai = sqrt(fxr*fxr+f1i*f1i+f2i*f2i+f3i*f3i);
          if(Float.isInfinite(phi)||Float.isNaN(phi)){
            ph32[i1] = 0f;
          } else { ph32[i1] = phi; }
          if(Float.isInfinite(pai)||Float.isNaN(pai)){
            pa32[i1] = 0f;
          } else { pa32[i1] = pai; }

        }
        for (int i1=1; i1<n1; i1++){
          float dpi = ph32[i1]-ph32[i1-1];
          if(dpi<-pi) {dpi+=2*pi;}
          if(dpi> pi) {dpi-=2*pi;}
          pf32[i1] = abs(dpi);
        }
      }
    }}); 
  }


  static public void applyForInsAmp(
    final float[][][] fx, final float[][][] pa){
    final int n3 = fx.length;
    final int n2 = fx[0].length; 
    final int n1 = fx[0][0].length; 
    final HilbertTransformFilter hbt = new HilbertTransformFilter();
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
      for (int i2=0; i2<n2; i2++){
        float[] pa32 = pa[i3][i2];
        float[] fx32 = fx[i3][i2];
        float[] fi = new float[n1];
        hbt.apply(n1,fx32,fi);
        for (int i1=0; i1<n1; i1++){
          float fxi = fi[i1];
          float fxr = fx32[i1];
          float pai = sqrt(fxr*fxr+fxi*fxi);
          if(Float.isInfinite(pai)||Float.isNaN(pai)){
            pa32[i1] = 0f;
          } else { pa32[i1] = pai; }
        }
      }
    }}); 
  }

  public void applyForU1(
    float sig1, float sig2, final float[][][] fx, final float[][][] u1){
    final int n3 = fx.length;
    final int n2 = fx[0].length; 
    final int n1 = fx[0][0].length; 
    float[][][] ut = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(sig1,sig2,sig2);
    lof.applyForNormal(fx,u1,ut,ut);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    zero(ut);
    float[][][] g1 = ut;
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    rgf.apply100(u1,g1);
    rgf.apply010(u1,g2);
    rgf.apply001(u1,g3);
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
      for (int i2=0; i2<n2; i2++){
        float[] g132 = g1[i3][i2];
        float[] g232 = g2[i3][i2];
        float[] g332 = g3[i3][i2];
        for (int i1=0; i1<n1; i1++){
          float g1i = g132[i1];
          float g2i = g232[i1];
          float g3i = g332[i1];
          u1[i3][i2][i1] = sqrt(g1i*g1i+g2i*g2i+g3i*g3i);
        }
      }
    }}); 
  }




  public float[][][] applyForTextureDensity(
    float sigma1, float sigma2, float ev, float[][][] ep, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    RecursiveGaussianFilter rgfGradient = new RecursiveGaussianFilter(1);
    float[][][] g1 = new float[n3][n2][n1];
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    rgfGradient.apply100(gx,g1);
    rgfGradient.apply010(gx,g2);
    rgfGradient.apply001(gx,g3);
    /*
    float[][][] g11 = g1;
    float[][][] g22 = g2;
    float[][][] g33 = g3;
    float[][][] g12 = new float[n3][n2][n1];
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g11,g22,g33,g12,g13,g23};
    RecursiveGaussianFilter rgfSmoother1 = new RecursiveGaussianFilter(sigma1);
    RecursiveGaussianFilter rgfSmoother2 = new RecursiveGaussianFilter(sigma2);
    RecursiveGaussianFilter rgfSmoother3 = new RecursiveGaussianFilter(sigma2);
    for (float[][][] g:gs) {
      rgfSmoother1.apply0XX(g,h);
      rgfSmoother2.applyX0X(h,g);
      rgfSmoother3.applyXX0(g,h);
      copy(h,g);
    }
    gs = new float[][][][]{g11,g22,g33,g12,g23,g13};
    */
    float[][][][] gs = new float[][][][]{gx,g1,g2,g3};
    return applyForDensity(ev,ep,gs);
  }

  private void computeGradientProducts(
    final float[][][] g1, final float[][][] g2, final float[][][] g3,
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33)
  {
    final int n1 = g1[0][0].length;
    final int n2 = g1[0].length;
    final int n3 = g1.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] g1i = g1[i3][i2];
          float[] g2i = g2[i3][i2];
          float[] g3i = g3[i3][i2];
          float[] g11i = g11[i3][i2];
          float[] g12i = g12[i3][i2];
          float[] g13i = g13[i3][i2];
          float[] g22i = g22[i3][i2];
          float[] g23i = g23[i3][i2];
          float[] g33i = g33[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float g1ii = g1i[i1];
            float g2ii = g2i[i1];
            float g3ii = g3i[i1];
            g11i[i1] = g1ii*g1ii;
            g22i[i1] = g2ii*g2ii;
            g33i[i1] = g3ii*g3ii;
            g12i[i1] = g1ii*g2ii;
            g13i[i1] = g1ii*g3ii;
            g23i[i1] = g2ii*g3ii;
          }
        }
      }
    });
  }



  public float[][][] applyForDensity(
    float ev, float[][][] ep, float[][][][] gs) {
    int n4 = gs.length;
    int n3 = gs[0].length;
    int n2 = gs[0][0].length;
    int n1 = gs[0][0][0].length;
    float[][][][] dps = new float[n4][n3][n2][n1];
    for (int i4=0; i4<n4; ++i4)
      dps[i4] = density(ev,ep,gs[i4]);
    if(n4>1) balanceDensity(dps);
    float[][][] dp = new float[n3][n2][n1];
    for (int i4=1; i4<n4; ++i4)
    for (int i3=0; i3<n3; ++i3)
    for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      dp[i3][i2][i1] += dps[i4][i3][i2][i1];
    return dp;
  }



  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _n1;
  private int _n2;
  private int _n3;
  private byte[][][] _phi;
  private int _gWindow;
  private float[][][] _gauss;
  private float _gaussThreshold;
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

  public int[][][] toGrayIntegers(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    int[][][] gx = new int[n3][n2][n1];
    fx = sub(fx,min(fx));
    float fmax = 255f/max(fx);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i3][i2][i1] = round(fx[i3][i2][i1]*fmax);
    }}}
    return gx;
  }

  private void balanceDensity(float[][][][] ds) {
    int n4 = ds.length;
    int n3 = ds[0].length;
    int n2 = ds[0][0].length;
    int n1 = ds[0][0][0].length;
    float[] sc = new float[n4];
    float[] sd = new float[n4];
    for (int i4=0; i4<n4; ++i4) {
      float[][][] d4 = ds[i4];
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if(d4[i3][i2][i1]<0f) {
          sc[i4] += 1f;
          sd[i4] += abs(d4[i3][i2][i1]);
        }
      }}}
    }
    div(sd,sc,sd);
    float sdmax = max(sd);
    for (int i4=0; i4<n4; ++i4) {
      float sci = sdmax/sd[i4];
      mul(ds[i4],sci,ds[i4]);
    }

  }

  public void densityX(
    float f, float[][][] fx, float[][][] gx) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] p1 = new float[256];
    float[] p2 = new float[256];
    int[][][] gg = toGrayIntegers(gx);
    density(f,fx,gg,p1,p2);
    for (int i3=1; i3<n3-1; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      int ggi = gg[i3][i2][i1];
      float di = log(p2[ggi])-log(p1[ggi]);
      if(Float.isNaN(di)) {continue;}
      if(Float.isInfinite(di)) {continue;}
      gx[i3][i2][i1] = di;
    }}}
  }


  public float[][][] density(
    float f, float[][][] fx, float[][][] gx) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] p1 = new float[256];
    float[] p2 = new float[256];
    float[][][] dp = new float[n3][n2][n1];
    int[][][] gg = toGrayIntegers(gx);
    density(f,fx,gg,p1,p2);
    for (int i3=1; i3<n3-1; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
    for (int i1=1; i1<n1-1; ++i1) {
      int ggi = gg[i3][i2][i1];
      float di = log(p2[ggi])-log(p1[ggi]);
      if(Float.isNaN(di)) {continue;}
      if(Float.isInfinite(di)) {continue;}
      dp[i3][i2][i1] = di;
    }}}
    return dp;
  }

  public void density(float f, float[][][] fx, 
    int[][][] gx, float[] p1, float[] p2) {
    double c1 = 0f;
    double c2 = 0f;
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    double[] p1s = new double[256];
    double[] p2s = new double[256];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int gxi = gx[i3][i2][i1];
      float fxi = fx[i3][i2][i1];
      if(fxi<f) {
        c1 += 1.0; p1s[gxi] += 1.0;
      } else { 
        c2 += 1.0; p2s[gxi] += 1.0;
      }
    }}}
    div(p1s,c1,p1s);
    div(p2s,c2,p2s);
    for (int ip=0; ip<256; ++ip) {
      p1[ip] = (float)p1s[ip];
      p2[ip] = (float)p2s[ip];
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2);
    rgf.apply0(p1,p1);
    rgf.apply0(p2,p2);
  }




  private void cycleOneX(final float[][][] fx) {
    for (int iter=0; iter<_speedIters; iter++) {
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
      if(fx3[i2][i1]<0f){switchIn(i3,louti,sp);}
    }
  }
  private void switchOutSlice3(int i3, Iterator<Sample>lini,float[][] fx3) {
    while(lini.hasNext()) {
      Sample sp = lini.next();
      int i1 = sp.p1;
      int i2 = sp.p2;
      if(fx3[i2][i1]>0f){switchOut(i3,lini,sp);}
    }
  }


  private void switchInSlice3(int i3, Iterator<Sample>louti) {
    while(louti.hasNext()) {
      Sample sp = louti.next();
      int i1 = sp.p1;
      int i2 = sp.p2;
      if(smoothSpeed(i1,i2,i3)>_gaussThreshold){switchIn(i3,louti,sp);}
    }
  }

  private void switchOutSlice3(int i3, Iterator<Sample>lini) {
    while(lini.hasNext()) {
      Sample sp = lini.next();
      int i1 = sp.p1;
      int i2 = sp.p2;
      if(smoothSpeed(i1,i2,i3)<_gaussThreshold){switchOut(i3,lini,sp);}
    }
  }


  private boolean cycleOne(float[][][] dp) {
    boolean converged = false;
    for (int iter=0; iter<_speedIters; iter++) {
      //outward evolution
      Iterator<Point> louti = _lout.iterator();
      while(louti.hasNext()) {
        Point pi = louti.next();
        int p1 = pi.p1;
        int p2 = pi.p2;
        int p3 = pi.p3;
        if(dp[p3][p2][p1]<0f) {switchIn(louti,pi);}
      }
      _lout.addAll(0,_loutAdd);
      _loutAdd.clear();
      cleanLin();
      //inward evolution
      Iterator<Point> lini = _lin.iterator();
      while(lini.hasNext()) {
        Point pi = lini.next();
        int p1 = pi.p1;
        int p2 = pi.p2;
        int p3 = pi.p3;
        if(dp[p3][p2][p1]>0f) {switchOut(lini,pi);}
      }
      _lin.addAll(0,_linAdd);
      _linAdd.clear();
      cleanLout();
      //converged = checkConvergence();
      //if(converged) {break;}
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
		float gfScale = 0;
    double sigmas = 1.0/(sigma*sigma);
    double sigmat = sigmas/sigma;
		// Rough heuristic: scale by number of elements in filter
    _gauss = new float[m3][m2][m1];
		// In theory could just calculate 1/8th and duplicate instead
		for (int i3=0; i3<m3; ++i3) {
      double d3  = i3-gw;
		for (int i2=0; i2<m2; ++i2) {
      double d2  = i2-gw;
			for (int i1=0; i1<m1; ++i1) {
        double d1  = i1-gw;
				double ds = d3*d3+d2*d2+d1*d1;
				float gf = (float)(sigmat*exp(-0.5*sigmas*ds));
        _gauss[i3][i2][i1]=gf;
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
		float f = 0;
		for(int d3=d3m; d3<d3p; ++d3) {
		for(int d2=d2m; d2<d2p; ++d2) {
		for(int d1=d1m; d1<d1p; ++d1) {
			if (_phi[p3+d3][p2+d2][p1+d1]<0) 
			  f += _gauss[gw+d3][gw+d2][gw+d1];
		}}}
		p.setSmoothSpeed(f);
	}

	private float smoothSpeed(int p1, int p2, int p3) {
		// Convolve neighbourhood of a point with a gaussian
    int gw = _gWindow;
		int d1m = max(-gw,-p1);
		int d2m = max(-gw,-p2);
		int d3m = max(-gw,-p3);
		int d1p = min(gw+1,_n1-p1);
		int d2p = min(gw+1,_n2-p2);
		int d3p = min(gw+1,_n3-p3);
		float f = 0f;
		for(int d3=d3m; d3<d3p; ++d3) {
		for(int d2=d2m; d2<d2p; ++d2) {
		for(int d1=d1m; d1<d1p; ++d1) {
			if (_phi[p3+d3][p2+d2][p1+d1]<0) 
			  f += _gauss[gw+d3][gw+d2][gw+d1];
		}}}
    return f;
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
    p.ct=0;
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
    p.ct=0;
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
    public float fe; //evolution speed
    public float fs; //smooth speed
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

    public void setEvolutionSpeed(float fe) {
      this.fe = fe;
    }
    public void setSmoothSpeed(float fs) {
      this.fs = fs;
    }
  }

  public class Sample {
	  public int p1; //1st coordinate
	  public int p2; //2nd coordinate
    public int ct=0;
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
