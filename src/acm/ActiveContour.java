/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package acm;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Active contour
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.18
 */

public class ActiveContour {

  public static class Snake {
    public float c1;
    public float c2;

    public int size=0;
    public ArrayList<Float> x1 = new ArrayList<Float>();
    public ArrayList<Float> u1 = new ArrayList<Float>();

    public ArrayList<Float> x2 = new ArrayList<Float>();
    public ArrayList<Float> u2 = new ArrayList<Float>();

    public float[] getArrayU1() {
      float[] ua = new float[size];
      for (int i=0; i<size; ++i) {
        ua[i] = u1.get(i);
      }
      return ua;
    }

    public float[] getArrayU2() {
      float[] ua = new float[size];
      for (int i=0; i<size; ++i) {
        ua[i] = u2.get(i);
      }
      return ua;
    }


    public float[] getArrayX1() {
      float[] xa = new float[size];
      for (int i=0; i<size; ++i) {
        xa[i] = x1.get(i);
      }
      return xa;
    }

    public float[] getArrayX2() {
      float[] xa = new float[size];
      for (int i=0; i<size; ++i) {
        xa[i] = x2.get(i);
      }
      return xa;
    }

  }

  public Snake getSnake() {
    return _snake;
  }

  public ActiveContour(Sampling s1, Sampling s2, float c1, float c2, float d) {
    _s1 = s1;
    _s2 = s2;
    _snake.c1=c1;
    _snake.c2=c2;
    _snake.size = 5;
    float[] x1 = new float[]{c1-d,c1  ,c1+d,c1  ,c1-d};
    float[] x2 = new float[]{c2  ,c2+d,c2  ,c2-d,c2  };
    _snake.x1.clear();
    _snake.x2.clear();
    _snake.u1.clear();
    _snake.u2.clear();
    for (float x1i:x1) {_snake.x1.add(x1i);}
    for (float x2i:x2) {_snake.x2.add(x2i);}
    resampleSnakeX();
  }


  public ActiveContour(
    Sampling s1, Sampling s2, float c1, float c2, float[] x1, float[] x2) 
  {
    _s1 = s1;
    _s2 = s2;
    int n = x1.length;
    _snake.c1=c1;
    _snake.c2=c2;
    _snake.size = n;
    _snake.x1.clear();
    _snake.x2.clear();
    _snake.u1.clear();
    _snake.u2.clear();
    for (int i=0; i<n; ++i) {
      _snake.x1.add(x1[i]);
      _snake.x2.add(x2[i]);
    }
    setNormals();
    //resampleSnakeX();
  }

  public void setIterations(int exiter) {
    _exiter = exiter;
  }

  public void setParameters(float kamma, float gamma, float alpha, float beta) {
    _kamma = kamma;
    _gamma = gamma;
    _alpha = alpha;
    _beta  = beta;
  }

  public void exForce(float sig1, float sig2,
    float[][] g, float[][][] fp, float[][][] fn) 
  {
    LocalOrientFilter lof = new LocalOrientFilter(sig1,sig2);
    lof.applyForNormal(g,fn[0],fn[1]);
    fp[0] = mul(fn[1],-1f);
    fp[1] = copy(fn[0]);
  }

  public void firstUpdate(float[][][] ft, float[][][] fn) {
    setIterations(1);
    releaseSnake(ft,fn,true);
    setIterations(_initer);
  }


  public Snake releaseSnake(
    float[][][] fp, float[][][] fn, boolean dc) 
  {
    float[] A = computeA(_alpha,_beta);
    for (int iter=0; iter<_exiter; ++iter) {
      int np = _snake.size;
      float[] b1 = new float[np];
      float[] r1 = new float[np];
      float[] b2 = new float[np];
      float[] r2 = new float[np];
      //float[] ba = fillfloat(beta,np);
      //float[] aa = fillfloat(alpha,np);
      Smoother smoother = new Smoother(2.0f);
      VecArrayFloat1 vb1 = new VecArrayFloat1(b1);
      VecArrayFloat1 vr1 = new VecArrayFloat1(r1);
      CgSolver cs1 = new CgSolver(_small,_initer);
      //A1 a11 = new A1(smoother,gamma,aa,ba);
      A1X a11 = new A1X(smoother,_gamma,A);
      //makeRhs1(s1,s2,kamma,gamma,_snake,f1,b1);
      makeRhs1(_s1,_s2,_kamma,_gamma,_snake,fp,fn,b1,dc);

      smoother.applyTranspose(b1);
      cs1.solve(a11,vb1,vr1);
      smoother.apply(r1);
      r1[np-1] = r1[0];

      VecArrayFloat1 vb2 = new VecArrayFloat1(b2);
      VecArrayFloat1 vr2 = new VecArrayFloat1(r2);
      CgSolver cs2 = new CgSolver(_small,_initer);
      //A1 a12 = new A1(smoother,gamma,aa,ba);
      A1X a12 = new A1X(smoother,_gamma,A);
      //makeRhs2(s1,s2,kamma,gamma,_snake,f2,b2);
      makeRhs2(_s1,_s2,_kamma,_gamma,_snake,fp,fn,b2,dc);
      smoother.applyTranspose(b2);
      cs2.solve(a12,vb2,vr2);
      smoother.apply(r2);
      r2[np-1] = r2[0];

      _snake.x1.clear();
      _snake.x2.clear();
      for (int ip=0; ip<np; ++ip) {
        _snake.x1.add(r1[ip]);
        _snake.x2.add(r2[ip]);
      }
      resampleSnakeX();
    }
    return _snake;
  }

  //get external forces on a snake contour for displaying
  public float[][][] exForceOnSnake(float[][][] fp, float[][][] fn) {
    int np = _snake.size;
    float[][] fps = new float[2][np];
    float[][] fns = new float[2][np];
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float cc1 = _snake.c1;
      float cc2 = _snake.c2;
      float x1i = _snake.x1.get(ip);
      float x2i = _snake.x2.get(ip);
      float u1i = _snake.u1.get(ip);
      float u2i = _snake.u2.get(ip);
      float fp1 = si.interpolate(_s1,_s2,fp[0],x1i,x2i);
      float fp2 = si.interpolate(_s1,_s2,fp[1],x1i,x2i);
      float fn1 = si.interpolate(_s1,_s2,fn[0],x1i,x2i);
      float fn2 = si.interpolate(_s1,_s2,fn[1],x1i,x2i);
      float xc1 = x1i-cc1;
      float xc2 = x2i-cc2;
      if((xc1*fp1+xc2*fp2)<0f){fp1=-fp1;fp2=-fp2;}
      if((u1i*fn1+u2i*fn2)<0f){fn1=-fn1;fn2=-fn2;}
      fps[0][ip] = fp1;
      fps[1][ip] = fp2;
      fns[0][ip] = fn1;
      fns[1][ip] = fn2;
    }
    return new float[][][]{fps,fns};
  }


  ////////////////////////////////////////////////////////////
  Sampling _s1,_s2;
  private float _kamma = 1.0f;
  private float _gamma = 0.5f;
  private float _alpha = 0.1f;
  private float _beta  = 0.1f;
  private int _exiter = 100;
  private int _initer = 200;
  private float _small = 0.01f;
  private Snake _snake = new Snake();

  private float[][] updateNormals(
    Sampling s1, Sampling s2, Snake sk, float[][] u1, float[][] u2) 
  {
    int np = sk.size;
    float[] g1 = new float[np];
    float[] g2 = new float[np];
    float[] x1 = sk.getArrayX1();
    float[] x2 = sk.getArrayX2();
    float[][] us = new float[2][np];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply1(x1,g1);
    rgf.apply1(x2,g2);
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float g1i = g1[ip];
      float g2i = g2[ip];
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      float u1i = si.interpolate(s1,s2,u1,x1i,x2i);
      float u2i = si.interpolate(s1,s2,u2,x1i,x2i);
      float usi = sqrt(u1i*u1i+u2i*u2i);
      u1i /= usi;
      u2i /= usi;
      if((u1i*g2i-u2i*g1i)<0f) {
        u1i = -u1i;
        u2i = -u2i;
      }
      us[0][ip] = u1i;
      us[1][ip] = u2i;
    }
    return us;
  }

  private float[] computeA(float alpha, float beta) {
    float[] A = new float[3];
    A[2] =                 beta;
    A[1] =     -alpha-4.0f*beta;
    A[0] = 2.0f*alpha+6.0f*beta;
    return A;
  }

  private static class A1X implements CgSolver.A {
    A1X(Smoother smoother, float gama, float[] A) {
      _A = A;
      _gama = gama;
      _smoother = smoother;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      v1y.zero();
      applyA(_gama,_A,z,y);
      _smoother.applyTranspose(y);
    }

    private float _gama;
    private float[] _A;
    private Smoother _smoother;
  }


  private static class A1 implements CgSolver.A {
    A1(Smoother smoother, float gama, float[] alpha, float[] beta) {
      _gama = gama;
      _alpha = alpha;
      _beta = beta;
      _smoother = smoother;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      v1y.zero();
      applyLhs(_gama,_alpha,_beta,z,y);
      _smoother.applyTranspose(y);
    }

    private float _gama;
    private float[] _beta;
    private float[] _alpha;
    private Smoother _smoother;
  }


    // Smoother used as a preconditioner. 
  private static class Smoother {
    public Smoother(float sigma) {
      _sigma = sigma;
    }

    public void apply(float[] x) {
      smooth(_sigma,x);
    }

    public void applyTranspose(float[] x) {
      smooth(_sigma,x);
    }

    private float _sigma;
  }

  private static void smooth(float sigma, float[] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
  }

  private static void applyA
    (float gamma, float[] A, float[] x, float[] y) 
  {
    int n = x.length;
    float a2 = A[2];
    float a1 = A[1];
    float a0 = A[0]+gamma;
    float[] xp = new float[n+4];
    xp[0  ] = x[n-2];
    xp[1  ] = x[n-1];
    xp[n+2] = x[0  ];
    xp[n+3] = x[1  ];
    for (int i=2; i<n+2; ++i)
      xp[i] = x[i-2];
    for (int i=2; i<n+2; ++i) {
      float xm2 = xp[i-2]; 
      float xm1 = xp[i-1]; 
      float xi0 = xp[i  ]; 
      float xp1 = xp[i+1]; 
      float xp2 = xp[i+2]; 
      y[i-2] = a0*xi0+a1*(xm1+xp1)+a2*(xm2+xp2);
    }
  }

  private static void applyLhs(
    float gama, float[] alpha, float[] beta, float[] x, float[] y) 
  {
    int n = x.length;
    float[] y1 = new float[n];
    float[] y2 = new float[n];
    float[] y3 = new float[n];
    float[] sx = mul(gama,x);
    zero(y);
    applyLhs1(alpha,x,y1);
    applyLhs2(beta,x,y2);
    applyLhs3(x,y3);
    //add(y1,y2,y);
    add(y3,y,y);
    add(sx,y,y);
  }

  private static void applyLhs3(
    float[] x, float[] y) 
  {
    float xi;
    int n = x.length;
    for (int i=1; i<n; ++i) {
      xi  = x[i  ];
      xi -= x[i-1];
      y[i-1] -= xi;
      y[i  ]  = xi;
    }
    mul(y,-1f,y);
  }


  private static void applyLhs1(
    float[] alpha, float[] x, float[] y) 
  {
    float xi;
    int n = x.length;
    for (int i=1; i<n; ++i) {
      xi  = x[i  ];
      xi -= x[i-1];
      xi *= alpha[i];
      y[i-1] -= xi;
      y[i  ]  = xi;
    }
    //mul(y,-1.0f,y);
  }

  private static void applyLhs2(
    float[] beta, float[] x, float[] y) 
  {
    float xi;
    int n = x.length;
    float[] t = new float[n];
    for (int i=1; i<n; ++i) {
      xi  = x[i  ];
      xi -= x[i-1];
      t[i-1] -= xi;
      t[i  ]  = xi;
    }
    float ti;
    mul(beta,t,t);
    for (int i=1; i<n; ++i) {
      ti  = t[i  ];
      ti -= t[i-1];
      y[i-1] -= ti;
      y[i  ]  = ti;
    }
  }

  private static void makeRhs1(Sampling s1, Sampling s2, 
    float kamma, float gamma, Snake sk, float[][][] fp, 
    float[][][] fn, float[] y,boolean dc) 
  {
    int np = sk.size;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float cc1 = sk.c1;
      float cc2 = sk.c2;
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      float u1i = sk.u1.get(ip);
      float u2i = sk.u2.get(ip);

      float fp1 = si.interpolate(s1,s2,fp[0],x1i,x2i);
      float fn1 = si.interpolate(s1,s2,fn[0],x1i,x2i);

      float fp2 = si.interpolate(s1,s2,fp[1],x1i,x2i);
      float fn2 = si.interpolate(s1,s2,fn[1],x1i,x2i);
      float xc1 = x1i-cc1;
      float xc2 = x2i-cc2;
      if((xc1*fp1+xc2*fp2)<0f) {fp1=-fp1;}
      if(dc) {
        if((u1i*fn1+u2i*fn2)<0f) {fn1=-fn1;}
      }
      y[ip] = kamma*(fp1+fn1)+x1i*gamma;
    }
  }

  private static void makeRhs2(Sampling s1, Sampling s2, 
    float kamma, float gamma, Snake sk, float[][][] fp, 
    float[][][] fn, float[] y, boolean dc) 
  {
    int np = sk.size;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float cc1 = sk.c1;
      float cc2 = sk.c2;
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      float u1i = sk.u1.get(ip);
      float u2i = sk.u2.get(ip);

      float fp1 = si.interpolate(s1,s2,fp[0],x1i,x2i);
      float fn1 = si.interpolate(s1,s2,fn[0],x1i,x2i);

      float fp2 = si.interpolate(s1,s2,fp[1],x1i,x2i);
      float fn2 = si.interpolate(s1,s2,fn[1],x1i,x2i);
      float xc1 = x1i-cc1;
      float xc2 = x2i-cc2;
      if((xc1*fp1+xc2*fp2)<0f) {fp2=-fp2;}
      if(dc) {
        if((u1i*fn1+u2i*fn2)<0f) {fn2=-fn2;}
      }
      y[ip] = kamma*(fp2+fn2)+x2i*gamma;
    }
  }


  private static void makeRhs1(Sampling s1, Sampling s2, 
    float kamma, float gamma, Snake sk, float[][] u1, float[] y) 
  {
    int np = sk.size;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      y[ip] = -kamma*si.interpolate(s1,s2,u1,x1i,x2i)+x1i*gamma;
    }
  }

  private static void makeRhs2(
    Sampling s1, Sampling s2, 
    float kamma, float gamma, Snake sk, float[][] u2,float[] y) 
  {
    int np = sk.size;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      y[ip] = -kamma*si.interpolate(s1,s2,u2,x1i,x2i)+x2i*gamma;
    }
  }

  private void resampleSnakeX() {
    resampleSnake();
    /*
    if(removeBadPoints()){
      resampleSnake();
    }
    */
  }


  private void resampleSnake() {
    int np = _snake.size;
    float[] ds = new float[np];
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    x1[0] = _snake.x1.get(0);
    x2[0] = _snake.x2.get(0);
    int k = 0;
    for (int ip=1; ip<np; ++ip) {
      float x1i = _snake.x1.get(ip);
      float x2i = _snake.x2.get(ip);
      float x1m = x1[ip-1];
      float x2m = x2[ip-1];
      float dx1 = x1i-x1m;
      float dx2 = x2i-x2m;
      float dsi = sqrt(dx1*dx1+dx2*dx2);
      if(dsi>0.0f) {
        k++;
        x1[k] = x1i;
        x2[k] = x2i;
        ds[k] = ds[ip-1]+dsi;
      }
    }
    System.out.println("k="+(k+1));
    x1 = copy(k+1,x1);
    x2 = copy(k+1,x2);
    ds = copy(k+1,ds);
    double d = 1.0;
    double f = 0.00;
    double l = ds[k];
    int n = (int)round((l-f)/d);
    Sampling ss = new Sampling(n,d,f);
    double[] sv = ss.getValues();
    CubicInterpolator cx1 = new CubicInterpolator(ds,x1);
    CubicInterpolator cx2 = new CubicInterpolator(ds,x2);
    _snake.size = n+1;
    _snake.x1.clear();
    _snake.x2.clear();
    for (int i=0; i<n; ++i) {
      float si = (float)sv[i];
      _snake.x1.add(cx1.interpolate(si));
      _snake.x2.add(cx2.interpolate(si));
    }
    _snake.x1.add(cx1.interpolate(0f));
    _snake.x2.add(cx2.interpolate(0f));

    _snake.u1.clear();
    _snake.u2.clear();
    float[] x1a = _snake.getArrayX1();
    float[] x2a = _snake.getArrayX2();
    for (int i=0; i<n; ++i) {
      float g1 = x1a[i+1]-x1a[i];
      float g2 = x2a[i+1]-x2a[i];
      float gs = sqrt(g1*g1+g2*g2);
      if(gs>0.0f){g1/=gs;g2/=gs;}
      _snake.u1.add(-g2);
      _snake.u2.add(g1);
    }
    float g1 = x1a[1]-x1a[0];
    float g2 = x2a[1]-x2a[0];
    float gs = sqrt(g1*g1+g2*g2);
    if(gs>0.0f){g1/=gs;g2/=gs;}
    _snake.u1.add(-g2);
    _snake.u2.add(g1);
  }

  private void setNormals() {
    _snake.u1.clear();
    _snake.u2.clear();
    float[] x1a = _snake.getArrayX1();
    float[] x2a = _snake.getArrayX2();
    int n = _snake.size;
    for (int i=0; i<n-1; ++i) {
      float g1 = x1a[i+1]-x1a[i];
      float g2 = x2a[i+1]-x2a[i];
      float gs = sqrt(g1*g1+g2*g2);
      if(gs>0.0f){g1/=gs;g2/=gs;}
      _snake.u1.add(-g2);
      _snake.u2.add(g1);
    }
    float g1 = x1a[1]-x1a[0];
    float g2 = x2a[1]-x2a[0];
    float gs = sqrt(g1*g1+g2*g2);
    if(gs>0.0f){g1/=gs;g2/=gs;}
    _snake.u1.add(-g2);
    _snake.u2.add(g1);

  }

  private boolean removeBadPoints() {
    int np = _snake.size;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int[] mp = new int[np];
    int[][] mk = new int[n1][n2];
    int[][] pk = new int[n1][n2];
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    for (int ip=0; ip<np; ++ip) {
      float x1i = _snake.x1.get(ip);
      float x2i = _snake.x2.get(ip);
      x1[ip] = x1i;
      x2[ip] = x2i;
      int i1 = round(x1i);
      int i2 = round(x2i);
      if(i1<0)    {mp[ip]=1;continue;}
      if(i2<0)    {mp[ip]=1;continue;}
      if(i1>n1-1) {mp[ip]=1;continue;}
      if(i2>n2-1) {mp[ip]=1;continue;}
      if(mk[i2][i1]>=1){
        int pp = pk[i2][i1];
        mp[pp] = 1;
        mp[ip] = 1;
      }
      pk[i2][i1] = ip;
      mk[i2][i1] += 1;
    }
    int k = 0;
    _snake.x1.clear();
    _snake.x2.clear();
    for (int ip=0; ip<np; ++ip) {
      if(mp[ip]==1) {continue;}
      _snake.x1.add(x1[ip]);
      _snake.x2.add(x2[ip]);
      k++;
    }
    _snake.size=k;
    if(sum(mp)>0){return true;}
    else{return false;}
  }

}
