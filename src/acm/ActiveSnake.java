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

public class ActiveSnake {

  public static class Snake {
    public int size=0;
    public ArrayList<Float> x1 = new ArrayList<Float>();
    public ArrayList<Float> x2 = new ArrayList<Float>();
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

  public ActiveSnake(float[] x1, float[] x2) {
    int n = x1.length;
    _snake.size = n;
    for (int i=0; i<n; ++i) {
      _snake.x1.add(x1[i]);
      _snake.x2.add(x2[i]);
    }
    resampleSnake();
  }

  public void setIterations(int exiter) {
    _exiter = exiter;
  }

  public Snake releaseSnake(Sampling s1, Sampling s2, 
    float kamma, float gamma, float  alpha, 
    float beta, float[][] f1, float[][] f2, float[][] p, float[][] el) 
  {
    for (int iter=0; iter<_exiter; ++iter) {
      int np = _snake.size;
      float[][] b = new float[2][np];
      float[][] r = new float[2][np];
      Smoother smoother = new Smoother(4.0f);
      float[][] ps = updateSlopes(s1,s2,_snake,p,el);
      VecArrayFloat2 vb = new VecArrayFloat2(b);
      VecArrayFloat2 vr = new VecArrayFloat2(r);
      CgSolver cs = new CgSolver(_small,_initer);
      A2 a2 = new A2(smoother,gamma,alpha,beta,ps);
      makeRhs(s1,s2,kamma,gamma,_snake,f1,f2,b);
      smoother.applyTranspose(b);
      cs.solve(a2,vb,vr);
      smoother.apply(r);
      r[0][np-1] = r[0][0];
      r[1][np-1] = r[1][0];

      _snake.x1.clear();
      _snake.x2.clear();
      for (int ip=0; ip<np; ++ip) {
        _snake.x1.add(r[0][ip]);
        _snake.x2.add(r[1][ip]);
      }
      resampleSnake();
    }
    return _snake;
  }

  ////////////////////////////////////////////////////////////
  private int _exiter = 100;
  private int _initer = 200;
  private float _small = 0.01f;
  private Snake _snake = new Snake();

  private float[][] updateSlopes(
    Sampling s1, Sampling s2, Snake sk, float[][] p, float[][] w) 
  {
    int np = sk.size;
    float[] g1 = new float[np];
    float[] g2 = new float[np];
    float[] x1 = sk.getArrayX1();
    float[] x2 = sk.getArrayX2();
    float[][] ps = new float[2][np];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply1(x1,g1);
    rgf.apply1(x2,g2);
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float g1i = g1[ip];
      float g2i = g2[ip];
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      float pi = si.interpolate(s1,s2,p,x1i,x2i);
      float wi = si.interpolate(s1,s2,w,x1i,x2i);
      ps[0][ip] = -pi;
      ps[1][ip] =  wi;
    }
    return ps;
  }

  /*

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

  */

  private static class A2 implements CgSolver.A {
    A2(Smoother smoother, float gama, float alpha, float beta, float[][] ps) {
      _ps = ps;
      _beta = beta;
      _gama = gama;
      _alpha = alpha;
      _smoother = smoother;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      applyLhs(_gama,_alpha,_beta,_ps,z,y);
      _smoother.applyTranspose(y);
    }

    private float[][] _ps;
    private float _gama;
    private float _beta;
    private float _alpha;
    private Smoother _smoother;
  }


    // Smoother used as a preconditioner. 
  private static class Smoother {
    public Smoother(float sigma) {
      _sigma = sigma;
    }

    public void apply(float[][] x) {
      smooth(_sigma,x[0]);
      smooth(_sigma,x[1]);
    }

    public void applyTranspose(float[][] x) {
      smooth(_sigma,x[0]);
      smooth(_sigma,x[1]);
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
    float gama, float alpha, float beta, float[][] ps, float[][] x, float[][] y) 
  {
    int n = x[0].length;
    float[][] y1 = new float[2][n];
    float[][] y2 = new float[2][n];
    float[][] sx = mul(gama,x);
    zero(y);
    applyLhs(alpha,ps,x,y1);
    applyLaplace(x,y2);
    y2 = mul(y2,beta);
    add(y1,y2,y);
    add(sx,y,y);
  }

  private static void applyLhs(
    float alpha, float[][] p, float[][] x, float[][] y) 
  {
    float x1i,x2i;
    zero(y);
    float[] x1 = x[0];
    float[] y1 = y[0];
    float[] x2 = x[1];
    float[] y2 = y[1];
    int n = x[0].length;
    for (int i=1; i<n; ++i) {
      float pi = p[0][i];
      float wi = p[1][i];
      float ps = pi*pi;
      float ws = wi*wi;
      float d11 = (pi+1f)*ws;
      float d22 = (pi+ps)*ws;
      x1i  = x1[i  ];
      x1i -= x1[i-1];
      x2i  = x2[i  ];
      x2i -= x2[i-1];

      //x1i *= (alpha+pi+1f);
      //x2i *= (alpha+pi+ps);
      float y1i = d11*x1i+d22*x2i+alpha*x1i; 
      float y2i = d11*x1i+d22*x2i+alpha*x2i; 

      y1[i-1] -= y1i;
      y1[i  ]  = y1i;
      y2[i-1] -= y2i;
      y2[i  ]  = y2i;
    }
  }

  private static void applyLaplace(
    float[][] x, float[][] y) 
  {
    float x1i,x2i;
    float[] x1 = x[0];
    float[] x2 = x[1];
    float[] y1 = y[0];
    float[] y2 = y[1];
    int n = x[0].length;
    float[] t1 = new float[n];
    float[] t2 = new float[n];
    for (int i=1; i<n; ++i) {
      x1i  = x1[i  ];
      x1i -= x1[i-1];
      x2i  = x2[i  ];
      x2i -= x2[i-1];

      t1[i-1] -= x1i;
      t1[i  ]  = x1i;
      t2[i-1] -= x2i;
      t2[i  ]  = x2i;

    }
    float t1i,t2i;
    for (int i=1; i<n; ++i) {
      t1i  = t1[i  ];
      t1i -= t1[i-1];
      t2i  = t2[i  ];
      t2i -= t2[i-1];

      y1[i-1] -= t1i;
      y1[i  ]  = t1i;
      y2[i-1] -= t2i;
      y2[i  ]  = t2i;
    }
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
    float beta, float[] x, float[] y) 
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

  private static void makeRhs(Sampling s1, Sampling s2, 
    float kamma, float gamma, Snake sk, float[][] u1, float[][] u2, float[][] y) 
  {
    int np = sk.size;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      y[0][ip] = -kamma*si.interpolate(s1,s2,u1,x1i,x2i)+x1i*gamma;
      y[1][ip] = -kamma*si.interpolate(s1,s2,u2,x1i,x2i)+x2i*gamma;
    }
  }


  private static void makeRhs1(Sampling s1, Sampling s2, 
    float kamma, float gamma, Snake sk, float[][] u1, float[] v, float[] y) 
  {
    int np = sk.size;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      y[ip] = -kamma*si.interpolate(s1,s2,u1,x1i,x2i)+x1i*gamma;
    }
    for (int ip=1; ip<np; ++ip) {
      y[ip] += v[ip]-v[ip-1];
    }
    y[0] += v[1]-v[0];
  }

  private static void makeRhs2(
    Sampling s1, Sampling s2, 
    float kamma, float gamma, Snake sk, float[][] u2,float[] v, float[] y) 
  {
    int np = sk.size;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      y[ip] = -kamma*si.interpolate(s1,s2,u2,x1i,x2i)+x2i*gamma;
    }
    for (int ip=1; ip<np; ++ip) {
      y[ip] += v[ip]-v[ip-1];
    }
    y[0] += v[1]-v[0];

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
    double d = 0.5;
    double f = 0.00;
    double l = ds[k];
    int n = (int)round((l-f)/d);
    Sampling ss = new Sampling(n,d,f);
    double[] sv = ss.getValues();
    CubicInterpolator cx1 = new CubicInterpolator(ds,x1);
    CubicInterpolator cx2 = new CubicInterpolator(ds,x2);
    _snake.size = n;
    _snake.x1.clear();
    _snake.x2.clear();
    for (int i=0; i<n; ++i) {
      float si = (float)sv[i];
      _snake.x1.add(cx1.interpolate(si));
      _snake.x2.add(cx2.interpolate(si));
    }
    _snake.x1.add(cx1.interpolate(0f));
    _snake.x2.add(cx2.interpolate(0f));
  }

}
