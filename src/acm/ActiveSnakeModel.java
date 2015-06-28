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
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class ActiveSnakeModel {

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

  public ActiveSnakeModel(float[] x1, float[] x2) {
    int n = x1.length;
    _snake.size = n;
    for (int i=0; i<n; ++i) {
      _snake.x1.add(x1[i]);
      _snake.x2.add(x2[i]);
    }
    resampleSnake();
  }

  public Snake releaseSnake(Sampling s1, Sampling s2, 
    float kamma, float gamma, float alpha, float beta, float[][] u1, float[][] u2) 
  {
    int np = _snake.size;
    double[] bw = new double[np];
    bw[0] = 2.0*alpha+6.0*beta; 
    bw[1] =    -alpha-4.0*beta; 
    bw[2] = beta;
    bw[np-2] = beta;
    bw[np-1] = -alpha-4.0*beta;
    double[][] A = new double[np][np];
    for (int i1=0; i1<np; ++i1) {
      for (int i2=0; i2<np; ++i2) {
        A[i2][i1] = bw[i2];
      }
      A[i1][i1] += gamma;
      bw = circleShift(bw);
    }
    DMatrix MA = new DMatrix(A);
    DMatrix MI = MA.inverse();
    double[][] ai = MI.get();
    for (int iter=0; iter<_exiter; ++iter) {
      double[] b1 = makeRhs1(s1,s2,kamma,gamma,_snake,u1);
      double[] b2 = makeRhs2(s1,s2,kamma,gamma,_snake,u2);
      float[] x1 = new float[np];
      float[] x2 = new float[np];
      for (int i1=0; i1<np; ++i1) {
      for (int i2=0; i2<np; ++i2) {
        x1[i1] += ai[i2][i1]*b1[i2];
        x2[i1] += ai[i2][i1]*b2[i2];
      }}
      _snake.x1.clear();
      _snake.x2.clear();
      for (int ip=0; ip<np; ++ip) {
        _snake.x1.add(x1[ip]);
        _snake.x2.add(x2[ip]);
      }
    }
    return _snake;
  }

  ////////////////////////////////////////////////////////////
  private int _exiter = 100;
  private int _initer = 200;
  private float _small = 0.01f;
  private Snake _snake = new Snake();

  private double[] circleShift(double[] bw) {
    int np = bw.length;
    double[] bs = new double[np];
    bs[0] = bw[np-1]; 
    for (int ip=1; ip<np; ++ip) {
      bs[ip] = bw[ip-1]; 
    }
    return bs;
  }

  private static double[] makeRhs1(
    Sampling s1, Sampling s2, float kamma, float gama, Snake sk, float[][] u1) 
  {
    int np = sk.size;
    double[] y = new double[np];
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      y[ip] = -kamma*si.interpolate(s1,s2,u1,x1i,x2i)+x1i*gama;
    }
    return y;
  }

  private static double[] makeRhs2(
    Sampling s1, Sampling s2, float kamma, float gama, Snake sk, float[][] u2) 
  {
    int np = sk.size;
    double[] y = new double[np];
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      float x1i = sk.x1.get(ip);
      float x2i = sk.x2.get(ip);
      y[ip] = -kamma*si.interpolate(s1,s2,u2,x1i,x2i)+x2i*gama;
    }
    return y;
  }


  private void resampleSnake() {
    int np = _snake.size;
    float[] ds = new float[np];
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    x1[0] = _snake.x1.get(0);
    x2[0] = _snake.x2.get(0);
    for (int ip=1; ip<np; ++ip) {
      float x1i = _snake.x1.get(ip);
      float x2i = _snake.x2.get(ip);
      x1[ip] = x1i;
      x2[ip] = x2i;
      float x1m = x1[ip-1];
      float x2m = x2[ip-1];
      float dx1 = x1i-x1m;
      float dx2 = x2i-x2m;
      ds[ip] = ds[ip-1]+sqrt(dx1*dx1+dx2*dx2);
    }
    double d = 2.0;
    double f = 0.00;
    double l = ds[np-1];
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

  }
}
