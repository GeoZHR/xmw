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

public class ActiveContour2 {

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

  public ActiveContour2(int n1, int n2, float c1, float c2, float d) {
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
    resampleSnakeX(n1,n2);
  }


  public ActiveContour2(
    float c1, float c2, float[] x1, float[] x2) 
  {
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

  public Snake updateSnake(int niter, float[][] v1, float[][] v2) {
    for (int it=0; it<niter; it++) {
      updateSnake(v1,v2);
    }
    return _snake;
  }

  public void updateSnake(float[][] v1, float[][] v2) {
    float[] x1 = _snake.getArrayX1();
    float[] x2 = _snake.getArrayX2();
    float[] u1 = _snake.getArrayU1();
    float[] u2 = _snake.getArrayU2();
    _snake.x1.clear();
    _snake.x2.clear();
    _snake.u1.clear();
    _snake.x2.clear();
    int nc = x1.length;
    int n2 = v1.length;
    int n1 = v1[0].length;
    for (int ic=0; ic<nc; ++ic) {
      int ip = ic+1; if(ip==nc) ip=0;
      int im = ic-1; if(im<0)   im=nc-1;
      float x1i = x1[ic];
      float x1p = x1[ip];
      float x1m = x1[im];
      float x2i = x2[ic];
      float x2p = x2[ip];
      float x2m = x2[im];
      float x1a = (x1i+x1p+x1m)/3f;
      float x2a = (x2i+x2p+x2m)/3f;
      float u1i = u1[ic];
      float u2i = u2[ic];
      float ai1 = x1a-x1i;
      float ai2 = x2a-x2i;
      float uas = ai1*u1i+ai2*u2i;
      float fn1 = uas*u1i;
      float fn2 = uas*u2i;
      float ft1 = ai1-fn1;
      float ft2 = ai2-fn2;
      int i1i = round(x1i);
      int i2i = round(x2i);
      i1i = min(i1i,n1-1);
      i2i = min(i2i,n2-1);
      i1i = max(i1i,0);
      i2i = max(i2i,0);
      float v1i = v1[i2i][i1i];
      float v2i = v2[i2i][i1i];
      float t1i = -v2i;
      float t2i =  v1i;
      float tui = t1i*u1i+t2i*u2i;
      if(tui<0f){t1i=-t1i;t2i=-t2i;}
      float x1u = x1i+fn1+ft1+4f*v1i+t1i;
      float x2u = x2i+fn2+ft2+4f*v2i+t2i;
      x1u = min(x1u,n1-1);
      x2u = min(x2u,n2-1);
      x1u = max(x1u,0);
      x2u = max(x2u,0);
      _snake.x1.add(x1u);
      _snake.x2.add(x2u);
    }
    resampleSnakeX(n1,n2);
    //setNormals();
  }

  private void resampleSnakeX(int n1, int n2) {
    resampleSnake(n1,n2);
    /*
    if(removeBadPoints()){
      resampleSnake();
    }
    */
  }


  private void resampleSnake(int n1, int n2) {
    int np = _snake.size;
    float[] ds = new float[np+1];
    float[] x1 = new float[np+1];
    float[] x2 = new float[np+1];
    x1[0] = _snake.x1.get(0);
    x2[0] = _snake.x2.get(0);
    int k = 0;
    for (int ip=1; ip<=np; ++ip) {
      float x1i,x2i;
      if(ip==np){
        x1i = _snake.x1.get(0);
        x2i = _snake.x2.get(0);
      } else {
        x1i = _snake.x1.get(ip);
        x2i = _snake.x2.get(ip);
      }
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
    _snake.size = n;
    _snake.x1.clear();
    _snake.x2.clear();
    for (int i=0; i<n; ++i) {
      float si = (float)sv[i];
      float x1i = cx1.interpolate(si);
      float x2i = cx2.interpolate(si);
      x1i = min(x1i,n1-1);
      x2i = min(x2i,n2-1);
      x1i = max(x1i,   0);
      x2i = max(x2i,   0);
      _snake.x1.add(x1i);
      _snake.x2.add(x2i);
    }
    _snake.u1.clear();
    _snake.u2.clear();
    float[] x1a = _snake.getArrayX1();
    float[] x2a = _snake.getArrayX2();
    for (int i=0; i<n; ++i) {
      int ip = i+1; if(ip==n) ip=0;
      int im = i-1; if(im<0)  im=n-1;
      float g1 = x1a[ip]-x1a[im];
      float g2 = x2a[ip]-x2a[im];
      float gs = sqrt(g1*g1+g2*g2);
      if(gs>0.0f){g1/=gs;g2/=gs;}
      _snake.u1.add(-g2);
      _snake.u2.add(g1);
    }
  }

  private void setNormals() {
    _snake.u1.clear();
    _snake.u2.clear();
    float[] x1a = _snake.getArrayX1();
    float[] x2a = _snake.getArrayX2();
    int n = _snake.size;
    for (int i=0; i<n; ++i) {
      int ip = i+1; if(ip==n) ip=0;
      int im = i-1; if(im<0)  im=n-1;
      float g1 = x1a[ip]-x1a[im];
      float g2 = x2a[ip]-x2a[im];
      float gs = sqrt(g1*g1+g2*g2);
      if(gs>0.0f){g1/=gs;g2/=gs;}
      _snake.u1.add(-g2);
      _snake.u2.add(g1);
    }
  }

  private Snake _snake = new Snake();


}
