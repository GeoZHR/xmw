package util;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;



/**
 * Optimal path picking
 * @author Xinming Wu 
 * @version 2016.03.23
 */

public class SmoothWithShaping {
 
  public SmoothWithShaping(float sigma) {
    _sigma1 = sigma;
  }

  public float[] smooth( float[] wx, float[] fx) {
    int n1 = fx.length;
    float[] b = new float[n1];
    float[] r = new float[n1];
    makeRhs(wx,fx,b);
    VecArrayFloat1 vb = new VecArrayFloat1(b);
    VecArrayFloat1 vr = new VecArrayFloat1(r);
    Smoother1 smoother1 = new Smoother1(_sigma1);
    A1 a1 = new A1(smoother1,wx);
    CgSolver cs = new CgSolver(0.001,200);
    smoother1.applyTranspose(b);
    cs.solve(a1,vb,vr);
    return r;
  }

  private void makeRhs(float[] wx, float[] fx, float[] b){
    int n1 = fx.length;
    for (int i1=0; i1<n1; ++i1) {
      float wxi = wx[i1];
      b[i1] = fx[i1]*wxi*wxi;
    }
  }

  // Conjugate-gradient operators.
  private static class A1 implements CgSolver.A {
    A1(Smoother1 s1, float[] wp) 
    {
      _s1 = s1;
      _wp = wp;
      float n1 = wp.length;
      _sc = sum(wp)/(n1);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      v1y.zero();
      _s1.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s1.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother1 _s1;
    private float[] _wp;
  }


  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] wp) 
    {
      _s2 = s2;
      _wp = wp;
      float n2 = wp.length;
      float n1 = wp[0].length;
      _sc = 4f*sum(wp)/(n1*n2);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      _s2.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s2.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother2 _s2;
    private float[][] _wp;
  }

  private static void applyLhs(float[] wp, float[] x, float[] y) {
    int n1 = wp.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += wp[i1]*wp[i1]*x[i1];
  }


  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    int n2 = wp.length;
    for (int i2=0; i2<n2; ++i2)
      applyLhs(wp[i2],x[i2],y[i2]);
  }

  private static void addAndScale(float sc, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      addAndScale(sc,x[i2],y[i2]);
    }
  }

  private static void addAndScale(float sc, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1) {
      y[i1] += sc*x[i1];
    }
  }


  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(float sigma1, float sigma2) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,x);
    }
    private float _sigma1,_sigma2;
  }


  private static class Smoother1 {
    public Smoother1(float sigma) {
      _sigma = sigma;
    }
    public void apply(float[] x) {
      smooth1(_sigma,x);
    }
    public void applyTranspose(float[] x) {
      smooth1(_sigma,x);
    }
    private float _sigma;
  }



  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
  }


  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }



  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1;
}
