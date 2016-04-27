/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package aii;

import vec.*;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/** 
 * Image-guided 2D acoustic impedance inversion with constraints from well logs
 * @Author: Xinming Wu and Dave Hale, Colorado School of Mines
 * @Version: 2015.10.10
 */

public class FitLogImpWithSeisRef {

  /**
   * Constructor.
   * @param sigma1 smoother half-width for 1st dimension.
   */
  public FitLogImpWithSeisRef(double sigma1) {
    _sigma1 = (float)sigma1;
  }

  /**
   * Sets parameters that control the number of solver iterations.
   * @param small stop iterations when error norm is reduced by this fraction.
   * @param niter maximum number of solver iterations.
   */
  public void setIterations(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  public void setSeisBalance(float sc) {
    _sc = sc;
  }

  public float[] fitImpedance(float[] px, float[] rx) {
    int n1 = rx.length;
    float[] r = copy(px);
    float[] b = new float[n1];
    VecArrayFloat1 vb = new VecArrayFloat1(b);
    VecArrayFloat1 vr = new VecArrayFloat1(r);
    CgSolver cs = new CgSolver(_small,_niter);
    Smoother1 smoother1 = new Smoother1(n1,_sigma1);
    makeRhs(rx,px,b);
    A1 a1 = new A1(smoother1);
    smoother1.applyTranspose(b);
    cs.solve(a1,vb,vr);
    smoother1.apply(r);
    return r;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _small = 0.001f; // stop CG iterations if residuals are small
  private int _niter = 200; // maximum number of inner CG iterations
  private static float _sc = 0.5f;

  private static class A1 implements CgSolver.A {
    A1(Smoother1 s1) { 
      _s1 = s1;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      _s1.apply(z);
      zero(y);
      applyLhs(z,y);
      _s1.applyTranspose(y);
    }
    private Smoother1 _s1;
  }

  private static void applyLhs( float[] x, float[] y) {
    int n1 = x.length;
    float sr =  _sc;
    float sp = 1f-sr;
    float wr=sr*sr; 
    float wp=sp*sp; 
    for (int i1=1; i1<n1; ++i1) {
      float xa=0.0f;
      xa  = x[i1  ];
      xa -= x[i1-1];
      xa *= wr;
      y[i1-1] -= xa;
      y[i1  ]  = xa;
    }
    for (int i1=0; i1<n1; ++i1)
      y[i1] += x[i1]*wp;
  }

  private static void makeRhs(
    float[] r, float[] p, float[] y) 
  {
    zero(y);
    int n1 = y.length;
    float sr =  _sc;
    float sp = 1f-sr;
    float wr=sr*sr; 
    float wp=sp*sp; 
    for (int i1=1; i1<n1; ++i1) {
      float ri = wr*r[i1];
      y[i1  ] += ri;
      y[i1-1] -= ri;
    }
    for (int i1=0; i1<n1; ++i1)
      y[i1  ] += wp*p[i1];
  }


  // Smoother used as a preconditioner.
  private static class Smoother1 {
    public Smoother1(int n1, float sigma1){
      _sigma1 = sigma1;
    }
    public void apply(float[] x) {
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[] x) {
      smooth1(_sigma1,x);
    }
    private float _sigma1;
  }


  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }


}
