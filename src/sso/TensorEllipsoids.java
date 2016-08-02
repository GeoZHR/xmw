/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package sso;


import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * make ellipsoids for displaying
 * @author Xinming Wu
 * @version 2016.07.30
 */
public class TensorEllipsoids extends Node {
  public TensorEllipsoids(
    Sampling s1, Sampling s2, Sampling s3,
    EigenTensors3 et, float[][][] ep, float[][] hz) 
  {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _ep = ep;
    _et = et;
    _hz = hz;
    _emax = findMaxEigenvalue();

  }

  /**
   * Sets the maximum size of the ellipsoids.
   * As this size is increased, the number of ellipsoids decreases.
   * @param size the maximum ellipsoid size, in samples.
   */
  public void setEllipsoidSize(int size) {
    _ellipsoidSize = size;
    dirtyDraw();
  }

  public void draw(DrawContext dc) {
    // Tensor sampling.
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    double d3 = _s3.getDelta();
    double d2 = _s2.getDelta();
    double d1 = _s1.getDelta();
    double f3 = _s3.getFirst();
    double f2 = _s2.getFirst();
    double f1 = _s1.getFirst();

    // Min/max (x,y,z) coordinates.
    double min1 = _s1.getFirst();
    double max1 = _s1.getLast();
    double min2 = _s2.getFirst();
    double max2 = _s2.getLast();
    double min3 = _s3.getFirst();
    double max3 = _s3.getLast();

    // Maximum length of eigenvectors u, v and w.
    float dmax = 0.5f*_ellipsoidSize;
    float d1max = (float)d1*dmax;
    float d2max = (float)d2*dmax;
    float d3max = (float)d3*dmax;

    // Distance between ellipsoid centers (in samples).
    int kec = (int)(2.0*dmax);

    // Scaling factor for the eigenvectors.
    float scale = dmax/sqrt(_emax);

    // Smallest eigenvalue permitted.
    float etiny = 0.0001f*_emax;

    for (int i3=0; i3<n3; i3+=3) {
    for (int i2=0; i2<n2; i2+=3) {
    //for (int i1=20; i1<n1-20; i1+=5) {

      int i1=round(_hz[i3][i2]);
      float epi = _ep[i3][i2][i1];
      //if(epi>0.2f) {
        float[] e = _et.getEigenvalues(i1,i2,i3);
        float[] u = _et.getEigenvectorU(i1,i2,i3);
        float[] v = _et.getEigenvectorV(i1,i2,i3);
        float[] w = _et.getEigenvectorW(i1,i2,i3);
        float eu = e[0], ev = e[1], ew = e[2];
        if (eu<=etiny) eu = etiny;
        if (ev<=etiny) ev = etiny;
        if (ew<=etiny) ew = etiny;
        float u1 = u[0], u2 = u[1], u3 = u[2];
        float v1 = v[0], v2 = v[1], v3 = v[2];
        float w1 = w[0], w2 = w[1], w3 = w[2];
        float su = scale*sqrt(eu);
        float sv = scale*sqrt(ev);
        float sw = scale*sqrt(ew);
        float c1 = (float)(i1*d1+f1);
        float c2 = (float)(i2*d2+f2);
        float c3 = (float)(i3*d3+f3);
        u3 *= su*d3; u2 *= su*d2; u1 *= su*d1;
        v3 *= sv*d3; v2 *= sv*d2; v1 *= sv*d1;
        w3 *= sw*d3; w2 *= sw*d2; w1 *= sw*d1;
        _ellipsoid.draw(c3,c2,c1,u3,u2,u1,v3,v2,v1,w3,w2,w1);
      //}
    }}
  }
  public BoundingSphere computeBoundingSphere(boolean finite) {
    return _bs;
  }

    /**
   * Finds the largest eigenvalue to be used for scaling.
   */
  private float findMaxEigenvalue() {
    int n1 = _et.getN1();
    int n2 = _et.getN2();
    int n3 = _et.getN3();
    float[] e = new float[3];
    float emax = 0.0f;
    for (int i3=0; i3<n3; i3+=3) {
    for (int i2=0; i2<n2; i2+=3) {
      int i1=round(_hz[i3][i2]);
      float epi = _ep[i3][i2][i1];
      //if(epi>0.2f) {
      _et.getEigenvalues(i1,i2,i3,e);
      float emaxi = max(e[0],e[1],e[2]);
      if (emax<emaxi)
        emax = emaxi;
      //}
    }}
    return emax;
  }

  private float _emax;
  private float[][][] _ep;
  private EigenTensors3 _et;
  private Sampling _s1,_s2,_s3;
  private float _ellipsoidSize = 5;
  private float[][] _hz;
  private EllipsoidGlyph _ellipsoid = new EllipsoidGlyph();
  private BoundingSphere _bs = 
    new BoundingSphere(new BoundingBox(-1,-1,-1,1,1,1));
}
