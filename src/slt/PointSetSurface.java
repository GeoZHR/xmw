/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package slt;

import ipf.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Surface reconstruction from oriented points. 
 * <p>
 * Based on the works by Gael Guennebaud and Markus Gross, 2007, 
 * Algebraic Point Set Surfaces; 
 * Gael Guennebaud, Marcel Germann and Markus Gross, 2008, 
 * Dynamic Sampling and Rendering of Algebraic Point Set Surfaces.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.13
 */
public class PointSetSurface {

  /**
   * Given oriented points, a signed function or scalar field 
   * is computed, and the zero-contour surface of this scalar field is the 
   * surface fitting the positions and normals of the oriented points.
   * @param n1 1st dimension of the scalar field to be computed
   * @param n2 2nd dimension of the scalar field to be computed
   * @param n3 3rd dimension of the scalar field to be computed
   * @param fc an array of points with positions and orientations
   */ 
  public float[][][] findScalarField(
    final int n1, final int n2, final int n3, float[][] xus) 
  {
    final float v = -30.f;
    final float[] x1s = xus[0];
    final float[] x2s = xus[1];
    final float[] x3s = xus[2];
    final float[] u1s = xus[3];
    final float[] u2s = xus[4];
    final float[] u3s = xus[5];
    final float[] eps = xus[6];
    final KdTree kt = new KdTree(new float[][]{x1s,x2s,x3s});
    final float[][][] sf = fillfloat(v,n1,n2,n3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xms = new float[3];
      float[] xps = new float[3];
      trace("i3="+i3);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        int d = 10;
        xms[0] = i1-d; xps[0] = i1+d;
        xms[1] = i2-d; xps[1] = i2+d;
        xms[2] = i3-d; xps[2] = i3+d;
        int[] ids = kt.findInRange(xms,xps);
        int nd = ids.length;
        while(nd<20 && d<=20) {
          d++;
          xms[0]--; xps[0]++;
          xms[1]--; xps[1]++;
          xms[2]--; xps[2]++;
          ids = kt.findInRange(xms,xps);
          nd = ids.length;
        }
        if(nd<10){continue;}
        float[] dxi = new float[nd];
        float[] x1i = new float[nd];
        float[] x2i = new float[nd];
        float[] x3i = new float[nd];
        float[] u1i = new float[nd];
        float[] u2i = new float[nd];
        float[] u3i = new float[nd];
        float[] epi = new float[nd];
        for (int ik=0; ik<nd; ++ik) {
          int ip = ids[ik];
          float x1 = x1s[ip];
          float x2 = x2s[ip];
          float x3 = x3s[ip];
          float d1 = x1-i1;
          float d2 = x2-i2;
          float d3 = x3-i3;
          x1i[ik] = x1;
          x2i[ik] = x2;
          x3i[ik] = x3;
          u1i[ik] = u1s[ip];
          u2i[ik] = u2s[ip];
          u3i[ik] = u3s[ip];
          epi[ik] = eps[ip];
          dxi[ik] = d1*d1+d2*d2+d3*d3;
        }
        int is = i1*i1+i2*i2+i3*i3;
        float[] ps = sphereFit(_beta,epi,dxi,x1i,x2i,x3i,u1i,u2i,u3i);
        sf[i3][i2][i1] = ps[0]+i1*ps[1]+i2*ps[2]+i3*ps[3]+is*ps[4];
      }}
    }});
    return sf;
  }

  public float[][][] findScalarField(
    int n1, int n2, int n3, FaultCell[] cells) 
  {
    int nc = cells.length;
    float[][] xus = new float[7][nc];
    for (int ic=0; ic<nc; ++ic) {
      xus[0][ic] = cells[ic].getX1();
      xus[1][ic] = cells[ic].getX2();
      xus[2][ic] = cells[ic].getX3();
      xus[3][ic] = cells[ic].getW1();
      xus[4][ic] = cells[ic].getW2();
      xus[5][ic] = cells[ic].getW3();
      xus[6][ic] = cells[ic].getFl();
    }
    return findScalarField(n1,n2,n3,xus);
  }


  // bigger bt will produce smoother surfaces
  private static float[] sphereFit(
    float bt, float[] ep, float[] dx,
    float[] x1, float[] x2, float[] x3, 
    float[] u1, float[] u2, float[] u3)
  {
    int np = x1.length;
    float[] wp = new float[np];
    float[] wn = new float[np];
    weightFunc(ep,dx,wp,wn);
    float c1 = 0.0f;
    float d1 = 0.0f;
    float b1 = 0.0f;
    float[] wpu = new float[3];
    float[] wnu = new float[3];
    float[] wpx = new float[3];
    float[] wnx = new float[3];
    for (int ip =0; ip<np; ++ip) {
      float wpi = wp[ip];
      if(wpi==0f){continue;}
      float wni = wn[ip];
      float xi1 = x1[ip];
      float xi2 = x2[ip];
      float xi3 = x3[ip];
      float ui1 = u1[ip];
      float ui2 = u2[ip];
      float ui3 = u3[ip];
      float xsi = xi1*xi1+xi2*xi2+xi3*xi3;
      b1 += wni*xsi;
      d1 += wpi*xsi;
      c1 += wpi*(xi1*ui1+xi2*ui2+xi3*ui3); 
      wnx[0] += wni*xi1;
      wnx[1] += wni*xi2;
      wnx[2] += wni*xi3;

      wpx[0] += wpi*xi1;
      wpx[1] += wpi*xi2;
      wpx[2] += wpi*xi3;

      wnu[0] += wni*ui1;
      wnu[1] += wni*ui2;
      wnu[2] += wni*ui3;

      wpu[0] += wpi*ui1;
      wpu[1] += wpi*ui2;
      wpu[2] += wpi*ui3;
    }
    float c2 = wnx[0]*wpu[0]+wnx[1]*wpu[1]+wnx[2]*wpu[2];
    float d2 = wnx[0]*wpx[0]+wnx[1]*wpx[1]+wnx[2]*wpx[2];
    float p4 = bt*(c1-c2)/(d1-d2);
    float p1 = wnu[0]-p4*wnx[0];
    float p2 = wnu[1]-p4*wnx[1];
    float p3 = wnu[2]-p4*wnx[2];
    p4 *= 0.5f;
    float p0 = -p1*wnx[0]-p2*wnx[1]-p3*wnx[2]-p4*b1;
    return new float[]{p0,p1,p2,p3,p4};
  }

  private static void weightFunc(
    float[] ep, float[] dx, float[] wp, float[] wn) 
  {
    float ws = 0.0f;
    int np = dx.length;
    // compute adaptve point density
    // Similar to Pauly et al 2003, 
    // Shape modeling with poing-sampled geometry.
    int k = np/2;
    float[] dt = copy(dx);
    quickPartialSort(k,dt);
    float hx = dt[k]*3f;
    // compute weights
    for (int ip=0; ip<np; ++ip) {
      float x = dx[ip]/hx;
      float wpi = 0.0001f;
      if (x<1f) {wpi = pow((1.0f-x),4)*ep[ip];}
      wp[ip] = wpi;
      ws += wpi;
    }
    for (int ip=0; ip<np; ++ip)
      wn[ip] = wp[ip]/ws;
  }

  private static void trace(String s) {
    System.out.println(s);
  }


  private float _beta = 1.0f;

}
