/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.09.16
 */
public class FaultSurfer {

  public void setTransform(boolean X2ToX1, boolean X3ToX1) {
    _X2ToX1 = X2ToX1;
    _X3ToX1 = X3ToX1;
  }

  public float[][][][] applyTransform(FaultCell[] fc, float[][][] g){
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    float[][][][] ts = new float[4][n3][n1][n2];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          ts[0][i3][i1][i2] = g[i3][i2][i1];
        }
      }
    }
    faultSlopeInterp(20.0f,fc,ts[1],ts[2],ts[3]);
    return ts;
  }

  public float[][] surfer(
    int n1, int n2, int n3, FaultSkin[] fs) 
  {
      float[][] kk = new float[3][1];
      FaultSkin[] fst = new FaultSkin[1];
      fst[0] = fs[0];
      kk[0][0] = FaultSkin.getCells(fst)[50].i2;
      kk[1][0] = FaultSkin.getCells(fst)[50].i1;
      kk[2][0] = FaultSkin.getCells(fst)[50].i2;
      FaultCell[] fc = FaultSkin.getCells(fs);
      System.out.println("k1="+kk[0][0]);
      float[][]   sf = new float[n3][n1];
      float[][][] p2 = new float[n3][n1][n2];
      float[][][] p3 = new float[n3][n1][n2];
      float[][][] wp = new float[n3][n1][n2];
      faultSlopeInterp(20.0f,fc,p2,p3,wp);
      surfer(kk,p2,p3,wp,sf);
      return sf;
  }

  private void surfer(float[][] ks, float[][][] p2, 
    float[][][] p3, float[][][] wp, float[][] sf) 
  {
    int n3 = p2.length;
    int n2 = p2[0].length;
    int n1 = p2[0][0].length;
    SurfaceExtractorC se = new SurfaceExtractorC();
    add(sf,ks[0][0],sf);
    //sf = se.surfaceInitialization(n2,n3,n1,ks[0],ks[1],ks[2]);
    se.surfaceUpdateFromSlopes(wp,p2,p3,ks[0],ks[1],ks[2],sf);
  }

  public void faultSlopeInterp(float pMax, FaultCell[] fc, 
    float[][][] p2, float[][][] p3, float[][][] wp) {
    float pMin = -pMax;
    int nc = fc.length; 
    int n3 = p2.length; 
    int n2 = p2[0].length; 
    int n1 = p2[0][0].length; 
    float[][] xf = new float[3][nc];
    if (_X2ToX1) {
      for (int ic=0; ic<nc; ++ic) {
        xf[0][ic] = fc[ic].i2;
        xf[1][ic] = fc[ic].i1;
        xf[2][ic] = fc[ic].i3;
      }
    }
    KdTree kt = new KdTree(xf);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          int ic = kt.findNearest(new float[]{i1,i2,i3});
          float u1 = fc[ic].w2;
          float u2 = fc[ic].w1;
          float u3 = fc[ic].w3;
          if (u1<0.0f) {
            u1 = -u1;
            u2 = -u2;
            u3 = -u3;
          }
          if (-u2<pMin*u1) {u2 = -pMin*u1;}
          if (-u2>pMax*u1) {u2 = -pMax*u1;}
          if (-u3<pMin*u1) {u3 = -pMin*u1;}
          if (-u3>pMax*u1) {u3 = -pMax*u1;}
          if(u1==0.0f){
            p2[i3][i2][i1] = (u2<0.0f)?pMax:pMin; 
            p3[i3][i2][i1] = (u3<0.0f)?pMax:pMin; 
          } else {
            p2[i3][i2][i1] = -u2/u1;
            p3[i3][i2][i1] = -u3/u1;
          }
          float d1 = i1-fc[ic].i1;
          float d2 = i2-fc[ic].i2;
          float d3 = i3-fc[ic].i3;
          float ds = sqrt(d1*d1+d2*d2+d3*d3);
          wp[i3][i2][i1] = fc[ic].fl;
          if(ds==0.0f) {
            wp[i3][i2][i1] = fc[ic].fl;
          } else {
            wp[i3][i2][i1] = fc[ic].fl/ds;
          }
        }
      }
    }
  }

  private boolean _X2ToX1 = true;
  private boolean _X3ToX1 = false;
}
