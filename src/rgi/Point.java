/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package rgi;

import java.io.Serializable;

import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * A point is a kown sample with a sample value and its corresponding 
 * coordinates in the original, unfaulted, and flattened spaces.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.07.26
 */

public class Point implements Serializable {
  private static final long serialVersionUID = 1L;

  public float[] getSamplesX() {
    return new float[]{this.fv,this.x1,this.x2,this.x3};
  }

  public float[] getSamplesW() {
    return new float[]{this.fv,this.w1,this.w2,this.w3};
  }

  public float[] getSamplesU() {
    return new float[]{this.fv,this.u1,this.u2,this.u3};
  }

  public void setCoordW (float w1, float w2, float w3) {
    this.w1 = w1;
    this.w2 = w2;
    this.w3 = w3;
  }

  public void setCoordU (float u1, float u2, float u3) {
    this.u1 = u1;
    this.u2 = u2;
    this.u3 = u3;
  }

  public void setU1(float u1) {
    this.u1 = u1;
  }
 
  /////////////////////////////////////////////////////////////////////////
  // package

  float fv;       // sample value
  float x1,x2,x3; // original coordinates
  float w1,w2,w3; // unfault coordinates
  float u1,u2,u3; // flattened coordinates
  boolean onFault;

  Point(float fv, float x1, float x2, float x3, boolean onFault) {
    set(fv,x1,x2,x3,onFault);
  }

  Point(float fv, float x1, float x2, float x3,
        float w1, float w2, float w3, boolean onFault) {
    set(fv,x1,x2,x3,w1,w2,w3,onFault);
  }


  /////////////////////////////////////////////////////////////////////////
  // private

  private void set(float fv, 
    float x1, float x2, float x3, boolean onFault) {
    this.x1 = x1; 
    this.x2 = x2; 
    this.x3 = x3;
    this.fv = fv; 
    this.onFault = onFault;
  }

  private void set(float fv,float x1, float x2, float x3,
    float w1, float w2, float w3, boolean onFault) {
    this.x1 = x1; 
    this.x2 = x2; 
    this.x3 = x3;
    this.w1 = w1; 
    this.w2 = w2; 
    this.w3 = w3;
    this.fv = fv; 
    this.onFault = onFault;
  }


  private static class FloatList {
    public int n = 0;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public void add(float[] f) {
      int m = f.length;
      for (int i=0; i<m; ++i)
        add(f[i]);
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }
}
