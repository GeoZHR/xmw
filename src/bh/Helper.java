/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package bh;

import static edu.mines.jtk.util.ArrayMath.*;

public class Helper {

  public float[][][] getTopBottomHorizons(float[][][] dx) {
    int n3 = dx.length;
    int n2 = dx[0].length;
    int n1 = dx[0][0].length;
    float[][][] hs = new float[2][n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      float dxm = dx[i3][i2][i1-1];
      float dxi = dx[i3][i2][i1  ];
      if(dxm==0f && dxi>0f)
        hs[0][i3][i2] = i1;
      if(dxm>0f && dxi==0f)
        hs[1][i3][i2] = i1;
    }}}
    return hs;
  }

  public float[][][] flattenByHorizon(float[][] hx, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    float[][][] gf = new float[n3][n2][n1];
    int hm = round(max(hx));
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int k1 = round(i1+hx[i3][i2]-hm);
      k1 = max(k1,0);
      k1 = min(k1,n1-1);
      gf[i3][i2][i1] = gx[i3][i2][k1];

    }}}
    return gf;
  }

}
