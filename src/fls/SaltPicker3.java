/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.*;
import pik.*;

/**
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.16.07
 */


public class SaltPicker3 {

  public float[][][] pick3(
    int d3, float[] c3, float[][] c2, float[][] c1, float[][][] pa, 
    float[][][] fss) 
  {
    int nc = c3.length;
    int n3 = pa.length;
    int n2 = pa[0].length;
    int n1 = pa[0][0].length;
    SaltPicker2 sp2 = new SaltPicker2();
    float[][][] pks = new float[n3][4][];
    for (int ic=0; ic<nc; ++ic) {
      int k3 = round(c3[ic]);
      float[][] xu = sp2.initialBoundary(1,c1[ic],c2[ic]);
      sp2.refine(50,1,10,1,xu,pa[k3]);
      float[][] xp = copy(xu);
      float[][] xm = copy(xu);
      float[][] xt = sp2.regridBoundary(1.0f,xu[0],xu[1]);
      pks[k3] = xt;
      float[][] wx = new float[n2][n1];
      pointsToImage(xt[1],xt[0],pa[k3],wx);
      signAsignmentH(wx,fss[k3]);
      int e3 = min(k3+d3,n3-1);
      int b3 = max(k3-d3,0);
      for (int i3=k3+1; i3<=e3; i3++) {
        xp = sp2.pickNext(5,1,3,1,xp[0],xp[1],pa[i3]);
        xt = sp2.regridBoundary(1.0f,xp[0],xp[1]);
        pks[i3] = xt;
        zero(wx);
        pointsToImage(xt[1],xt[0],pa[k3],wx);
        signAsignmentH(wx,fss[i3]);

      }
      for (int i3=k3-1; i3>=b3; i3--) {
        xm = sp2.pickNext(5,1,3,1,xm[0],xm[1],pa[i3]);
        xt = sp2.regridBoundary(1.0f,xm[0],xm[1]);
        pks[i3] = xt;
        zero(wx);
        pointsToImage(xt[1],xt[0],pa[k3],wx);
        signAsignmentH(wx,fss[i3]);
      }
    }
    return pks;
  }

  public void pointsToImage(
    float[] ys, float[] zs, float[][] pa, float[][] wx) 
  {
    int n2 = pa.length;
    int n1 = pa[0].length;
    int np = ys.length;
    for (int ip=0; ip<np; ++ip) {
      int i2 = round(ys[ip]);
      int i1 = round(zs[ip]);
      i2 = max(i2,0);
      i2 = min(i2,n2-1);
      i1 = max(i1,0);
      i1 = min(i1,n1-1);
      wx[i2][i1] = pa[i2][i1];
    }
  }

  public void signAsignmentH(
    float[][] wx, float[][] fx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    for (int i1=0; i1<n1; ++i1) {
      float mk = 1f;
    for (int i2=0; i2<n2; ++i2) {
      if(wx[i2][i1]!=0f) {
        fx[i2][i1] = 0f;
        if(i2==0) mk*=-1f;
        if(i2>0&&wx[i2-1][i1]==0) mk *= -1f;
      } else {
        fx[i2][i1] = mk;
      }
    }}
    //despike
    float[][] fs = zerofloat(n1,n2);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(10);
    ref.apply(fx,fs);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(fx[i2][i1]!=0) {
        if(fs[i2][i1]>0) fx[i2][i1] =  1f;
        if(fs[i2][i1]<0) fx[i2][i1] = -1f;
      }
    }}
  }


}
