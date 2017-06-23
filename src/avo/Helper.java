/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package avo;

import java.io.*;
import java.util.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Helper {

  public float[][][] sortLogs(float[][] rv1, float[][] rv2, 
    float[][] rv3, float[][] rv4) {
    int[] nps = new int[4];
    int np1 = rv1[0].length;
    int np2 = rv2[0].length;
    int np3 = rv3[0].length;
    int np4 = rv4[0].length;
    nps[0] = np1;
    nps[1] = np2;
    nps[2] = np3;
    nps[3] = np4;
    int npm = max(nps);
    //float[][][] las = new float[3][4][npm];
    float[][][] las = fillfloat(-999.25f,npm,4,3);
    for (int ip=0; ip<np1; ++ip) { 
      las[0][0][ip] = rv1[0][ip];
      las[1][0][ip] = rv1[1][ip];
      las[2][0][ip] = rv1[2][ip];
    }
    for (int ip=0; ip<np2; ++ip) { 
      las[0][1][ip] = rv2[0][ip];
      las[1][1][ip] = rv2[1][ip];
      las[2][1][ip] = rv2[2][ip];
    }
    for (int ip=0; ip<np3; ++ip) { 
      las[0][2][ip] = rv3[0][ip];
      las[1][2][ip] = rv3[1][ip];
      las[2][2][ip] = rv3[2][ip];
    }
    for (int ip=0; ip<np4; ++ip) { 
      las[0][3][ip] = rv4[0][ip];
      las[1][3][ip] = rv4[1][ip];
      las[2][3][ip] = rv4[2][ip];
    }
    return las;
  }

  public float[][][] resampleLogs(
    int m1, int d1, float[][][] las) 
  {
    int n3 = las.length;
    int n2 = las[0].length;
    float[][][] lar = new float[n3][n2][m1];
    MedianFinder mf = new MedianFinder(7);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<m1; ++i1) {
      int b1 = i1*d1;
      float[] rht = copy(d1,b1,las[0][i2]);
      float[] vpt = copy(d1,b1,las[1][i2]);
      float[] vst = copy(d1,b1,las[2][i2]);
      float vsm = mf.findMedian(vst);
      int[] id = new int[1];
      float[] vsa = abs(sub(vst,vsm));
      min(vsa,id);
      lar[0][i2][i1] = rht[id[0]];
      lar[1][i2][i1] = vpt[id[0]];
      lar[2][i2][i1] = vst[id[0]];
    }}
    return lar;
  }
}
