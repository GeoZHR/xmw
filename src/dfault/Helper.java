/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dfault;

import util.*;
import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Generates 2D fake training data for DeepFault.
 * <em>
 * Jacobians of functions used in folding and faulting have been implemented
 * but not tested. Therefore, beware of errors in calculated slopes p2 and p3.
 * </em>
 * @author Xinming Wu, University of Texas at Austin.
 * @version 2017.08.01
 */
public class Helper {

  public float[][] thin(float[][] ep) {
    int n2 = ep.length;
    int n1 = ep[0].length;
    float[][] ept = new float[n2][n1];
    float[][] eps = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
    rgf.applyX0(ep,eps);
    for (int i2=1; i2<n2-1; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float epi = eps[i2  ][i1];
      float epm = eps[i2-1][i1];
      float epp = eps[i2+1][i1];
      if(epi>epm&&epi>epp) 
        ept[i2][i1] = ep[i2][i1];
    }}
    return ept;
  }

  public float[][] pickSeeds(
    int d, float fm, float[][] ft) {
    final int n2 = ft.length;
    final int n1 = ft[0].length;
    final ArrayList<int[]> cs = new ArrayList<int[]>();
    final ArrayList<Float> fa = new ArrayList<Float>();
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float fti = ft[i2][i1];
      if(fti>fm) {
        fa.add(fti);
        cs.add(new int[]{i1,i2});
      }
    }}
    int np = cs.size();
    int[] is = new int[np];
    float[] fs = new float[np];
    for (int ip=0; ip<np; ++ip) {
      is[ip] = ip;
      fs[ip] = fa.get(ip);
    }
    quickIndexSort(fs,is);
    int[][] mark = new int[n2][n1];
    ArrayList<float[]> seeds = new ArrayList<float[]>();
    for (int ip=np-1; ip>=0; --ip) {
      int[] id = cs.get(is[ip]);
      int i1 = id[0];
      int i2 = id[1];
      int b1 = i1-d; b1=max(b1,0);
      int b2 = i2-d; b2=max(b2,0);
      int e1 = i1+d; e1=min(e1,n1-1);
      int e2 = i2+d; e2=min(e2,n2-1);
      boolean ok = true;
      for (int k2=b2;k2<=e2;k2++) {
      for (int k1=b1;k1<=e1;k1++) {
        if(mark[k2][k1]==1) {
          ok=false;
          break;
        }
      }}
      if(ok) {
        seeds.add(new float[]{i1,i2});
        mark[i2][i1] = 1;
      }
    }
    return seeds.toArray(new float[0][]);
  }

  public float[][][] setTestImage(int m1, int m2, float[][] seeds, float[][] gx) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    int h1 = round(m1*0.5f);
    int h2 = round(m2*0.5f);
    int np = seeds.length;
    float[][][] gs = new float[np][m2][m1];
    for (int ip=0; ip<np; ++ip) {
      int i1 = (int)seeds[ip][0];
      int i2 = (int)seeds[ip][1];
      for(int k2=-h2; k2<h2; k2++) {
      for(int k1=-h1; k1<h1; k1++) {
        int j1 = i1+k1;
        int j2 = i2+k2;
        j1 = min(max(j1,0),n1-1);
        j2 = min(max(j2,0),n2-1);
        gs[ip][k2+h2][k1+h1] = gx[j2][j1];
      }}
    }
    return gs;
  }


  public float[][][] setTestImage(int m1, int m2, float[][] gx) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    int h1 = round(m1*0.5f);
    int h2 = round(m2*0.5f);
    float[][][] gs = new float[n1*n2][m2][m1];
    int ct = 0;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      for(int k2=-h2; k2<h2; k2++) {
      for(int k1=-h1; k1<h1; k1++) {
        int j1 = i1+k1;
        int j2 = i2+k2;
        j1 = min(max(j1,0),n1-1);
        j2 = min(max(j2,0),n2-1);
        gs[ct][k2+h2][k1+h1] = gx[j2][j1];
      }}
      ct++;
    }}
    return gs;
  }

}
