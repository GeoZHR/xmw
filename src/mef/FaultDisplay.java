/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.util.Stopwatch;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import util.*;
import ipfx.FaultCell;
import static ipfx.FaultGeometry.*;

/**
 * 3D tensor voting. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.02.11
 */
public class FaultDisplay {


  public float[][][] setSlices(int k1, int k2, int k3, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] ft = copy(fl);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
      int i2m = i2-1;
      int i2p = i2+1;
      float fli = fl[i3][i2 ][k1];
      float flm = fl[i3][i2m][k1];
      float flp = fl[i3][i2p][k1];
      if(fli>0f&&flm<=0f) {ft[i3][i2m][k1]=fli;}
      if(fli>0f&&flp<=0f) {ft[i3][i2p][k1]=fli;}
    }}
    for (int i2=0; i2<n2; ++i2) {
    for (int i3=1; i3<n3-1; ++i3) {
      int i3m = i3-1;
      int i3p = i3+1;
      float fli = fl[i3 ][i2][k1];
      float flm = fl[i3m][i2][k1];
      float flp = fl[i3p][i2][k1];
      if(fli>0f&&flm<=0f) {ft[i3m][i2][k1]=fli;}
      if(fli>0f&&flp<=0f) {ft[i3p][i2][k1]=fli;}
    }}
    float[][][] fp = copy(ft);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
      int i2m = i2-1;
      float fli = ft[i3][i2 ][k1];
      float flm = ft[i3][i2m][k1];
      if(fli>0f&&flm<=0f) {fp[i3][i2m][k1]=fli;}
    }}
    for (int i1=0; i1<n1; ++i1) {
    for (int i2=1; i2<n2-1; ++i2) {
      int i2m = i2-1;
      int i2p = i2+1;
      float fli = fl[k3][i2 ][i1];
      float flm = fl[k3][i2m][i1];
      float flp = fl[k3][i2p][i1];
      if(fli>0f&&flm<=0f) {fp[k3][i2m][i1]=fli;}
      if(fli>0f&&flp<=0f) {fp[k3][i2p][i1]=fli;}
    }}
    for (int i1=0; i1<n1; ++i1) {
    for (int i3=1; i3<n3-1; ++i3) {
      int i3m = i3-1;
      int i3p = i3+1;
      float fli = fl[i3 ][k2][i1];
      float flm = fl[i3m][k2][i1];
      float flp = fl[i3p][k2][i1];
      if(fli>0.0f&&flm<=0f) {fp[i3m][k2][i1]=fli;}
      if(fli>0.0f&&flp<=0f) {fp[i3p][k2][i1]=fli;}
    }}
    return fp;
  }

  private FaultCell[][][] getCellArray(
    int n1, int n2, int n3, Sampling sp, Sampling st, FaultCell[] fcs) 
  {
    int np = sp.getCount();
    int nt = st.getCount();
    int[][] ct = new int[np][nt];
    FaultCell[][][] fca = new FaultCell[np][nt][];
    FaultCell[][][] fct = new FaultCell[np][nt][n1*max(n2,n3)*5];
    for (FaultCell fci:fcs) {
      int ip = sp.indexOfNearest(fci.getFp());
      int it = st.indexOfNearest(fci.getFt());
      int ic = ct[ip][it];
      fct[ip][it][ic] = fci;
      ct[ip][it] += 1;
    }
    for (int ip=0; ip<np; ++ip) {
    for (int it=0; it<nt; ++it) {
      int nc = ct[ip][it];
      fca[ip][it] = new FaultCell[nc];
      for (int ic=0; ic<nc; ++ic) {
        fca[ip][it][ic] = fct[ip][it][ic];
      }
    }}
    fct = null;
    return fca;
  }

  public FaultCell[] getFaultCells(int n1, int n2, int n3, FaultCell[] fcs) {
    FaultCell[][][] fcg = new FaultCell[n3][n2][n1];
    ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    for (FaultCell fc:fcs) {
      int i1 = fc.getI1();
      int i2 = fc.getI2();
      int i3 = fc.getI3();
      fcg[i3][i2][i1] = fc;
    }

    for (int i3=0; i3<n3; i3+=1) {
    for (int i2=0; i2<n2; i2+=2) {
    for (int i1=0; i1<n1; i1+=2) {
      FaultCell fc = fcg[i3][i2][i1];
      if(fc!=null) cellList.add(fc);
    }}}
    return cellList.toArray(new FaultCell[0]);
  }

  public void getFlImage(FaultCell[] fcs, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    for (FaultCell fc:fcs) {
      int i1 = fc.getI1();
      int i2 = fc.getI2();
      int i3 = fc.getI3();
      i1 = min(i1,n1-1); i1 = max(i1,0);
      i2 = min(i2,n2-1); i2 = max(i2,0);
      i3 = min(i3,n3-1); i3 = max(i3,0);
      int m1 = fc.getM1();
      int m2 = fc.getM2();
      int m3 = fc.getM3();
      int p1 = fc.getP1();
      int p2 = fc.getP2();
      int p3 = fc.getP3();
      m1 = min(m1,n1-1); m1 = max(m1,0);
      m2 = min(m2,n2-1); m2 = max(m2,0);
      m3 = min(m3,n3-1); m3 = max(m3,0);
      p1 = min(p1,n1-1); p1 = max(p1,0);
      p2 = min(p2,n2-1); p2 = max(p2,0);
      p3 = min(p3,n3-1); p3 = max(p3,0);
      fl[m3][m2][m1] = fc.getFl();
      fl[p3][p2][p1] = fc.getFl();
      fl[i3][i2][i1] = fc.getFl();
    }
  }

  public void getFaultImages(
    FaultCell[] fcs, float[][][] fl, float[][][] fp, float[][][] ft) 
  {
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    for (FaultCell fc:fcs) {
      int i1 = fc.getI1();
      int i2 = fc.getI2();
      int i3 = fc.getI3();
      i1 = min(i1,n1-1); i1 = max(i1,0);
      i2 = min(i2,n2-1); i2 = max(i2,0);
      i3 = min(i3,n3-1); i3 = max(i3,0);
      fl[i3][i2][i1] = fc.getFl();
      fp[i3][i2][i1] = fc.getFp();
      ft[i3][i2][i1] = fc.getFt();
      /*
      int m1 = fc.getM1();
      int m2 = fc.getM2();
      int m3 = fc.getM3();
      int p1 = fc.getP1();
      int p2 = fc.getP2();
      int p3 = fc.getP3();
      m1 = min(m1,n1-1); m1 = max(m1,0);
      m2 = min(m2,n2-1); m2 = max(m2,0);
      m3 = min(m3,n3-1); m3 = max(m3,0);
      p1 = min(p1,n1-1); p1 = max(p1,0);
      p2 = min(p2,n2-1); p2 = max(p2,0);
      p3 = min(p3,n3-1); p3 = max(p3,0);
      fl[m3][m2][m1] = fc.getFl();
      fl[p3][p2][p1] = fc.getFl();
      fp[m3][m2][m1] = fc.getFp();
      fp[p3][p2][p1] = fc.getFp();
      ft[m3][m2][m1] = fc.getFt();
      ft[p3][p2][p1] = fc.getFt();
      */
    }
  }
}
