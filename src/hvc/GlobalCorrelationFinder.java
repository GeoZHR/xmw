/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hvc;

import java.util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates local slopes of features in 2D and 3D images.
 * 
 * For a 2D image f(x1,x2), slope p2(x1,x2) is the ratio dx1/dx2 of 
 * linear features nearest the image sample at point (x1,x2). An 
 * estimate of the linearity (a number in [0,1]) of features nearest 
 * that sample may also be computed.
 * 
 * Likewise, for a 3D image f(x1,x2,x3), slopes p2(x1,x2,x3) and
 * p3(x1,x2,x3) are the ratios dx1/dx2 and dx1/dx3, respectively, 
 * of planar features nearest the image sample at point (x1,x2,x3). 
 * An estimate of the planarity (a number in [0,1]) of features 
 * nearest that sample may also be computed.
 *
 * All slopes are measured in unitless samples (dx1) per sample (dx2).
 * Minimum and maximum slopes may be specified, and default min and 
 * max bounds are -100 and 100 samples per sample, respectively.
 * Estimated slopes will be within the default or specified bounds.
 *
 * @author Xinming Wu, University of Texas at Austin.
 * @version 2018.08.03
 */
public class GlobalCorrelationFinder {

  /**
   * Constructs a global correlation finder with shift bounds.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public GlobalCorrelationFinder(int shiftMin, int shiftMax) {
    _shiftMin = shiftMin;
    _shiftMax = shiftMax;
  }

  public void setStrainMax(double strainMax) {
    _strainMax = strainMax;
  }
 
  /**
   * Finds slopes of features in the specified 2D image.
   * @param f array[n2][n1] of input image samples.
   * @param p2 array[n2][n1] of output slopes.
   */
  public float[][] findCorrelations(int dm, int[] ps, float[][] fx) {
    int np = ps.length;
    int n1 = fx[0].length;
    Sampling s1 = new Sampling(n1);
    DynamicWarping dw = new DynamicWarping(_shiftMin,_shiftMax);
    dw.setStrainMax(_strainMax);
    InverseInterpolator ii = new InverseInterpolator(s1,s1);
    ArrayList<float[]> dts = new ArrayList<float[]>();
    for (int ip=0; ip<np; ++ip) {
    for (int jp=ip+1; jp<np; ++jp) {
      int p2 = ps[ip];
      int m2 = ps[jp];
      if(abs(p2-m2)>dm) continue;
      float[] fp = fx[p2];
      float[] fm = fx[m2];
      float[] dt = new float[n1+3];
      float[] t = new float[n1];
      float[] u = dw.findShifts(fp,fm);
      for (int i1=0; i1<n1; ++i1)
        u[i1] += i1;
      dt[0] = p2;
      dt[1] = m2;
      dt[2] = round(max(u[0],0));
      ii.invert(u,t);
      for (int i1=0; i1<n1; ++i1) {
        dt[i1+3] = t[i1];
      }
      dts.add(dt);
    }}
    return dts.toArray(new float[0][]);
  }

  public void unpack(float[][] dks, int[][] ks, float[][] ds) {
    int ns = dks.length;
    int n1 = ds[0].length;
    for (int is=0; is<ns; ++is) {
      ks[is][0] = (int)dks[is][0];
      ks[is][1] = (int)dks[is][1];
      ks[is][2] = (int)dks[is][2];
    }
    for (int is=0; is<ns; ++is) {
      ds[is] = copy(n1,3,dks[is]);
    }
  }

 
  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _shiftMin,_shiftMax; 
  private double _strainMax;
}
