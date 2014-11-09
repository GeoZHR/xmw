/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wsi;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Migrated shot gather image flattening before stacking.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.11.08
 */
public class FlattenAndStack {

  public void setForSlopes(float sigma1, float sigma2, float pmax) {
    _pmax = pmax;
    _pSigma1 = sigma1;
    _pSigma2 = sigma2;
  }

  public void setFlattenSmoothings(float sigma1, float sigma2) {

  }

  public float[][] apply(
    final Sampling s1, final Sampling s2, final Sampling s3, final float[][][] x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    final float[][] g = zerofloat(n1,n2);
    final float[][] s = stack(x);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] f = zerofloat(n1,n3);
      System.out.println("i2="+i2);
      for (int i3=0; i3<n3; ++i3)
        f[i3] = x[i3][i2];
      float[][] p2 = new float[n3][n1];
      float[][] wp = new float[n3][n1];
      float[][] u1 = new float[n3][n1];
      LocalSlopeFinder lsf = new LocalSlopeFinder(_pSigma1,_pSigma2,_pmax);
      Flattener2C flc = new Flattener2C();
      flc.setWeight1(_weight1);
      flc.setSmoothings(_fSigma1,_fSigma2);
      flc.setIterations(_small,_niter);
      lsf.findSlopes(f,p2,wp); 
      cleanImage(p2,f);
      //int[] c = referenceTrace(mul(wp,u1));
      int[] c = referenceTrace(s[i2],f);
      wp = pow(wp,20.0f);
      Flattener2C.Mappings mp = flc.getMappingsFromSlopes(s1,s3,p2,wp,c);
      stack(g[i2],mp.flatten(f));
    }});
    return g;
  }

  private float _pSigma1 = 4.0f;
  private float _pSigma2 = 2.0f;
  private float _pmax = 6.0f;

  private int _niter = 200;
  private float _small = 0.01f;
  private float _fSigma1 = 4.0f;
  private float _fSigma2 = 8.0f;
  private float _weight1 = 0.5f;

  private void cleanImage(float[][] p2, float[][] x) {
    int n3 = x.length;
    int n1 = x[0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        if (p2[i3][i1]>=_pmax) {
          x[i3][i1] = 0.0f;
        }
      }
    }
  }

  private float[][] stack(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][] xs = new float[n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      add(x[i3],xs,xs);
    }
    return xs;
  }

  private void stack(float[] g, float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        g[i1] += f[i2][i1];
  }

  private int[] referenceTrace(float[] xsi, float[][] x2i) {
    int[] c = new int[1];
    int n3 = x2i.length;
    float sum = 50000000000.0f;
    for (int i3=0; i3<n3; ++i3) {
      float smi = sum(abs(sub(x2i[i3],xsi)));
      if (smi<sum) {
        sum = smi;
        c[0] = i3;
      }
    }
    return c;
  }


  private int[] referenceTrace(float[][] f) {
    int n2 = f.length;
    float sum = 0.0f;
    int[] c = new int[1];
    for (int i2=0; i2<n2; ++i2) {
      float smi = sum(abs(f[i2]));
      if (smi>sum) {
        c[0] = i2;
        sum = smi;
      }
    }
    System.out.println("ci="+c[0]);
    return c;
  }

}
