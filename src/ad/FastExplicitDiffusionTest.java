/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ad;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import util.*;

/**
 * Fast explicit diffusion filter. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.05
 */

public class FastExplicitDiffusionTest {

  /**
   * Construction.
   * @param m number of cycles
   * @param d stability limit
   */
  public void setCycles(int m, float d) {
    _m = m;
    _d = d;
  }

  public float[][][] apply(
    float sigma, EigenTensors et,float[][][] fx) 
  {
    float[][][] gx = copy(fx);
    float t = sigma*sigma*0.5f;
    FedStep fs = new FedStep(t,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    Stopwatch sw = new Stopwatch();
    sw.start();
    int ik = 0;
    int nk = _m*nc;
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      if (ik>0) {
        double timeUsed = sw.time();
        double timeLeft = ((double)nk/(double)ik-1.0)*timeUsed;
        int timeLeftSec = 1+(int)timeLeft;
        trace("Linear diffusion: done in "+timeLeftSec+" seconds");
      }
      applyLaplacian(et,-ts[ic],copy(gx),gx);
      ik++;
    }}
    sw.stop();
    trace("Linear diffusion: done");
    return gx;
  }

  private void applyLaplacian(final EigenTensors d, final float s, 
    final float[][][] x, final float[][][] y) 
  {
    int i3start = 1; 
    final int i3step = 2; 
    final int i3stop = x.length;
    for (int i3pass=0; i3pass<i3step; ++i3pass,++i3start) {
      Parallel.loop(i3start,i3stop,i3step,new Parallel.LoopInt() {
        public void compute(int i3) {apply22(i3,d,s,x,y);}});
    }
  }

  private void apply22(int i3, EigenTensors d, float s, 
    float[][][] x, float[][][] y) 
  {
    s *= 0.0625f;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float[] di = new float[6];
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x0m = x[i3  ][i2-1];
      float[] xm0 = x[i3-1][i2  ];
      float[] xmm = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y0m = y[i3  ][i2-1];
      float[] ym0 = y[i3-1][i2  ];
      float[] ymm = y[i3-1][i2-1];
      for (int i1=1,m1=0; i1<n1; ++i1,++m1) {
        d.getTensor(i1,i2,i3,di);
        float d11 = di[0];
        float d12 = di[1];
        float d13 = di[2];
        float d22 = di[3];
        float d23 = di[4];
        float d33 = di[5];
        float xa = x00[i1]-xmm[m1];
        float xb = x00[m1]-xmm[i1];
        float xc = x0m[i1]-xm0[m1];
        float xd = xm0[i1]-x0m[m1];
        float x1 = xa-xb+xc+xd;
        float x2 = xa+xb-xc+xd;
        float x3 = xa+xb+xc-xd;
        float y1 = (d11*x1+d12*x2+d13*x3)*s;
        float y2 = (d12*x1+d22*x2+d23*x3)*s;
        float y3 = (d13*x1+d23*x2+d33*x3)*s;
        float ya = y1+y2+y3; y00[i1] += ya; ymm[m1] -= ya;
        float yb = y1-y2+y3; y0m[i1] += yb; ym0[m1] -= yb;
        float yc = y1+y2-y3; ym0[i1] += yc; y0m[m1] -= yc;
        float yd = y1-y2-y3; ymm[i1] += yd; y00[m1] -= yd;
      }
    }
  }


    // Copys array x to array y.
  private static void scopy(float[] x, float[] y) {
    copy(x,y);
  }
  private static void scopy(float[][] x, float[][] y) {
    copy(x,y);
  }
  private static void scopy(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        scopy(x[i3],y[i3]);
      }
    });
  }



  private static void trace(String s) {
    System.out.println(s);
  }

  private int _m = 5; //number of cycles
  private float _d = 0.5f; //stability limit
}
