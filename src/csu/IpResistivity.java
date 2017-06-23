/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package csu;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class IpResistivity {

  public float[] stack(int nc, int np, float[] Ix) {
    float[] Is = new float[np];
    float scale = 1.0f/nc;
    for (int ic=0; ic<nc; ++ic) {
      int ns = ic*np;
      for (int ip=0; ip<np; ++ip) {
        Is[ip] += Ix[ns+ip]*scale;
      }
    }
    return Is;
  }

  public float[][] applyForComplex(float[] fx) {
    int nx = fx.length;
    Fft fft = new Fft(nx);
    Sampling sk = fft.getFrequencySampling1();
    int nk = sk.getCount();
    float[] gx = fft.applyForward(fx);
    float[][] gc = new float[3][nk];
    System.out.println("sf="+sk.getFirst());
    System.out.println("sl="+sk.getLast());
    System.out.println("sd="+sk.getDelta());
    System.out.println("nk="+sk.getCount());
    for (int kk=0, kr=0, ki=kr+1; kk<nk; ++kk, kr+=2, ki+=2) {
      float gr = gx[kr];
      float gi = gx[ki];
      gc[0][kk] = gr;
      gc[1][kk] = gi;
      gc[2][kk] = sqrt(gr*gr+gi*gi);
    }
    return gc;
  }

  public float[][] applyForResistivityAndPhase(float[][] Ic, float[][] Uc) {
    int n1 = Ic[0].length;
    float[][] rp = new float[2][n1];
    for (int i1=0; i1<n1; ++i1) {
      float ir = Ic[0][i1];
      float ii = Ic[1][i1];
      float ur = Uc[0][i1];
      float ui = Uc[1][i1];
      float is = ir*ir+ii*ii;
      float rr = ur*ir+ui*ii;
      float ri = ui*ir-ur*ii;
      rr /= is;
      ri /= is;
      rp[0][i1] = sqrt(rr*rr+ri*ri);
      rp[1][i1] = atan(ri/rr);
    }
    return rp;
  }

}
