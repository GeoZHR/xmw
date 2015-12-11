/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ad;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import ipfx.*;
import ipfx.FaultCell;
import static ipfx.FaultGeometry.*;

/**
 * Automatic drawing. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.04
 */
public class AutoDraw {


  public float[][] apply(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] v1 = new float[n2][n1];
    float[][] v2 = new float[n2][n1];
    float dt = 2.5f;
    float[][] fs = copy(fx);
    float sigma = 30f;
    float sigs = 1f/(sigma*sigma);
    updateDirections(fs,v1,v2);
    float[][] fp = copy(fs);
    float[][] ft = copy(fs);
    for (int iter=1; iter<100; ++iter) {
      for (int i2=1,i2m=0,i2p=2; i2<n2-1; ++i2,++i2m,++i2p) {
      for (int i1=1,i1m=0,i1p=2; i1<n1-1; ++i1,++i1m,++i1p) {
        g1[i2][i1] = (fs[i2][i1p]-fs[i2][i1m])*0.5f;
        g2[i2][i1] = (fs[i2p][i1]-fs[i2m][i1])*0.5f;
      }}
      for (int i2=1,i2m=0,i2p=2; i2<n2-1; ++i2,++i2m,++i2p) {
      for (int i1=1,i1m=0,i1p=2; i1<n1-1; ++i1,++i1m,++i1p) {
        float v1i = v1[i2][i1];
        float v2i = v2[i2][i1];
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        float eta = abs(v1i*g1i+v2i*g2i);
        eta = pow(eta,4f)*exp(-iter*iter*sigs);
        float fxii = ft[i2][i1];
        float fxpi = ft[i2][i1p];
        float fxmi = ft[i2][i1m];
        float fxip = ft[i2p][i1];
        float fxim = ft[i2m][i1];
        float fxd1 = max(fxpi-fxii,fxmi-fxii);
        float fxd2 = max(fxip-fxii,fxim-fxii);
        fxd1 = max(0,fxd1);
        fxd2 = max(0,fxd2);
        float fxdi = eta*dt*sqrt(fxd1*fxd1+fxd2*fxd2);
        fs[i2][i1] += fxdi; 
      }}
      if(iter>1) ft = sub(fs,fp);
      fp = copy(fs);
      updateDirections(fs,v1,v2);
    }
    return fs;
  }

  private void updateDirections(
    float[][] fx, float[][] v1, float[][] v2)
  {
    LocalOrientFilter lof = new LocalOrientFilter(4.0,4.0);
    lof.apply(fx,null,null,null,v1,v2,null,null,null);
  }
}
