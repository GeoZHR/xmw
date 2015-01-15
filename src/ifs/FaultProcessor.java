/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;


import static edu.mines.jtk.util.ArrayMath.*;

import static ifs.FaultScanner.*;

/**
  *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.12.29
 */
public class FaultProcessor {

  public FaultSkin[] computeFaultSkins(
    int minSkinSize,int dh, int dv, float[][][] f) {
    int nx = f.length;
    int ny = f[0].length;
    int nz = f[0][0].length;
    if(dh==1 && dv==1) {
      System.out.println("Use full grid method");
      FaultSkin[] sks = daveSkinning(f);
      FaultCell[] fcs = FaultSkin.getCells(sks);
      FaultSurfer fsr = new FaultSurfer(nz,ny,nx,fcs);
      return fsr.applySurferM(minSkinSize);
    } else {
      /*
      System.out.println("Use down grid method");
      float[][][] fh = sparseSample(dh,dv,f);
      int nxd = fh.length;
      int nyd = fh[0].length;
      int nzd = fh[0][0].length;
      Sampling sxu = new Sampling(nx);
      Sampling syu = new Sampling(ny);
      Sampling szu = new Sampling(nz);
      Sampling sxd = new Sampling(nxd,dh,0);
      Sampling syd = new Sampling(nyd,dh,0);
      Sampling szd = new Sampling(nzd,dv,0);
      FaultSkin[] sks = daveSkinning(fh);
      FaultCell[] fcs = FaultSkin.getCells(sks);
      FaultSurferM sfc = new SkinsFromCells(szd,syd,sxd,szu,syu,sxu,fcs);
      return sfc.applyForSkins();
      */
      return null;
    }
  }

  private FaultSkin[] daveSkinning(float[][][] f) {
    int minP = 0;
    int maxP = 360;
    int minT = 65;
    int maxT = 85;
    double sigmaP = 4.0;
    double sigmaT = 20.0;
    FaultScanner fsc = new FaultScanner(sigmaP,sigmaT);
    float[][][] ft = taper(10,0,0,f);
    float[][][][] ps = computeSlopes(ft);
    float[][][][] flpt = fsc.scan(minP,maxP,minT,maxT,ps[0],ps[1],ft,null);
    int n3 = flpt[0].length;
    int n2 = flpt[0][0].length;
    int n1 = flpt[0][0][0].length;
    int minSkinSize = round(max(n2,n3)*n1*0.2f);
    double lowerFl = 0.2;
    double upperFl = 0.5;
    FaultSkinner fsk = new FaultSkinner();
    fsk.setMinSkinSize(minSkinSize);
    fsk.setGrowLikelihoods(lowerFl,upperFl);
    FaultCell[] fcs = fsk.findCells(flpt);
    FaultSkin[] fss = fsk.findSkins(fcs);
    return fss;

  }

  private float[][][][] computeSlopes(float[][][] f) {
    double pmax = 5.0f;
    double sigma3 = 1.0f;
    double sigma2 = 1.0f;
    double sigma1 = 16.f;
    return slopes(sigma1,sigma2,sigma3,pmax,f);
  }

  private float[][][] sparseSample(int dh, int dv, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    int n3h = reSampling(dh,n3);
    int n2h = reSampling(dh,n2);
    int n1h = reSampling(dv,n1);
    float[][][] y = new float[n3h][n2h][n1h];
    for (int i3=0; i3<n3h; i3++) 
    for (int i2=0; i2<n2h; i2++) 
    for (int i1=0; i1<n1h; i1++) {
      int i1i = i1*dv;
      int i2i = i2*dh;
      int i3i = i3*dh;
      y[i3][i2][i1] = x[i3i][i2i][i1i];
    }
    return y;
  }

  private int reSampling(int d, int n) {
    int ns = 0; 
    for (int i=0; i<n; i+=d)
      ns += 1;
    return ns;
  }

}
