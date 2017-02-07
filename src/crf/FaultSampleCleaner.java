package crf;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

import mef.*;


/**
 * Recompute fault likelihoods at fault samples from a smoothed 
 * seismic image, and discard fault samples with relatively small 
 * fault likelihoods.
 * The seismic image is smoothed by structure-oriented smoothing with 
 * fault preserving.
 * @author Xinming Wu, BEG, University of Texas at Austin.
 * @version 2017.02.06
 */

public class FaultSampleCleaner {

  public void recomputeLikelihoods(
    final FaultSkin[] skins, final float[][][] fx) {
    final int ns = skins.length;
    final float[][][][] snd = semblanceNumDen(fx);
    loop(ns,new LoopInt() {
      public void compute(int is) {
        recomputeLikelihoods(skins[is],snd);
      }
    });
  }

  public void recomputeLikelihoods(FaultSkin skin, float[][][][] snd) {
    int nsmooth = 50;
    setNumAndDen(skin,snd);
    smoothNumAndDen(nsmooth,skin);
    computeLikelihoods(skin);
  }

  private void setNumAndDen(FaultSkin skin, float[][][][] snd) {
    float[][][] sn = snd[0];
    float[][][] sd = snd[1];
    for (FaultCell cell:skin) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      cell.setNum(sn[i3][i2][i1]);
      cell.setDen(sd[i3][i2][i1]);
    }
  }

  private void computeLikelihoods(FaultSkin skin) {
    for (FaultCell cell:skin) {
      float sm = cell.num/cell.den;
      sm *= sm;
      sm *= sm;
      sm *= sm;
      float fl = 1f-sm;
      cell.setFl(fl);
    }

  }

  private void smoothNumAndDen(int nsmooth, FaultSkin skin) {
    FaultCell.GetN getter = new FaultCell.GetN() {
      public float[] get(FaultCell cell) {
        return new float[]{cell.num,cell.den};
      }
    };
    FaultCell.SetN setter = new FaultCell.SetN() {
      public void set(FaultCell cell, float[] nds) {
        float num = nds[0]; 
        float den = nds[1]; 
        cell.setNum(num);
        cell.setDen(den);
      }
    };
    for (int ismooth=0; ismooth<nsmooth; ++ismooth)
      smoothN(getter,setter,skin);
  }

  void smoothN(FaultCell.GetN getter, FaultCell.SetN setter, FaultSkin skin) {
    int nval = 2;
    FaultCell[] cells = skin.getCells();
    int nc = cells.length;
    float[][] vals = new float[nc][nval];
    FaultCell[] cellNabors = new FaultCell[4];
    for (int ic=0; ic<nc; ++ic) {
      FaultCell cell = cells[ic];
      vals[ic] = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      cellNabors[2] = cell.cl;
      cellNabors[3] = cell.cr;
      float cnt = 1f;
      for (FaultCell cellNabor:cellNabors) {
        if (cellNabor!=null) {
          cnt++;
          float[] valsNabor = getter.get(cellNabor);
          for (int ival=0; ival<nval; ++ival)
            vals[ic][ival] += valsNabor[ival];
        }
      }
      if (cnt>1) {
        vals[ic][0] /= cnt;
        vals[ic][1] /= cnt;
      }
    }
    for (int ic=0; ic<nc; ++ic) {
      setter.set(cells[ic],vals[ic]);
    }
  }


    // Computes fault semblance numerators and denominators.
  private static float[][][][] semblanceNumDen(final float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] sn = new float[n3][n2][n1];
    float[][][] sd = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8,2,2);
    EigenTensors3 ets = lof.applyForTensors(f);
    ets.setEigenvalues(0.05f,1.0f,1.0f);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,8,f,sn);
    lsf.apply(ets,8,mul(f,f),sd);
    mul(sn,sn,sn);
    return new float[][][][]{sn,sd};
  }


  private void computeNumAndDen(
    float[][][] fx, float[][][] u1, float[][][] u2, float[][][] u3, 
    FaultSkin skin) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    SincInterpolator si = new SincInterpolator();
    float offset = 2f;
    for (FaultCell cell:skin) {
      float x1 = cell.getX1();
      float x2 = cell.getX2();
      float x3 = cell.getX3();
      float w1 = cell.getW1();
      float w2 = cell.getW2();
      float w3 = cell.getW3();
      float d2 =  cell.getV3()*offset;
      float d3 = -cell.getV2()*offset;
      float x2p = x2+d2;
      float x3p = x3+d3;
      float x2m = x2-d2;
      float x3m = x3-d3;
      float u1m = imageValueAt(x1,x2m,x3m,u1);
      float u2m = imageValueAt(x1,x2m,x3m,u2);
      float u3m = imageValueAt(x1,x2m,x3m,u3);
      float u1p = imageValueAt(x1,x2p,x3p,u1);
      float u2p = imageValueAt(x1,x2p,x3p,u2);
      float u3p = imageValueAt(x1,x2p,x3p,u3);
      float num = 0f;
      float den = 0f;
      for (float d=1; d<3; d++) {
        float v1p = d*w1;
        float v2p = d*w2;
        float v3p = d*w3;
        float v1m = -v1p;
        float v2m = -v2p;
        float v3m = -v3p;
        float vup = v1p*u1p+v2p*u2p+v3p*u3p;
        float vum = v1m*u1m+v2m*u2m+v3m*u3m;
        float p1p = v1p-vup*u1p+x1;
        float p2p = v2p-vup*u2p+x2p;
        float p3p = v3p-vup*u3p+x3p;
        float p1m = v1m-vum*u1m+x1;
        float p2m = v2m-vum*u2m+x2m;
        float p3m = v3m-vum*u3m+x3m;
        float fxp = si.interpolate(s1,s2,s3,fx,p1p,p2p,p3p);
        float fxm = si.interpolate(s1,s2,s3,fx,p1m,p2m,p3m);
        num += (fxp+fxm);
        den += (fxp*fxp+fxm*fxm);
      }
      cell.setDen(den);
      cell.setNum(num*num);
    }
  }

  private static float imageValueAt(
    float p1, float p2, float p3, float[][][]f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    int i1 = max(0,min(n1-1,round(p1)));
    int i2 = max(0,min(n2-1,round(p2)));
    int i3 = max(0,min(n3-1,round(p3)));
    return f[i3][i2][i1];
  }


}


