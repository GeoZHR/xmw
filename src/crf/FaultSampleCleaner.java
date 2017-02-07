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
    final FaultSkin[] skins, final float[][][] fx, 
    final float[][][] p2, final float[][][] p3) {
    final int ns = skins.length;
    final float[][][][] snd = semblanceNumDen(fx,p2,p3);
    loop(ns,new LoopInt() {
      public void compute(int is) {
        recomputeLikelihoods(skins[is],snd);
      }
    });
  }

  public void recomputeLikelihoods(FaultSkin skin, float[][][][] snd) {
    int nsmooth = 100;
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

    LocalOrientFilter lof = new LocalOrientFilter(8,2,2);
    EigenTensors ets = lof.applyForTensors(f);
    float[][][] sn = new float[n3][n2][n1];
    float[][][] sd = new float[n3][n2][n1];

    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      float[] xmm = new float[n1];
      float[] xm0 = new float[n1];
      float[] xmp = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xpm = new float[n1];
      float[] xp0 = new float[n1];
      float[] xpp = new float[n1];
      float[] gmm = new float[n1];
      float[] gm0 = new float[n1];
      float[] gmp = new float[n1];
      float[] g0m = new float[n1];
      float[] g0p = new float[n1];
      float[] gpm = new float[n1];
      float[] gp0 = new float[n1];
      float[] gpp = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p2mm = p2[i3m][i2m];
        float[] p2mp = p2[i3m][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p20p = p2[i3 ][i2p];
        float[] p2pm = p2[i3p][i2m];
        float[] p2pp = p2[i3p][i2p];
        float[] p3mm = p3[i3m][i2m];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3mp = p3[i3m][i2p];
        float[] p3pm = p3[i3p][i2m];
        float[] p3p0 = p3[i3p][i2 ];
        float[] p3pp = p3[i3p][i2p];
        float[] sn32 = sn[i3][i2];
        float[] sd32 = sd[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          xmm[i1] = i1-p3mm[i1]-p2mm[i1];
          xm0[i1] = i1-p3m0[i1]         ;
          xmp[i1] = i1-p3mp[i1]+p2mp[i1];
          x0m[i1] = i1         -p20m[i1];
          x0p[i1] = i1         +p20p[i1];
          xpm[i1] = i1+p3pm[i1]-p2pm[i1];
          xp0[i1] = i1+p3p0[i1]         ;
          xpp[i1] = i1+p3pp[i1]+p2pp[i1];
        }
        si.interpolate(n1,1.0,0.0,fmm,n1,xmm,gmm);
        si.interpolate(n1,1.0,0.0,fm0,n1,xm0,gm0);
        si.interpolate(n1,1.0,0.0,fmp,n1,xmp,gmp);
        si.interpolate(n1,1.0,0.0,f0m,n1,x0m,g0m);
        si.interpolate(n1,1.0,0.0,f0p,n1,x0p,g0p);
        si.interpolate(n1,1.0,0.0,fpm,n1,xpm,gpm);
        si.interpolate(n1,1.0,0.0,fp0,n1,xp0,gp0);
        si.interpolate(n1,1.0,0.0,fpp,n1,xpp,gpp);
        float[] hmm = gmm, hm0 = gm0, hmp = gmp;
        float[] h0m = g0m, h00 = f00, h0p = g0p;
        float[] hpm = gpm, hp0 = gp0, hpp = gpp;
        if (            i2==0   ) h0m = h00;
        if (            i2==n2-1) h0p = h00;
        if (i3==0               ) hm0 = h00;
        if (i3==n3-1            ) hp0 = h00;
        if (i3==0    && i2==0   ) hmm = h00;
        if (i3==0    && i2==n2-1) hmp = h00;
        if (i3==n3-1 && i2==0   ) hpm = h00;
        if (i3==n3-1 && i2==n2-1) hpp = h00;
        for (int i1=0; i1<n1; ++i1) {
          float hmmi = hmm[i1];
          float hm0i = hm0[i1];
          float hmpi = hmp[i1];
          float h0mi = h0m[i1];
          float h00i = h00[i1];
          float h0pi = h0p[i1];
          float hpmi = hpm[i1];
          float hp0i = hp0[i1];
          float hppi = hpp[i1];
          float sumn = hmmi+hm0i+hmpi+
                       h0mi+h00i+h0pi+
                       hpmi+hp0i+hppi;
          float sumd = hmmi*hmmi+hm0i*hm0i+hmpi*hmpi+
                       h0mi*h0mi+h00i*h00i+h0pi*h0pi+
                       hpmi*hpmi+hp0i*hp0i+hppi*hppi;
          sn32[i1] = sumn*sumn;
          sd32[i1] = 9.0f*sumd;
        }
      }
    }});
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


