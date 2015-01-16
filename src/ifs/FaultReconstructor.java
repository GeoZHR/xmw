/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import static ifs.FaultGeometry.*;

/**
 * Reconstruct fault images from fault cells.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.15
 */

public class FaultReconstructor {

  public FaultReconstructor(int n1, int n2, int n3, FaultCell[] fcs){
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fcs = fcs;
  }

  public float[][][][][] applyForFaultImages(int minSkinSize) {
    float[][][] fl1 = new float[_n3][_n2][_n1];
    float[][][] fp1 = new float[_n3][_n2][_n1];
    float[][][] ft1 = new float[_n3][_n2][_n1];
    float[][][] fl2 = new float[_n3][_n2][_n1];
    float[][][] fp2 = new float[_n3][_n2][_n1];
    float[][][] ft2 = new float[_n3][_n2][_n1];
    HashSet<FaultCell[]> hsf = sortFaultCells();
    int is=0;
    for (FaultCell[] fcs:hsf) {
      int nc = fcs.length;
      if(nc<minSkinSize){continue;}
      System.out.println("is="+is); is++;
      float[][][] fl = accumulateGaussians(fcs);
      float[][][][] flpt = smooth(fl);
      for (int i3=0; i3<_n3; ++i3) {
        for (int i2=0; i2<_n2; ++i2) {
          for (int i1=0; i1<_n1; ++i1) {
            float flt = fl1[i3][i2][i1];
            float fpt = fp1[i3][i2][i1];
            float ftt = ft1[i3][i2][i1];
            float fli = flpt[0][i3][i2][i1];
            float fpi = flpt[1][i3][i2][i1];
            float fti = flpt[2][i3][i2][i1];
            if(ftt==0.0f) {
              fl1[i3][i2][i1] = fli;
              fp1[i3][i2][i1] = fpi;
              ft1[i3][i2][i1] = fti;
            } else if(fli>flt && strikeDif(fpi,fpt)<30f) {
              fl1[i3][i2][i1] = fli;
              fp1[i3][i2][i1] = fpi;
              ft1[i3][i2][i1] = fti;
            } else if(fli>fl2[i3][i2][i1] && strikeDif(fpi,fpt)>=30f) {
              fl2[i3][i2][i1] = fli;
              fp2[i3][i2][i1] = fpi;
              ft2[i3][i2][i1] = fti;
            }
          }
        }
      }
    }
    float[][][][] fpt1 = new float[][][][]{fl1,fp1,ft1};
    float[][][][] fpt2 = new float[][][][]{fl2,fp2,ft2};
    return new float[][][][][]{fpt1,fpt2}; 
  }

  private float strikeDif(float fp1, float fp2) {
    float dfp = abs(fp1-fp2);
    if(dfp>(360-30)){dfp=abs(360-dfp);}
    return dfp;
  }

  private float[][][][] smooth(float[][][] fl) {
    float sigma1 = 8.0f;
    float sigma2 = 2.0f;
    float sigma  = 20.f;
    float c = sigma*sigma*0.5f;
    float[][][] u1  = new float[_n3][_n2][_n1];
    float[][][] u2  = new float[_n3][_n2][_n1];
    float[][][] u3  = new float[_n3][_n2][_n1];
    float[][][] fls = new float[_n3][_n2][_n1];
    float[][][] fps = new float[_n3][_n2][_n1];
    float[][][] fts = new float[_n3][_n2][_n1];
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    EigenTensors3 d = lof.applyForTensors(fl);
    d.setEigenvalues(0.001f,1.0f,1.0f);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.1,20);
    lsf.apply(d,c,fl,fls);
    div(fls,max(fls),fls); 
    lof.applyForNormal(fls,u1,u2,u3);
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          float fli = fls[i3][i2][i1];
          if(u2i!=0f && u3i!=0f && fli>0.1f) {
            fts[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
            fps[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
          }
        }
      }
    }
    return new float[][][][]{fls,fps,fts};
  }



  private float[][][] accumulateGaussians(FaultCell[] fcs) {
    float[][][] fg = new float[_n3][_n2][_n1];
    for (FaultCell fc:fcs) {
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      fg[i3][i2][i1] = fc.fl;
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    rgf.apply0XX(fg,fg);
    rgf.applyX0X(fg,fg);
    rgf.applyXX0(fg,fg);
    return div(fg,max(fg));
  }

  private HashSet<FaultCell[]> sortFaultCells() {
    HashSet<FaultCell> hfc = new HashSet<FaultCell>();
    for (FaultCell fc:_fcs){hfc.add(fc);}
    int nc = hfc.size();
    while(nc>100) {

    }

    final KdTree kt = setStrikeKdTree(hfc,id);
    final HashSet<Inte> hsf = new HashSet<FaultCell[]>();
    Parallel.loop(0,360,1,new Parallel.LoopInt() {
    public void compute(int fp) {
      FaultCell[] fcs = findStrikeInWindow(fp,kt);
      if(fcs!=null){hsf.add(fcs);}
    }});
    return hsf;
    //return hsf.toArray(new FaultCell[hsf.size()][]);
  }

  private KdTree setStrikeKdTree() {
    int nc = _fcs.length;
    float[][] xc = new float[1][nc];
    for (int ic=0; ic<nc; ++ic) {
      xc[0][ic] = _fcs[ic].fp;
    }
    return new KdTree(xc);
  }

  private FaultCell[] cellsInStrikeWd(HashSet<Integer> hsc) {
    int dp = 8;
    int maxN = 0;
    int nc = hsc.size();
    int[] dc = new int[nc];
    float[] pm1 = new float[1];
    float[] pp1 = new float[1];
    float[] pm2 = new float[1];
    float[] pp2 = new float[1];
    float[] pm3 = new float[1];
    float[] pp3 = new float[1];
    float[][] xc = setKdTreeStrike(hsc,dc);
    KdTree kt = new KdTree(xc);
    HashSet<Integer> hst = new HashSet<Integer>();
    for (int pt=0; pt<=360; ++pt) {
      float pmi = pt-dp;
      float ppi = pt+dp;
      int nd1=0, nd2=0, nd3=0;
      pm3[0] = pmi; pp3[0] = ppi;
      int[] id1=null, id2=null, id3=null;
      if (pmi<0.0f) {
        pm3[0] = 0.0f;
        pp1[0] = 360f;
        pm1[0] = 360f+pmi;
        id1 = kt.findInRange(pm1,pp1);
      }
      if(ppi>360f) {
        pp3[0] = 360f;
        pm2[0] = 0.0f;
        pp2[0] = ppi-360f;
        id2 = kt.findInRange(pm2,pp2);
      }
      id3 = kt.findInRange(pm3,pp3);
      if(id1!=null){nd1=id1.length;}
      if(id2!=null){nd2=id2.length;}
      if(id3!=null){nd3=id3.length;}
      int nd = nd1+nd2+nd3;
      if(nd>maxN) {
        maxN = nd;
        hst.clear();
        if(id1!=null) {
          for (int ik=0; ik<nd1; ++ik){
            int ip = id1[ik];
            hst.add(dc[ip]);
          }
        }
        if(id2!=null) {
          for (int ik=0; ik<nd2; ++ik){
            int ip = id2[ik];
            hst.add(dc[ip]);
          }
        }
        if(id3!=null) {
          for (int ik=0; ik<nd3; ++ik){
            int ip = id3[ik];
            hst.add(dc[ip]);
          }
        }
      }
    }
    int ik=-1;
    FaultCell[] fc = new FaultCell[maxN];
    for (int ic:hst) {
      ik++;
      fc[ik] = _fc[ic];
    }
    return fc;
  }

  private KdTree setKdTree(FaultCell[] fcs) {
    int nc = fcs.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
    }
    return new KdTree(xc);
  }

  private int _n1, _n2, _n3;
  private FaultCell[] _fcs;
}


