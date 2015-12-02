/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipfx;

import java.util.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

import static ipfx.FaultGeometry.*;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes fault skins from images of fault likelihoods, strikes and dips. A
 * fault skin is a linked list of fault cells. Each fault cell is an oriented
 * point located on a ridge in an image of fault likelihood. Each image sample
 * corresponds to either no cell or one cell.
 * <p>
 * A cell has up to four neighbors ("nabors") that lie above, below, left and
 * right of the cell when viewed from above the fault, that is, when looking
 * from the hanging wall toward the footwall. Links to nabors enables cells to
 * form a skin of connected cells, which represents a fault.
 * <p>
 * Links to left and right cell nabors can be used to iterate over all cells
 * along a fault trace, a path of constant depth that is everywhere tangent to
 * fault strike. Likewise, links to cell nabors above and below a cell can be
 * used to iterate up or down a fault. However, this simple up or down
 * iteration typically does not coincide with a fault curve that is everywhere
 * tangent to fault dip.
 * <p>
 * A fault skin is grown by linking nearby cells having similar properties,
 * beginning with a seed cell that has sufficiently high fault likelihood.
 * Several methods in this class set parameters that control this growing
 * process.
 *
 * @author Xinming Wu and Dave Hale, Colorado School of Mines
 * @version 2015.01.22
 */
public class FaultSkinnerX {

  /**
   * Constructs a fault skinner with default parameters.
   */
  public FaultSkinnerX() {
    _fs1min = -Float.MAX_VALUE;
    _fs1max =  Float.MAX_VALUE;
    _fllo = 0.2f;
    _flhi = 0.8f;
    _dflmax = 0.2f;
    _dfpmax = 10.0f;
    _dftmax = 10.0f;
    _dnpmax = 0.5f;
    _ds1max = 1.0f;
    _ncsmin = 400;
  }

  public void setParameters(float dfp, float dft, float dnp) {
    _dfpmax1 = dfp;
    _dftmax1 = dft;
    _dnpmax1 = dnp;
  }

  public void setSkinning(boolean useOldCells) {
    _useOldCells = useOldCells;
  }

  public FaultCell[] createCells(
    float x1, float x2, float x3,
    float fl, float fp, float ft,
    FaultCell[] fcs) 
  {
    int ic = 0;
    int nc = fcs.length+1;
    FaultCell cell = new FaultCell(x1,x2,x3,fl,fp,ft);
    FaultCell[] cells = new FaultCell[nc];
    for (FaultCell fc:fcs){cells[ic]=fc;ic++;}
    cells[nc-1] = cell;
    return cells;
  }

  public void resetCells(FaultCell[] cells) {
    int nc = cells.length;
    for (int ic=0; ic<nc; ++ic) {
      cells[ic].ca = null;
      cells[ic].cb = null;
      cells[ic].cl = null;
      cells[ic].cr = null;
      cells[ic].skin=null;
      cells[ic].notUsed=true;
    }
  }

  public FaultCell[] resetCells(FaultSkin[] skins) {
    ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
    for (FaultSkin skin:skins) {
      int ic = -1;
      int nc = skin.size();
      float[] fl = new float[nc];
      for (FaultCell fc:skin){ic++;fl[ic]=fc.fl;}
      breakIntersects(skin);
      ic = -1;
      for (FaultCell cell:skin) {
        ic++;
        if (cell.fl==0f){continue;}
        cell.ca=null;
        cell.cb=null;
        cell.cl=null;
        cell.cr=null;
        cell.skin=null;
        cell.fl = fl[ic];
        cell.notUsed=true;
        fcl.add(cell);
      }
    }
    return fcl.toArray(new FaultCell[0]);
  }


  private void breakIntersects(FaultSkin skin) {
    int ic = -1;
    int nc = skin.size();
    float[][] xc = new float[3][nc];
    FaultCell[] fcs = new FaultCell[nc];
    for (FaultCell fc:skin) {
      ic++;
      fcs[ic] = fc;
      xc[0][ic] = fc.x1;
      xc[1][ic] = fc.x2;
      xc[2][ic] = fc.x3;
    }
    KdTree kt = new KdTree(xc);
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (FaultCell cell:skin){hsc.add(cell);}
    ArrayList<FaultCell[]> fca = new ArrayList<FaultCell[]>();
    while(hsc.size()>0) {
      FaultCell[] cells = findInFpWind(hsc);
      fca.add(cells);
    }
    FaultCell[][] cells = fca.toArray(new FaultCell[0][0]);
    int ns = cells.length;
    boolean needBreak = false;
    for(int is=0; is<ns; ++is) {
      int nk = cells[is].length;
      float den = (float)nk/(float)nc;
      if(den>=0.1f&&den<0.9f){
        trace("test");
        needBreak=true;
        for (int ik=0; ik<nk; ++ik){
          cells[is][ik].fl=is+1;
        }
      }
    }
    if(!needBreak) {return;}
    for (FaultCell fc:skin){checkIntersect(fc);}
    for (FaultCell fc:skin) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      if(fc.intersect) {
        xmin[0] = fc.x1-10;
        xmax[0] = fc.x1+10;
        xmin[1] = fc.x2-10;
        xmax[1] = fc.x2+10;
        xmin[2] = fc.x3-10;
        xmax[2] = fc.x3+10;
        int[] id = kt.findInRange(xmin,xmax);
        int nd = id.length;
        if(nd<1){continue;}
        for (int ik=0; ik<nd; ++ik) {
          int ip = id[ik];
          float dfp = absDeltaFp(fc,fcs[ip]);
          if(dfp<10f){fcs[ip].fl=0f;}
        }
      }
    }
  }

  private void checkIntersect(FaultCell cell) {
    float flc = cell.fl;
    FaultCell ca = cell.ca;
    FaultCell cb = cell.cb;
    FaultCell cl = cell.cl;
    FaultCell cr = cell.cr;
    if(flc<1f){return;}
    if(ca!=null) {
      float fla = ca.fl;
      if(flc!=fla&&fla>1f){cell.intersect=true;}
    }
    if(cb!=null) {
      float flb = cb.fl;
      if(flc!=flb&&flb>1f){cell.intersect=true;}
    }
    if(cl!=null) {
      float fll = cl.fl;
      if(flc!=fll&&fll>1f){cell.intersect=true;}
    }
    if(cr!=null) {
      float flr = cr.fl;
      if(flc!=flr&&flr>1f){cell.intersect=true;}
    }
  }

  private FaultCell[] findInFpWind(HashSet<FaultCell> hsc) {
    int maxN = 0;
    float dp = 20f;
    int nc = hsc.size();
    float[] fpms = new float[2];
    float[] fpm1 = new float[1];
    float[] fpp1 = new float[1];
    float[] fpm2 = new float[1];
    float[] fpp2 = new float[1];
    float[] fpm3 = new float[1];
    float[] fpp3 = new float[1];
    FaultCell[] cells = new FaultCell[nc]; 
    KdTree kt = setFpKdTree(hsc,fpms,cells);
    HashSet<Integer> hsi = new HashSet<Integer>();
    float fpMin = fpms[0], fpMax=fpms[1];
    for (float pt=fpMin; pt<=fpMax; pt+=1f) {
      float fpm = pt-dp;
      float fpp = pt+dp;
      int nd1=0, nd2=0, nd3=0;
      fpm3[0] = fpm; fpp3[0] = fpp;
      int[] id1=null, id2=null, id3=null;
      if (fpm<0.0f) {
        fpm3[0] = 0.0f;
        fpp1[0] = 360f;
        fpm1[0] = 360f+fpm;
        id1 = kt.findInRange(fpm1,fpp1);
      }
      if(fpp>360f) {
        fpp3[0] = 360f;
        fpm2[0] = 0.0f;
        fpp2[0] = fpp-360f;
        id2 = kt.findInRange(fpm2,fpp2);
      }
      id3 = kt.findInRange(fpm3,fpp3);
      if(id1!=null){nd1=id1.length;}
      if(id2!=null){nd2=id2.length;}
      if(id3!=null){nd3=id3.length;}
      int nd = nd1+nd2+nd3;
      if(nd>maxN) {
        maxN = nd;
        hsi.clear();
        if(id1!=null) {
          for (int ik=0; ik<nd1; ++ik){
            hsi.add(id1[ik]);
          }
        }
        if(id2!=null) {
          for (int ik=0; ik<nd2; ++ik){
            hsi.add(id2[ik]);
          }
        }
        if(id3!=null) {
          for (int ik=0; ik<nd3; ++ik){
            hsi.add(id3[ik]);
          }
        }
      }
    }
    int ic = -1;
    int ni = hsi.size();
    FaultCell[] fcs = new FaultCell[ni];
    for (int ik:hsi) {
      ic ++;
      FaultCell cell = cells[ik];
      fcs[ic] = cell;
      hsc.remove(cell);
    }
    return fcs;
  }


  private KdTree setFpKdTree(
    HashSet<FaultCell> hsc, float[] fpms, FaultCell[] cells) 
  {
    int ic = -1;
    int nc = hsc.size();
    float[][] xc = new float[1][nc];
    float fpMin = 360f;
    float fpMax = 0.0f;
    for (FaultCell cell:hsc) {
      ic ++;
      float fp = cell.fp;
      xc[0][ic] = fp;
      if(fp>fpMax){fpMax=fp;}
      if(fp<fpMin){fpMin=fp;}
      cells[ic] = cell;
    }
    fpms[0]=fpMin;
    fpms[1]=fpMax;
    return new KdTree(xc);
  }


  public void getFl(FaultSkin[] sks, float[][][] fl) {
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        int i1i = fc.i1;
        int i2i = fc.i2;
        int i3i = fc.i3;
        fl[i3i][i2i][i1i] = fc.fl;
      }
    }
  }

  public void getFls(FaultSkin[] sks, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        int i1m = fc.i1;
        int i2m = fc.i2m;
        int i3m = fc.i3m;
        int i1p = fc.i1;
        int i2p = fc.i2p;
        int i3p = fc.i3p;
        if(i2m>=n2){continue;}
        if(i3m>=n3){continue;}
        fl[i3m][i2m][i1m] = fc.fl;
        fl[i3p][i2p][i1p] = fc.fl;
      }
    }
  }

  /**
   * Sets the minimum number of cells in a skin. Skins smaller than this will
   * be discarded.
   * <p>
   * The default minimum skin size is 400.
   */
  public void setMinSkinSize(int minSize) {
    _ncsmin = minSize;
  }

  /**
   * Sets lower and upper bounds on fault throw for cells in a skin.
   * These bounds should be set only after fault dip slips have been 
   * computed for all cells used to grow skins.
   * <p>
   * The default bounds are huge, so that throws are unrestricted.
   * @param minThrow the lower bound.
   * @param maxThrow the upper bound.
   */
  public void setMinMaxThrow(double minThrow, double maxThrow) {
    Check.argument(minThrow<=maxThrow,"minThrow does not exceed maxThrow");
    _fs1min = (float)minThrow;
    _fs1max = (float)maxThrow;
  }

  /**
   * Sets fault likelihood thresholds used to grow skins. Cells in a skin
   * should have, or be connected to cells that have, high fault likelihoods.
   * All cells in a skin will have fault likelihoods not less than the lower
   * threshold. At least one cell in a skin will have a fault likelihood not
   * less than the upper threshold. 
   * <p>
   * The default thresholds are 0.2 and 0.8, respectively.
   * @param lowerLikelihood lower threshold for fault likelihood.
   * @param upperLikelihood upper threshold for fault likelihood.
   */
  public void setGrowLikelihoods(
      double lowerLikelihood, double upperLikelihood) {
    Check.argument(lowerLikelihood<=upperLikelihood,
        "lowerLikelihood does not exceed upperLikelihood");
    _fllo = (float)lowerLikelihood;
    _flhi = (float)upperLikelihood;
  }

  /**
   * Sets the maximum difference in fault likelihood for a cell and its nabors.
   * <p>
   * The default maximum difference is 0.2.
   * @param maxDeltaLikelihood upper bound on difference in fault likelihood.
   */
  public void setMaxDeltaLikelihood(double maxDeltaLikelihood) {
    _dflmax = (float)maxDeltaLikelihood;
  }

  /**
   * Sets the maximum difference in fault strike for a cell and its nabors.
   * @param maxDeltaStrike upper bound on difference in fault strike.
   * <p>
   * The default maximum difference is 30 degrees.
   */
  public void setMaxDeltaStrike(double maxDeltaStrike) {
    _dfpmax = (float)maxDeltaStrike;
  }

  /**
   * Sets the maximum difference in fault dip for a cell and its nabors.
   * @param maxDeltaDip upper bound on difference in fault dip.
   * <p>
   * The default maximum difference is 10 degrees.
   */
  public void setMaxDeltaDip(double maxDeltaDip) {
    _dftmax = (float)maxDeltaDip;
  }

  /**
   * Sets the maximum difference in fault throw for a cell and its nabors.
   * @param maxDeltaThrow upper bound on difference in fault throw.
   * <p>
   * The default maximum difference is 1.0 samples.
   */
  public void setMaxDeltaThrow(double maxDeltaThrow) {
    _ds1max = (float)maxDeltaThrow;
  }

  /**
   * Sets the threshold planar distance. A cell should lie near the planes of
   * nabor cells. The specified threshold is the maximum distance to nabor
   * planes.
   * <p>
   * The default maximim planar distance is 0.5 samples.
   * @param maxPlanarDistance upper bound on planar distance, in samples.
   */
  public void setMaxPlanarDistance(double maxPlanarDistance) {
    _dnpmax = (float)maxPlanarDistance;
  }

  public FaultSkin[] findSkinsXX(FaultCell[] cells, float[][][] fl) {
    FaultSkin[] sks = skinsXX(cells,fl);
    //return sks;
    return reskin(sks,fl);
    //return resetSkins(sks,fl);
  }

  public FaultSkin[] reskin(FaultSkin[] sks, float[][][] fl) {
    int ns = sks.length;
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    ArrayList<FaultSkin> skl = new ArrayList<FaultSkin>();
    for (int is=0; is<ns; ++is) {
      System.out.println("skin="+is);
      FaultCell[] fcs = sks[is].getCells();
      float[][][][] flpt = new float[3][n3][n2][n1];
      faultImagesFromCells(30,fcs,flpt);
      FaultSkinner fs = new FaultSkinner();
      fs.setMaxDeltaStrike(10);
      //fs.setGrowLikelihoods(0.2,0.5);
      fs.setGrowLikelihoods(0.1,0.2);
      fs.setMinSkinSize(round(fcs.length/2));
      FaultCell[] cells = fs.findCells(flpt);
      System.out.println("cells="+cells.length);
      FaultSkin[] skins = fs.findSkins(cells);
      if(skins.length<1){continue;}
      for (FaultCell fc:skins[0]) {
        float x1 = fc.x1;
        float x2 = fc.x2;
        float x3 = fc.x3;
        int i1i = round(x1);
        int i2i = round(x2);
        int i3i = round(x3);
        int i2n = i2i;
        int i3n = i3i;
        if(i2i>x2) {i2i -= 1;}
        else       {i2n += 1;}
        if(i3i>x3) {i3i -= 1;}
        else       {i3n += 1;}
        if(i2i<0||i2n>=n2) {fc.fl=fl[i3i][i2i][i1i]; continue;}
        if(i3i<0||i3n>=n3) {fc.fl=fl[i3i][i2i][i1i]; continue;}
        x2 -= i2i; x3 -= i3i;
        float fla = fl[i3i][i2i][i1i];
        float flb = fl[i3i][i2n][i1i];
        float flc = fl[i3n][i2i][i1i];
        float fld = fl[i3n][i2n][i1i];
        fc.fl=fla*(1f-x2)*(1f-x3)+flb*(1f-x3)*x2+flc*(1f-x2)*x3+fld*x2*x3;
      }
      skl.add(skins[0]);
    }
    return skl.toArray(new FaultSkin[0]);
  }

  public float[][][][] faultImagesFromCellsX(
    final int n1, final int n2, final int n3, final int dt, final FaultCell[] fc)
  {
    int nc = fc.length;
    final float[][][] fl = new float[n3][n2][n1];
    final float[][][] fp = new float[n3][n2][n1];
    final float[][][] ft = new float[n3][n2][n1];
    float sigmaw = 4.0f;
    float sigmav = dt*100f;
    float sigmau = dt*100f;
    final float mark = -360f;
    final float[][][] fpt = fillfloat(mark,n1,n2,n3);
    final float[][] xc = new float[3][nc];
    setKdTreeNodes(fc,xc,fpt);
    final int[] bs1 = setBounds(n1,xc[0]);
    final int[] bs2 = setBounds(n2,xc[1]);
    final int[] bs3 = setBounds(n3,xc[2]);
    final KdTree kt = new KdTree(xc);
    final float sv = 1.f/(sigmav*sigmav); 
    final float su = 1.f/(sigmau*sigmau); 
    final float sw = 1.f/(sigmaw*sigmaw); 
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      System.out.println("i3="+i3);
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          if(fpt[i3][i2][i1]==mark){continue;}
          xmin[1] = i2-dt; xmax[1] = i2+dt;
          xmin[2] = i3-dt; xmax[2] = i3+dt;
          xmin[0] = i1-dt/4; xmax[0] = i1+dt/4;
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<10){continue;}
          ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            FaultCell fci = fc[ip];
            fcl.add(fci);
          }
          FaultCell[] fcs = fcl.toArray(new FaultCell[0]);
          if(fcs==null){continue;}
          float wps = 0.0f;
          for (FaultCell fci:fcs) {
            float x1i = fci.i1;
            float x2i = fci.i2;
            float x3i = fci.i3;
            float dx1 = x1i-i1;
            float dx2 = x2i-i2;
            float dx3 = x3i-i3;

            float d11 = dx1*dx1;
            float d22 = dx2*dx2;
            float d33 = dx3*dx3;
            float d12 = dx1*dx2;
            float d13 = dx1*dx3;
            float d23 = dx2*dx3;

            float w11 = fci.w11;
            float w22 = fci.w22;
            float w33 = fci.w33;
            float w12 = fci.w12;
            float w13 = fci.w13;
            float w23 = fci.w23;

            float u11 = fci.u11;
            float u22 = fci.u22;
            float u33 = fci.u33;
            float u12 = fci.u12;
            float u13 = fci.u13;
            float u23 = fci.u23;

            float v11 = fci.v11;
            float v22 = fci.v22;
            float v33 = fci.v33;
            float v12 = fci.v12;
            float v13 = fci.v13;
            float v23 = fci.v23;

            float wd1 = w12*d12*2.0f;
            float wd2 = w13*d13*2.0f;
            float wd3 = w23*d23*2.0f;

            float ud1 = u12*d12*2.0f;
            float ud2 = u13*d13*2.0f;
            float ud3 = u23*d23*2.0f;

            float vd1 = v12*d12*2.0f;
            float vd2 = v13*d13*2.0f;
            float vd3 = v23*d23*2.0f;

            float wds = w11*d11+w22*d22+w33*d33;
            float uds = u11*d11+u22*d22+u33*d33;
            float vds = v11*d11+v22*d22+v33*d33;

            float gss = 0.0f;
            float flc = fci.fl;
            float wpi = flc;//pow(flc,2.f);
            gss += (wd1+wd2+wd3+wds)*sw;
            gss += (ud1+ud2+ud3+uds)*su;
            gss += (vd1+vd2+vd3+vds)*sv;
            float fli = exp(-gss)*wpi;
            fl[i3][i2][i1] += fli;
            wps += wpi;
          }
          if(wps!=0f) {
            fl[i3][i2][i1] /= wps;
          }
        }
      }
    }});
    computeStrikeDip(fl,fp,ft);
    return new float[][][][]{fl,fp,ft};
  }


  public void faultImagesFromCells(
    final int dt, final FaultCell[] fc, final float[][][][] flpt)
  {
    final float[][][] fl = flpt[0];
    final float[][][] fp = flpt[1];
    final float[][][] ft = flpt[2];
    int nc = fc.length;
    final int n3 = fl.length;
    final int n2 = fl[0].length;
    final int n1 = fl[0][0].length;
    float sigmaNor = 4.0f;
    final float mark = -360f;
    final float[][][] fpt = fillfloat(mark,n1,n2,n3);
    final float[][] xc = new float[3][nc];
    setKdTreeNodes(fc,xc,fpt);
    final int[] bs1 = setBounds(n1,xc[0]);
    final int[] bs2 = setBounds(n2,xc[1]);
    final int[] bs3 = setBounds(n3,xc[2]);
    final KdTree kt = new KdTree(xc);
    final float sv = 0.25f/(dt*dt); 
    final float su = 0.25f/(dt*dt); 
    final float sw = 1.0f/(sigmaNor*sigmaNor); 
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          if(fpt[i3][i2][i1]==mark){continue;}
          xmin[1] = i2-dt; xmax[1] = i2+dt;
          xmin[2] = i3-dt; xmax[2] = i3+dt;
          xmin[0] = i1-dt/4; xmax[0] = i1+dt/4;
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<10){continue;}
          ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            FaultCell fci = fc[ip];
            float del = fci.fp-fpt[i3][i2][i1];
            float dfp = min(abs(del),abs(del+360.0f),abs(del-360.0f));
            if(dfp<=90f) {fcl.add(fci);}
          }
          FaultCell[] cells = fcl.toArray(new FaultCell[0]);
          FaultCell[] fcs = nabors(i1,i2,i3,cells);
          if(fcs==null){continue;}
          float wps = 0.0f;
          for (FaultCell fci:fcs) {
            float x1i = fci.i1;
            float x2i = fci.i2;
            float x3i = fci.i3;

            float dx1 = x1i-i1;
            float dx2 = x2i-i2;
            float dx3 = x3i-i3;
            float d11 = dx1*dx1;
            float d22 = dx2*dx2;
            float d33 = dx3*dx3;
            float d12 = dx1*dx2;
            float d13 = dx1*dx3;
            float d23 = dx2*dx3;

            float w11 = fci.w11;
            float w22 = fci.w22;
            float w33 = fci.w33;
            float w12 = fci.w12;
            float w13 = fci.w13;
            float w23 = fci.w23;

            float u11 = fci.u11;
            float u22 = fci.u22;
            float u33 = fci.u33;
            float u12 = fci.u12;
            float u13 = fci.u13;
            float u23 = fci.u23;

            float v11 = fci.v11;
            float v22 = fci.v22;
            float v33 = fci.v33;
            float v12 = fci.v12;
            float v13 = fci.v13;
            float v23 = fci.v23;

            float wd1 = w12*d12*2.0f;
            float wd2 = w13*d13*2.0f;
            float wd3 = w23*d23*2.0f;

            float ud1 = u12*d12*2.0f;
            float ud2 = u13*d13*2.0f;
            float ud3 = u23*d23*2.0f;

            float vd1 = v12*d12*2.0f;
            float vd2 = v13*d13*2.0f;
            float vd3 = v23*d23*2.0f;

            float wds = w11*d11+w22*d22+w33*d33;
            float uds = u11*d11+u22*d22+u33*d33;
            float vds = v11*d11+v22*d22+v33*d33;

            float gss = 0.0f;
            float flc = fci.fl;
            float wpi = pow(flc,2.f);
            gss += (wd1+wd2+wd3+wds)*sw;
            gss += (ud1+ud2+ud3+uds)*su;
            gss += (vd1+vd2+vd3+vds)*sv;
            float fli = exp(-gss)*wpi;
            fl[i3][i2][i1] += fli;
            wps += wpi;
          }
          if(wps!=0f) {
            fl[i3][i2][i1] /= wps;
          }
        }
      }
    }});
    computeStrikeDip(fl,fp,ft);
  }


  private void computeStrikeDip(
    float[][][] fl, float[][][] fp, float[][][] ft) 
  {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(8,4);
    lof.applyForNormal(fl,u1,u2,u3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fli = fl[i3][i2][i1];
      if(fli>0f) {
        int k1 = i1;
        int k2 = i2;
        int k3 = i3;
        if(k1==0) {k1=1;}
        if(k2==0) {k2=1;}
        if(k3==0) {k3=1;}
        if(k1==n1-1) {k1=n1-2;}
        if(k2==n2-1) {k2=n2-2;}
        if(k3==n3-1) {k3=n3-2;}
        float u1i = -u1[k3][k2][k1];
        float u2i = -u2[k3][k2][k1];
        float u3i = -u3[k3][k2][k1];
        if(u2i!=0.0f && u3i!=0.0f) {
          ft[i3][i2][i1] = faultDipFromNormalVector(u1i,u2i,u3i);
          fp[i3][i2][i1] = faultStrikeFromNormalVector(u1i,u2i,u3i);
        }
      }
    }}}
  }

  public FaultCell[] nabors(
    int i1, int i2, int i3, FaultCell[] cells) 
  {
    ArrayList<Float> dsL = new ArrayList<Float>();
    ArrayList<Float> dsR = new ArrayList<Float>();
    ArrayList<FaultCell> fcL = new ArrayList<FaultCell>();
    ArrayList<FaultCell> fcR = new ArrayList<FaultCell>();
    int nbc=0, nac=0;
    for (FaultCell fci:cells) {
      float x1 = fci.x1;
      float x2 = fci.x2;
      float x3 = fci.x3;
      float d1 = x1-i1;
      float d2 = x2-i2;
      float d3 = x3-i3;
      float v2 =  x3;
      float v3 = -x2;
      float dd = d2*v2+d3*v3;
      float ds = d1*d1+d2*d2+d3*d3;
      if(abs(d1)>1f) {ds *= abs(d1);}
      if(dd<=0f) {dsL.add(ds);fcL.add(fci);} 
      else if(dd>=0f) {dsR.add(ds);fcR.add(fci);}
      if(d1>0f){nac++;}
      if(d1<0f){nbc++;}
    }
    int nb = 20;
    int nl = dsL.size();
    int nr = dsR.size();
    if(nl<10||nr<10) {return null;}
    if(nac<5||nbc<5) {return null;}
    //if(nl<nb||nr<nb) {nb = min(nl,nr);}
    ArrayList<FaultCell> nbs = new ArrayList<FaultCell>();
    if(nl<nb) {
      for (FaultCell fc:fcL) {nbs.add(fc);}
    } else {
      int il = 0;
      int[] idl = new int[nl];
      float[] dsl = new float[nl];
      FaultCell[] fcl = fcL.toArray(new FaultCell[0]);
      for (float dsi:dsL) {dsl[il]=dsi; idl[il]=il; il++;}
      quickIndexSort(dsl,idl);
      for (int ik=0; ik<nb; ++ik) {
        int ic = idl[ik];
        nbs.add(fcl[ic]);
      }
    }
    if(nr<nb) {
      for (FaultCell fc:fcR) {nbs.add(fc);}
    } else {
      int ir = 0;
      int[] idr = new int[nr];
      float[] dsr = new float[nr];
      FaultCell[] fcr = fcR.toArray(new FaultCell[0]);
      for (float dsi:dsR) {dsr[ir]=dsi; idr[ir]=ir; ir++;}
      quickIndexSort(dsr,idr);
      for (int ik=0; ik<nb; ++ik) {
        int ic = idr[ik];
        nbs.add(fcr[ic]);
      }
    }
    return nbs.toArray(new FaultCell[0]);
  }


  private float[] computeStrikeAndDip(
    float g11, float g12, float g13,
    float g22, float g23, float g33)
  {
    float[] us = solveEigenproblems(g11,g12,g13,g22,g23,g33);
    float u1i = us[0];
    float u2i = us[1];
    float u3i = us[2];
    if(u2i!=0.0f && u3i!=0.0f) {
      float ft = faultDipFromNormalVector(u1i,u2i,u3i);
      float fp = faultStrikeFromNormalVector(u1i,u2i,u3i);
      return new float[]{fp,ft};
    } else {return null;}
  }

  private float[] solveEigenproblems(
    float g11, float g12, float g13,
    float g22, float g23, float g33)
  {
    double[] e = new double[3];
    double[][] a = new double[3][3];
    double[][] z = new double[3][3];
    a[0][0] = g11;
    a[0][1] = g12;
    a[0][2] = g13;
    a[1][0] = g12;
    a[1][1] = g22;
    a[1][2] = g23;
    a[2][0] = g13;
    a[2][1] = g23;
    a[2][2] = g33;
    Eigen.solveSymmetric33(a,z,e);
    float u1i = (float)z[0][0];
    float u2i = (float)z[0][1];
    float u3i = (float)z[0][2];
    if (u1i>0.0f) {
      u1i = -u1i;
      u2i = -u2i;
      u3i = -u3i;
    }
    return new float[]{u1i,u2i,u3i};
  }


  private int[] setBounds(int n, float[] x) {
    int[] bs = new int[2];
    int n1m = (int)min(x)-5; 
    int n1p = (int)max(x)+5; 
    if(n1m<0){n1m=0;}
    if(n1p>n){n1p=n;}
    bs[0] = n1m;
    bs[1] = n1p;
    return bs;
  }

  public void cpt(FaultCell[] fc, float[][][] fp) {
    int nc = fc.length;
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    float[][][] dt = new float[n3][n2][n1];
    float[][][] pt = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    for (int ic=0; ic<nc; ic++) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      pt[i3][i2][i1] = 1.0f;
      w1[i3][i2][i1] = fc[ic].w1;
      w2[i3][i2][i1] = fc[ic].w2;
      w3[i3][i2][i1] = fc[ic].w3;
    }
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(0.0f,pt,dt,k1,k2,k3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float di = dt[i3][i2][i1];
      if(di<=30f) {
        int j1 = k1[i3][i2][i1];
        int j2 = k2[i3][i2][i1];
        int j3 = k3[i3][i2][i1];
        float d1 = i1-j1;
        float d2 = i2-j2;
        float d3 = i3-j3;
        float u1 = w1[j3][j2][j1];
        float u2 = w2[j3][j2][j1];
        float u3 = w3[j3][j2][j1];
        float ds = abs(u1*d1+u2*d2+u3*d3);
        if(ds<=2.0f) {
          fp[i3][i2][i1] = 1.0f;
        }
      }
    }}}
  }


  private void setKdTreeNodes(
    FaultCell[] fc, float[][] xc, float[][][] fp) {
    float mark = -360f;
    int nc = fc.length;
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    float[][][] dt = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    float[][][] pt = fillfloat(mark,n1,n2,n3);
    for (int ic=0; ic<nc; ic++) {
      int i1 = fc[ic].i1;
      int i2 = fc[ic].i2;
      int i3 = fc[ic].i3;
      xc[0][ic] = i1;
      xc[1][ic] = i2;
      xc[2][ic] = i3;
      pt[i3][i2][i1] = fc[ic].fp;
      w1[i3][i2][i1] = fc[ic].w1;
      w2[i3][i2][i1] = fc[ic].w2;
      w3[i3][i2][i1] = fc[ic].w3;
    }
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(mark,pt,dt,k1,k2,k3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float di = dt[i3][i2][i1];
      if(di<=30f) {
        int j1 = k1[i3][i2][i1];
        int j2 = k2[i3][i2][i1];
        int j3 = k3[i3][i2][i1];
        float d1 = i1-j1;
        float d2 = i2-j2;
        float d3 = i3-j3;
        float u1 = w1[j3][j2][j1];
        float u2 = w2[j3][j2][j1];
        float u3 = w3[j3][j2][j1];
        float ds = abs(u1*d1+u2*d2+u3*d3);
        if(ds<=2.0f) {
          fp[i3][i2][i1] = pt[j3][j2][j1];
        }
      }
    }}}
  }

  private FaultSkin[] resetSkins(FaultSkin[] sks, float[][][] fl) {
    FaultSkinner fsk = new FaultSkinner();
    fsk.setMaxPlanarDistance(1.0);
    fsk.setGrowLikelihoods(0.1,0.2);
    ArrayList<FaultSkin> skinList = new ArrayList<FaultSkin>();
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        int i1 = fc.i1;
        int i2 = fc.i2;
        int i3 = fc.i3;
        fc.ca = null;
        fc.cb = null;
        fc.cl = null;
        fc.cr = null;
        fc.skin = null;
        fc.fl = fl[i3][i2][i1];
      }
      FaultCell[] fcs = sk.getCells();
      FaultSkin[] skin = fsk.findSkins(fcs);
      skinList.add(skin[0]);
    }
    return skinList.toArray(new FaultSkin[0]);
  }

  /**
   * Returns an array of new skins with cells from specified skins. The
   * returned skins may differ from those specified if either cell properties
   * or parameters for this skinner have changed since the cells were last
   * skinned.
   * @param skins array of skins.
   * @return array of new skins.
   */
  /*
  public FaultSkin[] reskin(FaultSkin[] skins) {
    FaultCell[] cells = FaultSkin.getCells(skins);
    for (FaultCell cell:cells)
      cell.skin = null;
    return findSkins(cells);
  }
  */

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors.
   * In these arrays, cells are represented by quads with specified size.
   * @param size the size (in samples) of the quads.
   * @param cells the cells for which to compute quads.
   */
  public static float[][] getXyzUvwRgb(float size, FaultCell[] cells) {
    FloatList xyz = new FloatList();
    FloatList uvw = new FloatList();
    FloatList fcl = new FloatList();
    size *= 0.5f;
    float[] qa = {0.0f,-size,-size};
    float[] qb = {0.0f, size,-size};
    float[] qc = {0.0f, size, size};
    float[] qd = {0.0f,-size, size};
    for (FaultCell cell:cells) {
      float x1 = cell.x1;
      float x2 = cell.x2;
      float x3 = cell.x3;
      float w1 = cell.w1;
      float w2 = cell.w2;
      float w3 = cell.w3;
      float fl = cell.fl;
      float fp = toRadians(cell.fp);
      float ft = toRadians(cell.ft);
      float cp = cos(fp);
      float sp = sin(fp);
      float ct = cos(ft);
      float st = sin(ft);
      float[] ra = rotatePoint(cp,sp,ct,st,qa);
      float[] rb = rotatePoint(cp,sp,ct,st,qb);
      float[] rc = rotatePoint(cp,sp,ct,st,qc);
      float[] rd = rotatePoint(cp,sp,ct,st,qd);
      float a1 = x1+ra[0], a2 = x2+ra[1], a3 = x3+ra[2];
      float b1 = x1+rb[0], b2 = x2+rb[1], b3 = x3+rb[2];
      float c1 = x1+rc[0], c2 = x2+rc[1], c3 = x3+rc[2];
      float d1 = x1+rd[0], d2 = x2+rd[1], d3 = x3+rd[2];
      xyz.add(a3); xyz.add(a2); xyz.add(a1); fcl.add(fl);
      xyz.add(b3); xyz.add(b2); xyz.add(b1); fcl.add(fl);
      xyz.add(c3); xyz.add(c2); xyz.add(c1); fcl.add(fl);
      xyz.add(d3); xyz.add(d2); xyz.add(d1); fcl.add(fl);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
    }
    float[] fc = fcl.trim();
    float fcmin = 0.0f;
    float fcmax = 1.0f;
    ColorMap cmap = new ColorMap(fcmin,fcmax,ColorMap.JET);
    float[] rgb = cmap.getRgbFloats(fc);
    return new float[][]{xyz.trim(),uvw.trim(),rgb};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _fllo; // lower threshold on fault likelihoods
  private float _flhi; // higher threshold on fault likelihoods
  private float _fs1min; // min fault throw
  private float _fs1max; // max fault throw
  private float _dflmax; // max difference between likelihoods of nabors
  private float _dfpmax; // max difference between strikes of nabors
  private float _dftmax; // max difference between dips of nabors
  private float _ds1max; // max difference between throws of nabors
  private float _dnpmax; // max distance to planes of nabors
  private int _ncsmin; // min number of cells that form a skin
  private float _dfpmax1;
  private float _dftmax1;
  private float _dnpmax1;
  private boolean _useOldCells;


  private void updateCells(
    FaultSkin fs, FaultCell[] cells, HashSet<FaultCell> hsc) 
  {
    int d = 3;
    FaultCell[] fcs = FaultSkin.getCells(new FaultSkin[]{fs});
    int nc = fcs.length;
    KdTree kt = setKdTree(cells);
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    for (int ic=0; ic<nc; ++ic) {
      FaultCell fc = fcs[ic];
      float x1 = fc.x1;
      float x2 = fc.x2;
      float x3 = fc.x3;
      float fp = fc.fp;
      xmin[0] = x1-d;
      xmin[1] = x2-d;
      xmin[2] = x3-d;
      xmax[0] = x1+d;
      xmax[1] = x2+d;
      xmax[2] = x3+d;
      int[] id = kt.findInRange(xmin,xmax);
      int nd = id.length;
      if(nd<1){continue;}
      for (int ik=0; ik<nd; ++ik) {
        int ip = id[ik];
        float dp = fp-cells[ip].fp;
        dp = min(abs(dp),abs(dp+360f),abs(dp-360f));
        if(dp<10f){cells[ip].notUsed=false; hsc.remove(cells[ip]);}
      }
    }
    _cells = new FaultCellGrid(hsc.toArray(new FaultCell[0]));
  }

  private KdTree setKdTree(FaultCell[] cells) {
    int nc = cells.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      xc[0][ic] = cells[ic].x1;
      xc[1][ic] = cells[ic].x2;
      xc[2][ic] = cells[ic].x3;
    }
    return new KdTree(xc);
  }

  // Returns skins constructed from specified cells.
  private FaultSkin[] skinsXX(FaultCell[] cells, float[][][] fl) {
    int sk = 0;
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    int ncell = cells.length;
    _cells = new FaultCellGrid(cells);
    HashSet<FaultCell> hsc = new HashSet<FaultCell>();
    for (FaultCell fc:cells) {hsc.add(fc);}

    // Empty list of skins.
    ArrayList<FaultSkin> skinList = new ArrayList<FaultSkin>();

    // Cell comparator for high-to-low ordering based on fault likelihoods.
    Comparator<FaultCell> flComparator = new Comparator<FaultCell>() {
      public int compare(FaultCell c1, FaultCell c2) {
        if (c1.fl<c2.fl)
          return 1;
        else if (c1.fl>c2.fl)
          return -1;
        else
          return 0;
      }
    };

    // Make a list of cells that might be seeds for new skins.
    ArrayList<FaultCell> seedList = new ArrayList<FaultCell>();
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = cells[icell];
      if (cell.fl>=_flhi && cell.s1>=_fs1min && cell.s1<=_fs1max)
        seedList.add(cell);
    }
    //int nseed = 1;
    int nseed = seedList.size();

    // Sort the list of seeds high-to-low by fault likelihood.
    FaultCell[] seeds = seedList.toArray(new FaultCell[0]);
    Arrays.sort(seeds,flComparator);
    seedList.clear();
    for (FaultCell seed:seeds)
      seedList.add(seed);
    //seedList.add(seeds[10000]);

    // While potential seeds remain, ...
    for (int kseed=0; kseed<nseed; ++kseed) {

      // Skip any potential seeds that are already in a skin.
      //while (kseed<nseed && seedList.get(kseed).skin!=null)
      while (kseed<nseed && !seedList.get(kseed).notUsed)
        ++kseed;

      // If we found a seed with which to construct a new skin, ...
      if (kseed<nseed) {
        System.out.println("skinNo="+sk);
        sk++;
        FaultCell seed = seedList.get(kseed);

        _fcg = new FaultCellGrow(getCells(hsc),fl);
        _fcg.setParameters(_dfpmax1,_dftmax1,_dnpmax1);
        _cellsX = new FaultCellGridX(n1,n2,n3);
        // Make a new empty skin.
        FaultSkin skin = new FaultSkin();

        // Make a priority queue of cells, initially with only the seed.
        PriorityQueue<FaultCell> growQueue = 
            new PriorityQueue<FaultCell>(1024,flComparator);
        growQueue.add(seed);

        // While the grow queue is not empty, ...
        int ct=0;
        while (!growQueue.isEmpty()) {
          if(ct%1000==0) {
            System.out.println("ct="+ct);
          }
          ct++;
          // Get and remove the cell with highest fault likelihood from the
          // grow queue. If not already in the skin, add them and link and
          // add any mutually best nabors to the grow queue.
          FaultCell cell = growQueue.poll();
          if (cell.skin==null) {
            for(FaultCell fc:skin) {
              fc.notInQ = true;
            }
            boolean missNabors = true;
            if(_useOldCells) {
              missNabors = findInOldCells(cell);
            }
            if(missNabors) {findInNewCells(cell);}
            FaultCell ca,cb,cl,cr;
            ca = findNaborAbove(cell);
            cb = findNaborBelow(ca);
            if (ca!=null && ca.skin==null && cb==cell) {
              linkAboveBelow(ca,cb);
              growQueue.add(ca);
            }

            cb = findNaborBelow(cell);
            ca = findNaborAbove(cb);
            if (cb!=null && cb.skin==null && ca==cell) {
              linkAboveBelow(ca,cb);
              growQueue.add(cb);
            }

            cl = findNaborLeft(cell);
            cr = findNaborRight(cl);
            if (cl!=null && cl.skin==null && cr==cell) {
              linkLeftRight(cl,cr);
              growQueue.add(cl);
            }

            cr = findNaborRight(cell);
            cl = findNaborLeft(cr);
            if (cr!=null && cr.skin==null && cl==cell) {
              linkLeftRight(cl,cr);
              growQueue.add(cr);
            }
            skin.add(cell);
          }
        }

        // Done growing. Add this skin to the list of skins. Here we include
        // skins that are too small. If we did not include them here, we would
        // need to put them in a list of small skins, so that we could later
        // remove all of their cells. (By not removing those cells now, we
        // prevent them from becoming parts of other skins.) Instead, we
        // simply put all skins in the list, and filter that list later.
        skinList.add(skin);
        updateCells(skin,cells,hsc);
      }
    }

    // Filter skins to include only those that are big enough. Remove all
    // cells from any skins that are too small.
    ArrayList<FaultSkin> bigSkinList = new ArrayList<FaultSkin>();
    for (FaultSkin skin:skinList) {
      if (skin.size()>=_ncsmin) {
        bigSkinList.add(skin);
      } else {
        for (FaultCell cell:skin) {
          cell.skin = null;
          cell.ca = null;
          cell.cb = null;
          cell.cl = null;
          cell.cr = null;
        }
      }
    }
    return bigSkinList.toArray(new FaultSkin[0]);
  }

  private boolean checkPlanar(FaultCell cell) {
    float nc = 0.0f;
    float dp = 0.0f;
    Queue<FaultCell> naborQueue = new LinkedList<FaultCell>();
    naborQueue.add(cell);
    while(!naborQueue.isEmpty()) {
      FaultCell fc = naborQueue.poll();
      FaultCell ca = fc.ca;
      FaultCell cb = fc.cb;
      FaultCell cl = fc.cl;
      FaultCell cr = fc.cr;
      if(ca!=null) {
        float[] ds = distance(ca,cell);
        if(ds[0]<10f) {
          nc += 1.0f;
          dp += ds[1];
          if(ca.notInQ) {
            ca.notInQ = false;
            naborQueue.add(ca);
          }
        }
      }

      if(cb!=null) {
        float[] ds = distance(cb,cell);
        if(ds[0]<10f) {
          nc += 1.0f;
          dp += ds[1];
          if(cb.notInQ) {
            cb.notInQ = false;
            naborQueue.add(cb);
          }
        }
      }

      if(cl!=null) {
        float[] ds = distance(cl,cell);
        if(ds[0]<10f) {
          nc += 1.0f;
          dp += ds[1];
          if(cl.notInQ) {
            cl.notInQ = false;
            naborQueue.add(cl);
          }
        }
      }

      if(cr!=null) {
        float[] ds = distance(cr,cell);
        if(ds[0]<10f) {
          nc += 1.0f;
          dp += ds[1];
          if(cr.notInQ) {
            cr.notInQ = false;
            naborQueue.add(cr);
          }
        }
      }
    }
    
    if(nc==0.0f)  {return true;}
    if(dp/nc<3f*_dnpmax){return true;}
    else          {return false;}
  }

  private float[] distance(FaultCell ca, FaultCell cb) {
    float wa1 = ca.w1;
    float wa2 = ca.w2;
    float wa3 = ca.w3;
    float dx1 = cb.x1-ca.x1;
    float dx2 = cb.x2-ca.x2;
    float dx3 = cb.x3-ca.x3;
    float dnp = abs(wa1*dx1+wa2*dx2+wa3*dx3);
    float dxs = sqrt(dx1*dx1+dx2*dx2+dx3*dx3);
    return new float[]{dxs,dnp};   
  }

  private boolean findInOldCells(FaultCell cell) {
    FaultCell ca = cell.ca;
    FaultCell cb = cell.cb;
    FaultCell cl = cell.cl;
    FaultCell cr = cell.cr;
    if(cl==null) {cl = findNaborLeftOld (cell);}
    if(cr==null) {cr = findNaborRightOld(cell);}
    if(ca==null) {ca = findNaborAboveOld(cell);}
    if(cb==null) {cb = findNaborBelowOld(cell);}
    if(cl==null) {return true;}
    if(cr==null) {return true;}
    if(ca==null) {return true;}
    if(cb==null) {return true;}
    _cellsX.setCell(ca);
    _cellsX.setCell(cb);
    _cellsX.setCell(cl);
    _cellsX.setCell(cr);
    _cellsX.setCell(cell);
    return false;
  }


  private void findInNewCells(FaultCell cell) {
    FaultCell[] fcs = _fcg.applyForCells(cell);
    if(fcs==null) { return;}
    FaultCell ca = cell.ca;
    FaultCell cb = cell.cb;
    FaultCell cl = cell.cl;
    FaultCell cr = cell.cr;
    FaultCellGrid cg = new FaultCellGrid(fcs);
    FaultCell cn = nearestCell(cell,fcs);
    if(cn!=null) {
      cell.x1 = cn.x1; cell.i1 = cn.i1;
      cell.x2 = cn.x2; cell.i2 = cn.i2;
      cell.x3 = cn.x3; cell.i3 = cn.i3;
    }
    _cellsX.setCell(cell);
    if(cl==null){
      cl=findNaborLeft (cg,cell);
      if(cl!=null) {cl.interp=true;}
    }
    if(cr==null){
      cr=findNaborRight(cg,cell);
      if(cr!=null) {cr.interp=true;}
    }
    if(ca==null){
      ca=findNaborAbove(cg,cell);
      if(ca!=null) {ca.interp=true;}
    }
    if(cb==null){
      cb=findNaborBelow(cg,cell);
      if(cb!=null) {cb.interp=true;}
    }
    if(ca!=null) _cellsX.setCell(ca);
    if(cb!=null) _cellsX.setCell(cb);
    if(cl!=null) _cellsX.setCell(cl);
    if(cr!=null) _cellsX.setCell(cr);
  }

  private FaultCell[] getCells(HashSet<FaultCell> hsc) {
    int ic = 0;
    int nc = hsc.size();
    FaultCell[] fcs = new FaultCell[nc];
    for (FaultCell fc:hsc) {
      fcs[ic] = fc;
      ic ++;
    }
    return fcs;
  }

  private FaultCell nearestCell(FaultCell fc, FaultCell[] fcs) {
    int i1 = fc.i1;
    FaultCell cell = null;
    float ds = Float.MAX_VALUE;
    for (FaultCell fci:fcs) {
      if(fci.i1!=i1){continue;}
      float d1 = fc.x1-fci.x1;
      float d2 = fc.x2-fci.x2;
      float d3 = fc.x3-fci.x3;
      float di = d1*d1+d2*d2+d3*d3;
      if(di<ds){ds=di;cell=fci;}
    }
    return cell;
  }


  // Returns true if the specified cells are nabors. This method assumes that
  // all links are mutual. For example, if c1 is the nabor above c2, then c2
  // must be the nabor below c1.
  private static boolean areNabors(FaultCell c1, FaultCell c2) {
    return c1.ca==c2 || c1.cb==c2 || c1.cl==c2 || c1.cr==c2;
  }

  // Methods to link mutually best nabors.
  private void linkAboveBelow(FaultCell ca, FaultCell cb) {
    cb.ca = ca;
    ca.cb = cb;
  }
  private void linkLeftRight(FaultCell cl, FaultCell cr) {
    cr.cl = cl;
    cl.cr = cr;
  }

  // Methods to find good nabors of a specified cell. These methods return
  // null if no nabor is good enough, based on various thresholds.
  private FaultCellGrow _fcg;
  private FaultCellGrid _cells;
  private FaultCellGridX _cellsX;

  private FaultCell findNaborAboveOld(FaultCell cell) {
    FaultCell ca = _cells.findCellAbove(cell);
    return canBeNaborsOld(cell,ca)?ca:null;
  }
  private FaultCell findNaborBelowOld(FaultCell cell) {
    FaultCell cb = _cells.findCellBelow(cell);
    return canBeNaborsOld(cell,cb)?cb:null;
  }
  private FaultCell findNaborLeftOld(FaultCell cell) {
    FaultCell cl = _cells.findCellLeft(cell);
    return canBeNaborsOld(cell,cl)?cl:null;
  }
  private FaultCell findNaborRightOld(FaultCell cell) {
    FaultCell cr = _cells.findCellRight(cell);
    return canBeNaborsOld(cell,cr)?cr:null;
  }


  private FaultCell findNaborAbove(FaultCell cell) {
    FaultCell ca = _cellsX.findCellAbove(cell);
    return canBeNaborsNew(cell,ca)?ca:null;
  }
  private FaultCell findNaborBelow(FaultCell cell) {
    FaultCell cb = _cellsX.findCellBelow(cell);
    return canBeNaborsNew(cell,cb)?cb:null;
  }
  private FaultCell findNaborLeft(FaultCell cell) {
    FaultCell cl = _cellsX.findCellLeft(cell);
    return canBeNaborsNew(cell,cl)?cl:null;
  }
  private FaultCell findNaborRight(FaultCell cell) {
    FaultCell cr = _cellsX.findCellRight(cell);
    return canBeNaborsNew(cell,cr)?cr:null;
  }

  private FaultCell findNaborAbove(FaultCellGrid cells, FaultCell cell) {
    FaultCell ca = cells.findCellAbove(cell);
    return canBeNaborsNew(cell,ca)?ca:null;
  }
  private FaultCell findNaborBelow(FaultCellGrid cells, FaultCell cell) {
    FaultCell cb = cells.findCellBelow(cell);
    return canBeNaborsNew(cell,cb)?cb:null;
  }
  private FaultCell findNaborLeft(FaultCellGrid cells, FaultCell cell) {
    FaultCell cl = cells.findCellLeft(cell);
    return canBeNaborsNew(cell,cl)?cl:null;
  }
  private FaultCell findNaborRight(FaultCellGrid cells, FaultCell cell) {
    FaultCell cr = cells.findCellRight(cell);
    return canBeNaborsNew(cell,cr)?cr:null;
  }


  // Returns true if two specified cells can be nabors. The two cells are
  // assumed to be within one sample of each other. This method uses other
  // attributes of the cells to determine whether or not they can be nabors.
  private boolean canBeNaborsOld(FaultCell ca, FaultCell cb) {
    boolean can = true;
    if (ca==null || cb==null) {
      can = false;
    } else if (minFl(ca,cb)<(_fllo+0.2f)) {
      can = false;
    } else if (absDeltaFl(ca,cb)>_dflmax) {
      can = false;
    } else if (absDeltaFp(ca,cb)>_dfpmax) {
      can = false;
    } else if (absDeltaFt(ca,cb)>_dftmax) {
      can = false;
    } else if (maxDistanceToPlane(ca,cb)>_dnpmax) {
      can = false;
    }
    return can;
  }

  private boolean canBeNaborsNew(FaultCell ca, FaultCell cb) {
    boolean can = true;
    if (ca==null || cb==null) {
      can = false;
    } else if (minFl(ca,cb)<_fllo) {
      can = false;
    } else if (minS1(ca,cb)<_fs1min) {
      can = false;
    } else if (maxS1(ca,cb)>_fs1max) {
      can = false;
    } else if (absDeltaFl(ca,cb)>_dflmax) {
      can = false;
    } else if (absDeltaFp(ca,cb)>_dfpmax) {
      can = false;
    } else if (absDeltaFt(ca,cb)>_dftmax) {
      can = false;
    } else if (absDeltaS1(ca,cb)>_ds1max) {
      can = false;
    } else if (maxDistanceToPlane(ca,cb)>_dnpmax) {
      can = false;
    }
    return can;
  }

  private static float minFl(FaultCell ca, FaultCell cb) {
    return min(ca.fl,cb.fl);
  }
  private static float minS1(FaultCell ca, FaultCell cb) {
    return min(ca.s1,cb.s1);
  }
  private static float maxS1(FaultCell ca, FaultCell cb) {
    return max(ca.s1,cb.s1);
  }
  private static float absDeltaFl(FaultCell ca, FaultCell cb) {
    return abs(ca.fl-cb.fl);
  }
  private static float absDeltaFp(FaultCell ca, FaultCell cb) {
    float del = ca.fp-cb.fp;
    return min(abs(del),abs(del+360.0f),abs(del-360.0f));
  }
  private static float absDeltaFt(FaultCell ca, FaultCell cb) {
    return abs(ca.ft-cb.ft);
  }
  private static float absDeltaS1(FaultCell ca, FaultCell cb) {
    return abs(ca.s1-cb.s1);
  }
  private static float maxDistanceToPlane(FaultCell ca, FaultCell cb) {
    float aw1 = ca.w1, aw2 = ca.w2, aw3 = ca.w3;
    float ax1 = ca.x1, ax2 = ca.x2, ax3 = ca.x3;
    float bw1 = cb.w1, bw2 = cb.w2, bw3 = cb.w3;
    float bx1 = cb.x1, bx2 = cb.x2, bx3 = cb.x3;
    float dx1 = ax1-bx1;
    float dx2 = ax2-bx2;
    float dx3 = ax3-bx3;
    float dab = abs(aw1*dx1+aw2*dx2+aw3*dx3);
    float dba = abs(bw1*dx1+bw2*dx2+bw3*dx3);
    return max(dab,dba);
  }

  // Rotates a specified point by strike (phi) and dip (theta) angles,
  // given specified cosines (cp and ct) and sines (sp and st) of those 
  // angles. The order of transformation is
  // (1) rotate around axis x3 by dip angle
  // (2) rotate around axis x1 by strike angle
  // Returns the coordinates of the rotated point.
  private static float[] rotatePoint(
      float cp, float sp, float ct, float st, float[] x) {
    float x1 = x[0], x2 = x[1], x3 = x[2];
    float y1 =     ct*x1+   st*x2;
    float y2 = -cp*st*x1+cp*ct*x2+sp*x3;
    float y3 =  sp*st*x1-sp*ct*x2+cp*x3;
    return new float[]{y1,y2,y3};
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static class FloatList {
    public int n;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public void add(float[] f) {
      int m = f.length;
      for (int i=0; i<m; ++i)
        add(f[i]);
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }
}
