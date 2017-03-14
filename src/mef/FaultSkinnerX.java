/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.util.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

import static mef.FaultGeometry.*;
import static edu.mines.jtk.util.Parallel.*;
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
    }
  }

  public void getFlpt(int minSize, FaultSkin[] sks, 
    float[][][] fl, float[][][] fp, float[][][] ft) {
    for (FaultSkin sk:sks) {
      if(sk.size()<minSize) {continue;}
      for (FaultCell fc:sk) {
        int i1i = fc.i1;
        int i2i = fc.i2;
        int i3i = fc.i3;
        fl[i3i][i2i][i1i] = fc.fl;
        fp[i3i][i2i][i1i] = fc.fp;
        ft[i3i][i2i][i1i] = fc.ft;
      }
    }
  }


  public void getFp(int minSize, FaultSkin[] sks, float[][][] fp) {
    for (FaultSkin sk:sks) {
      if(sk.size()<minSize) {continue;}
      for (FaultCell fc:sk) {
        int i1i = fc.i1;
        int i2i = fc.i2;
        int i3i = fc.i3;
        fp[i3i][i2i][i1i] = fc.fp;
      }
    }
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

  public void getFaultMask(FaultSkin[] sks, float[][][] fm) {
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        int i1i = fc.i1;
        int i2i = fc.i2;
        int i3i = fc.i3;
        fm[i3i][i2i][i1i] = 1f;
      }
    }
  }



  public void getFl(int minSize, FaultSkin[] sks, float[][][] fl) {
    for (FaultSkin sk:sks) {
      if(sk.size()<minSize) {continue;}
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
        if(i2m>=n2||i2m<0){continue;}
        if(i3m>=n3||i3m<0){continue;}
        if(i2p>=n2||i2p<0){continue;}
        if(i3p>=n3||i3p<0){continue;}
        fl[i3m][i2m][i1m] = fc.fl;
        fl[i3p][i2p][i1p] = fc.fl;
      }
    }
  }

  /**
   * Returns array of cells in ridge surfaces of fault likelihood.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   * @return array of cells.
   */
  public FaultCell[] findCells(float[][][][] flpt) {
    return cells(flpt);
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

  public FaultSkin[] findSkins(int n1, int n2, int n3, FaultCell[] cells) {
    return skins(n1,n2,n3,cells);
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
  private int _w1=10,_w2=_w1,_w3=_w1;
  private float[][][][] _gw;
  private Sampling _sp, _st;

  // Uses fault images to find cells, oriented points located on ridges.
  private FaultCell[] cells(float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

    // Vertical image boundaries are discontinuities that may look like
    // faults. If a fault appears to be near and nearly parallel to image
    // boundaries, then assume it is a boundary artifact and not truly a
    // fault.
    int imax = 5; // max number of samples considered to be near boundary
    float wwmax = 0.75f; // cosine of 30 degrees, squared

    // Loop over all samples. Construct cells for samples nearest to ridges.
    ArrayList<FaultCell> cellList = new ArrayList<FaultCell>();
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmi = f[i3m][i2 ];
        float[] fim = f[i3 ][i2m];
        float[] fip = f[i3 ][i2p];
        float[] fpi = f[i3p][i2 ];
        float[] fmm = f[i3m][i2m];
        float[] fpp = f[i3p][i2p];
        float[] fmp = f[i3m][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fii = f[i3 ][i2 ];
        float[] pii = p[i3 ][i2 ];
        float[] tii = t[i3 ][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          float fmii = fmi[i1 ];
          float fimi = fim[i1 ];
          float fipi = fip[i1 ];
          float fpii = fpi[i1 ];
          float fmmi = fmm[i1 ];
          float fppi = fpp[i1 ];
          float fmpi = fmp[i1 ];
          float fpmi = fpm[i1 ];
          float fiii = fii[i1 ];
          float piii = pii[i1 ];
          float tiii = tii[i1 ];

          // Most image samples will not have a fault cell.
          FaultCell cell = null;

          // Accumulators for ridge likelihoods and locations. Depending on
          // the limits on fault strike used below, we may find more than one
          // ridge.
          float nr = 0;
          float fl = 0.0f;
          float d2 = 0.0f;
          float d3 = 0.0f;

          // If S-N ridge, ...
          if ((fipi<fiii && fimi<fiii) &&
              ((337.5f<=piii || piii<= 22.5f) || 
               (157.5f<=piii && piii<=202.5f))) {
            float f1 = 0.5f*(fipi-fimi); // 1st derivative
            float f2 = fipi-2.0f*fiii+fimi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=_fllo) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w2 = w[1];
              if (imax<=i2 && i2<n2-imax || w2*w2<=wwmax) {
                fl += fr;
                d2 += dr;
                nr += 1;
              }
            }
          }

          // If SW-NE ridge, ...
          if ((fmpi<fiii && fpmi<fiii) &&
              (( 22.5f<=piii && piii<= 67.5f) || 
               (202.5f<=piii && piii<=247.5f))) {
            float f1 = 0.5f*(fmpi-fpmi); // 1st derivative
            float f2 = fmpi-2.0f*fiii+fpmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=_fllo) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w2 = w[1], w3 = w[2];
              if ((imax<=i2 && i2<n2-imax || w2*w2<=wwmax) &&
                  (imax<=i3 && i3<n3-imax || w3*w3<=wwmax)) {
                fl += fr;
                d2 += dr;
                d3 -= dr;
                nr += 1;
              }
            }
          }

          // If W-E ridge, ...
          if ((fpii<fiii && fmii<fiii) &&
              (( 67.5f<=piii && piii<=112.5f) ||
               (247.5f<=piii && piii<=292.5f))) {
            float f1 = 0.5f*(fpii-fmii); // 1st derivative
            float f2 = fmii-2.0f*fiii+fpii; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=_fllo) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w3 = w[2];
              if (imax<=i3 && i3<n3-imax || w3*w3<=wwmax) {
                fl += fr;
                d3 += dr;
                nr += 1;
              }
            }
          }

          // If NW-SE ridge, ...
          if ((fppi<fiii && fmmi<fiii) &&
              ((112.5f<=piii && piii<=157.5f) || 
               (292.5f<=piii && piii<=337.5f))) {
            float f1 = 0.5f*(fppi-fmmi); // 1st derivative
            float f2 = fppi-2.0f*fiii+fmmi; // 2nd derivative
            float dr = -f1/f2; // signed distance to ridge
            float fr = fiii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
            if (fr>=_fllo) {
              float[] w = faultNormalVectorFromStrikeAndDip(piii,tiii);
              float w2 = w[1], w3 = w[2];
              if ((imax<=i2 && i2<n2-imax || w2*w2<=wwmax) &&
                  (imax<=i3 && i3<n3-imax || w3*w3<=wwmax)) {
                fl += fr;
                d2 += dr;
                d3 += dr;
                nr += 1;
              }
            }
          }

          // If at least one ridge, construct a cell and add to list.
          if (nr>0) {
            fl /= nr;
            d2 /= nr;
            d3 /= nr;
            cell = new FaultCell(i1,i2+d2,i3+d3,fl,piii,tiii);
            cellList.add(cell);
          }
        }
      }
    }
    return cellList.toArray(new FaultCell[0]);
  }

  private void updateCells(int n1, int n2, int n3, FaultSkin fs) {
    int d = 3;
    FaultCell[] fcs = FaultSkin.getCells(new FaultSkin[]{fs});
    int nc = fcs.length;
    for (int ic=0; ic<nc; ++ic) {
      FaultCell fc = fcs[ic];
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      int b1 = max(i1-d,0), e1 = min(i1+d,n1-1);
      int b2 = max(i2-d,0), e2 = min(i2+d,n2-1);
      int b3 = max(i3-d,0), e3 = min(i3+d,n3-1);
      for (int k3=b3; k3<=e3; k3++) {
      for (int k2=b2; k2<=e2; k2++) {
      for (int k1=b1; k1<=e1; k1++) {
        FaultCell fcx = _cells[k3][k2][k1];
        if (fcx!=null) {
          float dp = fc.fp-fcx.fp;
          dp = min(abs(dp),abs(dp+360f),abs(dp-360f));
          if(dp<10f){fcx.skin=fs;}
        }
      }}}
    }
  }

  // Returns skins constructed from specified cells.
  private FaultSkin[] skins(int n1, int n2, int n3, FaultCell[] cells) {
    int sk = 0;
    int ncell = cells.length;
    _cells = new FaultCell[n3][n2][n1];
    for (FaultCell cell:cells) {
      int i1 = cell.i1;
      int i2 = cell.i2;
      int i3 = cell.i3;
      _cells[i3][i2][i1] = cell;
    }
    _mask = new short[n3][n2][n1];
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

    // While potential seeds remain, ...
    for (int kseed=0; kseed<nseed; ++kseed) {

      // Skip any potential seeds that are already in a skin.
      while (kseed<nseed && seedList.get(kseed).skin!=null)
        ++kseed;

      // If we found a seed with which to construct a new skin, ...
      if (kseed<nseed) {
        System.out.println("skinNo="+sk);
        sk++;
        FaultCell seed = seedList.get(kseed);

        // Make a new empty skin.
        FaultSkin skin = new FaultSkin();

        // Make a priority queue of cells, initially with only the seed.
        PriorityQueue<FaultCell> growQueue = 
            new PriorityQueue<FaultCell>(1024,flComparator);
        growQueue.add(seed);

        // While the grow queue is not empty, ...
        int ct=0;
        setMask(sk,seed);
        skin.add(seed);
        while (!growQueue.isEmpty()) {
          if(ct%2000==0) {
            System.out.println("ct="+ct);
          }
          ct++;
          // Get and remove the cell with highest fault likelihood from the
          // grow queue. If not already in the skin, add them and link and
          // add any mutually best nabors to the grow queue.
          FaultCell cell = growQueue.poll();
          FaultCell[] fcs = findNabors(n1,n2,n3,cell);
          if(fcs!=null) {
            FaultCellGrid cg = new FaultCellGrid(fcs);
            nearestCell(cell,fcs);
            if(cell!=null){
              FaultCell ca = cell.ca;
              FaultCell cb = cell.cb;
              FaultCell cl = cell.cl;
              FaultCell cr = cell.cr;
              if(ca==null){ 
                ca=findNaborAbove(cg,cell);
                if(ca!=null&&_mask[ca.i3][ca.i2][ca.i1]!=sk) {
                  linkAboveBelow(ca,cell);
                  growQueue.add(ca);
                  setMask(sk,ca);
                  skin.add(ca);
                }
              }
              if(cb==null){ 
                cb=findNaborBelow(cg,cell);
                if(cb!=null&&_mask[cb.i3][cb.i2][cb.i1]!=sk) {
                  linkAboveBelow(cell,cb);
                  growQueue.add(cb);
                  setMask(sk,cb);
                  skin.add(cb);
                }
              }
              if(cl==null){ 
                cl= findNaborLeft(cg,cell);
                if(cl!=null&&_mask[cl.i3][cl.i2][cl.i1]!=sk) {
                  linkLeftRight(cl,cell);
                  growQueue.add(cl);
                  setMask(sk,cl);
                  skin.add(cl);
                }
              }
              if(cr==null){ 
                cr=findNaborRight(cg,cell);
                if(cr!=null&&_mask[cr.i3][cr.i2][cr.i1]!=sk) {
                  linkLeftRight(cell,cr);
                  growQueue.add(cr);
                  setMask(sk,cr);
                  skin.add(cr);
                }
              }
            }
          }
        }

        // Done growing. Add this skin to the list of skins. Here we include
        // skins that are too small. If we did not include them here, we would
        // need to put them in a list of small skins, so that we could later
        // remove all of their cells. (By not removing those cells now, we
        // prevent them from becoming parts of other skins.) Instead, we
        // simply put all skins in the list, and filter that list later.
        skinList.add(skin);
        updateCells(n1,n2,n3,skin);
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


  private void nearestCell(FaultCell cell, FaultCell[] fcs) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    FaultCell cn = null;
    for (FaultCell fci:fcs) {
      if(fci.i1==i1&&fci.i2==i2&&fci.i3==i3){cn=fci;}
    }
    if (cn!=null) {
      cell.x1 = cn.x1; 
      cell.i1 = cn.i1;
      cell.x2 = cn.x2; 
      cell.i2 = cn.i2;
      cell.x3 = cn.x3; 
      cell.i3 = cn.i3;
      cell.fl = cn.fl;
      cell.fp = cn.fp;
      cell.ft = cn.ft;
    } else {
      cell = null;
    }
  }


  private void setMask(int sk, FaultCell cell) {
    int n3 = _mask.length;
    int n2 = _mask[0].length;
    int n1 = _mask[0][0].length;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float w1 = cell.w1;
    float w2 = cell.w2;
    float w3 = cell.w3;
    for (int d=0; d<5; d++) {
      float d1 = d*w1;
      float d2 = d*w2;
      float d3 = d*w3;
      int m1 = round(x1-d1); if(m1<0) {m1=0;} if(m1>=n1){m1=n1-1;} 
      int m2 = round(x2-d2); if(m2<0) {m2=0;} if(m2>=n2){m2=n2-1;} 
      int m3 = round(x3-d3); if(m3<0) {m3=0;} if(m3>=n3){m3=n3-1;} 
      int p1 = round(x1+d1); if(p1<0) {p1=0;} if(p1>=n1){p1=n1-1;} 
      int p2 = round(x2+d2); if(p2<0) {p2=0;} if(p2>=n2){p2=n2-1;} 
      int p3 = round(x3+d3); if(p3<0) {p3=0;} if(p3>=n3){p3=n3-1;} 
      _mask[m3][m2][m1] = (byte)sk;
      _mask[p3][p2][p1] = (byte)sk;
    }

  }

  public float[][][][] setGaussWeights(
    float sigu, float sigw, Sampling sp, Sampling st) 
  {
    _sp = sp;
    _st = st;
    _gw = gaussWeights(sigu,sigw);
    return _gw;
  }


  public float[][][][] setGaussWeights(Sampling sp, Sampling st) {
    _sp = sp;
    _st = st;
    _gw = gaussWeights(10,2);
    return _gw;
  }

  public FaultCell[] findNaborsX(
    final int n1, final int n2, final int n3, FaultCell cell) 
  { 
    float fp = cell.fp;
    float ft = cell.ft;
    final int i1 = cell.i1;
    final int i2 = cell.i2;
    final int i3 = cell.i3;
    float sp2 = abs(cell.w2/cell.w1);
    float sp3 = abs(cell.w3/cell.w1);
    int w2 = min(round(_w1*sp3),_w2);
    int w3 = min(round(_w1*sp2),_w3);
    int b1 = max(i1-_w1,0), e1 = min(i1+_w1,n1-1);
    int b2 = max(i2- w2,0), e2 = min(i2+ w2,n2-1);
    int b3 = max(i3- w3,0), e3 = min(i3+ w3,n3-1);
    final FaultCell[] fcs = 
      findCellsInBoxX(b1,e1,b2,e2,b3,e3,fp,ft);
    final int nc = fcs.length;
    if (nc<10) {return null;}
    final int m1=5;
    int t2=min(max(round(m1*sp2),m1),_w2*2);
    int t3=min(max(round(m1*sp3),m1),_w3*2);
    if(t2%2==0) {t2+=1;};
    if(t3%2==0) {t3+=1;};
    final int m2=t2;
    final int m3=t3;
    final int s1=i1-(m1-1)/2;
    final int s2=i2-(m2-1)/2;
    final int s3=i3-(m3-1)/2;
    float sig = 10f;
    final float cs = 0.5f/(sig*sig);
    final float[][][] fls = new float[m3][m2][m1];
    final float[][][][] flpt = new float[3][m3][m2][m1];
    loop(m3,new Parallel.LoopInt(){
    public void compute(int k3) {
      int p3 = k3+s3;
      float[][] fs3 = fls[k3];
      float[][] fl3 = flpt[0][k3];
      float[][] fp3 = flpt[1][k3];
      float[][] ft3 = flpt[2][k3];
      for (int k2=0; k2<m2; ++k2) {
        int p2 = k2+s2;
      for (int k1=0; k1<m1; ++k1) {
        int p1 = k1+s1;
        float fps = 0f;
        float fts = 0f;
        float gws = 0f;
        float gls = 0f;
        for (int ic=0; ic<nc; ++ic) {
          FaultCell fc = fcs[ic];
          float w1 = fc.w1;
          float w2 = fc.w2;
          float w3 = fc.w3;
          float r1 = p1-fc.x1;
          float r2 = p2-fc.x2;
          float r3 = p3-fc.x3;
          float rs = r1*r1+r2*r2+r3*r3;
          float gw = exp(-rs*cs);
          if(rs!=0f) {
            rs = 1f/sqrt(rs);
            r1 *= rs; r2 *= rs; r3 *= rs;
            float wr = w1*r1+w2*r2+w3*r3;
            wr = 1-wr*wr;
            wr *= wr;
            wr *= wr;
            wr *= wr;
            gw *= wr;
          }
          float gl = gw*fc.fl;
          gws += gw;
          gls += gl;
          fps += gl*fc.fp; // fix this
          fts += gl*fc.ft;
        }
        float fpi = fps/gls; 
        fl3[k2][k1] = gls;
        fp3[k2][k1] = fpi;
        ft3[k2][k1] = fts/gls;
        fs3[k2][k1] = gls/gws;
      }}
    }});
    sub(flpt[0],min(flpt[0]),flpt[0]);
    div(flpt[0],max(flpt[0]),flpt[0]);
    FaultCell[] cells = cells(flpt);
    ArrayList<FaultCell> fca = new ArrayList<FaultCell>();
    for (FaultCell fc:cells) {
      int k1 = fc.i1;
      int k2 = fc.i2;
      int k3 = fc.i3;
      float x1 = fc.x1+s1;
      float x2 = fc.x2+s2;
      float x3 = fc.x3+s3;
      int p1 = round(x1);
      int p2 = round(x2);
      int p3 = round(x3);
      if(k1<0||k1>=m1) {continue;}
      if(k2<0||k2>=m2) {continue;}
      if(k3<0||k3>=m3) {continue;}
      if(p1<0||p1>=n1) {continue;}
      if(p2<0||p2>=n2) {continue;}
      if(p3<0||p3>=n3) {continue;}
      float fpi = fc.fp;
      float fti = fc.ft;
      float fli = fls[k3][k2][k1];
      if(fli>_fllo)
        fca.add(new FaultCell(x1,x2,x3,fli,fpi,fti));
    }
    return fca.toArray(new FaultCell[0]);
  }

  private FaultCell[] findCellsInBoxX(
    int b1, int e1, int b2, int e2, int b3, int e3, float fp, float ft)
  {
    float dpm = (float)(3*_sp.getDelta());
    float dtm = (float)(3*_st.getDelta());
    ArrayList<FaultCell> fca = new ArrayList<FaultCell>();
    for (int i3=b3; i3<=e3; ++i3) {
    for (int i2=b2; i2<=e2; ++i2) {
    for (int i1=b1; i1<=e1; ++i1) {
      FaultCell fc = _cells[i3][i2][i1];
      if (fc!=null && fc.skin==null) {
        float fpi = fc.fp;
        float dti = abs(fc.ft-ft);
        if(fpi-fp>360-dpm) {fpi = 360-fpi;}
        if(fpi-fp<dpm-360) {fpi = 360+fpi;}
        float dpi = abs(fpi-fp);
        if (dpi<dpm && dti<dtm){fca.add(fc);}
      }
    }}}
    return fca.toArray(new FaultCell[0]);
  }



  public FaultCell[] findNabors(
    final int n1, final int n2, final int n3, FaultCell cell) 
  { 
    float fp = cell.fp;
    float ft = cell.ft;
    final int i1 = cell.i1;
    final int i2 = cell.i2;
    final int i3 = cell.i3;
    float sp2 = abs(cell.w2/cell.w1);
    float sp3 = abs(cell.w3/cell.w1);
    int w2 = min(round(_w1*sp3),_w2);
    int w3 = min(round(_w1*sp2),_w3);
    int b1 = max(i1-_w1,0), e1 = min(i1+_w1,n1-1);
    int b2 = max(i2- w2,0), e2 = min(i2+ w2,n2-1);
    int b3 = max(i3- w3,0), e3 = min(i3+ w3,n3-1);
    final int[][] xas = new int[3][];
    final float[][] fas = new float[3][];
    final int[] ida = 
      findCellsInBox(b1,e1,b2,e2,b3,e3,fp,ft,fas,xas);
    final int[] x1a = xas[0];
    final int[] x2a = xas[1];
    final int[] x3a = xas[2];
    final float[] fla = fas[0];
    final float[] fpa = fas[1];
    final float[] fta = fas[2];
    final int nc = ida.length;
    if (nc<10) {return null;}
    final int m1=5;
    int t2=min(max(round(m1*sp2),m1),_w2*2);
    int t3=min(max(round(m1*sp3),m1),_w3*2);
    if(t2%2==0) {t2+=1;};
    if(t3%2==0) {t3+=1;};
    final int m2=t2;
    final int m3=t3;
    final int s1=i1-(m1-1)/2;
    final int s2=i2-(m2-1)/2;
    final int s3=i3-(m3-1)/2;
    final float[][][] fls = new float[m3][m2][m1];
    final float[][][][] flpt = new float[3][m3][m2][m1];
    loop(m3,new Parallel.LoopInt(){
    public void compute(int k3) {
      int p3 = k3+s3;
      float[][] fs3 = fls[k3];
      float[][] fl3 = flpt[0][k3];
      float[][] fp3 = flpt[1][k3];
      float[][] ft3 = flpt[2][k3];
      for (int k2=0; k2<m2; ++k2) {
        int p2 = k2+s2;
      for (int k1=0; k1<m1; ++k1) {
        int p1 = k1+s1;
        float fps = 0f;
        float fts = 0f;
        float gws = 0f;
        float gls = 0f;
        for (int ic=0; ic<nc; ++ic) {
          int r1 = p1-x1a[ic];
          int r2 = p2-x2a[ic];
          int r3 = p3-x3a[ic];
          float gw = _gw[ida[ic]][r3][r2][r1];
          float gl = gw*fla[ic];
          gws += gw;
          gls += gl;
          fps += gl*fpa[ic];
          fts += gl*fta[ic];
        }
        float fpi = fps/gls; 
        if(fpi>360f){fpi-=360f;}
        fl3[k2][k1] = gls;
        fp3[k2][k1] = fpi;
        ft3[k2][k1] = fts/gls;
        fs3[k2][k1] = gls/gws;
      }}
    }});
    sub(flpt[0],min(flpt[0]),flpt[0]);
    div(flpt[0],max(flpt[0]),flpt[0]);
    FaultCell[] cells = cells(flpt);
    ArrayList<FaultCell> fca = new ArrayList<FaultCell>();
    for (FaultCell fc:cells) {
      int k1 = fc.i1;
      int k2 = fc.i2;
      int k3 = fc.i3;
      float x1 = fc.x1+s1;
      float x2 = fc.x2+s2;
      float x3 = fc.x3+s3;
      int p1 = round(x1);
      int p2 = round(x2);
      int p3 = round(x3);
      if(k1<0||k1>=m1) {continue;}
      if(k2<0||k2>=m2) {continue;}
      if(k3<0||k3>=m3) {continue;}
      if(p1<0||p1>=n1) {continue;}
      if(p2<0||p2>=n2) {continue;}
      if(p3<0||p3>=n3) {continue;}
      float fpi = fc.fp;
      float fti = fc.ft;
      float fli = fls[k3][k2][k1];
      if(fli>_fllo)
        fca.add(new FaultCell(x1,x2,x3,fli,fpi,fti));
    }
    return fca.toArray(new FaultCell[0]);
  }



  private int[] findCellsInBox(
    final int b1, final int e1, final int b2, 
    final int e2, final int b3, final int e3, 
    final float fp, final float ft, final float[][] fas, final int[][] xas)
  {
    float dpm = (float)(3*_sp.getDelta());
    float dtm = (float)(3*_st.getDelta());
    final ArrayList<Float> fpa = new ArrayList<Float>();
    final ArrayList<FaultCell> fca = new ArrayList<FaultCell>();
    for (int i3=b3; i3<=e3; ++i3) {
    for (int i2=b2; i2<=e2; ++i2) {
    for (int i1=b1; i1<=e1; ++i1) {
      FaultCell fc = _cells[i3][i2][i1];
      if (fc!=null && fc.skin==null) {
        float fpi = fc.fp;
        float fti = fc.ft;
        float dti = abs(fti-ft);
        if(fpi-fp>360-dpm) {fpi = 360-fpi;}
        if(fpi-fp<dpm-360) {fpi = 360+fpi;}
        float dpi = abs(fpi-fp);
        if (dpi<dpm && dti<dtm){fca.add(fc);fpa.add(fpi);}
      }
    }}}
    final int d1 = _w1*3;
    final int d2 = _w2*3;
    final int d3 = _w3*3;
    final int nt = _st.getCount();
    final int nc = fca.size();
    xas[0] = new int[nc];
    xas[1] = new int[nc];
    xas[2] = new int[nc];
    fas[0] = new float[nc];
    fas[1] = new float[nc];
    fas[2] = new float[nc];
    final int[] ids = new int[nc];
    if (nc>0) {
      loop(nc,new Parallel.LoopInt(){
      public void compute(int ic) {
        FaultCell fc = fca.get(ic);
        float fti = fc.ft;
        float fli = fc.fl;
        float fpi = fpa.get(ic);
        int it = _st.indexOfNearest(fti);
        int ip = _sp.indexOfNearest(fc.fp);
        ids[ic] = ip*nt+it;
        fas[0][ic] = fli;
        fas[2][ic] = fti;
        fas[1][ic] = fpi;
        xas[0][ic] = fc.i1-d1;
        xas[1][ic] = fc.i2-d2;
        xas[2][ic] = fc.i3-d3;
      }});
    }
    return ids;
  }


  public float[][][][] gaussWeights(
    float sigu, float sigw) 
  {
    final int np = _sp.getCount();
    final int nt = _st.getCount();
    final int c1 = _w1*3;
    final int c2 = _w2*3;
    final int c3 = _w3*3;
    final int n1 = _w1*6+1;
    final int n2 = _w2*6+1;
    final int n3 = _w3*6+1;
    final float sw = 0.25f/(sigw*sigw);
    final float su = 0.25f/(sigu*sigu);
    final float sv = su;
    final float[][][][] gws = new float[np*nt][n3][n2][n1];
    for (int ip=0; ip<np; ++ip) {
    for (int it=0; it<nt; ++it) {
      float fp = (float)_sp.getValue(ip);
      float ft = (float)_st.getValue(it);
      final float[][][] gwi = gws[ip*nt+it];
      float[] u = faultDipVectorFromStrikeAndDip(fp,ft);
      float[] v = faultStrikeVectorFromStrikeAndDip(fp,ft);
      float[] w = faultNormalVectorFromStrikeAndDip(fp,ft);
      float u1 = u[0];
      float u2 = u[1];
      float u3 = u[2];
      float v1 = v[0];
      float v2 = v[1];
      float v3 = v[2];
      float w1 = w[0];
      float w2 = w[1];
      float w3 = w[2];
      final float w11 = w1*w1;
      final float w12 = w1*w2;
      final float w13 = w1*w3;
      final float w22 = w2*w2;
      final float w23 = w2*w3;
      final float w33 = w3*w3;
      final float v11 = v1*v1;
      final float v12 = v1*v2;
      final float v13 = v1*v3;
      final float v22 = v2*v2;
      final float v23 = v2*v3;
      final float v33 = v3*v3;
      final float u11 = u1*u1;
      final float u12 = u1*u2;
      final float u13 = u1*u3;
      final float u22 = u2*u2;
      final float u23 = u2*u3;
      final float u33 = u3*u3;
      loop(n3,new LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float dx1 = i1-c1;
          float dx2 = i2-c2;
          float dx3 = i3-c3;
          float d11 = dx1*dx1;
          float d22 = dx2*dx2;
          float d33 = dx3*dx3;
          float d12 = dx1*dx2;
          float d13 = dx1*dx3;
          float d23 = dx2*dx3;

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
          gss += (wd1+wd2+wd3+wds)*sw;
          gss += (ud1+ud2+ud3+uds)*su;
          gss += (vd1+vd2+vd3+vds)*sv;
          gwi[i3][i2][i1] = exp(-gss);
        }}
      }});
    }}
    return gws;
  }


  // Returns true if the specified cells are nabors. This method assumes that
  // all links are mutual. For example, if c1 is the nabor above c2, then c2
  // must be the nabor below c1.

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
  private FaultCell[][][] _cells;
  private short[][][] _mask;


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
