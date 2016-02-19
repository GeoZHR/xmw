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

import util.*;

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
   * Returns array of cells in ridge surfaces of fault likelihood.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   * @return array of cells.
   */
  public FaultCell[] findCells(float[][][][] flpt) {
    return cells(0,0,0,flpt);
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
  private int _w1=10,_w2=10,_w3=10;
  private float[][][][][] _gw;
  private Sampling _sp, _st;

  // Uses fault images to find cells, oriented points located on ridges.
  private FaultCell[] cells(int s1, int s2, int s3, float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

    // Smooth fault likelihoods in 2nd and 3rd dimensions. This helps to
    // eliminate spurious ridges, and improves the accuracy of 2nd-order
    // finite-difference approximations (parabolic interpolation) used to
    // locate ridges.
    /*
    float[][][] fs = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
    rgf.applyX0X(f,fs);
    rgf.applyXX0(fs,fs);
    f = fs;
    */

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
            cell = new FaultCell(i1+s1,i2+d2+s2,i3+d3+s3,fl,piii,tiii);
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
        FaultCell fcx = _cells.get(k1,k2,k3);
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
    _cells = new FaultCellGrid(n1,n2,n3,cells);

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
      while (kseed<nseed && seedList.get(kseed).skin!=null)
        ++kseed;

      // If we found a seed with which to construct a new skin, ...
      if (kseed<nseed) {
        System.out.println("skinNo="+sk);
        sk++;
        FaultCell seed = seedList.get(kseed);

        _cellsX = new FaultCellGrid(n1,n2,n3);
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
            findInNewCells(n1,n2,n3,cell);
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
        System.out.println("skinSize="+skin.size());
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


  private void findInNewCells(int n1, int n2, int n3, FaultCell cell) {
    FaultCell[] fcs = null;//findNabors(n1,n2,n3,cell);
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
      cell.fl = cn.fl;
      cell.fp = cn.fp;
      cell.ft = cn.ft;
    }
    _cellsX.set(cell);
    if(cl==null){
      cl=findNaborLeft(cg,cell);
    }
    if(cr==null){
      cr=findNaborRight(cg,cell);
    }
    if(ca==null){
      ca=findNaborAbove(cg,cell);
    }
    if(cb==null){
      cb=findNaborBelow(cg,cell);
    }
    if(ca!=null) _cellsX.set(ca);
    if(cb!=null) _cellsX.set(cb);
    if(cl!=null) _cellsX.set(cl);
    if(cr!=null) _cellsX.set(cr);
  }

  public void setCellGrid(FaultCell[] cells) {
    _cells = new FaultCellGrid(cells);
  }

  public float[][][][][] setGaussWeights(Sampling sp, Sampling st) {
    _sp = sp;
    _st = st;
    _gw = gaussWeights(20,1);
    return _gw;
  }

  public FaultCell[] findNabors(
    int n1, int n2, int n3, final float[][][] fls, FaultCell cell) 
  {
    float fp = cell.fp;
    float ft = cell.ft;
    final int i1 = cell.i1;
    final int i2 = cell.i2;
    final int i3 = cell.i3;
    int b1 = max(i1-_w1,0), e1 = min(i1+_w1,n1-1);
    int b2 = max(i2-_w2,0), e2 = min(i2+_w2,n2-1);
    int b3 = max(i3-_w3,0), e3 = min(i3+_w3,n3-1);
    final FaultCell[] fcs = findCellsInBox(b1,e1,b2,e2,b3,e3,fp,ft);
    final int nc = fcs.length;
    if (nc<5) {return null;}
    final int m1=25;
    final int m2=25;//min(max(round(m1*abs(cell.w2/cell.w1)),m1),_w2*3);
    final int m3=25;//min(max(round(m1*abs(cell.w3/cell.w1)),m1),_w3*3);
    final int s1=i1-(m1-1)/2;
    final int s2=i2-(m2-1)/2;
    final int s3=i3-(m3-1)/2;
    final float[][][][] flpt = new float[3][m3][m2][m1];
    loop(m3,new Parallel.LoopInt(){
    public void compute(int k3) {
      int p3 = k3+s3;
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
        for (FaultCell fc:fcs) {
          float fp = fc.fp;
          float ft = fc.ft;
          int r1 = p1-fc.i1+_w1*3;
          int r2 = p2-fc.i2+_w2*3;
          int r3 = p3-fc.i3+_w3*3;
          int ip = _sp.indexOfNearest(fp);
          int it = _st.indexOfNearest(ft);
          float gw = _gw[ip][it][r3][r2][r1]*fc.fl;
          gws += gw;
          fps += fp*gw;
          fts += ft*gw;
        }
        fl3[k2][k1] = gws;
        ft3[k2][k1] = fts/gws;
        fp3[k2][k1] = fps/gws;
        fls[p3][p2][p1] = gws;
      }}
    }});
    sub(flpt[0],min(flpt[0]),flpt[0]);
    div(flpt[0],max(flpt[0]),flpt[0]);
    sub(fls,min(fls),fls);
    div(fls,max(fls),fls);

    return cells(s1,s2,s3,flpt);
  }


  private FaultCell[] findCellsInBox(
    int b1, int e1, int b2, int e2, int b3, int e3, float fp, float ft)
  {
    ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
    for (int i3=b3; i3<=e3; ++i3) {
    for (int i2=b2; i2<=e2; ++i2) {
    for (int i1=b1; i1<=e1; ++i1) {
      FaultCell cell = _cells.get(i1,i2,i3);
      if (cell!=null) {
        float dp = cell.fp-fp;
        float dt = abs(cell.ft-ft);
        dp = min(abs(dp),abs(dp+360.0f),abs(dp-360.0f));
        if (dp<20f && dt<15f) {fcl.add(cell);}
      }
    }}}
    return fcl.toArray(new FaultCell[0]);
  }



  public float[][][][][] gaussWeights(
    float sigu, float sigw) 
  {
    int n1 = _w1*6+1;
    int n2 = _w2*6+1;
    int n3 = _w3*6+1;
    int np = _sp.getCount();
    int nt = _st.getCount();
    float[][][] gw = new float[n1][n2][n3];
    final float[][][][][] gws = new float[np][nt][n3][n2][n1];
    gw[_w3*3][_w2*3][_w1*3] = 1f;
    RecursiveGaussianFilterP rgfu = new RecursiveGaussianFilterP(sigu);
    RecursiveGaussianFilterP rgfw = new RecursiveGaussianFilterP(sigw);
    rgfu.apply0XX(gw,gw);
    rgfw.applyX0X(gw,gw);
    rgfu.applyXX0(gw,gw);
    for (int ip=0; ip<np; ++ip) {
      final float phi = (float)_sp.getValue(ip);
      Rotator r = new Rotator(phi,n1,n2,n3);
      float[][][] rgw = r.rotate(gw);
      shearTheta(_st,rgw,gws[ip]);
    }
    return gws;
  }

  private void shearTheta(
    Sampling thetaSampling, final float[][][] fx, final float[][][][] gwp) 
  {
    final int n3 = n3(fx);
    final int n2 = n2(fx);
    final int n1 = n1(fx);
    final int m3 = gwp[0].length;
    final int m2 = gwp[0][0].length;
    final int m1 = gwp[0][0][0].length;
    final float[][][] ft = new float[n3][n2][n1];
    final Sampling st = thetaSampling;
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    int nt = st.getCount();
    for (int it=0; it<nt; ++it) {
      float ti = (float)st.getValue(it);
      float theta = toRadians(ti);
      final float shear = -1.0f/tan(theta);
      loop(n2,new LoopInt() {
      public void compute(int i2) {
        float[][] fx2 = extractSlice2(i2,fx);
        if (fx2==null)return;
        int n3 = fx2.length;
        float[][] sns = shear(si,shear,fx2);
        for (int i3=0,j3=i3lo(i2,fx); i3<n3; ++i3,++j3) {
          ft[j3][i2] = sns[i3];
        }
      }});
      int[] id = new int[3];
      max(ft,id);
      float[][][] gw = gwp[it];
      for (int i3=0,k3=id[2]-_w3*3; i3<m3; ++i3,++k3) {
        if(k3>=n3||k3<0) {continue;}
        float[][] ft3 = ft[k3];
        float[][] gw3 = gw[i3];
      for (int i2=0,k2=id[1]-_w2*3; i2<m2; ++i2,++k2) {
        if(k2>=n2||k2<0) {continue;}
        float[] ft32 = ft3[k2];
        float[] gw32 = gw3[i2];
      for (int i1=0,k1=id[0]-_w1*3; i1<m1; ++i1,++k1) {
        if(k1>=n1||k1<0) {continue;}
        gw32[i1] = ft32[k1];
      }}}
    }
  }

    // Makes an array like that specified, including any null arrays.
  private float[][][] like(float[][][] p) {
    int n1 = n1(p);
    int n2 = n2(p);
    int n3 = n3(p);
    int np = p.length;
    float[][][] q = new float[n3][n2][];
    for (int ip=0; ip<np; ++ip) {
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          q[i3][i2] = (p[i3][i2]!=null)?new float[n1]:null;
        }
      }
    }
    return q;
  }


    // Numbers of samples in 3D arrays (arrays of arrays of arrays),
  // which after rotation may contain some null arrays.
  private static int n1(float[][][] f) {
    int n1 = 0;
    int n2 = f[0].length;
    int n3 = f.length;
    for (int i3=0; i3<n3 && n1==0; ++i3) {
      for (int i2=0; i2<n2 && n1==0; ++i2) {
        if (f[i3][i2]!=null)
          n1 = f[i3][i2].length;
      }
    }
    return n1;
  }
  private static int n2(float[][][] f) {
    return f[0].length;
  }
  private static int n3(float[][][] f) {
    return f.length;
  }
  private static int n1(float[][][][] f) {
    return n1(f[0]);
  }
  private static int n2(float[][][][] f) {
    return n2(f[0]);
  }
  private static int n3(float[][][][] f) {
    return n3(f[0]);
  }

  // Get/set non-null slices of rotated 3D arrays
  private static float[][] extractSlice2(int i2, float[][][] x) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i3lo = i3lo(i2,x);
    int i3hi = i3hi(i2,x);
    int m3 = 1+i3hi-i3lo;
    float[][] x2 = (m3>0)?new float[m3][n1]:null;
    for (int i3=0; i3<m3; ++i3)
      copy(x[i3+i3lo][i2],x2[i3]);
    return x2;
  }
  private static void restoreSlice2(int i2, float[][][] x, float[][] x2) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i3lo = i3lo(i2,x);
    int i3hi = i3hi(i2,x);
    int m3 = 1+i3hi-i3lo;
    assert x2.length==m3:"x2 length is correct";
    for (int i3=0; i3<m3; ++i3)
      copy(x2[i3],x[i3+i3lo][i2]);
  }


  private static int i2lo(int i3, float[][][] x) {
    int n2 = x[0].length;
    int i2lo = 0;
    while (i2lo<n2 && x[i3][i2lo]==null)
      ++i2lo;
    return i2lo;
  }
  private static int i2hi(int i3, float[][][] x) {
    int n2 = x[0].length;
    int i2hi = n2-1;
    while (i2hi>=0 && x[i3][i2hi]==null)
      --i2hi;
    return i2hi;
  }
  private static int i3lo(int i2, float[][][] x) {
    int n3 = x.length;
    int i3lo = 0;
    while (i3lo<n3 && x[i3lo][i2]==null)
      ++i3lo;
    return i3lo;
  }
  private static int i3hi(int i2, float[][][] x) {
    int n3 = x.length;
    int i3hi = n3-1;
    while (i3hi>=0 && x[i3hi][i2]==null)
      --i3hi;
    return i3hi;
  }

  // Shear horizontally such that q(i1,i2) = p(i1,i2+s*i1).
  private static float[][] shear(
    SincInterpolator si, double s, float[][] p)
  {
    int n1 = p[0].length;
    int n2p = p.length;
    int n2q = n2p+(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] q = new float[n2q][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2p; ++i2)
        pp[i2] = p[i2][i1];
      double f2q = (s<0.0f)?s*i1:s*i1-dqp;
      si.interpolate(n2p,1.0,0.0,pp,n2q,1.0f,f2q,qq);
      for (int i2=0; i2<n2q; ++i2)
        q[i2][i1] = qq[i2];
    }
    return q;
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
  private FaultCellGrid _cells;
  private FaultCellGrid _cellsX;


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




    ///////////////////////////////////////////////////////////////////////////
  // image rotator

  private static class Rotator {

    Rotator(double phi, int n1, int n2, int n3) {
      _n1 = n1;

      // angle phi in radians, cosine and sine
      _phir = toRadians(phi);
      _cosp = cos(_phir);
      _sinp = sin(_phir);

      // center of rotation
      _x2c = 0.5*(n2-1.0);
      _x3c = 0.5*(n3-1.0);

      // input sampling
      _s2p = new Sampling(n2,1.0,0.0);
      _s3p = new Sampling(n3,1.0,0.0);

      // corners of input sampling rectangle
      double[] x2s = { 0.0, 0.0,n2-1,n2-1};
      double[] x3s = { 0.0,n3-1,n3-1, 0.0};

      // bounds after rotation
      double x2min =  Double.MAX_VALUE;
      double x3min =  Double.MAX_VALUE;
      double x2max = -Double.MAX_VALUE;
      double x3max = -Double.MAX_VALUE;
      for (int i=0; i<4; ++i) {
        double x2q = x2q(x2s[i],x3s[i]);
        double x3q = x3q(x2s[i],x3s[i]);
        if (x2q<x2min) x2min = x2q;
        if (x2q>x2max) x2max = x2q;
        if (x3q<x3min) x3min = x3q;
        if (x3q>x3max) x3max = x3q;
      }
      x2min = floor(x2min);
      x2max = ceil(x2max);
      x3min = floor(x3min);
      x3max = ceil(x3max);

      // sampling after rotation
      int n2q = max(2,1+(int)(x2max-x2min+0.5));
      int n3q = max(2,1+(int)(x3max-x3min+0.5));
      double d2q = 1.0;
      double d3q = 1.0;
      double f2q = x2min;
      double f3q = x3min;
      _s2q = new Sampling(n2q,d2q,f2q);
      _s3q = new Sampling(n3q,d3q,f3q);
      //trace("s2p: n2p="+n2);
      //trace("s3p: n3p="+n3);
      //trace("s2q: n2q="+n2q+" d2q="+d2q+" f2q="+f2q);
      //trace("s3q: n3q="+n3q+" d3q="+d3q+" f3q="+f3q);
    }

    float[][][] rotate(float[][][] p) {
      final float[][][] fp = p;
      final float[][] siTable = _siTable;
      final int nsinc = siTable.length;
      final int lsinc = siTable[0].length;
      final Sampling s2q = _s2q;
      final Sampling s3q = _s3q;
      final int n1 = _n1;
      final int n2p = _s2p.getCount();
      final int n3p = _s3p.getCount();
      final int n2q = _s2q.getCount();
      final int n3q = _s3q.getCount();
      final float[][][] q = new float[n3q][n2q][];
      loop(n3q,new LoopInt() {
        public void compute(int i3) {
          double x3q = s3q.getValue(i3);
          for (int i2=0; i2<n2q; ++i2) {
            double x2q = s2q.getValue(i2);
            double x2p = x2p(x2q,x3q);
            double x3p = x3p(x2q,x3q);
            if (inBounds(x2p,x3p)) {
              float[] q32 = q[i3][i2] = new float[n1];
              int i2p = (int)floor(x2p);
              int i3p = (int)floor(x3p);
              double f2p = x2p-i2p;
              double f3p = x3p-i3p;
              int k2p = (int)(f2p*(nsinc-1)+0.5);
              int k3p = (int)(f3p*(nsinc-1)+0.5);
              for (int k3s=0; k3s<lsinc; ++k3s) {
                float s3 = siTable[k3p][k3s];
                int j3p = i3p+k3s-lsinc/2+1;
                if (j3p<   0) j3p = 0;
                if (j3p>=n3p) j3p = n3p-1;
                for (int k2s=0; k2s<lsinc; ++k2s) {
                  float s2 = siTable[k2p][k2s];
                  int j2p = i2p+k2s-lsinc/2+1;
                  if (j2p<   0) j2p = 0;
                  if (j2p>=n2p) j2p = n2p-1;
                  float[] p32 = fp[j3p][j2p];
                  float s32 = s3*s2;
                  for (int i1=0; i1<n1; ++i1)
                    q32[i1] += p32[i1]*s32;
                }
              }
            }
          }
        }
      });
      return q;
    }

    /////////////////////////////////////////////////////////////////////////
    // private

    private int _n1; // number of samples in 1st dimension
    private double _phir,_cosp,_sinp; // angle phi in radians, cosine, sine
    private double _x2c,_x3c; // coordinates of center of rotation
    private Sampling _s2p,_s3p; // samplings in original coordinates
    private Sampling _s2q,_s3q; // samplings in rotated coordinates
    private static float[][] _siTable; // sinc interpolation coefficients
    private static int HALF_LSINC; // half length of sinc interpolator
    static {
      SincInterpolator si = new SincInterpolator();
      _siTable = si.getTable();
      HALF_LSINC = _siTable[0].length/2;
    }
    private double x2p(double x2q, double x3q) {
      return _x2c+(x2q-_x2c)*_sinp-(x3q-_x3c)*_cosp;
    }
    private double x3p(double x2q, double x3q) {
      return _x3c+(x2q-_x2c)*_cosp+(x3q-_x3c)*_sinp;
    }
    private double x2q(double x2p, double x3p) {
      return _x2c+(x2p-_x2c)*_sinp+(x3p-_x3c)*_cosp;
    }
    private double x3q(double x2p, double x3p) {
      return _x3c-(x2p-_x2c)*_cosp+(x3p-_x3c)*_sinp;
    }
    /*
    private double x2p(double x2q, double x3q) {
      return _x2c+(x2q-_x2c)*_cosp-(x3q-_x3c)*_sinp;
    }
    private double x3p(double x2q, double x3q) {
      return _x3c+(x2q-_x2c)*_sinp+(x3q-_x3c)*_cosp;
    }
    private double x2q(double x2p, double x3p) {
      return _x2c+(x2p-_x2c)*_cosp+(x3p-_x3c)*_sinp;
    }
    private double x3q(double x2p, double x3p) {
      return _x3c-(x2p-_x2c)*_sinp+(x3p-_x3c)*_cosp;
    }
    */
    private boolean inBounds(double x2p, double x3p) {
      return _s2p.getFirst()-HALF_LSINC<=x2p && 
             _s3p.getFirst()-HALF_LSINC<=x3p && 
              x2p<=_s2p.getLast()+HALF_LSINC &&
              x3p<=_s3p.getLast()+HALF_LSINC;
    }
  }

}
