/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import java.util.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;
import static ifs.FaultGeometry.*;

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
    return skinsXX(cells,fl);
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
            if(checkPlanar(cell)) {
              boolean missNabors = true;
              missNabors = findInOldCells(cell);
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
    if(dp/nc<1.0f){return true;}
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
    _cellsX.set(ca);
    _cellsX.set(cb);
    _cellsX.set(cl);
    _cellsX.set(cr);
    _cellsX.set(cell);
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
