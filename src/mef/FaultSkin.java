/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.io.*;
import java.util.*;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.io.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A linked list of fault cells that may be used to analyze faults.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */

public class FaultSkin implements Iterable<FaultCell>,Serializable {
  private static final long serialVersionUID = 1L;

  /**
   * Gets the cell that was the seed used to grow this skin.
   * @return the seed cell.
   */
  public FaultCell getSeed() {
    return _seed;
  }

  /**
   * Returns the number of cells in this skin.
   */
  public int size() {
    return _cellList.size();
  }

  /**
   * Gets an array of cells in this skin.
   * @return array of cells.
   */
  public FaultCell[] getCells() {
    return _cellList.toArray(new FaultCell[0]);
  }

  /**
   * Gets all cells in the specified skins.
   * @param skins array of skins for which to get cells.
   * @return array of cells.
   */
  public static FaultCell[] getCells(FaultSkin[] skins) {
    int ncell = countCells(skins);
    FaultCell[] cells = new FaultCell[ncell];
    int icell = 0;
    for (FaultSkin skin:skins)
      for (FaultCell cell:skin)
        cells[icell++] = cell;
    return cells;
  }

  public static FaultCell[] getCellsX(FaultSkin[] skins) {
    int ncell = countCells(skins);
    FaultCell[] cells = new FaultCell[ncell];
    int icell = 0;
    for (FaultSkin skin:skins)
      for (FaultCell cell:skin){
        cell.skin = null;
        cells[icell++] = cell;
      }
    return cells;
  }


  public static void setCells(FaultSkin[] skins, float[][][] fl) {
    for (FaultSkin skin:skins){
    for (FaultCell cell:skin) {
      int i1 = cell.i1;
      int i2 = cell.i2;
      int i3 = cell.i3;
      cell.fl = fl[i3][i2][i1];
      //float fli = fl[i3][i2][i1];
      //cell.fl = exp(fli)/exp(1f);
    }}
  }
  public static void getLikelihood(FaultSkin[] skins, float[][][] fl) {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        fl[i3][i2][i1] = cell.fl;
      }
    }
  }

  public static void setValueOnFaults(float v, 
    FaultSkin[] skins, float[][][] fl) 
  {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        fl[i3][i2][i1] = v;
      }
    }
  }

  public static void setValuesOnFaults(float v, 
    FaultSkin[] skins, float[][][] fl) 
  {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2m = cell.i2m;
        int i3m = cell.i3m;
        int i2p = cell.i2p;
        int i3p = cell.i3p;
        fl[i3m][i2m][i1] = v;
        fl[i3p][i2p][i1] = v;
      }
    }
  }

  public static void setValuesAwayFromFaults(
    FaultSkin[] skins, float[][][] f, float[][][] g) 
  {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2m = cell.i2m;
        int i3m = cell.i3m;
        int i2p = cell.i2p;
        int i3p = cell.i3p;
        if(i2m<0){i2m=0;}if(i2m>=n2){i2m=n2-1;}
        if(i2p<0){i2p=0;}if(i2p>=n2){i2p=n2-1;}
        if(i3m<0){i3m=0;}if(i3m>=n3){i3m=n3-1;}
        if(i3p<0){i3p=0;}if(i3p>=n3){i3p=n3-1;}
        g[i3m][i2m][i1] = f[i3m][i2m][i1];
        g[i3p][i2p][i1] = f[i3p][i2p][i1];
      }
    }
  }


  public static void getLikelihoods(FaultSkin[] skins, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] sc = new float[n3][n2][n1];
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1i = cell.i1;
        int i2m = cell.i2m;
        int i3m = cell.i3m;
        int i2p = cell.i2p;
        int i3p = cell.i3p;
        if(i2m<0){i2m=0;}if(i2m>=n2){i2m=n2-1;}
        if(i2p<0){i2p=0;}if(i2p>=n2){i2p=n2-1;}
        if(i3m<0){i3m=0;}if(i3m>=n3){i3m=n3-1;}
        if(i3p<0){i3p=0;}if(i3p>=n3){i3p=n3-1;}
        sc[i3m][i2m][i1i] += 1f;
        sc[i3p][i2p][i1i] += 1f;
        fl[i3m][i2m][i1i] += cell.fl;
        fl[i3p][i2p][i1i] += cell.fl;
      }
    }
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float sci = sc[i3][i2][i1];
      if(sci>1f){fl[i3][i2][i1] /= sci;}
    }}}
  }

  public static void getThrow(float mark, FaultSkin[] skins, float[][][] fs) {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        float s1 = cell.s1;
        float si = fs[i3][i2][i1];
        if(si==mark) {
          fs[i3][i2][i1] = s1;
        } else if(abs(s1)>abs(si)) {
          fs[i3][i2][i1] = s1;
        }
      }
    }
  }

  public static void getThrowThick(
    float mark, FaultSkin[] skins, float[][][] fs) {
    int n3 = fs.length;
    int n2 = fs[0].length;
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        int i2m = cell.i2m;
        int i2p = cell.i2p;
        int i3m = cell.i3m;
        int i3p = cell.i3p;
        if(i2m<0){i2m=0;}
        if(i3m<0){i3m=0;}
        if(i2p<0){i2p=0;}
        if(i3p<0){i3p=0;}
        if(i2p>=n2){i2p=n2-1;}
        if(i3p>=n3){i3p=n3-1;}
        if(i2m>=n2){i2m=n2-1;}
        if(i3m>=n3){i3m=n3-1;}
        float s1 = cell.s1;
        float si = fs[i3][i2][i1];
        if(si==mark) {
          fs[i3m][i2m][i1] = s1;
          fs[i3p][i2p][i1] = s1;
        } else if(abs(s1)>abs(si)) {
          fs[i3m][i2m][i1] = s1;
          fs[i3p][i2p][i1] = s1;
        }
      }
    }
  }


  /**
   * Returns the total number of cells in the specified skins.
   * @param skins array of skins for which to count cells.
   * @return the total number of cells.
   */
  public static int countCells(FaultSkin[] skins) {
    int ncell = 0;
    for (FaultSkin skin:skins)
      ncell += skin.size();
    return ncell;
  }

  /**
   * Returns an iterator for the cells in this skin. 
   * @return cell iterator.
   */
  public Iterator<FaultCell> iterator() {
    return _cellList.iterator();
  }

  /**
   * Returns array of arrays of cells linked above and below.
   * @return array of arrays of linked cells; by reference, not by copy.
   */
  public FaultCell[][] getCellsAB() {
    if (_cellsAB!=null)
      return _cellsAB;

    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();

    // For all cells in this skin, ...
    for (FaultCell cell:_cellList) {

      // If the cell is not already in an array, ...
      if (!cellSet.contains(cell)) {

        // Search above for the top cell.
        FaultCell c = cell;
        for (FaultCell ca=c.ca; ca!=null; ca=c.ca)
          c = ca;

        // Add the top cell and all cells below it to a new list.
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        for (; c!=null; c=c.cb) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the list of arrays.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    assert _cellList.size()==cellSet.size();

    // Convert the list of arrays to the array of arrays to be returned.
    _cellsAB = cellsList.toArray(new FaultCell[0][]);
    return _cellsAB;
  }

  /**
   * Returns array of arrays of cells linked left and right.
   * @return array of arrays of linked cells; by reference, not by copy.
   */
  public FaultCell[][] getCellsLR() {
    if (_cellsLR!=null)
      return _cellsLR;

    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();

    // For all cells in this skin, ...
    for (FaultCell cell:_cellList) {

      // If the cell is not already in an array, ...
      if (!cellSet.contains(cell)) {

        // Search left until we find no left nabor or until that left
        // nabor is the cell with which we began the search.
        FaultCell c = cell;
        for (FaultCell cl=c.cl; cl!=null && cl!=cell; cl=c.cl)
          c = cl;

        // Remember the leftmost cell found and add it to a new list.
        FaultCell cLeft = c;
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        cList.add(c);
        cellSet.add(c);

        // Add cells found to the right. Again beware of cycles.
        for (c=c.cr; c!=null && c!=cLeft; c=c.cr) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the list of arrays.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    assert _cellList.size()==cellSet.size();

    // Convert the list of arrays to the array of arrays to be returned.
    _cellsLR = cellsList.toArray(new FaultCell[0][]);
    checkCellArrays(_cellsLR);
    return _cellsLR;
  }

  /**
   * Smooths the normal vectors of cells in this skin.
   * @param nsmooth the number of smoothings.
   */
  public void smoothCellNormals(int nsmooth) {
    FaultCell.GetN getter = new FaultCell.GetN() {
      public float[] get(FaultCell cell) {
        return new float[]{cell.w1,cell.w2,cell.w3};
      }
    };
    FaultCell.SetN setter = new FaultCell.SetN() {
      public void set(FaultCell cell, float[] w) {
        float w1 = w[0]; 
        float w2 = w[1]; 
        float w3 = w[2];
        float ws = 1.0f/sqrt(w1*w1+w2*w2+w3*w3);
        w1 *= ws;
        w2 *= ws;
        w3 *= ws;
        cell.setNormalVector(w1,w2,w3);
      }
    };
    for (int ismooth=0; ismooth<nsmooth; ++ismooth)
      smoothN(getter,setter);
  }

  /**
   * Gets a cell nearest the centroid of this skin.
   * In illustrations, this cell is often a good representative.
   * @return the cell nearest the centroid.
   */
  public FaultCell getCellNearestCentroid() {
    float c1 = 0.0f;
    float c2 = 0.0f;
    float c3 = 0.0f;
    float cs = 0.0f;
    for (FaultCell c:_cellList) {
      c1 += c.fl*c.x1;
      c2 += c.fl*c.x2;
      c3 += c.fl*c.x3;
      cs += c.fl;
    }
    c1 /= cs;
    c2 /= cs;
    c3 /= cs;
    float dmin = Float.MAX_VALUE;
    FaultCell cmin = null;
    for (FaultCell c:_cellList) {
      float d = c.distanceSquaredTo(c1,c2,c3);
      if (d<dmin) {
        cmin = c;
        dmin = d;
      }
    }
    return cmin;
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, cells in this skin are represented by quads with specified
   * size, and colors corresponding to fault likelihoods.
   * @param size size (in samples) of the quads.
   * @param cmap colormap used to compute rgb colors from cell properties.
   * @param lhc true, if left-handed coordinate system; false, otherwise.
   */
  public float[][] getCellXyzUvwRgbForLikelihood(
      float size, ColorMap cmap, boolean lhc) {
    return FaultCell.getXyzUvwRgbForLikelihood(size,cmap,getCells(),lhc);
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, cells in this skin are represented by quads with specified
   * size, and colors corresponding to fault throws.
   * @param size size (in samples) of the quads.
   * @param cmap colormap used to compute rgb colors from cell properties.
   * @param lhc true, if left-handed coordinate system; false, otherwise.
   */
  public float[][] getCellXyzUvwRgbForThrow(
      float size, ColorMap cmap, boolean lhc) {
    return FaultCell.getXyzUvwRgbForThrow(size,cmap,getCells(),lhc);
  }

  /**
   * Gets arrays of packed cell coordinates for cell links.
   * Each returned array contains packed (x,y,z) cell coordinates for
   * exactly one above-below or left-right linked list of cells.
   * @return array of arrays of packed xyz cell coordinates.
   */
  public float[][] getCellLinksXyz() {
    FaultCell[][] cellsAB = getCellsAB();
    FaultCell[][] cellsLR = getCellsLR();
    int nsAB = cellsAB.length; // number of segments for AB links
    int nsLR = cellsLR.length; // number of segments for LR links
    int ns = nsAB+nsLR; // total number of segments
    float[][] xyz = new float[ns][];
    float[][] rgb = new float[ns][];
    for (int is=0; is<ns; ++is) { // for all segments, ...
      FaultCell[] cells = (is<nsAB)?cellsAB[is]:cellsLR[is-nsAB]; // the cells
      int ncell = cells.length; // number of cells in this segment
      int np = ncell; // number of points in this segment
      if (is>=nsAB && cells[0].cl==cells[ncell-1]) // if a LR cycle, ...
        ++np; // then add one more so we end with the starting point
      float[] xyzi = new float[3*np]; // xyz for this segment
      for (int ip=0,ic=0; ip<np; ++ip) {
        FaultCell cell = cells[ip%ncell]; // % to handle any LR cycle
        xyzi[ic++] = cell.x3;
        xyzi[ic++] = cell.x2;
        xyzi[ic++] = cell.x1;
      }
      xyz[is] = xyzi;
    }
    return xyz;
  }

  public float[][] getCellLinksRgb(float r, float g, float b, float[][] xyz) {
    int ns = xyz.length;
    float[][] rgb = new float[ns][];
    for (int is=0; is<ns; ++is) {
      int np = xyz[is].length/3;
      rgb[is] = new float[np*3];
      for (int ip=0; ip<np; ++ip) {
        int ik = 3*ip;
        rgb[is][ik+0] = r;
        rgb[is][ik+1] = g;
        rgb[is][ik+2] = b;
      }
    }
    return rgb;
  }

  public float[] getCellRgb(float r, float g, float b, float[] xyz) {
    int np = xyz.length/3;
    float[] rgb = new float[np*3];
    for (int ip=0; ip<np; ++ip) {
      int ik = 3*ip;
      rgb[ik+0] = r;
      rgb[ik+1] = g;
      rgb[ik+2] = b;
    }
    return rgb;
  }


  /**
   * Returns a fault skin read from a file with specified name.
   * @param fileName the fault skin file name.
   * @return the fault skin.
   */
  public static FaultSkin readFromFile(String fileName) {
    FaultSkin skin = new FaultSkin();
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      int ncell = ais.readInt();
      int i1seed = ais.readInt();
      int i2seed = ais.readInt();
      int i3seed = ais.readInt();
      ArrayList<FaultCell> cellList = new ArrayList<FaultCell>(ncell);
      for (int icell=0; icell<ncell; ++icell) {
        float x1 = ais.readFloat();
        float x2 = ais.readFloat();
        float x3 = ais.readFloat();
        float fl = ais.readFloat();
        float fp = ais.readFloat();
        float ft = ais.readFloat();
        FaultCell cell = new FaultCell(x1,x2,x3,fl,fp,ft);
        cell.skin = skin;
        cellList.add(cell);
        cell.s1 = ais.readFloat();
        cell.s2 = ais.readFloat();
        cell.s3 = ais.readFloat();
        cell.r1 = ais.readFloat();
        cell.r2 = ais.readFloat();
        cell.r3 = ais.readFloat();
      }
      FaultCell[] cells = cellList.toArray(new FaultCell[0]);
      FaultCellGrid fcg = new FaultCellGrid(cells);
      for (FaultCell cell:cellList) {
        int i1a = ais.readInt();
        int i2a = ais.readInt();
        int i3a = ais.readInt();
        int i1b = ais.readInt();
        int i2b = ais.readInt();
        int i3b = ais.readInt();
        int i1l = ais.readInt();
        int i2l = ais.readInt();
        int i3l = ais.readInt();
        int i1r = ais.readInt();
        int i2r = ais.readInt();
        int i3r = ais.readInt();
        if (i1a!=INULL) 
          cell.ca = fcg.get(i1a,i2a,i3a);
        if (i1b!=INULL) 
          cell.cb = fcg.get(i1b,i2b,i3b);
        if (i1l!=INULL) 
          cell.cl = fcg.get(i1l,i2l,i3l);
        if (i1r!=INULL) 
          cell.cr = fcg.get(i1r,i2r,i3r);
      }
      ais.close();
      skin._seed = fcg.get(i1seed,i2seed,i3seed);
      skin._cellList = cellList;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return skin;
  }

  public void smooth(int n) { 
    for(int i=0; i<n; ++i) {
      for (FaultCell cell:_cellList) {
        FaultCell ca = cell.ca;
        FaultCell cb = cell.cb;
        FaultCell cl = cell.cl;
        FaultCell cr = cell.cr;
        float cs = 1f;
        float x1 = cell.x1;
        float x2 = cell.x2;
        float x3 = cell.x3;
        float fl = cell.fl;
        float s1 = cell.s1;
        if(ca!=null) {
          cs += 1f;
          x1 += ca.x1;
          x2 += ca.x2;
          x3 += ca.x3;
          fl += ca.fl;
          s1 += ca.s1;
        }
        if(cb!=null) {
          cs += 1f;
          x1 += cb.x1;
          x2 += cb.x2;
          x3 += cb.x3;
          fl += cb.fl;
          s1 += cb.s1;
        }
        if(cl!=null) {
          cs += 1f;
          x1 += cl.x1;
          x2 += cl.x2;
          x3 += cl.x3;
          fl += cl.fl;
          s1 += cl.s1;
        }
        if(cr!=null) {
          cs += 1f;
          x1 += cr.x1;
          x2 += cr.x2;
          x3 += cr.x3;
          fl += cr.fl;
          s1 += cr.s1;
        }
        cs = 1f/cs;
        cell.x1 = x1*cs;
        cell.x2 = x2*cs;
        cell.x3 = x3*cs;
        cell.fl = fl*cs;
        cell.s1 = s1*cs;
      }
    }

  }


  /**
   * Writes a fault skin to a file with specified name.
   * @param fileName the fault skin file name.
   * @param skin the fault skin.
   */
  public static void writeToFile(String fileName, FaultSkin skin) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      aos.writeInt(skin.size());
      aos.writeInt(skin._seed.i1);
      aos.writeInt(skin._seed.i2);
      aos.writeInt(skin._seed.i3);
      for (FaultCell cell:skin) {
        aos.writeFloat(cell.x1); 
        aos.writeFloat(cell.x2); 
        aos.writeFloat(cell.x3);
        aos.writeFloat(cell.fl); 
        aos.writeFloat(cell.fp); 
        aos.writeFloat(cell.ft);
        aos.writeFloat(cell.s1); 
        aos.writeFloat(cell.s2); 
        aos.writeFloat(cell.s3);
        aos.writeFloat(cell.r1); 
        aos.writeFloat(cell.r2); 
        aos.writeFloat(cell.r3);
      }
      for (FaultCell cell:skin) {
        FaultCell[] nabors = new FaultCell[]{cell.ca,cell.cb,cell.cl,cell.cr};
        for (FaultCell nabor:nabors) {
          if (nabor!=null) {
            aos.writeInt(nabor.i1);
            aos.writeInt(nabor.i2);
            aos.writeInt(nabor.i3);
          } else {
            aos.writeInt(INULL);
            aos.writeInt(INULL);
            aos.writeInt(INULL);
          }
        }
      }
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public static FaultSkin readFromFileSlow(String fileName) {
    try {
      FileInputStream fis = new FileInputStream(fileName);
      ObjectInputStream ois = new ObjectInputStream(fis);
      FaultSkin skin = (FaultSkin)ois.readObject();
      ois.close();
      return skin;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public static void writeToFileSlow(String fileName, FaultSkin skin) {
    try {
      FileOutputStream fos = new FileOutputStream(fileName);
      ObjectOutputStream oos = new ObjectOutputStream(fos);
      oos.writeObject(skin);
      oos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }


  /////////////////////////////////////////////////////////////////////////
  // package

  /**
   * Constructs an empty skin.
   */
  FaultSkin() {
    _cellList = new ArrayList<FaultCell>();
  }

  /**
   * Adds the specified skinless cell to this skin.
   * @param cell the cell to be added.
   */
  void add(FaultCell cell) {
    assert cell.skin==null;
    cell.skin = this;
    if (_seed==null)
      _seed = cell;
    _cellList.add(cell);
    _cellsAB = null;
    _cellsLR = null;
  }

  /**
   * Smooths one value stored in the cells of this skin. The value smoothed is
   * that accessed by the specified getter and setter. Each smoothed value is
   * an average of the values in a cell and its cell nabors. 
   */
  void smooth1(FaultCell.Get1 getter, FaultCell.Set1 setter) {
    int ncell = size();
    float[] vals = new float[ncell];
    float[] cnts = new float[ncell];
    FaultCell[] cellNabors = new FaultCell[4];
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float valCell = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      cellNabors[2] = cell.cl;
      cellNabors[3] = cell.cr;
      for (FaultCell cellNabor:cellNabors) {
        if (cellNabor!=null) {
          float valNabor = getter.get(cellNabor);
          vals[icell] += valCell+valNabor;
          cnts[icell] += 2.0f;
        }
      }
    }
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float cnti = cnts[icell];
      float vali = vals[icell]/(cnti>0.0f?cnti:1.0f);
      setter.set(cell,vali);
    }
  }

  /**
   * Smooths multiple values stored in the cells of this skin. The values
   * smoothed are those accessed by the specified getter and setter. Each
   * smoothed value is an average of the values in a cell and its cell nabors. 
   */
  void smoothN(FaultCell.GetN getter, FaultCell.SetN setter) {
    int ncell = size();
    int nval = getter.get(_seed).length;
    float[][] vals = new float[ncell][nval];
    float[] cnts = new float[ncell];
    FaultCell[] cellNabors = new FaultCell[4];
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float[] valsCell = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      cellNabors[2] = cell.cl;
      cellNabors[3] = cell.cr;
      for (FaultCell cellNabor:cellNabors) {
        if (cellNabor!=null) {
          float[] valsNabor = getter.get(cellNabor);
          for (int ival=0; ival<nval; ++ival)
            vals[icell][ival] += valsCell[ival]+valsNabor[ival];
          cnts[icell] += 2.0f;
        }
      }
    }
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float cnti = cnts[icell];
      float scli = 1.0f/(cnti>0.0f?cnti:1.0f);
      for (int ival=0; ival<nval; ++ival)
        vals[icell][ival] *= scli;
      setter.set(cell,vals[icell]);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // private

  private static final int INULL = -Integer.MAX_VALUE; // null index

  private FaultCell _seed; // cell in this skin with highest fl; null, if empty
  private ArrayList<FaultCell> _cellList; // list of cells in this skin
  private FaultCell[][] _cellsAB; // arrays of cells from above to below
  private FaultCell[][] _cellsLR; // arrays of cells from left to right

  private void checkCellArrays() {
    if (_cellsAB!=null)
      checkCellArrays(_cellsAB);
    if (_cellsLR!=null)
      checkCellArrays(_cellsLR);
  }

  private static void checkCellArrays(FaultCell[][] cells) {
    HashSet<FaultCell> cellSet = new HashSet<FaultCell>();
    int ncell = cells.length;
    for (int icell=0; icell<ncell; ++icell) {
      int mcell = cells[icell].length;
      for (int jcell=0; jcell<mcell; ++jcell) {
        FaultCell c = cells[icell][jcell];
        assert cellSet.add(c);
      }
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
