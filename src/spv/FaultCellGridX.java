/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package spv;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Fault cells in a 3D sampling grid. Each grid sample indexed by (i1,i2,i3)
 * contains either one fault cell or null. The grid facilitates searches for
 * cell nabors in skins and fast iterations along fault traces tangent to
 * fault strike and fault curves tangent to fault dip.
 * <p> 
 * Grid indices need not (and typically do not) begin at zero. Index bounds
 * for a fault cell grid are determined by the minima and maxima of indices of
 * cells used to construct the grid.
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.19
 */
public class FaultCellGridX {

  /**
   * Constructs a fault grid for specified cells. Grid index bounds are
   * determined by the minimum and maximum indices of the specified cells.
   * @param cells array of cells to be included in the grid.
   */
  public FaultCellGridX(int n1, int n2, int n3) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _cells = new FaultCell[n3][n2][n1];
  }


  /**
   * Gets the fault cell with specified indices, if any.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   * @return the fault cell; null, if none or if indices are out of bounds.
   */
  public FaultCell get(int i1, int i2, int i3) {
    return _cells[i3][i2][i1];
  }

  /**
   * Sets the specified fault cell. Uses the cell's {x1,x2,x3} coordinates to
   * determine the indices of the cell in this grid.
   * @param cell the fault cell.
   */
  public void set(FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    _cells[i3][i2][i1] = cell;
  }

  public void set(FaultSkin skin) {
    for (FaultCell cell:skin)
      set(cell);
  }

  public void set(FaultCell[] cells) {
    for (FaultCell cell:cells)
      set(cell);
  }
  public void setCellsInBox(
    FaultCell cell, int d1, int d2, int d3) {
    int c1 = cell.i1;
    int c2 = cell.i2;
    int c3 = cell.i3;
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      if(ci!=null) ci.used = true;
    }}}
  }

  public FaultCell findCellInBox(
    int c1, int c2, int c3, int d1, int d2, int d3) {
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    float dx = 1000f;
    FaultCell cell = null;
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      float di = (c1-i1)*(c1-i1)+(c2-i2)*(c2-i2)+(c3-i3)*(c3-i3);
      if(ci.skin==null&& di<dx) {
        dx = di;
        cell = ci;
      }
    }}}
    return cell;
  }

  public boolean findCellsInBox(
    float x1, float x2, float x3, float fp, int d1, int d2, int d3) {
    int c1 = round(x1);
    int c2 = round(x2);
    int c3 = round(x3);
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    if(fp>180) fp = 360-fp;
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      if(ci!=null) {
        float fpi = ci.fp;
        if(fpi>180) fpi = 360-fpi;
        float dp = abs(fp-fpi);
        if(dp<40f) return true;
      }
    }}}
    return false;
  }


  public boolean findCellsInBox(
    FaultCell cell, int d1, int d2, int d3) {
    int c1 = cell.i1;
    int c2 = cell.i2;
    int c3 = cell.i3;
    float fp = cell.fp;
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      if(ci!=null) {
        float fpi = ci.fp;
        float dp = abs(fp-fpi);
        dp = min(dp,360-dp);
        if(dp<40f) return true;
      }
    }}}
    return false;
  }





  ///////////////////////////////////////////////////////////////////////////
  // private

  private FaultCell[][][] _cells; // array of cells
  private int _n1,_n2,_n3;


}
