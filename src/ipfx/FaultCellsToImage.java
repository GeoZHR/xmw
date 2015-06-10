/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipfx;


import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A fault cell is an oriented point located on a fault. Each fault cell 
 * contains three attributes of fault likelihood, strike and dip.
 * Each fault cell corresponds to one and only one seismic image sample.
 *
 * Methods are provided here to conveniently convert fault cells into 
 * fault images used as weights for other image processing like flattening.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.07
 */

public class FaultCellsToImage {

  /**
   * Constructs with specified fault cells.
   * @param cells input fault cells.
   */
  public FaultCellsToImage(FaultCell[] cells) {
    _cells = cells;
  }

  public float[][][] getLikelihoodThin(float flmin, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] fli = new float[n3][n2][n1];
    for (FaultCell cell:_cells) {
      int i1 = cell.i1;
      int i2 = cell.i2;
      int i3 = cell.i3;
      if(cell.fl<flmin) {continue;}
      if(i1<0) {i1=0;} if(i1>=n1){i1=n1-1;}
      if(i2<0) {i2=0;} if(i2>=n2){i2=n2-1;}
      if(i3<0) {i3=0;} if(i3>=n3){i3=n3-1;}
      fli[i3][i2][i1] = fl[i3][i2][i1];
    }
    return fli;
  }

  public float[][][] getLikelihoodThick(float flmin, float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] fli = new float[n3][n2][n1];
    for (FaultCell cell:_cells) {
      int i1i = cell.i1;
      int i2i = cell.i2;
      int i3i = cell.i3;
      int i2m = cell.i2m;
      int i3m = cell.i3m;
      int i2p = cell.i2p;
      int i3p = cell.i3p;
      if(cell.fl<flmin) {continue;}
      if(i1i<0) {i1i=0;} if(i1i>=n1){i1i=n1-1;}
      if(i2i<0) {i2i=0;} if(i2i>=n2){i2i=n2-1;}
      if(i3i<0) {i3i=0;} if(i3i>=n3){i3i=n3-1;}
      if(i2m<0) {i2m=0;} if(i2p>=n2){i2p=n2-1;}
      if(i3m<0) {i3m=0;} if(i3p>=n3){i3p=n3-1;}
      fli[i3i][i2i][i1i] = fl[i3i][i2i][i1i];
      fli[i3m][i2m][i1i] = fl[i3m][i2m][i1i];
      fli[i3p][i2p][i1i] = fl[i3p][i2p][i1i];
    }
    return fli;
  }


  private FaultCell[] _cells = null;
}

