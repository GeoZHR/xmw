/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.io.Serializable;

import edu.mines.jtk.io.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

import static mef.FaultGeometry.*;

/**
 * A fault cell is an oriented point located on a fault. Fault cells can
 * be linked to form fault skins, which may be used to analyze faults.
 * <p>
 * Fault cells are computed from images of fault likelihoods, strikes and
 * dips. Each fault cell is an oriented point located on a ridge in an image
 * of fault likelihood. Fault cells have indices (i1,i2,i3) that indicate
 * which image sample is nearest to this ridge. In this way an image sample
 * is associated with either no cell or one cell.
 * <p>
 * A fault cell has up to four neighbors ("nabors") that lie above, below,
 * left and right of the cell when viewed from above the fault, that is, when
 * looking from the hanging wall toward the footwall. Links to nabors enables
 * cells to form a skin of connected cells, which represents a fault.
 * <p>
 * Links to left and right cell nabors can be used to iterate over all cells
 * along a fault trace, a path of constant depth that is everywhere tangent to
 * fault strike. Likewise, links to cell nabors above and below a cell can be
 * used to iterate up or down a fault. However, this simple up or down
 * iteration typically does not coincide with a fault curve that is everywhere
 * tangent to fault dip.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */

public class FaultCell implements Serializable {
  private static final long serialVersionUID = 1L;

  public void setInterp(boolean needInterp) {
    this.needInterp = needInterp;
  }

  public void setNormal(float w1, float w2, float w3) {
    setNormalVector(w1,w2,w3);
  }

  public FaultCell getCa() {
    return this.ca;
  }

  public FaultCell getCb() {
    return this.cb;
  }

  public void setFl(float fl) {
    this.fl = fl;
  }

  /**
   * Gets the fault likelihood for this cell.
   * @return the fault likelihood.
   */
  public float getFl() {
    return fl;
  }

  public float getFp() {
    return fp;
  }

  public float getFt() {
    return ft;
  }

  public float getSmp() {
    return smp;
  }

  public void shift(float d1, float d2, float d3) {
    this.x1 += d1;
    this.x2 += d2;
    this.x3 += d3;
    this.i1 = round(this.x1);
    this.i2 = round(this.x2);
    this.i3 = round(this.x3);
  }


  public int[] getI() {
    return new int[]{i1,i2,i3};
  }

  public int[] getIm() {
    return new int[]{i1,i2m,i3m};
  }

  public int[] getIp() {
    return new int[]{i1,i2p,i3p};
  }

  public float[] getU() {
    return new float[]{u1,u2,u3};
  }
  public float[] getV() {
    return new float[]{v1,v2,v3};
  }

  public void setNum(float num) {
    this.num = num;
  }

  public void setDen(float den) {
    this.den = den;
  }


  public int getM1() {
    return i1;
  }
  public int getM2() {
    return i2m;
  }
  public int getM3() {
    return i3m;
  }

  public int getP1() {
    return i1;
  }
  public int getP2() {
    return i2p;
  }
  public int getP3() {
    return i3p;
  }

  public int getI1() {
    return i1;
  }
  public int getI2() {
    return i2;
  }
  public int getI3() {
    return i3;
  }
  public float getT1() {
    return t1;
  }
  public float getT2() {
    return t2;
  }
  public float getT3() {
    return t3;
  }
  public void setT1(float t1) {
    this.t1 = t1;
  }
  public void setT2(float t2) {
    this.t2 = t2;
  }
  public void setT3(float t3) {
    this.t3 = t3;
  }

  public void setUnfaultShifts(float[] ts) {
    this.t1=ts[0];
    this.t2=ts[1];
    this.t3=ts[2];
  }

  public static void getFlThick(
    float mark, FaultCell[] cells, float[][][] fs) {
    int n3 = fs.length;
    int n2 = fs[0].length;
    for (FaultCell cell:cells) {
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
      float fl = cell.fl;
      float fi = fs[i3][i2][i1];
      if(fi==mark) {
        fs[i3m][i2m][i1] = fl;
        fs[i3p][i2p][i1] = fl;
      } else if(abs(fl)>abs(fi)) {
        fs[i3m][i2m][i1] = fl;
        fs[i3p][i2p][i1] = fl;
      }
    }
  }

  public static void getFpThick(
    float mark, FaultCell[] cells, float[][][] fs) {
    int n3 = fs.length;
    int n2 = fs[0].length;
    for (FaultCell cell:cells) {
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
      float fp = cell.fp;
      float fi = fs[i3][i2][i1];
      if(fi==mark) {
        fs[i3m][i2m][i1] = fp;
        fs[i3p][i2p][i1] = fp;
      } else if(abs(fp)>abs(fi)) {
        fs[i3m][i2m][i1] = fp;
        fs[i3p][i2p][i1] = fp;
      }
    }
  }

  public static void getFtThick(
    float mark, FaultCell[] cells, float[][][] fs) {
    int n3 = fs.length;
    int n2 = fs[0].length;
    for (FaultCell cell:cells) {
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
      float ft = cell.ft;
      float fi = fs[i3][i2][i1];
      if(fi==mark) {
        fs[i3m][i2m][i1] = ft;
        fs[i3p][i2p][i1] = ft;
      } else if(abs(ft)>abs(fi)) {
        fs[i3m][i2m][i1] = ft;
        fs[i3p][i2p][i1] = ft;
      }
    }
  }


    /**
   * Gets the slip vectors (s1,s2,s3) of this cell.
   * @return array {s1,s2,s3} of coordinates.
   */
  public float[] getS() {
    return new float[]{s1,s2,s3};
  }
  /**
   * Gets the 1st component of the slip vector for this cell.
   * @return the 1st coordinate.
   */
  public float getS1() {
    return s1;
  }

  /**
   * Gets the 2nd component of the slip vector for this cell.
   * @return the 2nd component.
   */
  public float getS2() {
    return s2;
  }

  /**
   * Gets the 3rd component of the slip vector for this cell.
   * @return the 3rd component.
   */
  public float getS3() {
    return s3;
  }

  public float getR1() {
    return r1;
  }
  public float getR2() {
    return r2;
  }
  public float getR3() {
    return r3;
  }


  /**
   * Gets the coordinates (x1,x2,x3) of the location for this cell.
   * @return array {x1,x2,x3} of coordinates.
   */
  public float[] getX() {
    return new float[]{x1,x2,x3};
  }

  /**
   * Gets the 1st coordinate of the location for this cell.
   * @return the 1st coordinate.
   */
  public float getX1() {
    return x1;
  }

  /**
   * Gets the 2nd coordinate of the location for this cell.
   * @return the 2nd coordinate.
   */
  public float getX2() {
    return x2;
  }

  /**
   * Gets the 3rd coordinate of the location for this cell.
   * @return the 3rd coordinate.
   */
  public float getX3() {
    return x3;
  }

  /**
   * Gets the components (w1,w2,w3) of the normal vector for this cell.
   * @return array {w1,w2,w3} of coordinates.
   */
  public float[] getW() {
    return new float[]{w1,w2,w3};
  }

  /**
   * Gets the 1st component of the normal vector for this cell.
   * @return the 1st component.
   */
  public float getW1() {
    return w1;
  }

  /**
   * Gets the 2nd component of the normal vector for this cell.
   * @return the 2nd component.
   */
  public float getW2() {
    return w2;
  }

  /**
   * Gets the 3rd component of the normal vector for this cell.
   * @return the 3rd component.
   */
  public float getW3() {
    return w3;
  }

  /**
   * Gets the 2nd component of the normal vector for this cell.
   * @return the 2nd component.
   */
  public float getV2() {
    return v2;
  }

  /**
   * Gets the 3rd component of the normal vector for this cell.
   * @return the 3rd component.
   */
  public float getV3() {
    return v3;
  }

  /**
   * Returns an array of packed (x,y,z) coordinates for a fault curve.
   * The fault curve is everywhere tangent to fault dip, and contains the
   * point for this cell. Returned coordinates are in above-to-below order.
   * @return array of packed (x,y,z) coordinates.
   */
  public float[] getFaultCurveXyz() {
    FloatList xyz = new FloatList();
    float[] p = new float[3];

    // Gather xyz for this cell and above this cell by walking up dip.
    FaultCell cell = this;
    p[0] = x1; p[1] = x2; p[2] = x3;
    for (int j1=cell.i1; j1==cell.i1; --j1) {
      xyz.add(p[2]); xyz.add(p[1]); xyz.add(p[0]);
      cell = cell.walkUpDipFrom(p);
    }

    // Remember the number of xyz gathered, including xyz for this cell.
    int na = xyz.n/3;

    // Gather xyz for cells below this one by walking down dip. We have
    // already gathered the xyz for this cell, so now we skip to the one
    // below it.
    cell = this;
    p[0] = x1; p[1] = x2; p[2] = x3;
    cell = cell.walkDownDipFrom(p); // skip this cell
    for (int j1=this.i1+1; j1==cell.i1; ++j1) {
      xyz.add(p[2]); xyz.add(p[1]); xyz.add(p[0]);
      cell = cell.walkDownDipFrom(p);
    }

    // Flip the order of all xyz gathered above this cell while walking up.
    float[] xyzs = xyz.trim();
    for (int ia=1,i=ia*3,ja=na-1,j=ja*3; ia<ja; ++ia,i+=3,--ja,j-=3) {
      float xi = xyzs[i  ];
      float yi = xyzs[i+1];
      float zi = xyzs[i+2];
      xyzs[i  ] = xyzs[j  ];
      xyzs[i+1] = xyzs[j+1];
      xyzs[i+2] = xyzs[j+2];
      xyzs[j  ] = xi;
      xyzs[j+1] = yi;
      xyzs[j+2] = zi;
    }
    return xyzs;
  }

  /**
   * Returns an array of packed (x,y,z) coordinates for a fault trace.
   * The fault trace is everywhere tangent to fault strike, and contains 
   * the point for this cell. Returned coordinates are left-to-right order.
   * @return array of packed (x,y,z) coordinates.
   */
  public float[] getFaultTraceXyz() {
    FloatList xyz = new FloatList();

    // First gather coordinates for this cell.
    xyz.add(x3); xyz.add(x2); xyz.add(x1);

    // Then gather coordinates for cells to the left of this cell. Take care
    // to handle the case in which this cell is found while walking left.
    FaultCell c;
    for (c=this.cl; c!=null && c!=this; c=c.cl) {
      xyz.add(c.x3); xyz.add(c.x2); xyz.add(c.x1);
    }

    // Remember number of xyz gathered, including xyz for this cell.
    int nl = xyz.n/3;

    // If we did not end at this cell, then gather coordinates of all
    // cells to the right of this cell.
    if (c!=this) {
      for (c=this.cr; c!=null; c=c.cr) {
        xyz.add(c.x3); xyz.add(c.x2); xyz.add(c.x1);
      }
    }

    // Flip the order of all xyz gathered left of this cell.
    float[] xyzs = xyz.trim();
    for (int il=1,i=il*3,jl=nl-1,j=jl*3; il<jl; ++il,i+=3,--jl,j-=3) {
      float xi = xyzs[i  ];
      float yi = xyzs[i+1];
      float zi = xyzs[i+2];
      xyzs[i  ] = xyzs[j  ];
      xyzs[i+1] = xyzs[j+1];
      xyzs[i+2] = xyzs[j+2];
      xyzs[j  ] = xi;
      xyzs[j+1] = yi;
      xyzs[j+2] = zi;
    }
    return xyzs;
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, the specified cells are represented by quads with specified
   * size, and colors corresponding to fault throws.
   * @param size the size (in samples) of the quads.
   * @param cmap the colormap used to compute rgb colors from floats.
   * @param cells the cells for which to return the arrays.
   * @param lhc true, if left-handed coordinates; false, otherwise.
   * @return arrays {xyz,uvw,rgb}.
   */
  public static float[][] getXyzUvwRgbForThrow(
      float size, ColorMap cmap, FaultCell[] cells, boolean lhc) {
    return getXyzUvwRgb(size,cmap,cells,new Get1() {
      public float get(FaultCell cell) { return cell.s1; }
    },lhc);
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, the specified cells are represented by quads with specified
   * size, and colors corresponding to fault likelihoods.
   * @param size the size (in samples) of the quads.
   * @param cmap the colormap used to compute rgb colors from floats.
   * @param cells the cells for which to return the arrays.
   * @param lhc true, if left-handed coordinates; false, otherwise.
   * @return arrays {xyz,uvw,rgb}.
   */
  public static float[][] getXyzUvwRgbForLikelihood(
      float size, ColorMap cmap, FaultCell[] cells, boolean lhc) {
    return getXyzUvwRgb(size,cmap,cells,new Get1() {
      public float get(FaultCell cell) { return cell.fl; }
    },lhc);
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  int i1,i2,i3; // cell indices
  float x1,x2,x3; // cell coordinates
  float fl,fp,ft; // likelihood, strike (phi) and dip (theta)
  float u1,u2,u3,us; // dip vector and scale factor = 1/sin(theta)
  float v1,v2,v3; // strike vector
  float w1,w2,w3; // normal vector
  float w11,w12,w13,w22,w23,w33;
  float u11,u12,u13,u22,u23,u33;
  float v11,v12,v13,v22,v23,v33;
  public FaultCell ca,cb,cl,cr; // nabors above, below, left and right
  FaultSkin skin; // if not null, the skin to which this cell belongs
  int i2m,i2p; // sample indices i2 for minus and plus sides of cell
  int i3m,i3p; // sample indices i3 for minus and plus sides of cell
  float[] emp; // array of minus-plus alignment errors
  float smp; // shift from minus side to plus side of cell
  float s1,s2,s3; // fault dip-slip vector
  float r1,r2,r3; // fault dip-slip vector
  float t1,t2,t3; // fault dip-slip vector
  public float num = 0f;
  public float den = 0f;
  public float flr = 0f;
  boolean interp;
  boolean notUsed;
  boolean intersect;
  boolean needInterp;

  public interface Get1 { public float get(FaultCell cell); }
  public interface GetN { public float[] get(FaultCell cell); }
  public interface Set1 { public void set(FaultCell cell, float value); }
  public interface SetN { public void set(FaultCell cell, float[] values); }

  public FaultCell(float x1, float x2, float x3, float fl, float fp, float ft) {
    set(x1,x2,x3,fl,fp,ft);
  }

  /**
   * Sets the normal vector for this cell. Also modifies all other cell
   * properties related to the normal vector, such as fault strike and dip.
   */
  void setNormalVector(float w1, float w2, float w3) {
    float fp = faultStrikeFromNormalVector(w1,w2,w3);
    float ft = faultDipFromNormalVector(w1,w2,w3);
    set(x1,x2,x3,fl,fp,ft);
  }

  /**
   * Returns the distance squared from this cell to a point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return distance squared.
   */
  float distanceTo(float p1, float p2, float p3) {
    return sqrt(distanceSquaredTo(p1,p2,p3));
  }

  /**
   * Returns the distance squared from this cell to a point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return distance squared.
   */
  float distanceSquaredTo(float p1, float p2, float p3) {
    float d1 = p1-x1;
    float d2 = p2-x2;
    float d3 = p3-x3;
    return d1*d1+d2*d2+d3*d3;
  }

  /**
   * Returns the signed distance from the plane of this cell to a point. The
   * distance is positive if the point lies on the side of the plane toward
   * which the normal vector points, negative if on the other side.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return signed distance.
   */
  float distanceFromPlaneTo(float p1, float p2, float p3) {
    return w1*(p1-x1)+w2*(p2-x2)+w3*(p3-x3);
  }

  /**
   * Gets the cell above this cell that is nearest to the specified point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return the cell above; null, if none.
   */
  FaultCell getCellAboveNearestTo(float p1, float p2, float p3) {
    FaultCell cla = (cl!=null)?cl.ca:null;
    FaultCell cra = (cr!=null)?cr.ca:null;
    return nearestCell(ca,cla,cra,p1,p2,p3);
  }

  /**
   * Gets the cell below this cell that is nearest to the specified point.
   * @param p1 1st coordinate of point.
   * @param p2 2nd coordinate of point.
   * @param p3 3rd coordinate of point.
   * @return the cell below; null, if none.
   */
  FaultCell getCellBelowNearestTo(float p1, float p2, float p3) {
    FaultCell clb = (cl!=null)?cl.cb:null;
    FaultCell crb = (cr!=null)?cr.cb:null;
    return nearestCell(cb,clb,crb,p1,p2,p3);
  }

  /**
   * Walks a point up the fault along a curve tangent to fault dip.
   * The input point is assumed to lie in the plane of this cell,
   * and the output point will lie in the plane of the returned cell.
   * That returned cell will be one immediately above this cell, if
   * sufficient nabors exist, or this cell, otherwise.
   * @param p input and output array {p1,p2,p3} of point coordinates.
   * @return the cell with a plane that contains the output point.
   */
  FaultCell walkUpDipFrom(float[] p) {
    FaultCell cell = this;
    float p1 = p[0];
    float p2 = p[1];
    float p3 = p[2];
    //assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;

    // Use dip vector of specified cell to walk the point up its plane.
    p1 -= 1.0f;
    p2 -= cell.us*cell.u2;
    p3 -= cell.us*cell.u3;

    // If a cell below is found, project point horizontally onto its plane.
    FaultCell ca = cell.getCellAboveNearestTo(p1,p2,p3);
    if (ca!=null) {
      cell = ca;
      float us = cell.us;
      float ws = us*us*cell.distanceFromPlaneTo(p1,p2,p3);
      p2 -= ws*cell.w2;
      p3 -= ws*cell.w3;
    }

    // Return updated point and cell.
    //assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;
    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    return cell;
  }

  /**
   * Walks a point down the fault along a curve tangent to fault dip.
   * The input point is assumed to lie in the plane of this cell,
   * and the output point will lie in the plane of the returned cell.
   * That returned cell will be one immediately below this cell, if
   * sufficient cell nabors exist, or this cell, otherwise.
   * @param p input and output array {p1,p2,p3} of point coordinates.
   * @return the cell with a plane that contains the output point.
   */
  FaultCell walkDownDipFrom(float[] p) {
    FaultCell cell = this;
    float p1 = p[0];
    float p2 = p[1];
    float p3 = p[2];
    //assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;

    // Use dip vector of specified cell to walk the point down its plane.
    p1 += 1.0f;
    p2 += cell.us*cell.u2;
    p3 += cell.us*cell.u3;

    // If a cell below is found, project point horizontally onto its plane.
    FaultCell cb = cell.getCellBelowNearestTo(p1,p2,p3);
    if (cb!=null) {
      cell = cb;
      float us = cell.us;
      float ws = us*us*cell.distanceFromPlaneTo(p1,p2,p3);
      p2 -= ws*cell.w2;
      p3 -= ws*cell.w3;
    }

    // Return updated point and cell.
    //assert abs(cell.distanceFromPlaneTo(p1,p2,p3))<0.01f;
    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    return cell;
  }

  /////////////////////////////////////////////////////////////////////////
  // private

  public void set(
      float x1, float x2, float x3, 
      float fl, float fp, float ft) {
    this.x1 = x1; 
    this.x2 = x2; 
    this.x3 = x3;
    this.fl = fl; 
    this.fp = fp; 
    this.ft = ft;
    i1 = round(x1);
    i2 = round(x2);
    i3 = round(x3);
    intersect = false;
    notUsed = true;
    needInterp = true;
    interp = false;
    float[] u = faultDipVectorFromStrikeAndDip(fp,ft);
    float[] v = faultStrikeVectorFromStrikeAndDip(fp,ft);
    float[] w = faultNormalVectorFromStrikeAndDip(fp,ft);
    u1 = u[0]; u2 = u[1]; u3 = u[2]; us = 1.0f/u1;
    v1 = v[0]; v2 = v[1]; v3 = v[2];
    w1 = w[0]; w2 = w[1]; w3 = w[2];

    w11 = w1*w1; w22 = w2*w2; w33 = w3*w3;
    w12 = w1*w2; w13 = w1*w3; w23 = w2*w3;

    u11 = u1*u1; u22 = u2*u2; u33 = u3*u3;
    u12 = u1*u2; u13 = u1*u3; u23 = u2*u3;

    v11 = v1*v1; v22 = v2*v2; v33 = v3*v3;
    v12 = v1*v2; v13 = v1*v3; v23 = v2*v3;


    // Indices (i2m,i2p) and (i3m,i3p) for minus-plus pairs of samples.
    // Cell normal vector w points from the minus side to the plus side.
    i2m = i2p = i2;
    i3m = i3p = i3;
    if (x2>i2) {
      ++i2p;
    } else if (x2<i2) {
      --i2m;
    }
    if (x3>i3) {
      ++i3p;
    } else if (x3<i3) {
      --i3m;
    }
    if ((i2p-i2m)*w2<0.0f) {
      int i2t = i2m; 
      i2m = i2p; 
      i2p = i2t;
    }
    if ((i3p-i3m)*w3<0.0f) {
      int i3t = i3m; 
      i3m = i3p; 
      i3p = i3t;
    }
  }

  private static FaultCell nearestCell(
      FaultCell c1, FaultCell c2, FaultCell c3, 
      float p1, float p2, float p3) 
  {
    float ds1 = distanceSquared(c1,p1,p2,p3);
    float ds2 = distanceSquared(c2,p1,p2,p3);
    float ds3 = distanceSquared(c3,p1,p2,p3);
    float dsm = min(ds1,ds2,ds3);
    if (dsm==ds1) { 
      return c1;
    } else if (dsm==ds2) {
      return c2;
    } else {
      return c3;
    }
  }

  private static float distanceSquared(
      FaultCell c, float p1, float p2, float p3) {
    return c!=null ? c.distanceSquaredTo(p1,p2,p3) : Float.MAX_VALUE;
  }

  /**
   * Rotates a specified point by strike (phi) and dip (theta) angles. Uses
   * specified cosines (cp and ct) and sines (sp and st) of those angles. The
   * order of transformation is
   * (1) rotate around axis x3 by dip angle
   * (2) rotate around axis x1 by strike angle
   * Returns the coordinates of the rotated point.
   */
  private static float[] rotatePoint(
      float cp, float sp, float ct, float st, float[] x) {
    float x1 = x[0], x2 = x[1], x3 = x[2];
    float y1 =     ct*x1+   st*x2;
    float y2 = -cp*st*x1+cp*ct*x2+sp*x3;
    float y3 =  sp*st*x1-sp*ct*x2+cp*x3;
    return new float[]{y1,y2,y3};
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, the specified cells are represented by quads with specified
   * size, and colors corresponding to a property gotten from each cell.
   * @param size the size (in samples) of the quads.
   * @param cmap the colormap used to compute rgb colors from floats.
   * @param cells the cells for which to compute quads.
   * @param get1 used to get the float from a cell to be colormapped.
   * @param lhc true, if left-handed coordinate system; false, otherwise.
   * @return arrays {xyz,uvw,rgb}.
   */
  private static float[][] getXyzUvwRgb(
      float size, ColorMap cmap, FaultCell[] cells, 
      Get1 get1, boolean lhc) {
    FloatList xyz = new FloatList();
    FloatList uvw = new FloatList();
    FloatList fcl = new FloatList();
    size *= 0.5f;
    float[] qa = {0.0f,-size,-size};
    float[] qb = {0.0f, size,-size};
    float[] qc = {0.0f, size, size};
    float[] qd = {0.0f,-size, size};
    if (lhc) {
      float[] qt = qb; qb = qc; qc = qt;
              qt = qa; qa = qd; qd = qt;
    }
    for (FaultCell cell:cells) {
      if(cell==null){continue;}
      float x1 = cell.x1;
      float x2 = cell.x2;
      float x3 = cell.x3;
      float w1 = cell.w1;
      float w2 = cell.w2;
      float w3 = cell.w3;
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
      xyz.add(a3); xyz.add(a2); xyz.add(a1);
      xyz.add(b3); xyz.add(b2); xyz.add(b1);
      xyz.add(c3); xyz.add(c2); xyz.add(c1);
      xyz.add(d3); xyz.add(d2); xyz.add(d1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      uvw.add(w3); uvw.add(w2); uvw.add(w1);
      float fc = get1.get(cell);
      fcl.add(fc);
      fcl.add(fc);
      fcl.add(fc);
      fcl.add(fc);
    }
    float[] fc = fcl.trim();
    float[] rgb = cmap.getRgbFloats(fc);
    return new float[][]{xyz.trim(),uvw.trim(),rgb};
  }

  private static class FloatList {
    public int n = 0;
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

    /**
   * Writes a fault cells to a file with specified name.
   * @param fileName the fault cell file name.
   * @param cells the array of cells.
   */
  public static void writeToFile(String fileName, FaultCell[] cells) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      aos.writeInt(cells.length);
      for (FaultCell cell:cells) {
        aos.writeFloat(cell.x1); 
        aos.writeFloat(cell.x2); 
        aos.writeFloat(cell.x3);
        aos.writeFloat(cell.fl); 
        aos.writeFloat(cell.fp); 
        aos.writeFloat(cell.ft);
      }
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }


  /**
   * Returns fault cells read from a file with specified name.
   * @param fileName the fault cell file name.
   * @return the array of fault cells.
   */
  public static FaultCell[] readFromFile(String fileName) {
    FaultCell[] cells;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      int ncell = ais.readInt();
      cells = new FaultCell[ncell];
      for (int icell=0; icell<ncell; ++icell) {
        float x1 = ais.readFloat();
        float x2 = ais.readFloat();
        float x3 = ais.readFloat();
        float fl = ais.readFloat();
        float fp = ais.readFloat();
        float ft = ais.readFloat();
        cells[icell] = new FaultCell(x1,x2,x3,fl,fp,ft);
      }
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return cells;
  }


}
