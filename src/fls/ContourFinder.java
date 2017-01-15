/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fls;

import java.util.*;
import edu.mines.jtk.dsp.Sampling;

///////////////////////////////////////////////////////////////////////////
// Comments in functions below refer to cells with indices (i1,i2). Each 
// cell has north, east, south, and west boundaries, defined as shown here:
// 
//              north
//  (i1,i2+1)	--------- (i1+1,i2+1)
//            | cell  |
//      west  | i1,i2	| east
//            |       |
//    (i1,i2) --------- (i1+1,i2)
//              south
///////////////////////////////////////////////////////////////////////////
public class ContourFinder {

  public static class Contour {
    float fc; // contoured function value
    int ns = 0; // number of segments
    ArrayList<float[]> x1 = new ArrayList<float[]>(); // x1[] per segment
    ArrayList<float[]> x2 = new ArrayList<float[]>(); // x2[] per segment
    Contour (float fc) {
      this.fc = fc;
    }
    public float[][][] getCoords() {
      int nc = x1.size();
      float[][][] xs = new float[nc][2][];
      for (int ic=0; ic<nc; ++ic) {
        xs[ic][0] = x1.get(ic);
        xs[ic][1] = x2.get(ic);
      }
      return xs;
    }
    void append(FloatList x1List, FloatList x2List) {
      ++this.ns;
      this.x1.add(x1List.trim());
      this.x2.add(x2List.trim());
    }
  }

  private static class FloatList {
    public int n;
    public float[] a = new float[64];
    public void add(double f) {
      if (n==a.length) {
        float[] t = new float[2*a.length];
        for (int i=0; i<n; ++i)
          t[i] = a[i];
        a = t;
      }
      a[n++] = (float)f;
    }
    public float[] trim() {
      float[] t = new float[n];
      for (int i=0; i<n; ++i)
        t[i] = a[i];
      return t;
    }
  }



  /**
   * Returns a new contour for one specified function value fc.
   */
  public Contour findContour(float fc, 
    Sampling s1, Sampling s2, float[][] f)
  {
    int n1 = s1.getCount();
    double d1 = s1.getDelta();
    double f1 = s1.getFirst();
    int n2 = s2.getCount();
    double d2 = s2.getDelta();
    double f2 = s2.getFirst();
    int n1m1 = n1-1;
    int n2m1 = n2-1;
    int i1,i2,is,i;
    
    // Mark and count intersections with west and south edges of cells.
    byte[][] flags = new byte[n2][n1];
    int ni = 0;
    for (i2=0; i2<n2; ++i2) {
      for (i1=0; i1<n1; ++i1) {
        if (i2<n2m1 && between(fc,f[i2][i1],f[i2+1][i1])) {
          setw(i1,i2,flags);
          ++ni;
        }
        if (i1<n1m1 && between(fc,f[i2][i1],f[i2][i1+1])) {
          sets(i1,i2,flags);
          ++ni;
        }
      }
    }

    // Construct new contour.
    Contour c = new Contour(fc);

    // Append contour segments intersecting north boundary of grid.
    i2 = n2m1;
    for (i1=0,is=i1+i2*n1; i1<n1m1 && ni>0; ++i1,is+=1) {
      if (sset(i1,i2,flags)) {
        float d = delta(fc,f[i2][i1],f[i2][i1+1]);
        FloatList x1 = new FloatList();
        FloatList x2 = new FloatList();
        x1.add(f1+(i1+d)*d1);
        x2.add(f2+(i2  )*d2);
        clrs(i1,i2,flags);
        for (i=is-n1; i>=0; i=connect(i,fc,n1,d1,f1,n2,d2,f2,f,flags,x1,x2))
          --ni;
        c.append(x1,x2);
      }
    }

    // Append contour segments intersecting east boundary of grid.
    i1 = n1m1;
    for (i2=0,is=i1+i2*n1; i2<n2m1 && ni>0; ++i2,is+=n1) {
      if (wset(i1,i2,flags)) {
        float d = delta(fc,f[i2][i1],f[i2+1][i1]);
        FloatList x1 = new FloatList();
        FloatList x2 = new FloatList();
        x1.add(f1+(i1  )*d1);
        x2.add(f2+(i2+d)*d2);
        clrw(i1,i2,flags);
        for (i=is-1; i>=0; i=connect(i,fc,n1,d1,f1,n2,d2,f2,f,flags,x1,x2))
          --ni;
        c.append(x1,x2);
      }
    }

    // Append contour segments intersecting south boundary of grid.
    i2 = 0;
    for (i1=0,is=i1+i2*n1; i1<n1m1 && ni>0; ++i1,is+=1) {
      if (sset(i1,i2,flags)) {
        float d = delta(fc,f[i2][i1],f[i2][i1+1]);
        FloatList x1 = new FloatList();
        FloatList x2 = new FloatList();
        x1.add(f1+(i1+d)*d1);
        x2.add(f2+(i2  )*d2);
        clrs(i1,i2,flags);
        for (i=is; i>=0; i=connect(i,fc,n1,d1,f1,n2,d2,f2,f,flags,x1,x2))
          --ni;
        c.append(x1,x2);
      }
    }

    // Append contour segments intersecting west boundary of grid.
    i1 = 0;
    for (i2=0,is=i1+i2*n1; i2<n2m1 && ni>0; ++i2,is+=n1) {
      if (wset(i1,i2,flags)) {
        float d = delta(fc,f[i2][i1],f[i2+1][i1]);
        FloatList x1 = new FloatList();
        FloatList x2 = new FloatList();
        x1.add(f1+(i1  )*d1);
        x2.add(f2+(i2+d)*d2);
        clrw(i1,i2,flags);
        for (i=is; i>=0; i=connect(i,fc,n1,d1,f1,n2,d2,f2,f,flags,x1,x2))
          --ni;
        c.append(x1,x2);
      }
    }

    // Append contour segments intersecting interior cells.
    for (i2=1; i2<n2m1 && ni>0; ++i2) {
      for (i1=0,is=i1+i2*n1; i1<n1m1 && ni>0; ++i1,++is) {
        if (sset(i1,i2,flags)) {
          float d = delta(fc,f[i2][i1],f[i2][i1+1]);
          FloatList x1 = new FloatList();
          FloatList x2 = new FloatList();
          x1.add(f1+(i1+d)*d1);
          x2.add(f2+(i2  )*d2);
          clrs(i1,i2,flags);
          for (i=is; i>=0; i=connect(i,fc,n1,d1,f1,n2,d2,f2,f,flags,x1,x2))
            --ni;
          // Close the contours...
          x1.add(x1.a[0]);
          x2.add(x2.a[0]);
          c.append(x1,x2);
        } 
      }
    }

    return c;
  }

  /**
   * Connects two intersections of a contour for a cell (i1,i2), if possible.
   * When called, the index = i1+i2*n1 points to a cell (i1,i2) for which
   * one intersection has already been found and cleared. This method looks
   * for another intersection. If another intersection is found, this method 
   * clears it, appends the intersection coordinates (x,y) to the specified 
   * lists, and returns a modified index that indicates an adjacent cell. If 
   * another intersection is not found, or if the found intersection is a grid 
   * boundary, then this method returns -1.
   */
  private static int connect(
    int index, float fc,
    int n1, double d1, double f1, int n2, double d2, double f2, float[][] f, 
    byte[][] flags, FloatList x1, FloatList x2)
  {
    int i1 = index%n1;
    int i2 = index/n1;

    // If exiting north, ...
    if (sset(i1,i2+1,flags)) {
      float d = delta(fc,f[i2+1][i1],f[i2+1][i1+1]);
      x1.add(f1+(i1+d)*d1);
      x2.add(f2+(i2+1)*d2);
      clrs(i1,++i2,flags);
      return (i2<n2-1)?index+n1:-1;
    } 
    
    // Else if exiting east, ...
    else if (wset(i1+1,i2,flags)) {
      float d = delta(fc,f[i2][i1+1],f[i2+1][i1+1]);
      x1.add(f1+(i1+1)*d1);
      x2.add(f2+(i2+d)*d2);
      clrw(++i1,i2,flags);
      return (i1<n1-1)?index+1:-1;
    }

    // Else if exiting south, ...
    else if (sset(i1,i2,flags)) {
      float d = delta(fc,f[i2][i1],f[i2][i1+1]);
      x1.add(f1+(i1+d)*d1);
      x2.add(f2+(i2  )*d2);
      clrs(i1,i2,flags);
      return (i2>0)?index-n1:-1;
    } 
    
    // Else if exiting west, ...
    else if (wset(i1,i2,flags)) {
      float d = delta(fc,f[i2][i1],f[i2+1][i1]);
      x1.add(f1+(i1  )*d1);
      x2.add(f2+(i2+d)*d2);
      clrw(i1,i2,flags);
      return (i1>0)?index-1:-1;
    }
    
    // Else if no intersection exists, ...
    else {
      return -1;
    }
  }

  // Contour exit-directions.
  private static final byte WEST = 1;
  private static final byte SOUTH = 2;
  private static final byte NOT_WEST = ~WEST;
  private static final byte NOT_SOUTH = ~SOUTH;

  private static void sets(int i1, int i2, byte[][] flags) {
    flags[i2][i1] |= SOUTH;
  }
   
  private static void setw(int i1, int i2, byte[][] flags) {
    flags[i2][i1] |= WEST;
  }
 
  private static void clrs(int i1, int i2, byte[][] flags) {
    flags[i2][i1] &= NOT_SOUTH;
  }
 
  private static void clrw(int i1, int i2, byte[][] flags) {
    flags[i2][i1] &= NOT_WEST;
  }

  private static boolean sset(int i1, int i2, byte[][] flags) {
    return (flags[i2][i1]&SOUTH)!=0;
  }

  private static boolean wset(int i1, int i2, byte[][] flags) {
    return (flags[i2][i1]&WEST)!=0;
  }

  private static boolean between(float fc, float f1, float f2) {
    return (f1<=f2) ? f1<=fc && fc<f2 : f2<=fc && fc<f1;
  }

  private static float delta(float fc, float f1, float f2) {
    return (f1!=f2)?(fc-f1)/(f2-f1):1.0f;
  }
}
