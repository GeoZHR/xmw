/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package sso;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

import java.util.*;
import util.*;


/**
 * Create meshes of surfaces for displaying.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.02.18
 */

public class HorizonDisplay {

  public float[][][] horizonCurves(
    int k2, int k3, float[][][] sfs) 
  {
    int n1s = sfs.length;
    int n3 = sfs[0].length;
    int n2 = sfs[0][0].length;
    float[][][] xbs = new float[n1s*200][2][];
    float[][] surfs = new float[n1s*200][max(n2,n3)*3];
    float rgtMin = 0;
    float rgtMax = n1s;
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);
    int k = 0;
    for (int i1s=4; i1s<n1s; ++i1s) {
      int p = 0;
      for (int i3=0; i3<n3; ++i3) {
        float x3 = i3;
        float x2 = k2;
        float x1 = sfs[i1s][i3][k2];
        surfs[k][p++] = x3;
        surfs[k][p++] = x2;
        surfs[k][p++] = x1;
      }
      float[] vs = fillfloat(i1s,p);
      xbs[k][0] = copy(p,surfs[k]);
      xbs[k][1] = cp.getRgbFloats(vs);
      if(p>6){k ++;}
      p = 0;
      for (int i2=0; i2<n2; ++i2) {
        float x3 = k3;
        float x2 = i2;
        float x1 = sfs[i1s][k3][i2];
        surfs[k][p++] = x3;
        surfs[k][p++] = x2;
        surfs[k][p++] = x1;
      }
      vs = fillfloat(i1s,p);
      xbs[k][0] = copy(p,surfs[k]);
      xbs[k][1] = cp.getRgbFloats(vs);
      if(p>6){k ++;}
      p = 0;
    }
    float[][][] hzs = new float[k][2][];
    for (int i=0; i<k; ++i) {
      int np = xbs[i][0].length;
      hzs[i][0] = copy(np,0,xbs[i][0]); 
      hzs[i][1] = copy(np,0,xbs[i][1]); 
    }
    return hzs;
  }


  public float[][][] slice12(int k3, Sampling s2, float[][][] hv) {
    int ns = hv.length;
    int n2 = hv[0][0].length;
    float f2 = (float)s2.getFirst();
    float d2 = (float)s2.getDelta();
    float[][][] cs = new float[ns][2][n2];
    for (int is=0; is<ns; ++is) {
      cs[is][0] = hv[is][k3];
      cs[is][1] = rampfloat(f2,d2,n2);
    }
    return cs;
  }

  public float[][][] slice13(int k2, Sampling s3, float[][][] hv) {
    int ns = hv.length;
    int n3 = hv[0].length;
    float f3 = (float)s3.getFirst();
    float d3 = (float)s3.getDelta();
    float[][][] cs = new float[ns][2][n3];
    for (int is=0; is<ns; ++is) {
      cs[is][1] = rampfloat(f3,d3,n3);
      for (int i3=0; i3<n3; ++i3) 
        cs[is][0][i3] = hv[is][i3][k2];
    }
    return cs;
  }

  public float[][][] slice23(
    int k1, Sampling s2, Sampling s3, float[][][] hv) 
  {
    int ns = hv.length;
    int n3 = hv[0].length;
    int n2 = hv[0][0].length;
    float[][][] cvs = new float[ns][2][];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    SincInterpolator si = new SincInterpolator();
    for (int is=0; is<ns; ++is) {
      float[][] g2 = new float[n3][n2];
      float[][] g3 = new float[n3][n2];
      rgf.apply10(hv[is],g2);
      rgf.apply01(hv[is],g3);
      float[][] hd = sub(hv[is],k1);
      ArrayList<Float> x2 = new ArrayList<Float>();
      ArrayList<Float> x3 = new ArrayList<Float>();
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float g2i = g2[i3][i2];
        float g3i = g3[i3][i2];
        float gsi = -1f/sqrt(g2i*g2i+g3i*g3i+1);
        g2i *= gsi;
        g3i *= gsi;
        float x2p = i2+g2i;
        float x3p = i3+g3i;
        float x2m = i2-g2i;
        float x3m = i3-g3i;
        float hdm = si.interpolate(n2,1,0,n3,1,0,hd,x2m,x3m);
        float hdp = si.interpolate(n2,1,0,n3,1,0,hd,x2p,x3p);
        if (hdm*hdp<0f) {
          x2.add((float)s2.getValue(i2));
          x3.add((float)s3.getValue(i3));
        }
      }}
      int np = x2.size();
      cvs[is][0] = new float[np];
      cvs[is][1] = new float[np];
      for (int ip=0; ip<np; ++ip) {
        cvs[is][0][ip] = x2.get(ip);
        cvs[is][1][ip] = x3.get(ip);
      }
    }
    return cvs;
  }

  public float[][][][] slice23X(
    int k1, Sampling s2, Sampling s3, float[][][] hv) 
  {
    int ns = hv.length;
    float[][][][] cvs = new float[ns][2][][];
    for (int is=0; is<ns; ++is) {
      float[][] hd = sub(hv[is],k1);
      Contour ct = makeContour(0f,s2,s3,hd);
      ArrayList<float[]> x2s = ct.x1;
      ArrayList<float[]> x3s = ct.x2;
      int nc = x2s.size();
      cvs[is][0] = new float[nc][];
      cvs[is][1] = new float[nc][];
      for (int ic=0; ic<nc; ++ic) {
        float[] x2a = x2s.get(ic);
        float[] x3a = x3s.get(ic);
        int np = x2a.length;
        cvs[is][0][ic] = new float[np];
        cvs[is][1][ic] = new float[np];
        for (int ip=0; ip<np; ++ip) {
          cvs[is][0][ic][ip] = x2a[ip];
          cvs[is][1][ic][ip] = x3a[ip];
        }
      }
    }
    return cvs;
  }



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

  /**
   * Returns a new contour for one specified function value fc.
   */
  private static Contour makeContour(float fc, 
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


    ///////////////////////////////////////////////////////////////////////////
  // private

  private static class Contour {
    float fc; // contoured function value
    int ns = 0; // number of segments
    ArrayList<float[]> x1 = new ArrayList<float[]>(); // x1[] per segment
    ArrayList<float[]> x2 = new ArrayList<float[]>(); // x2[] per segment
    Contour (float fc) {
      this.fc = fc;
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
