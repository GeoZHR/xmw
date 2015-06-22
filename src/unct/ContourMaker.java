/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package unct;

import java.util.ArrayList;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.util.AxisTics;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Clips;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A view of a sampled function f(x1,x2), displayed with contour lines.
 * <p>
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.06.20
 */
public class ContourMaker {

  public ContourMaker(float[][] f) {
    set(f);
  }

  /**
   * Constructs contours of the specified sampled function f(x1,x2).
   * @param s1 the sampling of the variable x1; must be uniform.
   * @param s2 the sampling of the variable x2; must be uniform.
   * @param f array[n2][n1] of sampled function values f(x1,x2), where
   *  n1 and n2 denote the number of samples in s1 and s2, respectively.
   */
  public ContourMaker(Sampling s1, Sampling s2, float[][] f) {
    set(s1,s2,f);
  }

  /**
   * Sets the sampled function f(x1,x2) for this view.
   * Assumes zero first sample values and unit sampling intervals.
   * @param f array[n2][n1] of sampled function values f(x1,x2), where 
   *  n1 = f[0].length and n2 = f.length.
   */
  public void set(float[][] f) {
    set(new Sampling(f[0].length),new Sampling(f.length),f);
  }

  /**
   * Sets the sampled function f(x1,x2) for this view.
   * @param s1 the sampling of the variable x1; must be uniform.
   * @param s2 the sampling of the variable x2; must be uniform.
   * @param f array[n2][n1] of sampled function values f(x1,x2), where
   *  n1 and n2 denote the number of samples in s1 and s2, respectively.
   */
  public void set(Sampling s1, Sampling s2, float[][] f) {
    Check.argument(s1.isUniform(),"s1 is uniform");
    Check.argument(s2.isUniform(),"s2 is uniform");
    Check.argument(isRegular(f),"f is regular");
    Check.argument(s1.getCount()==f[0].length,"s1 consistent with f"); 
    Check.argument(s2.getCount()==f.length,"s2 consistent with f");
    _s1 = s1;
    _s2 = s2;
    _f = copy(f);
    _clips = new Clips(f);
    _cs = null;
    _cl = null;
  }

 
  /**
   * Enables or disables automatically computed readable contour values. 
   * Here, readable values are multiples of 1, 2, and 5 times some 
   * power of ten. If enabled, then any specified number of contours
   * serves as an upper bound on the number of contour values.
   * Readable contour values are the default.
   * @param readableContours true, for readable contours; false, otherwise.
   */
  public void setReadableContours(boolean readableContours) {
    if (_readableContours!=readableContours) {
      _readableContours = readableContours;
      _cs = null;
      _cl = null;
    }
  }

  /**
   * Sets the number of contour values. 
   * If readable contours are enabled, then this number is an upper bound,
   * and the actual number of contours may be less than this number.
   * Otherwise, if readable contours are disabled, then contour values 
   * will be evenly spaced between, but not including, the minimum and 
   * maximum clip values. The default number of contour values is 25.
   * @param n the number of contour values.
   */
  public void setContours(int n) {
    _nc = n;
    _cs = null;
    _cl = null;
  }

  /**
   * Sets the contour values to those in the specified array.
   * If this method is called, then clips (or percentiles) are not used 
   * to determine contour values, and readable contours are disabled.
   * @param c the array of contour values.
   */
  public void setContours(float[] c) {
    double[] cd = new double[c.length];
    for (int i=0; i<c.length; ++i) 
      cd[i] = (double)c[i];
    setContours(new Sampling(cd));
  }

  /**
   * Sets the contour values to the specified sampling.
   * If this method is called, then clips (or percentiles) are not used 
   * to determine contour values, and readable contours are disabled.
   * @param cs the contour sampling.
   */
  public void setContours(Sampling cs) {
    _readableContours = false;
    _cs = cs;
    _cl = null;
  }

  /**
   * Gets the contour values.  
   * @return array of contour values.
   */
  public float[] getContours() {
    updateContourSampling();
    float[] values = new float[_cs.getCount()];
    for (int n=0; n<values.length; n++)
      values[n] = _cl.get(n).fc;
    return values;
  }

  public Contour getContour(float fc) {
    return makeContour(fc,_s1,_s2,_f); 
  }

  /**
   * Sets the clips for this view. These values limit the range used
   * to determine contour values. Function values f(x1,x2) less than
   * clipMin and greater than clipMax are ignored.
   * <p>
   * Calling this method disables the computation of clips from percentiles.
   * Previous clip values will be forgotten.
   * @param clipMin the lower bound on contour values.
   * @param clipMax the upper bound on contour values.
   */
  public void setClips(float clipMin, float clipMax) {
    _clips.setClips(clipMin,clipMax);
    _cs = null;
    _cl = null;
  }

  /**
   * Gets the minimum clip value.
   * @return the minimum clip value.
   */
  public float getClipMin() {
    return _clips.getClipMin();
  }

  /**
   * Gets the maximum clip value.
   * @return the maximum clip value.
   */
  public float getClipMax() {
    return _clips.getClipMax();
  }

  /**
   * Sets the percentiles used to compute clips for this view. The default 
   * percentiles are 0 and 100, which correspond to the minimum and maximum 
   * values of the sampled function f(x1,x2).
   * <p>
   * Calling this method enables the computation of clips from percentiles.
   * Any clip values specified or computed previously will be forgotten.
   * @param percMin the percentile corresponding to clipMin.
   * @param percMax the percentile corresponding to clipMax.
   */
  public void setPercentiles(float percMin, float percMax) {
    _clips.setPercentiles(percMin,percMax);
    _cs = null;
    _cl = null;
  }

  /**
   * Gets the minimum percentile.
   * @return the minimum percentile.
   */
  public float getPercentileMin() {
    return _clips.getPercentileMin();
  }

  /**
   * Gets the maximum percentile.
   * @return the maximum percentile.
   */
  public float getPercentileMax() {
    return _clips.getPercentileMax();
  }

  
  ///////////////////////////////////////////////////////////////////////////
  // private

  public class Contour {
    public float fc; // contoured function value
    public int ns = 0; // number of segments
    public ArrayList<float[]> x1 = new ArrayList<float[]>(); // x1[] per segment
    public ArrayList<float[]> x2 = new ArrayList<float[]>(); // x2[] per segment
    Contour (float fc) {
      this.fc = fc;
    }
    void append(FloatList x1List, FloatList x2List) {
      System.out.println("test!!!!!!");
      ++this.ns;
      this.x1.add(x1List.trim());
      this.x2.add(x2List.trim());
    }
  }

  private class FloatList {
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

  // Contour line attributes.
  private float _lineWidth = 0.0f;

  // The sampled floats.
  private Sampling _s1; // sampling of 1st dimension
  private Sampling _s2; // sampling of 2nd dimension
  private float[][] _f; // copy of array of floats

  // Sampling of the function f(x1,x2) in the pixel (x,y) coordinate system. 
  private boolean _transposed; // true, if (x,y) <=> (x2,x1)
  private boolean _xflipped; // true, if axis decreases with increasing x
  private boolean _yflipped; // true, if axis decreases with increasing y
  private int _nx; // number of samples for x axis
  private double _dx; // sampling interval for x axis
  private double _fx; // first sample for x axis
  private int _ny; // number of samples for y axis
  private double _dy; // sampling interval for y axis
  private double _fy; // first sample for y axis

  // Clipping
  private Clips _clips;
  private float _clipMin;
  private float _clipMax;

  // Contour sampling and list of contours.
  private int _nc = 25; // number of contours; maybe less if readable contours
  private boolean _readableContours = true; // true, for readable contour vals
  private Sampling _cs; // contour sampling
  private ArrayList<Contour> _cl; // list of contours
  
  /**
   * Update the clips if necessary.
   */
  private void updateClips() {
    float clipMin = _clips.getClipMin();
    float clipMax = _clips.getClipMax();
    if (_clipMin!=clipMin || _clipMax!=clipMax) {
      _clipMin = clipMin;
      _clipMax = clipMax;
    }
  }

  /**
   * Updates the number of contours. This is based on the contour sampling,
   * which is determined by the clipping values of the data range.  Since
   * the defaults are 0 and 100, the entire range of data will be used
   * to determine the contour values.  If the clips are set, the contours
   * will be set to the new range.
   */
  private void updateContourSampling() {
    if (_cs==null) {
      updateClips();
      int nc;
      double dc,fc;
      if (_readableContours) {
        AxisTics at = new AxisTics(_clipMin,_clipMax,_nc);
        nc = at.getCountMajor();
        dc = at.getDeltaMajor();
        fc = at.getFirstMajor();
      } else {
        nc = _nc;
        dc = (_clipMax-_clipMin)/(nc+1);
        fc = _clipMin;
      }
      double[] cstep = new double[nc];
      cstep[0] = fc+dc;
      int count = 1;
      while (count<nc) {
        cstep[count] = cstep[count-1]+dc;
        count++;
      }
      _cs = new Sampling(cstep);
    }
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
  private Contour makeContour(float fc, 
    Sampling s1, Sampling s2, float[][] f)
  {
    System.out.println("test!!!!!!");
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

    System.out.println("test!!!!!!");
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
  private int connect(
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
  private final byte WEST = 1;
  private final byte SOUTH = 2;
  private final byte NOT_WEST = ~WEST;
  private final byte NOT_SOUTH = ~SOUTH;

  private void sets(int i1, int i2, byte[][] flags) {
    flags[i2][i1] |= SOUTH;
  }
   
  private void setw(int i1, int i2, byte[][] flags) {
    flags[i2][i1] |= WEST;
  }
 
  private void clrs(int i1, int i2, byte[][] flags) {
    flags[i2][i1] &= NOT_SOUTH;
  }
 
  private void clrw(int i1, int i2, byte[][] flags) {
    flags[i2][i1] &= NOT_WEST;
  }

  private boolean sset(int i1, int i2, byte[][] flags) {
    return (flags[i2][i1]&SOUTH)!=0;
  }

  private boolean wset(int i1, int i2, byte[][] flags) {
    return (flags[i2][i1]&WEST)!=0;
  }

  private boolean between(float fc, float f1, float f2) {
    return (f1<=f2) ? f1<=fc && fc<f2 : f2<=fc && fc<f1;
  }

  private float delta(float fc, float f1, float f2) {
    return (f1!=f2)?(fc-f1)/(f2-f1):1.0f;
  }
}
