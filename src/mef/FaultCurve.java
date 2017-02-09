package mef;

import java.io.*;
import java.util.*;

import edu.mines.jtk.io.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A linked list of fault points that may be used to analyze faults.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.21
 */


public class FaultCurve implements Iterable<FaultPoint>,Serializable {
  private static final long serialVersionUID = 1L;


    /**
   * Gets the point that was the seed used to grow this curve.
   * @return the seed cell.
   */
  public FaultPoint getSeed() {
    return _seed;
  }

  /**
   * Returns the number of points in this curve.
   */
  public int size() {
    return _pointList.size();
  }

  /**
   * Gets an array of points in this curve.
   * @return array of points.
   */
  public FaultPoint[] getPoints() {
    return _pointList.toArray(new FaultPoint[0]);
  }

  /**
   * Gets all cells in the specified curves.
   * @param curves array of curves for which to get points.
   * @return array of cells.
   */
  public static FaultPoint[] getPoints(FaultCurve[] curves) {
    int npoint = countPoints(curves);
    FaultPoint[] points = new FaultPoint[npoint];
    int ipoint = 0;
    for (FaultCurve curve:curves)
      for (FaultPoint point:curve)
        points[ipoint++] = point;
    return points;
  }

  public static void getFlImage(FaultCurve[] curves, float[][] fl) {
    for (FaultCurve curve:curves) {
    for (FaultPoint point:curve) {
      int i1 = point.i1;
      int i2 = point.i2;
      fl[i2][i1] = point.fl;
    }}
  }

  public static void getFlsImage(FaultCurve[] curves, float[][] fl) {
    for (FaultCurve curve:curves) {
    for (FaultPoint point:curve) {
      int i1 = point.i1;
      int i2 = point.i2;
      int i2p = point.i2p;
      fl[i2][i1] = point.fl;
      fl[i2p][i1] = point.fl;
    }}
  }


  public static void getFtImage(FaultCurve[] curves, float[][] ft) {
    for (FaultCurve curve:curves) {
    for (FaultPoint point:curve) {
      int i1 = point.i1;
      int i2 = point.i2;
      ft[i2][i1] = point.ft;
    }}
  }

  public static void getFsImage(FaultCurve[] curves, float[][] fs) {
    for (FaultCurve curve:curves) {
    for (FaultPoint point:curve) {
      int i1 = point.i1;
      int i2m = point.i2m;
      int i2p = point.i2p;
      fs[i2m][i1] = point.smp;
      fs[i2p][i1] = point.smp;
    }}
  }

  public static void getFsImage(
    FaultCurve[] curves, float[][] fs1, float[][] fs2) 
  {
    for (FaultCurve curve:curves) {
    for (FaultPoint point:curve) {
      int i1 = point.i1;
      int i2 = point.i2;
      fs1[i2][i1] = point.s1;
      fs2[i2][i1] = point.s2;
    }}
  }


  public static void getSlipImage(FaultCurve[] curves, float[][] fd) {
    for (FaultCurve curve:curves) {
    for (FaultPoint point:curve) {
      int i1 = point.i1;
      int i2 = point.i2;
      fd[i2][i1]  = point.s1/point.s2;
    }}
  }




  /**
   * Returns the total number of points in the specified curves.
   * @param curves array of curves for which to count points.
   * @return the total number of points.
   */
  public static int countPoints(FaultCurve[] curves) {
    int npoint = 0;
    for (FaultCurve curve:curves)
      npoint += curve.size();
    return npoint;
  }

    /**
   * Returns array of arrays of cells linked above and below.
   * @return array of arrays of linked cells; by reference, not by copy.
   */
  public FaultPoint[][] getCellsAB() {
    if (_pointsAB!=null)
      return _pointsAB;

    HashSet<FaultPoint> cellSet = new HashSet<FaultPoint>(size());
    ArrayList<FaultPoint[]> cellsList = new ArrayList<FaultPoint[]>();

    // For all cells in this skin, ...
    for (FaultPoint cell:_pointList) {

      // If the cell is not already in an array, ...
      if (!cellSet.contains(cell)) {

        // Search above for the top cell.
        FaultPoint c = cell;
        for (FaultPoint ca=c.ca; ca!=null; ca=c.ca)
          c = ca;

        // Add the top cell and all cells below it to a new list.
        ArrayList<FaultPoint> cList = new ArrayList<FaultPoint>();
        for (; c!=null; c=c.cb) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the list of arrays.
        cellsList.add(cList.toArray(new FaultPoint[0]));
      }
    }
    assert _pointList.size()==cellSet.size();

    // Convert the list of arrays to the array of arrays to be returned.
    _pointsAB = cellsList.toArray(new FaultPoint[0][]);
    return _pointsAB;
  }



  /**
   * Returns an iterator for the points in this curve. 
   * @return cell iterator.
   */
  public Iterator<FaultPoint> iterator() {
    return _pointList.iterator();
  }

    /////////////////////////////////////////////////////////////////////////
  // package

  /**
   * Constructs an empty curve.
   */
  FaultCurve() {
    _pointList = new ArrayList<FaultPoint>();
  }

  /**
   * Adds the specified fault point to this fault curve.
   * @param point the fault point to be added.
   */
  void add(FaultPoint point) {
    assert point.curve==null;
    point.curve = this;
    if (_seed==null)
      _seed = point;
    _pointList.add(point);
    _pointsAB = null;
  }


    /**
   * Smooths one value stored in the cells of this skin. The value smoothed is
   * that accessed by the specified getter and setter. Each smoothed value is
   * an average of the values in a cell and its cell nabors. 
   */
  void smooth1(FaultPoint.Get1 getter, FaultPoint.Set1 setter) {
    int ncell = size();
    float[] vals = new float[ncell];
    float[] cnts = new float[ncell];
    FaultPoint[] cellNabors = new FaultPoint[2];
    for (int icell=0; icell<ncell; ++icell) {
      FaultPoint cell = _pointList.get(icell);
      float valCell = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      for (FaultPoint cellNabor:cellNabors) {
        if (cellNabor!=null) {
          float valNabor = getter.get(cellNabor);
          vals[icell] += valCell+valNabor;
          cnts[icell] += 2.0f;
        }
      }
    }
    for (int icell=0; icell<ncell; ++icell) {
      FaultPoint cell = _pointList.get(icell);
      float cnti = cnts[icell];
      float vali = vals[icell]/(cnti>0.0f?cnti:1.0f);
      setter.set(cell,vali);
    }
  }



  /////////////////////////////////////////////////////////////////////////
  // private

  private static final int INULL = -Integer.MAX_VALUE; // null index

  private FaultPoint _seed; // cell in this curve with highest fl; null, if empty
  private ArrayList<FaultPoint> _pointList; // list of points in this curve
  private FaultPoint[][] _pointsAB; // arrays of points from above to below

}


