package mef;

import java.util.*;

import edu.mines.jtk.util.*;

import util.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Computes fault curves from images of fault likelihoods and dips. 
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.07.21
 */
public class FaultCurver {

  /**
   * Constructs a fault curver with default parameters.
   */
  public FaultCurver() {
    _fs1min = -Float.MAX_VALUE;
    _fs1max =  Float.MAX_VALUE;
    _fllo = 0.2f;
    _flhi = 0.8f;
    _dflmax = 0.2f;
    _dftmax = 10.0f;
    _ds1max = 1.0f;
    _ncsmin = 400;
  }


    /**
   * Sets the minimum number of points in a curve. curves smaller than this will
   * be discarded.
   * <p>
   * The default minimum curve size is 400.
   */
  public void setMinCurveSize(int minSize) {
    _ncsmin = minSize;
  }

  /**
   * Sets lower and upper bounds on fault throw for points in a curve.
   * These bounds should be set only after fault dip slips have been 
   * computed for all points used to grow curves.
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
   * Sets fault likelihood thresholds used to grow curves. Points in a curve
   * should have, or be connected to points that have, high fault likelihoods.
   * All points in a curve will have fault likelihoods not less than the lower
   * threshold. At least one point in a curve will have a fault likelihood not
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
   * Sets the maximum difference in fault likelihood for a point and its nabors.
   * <p>
   * The default maximum difference is 0.2.
   * @param maxDeltaLikelihood upper bound on difference in fault likelihood.
   */
  public void setMaxDeltaLikelihood(double maxDeltaLikelihood) {
    _dflmax = (float)maxDeltaLikelihood;
  }

  /**
   * Sets the maximum difference in fault dip for a point and its nabors.
   * @param maxDeltaDip upper bound on difference in fault dip.
   * <p>
   * The default maximum difference is 10 degrees.
   */
  public void setMaxDeltaDip(double maxDeltaDip) {
    _dftmax = (float)maxDeltaDip;
  }

  /**
   * Sets the maximum difference in fault throw for a point and its nabors.
   * @param maxDeltaThrow upper bound on difference in fault throw.
   * <p>
   * The default maximum difference is 1.0 samples.
   */
  public void setMaxDeltaThrow(double maxDeltaThrow) {
    _ds1max = (float)maxDeltaThrow;
  }


  /**
   * Returns array of points in ridge surfaces of fault likelihood.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   * @return array of points.
   */
  public FaultPoint[] findPoints(float[][][] flt) {
    return points(flt);
  }


  /**
   * Returns an array of curves comprised of specified points. Some points may be
   * unused. For example, points with fault likelihoods less than the lower
   * threshold for growing curves will be unused. Likewise, points that do not
   * form a curve with sufficient size will be unused.
   * @param points array of points from which to grow curves.
   * @return array of curves.
   */
  public FaultCurve[] findCurves(FaultPoint[] points) {
    return curves(points);
  }





  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _fllo; // lower threshold on fault likelihoods
  private float _flhi; // higher threshold on fault likelihoods
  private float _fs1min; // min fault throw
  private float _fs1max; // max fault throw
  private float _dflmax; // max difference between likelihoods of nabors
  private float _dftmax; // max difference between dips of nabors
  private float _ds1max; // max difference between throws of nabors
  private int _ncsmin; // min number of points that form a curve


    // Uses fault images to find points, oriented points located on ridges.
  private FaultPoint[] points(float[][][] flt) {
    float[][] f = flt[0];
    float[][] t = flt[1];
    int n1 = f[0].length;
    int n2 = f.length;

    // Smooth fault likelihoods in 2nd and 3rd dimensions. This helps to
    // eliminate spurious ridges, and improves the accuracy of 2nd-order
    // finite-difference approximations (parabolic interpolation) used to
    // locate ridges.
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
    rgf.applyX0(f,fs);
    f = fs;

    // Vertical image boundaries are discontinuities that may look like
    // faults. If a fault appears to be near and nearly parallel to image
    // boundaries, then assume it is a boundary artifact and not truly a
    // fault.
    int imax = 5; // max number of samples considered to be near boundary

    // Loop over all samples. Construct points for samples nearest to ridges.
    ArrayList<FaultPoint> pointList = new ArrayList<FaultPoint>();
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float[] fi = f[i2 ];
      float[] fm = f[i2m];
      float[] fp = f[i2p];
      float[] ti = t[i2 ];
      for (int i1=0; i1<n1; ++i1) {
        float fii = fi[i1 ];
        float fmi = fm[i1 ];
        float fpi = fp[i1 ];
        float tii = ti[i1 ];

        // Most image samples will not have a fault point.
        FaultPoint point = null;

        if (fpi<fii && fmi<fii) {
          float f1 = 0.5f*(fpi-fmi); // 1st derivative
          float f2 = fpi-2.0f*fii+fmi; // 2nd derivative
          float dr = -f1/f2; // signed distance to ridge
          float fr = fii+f1*dr+0.5f*f2*dr*dr; // fault likelihood
          if (fr>=_fllo) {
            if (imax<=i2 && i2<n2-imax) {
              point = new FaultPoint(i1,i2+dr,fr,tii);
              pointList.add(point);
            }
          }
        }
      }
    }
    return pointList.toArray(new FaultPoint[0]);
  }


    // Returns curves constructed from specified points.
  private FaultCurve[] curves(FaultPoint[] points) {
    int npoint = points.length;

    // Grid of points used to quickly find point nabors.
    FaultPointGrid pointGrid = new FaultPointGrid(points);

    // Empty list of curves.
    ArrayList<FaultCurve> curveList = new ArrayList<FaultCurve>();

    // point comparator for high-to-low ordering based on fault likelihoods.
    Comparator<FaultPoint> flComparator = new Comparator<FaultPoint>() {
      public int compare(FaultPoint c1, FaultPoint c2) {
        if (c1.fl<c2.fl)
          return 1;
        else if (c1.fl>c2.fl)
          return -1;
        else
          return 0;
      }
    };

    // Make a list of points that might be seeds for new curves.
    ArrayList<FaultPoint> seedList = new ArrayList<FaultPoint>();
    for (int ipoint=0; ipoint<npoint; ++ipoint) {
      FaultPoint point = points[ipoint];
      if (point.fl>=_flhi && point.s1>=_fs1min && point.s1<=_fs1max)
        seedList.add(point);
    }
    int nseed = seedList.size();

    // Sort the list of seeds high-to-low by fault likelihood.
    FaultPoint[] seeds = seedList.toArray(new FaultPoint[0]);
    Arrays.sort(seeds,flComparator);
    seedList.clear();
    for (FaultPoint seed:seeds)
      seedList.add(seed);

    // While potential seeds remain, ...
    for (int kseed=0; kseed<nseed; ++kseed) {

      // Skip any potential seeds that are already in a curve.
      while (kseed<nseed && seedList.get(kseed).curve!=null)
        ++kseed;

      // If we found a seed with which to construct a new curve, ...
      if (kseed<nseed) {
        FaultPoint seed = seedList.get(kseed);

        // Make a new empty curve.
        FaultCurve curve = new FaultCurve();

        // Make a priority queue of points, initially with only the seed.
        PriorityQueue<FaultPoint> growQueue = 
            new PriorityQueue<FaultPoint>(1024,flComparator);
        growQueue.add(seed);

        // While the grow queue is not empty, ...
        while (!growQueue.isEmpty()) {

          // Get and remove the point with highest fault likelihood from the
          // grow queue. If not already in the curve, add them and link and
          // add any mutually best nabors to the grow queue.
          FaultPoint point = growQueue.poll();
          if (point.curve==null) {
            curve.add(point);
            FaultPoint ca,cb;
            ca = findNaborAbove(pointGrid,point);
            cb = findNaborBelow(pointGrid,ca);
            if (ca!=null && ca.curve==null && cb==point) {
              linkAboveBelow(ca,cb);
              growQueue.add(ca);
            }
            cb = findNaborBelow(pointGrid,point);
            ca = findNaborAbove(pointGrid,cb);
            if (cb!=null && cb.curve==null && ca==point) {
              linkAboveBelow(ca,cb);
              growQueue.add(cb);
            }
          }
        }

        // Done growing. Add this curve to the list of curves. Here we include
        // curves that are too small. If we did not include them here, we would
        // need to put them in a list of small curves, so that we could later
        // remove all of their points. (By not removing those points now, we
        // prevent them from becoming parts of other curves.) Instead, we
        // simply put all curves in the list, and filter that list later.
        curveList.add(curve);
      }
    }

    // Filter curves to include only those that are big enough. Remove all
    // points from any curves that are too small.
    ArrayList<FaultCurve> bigCurveList = new ArrayList<FaultCurve>();
    for (FaultCurve curve:curveList) {
      if (curve.size()>=_ncsmin) {
        bigCurveList.add(curve);
      } else {
        for (FaultPoint point:curve) {
          point.curve = null;
          point.ca = null;
          point.cb = null;
        }
      }
    }
    return bigCurveList.toArray(new FaultCurve[0]);
  }


    // Methods to find good nabors of a specified point. These methods return
  // null if no nabor is good enough, based on various thresholds.
  private FaultPoint findNaborAbove(FaultPointGrid points, FaultPoint point) {
    FaultPoint pa = points.findPointAbove(point);
    return canBeNabors(point,pa)?pa:null;
  }
  private FaultPoint findNaborBelow(FaultPointGrid points, FaultPoint point) {
    FaultPoint cb = points.findPointBelow(point);
    return canBeNabors(point,cb)?cb:null;
  }

  // Returns true if two specified points can be nabors. The two points are
  // assumed to be within one sample of each other. This method uses other
  // attributes of the points to determine whether or not they can be nabors.
  private boolean canBeNabors(FaultPoint pa, FaultPoint pb) {
    boolean can = true;
    if (pa==null || pb==null) {
      can = false;
    } else if (minFl(pa,pb)<_fllo) {
      can = false;
    } else if (minS1(pa,pb)<_fs1min) {
      can = false;
    } else if (maxS1(pa,pb)>_fs1max) {
      can = false;
    } else if (absDeltaFl(pa,pb)>_dflmax) {
      can = false;
    } else if (absDeltaFt(pa,pb)>_dftmax) {
      can = false;
    } else if (absDeltaS1(pa,pb)>_ds1max) {
      can = false;
    }
    return can;
  }
  private static float minFl(FaultPoint pa, FaultPoint pb) {
    return min(pa.fl,pb.fl);
  }
  private static float minS1(FaultPoint ca, FaultPoint cb) {
    return min(ca.s1,cb.s1);
  }
  private static float maxS1(FaultPoint ca, FaultPoint cb) {
    return max(ca.s1,cb.s1);
  }
  private static float absDeltaFl(FaultPoint ca, FaultPoint cb) {
    return abs(ca.fl-cb.fl);
  }
  private static float absDeltaFt(FaultPoint ca, FaultPoint cb) {
    return abs(ca.ft-cb.ft);
  }
  private static float absDeltaS1(FaultPoint ca, FaultPoint cb) {
    return abs(ca.s1-cb.s1);
  }


    // Methods to link mutually best nabors.
  private void linkAboveBelow(FaultPoint ca, FaultPoint cb) {
    cb.ca = ca;
    ca.cb = cb;
  }




}


