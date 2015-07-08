package rgi;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import java.util.*;

/**
 * Interpolation in space-RGT domain.
 * @author Xinming Wu
 * @version 2015.06.13
 */
public class RgtInterp3 {

  public RgtInterp3() {
  }

  public RgtInterp3(float pnull, float[][][] p) {
    float[][] fxs = getPoints(pnull,p);
    _fx = copy(fxs[0]);
    _x1 = copy(fxs[1]);
    _x2 = copy(fxs[2]);
    _x3 = copy(fxs[3]);
  }
  /**
   * Constructs an interpolator.
   * @param fx know values at the known points.
   * @param x1 1st coordinates of known points.
   * @param x2 2nd coordinates of known points.
   * @param x3 3rd coordinates of known points.
   */
  public RgtInterp3(
    float[] fx, float[] x1, float[] x2, float[] x3) 
  {
    _fx = copy(fx);
    _x1 = copy(x1);
    _x2 = copy(x2);
    _x3 = copy(x3);
  }

  public void setRgt(float[][][] u1) {
    _u1 = u1;
  }

  public void setScales(float sv, float sh) {
    _sv = sv;
    _sh = sh;
  }

  public float[][] getPoints(float pnull, float[][][] p) {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    ArrayList<Float> fx = new ArrayList<Float>();
    ArrayList<Float> x1 = new ArrayList<Float>();
    ArrayList<Float> x2 = new ArrayList<Float>();
    ArrayList<Float> x3 = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float pi = p[i3][i2][i1];
      if (pi!=pnull) {
        fx.add(pi);
        x1.add((float)i1);
        x2.add((float)i2);
        x3.add((float)i3);
      }
    }}}
    int np = fx.size();
    float[][] fxs = new float[4][np];
    for (int ip=0; ip<np; ++ip) {
      fxs[0][ip] = fx.get(ip);
      fxs[1][ip] = x1.get(ip);
      fxs[2][ip] = x2.get(ip);
      fxs[3][ip] = x3.get(ip);
    }
    return fxs;
  }

  /**
   * Computes gridded values using nearest neighbors.
   * Gridded values in the array p are computed for only unknown 
   * samples with value equal to the specified null value. Any
   * known (non-null) sample values in the array p are not changed.
   * <p>
   * This method also computes and returns an array of times to
   * nearest-neighbor samples. Times are zero for known samples 
   * and are positive for gridded samples.
   * @param pnull the null value representing unknown samples.
   * @param p array of sample values to be gridded.
   * @return array of times to nearest known samples.
   */
  public float[][][] gridNearest(float pnull, float[][][] p) {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    float[][][] ds = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    System.out.println("sh="+1.0f/_sh);
    System.out.println("sv="+1.0f/_sv);
    ClosestPointTransform cpt = 
      new ClosestPointTransform(1.0f/_sv,1.0f/_sh,1.0f/_sh);
    cpt.apply(pnull,p,ds,k1,k2,k3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      int p1 = k1[i3][i2][i1];
      int p2 = k2[i3][i2][i1];
      int p3 = k3[i3][i2][i1];
      p[i3][i2][i1] = p[p3][p2][p1];
    }}}
    return ds;
  }

  public float[][][] gridNearestK(float pnull, float[][][] p) {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    float[][][] ds = new float[n3][n2][n1];
    ArrayList<Float> fx = new ArrayList<Float>();
    ArrayList<Float> x1 = new ArrayList<Float>();
    ArrayList<Float> x2 = new ArrayList<Float>();
    ArrayList<Float> x3 = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float pi = p[i3][i2][i1];
      if(pi!=pnull) {
        fx.add(pi);
        x1.add((float)i1);
        x2.add((float)i2);
        x3.add((float)i3);
      }
    }}}
    int np = fx.size();
    float[] f = new float[np];
    float[][] x = new float[3][np];
    for (int ip=0; ip<np; ++ip) {
      f[ip] = fx.get(ip);
      x[0][ip] = x1.get(ip);
      x[1][ip] = x2.get(ip);
      x[2][ip] = x3.get(ip);
    }
    float[] scales = new float[3];
    scales[0] = 100.0f;
    scales[1] = 0.001f;
    scales[2] = 0.001f;
    KdTree kt = new KdTree(x);
    kt.setScales(scales);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float pi = p[i3][i2][i1];
      if(pi!=pnull){continue;}
      float[] xi = new float[]{i1,i2,i3};
      int k = kt.findNearest(xi);
      p[i3][i2][i1] = f[k];
      float d1 = (x[0][k]-i1)*scales[0];
      float d2 = (x[1][k]-i2)*scales[1];
      float d3 = (x[2][k]-i3)*scales[2];
      ds[i3][i2][i1] = d1*d1+d2*d2+d3*d3;
    }}}
    return ds;
  }


  public float[][][][] grid(Sampling s1, Sampling s2, Sampling s3) {
    Check.argument(s1.isUniform(),"s1 is uniform");
    Check.argument(s2.isUniform(),"s2 is uniform");
    Check.argument(s3.isUniform(),"s2 is uniform");
    Check.state(_fx!=null,"scattered samples have been set");
    Check.state(_x1!=null,"scattered samples have been set"); 
    Check.state(_x2!=null,"scattered samples have been set");
    Check.state(_x3!=null,"scattered samples have been set");
    convertX1(s1,s2,s3);
    int n3 = s3.getCount();
    int n2 = s2.getCount();
    double fu1 = min(_u1);
    double du1 = s1.getDelta()*0.5;
    double lu1 = max(_u1);
    int nu1 = (int)((lu1-fu1)/du1)+1;
    Sampling su1 = new Sampling(nu1,du1,fu1);
    float pnull = -FLT_MAX;
    float tnull = -FLT_MAX;
    SimpleGridder3 sg = new SimpleGridder3(_fx,_x1,_x2,_x3);
    sg.setNullValue(pnull);
    float[][][] p = sg.grid(su1,s2,s3);
    checkPoints(pnull,p);
    float[][][] t = new float[n3][n2][nu1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<nu1; ++i1) {
        t[i3][i2][i1] = (p[i3][i2][i1]!=pnull)?0.0f:tnull;
      }
    }}
    t = gridNearest(pnull,p);
    t = mul(sqrt(t),2f);
    float[][][] q = p;
    if (_blending) {
      q = new float[n3][n2][nu1];
      gridBlended(t,p,q);
    }
    float[][][] ti = unflatten(su1,t);
    float[][][] pi = unflatten(su1,p);
    float[][][] qi = unflatten(su1,q);
    return new float[][][][]{ti,pi,qi};
  }

    /**
   * Computes gridded values using blended neighbors. 
   * Note that blended-neighbor gridding can be performed only 
   * after nearest-neighbor gridding. Blending does not change
   * the values of known samples for which times are zero.
   * @param t array of times to nearest known samples.
   * @param p array of nearest-neighbor gridded values.
   * @param q array of blended-neighbor gridded values.
   */
  public void gridBlended(float[][][] t, float[][][] p, float[][][] q) {
    int n1 = t[0][0].length;
    int n2 = t[0].length;
    int n3 = t.length;

    // Compute time squared. If necessary, shift to account for the shift 
    // in the finite-difference stencil used in the local diffusion kernel.
    float[][][] s = mul(t,t);
    if (_ldk.getStencil()!=LocalDiffusionKernel.Stencil.D21) {
      for (int i3=n3-1; i3>0; --i3) {
        for (int i2=n2-1; i2>0; --i2) {
          for (int i1=n1-1; i1>0; --i1) {
            s[i3][i2][i1] = 0.125f*(s[i3  ][i2  ][i1  ] +
                                    s[i3  ][i2  ][i1-1] +
                                    s[i3  ][i2-1][i1  ] +
                                    s[i3  ][i2-1][i1-1] +
                                    s[i3-1][i2  ][i1  ] +
                                    s[i3-1][i2  ][i1-1] +
                                    s[i3-1][i2-1][i1  ] +
                                    s[i3-1][i2-1][i1-1]);
          }
        }
      }
    }

    // Construct and apply a local smoothing filter.
    LocalSmoothingFilter lsf = new LocalSmoothingFilter(0.01,10000,_ldk);
    lsf.setPreconditioner(true);
    float pavg = sum(p)/n1/n2/n3;
    float[][][] r = sub(p,pavg);
    // Smoothing attenuates finite-difference errors near Nyquist.
    // The problem with this smoothing is that it makes q != p when
    // known samples are adjacent, as when interpolating well logs.
    // This smoothing should be unnecessary for Stencil.D21.
    if (_ldk.getStencil()!=LocalDiffusionKernel.Stencil.D21)
      lsf.applySmoothS(r,r);
    setTensors(n1,n2,n3);
    lsf.apply(_tensors,_c,s,r,q);
    add(q,pavg,q);

    // Restore the known sample values. Due to errors in finite-difference
    // approximations, these values may have changed during smoothing.
    // Even with time adjustments, this restoration is still necessary
    // if we used applySmoothS above. Best to just do this in any case.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (t[i3][i2][i1]==0.0f) {
            q[i3][i2][i1] = p[i3][i2][i1];
          }
        }
      }
    }
  }

  private float _sv = 1.0f; // scale for distance in vertical direction.
  private float _sh = 1.0f; // scale for distance in horizontal direction.
  private float[] _fx = null; // known values at scattered points.
  private float[] _x1 = null; // 1st coordinates of scattered points.
  private float[] _x2 = null; // 2nd coordinates of scattered points.
  private float[] _x3 = null; // 3rd coordinates of scattered points.
  private float[][][] _u1 = null; // array of RGT values.

  private boolean _blending = true;
  private Tensors3 _tensors;
  private float _c = 0.7f;
  private LocalDiffusionKernel _ldk =
    new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D22);

  private void convertX1(Sampling s1, Sampling s2, Sampling s3) {
    int np = _x1.length;
    SincInterpolator si = new SincInterpolator();
    for (int ip=0; ip<np; ++ip) {
      _x1[ip]=si.interpolate(s1,s2,s3,_u1,_x1[ip],_x2[ip],_x3[ip]);
    }
  }

  private void checkPoints(float pnull, float[][][] p) {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    float[] x1s = new float[n1];
    float[] fxs = new float[n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k = 0;
      int i1e=0;
      int i1b=n1;
      for (int i1=0; i1<n1; ++i1) {
        float pi = p[i3][i2][i1];
        if (pi!=pnull) {
          x1s[k] = i1;
          fxs[k] = pi;
          if (i1e<i1) {i1e=i1;}
          if (i1b>i1) {i1b=i1;}
          k++;
        }
      }
      if(k<=1) {continue;}
      float[][] x = new float[1][k];
      copy(k,x1s,x[0]);
      KdTree kt = new KdTree(x);
      for (int i1=i1b; i1<=i1e; ++i1) {
        int xk = kt.findNearest(new float[]{i1});
        p[i3][i2][i1] = fxs[xk];
      }
    }}
  }

  private float[][][] unflatten(Sampling su1, float[][][] f) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    int nu1 = su1.getCount();
    double du1 = su1.getDelta();
    double fu1 = su1.getFirst();
    SincInterpolator si = new SincInterpolator();
    float[][][] g = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      si.interpolate(nu1,du1,fu1,f[i3][i2],n1,_u1[i3][i2],g[i3][i2]);
    }}
    return g;
  }

  private void setTensors(int n1, int n2, int n3) {
    float[][][] u1 = fillfloat(1.0f,n1,n2,n3);
    float[][][] u2 = fillfloat(0.0f,n1,n2,n3);
    float[][][] w1 = fillfloat(0.0f,n1,n2,n3);
    float[][][] w2 = fillfloat(0.0f,n1,n2,n3);
    float[][][] av = fillfloat(1.0f,n1,n2,n3);
    float[][][] aw = fillfloat(1.0f,n1,n2,n3);
    float[][][] au = fillfloat(_sv/_sh,n1,n2,n3);
    _tensors = new EigenTensors3(u1,u2,w1,w2,au,av,aw,true);
  }
}
