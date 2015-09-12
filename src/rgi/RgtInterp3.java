package rgi;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;
import java.util.*;

/**
 * Interpolation in space-RGT domain.
 * @author Xinming Wu
 * @version 2015.06.13
 */

public class RgtInterp3 {

  public RgtInterp3() {
  }

  public RgtInterp3(Point[][] ps) {
    _ps = ps;
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

  /*
  public float[][][] gridXX(int m1,
    Sampling s2, Sampling s3, float[][][] x1s, float[][][] g) 
  {
    int n3 = s3.getCount();
    int n2 = s2.getCount();
    int n1 = g[0][0].length;
    Sampling s1 = new Sampling(n1);
    float[][] fxs = getSamples(_ps);
    _fx = fxs[0];
    _x1 = fxs[1];
    _x2 = fxs[2];
    _x3 = fxs[3];
    float pnull = -FLT_MAX;
    float tnull = -FLT_MAX;
    SimpleGridder3 sg = new SimpleGridder3(_fx,_x1,_x2,_x3);
    sg.setNullValue(pnull);
    float[][][] p1 = sg.grid(s1,s2,s3);
    float[][][] p2 = copy(p1);
    float[][][] t = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        t[i3][i2][i1] = (p1[i3][i2][i1]!=pnull)?0.0f:tnull;
      }
    }}
    t = gridNearest(pnull,p1);
    t = mul(sqrt(t),2f);
    float[][][] s = coherence(g);
    float[][] ti1  = new float[n3][n2];
    float[][] gi1  = new float[n3][n2];
    float[][] si1  = new float[n3][n2];
    float[][] p1i1 = new float[n3][n2];
    float[][] p2i1 = new float[n3][n2];
    for (int i1=0; i1<n1; ++i1) {
      System.out.println("i1="+i1);
      if(slice(i1,n2,n3,pnull,p2,p2i1)) {
        slice(i1,n2,n3,g,gi1);
        slice(i1,n2,n3,s,si1);
        float[][] q = gridSlice(pnull,p2i1,gi1,si1);
        setValue(i1,n2,n3,p2,q);
      } else {
        slice(i1,n2,n3,s,si1);
        slice(i1,n2,n3,g,gi1);
        slice(i1,n2,n3,t,ti1);
        slice(i1,n2,n3,p1,p1i1);
        float[][] q = gridSlice(ti1,p1i1,gi1,si1);
        setValue(i1,n2,n3,p2,q);
      }
    }
    float[][][] q = unflatten(m1,x1s,p2);
    return q;
  }
  */

  private float[][][] unflatten(int m1, float[][][] x1s, float[][][] g) {
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    float[][][] f = new float[n3][n2][m1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float x1m = x1s[i3][i2][0];
      float yxm =   g[i3][i2][0];
      FloatList xl = new FloatList();
      FloatList yl = new FloatList();
      xl.add(x1m);
      yl.add(yxm);
      for (int i1=1; i1<n1; ++i1) {
        float x1i = x1s[i3][i2][i1];
        float yxi =   g[i3][i2][i1];
        if(x1i>x1m) {
          xl.add(x1i);
          yl.add(yxi);
          x1m = x1i;
          yxm = yxi;
        }
      }
      float[] x = xl.trim();
      float[] y = yl.trim();
      CubicInterpolator ci = new CubicInterpolator(x,y);
      for (int i1=0; i1<m1; ++i1) {
        f[i3][i2][i1] = ci.interpolate(i1);
      }
    }}
    return f;
  }


  public float[][][][] gridX(
    Sampling s1, Sampling s2, Sampling s3, float[][][] f) 
  {
    Check.argument(s1.isUniform(),"s1 is uniform");
    Check.argument(s2.isUniform(),"s2 is uniform");
    Check.argument(s3.isUniform(),"s2 is uniform");
    //convertX1(s1,s2,s3);
    int n3 = s3.getCount();
    int n2 = s2.getCount();
    double fu1 = 0.0;
    double lu1 = max(_u1);
    double du1 = s1.getDelta();
    int nu1 = (int)((lu1-fu1)/du1)+1;
    Sampling su1 = new Sampling(nu1,du1,fu1);
    float[][] fxs = getSamples(_ps);
    _fx = fxs[0];
    _x1 = fxs[1];
    _x2 = fxs[2];
    _x3 = fxs[3];
    float pnull = -FLT_MAX;
    float tnull = -FLT_MAX;
    SimpleGridder3 sg = new SimpleGridder3(_fx,_x1,_x2,_x3);
    sg.setNullValue(pnull);
    float[][][] p1 = sg.grid(su1,s2,s3);
    float[][][] p2 = copy(p1);
    float[][][] t = new float[n3][n2][nu1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<nu1; ++i1) {
        t[i3][i2][i1] = (p1[i3][i2][i1]!=pnull)?0.0f:tnull;
      }
    }}
    t = gridNearest(pnull,p1);
    t = mul(sqrt(t),2f);
    float[][][] g = flatten(s1,su1,f);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float g1 = g[i3][i2][110];
      float g2 = g[i3][i2][111];
      float g3 = g[i3][i2][109];
      float ga = (g1+g2+g3)/3.f;
      g[i3][i2][111]=ga;
      g[i3][i2][110]=ga;
      g[i3][i2][109]=ga;
      g[i3][i2][108]=ga;
      g[i3][i2][107]=ga;
      g[i3][i2][106]=ga;
      g[i3][i2][105]=ga;
    }}

    //float[][][] s = coherence(g);
    float[][] ti1  = new float[n3][n2];
    float[][] gi1  = new float[n3][n2];
    float[][] p1i1 = new float[n3][n2];
    float[][] p2i1 = new float[n3][n2];
    //float[][] si1  = new float[n3][n2];
    float[][] si1 = fillfloat(1.0f,n2,n3);
    for (int i1=0; i1<nu1; ++i1) {
      System.out.println("i1="+i1);
      if(slice(i1,n2,n3,pnull,p2,p2i1)) {
        slice(i1,n2,n3,g,gi1);
        //slice(i1,n2,n3,s,si1);
        //float[][] q = gridSlice(pnull,p2i1,gi1,si1);
        float[][] q = gridSlice(pnull,p2i1,gi1);
        setValue(i1,n2,n3,p2,q);
      } else {
        //slice(i1,n2,n3,s,si1);
        slice(i1,n2,n3,g,gi1);
        slice(i1,n2,n3,t,ti1);
        slice(i1,n2,n3,p1,p1i1);
        float[][] q = gridSlice(ti1,p1i1,gi1,si1);
        setValue(i1,n2,n3,p2,q);
      }
    }
    float[][][] q = unflatten(su1,p2);
    return new float[][][][]{p2,q};
  }

  public void smoothOnFaults(FaultSkin[] skins, float[][][] g) {
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
    float[][][] s = fillfloat(1.0f,n1,n2,n3);
    smooth2(4.0f,s,g);
    smooth3(4.0f,s,g);
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
    double du1 = s1.getDelta()*1.0;
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
  private Point[][] _ps=null;

  private boolean _blending = true;
  private Tensors3 _tensors;
  private float _c = 0.7f;
  private LocalDiffusionKernel _ldk =
    new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D22);

  private float[][] getSamples(Point[][] ps) {
    int ns = ps.length;
    FloatList fvl = new FloatList();
    FloatList x1l = new FloatList();
    FloatList x2l = new FloatList();
    FloatList x3l = new FloatList();
    for (int is=0; is<ns; ++is) {
    int np = ps[is].length;
    for (int ip=0; ip<np; ++ip) {
      Point pi = ps[is][ip];
      fvl.add(pi.fv);
      x1l.add(pi.u1);
      x2l.add(pi.u2);
      x3l.add(pi.u3);
    }}
    float[] fv = fvl.trim();
    float[] x1 = x1l.trim();
    float[] x2 = x2l.trim();
    float[] x3 = x3l.trim();
    return new float[][]{fv,x1,x2,x3};
  }
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

  public Point[][] checkPoints(Sampling s1, Point[][] ps) {
    int ns = ps.length;
    int n1 = s1.getCount();
    double f1 = s1.getFirst();
    double d1 = s1.getDelta();
    Point[][] pp = new Point[ns][];
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
      float[][] x = new float[1][np];
      for (int ip=0; ip<np; ++ip) {
        x[0][ip] = ps[is][ip].u1; 
      }
      float xb = min(x);
      float xe = max(x);
      KdTree kt = new KdTree(x);
      ArrayList<Point> pa = new ArrayList<Point>();
      for (int i1=0; i1<n1; ++i1) {
        float u1 = (float)(f1+d1*i1);
        int xk = kt.findNearest(new float[]{u1});
        float xi = x[0][xk];
        Point pi = ps[is][xk];
        float di = abs(round(xi-u1));
        if((xi==xb||xi==xe||pi.onFault)&&di>0.0f){continue;}
        Point pk = new Point(pi.fv,pi.x1,pi.x2,pi.x3,pi.onFault);
        pk.setCoordW(pi.w1,pi.w2,pi.w3);
        pk.setCoordU(u1,   pi.u2,pi.u3);
        pa.add(pk);
      }
      int ik = 0;
      int npp = pa.size();
      pp[is] = new Point[npp];
      for (Point pi:pa) {
        pp[is][ik] = pi;
        ik++;
      }
    }
    return pp;
  }


  private float[][][] unflatten(Sampling su1, float[][][] f) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    int nu1 = su1.getCount();
    double du1 = su1.getDelta();
    double fu1 = su1.getFirst();
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    float[][][] g = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      si.interpolate(nu1,du1,fu1,f[i3][i2],n1,_u1[i3][i2],g[i3][i2]);
    }}
    return g;
  }

  private float[][][] unflattenX(Sampling su1, float[][][] f) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    int nu1 = su1.getCount();
    double du1 = su1.getDelta();
    double fu1 = su1.getFirst();
    float[][][] g = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      double u1i = (double)_u1[i3][i2][i1];
      int k1 = (int)round((u1i-fu1)/du1);
      if(k1>=nu1){k1=nu1-1;}
      g[i3][i2][i1] = f[i3][i2][k1];
    }}}
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

  private float[][][] flatten(
    Sampling sx1, Sampling su1, final float[][][] f) 
  {
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int nx1 = sx1.getCount();
    final int nu1 = su1.getCount();
    final float[][][] x1 = new float[n3][n2][nu1];
    final InverseInterpolator ii = new InverseInterpolator(sx1,su1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(_u1[i3][i2],x1[i3][i2]);
    }});

    final double d1 = sx1.getDelta();
    final double f1 = sx1.getFirst();
    final float[][][] g = new float[n3][n2][nu1];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        si.interpolate(nx1,d1,f1,f[i3][i2],nu1,x1[i3][i2],g[i3][i2]);
    }});
    return g;
  }


  public float[][][] coherence(float[][][] f) {
    int sig1 = 2;
    int sig2 = sig1*6;
    float sigma1 = 2.0f;
    float sigma2 = 4.0f;
    float sigma3 = 4.0f;
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2,sigma3);
    EigenTensors3 et = lof.applyForTensors(f);
    LocalSemblanceFilter lsf = new LocalSemblanceFilter(sig1,sig2);
    float[][][] s = lsf.semblance(LocalSemblanceFilter.Direction3.V,et,f);
    return pow(s,10f);
  }


  private float[][] gridSlice(
    float pnull, float[][] p, float[][] f) 
  {
    EigenTensors2 et = makeImageTensors(f);
    BlendedGridder2 bg = new BlendedGridder2(et);
    float[][] t = bg.gridNearest(pnull,p);
    float[][] q = copy(p);
    bg.gridBlended(t,p,q);
    return q;
  }

  public EigenTensors2 makeImageTensors(float[][]f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float sigmad  = 3.0f;
    float sigmas = 8.0f;
    LocalOrientFilter lofs= new LocalOrientFilter(sigmas);
    EigenTensors2 ets = lofs.applyForTensors(f);
    float[][] ed = edge(f);
    float[][] sc = coherence(ets,ed);
    LocalOrientFilter lofd = new LocalOrientFilter(sigmad);
    EigenTensors2 et  = lofd.applyForTensors(ed);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,n2*4,sc,sc);
    sc = pow(sc,10);
    sc = sub(sc,min(sc));
    sc = div(sc,max(sc));
    float[][] eu = new float[n2][n1];
    float[][] ev = new float[n2][n1];
    et.getEigenvalues(eu,ev);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float sci = sc[i2][i1];
      if(sci<0.15f) {
        eu[i2][i1] *=0.0001f;
        ev[i2][i1] *=0.0001f;
      } else {
        eu[i2][i1] *=0.002f;
        ev[i2][i1] *=0.0001f;
      }
    }}
    et.setEigenvalues(eu,ev);
    et.invertStructure(1.0,1.0);
    return et;
  }

  public float[][] edge(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] u1 = new float[n2][n1]; 
    float[][] u2 = new float[n2][n1]; 
    float[][] g1 = new float[n2][n1]; 
    float[][] g2 = new float[n2][n1]; 
    LocalOrientFilter lof = new LocalOrientFilter(8.0);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(3.0);
    lof.applyForNormal(f,u1,u2);
    rgf.apply10(f,g1);
    rgf.apply01(f,g2);
    float[][] gu = add(mul(u1,g1),mul(u2,g2));
    gu = abs(gu);
    gu = sub(gu,min(gu));
    gu = div(gu,max(gu));
    return gu;
  }

  public float[][] coherence(EigenTensors2 et, float[][] f) {
    int sig2 = 2;
    int sig1 = 32;
    LocalSemblanceFilter lsf = new LocalSemblanceFilter(sig1,sig2);
    return lsf.semblance(LocalSemblanceFilter.Direction2.V,et,f);
  }


  private float[][] gridSlice(
    float pnull, float[][] p, float[][] f, float[][] s) 
  {
    EigenTensors2 et = makeImageTensors(f,s);
    BlendedGridder2 bg = new BlendedGridder2(et);
    float[][] t = bg.gridNearest(pnull,p);
    float[][] q = copy(p);
    bg.gridBlended(t,p,q);
    return q;
  }

  private EigenTensors2 makeImageTensors(float[][]f, float[][] s) {
    float sigma1 = 4.0f;
    float sigma2 = 4.0f;
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    EigenTensors2 et = lof.applyForTensors(f);
    s = sub(1.0f,s);
    s = clip(0.0001f,1.0f,s);
    et.scale(s);
    et.invertStructure(1.0,1.0);
    return et;
  }


  private float[][] gridSlice(
    float[][] t, float[][] p, float[][] f, float[][] s) 
  {
    EigenTensors2 et = makeImageTensors(f,s);
    BlendedGridder2 bg = new BlendedGridder2(et);
    float[][] q = copy(p);
    bg.gridBlended(t,p,q);
    return q;
  }

  private boolean slice(
    int i1, int n2, int n3, float vnull, float[][][] f, float[][] f1) {
    boolean valid = false;
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
      float fi = f[i3][i2][i1];
      if(fi!=vnull){valid=true;}
      f1[i3][i2] = fi;
    }}
    return valid;
  }

  private void slice(
    int i1, int n2, int n3, float[][][] f, float[][] f1) {
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
      f1[i3][i2] = f[i3][i2][i1];
    }}
  }

  private void setValue(
    int i1, int n2, int n3, float[][][] f, float[][] f1) {
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
      f[i3][i2][i1]=f1[i3][i2];
    }}
  }

  private class FloatList {
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

    // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] st = fillfloat(1.0f,n2);
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      if (s!=null) {
        for (int i2=0; i2<n2; ++i2)
          st[i2] = s[i2][i1];
      }
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }
  private static void smooth2(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] s3 = (s!=null)?s[i3]:null;
      float[][] x3 = x[i3];
      smooth2(sigma,s3,x3);
    }});
  }

  // Smoothing for dimension 3.
  private static void smooth3(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] s2 = (s!=null)?new float[n3][]:null;
      float[][] x2 = new float[n3][];
      for (int i3=0; i3<n3; ++i3) {
        if (s!=null)
          s2[i3] = s[i3][i2];
        x2[i3] = x[i3][i2];
      }
      smooth2(sigma,s2,x2);
    }});
  }




}
