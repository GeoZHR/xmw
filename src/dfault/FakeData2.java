/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dfault;

import java.util.Random;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Generates 2D fake training data for DeepFault.
 * <em>
 * Jacobians of functions used in folding and faulting have been implemented
 * but not tested. Therefore, beware of errors in calculated slopes p2 and p3.
 * </em>
 * @author Xinming Wu, University of Texas at Austin.
 * @version 2017.08.01
 */
public class FakeData2 {

  public FakeData2(float fd, float dd, float ld) {
    float pi = (float)Math.PI;
    float dt = 1; //interval of fault throws
    float da = 1; //interval of vertical shifts
    float db = pi/18f; //interval of phase shifts
    float dc = 0.2f; //interval of phase scales
    float ds = 0.02f; //interval of vertical shears
    float dn = 0.8f; //interval of noise

    float ft = -30.f;  //first fault throws
    float fa = 0.00f;    //first vertical shifts
    float fb = 0.00f;   //first phase shifts
    float fc = -5.0f;  //first phase scales
    float fs = -0.5f;  //first vertical shears
    float fn = 0.00f;    //first noise

    float lt = 30.0f;   //last fault throws
    float la = 20.0f;   //last vertical shifts
    float lb = pi*2f; //last phase shifts
    float lc = 5.00f;  //last phase scales
    float ls = 0.50f;  //last vertical shear
    float ln = 0.80f;  //last noise

    _tx=count(ft,dt,lt); //number of fault throws
    _ax=count(fa,da,la); //number of vertical shifts
    _bx=count(fb,db,lb); //number of phase shifts
    _cx=count(fc,dc,lc); //number of phase scales
    _sx=count(fs,ds,ls); //number of vertical shears
    _nx=count(fn,dn,ln); //number of noise
    _hx=new float[]{-12,-11,-10,-9,-8,-7,-6,-5,0,0,0,1,0,0,1,0,0,0,0,5,6,7,8,9,10,11,12};
    _nh = _hx.length;
    _nt = _tx.length;
    _na = _ax.length;
    _nb = _bx.length;
    _nc = _cx.length;
    _ns = _sx.length;
    _nn = _nx.length;
    _tx[0] = 0;
    _tx[_nt-1] = 0;
    _tx[6] = 0;
    _tx[_nt-5] = 0;
    int nd = 0;
    for (float x= fd; x<=ld; x+=dd)
      nd++;
    nd *=2;
    _nd = nd;
    _dx = new float[_nd];
    for (int id=0; id<nd/2; ++id) {
      _dx[id] = -ld+dd*id;
    }
    for (int id=nd/2; id<nd; ++id) {
      _dx[id] = fd+dd*(id-nd/2);
    }
  }

  /**
   * Returns a fake 2D seismic image with slopes.
   * The fake seismic image contains sinusoidal folding, horizontal and
   * dipping layers, two unconformities, and two intersecting faults with
   * throws that increase linearly with depth. While the image may have
   * a specified amount of additive noise, the slopes are noise-free.
   * @param noise rms of noise (relative to signal) added to the image.
   */
  public static float[][][] seismicAndSlopes2d2014A(double noise) {
    int n1 = 300;
    int n2 = 200;
    float[][][] p = makeReflectivityWithNormals(n1,n2);
    float pi = (float)Math.PI;
    Linear1 throw1 = new Linear1(10f,0f);
    LinearFault2 fault1 = new LinearFault2(n1*0.5f,n2*0.5f, 15.0f,throw1);
    Sinusoidal2 fold = new Sinusoidal2(20.0f,pi,0*pi/n2);
    VerticalShear2 shear = new VerticalShear2(new Linear1(0.0f,-1.0f));
    p = apply(fold,p);
    p = apply(shear,p);
    p = apply(fault1,p);
    p = addWavelet(0.1,p);
    p[0] = addNoise(noise,p[0]);
    p[1] = neg(div(p[2],p[1]));
    return new float[][][]{p[0],p[1]};
  }

  public float[][] getTrainDataAndLabels(float[][][] sx) {
    final int n3 = sx.length;
    final int n2 = sx[0].length;
    final int n1 = sx[0][0].length;
    final int m1 = 300;
    final int m2 = 200;
    final float[][] label = new float[n3][_nd];
    //final float[] label = new float[n3];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      int f1 = round((m1-n1)*0.5f);
      int f2 = round((m2-n2)*0.5f);
      Random r = new Random();
      int ia = r.nextInt(_na);
      int ib = r.nextInt(_nb);
      int ic = r.nextInt(_nc);
      int id = r.nextInt(_nd);
      int is = r.nextInt(_ns);
      int it = r.nextInt(_nt);
      int ih = r.nextInt(_nh);
      int in = r.nextInt(_nn);
      float axi = _ax[ia];
      float bxi = _bx[ib];
      float cxi = _cx[ic];
      float sxi = _sx[is];
      float txi = _tx[it];
      float dxi = _dx[id];
      float hxi = _hx[ih];
      float nxi = 0.0f;//_nx[in];
      float pi = (float)Math.PI;
      f2 += hxi;
      float[][][] p = makeReflectivityWithNormals(m1,m2);
      Linear1 throw1 = new Linear1(txi,0f);
      LinearFault2 fault1 = new LinearFault2(m1*0.5f,m2*0.5f,dxi,throw1);
      Sinusoidal2 fold = new Sinusoidal2(axi,bxi,cxi*pi/m2);
      VerticalShear2 shear = new VerticalShear2(new Linear1(0.0f,sxi));
      p = apply(fold,p);
      p = apply(shear,p);
      p = apply(fault1,p);
      p = addWavelet(0.1,p);
      //float[][] u1 = new float[m2][m1];
      //float[][] u2 = new float[m2][m1];
      //float[][] el = new float[m2][m1];
      //LocalOrientFilter lof = new LocalOrientFilter(2,1);
      //lof.applyForNormalLinear(p[0],u1,u2,el);
      //p[0] = addNoise(nxi,p[0]);
      sx[i3] = copy(n1,n2,f1,f2,p[0]);
      //ex[i3] = copy(n1,n2,f1,f2,el);
      if(abs(txi)>0f&&abs(hxi)<=1f) {
        label[i3][id] = 1;
        //label[i3] = 1;
      }
      /*
      for (float i2=0; i2<n2; i2+=0.1f) {
        float theta = pi*(0.5f-dxi/180f);
        float pxi = tan(theta);
        int i1 = round(pxi*i2+(n1-n2*pxi)*0.5f);
        int k2 = round(i2);
        if(i1>=0&&i1<n1) 
          data[i3][k2][i1] = 1f;
      }
      */
    }});
    return label;
  }

  public float[][][] dipToFaultImage(int n1, int n2, float[][] ps) {
    int np = ps.length;
    float pi = (float)Math.PI;
    float[][][] gs = new float[np][n2][n1];
    for (int ip=0; ip<np; ip+=1) {
      int[] id = new int[1];
      float pm = max(ps[ip],id);
      float di = _dx[id[0]];
      System.out.println("pm="+pm);
      if(pm==0f) continue;
      for (float i2=0; i2<n2; i2+=0.1f) {
        float theta = pi*(0.5f-di/180f);
        float pxi = tan(theta);
        int i1 = round(pxi*i2+(n1-n2*pxi)*0.5f);
        int k2 = round(i2);
        if(i1>=0&&i1<n1) 
          gs[ip][k2][i1] = pm;
      }
    }
    return gs;
  }

  private int _na,_nb,_nc,_ns,_nn,_nt,_nd,_nh;
  private float[] _ax,_bx,_cx,_sx,_nx,_tx,_dx,_hx;


  private float[] count(float f, float d, float l) {
    int n = 0;
    for (float x=f; x<=l; x+=d)
      n++;
    float[] y = new float[n];
    int k = 0;
    for (float x=f; x<=l; x+=d)
      y[k++] = x;
    return y;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static SincInterpolator _si = new SincInterpolator();
  static {
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

  private static float[][] addNoise(double nrms, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    Random r = new Random(1);
    float[][] g = mul(2.0f,sub(randfloat(r,n1,n2),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    rgf.apply10(g,g); // 1st derivative enhances high-frequencies
    g = mul(g,(float)nrms*rms(f)/rms(g));
    return add(f,g);
  }

  /**
   * Coordinates, Jacobians, and transforms.
   */
  private static class C2 {
    float c1,c2;
    C2(float c1, float c2) { 
      this.c1 = c1; 
      this.c2 = c2; 
    }
  }

  private static class D2 {
    float d11,d12,
          d21,d22;
    D2(float d11, float d12,
       float d21, float d22) { 
      this.d11 = d11;  this.d12 = d12; 
      this.d21 = d21;  this.d22 = d22; 
    }
  }

  private interface T1 {
    float f(float x);
    float df(float x);
  }
  private interface T2 {
    C2 f(float x1, float x2);
    D2 df(float x1, float x2);
  }
  /**
   * Applies a 2D coordinate transform to an image.
   * @param t coordinate transform f(x).
   * @param p input image p(x).
   * @return transformed image q(x) = p(f(x)).
   */
  private static float[][] apply(T2 t, float[][] p) {
    int n1 = p[0].length;
    int n2 = p.length;
    float[][] q = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        C2 f = t.f(i1,i2);
        float f1 = f.c1;
        float f2 = f.c2;
        q[i2][i1] = _si.interpolate(n1,1.0,0.0,n2,1.0,0.0,p,f1,f2);
      }
    }
    return q;
  }

  /**
   * Applies a 2D coordinate transform to an image and normal vectors.
   * The input array {p0,p1,p2} contains the input image p0 and 1st and 2nd
   * components of normal vectors, p1 and p2. The returned array {q0,q1,q2}
   * contains the corresponding transformed image and normal vectors. All
   * normal vectors are unit vectors.
   * @param t coordinate transform f(x).
   * @param p input image p(x) and normal vectors.
   * @return transformed image q(x) = p(f(x) and normal vectors.
   */
  private static float[][][] apply(T2 t, float[][][] p) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] q0 = apply(t,p[0]);
    float[][] q1 = apply(t,p[1]);
    float[][] q2 = apply(t,p[2]);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        D2 d = t.df(i1,i2);
        float q1i = d.d11*q1[i2][i1]+d.d21*q2[i2][i1];
        float q2i = d.d12*q1[i2][i1]+d.d22*q2[i2][i1];
        float qsi = 1.0f/sqrt(q1i*q1i+q2i*q2i);
        q1[i2][i1] = q1i*qsi;
        q2[i2][i1] = q2i*qsi;
      }
    }
    return new float[][][]{q0,q1,q2};
  }

  /**
   * A linear 1D coordinate mapping f(x) = a0+a1*x.
   */
  private static class Linear1 implements T1 {
    public Linear1(float a0, float a1) { 
      _a0 = a0; _a1 = a1;
    }
    public float f(float x) { 
      return _a0+_a1*x; 
    }
    public float df(float x) { 
      return _a1; 
    }
    private float _a0,_a1;
  }

  /**
   * A linear fault in a 2D image.
   */
  private static class LinearFault2 implements T2 {

    /**
     * Constructs a linear fault.
     * @param fx1 coordinate x1 of a reference point on the fault.
     * @param fx2 coordinate x2 of a reference point on the fault.
     * @param ftheta fault dip, measured in degrees from vertical.
     * @param fthrow fault throw, a function of coordinate x1.
     */
    public LinearFault2(float fx1, float fx2, float ftheta, T1 fthrow) {

      // Reference point (the origin in fault-line coordinates).
      _r1 = fx1;
      _r2 = fx2;

      // Tangent of fault dip.
      float rtheta = toRadians(ftheta);
      _ttheta = tan(rtheta);

      // Fault normal vector.
      float ctheta = cos(rtheta);
      float stheta = sin(rtheta);
      _u1 = -stheta;
      _u2 =  ctheta;

      // Ensure vertical component of normal vector is non-positive.
      if (_u1>0.0f) {
        _u1 = -_u1;
        _u2 = -_u2;
      }

      // Constant needed to locate points with respect to plane.
      _u0 = -(fx1*_u1+fx2*_u2);

      // Fault throw.
      _t1 = fthrow;
    }
    public C2 f(float x1, float x2) {
      if (faulted(x1,x2)) {
        x1 -= _r1;
        x2 -= _r2;
        float t = _t1.f(x1);
        x1 -= t;
        x2 -= t*_ttheta;
        x1 += _r1;
        x2 += _r2;
      }
      return new C2(x1,x2);
    }
    public D2 df(float x1, float x2) {
      float d11 = 1.0f, d12 = 0.0f,
            d21 = 0.0f, d22 = 1.0f;
      if (faulted(x1,x2)) {
        x1 -= _r1;
        x2 -= _r2;
        float dt = _t1.df(x1);
        d11 -= dt;
        d21 -= dt*_ttheta;
      }
      return new D2(d11,d12,
                    d21,d22);
    }
    private float _r1,_r2;
    private float _ttheta;
    private float _u0,_u1,_u2;
    private T1 _t1;
    private boolean faulted(float x1, float x2) {
      return _u0+_u1*x1+_u2*x2>=0.0f;
    }
  }


  /**
   * Vertical shear of a 2D image.
   */
  private static class VerticalShear2 implements T2 {
    public VerticalShear2(T1 s1) {
      _s1 = s1;
    }
    public C2 f(float x1, float x2) {
      x1 -= _s1.f(x2);
      return new C2(x1,x2);
    }
    public D2 df(float x1, float x2) {
      float d12 = -_s1.df(x2);
      return new D2(1.0f,  d12,
                    0.0f, 1.0f);
    }
    private T1 _s1;
  }

  /**
   * Sinusoidal folding in a 2D image.
   */
  private static class Sinusoidal2 implements T2 {
    public Sinusoidal2(float a, float b, float c) {
      _a = a;
      _b = b;
      _c = c;
    }
    public C2 f(float x1, float x2) {
      float f1 = x1-_a*sin(_b+_c*x2);
      float f2 = x2;
      return new C2(f1,f2);
    }
    public D2 df(float x1, float x2) {
      float d11 = 1.0f;
      float d12 = -_a*_c*cos(_b+_c*x2);

      float d21 = 0.0f;
      float d22 = 1.0f;
      return new D2(d11,d12,d21,d22);
    }
    private float _a,_b,_c;
  }

  private static float[][][] makeReflectivityWithNormals(int n1, int n2) {
    Random random = new Random(31);
    float[] r = pow(mul(2.0f,sub(randfloat(random,n1),0.5f)),5.0f);
    float[][][] p = new float[3][n2][n1];
    for (int i2=0; i2<n2; ++i2)
      copy(r,p[0][i2]);
    p[1] = fillfloat(1.0f,n1,n2);
    p[2] = fillfloat(0.0f,n1,n2);
    return p;
  }

  private static float[][][] addWavelet(double fpeak, float[][][] p) {
    double sigma = max(1.0,1.0/(2.0*PI*fpeak));
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] p0 = p[0];
    float[][] p1 = p[1];
    float[][] p2 = p[2];
    float[][] q = copy(p0);
    float[][] q1 = new float[n2][n1];
    float[][] q2 = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    for (int id=0; id<2; ++id) { // 2nd directional derivative of Gaussian
      rgf.apply10(q,q1);
      rgf.apply01(q,q2);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          q[i2][i1] = p1[i2][i1]*q1[i2][i1]+p2[i2][i1]*q2[i2][i1];
        }
      }
    }
    q = mul(q,-1.0f/rms(q)); // negate for Ricker wavelet
    return new float[][][]{q,p1,p2};
  }

  private static float rms(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    double sum = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fi = f[i2][i1];
        sum += fi*fi;
      }
    }
    return (float)sqrt(sum/n1/n2);
  }

}
