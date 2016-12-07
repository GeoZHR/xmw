package med;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Transform an image from xyz to polar coordinates
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.10.08
 */

public class PolarCoordinates {

  public PolarCoordinates(float c1, float c2, float rmax) {
    _c1 = c1;
    _c2 = c2;
    _rmax = rmax;
    _st = angleSampling(rmax,0,360);
    _sp = angleSampling(rmax,0,360);
  }

  public PolarCoordinates(float c1, float c2, float c3, float rmin, float rmax) {
    _c1 = c1;
    _c2 = c2;
    _c3 = c3;
    _rmax = rmax;
    _st = angleSampling(rmax,0,360);
    _sp = angleSampling(rmax,0,360);
    int nr = round(rmax-rmin);
    _sr = new Sampling(nr,1.0,rmin);
  }


  public float[][] forwardTransform(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    int nr = round(_rmax);
    int np = _sp.getCount();
    float[][] fr = new float[np][nr];
    SincInterpolator si = new SincInterpolator();
    for (int ir=0; ir<nr; ++ir) {
    for (int ip=0; ip<np; ++ip) {
      float phi = (float)toRadians(_sp.getValue(ip));
      float rx1 = ir*sin(phi);
      float rx2 = ir*cos(phi);
      float x1i = _c1-rx1;
      float x2i = _c2-rx2;
      fr[ip][ir] = si.interpolate(s1,s2,fx,x1i,x2i);
    }}
    return fr;
  }

  public float[][][] forwardTransform(final float[][][] fx) {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final int nr = _sr.getCount();
    final int np = _sp.getCount();
    final int nt = _st.getCount();
    final float[][][] fr = new float[nt][np][nr];
    SincInterpolator si = new SincInterpolator();
    loop(nt,new LoopInt() {
    public void compute(int it) {
      float theta = (float)toRadians(_sp.getValue(it));
      float sti = sin(theta);
      float cti = cos(theta);
      for (int ip=0; ip<np; ++ip) {
      for (int ir=0; ir<nr; ++ir) {
        float ri = (float)_sr.getValue(ir); 
        float phi = (float)toRadians(_sp.getValue(ip));
        float rx1 = ri*cti;
        float rx2 = ri*sti*sin(phi);
        float rx3 = ri*sti*cos(phi);
        float x1i = _c1-rx1;
        float x2i = _c2-rx2;
        float x3i = _c3-rx3;
        fr[it][ip][ir] = si.interpolate(s1,s2,s3,fx,x1i,x2i,x3i);
      }}
    }});
    return fr;
  }

  public float[][] reverseTransform(float[] r) {
    int np = _sp.getCount();
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    for (int ip=0; ip<np; ++ip) {
      float ri = r[ip];
      float phi = (float)toRadians(_sp.getValue(ip));
      float rx1 = ri*sin(phi);
      float rx2 = ri*cos(phi);
      x1[ip] = _c1-rx1;
      x2[ip] = _c2-rx2;
    }
    return new float[][]{x1,x2};
  }

  public float[] reverseTransform(float[][] r) {
    float[] tpr = buildTrigs(_st,_sp,r);
    int nc = tpr.length;
    float[] xyz = new float[nc];
    float dr = (float)_sr.getDelta();
    float fr = (float)_sr.getFirst();
    for (int ic=0; ic<nc; ic+=3) {
      float it = tpr[ic  ];
      float ip = tpr[ic+1];
      float ir = tpr[ic+2]*dr+fr;
      float phi = (float)toRadians(ip);
      float theta = (float)toRadians(it);
      float sti = sin(theta);
      float cti = cos(theta);
      float rx1 = ir*cti;
      float rx2 = ir*sti*sin(phi);
      float rx3 = ir*sti*cos(phi);
      xyz[ic  ] = _c3-rx3;
      xyz[ic+1] = _c2-rx2;
      xyz[ic+2] = _c1-rx1;
    }
    return xyz;
  }


  public float[] reverseTransform(int n1, int n2, float[] r) {
    float[] x1 = fillfloat(n1-1,n2);
    int np = _sp.getCount();
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    float ri = r[0];
    float phi = (float)toRadians(_sp.getValue(0));
    float rx1 = ri*sin(phi);
    float rx2 = ri*cos(phi);
    float x1i = _c1-rx1;
    float x2i = _c2-rx2;
    float x2t = x2i;
    x1a.add(x1i);
    x2a.add(x2i);
    for (int ip=1; ip<np; ++ip) {
      ri = r[ip];
      phi = (float)toRadians(_sp.getValue(ip));
      rx1 = ri*sin(phi);
      rx2 = ri*cos(phi);
      x1i = _c1-rx1;
      x2i = _c2-rx2;
      if (x2i>x2t) {
        x2t = x2i;
        x1a.add(x1i);
        x2a.add(x2i);
      }
    }
    int mp = x1a.size(); 
    float[] x1p = new float[mp];
    float[] x2p = new float[mp];
    for (int ip=0; ip<mp; ++ip) {
      x1p[ip] = x1a.get(ip);
      x2p[ip] = x2a.get(ip);
    }
    int b2 = round(min(x2p));
    int e2 = round(max(x2p));
    b2 = max(b2,0);
    e2 = min(e2,n2-1);
    CubicInterpolator ci = new CubicInterpolator(x2p,x1p);
    for (int i2=b2; i2<=e2; ++i2) {
      x1[i2] = ci.interpolate(i2);
    }
    return x1;
  }

  private float[] buildTrigs(Sampling sx, Sampling sy, float[][] z){
    int i = 0;
    int nx = z.length;
    int ny = z[0].length;
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=0;ix<nx-1; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      //if(x0>180f) {continue;}
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny-1; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        if(y0<0f||y0>360f) {continue;}
        float y1 = (float)sy.getValue(iy+1);
        xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
      }
    }
    return copy(i,0,xyz);
  }


  private static Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(2.0/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  private float _c1,_c2,_c3,_rmax;
  private Sampling _st,_sp,_sr;

}
