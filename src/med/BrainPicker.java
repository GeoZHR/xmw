package med;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Transform an image from xyz to polar coordinates
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.10.08
 */

public class BrainPicker {

  public BrainPicker(
    float x1, float x2, float rx, float minPhi, float maxPhi) 
  {
    _x1 = x1;
    _x2 = x2;
    _rx = rx;
    _minPhi = minPhi;
    _maxPhi = maxPhi;
  }

  public float[][] applyTransform(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling sp = angleSampling(_rx,_minPhi,_maxPhi);
    int nr = round(_rx);
    int np = sp.getCount();
    float[][] fr = new float[np][nr];
    SincInterpolator si = new SincInterpolator();
    for (int ir=0; ir<nr; ++ir) {
    for (int ip=0; ip<np; ++ip) {
      float phi = (float)toRadians(sp.getValue(ip));
      float rx1 = ir*sin(phi);
      float rx2 = ir*cos(phi);
      float x1i = _x1-rx1;
      float x2i = _x2-rx2;
      fr[ip][ir] = si.interpolate(s1,s2,fx,x1i,x2i);
    }}
    return fr;
  }

  public float[][] reverseTransform(float[] r) {
    Sampling sp = angleSampling(_rx,_minPhi,_maxPhi);
    int np = sp.getCount();
    float[] x1 = new float[np];
    float[] x2 = new float[np];
    for (int ip=0; ip<np; ++ip) {
      float ri = r[ip];
      float phi = (float)toRadians(sp.getValue(ip));
      float rx1 = ri*sin(phi);
      float rx2 = ri*cos(phi);
      x1[ip] = _x1-rx1;
      x2[ip] = _x2-rx2;
    }
    return new float[][]{x1,x2};
  }


  public float[] reverseTransform(int n1, int n2, float[] r) {
    float[] x1 = fillfloat(n1-1,n2);
    Sampling sp = angleSampling(_rx,_minPhi,_maxPhi);
    int np = sp.getCount();
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    float ri = r[0];
    float phi = (float)toRadians(sp.getValue(0));
    float rx1 = ri*sin(phi);
    float rx2 = ri*cos(phi);
    float x1i = _x1-rx1;
    float x2i = _x2-rx2;
    float x2t = x2i;
    x1a.add(x1i);
    x2a.add(x2i);
    for (int ip=1; ip<np; ++ip) {
      ri = r[ip];
      phi = (float)toRadians(sp.getValue(ip));
      rx1 = ri*sin(phi);
      rx2 = ri*cos(phi);
      x1i = _x1-rx1;
      x2i = _x2-rx2;
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

  private static Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  private float _x1,_x2,_rx,_minPhi,_maxPhi;

}
