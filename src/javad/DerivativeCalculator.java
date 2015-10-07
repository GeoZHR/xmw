package javad;

import fit.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class DerivativeCalculator {

  public float[][] upSample(float d1, float d2, float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    DoubleList dl1 = new DoubleList();
    DoubleList dl2 = new DoubleList();
    for (float k=0.0f; k<n1; k+=d1)
      dl1.add(k);
    for (float k=0.0f; k<n2; k+=d2)
      dl2.add(k);
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling ss1 = new Sampling(dl1.trim());
    Sampling ss2 = new Sampling(dl2.trim());
    int ns1 = ss1.getCount();
    int ns2 = ss2.getCount();
    float[][] fs = new float[ns2][ns1];
    SincInterpolator si = new SincInterpolator();
    for (int i2=0; i2<ns2; i2++) {
    for (int i1=0; i1<ns1; i1++) {
      double x1 = ss1.getValue(i1);
      double x2 = ss2.getValue(i2);
      fs[i2][i1] = si.interpolate(s1,s2,f,x1,x2);
    }}
    return fs;
  }

  public float[] regrid(double rho, float[] cv) {
    int n = cv.length;
    float[] sv = new float[n];
    DoubleList xc = new DoubleList();
    DoubleList yc = new DoubleList();
    for (int i=0; i<n; i+=10) {
      xc.add(i); 
      yc.add(cv[i]);
    }
    double[] ws = filldouble(1.0,n);
    SmoothSplines ss = new SmoothSplines(rho,ws,xc.trim(),yc.trim());
    for (int i=0; i<n; i++) {
      sv[i] = (float)ss.interpolate0(i);
    }
    return sv;
  }
  public float[] regrid0(float[] cv) {
    int n = cv.length;
    float[] sv = new float[n];
    FloatList xc1 = new FloatList();
    FloatList yc1 = new FloatList();
    FloatList xc2 = new FloatList();
    FloatList yc2 = new FloatList();
    for (int i=0; i<n; i+=20) {
      xc1.add(i); yc1.add(cv[i]);
    }
    for (int i=10; i<n; i+=20) {
      xc2.add(i); yc2.add(cv[i]);
    }
    CubicInterpolator c1 = 
      new CubicInterpolator(CubicInterpolator.Method.SPLINE,xc1.trim(),yc1.trim());
    CubicInterpolator c2 = 
      new CubicInterpolator(CubicInterpolator.Method.SPLINE,xc2.trim(),yc2.trim());
    for (int i=0; i<n; i++) {
      sv[i] = 0.5f*(c1.interpolate(i)+c2.interpolate(i));
    }
    return sv;
  }

  public float[] regrid1(float[] cv) {
    int n = cv.length;
    float[] sv = new float[n];
    float[] vt = regrid0(cv);
    FloatList xc1 = new FloatList();
    FloatList yc1 = new FloatList();
    FloatList xc2 = new FloatList();
    FloatList yc2 = new FloatList();
    for (int i=0; i<n; i+=20) {
      xc1.add(i); yc1.add(vt[i]);
    }
    for (int i=10; i<n; i+=20) {
      xc2.add(i); yc2.add(vt[i]);
    }
    CubicInterpolator c1 = 
      new CubicInterpolator(CubicInterpolator.Method.SPLINE,xc1.trim(),yc1.trim());
    CubicInterpolator c2 = 
      new CubicInterpolator(CubicInterpolator.Method.SPLINE,xc2.trim(),yc2.trim());
    for (int i=0; i<n; i++) {
      sv[i] = 0.5f*(c1.interpolate1(i)+c2.interpolate1(i));
    }
    return sv;
  }


  private class DoubleList {
  public int n;
  public double[] a = new double[1024];
  public void add(double f) {
    if (n==a.length) {
      double[] t = new double[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = f;
  }
  public double[] trim() {
    if (n==0)
      return null;
    double[] t = new double[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
  }

  private class FloatList {
  public int n;
  public float[] a = new float[1024];
  public void add(float f) {
    if (n==a.length) {
      float[] t = new float[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = f;
  }
  public float[] trim() {
    if (n==0)
      return null;
    float[] t = new float[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
  }

}
