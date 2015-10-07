package fit;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Fit a 1D curve with smoothing cubic splines.
 * <p>
 * Reference:
 * <ul><li>
 * D.S.G. Pollock, 1993, Smoothing with cubic splines
 * http://www.physics.muni.cz/~jancely/NM/Texty/Numerika/CubicSmoothingSpline.pdf
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.10.28
 */

public class SmoothSplines {

  public SmoothSplines(double rho, double[] w, double[] x,double[] y) {
    _x = x;
    _y = y;
    _rho = rho;
    int np = w.length;
    _w = filldouble(1.0,np);
    if (w!=null) {_w = w;}
    _c = new double[np][4];
    initSplines();
  }

  public SmoothSplines(double rho, double[] x, double[] y)  {
     this(rho, null, x, y);
  }

  public double interpolate0 (double x) {
    int i = index(x);
    double dx = x-_x[i];
    return interpolate0(dx,_c[i]);
  }

  private int index(double x) {
    int index = binarySearch(_x,x,_index);
    if (index<0) 
      index = (index<-1)?-2-index:0;
    _index = index;
    return index;
  }

  private static double interpolate0(double dx, double[] c) {
    int nc = c.length;
    double ci = c[nc-1];
    for (int i=nc-2; i>=0; i--)
      ci = c[i] + dx*ci;
    return ci;
  }


  private static double interpolate1(double dx, double[] c) {
    int nc = c.length;
    double ci = computeC(nc-1,c);
    for (int k=nc-2; k>=nc; k--)
      ci = computeC(k,c)+dx*ci;
    return ci;
  }

  private static double computeC(int i, double[] c) {
   int n = c.length;
   double ci = c[i];
   for (int k=i; k>i-n; k--)
     ci *= k;
   return ci;
  }

  private void initSplines () {
    int n = _x.length;
    double[] h = new double[n];
    double[] r = new double[n];
    double[] f = new double[n];
    double[] u = new double[n];
    double[] v = new double[n];
    double[] w = new double[n];
    double[] s = new double[n];
    double[] q = new double[n+1];
    for (int i=0; i<n; ++i) {
      double wi = _w[i];
      s[i] = (wi<=0.0)?1.0e200:1.0/sqrt(wi);
    }
    double mu=(_rho<=0)?1.0e200:2.0*(1.0-_rho)/(3.0*_rho);
    h[0] = _x[1] - _x[0];
    r[0] = 3.0/h[0];

    for (int i=1; i<n-1; i++) {
      h[i] = _x[i+1]-_x[i];
      r[i] = 3.0/h[i];
      q[i] = 3.0*(_y[i+1]-_y[i])/h[i]-3.0*(_y[i]-_y[i-1])/h[i-1];
    }

    for (int i=1;i<n-1; i++) {
      double ri =  r[i  ];
      double rm =  r[i-1];
      double rp =  r[i+1];
      double si =  s[i  ];
      double sm =  s[i-1];
      double sp =  s[i+1];
      double xm = _x[i-1];
      double xp = _x[i+1];
      double fi = -rm-ri;
      f[i] = fi;
      u[i] = rm*rm*sm+fi*fi*si+ri*ri*sp;
      u[i] = mu*u[i]+2.0*(xp-xm);
      v[i] = fi*ri*si-ri*(ri+rp)*sp;
      v[i] = mu*v[i]+h[i];
      w[i] = mu*ri*rp*sp;
    }
    q = quincunx(u, v, w, q);

    // spline parameters
    double d3 = 1.0/3.0;
    double r0 = r[0], s0 = s[0], h0 = h[0];
    double q1 = q[1], f1 = f[1], r1 = r[1];
    double q2 = q[2], y1 =_y[1], y0 =_y[0];
    double c00 = y0-mu*r0*q1*s0;
    double c10 = y1-mu*(f1*q1+r1*q2)*s0;
    double c01 = (c10-c00)/h0-q1*h0*d3;
    double c02 = 0.0;
    double c03 = d3*q1/h0;
    _c[0][0] = c00; _c[0][1] = c01;
    _c[0][2] = c02; _c[0][3] = c03;
    _c[1][0] = c10;  r[0]  = 0;
    for (int i=1; i<n-1; ++i) {
      double qi = q[i  ];
      double ri = r[i  ];
      double hi = h[i  ];
      double fi = f[i  ];
      double si = s[i  ];
      double qm = q[i-1];
      double qp = q[i+1];
      double rm = r[i-1];
      double hm = h[i-1];
      double yi = _y[i ];
      double cm1 = _c[i-1][1];
      _c[i][3] = d3*(qp-qi)/hi;
      _c[i][2] = qi;
      _c[i][1] = (qi+qm)*hm+cm1;
      _c[i][0] = yi-mu*(rm*qm+fi*qi+ri*qm)*si;
    }

    _c[n-1][3] = 0;
    _c[n-1][2] = 0;
    _c[n-1][1] = q[n-2]*h[n-2]+_c[n-2][1];
    _c[n-1][0] = _y[n-1]-mu*(r[n-2]*q[n-2])*s[n-1];

   }


   private static double[] quincunx(
     double[] u, double[] v, double[] w, double[] q) 
   {
     int n = u.length;
     u[0] = 0;
     v[1] /= u[1];
     w[1] /= u[1];
     for (int i=2;i<n;i++) {
       u[i] -= u[i-2]*w[i-2]*w[i-2]+u[i-1]*v[i-1]*v[i-1];
       v[i] -= u[i-1]*v[i-1]*w[i-1];
       v[i] /= u[i];
       w[i] /= u[i];
     }

     q[1] = q[1]-v[0]*q[0];
     for (int i=2;i<n;i++)
       q[i] -= v[i-1]*q[i-1]+w[i-2]*q[i-2];
     for (int i=1;i<n; i++)
       q[i] /= u[i];

     q[n-1]  = 0.0;
     q[n-2] -= v[n-2]*q[n-1];
     for (int i=n-3; i>0; i--)
       q[i] -= v[i]*q[i+1]+w[i]*q[i+2];
     return q;
   }

   private double[] _x;
   private double[] _y;
   private double[] _w=null;
   private double _rho=0.50;
   private double[][] _c=null;
   private int _index;
}
