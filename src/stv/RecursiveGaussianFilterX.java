package stv;

import static edu.mines.jtk.util.ArrayMath.*;

public class RecursiveGaussianFilterX {

  public RecursiveGaussianFilterX(float sigma) {
    makeBs(sigma);
  }

  public float[] gauss(int n, float sigma, float mu) {
    float[] y = new float[n];
    float pi = (float)Math.PI;
    for (int i=0; i<n; ++i) {
      y[i] = 1.0f/(sigma*sqrt(2*pi))*exp(-(i-mu)*(i-mu)/(2*sigma*sigma));
    }
    return y;
  }

  public void apply00(float[][] x, float[][] y) {
    apply0X(x,y);
    applyX0(y,y);
  }

  public void applyX0(float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[] x1 = new float[n2];
    float[] y1 = new float[n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2)
        x1[i2] = x[i2][i1];
      apply0(x1,y1);
      for (int i2=0; i2<n2; ++i2)
        y[i2][i1] = y1[i2];
    }
  }

  public void apply0X(float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      apply0(x[i2],y[i2]);
  }

  public void apply0(float[] x, float[] y) {
    int n1 = x.length;
    double pd = _sc*x[0];
    applyForward(pd,x,y);
    double up = x[n1-1]*_sc;
    double vp = up*_sc;
    applyBackward(up,vp,y,y);
  }

  public void applyT(float[] x, float[] y) {
    applyForwardX(x,y);
    applyBackwardX(y,y);
  }


  public void applyForward(double pd, float[] x, float[] y) {
    int n = x.length;
    double yim1 = pd;
    double yim2 = pd;
    double yim3 = pd;
    for (int i=0; i<n; ++i) {
      double xi = x[i];
      double yi = xi+_b1*yim1+_b2*yim2+_b3*yim3;
      y[i] = (float)yi;
      yim3 = yim2;
      yim2 = yim1;
      yim1 = yi;
    }
  }

  public void applyBackward(double up, double vp, float[] x, float[] y) {
    int n = x.length;
    double u0 = x[n-1]-up;
    double u1 = x[n-2]-up;
    double u2 = x[n-3]-up;
    double a0 = _a0*_a0;
    double mu0 = _m[0][0]*u0+_m[0][1]*u1+_m[0][2]*u2;
    double mu1 = _m[1][0]*u0+_m[1][1]*u1+_m[1][2]*u2;
    double mu2 = _m[2][0]*u0+_m[2][1]*u1+_m[2][2]*u2;
    double yip1 = (mu0+vp)*a0;
    double yip2 = (mu1+vp)*a0;
    double yip3 = (mu2+vp)*a0;
    y[n-1] = (float)yip1;
    double b1 = _b1;
    double b2 = _b2;
    double b3 = _b3;
    for (int i=n-2; i>=0; --i) {
      double xi = x[i];
      double yi = a0*xi+b1*yip1+b2*yip2+b3*yip3;
      y[i] = (float)yi;
      yip3 = yip2;
      yip2 = yip1;
      yip1 = yi;
    }
  }

  public void applyForwardX(float[] x, float[] y) {
    int n = x.length;
    double yim1 = 0;
    double yim2 = 0;
    double yim3 = 0;
    for (int i=0; i<n; ++i) {
      double xi = x[i];
      double yi = _a0*xi+_b1*yim1+_b2*yim2+_b3*yim3;
      y[i] = (float)yi;
      yim3 = yim2;
      yim2 = yim1;
      yim1 = yi;
    }
  }


  public void applyBackwardX(float[] x, float[] y) {
    int n = x.length;
    double yip1 = 0.0;
    double yip2 = 0.0;
    double yip3 = 0.0;
    for (int i=n-1; i>=0; --i) {
      double xi = x[i];
      double yi = _a0*xi+_b1*yip1+_b2*yip2+_b3*yip3;
      y[i] = (float)yi;
      yip3 = yip2;
      yip2 = yip1;
      yip1 = yi;
    }
  }

  private void makeBs(float sigma) {
    double q;
    //sigma = ceil(sigma);
    if (sigma<=2.5) {
      q = 3.97156-4.14554*sqrt(1.0-0.26891*sigma);
    } else {
      q = 0.98711*sigma-0.96330;
    }
    double q2 = q*q;
    double q3 = q2*q;
    double t0  = 1.57825+2.444134*q+1.4281*q2+0.422205*q3; 
    double t1 = 2.44413*q+2.85619*q2+1.26661*q3;
    double t2 = -1.4281*q2-1.26661*q3;
    double t3 = 0.422205*q3; 
    double b1 = t1/t0;
    double b2 = t2/t0;
    double b3 = t3/t0;
    _b1 = b1;
    _b2 = b2;
    _b3 = b3;
    _a0 = 1.0-(b1+b2+b3);
    _m[0][0] = -b3*b1+1.0-b3*b3-b2;
    _m[0][1] = (b3+b1)*(b2+b3*b1);
    _m[0][2] = b3*(b1+b3*b2);
    _m[1][0] = b1+b3*b2;
    _m[1][1] = -(b2-1.0)*(b2+b3*b1);
    _m[1][2] = -(b3*b1+b3*b3+b2-1.0)*b3;
    _m[2][0] = b3*b1+b2+b1*b1-b2*b2;
    _m[2][1] = b1*b2+b3*b2*b2-b1*b3*b3-b3*b3*b3-b3*b2+b3;
    _m[2][2] = b3*(b1+b3*b2);
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      _m[i][j] *= (1.0+b2+(b1-b3)*b3);
      _m[i][j] /= (1.0+b1-b2+b3)*(1.0-b1-b2-b3);
    }
    _sc = 1.0/(1.0-_b1-_b2-_b3);
  }

  private double _a0,_b1,_b2,_b3,_sc;
  private double[][] _m = new double[3][3];
}
