package scl;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

import util.*;

public class ChannelScanner {

  /**
   * Constructs a channel scanner with specified parameters.
   * @param sigmaMin minimum scale of channels.
   * @param sigmaMax maximum scale of channels.
   */
  public ChannelScanner(double sigmaMin, double sigmaMax) {
    _sigmaMin = sigmaMin;
    _sigmaMax = sigmaMax;
  }

  public float[][] initialOrient(
    float[][] u1, float[][] u2) {
    int n2 = u1.length;
    int n1 = u1[0].length;
    float pi = (float)Math.PI;
    float[][] os = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      os[i2][i1] = 0.5f*pi+atan2(u1i,u2i);
    }}
    return os;
  }


  public float[][] smooth(
    float sigma, float[][] u1, float[][] u2, float[][] cl) 
  {
    int n2 = cl.length;
    int n1 = cl[0].length;
    float[][] cs = new float[n2][n1];
    float[][] eu = fillfloat(0.0001f,n1,n2);
    float[][] ev = fillfloat(1.0000f,n1,n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    EigenTensors2 ets = new EigenTensors2(u1,u2,eu,ev);
    lsf.apply(ets,sigma,cl,cs);
    return cs;
  }


  public float[][][] scan(float gamma, float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] cl = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] u1k = new float[n2][n1];
    float[][] u2k = new float[n2][n1];
    float sigMin = (float)_sigmaMin;
    float sigMax = (float)_sigmaMax;
    for (float sigma=sigMin; sigma<=sigMax; sigma = sigma+1.0f) {
      float[][] clk = scan(sigma,gamma,x,u1k,u2k);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float clt = clk[i2][i1];
        if(clt>cl[i2][i1]) {
          cl[i2][i1] = clt;
          u1[i2][i1] = u1k[i2][i1];
          u2[i2][i1] = u2k[i2][i1];
        }
      }}
    }
    return new float[][][]{u1,u2,cl};
  }

  public float[][][] scan(float gamma, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] cl = new float[n3][n2][n1];
    float sigMin = (float)_sigmaMin;
    float sigMax = (float)_sigmaMax;
    for (float sigma=sigMin; sigma<=sigMax; sigma = sigma+1.0f) {
      float[][][] ct = scan(sigma,gamma,x);
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float cti = ct[i3][i2][i1];
        if(cti>cl[i3][i2][i1]) {
          cl[i3][i2][i1] = cti;
        }
      }}}
    }
    return cl;
  }


    /**
   * Applies this filter for the specified image and outputs. All
   * outputs are optional and are computed for only non-null arrays.
   * @param x input array for 2-D image
   * @param u1 1st component of 1st eigenvector.
   * @param u2 2nd component of 1st eigenvector.
   * @param v1 1st component of 2nd eigenvector.
   * @param v2 2nd component of 2nd eigenvector.
   * @param eu largest eigenvalue corresponding to the eigenvector u.
   * @param ev smallest eigenvalue corresponding to the eigenvector v.
   */
  public float[][] scan(float sigma, float gamma, 
    float[][] x, float[][] u1, float[][] u2) 
  {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] h11 = new float[n2][n1];
    float[][] h12 = new float[n2][n1];
    float[][] h22 = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    rgf.apply2X(x,h11);
    rgf.apply1X(x,h12);
    rgf.applyX1(h12,h12);
    rgf.applyX2(x,h22);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    LocalOrientFilterP lof = new LocalOrientFilterP(4f,4f);
    EigenTensors2 ets = lof.applyForTensors(x);
    ets.setEigenvalues(0.0001f,1.0f);
    lsf.apply(ets,64,h11,h11);
    lsf.apply(ets,64,h12,h12);
    lsf.apply(ets,64,h22,h22);

    mul(h11,sigma*sigma,h11);
    mul(h12,sigma*sigma,h12);
    mul(h22,sigma*sigma,h22);

    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    float[][] cl = new float[n2][n1];
    float beta = 0.5f;
    float gammas = 0.5f/(gamma*gamma);
    float betas = 0.5f/(beta*beta);

    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = h11[i2][i1];
        a[0][1] = h12[i2][i1];
        a[1][0] = h12[i2][i1];
        a[1][1] = h22[i2][i1];
        solveSymmetric22(a,z,e);
        float eui = e[0];
        float evi = e[1];
        float u1i = z[0][0];
        float u2i = z[0][1];
        if (u1i<0.0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        if (eui>=0f) {continue;}
        float b = evi*evi/(eui*eui);
        float s = eui*eui+evi*evi;
        u1[i2][i1] = u1i;
        u2[i2][i1] = u2i;
        cl[i2][i1] = exp(-b*betas)*(1f-exp(-s*gammas));
      }
    }
    return cl;
  }

  public float[][][] scan(float sigma, float gamma, float[][][] x) {
    final int n3 = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    float[][][] g1 = new float[n3][n2][n1];
    final float[][][] h11 = new float[n3][n2][n1];
    final float[][][] h12 = new float[n3][n2][n1];
    final float[][][] h13 = new float[n3][n2][n1];
    final float[][][] h22 = new float[n3][n2][n1];
    final float[][][] h23 = new float[n3][n2][n1];
    final float[][][] h33 = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    rgf.apply2XX(x,h11);
    rgf.applyX2X(x,h22);
    rgf.applyXX2(x,h33);
    rgf.apply1XX(x,g1);
    rgf.applyX1X(g1,h12);
    rgf.applyXX1(g1,h13);
    rgf.applyX1X(x,h23);
    rgf.applyXX1(h23,h23);

    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    LocalOrientFilterP lof = new LocalOrientFilterP(2.0,6.0,6.0);
    EigenTensors3 ets = lof.applyForTensors(x);
    ets.setEigenvalues(0.0001f,0.0001f,1.0f);
    lsf.apply(ets,64,h11,h11);
    lsf.apply(ets,64,h12,h12);
    lsf.apply(ets,64,h13,h13);
    lsf.apply(ets,64,h22,h22);
    lsf.apply(ets,64,h23,h23);
    lsf.apply(ets,64,h33,h33);

    float alpha  = 0.5f;
    final float sigmas = sigma*sigma;
    final float alphas = 0.5f/(alpha*alpha);
    final float gammas = 0.5f/(gamma*gamma);
    final float betas  = alphas;
    final float[][][] cl = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      double[] e = new double[3];
      double[][] a = new double[3][3];
      double[][] z = new double[3][3];
      for (int i2=0; i2<n2; ++i2) {
        float[] h11i = h11[i3][i2];
        float[] h12i = h12[i3][i2];
        float[] h13i = h13[i3][i2];
        float[] h22i = h22[i3][i2];
        float[] h23i = h23[i3][i2];
        float[] h33i = h33[i3][i2];
        float[] cli = cl[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          a[0][0] = h11i[i1]*sigmas;
          a[0][1] = h12i[i1]*sigmas;
          a[0][2] = h13i[i1]*sigmas;
          a[1][0] = h12i[i1]*sigmas;
          a[1][1] = h22i[i1]*sigmas;
          a[1][2] = h23i[i1]*sigmas;
          a[2][0] = h13i[i1]*sigmas;
          a[2][1] = h23i[i1]*sigmas;
          a[2][2] = h33i[i1]*sigmas;
          solveSymmetric33Jacobi(a,z,e);
          float e1 = (float)e[2];
          float e2 = (float)e[1];
          float e3 = (float)e[0];
          float ra = e2*e2/(e3*e3);
          float rb = e1*e1/abs(e2*e3);
          float rc = e1*e1+e2*e2+e3*e3;
          if(e2>=0||e3>=0) {continue;}
          cli[i1] = (1-exp(-ra*alphas))*exp(-rb*betas)*(1-exp(-rc*gammas));
        }
      }
     }});
    return cl;
  }


  /**
   * Computes eigenvalues and eigenvectors for a symmetric 2x2 matrix A.
   * If the eigenvectors are placed in columns in a matrix V, and the 
   * eigenvalues are placed in corresponding columns of a diagonal 
   * matrix D, then AV = VD.
   * @param a the symmetric matrix A.
   * @param v the array of eigenvectors v[0] and v[1].
   * @param d the array of eigenvalues d[0] and d[1].
   */
  public static void solveSymmetric22(float[][] a, float[][] v, float[] d) {

    // Copy matrix to local variables.
    float a00 = a[0][0];
    float a01 = a[0][1],  a11 = a[1][1];

    // Initial eigenvectors. 
    float v00 = 1.0f,     v01 = 0.0f;
    float v10 = 0.0f,     v11 = 1.0f;

    // If off-diagonal element is non-zero, zero it with a Jacobi rotation.
    if (a01!=0.0f) {
      float tiny = 0.1f*sqrt(FLT_EPSILON); // avoid overflow in r*r below
      float c,r,s,t,u,vpr,vqr;
      u = a11-a00;
      if (abs(a01)<tiny*abs(u)) {
        t = a01/u;
      } else {
        r = 0.5f*u/a01;
        t = (r>=0.0f)?1.0f/(r+sqrt(1.0f+r*r)):1.0f/(r-sqrt(1.0f+r*r));
      }
      c = 1.0f/sqrt(1.0f+t*t);
      s = t*c;
      u = s/(1.0f+c);
      r = t*a01;
      a00 -= r;
      a11 += r;
      //a01 = 0.0f;
      vpr = v00;
      vqr = v10;
      v00 = vpr-s*(vqr+vpr*u);
      v10 = vqr+s*(vpr-vqr*u);
      vpr = v01;
      vqr = v11;
      v01 = vpr-s*(vqr+vpr*u);
      v11 = vqr+s*(vpr-vqr*u);
    }

    // Copy eigenvalues and eigenvectors to output arrays.
    d[0] = a00;
    d[1] = a11;
    v[0][0] = v00;  v[0][1] = v01;
    v[1][0] = v10;  v[1][1] = v11;

    // Sort eigenvalues (and eigenvectors) in descending order.
    if (abs(d[0])<abs(d[1])) {
      float dt = d[1];
      d[1] = d[0];
      d[0] = dt;
      float[] vt = v[1];
      v[1] = v[0];
      v[0] = vt;
    }
  }


  private static void solveSymmetric33Jacobi(
    double[][] a, double[][] v, double[] d) 
  {

    // Copy matrix to local variables.
    double a00 = a[0][0],
           a01 = a[0][1],  a11 = a[1][1],
           a02 = a[0][2],  a12 = a[1][2],  a22 = a[2][2];

    // Initial eigenvectors. 
    double v00 = 1.0,  v01 = 0.0,  v02 = 0.0,
           v10 = 0.0,  v11 = 1.0,  v12 = 0.0,
           v20 = 0.0,  v21 = 0.0,  v22 = 1.0;

    // Tiny constant to avoid overflow of r*r (in computation of t) below.
    double tiny = 0.1*sqrt(DBL_EPSILON);
    
    // Absolute values of off-diagonal elements.
    double aa01 = abs(a01);
    double aa02 = abs(a02);
    double aa12 = abs(a12);

    // Apply Jacobi rotations until all off-diagonal elements are zero.
    // Count rotations, just in case this does not converge.
    for (int nrot=0; aa01+aa02+aa12>0.0; ++nrot) {
      Check.state(nrot<100,"number of Jacobi rotations is less than 100");
      double c,r,s,t,u,vpr,vqr,apr,aqr;

      // If a01 is the largest off-diagonal element, ...
      if (aa01>=aa02 && aa01>=aa12) {
        u = a11-a00;
        if (abs(a01)<tiny*abs(u)) {
          t = a01/u;
        } else {
          r = 0.5*u/a01;
          t = (r>=0.0)?1.0/(r+sqrt(1.0+r*r)):1.0/(r-sqrt(1.0+r*r));
        }
        c = 1.0/sqrt(1.0+t*t);
        s = t*c;
        u = s/(1.0+c);
        r = t*a01;
        a00 -= r;
        a11 += r;
        a01 = 0.0;
        apr = a02;
        aqr = a12;
        a02 = apr-s*(aqr+apr*u);
        a12 = aqr+s*(apr-aqr*u);
        vpr = v00;
        vqr = v10;
        v00 = vpr-s*(vqr+vpr*u);
        v10 = vqr+s*(vpr-vqr*u);
        vpr = v01;
        vqr = v11;
        v01 = vpr-s*(vqr+vpr*u);
        v11 = vqr+s*(vpr-vqr*u);
        vpr = v02;
        vqr = v12;
        v02 = vpr-s*(vqr+vpr*u);
        v12 = vqr+s*(vpr-vqr*u);
      } 
      
      // Else if a02 is the largest off-diagonal element, ...
      else if (aa02>=aa01 && aa02>=aa12) {
        u = a22-a00;
        if (abs(a02)<tiny*abs(u)) {
          t = a02/u;
        } else {
          r = 0.5*u/a02;
          t = (r>=0.0)?1.0/(r+sqrt(1.0+r*r)):1.0/(r-sqrt(1.0+r*r));
        }
        c = 1.0/sqrt(1.0+t*t);
        s = t*c;
        u = s/(1.0+c);
        r = t*a02;
        a00 -= r;
        a22 += r;
        a02 = 0.0;
        apr = a01;
        aqr = a12;
        a01 = apr-s*(aqr+apr*u);
        a12 = aqr+s*(apr-aqr*u);
        vpr = v00;
        vqr = v20;
        v00 = vpr-s*(vqr+vpr*u);
        v20 = vqr+s*(vpr-vqr*u);
        vpr = v01;
        vqr = v21;
        v01 = vpr-s*(vqr+vpr*u);
        v21 = vqr+s*(vpr-vqr*u);
        vpr = v02;
        vqr = v22;
        v02 = vpr-s*(vqr+vpr*u);
        v22 = vqr+s*(vpr-vqr*u);
      } 

      // Else if a12 is the largest off-diagonal element, ...
      else {
        u = a22-a11;
        if (abs(a12)<tiny*abs(u)) {
          t = a12/u;
        } else {
          r = 0.5*u/a12;
          t = (r>=0.0)?1.0/(r+sqrt(1.0+r*r)):1.0/(r-sqrt(1.0+r*r));
        }
        c = 1.0/sqrt(1.0+t*t);
        s = t*c;
        u = s/(1.0+c);
        r = t*a12;
        a11 -= r;
        a22 += r;
        a12 = 0.0;
        apr = a01;
        aqr = a02;
        a01 = apr-s*(aqr+apr*u);
        a02 = aqr+s*(apr-aqr*u);
        vpr = v10;
        vqr = v20;
        v10 = vpr-s*(vqr+vpr*u);
        v20 = vqr+s*(vpr-vqr*u);
        vpr = v11;
        vqr = v21;
        v11 = vpr-s*(vqr+vpr*u);
        v21 = vqr+s*(vpr-vqr*u);
        vpr = v12;
        vqr = v22;
        v12 = vpr-s*(vqr+vpr*u);
        v22 = vqr+s*(vpr-vqr*u);
      }

      // Update absolute values of all off-diagonal elements.
      aa01 = abs(a01);
      aa02 = abs(a02);
      aa12 = abs(a12);
    }

    // Copy eigenvalues and eigenvectors to output arrays.
    d[0] = a00;
    d[1] = a11;
    d[2] = a22;
    v[0][0] = v00;  v[0][1] = v01;  v[0][2] = v02;
    v[1][0] = v10;  v[1][1] = v11;  v[1][2] = v12;
    v[2][0] = v20;  v[2][1] = v21;  v[2][2] = v22;

    // Sort eigenvalues (and eigenvectors) in descending order.
    sortDescending33(v,d);
  }


    /**
   * Sorts eigenvalues d and eigenvectors v in descending order.
   */
  private static void sortDescending33(double[][] v, double[] d) {
    for (int i=0; i<3; ++i) {
      for (int j=i; j>0 && abs(d[j-1])<abs(d[j]); --j) {
        double dj = d[j];
        d[j] = d[j-1];
        d[j-1] = dj;
        double[] vj = v[j];
        v[j] = v[j-1];
        v[j-1] = vj;
      }
    }
  }

  private double _sigmaMin, _sigmaMax;
  
}
