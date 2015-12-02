package stv;

import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class SteerableFilter2 {

  public SteerableFilter2 (int nd, float sigma) {
    _nd = nd;
    _sigma = sigma;
    _se = new int[nd];
    for (int id=0; id<nd; ++id) 
      _se[id] = 2*id;
  }

  public float[][] applyFilterX(float[][] fx, float[][] u1, float[][] u2) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int k1 = max((int)(_sigma*4),20);
    int k2 = max((int)(_sigma*4),20);
    //int k1 = round(_sigma);
    //int k2 = round(_sigma);
    float hpi = (float)(Math.PI*0.5);
    float[][] fk = new float[n2][n1];
    float[][] fsr = new float[n2][n1];
    float[][] fsi = new float[n2][n1];
    float[][] gsr = new float[n2][n1];
    float[][] gsi = new float[n2][n1];
    float[][][][] bs = filterBasis(k1,k2,_sigma);
    for (int id=0; id<_nd; ++id) {
      int ms = _se[id];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float phi = hpi+atan2(u1i,u2i);
        float scr =  cos(phi*ms);
        float sci = -sin(phi*ms);
        fsr[i2][i1] = scr*fx[i2][i1];
        fsi[i2][i1] = sci*fx[i2][i1];
      }}
      conv(n1,n2,0,0,fsr,k1,k2,0,0,bs[0][id],n1,n2,k1/2,k2/2,gsr);
      conv(n1,n2,0,0,fsi,k1,k2,0,0,bs[1][id],n1,n2,k1/2,k2/2,gsi);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float gri = gsr[i2][i1];
        float gii = gsi[i2][i1];
        fk[i2][i1] += gri-gii;
      }}
    }
    return fk;
  }


  public float[][] applyFilter(float[][] fx, float[][] u1, float[][] u2) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int k1 = max((int)(_sigma+2),20);
    int k2 = max((int)(_sigma+2),20);
    float hpi = (float)(Math.PI*0.5);
    float[][] fk = new float[n2][n1];
    float[][] fsr = new float[n2][n1];
    float[][] fsi = new float[n2][n1];
    float[][][][] bs = filterBasis(k1,k2,_sigma);
    for (int id=0; id<_nd; ++id) {
      int ms = _se[id];
      conv(n1,n2,0,0,fx,k1,k2,0,0,bs[0][id],n1,n2,k1/2,k2/2,fsr);
      conv(n1,n2,0,0,fx,k1,k2,0,0,bs[1][id],n1,n2,k1/2,k2/2,fsi);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float phi = hpi+atan2(u1i,u2i);
        float scr =  cos(phi*ms);
        float sci = -sin(phi*ms);
        float bsr = fsr[i2][i1];
        float bsi = fsi[i2][i1];
        fk[i2][i1] += bsr*scr-bsi*sci;
      }}
    }
    return fk;
  }


  private float[][][][] filterBasis(int n1, int n2, float sig) {
    float c1 = ceil(n1/2f);
    float c2 = ceil(n2/2f);
    float sg = -0.5f/(sig*sig);
    float[] xc = new float[2];
    float[][][][] bs = new float[2][_nd][n2][n1];
    double[] ss = new double[_nd];
    int nd = 50;
    int nf = nd/2;
    for (int is=0; is<_nd; ++is) {
      ss[is] = (factorial(nd-nf+1,nd,1)/factorial(2,nf,1));
      if (is!=0) {ss[is]*=2.0;}
      nf ++;
    }
    div(ss,sum(ss)*sig,ss);
    for (int is=0; is<_nd; ++is) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float x1i = i1-c1;
      float x2i = i2-c2;
      float xsi = x1i*x1i+x2i*x2i;
      double exi = exp(xsi*sg)*ss[is];
      if (is==0) {
        bs[is][0][i2][i1] = (float)exi;
      } else {
        if(xsi==0.0f){
          bs[0][is][i2][i1] = (float)(ss[is]);
          bs[1][is][i2][i1] = 0.0f;
        } else {
          xc[0] = x2i/sqrt(xsi);
          xc[1] = x1i/sqrt(xsi);
          xc = cpow(xc,_se[is]);
          bs[0][is][i2][i1] = (float)(exi*xc[0]);
          bs[1][is][i2][i1] = (float)(exi*xc[1]);
        }
      }
    }}}
    return bs;
  }

  private double factorial(int nb, int ne, double fs) {
    for (int k=nb; k<=ne; ++k) 
      fs *= k;
    return fs;
  }

  private int _nd;
  private int[] _se;
  private float _sigma = 10f;

}
