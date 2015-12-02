package stv;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class HermiteGaussFilter2 {

  public HermiteGaussFilter2 (int order, float lambda) {
    _nd = order;
    _lambda = lambda;
    float lbs = _lambda*_lambda;
    _beta  = (1f/lbs-lbs)*0.25f;
    _alpha = (1f/lbs+lbs)*0.25f;
  }


  public float[][] applyFilter(float sig, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int nfft1 = n1;//FftComplex.nfftSmall(n1);
    int nfft2 = n2;//FftComplex.nfftSmall(n2);
    /*
    FftComplex fft1 = new FftComplex(nfft1);
    FftComplex fft2 = new FftComplex(nfft2);

    float[][] fc = czerofloat(nfft1,nfft2);
    complexImage(fx,fc);
    fft1.complexToComplex1(-1,nfft2,fc,fc);
    fft2.complexToComplex2(-1,nfft1,fc,fc);
    float[][][] bs = filterBasis(nfft1,nfft2);
    int ns = bs.length;
    for (int is=0; is<ns; ++is) {
      fft1.complexToComplex1(-1,nfft2,bs[is],bs[is]);
      fft2.complexToComplex2(-1,nfft1,bs[is],bs[is]);
    }
    float[][] u2 = czerofloat(nfft1,nfft2);
    for (int is=0; is<ns; ++is){
    for (int i2=0; i2<nfft2; ++i2){
    for (int i1=0,ir=0,ii=1; i1<nfft1; ++i1,ir+=2,ii+=2){
      float fcr = fc[i2][ir];
      float fci = fc[i2][ii];
      float bsr = bs[is][i2][ir];
      float bsi = bs[is][i2][ii];
      u2[i2][ir] = fcr*bsr-fci*bsi;
      u2[i2][ii] = fcr*bsi+fci*bsr;
    }}}
    fft1.complexToComplex1(1,nfft2,u2,u2);
    fft2.complexToComplex2(1,nfft1,u2,u2);
    fft1.scale(nfft1,nfft2,u2);
    fft2.scale(nfft1,nfft2,u2);

    float[][] cu2 = cabs(u2);
    float[][] au2 = zerofloat(nfft1,nfft2);
    int j1 = nfft1/2;
    int j2 = nfft2/2;
    copy(nfft1-j1,nfft2-j2,0,0,cu2,j1,j2,au2);
    copy(j1,j2,nfft1-j1,nfft2-j2,cu2,0,0,au2);
    copy(nfft1-j1,j2,0,nfft2-j2,cu2,j1,0,au2);
    copy(j1,nfft2-j2,nfft1-j1,0,cu2,0,j2,au2);
    */

    float[][][][] bs = filterBasisX(nfft1,nfft2,sig);
    int ns = bs[0].length;
    System.out.println("ns="+ns);
    float[][] fc = new float[n2][n1];
    float phi = -(float)(Math.PI/180);
    phi *=20f;
    float[] sc = new float[2];
    for (int is=0; is<ns; ++is) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float bsr = bs[0][is][i2][i1];
      float bsi = bs[1][is][i2][i1];
      sc[0] = cos(phi);
      sc[1] = sin(phi);
      sc = cpow(sc,is*2);
      float fcr = bsr*sc[0]-bsi*sc[1];
      fc[i2][i1] += fcr;
    }}}
    return fc;
  }

  public float[][] applyFilterX(
    float sig, float[][] fx, float[][] u1, float[][] u2) {
    int k1 = 50;
    int k2 = 50;
    int n2 = fx.length;
    int n1 = fx[0].length;
    float hpi = (float)(Math.PI*0.5);
    float[][] fk = new float[n2][n1];
    float[][] fsr = new float[n2][n1];
    float[][] fsi = new float[n2][n1];
    float[][] gsr = new float[n2][n1];
    float[][] gsi = new float[n2][n1];
    float[][][][] bs = filterBasisX(k1,k2,sig);
    for (int id=0; id<_nd; ++id) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float phi = atan2(u1i,u2i);
        float scr =  cos(phi*id*2);
        float sci = -sin(phi*id*2);
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

  public float[][] applyFilter10(
    final float sig, final float[][] fx, final float[][] u1, final float[][] u2) {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final int nfft1 = FftComplex.nfftSmall(n1);
    final int nfft2 = FftComplex.nfftSmall(n2);
    final float[][] gx = zerofloat(nfft1,nfft2);
    final float[][] fw = czerofloat(nfft1,nfft2);
    final FftComplex fft1 = new FftComplex(nfft1);
    final FftComplex fft2 = new FftComplex(nfft2);
    for (int id=0; id<_nd; ++id) {
      final float[][] fc = czerofloat(nfft1,nfft2);
      final float[][] wc = czerofloat(nfft1,nfft2);
      filterKernal10(id,sig,wc);
      imageFeatures(id,fx,u1,u2,fc);
      fft1.complexToComplex1(-1,nfft2,fc,fc);
      fft2.complexToComplex2(-1,nfft1,fc,fc);
      fft1.complexToComplex1(-1,nfft2,wc,wc);
      fft2.complexToComplex2(-1,nfft1,wc,wc);
      for (int i2=0; i2<nfft2; ++i2){
      for (int i1=0,ir=0,ii=1; i1<nfft1; ++i1,ir+=2,ii+=2) {
        float cri = fc[i2][ir];
        float cii = fc[i2][ii];
        float wri = wc[i2][ir];
        float wii = wc[i2][ii];
        fw[i2][ir] += (cri*wri-cii*wii);
        fw[i2][ii] += (cri*wii+cii*wri);
      }}
    }
    fft1.complexToComplex1(1,nfft2,fw,fw);
    fft2.complexToComplex2(1,nfft1,fw,fw);
    fft1.scale(nfft1,nfft2,fw);
    fft2.scale(nfft1,nfft2,fw);
    int j1 = nfft1/2;
    int j2 = nfft2/2;
    float[][] ga = zerofloat(nfft1,nfft2);
    for (int i2=0; i2<nfft2; ++i2){
    for (int i1=0,ir=0; i1<nfft1; ++i1,ir+=2) {
      ga[i2][i1] = fw[i2][ir];
    }}
    copy(nfft1-j1,nfft2-j2,0,0,ga,j1,j2,gx);
    copy(j1,j2,nfft1-j1,nfft2-j2,ga,0,0,gx);
    copy(nfft1-j1,j2,0,nfft2-j2,ga,j1,0,gx);
    copy(j1,nfft2-j2,nfft1-j1,0,ga,0,j2,gx);
    return copy(n1,n2,gx);
  }


  public float[][] applyFilter00(
    final float sig, final float[][] fx, final float[][] u1, final float[][] u2) {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final int nfft1 = FftComplex.nfftSmall(n1);
    final int nfft2 = FftComplex.nfftSmall(n2);
    final float[][] gx = zerofloat(nfft1,nfft2);
    final float[][] fw = czerofloat(nfft1,nfft2);
    final FftComplex fft1 = new FftComplex(nfft1);
    final FftComplex fft2 = new FftComplex(nfft2);
    for (int id=0; id<_nd; ++id) {
      final float[][] fc = czerofloat(nfft1,nfft2);
      final float[][] wc = czerofloat(nfft1,nfft2);
      filterKernal00(id,sig,wc);
      imageFeatures(id,fx,u1,u2,fc);
      fft1.complexToComplex1(-1,nfft2,fc,fc);
      fft2.complexToComplex2(-1,nfft1,fc,fc);
      fft1.complexToComplex1(-1,nfft2,wc,wc);
      fft2.complexToComplex2(-1,nfft1,wc,wc);
      for (int i2=0; i2<nfft2; ++i2){
      for (int i1=0,ir=0,ii=1; i1<nfft1; ++i1,ir+=2,ii+=2) {
        float cri = fc[i2][ir];
        float cii = fc[i2][ii];
        float wri = wc[i2][ir];
        float wii = wc[i2][ii];
        fw[i2][ir] += (cri*wri-cii*wii);
        fw[i2][ii] += (cri*wii+cii*wri);
      }}
    }
    fft1.complexToComplex1(1,nfft2,fw,fw);
    fft2.complexToComplex2(1,nfft1,fw,fw);
    fft1.scale(nfft1,nfft2,fw);
    fft2.scale(nfft1,nfft2,fw);
    int j1 = nfft1/2;
    int j2 = nfft2/2;
    float[][] ga = zerofloat(nfft1,nfft2);
    for (int i2=0; i2<nfft2; ++i2){
    for (int i1=0,ir=0; i1<nfft1; ++i1,ir+=2) {
      ga[i2][i1] = fw[i2][ir];
    }}
    copy(nfft1-j1,nfft2-j2,0,0,ga,j1,j2,gx);
    copy(j1,j2,nfft1-j1,nfft2-j2,ga,0,0,gx);
    copy(nfft1-j1,j2,0,nfft2-j2,ga,j1,0,gx);
    copy(j1,nfft2-j2,nfft1-j1,0,ga,0,j2,gx);
    return copy(n1,n2,gx);
  }

  private void imageFeatures(
    int s, float[][] fx, float[][] u1, float[][] u2, float[][] fc) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0,ir=0,ii=1; i1<n1; ++i1,ir+=2,ii+=2) {
      float fxi = fx[i2][i1];
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      float phi = atan2(u1i,u2i);
      float scr =  cos(phi*s*2);
      float sci = -sin(phi*s*2);
      fc[i2][ir] = scr*fxi;
      fc[i2][ii] = sci*fxi;
    }}
  }



  public float[][] checkKernal(int n1, int n2, double sig) {
    float[][] fb = new float[n2][n1];
    float[][] fk = czerofloat(n1,n2);
    for (int id=0; id<_nd; ++id) {
      filterKernal10(id,sig,fk);
      for (int i2=0; i2<n2; ++i2){
      for (int i1=0,ir=0,ii=1; i1<n1; ++i1,ir+=2,ii+=2) {
        float fkr = fk[i2][ir];
        float fki = fk[i2][ii];
        fb[i2][i1] += fkr;
      }}
    }
    return fb;
  }
  private void filterKernal00(int s, double sig, float[][] fk) {
    int n2 = fk.length;
    int n1 = fk[0].length/2;
    double c1 = ceil(n1/2f);
    double c2 = ceil(n2/2f);
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0,ir=0,ii=1; i1<n1; ++i1,ir+=2,ii+=2) {
      double x1i = i1-c1;
      double x2i = i2-c2;
      double[] xc = new double[2];
      double xsi = x1i*x1i+x2i*x2i;
      if(xsi==0f) {
        xc[0] = 1f;
        xc[1] = 0f;
      } else {
        xc[0] = x2i/sqrt(xsi);
        xc[1] = x1i/sqrt(xsi);
      }
      xsi /= sig;
      double bx = _beta*xsi;
      double ax = _alpha*xsi;
      if (bx>700||ax>700){continue;}
      double exi = exp(-ax);
      double bsi = bessi(bx,s);
      double ssi = exi*bsi;
      xc = cpow(xc,(s*2));
      fk[i2][ir] = (float)(ssi*xc[0]);
      fk[i2][ii] = (float)(ssi*xc[1]);
    }}
  }

  private void filterKernal10(int s, double sig, float[][] fk) {
    int n2 = fk.length;
    int n1 = fk[0].length/2;
    double c1 = ceil(n1/2f);
    double c2 = ceil(n2/2f);
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0,ir=0,ii=1; i1<n1; ++i1,ir+=2,ii+=2) {
      double x1i = i1-c1;
      double x2i = i2-c2;
      double[] xc = new double[2];
      double xsi = x1i*x1i+x2i*x2i;
      double xsr = sqrt(xsi);
      if(xsi==0f) {
        xc[0] = 0f;
        xc[1] = 0f;
      } else {
        xc[0] = x2i/xsr;
        xc[1] = x1i/xsr;
      }
      xsi /= sig;
      double bx = _beta*xsi;
      double ax = _alpha*xsi;
      if (bx>700||ax>700){continue;}
      double exi = exp(-ax);
      double bsi = bessi(bx,s);
      double bsp = bessi(bx,s+1);
      double ssi = _lambda*xsr*exi*(bsp+bsi);
      xc = cpow(xc,(s*2+1));
      fk[i2][ir] = (float)(ssi*xc[0]);
      fk[i2][ii] = (float)(ssi*xc[1]);
    }}
  }


  private void filterKernal20(int s, double sig, float[][] fk) {
    int n2 = fk.length;
    int n1 = fk[0].length/2;
    double c1 = ceil(n1/2f);
    double c2 = ceil(n2/2f);
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0,ir=0,ii=1; i1<n1; ++i1,ir+=2,ii+=2) {
      double x1i = i1-c1;
      double x2i = i2-c2;
      double[] xc = new double[2];
      double xsi = x1i*x1i+x2i*x2i;
      double xsr = sqrt(xsi);
      if(xsi==0f) {
        xc[0] = 0f;
        xc[1] = 0f;
      } else {
        xc[0] = x2i/xsr;
        xc[1] = x1i/xsr;
      }
      xsi /= sig;
      double bx = _beta*xsi;
      double ax = _alpha*xsi;
      if (bx>700||ax>700){continue;}
      double exi = exp(-ax);
      double bsi = bessi(bx,s);
      double bsm = bessi(bx,s-1);
      double bsp = bessi(bx,s+1);
      double sci = _lambda*_lambda*xsi;
      double ssi = exi*(2f*(sci-1)*bsi+sci*(bsp+bsm));
      xc = cpow(xc,(s*2));
      fk[i2][ir] = (float)(ssi*xc[0]);
      fk[i2][ii] = (float)(ssi*xc[1]);
    }}
  }



  private float[][][][] filterBasisX(int n1, int n2, double sig) {
    double c1 = ceil(n1/2f);
    double c2 = ceil(n2/2f);
    float[][][][] bs = new float[2][_nd][n2][n1];
    for (int id=0; id<_nd; ++id) {
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0; i1<n1; ++i1) {
      double x1i = i1-c1;
      double x2i = i2-c2;
      double[] xc = new double[2];
      double xsi = x1i*x1i+x2i*x2i;
      if(xsi==0f) {
        xc[0] = 1f;
        xc[1] = 0f;
      } else {
        xc[0] = x2i/sqrt(xsi);
        xc[1] = x1i/sqrt(xsi);
      }
      xsi /= sig;
      double bx = _beta*xsi;
      double ax = _alpha*xsi;
      if (bx>700||ax>700){continue;}
      double exi = exp(-ax);
      double bsi = bessi(bx,id);
      double ssi = exi*bsi;
      xc = cpow(xc,(id*2));
      bs[0][id][i2][i1] = (float)(ssi*xc[0]);
      bs[1][id][i2][i1] = (float)(ssi*xc[1]);
    }}}
    return bs;
  }

  public static final int ACC = 40;
  public static final double BIGNO = 1.0e10;
  public static final double BIGNI = 1.0e-10;

  public double bessi(double x, int n) {
    if (n == 0) return Bessel.i0(x);
    if (n == 1) return Bessel.i1(x);
    int j;
    double bi, bim, bip, tox, ans;

    if (x == 0.0)
      return 0.0;
    else {
      tox = 2.0 / Math.abs(x);
      bip = ans = 0.0;
      bi = 1.0;
      for (j = 2 * (n + (int)Math.sqrt(ACC * n)); j > 0; j--) {
        bim = bip + j * tox * bi;
        bip = bi;
        bi = bim;
        if (Math.abs(bi) > BIGNO) {
          ans *= BIGNI;
          bi *= BIGNI;
          bip *= BIGNI;
        }
        if (j == n) ans = bip;
      }
      ans *= Bessel.i0(x) / bi;
      return (x < 0.0 && ((n & 1) != 0)) ? -ans : ans;
    }
  }

  private int _nd = 3;
  private float _lambda = 2;
  private float _beta, _alpha;
}
