/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package stv;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * 2D tensor voting using steerable filters.
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.11.03
 */
public class TensorVoting2 {

  public TensorVoting2 (float sigma) {
    _sigma = sigma;
  }

  public float[][][] applyVoting(float[][] ss, float[][] os) {
    int n2 = ss.length;
    int n1 = ss[0].length;
    int nfft1 = FftComplex.nfftSmall(n1);
    int nfft2 = FftComplex.nfftSmall(n2);
    float[][] c0 = czerofloat(nfft1,nfft2);
    float[][] c2 = czerofloat(nfft1,nfft2);
    float[][] c4 = czerofloat(nfft1,nfft2);
    float[][] c6 = czerofloat(nfft1,nfft2);
    float[][] w0 = czerofloat(nfft1,nfft2);
    float[][] w2 = czerofloat(nfft1,nfft2);
    float[][] w4 = czerofloat(nfft1,nfft2);
    float[][] w6 = czerofloat(nfft1,nfft2);
    float[][] w8 = czerofloat(nfft1,nfft2);

    imageFeatures(0,ss,os,c0);
    imageFeatures(2,ss,os,c2);
    imageFeatures(4,ss,os,c4);
    imageFeatures(6,ss,os,c6);
    float[][] cb = cconj(c2);
    filterKernals(0,_sigma,w0);
    filterKernals(2,_sigma,w2);
    filterKernals(4,_sigma,w4);
    filterKernals(6,_sigma,w6);
    filterKernals(8,_sigma,w8);

    FftComplex fft1 = new FftComplex(nfft1);
    FftComplex fft2 = new FftComplex(nfft2);
    fft1.complexToComplex1(-1,nfft2,c0,c0);
    fft2.complexToComplex2(-1,nfft1,c0,c0);
    fft1.complexToComplex1(-1,nfft2,c2,c2);
    fft2.complexToComplex2(-1,nfft1,c2,c2);
    fft1.complexToComplex1(-1,nfft2,c4,c4);
    fft2.complexToComplex2(-1,nfft1,c4,c4);
    fft1.complexToComplex1(-1,nfft2,c6,c6);
    fft2.complexToComplex2(-1,nfft1,c6,c6);
    fft1.complexToComplex1(-1,nfft2,cb,cb);
    fft2.complexToComplex2(-1,nfft1,cb,cb);

    fft1.complexToComplex1(-1,nfft2,w0,w0);
    fft2.complexToComplex2(-1,nfft1,w0,w0);
    fft1.complexToComplex1(-1,nfft2,w2,w2);
    fft2.complexToComplex2(-1,nfft1,w2,w2);
    fft1.complexToComplex1(-1,nfft2,w4,w4);
    fft2.complexToComplex2(-1,nfft1,w4,w4);
    fft1.complexToComplex1(-1,nfft2,w6,w6);
    fft2.complexToComplex2(-1,nfft1,w6,w6);
    fft1.complexToComplex1(-1,nfft2,w8,w8);
    fft2.complexToComplex2(-1,nfft1,w8,w8);

    float[][] u2 = czerofloat(nfft1,nfft2);
    for (int i2=0; i2<nfft2; ++i2){
    for (int i1=0,ir=0,ii=1; i1<nfft1; ++i1,ir+=2,ii+=2) {
      float c0r = c0[i2][ir];
      float c0i = c0[i2][ii];
      float c2r = c2[i2][ir];
      float c2i = c2[i2][ii];
      float c4r = c4[i2][ir];
      float c4i = c4[i2][ii];
      float c6r = c6[i2][ir];
      float c6i = c6[i2][ii];
      float cbr = cb[i2][ir];
      float cbi = cb[i2][ii];
      float w0r = w0[i2][ir];
      float w0i = w0[i2][ii];
      float w2r = w2[i2][ir];
      float w2i = w2[i2][ii];
      float w4r = w4[i2][ir];
      float w4i = w4[i2][ii];
      float w6r = w6[i2][ir];
      float w6i = w6[i2][ii];
      float w8r = w8[i2][ir];
      float w8i = w8[i2][ii];
      float wc0br = w0r*cbr-w0i*cbi;
      float wc0bi = w0r*cbi+w0i*cbr;
      float wc20r = w2r*c0r-w2i*c0i;
      float wc20i = w2r*c0i+w2i*c0r;
      float wc42r = w4r*c2r-w4i*c2i;
      float wc42i = w4r*c2i+w4i*c2r;
      float wc64r = w6r*c4r-w6i*c4i;
      float wc64i = w6r*c4i+w6i*c4r;
      float wc86r = w8r*c6r-w8i*c6i;
      float wc86i = w8r*c6i+w8i*c6r;
      //float wc00r = w0r*c0r-w0i*c0i;
      //float wc00i = w0r*c0i+w0r*c0i;
      //float wc22r = w2r*c2r-w2i*c2i;
      //float wc22i = w2r*c2i+w2r*c2i;
      //float wc44r = w4r*c4r-w4i*c4i;
      //float wc44i = w4r*c4i+w4r*c4i;
      //u0[i2][ir] = 6f*wc00r+8f*wc22r+2f*wc44r;
      //u0[i2][ii] = 6f*wc00i+8f*wc22i+2f*wc44i;
      u2[i2][ir] = wc0br+4f*wc20r+6f*wc42r+4*wc64r+wc86r;
      u2[i2][ii] = wc0bi+4f*wc20i+6f*wc42i+4*wc64i+wc86i;
    }}
    fft1.complexToComplex1(1,nfft2,u2,u2);
    fft2.complexToComplex2(1,nfft1,u2,u2);
    fft1.scale(nfft1,nfft2,u2);
    fft2.scale(nfft1,nfft2,u2);

    float[][] cu2 = cabs(u2);
    float[][] cu0 = zerofloat(nfft1,nfft2);
    float[][] au2 = zerofloat(nfft1,nfft2);
    float[][] au0 = zerofloat(nfft1,nfft2);
    int j1 = nfft1/2;
    int j2 = nfft2/2;
    copy(nfft1-j1,nfft2-j2,0,0,cu2,j1,j2,au2);
    copy(j1,j2,nfft1-j1,nfft2-j2,cu2,0,0,au2);
    copy(nfft1-j1,j2,0,nfft2-j2,cu2,j1,0,au2);
    copy(j1,nfft2-j2,nfft1-j1,0,cu2,0,j2,au2);
    for (int i2=0; i2<nfft2; ++i2){
    for (int i1=0,ir=0,ii=1; i1<nfft1; ++i1,ir+=2,ii+=2) {
      float u2r = u2[i2][ir];
      float u2i = u2[i2][ii];
      cu0[i2][i1] = 0.5f*atan2(u2i,u2r);
    }}
    copy(nfft1-j1,nfft2-j2,0,0,cu0,j1,j2,au0);
    copy(j1,j2,nfft1-j1,nfft2-j2,cu0,0,0,au0);
    copy(nfft1-j1,j2,0,nfft2-j2,cu0,j1,0,au0);
    copy(j1,nfft2-j2,nfft1-j1,0,cu0,0,j2,au0);
    return new float[][][]{au2,au0};
  }

  public float[][][] initialTensorField(
    float sigma, float sigma1, float sigma2, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float pi = (float)Math.PI;
    float[][] ss = new float[n2][n1];
    float[][] os = new float[n2][n1];
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] el = new float[n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
    rgf.apply10(fx,g1);
    rgf.apply01(fx,g2);
    lof.applyForNormalLinear(fx,u1,u2,el);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      ss[i2][i1] = sqrt(g1i*g1i+g2i*g2i);
      //ss[i2][i1] = abs(g1i*u1i+g2i*u2i);
      os[i2][i1] = 0.5f*pi+atan2(u1i,u2i);
    }}
    return new float[][][]{ss,os};
  }

  private void filterKernals(
    float m, float sigma, float[][] fk) {
    int n2 = fk.length;
    int n1 = fk[0].length/2;
    float c1 = ceil(n1/2f);
    float c2 = ceil(n2/2f);
    float sigmar = -0.5f/(sigma*sigma);
    float[] xc = new float[2];
    for (int i2=0; i2<n2; ++i2){
    for (int i1=0,ir=0,ii=1; i1<n1; ++i1,ir+=2,ii+=2) {
      float x1i = i1-c1;
      float x2i = i2-c2;
      if(x1i==0f && x2i==0f){
        fk[i2][ir] = 1f;
        fk[i2][ii] = 0f;
        continue;
      }
      float xsi = x1i*x1i+x2i*x2i;
      float exi = exp(xsi*sigmar);
      xc[0] = x2i/sqrt(xsi);
      xc[1] = x1i/sqrt(xsi);
      xc = cpow(xc,m);
      fk[i2][ir] = exi*xc[0];
      fk[i2][ii] = exi*xc[1];
    }}
  }

  private void imageFeatures(
    float m, float[][] ss, float[][] os, float[][] cs) {
    int n2 = ss.length;
    int n1 = ss[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0,ir=0,ii=1; i1<n1; ++i1,ir+=2,ii+=2) {
      float ssi = ss[i2][i1];
      float osi = os[i2][i1];
      cs[i2][ir] = ssi*cos(-m*osi);
      cs[i2][ii] = ssi*sin(-m*osi);
    }}
  }

  public float[][] findRidges(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] rs = new float[n2][n1];
    SincInterpolator si  =new SincInterpolator();
    LocalOrientFilter lof = new LocalOrientFilter(2.0,2.0);
    lof.applyForNormal(fx,u1,u2);
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      float x1p = i1+u1i;
      float x2p = i2+u2i;
      float x1m = i1-u1i;
      float x2m = i2-u2i;
      float fxi = fx[i2][i1];
      float fxm = si.interpolate(s1,s2,fx,x1m,x2m);
      float fxp = si.interpolate(s1,s2,fx,x1p,x2p);
      if(fxi>fxm&&fxi>fxp){
        rs[i2][i1]=fx[i2][i1];
        if(abs(u1i)>abs(u2i)) {
          int i1m = max(i1-1,0);
          int i1p = min(i1+1,n1-1);
          rs[i2][i1m]=fx[i2][i1m];
          rs[i2][i1p]=fx[i2][i1p];
        } else {
          int i2m = max(i2-1,0);
          int i2p = min(i2+1,n2-1);
          rs[i2m][i1]=fx[i2m][i1];
          rs[i2p][i1]=fx[i2p][i1];
        }
      }
    }}
    return rs;
  }
  private float _sigma=8.0f;


}
