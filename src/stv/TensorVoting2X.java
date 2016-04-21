/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package stv;

import util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

/**
 * 2D tensor voting using steerable filters.
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.11.03
 */
public class TensorVoting2X {

  public TensorVoting2X (int nd, float sigma) {
    _nd = nd;
    _sigma = sigma;
    setScales(nd);
  }

  private void setScales(int nd) {
    int nf = nd/2;
    double small = 1.0/pow(10.0,nd);
    _scs = filldouble(1.0,nd+1);
    double fn = factorial(1,nd,1)*small;
    for (int id=1; id<nf; ++id) {
      double fd = factorial(1,id,1)*factorial(1,nd-id,1)*small;
      _scs[id   ] = fn/fd;
      _scs[nd-id] = fn/fd;
    }
    double fd = factorial(1,nf,1)*factorial(1,nd-nf,1)*small;
    _scs[nf] = fn/fd;
    div(_scs,sum(_scs),_scs);
  }

  private float factorial(int nb, int ne, float fs) {
    for (int k=nb; k<=ne; ++k) 
      fs *= k;
    return fs;
  }


  public float[][][] applyVoting(final float[][] ss, final float[][] os) {
    final int ns = _nd+1;
    final int n2 = ss.length;
    final int n1 = ss[0].length;
    final int nfft1 = FftComplex.nfftSmall(n1);
    final int nfft2 = FftComplex.nfftSmall(n2);
    final float[][] u2 = czerofloat(nfft1,nfft2);
    final FftComplex fft1 = new FftComplex(nfft1);
    final FftComplex fft2 = new FftComplex(nfft2);
    loop(ns, new LoopInt() {
    public void compute(int is) {
      float sci = (float)_scs[is];
      float msi = _nd-(is+1)*2;
      float[][] ci = czerofloat(nfft1,nfft2);
      float[][] wi = czerofloat(nfft1,nfft2);
      imageFeatures( msi,ss,os,ci);
      filterKernals(-msi+2,_sigma,wi);
      fft1.complexToComplex1(-1,nfft2,ci,ci);
      fft2.complexToComplex2(-1,nfft1,ci,ci);
      fft1.complexToComplex1(-1,nfft2,wi,wi);
      fft2.complexToComplex2(-1,nfft1,wi,wi);
      for (int i2=0; i2<nfft2; ++i2){
      for (int i1=0,ir=0,ii=1; i1<nfft1; ++i1,ir+=2,ii+=2) {
        float cri = ci[i2][ir];
        float cii = ci[i2][ii];
        float wri = wi[i2][ir];
        float wii = wi[i2][ii];
        u2[i2][ir] += (cri*wri-cii*wii)*sci;
        u2[i2][ii] += (cri*wii+cii*wri)*sci;
      }}
    }});
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

  public float[][][] smoothEdge(
    float sigma1, float sigma2, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float pi = (float)Math.PI;
    float[][] ss = new float[n2][n1];
    float[][] os = new float[n2][n1];
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] gx = new float[n2][n1];
    LocalOrientFilterP lof = new LocalOrientFilterP(sigma1,sigma2);
    lof.applyForNormal(fx,u1,u2);
    EigenTensors2 ets = lof.applyForTensors(fx);
    ets.setEigenvalues(0.0001f,1.0f);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(ets,10,fx,gx);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
    rgf.apply10(gx,g1);
    rgf.apply01(gx,g2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float g1i = g1[i2][i1];
      float g2i = g2[i2][i1];
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      ss[i2][i1] = abs(g1i*u1i+g2i*u2i);
      os[i2][i1] = 0.5f*pi+atan2(u1i,u2i);
    }}
    ss = sub(ss,min(ss));
    ss = div(ss,max(ss));
    return new float[][][]{ss,gx};
  }

 
  public float[][] initialOrient(
    float sigma1, float sigma2, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float pi = (float)Math.PI;
    float[][] os = new float[n2][n1];
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] el = new float[n2][n1];
    LocalOrientFilterP lof = new LocalOrientFilterP(sigma1,sigma2);
    lof.applyForNormalLinear(fx,u1,u2,el);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i2][i1];
      float u2i = u2[i2][i1];
      os[i2][i1] = 0.5f*pi+atan2(u1i,u2i);
    }}
    return os;
  }

  public float[][][] initialTensorField(
    float gmin, float gmax,
    float sigma1, float sigma2, float[][] fx) {
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
    LocalOrientFilterP lof = new LocalOrientFilterP(sigma1,sigma2);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
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
    ss = sub(ss,min(ss));
    ss = div(ss,max(ss));
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float ssi = ss[i2][i1];
      if (ssi<gmin || ssi>gmax)
        ss[i2][i1] = 0.0f;
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
        fk[i2][ir] = 0f;
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
      cs[i2][ir] = ssi*cos(m*osi);
      cs[i2][ii] = ssi*sin(m*osi);
    }}
  }

  public float[][] findRidges(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] rs = new float[n2][n1];
    SincInterpolator si  =new SincInterpolator();
    LocalOrientFilterP lof = new LocalOrientFilterP(2.0,2.0);
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
        /*
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
        */
      }
    }}
    return rs;
  }

  private int _nd = 4;
  private float _sigma=8.0f;
  private double[] _scs;


}
