package sso;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Testing code for seismic carbonate interpretation
 * @author Xinming Wu
 * @version 2016.05.12
 */

public class HorizonVolume {
 
  // scale the curvature term 
  public void setWeights(float w){
    _weight = w;
  }
  
  public void setSmoothings(float sigma1, float sigma2){
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }
  
  public void setExternalIterations(int exniter){
    _exniter = exniter;
  }


  // Interpolate an initial surface passing through control points
  public float[][][] applyForInitial(int n2, int n3, float lmt, 
    float[] k1, int k2, int k3) 
  {
    int ns = k1.length;
    float[][][] hv = new float[ns][n3][n2];
    for (int is=0; is<ns; ++is) {
      hv[is] = fillfloat(k1[is],n2,n3);
    }
    return hv;
  }

  
  public void applyForInitial(int k2, int k3, 
    float lmt, float[][] hz, float[] k1, float[][][] hv)
  {
    int ns = hv.length;
    int n3 = hv[0].length;
    int n2 = hv[0][0].length;
    hv[0] = copy(hz);
    k1[0] = round(hz[k3][k2]);
    for (int is=1; is<ns; ++is) {
      k1[is] = k1[is-1]+1f;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float x1 = hv[is-1][i3][i2]+1f;
      if(x1>lmt) {x1=lmt;}
      hv[is][i3][i2]=x1;
    }}}
  }

  public float[][][] applyForHorizonVolume(
    float[] k1, float[] k2, float[] k3, 
    float[][][] ep, float[][][] p, float[][][] q) 
  {
    int ns = k1.length;
    int n3 = ep.length;
    int n2 = ep[0].length;
    int n1 = ep[0][0].length;
    float[][][] hv = new float[ns][][];
    //initial horizon volume with horizontal surfaces
    for (int is=0; is<ns; ++is)
      hv[is] = fillfloat(k1[is],n2,n3); 
    float lmt = (float)n1-1.f;
    float[][][] b  = new float[ns][n3][n2]; 
    float[][][] pi = new float[ns][n3][n2]; 
    float[][][] qi = new float[ns][n3][n2]; 
    float[][][] wi = new float[ns][n3][n2]; 
    //checkControlPoints(k2, k3, hv); 
    float[][] k1s = new float[][]{k1};
    float[][] k2s = new float[][]{k2};
    float[][] k3s = new float[][]{k3};
    float[][][] ks = new float[][][]{k1s,k2s,k3s};
    //float[][][] ks = updateConstraints(k1,k2,k3,p,q,ep,hv);
    float[][][] hvt = copy(hv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat3 vb    = new VecArrayFloat3(b);
      VecArrayFloat3 vhv = new VecArrayFloat3(hv);
      updateSlopesAndWeights(p,q,ep,hv,pi,qi,wi);
      A3 a3 = new A3(wi,_weight);
      M3 m3 = new M3(_sigma1,_sigma2,wi,ks[1],ks[2]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi,pi,qi,b);
      cs.solve(a3,m3,vb,vhv);
      hv = vhv.getArray();
      horizonCorrectionX(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(n3*n2*ns); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
    return hv;
  }

  
  public float[][][] applyForHorizonVolume(
    float[] k1, int k2, int k3, 
    float[][][] ep, float[][][] p, float[][][] q) 
  {
    int ns = k1.length;
    int n3 = ep.length;
    int n2 = ep[0].length;
    int n1 = ep[0][0].length;
    float[][][] hv = new float[ns][][];
    //initial horizon volume with horizontal surfaces
    for (int is=0; is<ns; ++is)
      hv[is] = fillfloat(k1[is],n2,n3); 
    float lmt = (float)n1-1.f;
    float[][][] b  = new float[ns][n3][n2]; 
    float[][][] pi = new float[ns][n3][n2]; 
    float[][][] qi = new float[ns][n3][n2]; 
    float[][][] wi = new float[ns][n3][n2]; 
    checkControlPoints(k2, k3, hv); 
    float[][][] ks = updateConstraints(k1,k2,k3,p,q,ep,hv);
    float[][][] hvt = copy(hv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat3 vb    = new VecArrayFloat3(b);
      VecArrayFloat3 vhv = new VecArrayFloat3(hv);
      updateSlopesAndWeights(p,q,ep,hv,pi,qi,wi);
      A3 a3 = new A3(wi,_weight);
      M3 m3 = new M3(_sigma1,_sigma2,wi,ks[1],ks[2]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi,pi,qi,b);
      cs.solve(a3,m3,vb,vhv);
      hv = vhv.getArray();
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(n3*n2*ns); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
    return hv;
  }

  // Updates surfaces using the seismic normal vectors and control points.
  public void horizonUpdateFromSlopes(
     float[][][] ep, float[][][] p ,float[][][] q,
     float[] k1, int k2, int k3, float[][] hz, float[][][] hv)
  {	
    int n3 = p.length; 
    int n2 = p[0].length; 
    int n1 = p[0][0].length; 
    int ns = hv.length;
    float lmt = (float)n1-1.f;
    float[][][] b  = new float[ns][n3][n2]; 
    float[][][] pi = new float[ns][n3][n2]; 
    float[][][] qi = new float[ns][n3][n2]; 
    float[][][] wi = new float[ns][n3][n2]; 
    checkControlPoints(k2, k3, hv); 
    float[][][] ks = updateConstraints(k1,k2,k3,p,q,ep,hv);
    float[][][] hvt = copy(hv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat3 vb    = new VecArrayFloat3(b);
      VecArrayFloat3 vhv = new VecArrayFloat3(hv);
      updateSlopesAndWeights(p,q,ep,hv,pi,qi,wi);
      A3 a3 = new A3(wi,_weight);
      M3 m3 = new M3(_sigma1,_sigma2,wi,ks[1],ks[2]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi,pi,qi,b);
      cs.solve(a3,m3,vb,vhv);
      hv = vhv.getArray();
      hv[0] = copy(hz);
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(n3*n2*ns); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
  }


  public float[][] updateWeights(float[][] surf, float[][][] w) {
    int n3 = w.length;
    int n2 = w[0].length;
    int n1 = w[0][0].length;
    float[][] wi = new float[n3][n2];
    SincInterpolator wsi = new SincInterpolator();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        double x2i = (double)i2;
        double x3i = (double)i3;
        double x1i = (double)surf[i3][i2];
	      wi[i3][i2] = 
          wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,w,x1i,x2i,x3i);
      }
    }
    return wi;
  }

  private static void updateSlopesAndWeights (
    final float[][][] p, final float[][][] q, final float[][][] ep,
    final float[][][] hv, final float[][][] pi, 
    final float[][][] qi, final float[][][] wi)
  {
    final int n3 = p.length;
    final int n2 = p[0].length;
    final int n1 = p[0][0].length;
    final int ns = hv.length;
    final SincInterpolator psi = new SincInterpolator();
    final SincInterpolator qsi = new SincInterpolator();
    final SincInterpolator wsi = new SincInterpolator();
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      float[][] hs = hv[is];
      float[][] ws = wi[is];
      float[][] ps = pi[is];
      float[][] qs = qi[is];
      for (int i3=0; i3<n3; i3++){
        for (int i2=0; i2<n2; i2++){
          double x2i = (double)i2;
          double x3i = (double)i3;
          double x1i = (double)hs[i3][i2];
	        float wi = 
            wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,ep,x1i,x2i,x3i);
          ws[i3][i2] = (wi>0.0005f)?wi:0.0005f;
          ps[i3][i2] = 
            psi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, p,x1i,x2i,x3i);
	        qs[i3][i2] = 
            qsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, q,x1i,x2i,x3i);
        }
      }
    }});
  }

  private void horizonCorrectionX(float lmt, float[][][] hv) {
    int ns = hv.length;
    int n3 = hv[0].length;
    int n2 = hv[0][0].length;
    for(int is=1; is<ns; ++is) {
    for(int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      if (hv[is][i3][i2]<0.f) hv[is][i3][i2]=0.f;
      if (hv[is][i3][i2]>lmt) hv[is][i3][i2]=lmt;
    }}}
  }

  private void horizonCorrection(float lmt, float[][][] hv) {
    int ns = hv.length;
    int n3 = hv[0].length;
    int n2 = hv[0][0].length;
    for(int is=1; is<ns; ++is) {
    for(int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      if (hv[is][i3][i2]<0.f) hv[is][i3][i2]=0.f;
      if (hv[is][i3][i2]<hv[is-1][i3][i2]) hv[is][i3][i2]=hv[is-1][i3][i2];
      if (hv[is][i3][i2]>lmt) hv[is][i3][i2]=lmt;
    }}}
  }

  private static class A3 implements CgSolver.A{
    A3(float[][][] wp, float w){
      _w  = w;
      _wp = wp;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat3 v3x = (VecArrayFloat3) vx;
      VecArrayFloat3 v3y = (VecArrayFloat3) vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      int ns = y.length;
      int n2 = y[0].length;
      int n1 = y[0][0].length; 
      float[][][] yy = new float[ns][n2][n1];
      float[][][] yt = new float[ns][n2][n1];
      VecArrayFloat3 v3yy = new VecArrayFloat3(yy);
      v3y.zero();
      v3yy.zero();
      applyLhs(_wp,x,y);
      if (_w>0.0f) {
        applyLhs(_wp,x,yt);
        applyLhs(_wp,yt,yy);
        v3y.add(1.f,v3yy,_w);
      }
    }
    private float[][][] _wp;
    private float _w;
  }

   // Preconditioner; includes smoothers and (optional) constraints.
  private static class M3 implements CgSolver.A {
    M3(float sigma1, float sigma2, float[][][] wp, float[][] k2, float[][] k3) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
      if (k2!=null && k3!=null) {
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      copy(x,y);
      constrain(_k2,_k3,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      constrain(_k2,_k3,y);
    }
    private float _sigma1,_sigma2;
    private float[][][] _wp;
    private float[][] _k2,_k3;
  }

  private float[][][] updateConstraints( float[] k1, int k2, int k3, 
    float[][][] p, float[][][] q, float[][][] w, float[][][] hv){
    int w2 = 10; 
    int w3 = w2;
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    int nc = k1.length;
    int[] cot = new int[nc];
    float[][][] cMark = new float[nc][n3][n2];
    float[][][] cConf = fillfloat(-1.0f,n2,n3,nc);
    float wh = max(w)*0.8f;
    int c3 = k3;
    int c2 = k2;
    int i3b = c3-w3; if(i3b<0) i3b=0;
    int i2b = c2-w2; if(i2b<0) i2b=0;
    int i3e = c3+w3; if(i3e>=n3) i3e=n3-1;
    int i2e = c2+w2; if(i2e>=n2) i2e=n2-1;
    int n2s = i2e-i2b+1;
    int n3s = i3e-i3b+1;
    float[][][] ws = copy(n1,n2s,n3s,0,i2b,i3b,w);
    float[][][] ps = copy(n1,n2s,n3s,0,i2b,i3b,p);
    float[][][] qs = copy(n1,n2s,n3s,0,i2b,i3b,q);
    float[][][] hvs = new float[nc][n3s][n2s];
    for (int ic=0; ic<nc; ++ic) {
      float c1 = k1[ic];
      hvs[ic]=fillfloat(c1,n2s,n3s);
    }
    subsetUpdate(ws,ps,qs,k2-i2b,k3-i3b,hvs);
    for (int ic=0; ic<nc; ++ic) {
      System.out.println("ic="+ic);
      for(int i3s=2; i3s<n3s-2; ++i3s) {
        for(int i2s=2; i2s<n2s-2; ++i2s) {
          int i2 = i2s+i2b;
          int i3 = i3s+i3b;
          int d2 = i2-c2;
          int d3 = i3-c3;
          float dsi = d2*d2+d3*d3;
          float cfi = cConf[ic][i3][i2];
          float k1i = hvs[ic][i3s][i2s];
          if(k1i<0) {continue;}
          float wpi = w[i3][i2][round(k1i)];
          if(dsi==0.0f) {
            cot[ic]++;
            cConf[ic][i3][i2] = dsi;
            cMark[ic][i3][i2] = k1i;
            continue;
          }
          if(wpi<wh && dsi>8.0f){continue;}
          if(cfi==-1.f) {
            cot[ic]++;
            cConf[ic][i3][i2] = dsi;
            cMark[ic][i3][i2] = k1i;
            continue;
          }
          if(dsi<cfi && abs(cMark[ic][i3][i2]-k1i)<1.0f) {
            cConf[ic][i3][i2] = dsi;
            cMark[ic][i3][i2] = k1i;
          }
        }
      }
    }
    float[][][] ks = new float[3][nc][];
    for (int ic=0; ic<nc; ++ic) {
      ks[0][ic] = new float[cot[ic]];
      ks[1][ic] = new float[cot[ic]];
      ks[2][ic] = new float[cot[ic]];
      int k = 0;
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          float cfi = cConf[ic][i3][i2];
          if(cfi>=0.0f) {
            ks[2][ic][k] = i3;
            ks[1][ic][k] = i2;
            ks[0][ic][k] = cMark[ic][i3][i2];
            hv[ic][i3][i2] = cMark[ic][i3][i2];
            k++;
          }
        }
      }
    }
    return ks;
  }

  // Updates the surface using the seismic normal vectors and control points.
  public void subsetUpdate
    (float[][][] w, float[][][] p ,float[][][] q,
     int k2, int k3,float[][][] hv)
  {	
    int n3 = p.length; 
    int n2 = p[0].length; 
    int n1 = p[0][0].length; 
    int ns = hv.length;
    float lmt = (float)n1-1.f;
    float[][][] hvt = copy(hv);
    float[][][] b  = new float[ns][n3][n2]; 
    float[][][] pi = new float[ns][n3][n2]; 
    float[][][] qi = new float[ns][n3][n2]; 
    float[][][] wi = new float[ns][n3][n2]; 
    float[][] k2s = fillfloat(k2,ns,1);
    float[][] k3s = fillfloat(k3,ns,1);
    for (int n=1; n<=30; n++){
      VecArrayFloat3 vb  = new VecArrayFloat3(b);
      VecArrayFloat3 vhv = new VecArrayFloat3(hv);
      updateSlopesAndWeights(p,q,w,hv,pi,qi,wi);
      A3 a3 = new A3(wi,_weight);
      M3 m3 = new M3(_sigma1,_sigma2,wi,k2s,k3s);
      CgSolver cs = new CgSolver(_small, _niter);
      vb.zero();
      makeRhs(wi,pi,qi,b);
      cs.solve(a3,m3,vb,vhv);
      hv = vhv.getArray();
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(ns*n3*n2); 
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
  }


  private static void checkControlPoints(
    int k2, int k3, float[][][] f) 
  {
    int ns = f.length;
    for (int is=0; is<ns; ++is) {
      System.out.println(" i2="+k2+" i3="+k3+" f1="+f[is][k3][k2]);
    }
  }

  private static void constrain(float[][] k2, float[][] k3, float[][][] x) {
    int ns = k2.length;
    for (int is=0; is<ns; ++is) {
      int np = k2[is].length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[is][ip]; 
        int i3 = (int)k3[is][ip]; 
        x[is][i3][i2] = 0.f;
      }
    }
  }

  private static void smooth1(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int ns = s.length;
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      smooth1(sigma,s[is],x[is]);
    }});
  }
  private static void smooth2(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int ns = s.length;
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      smooth2(sigma,s[is],x[is]);
    }});
  }

  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    float[] st = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        xt[i1] = x[i2][i1];
        st[i1] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    float[] st = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        xt[i2] = x[i2][i1];
        st[i2] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }

  private static void applyLhs(
    final float[][][] wp, final float[][][] x, final float[][][] y) 
  {
    zero(y);
    final int ns = x.length;
    final int n2 = x[0].length;
    final int n1 = x[0][0].length;
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      float[][] xs = x[is];
      float[][] ys = y[is];
      float[][] ws = wp[is];
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float wpi = ws[i2][i1];
          if(wpi<0.05f) {wpi=0.05f;}
          float wps = wpi*wpi;
          float d11 = wps;
          float d22 = wps;
          float xa = 0.0f;
          float xb = 0.0f;
          xa += xs[i2  ][i1  ];
          xb -= xs[i2  ][i1-1];
          xb += xs[i2-1][i1  ];
          xa -= xs[i2-1][i1-1];
          float x1 = 0.5f*(xa+xb);
          float x2 = 0.5f*(xa-xb);
          float y1 = d11*x1;
          float y2 = d22*x2;
          float ya = 0.5f*(y1+y2);
          float yb = 0.5f*(y1-y2);
          ys[i2  ][i1  ] += ya;
          ys[i2  ][i1-1] -= yb;
          ys[i2-1][i1  ] += yb;
          ys[i2-1][i1-1] -= ya;
        }
      }
    }});
  }
  
  //private static void makeRhs
  private static void makeRhs(
    final float[][][] wp, final float[][][] p2, 
    final float[][][] p3, final float[][][] y) 
  {
    zero(y);
    final int ns = y.length;
    final int n2 = y[0].length;
    final int n1 = y[0][0].length;
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      float[][] ys  = y[is];
      float[][] p2s = p2[is];
      float[][] p3s = p3[is];
      float[][] wps = wp[is];
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float wpi = wps[i2][i1];
          if(wpi<0.05f) {wpi=0.05f;}
          float p2i = p2s[i2][i1];
          float p3i = p3s[i2][i1];
          float b11 = wpi;
          float b22 = wpi;
          float x1 = wpi*p2i;
          float x2 = wpi*p3i;
          float y1 = b11*x1;
          float y2 = b22*x2;
          float ya = 0.5f*(y1+y2);
          float yb = 0.5f*(y1-y2);
          ys[i2  ][i1  ] += ya;
          ys[i2  ][i1-1] -= yb;
          ys[i2-1][i1  ] += yb;
          ys[i2-1][i1-1] -= ya;
        }
      }
    }});
  }
 
  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 100; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
}
