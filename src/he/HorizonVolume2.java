package he;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Testing code for seismic carbonate interpretation
 * @author Xinming Wu
 * @version 2016.08.05
 */

public class HorizonVolume2 {
 
  // scale the curvature term 
  public void setWeights(float w){
    _weight = w;
  }
  
  public void setSmoothings(float sigma){
    _sigma = sigma;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }
  
  public void setExternalIterations(int exniter){
    _exniter = exniter;
  }

    // Interpolate an initial surface passing through control points
  public float[] curveInitialization
    (int n2, float lmt, float[] k1, float[] k2) 
  {
    if (k1.length==1) {
      float[] y = zerofloat(n2);
      add(y,k1[0],y);
      return y; 
    } else {
      CubicInterpolator ci = new CubicInterpolator(k2,k1);
      float[] x = new float[n2];
      float[] y = new float[n2];
      for (int i2=0; i2<n2; ++i2)
        x[i2] = i2;
      ci.interpolate(x,y);
      return y;
    }
  }


  // Interpolate an initial surface passing through control points
  public float[][] applyForInitial(int n2, float lmt, 
    float[] k1, int k2) 
  {
    int ns = k1.length;
    float[][] hv = new float[ns][n2];
    for (int is=0; is<ns; ++is) {
      hv[is] = fillfloat(k1[is],n2);
    }
    return hv;
  }

  
  public void applyForInitial(int k2,
    float lmt, float[] hz, float[] k1, float[][] hv)
  {
    int ns = hv.length;
    int n2 = hv[0].length;
    hv[0] = copy(hz);
    k1[0] = round(hz[k2]);
    for (int is=1; is<ns; ++is) {
      k1[is] = k1[is-1]+1f;
    for (int i2=0; i2<n2; ++i2) {
      float x1 = hv[is-1][i2]+1f;
      if(x1>lmt) {x1=lmt;}
      hv[is][i2]=x1;
    }}
  }

  public float[][] applyForHorizonVolume(
    float[][] k1, float[][] k2, float[][] ep, float[][] p) 
  {
    int ns = k1.length;
    int n2 = ep.length;
    int n1 = ep[0].length;
    float[][] hv = new float[ns][];
    //initial horizon volume with horizontal surfaces
    for (int is=0; is<ns; ++is) {
      hv[is] = curveInitialization(n2,n1-1,k1[is],k2[is]); 
    }
    float lmt = (float)n1-1.f;
    float[][] b  = new float[ns][n2]; 
    float[][] pi = new float[ns][n2]; 
    float[][] wi = new float[ns][n2]; 
    float[][][] ks = new float[][][]{k1,k2};
    float[][] hvt = copy(hv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat2 vb    = new VecArrayFloat2(b);
      VecArrayFloat2 vhv = new VecArrayFloat2(hv);
      updateSlopesAndWeights(p,ep,hv,pi,wi);
      A2 a2 = new A2(wi,_weight);
      M2 m2 = new M2(_sigma,wi,ks[1]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi,pi,b);
      cs.solve(a2,m2,vb,vhv);
      hv = vhv.getArray();
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(n2*ns); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
    return hv;
  }


  public float[][] applyForHorizonVolume(
    float[] k1, float[] k2, float[][] ep, float[][] p) 
  {
    int ns = k1.length;
    int n2 = ep.length;
    int n1 = ep[0].length;
    float[][] hv = new float[ns][];
    //initial horizon volume with horizontal surfaces
    for (int is=0; is<ns; ++is)
      hv[is] = fillfloat(k1[is],n2); 
    float lmt = (float)n1-1.f;
    float[][] b  = new float[ns][n2]; 
    float[][] pi = new float[ns][n2]; 
    float[][] wi = new float[ns][n2]; 
    //checkControlPoints(k2, hv); 
    //float[][][] ks = updateConstraints(k1,k2,p,ep,hv);
    int np = k1.length;
    float[][][] ks = new float[2][np][1];
    for (int ip=0; ip<np; ++ip) {
      ks[0][ip][0] = k1[ip];
      ks[1][ip][0] = k2[ip];
    }
    /*
    ks[0][0] = new float[2];
    ks[1][0] = new float[2];
    ks[0][np-1] = new float[2];
    ks[1][np-1] = new float[2];
    ks[0][0][0] = k1[0];
    ks[1][0][0] = k2[0];
    ks[0][0][1] = 40;
    ks[1][0][1] = 120;
    hv[0][120] = 40;
    ks[0][np-1][0] = k1[np-1];
    ks[1][np-1][0] = k2[np-1];
    ks[0][np-1][1] = 81;
    ks[1][np-1][1] = 50;
    hv[np-1][50] = 81;
    */

    float[][] hvt = copy(hv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat2 vb    = new VecArrayFloat2(b);
      VecArrayFloat2 vhv = new VecArrayFloat2(hv);
      updateSlopesAndWeights(p,ep,hv,pi,wi);
      A2 a2 = new A2(wi,_weight);
      M2 m2 = new M2(_sigma,wi,ks[1]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi,pi,b);
      cs.solve(a2,m2,vb,vhv);
      hv = vhv.getArray();
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(n2*ns); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
    return hv;
  }

  
  public float[][] applyForHorizonVolume(
    float[] k1, int k2, float[][] ep, float[][] p) 
  {
    int ns = k1.length;
    int n2 = ep.length;
    int n1 = ep[0].length;
    float[][] hv = new float[ns][];
    //initial horizon volume with horizontal surfaces
    for (int is=0; is<ns; ++is)
      hv[is] = fillfloat(k1[is],n2); 
    float lmt = (float)n1-1.f;
    float[][] b  = new float[ns][n2]; 
    float[][] pi = new float[ns][n2]; 
    float[][] wi = new float[ns][n2]; 
    checkControlPoints(k2, hv); 
    float[][][] ks = updateConstraints(k1,k2,p,ep,hv);
    float[][] hvt = copy(hv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat2 vb    = new VecArrayFloat2(b);
      VecArrayFloat2 vhv = new VecArrayFloat2(hv);
      updateSlopesAndWeights(p,ep,hv,pi,wi);
      A2 a2 = new A2(wi,_weight);
      M2 m2 = new M2(_sigma,wi,ks[1]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi,pi,b);
      cs.solve(a2,m2,vb,vhv);
      hv = vhv.getArray();
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(n2*ns); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
    return hv;
  }

  public float[][] rgtFromHorizonVolume(int n1, int dv, float[][] hx) {
    final int m1 = hx.length;
    final int n2 = hx[0].length;
    final float[][] ux = new float[n2][n1];
    System.out.println("n2="+n2);
    Parallel.loop(n2, new Parallel.LoopInt() {
    public void compute(int i2) {
      float[] x1s = new float[m1];
      float[] u1s = new float[m1];
      x1s[0] = hx[0][i2];
      for (int k1=1; k1<m1; ++k1) {
        u1s[k1] = k1*dv;
        x1s[k1] = hx[k1][i2];
        if(x1s[k1]<=x1s[k1-1]) {
          x1s[k1] = x1s[k1-1]+0.001f;
        }
      }
      CubicInterpolator ci  = new CubicInterpolator(x1s,u1s);
      int b1 = round(hx[0][i2]); b1 = max(b1,0);
      int e1 = round(hx[m1-1][i2]); e1 = min(e1,n1);
      for (int k1=b1; k1<e1; k1++) {
        ux[i2][k1] = ci.interpolate(k1);
      }
    }});
    return ux;
  }

  // Updates surfaces using the seismic normal vectors and control points.
  public void horizonUpdateFromSlopes(
     float[][] ep, float[][] p ,float[][][] q,
     float[] k1, int k2, float[] hz, float[][] hv)
  {	
    int ns = hv.length;
    int n2 = p.length; 
    int n1 = p[0].length; 
    float lmt = (float)n1-1.f;
    float[][] b  = new float[ns][n2]; 
    float[][] pi = new float[ns][n2]; 
    float[][] wi = new float[ns][n2]; 
    checkControlPoints(k2,hv); 
    float[][][] ks = updateConstraints(k1,k2,p,ep,hv);
    float[][] hvt = copy(hv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat2 vb    = new VecArrayFloat2(b);
      VecArrayFloat2 vhv = new VecArrayFloat2(hv);
      updateSlopesAndWeights(p,ep,hv,pi,wi);
      A2 a2 = new A2(wi,_weight);
      M2 m2 = new M2(_sigma,wi,ks[1]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi,pi,b);
      cs.solve(a2,m2,vb,vhv);
      hv = vhv.getArray();
      hv[0] = copy(hz);
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(n2*ns); 
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
    final float[][] p, final float[][] ep,
    final float[][] hv, final float[][] pi, 
    final float[][] wi)
  {
    final int n2 = p.length;
    final int n1 = p[0].length;
    final int ns = hv.length;
    final SincInterpolator psi = new SincInterpolator();
    final SincInterpolator wsi = new SincInterpolator();
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      float[] hs = hv[is];
      float[] ws = wi[is];
      float[] ps = pi[is];
      for (int i2=0; i2<n2; i2++){
        double x2i = (double)i2;
        double x1i = (double)hs[i2];
	      float wi = wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,ep,x1i,x2i);
        ps[i2]   = psi.interpolate(n1,1.0,0.0,n2,1.0,0.0, p,x1i,x2i);
        ws[i2] = (wi>0.0005f)?wi:0.0005f;
      }
    }});
  }

  private void horizonCorrection(float lmt, float[][] hv) {
    int ns = hv.length;
    int n2 = hv[0].length;
    for(int is=1; is<ns; ++is) {
    for (int i2=0; i2<n2; ++i2) {
      if (hv[is][i2]<0.f) hv[is][i2]=0.f;
      //if (hv[is][i2]<hv[is-1][i2]) hv[is][i2]=hv[is-1][i2];
      if (hv[is][i2]>lmt) hv[is][i2]=lmt;
    }}
  }

  private static class A2 implements CgSolver.A{
    A2(float[][] wp, float w){
      _w  = w;
      _wp = wp;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      int ns = y.length;
      int n2 = y[0].length;
      float[][] yy = new float[ns][n2];
      float[][] yt = new float[ns][n2];
      VecArrayFloat2 v2yy = new VecArrayFloat2(yy);
      v2y.zero();
      v2yy.zero();
      applyLhs(_wp,x,y);
      if (_w>0.0f) {
        applyLhs(_wp,x,yt);
        applyLhs(_wp,yt,yy);
        v2y.add(1.f,v2yy,_w);
      }
    }
    private float[][] _wp;
    private float _w;
  }

   // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma, float[][] wp, float[][] k2) {
      _sigma = sigma;
      _wp = wp;
      if (k2!=null) {
        _k2 = copy(k2);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k2,y);
      smooth1(_sigma,_wp,y);
      constrain(_k2,y);
    }
    private float _sigma;
    private float[][] _wp;
    private float[][] _k2;
  }


  private float[][][] updateConstraints( float[] k1, int k2,
    float[][] p, float[][] w, float[][] hv){
    int w2 = 10; 
    int n2 = p.length;
    int n1 = p[0].length;
    int nc = k1.length;
    int[] cot = new int[nc];
    float[][] cMark = new float[nc][n2];
    float[][] cConf = fillfloat(-1.0f,n2,nc);
    float wh = max(w)*0.8f;
    int c2 = k2;
    int i2b = c2-w2; if(i2b<0) i2b=0;
    int i2e = c2+w2; if(i2e>=n2) i2e=n2-1;
    int n2s = i2e-i2b+1;
    float[][] ws = copy(n1,n2s,0,i2b,w);
    float[][] ps = copy(n1,n2s,0,i2b,p);
    float[][] hvs = new float[nc][n2s];
    for (int ic=0; ic<nc; ++ic) {
      float c1 = k1[ic];
      hvs[ic]=fillfloat(c1,n2s);
    }
    subsetUpdate(ws,ps,k2-i2b,hvs);
    for (int ic=0; ic<nc; ++ic) {
      System.out.println("ic="+ic);
      for(int i2s=2; i2s<n2s-2; ++i2s) {
        int i2 = i2s+i2b;
        int d2 = i2-c2;
        float dsi = d2*d2;
        float cfi = cConf[ic][i2];
        float k1i = hvs[ic][i2s];
        if(k1i<0) {continue;}
        float wpi = w[i2][round(k1i)];
        if(dsi==0.0f) {
          cot[ic]++;
          cConf[ic][i2] = dsi;
          cMark[ic][i2] = k1i;
          continue;
        }
        if(wpi<wh && dsi>8.0f){continue;}
        if(cfi==-1.f) {
          cot[ic]++;
          cConf[ic][i2] = dsi;
          cMark[ic][i2] = k1i;
          continue;
        }
        if(dsi<cfi && abs(cMark[ic][i2]-k1i)<1.0f) {
          cConf[ic][i2] = dsi;
          cMark[ic][i2] = k1i;
        }
      }
    }
    float[][][] ks = new float[2][nc][];
    for (int ic=0; ic<nc; ++ic) {
      ks[0][ic] = new float[cot[ic]];
      ks[1][ic] = new float[cot[ic]];
      int k = 0;
      for (int i2=0; i2<n2; ++i2) {
        float cfi = cConf[ic][i2];
        if(cfi>=0.0f) {
          ks[1][ic][k] = i2;
          ks[0][ic][k] = cMark[ic][i2];
          hv[ic][i2] = cMark[ic][i2];
          k++;
        }
      }
    }
    return ks;
  }

  // Updates the surface using the seismic normal vectors and control points.
  public void subsetUpdate
    (float[][] w, float[][] p,
     int k2, float[][] hv)
  {	
    int n2 = p.length; 
    int n1 = p[0].length; 
    int ns = hv.length;
    float lmt = (float)n1-1.f;
    float[][] hvt = copy(hv);
    float[][] b  = new float[ns][n2]; 
    float[][] pi = new float[ns][n2]; 
    float[][] wi = new float[ns][n2]; 
    float[][] k2s = fillfloat(k2,ns,1);
    for (int n=1; n<=30; n++){
      VecArrayFloat2 vb  = new VecArrayFloat2(b);
      VecArrayFloat2 vhv = new VecArrayFloat2(hv);
      updateSlopesAndWeights(p,w,hv,pi,wi);
      A2 a2 = new A2(wi,_weight);
      M2 m2 = new M2(_sigma,wi,k2s);
      CgSolver cs = new CgSolver(_small, _niter);
      vb.zero();
      makeRhs(wi,pi,b);
      cs.solve(a2,m2,vb,vhv);
      hv = vhv.getArray();
      horizonCorrection(lmt,hv);
      float ad = sum(abs(sub(hvt,hv)))/(ns*n2); 
      if (ad<0.02f) break;
      hvt = copy(hv);
    }
  }


  private static void checkControlPoints(
    int k2, float[][] f) 
  {
    int ns = f.length;
    for (int is=0; is<ns; ++is) {
      System.out.println(" i2="+k2+" f1="+f[is][k2]);
    }
  }

  private static void constrain(float[][] k2, float[][] x) {
    int ns = k2.length;
    for (int is=0; is<ns; ++is) {
      int np = k2[is].length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[is][ip]; 
        x[is][i2] = 0.f;
      }
    }
  }

  private static void smooth1(
    final float sigma, final float[][] s, final float[][] x) 
  {
    final int ns = s.length;
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      smooth1(sigma,s[is],x[is]);
    }});
  }

  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[] s, float[] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n2 = x.length;
    float[] xt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(c,s,x,xt);
    copy(xt,x);
  }

  private static void applyLhs(
    final float[][] wp, final float[][] x, final float[][] y) 
  {
    zero(y);
    final int ns = x.length;
    final int n2 = x[0].length;
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      float[] xs = x[is];
      float[] ys = y[is];
      float[] ws = wp[is];
      for (int i2=1; i2<n2; ++i2) {
        float wi=ws[i2];
        float xa=0.0f;
        xa  = xs[i2  ];
        xa -= xs[i2-1];
        xa *= wi*wi;
        ys[i2-1] -= xa;
        ys[i2  ]  = xa;
      }
    }});
  }
  
  //private static void makeRhs
  private static void makeRhs(
    final float[][] wp, final float[][] p, final float[][] y) 
  {
    zero(y);
    final int ns = y.length;
    final int n2 = y[0].length;
    Parallel.loop(ns, new Parallel.LoopInt() {
    public void compute(int is) {
      float[] ys = y[is];
      float[] ps = p[is];

      float[] ws = wp[is];
      for (int i2=1; i2<n2; ++i2) {
        float wi = ws[i2];
        float pi = wi*wi*ps[i2];
        ys[i2  ] += pi;
        ys[i2-1] -= pi;
      }
    }});
  }
 
  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma = 6.0f; // precon smoothing extent for 1st dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 100; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
}
