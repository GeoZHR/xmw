package mhe;

import java.awt.*;
import java.util.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Extract a single seismic horizon surface with control points
 * the constraints derived from control points are incorporated in a 
 * preconditioner in the conjugate gradient method used to solve the 
 * linear system for horizon extracting.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.11
 */
public class GlobalHorizon3 {
 
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
  public float[][] surfaceInitialization
    (int n2, int n3, float lmt, float[] k1, float[] k2, float[] k3) 
  {
    if (k1.length==1) {
      float[][] surf = zerofloat(n2,n3);
      add(surf,k1[0],surf);
      return surf; 
    } else {
      Sampling s2 = new Sampling(n2,1.0f,0.0f);
      Sampling s3 = new Sampling(n3,1.0f,0.0f);
      RadialInterpolator2.Biharmonic bs = new RadialInterpolator2.Biharmonic();
      RadialGridder2 rg = new RadialGridder2(bs,k1,k2,k3);    
      float[][] surf = rg.grid(s2,s3);
      surfaceCorrect(surf,lmt);
      checkControlPoints(k2, k3, surf); 
      return surf;
    }
  }

  // Interpolate an initial surface passing through control points
  public float[][] surfaceInitializationFast
    (int n2, int n3, float lmt, float[] k1, float[] k2, float[] k3) 
  {
    if (k1.length==1) {
      float[][] surf = zerofloat(n2,n3);
      add(surf,k1[0],surf);
      return surf; 
    } else {
      Sampling s2 = new Sampling(n2,1.0f,0.0f);
      Sampling s3 = new Sampling(n3,1.0f,0.0f);
      BlendedGridder2 bg = new BlendedGridder2(k1,k2,k3);
      float[][] surf = bg.grid(s2,s3);
      surfaceCorrect(surf,lmt);
      checkControlPoints(k1,k2,k3,surf); 
      return surf;
    }
  }

  // Updates the surface using the seismic normal vectors and control points.
  public float[][] surfaceUpdateFromSlopesAndCorrelations(
    int um, int wd, float[][][] fx,
    float[][][] ep, float[][][] p ,float[][][] q, 
    float[] k2, float[] k3, float[][] surf)
  {	
    int n3 = p.length; 
    int n2 = p[0].length; 
    int n1 = p[0][0].length; 
    float lmt = (float)n1-1.f;
    float[][] surft = copy(surf);
    float[][] b   = new float[n3][n2]; 
    float[][] pi1 = new float[n3][n2]; 
    float[][] qi1 = new float[n3][n2]; 
    float[][] wi1 = new float[n3][n2]; 
    checkControlPoints(k2, k3, surf); 
    VecArrayFloat2 vb    = new VecArrayFloat2(b);
    VecArrayFloat2 vsurf = new VecArrayFloat2(surf);
    GlobalCorrelationFinder gcf = new GlobalCorrelationFinder(-10,10);
    int dm = 10;
    int dc = 10;
    float[] ui1 = null;
    float[][] pc = null;
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      dm = min(dm,100);
      if(n>1) {
        pc = gcf.getTraceIndexes(10,10,dm,dc,k2,k3,n2,n3,0.01f);
        int ns = pc.length;
        ui1 = new float[ns];
      }
      updateSlopesAndWeights(um,wd,fx,p,q,ep,pc,surf,pi1,qi1,wi1,ui1);
      A2 a2 = new A2(_weight,wi1,pc);
      M2 m2 = new M2(_sigma1,_sigma2,wi1,k2,k3);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi1,pi1,qi1,pc,ui1,b);
      cs.solve(a2,m2,vb,vsurf);
      checkControlPoints(k2, k3, surf); 
      surf = vsurf.getArray();
      surfaceCorrect(surf,lmt);
      float ad = sum(abs(sub(surft,surf)))/(n3*n2); 
      System.out.println(" Average adjustments per sample = "+ad);
      surft = copy(surf);
      dm += 20;
      dc += 20;
    }
    return vsurf.getArray();
  }



  // Updates the surface using the seismic normal vectors and control points.
  public float[][] surfaceUpdateFromSlopesAndCorrelations(
    float[][][] ep, float[][][] p ,float[][][] q,
    float[] k2, float[] k3, float[][] surf)
  {	
    int n3 = p.length; 
    int n2 = p[0].length; 
    int n1 = p[0][0].length; 
    float lmt = (float)n1-1.f;
    float[][] surft = copy(surf);
    float[][] b   = new float[n3][n2]; 
    float[][] pi1 = new float[n3][n2]; 
    float[][] qi1 = new float[n3][n2]; 
    float[][] wi1 = new float[n3][n2]; 
    checkControlPoints(k2, k3, surf); 
    VecArrayFloat2 vb    = new VecArrayFloat2(b);
    VecArrayFloat2 vsurf = new VecArrayFloat2(surf);
    int niter = 100;
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      updateSlopesAndWeights(p,q,ep,surf,pi1,qi1,wi1);
      A2 a2 = new A2(_weight,wi1,null);
      M2 m2 = new M2(_sigma1,_sigma2,wi1,k2,k3);
      if(n>5) {niter=_niter;}
      CgSolver cs = new CgSolver(_small,niter);
      vb.zero();
      makeRhs(wi1,pi1,qi1,b);
      cs.solve(a2,m2,vb,vsurf);
      checkControlPoints(k2, k3, surf); 
      surf = vsurf.getArray();
      surfaceCorrect(surf,lmt);
      float ad = sum(abs(sub(surft,surf)))/(n3*n2); 
      System.out.println(" Average adjustments per sample = "+ad);
      //if (ad<0.02f) break;
      surft = copy(surf);
    }
    return vsurf.getArray();
    // show the constraint force for each control point,
    // the points with little constraint force are not important, 
    // and hence can be removed from the constraints.
    //checkConstraintForce(k2,k3,cf);
  }

  public float[][][] heightRgb(
    ColorMap mp, float[][] sf) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = sf[i3][i2];
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }

  public float[][][] amplitudeRgb(
    ColorMap mp, float[][][] fx, float[][] sf) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = si.interpolate(n1,1,0,fx[i3][i2],sf[i3][i2]);
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
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
    float[][][] p, float[][][] q, float[][][] ep,
    float[][] surf, float[][] pi1, float[][] qi1,float[][] wi1)
  {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    SincInterpolator psi = new SincInterpolator();
    SincInterpolator qsi = new SincInterpolator();
    SincInterpolator wsi = new SincInterpolator();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        double x2i = (double)i2;
        double x3i = (double)i3;
        double x1i = (double)surf[i3][i2];
	      float wi = 
          wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,ep,x1i,x2i,x3i);
        wi1[i3][i2] = (wi>0.0005f)?wi:0.0005f;
        pi1[i3][i2] = 
          psi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, p,x1i,x2i,x3i);
	      qi1[i3][i2] = 
          qsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, q,x1i,x2i,x3i);
      }
    }
  }

  private static void updateSlopesAndWeights (
    int um, int dc, float[][][] fx,
    float[][][] p, float[][][] q, float[][][] ep, float[][] pc,
    float[][] surf, float[][] pi1, float[][] qi1,float[][] wi1, float[] ui1)
  {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    SincInterpolator psi = new SincInterpolator();
    SincInterpolator qsi = new SincInterpolator();
    SincInterpolator wsi = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++){
        double x2i = (double)i2;
        double x3i = (double)i3;
        double x1i = (double)surf[i3][i2];
	      float wi = 
          wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,ep,x1i,x2i,x3i);
        wi1[i3][i2] = (wi>0.0005f)?wi:0.0005f;
        pi1[i3][i2] = 
          psi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, p,x1i,x2i,x3i);
	      qi1[i3][i2] = 
          qsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, q,x1i,x2i,x3i);
      }
    }});
    if(pc!=null) updateCorrelations(um,dc,fx,pc,surf,ui1);
  }


  private static void updateCorrelations(
    int um, int dc, float[][][] fx, float[][] pc, float[][] surf, float[] ui1) {
    final int ns = ui1.length;
    int n1 = fx[0][0].length;
    final SincInterpolator usi = new SincInterpolator();
    final DynamicWarping dw = new DynamicWarping(-um,um);
    dw.setStrainMax(0.25);
    dw.setErrorSmoothing(2);
    dw.setShiftSmoothing(4);
    final Random rd = new Random();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      int p2 = (int)pc[is][0];
      int p3 = (int)pc[is][1];
      int m2 = (int)pc[is][2];
      int m3 = (int)pc[is][3];
      int cp = round(surf[p3][p2]);
      int cm = round(surf[m3][m2]);
      int nc = dc*2+1;
      float[] fp = new float[nc];
      float[] fm = new float[nc];
      for (int i1 = cp-dc, p1=0; i1<=cp+dc; i1++, p1++) {
        if(i1<0||i1>=n1) {
          int k1 = rd.nextInt(n1);
          fp[p1] = fx[p3][p2][k1];
        } else {
          fp[p1] = fx[p3][p2][i1];
        }
      }
      for (int i1 = cm-dc, m1=0; i1<=cm+dc; i1++, m1++) {
        if(i1<0||i1>=n1) {
          int k1 = rd.nextInt(n1);
          fm[m1] = fx[m3][m2][k1];
        } else {
          fm[m1] = fx[m3][m2][i1];
        }
      }
      float[] ui = dw.findShifts(fp,fm);
      int dm = min(cm,abs(n1-cm));
      int dp = min(cp,abs(n1-cp));
      if(min(dm,dp)<2) pc[is][4] = 0f;
      ui1[is] = usi.interpolate(nc,1.0,0.0,ui,surf[m3][m2]-cm+dc)+cm-cp;
    }});
  }


  private static float correlate(float[] f, float[] g) {
    int n1 = f.length;
    float ff = 0f;
    float gg = 0f;
    float fg = 0f;
    float fa = 0f;
    float ga = 0f;
    int dc = 10;
    int c1 = (n1-1)/2;
    int b1 = c1-dc;
    int e1 = c1+dc;
    int nd = dc*2+1;
    for (int i1=b1; i1<=e1; ++i1) {
      fa += f[i1];
      ga += g[i1];
    }
    fa /= nd;
    ga /= nd;
    for (int i1=b1; i1<=e1; ++i1) {
      float fi = f[i1]-fa;
      float gi = g[i1]-ga;
      fg += fi*gi;
      ff += fi*fi;
      gg += gi*gi;
    }
    return (fg*fg)/(ff*gg);
  }



  private void surfaceCorrect(float[][] surf, float lmt) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    for(int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (surf[i2][i1]<0.f) surf[i2][i1]=0.f;
        if (surf[i2][i1]>lmt) surf[i2][i1]=lmt;
      }
    }
  }

  private static class A2 implements CgSolver.A{
    A2(float w1, float[][] wp, float[][] pc){
      _w1 = w1;
      _wp = wp;
      _pc = pc;
      testSpd();
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      int n1 = y[0].length; int n2 = y.length;
      float[][] yy = new float[n2][n1];
      VecArrayFloat2 v2yy = new VecArrayFloat2(yy);
      v2y.zero();
      v2yy.zero();
      applyLhs(_wp,x,y);
      if(_pc!=null) screenLhs(_pc,_wp,x,y);
      if (_w1>0.0f) {
        applyLhs(_wp,y,yy);
        v2y.add(1.f,v2yy,_w1);
      }
    }
    private float _w1;
    private float[][] _wp;
    private float[][] _pc;
    public void testSpd() {
    // symmetric: y'Ax = x'(A'y) = x'Ay
    // positive-semidefinite: x'Ax >= 0
      int n2 = _wp.length;
      int n1 = _wp[0].length;
      float[][] x = sub(randfloat(n1,n2),0.5f);
      float[][] y = sub(randfloat(n1,n2),0.5f);
      float[][] ax = zerofloat(n1,n2);
      float[][] ay = zerofloat(n1,n2);
      VecArrayFloat2 vx = new VecArrayFloat2(x);
      VecArrayFloat2 vy = new VecArrayFloat2(y);
      VecArrayFloat2 vax = new VecArrayFloat2(ax);
      VecArrayFloat2 vay = new VecArrayFloat2(ay);
      apply(vx,vax);
      apply(vy,vay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("A3: yax="+yax+" xay="+xay);
      System.out.println("A3: xax="+xax+" yay="+yay);
    }
  }

  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[][] wp, float[] k2, float[] k3) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
      if (k2!=null && k3!=null) {
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k2,_k3,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      constrain(_k2,_k3,y);
    }
    private float _sigma1,_sigma2;
    private float[][] _wp;
    private float[] _k2,_k3;
  }

  private static void checkControlPoints(float[] k2, float[] k3, float[][] f) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
      }
    }
  }

  private static void checkControlPoints(
    float[] k1, float[] k2, float[] k3, float[][] f) 
  {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
        f[i3][i2] = k1[ip];
      }
    }
  }


  private static void constrain(float[] k2, float[] k3, float[][] x) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2] = 0.f;
      }
    }
  }

  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    int n2 = x.length;
    int n1 = x[0].length;
    float c = 0.5f*sigma*sigma;
    final LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[] xt = zerofloat(n1);
      float[] yt = zerofloat(n1);
      float[] st = zerofloat(n1);
      for (int i1=0; i1<n1; ++i1) {
        xt[i1] = x[i2][i1];
        st[i1] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }});
  }

  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    final int n2 = x.length;
    final int n1 = x[0].length;
    final float c = 0.5f*sigma*sigma;
    final LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[] xt = zerofloat(n2);
      float[] yt = zerofloat(n2);
      float[] st = zerofloat(n2);
      for (int i2=0; i2<n2; ++i2) {
        xt[i2] = x[i2][i1];
        st[i2] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }});
  }

  /*
  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    final int n2 = x.length;
    Parallel.loop(1,n2,2,new Parallel.LoopInt() { // i2 = 1, 3, 5, ...
    public void compute(int i2) {
      applyLhsSlice2(i2,wp,x,y);
    }});
    Parallel.loop(2,n2,2,new Parallel.LoopInt() { // i2 = 2, 4, 6, ...
    public void compute(int i2) {
      applyLhsSlice2(i2,wp,x,y);
    }});
  }
  */

  private static void applyLhs(
    float[][] wp, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
    for (int i1=1; i1<n1; ++i1) {
      float wpi = (wp!=null)?wp[i2][i1]:1.000f;
      if(wpi<0.05f) {wpi=0.05f;}
      float ws = wpi*wpi*0.25f;
      float xa = 0.0f;
      float xb = 0.0f;
      xa += x[i2  ][i1  ];
      xb -= x[i2  ][i1-1];
      xb += x[i2-1][i1  ];
      xa -= x[i2-1][i1-1];
      float x1 = xa+xb;
      float x2 = xa-xb;
      float y1 = ws*x1;
      float y2 = ws*x2;
      float ya = y1+y2;
      float yb = y1-y2;
      y[i2  ][i1  ] += ya;
      y[i2  ][i1-1] -= yb;
      y[i2-1][i1  ] += yb;
      y[i2-1][i1-1] -= ya;
    }}
  }

  private static void screenLhs(
    float [][] ks, float[][] w, float[][] x, float[][] y) {
    int ns = ks.length;
    for (int is=0; is<ns; ++is) {
      int p2 = (int)ks[is][0];
      int p3 = (int)ks[is][1];
      int m2 = (int)ks[is][2];
      int m3 = (int)ks[is][3];
      float scale = ks[is][4];
      float wi = (w[p3][p2]+w[m3][m2])*0.5f;
      float ws = wi*wi; 
      float dx = 0.0f;
      dx += x[p3][p2];
      dx -= x[m3][m2];
      dx *= scale*ws;
      y[m3][m2] -= dx;
      y[p3][p2] += dx;
    }
  }

  //private static void makeRhs
  private static void makeRhs(
    float[][] wp, float[][] p2, float[][] p3, float[][] y) {
    zero(y);
    int n2 = y.length;
    int n1 = y[0].length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f) {wpi=0.05f;}
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float ws = wpi*wpi*0.5f;
        float y1 = ws*p2i;
        float y2 = ws*p3i;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private static void makeRhs(
    float[][] wp, float[][] p2, float[][] p3, 
    float[][] pc, float[] us, float[][] y) {
    zero(y);
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f) {wpi=0.05f;}
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float wps = wpi*wpi*0.5f;
        float y1 = wps*p2i;
        float y2 = wps*p3i;
        float ya = (y1+y2);
        float yb = (y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
    if(pc!=null) screenRhs(wp,pc,us,y);
  }


  private static void screenRhs(
    float[][] w, float[][] ks, float[] us, float[][] y) {
    int ns = ks.length;
    for (int is=0; is<ns; ++is) {
      int p2 = (int)ks[is][0];
      int p3 = (int)ks[is][1];
      int m2 = (int)ks[is][2];
      int m3 = (int)ks[is][3];
      float scale = ks[is][4];
      float wi = (w[p3][p2]+w[m3][m2])*0.5f;
      float ws = wi*wi; 
      float ui = -us[is];
      ui *= scale*ws;
      y[m3][m2] -= ui;
      y[p3][p2] += ui;
    }
  }


  private void checkConstraintForce(float[] k2, float[] k3, float[] cf) {
    int np = k2.length;
    div(cf,max(cf),cf);
    for (int ip=0; ip<np; ++ip) {
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      System.out.println(" i2="+i2+" i3="+i3+" ad="+cf[ip]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
}
