package mhe;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import vec.*;
import util.*;

/**
 * Extract a single seismic horizon curve with control points,
 * the constraints derived from control points are incorporated in a 
 * preconditioner in the conjugate gradient method used to solve the 
 * linear system for horizon extracting.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2018.08.08
 */

public class MultigridHorizon2{
 
  // scale the curvature term 
  public void setWeight(float w){
    _weight = w;
  }
  public void setSmoothing(float sigma1){
    _sigma1 = sigma1;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }
  
  public void setExternalIterations(int exniter){
    _exniter = exniter;
  }

  public void setConstrain(int d2, int dm, float wc) {
    _d2 = d2;
    _dm = dm;
    _wc = wc;
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

  // Updates the surface using the seismic normal vectors and control points.
  public float[] curveUpdateFromSlopes
    (float[][] ep, float[][] p ,
     float[] k1, float[] k2, float[] cv)
  {	
    int n2 = p.length; 
    int n1 = p[0].length; 
    float lmt = (float)n1-1.f;
    float[] cvt = copy(cv);
    float[] b   = new float[n2]; 
    float[] pi1 = new float[n2]; 
    float[] wi1 = new float[n2]; 
    checkControlPoints(k2,cv); 
    //float[][] ks = updateConstraints(k1,k2,p,ep,cv);
    float[][] ks = new float[][]{k1,k2};
    //int np = k1.length;
    //float[] cf = new float[np]; // force of constraints
    VecArrayFloat1 vb  = new VecArrayFloat1(b);
    VecArrayFloat1 vcv = new VecArrayFloat1(cv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      updateSlopesAndWeights(p,ep,cv,pi1,wi1);
      A1 a1 = new A1(_weight,wi1,null);
      M1 m1 = new M1(_sigma1,wi1,ks[1]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi1,pi1,null,null,b);
      cs.solve(a1,m1,vb,vcv);
      cv = vcv.getArray();
      cv = clip(0f,lmt,cv);
      float ad = sum(abs(sub(cvt,cv)))/n2; 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.01f) break;
      cvt = copy(cv);
    }
    return vcv.getArray();
    // show the constraint force for each control point,
    // the points with little constraint force are not important, 
    // and hence can be removed from the constraints.
    //checkConstraintForce(k2,k3,cf);
  }

  public float[] curveUpdateFromSlopes (
     int um, int wd, float[][] f,
     float[][] ep, float[][] p2, 
     float[] k1, float[] k2, float[] cv)
  {	
    int n2 = p2.length; 
    int n1 = p2[0].length; 
    float lmt = (float)n1-1.f;
    float[] cvt = copy(cv);
    float[] b   = new float[n2]; 
    float[] pi1 = new float[n2]; 
    float[] wi1 = new float[n2]; 
    VecArrayFloat1 vb  = new VecArrayFloat1(b);
    VecArrayFloat1 vcv = new VecArrayFloat1(cv);

    checkControlPoints(k2,cv); 
    //float[][] ks = updateConstraints(k1,k2,p,ep,cv);
    float[][] ks = new float[][]{k1,k2};
    //int np = k1.length;
    //float[] cf = new float[np]; // force of constraints
    float[][] cvs = new float[_exniter+1][n2];
    cvs[0] = copy(cv);
    int dm = 10;
    int dc = 10;
    GlobalCorrelationFinder gcf = new GlobalCorrelationFinder(-1,1);
    float[] ui1 = null;
    float[][] pc = null;
    for (int n=1; n<=_exniter; n++){
      dm = min(dm,60);
      System.out.println(" Iteration "+n+"......");
      if (n>1) {
         pc = gcf.getTraceIndexes(5,dm,dc,n2,k2,0.01f);
         int ns = pc.length;
         ui1 = new float[ns];
      }
      updateSlopesAndWeights(um,wd,f,p2,ep,pc,cv,pi1,wi1,ui1);
      A1 a1 = new A1(_weight,wi1,pc);
      M1 m1 = new M1(_sigma1,wi1,ks[1]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi1,pi1,pc,ui1,b);
      cs.solve(a1,m1,vb,vcv);
      cv = vcv.getArray();
      cv = clip(0f,lmt,cv);
      float ad = sum(abs(sub(cvt,cv)))/n2; 
      System.out.println(" Average adjustments per sample = "+ad);
      cvt = copy(cv);
      cvs[n-1] = copy(cv);
      dm += 20;
      dc += 20;
    }
    return vcv.getArray();
  
    // show the constraint force for each control point,
    // the points with little constraint force are not important, 
    // and hence can be removed from the constraints.
    //checkConstraintForce(k2,k3,cf);
  }


  /*
  public float[][] curveUpdateFromSlopes (
     int um, float[][] f,
     float[][] ep, float[][] p2, 
     float[][] pc, float[] k1, float[] k2, float[] cv)
  {	
    int ns = pc.length;
    int n2 = p2.length; 
    int n1 = p2[0].length; 
    float lmt = (float)n1-1.f;
    float[] cvt = copy(cv);
    float[] b   = new float[n2]; 
    float[] pi1 = new float[n2]; 
    float[] wi1 = new float[n2]; 
    float[] ui1 = new float[ns]; 
    VecArrayFloat1 vb  = new VecArrayFloat1(b);
    VecArrayFloat1 vcv = new VecArrayFloat1(cv);

    checkControlPoints(k2,cv); 
    //float[][] ks = updateConstraints(k1,k2,p,ep,cv);
    float[][] ks = new float[][]{k1,k2};
    //int np = k1.length;
    //float[] cf = new float[np]; // force of constraints
    float du = (float)um/(float)_exniter;
    int wd = um;
    float[][] cvs = new float[_exniter+1][n2];
    cvs[0] = copy(cv);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      if(n<=2) {
        updateSlopesAndWeights(p2,ep,cv,pi1,wi1);
        A1 a1 = new A1(_weight,wi1,null);
        M1 m1 = new M1(_sigma1,wi1,ks[1]);
        CgSolver cs = new CgSolver(_small,_niter);
        vb.zero();
        makeRhs(wi1,pi1,null,null,b);
        cs.solve(a1,m1,vb,vcv);
      } else {
        updateSlopesAndWeights(wd,f,p2,ep,pc,cv,pi1,wi1,ui1);
        A1 a1 = new A1(_weight,wi1,pc);
        M1 m1 = new M1(_sigma1,wi1,ks[1]);
        CgSolver cs = new CgSolver(_small,_niter);
        vb.zero();
        makeRhs(wi1,pi1,pc,ui1,b);
        cs.solve(a1,m1,vb,vcv);
        wd = round(um-(n-2)*du);
        wd = max(wd,5);
      }
      cv = vcv.getArray();
      cv = clip(0f,lmt,cv);
      float ad = sum(abs(sub(cvt,cv)))/n2; 
      System.out.println(" Average adjustments per sample = "+ad);
      cvt = copy(cv);
      cvs[n-1] = copy(cv);
      if (ad<0.01f&&n>2) break;
    }
    return cvs;
  
    // show the constraint force for each control point,
    // the points with little constraint force are not important, 
    // and hence can be removed from the constraints.
    //checkConstraintForce(k2,k3,cf);
  }
  */

  public float[] regrid1(float[] cv) {
    int n = cv.length;
    float[] sv = new float[n];
    FloatList xc = new FloatList();
    FloatList yc = new FloatList();
    for (int i=0; i<n; i+=16) {
      xc.add(i); yc.add(cv[i]);
    }
    CubicInterpolator ci = 
      new CubicInterpolator(CubicInterpolator.Method.SPLINE,xc.trim(),yc.trim());
    for (int i=0; i<n; i++) {
      sv[i] = ci.interpolate(i);
    }
    FloatList xs = new FloatList();
    FloatList ys = new FloatList();
    for (int i=8; i<n; i+=16) {
      xs.add(i);ys.add(sv[i]);
    }
    CubicInterpolator si = 
      new CubicInterpolator(CubicInterpolator.Method.SPLINE,xs.trim(),ys.trim());
    for (int i=0; i<n; i++) {
      sv[i] = si.interpolate1(i);
    }
    return sv;
  }

  public float[] regrid2(float[] cv) {
    int n = cv.length;
    float[] sv = new float[n];
    FloatList xc = new FloatList();
    FloatList yc = new FloatList();
    for (int i=0; i<n; i+=16) {
      xc.add(i); yc.add(cv[i]);
    }
    CubicInterpolator ci = 
      new CubicInterpolator(CubicInterpolator.Method.SPLINE,xc.trim(),yc.trim());
    for (int i=0; i<n; i++) {
      sv[i] = ci.interpolate1(i);
    }
    return sv;
  }


  private static void updateSlopesAndWeights (
    float[][]p, float[][]ep,
    float[] cv, float[] pi1, float[] wi1)
  {
    int n2 = p.length;
    int n1 = p[0].length;
    SincInterpolator psi = new SincInterpolator();
    SincInterpolator wsi = new SincInterpolator();
    for (int i2=0; i2<n2; i2++){
      double x2i = (double)i2;
      double x1i = (double)cv[i2];
	    float wi = wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,ep,x1i,x2i);
      wi1[i2] = (wi>0.0005f)?wi:0.0005f;
      pi1[i2] = psi.interpolate(n1,1.0,0.0,n2,1.0,0.0,p,x1i,x2i);
    }
  }

  private static void updateSlopesAndWeights (
    int um, int dc,
    float[][] fx, float[][]p, float[][]ep, float[][] pc,
    float[] cv, float[] pi1, float[] wi1, float[] ui1)
  {
    int n2 = p.length;
    int n1 = p[0].length;
    SincInterpolator psi = new SincInterpolator();
    SincInterpolator wsi = new SincInterpolator();
    for (int i2=0; i2<n2; i2++){
      double x2i = (double)i2;
      double x1i = (double)cv[i2];
	    float wi = wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,ep,x1i,x2i);
      wi1[i2] = (wi>0.0005f)?wi:0.0005f;
      pi1[i2] = psi.interpolate(n1,1.0,0.0,n2,1.0,0.0,p,x1i,x2i);
    }
    if(pc!=null) updateCorrelations(um,dc,fx,pc,cv,ui1);
  }


  private static void updateCorrelations(
    int um, int dc, float[][] fx, float[][] pc, float[] cv, float[] ui1) {
    int n1 = fx[0].length;
    int ns = ui1.length;
    final SincInterpolator usi = new SincInterpolator();
    final DynamicWarping dw = new DynamicWarping(-um,um);
    dw.setStrainMax(0.25);
    final Random rd = new Random();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      int p2 = (int)pc[is][0];
      int m2 = (int)pc[is][1];
      int cp = round(cv[p2]);
      int cm = round(cv[m2]);
      int nc = dc*2+1;
      float[] fp = new float[nc];
      float[] fm = new float[nc];
      for (int i1 = cp-dc, p1=0; i1<=cp+dc; i1++, p1++) {
        if(i1<0||i1>=n1) {
          int k1 = rd.nextInt(n1);
          fp[p1] = fx[p2][k1];
        } else {
          fp[p1] = fx[p2][i1];
        }
      }
      for (int i1 = cm-dc, m1=0; i1<=cm+dc; i1++, m1++) {
        if(i1<0||i1>=n1) {
          int k1 = rd.nextInt(n1);
          fm[m1] = fx[m2][k1];
        } else {
          fm[m1] = fx[m2][i1];
        }
      }
      float[] ui = dw.findShifts(fp,fm);
      ui1[is] = usi.interpolate(nc,1.0,0.0,ui,cv[m2]-cm+dc)+cm-cp;
    }});
  }



  private static class A1 implements CgSolver.A {
    A1(float w1, float[] wp, float[][] pc){
      _w1 = w1;
      _wp = wp;
      _pc = pc;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat1 v1x = (VecArrayFloat1) vx;
      VecArrayFloat1 v1y = (VecArrayFloat1) vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      int n1 = y.length;
      float[] yy = new float[n1];
      float[] yt = new float[n1];
      VecArrayFloat1 v1yy = new VecArrayFloat1(yy);
      v1y.zero();
      v1yy.zero();
      applyLhs(_wp,x,y);
      if(_pc!=null) screenLhs(_pc,_wp,x,y);
      if (_w1>0.0f) {
        applyLhs(_wp,x,yt);
        applyLhs(_wp,yt,yy);
        v1y.add(1.f,v1yy,_w1);
      }
    }
    private float _w1;
    private float[] _wp;
    private float[][] _pc = null;
  }

   // Preconditioner; includes smoothers and (optional) constraints.
  private static class M1 implements CgSolver.A {
    M1(float sigma1, float[] wp, float[] k2) {
      _sig1 = sigma1;
      _wp = wp;
      if (k2!=null) {
        _k2 = copy(k2);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      copy(x,y);
      constrain(_k2,y);
      //smooth1(_sig1,_wp,y);
      smooth1(_sig1,y);
      constrain(_k2,y);
    }
    private float _sig1;
    private float[] _wp;
    private float[] _k2;
  }

  private static void checkControlPoints(float[] k2, float[] f) {
    if (k2!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        System.out.println(" i2="+i2+" f1="+f[i2]);
      }
    }
  }

  private static void constrain(float[] k2, float[] x) {
    if (k2!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        x[i2] = 0.0f;
      }
    }
  }

  private static void smooth1(float sigma, float[] x){
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
    ref.apply(x,x);
  }

  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[] s, float[] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x.length;
    float[] yt = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(c,s,x,yt);
    for (int i1=0; i1<n1; ++i1)
      x[i1] = yt[i1];
  }


  private static void applyLhs(float[] w, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=1; i1<n1; ++i1) {
      float xa=0.0f;
      float wi=w[i1];
      float ws=wi*wi; 
      xa  = x[i1  ];
      xa -= x[i1-1];
      xa *= ws;
      y[i1-1] -= xa;
      y[i1  ]  = xa;
    }
  }

  private static void screenLhs(
    float [][] ks, float[] w, float[] x, float[] y) {
    int ns = ks.length;
    for (int is=0; is<ns; ++is) {
      int p1 = (int)ks[is][0];
      int m1 = (int)ks[is][1];
      float scale = ks[is][2];
      float wi = (w[p1]+w[m1])*0.5f;
      float ws=wi*wi; 
      ws *= ws;
      float dx = 0.0f;
      dx += x[p1];
      dx -= x[m1];
      dx *= scale*ws;
      y[m1] -= dx;
      y[p1] += dx;
    }
  }


  //private static void makeRhs
  private static void makeRhs(
    float[] w, float[] p, float[][] ks, float[] us, float[] y) 
  {
    zero(y);
    int n1 = y.length;
    for (int i1=1; i1<n1; ++i1) {
      float wi = w[i1];
      float ws = wi*wi;
      float pi = ws*p[i1];
      y[i1  ] += pi;
      y[i1-1] -= pi;
    }
    if(ks!=null) screenRhs(w,ks,us,y);
  }

  private static void screenRhs(
    float[] w,
    float[][] ks, float[] us, float[] y) {
    int ns = ks.length;
    for (int is=0; is<ns; ++is) {
      int p1 = (int)ks[is][0];
      int m1 = (int)ks[is][1];
      float scale = ks[is][2];
      float wi = (w[p1]+w[m1])*0.5f;
      float ws=wi*wi; 
      ws *= ws;
      float ui = -us[is];
      ui *= scale*ws;
      y[m1] -= ui;
      y[p1] += ui;
    }
  }

 
  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
  private float _wc = 1f;
  private int _d2 = 2;
  private int _dm = 20;


  private class FloatList {
  public int n;
  public float[] a = new float[1024];
  public void add(float f) {
    if (n==a.length) {
      float[] t = new float[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = f;
  }
  public float[] trim() {
    if (n==0)
      return null;
    float[] t = new float[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
 }
}
