package he;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Extract a single seismic horizon curve with control points,
 * the constraints derived from control points are incorporated in a 
 * preconditioner in the conjugate gradient method used to solve the 
 * linear system for horizon extracting.
 * @author Xinming Wu and Dave Hale
 * @version 2015.10.01
 */

public class HorizonExtractor2{
 
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

  // find the peak or trough nearest to each control point
  public float[] refineConstraints(
    float[] k1, float[] k2, float[] k3, float[][][] u) 
  {
    int np = k1.length;
    int n1 = u[0][0].length;
    for (int ip=0; ip<np; ++ip) {
      float k1i = k1[ip];
      int i1 = (int)k1i;
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      int i1m = i1-1; 
      int i1p = i1+1; 
      if(i1m<0  ){i1m=0;   i1=1;   i1p=2;}
      if(i1p>=n1){i1m=n1-3;i1=n1-2;i1p=n1-1;}
      float um = u[i3][i2][i1m];
      float ui = u[i3][i2][i1 ];
      float up = u[i3][i2][i1p];
      float kp = parabolicPeak(k1i,um,ui,up);
      if(abs(kp-k1i)<2.0 && k1i*kp>=0.0f){k1[ip]=kp;}
    }
    return k1;
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
  public void curveUpdateFromSlopes
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
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat1 vb  = new VecArrayFloat1(b);
      VecArrayFloat1 vcv = new VecArrayFloat1(cv);
      updateSlopesAndWeights(p,ep,cv,pi1,wi1);
      A1 a1 = new A1(_weight,wi1);
      M1 m1 = new M1(_sigma1,wi1,ks[1]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(wi1,pi1,b);
      cs.solve(a1,m1,vb,vcv);
      cv = vcv.getArray();
      clip(0f,lmt,cv);
      float ad = sum(abs(sub(cvt,cv)))/n2; 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.01f) break;
      cvt = copy(cv);
    }
  
    // show the constraint force for each control point,
    // the points with little constraint force are not important, 
    // and hence can be removed from the constraints.
    //checkConstraintForce(k2,k3,cf);
  }

  public void surfaceRefine(float[][] surf, float[][][] u) {
    int n3 = u.length;
    int n2 = u[0].length;
    float[] k1 = new float[1];
    float[] k2 = new float[1];
    float[] k3 = new float[1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        k2[0] = (float)i2;
        k3[0] = (float)i3;
        k1[0] = surf[i3][i2];
        surf[i3][i2]=refineConstraints(k1,k2,k3,u)[0];
      }
    }
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

  private static class A1 implements CgSolver.A {
    A1(float w, float[] wp){
      _w  = w;
      _wp = wp;
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
      if (_w>0.0f) {
        applyLhs(_wp,x,yt);
        applyLhs(_wp,yt,yy);
        v1y.add(1.f,v1yy,_w);
      }
    }
    private float _w;
    private float[] _wp;
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

  private float[][] updateConstraints(float[] k1, float[] k2,
    float[][] p, float[][] w, float[] cv)
  {
    int cot = 0;
    int w2 = 10; 
    int n2 = p.length;
    int n1 = p[0].length;
    float[] cMark = new float[n2];
    float[] cConf = fillfloat(-1.0f,n2);
    int nc = k1.length;
    float wh = max(w)*0.8f;
    for (int ic=0; ic<nc; ++ic) {
      System.out.println("ic="+ic);
      float c1 = k1[ic];
      int c2 = (int)k2[ic];
      int i2b = c2-w2; if(i2b<0) i2b=0;
      int i2e = c2+w2; if(i2e>=n2) i2e=n2-1;
      int n2s = i2e-i2b+1;
      float[][] ws = copy(n1,n2s,0,i2b,w);
      float[][] ps = copy(n1,n2s,0,i2b,p);
      float[] k2s = new float[]{c2-i2b};
      float[] cvs = fillfloat(c1,n2s);
      subsetUpdate(ws,ps,k2s,cvs);
      for(int i2s=2; i2s<n2s-2; ++i2s) {
        int i2 = i2s+i2b;
        int d2 = i2-c2;
        float dsi = d2*d2;
        float cfi = cConf[i2];
        float k1i = cvs[i2s];
        float wpi = w[i2][round(k1i)];
        if(dsi==0.0f) {
          cot++;
          cConf[i2] = dsi;
          cMark[i2] = k1i;
          continue;
        }
        if(wpi<wh && dsi>8.0f){continue;}
        if(cfi==-1.f) {
          cot++;
          cConf[i2] = dsi;
          cMark[i2] = k1i;
          continue;
        }
        if(dsi<cfi && abs(cMark[i2]-k1i)<1.0f) {
          cConf[i2] = dsi;
          cMark[i2] = k1i;
        }
      }
    }
    int k = 0;
    float[][] ks = new float[2][cot];
    for (int i2=0; i2<n2; ++i2) {
      float cfi = cConf[i2];
      if(cfi>=0.0f) {
        ks[1][k] = i2;
        ks[0][k] = cMark[i2];
        cv[i2] = cMark[i2];
        k++;
      }
    }
    return ks;
  }

  // Updates the surface using the seismic normal vectors and control points.
  public void subsetUpdate(float[][]w, float[][] p,float[] k2, float[] cv)
  {	
    int n2 = p.length; 
    int n1 = p[0].length; 
    float lmt = (float)n1-1.f;
    float[] cvt = copy(cv);
    float[] b   = new float[n2]; 
    float[] pi1 = new float[n2]; 
    float[] wi1 = new float[n2]; 
    for (int n=1; n<=30; n++){
      VecArrayFloat1 vb  = new VecArrayFloat1(b);
      VecArrayFloat1 vcv = new VecArrayFloat1(cv);
      updateSlopesAndWeights(p,w,cv,pi1,wi1);
      A1 a1 = new A1(_weight,wi1);
      M1 m1 = new M1(_sigma1,wi1,k2);
      CgSolver cs = new CgSolver(_small, _niter);
      vb.zero();
      makeRhs(wi1,pi1,b);
      cs.solve(a1,m1,vb,vcv);
      cv = vcv.getArray();
      clip(0f,lmt,cv);
      float ad = sum(abs(sub(cvt,cv)))/n2; 
      if (ad<0.01f) break;
      cvt = copy(cv);
    }
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
    ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE);
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

  //private static void makeRhs
  private static void makeRhs(
    float[] w, float[] p, float[] y) 
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
  }
 
  // use 3 points to fit a parabolic curve and find its peak
  private float parabolicPeak(float z, float um, float ui, float up) {
    float a = um-up;
    float b = 2.0f*(um+up)-4.0f*ui;
    return (z+a/b);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
}
